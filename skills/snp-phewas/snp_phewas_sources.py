"""Per-dataset readers: one variant x the dataset's outcomes -> one raw record per outcome.

A raw record is source-native (source alleles, source sign); orientation onto the Open
Targets effect allele happens once, in `snp_phewas_lookup`, for every dataset alike.
`record_from_fields` turns parsed file lines into a record given that file's header, and
`PanUkbbReader.read_outcome` / `FinnGenReader.read_outcome_tabix` read one outcome file.
"""

from __future__ import annotations

import logging
import math
import os
import re
import struct
import time
import zlib
from concurrent.futures import ThreadPoolExecutor
from dataclasses import dataclass, field
from functools import lru_cache
from typing import Any, Callable

import snp_phewas_schema as S
from snp_phewas_catalog import DatasetSpec, OutcomeSpec
from snp_phewas_pheweb import FinnGenPhewebClient
from snp_phewas_tabix import (
    IndexFormatError,
    SourceFileMissing,
    SourceUnreachable,
    default_index_cache,
    ensure_index,
    fetch_position,
    read_header,
)
from snp_phewas_variant import Placement

log = logging.getLogger(__name__)

DEFAULT_WORKERS = 16


@dataclass
class RawRecord:
    status: str
    src_ref: str | None = None
    src_alt: str | None = None
    beta: float | None = None
    se: float | None = None
    neglog10_p: float | None = None
    eaf: float | None = None
    low_confidence: bool = False
    variant_id_native: str | None = None
    route: str | None = None
    detail: str | None = None
    extra: dict[str, Any] = field(default_factory=dict)


def _f(v: Any) -> float | None:
    if v is None or v == "" or v == "NA":
        return None
    try:
        x = float(v)
    except (TypeError, ValueError):
        return None
    return None if math.isnan(x) else x


# --- column contracts (selected by NAME; Pan-UKBB column order varies with the strata) ---

@dataclass(frozen=True)
class Columns:
    chrom: str
    pos: str
    ref: str
    alt: str
    beta: str
    se: str
    neglog10_p: str
    eaf: str
    low_confidence: str | None = None
    eaf_fallback: str | None = None


def pan_ukbb_columns(ancestry: str) -> Columns:
    """Pan-UKBB columns for one ancestry. Continuous outcomes carry `af_<POP>`; case/control
    outcomes carry `af_cases_<POP>` / `af_controls_<POP>` instead, and the controls'
    frequency is used as the EAF (the population-frequency proxy the file offers)."""
    return Columns("chr", "pos", "ref", "alt", f"beta_{ancestry}", f"se_{ancestry}",
                   f"neglog10_pval_{ancestry}", f"af_{ancestry}", f"low_confidence_{ancestry}",
                   eaf_fallback=f"af_controls_{ancestry}")


FINNGEN_COLUMNS = Columns("chrom", "pos", "ref", "alt", "beta", "sebeta", "mlogp", "af_alt")


def pick_allele_match(records: list[dict[str, str]], cols: Columns, tgt: Placement
                      ) -> dict[str, str] | None:
    """The record at the position whose alleles are the target's, in either order."""
    want = {(tgt.ref, tgt.alt), (tgt.alt, tgt.ref)}
    for r in records:
        if (r[cols.ref].upper(), r[cols.alt].upper()) in want:
            return r
    return None


def record_from_fields(records: list[dict[str, str]], cols: Columns, tgt: Placement,
                       route: str) -> RawRecord:
    """Classify the records found at the target position for one outcome file."""
    if not records:
        return RawRecord(S.ST_ABSENT, route=route)
    r = pick_allele_match(records, cols, tgt)
    if r is None:
        found = ", ".join(f"{x[cols.ref]}>{x[cols.alt]}" for x in records)
        first = records[0]
        return RawRecord(S.ST_ALLELE_MISMATCH, src_ref=first[cols.ref], src_alt=first[cols.alt],
                         variant_id_native=f"{first[cols.chrom]}_{first[cols.pos]}_{first[cols.ref]}_{first[cols.alt]}",
                         route=route, detail=f"alleles at this position: {found}")
    vid = f"{r[cols.chrom]}_{r[cols.pos]}_{r[cols.ref]}_{r[cols.alt]}"
    if cols.beta not in r:
        return RawRecord(S.ST_NOT_TESTED, src_ref=r[cols.ref], src_alt=r[cols.alt],
                         variant_id_native=vid, route=route,
                         detail="the file has no result column for this ancestry group")
    beta, se, nlp = _f(r.get(cols.beta)), _f(r.get(cols.se)), _f(r.get(cols.neglog10_p))
    if beta is None or se is None or nlp is None:
        return RawRecord(S.ST_NOT_TESTED, src_ref=r[cols.ref], src_alt=r[cols.alt],
                         variant_id_native=vid, route=route)
    low = bool(cols.low_confidence) and (r.get(cols.low_confidence) or "").lower() == "true"
    return RawRecord(S.ST_LOW_CONFIDENCE if low else S.ST_OK,
                     src_ref=r[cols.ref], src_alt=r[cols.alt], beta=beta, se=se,
                     neglog10_p=nlp, low_confidence=low,
                     eaf=_f(r.get(cols.eaf)) if cols.eaf in r else _f(r.get(cols.eaf_fallback or "")),
                     variant_id_native=vid, route=route)


@lru_cache(maxsize=8192)
def _header(data_path: str) -> tuple[str, ...]:
    return tuple(read_header(data_path))


def release_index_cache(spec: DatasetSpec, base: str | None) -> str:
    """The index cache for ONE dataset release. A release bump reads fresh indexes rather
    than a stale one whose offsets point into a re-uploaded file (which would read as
    absent, not as an error)."""
    tag = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{spec.dataset_id}-{spec.release}")
    return os.path.join(base or default_index_cache(), tag)


def read_outcome_file(data_path: str, index_path: str, cols: Columns, tgt: Placement, *,
                      index_cache: str | None, route: str = "tabix") -> RawRecord:
    """One outcome file -> one raw record. Transport failures become outage records."""
    try:
        idx = ensure_index(index_path, index_cache)
        header = _header(data_path)
        hits = fetch_position(data_path, idx, tgt.chrom, tgt.pos)
    except SourceFileMissing as e:
        return RawRecord(S.ST_FILE_MISSING, route=route, detail=str(e))
    except SourceUnreachable as e:
        return RawRecord(S.ST_UNREACHABLE, route=route, detail=str(e))
    except (zlib.error, struct.error, IndexFormatError, UnicodeDecodeError, ValueError) as e:
        # One corrupt file is one unreadable outcome, not the end of a 7,160-file lookup.
        return RawRecord(S.ST_RECORD_UNPARSABLE, route=route, detail=f"{type(e).__name__}: {e}")
    except OSError as e:
        # A local failure (the index cache's disk is full, a permission error) is a read
        # that did not happen: one unreadable outcome, never a crashed lookup.
        return RawRecord(S.ST_UNREACHABLE, route=route,
                         detail=f"local read failed ({type(e).__name__}: {e.strerror or e})")
    records = [dict(zip(header, f)) for f in hits if len(f) == len(header)]
    if hits and not records:
        return RawRecord(S.ST_RECORD_UNPARSABLE, route=route,
                         detail=f"{len(hits)} record(s) at the position with a field count "
                                f"other than the header's {len(header)}")
    try:
        return record_from_fields(records, cols, tgt, route)
    except KeyError as e:
        # A file whose header lacks a contract column is one unreadable outcome, not a
        # crashed lookup.
        return RawRecord(S.ST_RECORD_UNPARSABLE, route=route,
                         detail=f"the file has no {e.args[0]!r} column")


def run_pool(fn: Callable[[OutcomeSpec], RawRecord], outcomes: list[OutcomeSpec],
             workers: int, label: str) -> dict[str, RawRecord]:
    """Map `fn` over outcomes with a thread pool, logging progress as it goes."""
    out: dict[str, RawRecord] = {}
    n = len(outcomes)
    step = max(1, n // 20)
    t0 = time.time()
    with ThreadPoolExecutor(max_workers=max(1, workers)) as ex:
        for i, (o, rec) in enumerate(zip(outcomes, ex.map(fn, outcomes)), 1):
            out[o.outcome_id] = rec
            if i % step == 0 or i == n:
                log.info("%s: %d/%d outcomes read (%.0fs)", label, i, n, time.time() - t0)
    return out


class PanUkbbReader:
    def __init__(self, spec: DatasetSpec, *, index_cache: str | None = None,
                 workers: int = DEFAULT_WORKERS, data_base: str | None = None,
                 index_base: str | None = None) -> None:
        self.spec = spec
        self.index_cache = release_index_cache(spec, index_cache)
        self.workers = workers
        self.data_base = data_base or spec.family_spec["data_base_url"]
        self.index_base = index_base or spec.family_spec["index_base_url"]
        self.cols = pan_ukbb_columns(spec.ancestry)

    def paths(self, o: OutcomeSpec) -> tuple[str, str]:
        return self.data_base + o.filename, self.index_base + o.filename_tabix

    def read_outcome(self, o: OutcomeSpec, tgt: Placement) -> RawRecord:
        data, idx = self.paths(o)
        return read_outcome_file(data, idx, self.cols, tgt, index_cache=self.index_cache)

    def read(self, outcomes: list[OutcomeSpec], tgt: Placement) -> tuple[dict[str, RawRecord], dict]:
        t0 = time.time()
        recs = run_pool(lambda o: self.read_outcome(o, tgt), outcomes, self.workers,
                        self.spec.dataset_id)
        return recs, {"route": "tabix", "n_files_read": len(outcomes),
                      "elapsed_s": round(time.time() - t0, 2), "workers": self.workers}


FINNGEN_ROUTES = ("hybrid", "api", "tabix")


class FinnGenReader:
    """FinnGen via the PheWeb variant API and/or per-endpoint tabix reads.

    `hybrid` (default): one API request, then a file read for every endpoint the API
    censored (it reports values only below p = 0.05; see `snp_phewas_pheweb`), so every
    tested endpoint carries its real beta / SE / p. If the API cannot be reached, every
    endpoint is read from the files instead. `api`: the one request only, censored
    endpoints labelled as such. `tabix`: a file read for every endpoint.
    """

    def __init__(self, spec: DatasetSpec, *, route: str = "hybrid",
                 client: FinnGenPhewebClient | None = None, index_cache: str | None = None,
                 workers: int = DEFAULT_WORKERS, data_base: str | None = None) -> None:
        if route not in FINNGEN_ROUTES:
            raise ValueError(f"FinnGen route must be one of {FINNGEN_ROUTES}")
        self.spec = spec
        self.route = route
        self.client = client or FinnGenPhewebClient(spec.family_spec["api_base_url"])
        self.index_cache = release_index_cache(spec, index_cache)
        self.workers = workers
        self.data_base = data_base or spec.family_spec["data_base_url"]

    def read_outcome_tabix(self, o: OutcomeSpec, tgt: Placement) -> RawRecord:
        data = self.data_base + o.filename
        return read_outcome_file(data, data + ".tbi", FINNGEN_COLUMNS, tgt,
                                 index_cache=self.index_cache)

    def _tabix(self, outcomes: list[OutcomeSpec], tgt: Placement) -> dict[str, RawRecord]:
        return run_pool(lambda o: self.read_outcome_tabix(o, tgt), outcomes, self.workers,
                        self.spec.dataset_id)

    def read(self, outcomes: list[OutcomeSpec], tgt: Placement,
             catalog_ids: set[str] | None = None,
             ) -> tuple[dict[str, RawRecord], dict, list[tuple[str, RawRecord, dict]]]:
        """(records by outcome_id, provenance, source results outside the catalog).

        `catalog_ids` is the dataset's full catalog (defaults to `outcomes`), so a result
        for an outcome that was merely not queried is not mistaken for an uncatalogued one.
        """
        t0 = time.time()
        if self.route == "tabix":
            recs = self._tabix(outcomes, tgt)
            return recs, {"route": "tabix", "n_files_read": len(outcomes),
                          "elapsed_s": round(time.time() - t0, 2)}, []
        cat = catalog_ids or {o.outcome_id for o in outcomes}
        recs, prov, uncat = self._read_api(outcomes, tgt, t0, cat)
        if self.route == "api":
            return recs, prov, uncat
        if prov.get("api_outage"):
            # The files sit on another host: an API outage is no reason to read nothing.
            recs = self._tabix(outcomes, tgt)
            prov.update(route="tabix", n_files_read=len(outcomes), api_outage_fallback=True,
                        elapsed_s=round(time.time() - t0, 2))
            return recs, prov, uncat
        if prov.get("api_absent"):
            # Confirm an API "no such variant" against one endpoint file before trusting it.
            probe = self.read_outcome_tabix(outcomes[0], tgt) if outcomes else None
            if probe is not None and probe.status not in (S.ST_ABSENT,):
                recs = self._tabix(outcomes, tgt)
                prov.update(route="tabix", n_files_read=len(outcomes),
                            api_absent_contradicted=True)
            else:
                prov["n_files_read"] = 1 if outcomes else 0
                prov["route"] = "api+tabix"
            prov["elapsed_s"] = round(time.time() - t0, 2)
            return recs, prov, uncat
        fill = [o for o in outcomes if recs[o.outcome_id].status in (S.ST_API_CENSORED, S.ST_ABSENT)]
        filled = self._tabix(fill, tgt)
        recs.update(filled)
        prov.update(route="api+tabix", n_files_read=len(fill),
                    elapsed_s=round(time.time() - t0, 2))
        return recs, prov, uncat

    def _read_api(self, outcomes: list[OutcomeSpec], tgt: Placement, t0: float,
                  catalog_ids: set[str]):
        chrom = "23" if tgt.chrom == "X" else tgt.chrom
        n_req = 0
        try:
            n_req += 1
            payload = self.client.fetch_variant(chrom, tgt.pos, tgt.ref, tgt.alt)
            src_ref, src_alt = tgt.ref, tgt.alt
            if payload is None:
                n_req += 1
                payload = self.client.fetch_variant(chrom, tgt.pos, tgt.alt, tgt.ref)
                src_ref, src_alt = tgt.alt, tgt.ref
        except SourceUnreachable as e:
            recs = {o.outcome_id: RawRecord(S.ST_UNREACHABLE, route="api", detail=str(e))
                    for o in outcomes}
            return recs, {"route": "api", "n_requests": n_req, "api_outage": True,
                          "api_error": str(e),
                          "elapsed_s": round(time.time() - t0, 2)}, []
        prov: dict[str, Any] = {"route": "api", "n_requests": n_req,
                                "elapsed_s": round(time.time() - t0, 2)}
        if payload is None:
            prov["api_absent"] = True
            return ({o.outcome_id: RawRecord(S.ST_ABSENT, route="api") for o in outcomes},
                    prov, [])
        vid = f"{tgt.chrom}_{tgt.pos}_{src_ref}_{src_alt}"
        by_code: dict[str, dict] = {}
        for res in payload["results"]:
            by_code[str(res.get("phenocode"))] = res
        recs: dict[str, RawRecord] = {}
        for o in outcomes:
            res = by_code.get(o.outcome_id)
            recs[o.outcome_id] = (self._api_record(res, src_ref, src_alt, vid) if res
                                  else RawRecord(S.ST_ABSENT, route="api",
                                                 detail="not in the variant's results"))
        uncat = [(code, self._api_record(res, src_ref, src_alt, vid), res)
                 for code, res in sorted(by_code.items()) if code not in catalog_ids]
        prov["n_results"] = len(by_code)
        prov["n_with_values"] = sum(1 for r in recs.values() if r.status == S.ST_OK)
        return recs, prov, uncat

    @staticmethod
    def _api_record(res: dict, src_ref: str, src_alt: str, vid: str) -> RawRecord:
        beta, se = _f(res.get("beta")), _f(res.get("sebeta"))
        nlp = _f(res.get("mlogp"))
        if nlp is None and _f(res.get("pval")):
            nlp = -math.log10(float(res["pval"]))
        if beta is None or se is None or nlp is None:
            return RawRecord(S.ST_API_CENSORED, src_ref=src_ref, src_alt=src_alt,
                             variant_id_native=vid, route="api")
        # The API's `maf` field is the files' `af_alt`, not a minor-allele frequency: equal on
        # the 2 variants checked 2026-10-17, including 1:109266488-G-A where alt is the major
        # allele (af_alt 0.751, API maf 0.751).
        return RawRecord(S.ST_OK, src_ref=src_ref, src_alt=src_alt, beta=beta, se=se,
                         neglog10_p=nlp, eaf=_f(res.get("maf")), variant_id_native=vid,
                         route="api",
                         extra={"n_case": res.get("n_case"), "n_control": res.get("n_control")})
