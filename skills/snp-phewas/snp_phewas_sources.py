"""Per-dataset readers: one variant x the dataset's outcomes -> one raw record per outcome.

A raw record is source-native (source alleles, source sign); orientation onto the Open
Targets effect allele happens once, in `snp_phewas_lookup`, for every dataset alike.
`record_from_fields` turns parsed file lines into a record given that file's header, and
`PanUkbbReader.read_outcome` / `FinnGenReader.read_outcome_tabix` / `MvpReader.read_outcome`
read one outcome file.
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
from statistics import NormalDist
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
    # A direction the SOURCE could not establish (GWAS Catalog `hm_code` 9 / 14-18); the
    # lookup's own orientation can still withhold one for its own reason.
    direction_withheld: str | None = None


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
    """Column names of one file layout. `beta` / `se` / `neglog10_p` may be None when the
    layout carries them in another form: an odds ratio with its 95% interval instead of
    beta / SE (`odds_ratio`, `ci_lower`, `ci_upper`), and a probability instead of -log10 p
    (`p_value`). `hm_code` is the GWAS Catalog harmonizer's per-record code."""

    chrom: str
    pos: str
    ref: str
    alt: str
    beta: str | None
    se: str | None
    neglog10_p: str | None
    eaf: str
    low_confidence: str | None = None
    eaf_fallback: str | None = None
    p_value: str | None = None
    odds_ratio: str | None = None
    ci_lower: str | None = None
    ci_upper: str | None = None
    hm_code: str | None = None


def pan_ukbb_columns(ancestry: str) -> Columns:
    """Pan-UKBB columns for one ancestry. Continuous outcomes carry `af_<POP>`; case/control
    outcomes carry `af_cases_<POP>` / `af_controls_<POP>` instead, and the controls'
    frequency is used as the EAF (the population-frequency proxy the file offers)."""
    return Columns("chr", "pos", "ref", "alt", f"beta_{ancestry}", f"se_{ancestry}",
                   f"neglog10_pval_{ancestry}", f"af_{ancestry}", f"low_confidence_{ancestry}",
                   eaf_fallback=f"af_controls_{ancestry}")


FINNGEN_COLUMNS = Columns("chrom", "pos", "ref", "alt", "beta", "sebeta", "mlogp", "af_alt")

# GWAS Catalog harmonized files (GWAS-SSF v1.0). The harmonizer orients `effect_allele` to the
# Ensembl ALT where it can, and its `variant_id` is chrom_pos_other_effect, so the file is
# reference-anchored with ref = `other_allele`, alt = `effect_allele`. Quantitative studies
# carry beta / SE; case/control studies carry an odds ratio and its 95% interval with
# `standard_error = NA`. `p_value` is a probability and underflows to 0 at strong hits.
MVP_COLUMNS_QUANT = Columns("chromosome", "base_pair_location", "other_allele", "effect_allele",
                            "beta", "standard_error", None, "effect_allele_frequency",
                            p_value="p_value", hm_code="hm_code")
MVP_COLUMNS_BINARY = Columns("chromosome", "base_pair_location", "other_allele",
                             "effect_allele", None, None, None, "effect_allele_frequency",
                             p_value="p_value", odds_ratio="odds_ratio", ci_lower="ci_lower",
                             ci_upper="ci_upper", hm_code="hm_code")
#: Every code the GWAS Catalog harmonizer writes to `hm_code`, as code -> (meaning,
#: direction withheld), copied from the Catalog's harmonization code table
#: (https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics, read 2026-10-24). 1-8 and
#: 10-13 are oriented records and keep their direction; 9 and 14-18 are records the
#: harmonizer could NOT orient to the forward strand, so their direction is withheld, never
#: guessed. The two MVP files read here carry only 5, 6 and 10-13 (see SKILL.md); the rule
#: covers the codes they do not exercise, such as 9 (palindromic mode `drop`). A code
#: outside this table, empty included, makes the record unparsable rather than silently kept.
HM_CODES: dict[str, tuple[str, bool]] = {
    "1": ("Palindromic; Infer strand; Forward strand; Alleles correct", False),
    "2": ("Palindromic; Infer strand; Forward strand; Flipped alleles", False),
    "3": ("Palindromic; Infer strand; Reverse strand; Alleles correct", False),
    "4": ("Palindromic; Infer strand; Reverse strand; Flipped alleles", False),
    "5": ("Palindromic; Assume forward strand; Alleles correct", False),
    "6": ("Palindromic; Assume forward strand; Flipped alleles", False),
    "7": ("Palindromic; Assume reverse strand; Alleles correct", False),
    "8": ("Palindromic; Assume reverse strand; Flipped alleles", False),
    "9": ("Palindromic; Drop palindromic; Not orientated", True),
    "10": ("Forward strand; Alleles correct", False),
    "11": ("Forward strand; Flipped alleles", False),
    "12": ("Reverse strand; Alleles correct", False),
    "13": ("Reverse strand; Flipped alleles", False),
    "14": ("Required fields are not known; Not orientated", True),
    "15": ("No matching variants in reference VCF; Not orientated", True),
    "16": ("Multiple matching variants in reference VCF; Not orientated", True),
    "17": ("Palindromic; Infer strand; EAF or reference VCF AF not known; Not orientated", True),
    "18": ("Palindromic; Infer strand; EAF < specified minor allele frequency threshold; "
           "Not orientated", True),
}
#: The codes whose direction is withheld (9 and 14-18), derived from the table.
HM_CODES_NOT_ORIENTATED = frozenset(c for c, (_, withheld) in HM_CODES.items() if withheld)
Z_95 = 1.959963984540054  # the two-sided 95% normal quantile an interval is built on
P_UNDERFLOW_DETAIL = "p underflowed to 0 in the source; -log10 p derived from beta/SE"


def mvp_columns_for(header: list[str] | tuple[str, ...]) -> Columns:
    """The MVP layout of one file, chosen by which effect column its header carries."""
    if "beta" in header:
        return MVP_COLUMNS_QUANT
    if "odds_ratio" in header:
        return MVP_COLUMNS_BINARY
    raise KeyError("beta")  # neither: one unreadable outcome (read_outcome_file catches it)


def neglog10_p_from_z(z: float) -> float:
    """-log10 of the two-sided normal p for `z`, in log space so it never underflows.

    p = erfc(|z| / sqrt 2). Below x = |z| / sqrt 2 = 25 `math.erfc` does not underflow
    (erfc(25) is about 8e-274); above it the asymptotic expansion
    ln erfc(x) = -x^2 - ln(x sqrt(pi)) + ln(1 - 1/(2x^2) + 3/(4x^4) - 15/(8x^6))
    is used. Against scipy's `norm.logsf` the absolute difference is below 2e-11 in
    -log10 p over |z| from 0.1 to 300 (measured 2026-10-19; the unit test re-checks it
    whenever scipy is installed).
    """
    x = abs(z) / math.sqrt(2.0)
    if x < 25.0:
        return -math.log10(math.erfc(x))
    x2 = x * x
    series = 1.0 - 1.0 / (2 * x2) + 3.0 / (4 * x2 * x2) - 15.0 / (8 * x2 ** 3)
    ln_p = -x2 - math.log(x * math.sqrt(math.pi)) + math.log(series)
    return -ln_p / math.log(10.0)


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
    effect_col = cols.beta or cols.odds_ratio
    if effect_col not in r:
        return RawRecord(S.ST_NOT_TESTED, src_ref=r[cols.ref], src_alt=r[cols.alt],
                         variant_id_native=vid, route=route,
                         detail="the file has no result column for this ancestry group")
    base = dict(src_ref=r[cols.ref], src_alt=r[cols.alt], variant_id_native=vid, route=route)
    hm_code = None
    if cols.hm_code:  # a harmonized file: every record carries a code from the Catalog's table
        if cols.hm_code not in r:
            return RawRecord(S.ST_RECORD_UNPARSABLE, **base,
                             detail="the harmonized file has no hm_code column")
        hm_code = (r.get(cols.hm_code) or "").strip()
        if hm_code not in HM_CODES:
            return RawRecord(S.ST_RECORD_UNPARSABLE, **base,
                             detail=f"hm_code {hm_code!r} is not in the GWAS Catalog "
                                    f"harmonizer's code table")
    detail = None
    if cols.beta:
        beta, se = _f(r.get(cols.beta)), _f(r.get(cols.se))
    else:  # an odds ratio: beta = ln OR; SE from the source's p, else from the 95% interval
        orv, lo, hi = (_f(r.get(c)) for c in (cols.odds_ratio, cols.ci_lower, cols.ci_upper))
        if orv is None:
            return RawRecord(S.ST_NOT_TESTED, **base)
        if orv <= 0:
            return RawRecord(S.ST_RECORD_UNPARSABLE, **base,
                             detail=f"odds ratio {orv} has no log")
        beta = math.log(orv)
        # the interval is the fallback only, so a record with an odds ratio and a p but no
        # interval is still a value
        se = ((math.log(hi) - math.log(lo)) / (2 * Z_95)
              if lo is not None and hi is not None and 0 < lo < hi else None)
        # The GWAS Catalog harmonized MVP files carry the odds ratio and interval rounded to
        # about three significant digits (1.015228... is 1 / 0.985), so an interval-width SE
        # disagrees with the file's own p on many records: measured on the vendored IHD slice,
        # 85 of 150 records differ by more than 1e-3 in p or 0.05 in log10 p. Where the file's
        # p is a usable probability the SE is the one that reproduces it (z from p, SE =
        # |beta| / z); the interval is used only when p cannot give it (p = 0, 1 or a null
        # effect).
        p_src = _f(r.get(cols.p_value or ""))
        if p_src is not None and 0.0 < p_src < 1.0 and beta != 0.0 and p_src / 2 > 0.0:
            z = -NormalDist().inv_cdf(p_src / 2)
            if z > 0:
                se = abs(beta) / z
    if cols.neglog10_p:
        nlp = _f(r.get(cols.neglog10_p))
    else:
        p = _f(r.get(cols.p_value or ""))
        nlp = None
        if p is not None:
            if 0.0 < p <= 1.0:
                nlp = -math.log10(p)
            elif p == 0.0:
                if beta is None or not se or se <= 0 or not math.isfinite(beta / se):
                    return RawRecord(S.ST_RECORD_UNPARSABLE, **base,
                                     detail="p underflowed to 0 and no usable beta/SE")
                nlp, detail = neglog10_p_from_z(beta / se), P_UNDERFLOW_DETAIL
            else:
                return RawRecord(S.ST_RECORD_UNPARSABLE, **base,
                                 detail=f"p value {p} is outside [0, 1]")
    if beta is None or se is None or nlp is None:
        return RawRecord(S.ST_NOT_TESTED, **base)
    if not (math.isfinite(beta) and math.isfinite(se)):
        return RawRecord(S.ST_RECORD_UNPARSABLE, **base, detail="beta or SE is not finite")
    low = bool(cols.low_confidence) and (r.get(cols.low_confidence) or "").lower() == "true"
    withheld = S.WITHHELD_SOURCE_NOT_ORIENTATED if hm_code in HM_CODES_NOT_ORIENTATED else None
    return RawRecord(S.ST_LOW_CONFIDENCE if low else S.ST_OK, beta=beta, se=se,
                     neglog10_p=nlp, low_confidence=low,
                     eaf=_f(r.get(cols.eaf)) if cols.eaf in r else _f(r.get(cols.eaf_fallback or "")),
                     detail=detail, direction_withheld=withheld, **base)


@lru_cache(maxsize=8192)
def _header(data_path: str) -> tuple[str, ...]:
    return tuple(read_header(data_path))


def release_index_cache(spec: DatasetSpec, base: str | None) -> str:
    """The index cache for ONE dataset release. A release bump reads fresh indexes rather
    than a stale one whose offsets point into a re-uploaded file (which would read as
    absent, not as an error)."""
    tag = re.sub(r"[^A-Za-z0-9._-]+", "_", f"{spec.dataset_id}-{spec.release}")
    return os.path.join(base or default_index_cache(), tag)


def read_outcome_file(data_path: str, index_path: str, cols: Columns | Callable[..., Columns],
                      tgt: Placement, *, index_cache: str | None, route: str = "tabix"
                      ) -> RawRecord:
    """One outcome file -> one raw record. Transport failures become outage records.

    `cols` is a layout, or a function choosing the layout from the file's header (the GWAS
    Catalog writes case/control and quantitative studies with different effect columns)."""
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
        if callable(cols):
            cols = cols(header)
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


def not_readable_record(o: OutcomeSpec) -> RawRecord | None:
    """`source_file_not_readable` for a catalog outcome the source holds with no
    position-indexed file (an MVP study the GWAS Catalog has not harmonized), with the
    reason first in `detail` (`raw_only_unindexed: ...`); None for a readable outcome. Known
    from the catalog, so it costs no network call and holds whatever the lookup's scope."""
    if o.access_status is None or (o.access_status == "harmonised" and o.filename):
        return None
    from snp_phewas_mvp_catalog import ACCESS_STATUS_LABELS

    label = ACCESS_STATUS_LABELS.get(o.access_status, "no position-indexed file")
    return RawRecord(S.ST_NOT_READABLE, route="catalog", detail=f"{o.access_status}: {label}")


class MvpReader:
    """Million Veteran Program studies through the GWAS Catalog's harmonized files.

    A catalog row the Catalog has not harmonized (`access_status` raw_only_unindexed or
    no_summary_statistics) has no position-indexed file: it is answered
    `source_file_not_readable` with the reason, without a network call."""

    def __init__(self, spec: DatasetSpec, *, index_cache: str | None = None,
                 workers: int = DEFAULT_WORKERS, data_base: str | None = None) -> None:
        self.spec = spec
        self.index_cache = release_index_cache(spec, index_cache)
        self.workers = workers
        self.data_base = data_base or spec.family_spec["data_base_url"]

    def paths(self, o: OutcomeSpec) -> tuple[str, str]:
        return self.data_base + o.filename, self.data_base + (o.filename_tabix or o.filename + ".tbi")

    def read_outcome(self, o: OutcomeSpec, tgt: Placement) -> RawRecord:
        unreadable = not_readable_record(o)
        if unreadable is not None:
            return unreadable
        data, idx = self.paths(o)
        return read_outcome_file(data, idx, mvp_columns_for, tgt, index_cache=self.index_cache)

    def read(self, outcomes: list[OutcomeSpec], tgt: Placement) -> tuple[dict[str, RawRecord], dict]:
        t0 = time.time()
        recs = run_pool(lambda o: self.read_outcome(o, tgt), outcomes, self.workers,
                        self.spec.dataset_id)
        n_files = sum(1 for o in outcomes if o.access_status == "harmonised" and o.filename)
        return recs, {"route": "tabix", "n_files_read": n_files,
                      "n_not_readable": len(outcomes) - n_files,
                      "elapsed_s": round(time.time() - t0, 2), "workers": self.workers}


FINNGEN_ROUTES = ("hybrid", "api", "tabix")


class FinnGenReader:
    """FinnGen via the PheWeb variant API and/or per-endpoint tabix reads.

    `hybrid` (default): one API request, then a file read for every endpoint the API
    censored (it reports values only below p = 0.05; see `snp_phewas_pheweb`), so every
    tested endpoint carries its real beta / SE / p. If the API cannot be reached, every
    endpoint is read from the files instead. `api`: the one request only, censored
    endpoints labeled as such. `tabix`: a file read for every endpoint.
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
        detail = None
        if nlp is None:
            # The same p handling as the file route: a p of 0 is an underflow, not a censored
            # value, so it is derived from beta/SE (a truthiness test read it as p >= 0.05).
            p = _f(res.get("pval"))
            if p is not None and 0.0 < p <= 1.0:
                nlp = -math.log10(p)
            elif (p == 0.0 and beta is not None and se is not None and se > 0
                  and math.isfinite(beta / se)):
                nlp, detail = neglog10_p_from_z(beta / se), P_UNDERFLOW_DETAIL
        if beta is None or se is None or nlp is None:
            return RawRecord(S.ST_API_CENSORED, src_ref=src_ref, src_alt=src_alt,
                             variant_id_native=vid, route="api")
        # The API's `maf` field is the files' `af_alt`, not a minor-allele frequency: equal on
        # the 2 variants checked 2026-10-17, including 1:109266488-G-A where alt is the major
        # allele (af_alt 0.751, API maf 0.751).
        return RawRecord(S.ST_OK, src_ref=src_ref, src_alt=src_alt, beta=beta, se=se,
                         neglog10_p=nlp, eaf=_f(res.get("maf")), variant_id_native=vid,
                         route="api", detail=detail,
                         extra={"n_case": res.get("n_case"), "n_control": res.get("n_control")})
