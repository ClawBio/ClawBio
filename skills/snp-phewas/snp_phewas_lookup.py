"""`lookup_variant`: one variant -> one row per (dataset, outcome), every outcome accounted.

Every value is read live from the sources. `readers=` injects a reader per dataset id
(the offline demo and tests use it; a reader that returns the same records from a
local mirror or a precomputed table plugs in here too).

Scope: by default a Pan-UKBB dataset is read only over the outcomes the cross-biobank
crosswalk pairs with a FinnGen endpoint (`scope="harmonized"`); every other Pan-UKBB
outcome is still a row, with status `not_queried`, and the result's notes say so in plain
words. `scope="all"` reads every outcome.
"""

from __future__ import annotations

import logging
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Callable, Iterable, Mapping

import snp_phewas_orient as O
import snp_phewas_schema as S
from snp_phewas_catalog import (
    DatasetSpec,
    OutcomeSpec,
    get_dataset,
    load_outcome_catalog,
    not_served,
    resolve_dataset_ids,
)
from snp_phewas_crosswalk import harmonized_finngen_endpoints, harmonized_pan_ukbb_outcomes
from snp_phewas_tabix import SourceUnreachable
from snp_phewas_sources import (
    DEFAULT_WORKERS,
    FinnGenReader,
    PanUkbbReader,
    RawRecord,
)
from snp_phewas_variant import (
    MapFn,
    Placement,
    ResolvedVariant,
    map_placement,
    resolve_variant,
)

log = logging.getLogger(__name__)

SCOPES = ("harmonized", "all")


@dataclass
class SnpPhewasResult:
    variant: ResolvedVariant
    dataset_ids: list[str]
    rows: list[S.OutcomeRow]
    uncatalogued: list[S.OutcomeRow]
    provenance: dict[str, dict[str, Any]]
    notes: list[str] = field(default_factory=list)
    scope: str = "all"

    def dataset_rows(self, dataset_id: str) -> list[S.OutcomeRow]:
        return [r for r in self.rows if r.dataset_id == dataset_id]

    @property
    def tested(self) -> list[S.OutcomeRow]:
        return [r for r in self.rows if r.status in S.TESTED_STATUSES]

    @property
    def excluded(self) -> list[S.OutcomeRow]:
        return [r for r in self.rows if r.status not in S.TESTED_STATUSES]

    def status_counts(self, dataset_id: str) -> Counter:
        return Counter(r.status for r in self.dataset_rows(dataset_id))

    def attributions(self) -> list[str]:
        return [get_dataset(d).attribution for d in self.dataset_ids]


def _row(spec: DatasetSpec, o: OutcomeSpec, v: ResolvedVariant, tgt: Placement | None,
         rec: RawRecord) -> S.OutcomeRow:
    g = v.grch38
    row = S.OutcomeRow(
        dataset_id=spec.dataset_id, dataset_release=spec.release, outcome_id=o.outcome_id,
        outcome_label=o.outcome_label, outcome_category_native=o.outcome_category_native,
        outcome_kind=o.outcome_kind, n_cases=o.n_cases, n_controls=o.n_controls,
        n_total=o.n_total, variant_id_grch38=g.variant_id,
        variant_id_native=rec.variant_id_native, rs_id=v.rs_id, effect_allele=g.alt,
        other_allele=g.ref, status=rec.status, source_route=rec.route,
        status_detail=rec.detail, low_confidence=rec.low_confidence,
    )
    if tgt is not None and rec.src_ref and rec.src_alt:
        ori = O.orient(rec.src_ref, rec.src_alt, tgt.ref, tgt.alt)
        row.palindromic = ori.palindromic
        if not ori.matches:
            row.status = S.ST_ALLELE_MISMATCH
            return row
        if rec.status in S.TESTED_STATUSES:
            beta, eaf = O.apply(ori, rec.beta, rec.eaf)
            row.beta, row.eaf = beta, eaf
            row.se, row.neglog10_p = rec.se, rec.neglog10_p
            row.allele_flip_applied = ori.flip
            row.direction_withheld_reason = ori.withheld_reason
    return row


def _target(spec: DatasetSpec, v: ResolvedVariant, map_fn: MapFn | None,
            cache: dict[str, tuple[Placement | None, str, str | None]]
            ) -> tuple[Placement | None, str, str | None]:
    """(placement on the dataset's build, note, failure status or None)."""
    if spec.build == "GRCh38":
        return v.grch38, "", None
    if spec.build not in cache:
        try:
            p, note = map_placement(v.grch38, spec.build, map_fn=map_fn)
            cache[spec.build] = (p, note, None if p else S.ST_BUILD_UNMAPPABLE)
        except SourceUnreachable as e:
            cache[spec.build] = (None, f"the {spec.build} position could not be looked up ({e})",
                                 S.ST_UNREACHABLE)
    return cache[spec.build]


def _reader_for(spec: DatasetSpec, *, finngen_route: str, workers: int,
                index_cache: str | None):
    if spec.family == "pan_ukbb":
        return PanUkbbReader(spec, index_cache=index_cache, workers=workers)
    if spec.family == "finngen":
        return FinnGenReader(spec, route=finngen_route, index_cache=index_cache, workers=workers)
    raise ValueError(f"no reader for {spec.family}")  # pragma: no cover


def lookup_variant(
    variant: str,
    *,
    datasets: str | list[str] | None = None,
    scope: str = "harmonized",
    build: str = "GRCh38",
    outcomes: Mapping[str, Iterable[str]] | None = None,
    finngen_route: str = "hybrid",
    workers: int = DEFAULT_WORKERS,
    index_cache: str | None = None,
    readers: Mapping[str, Any] | None = None,
    rsid_resolver: Callable[[str], str] | None = None,
    map_fn: MapFn | None = None,
    harmonize: bool = True,
) -> SnpPhewasResult:
    """Look one variant up in every outcome of each curated dataset.

    `outcomes` optionally restricts the read per dataset ({dataset_id: outcome ids}) and
    overrides `scope` for that dataset; the other catalog outcomes still get a row, with
    status `not_queried`. `readers` injects a reader per dataset_id (demo, tests, local
    mirrors).
    """
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES}")
    ids = resolve_dataset_ids(datasets)
    v = resolve_variant(variant, build=build, rsid_resolver=rsid_resolver, map_fn=map_fn)
    notes = [v.resolution_note + "."]
    rows: list[S.OutcomeRow] = []
    uncat_rows: list[S.OutcomeRow] = []
    provenance: dict[str, dict[str, Any]] = {}
    build_cache: dict[str, tuple[Placement | None, str, str | None]] = {}

    for dsid in ids:
        spec = get_dataset(dsid)
        catalog = load_outcome_catalog(dsid)
        cat_ids = {o.outcome_id for o in catalog}
        want = set(outcomes[dsid]) if outcomes and dsid in outcomes else None
        if want is not None and want - cat_ids:
            raise ValueError(f"{dsid}: unknown outcome ids {sorted(want - cat_ids)[:5]}")
        # Two views: the outcomes matched across the two biobanks (the default; what the
        # replication view compares, on both sides) or every outcome (--all-outcomes). The cut
        # exists for the replication view, so a lookup of one biobank reads its whole catalog.
        other = "finngen" if spec.family == "pan_ukbb" else "pan_ukbb"
        scoped = (want is None and scope == "harmonized"
                  and spec.family in ("pan_ukbb", "finngen")
                  and any(get_dataset(d).family == other for d in ids))
        if scoped:
            matched = (harmonized_pan_ukbb_outcomes() if spec.family == "pan_ukbb"
                       else harmonized_finngen_endpoints())
            want = matched & cat_ids
        queried = [o for o in catalog if want is None or o.outcome_id in want]
        prov: dict[str, Any] = {"release": spec.release, "build": spec.build,
                                "attribution": spec.attribution,
                                "n_catalog": len(catalog), "n_queried": len(queried),
                                "scope": ("harmonized" if scoped else
                                          "requested" if want is not None else "all")}
        uncat: list = []
        tgt, map_note, fail = _target(spec, v, map_fn, build_cache)
        if map_note:
            prov["build_note"] = map_note
        recs: dict[str, RawRecord] = {}
        if fail is not None:
            recs = {o.outcome_id: RawRecord(fail, detail=map_note) for o in queried}
        else:
            reader = (readers or {}).get(dsid) or _reader_for(
                spec, finngen_route=finngen_route, workers=workers, index_cache=index_cache)
            if spec.family == "finngen":
                recs, rprov, uncat = reader.read(queried, tgt, cat_ids)
            else:
                recs, rprov = reader.read(queried, tgt)
            prov.update(rprov)
        lost = {o.outcome_id for o in queried} - set(recs)
        if lost:
            raise S.PartitionError(f"{dsid}: the reader returned no record for "
                                   f"{len(lost)} queried outcomes (e.g. {sorted(lost)[:3]})")
        ds_rows = [_row(spec, o, v, tgt, recs.get(o.outcome_id, RawRecord(S.ST_NOT_QUERIED)))
                   for o in catalog]
        S.assert_partition(dsid, ds_rows, cat_ids)
        for code, rec, res in uncat:
            o = OutcomeSpec(dsid, code, str(res.get("phenostring") or code),
                            str(res.get("category") or ""), "binary", None,
                            res.get("n_case"), res.get("n_control"), None)
            r = _row(spec, o, v, tgt, rec)
            r.status_detail = (f"source value status was {r.status}; the outcome is not in "
                               f"the published catalog")
            r.status = S.ST_NOT_IN_CATALOG
            uncat_rows.append(r)
        prov["status_counts"] = dict(Counter(r.status for r in ds_rows))
        prov["n_uncatalogued"] = len(uncat)
        provenance[dsid] = prov
        rows.extend(ds_rows)
        notes.extend(_dataset_notes(spec, ds_rows, prov, len(uncat)))

    others = not_served()
    if others:
        notes.append(" ".join(others.values()))
    result = SnpPhewasResult(variant=v, dataset_ids=ids, rows=rows, uncatalogued=uncat_rows,
                             provenance=provenance, notes=notes, scope=scope)
    if harmonize:
        from snp_phewas_harmonize import attach_harmonization

        attach_harmonization(result)
    return result


def _n(k: int, one: str, many: str) -> str:
    return f"{k:,} {one if k == 1 else many}"


def _dataset_notes(spec: DatasetSpec, rows: list[S.OutcomeRow], prov: dict, n_uncat: int
                   ) -> list[str]:
    c = Counter(r.status for r in rows)
    n = len(rows)
    tested = sum(c[s] for s in S.TESTED_STATUSES)
    out = [f"{spec.label}: tested in {tested:,} of {n:,} outcomes."]
    untested = {s: k for s, k in c.items() if s not in S.TESTED_STATUSES}
    if untested:
        top, k = max(untested.items(), key=lambda kv: kv[1])
        m = n - tested
        if m == k:
            out.append(f"{spec.label}: {_n(m, 'outcome has', 'outcomes have')} no value: "
                       f"{S.STATUS_LABELS[top]}.")
        else:
            out.append(f"{spec.label}: {_n(m, 'outcome has', 'outcomes have')} no value; the "
                       f"most common reason ({k:,}) is: {S.STATUS_LABELS[top]}.")
    if c[S.ST_UNREACHABLE]:
        out.append(f"{spec.label}: {_n(c[S.ST_UNREACHABLE], 'outcome', 'outcomes')} could not be read because "
                   f"the source was unreachable; they are not evidence of no association.")
    if c[S.ST_LOW_CONFIDENCE]:
        out.append(f"{spec.label}: {_n(c[S.ST_LOW_CONFIDENCE], 'value is', 'values are')} flagged low confidence "
                   f"by the source (too few carriers) and are kept out of the replication view.")
    withheld = sum(1 for r in rows if r.direction_withheld_reason)
    if withheld:
        out.append(f"{spec.label}: no direction shown for {_n(withheld, 'outcome', 'outcomes')}: "
                   + S.DIRECTION_WITHHELD_LABELS[S.WITHHELD_PALINDROMIC_ORDER].split(": ", 1)[1] + ".")
    flipped = sum(1 for r in rows if r.allele_flip_applied)
    if flipped:
        out.append(f"{spec.label}: the source lists the alleles in the opposite order, so "
                   f"{flipped:,} effects were re-signed onto the {rows[0].effect_allele} allele.")
    if prov.get("api_outage_fallback"):
        if any(r.status not in S.TRANSIENT_STATUSES for r in rows
               if r.status != S.ST_NOT_QUERIED):
            out.append(f"{spec.label}: the web API could not be reached, so every endpoint "
                       f"was read from the files instead.")
        else:
            out.append(f"{spec.label}: the web API could not be reached, and the files could "
                       f"not be read either; no value is shown, and this is an outage, not "
                       f"an absence.")
    if prov.get("route") in ("tabix", "api+tabix") and prov.get("elapsed_s") is not None:
        api = " plus one web-API request" if prov["route"] == "api+tabix" else ""
        out.append(f"{spec.label}: read from {prov['n_files_read']:,} files{api} in "
                   f"{prov['elapsed_s']:.0f} s.")
    if prov.get("scope") == "harmonized" and prov["n_queried"] < prov["n_catalog"]:
        out.append(f"{spec.label}: read over the {prov['n_queried']:,} outcomes matched to the "
                   f"other biobank, not all {prov['n_catalog']:,} in the catalog; the other "
                   f"{prov['n_catalog'] - prov['n_queried']:,} are listed as not part of this "
                   f"lookup. Use --all-outcomes to read every outcome in both biobanks.")
    elif prov.get("n_queried") is not None and prov["n_queried"] < prov["n_catalog"]:
        out.append(f"{spec.label}: only {prov['n_queried']:,} of {prov['n_catalog']:,} outcomes "
                   f"were requested; the rest are listed as not part of this lookup.")
    if n_uncat:
        out.append(f"{spec.label}: the source also returned {n_uncat} outcome(s) its published "
                   f"catalog does not list; they are kept separately, not in the counts above.")
    if spec.build != "GRCh38" and prov.get("build_note"):
        out.append(f"{spec.label} uses {spec.build}: {prov['build_note']}.")
    return out
