"""`lookup_variant`: one variant -> one row per (dataset, outcome), every outcome accounted.

Every value is read live from the sources. `readers=` injects a reader per dataset id
(the offline demo and tests use it; a reader that returns the same records from a
local mirror or a precomputed table plugs in here too).

Scope: with two or more biobanks in the lookup, each biobank is read by default only over
its outcomes with a readable counterpart in ANOTHER biobank of the same lookup
(`scope="harmonized"`, `matched_outcome_ids`): one rule for every biobank, never more of one
because its data is quicker to get. Every other outcome is still a row, with status
`not_queried`, and the result's notes say so in plain words. `scope="all"` reads every
outcome in every biobank; a lookup of one biobank alone reads its whole catalog.
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
from snp_phewas_crosswalk import is_comparable, load_crosswalk
from snp_phewas_tabix import SourceUnreachable
from snp_phewas_sources import (
    DEFAULT_WORKERS,
    FinnGenReader,
    MvpReader,
    PanUkbbReader,
    RawRecord,
    not_readable_record,
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


def _family_names(families: list[str] | None) -> str:
    from snp_phewas_catalog import family_label

    names = [family_label(f) for f in (families or [])]
    return " or ".join(names) if names else "another biobank"


def matched_outcome_ids(dataset_id: str, dataset_ids: list[str]) -> frozenset[str]:
    """Outcome ids of `dataset_id` with a counterpart the replication view can actually
    compare: an outcome of ANOTHER family's dataset in this same lookup that shares its
    harmonized id and has a readable file. A study the source holds with no readable file is
    no counterpart (it can never be compared), and neither is an outcome of a dataset outside
    the lookup (another MVP ancestry group, say)."""
    cw = load_crosswalk()
    fam = get_dataset(dataset_id).family

    def keys(dsid: str) -> dict[str, str]:
        f = get_dataset(dsid).family
        out = {}
        for o in load_outcome_catalog(dsid):
            r = cw.get((f, o.outcome_id))
            if r is not None and is_comparable(r) and not_readable_record(o) is None:
                out[o.outcome_id] = r.harmonized_outcome_id
        return out

    others: set[str] = set()
    for d in dataset_ids:
        if get_dataset(d).family != fam:
            others.update(keys(d).values())
    return frozenset(oid for oid, h in keys(dataset_id).items() if h in others)


def matched_scope_sentence(label: str, rows: list[S.OutcomeRow], prov: dict[str, Any]) -> str:
    """What a matched-outcome read covered, in plain language, counted from the rows (a count
    taken from the request would call unreadable studies "read")."""
    n_cat = prov["n_catalog"]
    read = sum(1 for r in rows if r.status not in (S.ST_NOT_QUERIED, S.ST_NOT_READABLE))
    not_looked = sum(1 for r in rows if r.status == S.ST_NOT_QUERIED)
    others = _family_names(prov.get("scope_families"))
    s = (f"{label}: read over the {read:,} outcomes with a readable counterpart in "
         f"{others} (the ones the replication view compares), not all {n_cat:,}; "
         f"{not_looked:,} were not looked up")
    n_nr = sum(1 for r in rows if r.status == S.ST_NOT_READABLE)
    if n_nr:
        s += f" and {n_nr:,} have no readable file"
    return (s + "; none of them is counted as a null. Run with --all-outcomes to read every "
            "outcome in every biobank.")


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
            # the lookup's own reason first; else one the source gave (it could not orient
            # the record to the forward strand). Either way beta and EAF are withheld, as
            # the palindromic rule withholds them; SE and p do not depend on orientation.
            withheld = ori.withheld_reason or rec.direction_withheld
            if withheld:
                beta, eaf = None, None
            row.beta, row.eaf = beta, eaf
            row.se, row.neglog10_p = rec.se, rec.neglog10_p
            row.allele_flip_applied = ori.flip and not withheld
            row.direction_withheld_reason = withheld
    return row


def _not_looked_up(o: OutcomeSpec) -> RawRecord:
    """The record of a catalog outcome this lookup did not read: `not_queried`, unless the
    source holds it with no readable file at all, which is said whatever was asked."""
    return not_readable_record(o) or RawRecord(S.ST_NOT_QUERIED)


def _target(spec: DatasetSpec, v: ResolvedVariant, map_fn: MapFn | None,
            cache: dict[str, tuple[Placement | None, str, str | None]]
            ) -> tuple[Placement | None, str, str | None]:
    """(placement on the dataset's build, note, failure status or None).

    A dataset that released autosomal results only (`contigs: autosomes`, MVP) answers an
    X / Y / MT variant `absent_in_dataset` for every outcome without reading a file."""
    if (spec.family_spec.get("contigs") == "autosomes"
            and v.grch38.chrom.removeprefix("chr") not in S.AUTOSOMES):
        return None, S.AUTOSOMES_ONLY_DETAIL, S.ST_ABSENT
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
    if spec.family == "mvp":
        return MvpReader(spec, index_cache=index_cache, workers=workers)
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
    rsid_resolver: Callable[[str], str | list[str]] | None = None,
    allele_weight: Callable[[str], tuple[int, float]] | None = None,
    map_fn: MapFn | None = None,
    harmonize: bool = True,
) -> SnpPhewasResult:
    """Look one variant up in every outcome of each curated dataset.

    `outcomes` optionally restricts the read per dataset ({dataset_id: outcome ids}) and
    overrides `scope` for that dataset; the other catalog outcomes still get a row, with
    status `not_queried`. `readers` injects a reader per dataset_id (demo, tests, local
    mirrors).

    `scope="harmonized"` (the default): when the lookup includes another biobank family,
    each dataset is read over the outcomes with a READABLE counterpart in another family's
    dataset of the same lookup (`matched_outcome_ids`), the ones the replication view
    compares. Counted on the committed crosswalk: Pan-UK Biobank EUR 589 and FinnGen 622
    with those two; with MVP EUR added (the default three), 1,156 / 622 / 727. The SAME rule
    in every biobank: FinnGen's one API request returns every endpoint, and the ones outside
    the set are not kept. Said in the notes. `scope="all"` reads every outcome.
    """
    if scope not in SCOPES:
        raise ValueError(f"scope must be one of {SCOPES}")
    ids = resolve_dataset_ids(datasets)
    v = resolve_variant(variant, build=build, rsid_resolver=rsid_resolver, map_fn=map_fn,
                        allele_weight=allele_weight)
    notes = [v.resolution_note + "."]
    rows: list[S.OutcomeRow] = []
    uncat_rows: list[S.OutcomeRow] = []
    provenance: dict[str, dict[str, Any]] = {}
    build_cache: dict[str, tuple[Placement | None, str, str | None]] = {}

    def one_dataset(dsid: str):
        """One dataset's rows, uncataloged rows, provenance and notes (no shared writes
        besides the build-placement cache, which holds one answer per build)."""
        ds_uncat: list[S.OutcomeRow] = []
        spec = get_dataset(dsid)
        catalog = load_outcome_catalog(dsid)
        cat_ids = {o.outcome_id for o in catalog}
        want = set(outcomes[dsid]) if outcomes and dsid in outcomes else None
        if want is not None and want - cat_ids:
            raise ValueError(f"{dsid}: unknown outcome ids {sorted(want - cat_ids)[:5]}")
        others = sorted({get_dataset(d).family for d in ids} - {spec.family})
        # The cut exists for the replication view, so a lookup of one biobank reads its
        # whole catalog.
        scoped = want is None and scope == "harmonized" and bool(others)
        if scoped:
            want = set(matched_outcome_ids(dsid, ids))
        queried = [o for o in catalog if want is None or o.outcome_id in want]
        prov: dict[str, Any] = {"release": spec.release, "build": spec.build,
                                "attribution": spec.attribution,
                                "n_catalog": len(catalog), "n_queried": len(queried),
                                "scope": ("harmonized" if scoped else
                                          "requested" if want is not None else "all")}
        if scoped:
            prov["scope_families"] = others
        uncat: list = []
        tgt, map_note, fail = _target(spec, v, map_fn, build_cache)
        if map_note == S.AUTOSOMES_ONLY_DETAIL:
            prov["autosomes_only"] = True
        elif map_note:
            prov["build_note"] = map_note
        recs: dict[str, RawRecord] = {}
        if fail is not None:
            recs = {o.outcome_id: not_readable_record(o) or RawRecord(fail, detail=map_note)
                    for o in queried}
        else:
            # A study the source holds with no readable file is answered from the catalog
            # here, before any reader, so no reader (injected or not) can call it absent.
            unread = {o.outcome_id: rec for o in queried
                      if (rec := not_readable_record(o)) is not None}
            to_read = [o for o in queried if o.outcome_id not in unread]
            reader = (readers or {}).get(dsid) or _reader_for(
                spec, finngen_route=finngen_route, workers=workers, index_cache=index_cache)
            if spec.family == "finngen":
                recs, rprov, uncat = reader.read(to_read, tgt, cat_ids)
            else:
                recs, rprov = reader.read(to_read, tgt)
            recs.update(unread)
            prov.update(rprov)
            if unread:
                prov["n_not_readable"] = len(unread)
        lost = {o.outcome_id for o in queried} - set(recs)
        if lost:
            raise S.PartitionError(f"{dsid}: the reader returned no record for "
                                   f"{len(lost)} queried outcomes (e.g. {sorted(lost)[:3]})")
        ds_rows = [_row(spec, o, v, tgt, recs.get(o.outcome_id) or _not_looked_up(o))
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
            ds_uncat.append(r)
        prov["status_counts"] = dict(Counter(r.status for r in ds_rows))
        prov["n_uncatalogued"] = len(uncat)
        return ds_rows, ds_uncat, prov, _dataset_notes(spec, ds_rows, prov, len(uncat))

    # The datasets are independent reads, so they run at once; results merge in dataset
    # order. Measured on rs12740374 (three-biobank matched read, empty index cache, 16
    # workers each): 159 s read at once, 308 s one after another (104 + 47 + 153 s).
    from concurrent.futures import ThreadPoolExecutor

    with ThreadPoolExecutor(max_workers=max(1, len(ids))) as ex:
        per_dataset = list(ex.map(one_dataset, ids))
    for dsid, (ds_rows, ds_uncat, prov, ds_notes) in zip(ids, per_dataset):
        rows.extend(ds_rows)
        uncat_rows.extend(ds_uncat)
        provenance[dsid] = prov
        notes.extend(ds_notes)

    others = not_served()
    if others:
        notes.append(" ".join(others.values()))
    result = SnpPhewasResult(variant=v, dataset_ids=ids, rows=rows, uncatalogued=uncat_rows,
                             provenance=provenance, notes=notes, scope=scope)
    if harmonize:
        from snp_phewas_harmonize import attach_harmonization

        attach_harmonization(result)
    return result


def not_readable_sentence(label: str, rows: list[S.OutcomeRow]) -> str:
    """How many of a dataset's outcomes could not be read because the source holds no
    position-indexed file, split by reason (counted from the rows, never recalled)."""
    from snp_phewas_mvp_catalog import ACCESS_NO_SUMSTATS, ACCESS_RAW_ONLY

    nr = [r for r in rows if r.status == S.ST_NOT_READABLE]
    by = Counter((r.status_detail or "").split(":", 1)[0] for r in nr)
    readable = sum(1 for r in rows if r.status != S.ST_NOT_READABLE)
    parts = []
    if by[ACCESS_RAW_ONLY]:
        parts.append(f"{by[ACCESS_RAW_ONLY]:,} not harmonized by the GWAS Catalog")
    if by[ACCESS_NO_SUMSTATS]:
        parts.append(f"{by[ACCESS_NO_SUMSTATS]:,} deposited without summary statistics")
    other = len(nr) - by[ACCESS_RAW_ONLY] - by[ACCESS_NO_SUMSTATS]
    if other:
        parts.append(f"{other:,} for another reason")
    return (f"{label}: {readable:,} of {len(rows):,} outcomes are readable; the other "
            f"{len(nr):,} cannot be read from the GWAS Catalog ({'; '.join(parts)}); they "
            f"are listed, not counted as nulls, and are outside the significance threshold's "
            f"denominator.")


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
    for reason, k in sorted(Counter(r.direction_withheld_reason for r in rows
                                    if r.direction_withheld_reason).items()):
        out.append(f"{spec.label}: no direction shown for {_n(k, 'outcome', 'outcomes')}: "
                   + S.DIRECTION_WITHHELD_LABELS[reason].split(": ", 1)[1] + ".")
    if c[S.ST_NOT_READABLE]:
        out.append(not_readable_sentence(spec.label, rows))
    if prov.get("autosomes_only"):
        out.append(f"{spec.label}: released autosomal results only, so this variant is "
                   f"absent from every one of its outcomes.")
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
    if prov.get("scope") == "harmonized":
        out.append(matched_scope_sentence(spec.label, rows, prov))
    elif prov.get("n_queried") is not None and prov["n_queried"] < prov["n_catalog"]:
        out.append(f"{spec.label}: only {prov['n_queried']:,} of {prov['n_catalog']:,} outcomes "
                   f"were requested; the rest are listed as not part of this lookup.")
    if n_uncat:
        out.append(f"{spec.label}: the source also returned {n_uncat} outcome(s) its published "
                   f"catalog does not list; they are kept separately, not in the counts above.")
    if spec.build != "GRCh38" and prov.get("build_note"):
        out.append(f"{spec.label} uses {spec.build}: {prov['build_note']}.")
    return out
