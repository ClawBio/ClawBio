"""Pre-fetched-input adapters for locuscompare-region-render.

Supports the `sumstats_path` / `ld.source: synthetic` / `gene_track.source: synthetic`
config blocks documented in INPUT_SCHEMA.md. Adapts user-supplied TSVs (or any
sibling skill's TSV output via the `examples/recipes/` harmonisation scripts)
into the duck-typed client surface the core orchestrator expects, so no
fetcher / LD / GENCODE network call is made.

The TSV column contract is the one in INPUT_SCHEMA.md (required columns:
variant_id, chromosome, position_bp, allele_a, allele_b, beta, se, p; optional:
maf, eaf, n, rsid, info, molecular_trait_id, study_id, ...).
"""
from __future__ import annotations

import gzip
import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Iterable

# Sibling-skill dataclass imports. These are added to sys.path by cli.py's
# bootstrap block before this module is imported, so a flat `from <module>`
# resolves correctly inside the fork's flat skill layout.
from eqtl_catalogue_region_fetch import (
    EQTLCatalogueDuplicateRowConflict,
    EQTLCatalogueRelease,
    collapse_duplicate_rsid_rows,
    RegionResult as EQTLRegionResult,
    RegionVariant as EQTLRegionVariant,
)
from gwas_catalog_region_fetch import (
    GWASCatalogRelease,
    RegionResult as GWASRegionResult,
    RegionVariant as GWASRegionVariant,
)


# ---------- Canonical INPUT_SCHEMA.md TSV parsing ----------

REQUIRED_COLUMNS = (
    "variant_id", "chromosome", "position_bp",
    "allele_a", "allele_b", "beta", "se", "p",
)


class PrefetchedSchemaError(ValueError):
    """Raised when a pre-fetched TSV does not satisfy INPUT_SCHEMA.md."""


def _open_text(path: Path):
    """Open a TSV transparently for .tsv / .tsv.gz / .tsv.bgz."""
    suffix = "".join(path.suffixes[-2:]).lower()
    if suffix in (".tsv.gz", ".tsv.bgz") or path.suffix.lower() in (".gz", ".bgz"):
        return gzip.open(path, mode="rt", encoding="utf-8")
    return path.open("rt", encoding="utf-8")


# Cell values read as missing. NaN counts: a NaN never equals itself, so two
# identical copies carrying one would otherwise read as disagreeing.
_MISSING_TOKENS = ("NA", "NAN")


def _parse_float_or_none(s: str) -> float | None:
    s = s.strip()
    if s == "" or s.upper() in _MISSING_TOKENS:
        return None
    return float(s)


def _parse_int_or_none(s: str) -> int | None:
    s = s.strip()
    if s == "" or s.upper() in _MISSING_TOKENS:
        return None
    return int(s)


def load_sumstats_tsv(path: Path) -> list[EQTLRegionVariant]:
    """Load a pre-fetched sumstats slice into the canonical RegionVariant list.

    Returns eQTL-Catalogue-flavoured RegionVariant objects (shape-compatible
    with the GWAS Catalog and UKB-PPP variants — the renderer joins on
    variant_id and consumes only the chrom/pos/ref/alt/beta/se/p columns).

    Per INPUT_SCHEMA.md `allele_a` is the non-effect allele and `allele_b` is
    the effect allele (`beta` per copy of `allele_b`). We surface them as
    `ref` / `alt` on RegionVariant so the downstream harmoniser sees the
    same OT GRCh38 ALT-effect convention the live fetchers emit.

    Raises PrefetchedSchemaError on missing required columns, malformed
    rows, or out-of-range values.
    """
    with _open_text(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if not ln.startswith("#") and ln.strip()]
    if not lines:
        raise PrefetchedSchemaError(f"{path}: empty file")
    header = lines[0].split("\t")
    header_idx = {col: i for i, col in enumerate(header)}
    missing = [c for c in REQUIRED_COLUMNS if c not in header_idx]
    if missing:
        raise PrefetchedSchemaError(
            f"{path}: missing required columns {missing} "
            f"(see INPUT_SCHEMA.md for the canonical TSV schema)"
        )

    variants: list[EQTLRegionVariant] = []
    for lineno, ln in enumerate(lines[1:], start=2):
        cols = ln.split("\t")
        if len(cols) < len(header):
            raise PrefetchedSchemaError(
                f"{path}:{lineno}: expected {len(header)} columns, got {len(cols)}"
            )
        try:
            variant_id = cols[header_idx["variant_id"]].strip()
            chrom = cols[header_idx["chromosome"]].strip()
            pos = int(cols[header_idx["position_bp"]])
            allele_a = cols[header_idx["allele_a"]].strip()
            allele_b = cols[header_idx["allele_b"]].strip()
            beta = _parse_float_or_none(cols[header_idx["beta"]])
            se = _parse_float_or_none(cols[header_idx["se"]])
            p = _parse_float_or_none(cols[header_idx["p"]])
        except ValueError as e:
            raise PrefetchedSchemaError(f"{path}:{lineno}: parse error: {e}") from e

        # Drop rows missing required beta/se/p per INPUT_SCHEMA.md §Missing data
        if beta is None or se is None or p is None:
            continue

        maf = (
            _parse_float_or_none(cols[header_idx["maf"]])
            if "maf" in header_idx else None
        )
        eaf = (
            _parse_float_or_none(cols[header_idx["eaf"]])
            if "eaf" in header_idx else None
        )

        variants.append(EQTLRegionVariant(
            variant_id=variant_id,
            chromosome=chrom,
            position=pos,
            ref=allele_a,
            alt=allele_b,
            beta=beta,
            se=se,
            p_value=p,
            maf=maf,
            effect_allele_frequency=eaf,
            raw={},
        ))
        # The optional `rsid` and `molecular_trait_id` columns ride on `raw`: a
        # row collapsed onto another copy of the same variant keeps its rsID in
        # `raw["rsids"]`, and the exposure client checks the file holds one trait.
        for col in ("rsid", "molecular_trait_id"):
            value = cols[header_idx[col]].strip() if col in header_idx else ""
            if value and value.upper() not in _MISSING_TOKENS:
                variants[-1].raw[col] = value
    return variants


def load_synthetic_ld(path: Path) -> dict[str, float]:
    """Load a two-column synthetic LD matrix TSV (`partner_variant_id`, `r2`).

    Returns a dict mapping partner variant_id to r2. The lead's own r2 entry
    is not required (downstream code seeds it to 1.0).
    """
    with _open_text(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if not ln.startswith("#") and ln.strip()]
    if not lines:
        raise PrefetchedSchemaError(f"{path}: empty file")
    header = lines[0].split("\t")
    header_idx = {col: i for i, col in enumerate(header)}
    for required in ("partner_variant_id", "r2"):
        if required not in header_idx:
            raise PrefetchedSchemaError(
                f"{path}: synthetic LD matrix missing required column {required!r}; "
                f"expected columns: partner_variant_id, r2"
            )
    out: dict[str, float] = {}
    for lineno, ln in enumerate(lines[1:], start=2):
        cols = ln.split("\t")
        if len(cols) < len(header):
            raise PrefetchedSchemaError(
                f"{path}:{lineno}: expected {len(header)} columns, got {len(cols)}"
            )
        try:
            partner = cols[header_idx["partner_variant_id"]].strip()
            r2 = float(cols[header_idx["r2"]])
        except ValueError as e:
            raise PrefetchedSchemaError(f"{path}:{lineno}: parse error: {e}") from e
        if not (0.0 <= r2 <= 1.0 + 1e-9):
            raise PrefetchedSchemaError(
                f"{path}:{lineno}: r2={r2} out of [0, 1]"
            )
        out[partner] = r2
    return out


def load_synthetic_gene_track(path: Path):
    """Load a synthetic gene track TSV. Required columns: gene_symbol, start,
    end, strand. Optional: biotype (default 'protein_coding').

    Returns a list of GeneTrackEntry (imported lazily to avoid pulling in
    matplotlib until render time).
    """
    from regional_plot import GeneTrackEntry

    with _open_text(path) as fh:
        lines = [ln.rstrip("\n") for ln in fh if not ln.startswith("#") and ln.strip()]
    if not lines:
        raise PrefetchedSchemaError(f"{path}: empty file")
    header = lines[0].split("\t")
    header_idx = {col: i for i, col in enumerate(header)}
    for required in ("gene_symbol", "start", "end", "strand"):
        if required not in header_idx:
            raise PrefetchedSchemaError(
                f"{path}: synthetic gene track missing required column {required!r}; "
                f"expected: gene_symbol, start, end, strand (biotype optional)"
            )
    entries = []
    for lineno, ln in enumerate(lines[1:], start=2):
        cols = ln.split("\t")
        try:
            entries.append(GeneTrackEntry(
                gene_symbol=cols[header_idx["gene_symbol"]].strip(),
                start=int(cols[header_idx["start"]]),
                end=int(cols[header_idx["end"]]),
                strand=cols[header_idx["strand"]].strip() or "+",
                exons=[],
                biotype=(
                    cols[header_idx["biotype"]].strip()
                    if "biotype" in header_idx else "protein_coding"
                ),
            ))
        except (ValueError, IndexError) as e:
            raise PrefetchedSchemaError(f"{path}:{lineno}: parse error: {e}") from e
    return entries


# ---------- Duck-typed stub clients ----------

def _now_utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


@dataclass
class PrefetchedEQTLClient:
    """Stub eQTL Catalogue client backed by pre-loaded variants. Implements
    the `fetch_region(...)` surface the core orchestrator calls; returns an
    EQTLRegionResult with a synthetic release block tagged 'prefetched'."""

    variants: list[EQTLRegionVariant]
    dataset_id: str = "prefetched"
    dataset_release: str = "prefetched"
    study_label: str | None = None
    tissue_label: str | None = None
    quant_method: str = "ge"
    notes: list[str] = field(default_factory=list)

    def fetch_region(
        self,
        *,
        dataset_id: str,
        chromosome: str,
        start_bp: int,
        end_bp: int,
        gene_id: str | None = None,
    ) -> EQTLRegionResult:
        # The plot draws ONE molecular trait. A file that mixes several
        # (transcripts or splicing events of one gene) would pair a variant once
        # per trait, so it stops here and names them.
        traits = sorted({v.raw["molecular_trait_id"] for v in self.variants
                         if v.raw.get("molecular_trait_id")})
        if len(traits) > 1:
            raise PrefetchedSchemaError(
                f"pre-fetched exposure file holds {len(traits)} molecular traits "
                f"(e.g. {', '.join(traits[:3])}); the plot draws one trait, so "
                f"filter the file to one molecular_trait_id first"
            )
        # One row per variant, as the live fetch returns. A slice cut from the
        # eQTL Catalogue carries one row per rsID, so a variant with several
        # rsIDs repeats with identical statistics: those copies collapse to the
        # first, and each dropped copy is kept on the result's `collapsed`.
        # Copies that DISAGREE cannot be merged, and picking one would be a
        # guess, so the render stops.
        try:
            collapse = collapse_duplicate_rsid_rows(
                self.variants, dataset_id=self.dataset_id or dataset_id,
            )
        except EQTLCatalogueDuplicateRowConflict as e:
            detail = str(e).split(" disagree on ", 1)[-1].split(";", 1)[0]
            vids = [v.variant_id for v in self.variants]
            repeated = sorted({vid for vid in vids if vids.count(vid) > 1})
            raise PrefetchedSchemaError(
                f"pre-fetched exposure file has rows of one variant that disagree "
                f"on {detail} (variants with repeated rows include "
                f"{', '.join(repeated[:3])}); repeated rows must be identical "
                f"copies, so no row was picked"
            ) from e
        variants = collapse.kept
        notes = list(self.notes) + ["sumstats loaded from pre-fetched TSV (no live fetch)"]
        if collapse.collapsed:
            notes.append(
                f"{len(collapse.collapsed)} repeated row(s) of a variant with identical "
                f"statistics (one row per rsID) collapsed to one row each"
            )
        return EQTLRegionResult(
            dataset_id=self.dataset_id or dataset_id,
            chromosome=chromosome,
            region_start_bp=start_bp,
            region_end_bp=end_bp,
            n_variants=len(variants),
            variants=variants,
            n_rows_read=len(self.variants),
            collapsed=collapse.collapsed,
            release=EQTLCatalogueRelease(
                api_version="prefetched",
                dataset_release=self.dataset_release,
                fetched_at_utc=_now_utc(),
                study_label=self.study_label,
                tissue_label=self.tissue_label,
                quant_method=self.quant_method,
            ),
            notes=notes,
        )


@dataclass
class PrefetchedGWASClient:
    """Stub GWAS Catalog client backed by pre-loaded variants."""

    variants: list[GWASRegionVariant]
    accession: str = "prefetched"
    notes: list[str] = field(default_factory=list)

    def fetch_region(
        self,
        *,
        accession: str,
        chromosome: str,
        start_bp: int,
        end_bp: int,
    ) -> GWASRegionResult:
        return GWASRegionResult(
            accession=self.accession or accession,
            chromosome=chromosome,
            region_start_bp=start_bp,
            region_end_bp=end_bp,
            n_variants=len(self.variants),
            variants=self.variants,
            release=GWASCatalogRelease(
                accession=self.accession or accession,
                harmonised_url="prefetched",
                fetched_at_utc=_now_utc(),
            ),
            notes=list(self.notes) + ["sumstats loaded from pre-fetched TSV (no live fetch)"],
        )


@dataclass
class PrefetchedLDClient:
    """Stub LD client backed by a pre-loaded r2 dict. The core orchestrator
    only consumes `pairs`, `plink_version`, `panel_id`, `panel_version` from
    the LD result; those are reported here as synthetic markers."""

    r2_by_partner: dict[str, float]
    panel_id: str = "synthetic"
    panel_version: str = "1.0"
    plink_version: str = "prefetched"
    super_pop: str = "synthetic"

    def r2_with_lead(
        self,
        lead: str,
        partners: Iterable[str],
        chromosome: str,
        window_bp: int,
    ):
        from ondemand_client import OnDemandLDPair, OnDemandLDResult

        partner_list = list(partners)
        pairs = []
        for p in partner_list:
            if p in self.r2_by_partner:
                pairs.append(OnDemandLDPair(partner_variant_id=p, r2=self.r2_by_partner[p]))
        return OnDemandLDResult(
            panel_id=self.panel_id,
            panel_version=self.panel_version,
            super_pop=self.super_pop,
            plink_version=self.plink_version,
            chromosome=chromosome,
            lead_variant_id=lead,
            window_bp=window_bp,
            n_partners_requested=len(partner_list),
            n_partners_returned=len(pairs),
            pairs=pairs,
            fetched_at_utc=_now_utc(),
            notes=["LD r² loaded from pre-fetched synthetic matrix (no plink call)"],
        )


def gwas_variants_from_eqtl(variants: list[EQTLRegionVariant]) -> list[GWASRegionVariant]:
    """Re-shape EQTL-flavoured RegionVariant rows into GWAS-flavoured ones.

    The outcome side of the orchestrator expects GWAS Catalog's RegionVariant
    (one extra optional `odds_ratio` field; otherwise identical). We load
    both sides via the canonical INPUT_SCHEMA TSV loader (which yields eQTL-
    flavoured rows) and reshape the outcome list here.
    """
    out: list[GWASRegionVariant] = []
    for v in variants:
        out.append(GWASRegionVariant(
            variant_id=v.variant_id,
            chromosome=v.chromosome,
            position=v.position,
            ref=v.ref,
            alt=v.alt,
            beta=v.beta,
            se=v.se,
            p_value=v.p_value,
            odds_ratio=None,
            effect_allele_frequency=v.effect_allele_frequency,
            raw={},
        ))
    return out
