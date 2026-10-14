"""Tier-2 LocusCompare orchestrator.

Wires the three execution skills (eQTL Catalogue, GWAS Catalog harmonised, LD
reference via plink + 1000G) + the wald_ratio harmoniser + the
render_full_locuscompare renderer + the manifest block builder, all driven by
a small per-row StudyIdMapping that resolves OT studyIds to the
upstream-source identifiers each fetcher needs.

Inputs:
- lead variant id (chr_pos_ref_alt, GRCh38)
- chromosome + lead_position_bp + window_bp
- StudyIdMapping (one per Tier-1 coloc row that gates Tier-2 rendering)
- pre-built clients (so callers can mock for tests)

Outputs:
- a PNG at the requested out_path (the 4-panel render_full_locuscompare figure)
- a `regional_locuscompare` manifest block documenting the four input
  panels, the lead variant, and the LD reference panel.

If `ld_client` is None (e.g. plink not installed in the current environment),
the renderer still produces all four panels but with r² coloring substituted
by a uniform grey; the manifest block records `ld_panel: "none"` and a caveat
in the caveats list.

When the eQTL Catalogue or GWAS Catalog source can't resolve the studyId, the
orchestrator raises Tier2NotAvailable so the caller can route to a credible-
set-only fallback.
"""

from __future__ import annotations

import time
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import sys
from pathlib import Path

import requests
import yaml

# Resolve sibling-skill imports for the fork's flat layout. This skill
# orchestrates eqtl-catalogue-region-fetch, gwas-catalog-region-fetch,
# and ld-1000g-region-compute, which live as sibling skill directories
# under skills/ in the same checkout. Inject each sibling-skill dir
# onto sys.path so the modules import cleanly.
_SKILL_DIR = Path(__file__).resolve().parent
_SKILLS_ROOT = _SKILL_DIR.parent
for _p in (
    _SKILL_DIR,
    _SKILLS_ROOT / "eqtl-catalogue-region-fetch",
    _SKILLS_ROOT / "gwas-catalog-region-fetch",
    _SKILLS_ROOT / "ld-1000g-region-compute",
    _SKILLS_ROOT / "ukb-ppp-region-fetch",
):
    if _p.exists() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

# Sibling-skill imports (one per upstream primitive in the suite).
from eqtl_catalogue_region_fetch import (
    EQTLCatalogueAPIError,
    EQTLCatalogueClient,
    RegionResult as EQTLRegionResult,
)
from gwas_catalog_region_fetch import (
    GWASCatalogClient,
    GWASCatalogFetchError,
    RegionResult as GWASRegionResult,
)
from ld_1000g_region_compute import (
    LDComputeError,
    SuperPop,
)
from ondemand_client import (
    OnDemand1000GLDClient,
)
from ukb_ppp_region_fetch import (
    RegionResult as UKBPPPRegionResult,
    UKBPPPAccessError,
    UKBPPPClient,
)
# The odds ratio -> log-odds conversion and the beta-source labels. The module is kept
# byte-identical wherever it is copied, so the conversion gives the same numbers in
# every copy; edit it only by re-syncing the whole file. It imports only the stdlib.
from _region_harmonise import (
    BETA_SOURCE_MIXED,
    BETA_SOURCE_NATIVE,
    BETA_SOURCE_OR_DERIVED,
    OR_NO_SE_REASON,
    beta_source_label,
    derive_beta_from_or,
    outcome_effect_scale_label,
)
from regional_plot import (
    GeneTrackEntry,
    LocusVariant,
    RegionalLocusCompareInput,
    harmonise_regions_for_locuscompare,
    render_full_locuscompare,
)

# GENCODE gene-track parser: bundled embedded helper, lives next to this
# script under _fetchers/. The on-demand Ensembl-REST fetcher avoids the
# multi-GB local-GTF dependency.
from _fetchers.gencode_ondemand import fetch_region_genes_remote


# Exposure kind dispatch (v1.3). The eQTL Catalogue fetcher serves
# eQTL + sQTL + sceQTL via the same code path; UKB-PPP serves pQTL.
# Unknown prefixes raise Tier2NotAvailable so the caller falls back to
# Tier-1 (credible-set-only) with the documented caption flag.
EXPOSURE_KIND_EQTL_CATALOGUE = "eqtl_catalogue"
EXPOSURE_KIND_UKB_PPP = "ukb_ppp"


def is_pqtl_study_id(ot_study_id: str) -> bool:
    """Return True iff an OT QTL studyId references a UKB-PPP pQTL study.

    UKB-PPP pQTL studyIds follow the `UKB_PPP_<ancestry>_<protein>` pattern.
    The prefix match is case-insensitive against the canonical OT export.
    """
    return ot_study_id.upper().startswith("UKB_PPP_")


def dispatch_exposure_kind(ot_study_id: str) -> str:
    """Map an OT QTL studyId to the exposure-fetcher kind.

    Returns one of EXPOSURE_KIND_EQTL_CATALOGUE or EXPOSURE_KIND_UKB_PPP.
    The default is eQTL Catalogue because that fetcher serves all four
    eQTL-Catalogue-hosted quant methods (ge / exon / tx / txrev /
    leafcutter) plus the sceQTL studies in v7+.
    """
    if is_pqtl_study_id(ot_study_id):
        return EXPOSURE_KIND_UKB_PPP
    return EXPOSURE_KIND_EQTL_CATALOGUE


@dataclass
class StudyIdMapping:
    """Per-row resolution of OT studyIds to upstream-source ids.

    The OT studyId of the QTL credible set (e.g.
    `quach_2016_ge_monocyte_iav_ensg00000115808`) maps to an eQTL Catalogue
    `dataset_id` (e.g. `QTD000110`). For UKB-PPP pQTL rows the studyId
    encodes the protein + ancestry (e.g. `UKB_PPP_EUR_SORT1`); the mapping
    additionally captures the canonical HGNC protein label + ancestry code
    for the fetcher's name -> Synapse-fileID resolver. The OT outcome
    studyId (e.g. `FINNGEN_R12_I9_HEARTFAIL`) maps to a GWAS Catalog
    `accession` (e.g. `GCST90475990`).

    The canonical mapping table is auto-resolvable for most rows; manual
    lookup is acceptable as a starting point. This dataclass is the row-
    level contract; the lookup mechanism (auto or manual) is the
    orchestrator's call.

    `outcome_trait_label` is human-readable text for the panel title (e.g.
    "hypertrophic cardiomyopathy"). The exposure side's tissue / condition /
    quant labels come from the eQTL Catalogue dataset metadata at fetch time
    so they don't need to be repeated in the YAML lookup. For pQTL rows
    the Olink panel + ancestry-label are read from UKB-PPP at fetch time
    too; only protein_label + ancestry are needed in the YAML.

    Per v1.3 dispatch, exactly one of `eqtl_catalogue_dataset_id` or
    (`ukb_ppp_protein_label`, `ukb_ppp_ancestry`) must be populated for
    the row to resolve. The contract is validated at fetch time, not at
    dataclass construction, so YAML round-trips don't break for in-
    progress rows.
    """
    ot_left_study_id: str
    ot_right_study_id: str
    gwas_catalog_accession: str
    eqtl_catalogue_dataset_id: str = ""
    ukb_ppp_protein_label: str = ""
    ukb_ppp_ancestry: str = ""
    ancestry_left: str = "EUR"
    ancestry_right: str = "EUR"
    outcome_trait_label: str = ""
    exposure_gene_symbol: str = ""
    notes: list[str] = field(default_factory=list)


class Tier2NotAvailable(Exception):
    """Raised when one of the upstream sources cannot resolve the requested
    studyId or region. Caller should fall back to a credible-set-only
    rendering path.
    """


@dataclass
class Tier2Result:
    plot_path: Path
    manifest_block: dict[str, Any]
    n_pairs: int
    n_palindromic_excluded: int
    notes: list[str]


@dataclass
class LocusCompareSpec:
    """Generic input to the locuscompare core. Decoupled from OT.

    Built either directly (non-OT entry vectors: gwas-lookup follow-up,
    fine-mapping chain, raw user-supplied harmonised TSVs) or by the OT-shim
    wrapper `render_tier2_for_lead` which translates a `StudyIdMapping` into
    this generic spec.

    `exposure_kind` selects the fetcher backend: EXPOSURE_KIND_EQTL_CATALOGUE
    consumes `eqtl_dataset_id` + `molecular_trait_id`; EXPOSURE_KIND_UKB_PPP
    consumes `pqtl_protein_label` + `pqtl_ancestry`. Adding a new molQTL
    backend means adding a new EXPOSURE_KIND_* + the corresponding fields
    here + the dispatch branch in `_render_for_spec`.
    """
    lead_variant_id: str
    chromosome: str
    lead_position_bp: int
    window_bp: int

    eqtl_dataset_id: str
    molecular_trait_id: str | None

    gwas_accession: str

    exposure_kind: str = EXPOSURE_KIND_EQTL_CATALOGUE
    pqtl_protein_label: str = ""
    pqtl_ancestry: str = ""

    exposure_gene_symbol: str = ""
    outcome_trait_label: str = ""

    # Optional rs-ID for the lead variant. Set by upstream agents
    # (e.g. ai_scientist's coloc_with_mr workflow resolves via OT / dbSNP);
    # propagates into the manifest + report only. The orchestrator joins on
    # variant_id (chr_pos_ref_alt), not rs_id; this field is human-readability
    # metadata only.
    lead_rs_id: str | None = None

    exposure_id_extra: str = ""
    outcome_id_extra: str = ""
    provenance_prefix: str = ""
    # The outcome's upstream study id when it differs from `gwas_accession` (e.g. an
    # OT `FINNGEN_R12_*` id mapped to a GWAS Catalog accession). Used only to name the
    # outcome effect's scale: FinnGen disease endpoints report log odds.
    outcome_source_study_id: str = ""

    release_tag: str = ""

    notes: list[str] = field(default_factory=list)

    super_pop: SuperPop = SuperPop.EUR
    intersected_pip_product: float | None = None
    extra_caveats: list[str] = field(default_factory=list)
    gencode_gtf_path: Path | None = None
    gene_biotypes: tuple[str, ...] | None = ("protein_coding",)
    # Optional pre-fetched gene track. When provided, bypasses the local-GTF
    # path entirely. Use this to inject genes resolved from a remote source
    # (e.g. locuscompare's Ensembl REST on-demand fetcher) without writing a
    # synthetic GTF to disk.
    prefetched_gene_track: list | None = None


def render_locuscompare_for_lead(
    spec: LocusCompareSpec,
    *,
    eqtl_client: EQTLCatalogueClient,
    gwas_client: GWASCatalogClient,
    ld_client: OnDemand1000GLDClient | None,
    out_path: Path,
    ukb_ppp_client: UKBPPPClient | None = None,
) -> Tier2Result:
    """Run the full locuscompare pipeline for one lead variant. Generic, no OT.

    Steps:
    1. Fetch exposure region (eQTL Catalogue or UKB-PPP, per `spec.exposure_kind`).
    2. Fetch outcome region from GWAS Catalog harmonised.
    3. Compute r² between the lead and every harmonised partner via plink.
    4. Join + harmonise per skills.knowledge.wald_ratio.harmonise_regions_for_locuscompare.
    5. Render the 4-panel figure.
    6. Build the manifest block.
    """
    return _render_for_spec(
        spec=spec,
        eqtl_client=eqtl_client,
        gwas_client=gwas_client,
        ld_client=ld_client,
        out_path=out_path,
        ukb_ppp_client=ukb_ppp_client,
    )


def render_tier2_for_lead(
    *,
    lead_variant_id: str,
    chromosome: str,
    lead_position_bp: int,
    window_bp: int,
    study_mapping: StudyIdMapping,
    eqtl_client: EQTLCatalogueClient,
    gwas_client: GWASCatalogClient,
    ld_client: OnDemand1000GLDClient | None,
    out_path: Path,
    ot_release: str,
    super_pop: SuperPop = SuperPop.EUR,
    intersected_pip_product: float | None = None,
    extra_caveats: list[str] | None = None,
    gencode_gtf_path: Path | None = None,
    gene_biotypes: tuple[str, ...] | None = ("protein_coding",),
    ukb_ppp_client: UKBPPPClient | None = None,
) -> Tier2Result:
    """OT-shim wrapper. Resolves an OT row into a generic LocusCompareSpec
    and delegates to `render_locuscompare_for_lead`.

    Existing callers (tier2_cli, tests, downstream orchestrators) continue
    to work unchanged. The OT-specific work is concentrated here:
    - Dispatch by OT studyId prefix: `UKB_PPP_*` -> UKB-PPP pQTL fetcher;
      everything else -> eQTL Catalogue fetcher (the existing
      ge/exon/tx/txrev/leafcutter/sceQTL path).
    - Extract ENSG from OT QTL studyId for the eQTL Cat molecular_trait_id.
    - Trigger the FinnGen ancestry caveat.
    - Format OT-flavored label / provenance strings.

    The dispatched exposure-fetcher must be passed in (`eqtl_client` for
    eQTL/sQTL/sceQTL rows, `ukb_ppp_client` for pQTL rows). When dispatch
    selects a backend whose client was not supplied, raises
    `Tier2NotAvailable` so the caller falls back to Tier-1 with the
    documented caption flag.
    """
    runtime_caveats = list(extra_caveats or [])
    if "finngen" in study_mapping.ot_right_study_id.lower():
        runtime_caveats.append(
            "FinnGen Finnish-EUR; 1000G EUR proxy used. "
            "Common-variant LD agrees within ~0.05 r² per Locke 2019."
        )

    exposure_kind = dispatch_exposure_kind(study_mapping.ot_left_study_id)

    if exposure_kind == EXPOSURE_KIND_UKB_PPP:
        if ukb_ppp_client is None:
            raise Tier2NotAvailable(
                f"OT studyId {study_mapping.ot_left_study_id!r} is a UKB-PPP "
                f"pQTL row but no ukb_ppp_client was supplied; fall back to "
                f"Tier-1 (CS-only) and caption 'UKB-PPP client not "
                f"configured for this render'."
            )
        if not study_mapping.ukb_ppp_protein_label:
            raise Tier2NotAvailable(
                f"OT studyId {study_mapping.ot_left_study_id!r} is a UKB-PPP "
                f"pQTL row but `ukb_ppp_protein_label` is empty in the "
                f"mapping; fix the YAML row or fall back to Tier-1."
            )
        spec = LocusCompareSpec(
            lead_variant_id=lead_variant_id,
            chromosome=chromosome,
            lead_position_bp=lead_position_bp,
            window_bp=window_bp,
            # eqtl_dataset_id + molecular_trait_id unused by the pQTL path
            # but the dataclass requires them; pass safe defaults.
            eqtl_dataset_id="",
            molecular_trait_id=None,
            gwas_accession=study_mapping.gwas_catalog_accession,
            exposure_kind=EXPOSURE_KIND_UKB_PPP,
            pqtl_protein_label=study_mapping.ukb_ppp_protein_label,
            pqtl_ancestry=study_mapping.ukb_ppp_ancestry or "EUR",
            exposure_gene_symbol=study_mapping.exposure_gene_symbol or
            study_mapping.ukb_ppp_protein_label,
            outcome_trait_label=study_mapping.outcome_trait_label,
            exposure_id_extra=f" (OT studyId {study_mapping.ot_left_study_id})",
            outcome_id_extra=f" (OT studyId {study_mapping.ot_right_study_id})",
            outcome_source_study_id=study_mapping.ot_right_study_id,
            provenance_prefix=f"OT release: {ot_release} | ",
            release_tag=ot_release,
            notes=list(study_mapping.notes),
            super_pop=super_pop,
            intersected_pip_product=intersected_pip_product,
            extra_caveats=runtime_caveats,
            gencode_gtf_path=gencode_gtf_path,
            gene_biotypes=gene_biotypes,
        )
    else:
        molecular_trait_id = _extract_ensg_from_ot_study_id(
            study_mapping.ot_left_study_id,
        )
        spec = LocusCompareSpec(
            lead_variant_id=lead_variant_id,
            chromosome=chromosome,
            lead_position_bp=lead_position_bp,
            window_bp=window_bp,
            eqtl_dataset_id=study_mapping.eqtl_catalogue_dataset_id,
            molecular_trait_id=molecular_trait_id,
            gwas_accession=study_mapping.gwas_catalog_accession,
            exposure_kind=EXPOSURE_KIND_EQTL_CATALOGUE,
            exposure_gene_symbol=study_mapping.exposure_gene_symbol,
            outcome_trait_label=study_mapping.outcome_trait_label,
            exposure_id_extra=f" (OT studyId {study_mapping.ot_left_study_id})",
            outcome_id_extra=f" (OT studyId {study_mapping.ot_right_study_id})",
            outcome_source_study_id=study_mapping.ot_right_study_id,
            provenance_prefix=f"OT release: {ot_release} | ",
            release_tag=ot_release,
            notes=list(study_mapping.notes),
            super_pop=super_pop,
            intersected_pip_product=intersected_pip_product,
            extra_caveats=runtime_caveats,
            gencode_gtf_path=gencode_gtf_path,
            gene_biotypes=gene_biotypes,
        )

    return _render_for_spec(
        spec=spec,
        eqtl_client=eqtl_client,
        gwas_client=gwas_client,
        ld_client=ld_client,
        out_path=out_path,
        ukb_ppp_client=ukb_ppp_client,
    )


def _render_for_spec(
    *,
    spec: LocusCompareSpec,
    eqtl_client: EQTLCatalogueClient,
    gwas_client: GWASCatalogClient,
    ld_client: OnDemand1000GLDClient | None,
    out_path: Path,
    ukb_ppp_client: UKBPPPClient | None = None,
) -> Tier2Result:
    """Core implementation shared by both entry points.

    Lifted from the original render_tier2_for_lead body, with all
    `study_mapping.X` references replaced by `spec.X` reads. v1.3 adds
    a UKB-PPP exposure-fetch branch on `spec.exposure_kind`.
    """
    out_path = Path(out_path)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    lead_variant_id = spec.lead_variant_id
    chromosome = spec.chromosome
    window_bp = spec.window_bp
    lead_position_bp = spec.lead_position_bp
    super_pop = spec.super_pop
    intersected_pip_product = spec.intersected_pip_product
    gencode_gtf_path = spec.gencode_gtf_path
    gene_biotypes = spec.gene_biotypes

    half = max(window_bp // 2, 1)
    start_bp = max(0, lead_position_bp - half)
    end_bp = lead_position_bp + half

    notes: list[str] = []
    extra_caveats = list(spec.extra_caveats)
    # v1.4 fetcher notes: each fetcher's RegionResult.notes propagates into a
    # separate manifest field, prefixed with the originating fetcher name.
    # Kept distinct from `extra_caveats` (curated, user-facing) so the noisy
    # data-source channel can be triaged independently.
    data_source_warnings: list[str] = []

    eqtl_dataset_id = spec.eqtl_dataset_id
    gwas_accession = spec.gwas_accession
    molecular_trait_id = spec.molecular_trait_id

    # 1. Exposure region. v1.3 dispatch by spec.exposure_kind. The two
    # backends return shape-compatible RegionResult objects (same
    # variant fields, equivalent release metadata); downstream join +
    # render code is agnostic to which produced the rows. `exposure`
    # plus `exposure_source`, `exposure_source_release`, and
    # `exposure_study_id` capture all the source-specific metadata that
    # the manifest block + caption builder consume.
    pqtl = None
    if spec.exposure_kind == EXPOSURE_KIND_UKB_PPP:
        if ukb_ppp_client is None:
            raise Tier2NotAvailable(
                "spec.exposure_kind == ukb_ppp but no ukb_ppp_client supplied"
            )
        try:
            pqtl = ukb_ppp_client.fetch_region(
                protein_label=spec.pqtl_protein_label,
                ancestry=spec.pqtl_ancestry or "EUR",
                chromosome=chromosome.lstrip("chr"),
                start_bp=start_bp,
                end_bp=end_bp,
            )
        except UKBPPPAccessError as e:
            raise Tier2NotAvailable(
                f"UKB-PPP cannot resolve {spec.pqtl_protein_label!r} "
                f"in ancestry {spec.pqtl_ancestry or 'EUR'}: {e!s}"
            ) from e
        if pqtl.n_variants == 0:
            raise Tier2NotAvailable(
                f"UKB-PPP returned zero variants for "
                f"{spec.pqtl_protein_label} ({spec.pqtl_ancestry or 'EUR'}) "
                f"at {chromosome}:{start_bp}-{end_bp}"
            )
        exposure_variants = pqtl.variants
        exposure_source = "ukb_ppp"
        exposure_source_release = pqtl.release.release_label
        exposure_study_id = pqtl.release.synapse_id
        exposure_protein_label = pqtl.release.protein_label
        exposure_ancestry_code = pqtl.release.ancestry
        exposure_ancestry_label = pqtl.release.ancestry_label
        # eQTL Catalogue path leaves these unset; mark for the
        # _build_short_labels switch.
        exposure_eqtl_release = None
        data_source_warnings.extend(f"ukb_ppp: {n}" for n in pqtl.notes)
    else:
        try:
            # Filter on the `gene_id` column, not `molecular_trait_id`.
            # spec.molecular_trait_id is the parent ENSG (extracted from
            # the OT studyId). For `ge` / `microarray` the molecular_trait_id
            # column equals gene_id, so both work; for splicing / exon /
            # transcript quant methods the molecular_trait_id column is a
            # cluster / exon / transcript id and would never match the ENSG.
            # `gene_id` is the portable filter across all quant methods.
            exposure: EQTLRegionResult = eqtl_client.fetch_region(
                dataset_id=eqtl_dataset_id,
                chromosome=chromosome.lstrip("chr"),
                start_bp=start_bp,
                end_bp=end_bp,
                gene_id=molecular_trait_id,
            )
        except EQTLCatalogueAPIError as e:
            raise Tier2NotAvailable(
                f"eQTL Catalogue cannot resolve dataset {eqtl_dataset_id}: {e!s}"
            ) from e
        if exposure.n_variants == 0:
            raise Tier2NotAvailable(
                f"eQTL Catalogue returned zero variants for "
                f"{eqtl_dataset_id} at {chromosome}:{start_bp}-{end_bp}"
            )
        exposure_variants = exposure.variants
        exposure_source = "eqtl_catalogue"
        exposure_source_release = exposure.release.dataset_release or "v7+"
        exposure_study_id = eqtl_dataset_id
        exposure_protein_label = ""
        exposure_ancestry_code = ""
        exposure_ancestry_label = ""
        exposure_eqtl_release = exposure.release
        data_source_warnings.extend(f"eqtl_catalogue: {n}" for n in exposure.notes)
        # Surface the file-class disclosure when the window came out of
        # .cc.tsv.gz rather than the full nominal pass. Per the catalogue's
        # docs/Columns_parquet.md, .cc keeps only the traits with permutation
        # FDR < 1% and a credible set, each with all its variants; for exon /
        # transcript / splicing data it also keeps one trait per fine-mapped
        # signal, so a trait that was not the top one of its signal is absent
        # from the file. The sentence is scoped by quant method so it stays
        # true if a ge .cc is ever rendered.
        #
        # Gated on the file class the fetch actually served, not on
        # quant_method. quant_method is a best-effort presentation label and is
        # None whenever the metadata it comes from is unavailable, so gating on
        # it skips the disclosure for a genuinely .cc-backed panel. It is also
        # wrong in the other direction: QTD000584 (aptamer) ships .all, so the
        # quant-method rule would state ".cc.tsv.gz" about a panel built from
        # the full nominal pass. The fetcher records the class of the file it
        # opened on the result, so that is the authority here.
        #
        # file_class is None only for a result that never read one (a
        # hand-built test double or a replayed fixture); the quant-method rule
        # is kept as the fallback for that case alone.
        _fc = (exposure.file_class or "").lower()
        _qm = (exposure.release.quant_method or "").lower()
        _credible_set_filtered = (
            _fc == "cc" if _fc else bool(_qm) and _qm not in {"ge", "microarray"}
        )
        if _credible_set_filtered:
            # The quant_method clause is dropped rather than printed empty when
            # the label is unavailable: the disclosure is about the file class.
            _qm_clause = f"; quant_method={_qm}" if _qm else ""
            if _qm in {"ge", "microarray", "aptamer"}:
                _rule = ("keeps only molecular traits with a significant signal "
                         "and a fine-mapped credible set, each with all its variants")
            else:
                _rule = ("keeps one molecular trait per fine-mapped signal, so a "
                         "trait that was not the top one of its signal is absent "
                         "from the file")
            extra_caveats.append(
                f"sumstats are credible-set-filtered (eQTL Catalogue "
                f".cc.tsv.gz{_qm_clause}); {_rule}"
            )

    # 2. Outcome region (GWAS Catalog harmonised).
    try:
        outcome: GWASRegionResult = gwas_client.fetch_region(
            accession=gwas_accession,
            chromosome=chromosome.lstrip("chr"),
            start_bp=start_bp,
            end_bp=end_bp,
        )
    except GWASCatalogFetchError as e:
        raise Tier2NotAvailable(
            f"GWAS Catalog cannot resolve accession {gwas_accession}: {e!s}"
        ) from e
    if outcome.n_variants == 0:
        raise Tier2NotAvailable(
            f"GWAS Catalog returned zero variants for "
            f"{gwas_accession} at {chromosome}:{start_bp}-{end_bp}"
        )
    data_source_warnings.extend(f"gwas_catalog: {n}" for n in outcome.notes)

    # 2b. Odds ratio -> log-odds effect size, on the FETCHED rows, before the allele
    # flip in `harmonise_regions_for_locuscompare`, so a flipped row negates the derived
    # beta exactly as it would a reported one. Without it, a case-control study that
    # publishes odds ratios and no beta drew bare effect-size axes.
    n_native_outcome = sum(1 for v in outcome.variants if v.beta is not None)
    outcome_rows, or_events = derive_beta_from_or(list(outcome.variants))
    if or_events and n_native_outcome:
        # Some rows reported a beta and others only an odds ratio. The derived-log-odds
        # label would be false for the reported points, so the window is labelled mixed
        # and no scale is claimed for it.
        outcome_beta_source: str | None = BETA_SOURCE_MIXED
    elif or_events:
        outcome_beta_source = BETA_SOURCE_OR_DERIVED
    elif n_native_outcome:
        outcome_beta_source = BETA_SOURCE_NATIVE
    else:
        outcome_beta_source = None  # the window carries no outcome effect size at all
    outcome_scale_label = (
        # Not a scale, a warning on the axis: the reported points' scale is unknown.
        "mixed scales: reported effect sizes and derived log odds ratios"
        if outcome_beta_source == BETA_SOURCE_MIXED
        else outcome_effect_scale_label(
            outcome_beta_source, spec.outcome_source_study_id or gwas_accession,
        )
    )
    n_or_unconverted = sum(
        1 for v in outcome_rows if v.beta is None and v.odds_ratio is not None
    )
    if or_events:
        notes.append(
            f"outcome effect sizes are log odds ratios derived from the reported odds "
            f"ratios ({len(or_events)} of {len(outcome_rows)} outcome variants converted)"
        )
    if or_events and n_native_outcome:
        notes.append(
            f"the outcome window mixes effect-size scales: {n_native_outcome} variants "
            f"report an effect size and {len(or_events)} carry a log odds ratio derived "
            f"from a reported odds ratio, so the effect-size axis names no single scale"
        )
    if or_events and n_or_unconverted:
        notes.append(
            f"{n_or_unconverted} outcome variants report an odds ratio from which "
            f"{OR_NO_SE_REASON}, so no effect size was derived for them; they appear on "
            f"the p-value panels only"
        )

    # 3. LD r² (optional; gracefully degrade when plink not installed).
    r2_by_variant: dict[str, float] = {lead_variant_id: 1.0}
    plink_version = ""
    panel_id = "none"
    panel_version = ""
    if ld_client is not None:
        partner_ids = [v.variant_id for v in exposure_variants if v.variant_id != lead_variant_id]
        try:
            ld = ld_client.r2_with_lead(
                lead=lead_variant_id,
                partners=partner_ids,
                chromosome=chromosome,
                window_bp=window_bp,
            )
        except LDComputeError as e:
            notes.append(f"LD computation failed: {e!s}; rendering without LD coloring")
            extra_caveats.append("LD r² unavailable (compute failed); points show as grey")
        else:
            for pair in ld.pairs:
                r2_by_variant[pair.partner_variant_id] = pair.r2
            plink_version = ld.plink_version
            panel_id = ld.panel_id
            panel_version = ld.panel_version
    else:
        notes.append("ld_client is None; rendering without LD coloring")
        extra_caveats.append("LD r² unavailable (no plink client provided); points show as grey")

    # 4. Harmonise + join.
    pairs = harmonise_regions_for_locuscompare(
        exposure_variants=exposure_variants,
        outcome_variants=outcome_rows,
        r2_by_variant=r2_by_variant,
        lead_variant_id=lead_variant_id,
    )
    if not pairs:
        raise Tier2NotAvailable(
            f"no joinable variants between exposure {exposure_study_id} "
            f"and outcome {gwas_accession} at "
            f"{chromosome}:{start_bp}-{end_bp}"
        )

    n_palindromic = sum(1 for p in pairs if p.palindromic_excluded)

    # Never bare axes: when the effect-size panel will have nothing to plot, say which
    # side lacks effect sizes and why, on the panel, in the notes and in the manifest.
    effect_size_reason = _effect_size_unavailable_reason(
        pairs=pairs,
        exposure_rows=exposure_variants,
        outcome_rows=outcome_rows,
        exposure_label=(
            f"UKB-PPP protein {exposure_protein_label}"
            if spec.exposure_kind == EXPOSURE_KIND_UKB_PPP
            else _study_label(exposure_study_id)
        ),
        outcome_label=_study_label(gwas_accession),
    )
    if effect_size_reason:
        notes.append(f"effect-size panel not drawn: {effect_size_reason}")

    # 5. Render.
    fetched_at = _now_utc()
    pip_label = (
        f" (intersected; PIP_L x PIP_R = {intersected_pip_product:.4f})"
        if intersected_pip_product is not None else ""
    )
    # Convert the raw exposure / outcome region rows to wald_ratio LocusVariant
    # so the renderer's per-side manhattan tracks can show ALL variants in the
    # source's region (LocusZoom-style "bottom of plot fills the window"),
    # not just the joined intersection. The eQTL Catalogue and UKB-PPP
    # RegionVariant shapes are field-compatible.
    exposure_track = [_eqtl_to_locus_variant(v) for v in exposure_variants]
    outcome_track = [_gwas_to_locus_variant(v) for v in outcome_rows]

    exposure_short_label, outcome_short_label = _build_short_labels_from_spec(
        spec,
        exposure_release=exposure_eqtl_release,
        pqtl_release=pqtl.release if spec.exposure_kind == EXPOSURE_KIND_UKB_PPP else None,
    )

    # Gene track. Two sources, in priority order:
    # 1. spec.prefetched_gene_track — caller-supplied list (e.g. from the
    #    locuscompare CLI's Ensembl REST on-demand fetcher).
    # 2. on-demand Ensembl REST fetch via fetch_region_genes_remote (fallback
    #    when no prefetched track was supplied; graceful-degrades to no track
    #    if Ensembl is offline).
    gene_track: list[GeneTrackEntry] = []
    if spec.prefetched_gene_track is not None:
        gene_track = list(spec.prefetched_gene_track)
    else:
        try:
            genes, release_meta, gt_notes = fetch_region_genes_remote(
                chromosome=chromosome.lstrip("chr"),
                start_bp=start_bp, end_bp=end_bp,
                biotypes=gene_biotypes,
            )
            gene_track = [
                GeneTrackEntry(
                    gene_symbol=g.gene_symbol,
                    start=g.start, end=g.end, strand=g.strand,
                    exons=[(e.start, e.end) for e in g.exons],
                    biotype=g.biotype,
                )
                for g in genes
            ]
            notes.extend(gt_notes)
        except (requests.RequestException, ValueError) as e:
            notes.append(f"gene track unavailable: {e!s}")
            extra_caveats.append("gene track unavailable (Ensembl REST error)")

    if spec.exposure_kind == EXPOSURE_KIND_UKB_PPP:
        exposure_label_long = (
            f"UKB-PPP pQTL ({exposure_source_release}); "
            f"protein {exposure_protein_label}; "
            f"ancestry {exposure_ancestry_label} ({exposure_ancestry_code})"
            f"{spec.exposure_id_extra}"
        )
    else:
        exposure_label_long = (
            f"eQTL Catalogue {exposure_source_release}; "
            f"study {exposure_study_id}{spec.exposure_id_extra}"
        )

    inp = RegionalLocusCompareInput(
        pairs=pairs,
        lead_variant_id=lead_variant_id,
        chromosome=chromosome,
        window_bp=window_bp,
        ld_panel_label=(
            f"{panel_id} ({super_pop.value}); plink {plink_version}"
            if ld_client is not None and panel_id != "none"
            else "no LD reference (plot rendered without LD coloring)"
        ),
        window_label=f"+/-{window_bp // 2000} kb of lead {lead_variant_id}{pip_label}",
        exposure_label=exposure_label_long,
        outcome_label=(
            f"GWAS Catalog harmonised; study {gwas_accession}{spec.outcome_id_extra}"
        ),
        provenance_label=f"{spec.provenance_prefix}Rendered: {fetched_at}",
        caveats=extra_caveats + spec.notes,
        exposure_track_variants=exposure_track,
        outcome_track_variants=outcome_track,
        r2_by_variant=dict(r2_by_variant),
        exposure_short_label=exposure_short_label,
        outcome_short_label=outcome_short_label,
        gene_track=gene_track,
        focal_gene_symbol=spec.exposure_gene_symbol or None,
        outcome_effect_scale_label=outcome_scale_label,
        effect_size_unavailable_reason=effect_size_reason,
    )
    render_full_locuscompare(inp, out_path)

    # 6. Manifest block. v1.3 adds three optional fields for pQTL renders
    # (exposure_protein_label, exposure_ancestry, exposure_ancestry_label)
    # so the caption builder + downstream tooling can distinguish a
    # `(plasma sortilin, EUR)` render from a `(GTEx liver, EUR)` render
    # purely from the manifest, without re-resolving the study ID.
    block = build_regional_locuscompare_block(
        ot_release=spec.release_tag,
        exposure_source=exposure_source,
        exposure_source_release=exposure_source_release,
        exposure_study_id=exposure_study_id,
        exposure_protein_label=exposure_protein_label,
        exposure_ancestry=exposure_ancestry_code,
        exposure_ancestry_label=exposure_ancestry_label,
        exposure_harmonisation_version="",  # neither source surfaces this
        outcome_source="gwas_catalog_harmonised",
        outcome_source_release=outcome.release.fetched_at_utc.split("T")[0],
        outcome_study_id=gwas_accession,
        outcome_harmonisation_version="",
        ld_panel=panel_id,
        ld_panel_super_pop=super_pop.value,
        ld_panel_version=panel_version,
        plink_version=plink_version,
        window_bp=window_bp,
        lead_variant_id=lead_variant_id,
        n_pairs=len(pairs),
        n_palindromic_excluded=n_palindromic,
        scatter_downsampled=len(pairs) > 5000,
        scatter_downsample_target=5000,
        ancestry_caveats=extra_caveats + spec.notes,
        data_source_warnings=data_source_warnings,
        plot_artifact=str(out_path.name),
        fetched_at=fetched_at,
    )
    # Where the outcome effect sizes came from, with the human label beside the code,
    # and why the effect-size panel is empty when it is.
    block["outcome_beta_source"] = outcome_beta_source
    block["outcome_beta_source_label"] = beta_source_label(outcome_beta_source)
    block["outcome_effect_scale_label"] = outcome_scale_label
    block["effect_size_panel_unavailable_reason"] = effect_size_reason

    return Tier2Result(
        plot_path=out_path,
        manifest_block=block,
        n_pairs=len(pairs),
        n_palindromic_excluded=n_palindromic,
        notes=notes,
    )


def _study_label(study_id: str) -> str:
    """A study id for a reader. The CLI's pre-fetched path has no accession and passes
    the placeholder `prefetched`, which names our input route, not a study."""
    return "supplied as a file" if study_id == "prefetched" else study_id


def _effect_size_unavailable_reason(
    *,
    pairs,
    exposure_rows,
    outcome_rows,
    exposure_label: str,
    outcome_label: str,
) -> str | None:
    """Why the effect-size panel has no point to draw, naming the side; None if it has.

    Mirrors the panel's own filter: palindromic pairs are left off it, and a point needs
    a beta on both sides. `outcome_rows` are the rows AFTER the odds-ratio conversion.
    """
    drawable = [p for p in pairs if not p.palindromic_excluded]
    if any(p.beta_exposure is not None and p.beta_outcome is not None for p in drawable):
        return None
    reasons: list[str] = []
    if not any(v.beta is not None for v in exposure_rows):
        reasons.append(
            f"the exposure study ({exposure_label}) publishes no per-variant effect "
            f"sizes in this window"
        )
    if not any(v.beta is not None for v in outcome_rows):
        if any(getattr(v, "odds_ratio", None) is not None for v in outcome_rows):
            # The conversion returns None when it cannot derive a standard error, and
            # does not fill a beta without one (OR_NO_SE_REASON).
            reasons.append(
                f"the outcome study ({outcome_label}) reports odds ratios from which "
                f"{OR_NO_SE_REASON}, so no log odds effect size could be derived"
            )
        else:
            reasons.append(
                f"the outcome study ({outcome_label}) publishes p-values only, with no "
                f"effect size or odds ratio, in this window"
            )
    if reasons:
        return "; ".join(reasons)
    if pairs and not drawable:
        return (
            "every variant joined across the two studies is strand-ambiguous "
            "(palindromic) and is left off this panel"
        )
    return (
        "no variant joined across the two studies carries an effect size on both sides"
    )


def _build_short_labels_from_spec(
    spec: LocusCompareSpec,
    *,
    exposure_release=None,
    pqtl_release=None,
) -> tuple[str, str]:
    """Construct front-and-center one-line panel titles for the manhattans.

    Branches by `spec.exposure_kind`:
    - eqtl_catalogue: reads study_label / sample_group / quant_method from
      the eQTL Catalogue release dataclass (existing behavior).
    - ukb_ppp: reads protein_label / ancestry_label from the UKB-PPP
      release dataclass; surfaces the protein + ancestry + sample size
      for the pQTL render's panel title.
    """
    # Outcome side (shared).
    if spec.outcome_trait_label:
        outcome_short = (
            f"Outcome (GWAS): {spec.outcome_trait_label} "
            f"({spec.gwas_accession})"
        )
    else:
        outcome_short = f"Outcome (GWAS): {spec.gwas_accession}"

    # Exposure side.
    if spec.exposure_kind == EXPOSURE_KIND_UKB_PPP and pqtl_release is not None:
        bits: list[str] = []
        if pqtl_release.protein_hgnc:
            bits.append(pqtl_release.protein_hgnc)
        bits.append("UKB-PPP")
        if pqtl_release.olink_panel:
            bits.append(f"Olink {pqtl_release.olink_panel}")
        if pqtl_release.ancestry_label:
            n_label = (
                f"; N = {pqtl_release.n_samples:,}"
                if pqtl_release.n_samples else ""
            )
            bits.append(f"{pqtl_release.ancestry_label}{n_label}")
        descriptor = " | ".join(bits) if bits else pqtl_release.protein_label
        exposure_short = (
            f"Exposure (pQTL): {descriptor} "
            f"({pqtl_release.olink_reagent_id})"
        )
        return exposure_short, outcome_short

    # eQTL Catalogue path (default).
    if exposure_release is None:
        # Should not happen for a well-formed spec; fall back to dataset_id.
        return (
            f"Exposure (eQTL): {spec.eqtl_dataset_id}",
            outcome_short,
        )
    bits = []
    if spec.exposure_gene_symbol:
        bits.append(spec.exposure_gene_symbol)
    if exposure_release.study_label:
        bits.append(exposure_release.study_label)
    if exposure_release.sample_group:
        bits.append(exposure_release.sample_group)
    elif exposure_release.tissue_label and exposure_release.condition_label:
        bits.append(f"{exposure_release.tissue_label} ({exposure_release.condition_label})")
    if exposure_release.quant_method:
        bits.append(exposure_release.quant_method)
    descriptor = " | ".join(bits) if bits else spec.eqtl_dataset_id
    exposure_short = (
        f"Exposure (eQTL): {descriptor} "
        f"({spec.eqtl_dataset_id})"
    )
    return exposure_short, outcome_short


def load_study_id_mappings(yaml_path: Path) -> dict[tuple[str, str], StudyIdMapping]:
    """Read the manual lookup table (YAML).

    Schema (v1.3, polymorphic by exposure kind):
      mappings:
        # eQTL / sQTL / sceQTL row (eQTL Catalogue backend):
        - ot_left_study_id: quach_2016_ge_monocyte_iav_ensg00000115808
          eqtl_catalogue_dataset_id: QTD000110
          ot_right_study_id: FINNGEN_R12_I9_HEARTFAIL
          gwas_catalog_accession: GCST90475990
          ancestry_left: EUR
          ancestry_right: EUR (Finnish)
          notes:
            - "Quach 2016 IAV-stimulated monocyte; verified live 2026-05-04"

        # pQTL row (UKB-PPP backend):
        - ot_left_study_id: UKB_PPP_EUR_SORT1
          ukb_ppp_protein_label: SORT1
          ukb_ppp_ancestry: EUR
          ot_right_study_id: GCST90269602
          gwas_catalog_accession: GCST90269602
          ancestry_left: EUR
          ancestry_right: EUR
          exposure_gene_symbol: SORT1
          outcome_trait_label: cholesterol VLDL
          notes:
            - "Sun 2023 UKB-PPP plasma cis-pQTL; verified live 2026-05-15"

    Exactly one of `eqtl_catalogue_dataset_id` or
    (`ukb_ppp_protein_label`, `ukb_ppp_ancestry`) must be populated per row;
    the orchestrator dispatches by the OT studyId prefix and reads the
    matching fields.

    Returns a dict keyed by (ot_left_study_id, ot_right_study_id) for fast
    per-row lookup.
    """
    yaml_path = Path(yaml_path)
    if not yaml_path.is_file():
        return {}
    data = yaml.safe_load(yaml_path.read_text()) or {}
    out: dict[tuple[str, str], StudyIdMapping] = {}
    for entry in data.get("mappings", []):
        m = StudyIdMapping(
            ot_left_study_id=entry["ot_left_study_id"],
            ot_right_study_id=entry["ot_right_study_id"],
            gwas_catalog_accession=entry["gwas_catalog_accession"],
            eqtl_catalogue_dataset_id=entry.get("eqtl_catalogue_dataset_id", ""),
            ukb_ppp_protein_label=entry.get("ukb_ppp_protein_label", ""),
            ukb_ppp_ancestry=entry.get("ukb_ppp_ancestry", ""),
            ancestry_left=entry.get("ancestry_left", "EUR"),
            ancestry_right=entry.get("ancestry_right", "EUR"),
            outcome_trait_label=entry.get("outcome_trait_label", ""),
            exposure_gene_symbol=entry.get("exposure_gene_symbol", ""),
            notes=list(entry.get("notes", [])),
        )
        out[(m.ot_left_study_id, m.ot_right_study_id)] = m
    return out


def _eqtl_to_locus_variant(v) -> LocusVariant:
    """Map an exposure-side RegionVariant (eQTL Catalogue or UKB-PPP; both
    expose the same field set) to a wald_ratio LocusVariant. The renderer's
    per-side manhattan reads the LocusVariant fields it shares with the
    wald_ratio harmoniser.
    """
    return LocusVariant(
        variant_id=v.variant_id,
        chromosome=v.chromosome,
        position=v.position,
        ref=v.ref,
        alt=v.alt,
        pip=None,
        beta=v.beta,
        se=v.se,
        p_value=v.p_value,
        is95=None,
        is99=None,
        r2_lead=None,
    )


def _gwas_to_locus_variant(v) -> LocusVariant:
    """Same as `_eqtl_to_locus_variant` for GWAS Catalog rows."""
    return LocusVariant(
        variant_id=v.variant_id,
        chromosome=v.chromosome,
        position=v.position,
        ref=v.ref,
        alt=v.alt,
        pip=None,
        beta=v.beta,
        se=v.se,
        p_value=v.p_value,
        is95=None,
        is99=None,
        r2_lead=None,
    )


def _now_utc() -> str:
    return time.strftime("%Y-%m-%dT%H:%M:%SZ", time.gmtime())


def _extract_ensg_from_ot_study_id(ot_study_id: str) -> str | None:
    """OT QTL studyIds end with `_ensg<digits>` (lowercase). Return the
    uppercased ENSG id (the format eQTL Catalogue's REST API expects) or
    None if the studyId does not match.
    """
    import re
    m = re.search(r"(ensg\d+)$", ot_study_id, flags=re.IGNORECASE)
    if m:
        return m.group(1).upper()
    return None


# ---- Manifest block builder (inlined from skills/decision/target_validation/
#      renderers/manifest.py for fork self-containment; keep byte-equivalent
#      if the canonical changes) ----

def build_regional_locuscompare_block(
    *,
    ot_release: str,
    exposure_source: str,
    exposure_source_release: str,
    exposure_study_id: str,
    outcome_source: str,
    outcome_source_release: str,
    outcome_study_id: str,
    ld_panel: str,
    ld_panel_super_pop: str,
    ld_panel_version: str,
    plink_version: str,
    window_bp: int,
    lead_variant_id: str,
    n_pairs: int,
    n_palindromic_excluded: int = 0,
    scatter_downsampled: bool = False,
    scatter_downsample_target: int = 5000,
    ancestry_caveats: list[str] | None = None,
    data_source_warnings: list[str] | None = None,
    exposure_harmonisation_version: str = "",
    outcome_harmonisation_version: str = "",
    plot_artifact: str | None = None,
    fetched_at: str | None = None,
    exposure_protein_label: str = "",
    exposure_ancestry: str = "",
    exposure_ancestry_label: str = "",
) -> dict[str, Any]:
    """Construct one `regional_locuscompare` manifest block.

    The orchestrator builds this dict per Tier-2 render and passes a list to
    `write_manifest(..., tier2_renders=[...])`.

    pQTL caption fields (exposure_protein_label, exposure_ancestry,
    exposure_ancestry_label) remain empty strings for eQTL / sQTL / sceQTL
    renders where the existing tissue / quant fields carry equivalent context;
    they're populated for UKB-PPP renders where plasma is implicit but
    protein identity isn't.

    `data_source_warnings` is an additive sibling to `ancestry_caveats`
    that carries fetcher-emitted warnings (schema drift, pagination notes,
    etc.) prefixed with the originating fetcher name. Kept separate from
    `ancestry_caveats` so the curated bucket stays clean and the noisy
    fetcher channel can be observed / triaged independently.
    """
    return {
        "ot_release": ot_release,
        "exposure_source": exposure_source,
        "exposure_source_release": exposure_source_release,
        "exposure_study_id": exposure_study_id,
        "exposure_protein_label": exposure_protein_label,
        "exposure_ancestry": exposure_ancestry,
        "exposure_ancestry_label": exposure_ancestry_label,
        "exposure_harmonisation_version": exposure_harmonisation_version,
        "outcome_source": outcome_source,
        "outcome_source_release": outcome_source_release,
        "outcome_study_id": outcome_study_id,
        "outcome_harmonisation_version": outcome_harmonisation_version,
        "ld_panel": ld_panel,
        "ld_panel_super_pop": ld_panel_super_pop,
        "ld_panel_version": ld_panel_version,
        "plink_version": plink_version,
        "window_bp": window_bp,
        "lead_variant_id": lead_variant_id,
        "n_pairs": n_pairs,
        "n_palindromic_excluded": n_palindromic_excluded,
        "scatter_downsampled": scatter_downsampled,
        "scatter_downsample_target": scatter_downsample_target,
        "ancestry_caveats": list(ancestry_caveats or []),
        "data_source_warnings": list(data_source_warnings or []),
        "plot_artifact": plot_artifact,
        "fetched_at": fetched_at or _now_utc(),
    }


def _now_utc() -> str:
    import time as _time
    return _time.strftime("%Y-%m-%dT%H:%M:%SZ", _time.gmtime())


__all__ = [
    "EXPOSURE_KIND_EQTL_CATALOGUE",
    "EXPOSURE_KIND_UKB_PPP",
    "LocusCompareSpec",
    "StudyIdMapping",
    "Tier2NotAvailable",
    "Tier2Result",
    "dispatch_exposure_kind",
    "is_pqtl_study_id",
    "load_study_id_mappings",
    "render_locuscompare_for_lead",
    "render_tier2_for_lead",
]
