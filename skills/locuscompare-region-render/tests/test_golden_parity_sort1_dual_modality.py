"""Golden-parity test: SORT1 dual-modality (eQTL + pQTL) render.

This single fixture exercises the cross-backend dispatch + shared-outcome
harmonisation that the v1.3 release added. Two render passes run side by
side, both reading from offline cassettes captured 2026-05-15:

- eQTL exposure: GTEx minor salivary gland ge-eQTL (eQTL Cat QTD000276)
  joined to GWAS Catalog harmonised GCST90269602 (cholesterol in medium
  VLDL).
- pQTL exposure: UKB-PPP EUR SORT1 plasma sortilin (bundled slice) joined
  to the SAME GCST90269602 outcome.

Both renders share the canonical SORT1 1p13.3 lead variant
(`1_109274968_G_T`, rs12740374, Musunuru 2010) and a +/-500 kb window.

The manifest_block fields locked in `expected.yaml` cover everything that
should not drift across code changes: study/dataset identifiers, lead
variant, window, harmonisation counts (n_pairs / n_palindromic_excluded),
and the protein label / ancestry that distinguish the pQTL row from the
eQTL row in caption text.

Fields explicitly NOT locked (see expected.yaml comments):
- `fetched_at`: timestamp emitted at render time
- `plot_artifact`: filename basename varies with out_path
- `ld_panel*` / `plink_version`: depend on whether plink 1.9 is installed
  on the runner; the fixture is captured with `ld_panel == "none"` and
  the LD block stays out of the assertion.

Refresh policy: a failing assertion means the orchestrator changed
behaviour. Inspect the diff. Do NOT refresh expected.yaml without
explicit rationale; the fixture is the source of truth for cross-
backend dispatch.
"""

from __future__ import annotations

import gzip
import json
import sys
from pathlib import Path
from typing import Any

import pytest
import yaml

# Inject sibling-skill dirs onto sys.path so the bare-name imports below
# resolve under ClawBio's flat skills/<kebab-name>/ layout. These four
# skills are the published primitives this orchestrator composes.
_SKILL_DIR = Path(__file__).resolve().parent.parent
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

# Sibling-skill imports. On a standalone fork checkout (this branch only,
# with the other four primitives not yet merged into upstream/main),
# collection would fail at import time. Skip the whole module instead so
# the rest of the suite still runs. The golden-parity assertions are
# meaningful only once the eqtl / gwas / ld / ukb-ppp primitives sit
# next to this skill in the same checkout.
try:
    from eqtl_catalogue_region_fetch import (  # noqa: E402
        EQTLCatalogueClient,
        EQTLCatalogueRelease,
    )
    from eqtl_catalogue_region_fetch import (  # noqa: E402
        RegionResult as EQTLRegionResult,
        RegionVariant as EQTLRegionVariant,
    )
    from gwas_catalog_region_fetch import (  # noqa: E402
        GWASCatalogClient,
        GWASCatalogRelease,
    )
    from gwas_catalog_region_fetch import (  # noqa: E402
        RegionResult as GWASRegionResult,
        RegionVariant as GWASRegionVariant,
    )
    from ld_1000g_region_compute import SuperPop  # noqa: E402
    from locuscompare_region_render import (  # noqa: E402
        EXPOSURE_KIND_EQTL_CATALOGUE,
        EXPOSURE_KIND_UKB_PPP,
        LocusCompareSpec,
        render_locuscompare_for_lead,
    )
    from ukb_ppp_region_fetch import (  # noqa: E402
        UKBPPPClient,
        UKBPPPRelease,
    )
    from ukb_ppp_region_fetch import (  # noqa: E402
        RegionResult as UKBPPPRegionResult,
        RegionVariant as UKBPPPRegionVariant,
    )
except ModuleNotFoundError as _e:
    pytest.skip(
        f"golden-parity replay requires the 4 sibling skills "
        f"(eqtl-catalogue-region-fetch, gwas-catalog-region-fetch, "
        f"ld-1000g-region-compute, ukb-ppp-region-fetch) checked out "
        f"alongside this one: {_e}",
        allow_module_level=True,
    )


FIXTURE_DIR = Path(__file__).resolve().parent / "fixtures" / "golden" / "sort1_dual_modality_eqtl_pqtl_ldl"
INPUTS_DIR = FIXTURE_DIR / "inputs"
EXPECTED_PATH = FIXTURE_DIR / "expected.yaml"


# ----------------------- cassette loaders -----------------------


def _load_eqtl_cassette() -> EQTLRegionResult:
    with gzip.open(
        INPUTS_DIR / "eqtl_exposure_QTD000276_chr1_109_274_968_window.json.gz",
        mode="rt",
    ) as f:
        data = json.load(f)
    variants = [
        EQTLRegionVariant(**{k: v for k, v in row.items() if k != "raw"} | {"raw": row.get("raw", {})})
        for row in data["variants"]
    ]
    release = EQTLCatalogueRelease(**data["release"])
    return EQTLRegionResult(
        dataset_id=data["dataset_id"],
        chromosome=data["chromosome"],
        region_start_bp=data["region_start_bp"],
        region_end_bp=data["region_end_bp"],
        n_variants=data["n_variants"],
        variants=variants,
        release=release,
        notes=list(data.get("notes", [])),
    )


def _load_gwas_cassette() -> GWASRegionResult:
    with gzip.open(
        INPUTS_DIR / "gwas_outcome_GCST90269602_chr1_109_274_968_window.json.gz",
        mode="rt",
    ) as f:
        data = json.load(f)
    variants = [
        GWASRegionVariant(**{k: v for k, v in row.items() if k != "raw"} | {"raw": row.get("raw", {})})
        for row in data["variants"]
    ]
    release = GWASCatalogRelease(**data["release"])
    return GWASRegionResult(
        accession=data["accession"],
        chromosome=data["chromosome"],
        region_start_bp=data["region_start_bp"],
        region_end_bp=data["region_end_bp"],
        n_variants=data["n_variants"],
        variants=variants,
        release=release,
        notes=list(data.get("notes", [])),
    )


def _load_pqtl_cassette() -> UKBPPPRegionResult:
    with gzip.open(
        INPUTS_DIR / "pqtl_exposure_SORT1_EUR.json.gz",
        mode="rt",
    ) as f:
        data = json.load(f)
    variants = [
        UKBPPPRegionVariant(**{k: v for k, v in row.items() if k != "raw"} | {"raw": row.get("raw", {})})
        for row in data["variants"]
    ]
    release = UKBPPPRelease(**data["release"])
    return UKBPPPRegionResult(
        protein_label_short=data["protein_label_short"],
        ancestry=data["ancestry"],
        chromosome=data["chromosome"],
        region_start_bp=data["region_start_bp"],
        region_end_bp=data["region_end_bp"],
        n_variants=data["n_variants"],
        variants=variants,
        release=release,
        notes=list(data.get("notes", [])),
    )


# ----------------------- mocked clients -----------------------


class _CassetteEQTLClient:
    """Stand-in for `EQTLCatalogueClient` that returns the offline cassette
    for any fetch_region call. Records the call so the test can confirm the
    orchestrator called us with the expected window."""

    def __init__(self, cassette: EQTLRegionResult) -> None:
        self._cassette = cassette
        self.calls: list[dict[str, Any]] = []

    def fetch_region(self, **kwargs: Any) -> EQTLRegionResult:
        self.calls.append(kwargs)
        return self._cassette


class _CassetteGWASClient:
    def __init__(self, cassette: GWASRegionResult) -> None:
        self._cassette = cassette
        self.calls: list[dict[str, Any]] = []

    def fetch_region(self, **kwargs: Any) -> GWASRegionResult:
        self.calls.append(kwargs)
        return self._cassette


class _CassetteUKBPPPClient:
    def __init__(self, cassette: UKBPPPRegionResult) -> None:
        self._cassette = cassette
        self.calls: list[dict[str, Any]] = []

    def fetch_region(self, **kwargs: Any) -> UKBPPPRegionResult:
        self.calls.append(kwargs)
        return self._cassette


# ----------------------- spec builders -----------------------


def _eqtl_spec() -> LocusCompareSpec:
    """eQTL Catalogue exposure (QTD000276, SORT1 ge in minor salivary gland)
    vs GCST90269602 cholesterol-VLDL outcome."""
    return LocusCompareSpec(
        lead_variant_id="1_109274968_G_T",
        chromosome="1",
        lead_position_bp=109_274_968,
        window_bp=1_000_000,
        eqtl_dataset_id="QTD000276",
        molecular_trait_id="ENSG00000134243",
        gwas_accession="GCST90269602",
        exposure_kind=EXPOSURE_KIND_EQTL_CATALOGUE,
        exposure_gene_symbol="SORT1",
        outcome_trait_label="cholesterol in medium VLDL",
        release_tag="26.03",
        super_pop=SuperPop.EUR,
        # Skip the network gene-track fetch; renderer surfaces a caveat,
        # which is fine for the manifest-block assertion.
        prefetched_gene_track=[],
    )


def _pqtl_spec() -> LocusCompareSpec:
    """UKB-PPP EUR exposure (SORT1 plasma sortilin) vs the same outcome."""
    return LocusCompareSpec(
        lead_variant_id="1_109274968_G_T",
        chromosome="1",
        lead_position_bp=109_274_968,
        window_bp=1_000_000,
        eqtl_dataset_id="",
        molecular_trait_id=None,
        gwas_accession="GCST90269602",
        exposure_kind=EXPOSURE_KIND_UKB_PPP,
        pqtl_protein_label="SORT1",
        pqtl_ancestry="EUR",
        exposure_gene_symbol="SORT1",
        outcome_trait_label="cholesterol in medium VLDL",
        release_tag="26.03",
        super_pop=SuperPop.EUR,
        prefetched_gene_track=[],
    )


# ----------------------- the actual replay test -----------------------


@pytest.fixture(scope="module")
def expected_blocks() -> dict[str, dict[str, Any]]:
    """Locked manifest_block fields for both renders (see expected.yaml)."""
    return yaml.safe_load(EXPECTED_PATH.read_text())


def _assert_block_matches_locked(
    actual: dict[str, Any],
    locked: dict[str, Any],
) -> None:
    """Assert each locked key in `locked` matches the same key in `actual`.

    The test asserts only the keys named in `locked` (a deliberate subset
    of `actual`). Fields like `fetched_at`, `plot_artifact`, and the LD
    block are intentionally absent from `locked` and not compared here.
    """
    diffs: list[str] = []
    for key, want in locked.items():
        got = actual.get(key)
        if got != want:
            diffs.append(f"  {key}: got {got!r}, expected {want!r}")
    assert not diffs, "manifest_block drifted from golden expected.yaml:\n" + "\n".join(diffs)


def test_sort1_dual_modality_eqtl_pqtl_ldl_replay(tmp_path, expected_blocks):
    """Replay both SORT1 exposure backends offline and verify both
    manifest_blocks match the locked expected.yaml. Fails loudly if either
    the eQTL or pQTL path drifts away from the 2026-05-15 capture."""
    eqtl_client = _CassetteEQTLClient(_load_eqtl_cassette())
    gwas_client = _CassetteGWASClient(_load_gwas_cassette())
    pqtl_client = _CassetteUKBPPPClient(_load_pqtl_cassette())

    # --- eQTL render ---
    eqtl_result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "eqtl" / "sort1_eqtl.png",
        ukb_ppp_client=None,
    )
    _assert_block_matches_locked(
        eqtl_result.manifest_block, expected_blocks["eqtl_render"]["manifest_block"]
    )
    assert eqtl_result.n_pairs == expected_blocks["eqtl_render"]["manifest_block"]["n_pairs"]
    assert (
        eqtl_result.n_palindromic_excluded
        == expected_blocks["eqtl_render"]["manifest_block"]["n_palindromic_excluded"]
    )

    # --- pQTL render (same outcome / lead / window, different exposure backend) ---
    pqtl_result = render_locuscompare_for_lead(
        spec=_pqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "pqtl" / "sort1_pqtl.png",
        ukb_ppp_client=pqtl_client,  # type: ignore[arg-type]
    )
    _assert_block_matches_locked(
        pqtl_result.manifest_block, expected_blocks["pqtl_render"]["manifest_block"]
    )
    assert pqtl_result.n_pairs == expected_blocks["pqtl_render"]["manifest_block"]["n_pairs"]
    assert (
        pqtl_result.n_palindromic_excluded
        == expected_blocks["pqtl_render"]["manifest_block"]["n_palindromic_excluded"]
    )

    # Verify the dispatch wiring fired the right backend per spec:
    # eqtl_client got 1 call (eqtl render), pqtl_client got 1 call (pqtl render).
    assert len(eqtl_client.calls) == 1, "eQTL client should be called exactly once (eQTL render)"
    assert len(pqtl_client.calls) == 1, "pQTL client should be called exactly once (pQTL render)"
    # Both renders share the same outcome accession (the dual-modality contract).
    assert len(gwas_client.calls) == 2, "GWAS client should be called twice (once per render)"
    assert all(c["accession"] == "GCST90269602" for c in gwas_client.calls)


# ----------------------- caveat for credible-set-filtered sumstats -----------------------


def _caveat_marker(qm: str) -> str:
    return f"credible-set-filtered (eQTL Catalogue .cc.tsv.gz; quant_method={qm})"


@pytest.mark.parametrize("quant_method", ["txrev", "leafcutter", "exon", "tx"])
def test_eqtl_render_adds_credible_set_caveat_for_non_ge_quant_methods(
    tmp_path, quant_method,
):
    """Per the v1.3 sQTL caption-disclosure update: when the eQTL exposure
    pulls from a .cc.tsv.gz file (any quant_method other than ge or
    microarray), the manifest must surface a user-facing caveat
    explaining the credible-set filter. Without it, a user could mistake
    the sparse rendered window for a full nominal-pass run."""
    cassette = _load_eqtl_cassette()
    cassette.release.quant_method = quant_method
    eqtl_client = _CassetteEQTLClient(cassette)
    gwas_client = _CassetteGWASClient(_load_gwas_cassette())
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / f"sort1_eqtl_{quant_method}.png",
        ukb_ppp_client=None,
    )
    caveats = result.manifest_block["ancestry_caveats"]
    assert any(_caveat_marker(quant_method) in c for c in caveats), (
        f"expected credible-set caveat for quant_method={quant_method!r}; "
        f"got caveats: {caveats}"
    )


@pytest.mark.parametrize("quant_method", ["ge", "microarray", ""])
def test_eqtl_render_omits_credible_set_caveat_for_full_nominal_pass(
    tmp_path, quant_method,
):
    """ge + microarray exposures pull from .all.tsv.gz (full nominal pass),
    so the credible-set caveat must NOT fire. Empty / unknown quant_method
    falls back to the legacy .all path and likewise should not fire."""
    cassette = _load_eqtl_cassette()
    cassette.release.quant_method = quant_method or None
    eqtl_client = _CassetteEQTLClient(cassette)
    gwas_client = _CassetteGWASClient(_load_gwas_cassette())
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / f"sort1_eqtl_{quant_method or 'unknown'}.png",
        ukb_ppp_client=None,
    )
    caveats = result.manifest_block["ancestry_caveats"]
    assert not any("credible-set-filtered" in c for c in caveats), (
        f"unexpected credible-set caveat for quant_method={quant_method!r}; "
        f"got caveats: {caveats}"
    )


def test_a_cc_backed_panel_discloses_the_filter_even_with_no_quant_method(tmp_path):
    """The disclosure is about the FILE CLASS, so losing the quant_method label
    must not lose it. `quant_method` is an optional presentation label and is None
    when the dataset metadata is unavailable; gated on it, a genuinely `.cc`-backed
    panel rendered with no disclosure that its summary statistics are credible-set
    filtered, and a reader not told the panel is filtered reads an absence as
    biology."""
    cassette = _load_eqtl_cassette()
    cassette.release.quant_method = None
    cassette.file_class = "cc"
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=_CassetteEQTLClient(cassette),  # type: ignore[arg-type]
        gwas_client=_CassetteGWASClient(_load_gwas_cassette()),  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "sort1_cc_no_qm.png",
        ukb_ppp_client=None,
    )
    caveats = result.manifest_block["ancestry_caveats"]
    assert any("credible-set-filtered" in c for c in caveats), (
        f"a .cc-backed panel rendered with no filter disclosure; caveats: {caveats}")
    # and it reads the same with or without the label, so it never prints a
    # dangling empty clause
    assert not any("quant_method=)" in c or "quant_method=;" in c for c in caveats)


def test_an_all_backed_panel_is_not_labelled_credible_set_filtered(tmp_path):
    """The other direction, and it is an affirmatively WRONG statement rather than
    a missing one. QTD000584 is an aptamer dataset that ships `.all` regardless of
    its quant_method (documented in eqtl-catalogue-region-fetch), so the
    quant-method rule told the reader a full-nominal-pass panel was credible-set
    filtered. The class the fetch actually served is the authority."""
    cassette = _load_eqtl_cassette()
    cassette.release.quant_method = "aptamer"
    cassette.file_class = "all"
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=_CassetteEQTLClient(cassette),  # type: ignore[arg-type]
        gwas_client=_CassetteGWASClient(_load_gwas_cassette()),  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "sort1_all_aptamer.png",
        ukb_ppp_client=None,
    )
    caveats = result.manifest_block["ancestry_caveats"]
    assert not any("credible-set-filtered" in c for c in caveats), (
        f"an .all-backed panel was labelled credible-set filtered; caveats: {caveats}")


# ----------------------- fetcher-notes propagation -----------------------
#
# Each fetcher's `RegionResult.notes` (schema-drift complaints, pagination
# warnings, etc.) propagates into the manifest's `data_source_warnings`
# field, prefixed with the originating fetcher name.


def test_eqtl_fetcher_notes_propagate_to_data_source_warnings(tmp_path):
    """A note emitted by the eQTL Catalogue fetcher must surface in the
    manifest's `data_source_warnings` list with the `eqtl_catalogue:`
    prefix. Models the schema-drift complaint the fetcher already emits
    today when a row's column count doesn't match the header."""
    eqtl_cassette = _load_eqtl_cassette()
    eqtl_cassette.notes = ["row column count 18 != header 19; skipping"]
    eqtl_client = _CassetteEQTLClient(eqtl_cassette)
    gwas_client = _CassetteGWASClient(_load_gwas_cassette())
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "eqtl_note_prop.png",
        ukb_ppp_client=None,
    )
    warnings = result.manifest_block["data_source_warnings"]
    assert (
        "eqtl_catalogue: row column count 18 != header 19; skipping" in warnings
    ), f"expected eqtl_catalogue-prefixed note in data_source_warnings; got {warnings}"


def test_gwas_fetcher_notes_propagate_to_data_source_warnings(tmp_path):
    """A note emitted by the GWAS Catalog fetcher must surface in the
    manifest's `data_source_warnings` list with the `gwas_catalog:`
    prefix."""
    eqtl_client = _CassetteEQTLClient(_load_eqtl_cassette())
    gwas_cassette = _load_gwas_cassette()
    gwas_cassette.notes = ["row column count 14 != header 15; skipping"]
    gwas_client = _CassetteGWASClient(gwas_cassette)
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "gwas_note_prop.png",
        ukb_ppp_client=None,
    )
    warnings = result.manifest_block["data_source_warnings"]
    assert (
        "gwas_catalog: row column count 14 != header 15; skipping" in warnings
    ), f"expected gwas_catalog-prefixed note in data_source_warnings; got {warnings}"


def test_ukb_ppp_fetcher_notes_propagate_to_data_source_warnings(tmp_path):
    """A note emitted by the UKB-PPP fetcher must surface in the manifest's
    `data_source_warnings` list with the `ukb_ppp:` prefix. Triggered on
    the pQTL render path."""
    eqtl_client = _CassetteEQTLClient(_load_eqtl_cassette())
    gwas_client = _CassetteGWASClient(_load_gwas_cassette())
    pqtl_cassette = _load_pqtl_cassette()
    pqtl_cassette.notes = ["row column count 12 != header 13; skipping"]
    pqtl_client = _CassetteUKBPPPClient(pqtl_cassette)
    result = render_locuscompare_for_lead(
        spec=_pqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "ukb_ppp_note_prop.png",
        ukb_ppp_client=pqtl_client,  # type: ignore[arg-type]
    )
    warnings = result.manifest_block["data_source_warnings"]
    assert (
        "ukb_ppp: row column count 12 != header 13; skipping" in warnings
    ), f"expected ukb_ppp-prefixed note in data_source_warnings; got {warnings}"


def test_data_source_warnings_does_not_pollute_ancestry_caveats(tmp_path):
    """Regression guard: fetcher notes go into `data_source_warnings`
    (separate manifest field), not into the curated `ancestry_caveats`
    bucket. A render with fetcher notes from BOTH sides must leave
    `ancestry_caveats` containing only the curated entries (here: the
    no-plink2 grey-LD caveat) without any prefixed fetcher messages."""
    eqtl_cassette = _load_eqtl_cassette()
    eqtl_cassette.notes = ["eqtl note from cassette"]
    eqtl_client = _CassetteEQTLClient(eqtl_cassette)
    gwas_cassette = _load_gwas_cassette()
    gwas_cassette.notes = ["gwas note from cassette"]
    gwas_client = _CassetteGWASClient(gwas_cassette)
    result = render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=eqtl_client,  # type: ignore[arg-type]
        gwas_client=gwas_client,  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "no_pollution.png",
        ukb_ppp_client=None,
    )
    caveats = result.manifest_block["ancestry_caveats"]
    warnings = result.manifest_block["data_source_warnings"]
    # Fetcher notes belong only to data_source_warnings.
    assert "eqtl_catalogue: eqtl note from cassette" in warnings
    assert "gwas_catalog: gwas note from cassette" in warnings
    # ancestry_caveats stays curated (the v1.3 no-plink2 caveat is the only
    # entry expected here; no fetcher prefixes should appear).
    assert not any(
        c.startswith(("eqtl_catalogue:", "gwas_catalog:", "ukb_ppp:")) for c in caveats
    ), f"fetcher notes leaked into ancestry_caveats: {caveats}"
