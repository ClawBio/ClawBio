"""Unit tests for `_prefetched.py` adapter clients.

Audit gap (PR #272): existing tests in `test_locuscompare_region_render.py`
cover the canonical-TSV loaders. These tests cover the duck-typed stub
clients (PrefetchedEQTLClient / PrefetchedGWASClient / PrefetchedLDClient)
that satisfy the core orchestrator's fetcher interface, plus the
`gwas_variants_from_eqtl` reshape used on the outcome side.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest


SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

# Inject sibling-skill directories so the prefetched module's sibling imports
# resolve in the flat fork layout.
_SKILLS_ROOT = SKILL_DIR.parent
for _p in (
    SKILL_DIR,
    _SKILLS_ROOT / "eqtl-catalogue-region-fetch",
    _SKILLS_ROOT / "gwas-catalog-region-fetch",
    _SKILLS_ROOT / "ld-1000g-region-compute",
):
    if _p.exists() and str(_p) not in sys.path:
        sys.path.insert(0, str(_p))


from _prefetched import (  # noqa: E402
    PrefetchedEQTLClient,
    PrefetchedGWASClient,
    PrefetchedLDClient,
    gwas_variants_from_eqtl,
)
from eqtl_catalogue_region_fetch import RegionVariant as EQTLRegionVariant  # noqa: E402
from gwas_catalog_region_fetch import RegionVariant as GWASRegionVariant  # noqa: E402


def _eqtl_variant(vid="1_100_A_G", ref="A", alt="G", beta=0.1, se=0.01, p=1e-6):
    return EQTLRegionVariant(
        variant_id=vid,
        chromosome="1",
        position=int(vid.split("_")[1]),
        ref=ref,
        alt=alt,
        beta=beta,
        se=se,
        p_value=p,
        maf=0.2,
        effect_allele_frequency=0.2,
        raw={},
    )


def test_prefetched_eqtl_client_fetch_region_returns_loaded_variants():
    variants = [_eqtl_variant("1_100_A_G"), _eqtl_variant("1_200_C_T", ref="C", alt="T")]
    client = PrefetchedEQTLClient(
        variants=variants,
        dataset_id="QTD000276",
        dataset_release="prefetched",
        study_label="SORT1 expression — minor salivary gland",
        tissue_label="minor salivary gland",
        quant_method="ge",
    )

    result = client.fetch_region(
        dataset_id="ignored-passthrough",
        chromosome="1",
        start_bp=1,
        end_bp=1_000_000,
    )

    assert result.dataset_id == "QTD000276"
    assert result.n_variants == 2
    assert result.region_start_bp == 1
    assert result.region_end_bp == 1_000_000
    assert result.variants[0].variant_id == "1_100_A_G"
    assert any("pre-fetched" in n for n in result.notes)
    assert result.release.dataset_release == "prefetched"
    assert result.release.tissue_label == "minor salivary gland"


def test_prefetched_eqtl_client_falls_back_to_caller_dataset_id_when_empty():
    """`dataset_id=""` on the client should defer to the caller's value."""
    client = PrefetchedEQTLClient(variants=[_eqtl_variant()], dataset_id="")
    result = client.fetch_region(
        dataset_id="caller-supplied",
        chromosome="1",
        start_bp=1,
        end_bp=1_000,
    )
    assert result.dataset_id == "caller-supplied"


def test_prefetched_gwas_client_fetch_region_returns_loaded_variants():
    gwas_variants = gwas_variants_from_eqtl([_eqtl_variant("1_100_A_G")])
    client = PrefetchedGWASClient(variants=gwas_variants, accession="GCST90269602")
    result = client.fetch_region(
        accession="ignored-passthrough",
        chromosome="1",
        start_bp=1,
        end_bp=1_000_000,
    )
    assert result.accession == "GCST90269602"
    assert result.n_variants == 1
    assert result.variants[0].variant_id == "1_100_A_G"
    assert any("pre-fetched" in n for n in result.notes)


def test_gwas_variants_from_eqtl_preserves_core_fields():
    eqtl = [
        _eqtl_variant("1_100_A_G", ref="A", alt="G", beta=0.5, se=0.05, p=1e-10),
        _eqtl_variant("1_200_C_T", ref="C", alt="T", beta=-0.3, se=0.04, p=1e-8),
    ]
    gwas = gwas_variants_from_eqtl(eqtl)
    assert len(gwas) == 2
    assert all(isinstance(v, GWASRegionVariant) for v in gwas)
    assert [v.variant_id for v in gwas] == ["1_100_A_G", "1_200_C_T"]
    assert gwas[0].beta == pytest.approx(0.5)
    assert gwas[1].beta == pytest.approx(-0.3)
    # GWAS-specific field gets a sensible default (odds_ratio absent on eQTL side).
    assert gwas[0].odds_ratio is None
    assert gwas[1].odds_ratio is None


def test_prefetched_ld_client_returns_only_known_partners():
    """The LD adapter must drop partner ids missing from the matrix and
    surface counts honestly (`n_partners_requested` vs `n_partners_returned`)."""
    client = PrefetchedLDClient(
        r2_by_partner={
            "1_100_A_G": 0.95,
            "1_200_C_T": 0.40,
        },
        panel_id="synthetic-test",
        panel_version="1.0",
        plink_version="prefetched",
        super_pop="synthetic",
    )

    result = client.r2_with_lead(
        lead="1_500_A_T",
        partners=["1_100_A_G", "1_999_T_C", "1_200_C_T"],
        chromosome="1",
        window_bp=1_000_000,
    )

    assert result.lead_variant_id == "1_500_A_T"
    assert result.window_bp == 1_000_000
    assert result.panel_id == "synthetic-test"
    assert result.plink_version == "prefetched"
    assert result.n_partners_requested == 3
    assert result.n_partners_returned == 2
    returned_ids = {p.partner_variant_id for p in result.pairs}
    assert returned_ids == {"1_100_A_G", "1_200_C_T"}
    r2_by_id = {p.partner_variant_id: p.r2 for p in result.pairs}
    assert r2_by_id["1_100_A_G"] == pytest.approx(0.95)
    assert r2_by_id["1_200_C_T"] == pytest.approx(0.40)
    assert any("pre-fetched synthetic matrix" in n for n in result.notes)


def test_prefetched_ld_client_handles_empty_partner_list():
    client = PrefetchedLDClient(r2_by_partner={"1_100_A_G": 0.9})
    result = client.r2_with_lead(
        lead="1_500_A_T", partners=[], chromosome="1", window_bp=1_000_000,
    )
    assert result.n_partners_requested == 0
    assert result.n_partners_returned == 0
    assert result.pairs == []


# ---------- One-row-per-rsID copies in a pre-fetched slice ----------

_TSV_HEADER = "variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\tp\trsid\n"


def _client_from_tsv(tmp_path: Path, rows: list[str]) -> PrefetchedEQTLClient:
    from _prefetched import load_sumstats_tsv

    tsv = tmp_path / "exposure.tsv"
    tsv.write_text(_TSV_HEADER + "".join(r + "\n" for r in rows))
    return PrefetchedEQTLClient(variants=load_sumstats_tsv(tsv), dataset_id="QTD000276")


def _fetch(client: PrefetchedEQTLClient):
    return client.fetch_region(dataset_id="x", chromosome="1", start_bp=1, end_bp=10_000)


def test_prefetched_eqtl_client_collapses_duplicate_rsid_rows(tmp_path: Path):
    """A slice cut from the eQTL Catalogue repeats a variant once per rsID with
    identical statistics: the render gets one row, keeping both rsIDs, and the
    dropped copy is recorded with its reason."""
    result = _fetch(_client_from_tsv(tmp_path, [
        "1_100_T_TA\t1\t100\tT\tTA\t0.2\t0.01\t1e-6\trs1",
        "1_100_T_TA\t1\t100\tT\tTA\t0.2\t0.01\t1e-6\trs2",
        "1_200_C_T\t1\t200\tC\tT\t0.1\t0.01\t1e-3\trs3",
    ]))
    assert [v.variant_id for v in result.variants] == ["1_100_T_TA", "1_200_C_T"]
    assert result.n_variants == 2
    assert result.n_rows_read == 3
    assert result.variants[0].raw["rsids"] == "rs1;rs2"
    assert [(c.variant_id, c.rsid, c.kept_rsid, c.reason) for c in result.collapsed] == [
        ("1_100_T_TA", "rs2", "rs1", "duplicate_rsid_row"),
    ]
    assert any("1 repeated row(s)" in n for n in result.notes)


def test_prefetched_eqtl_client_refuses_copies_that_disagree(tmp_path: Path):
    """Rows of one variant that disagree cannot be merged, and picking one would
    be a guess: the render stops. The error is the pre-fetched file's own, so
    the render does not report a local file as an eQTL Catalogue failure."""
    from _prefetched import PrefetchedSchemaError
    from eqtl_catalogue_region_fetch import EQTLCatalogueAPIError

    client = _client_from_tsv(tmp_path, [
        "1_100_A_G\t1\t100\tA\tG\t0.2\t0.01\t1e-6\trs1",
        "1_100_A_G\t1\t100\tA\tG\t-0.3\t0.02\t1e-4\trs1",
    ])
    with pytest.raises(PrefetchedSchemaError, match=r"disagree on beta \(0.2 vs -0.3\).*1_100_A_G") as exc:
        _fetch(client)
    assert not isinstance(exc.value, EQTLCatalogueAPIError)
    assert "catalogue" not in str(exc.value).lower()


def test_prefetched_eqtl_client_refuses_a_file_with_several_traits(tmp_path: Path):
    """The plot draws one molecular trait; a file mixing several would pair a
    variant once per trait, so it stops and names them."""
    from _prefetched import PrefetchedSchemaError, load_sumstats_tsv

    tsv = tmp_path / "exposure.tsv"
    tsv.write_text(
        "variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\tp\tmolecular_trait_id\n"
        "1_100_A_G\t1\t100\tA\tG\t0.2\t0.01\t1e-6\tENST1\n"
        "1_100_A_G\t1\t100\tA\tG\t-0.3\t0.02\t1e-4\tENST2\n"
    )
    client = PrefetchedEQTLClient(variants=load_sumstats_tsv(tsv), dataset_id="QTD000276")
    with pytest.raises(PrefetchedSchemaError, match="2 molecular traits.*ENST1, ENST2"):
        _fetch(client)


def test_prefetched_eqtl_client_collapses_copies_with_nan_and_na(tmp_path: Path):
    """`NaN` reads as missing like `NA`, so identical copies carrying it still
    collapse; an `NA` rsID is not stored as an rsID."""
    from _prefetched import load_sumstats_tsv

    tsv = tmp_path / "exposure.tsv"
    tsv.write_text(
        "variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\tp\tmaf\trsid\n"
        "1_100_T_TA\t1\t100\tT\tTA\t0.2\t0.01\t1e-6\tnan\tNA\n"
        "1_100_T_TA\t1\t100\tT\tTA\t0.2\t0.01\t1e-6\tNaN\trs2\n"
    )
    result = _fetch(PrefetchedEQTLClient(variants=load_sumstats_tsv(tsv), dataset_id="QTD000276"))
    assert result.n_variants == 1
    assert "rsid" not in result.variants[0].raw
    assert result.variants[0].raw["rsids"] == "rs2"


def test_prefetched_eqtl_client_without_repeats_is_unchanged(tmp_path: Path):
    result = _fetch(_client_from_tsv(tmp_path, [
        "1_100_A_G\t1\t100\tA\tG\t0.2\t0.01\t1e-6\trs1",
        "1_200_C_T\t1\t200\tC\tT\t0.1\t0.01\t1e-3\t",
    ]))
    assert result.n_variants == 2
    assert result.collapsed == []
    assert result.variants[1].raw.get("rsid") is None
    assert result.notes == ["sumstats loaded from pre-fetched TSV (no live fetch)"]
