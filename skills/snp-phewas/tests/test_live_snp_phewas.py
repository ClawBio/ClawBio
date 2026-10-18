"""Live smoke (network): the response SHAPES the offline fixtures assume, the golden values
on the live sources, and the two source behaviours this skill exists to handle.

Run with: RUN_LIVE_TESTS=1 python -m pytest skills/snp-phewas/tests/test_live_snp_phewas.py
"""

from __future__ import annotations

import os

import pytest

pytestmark = [
    pytest.mark.live,
    pytest.mark.network,
    pytest.mark.skipif(os.environ.get("RUN_LIVE_TESTS") != "1", reason="set RUN_LIVE_TESTS=1"),
]

PAN = "https://pan-ukb-us-east-1.s3.amazonaws.com/"
LDL = "biomarkers-30780-both_sexes-irnt.tsv.bgz"
FG = "https://storage.googleapis.com/finngen-public-data-r12/summary_stats/release/finngen_R12_I9_CHD.gz"


@pytest.fixture(scope="module")
def cache(tmp_path_factory):
    return str(tmp_path_factory.mktemp("tbi"))


def test_finngen_api_shape_golden_and_p_threshold():
    from snp_phewas_pheweb import FinnGenPhewebClient

    payload = FinnGenPhewebClient().fetch_variant("1", 109274968, "G", "T")
    res = {r["phenocode"]: r for r in payload["results"]}
    assert len(res) >= 2469
    chd = res["I9_CHD"]
    assert {"beta", "sebeta", "pval", "mlogp", "maf", "n_case", "n_control"} <= set(chd)
    assert chd["beta"] == pytest.approx(-0.0801145) and chd["sebeta"] == pytest.approx(0.00864016)
    assert chd["pval"] == pytest.approx(1.82e-20, rel=1e-2)
    # The API fills values only below p = 0.05; the rest come back null (snp_phewas_pheweb).
    valued = [r for r in payload["results"] if r.get("pval") is not None]
    assert valued and all(r["pval"] < 0.05 for r in valued)
    assert any(r.get("pval") is None for r in payload["results"])
    assert FinnGenPhewebClient().fetch_variant("1", 109274968, "G", "A") is None


def test_finngen_tabix_equals_api_live(cache):
    from snp_phewas_tabix import ensure_index, fetch_position, read_header

    hdr = read_header(FG)
    (row,) = fetch_position(FG, ensure_index(FG + ".tbi", cache), "1", 109274968)
    rec = dict(zip(hdr, row))
    assert (rec["ref"], rec["alt"]) == ("G", "T")
    assert float(rec["beta"]) == pytest.approx(-0.0801145)
    assert float(rec["sebeta"]) == pytest.approx(0.00864016)
    assert float(rec["af_alt"]) == pytest.approx(0.214503)


def test_pan_ukbb_ldl_golden_and_the_index_defect(cache):
    from snp_phewas_tabix import ensure_index, fetch_position, read_header

    idx = ensure_index(PAN + "sumstats_flat_files_tabix/" + LDL + ".tbi", cache)
    hdr = read_header(PAN + "sumstats_flat_files/" + LDL)
    assert {"chr", "pos", "ref", "alt", "beta_EUR", "se_EUR", "neglog10_pval_EUR", "af_EUR",
            "low_confidence_EUR"} <= set(hdr)
    (row,) = fetch_position(PAN + "sumstats_flat_files/" + LDL, idx, "1", 109817590)
    rec = dict(zip(hdr, row))
    assert float(rec["beta_EUR"]) == pytest.approx(-0.1345)
    assert float(rec["se_EUR"]) == pytest.approx(0.002913)
    assert float(rec["neglog10_pval_EUR"]) == pytest.approx(464.8)
    # Positions an exact-position htslib query on this file missed (measured 2026-10-17); the
    # reader must find every one of them.
    for chrom, pos in (("1", 109817590), ("5", 51097763), ("5", 51401766), ("5", 50327893)):
        assert fetch_position(PAN + "sumstats_flat_files/" + LDL, idx, chrom, pos), (chrom, pos)


def test_ensembl_grch37_mapping_live():
    from snp_phewas_variant import Placement, map_placement

    p37, note = map_placement(Placement("GRCh38", "1", 109274968, "G", "T"), "GRCh37",
                              map_fn=None)
    assert (p37.pos, p37.ref, p37.alt, p37.strand_vs_grch38) == (109817590, "G", "T", 1)


def test_rsid_resolution_live():
    from snp_phewas_variant import resolve_variant

    v = resolve_variant("rs12740374")
    assert v.variant_id_grch38 == "1_109274968_G_T" and v.rs_id == "rs12740374"


def test_open_targets_resolver_returns_ranked_ids_live():
    from snp_phewas_ot_resolver import search_variant_ids

    assert search_variant_ids("rs12740374")[0] == "1_109274968_G_T"
    assert search_variant_ids("rs0") == []


def test_quick_mode_end_to_end_live(tmp_path):
    import snp_phewas

    assert snp_phewas.main(["--variant", "rs12740374", "--quick", "--output", str(tmp_path)]) == 0
    report = (tmp_path / "report.md").read_text()
    assert "Quick mode: FinnGen web API only" in report
    rows = (tmp_path / "tables" / "snp_phewas_rows.tsv").read_text().splitlines()
    assert len(rows) == 1 + 2469
    assert any("\tI9_CHD\t" in r and "\tok\t" in r for r in rows)
