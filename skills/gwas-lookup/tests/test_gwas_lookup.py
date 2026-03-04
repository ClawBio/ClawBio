"""
test_gwas_lookup.py — Automated test suite for GWAS Lookup skill.
Run with: pytest skills/gwas-lookup/tests/test_gwas_lookup.py -v

Uses pre-fetched JSON fixtures for all tests — no network required.
"""

import json
import sys
from pathlib import Path

# Add parent dir to path so we can import the skill modules
SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

FIXTURES = Path(__file__).parent / "fixtures"
DEMO_DATA = SKILL_DIR / "data" / "demo_rs3798220.json"


def load_fixture(name: str) -> dict:
    path = FIXTURES / f"{name}.json"
    return json.loads(path.read_text())


def load_demo_data() -> dict:
    return json.loads(DEMO_DATA.read_text())


# ── Normalisation ─────────────────────────────────────────────────────────────


def test_merge_gwas_sorts_by_pval():
    """GWAS associations should be sorted by p-value ascending."""
    from core.normalise import merge_gwas

    gwas_catalog = load_fixture("gwas_catalog")
    credsets = load_fixture("open_targets_credsets")
    merged = merge_gwas(gwas_catalog, credsets)

    assert len(merged) > 0
    pvals = [a["pval"] for a in merged if a["pval"] is not None]
    assert pvals == sorted(pvals), "GWAS associations should be sorted by p-value"


def test_merge_gwas_includes_both_sources():
    """Merged GWAS should include entries from both GWAS Catalog and Open Targets."""
    from core.normalise import merge_gwas

    gwas_catalog = load_fixture("gwas_catalog")
    credsets = load_fixture("open_targets_credsets")
    merged = merge_gwas(gwas_catalog, credsets)

    sources = {a["source"] for a in merged}
    assert "gwas_catalog" in sources
    assert "open_targets" in sources


def test_merge_gwas_flags_significant():
    """Genome-wide significant hits (p < 5e-8) should be flagged."""
    from core.normalise import merge_gwas

    gwas_catalog = load_fixture("gwas_catalog")
    credsets = load_fixture("open_targets_credsets")
    merged = merge_gwas(gwas_catalog, credsets)

    significant = [a for a in merged if a.get("genome_wide_significant")]
    assert len(significant) > 0, "Should have genome-wide significant hits"

    for a in significant:
        assert a["pval"] < 5e-8


def test_merge_phewas_structure():
    """PheWAS merge should return dict with ukb, finngen, bbj keys."""
    from core.normalise import merge_phewas

    ukb = load_fixture("pheweb_ukb")
    finngen = {"source": "finngen", "status": "ok", "associations": []}
    bbj = {"source": "pheweb_bbj", "status": "ok", "associations": []}

    result = merge_phewas(ukb, finngen, bbj)
    assert "ukb" in result
    assert "finngen" in result
    assert "bbj" in result
    assert len(result["ukb"]) == 2  # 2 UKB associations in fixture


def test_merge_eqtls():
    """eQTL merge should combine GTEx and eQTL Catalogue results."""
    from core.normalise import merge_eqtls

    gtex = load_fixture("gtex")
    eqtl_cat = load_fixture("eqtl_catalogue")
    merged = merge_eqtls(gtex, eqtl_cat)

    assert len(merged) == 3  # 2 GTEx + 1 eQTL Catalogue
    sources = {e["source"] for e in merged}
    assert "gtex" in sources
    assert "eqtl_catalogue" in sources


def test_merge_all_structure():
    """merge_all should produce the expected top-level keys."""
    from core.normalise import merge_all

    api_results = {
        "gwas_catalog": load_fixture("gwas_catalog"),
        "open_targets_credsets": load_fixture("open_targets_credsets"),
        "pheweb_ukb": load_fixture("pheweb_ukb"),
        "finngen": {"source": "finngen", "status": "ok", "associations": []},
        "pheweb_bbj": {"source": "pheweb_bbj", "status": "ok", "associations": []},
        "gtex": load_fixture("gtex"),
        "eqtl_catalogue": load_fixture("eqtl_catalogue"),
    }

    merged = merge_all(api_results)
    assert "gwas_associations" in merged
    assert "phewas" in merged
    assert "eqtl_associations" in merged
    assert "credible_sets" in merged
    assert "data_sources" in merged
    assert "summary" in merged

    summary = merged["summary"]
    assert summary["total_gwas"] > 0
    assert summary["total_eqtls"] > 0


# ── Graceful degradation ─────────────────────────────────────────────────────


def test_merge_with_error_api():
    """Report should still generate when one API returns an error."""
    from core.normalise import merge_all

    api_results = {
        "gwas_catalog": load_fixture("gwas_catalog"),
        "open_targets_credsets": load_fixture("open_targets_credsets"),
        "pheweb_ukb": load_fixture("pheweb_ukb"),
        "finngen": load_fixture("error_api"),  # error
        "pheweb_bbj": {"source": "pheweb_bbj", "status": "error", "message": "404"},
        "gtex": load_fixture("gtex"),
        "eqtl_catalogue": load_fixture("eqtl_catalogue"),
    }

    merged = merge_all(api_results)
    # Should still have GWAS and eQTL results despite PheWAS errors
    assert merged["summary"]["total_gwas"] > 0
    assert merged["summary"]["total_eqtls"] > 0
    # FinnGen and BBJ should show as error
    assert merged["data_sources"]["finngen"]["status"] == "error"
    assert merged["data_sources"]["pheweb_bbj"]["status"] == "error"


def test_merge_with_all_errors():
    """merge_all should produce a valid structure even if all APIs fail."""
    from core.normalise import merge_all

    api_results = {
        "gwas_catalog": {"source": "gwas_catalog", "status": "error", "message": "timeout"},
        "open_targets_credsets": {"source": "open_targets_credsets", "status": "error", "message": "timeout"},
        "pheweb_ukb": {"source": "pheweb_ukb", "status": "error", "message": "timeout"},
        "finngen": {"source": "finngen", "status": "error", "message": "timeout"},
        "pheweb_bbj": {"source": "pheweb_bbj", "status": "error", "message": "timeout"},
        "gtex": {"source": "gtex", "status": "error", "message": "timeout"},
        "eqtl_catalogue": {"source": "eqtl_catalogue", "status": "error", "message": "timeout"},
    }

    merged = merge_all(api_results)
    assert merged["summary"]["total_gwas"] == 0
    assert merged["summary"]["total_eqtls"] == 0


# ── Report generation ─────────────────────────────────────────────────────────


def test_report_includes_disclaimer():
    """Report markdown should include the ClawBio disclaimer."""
    from core.report import generate_markdown

    demo = load_demo_data()
    from core.normalise import merge_all
    merged = merge_all(demo["api_results"])

    report = generate_markdown(demo["variant"], merged)
    assert "research and educational tool" in report
    assert "not a medical device" in report


def test_report_includes_variant_info():
    """Report should include variant rsID, coordinates, and consequence."""
    from core.report import generate_markdown
    from core.normalise import merge_all

    demo = load_demo_data()
    merged = merge_all(demo["api_results"])

    report = generate_markdown(demo["variant"], merged)
    assert "rs3798220" in report
    assert "160540105" in report
    assert "missense_variant" in report


def test_report_includes_gwas_table():
    """Report should include a GWAS associations table."""
    from core.report import generate_markdown
    from core.normalise import merge_all

    demo = load_demo_data()
    merged = merge_all(demo["api_results"])

    report = generate_markdown(demo["variant"], merged)
    assert "GWAS Associations" in report
    assert "Lipoprotein" in report


def test_write_tables(tmp_path):
    """CSV tables should be written to the output directory."""
    from core.report import write_tables
    from core.normalise import merge_all

    demo = load_demo_data()
    merged = merge_all(demo["api_results"])

    write_tables(tmp_path, merged)
    tables_dir = tmp_path / "tables"
    assert tables_dir.exists()
    assert (tables_dir / "gwas_associations.csv").exists()


# ── Demo mode ─────────────────────────────────────────────────────────────────


def test_demo_data_loads():
    """Demo data file should load and contain expected structure."""
    demo = load_demo_data()
    assert "variant" in demo
    assert "api_results" in demo
    assert demo["variant"]["rsid"] == "rs3798220"
    assert demo["variant"]["chr"] == "6"
    assert "gwas_catalog" in demo["api_results"]


def test_demo_data_has_all_sources():
    """Demo data should include results from all 8 API modules."""
    demo = load_demo_data()
    expected = [
        "gwas_catalog", "open_targets", "open_targets_credsets",
        "pheweb_ukb", "finngen", "pheweb_bbj", "gtex", "eqtl_catalogue",
    ]
    for src in expected:
        assert src in demo["api_results"], f"Missing demo data for {src}"
        assert demo["api_results"][src]["status"] == "ok", f"{src} should be ok in demo"


def test_demo_full_pipeline(tmp_path):
    """Full pipeline should run with demo data and produce report.md."""
    from core.normalise import merge_all
    from core.report import generate_markdown, write_tables, write_reproducibility

    demo = load_demo_data()
    variant = demo["variant"]
    merged = merge_all(demo["api_results"])

    # Write report
    report = generate_markdown(variant, merged)
    (tmp_path / "report.md").write_text(report)
    assert (tmp_path / "report.md").exists()

    # Write tables
    write_tables(tmp_path, merged)
    assert (tmp_path / "tables" / "gwas_associations.csv").exists()

    # Write reproducibility
    write_reproducibility(tmp_path, variant, [])
    assert (tmp_path / "reproducibility" / "commands.sh").exists()
    assert (tmp_path / "reproducibility" / "api_versions.json").exists()
