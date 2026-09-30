"""
test_gwas_lookup.py — Automated test suite for GWAS Lookup skill.
Run with: pytest skills/gwas-lookup/tests/test_gwas_lookup.py -v

Uses pre-fetched JSON fixtures for all tests — no network required.
"""

import json
import sys
from pathlib import Path

import pytest

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
    from gwas_lookup_core.normalise import merge_gwas

    gwas_catalog = load_fixture("gwas_catalog")
    credsets = load_fixture("open_targets_credsets")
    merged = merge_gwas(gwas_catalog, credsets)

    assert len(merged) > 0
    pvals = [a["pval"] for a in merged if a["pval"] is not None]
    assert pvals == sorted(pvals), "GWAS associations should be sorted by p-value"


def test_merge_gwas_includes_both_sources():
    """Merged GWAS should include entries from both GWAS Catalog and Open Targets."""
    from gwas_lookup_core.normalise import merge_gwas

    gwas_catalog = load_fixture("gwas_catalog")
    credsets = load_fixture("open_targets_credsets")
    merged = merge_gwas(gwas_catalog, credsets)

    sources = {a["source"] for a in merged}
    assert "gwas_catalog" in sources
    assert "open_targets" in sources


def test_merge_gwas_flags_significant():
    """Genome-wide significant hits (p < 5e-8) should be flagged."""
    from gwas_lookup_core.normalise import merge_gwas

    gwas_catalog = load_fixture("gwas_catalog")
    credsets = load_fixture("open_targets_credsets")
    merged = merge_gwas(gwas_catalog, credsets)

    significant = [a for a in merged if a.get("genome_wide_significant")]
    assert len(significant) > 0, "Should have genome-wide significant hits"

    for a in significant:
        assert a["pval"] < 5e-8


def test_merge_phewas_structure():
    """PheWAS merge should return dict with ukb, finngen, bbj keys."""
    from gwas_lookup_core.normalise import merge_phewas

    ukb = load_fixture("pheweb_ukb")
    finngen = {"source": "finngen", "status": "ok", "associations": []}
    bbj = {"source": "pheweb_bbj", "status": "ok", "associations": []}

    result = merge_phewas(ukb, finngen, bbj)
    assert "ukb" in result
    assert "finngen" in result
    assert "bbj" in result
    assert len(result["ukb"]) == 2  # 2 UKB associations in fixture


def test_merge_eqtls():
    """eQTL merge should normalise GTEx results."""
    from gwas_lookup_core.normalise import merge_eqtls

    merged = merge_eqtls(load_fixture("gtex"))

    assert len(merged) == 2
    assert {e["source"] for e in merged} == {"gtex"}


def test_merge_all_structure():
    """merge_all should produce the expected top-level keys."""
    from gwas_lookup_core.normalise import merge_all

    api_results = {
        "gwas_catalog": load_fixture("gwas_catalog"),
        "open_targets_credsets": load_fixture("open_targets_credsets"),
        "pheweb_ukb": load_fixture("pheweb_ukb"),
        "finngen": {"source": "finngen", "status": "ok", "associations": []},
        "pheweb_bbj": {"source": "pheweb_bbj", "status": "ok", "associations": []},
        "gtex": load_fixture("gtex"),
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


def test_resolve_variant_imports_api_modules(monkeypatch):
    """resolve_variant should use the sibling api package without package-parent assumptions."""
    from gwas_lookup_core.resolve import resolve_variant
    from gwas_lookup_api import ensembl, portaldev

    monkeypatch.setattr(
        ensembl,
        "get_variant_info",
        lambda rsid, cache_dir=None, use_cache=True: {
            "status": "ok",
            "chr": "6",
            "pos_grch38": 160540105,
            "pos_grch37": 161005610,
            "ref_allele": "T",
            "alt_alleles": ["C"],
            "allele_string": "T/C",
            "var_class": "SNV",
            "most_severe_consequence": "missense_variant",
            "minor_allele": "C",
            "maf": 0.12,
            "populations": [],
        },
    )
    monkeypatch.setattr(
        ensembl,
        "get_vep_annotation",
        lambda rsid, cache_dir=None, use_cache=True: {
            "status": "ok",
            "consequences": [],
        },
    )
    monkeypatch.setattr(
        portaldev,
        "resolve_rsid",
        lambda rsid, cache_dir=None, use_cache=True: {"status": "error"},
    )

    resolved = resolve_variant("rs3798220")

    assert resolved["rsid"] == "rs3798220"
    assert resolved["variant_ids"]["open_targets"] == "6_160540105_T_C"
    assert resolved["variant_ids"]["bbj"] == "6:161005610-T-C"


# ── Graceful degradation ─────────────────────────────────────────────────────


def test_merge_with_error_api():
    """Report should still generate when one API returns an error."""
    from gwas_lookup_core.normalise import merge_all

    api_results = {
        "gwas_catalog": load_fixture("gwas_catalog"),
        "open_targets_credsets": load_fixture("open_targets_credsets"),
        "pheweb_ukb": load_fixture("pheweb_ukb"),
        "finngen": load_fixture("error_api"),  # error
        "pheweb_bbj": {"source": "pheweb_bbj", "status": "error", "message": "404"},
        "gtex": load_fixture("gtex"),
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
    from gwas_lookup_core.normalise import merge_all

    api_results = {
        "gwas_catalog": {"source": "gwas_catalog", "status": "error", "message": "timeout"},
        "open_targets_credsets": {"source": "open_targets_credsets", "status": "error", "message": "timeout"},
        "pheweb_ukb": {"source": "pheweb_ukb", "status": "error", "message": "timeout"},
        "finngen": {"source": "finngen", "status": "error", "message": "timeout"},
        "pheweb_bbj": {"source": "pheweb_bbj", "status": "error", "message": "timeout"},
        "gtex": {"source": "gtex", "status": "error", "message": "timeout"},
    }

    merged = merge_all(api_results)
    assert merged["summary"]["total_gwas"] == 0
    assert merged["summary"]["total_eqtls"] == 0


# ── Report generation ─────────────────────────────────────────────────────────


def test_report_includes_disclaimer():
    """Report markdown should include the ClawBio disclaimer."""
    from gwas_lookup_core.report import generate_markdown

    demo = load_demo_data()
    from gwas_lookup_core.normalise import merge_all
    merged = merge_all(demo["api_results"])

    report = generate_markdown(demo["variant"], merged)
    assert "research and educational tool" in report
    assert "not a medical device" in report


def test_report_includes_variant_info():
    """Report should include variant rsID, coordinates, and consequence."""
    from gwas_lookup_core.report import generate_markdown
    from gwas_lookup_core.normalise import merge_all

    demo = load_demo_data()
    merged = merge_all(demo["api_results"])

    report = generate_markdown(demo["variant"], merged)
    assert "rs3798220" in report
    assert "160540105" in report
    assert "missense_variant" in report


def test_report_includes_gwas_table():
    """Report should include a GWAS associations table."""
    from gwas_lookup_core.report import generate_markdown
    from gwas_lookup_core.normalise import merge_all

    demo = load_demo_data()
    merged = merge_all(demo["api_results"])

    report = generate_markdown(demo["variant"], merged)
    assert "GWAS Associations" in report
    assert "Lipoprotein" in report


def test_write_tables(tmp_path):
    """CSV tables should be written to the output directory."""
    from gwas_lookup_core.report import write_tables
    from gwas_lookup_core.normalise import merge_all

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
    """Demo data should include results from every queried API."""
    demo = load_demo_data()
    expected = [
        "gwas_catalog", "open_targets_credsets",
        "pheweb_ukb", "finngen", "pheweb_bbj", "gtex",
    ]
    for src in expected:
        assert src in demo["api_results"], f"Missing demo data for {src}"
        assert demo["api_results"][src]["status"] == "ok", f"{src} should be ok in demo"


def test_demo_full_pipeline(tmp_path):
    """Full pipeline should run with demo data and produce report.md."""
    from gwas_lookup_core.normalise import merge_all
    from gwas_lookup_core.report import generate_markdown, write_tables, write_reproducibility

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


def test_gwas_catalog_forwards_max_hits_as_page_size(monkeypatch):
    """GWAS Catalog pages at 20 unless size is sent with the request."""
    from gwas_lookup_api import gwas_catalog

    captured = {}

    class FakeClient:
        def get(self, endpoint, params=None):
            captured["endpoint"] = endpoint
            captured["params"] = params or {}
            page_size = captured["params"].get("size", 20)
            return {
                "_embedded": {
                    "associations": [
                        {"pvalue": i, "efoTraits": [], "riskAlleles": []}
                        for i in range(30)
                    ][:page_size]
                }
            }

    monkeypatch.setattr(gwas_catalog, "_make_client", lambda *args, **kwargs: FakeClient())
    result = gwas_catalog.get_associations("rs1", max_hits=25)
    assert captured["endpoint"] == "singleNucleotidePolymorphisms/rs1/associations"
    assert captured["params"] == {"size": 25}
    assert result["status"] == "ok"
    assert result["total_associations"] == 25
    assert len(result["associations"]) == 25


# ── Open Targets Platform schema ──────────────────────────────────────────────


def test_open_targets_credsets_parse_platform_schema(monkeypatch):
    """Credible sets come from variant.credibleSets.rows with this variant's locus stats."""
    from gwas_lookup_api import open_targets

    class FakeClient:
        def post(self, endpoint, json_body, params=None):
            return load_fixture("open_targets_platform_credsets")

    monkeypatch.setattr(open_targets, "_make_client", lambda *a, **k: FakeClient())
    result = open_targets.get_credible_sets("6", 160540105, "T", "C")

    assert result["status"] == "ok"
    assert result["total_credible_sets"] == 226  # count from the API, not the page
    first = result["credible_sets"][0]
    assert first["study_id"] == "GCST90498998"
    assert first["trait"] == "Total lipids in medium VLDL"
    assert first["pval"] == pytest.approx(6.561e-127, rel=1e-3)
    assert first["posterior_probability"] == pytest.approx(0.8493, rel=1e-3)
    assert first["beta"] == pytest.approx(-0.18998)
    assert first["is_95_credible"] is True


def test_open_targets_graphql_errors_are_errors_not_empty(monkeypatch):
    """GraphQL reports schema errors in the body; they must not read as 'no data'."""
    from gwas_lookup_api import open_targets

    class FakeClient:
        def post(self, endpoint, json_body, params=None):
            return {"data": None, "errors": [{"message": "Cannot query field 'rsId' on type 'Variant'."}]}

    monkeypatch.setattr(open_targets, "_make_client", lambda *a, **k: FakeClient())
    result = open_targets.get_credible_sets("6", 160540105, "T", "C")

    assert result["status"] == "error"
    assert "Cannot query field" in result["message"]


# ── Failed sources are surfaced ───────────────────────────────────────────────


def test_failed_sources_listed_in_result_json(tmp_path):
    """A source that errors must be named in result.json, not only printed."""
    import gwas_lookup

    demo = load_demo_data()
    demo["api_results"]["gtex"] = {"source": "gtex", "status": "error", "message": "HTTP 410 Gone"}
    gwas_lookup.run_lookup("rs3798220", tmp_path, make_figures=False, demo_data=demo)

    summary = json.loads((tmp_path / "result.json").read_text())["summary"]
    assert summary["apis_failed"] == {"gtex": "HTTP 410 Gone"}


def test_credible_set_total_reports_api_count_not_page():
    """The report must not present the first page as the whole set."""
    from gwas_lookup_core.normalise import merge_all

    credsets = load_fixture("open_targets_credsets")
    credsets["total_credible_sets"] = 226
    merged = merge_all({"open_targets_credsets": credsets})
    assert merged["summary"]["total_credible_sets"] == 226
