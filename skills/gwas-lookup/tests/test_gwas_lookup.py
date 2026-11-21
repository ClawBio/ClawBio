"""
test_gwas_lookup.py — Automated test suite for GWAS Lookup skill.
Run with: pytest skills/gwas-lookup/tests/test_gwas_lookup.py -v

Uses pre-fetched JSON fixtures for all tests — no network required.
"""

import json
import sys
from pathlib import Path

import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

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


def test_report_labels_effect_estimates():
    """Each effect in the GWAS table must say which measure it is reported in."""
    from gwas_lookup_core.report import generate_markdown
    from gwas_lookup_core.normalise import merge_all

    demo = load_demo_data()
    merged = merge_all(demo["api_results"])

    report = generate_markdown(demo["variant"], merged)

    assert "| OR/Beta |" not in report
    assert "(OR)" in report
    assert "fine-mapped" in report


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
    write_reproducibility(tmp_path, "python skills/gwas-lookup/gwas_lookup.py --demo", [])
    assert (tmp_path / "reproducibility" / "commands.sh").exists()
    assert (tmp_path / "reproducibility" / "api_versions.json").exists()


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


def test_demo_writes_verifiable_reproducibility_bundle(tmp_path):
    """--demo writes commands.sh, environment.yml and checksums that all verify."""
    import hashlib
    import subprocess

    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(SKILL_DIR / "gwas_lookup.py"), "--demo", "--output", str(out)],
        check=True, capture_output=True,
    )

    repro = out / "reproducibility"
    commands = (repro / "commands.sh").read_text()
    assert "--demo" in commands
    assert "requests" in (repro / "environment.yml").read_text()

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "result.json", "reproducibility/environment.yml"} <= labels
    assert "tables/gwas_associations.csv" in labels
    for line in lines:
        digest, label = line.split("  ", 1)
        target = out / label
        assert target.resolve().is_relative_to(out.resolve()), label
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, label


def test_credible_set_total_reports_api_count_not_page():
    """The report must not present the first page as the whole set."""
    from gwas_lookup_core.normalise import merge_all

    credsets = load_fixture("open_targets_credsets")
    credsets["total_credible_sets"] = 226
    merged = merge_all({"open_targets_credsets": credsets})
    assert merged["summary"]["total_credible_sets"] == 226


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/gwas-lookup/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("gwas-lookup", tmp_path) == snapshot({'raw_results.json': """\
{
  "variant": {
    "rsid": "rs3798220",
    "chr": "6",
    "pos_grch38": 160540105,
    "pos_grch37": 160531442,
    "ref": "T",
    "alt": "C",
    "alt_alleles": [
      "C"
    ],
    "allele_string": "T/C",
    "var_class": "SNP",
    "most_severe_consequence": "missense_variant",
    "minor_allele": "C",
    "maf": 0.0178,
    "populations": [
      {
        "population": "1000GENOMES:phase_3:EUR",
        "allele": "C",
        "frequency": 0.0119
      },
      {
        "population": "1000GENOMES:phase_3:AFR",
        "allele": "C",
        "frequency": 0.0015
      },
      {
        "population": "1000GENOMES:phase_3:EAS",
        "allele": "C",
        "frequency": 0.0
      },
      {
        "population": "1000GENOMES:phase_3:AMR",
        "allele": "C",
        "frequency": 0.2176
      },
      {
        "population": "1000GENOMES:phase_3:SAS",
        "allele": "C",
        "frequency": 0.0041
      },
      {
        "population": "gnomAD:nfe",
        "allele": "C",
        "frequency": 0.0131
      },
      {
        "population": "gnomAD:afr",
        "allele": "C",
        "frequency": 0.0012
      },
      {
        "population": "gnomAD:eas",
        "allele": "C",
        "frequency": 0.0001
      },
      {
        "population": "gnomAD:amr",
        "allele": "C",
        "frequency": 0.0398
      },
      {
        "population": "gnomAD:fin",
        "allele": "C",
        "frequency": 0.0089
      }
    ],
    "variant_ids": {
      "open_targets": "6_160540105_T_C",
      "gtex": "chr6_160540105_T_C_b38",
      "pheweb": "6:160540105-T-C",
      "bbj": "6:160531442-T-C"
    },
    "ensembl_info": {
      "source": "ensembl_variation",
      "status": "ok",
      "rsid": "rs3798220"
    },
    "vep_info": {
      "source": "ensembl_vep",
      "status": "ok",
      "rsid": "rs3798220",
      "most_severe_consequence": "missense_variant",
      "consequences": [
        {
          "gene_symbol": "LPA",
          "gene_id": "ENSG00000198670",
          "consequence_terms": [
            "missense_variant"
          ],
          "impact": "MODERATE",
          "biotype": "protein_coding",
          "sift": "tolerated",
          "polyphen": "benign"
        },
        {
          "gene_symbol": "LPA",
          "gene_id": "ENSG00000198670",
          "consequence_terms": [
            "3_prime_UTR_variant"
          ],
          "impact": "MODIFIER",
          "biotype": "protein_coding",
          "sift": "",
          "polyphen": ""
        }
      ]
    }
  },
  "api_results": {
    "gwas_catalog": {
      "source": "gwas_catalog",
      "status": "ok",
      "rsid": "rs3798220",
      "total_associations": 8,
      "associations": [
        {
          "pvalue": 1e-45,
          "pvalue_mlog": 1,
          "pvalue_exponent": -45,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.02",
          "or_beta": 1.51,
          "beta_num": null,
          "beta_direction": null,
          "beta_unit": null,
          "ci": "[1.38-1.66]",
          "traits": [
            "Lipoprotein(a) levels"
          ],
          "study_accession": "GCST005140"
        },
        {
          "pvalue": 2e-38,
          "pvalue_mlog": 2,
          "pvalue_exponent": -38,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.015",
          "or_beta": 1.92,
          "beta_num": null,
          "beta_direction": null,
          "beta_unit": null,
          "ci": "[1.74-2.12]",
          "traits": [
            "Coronary heart disease"
          ],
          "study_accession": "GCST003116"
        },
        {
          "pvalue": 5e-22,
          "pvalue_mlog": 5,
          "pvalue_exponent": -22,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.013",
          "or_beta": null,
          "beta_num": 0.84,
          "beta_direction": "increase",
          "beta_unit": "unit increase",
          "ci": "[0.71-0.97]",
          "traits": [
            "Lipoprotein(a) measurement"
          ],
          "study_accession": "GCST90002370"
        },
        {
          "pvalue": 3e-18,
          "pvalue_mlog": 3,
          "pvalue_exponent": -18,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.02",
          "or_beta": 1.36,
          "beta_num": null,
          "beta_direction": null,
          "beta_unit": null,
          "ci": "[1.24-1.49]",
          "traits": [
            "Aortic valve stenosis"
          ],
          "study_accession": "GCST90038614"
        },
        {
          "pvalue": 8e-15,
          "pvalue_mlog": 8,
          "pvalue_exponent": -15,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.014",
          "or_beta": 1.27,
          "beta_num": null,
          "beta_direction": null,
          "beta_unit": null,
          "ci": "[1.19-1.36]",
          "traits": [
            "Peripheral artery disease"
          ],
          "study_accession": "GCST90018890"
        },
        {
          "pvalue": 1e-12,
          "pvalue_mlog": 1,
          "pvalue_exponent": -12,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.018",
          "or_beta": 1.22,
          "beta_num": null,
          "beta_direction": null,
          "beta_unit": null,
          "ci": "[1.15-1.29]",
          "traits": [
            "Ischemic stroke"
          ],
          "study_accession": "GCST006906"
        },
        {
          "pvalue": 4e-10,
          "pvalue_mlog": 4,
          "pvalue_exponent": -10,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.01",
          "or_beta": 1.18,
          "beta_num": null,
          "beta_direction": null,
          "beta_unit": null,
          "ci": "[1.12-1.25]",
          "traits": [
            "Heart failure"
          ],
          "study_accession": "GCST90043673"
        },
        {
          "pvalue": 2e-08,
          "pvalue_mlog": 2,
          "pvalue_exponent": -8,
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.012",
          "or_beta": 1.15,
          "beta_num": null,
          "beta_direction": null,
          "beta_unit": null,
          "ci": "[1.09-1.22]",
          "traits": [
            "Venous thromboembolism"
          ],
          "study_accession": "GCST90018997"
        }
      ]
    },
    "open_targets_credsets": {
      "source": "open_targets_credsets",
      "status": "ok",
      "variant_id": "6_160540105_T_C",
      "credible_sets": [
        {
          "study_id": "GCST005140",
          "trait": "Lipoprotein(a) levels",
          "posterior_probability": 0.92,
          "pval": 1e-45,
          "beta": 0.65,
          "is_95_credible": true,
          "is_99_credible": true
        },
        {
          "study_id": "GCST003116",
          "trait": "Coronary heart disease",
          "posterior_probability": 0.78,
          "pval": 2e-38,
          "beta": 0.31,
          "is_95_credible": true,
          "is_99_credible": true
        },
        {
          "study_id": "GCST90038614",
          "trait": "Aortic valve stenosis",
          "posterior_probability": 0.45,
          "pval": 3e-18,
          "beta": 0.22,
          "is_95_credible": false,
          "is_99_credible": true
        }
      ]
    },
    "pheweb_ukb": {
      "source": "pheweb_ukb",
      "status": "ok",
      "variant": "6:160540105-T-C",
      "total_associations": 5,
      "associations": [
        {
          "phenocode": "411.4",
          "phenostring": "Coronary atherosclerosis",
          "category": "circulatory system",
          "pval": 3.2e-25,
          "beta": 0.28,
          "sebeta": 0.027,
          "maf": 0.012,
          "num_cases": 34210,
          "num_controls": 386290
        },
        {
          "phenocode": "411",
          "phenostring": "Ischemic heart disease",
          "category": "circulatory system",
          "pval": 1.5e-20,
          "beta": 0.24,
          "sebeta": 0.026,
          "maf": 0.012,
          "num_cases": 41280,
          "num_controls": 379220
        },
        {
          "phenocode": "424.1",
          "phenostring": "Aortic valve disease",
          "category": "circulatory system",
          "pval": 8.7e-12,
          "beta": 0.35,
          "sebeta": 0.051,
          "maf": 0.012,
          "num_cases": 8420,
          "num_controls": 412080
        },
        {
          "phenocode": "443.9",
          "phenostring": "Peripheral vascular disease",
          "category": "circulatory system",
          "pval": 2.1e-08,
          "beta": 0.21,
          "sebeta": 0.037,
          "maf": 0.012,
          "num_cases": 12350,
          "num_controls": 408150
        },
        {
          "phenocode": "LPA",
          "phenostring": "Lipoprotein(a)",
          "category": "biomarkers",
          "pval": 1e-180,
          "beta": 1.85,
          "sebeta": 0.055,
          "maf": 0.012,
          "num_cases": null,
          "num_controls": null
        }
      ]
    },
    "finngen": {
      "source": "finngen",
      "status": "ok",
      "variant": "6:160540105-T-C",
      "total_associations": 3,
      "associations": [
        {
          "phenocode": "I9_CHD",
          "phenostring": "Coronary heart disease",
          "category": "Circulatory system",
          "pval": 4.5e-15,
          "beta": 0.31,
          "sebeta": 0.039,
          "maf": 0.009,
          "num_cases": 25400,
          "num_controls": 334100
        },
        {
          "phenocode": "I9_AORVALam",
          "phenostring": "Aortic valve stenosis/insufficiency",
          "category": "Circulatory system",
          "pval": 6.2e-10,
          "beta": 0.42,
          "sebeta": 0.068,
          "maf": 0.009,
          "num_cases": 5800,
          "num_controls": 353700
        },
        {
          "phenocode": "I9_PAD",
          "phenostring": "Peripheral artery disease",
          "category": "Circulatory system",
          "pval": 3.8e-06,
          "beta": 0.25,
          "sebeta": 0.054,
          "maf": 0.009,
          "num_cases": 7200,
          "num_controls": 352300
        }
      ]
    },
    "pheweb_bbj": {
      "source": "pheweb_bbj",
      "status": "ok",
      "variant": "6:160531442-T-C",
      "coordinate_system": "GRCh37",
      "total_associations": 1,
      "associations": [
        {
          "phenocode": "Angina_pectoris",
          "phenostring": "Angina pectoris",
          "category": "Circulatory",
          "pval": 0.032,
          "beta": 0.45,
          "sebeta": 0.21,
          "maf": 0.001,
          "num_cases": 4580,
          "num_controls": 176800
        }
      ]
    },
    "gtex": {
      "source": "gtex",
      "status": "ok",
      "variant_id": "chr6_160540105_T_C_b38",
      "total_eqtls": 3,
      "eqtls": [
        {
          "gene_symbol": "LPA",
          "gene_id": "ENSG00000198670.11",
          "tissue": "Liver",
          "pval": 2.3e-08,
          "nes": -0.82,
          "tissue_name": "Liver"
        },
        {
          "gene_symbol": "LPA",
          "gene_id": "ENSG00000198670.11",
          "tissue": "Adipose_Subcutaneous",
          "pval": 4.1e-05,
          "nes": -0.45,
          "tissue_name": "Adipose - Subcutaneous"
        },
        {
          "gene_symbol": "SLC22A3",
          "gene_id": "ENSG00000146477.13",
          "tissue": "Liver",
          "pval": 0.0012,
          "nes": 0.31,
          "tissue_name": "Liver"
        }
      ]
    }
  },
  "merged": {
    "gwas_associations": [
      {
        "source": "gwas_catalog",
        "trait": "Lipoprotein(a) levels",
        "pval": 1e-45,
        "or_beta": 1.51,
        "beta": null,
        "beta_direction": null,
        "ci": "[1.38-1.66]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.02",
        "study": "GCST005140",
        "genome_wide_significant": true
      },
      {
        "source": "open_targets",
        "trait": "Lipoprotein(a) levels",
        "pval": 1e-45,
        "or_beta": null,
        "beta": 0.65,
        "beta_direction": "",
        "ci": "",
        "risk_allele": "",
        "risk_frequency": "",
        "study": "GCST005140",
        "genome_wide_significant": true,
        "posterior_probability": 0.92,
        "is_95_credible": true
      },
      {
        "source": "gwas_catalog",
        "trait": "Coronary heart disease",
        "pval": 2e-38,
        "or_beta": 1.92,
        "beta": null,
        "beta_direction": null,
        "ci": "[1.74-2.12]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.015",
        "study": "GCST003116",
        "genome_wide_significant": true
      },
      {
        "source": "open_targets",
        "trait": "Coronary heart disease",
        "pval": 2e-38,
        "or_beta": null,
        "beta": 0.31,
        "beta_direction": "",
        "ci": "",
        "risk_allele": "",
        "risk_frequency": "",
        "study": "GCST003116",
        "genome_wide_significant": true,
        "posterior_probability": 0.78,
        "is_95_credible": true
      },
      {
        "source": "gwas_catalog",
        "trait": "Lipoprotein(a) measurement",
        "pval": 5e-22,
        "or_beta": null,
        "beta": 0.84,
        "beta_direction": "increase",
        "ci": "[0.71-0.97]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.013",
        "study": "GCST90002370",
        "genome_wide_significant": true
      },
      {
        "source": "gwas_catalog",
        "trait": "Aortic valve stenosis",
        "pval": 3e-18,
        "or_beta": 1.36,
        "beta": null,
        "beta_direction": null,
        "ci": "[1.24-1.49]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.02",
        "study": "GCST90038614",
        "genome_wide_significant": true
      },
      {
        "source": "open_targets",
        "trait": "Aortic valve stenosis",
        "pval": 3e-18,
        "or_beta": null,
        "beta": 0.22,
        "beta_direction": "",
        "ci": "",
        "risk_allele": "",
        "risk_frequency": "",
        "study": "GCST90038614",
        "genome_wide_significant": true,
        "posterior_probability": 0.45,
        "is_95_credible": false
      },
      {
        "source": "gwas_catalog",
        "trait": "Peripheral artery disease",
        "pval": 8e-15,
        "or_beta": 1.27,
        "beta": null,
        "beta_direction": null,
        "ci": "[1.19-1.36]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.014",
        "study": "GCST90018890",
        "genome_wide_significant": true
      },
      {
        "source": "gwas_catalog",
        "trait": "Ischemic stroke",
        "pval": 1e-12,
        "or_beta": 1.22,
        "beta": null,
        "beta_direction": null,
        "ci": "[1.15-1.29]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.018",
        "study": "GCST006906",
        "genome_wide_significant": true
      },
      {
        "source": "gwas_catalog",
        "trait": "Heart failure",
        "pval": 4e-10,
        "or_beta": 1.18,
        "beta": null,
        "beta_direction": null,
        "ci": "[1.12-1.25]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.01",
        "study": "GCST90043673",
        "genome_wide_significant": true
      },
      {
        "source": "gwas_catalog",
        "trait": "Venous thromboembolism",
        "pval": 2e-08,
        "or_beta": 1.15,
        "beta": null,
        "beta_direction": null,
        "ci": "[1.09-1.22]",
        "risk_allele": "rs3798220-C",
        "risk_frequency": "0.012",
        "study": "GCST90018997",
        "genome_wide_significant": true
      }
    ],
    "phewas": {
      "ukb": [
        {
          "source": "ukb",
          "phenocode": "LPA",
          "phenostring": "Lipoprotein(a)",
          "category": "biomarkers",
          "pval": 1e-180,
          "beta": 1.85,
          "se": 0.055,
          "maf": 0.012,
          "num_cases": null,
          "num_controls": null,
          "genome_wide_significant": true
        },
        {
          "source": "ukb",
          "phenocode": "411.4",
          "phenostring": "Coronary atherosclerosis",
          "category": "circulatory system",
          "pval": 3.2e-25,
          "beta": 0.28,
          "se": 0.027,
          "maf": 0.012,
          "num_cases": 34210,
          "num_controls": 386290,
          "genome_wide_significant": true
        },
        {
          "source": "ukb",
          "phenocode": "411",
          "phenostring": "Ischemic heart disease",
          "category": "circulatory system",
          "pval": 1.5e-20,
          "beta": 0.24,
          "se": 0.026,
          "maf": 0.012,
          "num_cases": 41280,
          "num_controls": 379220,
          "genome_wide_significant": true
        },
        {
          "source": "ukb",
          "phenocode": "424.1",
          "phenostring": "Aortic valve disease",
          "category": "circulatory system",
          "pval": 8.7e-12,
          "beta": 0.35,
          "se": 0.051,
          "maf": 0.012,
          "num_cases": 8420,
          "num_controls": 412080,
          "genome_wide_significant": true
        },
        {
          "source": "ukb",
          "phenocode": "443.9",
          "phenostring": "Peripheral vascular disease",
          "category": "circulatory system",
          "pval": 2.1e-08,
          "beta": 0.21,
          "se": 0.037,
          "maf": 0.012,
          "num_cases": 12350,
          "num_controls": 408150,
          "genome_wide_significant": true
        }
      ],
      "finngen": [
        {
          "source": "finngen",
          "phenocode": "I9_CHD",
          "phenostring": "Coronary heart disease",
          "category": "Circulatory system",
          "pval": 4.5e-15,
          "beta": 0.31,
          "se": 0.039,
          "maf": 0.009,
          "num_cases": 25400,
          "num_controls": 334100,
          "genome_wide_significant": true
        },
        {
          "source": "finngen",
          "phenocode": "I9_AORVALam",
          "phenostring": "Aortic valve stenosis/insufficiency",
          "category": "Circulatory system",
          "pval": 6.2e-10,
          "beta": 0.42,
          "se": 0.068,
          "maf": 0.009,
          "num_cases": 5800,
          "num_controls": 353700,
          "genome_wide_significant": true
        },
        {
          "source": "finngen",
          "phenocode": "I9_PAD",
          "phenostring": "Peripheral artery disease",
          "category": "Circulatory system",
          "pval": 3.8e-06,
          "beta": 0.25,
          "se": 0.054,
          "maf": 0.009,
          "num_cases": 7200,
          "num_controls": 352300,
          "genome_wide_significant": false
        }
      ],
      "bbj": [
        {
          "source": "bbj",
          "phenocode": "Angina_pectoris",
          "phenostring": "Angina pectoris",
          "category": "Circulatory",
          "pval": 0.032,
          "beta": 0.45,
          "se": 0.21,
          "maf": 0.001,
          "num_cases": 4580,
          "num_controls": 176800,
          "genome_wide_significant": false
        }
      ]
    },
    "eqtl_associations": [
      {
        "source": "gtex",
        "gene": "LPA",
        "gene_id": "ENSG00000198670.11",
        "tissue": "Liver",
        "tissue_name": "Liver",
        "pval": 2.3e-08,
        "effect_size": -0.82,
        "genome_wide_significant": true
      },
      {
        "source": "gtex",
        "gene": "LPA",
        "gene_id": "ENSG00000198670.11",
        "tissue": "Adipose_Subcutaneous",
        "tissue_name": "Adipose - Subcutaneous",
        "pval": 4.1e-05,
        "effect_size": -0.45,
        "genome_wide_significant": false
      },
      {
        "source": "gtex",
        "gene": "SLC22A3",
        "gene_id": "ENSG00000146477.13",
        "tissue": "Liver",
        "tissue_name": "Liver",
        "pval": 0.0012,
        "effect_size": 0.31,
        "genome_wide_significant": false
      }
    ],
    "credible_sets": [
      {
        "study_id": "GCST005140",
        "trait": "Lipoprotein(a) levels",
        "posterior_probability": 0.92,
        "pval": 1e-45,
        "beta": 0.65,
        "is_95_credible": true,
        "is_99_credible": true
      },
      {
        "study_id": "GCST003116",
        "trait": "Coronary heart disease",
        "posterior_probability": 0.78,
        "pval": 2e-38,
        "beta": 0.31,
        "is_95_credible": true,
        "is_99_credible": true
      },
      {
        "study_id": "GCST90038614",
        "trait": "Aortic valve stenosis",
        "posterior_probability": 0.45,
        "pval": 3e-18,
        "beta": 0.22,
        "is_95_credible": false,
        "is_99_credible": true
      }
    ],
    "data_sources": {
      "gwas_catalog": {
        "status": "ok",
        "message": ""
      },
      "open_targets_credsets": {
        "status": "ok",
        "message": ""
      },
      "pheweb_ukb": {
        "status": "ok",
        "message": ""
      },
      "finngen": {
        "status": "ok",
        "message": ""
      },
      "pheweb_bbj": {
        "status": "ok",
        "message": ""
      },
      "gtex": {
        "status": "ok",
        "message": ""
      }
    },
    "summary": {
      "total_gwas": 11,
      "total_gwas_significant": 11,
      "total_phewas_ukb": 5,
      "total_phewas_finngen": 3,
      "total_phewas_bbj": 1,
      "total_eqtls": 3,
      "total_credible_sets": 3
    }
  }
}\
""", 'report.md': """\
# GWAS Lookup Report: rs3798220

**Date**: <TIMESTAMP> UTC
**Variant**: rs3798220
**Location (GRCh38)**: chr6:160540105
**Location (GRCh37)**: chr6:160531442
**Alleles**: T/C
**Class**: SNP
**Consequence**: missense_variant
**Minor allele**: C (MAF: 0.0178)

---

## Variant Effect Prediction (VEP)

| Gene | Consequence | Impact | SIFT | PolyPhen |
|------|------------|--------|------|----------|
| LPA | missense_variant | MODERATE | tolerated | benign |
| LPA | 3_prime_UTR_variant | MODIFIER |  |  |

## GWAS Associations (11 total, 11 genome-wide significant)

| Trait | P-value | Effect (type) | Risk Allele | Source | Study |
|-------|---------|---------------|-------------|--------|-------|
| Lipoprotein(a) levels ** | 1.00e-45 | 1.510 (OR) | rs3798220-C | gwas_catalog | GCST005140 |
| Lipoprotein(a) levels ** | 1.00e-45 | 0.650 (beta, fine-mapped) |  | open_targets | GCST005140 |
| Coronary heart disease ** | 2.00e-38 | 1.920 (OR) | rs3798220-C | gwas_catalog | GCST003116 |
| Coronary heart disease ** | 2.00e-38 | 0.310 (beta, fine-mapped) |  | open_targets | GCST003116 |
| Lipoprotein(a) measurement ** | 5.00e-22 | 0.840 (beta) | rs3798220-C | gwas_catalog | GCST90002370 |
| Aortic valve stenosis ** | 3.00e-18 | 1.360 (OR) | rs3798220-C | gwas_catalog | GCST90038614 |
| Aortic valve stenosis ** | 3.00e-18 | 0.220 (beta, fine-mapped) |  | open_targets | GCST90038614 |
| Peripheral artery disease ** | 8.00e-15 | 1.270 (OR) | rs3798220-C | gwas_catalog | GCST90018890 |
| Ischemic stroke ** | 1.00e-12 | 1.220 (OR) | rs3798220-C | gwas_catalog | GCST006906 |
| Heart failure ** | 4.00e-10 | 1.180 (OR) | rs3798220-C | gwas_catalog | GCST90043673 |
| Venous thromboembolism ** | 2.00e-08 | 1.150 (OR) | rs3798220-C | gwas_catalog | GCST90018997 |

*Open Targets rows are fine-mapped credible sets. Their beta is a posterior estimate for this variant, not a GWAS Catalog effect size, and the query does not request an effect allele.*

## PheWAS: UKB-TOPMed (5 phenotypes, 5 significant)

| Phenotype | Category | P-value | Beta | MAF |
|-----------|----------|---------|------|-----|
| Lipoprotein(a) ** | biomarkers | 1.00e-180 | 1.850 | 0.0120 |
| Coronary atherosclerosis ** | circulatory system | 3.20e-25 | 0.280 | 0.0120 |
| Ischemic heart disease ** | circulatory system | 1.50e-20 | 0.240 | 0.0120 |
| Aortic valve disease ** | circulatory system | 8.70e-12 | 0.350 | 0.0120 |
| Peripheral vascular disease ** | circulatory system | 2.10e-08 | 0.210 | 0.0120 |

## PheWAS: FinnGen (3 phenotypes, 2 significant)

| Phenotype | Category | P-value | Beta | MAF |
|-----------|----------|---------|------|-----|
| Coronary heart disease ** | Circulatory system | 4.50e-15 | 0.310 | 0.0090 |
| Aortic valve stenosis/insufficiency ** | Circulatory system | 6.20e-10 | 0.420 | 0.0090 |
| Peripheral artery disease | Circulatory system | 3.80e-06 | 0.250 | 0.0090 |

## PheWAS: Biobank Japan (1 phenotypes, 0 significant)

| Phenotype | Category | P-value | Beta | MAF |
|-----------|----------|---------|------|-----|
| Angina pectoris | Circulatory | 3.20e-02 | 0.450 | 0.0010 |

## eQTL Associations (3 total)

| Gene | Tissue | P-value | Effect Size | Source |
|------|--------|---------|-------------|--------|
| LPA | Liver | 2.30e-08 | -0.820 | gtex |
| LPA | Adipose - Subcutaneous | 4.10e-05 | -0.450 | gtex |
| SLC22A3 | Liver | 1.20e-03 | 0.310 | gtex |

## Fine-Mapping / Credible Sets (3 sets)

| Trait | Study | Post. Prob. | P-value | 95% CS | 99% CS |
|-------|-------|-------------|---------|--------|--------|
| Lipoprotein(a) levels | GCST005140 | 0.920 | 1.00e-45 | Yes | Yes |
| Coronary heart disease | GCST003116 | 0.780 | 2.00e-38 | Yes | Yes |
| Aortic valve stenosis | GCST90038614 | 0.450 | 3.00e-18 | No | Yes |

## Data Sources

| Source | Status |
|--------|--------|
| gwas_catalog | OK |
| open_targets_credsets | OK |
| pheweb_ukb | OK |
| finngen | OK |
| pheweb_bbj | OK |
| gtex | OK |

## Methods

This report was generated by the ClawBio GWAS Lookup skill, which queries
8 genomic databases in parallel: Ensembl, GWAS Catalog, Open Targets,
UKB-TOPMed PheWeb, FinnGen, Biobank Japan PheWeb, GTEx, and LocusZoom PortalDev.

- **Genome-wide significance threshold**: p < 5 x 10^-8
- **Coordinate system**: GRCh38 (primary), GRCh37 for BBJ
- **Cache**: 24-hour local file cache

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
""", 'result.json': """\
{
  "skill": "gwas-lookup",
  "version": "0.2.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "",
  "datasets": {},
  "summary": {
    "rsid": "rs3798220",
    "chr": "6",
    "pos_grch38": 160540105,
    "total_gwas": 11,
    "total_gwas_significant": 11,
    "total_phewas_ukb": 5,
    "total_phewas_finngen": 3,
    "total_phewas_bbj": 1,
    "total_eqtls": 3,
    "total_credible_sets": 3,
    "apis_queried": 6,
    "apis_skipped": 0,
    "apis_failed": {}
  },
  "data": {
    "variant": {
      "rsid": "rs3798220",
      "chr": "6",
      "pos_grch38": 160540105,
      "pos_grch37": 160531442,
      "ref": "T",
      "alt": "C",
      "alt_alleles": [
        "C"
      ],
      "allele_string": "T/C",
      "var_class": "SNP",
      "most_severe_consequence": "missense_variant",
      "minor_allele": "C",
      "maf": 0.0178,
      "populations": [
        {
          "population": "1000GENOMES:phase_3:EUR",
          "allele": "C",
          "frequency": 0.0119
        },
        {
          "population": "1000GENOMES:phase_3:AFR",
          "allele": "C",
          "frequency": 0.0015
        },
        {
          "population": "1000GENOMES:phase_3:EAS",
          "allele": "C",
          "frequency": 0.0
        },
        {
          "population": "1000GENOMES:phase_3:AMR",
          "allele": "C",
          "frequency": 0.2176
        },
        {
          "population": "1000GENOMES:phase_3:SAS",
          "allele": "C",
          "frequency": 0.0041
        },
        {
          "population": "gnomAD:nfe",
          "allele": "C",
          "frequency": 0.0131
        },
        {
          "population": "gnomAD:afr",
          "allele": "C",
          "frequency": 0.0012
        },
        {
          "population": "gnomAD:eas",
          "allele": "C",
          "frequency": 0.0001
        },
        {
          "population": "gnomAD:amr",
          "allele": "C",
          "frequency": 0.0398
        },
        {
          "population": "gnomAD:fin",
          "allele": "C",
          "frequency": 0.0089
        }
      ],
      "variant_ids": {
        "open_targets": "6_160540105_T_C",
        "gtex": "chr6_160540105_T_C_b38",
        "pheweb": "6:160540105-T-C",
        "bbj": "6:160531442-T-C"
      },
      "ensembl_info": {
        "source": "ensembl_variation",
        "status": "ok",
        "rsid": "rs3798220"
      },
      "vep_info": {
        "source": "ensembl_vep",
        "status": "ok",
        "rsid": "rs3798220",
        "most_severe_consequence": "missense_variant",
        "consequences": [
          {
            "gene_symbol": "LPA",
            "gene_id": "ENSG00000198670",
            "consequence_terms": [
              "missense_variant"
            ],
            "impact": "MODERATE",
            "biotype": "protein_coding",
            "sift": "tolerated",
            "polyphen": "benign"
          },
          {
            "gene_symbol": "LPA",
            "gene_id": "ENSG00000198670",
            "consequence_terms": [
              "3_prime_UTR_variant"
            ],
            "impact": "MODIFIER",
            "biotype": "protein_coding",
            "sift": "",
            "polyphen": ""
          }
        ]
      }
    },
    "merged": {
      "gwas_associations": [
        {
          "source": "gwas_catalog",
          "trait": "Lipoprotein(a) levels",
          "pval": 1e-45,
          "or_beta": 1.51,
          "beta": null,
          "beta_direction": null,
          "ci": "[1.38-1.66]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.02",
          "study": "GCST005140",
          "genome_wide_significant": true
        },
        {
          "source": "open_targets",
          "trait": "Lipoprotein(a) levels",
          "pval": 1e-45,
          "or_beta": null,
          "beta": 0.65,
          "beta_direction": "",
          "ci": "",
          "risk_allele": "",
          "risk_frequency": "",
          "study": "GCST005140",
          "genome_wide_significant": true,
          "posterior_probability": 0.92,
          "is_95_credible": true
        },
        {
          "source": "gwas_catalog",
          "trait": "Coronary heart disease",
          "pval": 2e-38,
          "or_beta": 1.92,
          "beta": null,
          "beta_direction": null,
          "ci": "[1.74-2.12]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.015",
          "study": "GCST003116",
          "genome_wide_significant": true
        },
        {
          "source": "open_targets",
          "trait": "Coronary heart disease",
          "pval": 2e-38,
          "or_beta": null,
          "beta": 0.31,
          "beta_direction": "",
          "ci": "",
          "risk_allele": "",
          "risk_frequency": "",
          "study": "GCST003116",
          "genome_wide_significant": true,
          "posterior_probability": 0.78,
          "is_95_credible": true
        },
        {
          "source": "gwas_catalog",
          "trait": "Lipoprotein(a) measurement",
          "pval": 5e-22,
          "or_beta": null,
          "beta": 0.84,
          "beta_direction": "increase",
          "ci": "[0.71-0.97]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.013",
          "study": "GCST90002370",
          "genome_wide_significant": true
        },
        {
          "source": "gwas_catalog",
          "trait": "Aortic valve stenosis",
          "pval": 3e-18,
          "or_beta": 1.36,
          "beta": null,
          "beta_direction": null,
          "ci": "[1.24-1.49]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.02",
          "study": "GCST90038614",
          "genome_wide_significant": true
        },
        {
          "source": "open_targets",
          "trait": "Aortic valve stenosis",
          "pval": 3e-18,
          "or_beta": null,
          "beta": 0.22,
          "beta_direction": "",
          "ci": "",
          "risk_allele": "",
          "risk_frequency": "",
          "study": "GCST90038614",
          "genome_wide_significant": true,
          "posterior_probability": 0.45,
          "is_95_credible": false
        },
        {
          "source": "gwas_catalog",
          "trait": "Peripheral artery disease",
          "pval": 8e-15,
          "or_beta": 1.27,
          "beta": null,
          "beta_direction": null,
          "ci": "[1.19-1.36]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.014",
          "study": "GCST90018890",
          "genome_wide_significant": true
        },
        {
          "source": "gwas_catalog",
          "trait": "Ischemic stroke",
          "pval": 1e-12,
          "or_beta": 1.22,
          "beta": null,
          "beta_direction": null,
          "ci": "[1.15-1.29]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.018",
          "study": "GCST006906",
          "genome_wide_significant": true
        },
        {
          "source": "gwas_catalog",
          "trait": "Heart failure",
          "pval": 4e-10,
          "or_beta": 1.18,
          "beta": null,
          "beta_direction": null,
          "ci": "[1.12-1.25]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.01",
          "study": "GCST90043673",
          "genome_wide_significant": true
        },
        {
          "source": "gwas_catalog",
          "trait": "Venous thromboembolism",
          "pval": 2e-08,
          "or_beta": 1.15,
          "beta": null,
          "beta_direction": null,
          "ci": "[1.09-1.22]",
          "risk_allele": "rs3798220-C",
          "risk_frequency": "0.012",
          "study": "GCST90018997",
          "genome_wide_significant": true
        }
      ],
      "phewas": {
        "ukb": [
          {
            "source": "ukb",
            "phenocode": "LPA",
            "phenostring": "Lipoprotein(a)",
            "category": "biomarkers",
            "pval": 1e-180,
            "beta": 1.85,
            "se": 0.055,
            "maf": 0.012,
            "num_cases": null,
            "num_controls": null,
            "genome_wide_significant": true
          },
          {
            "source": "ukb",
            "phenocode": "411.4",
            "phenostring": "Coronary atherosclerosis",
            "category": "circulatory system",
            "pval": 3.2e-25,
            "beta": 0.28,
            "se": 0.027,
            "maf": 0.012,
            "num_cases": 34210,
            "num_controls": 386290,
            "genome_wide_significant": true
          },
          {
            "source": "ukb",
            "phenocode": "411",
            "phenostring": "Ischemic heart disease",
            "category": "circulatory system",
            "pval": 1.5e-20,
            "beta": 0.24,
            "se": 0.026,
            "maf": 0.012,
            "num_cases": 41280,
            "num_controls": 379220,
            "genome_wide_significant": true
          },
          {
            "source": "ukb",
            "phenocode": "424.1",
            "phenostring": "Aortic valve disease",
            "category": "circulatory system",
            "pval": 8.7e-12,
            "beta": 0.35,
            "se": 0.051,
            "maf": 0.012,
            "num_cases": 8420,
            "num_controls": 412080,
            "genome_wide_significant": true
          },
          {
            "source": "ukb",
            "phenocode": "443.9",
            "phenostring": "Peripheral vascular disease",
            "category": "circulatory system",
            "pval": 2.1e-08,
            "beta": 0.21,
            "se": 0.037,
            "maf": 0.012,
            "num_cases": 12350,
            "num_controls": 408150,
            "genome_wide_significant": true
          }
        ],
        "finngen": [
          {
            "source": "finngen",
            "phenocode": "I9_CHD",
            "phenostring": "Coronary heart disease",
            "category": "Circulatory system",
            "pval": 4.5e-15,
            "beta": 0.31,
            "se": 0.039,
            "maf": 0.009,
            "num_cases": 25400,
            "num_controls": 334100,
            "genome_wide_significant": true
          },
          {
            "source": "finngen",
            "phenocode": "I9_AORVALam",
            "phenostring": "Aortic valve stenosis/insufficiency",
            "category": "Circulatory system",
            "pval": 6.2e-10,
            "beta": 0.42,
            "se": 0.068,
            "maf": 0.009,
            "num_cases": 5800,
            "num_controls": 353700,
            "genome_wide_significant": true
          },
          {
            "source": "finngen",
            "phenocode": "I9_PAD",
            "phenostring": "Peripheral artery disease",
            "category": "Circulatory system",
            "pval": 3.8e-06,
            "beta": 0.25,
            "se": 0.054,
            "maf": 0.009,
            "num_cases": 7200,
            "num_controls": 352300,
            "genome_wide_significant": false
          }
        ],
        "bbj": [
          {
            "source": "bbj",
            "phenocode": "Angina_pectoris",
            "phenostring": "Angina pectoris",
            "category": "Circulatory",
            "pval": 0.032,
            "beta": 0.45,
            "se": 0.21,
            "maf": 0.001,
            "num_cases": 4580,
            "num_controls": 176800,
            "genome_wide_significant": false
          }
        ]
      },
      "eqtl_associations": [
        {
          "source": "gtex",
          "gene": "LPA",
          "gene_id": "ENSG00000198670.11",
          "tissue": "Liver",
          "tissue_name": "Liver",
          "pval": 2.3e-08,
          "effect_size": -0.82,
          "genome_wide_significant": true
        },
        {
          "source": "gtex",
          "gene": "LPA",
          "gene_id": "ENSG00000198670.11",
          "tissue": "Adipose_Subcutaneous",
          "tissue_name": "Adipose - Subcutaneous",
          "pval": 4.1e-05,
          "effect_size": -0.45,
          "genome_wide_significant": false
        },
        {
          "source": "gtex",
          "gene": "SLC22A3",
          "gene_id": "ENSG00000146477.13",
          "tissue": "Liver",
          "tissue_name": "Liver",
          "pval": 0.0012,
          "effect_size": 0.31,
          "genome_wide_significant": false
        }
      ],
      "credible_sets": [
        {
          "study_id": "GCST005140",
          "trait": "Lipoprotein(a) levels",
          "posterior_probability": 0.92,
          "pval": 1e-45,
          "beta": 0.65,
          "is_95_credible": true,
          "is_99_credible": true
        },
        {
          "study_id": "GCST003116",
          "trait": "Coronary heart disease",
          "posterior_probability": 0.78,
          "pval": 2e-38,
          "beta": 0.31,
          "is_95_credible": true,
          "is_99_credible": true
        },
        {
          "study_id": "GCST90038614",
          "trait": "Aortic valve stenosis",
          "posterior_probability": 0.45,
          "pval": 3e-18,
          "beta": 0.22,
          "is_95_credible": false,
          "is_99_credible": true
        }
      ],
      "data_sources": {
        "gwas_catalog": {
          "status": "ok",
          "message": ""
        },
        "open_targets_credsets": {
          "status": "ok",
          "message": ""
        },
        "pheweb_ukb": {
          "status": "ok",
          "message": ""
        },
        "finngen": {
          "status": "ok",
          "message": ""
        },
        "pheweb_bbj": {
          "status": "ok",
          "message": ""
        },
        "gtex": {
          "status": "ok",
          "message": ""
        }
      },
      "summary": {
        "total_gwas": 11,
        "total_gwas_significant": 11,
        "total_phewas_ukb": 5,
        "total_phewas_finngen": 3,
        "total_phewas_bbj": 1,
        "total_eqtls": 3,
        "total_credible_sets": 3
      }
    }
  }
}\
""", 'tables/credible_sets.csv': """\
study_id,trait,posterior_probability,pval,beta,is_95_credible,is_99_credible
GCST005140,Lipoprotein(a) levels,0.92,1e-45,0.65,True,True
GCST003116,Coronary heart disease,0.78,2e-38,0.31,True,True
GCST90038614,Aortic valve stenosis,0.45,3e-18,0.22,False,True
""", 'tables/eqtl_associations.csv': """\
source,gene,gene_id,tissue,tissue_name,pval,effect_size,genome_wide_significant
gtex,LPA,ENSG00000198670.11,Liver,Liver,2.3e-08,-0.82,True
gtex,LPA,ENSG00000198670.11,Adipose_Subcutaneous,Adipose - Subcutaneous,4.1e-05,-0.45,False
gtex,SLC22A3,ENSG00000146477.13,Liver,Liver,0.0012,0.31,False
""", 'tables/gwas_associations.csv': """\
source,trait,pval,or_beta,beta,beta_direction,ci,risk_allele,risk_frequency,study,genome_wide_significant,posterior_probability,is_95_credible
gwas_catalog,Lipoprotein(a) levels,1e-45,1.51,,,[1.38-1.66],rs3798220-C,0.02,GCST005140,True,,
open_targets,Lipoprotein(a) levels,1e-45,,0.65,,,,,GCST005140,True,0.92,True
gwas_catalog,Coronary heart disease,2e-38,1.92,,,[1.74-2.12],rs3798220-C,0.015,GCST003116,True,,
open_targets,Coronary heart disease,2e-38,,0.31,,,,,GCST003116,True,0.78,True
gwas_catalog,Lipoprotein(a) measurement,5e-22,,0.84,increase,[0.71-0.97],rs3798220-C,0.013,GCST90002370,True,,
gwas_catalog,Aortic valve stenosis,3e-18,1.36,,,[1.24-1.49],rs3798220-C,0.02,GCST90038614,True,,
open_targets,Aortic valve stenosis,3e-18,,0.22,,,,,GCST90038614,True,0.45,False
gwas_catalog,Peripheral artery disease,8e-15,1.27,,,[1.19-1.36],rs3798220-C,0.014,GCST90018890,True,,
gwas_catalog,Ischemic stroke,1e-12,1.22,,,[1.15-1.29],rs3798220-C,0.018,GCST006906,True,,
gwas_catalog,Heart failure,4e-10,1.18,,,[1.12-1.25],rs3798220-C,0.01,GCST90043673,True,,
gwas_catalog,Venous thromboembolism,2e-08,1.15,,,[1.09-1.22],rs3798220-C,0.012,GCST90018997,True,,
""", 'tables/phewas_bbj.csv': """\
source,phenocode,phenostring,category,pval,beta,se,maf,num_cases,num_controls,genome_wide_significant
bbj,Angina_pectoris,Angina pectoris,Circulatory,0.032,0.45,0.21,0.001,4580,176800,False
""", 'tables/phewas_finngen.csv': """\
source,phenocode,phenostring,category,pval,beta,se,maf,num_cases,num_controls,genome_wide_significant
finngen,I9_CHD,Coronary heart disease,Circulatory system,4.5e-15,0.31,0.039,0.009,25400,334100,True
finngen,I9_AORVALam,Aortic valve stenosis/insufficiency,Circulatory system,6.2e-10,0.42,0.068,0.009,5800,353700,True
finngen,I9_PAD,Peripheral artery disease,Circulatory system,3.8e-06,0.25,0.054,0.009,7200,352300,False
""", 'tables/phewas_ukb.csv': """\
source,phenocode,phenostring,category,pval,beta,se,maf,num_cases,num_controls,genome_wide_significant
ukb,LPA,Lipoprotein(a),biomarkers,1e-180,1.85,0.055,0.012,,,True
ukb,411.4,Coronary atherosclerosis,circulatory system,3.2e-25,0.28,0.027,0.012,34210,386290,True
ukb,411,Ischemic heart disease,circulatory system,1.5e-20,0.24,0.026,0.012,41280,379220,True
ukb,424.1,Aortic valve disease,circulatory system,8.7e-12,0.35,0.051,0.012,8420,412080,True
ukb,443.9,Peripheral vascular disease,circulatory system,2.1e-08,0.21,0.037,0.012,12350,408150,True
"""})
