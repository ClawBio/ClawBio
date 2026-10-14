"""Safety regressions use synthetic marker calls, never private patient data."""

import html
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest
from hypothesis import given, strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))
import pharmgx_reporter as reporter
from clawbio.cli import format_pharmgx_preview
from clawbio.common.checksums import sha256_file


def pgx(calls):
    return {
        rsid: {"genotype": genotype, **reporter.PGX_SNPS[rsid]}
        for rsid, genotype in calls.items()
    }


def normal_profiles():
    return {
        gene: {"diplotype": "*1/*1", "phenotype": "Normal Metabolizer"}
        for gene in reporter.GENE_DEFS
    }


@pytest.mark.parametrize("genotype", ["II", "DD", "NN", "A", "AX", ""])
def test_unsupported_genotype_never_becomes_reference(genotype):
    diplotype = reporter.call_diplotype("CYP2D6", pgx({"rs5030655": genotype}))
    assert "Indeterminate" in reporter.call_phenotype("CYP2D6", diplotype)
    assert "*1/*1" not in diplotype


def test_detected_variant_does_not_hide_missing_coverage():
    diplotype = reporter.call_diplotype("CYP2D6", pgx({"rs3892097": "CT"}))
    assert "Indeterminate" in reporter.call_phenotype("CYP2D6", diplotype)
    assert "1/5" in diplotype


@given(st.sampled_from(["AA", "AG", "GA", "GG"]))
def test_verified_cyp3a5_strand_conversion_is_invariant(genotype):
    coding = pgx({"rs776746": genotype, "rs10264272": "GG"})
    opposite = pgx({
        "rs776746": reporter.flip_genotype(genotype),
        "rs10264272": "GG",
    })
    assert reporter.normalize_marker_genotype("rs776746", genotype) == reporter.normalize_marker_genotype(
        "rs776746", reporter.flip_genotype(genotype)
    )
    assert reporter.call_diplotype("CYP3A5", coding) == reporter.call_diplotype("CYP3A5", opposite)


@pytest.mark.parametrize("genotype", ["TT", "II", "DD"])
def test_ugt_repeat_proxy_cannot_establish_normal_function(genotype):
    diplotype = reporter.call_diplotype(
        "UGT1A1", pgx({"rs8175347": genotype, "rs4148323": "GG"})
    )
    assert "Indeterminate" in reporter.call_phenotype("UGT1A1", diplotype)
    assert "repeat" in diplotype.lower()


def test_cyp1a2_does_not_claim_clinical_ultrarapid_phenotype():
    diplotype = reporter.call_diplotype(
        "CYP1A2", pgx({"rs762551": "AA", "rs2069514": "GG"})
    )
    phenotype = reporter.call_phenotype("CYP1A2", diplotype)
    assert "Indeterminate" in phenotype
    assert "Ultrarapid" not in phenotype


def test_missing_mthfr_marker_cannot_establish_normal_activity():
    diplotype = reporter.call_diplotype("MTHFR", pgx({"rs1801131": "TT"}))
    assert "Indeterminate" in reporter.call_phenotype("MTHFR", diplotype)


@pytest.mark.parametrize("phenotype", ["", "Normal mystery", "Unknown", "Indeterminate"])
def test_unrecognised_phenotype_never_defaults_to_normal(phenotype):
    assert reporter.phenotype_to_key(phenotype) == "indeterminate"


def test_nudt15_marker_has_correct_star_allele():
    # CPIC 2025-update supplement: rs147390019 defines *4, not *2.
    assert reporter.PGX_SNPS["rs147390019"]["allele"] == "*4"
    assert reporter.GENE_DEFS["NUDT15"]["variants"]["rs147390019"]["allele"] == "*4"


@pytest.mark.parametrize("drug", ["Azathioprine", "Mercaptopurine", "Thioguanine"])
@pytest.mark.parametrize("nudt15", ["Normal Metabolizer", "Intermediate Metabolizer",
                                  "Poor Metabolizer", "Indeterminate (incomplete coverage)"])
def test_thiopurines_consider_both_genes_and_match_single_lookup(drug, nudt15):
    profiles = normal_profiles()
    profiles["NUDT15"]["phenotype"] = nudt15
    full = next(d for bucket in reporter.lookup_drugs(profiles).values()
                for d in bucket if d["drug"] == drug)
    single = reporter.lookup_single_drug(drug, profiles)
    expected = {
        "Normal Metabolizer": "standard",
        "Intermediate Metabolizer": "caution",
        "Poor Metabolizer": "caution",
        "Indeterminate (incomplete coverage)": "indeterminate",
    }[nudt15]
    assert full["classification"] == single["classification"] == expected
    assert "NUDT15" in full["gene"]
    assert full.get("note")


def test_thiopurine_risk_is_not_hidden_by_other_gene_uncertainty():
    profiles = normal_profiles()
    profiles["TPMT"]["phenotype"] = "Indeterminate"
    profiles["NUDT15"]["phenotype"] = "Poor Metabolizer"
    result = reporter.lookup_single_drug("Azathioprine", profiles)
    assert result["classification"] == "caution"
    assert "TPMT" in result["note"] and "NUDT15" in result["note"]


@pytest.mark.parametrize("drug", ["Clozapine", "Methotrexate"])
def test_no_prescribing_evidence_is_not_a_genotype_dose_alert(drug):
    result = reporter.lookup_single_drug(drug, normal_profiles())
    assert result["classification"] == "indeterminate"
    assert result["reason"] == "no_prescribing_recommendation"
    assert "genotype" in result["note"].lower()
    assert result["source"]
    card = reporter.format_dosage_card(result, visible_dose="50mg")
    assert "50mg may need adjustment" not in card


@pytest.mark.parametrize("phenotype", ["Normal Metabolizer", "Intermediate Metabolizer",
                                      "Poor Metabolizer"])
def test_methadone_uses_standard_prescribing_not_genotype_caution(phenotype):
    profiles = normal_profiles()
    profiles["CYP2B6"]["phenotype"] = phenotype
    result = reporter.lookup_single_drug("Methadone", profiles)
    assert result["classification"] == "standard"
    assert "standard prescribing" in result["note"].lower()
    assert result["source"]


def test_unresolved_phenotype_has_no_extracted_patient_recommendation():
    enrichment = {"codeine": {"_guidelines": [{"textMarkdown": "unused"}]}}
    results = {"indeterminate": [{"drug": "Codeine", "gene": "CYP2D6"}]}
    profiles = {"CYP2D6": {"phenotype": "Indeterminate (incomplete coverage)"}}
    from unittest.mock import patch
    with patch("clawbio.common.rec_shortener.extract_all_source_recs") as extract:
        reporter.extract_phenotype_recs(enrichment, results, profiles)
    extract.assert_not_called()


def test_terminal_preserves_phenotypes_and_input_identity(capsys):
    report = """**Input**: `synthetic.txt`
| Gene | Full Name | Diplotype | Phenotype |
|------|-----------|-----------|-----------|
| CYP2C19 | Enzyme | *1/*17 | Rapid Metabolizer |
| CYP2D6 | Enzyme | Indeterminate | Indeterminate (incomplete coverage) |
| HLA-B | HLA | Not assessed | Indeterminate (not in panel) |
"""
    format_pharmgx_preview(report, "report.md")
    output = capsys.readouterr().out
    assert "Rapid Metabolizer" in output
    assert "Indeterminate (incomplete coverage)" in output
    assert "synthetic.txt" in output
    assert "Corpasome" not in output


def test_terminal_does_not_invent_warfarin_genotypes_or_doses(capsys):
    report = """### Actionable Alerts
**AVOID / USE ALTERNATIVE:**
- **Warfarin** (Coumadin) [CYP2C9+VKORC1]: Clinical review required.
"""
    format_pharmgx_preview(report, "report.md")
    output = capsys.readouterr().out
    assert "50–80%" not in output
    assert "Both copies" not in output
    assert "*1/*2" not in output


def test_reports_disclose_missing_markers_and_limitations_consistently(tmp_path):
    source = tmp_path / "synthetic.txt"
    source.write_text("# rsid\tchromosome\tposition\tgenotype\nrs3892097\t22\t42524947\tCT\n")
    markers = pgx({"rs3892097": "CT"})
    profiles = normal_profiles()
    profiles["CYP2D6"] = {
        "diplotype": reporter.call_diplotype("CYP2D6", markers),
        "phenotype": "Indeterminate (incomplete coverage)",
    }
    results = reporter.lookup_drugs(profiles)
    markdown = reporter.generate_report(source, "23andme", 1, markers, profiles, results)
    webpage = html.unescape(reporter.generate_html_report(
        source, "23andme", 1, markers, profiles, results
    ))
    for text in [markdown, webpage]:
        assert "rs5030655" in text
        assert "rs28371725" in text
        assert "Copy number" in text
        assert "UGT1A1" in text and "repeat" in text.lower()
        assert "HLA-B" in text and "MT-RNR1" in text
        assert "unmapped diplotypes" not in text
        assert "No high-risk gene-drug interactions detected" not in text
        assert "not proof" in text.lower()
    unassessed = markdown.split("## Genes Not Assessed by This Panel")[1].split("##")[0]
    assert "CYP2B6" not in unassessed
    assert f"{len(reporter.PGX_SNPS)} pharmacogenomic variants across {len(reporter.GENE_DEFS)} genes" in markdown
    assert "Variant allele" in markdown and "Variant allele" in webpage


def test_cli_api_and_json_have_same_calls_and_coverage(tmp_path):
    source = tmp_path / "synthetic.txt"
    source.write_text("# rsid\tchromosome\tposition\tgenotype\nrs4244285\t10\t96541616\tGG\n"
                      "rs4148323\t2\t234669144\tGG\n")
    output = tmp_path / "out"
    completed = subprocess.run([
        sys.executable, str(SKILL_DIR / "pharmgx_reporter.py"), "--input", str(source),
        "--output", str(output), "--no-enrich",
    ], capture_output=True, text=True)
    assert completed.returncode == 0, completed.stderr
    result = json.loads((output / "result.json").read_text())
    spec = importlib.util.spec_from_file_location("pharmgx_safety_api", SKILL_DIR / "api.py")
    api = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(api)
    api_result = api.run({"rs4244285": "GG", "rs4148323": "GG"})
    assert result["data"]["gene_profiles"] == api_result["gene_profiles"]
    assert "marker_coverage" in result["data"]
    assert result["data"]["marker_coverage"] == api_result["marker_coverage"]
    assert result["input_checksum"] == "sha256:" + sha256_file(source)


def test_standard_dosage_card_does_not_promise_efficacy():
    result = reporter.lookup_single_drug("Clopidogrel", normal_profiles())
    card = reporter.format_dosage_card(result, visible_dose="75mg")
    assert "expected to be effective" not in card
    assert "supports Clopidogrel 75mg" not in card


def test_incomplete_qualifier_cannot_be_stripped_to_a_normal_phenotype():
    assert reporter.phenotype_to_key("Normal Metabolizer (incomplete coverage)") == "indeterminate"


def test_warfarin_sensitivity_is_a_dosing_review_not_a_blanket_contraindication():
    profiles = {
        "CYP2C9": {"phenotype": "Normal Metabolizer"},
        "VKORC1": {"phenotype": "High Warfarin Sensitivity"},
    }
    classification, note = reporter.get_warfarin_rec(profiles)
    assert classification == "caution"
    assert "INR" in note and "algorithm" in note.lower()


def test_reproduction_command_preserves_paths_and_offline_option(tmp_path):
    source = tmp_path / "synthetic input.txt"
    source.write_text("synthetic data\n")
    output = tmp_path / "report"
    path = reporter.write_commands_sh(output, source, no_enrich=True)
    script = path.read_text()
    assert str(source.resolve()) in script
    assert str(SKILL_DIR / "pharmgx_reporter.py") in script
    assert "--no-enrich" in script
    assert "already exists" in script
    assert subprocess.run(["bash", "-n", str(path)], capture_output=True).returncode == 0
