"""Synthetic evidence contracts; no patient data or clinical validation."""
import copy
import importlib.util
import itertools
import json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parent.parent
spec = importlib.util.spec_from_file_location("severity_under_test", ROOT / "severity_evidence.py")
severity = importlib.util.module_from_spec(spec)
spec.loader.exec_module(severity)

def bundle():
    return {"schema_version": "1.0", "hpo_release": "2025-01-16",
            "upstream_commit": "03bc5f6be6456a3ca2c5206e399f0a38f879ce57",
            "records": [{"gene": "SYNTHETIC", "disease_id": "MONDO:0000001",
                         "disease": "Synthetic condition", "hpo_id": "HP:0000001",
                         "frequency": "30%", "qualifier": "", "tier": 1,
                         "citation": "https://example.org/synthetic-evidence",
                         "source_span": "Synthetic Tier 1 evidence, not a clinical assertion.",
                         "source_version": "fixture-v1", "hpo_release": "2025-01-16"}]}

def evaluate(b):
    return severity.evaluate(b, "SYNTHETIC", "MONDO:0000001")

def test_valid_and_provenance():
    b = bundle()
    r = evaluate(b)
    assert r["severity"] == "Severe"
    assert r["hpo_release"] == b["hpo_release"]
    assert len(r["evidence_sha256"]) == 64
    assert r["upstream_commit"] == b["upstream_commit"]
    assert r["evidence"][0]["source_span"] == b["records"][0]["source_span"]
    assert b == bundle()

@pytest.mark.parametrize("field", ["gene", "disease_id", "disease", "hpo_id",
 "frequency", "qualifier", "tier", "citation", "source_span", "source_version", "hpo_release"])
def test_missing_record_field_abstains(field):
    b = bundle()
    del b["records"][0][field]
    assert evaluate(b)["severity"] == "severity_unknown"

@pytest.mark.parametrize("field", ["schema_version", "hpo_release", "upstream_commit", "records"])
def test_missing_envelope_abstains(field):
    b = bundle()
    del b[field]
    assert evaluate(b)["severity"] == "severity_unknown"

@pytest.mark.parametrize("frequency", [None, "", "NaN", float("nan"), True, -1, 1.1,
                                       "8-Mar", "3/8/2024", "4/0", "HP:9999999"])
def test_ambiguous_frequency_abstains(frequency):
    b = bundle()
    b["records"][0]["frequency"] = frequency
    assert evaluate(b)["severity"] == "severity_unknown"

@pytest.mark.parametrize("frequency,expected", [
    (0.3, "Severe"), ("3/10", "Severe"), ("30%", "Severe"),
    ("HP:0040280", "Severe"), ("HP:0040281", "Severe"), ("HP:0040282", "Severe"),
    (0.299, "severity_unknown"), ("29%", "severity_unknown"),
    ("HP:0040283", "severity_unknown"), ("HP:0040284", "severity_unknown"),
    ("HP:0040285", "severity_unknown")])
def test_frequency_filter(frequency, expected):
    b = bundle()
    b["records"][0]["frequency"] = frequency
    assert evaluate(b)["severity"] == expected

def test_negation_and_duplicates():
    b = bundle()
    b["records"][0]["qualifier"] = " NOT "
    assert evaluate(b)["severity"] == "severity_unknown"
    b = bundle()
    b["records"].append(copy.deepcopy(b["records"][0]))
    assert evaluate(b)["severity"] == "severity_unknown"

@pytest.mark.parametrize("field,value", [("tier", True), ("tier", 0), ("tier", 5),
 ("tier", "1"), ("hpo_id", "x"), ("disease_id", "x"), ("source_version", "latest"),
 ("hpo_release", "2024-01-01"), ("qualifier", "maybe"), ("citation", "made up")])
def test_invalid_record_abstains(field, value):
    b = bundle()
    b["records"][0][field] = value
    assert evaluate(b)["severity"] == "severity_unknown"

@pytest.mark.parametrize("value", [None, [], "bad", 2, {"records": [None]}])
def test_malformed_bundle(value):
    assert evaluate(value)["severity"] == "severity_unknown"

def test_incompatible_versions():
    for field, value in [("schema_version", "2"), ("upstream_commit", "main"),
                         ("hpo_release", "latest"), ("hpo_release", "2025-99-99")]:
        b = bundle()
        b[field] = value
        assert evaluate(b)["severity"] == "severity_unknown"

def test_no_gene_only_join():
    assert severity.evaluate(bundle(), "SYNTHETIC", None)["severity"] == "severity_unknown"
    assert severity.evaluate(bundle(), "OTHER", "MONDO:0000001")["severity"] == "severity_unknown"

def test_all_filtered_is_unknown_not_mild():
    b = bundle()
    b["records"][0]["tier"] = 4
    assert evaluate(b)["severity"] == "Mild"
    b["records"][0]["frequency"] = 0.1
    assert evaluate(b)["severity"] == "severity_unknown"

def test_upstream_parity():
    spec = importlib.util.spec_from_file_location("upstream_reference", ROOT / "tests/fixtures/upstream_aggregation.py")
    ref = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(ref)
    for counts in itertools.product(range(6), repeat=3):
        row = dict(zip(("1", "2", "3"), counts))
        assert severity.aggregate(row) == ref.calculate_severity(row)

def test_api_optional_and_noninterfering(tmp_path):
    import sys
    sys.path.insert(0, str(ROOT))
    api_spec = importlib.util.spec_from_file_location("severity_api_under_test", ROOT / "api.py")
    api_module = importlib.util.module_from_spec(api_spec)
    api_spec.loader.exec_module(api_module)
    run = api_module.run
    entry = {"id": "synthetic-variant", "gene": "SYNTHETIC", "disease_id": "MONDO:0000001",
             "condition": "Synthetic condition", "condition_severity": "mild",
             "ref_allele": "A", "alt_allele": "G", "inheritance": "ar",
             "clinvar_significance": "Pathogenic", "panels": ["synthetic"]}
    panel = tmp_path / "panel.json"
    panel.write_text(json.dumps([entry]))
    base = run({"synthetic-variant": "AG"}, {"panel_path": panel})
    assert "severity_evidence" not in base
    changed = run({"synthetic-variant": "AG"}, {"panel_path": panel, "severity_evidence": bundle()})
    assert changed["findings"][0]["disease_severity"]["severity"] == "Severe"
    assert "severity_disagreement" in changed["findings"][0]["disease_severity"]["warnings"]
    clean = copy.deepcopy(changed)
    del clean["severity_evidence"]
    del clean["findings"][0]["disease_severity"]
    assert clean == base

def test_duplicate_panel_ids_abstain():
    rows = [{"id": "same", "gene": "SYNTHETIC", "disease_id": "MONDO:0000001"},
            {"id": "same", "gene": "OTHER", "disease_id": "MONDO:0000002"}]
    findings = [{"id": "same", "gene": "SYNTHETIC"}, {"id": "same", "gene": "OTHER"}]
    severity.annotate(findings, rows, bundle())
    assert all(f["disease_severity"]["severity"] == "severity_unknown" for f in findings)
    assert all("ambiguous_panel_identity" in f["disease_severity"]["reasons"] for f in findings)
