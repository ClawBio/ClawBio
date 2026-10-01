"""Offline contract tests. All positive evidence below is synthetic, not clinical validation."""
import copy
import hashlib
import importlib.util
import json
import os
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[3]
SKILL = ROOT / "skills/gwas-prs"


def validator():
    spec = importlib.util.spec_from_file_location("prs_evidence", SKILL / "prs_evidence.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def case(tmp_path):
    score = tmp_path / "score.txt"
    score.write_text("rsID\teffect_allele\tother_allele\teffect_weight\nrs1\tA\tG\t0.2\n")
    inp = tmp_path / "input.txt"
    inp.write_text("# synthetic genotype\nrs1\t1\t1\tAG\n")
    evidence = {
        "schema_version": 1,
        "input": {"sha256": hashlib.sha256(inp.read_bytes()).hexdigest(),
                  "build": "GRCh37", "strand": "forward", "population": "synthetic-population",
                  "context": "synthetic-adults", "source": "synthetic input manifest"},
        "scores": {"TEST-SCORE": {
            "sha256": hashlib.sha256(score.read_bytes()).hexdigest(),
            "version": "1", "source": "synthetic score manifest", "build": "GRCh37",
            "variant_count": 1,
            "validation": {"population": "synthetic-population", "context": "synthetic-adults",
                "independent": True, "n": 1000, "source": "synthetic validation",
                "metric": "r2", "estimate": 0.03, "ci_lower": 0.01, "ci_upper": 0.05},
            "reference": {"population": "synthetic-population", "context": "synthetic-adults",
                "source": "synthetic reference", "n": 1000, "mean": 0.2, "sd": 0.1,
                "distribution": "normal", "missingness": "complete_only"},
            "heritability": {"estimate": 0.1, "source": "synthetic"},
            "causal_evidence": {"method": "MR", "source": "synthetic"}
        }}
    }
    return dict(score_id="TEST-SCORE", scoring_path=score, input_path=inp,
                build="GRCh37", variants=[{"rsid": "rs1", "effect_allele": "A",
                                          "other_allele": "G", "effect_weight": 0.2}],
                genotypes={"rs1": "AG"}, raw_score=0.2, evidence=evidence, curated=False)


def test_supported_percentile_is_not_disease_risk(case):
    result = validator().assess(**case)
    assert result["status"] == "supported"
    assert result["percentile"] == pytest.approx(50)
    assert result["risk_category"] is None
    assert result["causal_claim_supported"] is False
    assert result["predictive_performance"]["estimate"] == 0.03
    assert result["heritability"]["estimate"] == 0.1


@pytest.mark.parametrize("mutation,code", [
    ("missing", "EVIDENCE_MISSING"), ("score_hash", "SCORE_IDENTITY_MISMATCH"),
    ("input_hash", "INPUT_IDENTITY_MISMATCH"), ("build", "BUILD_MISMATCH"),
    ("population", "VALIDATION_CONTEXT_MISMATCH"), ("context", "VALIDATION_CONTEXT_MISMATCH"),
    ("independence", "VALIDATION_INCOMPLETE"), ("reference", "REFERENCE_INCOMPLETE"),
    ("reference_population", "REFERENCE_CONTEXT_MISMATCH"),
    ("nan", "REFERENCE_INCOMPLETE"), ("boolean", "REFERENCE_INCOMPLETE"),
    ("ci", "VALIDATION_INCOMPLETE"), ("no_calls", "ALLELE_INCOMPATIBLE"),
    ("allele", "ALLELE_INCOMPATIBLE"), ("missing_allele", "ALLELE_INCOMPATIBLE"),
    ("missing_variant", "INCOMPLETE_VARIANTS"), ("duplicate", "DUPLICATE_VARIANTS"),
    ("count", "VARIANT_COUNT_MISMATCH"), ("curated", "ILLUSTRATIVE_PANEL"),
])
def test_withholding_has_machine_readable_reasons(case, mutation, code):
    card = case["evidence"]["scores"]["TEST-SCORE"]
    if mutation == "missing": case["evidence"] = None
    elif mutation == "score_hash": card["sha256"] = "0" * 64
    elif mutation == "input_hash": case["evidence"]["input"]["sha256"] = "0" * 64
    elif mutation == "build": case["evidence"]["input"]["build"] = "GRCh38"
    elif mutation == "population": card["validation"]["population"] = "another population"
    elif mutation == "context": card["validation"]["context"] = "children"
    elif mutation == "independence": card["validation"]["independent"] = "true"
    elif mutation == "reference": del card["reference"]
    elif mutation == "reference_population": card["reference"]["population"] = "another population"
    elif mutation == "nan": card["reference"]["sd"] = float("nan")
    elif mutation == "boolean": card["reference"]["sd"] = True
    elif mutation == "ci": card["validation"]["ci_lower"] = 0.5
    elif mutation == "no_calls": case["genotypes"]["rs1"] = "00"
    elif mutation == "allele": case["genotypes"]["rs1"] = "CC"
    elif mutation == "missing_allele": del case["variants"][0]["other_allele"]
    elif mutation == "missing_variant": case["genotypes"] = {}
    elif mutation == "duplicate": case["variants"] *= 2
    elif mutation == "count": card["variant_count"] = 2
    elif mutation == "curated": case["curated"] = True
    result = validator().assess(**case)
    assert result["status"] != "supported"
    assert code in [r["code"] for r in result["reasons"]]
    assert result["percentile"] is None
    assert result["z_score"] is None
    assert result["risk_category"] is None


def test_deterministic_and_does_not_mutate(case):
    original = copy.deepcopy(case)
    module = validator()
    assert module.assess(**case) == module.assess(**case)
    assert case == original


@pytest.mark.parametrize("value", [[], "bad", {"schema_version": 2}, {"schema_version": 1, "scores": []}])
def test_malformed_evidence_fails_closed(case, value):
    case["evidence"] = value
    assert validator().assess(**case)["status"] != "supported"


def test_cli_evidence_survives_every_output_and_replay(tmp_path):
    evidence = tmp_path / "evidence.json"
    evidence.write_text('{"schema_version": 1, "input": {}, "scores": {}}')
    out = tmp_path / "report"
    run = subprocess.run([sys.executable, str(SKILL / "gwas_prs.py"), "--demo",
                          "--evidence-json", str(evidence), "--output", str(out)],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stderr + run.stdout
    rows = json.loads((out / "prs_results.json").read_text())
    assert rows and all(r["percentile"] is None for r in rows)
    assert all(r["evidence_assessment"]["status"] != "supported" for r in rows)
    envelope = json.loads((out / "result.json").read_text())
    assert envelope["summary"]["percentile"] is None
    assert "Evidence assessment" in (out / "prs_report.md").read_text()
    assert "--evidence-json" in (out / "reproducibility/commands.sh").read_text()
    provenance = json.loads((out / "reproducibility/provenance.json").read_text())
    assert provenance["parameters"]["evidence_sha256"] == hashlib.sha256(evidence.read_bytes()).hexdigest()
    assert str(evidence) not in (out / "reproducibility/provenance.json").read_text()

    replay_dir = tmp_path / "replayed"
    env = dict(os.environ, CLAWBIO_ROOT=str(ROOT),
               PRS_EVIDENCE_FILE=str(evidence), PRS_REPLAY_OUTPUT=str(replay_dir))
    env["PATH"] = str(Path(sys.executable).parent) + os.pathsep + env.get("PATH", "")
    replay = subprocess.run(["bash", str(out / "reproducibility/commands.sh")],
                            env=env, capture_output=True, text=True)
    assert replay.returncode == 0, replay.stdout + replay.stderr
    assert (replay_dir / "prs_results.json").read_bytes() == (out / "prs_results.json").read_bytes()

def test_raw_score_must_match_compatible_dosages(case):
    case["raw_score"] = 999
    result = validator().assess(**case)
    assert "RAW_SCORE_MISMATCH" in [r["code"] for r in result["reasons"]]
    assert result["percentile"] is None


def test_strand_must_be_declared(case):
    del case["evidence"]["input"]["strand"]
    result = validator().assess(**case)
    assert "STRAND_UNKNOWN" in [r["code"] for r in result["reasons"]]


@pytest.mark.parametrize("trait", ["Synthetic personality trait", "Synthetic metabolic trait"])
@pytest.mark.parametrize("provide_evidence", [True, False])
def test_real_input_path_supported_or_withheld(case, monkeypatch, tmp_path, trait, provide_evidence):
    sys.path.insert(0, str(SKILL))
    import gwas_prs
    score_id = "PGS999999"  # synthetic identifier; no API requests
    (tmp_path / (score_id + "_hmPOS_GRCh37.txt")).write_bytes(case["scoring_path"].read_bytes())
    inp = case["input_path"]
    inp.write_text("# This data file generated by 23andMe\n# rsid\tchromosome\tposition\tgenotype\nrs1\t1\t1\tAG\n")
    ev = case["evidence"]
    ev["input"]["sha256"] = hashlib.sha256(inp.read_bytes()).hexdigest()
    ev["scores"][score_id] = ev["scores"].pop("TEST-SCORE")
    evfile = tmp_path / "manifest.json"
    evfile.write_text(json.dumps(ev))
    monkeypatch.setattr(gwas_prs, "DATA_DIR", tmp_path)
    monkeypatch.setattr(gwas_prs.PGSCatalogClient, "get_score_metadata",
                        lambda *a, **kw: {"trait_reported": [trait], "variants_number": 1})
    out = tmp_path / "result"
    argv = ["gwas_prs.py", "--input", str(inp), "--pgs-id", score_id,
            "--output", str(out), "--cache-dir", str(tmp_path / "cache")]
    if provide_evidence:
        argv += ["--evidence-json", str(evfile)]
    monkeypatch.setattr(sys, "argv", argv)
    gwas_prs.main()
    rows = json.loads((out / "prs_results.json").read_text())
    assert rows[0]["trait"] == trait
    assert rows[0]["percentile"] == (50 if provide_evidence else None)
    assert rows[0]["risk_category"] is None
    assert rows[0]["evidence_assessment"]["status"] == ("supported" if provide_evidence else "not_established")

def test_offline_benchmark_runs_both_traits_and_exercises_abstention(tmp_path):
    run = subprocess.run([sys.executable, str(SKILL / "benchmark_evidence.py"),
                          "--demo", "--output", str(tmp_path / "benchmark")],
                         capture_output=True, text=True)
    assert run.returncode == 0, run.stdout + run.stderr
    result = json.loads((tmp_path / "benchmark/result.json").read_text())
    assert result["all_passed"]
    assert result["supported_cases"] >= 2
    assert result["withheld_cases"] >= 8
    assert result["false_supported"] == 0
    assert result["false_withheld"] == 0
    assert result["evidence_type"] == "synthetic_software_contract"
