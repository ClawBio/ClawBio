"""Tests for drug-repurposing-screen skill."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest
from inline_snapshot import snapshot
from clawbio.common import demo_runner

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from screen_engine import (
    UnsafeSampleQueryError,
    apply_sample_query,
    bh_fdr,
    classify_selectivity,
    fit_dose_response,
    run_pipeline,
    ssmd,
)
from drug_repurposing_screen import main as cli_main

needs_pyarrow = pytest.mark.skipif(importlib.util.find_spec("pyarrow") is None, reason="needs pyarrow (uv sync --extra drug-repurposing-screen)")


@pytest.fixture(scope="module")
def demo_bundle(tmp_path_factory):
    out = tmp_path_factory.mktemp("demo_regen") / "demo"
    import subprocess
    subprocess.run(
        [sys.executable, str(SKILL_DIR / "generate_demo_bundle.py"), "--output-dir", str(out)],
        check=True, cwd=str(SKILL_DIR),
    )
    return out


class TestSafeSampleQuery:
    @pytest.fixture
    def sample_frame(self) -> pd.DataFrame:
        return pd.DataFrame({
            "sample_id": ["S1", "S2", "S3"],
            "context": ["context_A", "context_B", "context_A"],
            "lineage": ["organoid", "fibroblast", "organoid"],
            "passage": [3, 5, 7],
        })

    def test_valid_equality_query(self, sample_frame):
        mask = apply_sample_query(sample_frame, "context == 'context_A'")
        assert mask.tolist() == [True, False, True]

    def test_valid_compound_query(self, sample_frame):
        mask = apply_sample_query(
            sample_frame,
            "(context == 'context_A') & (lineage == 'organoid')",
        )
        assert mask.tolist() == [True, False, True]

    def test_valid_numeric_query(self, sample_frame):
        mask = apply_sample_query(sample_frame, "passage <= 5")
        assert mask.tolist() == [True, True, False]

    def test_empty_query_selects_all(self, sample_frame):
        mask = apply_sample_query(sample_frame, None)
        assert mask.all()

    def test_unknown_column_rejected(self, sample_frame):
        with pytest.raises(UnsafeSampleQueryError, match="Unknown sample_info column"):
            apply_sample_query(sample_frame, "missing_col == 'x'")

    def test_call_expression_rejected(self, sample_frame):
        with pytest.raises(UnsafeSampleQueryError, match="forbidden token"):
            apply_sample_query(sample_frame, "__import__('os')")

    def test_attribute_access_rejected(self, sample_frame):
        with pytest.raises(UnsafeSampleQueryError, match="Attribute access"):
            apply_sample_query(sample_frame, "context.real")

    def test_subscript_rejected(self, sample_frame):
        with pytest.raises(UnsafeSampleQueryError, match="Subscripting"):
            apply_sample_query(sample_frame, "context[0] == 'A'")


class TestStatistics:
    def test_ssmd_good_separation(self):
        neg = np.array([100, 102, 98, 101, 99], float)
        pos = np.array([20, 22, 18, 21, 19], float)
        assert ssmd(neg, pos) > 2.0

    def test_bh_fdr_monotone(self):
        p = np.array([0.001, 0.01, 0.05, 0.2])
        q = bh_fdr(p)
        assert q[0] <= q[1] <= q[2]

    def test_classify_selective(self):
        assert classify_selectivity(0.4, 0.3, 0.6) == "context_selective"

    def test_dose_response_good_curve(self):
        doses = np.array([0.01, 0.1, 1, 10], float)
        v = np.array([0.95, 0.8, 0.3, 0.1], float)
        fit = fit_dose_response(doses, v)
        assert fit["converged"]
        assert fit["quality"] in {"good", "low_r2", "right_censored"}


class TestDemoPipeline:
    def test_pipeline_finds_hits(self, demo_bundle):
        result = run_pipeline(
            demo_bundle,
            demo_bundle / "schema.yaml",
            demo_bundle / "objective.yaml",
        )
        assert result.summary["n_hits"] >= 1
        assert result.summary["n_samples"] == 10

    def test_selective_compounds_ranked(self, demo_bundle):
        result = run_pipeline(
            demo_bundle,
            demo_bundle / "schema.yaml",
            demo_bundle / "objective.yaml",
        )
        top_ids = set(result.priority.head(5)["compound_id"])
        assert "BRD-0003" in top_ids or "BRD-0007" in top_ids or "BRD-0015" in top_ids

    @needs_pyarrow
    def test_cli_demo_end_to_end(self, tmp_path):
        out = tmp_path / "demo_out"
        rc = cli_main(["--demo", "--output", str(out)])
        assert rc == 0
        assert (out / "report.md").exists()
        assert (out / "tables/priority_table.csv").exists()
        assert (out / "result.json").exists()
        summary = json.loads((out / "result.json").read_text())["summary"]
        assert summary["n_hits"] >= 1


class TestCLIBehaviour:
    @needs_pyarrow
    def test_cli_bundle_mode_round_trip(self, demo_bundle, tmp_path):
        out = tmp_path / "bundle_out"
        rc = cli_main([
            "--bundle", str(demo_bundle),
            "--schema", str(demo_bundle / "schema.yaml"),
            "--objective", str(demo_bundle / "objective.yaml"),
            "--output", str(out),
        ])
        assert rc == 0
        assert (out / "tables/priority_table.csv").exists()
        assert (out / "cache/priority.parquet").exists()
        assert (out / "reproducibility/schema.yaml").exists()
        assert (out / "reproducibility/objective.yaml").exists()
        assert (out / "reproducibility/commands.sh").exists()

    def test_cli_missing_required_flags_errors(self, tmp_path):
        out = tmp_path / "fail_out"
        with pytest.raises(SystemExit):
            cli_main(["--output", str(out)])

    def test_cli_missing_schema_returns_nonzero(self, demo_bundle, tmp_path):
        out = tmp_path / "missing_schema"
        rc = cli_main([
            "--bundle", str(demo_bundle),
            "--schema", str(demo_bundle / "does_not_exist.yaml"),
            "--objective", str(demo_bundle / "objective.yaml"),
            "--output", str(out),
        ])
        assert rc != 0

    def test_cli_missing_objective_returns_nonzero(self, demo_bundle, tmp_path):
        out = tmp_path / "missing_obj"
        rc = cli_main([
            "--bundle", str(demo_bundle),
            "--schema", str(demo_bundle / "schema.yaml"),
            "--objective", str(demo_bundle / "no_such.yaml"),
            "--output", str(out),
        ])
        assert rc != 0

    @needs_pyarrow
    def test_report_contains_disclaimer(self, tmp_path):
        out = tmp_path / "disc"
        rc = cli_main(["--demo", "--output", str(out)])
        assert rc == 0
        report = (out / "report.md").read_text()
        assert "ClawBio is a research and educational tool" in report
        assert "not a medical device" in report

    @needs_pyarrow
    def test_resume_short_circuits_when_cache_present(self, tmp_path):
        out = tmp_path / "resume"
        rc1 = cli_main(["--demo", "--output", str(out)])
        assert rc1 == 0
        priority_mtime_before = (out / "cache" / "priority.parquet").stat().st_mtime
        rc2 = cli_main(["--demo", "--output", str(out), "--resume"])
        assert rc2 == 0
        priority_mtime_after = (out / "cache" / "priority.parquet").stat().st_mtime
        assert priority_mtime_after == priority_mtime_before


class TestEngineFallbacks:
    def test_biomarker_sweep_handles_missing_features_dir(self, demo_bundle, tmp_path):
        from screen_engine import biomarker_sweep
        from screen_engine import load_yaml as _load_yaml
        schema = _load_yaml(demo_bundle / "schema.yaml")
        objective = _load_yaml(demo_bundle / "objective.yaml")
        empty_sel = pd.DataFrame({"compound_id": ["BRD-0003"], "selectivity_class": ["context_selective"]})
        out = biomarker_sweep(empty_sel, pd.DataFrame({"sample_id": []}), None, objective, schema)
        assert list(out.columns) == ["compound_id", "feature", "feature_type", "rho", "pvalue", "qvalue"]
        assert len(out) == 0

    def test_classify_selectivity_buckets(self):
        assert classify_selectivity(0.01, 0.95, 0.2) == "inactive"
        assert classify_selectivity(0.85, 0.5, 0.3) == "broadly_active"
        assert classify_selectivity(0.3, 0.4, 0.8) == "context_selective"
        assert classify_selectivity(0.4, 0.4, 0.2) == "other"


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/drug-repurposing-screen/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    # 12 significant digits: scipy 1.18 moved bimodality_coef in its 17th digit
    # (a last-bit rounding change in its moment code), which is not a real change.
    outputs = demo_runner.demo_outputs("drug-repurposing-screen", tmp_path, float_digits=12)
    assert outputs == snapshot({'report.html': """\
<html><body><pre># Drug Repurposing Screen Report

**Objective:** Approved compounds selective in IBD organoid context
**Generated:** <TIMESTAMP> UTC

## Summary

- Samples screened: 10
- Compounds tested: 22
- Primary hits: 3
- Context-selective compounds: 0
- Top candidate: `BRD-0003`

## Top prioritised candidates

| rank | compound_id | compound_name | selectivity_class | priority | feature | feature_type | clinical_phase |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | BRD-0003 | Drug_0003 | other | 0.861599124287 | MT1A | expression | Launched |
| 2 | BRD-0015 | Drug_0015 | other | 0.831540557481 | MT1A | expression | Launched |
| 3 | BRD-0007 | Drug_0007 | other | 0.799129688821 | cg_00015 | methylation | Phase 2 |
| 4 | BRD-0012 | Drug_0012 | nan | 0.325 | nan | nan | Launched |
| 5 | BRD-0018 | Drug_0018 | nan | 0.325 | nan | nan | Launched |
| 6 | BRD-0009 | Drug_0009 | nan | 0.325 | nan | nan | Launched |
| 7 | BRD-0006 | Drug_0006 | nan | 0.325 | nan | nan | Launched |
| 8 | BRD-0001 | Drug_0001 | nan | 0.265 | nan | nan | Phase 2 |
| 9 | BRD-0004 | Drug_0004 | nan | 0.265 | nan | nan | Phase 2 |
| 10 | BRD-0002 | Drug_0002 | nan | 0.265 | nan | nan | Phase 2 |

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.</pre></body></html>\
""", 'report.md': """\
# Drug Repurposing Screen Report

**Objective:** Approved compounds selective in IBD organoid context
**Generated:** <TIMESTAMP> UTC

## Summary

- Samples screened: 10
- Compounds tested: 22
- Primary hits: 3
- Context-selective compounds: 0
- Top candidate: `BRD-0003`

## Top prioritised candidates

| rank | compound_id | compound_name | selectivity_class | priority | feature | feature_type | clinical_phase |
| --- | --- | --- | --- | --- | --- | --- | --- |
| 1 | BRD-0003 | Drug_0003 | other | 0.861599124287 | MT1A | expression | Launched |
| 2 | BRD-0015 | Drug_0015 | other | 0.831540557481 | MT1A | expression | Launched |
| 3 | BRD-0007 | Drug_0007 | other | 0.799129688821 | cg_00015 | methylation | Phase 2 |
| 4 | BRD-0012 | Drug_0012 | nan | 0.325 | nan | nan | Launched |
| 5 | BRD-0018 | Drug_0018 | nan | 0.325 | nan | nan | Launched |
| 6 | BRD-0009 | Drug_0009 | nan | 0.325 | nan | nan | Launched |
| 7 | BRD-0006 | Drug_0006 | nan | 0.325 | nan | nan | Launched |
| 8 | BRD-0001 | Drug_0001 | nan | 0.265 | nan | nan | Phase 2 |
| 9 | BRD-0004 | Drug_0004 | nan | 0.265 | nan | nan | Phase 2 |
| 10 | BRD-0002 | Drug_0002 | nan | 0.265 | nan | nan | Phase 2 |

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.\
""", 'result.json': """\
{
  "summary": {
    "objective_name": "Approved compounds selective in IBD organoid context",
    "n_samples": 10,
    "n_compounds": 22,
    "n_hits": 3,
    "n_context_selective": 0,
    "top_compound": "BRD-0003"
  },
  "priority_table": [
    {
      "compound_id": "BRD-0003",
      "median_viability": 0.597968578575,
      "n_samples_hit": 10,
      "is_hit": true,
      "selectivity_class": "other",
      "kill_rate_all": 0.5,
      "median_viability_all": 0.612311388908,
      "bimodality_coef": 0.421163535019,
      "target_kill_rate": 1.0,
      "off_target_kill_rate": 0.0,
      "context_selectivity_score": 1.0,
      "compound_name": "Drug_0003",
      "clinical_phase": "Launched",
      "moa": "kinase inhibitor",
      "feature": "MT1A",
      "feature_type": "expression",
      "rho": -0.90303030303,
      "qvalue": 0.0536035028511,
      "priority": 0.861599124287,
      "rank": 1
    },
    {
      "compound_id": "BRD-0015",
      "median_viability": 0.627582464004,
      "n_samples_hit": 10,
      "is_hit": true,
      "selectivity_class": "other",
      "kill_rate_all": 0.5,
      "median_viability_all": 0.635056282656,
      "bimodality_coef": 0.416911917892,
      "target_kill_rate": 1.0,
      "off_target_kill_rate": 0.0,
      "context_selectivity_score": 1.0,
      "compound_name": "Drug_0015",
      "clinical_phase": "Launched",
      "moa": "kinase inhibitor",
      "feature": "MT1A",
      "feature_type": "expression",
      "rho": -0.757575757576,
      "qvalue": 0.173837770075,
      "priority": 0.831540557481,
      "rank": 2
    },
    {
      "compound_id": "BRD-0007",
      "median_viability": 0.605601528096,
      "n_samples_hit": 10,
      "is_hit": true,
      "selectivity_class": "other",
      "kill_rate_all": 0.5,
      "median_viability_all": 0.610963015305,
      "bimodality_coef": 0.416791028671,
      "target_kill_rate": 1.0,
      "off_target_kill_rate": 0.0,
      "context_selectivity_score": 1.0,
      "compound_name": "Drug_0007",
      "clinical_phase": "Phase 2",
      "moa": "kinase inhibitor",
      "feature": "cg_00015",
      "feature_type": "methylation",
      "rho": 0.878787878788,
      "qvalue": 0.0634812447151,
      "priority": 0.799129688821,
      "rank": 3
    },
    {
      "compound_id": "BRD-0012",
      "median_viability": 0.989409640822,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0012",
      "clinical_phase": "Launched",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.325,
      "rank": 4
    },
    {
      "compound_id": "BRD-0018",
      "median_viability": 0.964979267554,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0018",
      "clinical_phase": "Launched",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.325,
      "rank": 5
    },
    {
      "compound_id": "BRD-0009",
      "median_viability": 0.94689128605,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0009",
      "clinical_phase": "Launched",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.325,
      "rank": 6
    },
    {
      "compound_id": "BRD-0006",
      "median_viability": 0.960719667256,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0006",
      "clinical_phase": "Launched",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.325,
      "rank": 7
    },
    {
      "compound_id": "BRD-0001",
      "median_viability": 0.980469371504,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0001",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 8
    },
    {
      "compound_id": "BRD-0004",
      "median_viability": 0.992297841254,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0004",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 9
    },
    {
      "compound_id": "BRD-0002",
      "median_viability": 0.970782116929,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0002",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 10
    },
    {
      "compound_id": "BRD-0008",
      "median_viability": 0.960709756296,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0008",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 11
    },
    {
      "compound_id": "BRD-0005",
      "median_viability": 0.989642127546,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0005",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 12
    },
    {
      "compound_id": "BRD-0011",
      "median_viability": 0.966581211326,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0011",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 13
    },
    {
      "compound_id": "BRD-0010",
      "median_viability": 0.941977305512,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0010",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 14
    },
    {
      "compound_id": "BRD-0014",
      "median_viability": 0.978400024325,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0014",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 15
    },
    {
      "compound_id": "BRD-0013",
      "median_viability": 1.00013728274,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0013",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 16
    },
    {
      "compound_id": "BRD-0016",
      "median_viability": 0.952907475381,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0016",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 17
    },
    {
      "compound_id": "BRD-0017",
      "median_viability": 0.935989344054,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0017",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 18
    },
    {
      "compound_id": "BRD-0019",
      "median_viability": 0.949237435242,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0019",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 19
    },
    {
      "compound_id": "BRD-0020",
      "median_viability": 0.950019487809,
      "n_samples_hit": 0,
      "is_hit": false,
      "selectivity_class": NaN,
      "kill_rate_all": NaN,
      "median_viability_all": NaN,
      "bimodality_coef": NaN,
      "target_kill_rate": NaN,
      "off_target_kill_rate": NaN,
      "context_selectivity_score": NaN,
      "compound_name": "Drug_0020",
      "clinical_phase": "Phase 2",
      "moa": "unknown",
      "feature": NaN,
      "feature_type": NaN,
      "rho": NaN,
      "qvalue": NaN,
      "priority": 0.265,
      "rank": 20
    }
  ],
  "objective": "Approved compounds selective in IBD organoid context"
}\
""", 'tables/biomarker_univariate_all_matrices.csv': """\
compound_id,feature,feature_type,rho,pvalue,qvalue
BRD-0003,MT1A,expression,-0.90303030303,0.000343612197763,0.0536035028511
BRD-0007,cg_00015,methylation,0.878787878788,0.000813862111732,0.0634812447151
BRD-0007,GENE_22,expression,-0.854545454545,0.00163680331599,0.0851137724313
BRD-0007,cg_00020,methylation,-0.830303030303,0.00294022702328,0.114668853908
BRD-0007,GENE_17,expression,0.806060606061,0.00486206110196,0.126413588651
BRD-0007,cg_context_A,methylation,-0.806060606061,0.00486206110196,0.126413588651
BRD-0003,cg_00015,methylation,0.793939393939,0.0060999233137,0.135941148134
BRD-0015,MT1A,expression,-0.757575757576,0.0111434467997,0.173837770075
BRD-0015,cg_00004,methylation,-0.757575757576,0.0111434467997,0.173837770075
BRD-0015,cg_context_A,methylation,-0.757575757576,0.0111434467997,0.173837770075
BRD-0015,cg_00020,methylation,-0.745454545455,0.0133301463154,0.189045711383
BRD-0007,GENE_1,expression,0.721212121212,0.0185731550895,0.241451016163
BRD-0015,GENE_18,expression,-0.709090909091,0.0216659233675,0.25999108041
BRD-0015,cg_00007,methylation,0.684848484848,0.0288827975067,0.281607275691
BRD-0003,cg_context_A,methylation,-0.684848484848,0.0288827975067,0.281607275691
BRD-0003,GENE_21,expression,-0.684848484848,0.0288827975067,0.281607275691
BRD-0007,cg_00003,methylation,-0.660606060606,0.0375883775714,0.325765938952
BRD-0007,MT1A,expression,-0.660606060606,0.0375883775714,0.325765938952
BRD-0007,GENE_26,expression,0.624242424242,0.0537177672172,0.441051141362
BRD-0007,cg_00011,methylation,0.612121212121,0.0599721424773,0.467782711323
BRD-0003,GENE_26,expression,0.6,0.066688,0.472878545455
BRD-0015,GENE_9,expression,0.6,0.066688,0.472878545455
BRD-0003,GENE_22,expression,-0.587878787879,0.0738777068887,0.501083577158
BRD-0003,cg_00020,methylation,-0.563636363636,0.0897240283171,0.559877936699
BRD-0015,cg_00015,methylation,0.563636363636,0.0897240283171,0.559877936699
BRD-0007,cg_00017,methylation,-0.551515151515,0.0984011766696,0.590407060018
BRD-0015,GENE_10,expression,0.515151515152,0.127552869703,0.710651702632
BRD-0015,GENE_4,expression,0.515151515152,0.127552869703,0.710651702632
BRD-0015,cg_00017,methylation,-0.50303030303,0.138333698394,0.719335231651
BRD-0003,GENE_10,expression,0.50303030303,0.138333698394,0.719335231651
BRD-0015,GENE_21,expression,-0.490909090909,0.149655673421,0.719930764878
BRD-0015,GENE_25,expression,-0.478787878788,0.161522928017,0.719930764878
BRD-0003,GENE_17,expression,0.478787878788,0.161522928017,0.719930764878
BRD-0015,GENE_1,expression,0.478787878788,0.161522928017,0.719930764878
BRD-0007,GENE_8,expression,-0.478787878788,0.161522928017,0.719930764878
BRD-0007,GENE_11,expression,-0.466666666667,0.173938597165,0.753733921049
BRD-0015,GENE_12,expression,0.442424242424,0.200422686712,0.762583881148
BRD-0015,cg_00010,methylation,-0.454545454545,0.186904810349,0.762583881148
BRD-0003,cg_00012,methylation,-0.442424242424,0.200422686712,0.762583881148
BRD-0003,cg_00011,methylation,0.442424242424,0.200422686712,0.762583881148
BRD-0003,cg_00004,methylation,-0.442424242424,0.200422686712,0.762583881148
BRD-0015,GENE_29,expression,0.418181818182,0.229112840983,0.79425784874
BRD-0007,GENE_21,expression,-0.418181818182,0.229112840983,0.79425784874
BRD-0003,GENE_2,expression,0.418181818182,0.229112840983,0.79425784874
BRD-0015,GENE_26,expression,0.430303030303,0.214492332533,0.79425784874
BRD-0003,GENE_30,expression,0.406060606061,0.244282294087,0.810809316543
BRD-0007,GENE_2,expression,0.406060606061,0.244282294087,0.810809316543
BRD-0007,GENE_7,expression,0.393939393939,0.259997766835,0.811193032525
BRD-0007,GENE_28,expression,-0.393939393939,0.259997766835,0.811193032525
BRD-0015,GENE_17,expression,0.393939393939,0.259997766835,0.811193032525
BRD-0007,GENE_29,expression,-0.260606060606,0.467089054386,0.847277819585
BRD-0007,GENE_13,expression,0.284848484848,0.425038154892,0.847277819585
BRD-0007,GENE_23,expression,0.345454545455,0.328226511471,0.847277819585
BRD-0007,GENE_30,expression,0.309090909091,0.384841230261,0.847277819585
BRD-0007,GENE_10,expression,0.345454545455,0.328226511471,0.847277819585
BRD-0007,GENE_16,expression,0.309090909091,0.384841230261,0.847277819585
BRD-0003,cg_00018,methylation,-0.309090909091,0.384841230261,0.847277819585
BRD-0003,cg_00017,methylation,-0.369696969697,0.293050075293,0.847277819585
BRD-0007,cg_00013,methylation,-0.272727272727,0.445838341543,0.847277819585
BRD-0007,cg_00009,methylation,0.272727272727,0.445838341543,0.847277819585
BRD-0003,cg_00010,methylation,-0.272727272727,0.445838341543,0.847277819585
BRD-0007,GENE_12,expression,-0.260606060606,0.467089054386,0.847277819585
BRD-0003,cg_00007,methylation,0.345454545455,0.328226511471,0.847277819585
BRD-0003,cg_00005,methylation,-0.345454545455,0.328226511471,0.847277819585
BRD-0007,cg_00001,methylation,0.260606060606,0.467089054386,0.847277819585
BRD-0015,GENE_23,expression,0.284848484848,0.425038154892,0.847277819585
BRD-0003,cg_00013,methylation,-0.369696969697,0.293050075293,0.847277819585
BRD-0003,cg_00003,methylation,-0.272727272727,0.445838341543,0.847277819585
BRD-0015,GENE_22,expression,-0.321212121212,0.365468310439,0.847277819585
BRD-0007,cg_00007,methylation,0.345454545455,0.328226511471,0.847277819585
BRD-0015,GENE_27,expression,-0.272727272727,0.445838341543,0.847277819585
BRD-0003,GENE_20,expression,-0.333333333333,0.346593507087,0.847277819585
BRD-0003,GENE_25,expression,-0.29696969697,0.40470167127,0.847277819585
BRD-0003,GENE_18,expression,-0.321212121212,0.365468310439,0.847277819585
BRD-0003,GENE_3,expression,0.333333333333,0.346593507087,0.847277819585
BRD-0003,GENE_1,expression,0.369696969697,0.293050075293,0.847277819585
BRD-0003,GENE_11,expression,-0.29696969697,0.40470167127,0.847277819585
BRD-0003,GENE_16,expression,0.309090909091,0.384841230261,0.847277819585
BRD-0003,GENE_7,expression,0.345454545455,0.328226511471,0.847277819585
BRD-0003,GENE_24,expression,-0.260606060606,0.467089054386,0.847277819585
BRD-0015,cg_00012,methylation,-0.29696969697,0.40470167127,0.847277819585
BRD-0015,cg_00002,methylation,-0.333333333333,0.346593507087,0.847277819585
BRD-0015,cg_00003,methylation,-0.284848484848,0.425038154892,0.847277819585
BRD-0015,cg_00009,methylation,-0.284848484848,0.425038154892,0.847277819585
BRD-0003,GENE_19,expression,-0.333333333333,0.346593507087,0.847277819585
BRD-0003,GENE_23,expression,0.309090909091,0.384841230261,0.847277819585
BRD-0003,GENE_9,expression,0.248484848485,0.488776304519,0.876426477069
BRD-0003,cg_00016,methylation,-0.236363636364,0.510885317515,0.895484376768
BRD-0007,cg_00006,methylation,0.236363636364,0.510885317515,0.895484376768
BRD-0007,GENE_4,expression,0.224242424242,0.533400561273,0.924560972872
BRD-0015,cg_00014,methylation,0.212121212121,0.556305775103,0.933158074366
BRD-0015,GENE_16,expression,0.212121212121,0.556305775103,0.933158074366
BRD-0015,GENE_8,expression,-0.212121212121,0.556305775103,0.933158074366
BRD-0007,cg_00004,methylation,-0.2,0.579584,0.951737936842
BRD-0007,cg_00018,methylation,-0.2,0.579584,0.951737936842
BRD-0003,cg_00019,methylation,0.187878787879,0.603217610029,0.970123166645
BRD-0003,GENE_4,expression,0.187878787879,0.603217610029,0.970123166645
BRD-0007,GENE_25,expression,0.175757575758,0.627188344776,0.978413817851
BRD-0007,cg_00019,methylation,-0.175757575758,0.627188344776,0.978413817851
BRD-0015,GENE_3,expression,0.175757575758,0.627188344776,0.978413817851
BRD-0003,GENE_14,expression,0.139393939394,0.70093188491,0.985093459874
BRD-0007,GENE_9,expression,0.139393939394,0.70093188491,0.985093459874
BRD-0003,cg_00002,methylation,-0.151515151515,0.676065175998,0.985093459874
BRD-0015,GENE_11,expression,-0.139393939394,0.70093188491,0.985093459874
BRD-0003,cg_00001,methylation,0.139393939394,0.70093188491,0.985093459874
BRD-0015,GENE_19,expression,-0.163636363636,0.651477342796,0.985093459874
BRD-0015,GENE_5,expression,0.139393939394,0.70093188491,0.985093459874
BRD-0015,GENE_15,expression,0.139393939394,0.70093188491,0.985093459874
BRD-0007,GENE_15,expression,0.139393939394,0.70093188491,0.985093459874
BRD-0015,cg_00001,methylation,0.163636363636,0.651477342796,0.985093459874
BRD-0007,GENE_19,expression,-0.151515151515,0.676065175998,0.985093459874
BRD-0003,GENE_6,expression,-0.0787878787879,0.828717394697,0.986742911195
BRD-0015,GENE_2,expression,-0.127272727273,0.726057014763,0.986742911195
BRD-0007,GENE_27,expression,-0.0787878787879,0.828717394697,0.986742911195
BRD-0007,GENE_20,expression,0.0181818181818,0.960240418129,0.986742911195
BRD-0007,GENE_24,expression,0.0181818181818,0.960240418129,0.986742911195
BRD-0007,GENE_18,expression,-0.030303030303,0.933772958094,0.986742911195
BRD-0007,GENE_5,expression,0.030303030303,0.933772958094,0.986742911195
BRD-0007,GENE_3,expression,0.127272727273,0.726057014763,0.986742911195
BRD-0007,GENE_14,expression,0.0787878787879,0.828717394697,0.986742911195
BRD-0007,GENE_6,expression,-0.0181818181818,0.960240418129,0.986742911195
BRD-0003,GENE_13,expression,0.0787878787879,0.828717394697,0.986742911195
BRD-0003,GENE_12,expression,0.0424242424242,0.907363817813,0.986742911195
BRD-0003,GENE_8,expression,-0.0909090909091,0.802771731207,0.986742911195
BRD-0003,GENE_15,expression,-0.10303030303,0.776998463444,0.986742911195
BRD-0003,GENE_5,expression,0.0666666666667,0.854813088249,0.986742911195
BRD-0003,GENE_28,expression,-0.10303030303,0.776998463444,0.986742911195
BRD-0003,GENE_29,expression,-0.127272727273,0.726057014763,0.986742911195
BRD-0007,cg_00014,methylation,-0.030303030303,0.933772958094,0.986742911195
BRD-0007,cg_00012,methylation,0.00606060606061,0.986742911195,0.986742911195
BRD-0007,cg_00008,methylation,-0.0424242424242,0.907363817813,0.986742911195
BRD-0007,cg_00010,methylation,0.0787878787879,0.828717394697,0.986742911195
BRD-0007,cg_00002,methylation,-0.115151515152,0.751419652326,0.986742911195
BRD-0007,cg_00005,methylation,-0.030303030303,0.933772958094,0.986742911195
BRD-0003,cg_00006,methylation,0.00606060606061,0.986742911195,0.986742911195
BRD-0003,cg_00014,methylation,-0.030303030303,0.933772958094,0.986742911195
BRD-0003,cg_00008,methylation,0.0545454545455,0.881036181162,0.986742911195
BRD-0003,cg_00009,methylation,-0.0181818181818,0.960240418129,0.986742911195
BRD-0015,GENE_6,expression,-0.115151515152,0.751419652326,0.986742911195
BRD-0015,GENE_7,expression,-0.0909090909091,0.802771731207,0.986742911195
BRD-0015,GENE_24,expression,0.00606060606061,0.986742911195,0.986742911195
BRD-0015,GENE_20,expression,-0.115151515152,0.751419652326,0.986742911195
BRD-0015,GENE_30,expression,0.0666666666667,0.854813088249,0.986742911195
BRD-0015,GENE_28,expression,-0.0666666666667,0.854813088249,0.986742911195
BRD-0015,GENE_13,expression,-0.030303030303,0.933772958094,0.986742911195
BRD-0015,GENE_14,expression,0.0424242424242,0.907363817813,0.986742911195
BRD-0003,GENE_27,expression,-0.00606060606061,0.986742911195,0.986742911195
BRD-0007,cg_00016,methylation,0.0545454545455,0.881036181162,0.986742911195
BRD-0015,cg_00008,methylation,0.00606060606061,0.986742911195,0.986742911195
BRD-0015,cg_00006,methylation,-0.030303030303,0.933772958094,0.986742911195
BRD-0015,cg_00005,methylation,0.0909090909091,0.802771731207,0.986742911195
BRD-0015,cg_00011,methylation,-0.0666666666667,0.854813088249,0.986742911195
BRD-0015,cg_00016,methylation,0.10303030303,0.776998463444,0.986742911195
BRD-0015,cg_00013,methylation,-0.0666666666667,0.854813088249,0.986742911195
BRD-0015,cg_00019,methylation,-0.0909090909091,0.802771731207,0.986742911195
BRD-0015,cg_00018,methylation,0.0181818181818,0.960240418129,0.986742911195
""", 'tables/priority_table.csv': """\
compound_id,median_viability,n_samples_hit,is_hit,selectivity_class,kill_rate_all,median_viability_all,bimodality_coef,target_kill_rate,off_target_kill_rate,context_selectivity_score,compound_name,clinical_phase,moa,feature,feature_type,rho,qvalue,priority,rank
BRD-0003,0.597968578575,10,True,other,0.5,0.612311388908,0.421163535019,1.0,0.0,1.0,Drug_0003,Launched,kinase inhibitor,MT1A,expression,-0.90303030303,0.0536035028511,0.861599124287,1
BRD-0015,0.627582464004,10,True,other,0.5,0.635056282656,0.416911917892,1.0,0.0,1.0,Drug_0015,Launched,kinase inhibitor,MT1A,expression,-0.757575757576,0.173837770075,0.831540557481,2
BRD-0007,0.605601528096,10,True,other,0.5,0.610963015305,0.416791028671,1.0,0.0,1.0,Drug_0007,Phase 2,kinase inhibitor,cg_00015,methylation,0.878787878788,0.0634812447151,0.799129688821,3
BRD-0012,0.989409640822,0,False,,,,,,,,Drug_0012,Launched,unknown,,,,,0.325,4
BRD-0018,0.964979267554,0,False,,,,,,,,Drug_0018,Launched,unknown,,,,,0.325,5
BRD-0009,0.94689128605,0,False,,,,,,,,Drug_0009,Launched,unknown,,,,,0.325,6
BRD-0006,0.960719667256,0,False,,,,,,,,Drug_0006,Launched,unknown,,,,,0.325,7
BRD-0001,0.980469371504,0,False,,,,,,,,Drug_0001,Phase 2,unknown,,,,,0.265,8
BRD-0004,0.992297841254,0,False,,,,,,,,Drug_0004,Phase 2,unknown,,,,,0.265,9
BRD-0002,0.970782116929,0,False,,,,,,,,Drug_0002,Phase 2,unknown,,,,,0.265,10
BRD-0008,0.960709756296,0,False,,,,,,,,Drug_0008,Phase 2,unknown,,,,,0.265,11
BRD-0005,0.989642127546,0,False,,,,,,,,Drug_0005,Phase 2,unknown,,,,,0.265,12
BRD-0011,0.966581211326,0,False,,,,,,,,Drug_0011,Phase 2,unknown,,,,,0.265,13
BRD-0010,0.941977305512,0,False,,,,,,,,Drug_0010,Phase 2,unknown,,,,,0.265,14
BRD-0014,0.978400024325,0,False,,,,,,,,Drug_0014,Phase 2,unknown,,,,,0.265,15
BRD-0013,1.00013728274,0,False,,,,,,,,Drug_0013,Phase 2,unknown,,,,,0.265,16
BRD-0016,0.952907475381,0,False,,,,,,,,Drug_0016,Phase 2,unknown,,,,,0.265,17
BRD-0017,0.935989344054,0,False,,,,,,,,Drug_0017,Phase 2,unknown,,,,,0.265,18
BRD-0019,0.949237435242,0,False,,,,,,,,Drug_0019,Phase 2,unknown,,,,,0.265,19
BRD-0020,0.950019487809,0,False,,,,,,,,Drug_0020,Phase 2,unknown,,,,,0.265,20
""", 'tables/selectivity.csv': """\
compound_id,selectivity_class,kill_rate_all,median_viability_all,bimodality_coef,target_kill_rate,off_target_kill_rate,context_selectivity_score
BRD-0003,other,0.5,0.612311388908,0.421163535019,1.0,0.0,1.0
BRD-0007,other,0.5,0.610963015305,0.416791028671,1.0,0.0,1.0
BRD-0015,other,0.5,0.635056282656,0.416911917892,1.0,0.0,1.0
"""})
