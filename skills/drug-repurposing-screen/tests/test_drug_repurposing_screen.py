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
    assert demo_runner.demo_outputs("drug-repurposing-screen", tmp_path) == snapshot({'report.html': """\
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
| 1 | BRD-0003 | Drug_0003 | other | 0.8615991242872321 | MT1A | expression | Launched |
| 2 | BRD-0015 | Drug_0015 | other | 0.8315405574811925 | MT1A | expression | Launched |
| 3 | BRD-0007 | Drug_0007 | other | 0.7991296888212219 | cg_00015 | methylation | Phase 2 |
| 4 | BRD-0012 | Drug_0012 | nan | 0.325 | nan | nan | Launched |
| 5 | BRD-0018 | Drug_0018 | nan | 0.325 | nan | nan | Launched |
| 6 | BRD-0009 | Drug_0009 | nan | 0.325 | nan | nan | Launched |
| 7 | BRD-0006 | Drug_0006 | nan | 0.325 | nan | nan | Launched |
| 8 | BRD-0001 | Drug_0001 | nan | 0.26499999999999996 | nan | nan | Phase 2 |
| 9 | BRD-0004 | Drug_0004 | nan | 0.26499999999999996 | nan | nan | Phase 2 |
| 10 | BRD-0002 | Drug_0002 | nan | 0.26499999999999996 | nan | nan | Phase 2 |

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
| 1 | BRD-0003 | Drug_0003 | other | 0.8615991242872321 | MT1A | expression | Launched |
| 2 | BRD-0015 | Drug_0015 | other | 0.8315405574811925 | MT1A | expression | Launched |
| 3 | BRD-0007 | Drug_0007 | other | 0.7991296888212219 | cg_00015 | methylation | Phase 2 |
| 4 | BRD-0012 | Drug_0012 | nan | 0.325 | nan | nan | Launched |
| 5 | BRD-0018 | Drug_0018 | nan | 0.325 | nan | nan | Launched |
| 6 | BRD-0009 | Drug_0009 | nan | 0.325 | nan | nan | Launched |
| 7 | BRD-0006 | Drug_0006 | nan | 0.325 | nan | nan | Launched |
| 8 | BRD-0001 | Drug_0001 | nan | 0.26499999999999996 | nan | nan | Phase 2 |
| 9 | BRD-0004 | Drug_0004 | nan | 0.26499999999999996 | nan | nan | Phase 2 |
| 10 | BRD-0002 | Drug_0002 | nan | 0.26499999999999996 | nan | nan | Phase 2 |

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
      "median_viability": 0.5979685785753599,
      "n_samples_hit": 10,
      "is_hit": true,
      "selectivity_class": "other",
      "kill_rate_all": 0.5,
      "median_viability_all": 0.6123113889084708,
      "bimodality_coef": 0.42116353501861975,
      "target_kill_rate": 1.0,
      "off_target_kill_rate": 0.0,
      "context_selectivity_score": 1.0,
      "compound_name": "Drug_0003",
      "clinical_phase": "Launched",
      "moa": "kinase inhibitor",
      "feature": "MT1A",
      "feature_type": "expression",
      "rho": -0.9030303030303028,
      "qvalue": 0.05360350285107206,
      "priority": 0.8615991242872321,
      "rank": 1
    },
    {
      "compound_id": "BRD-0015",
      "median_viability": 0.6275824640041725,
      "n_samples_hit": 10,
      "is_hit": true,
      "selectivity_class": "other",
      "kill_rate_all": 0.5,
      "median_viability_all": 0.635056282655958,
      "bimodality_coef": 0.41691191789245413,
      "target_kill_rate": 1.0,
      "off_target_kill_rate": 0.0,
      "context_selectivity_score": 1.0,
      "compound_name": "Drug_0015",
      "clinical_phase": "Launched",
      "moa": "kinase inhibitor",
      "feature": "MT1A",
      "feature_type": "expression",
      "rho": -0.7575757575757575,
      "qvalue": 0.17383777007522966,
      "priority": 0.8315405574811925,
      "rank": 2
    },
    {
      "compound_id": "BRD-0007",
      "median_viability": 0.6056015280961387,
      "n_samples_hit": 10,
      "is_hit": true,
      "selectivity_class": "other",
      "kill_rate_all": 0.5,
      "median_viability_all": 0.6109630153051955,
      "bimodality_coef": 0.41679102867139545,
      "target_kill_rate": 1.0,
      "off_target_kill_rate": 0.0,
      "context_selectivity_score": 1.0,
      "compound_name": "Drug_0007",
      "clinical_phase": "Phase 2",
      "moa": "kinase inhibitor",
      "feature": "cg_00015",
      "feature_type": "methylation",
      "rho": 0.8787878787878788,
      "qvalue": 0.06348124471511239,
      "priority": 0.7991296888212219,
      "rank": 3
    },
    {
      "compound_id": "BRD-0012",
      "median_viability": 0.9894096408221311,
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
      "median_viability": 0.9649792675542913,
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
      "median_viability": 0.9468912860499223,
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
      "median_viability": 0.9607196672556051,
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
      "median_viability": 0.9804693715037797,
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
      "priority": 0.26499999999999996,
      "rank": 8
    },
    {
      "compound_id": "BRD-0004",
      "median_viability": 0.9922978412544275,
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
      "priority": 0.26499999999999996,
      "rank": 9
    },
    {
      "compound_id": "BRD-0002",
      "median_viability": 0.9707821169293829,
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
      "priority": 0.26499999999999996,
      "rank": 10
    },
    {
      "compound_id": "BRD-0008",
      "median_viability": 0.9607097562961999,
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
      "priority": 0.26499999999999996,
      "rank": 11
    },
    {
      "compound_id": "BRD-0005",
      "median_viability": 0.9896421275456887,
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
      "priority": 0.26499999999999996,
      "rank": 12
    },
    {
      "compound_id": "BRD-0011",
      "median_viability": 0.9665812113256658,
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
      "priority": 0.26499999999999996,
      "rank": 13
    },
    {
      "compound_id": "BRD-0010",
      "median_viability": 0.9419773055122245,
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
      "priority": 0.26499999999999996,
      "rank": 14
    },
    {
      "compound_id": "BRD-0014",
      "median_viability": 0.9784000243252334,
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
      "priority": 0.26499999999999996,
      "rank": 15
    },
    {
      "compound_id": "BRD-0013",
      "median_viability": 1.000137282743537,
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
      "priority": 0.26499999999999996,
      "rank": 16
    },
    {
      "compound_id": "BRD-0016",
      "median_viability": 0.9529074753805287,
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
      "priority": 0.26499999999999996,
      "rank": 17
    },
    {
      "compound_id": "BRD-0017",
      "median_viability": 0.935989344054436,
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
      "priority": 0.26499999999999996,
      "rank": 18
    },
    {
      "compound_id": "BRD-0019",
      "median_viability": 0.949237435241945,
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
      "priority": 0.26499999999999996,
      "rank": 19
    },
    {
      "compound_id": "BRD-0020",
      "median_viability": 0.9500194878088266,
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
      "priority": 0.26499999999999996,
      "rank": 20
    }
  ],
  "objective": "Approved compounds selective in IBD organoid context"
}\
""", 'tables/biomarker_univariate_all_matrices.csv': """\
compound_id,feature,feature_type,rho,pvalue,qvalue
BRD-0003,MT1A,expression,-0.9030303030303028,0.0003436121977632824,0.05360350285107206
BRD-0007,cg_00015,methylation,0.8787878787878788,0.0008138621117322101,0.06348124471511239
BRD-0007,GENE_22,expression,-0.8545454545454544,0.0016368033159867156,0.0851137724313092
BRD-0007,cg_00020,methylation,-0.8303030303030302,0.002940227023279509,0.11466885390790085
BRD-0007,GENE_17,expression,0.8060606060606059,0.004862061101964607,0.1264135886510798
BRD-0007,cg_context_A,methylation,-0.8060606060606059,0.004862061101964607,0.1264135886510798
BRD-0003,cg_00015,methylation,0.7939393939393938,0.006099923313696911,0.13594114813381686
BRD-0015,MT1A,expression,-0.7575757575757575,0.011143446799694208,0.17383777007522966
BRD-0015,cg_00004,methylation,-0.7575757575757575,0.011143446799694208,0.17383777007522966
BRD-0015,cg_context_A,methylation,-0.7575757575757575,0.011143446799694208,0.17383777007522966
BRD-0015,cg_00020,methylation,-0.7454545454545454,0.013330146315440054,0.1890457113826044
BRD-0007,GENE_1,expression,0.7212121212121211,0.01857315508946023,0.24145101616298295
BRD-0015,GENE_18,expression,-0.709090909090909,0.02166592336753244,0.25999108041038926
BRD-0015,cg_00007,methylation,0.6848484848484848,0.028882797506732765,0.2816072756906445
BRD-0003,cg_context_A,methylation,-0.6848484848484848,0.028882797506732765,0.2816072756906445
BRD-0003,GENE_21,expression,-0.6848484848484848,0.028882797506732765,0.2816072756906445
BRD-0007,cg_00003,methylation,-0.6606060606060605,0.03758837757140934,0.3257659389522143
BRD-0007,MT1A,expression,-0.6606060606060605,0.03758837757140934,0.3257659389522143
BRD-0007,GENE_26,expression,0.6242424242424242,0.053717767217167436,0.4410511413620063
BRD-0007,cg_00011,methylation,0.6121212121212121,0.05997214247725527,0.4677827113225911
BRD-0003,GENE_26,expression,0.6,0.06668799999999998,0.4728785454545453
BRD-0015,GENE_9,expression,0.6,0.06668799999999998,0.4728785454545453
BRD-0003,GENE_22,expression,-0.5878787878787878,0.07387770688865808,0.5010835771578548
BRD-0003,cg_00020,methylation,-0.5636363636363636,0.0897240283170912,0.559877936698649
BRD-0015,cg_00015,methylation,0.5636363636363636,0.0897240283170912,0.559877936698649
BRD-0007,cg_00017,methylation,-0.5515151515151515,0.09840117666963498,0.59040706001781
BRD-0015,GENE_10,expression,0.5151515151515151,0.1275528697032343,0.7106517026323054
BRD-0015,GENE_4,expression,0.5151515151515151,0.1275528697032343,0.7106517026323054
BRD-0015,cg_00017,methylation,-0.5030303030303029,0.13833369839449194,0.7193352316513582
BRD-0003,GENE_10,expression,0.5030303030303029,0.13833369839449194,0.7193352316513582
BRD-0015,GENE_21,expression,-0.49090909090909085,0.1496556734208368,0.7199307648778035
BRD-0015,GENE_25,expression,-0.4787878787878787,0.16152292801745594,0.7199307648778035
BRD-0003,GENE_17,expression,0.4787878787878787,0.16152292801745594,0.7199307648778035
BRD-0015,GENE_1,expression,0.4787878787878787,0.16152292801745594,0.7199307648778035
BRD-0007,GENE_8,expression,-0.4787878787878787,0.16152292801745594,0.7199307648778035
BRD-0007,GENE_11,expression,-0.4666666666666666,0.17393859716506654,0.7537339210486217
BRD-0015,GENE_12,expression,0.4424242424242424,0.20042268671194213,0.7625838811478773
BRD-0015,cg_00010,methylation,-0.4545454545454545,0.18690481034933196,0.7625838811478773
BRD-0003,cg_00012,methylation,-0.4424242424242424,0.20042268671194213,0.7625838811478773
BRD-0003,cg_00011,methylation,0.4424242424242424,0.20042268671194213,0.7625838811478773
BRD-0003,cg_00004,methylation,-0.4424242424242424,0.20042268671194213,0.7625838811478773
BRD-0015,GENE_29,expression,0.41818181818181815,0.22911284098281892,0.794257848740439
BRD-0007,GENE_21,expression,-0.41818181818181815,0.22911284098281892,0.794257848740439
BRD-0003,GENE_2,expression,0.41818181818181815,0.22911284098281892,0.794257848740439
BRD-0015,GENE_26,expression,0.43030303030303024,0.21449233253280775,0.794257848740439
BRD-0003,GENE_30,expression,0.406060606060606,0.24428229408662647,0.8108093165428454
BRD-0007,GENE_2,expression,0.406060606060606,0.24428229408662647,0.8108093165428454
BRD-0007,GENE_7,expression,0.3939393939393939,0.25999776683488757,0.8111930325248493
BRD-0007,GENE_28,expression,-0.3939393939393939,0.25999776683488757,0.8111930325248493
BRD-0015,GENE_17,expression,0.3939393939393939,0.25999776683488757,0.8111930325248493
BRD-0007,GENE_29,expression,-0.2606060606060606,0.46708905438634024,0.8472778195845242
BRD-0007,GENE_13,expression,0.28484848484848485,0.4250381548921454,0.8472778195845242
BRD-0007,GENE_23,expression,0.34545454545454546,0.3282265114713672,0.8472778195845242
BRD-0007,GENE_30,expression,0.309090909090909,0.38484123026128353,0.8472778195845242
BRD-0007,GENE_10,expression,0.34545454545454546,0.3282265114713672,0.8472778195845242
BRD-0007,GENE_16,expression,0.309090909090909,0.38484123026128353,0.8472778195845242
BRD-0003,cg_00018,methylation,-0.309090909090909,0.38484123026128353,0.8472778195845242
BRD-0003,cg_00017,methylation,-0.3696969696969697,0.29305007529282234,0.8472778195845242
BRD-0007,cg_00013,methylation,-0.2727272727272727,0.4458383415427516,0.8472778195845242
BRD-0007,cg_00009,methylation,0.2727272727272727,0.4458383415427516,0.8472778195845242
BRD-0003,cg_00010,methylation,-0.2727272727272727,0.4458383415427516,0.8472778195845242
BRD-0007,GENE_12,expression,-0.2606060606060606,0.46708905438634024,0.8472778195845242
BRD-0003,cg_00007,methylation,0.34545454545454546,0.3282265114713672,0.8472778195845242
BRD-0003,cg_00005,methylation,-0.34545454545454546,0.3282265114713672,0.8472778195845242
BRD-0007,cg_00001,methylation,0.2606060606060606,0.46708905438634024,0.8472778195845242
BRD-0015,GENE_23,expression,0.28484848484848485,0.4250381548921454,0.8472778195845242
BRD-0003,cg_00013,methylation,-0.3696969696969697,0.29305007529282234,0.8472778195845242
BRD-0003,cg_00003,methylation,-0.2727272727272727,0.4458383415427516,0.8472778195845242
BRD-0015,GENE_22,expression,-0.32121212121212117,0.36546831043867045,0.8472778195845242
BRD-0007,cg_00007,methylation,0.34545454545454546,0.3282265114713672,0.8472778195845242
BRD-0015,GENE_27,expression,-0.2727272727272727,0.4458383415427516,0.8472778195845242
BRD-0003,GENE_20,expression,-0.33333333333333326,0.3465935070873343,0.8472778195845242
BRD-0003,GENE_25,expression,-0.29696969696969694,0.40470167127015705,0.8472778195845242
BRD-0003,GENE_18,expression,-0.32121212121212117,0.36546831043867045,0.8472778195845242
BRD-0003,GENE_3,expression,0.33333333333333326,0.3465935070873343,0.8472778195845242
BRD-0003,GENE_1,expression,0.3696969696969697,0.29305007529282234,0.8472778195845242
BRD-0003,GENE_11,expression,-0.29696969696969694,0.40470167127015705,0.8472778195845242
BRD-0003,GENE_16,expression,0.309090909090909,0.38484123026128353,0.8472778195845242
BRD-0003,GENE_7,expression,0.34545454545454546,0.3282265114713672,0.8472778195845242
BRD-0003,GENE_24,expression,-0.2606060606060606,0.46708905438634024,0.8472778195845242
BRD-0015,cg_00012,methylation,-0.29696969696969694,0.40470167127015705,0.8472778195845242
BRD-0015,cg_00002,methylation,-0.33333333333333326,0.3465935070873343,0.8472778195845242
BRD-0015,cg_00003,methylation,-0.28484848484848485,0.4250381548921454,0.8472778195845242
BRD-0015,cg_00009,methylation,-0.28484848484848485,0.4250381548921454,0.8472778195845242
BRD-0003,GENE_19,expression,-0.33333333333333326,0.3465935070873343,0.8472778195845242
BRD-0003,GENE_23,expression,0.309090909090909,0.38484123026128353,0.8472778195845242
BRD-0003,GENE_9,expression,0.24848484848484845,0.4887763045192432,0.8764264770689878
BRD-0003,cg_00016,methylation,-0.23636363636363633,0.5108853175152003,0.8954843767682162
BRD-0007,cg_00006,methylation,0.23636363636363633,0.5108853175152003,0.8954843767682162
BRD-0007,GENE_4,expression,0.22424242424242422,0.5334005612725947,0.9245609728724975
BRD-0015,cg_00014,methylation,0.2121212121212121,0.5563057751029299,0.9331580743662049
BRD-0015,GENE_16,expression,0.2121212121212121,0.5563057751029299,0.9331580743662049
BRD-0015,GENE_8,expression,-0.2121212121212121,0.5563057751029299,0.9331580743662049
BRD-0007,cg_00004,methylation,-0.19999999999999998,0.5795840000000002,0.9517379368421055
BRD-0007,cg_00018,methylation,-0.19999999999999998,0.5795840000000002,0.9517379368421055
BRD-0003,cg_00019,methylation,0.18787878787878787,0.6032176100292089,0.9701231666449133
BRD-0003,GENE_4,expression,0.18787878787878787,0.6032176100292089,0.9701231666449133
BRD-0007,GENE_25,expression,0.17575757575757575,0.6271883447764846,0.9784138178513159
BRD-0007,cg_00019,methylation,-0.17575757575757575,0.6271883447764846,0.9784138178513159
BRD-0015,GENE_3,expression,0.17575757575757575,0.6271883447764846,0.9784138178513159
BRD-0003,GENE_14,expression,0.13939393939393938,0.7009318849100584,0.9850934598735955
BRD-0007,GENE_9,expression,0.13939393939393938,0.7009318849100584,0.9850934598735955
BRD-0003,cg_00002,methylation,-0.1515151515151515,0.6760651759978539,0.9850934598735955
BRD-0015,GENE_11,expression,-0.13939393939393938,0.7009318849100584,0.9850934598735955
BRD-0003,cg_00001,methylation,0.13939393939393938,0.7009318849100584,0.9850934598735955
BRD-0015,GENE_19,expression,-0.16363636363636364,0.6514773427962428,0.9850934598735955
BRD-0015,GENE_5,expression,0.13939393939393938,0.7009318849100584,0.9850934598735955
BRD-0015,GENE_15,expression,0.13939393939393938,0.7009318849100584,0.9850934598735955
BRD-0007,GENE_15,expression,0.13939393939393938,0.7009318849100584,0.9850934598735955
BRD-0015,cg_00001,methylation,0.16363636363636364,0.6514773427962428,0.9850934598735955
BRD-0007,GENE_19,expression,-0.1515151515151515,0.6760651759978539,0.9850934598735955
BRD-0003,GENE_6,expression,-0.07878787878787878,0.8287173946974606,0.9867429111949892
BRD-0015,GENE_2,expression,-0.12727272727272726,0.7260570147627894,0.9867429111949892
BRD-0007,GENE_27,expression,-0.07878787878787878,0.8287173946974606,0.9867429111949892
BRD-0007,GENE_20,expression,0.01818181818181818,0.9602404181286242,0.9867429111949892
BRD-0007,GENE_24,expression,0.01818181818181818,0.9602404181286242,0.9867429111949892
BRD-0007,GENE_18,expression,-0.0303030303030303,0.9337729580941466,0.9867429111949892
BRD-0007,GENE_5,expression,0.0303030303030303,0.9337729580941466,0.9867429111949892
BRD-0007,GENE_3,expression,0.12727272727272726,0.7260570147627894,0.9867429111949892
BRD-0007,GENE_14,expression,0.07878787878787878,0.8287173946974606,0.9867429111949892
BRD-0007,GENE_6,expression,-0.01818181818181818,0.9602404181286242,0.9867429111949892
BRD-0003,GENE_13,expression,0.07878787878787878,0.8287173946974606,0.9867429111949892
BRD-0003,GENE_12,expression,0.04242424242424241,0.9073638178128159,0.9867429111949892
BRD-0003,GENE_8,expression,-0.09090909090909088,0.8027717312071619,0.9867429111949892
BRD-0003,GENE_15,expression,-0.10303030303030303,0.7769984634438898,0.9867429111949892
BRD-0003,GENE_5,expression,0.06666666666666665,0.8548130882487427,0.9867429111949892
BRD-0003,GENE_28,expression,-0.10303030303030303,0.7769984634438898,0.9867429111949892
BRD-0003,GENE_29,expression,-0.12727272727272726,0.7260570147627894,0.9867429111949892
BRD-0007,cg_00014,methylation,-0.0303030303030303,0.9337729580941466,0.9867429111949892
BRD-0007,cg_00012,methylation,0.006060606060606061,0.9867429111949892,0.9867429111949892
BRD-0007,cg_00008,methylation,-0.04242424242424241,0.9073638178128159,0.9867429111949892
BRD-0007,cg_00010,methylation,0.07878787878787878,0.8287173946974606,0.9867429111949892
BRD-0007,cg_00002,methylation,-0.11515151515151514,0.7514196523258485,0.9867429111949892
BRD-0007,cg_00005,methylation,-0.0303030303030303,0.9337729580941466,0.9867429111949892
BRD-0003,cg_00006,methylation,0.006060606060606061,0.9867429111949892,0.9867429111949892
BRD-0003,cg_00014,methylation,-0.0303030303030303,0.9337729580941466,0.9867429111949892
BRD-0003,cg_00008,methylation,0.05454545454545454,0.8810361811618526,0.9867429111949892
BRD-0003,cg_00009,methylation,-0.01818181818181818,0.9602404181286242,0.9867429111949892
BRD-0015,GENE_6,expression,-0.11515151515151514,0.7514196523258485,0.9867429111949892
BRD-0015,GENE_7,expression,-0.09090909090909088,0.8027717312071619,0.9867429111949892
BRD-0015,GENE_24,expression,0.006060606060606061,0.9867429111949892,0.9867429111949892
BRD-0015,GENE_20,expression,-0.11515151515151514,0.7514196523258485,0.9867429111949892
BRD-0015,GENE_30,expression,0.06666666666666665,0.8548130882487427,0.9867429111949892
BRD-0015,GENE_28,expression,-0.06666666666666665,0.8548130882487427,0.9867429111949892
BRD-0015,GENE_13,expression,-0.0303030303030303,0.9337729580941466,0.9867429111949892
BRD-0015,GENE_14,expression,0.04242424242424241,0.9073638178128159,0.9867429111949892
BRD-0003,GENE_27,expression,-0.006060606060606061,0.9867429111949892,0.9867429111949892
BRD-0007,cg_00016,methylation,0.05454545454545454,0.8810361811618526,0.9867429111949892
BRD-0015,cg_00008,methylation,0.006060606060606061,0.9867429111949892,0.9867429111949892
BRD-0015,cg_00006,methylation,-0.0303030303030303,0.9337729580941466,0.9867429111949892
BRD-0015,cg_00005,methylation,0.09090909090909088,0.8027717312071619,0.9867429111949892
BRD-0015,cg_00011,methylation,-0.06666666666666665,0.8548130882487427,0.9867429111949892
BRD-0015,cg_00016,methylation,0.10303030303030303,0.7769984634438898,0.9867429111949892
BRD-0015,cg_00013,methylation,-0.06666666666666665,0.8548130882487427,0.9867429111949892
BRD-0015,cg_00019,methylation,-0.09090909090909088,0.8027717312071619,0.9867429111949892
BRD-0015,cg_00018,methylation,0.01818181818181818,0.9602404181286242,0.9867429111949892
""", 'tables/priority_table.csv': """\
compound_id,median_viability,n_samples_hit,is_hit,selectivity_class,kill_rate_all,median_viability_all,bimodality_coef,target_kill_rate,off_target_kill_rate,context_selectivity_score,compound_name,clinical_phase,moa,feature,feature_type,rho,qvalue,priority,rank
BRD-0003,0.5979685785753599,10,True,other,0.5,0.6123113889084708,0.42116353501861975,1.0,0.0,1.0,Drug_0003,Launched,kinase inhibitor,MT1A,expression,-0.9030303030303028,0.05360350285107206,0.8615991242872321,1
BRD-0015,0.6275824640041725,10,True,other,0.5,0.635056282655958,0.41691191789245413,1.0,0.0,1.0,Drug_0015,Launched,kinase inhibitor,MT1A,expression,-0.7575757575757575,0.17383777007522966,0.8315405574811925,2
BRD-0007,0.6056015280961387,10,True,other,0.5,0.6109630153051955,0.41679102867139545,1.0,0.0,1.0,Drug_0007,Phase 2,kinase inhibitor,cg_00015,methylation,0.8787878787878788,0.06348124471511239,0.7991296888212219,3
BRD-0012,0.9894096408221311,0,False,,,,,,,,Drug_0012,Launched,unknown,,,,,0.325,4
BRD-0018,0.9649792675542913,0,False,,,,,,,,Drug_0018,Launched,unknown,,,,,0.325,5
BRD-0009,0.9468912860499223,0,False,,,,,,,,Drug_0009,Launched,unknown,,,,,0.325,6
BRD-0006,0.9607196672556051,0,False,,,,,,,,Drug_0006,Launched,unknown,,,,,0.325,7
BRD-0001,0.9804693715037797,0,False,,,,,,,,Drug_0001,Phase 2,unknown,,,,,0.26499999999999996,8
BRD-0004,0.9922978412544275,0,False,,,,,,,,Drug_0004,Phase 2,unknown,,,,,0.26499999999999996,9
BRD-0002,0.9707821169293829,0,False,,,,,,,,Drug_0002,Phase 2,unknown,,,,,0.26499999999999996,10
BRD-0008,0.9607097562961999,0,False,,,,,,,,Drug_0008,Phase 2,unknown,,,,,0.26499999999999996,11
BRD-0005,0.9896421275456887,0,False,,,,,,,,Drug_0005,Phase 2,unknown,,,,,0.26499999999999996,12
BRD-0011,0.9665812113256658,0,False,,,,,,,,Drug_0011,Phase 2,unknown,,,,,0.26499999999999996,13
BRD-0010,0.9419773055122245,0,False,,,,,,,,Drug_0010,Phase 2,unknown,,,,,0.26499999999999996,14
BRD-0014,0.9784000243252334,0,False,,,,,,,,Drug_0014,Phase 2,unknown,,,,,0.26499999999999996,15
BRD-0013,1.000137282743537,0,False,,,,,,,,Drug_0013,Phase 2,unknown,,,,,0.26499999999999996,16
BRD-0016,0.9529074753805287,0,False,,,,,,,,Drug_0016,Phase 2,unknown,,,,,0.26499999999999996,17
BRD-0017,0.935989344054436,0,False,,,,,,,,Drug_0017,Phase 2,unknown,,,,,0.26499999999999996,18
BRD-0019,0.949237435241945,0,False,,,,,,,,Drug_0019,Phase 2,unknown,,,,,0.26499999999999996,19
BRD-0020,0.9500194878088266,0,False,,,,,,,,Drug_0020,Phase 2,unknown,,,,,0.26499999999999996,20
""", 'tables/selectivity.csv': """\
compound_id,selectivity_class,kill_rate_all,median_viability_all,bimodality_coef,target_kill_rate,off_target_kill_rate,context_selectivity_score
BRD-0003,other,0.5,0.6123113889084708,0.42116353501861975,1.0,0.0,1.0
BRD-0007,other,0.5,0.6109630153051955,0.41679102867139545,1.0,0.0,1.0
BRD-0015,other,0.5,0.635056282655958,0.41691191789245413,1.0,0.0,1.0
"""})
