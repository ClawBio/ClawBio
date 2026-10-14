"""Tests for drug-repurposing-screen skill."""

from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

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
    out = tmp_path_factory.mktemp("demo_regen")
    import subprocess
    subprocess.run(
        [sys.executable, str(SKILL_DIR / "generate_demo_bundle.py")],
        check=True, cwd=str(SKILL_DIR),
    )
    return SKILL_DIR / "demo"


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
