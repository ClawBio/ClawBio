"""Tests for the scientific-audit baseline regression gate.

Provenance: #342. The bench was invoked as `clawbio-bench --smoke || true`, so
the job reported SUCCESS on every PR while the bench itself returned
`pass: false` with seven blocking harnesses. A gate that cannot fail is worse
than no gate, because the green tick stops anyone looking.

The gate deliberately does NOT require `overall.pass`. Seven harnesses are
below 100% today; demanding a clean bench would block every PR on pre-existing
debt and the `|| true` would go straight back on. It fails on *regression*
against a committed baseline instead.
"""
import importlib.util
import json
from pathlib import Path

import pytest


def _load():
    module_path = Path(__file__).resolve().parents[1] / "scripts" / "check_bench_baseline.py"
    spec = importlib.util.spec_from_file_location("check_bench_baseline", module_path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


MOD = _load()

BASELINE = {"harnesses": {"alpha": 90.0, "beta": 0.0}}


def _report(rates, harness_errors=None):
    harness_errors = harness_errors or {}
    return {
        "harnesses": {
            name: {"pass_rate": rate, "harness_errors": harness_errors.get(name, 0)}
            for name, rate in rates.items()
        }
    }


def test_unchanged_rates_pass():
    failures = MOD.compare(_report({"alpha": 90.0, "beta": 0.0}), BASELINE)
    assert failures == []


def test_improvement_passes():
    failures = MOD.compare(_report({"alpha": 97.5, "beta": 40.0}), BASELINE)
    assert failures == []


def test_regression_fails_and_names_the_harness():
    failures = MOD.compare(_report({"alpha": 80.0, "beta": 0.0}), BASELINE)
    assert len(failures) == 1
    assert "alpha" in failures[0]
    assert "90.0" in failures[0] and "80.0" in failures[0]


def test_regression_below_a_zero_baseline_is_impossible_but_beta_cannot_silently_vanish():
    """A harness present in the baseline and absent from the report is a
    regression, not a pass. Deleting a failing harness must not go green."""
    failures = MOD.compare(_report({"alpha": 90.0}), BASELINE)
    assert len(failures) == 1
    assert "beta" in failures[0]
    assert "missing" in failures[0].lower()


def test_new_harness_without_baseline_does_not_fail():
    """A newly added harness has no recorded rate; record it, do not block."""
    failures = MOD.compare(_report({"alpha": 90.0, "beta": 0.0, "gamma": 55.0}), BASELINE)
    assert failures == []


def test_float_noise_within_tolerance_passes():
    failures = MOD.compare(_report({"alpha": 89.995, "beta": 0.0}), BASELINE)
    assert failures == []


def test_drop_beyond_tolerance_fails():
    failures = MOD.compare(_report({"alpha": 89.9, "beta": 0.0}), BASELINE)
    assert len(failures) == 1


def test_new_harness_error_is_a_regression():
    """A harness that starts erroring reports pass_rate 0.0 or vanishes; catch
    the explicit error count too, since an errored harness proves nothing."""
    report = _report({"alpha": 90.0, "beta": 0.0}, harness_errors={"alpha": 2})
    failures = MOD.compare(report, BASELINE, baseline_errors={"alpha": 0, "beta": 0})
    assert len(failures) == 1
    assert "error" in failures[0].lower()


def test_committed_baseline_matches_schema():
    """The checked-in baseline must be loadable and non-empty, or the gate is
    silently vacuous."""
    baseline = MOD.load_baseline(
        Path(__file__).resolve().parents[1] / "bench_baseline.json"
    )
    assert baseline["harnesses"]
    assert "nutrigx-advisor" in baseline["harnesses"]
    for name, rate in baseline["harnesses"].items():
        assert isinstance(rate, (int, float)), name


def test_main_exits_nonzero_on_regression(tmp_path):
    report = tmp_path / "aggregate_report.json"
    report.write_text(json.dumps(_report({"alpha": 10.0, "beta": 0.0})))
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps(BASELINE))
    assert MOD.main([str(report), "--baseline", str(baseline)]) == 1


def test_main_exits_zero_when_clean(tmp_path):
    report = tmp_path / "aggregate_report.json"
    report.write_text(json.dumps(_report({"alpha": 90.0, "beta": 0.0})))
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps(BASELINE))
    assert MOD.main([str(report), "--baseline", str(baseline)]) == 0


def test_missing_report_is_an_error_not_a_pass(tmp_path):
    """The bench failing to produce a report at all must not read as success.
    That is the #342 failure mode in a different costume."""
    baseline = tmp_path / "baseline.json"
    baseline.write_text(json.dumps(BASELINE))
    assert MOD.main([str(tmp_path / "nope.json"), "--baseline", str(baseline)]) == 2
