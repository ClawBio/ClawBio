"""Tests for mendelian-randomisation skill.

Validates MR estimators (IVW, Egger, weighted median/mode), sensitivity
tests, instrument diagnostics, and end-to-end demo mode.
"""

from __future__ import annotations

import json
import math
import random
import sys
from pathlib import Path

import numpy as np
from scipy import stats
import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from mendelian_randomisation import (
    Instrument,
    MIN_SENSITIVITY_INSTRUMENTS,
    MREstimate,
    NoInstrumentsError,
    cochran_q,
    egger_weighted_dispersion,
    compute_i_squared_gx,
    generate_report,
    ivw,
    leave_one_out,
    load_demo_instruments,
    mr_egger,
    run_pipeline,
    run_sensitivity,
    steiger_test,
    weighted_median,
    weighted_mode,
)


@pytest.fixture
def demo_instruments():
    instruments, _, _ = load_demo_instruments()
    return instruments


@pytest.fixture
def known_causal():
    """Create instruments with a known causal effect of 0.5."""
    rng = np.random.default_rng(42)
    instruments = []
    for i in range(20):
        bx = rng.uniform(0.03, 0.07) * rng.choice([1, -1])
        se_x = abs(bx) / rng.uniform(5, 10)
        by = bx * 0.5 + rng.normal(0, se_x * 0.2)
        se_y = se_x * rng.uniform(1.5, 2.0)
        instruments.append(Instrument(
            snp=f"rs{i}", effect_allele="A", other_allele="G",
            eaf=rng.uniform(0.2, 0.8),
            beta_exposure=float(bx), se_exposure=float(se_x),
            pval_exposure=1e-10,
            beta_outcome=float(by), se_outcome=float(se_y),
            pval_outcome=0.01,
            f_statistic=float((bx / se_x) ** 2),
        ))
    return instruments


# ---------------------------------------------------------------------------
# Unit tests — estimators
# ---------------------------------------------------------------------------
class TestIVW:
    def test_returns_mr_estimate(self, demo_instruments):
        est = ivw(demo_instruments)
        assert isinstance(est, MREstimate)
        assert est.method == "IVW"
        assert est.n_snps == 30

    def test_recovers_true_effect(self, demo_instruments):
        est = ivw(demo_instruments)
        assert 0.3 < est.estimate < 0.9

    def test_significant_pvalue(self, demo_instruments):
        est = ivw(demo_instruments)
        assert est.pvalue < 0.05

    def test_ci_contains_estimate(self, demo_instruments):
        est = ivw(demo_instruments)
        assert est.ci_lower < est.estimate < est.ci_upper

    def test_known_causal_effect(self, known_causal):
        est = ivw(known_causal)
        assert abs(est.estimate - 0.5) < 0.3


class TestMREgger:
    def test_returns_estimate_and_intercept(self, demo_instruments):
        est, intercept, se_int, p_int = mr_egger(demo_instruments)
        assert est.method == "MR-Egger"
        assert isinstance(intercept, float)
        assert isinstance(p_int, float)

    def test_intercept_near_zero_when_no_pleiotropy(self, demo_instruments):
        _, intercept, _, p_int = mr_egger(demo_instruments)
        assert abs(intercept) < 0.1
        assert p_int > 0.05

    def test_matches_twosamplemr_mr_egger_regression_on_the_demo(self, demo_instruments):
        """Cross-implementation parity. TwoSampleMR's mr_egger_regression (R/mr.R at
        commit c14776b8, the function loaded verbatim into R 4.6.0) on these 30
        instruments returns slope 0.602231, se 0.081583, intercept -0.000187,
        intercept se 0.003121, slope p 4.874e-08, intercept p 0.9526 (both on a t
        reference with n - 2 df). Before instruments were oriented, this code
        returned 0.598932 and 0.039122 for the slope; before the t reference, it
        returned slope p 1.56e-13 on a normal."""
        est, intercept, se_int, p_int = mr_egger(demo_instruments)
        assert est.estimate == pytest.approx(0.602231, abs=5e-7)
        assert est.se == pytest.approx(0.081583, abs=5e-7)
        assert intercept == pytest.approx(-0.000187, abs=5e-7)
        assert se_int == pytest.approx(0.003121, abs=5e-7)
        assert est.pvalue == pytest.approx(4.874e-8, rel=1e-3)
        assert p_int == pytest.approx(0.9526, abs=5e-5)

    def test_egger_p_value_is_a_t_on_n_minus_2_df(self, demo_instruments):
        """On three instruments Egger has one residual degree of freedom. A normal
        reference called the slope significant (p = 0.006); the t on 1 df gives
        0.22, which is what TwoSampleMR reports. Recomputed here from the returned
        estimate and SE so the test pins the reference, not a number."""
        from scipy import stats as _st
        est, intercept, se_int, p_int = mr_egger(demo_instruments[:3])
        assert est.pvalue == pytest.approx(2 * _st.t.sf(abs(est.estimate / est.se), df=1), rel=1e-9)
        assert p_int == pytest.approx(2 * _st.t.sf(abs(intercept / se_int), df=1), rel=1e-9)
        assert est.pvalue > 0.1  # the normal reference gave 0.006 here

    def test_intercept_does_not_depend_on_allele_coding(self, demo_instruments):
        """Re-coding an instrument to its other allele negates both effects. The
        slope never depended on that; the intercept did, before instruments were
        oriented to a positive exposure effect (as TwoSampleMR does). Every other
        instrument is re-coded, so 18 of 30 demo rows are negative before and half
        of those flip again here."""
        est0, int0, se0, p0 = mr_egger(demo_instruments)
        recoded = [
            Instrument(**{**vars(i), "beta_exposure": -i.beta_exposure, "beta_outcome": -i.beta_outcome})
            if k % 2 else i
            for k, i in enumerate(demo_instruments)
        ]
        est1, int1, se1, p1 = mr_egger(recoded)
        assert est1.estimate == pytest.approx(est0.estimate, rel=1e-12)
        assert est1.se == pytest.approx(est0.se, rel=1e-12)
        assert int1 == pytest.approx(int0, rel=1e-12)
        assert se1 == pytest.approx(se0, rel=1e-12)
        assert p1 == pytest.approx(p0, rel=1e-12)


class TestWeightedMedian:
    def test_returns_estimate(self, demo_instruments):
        est = weighted_median(demo_instruments)
        assert est.method == "Weighted Median"
        assert est.n_snps == 30

    def test_consistent_with_ivw(self, demo_instruments):
        ivw_est = ivw(demo_instruments)
        wm_est = weighted_median(demo_instruments)
        assert abs(ivw_est.estimate - wm_est.estimate) < 0.3


class TestWeightedMode:
    def test_returns_estimate(self, demo_instruments):
        est = weighted_mode(demo_instruments)
        assert est.method == "Weighted Mode"


# ---------------------------------------------------------------------------
# Unit tests — sensitivity
# ---------------------------------------------------------------------------
class TestSensitivity:
    def test_cochran_q_no_heterogeneity(self, demo_instruments):
        ivw_est = ivw(demo_instruments)
        q, p, df = cochran_q(demo_instruments, ivw_est)
        assert df == 29
        assert p > 0.05

    def test_steiger_correct_direction(self, demo_instruments):
        # UPDATED: the test now returns (correct, p_or_None, note). The p is None here
        # because the demo instruments carry no sample sizes, and without them there is
        # no basis for one -- see `test_steiger_verdict_is_invariant_to_the_units`.
        correct, p, note = steiger_test(demo_instruments)
        assert correct is True
        assert p is None
        assert "no sample sizes" in note

    def test_i_squared_gx(self, demo_instruments):
        i2 = compute_i_squared_gx(demo_instruments)
        assert 0 <= i2 <= 1

    def test_leave_one_out_length(self, demo_instruments):
        loo = leave_one_out(demo_instruments)
        assert len(loo) == 30

    def test_run_sensitivity_returns_all_fields(self, demo_instruments):
        ivw_est = ivw(demo_instruments)
        sens = run_sensitivity(demo_instruments, ivw_est)
        assert sens.mean_f_statistic > 0
        assert sens.cochran_q >= 0


# ---------------------------------------------------------------------------
# Unit tests — instrument properties
# ---------------------------------------------------------------------------
class TestInstrumentProperties:
    def test_palindromic_AT(self):
        i = Instrument("rs1", "A", "T", 0.5, 0.05, 0.01, 1e-10, 0.03, 0.02, 0.1, 25.0)
        assert i.is_palindromic is True
        assert i.palindromic_ambiguous is True

    def test_palindromic_not_ambiguous(self):
        i = Instrument("rs1", "A", "T", 0.2, 0.05, 0.01, 1e-10, 0.03, 0.02, 0.1, 25.0)
        assert i.is_palindromic is True
        assert i.palindromic_ambiguous is False

    def test_not_palindromic(self):
        i = Instrument("rs1", "A", "C", 0.5, 0.05, 0.01, 1e-10, 0.03, 0.02, 0.1, 25.0)
        assert i.is_palindromic is False

    def test_weak_instrument(self):
        i = Instrument("rs1", "A", "G", 0.3, 0.05, 0.02, 1e-10, 0.03, 0.02, 0.1, 6.0)
        assert i.weak_instrument is True

    def test_strong_instrument(self):
        i = Instrument("rs1", "A", "G", 0.3, 0.05, 0.005, 1e-10, 0.03, 0.02, 0.1, 100.0)
        assert i.weak_instrument is False


# ---------------------------------------------------------------------------
# Integration tests — demo pipeline
# ---------------------------------------------------------------------------
class TestDemoPipeline:
    def test_demo_end_to_end(self, tmp_path):
        instruments, exposure, outcome = load_demo_instruments()
        summary = run_pipeline(instruments, exposure, outcome, tmp_path, demo=True)

        assert summary["n_instruments"] == 30
        assert 0.3 < summary["ivw_estimate"] < 0.9
        assert summary["ivw_pvalue"] < 0.05
        assert summary["n_weak"] == 0

        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()
        assert (tmp_path / "tables" / "mr_results.tsv").exists()
        assert (tmp_path / "tables" / "sensitivity.tsv").exists()
        assert (tmp_path / "tables" / "harmonised_instruments.tsv").exists()
        assert (tmp_path / "figures" / "scatter.png").exists()
        assert (tmp_path / "figures" / "forest.png").exists()
        assert (tmp_path / "figures" / "funnel.png").exists()
        assert (tmp_path / "figures" / "leave_one_out.png").exists()

    def test_demo_report_content(self, tmp_path):
        instruments, exposure, outcome = load_demo_instruments()
        run_pipeline(instruments, exposure, outcome, tmp_path, demo=True)

        report = (tmp_path / "report.md").read_text()
        assert "Mendelian Randomisation" in report
        assert "BMI" in report
        assert "T2D" in report
        assert "ClawBio is a research" in report
        assert "IVW" in report

        result = json.loads((tmp_path / "result.json").read_text())
        assert result["exposure"] == "Body mass index (BMI)"
        assert len(result["estimates"]) == 4
        assert result["sensitivity"]["n_weak"] == 0

    def test_all_methods_consistent(self, tmp_path):
        instruments, exposure, outcome = load_demo_instruments()
        run_pipeline(instruments, exposure, outcome, tmp_path, demo=True)

        result = json.loads((tmp_path / "result.json").read_text())
        estimates = {e["method"]: e["estimate"] for e in result["estimates"]}
        ivw_est = estimates["IVW"]
        for method, est in estimates.items():
            assert abs(est - ivw_est) < 0.2, f"{method} ({est:.3f}) diverges from IVW ({ivw_est:.3f})"


# ---------------------------------------------------------------------------
# MR-Egger applicability: two independent conditions, and a trap in the obvious
# test suite for them
# ---------------------------------------------------------------------------
#
# MR-Egger fits a slope AND an intercept. Two things have to hold, and neither implies
# the other:
#
#   STATISTICAL  n >= 3, or there is no residual degree of freedom at all.
#   NUMERICAL    at least two distinct beta_exposure values, which is not an n < 3
#                problem -- it is reachable at any n from a merge that duplicates one
#                SNP across rows, tissues or proxies.
#
# THE TRAP, and why the cases below are shaped the way they are: the n < 3 return comes
# FIRST, so every n=1 and n=2 input short-circuits before it ever reaches the
# conditioning test. A suite made only of n=1 and n=2 cases therefore passes whether or
# not the centered denominator and the relative threshold work at all -- which is the
# "a gate never observed failing is not known to work" failure, committed by the very
# tests written to prevent it. The n>=3 near-collinear case is the only input that
# reaches the second guard, so it is the load-bearing one here.


def _egger_input(bxs, by_scale=0.5, se_y=0.02):
    """Instruments with the given exposure effects and a clean linear outcome."""
    return [
        Instrument(snp=f"rs{k}", effect_allele="A", other_allele="G", eaf=0.3,
                   beta_exposure=float(bx), se_exposure=abs(float(bx)) / 8 or 0.01,
                   pval_exposure=1e-10,
                   beta_outcome=float(bx) * by_scale, se_outcome=se_y,
                   pval_outcome=0.01, f_statistic=64.0)
        for k, bx in enumerate(bxs)
    ]


@pytest.mark.parametrize("n", [1, 2])
def test_egger_below_three_instruments_is_not_applicable(n):
    """The statistical condition. Previously n=1 was a three-way coin flip -- guard,
    crash, or a finite slope with a median SE of 2.1e+07 -- and n=2 exited 0 while
    writing an infinite standard error."""
    est, intercept, int_se, int_p = mr_egger(_egger_input([0.2, 0.5][:n]))

    assert est.applicable is False
    assert est.n_snps == n
    assert "at least 3" in est.reason
    assert math.isnan(est.estimate) and math.isnan(est.se)
    # the intercept (the pleiotropy test) is undefined on exactly the same inputs
    assert math.isnan(intercept) and math.isnan(int_se) and math.isnan(int_p)


def test_egger_with_near_identical_exposure_effects_is_not_applicable():
    """The NUMERICAL condition, and the only input in this file that reaches it.

    n=5, so the instrument-count return does not fire; the exposure effects agree to
    ~1e-11 relative, so the slope is unidentified. The old code returned a number here:
    the same shape gave slope=3960 with se=2.64e+06 and no warning of any kind.
    """
    base = 0.4
    est, intercept, _, _ = mr_egger(
        _egger_input([base + k * 1e-11 for k in range(5)]))

    assert est.applicable is False
    assert est.n_snps == 5
    assert "too close to identical" in est.reason
    assert math.isnan(est.estimate)
    assert math.isnan(intercept)


def test_the_centered_dispersion_is_non_negative_where_the_expanded_form_is_not():
    """The sign, which is the root cause rather than the threshold.

    Tested on the denominator DIRECTLY rather than through `mr_egger`, and the reason
    is worth stating: once the conditioning check is in place it refuses these inputs
    anyway, so the public path behaves identically with either form. Centering is a
    STRUCTURAL property, not a behavioural one, and a test routed through the public
    API would pass without it -- which is the failure this whole file is about.

    The construction matters. Equally-spaced exposure effects with equal weights do NOT
    reproduce the defect: the cancellation stays non-negative by symmetry. It needs
    effects that are near-identical but irregularly spaced, with the varying weights
    real standard errors produce. Over the 200 draws below the expanded form goes
    negative for roughly a third of them, matching the ~35% measured at n=5 across
    4,000 draws; the centered form was non-negative in all 64,000 draws of that sweep.
    """
    rng = random.Random(20260829)
    saw_expanded_go_negative = 0
    for _ in range(200):
        base = rng.uniform(0.05, 0.9)
        bx = np.array([base * (1 + rng.uniform(-1e-12, 1e-12)) for _ in range(5)])
        se = np.array([rng.uniform(0.005, 0.05) for _ in range(5)])
        w = 1.0 / se ** 2

        centered = egger_weighted_dispersion(w, bx)
        assert centered >= 0.0, f"centered dispersion went negative: {centered!r}"
        assert math.isfinite(centered)

        expanded = np.sum(w) * np.sum(w * bx ** 2) - np.sum(w * bx) ** 2
        if expanded < 0:
            saw_expanded_go_negative += 1

    # exactly zero when the exposure effects are identical, which is the algebraic
    # statement that the slope is unidentified
    assert egger_weighted_dispersion(np.full(4, 2500.0), np.full(4, 0.4)) == 0.0

    # and the failure the centering removes is observed on this very sweep, not
    # asserted from a story about it
    assert saw_expanded_go_negative > 20, (
        f"the expanded form went negative only {saw_expanded_go_negative} times in 200 "
        "draws, so this test is not observing the defect it exists for")


@pytest.mark.parametrize("n", [3, 5])
def test_a_well_conditioned_fit_still_produces_an_estimate(n):
    """The false-positive direction: neither guard may fire on ordinary data."""
    est, intercept, int_se, int_p = mr_egger(
        _egger_input([0.1 * (k + 1) for k in range(n)]))

    assert est.applicable is True
    assert est.reason == ""
    assert math.isfinite(est.estimate) and math.isfinite(est.se)
    assert est.se < 1.0                      # not the inflated-SE branch
    assert est.estimate == pytest.approx(0.5, abs=0.05)
    assert all(math.isfinite(v) for v in (intercept, int_se, int_p))


def test_a_two_instrument_run_writes_parseable_json_and_a_report_that_says_why(tmp_path):
    """End to end, because every symptom of this bug was in the ARTIFACTS.

    A two-instrument run used to exit 0 having written `"se": Infinity` into
    result.json -- which RFC 8259 has no token for, so a strict parser rejects the whole
    document -- and a report.md asserting the result was "robust", because the
    consistency check compared point estimates and never looked at the standard errors.
    """
    insts = _egger_input([0.2, 0.5])
    out = tmp_path / "run"
    out.mkdir()
    run_pipeline(insts, exposure="X", outcome="Y", output_dir=out, demo=True)

    # json.load rejects bare Infinity/NaN by default, so this parses only if nothing
    # non-finite was written.
    result = json.loads((out / "result.json").read_text())
    egger = next(e for e in result["estimates"] if e["method"] == "MR-Egger")
    assert egger["applicable"] is False
    assert "at least 3" in egger["reason"]
    assert "estimate" not in egger and "se" not in egger
    assert result["sensitivity"]["egger_intercept"] is None

    # Below three instruments none of the three sensitivity estimators applies (the
    # weighted median and mode are order statistics of the ratios; TwoSampleMR returns
    # NA for all three below 3 SNPs), so every one of them is declared, not numbered.
    for method in ("Weighted Median", "Weighted Mode"):
        row = next(e for e in result["estimates"] if e["method"] == method)
        assert row["applicable"] is False and "estimate" not in row, row

    report = (out / "report.md").read_text()
    assert "MR-Egger was not computed" in report
    # The POSITIVE assertion, not merely the absence of the old sentence. With nothing
    # to compare against, the report must say the IVW stands alone -- not certify it as
    # "consistent across methods" on the strength of estimators that never ran, which
    # is what a two-instrument run used to do.
    assert ("No sensitivity estimator applies to this instrument set, so the IVW "
            "estimate stands alone and is not corroborated.") in report
    assert "consistent estimates across" not in report
    assert "| Egger intercept | not computed |" in report

    table = (out / "tables" / "mr_results.tsv").read_text()
    egger_row = next(l for l in table.splitlines() if l.startswith("MR-Egger"))
    assert "not_applicable" in egger_row
    assert "inf" not in egger_row.lower()


def test_the_json_writer_refuses_any_non_finite_value(tmp_path):
    """`allow_nan=False` is a BACKSTOP, and the fix removed every input that trips it.

    That is the point of it -- the pipeline no longer produces a non-finite -- but it
    also means no end-to-end test can observe it working. So drive the writer directly
    with a value it must refuse. Without the flag, Python writes a bare `Infinity`,
    which RFC 8259 has no token for, so a strict parser rejects the WHOLE document
    rather than one field. That is what a two-instrument run used to ship, exit code 0.
    """
    from mendelian_randomisation import SensitivityResults, _write_result_json

    broken = MREstimate("IVW", float("inf"), 1.0, 0.0, 2.0, 0.5, 4)
    with pytest.raises(ValueError):
        _write_result_json([broken], SensitivityResults(), 0.1, 0.5,
                           "X", "Y", tmp_path, "2026-01-01T00:00:00Z", True)


def test_the_scatter_plot_does_not_advertise_an_estimator_that_did_not_run(tmp_path):
    """A NaN slope draws no visible line and still claims a legend entry.

    matplotlib silently omits a line whose y values are NaN, so the plot LOOKS right --
    but the legend still gains an "MR-Egger (nan)" row, and a reader takes a legend
    entry as a statement that the method was fitted. The absence has to be an absence
    everywhere, not only in the tables.

    Recorded at the point the label is REQUESTED rather than read off the finished
    figure: `scatter_plot` closes its figure, so inspecting the current one afterwards
    reads an empty canvas and the test passes for the wrong reason.
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.axes
    from mendelian_randomisation import scatter_plot

    requested = []
    real_plot = matplotlib.axes.Axes.plot

    def spy(self, *a, **kw):
        if "label" in kw:
            requested.append(kw["label"])
        return real_plot(self, *a, **kw)

    insts = _egger_input([0.2, 0.5])
    estimates = [
        ivw(insts),
        MREstimate.not_applicable("MR-Egger", 2, "requires at least 3 instruments"),
        weighted_median(insts),
    ]
    matplotlib.axes.Axes.plot = spy
    try:
        scatter_plot(insts, estimates, tmp_path / "scatter.png")
    finally:
        matplotlib.axes.Axes.plot = real_plot

    assert (tmp_path / "scatter.png").exists()
    assert requested, "no labelled line was drawn at all; the spy is not wired"
    assert not any("MR-Egger" in lab for lab in requested), requested
    assert not any("nan" in lab.lower() for lab in requested), requested
    assert any("IVW" in lab for lab in requested), requested


def test_a_zero_outcome_standard_error_is_refused_rather_than_passed_through():
    """The conditioning check is written `not (rel > t)`, and this is why.

    An instrument with `se_outcome = 0` gives it infinite weight, which makes both the
    dispersion and its scale infinite and their ratio NaN. Every comparison against NaN
    is False, so the natural-reading `rel < threshold` would be False and the estimator
    would proceed on a quantity that is not a number. The negated form turns that same
    False into a refusal.

    Reachable from a malformed input file rather than only in principle, which is why it
    is tested instead of merely commented.
    """
    insts = [
        Instrument(snp=f"rs{k}", effect_allele="A", other_allele="G", eaf=0.3,
                   beta_exposure=0.1 * (k + 1), se_exposure=0.01, pval_exposure=1e-10,
                   beta_outcome=0.05 * (k + 1), se_outcome=se, pval_outcome=0.01,
                   f_statistic=64.0)
        for k, se in enumerate([0.02, 0.0, 0.02])
    ]

    # The infinite weight is the point of the test; numpy's divide-by-zero and
    # invalid-value warnings on the way to the NaN ratio are expected, not a defect.
    with np.errstate(divide="ignore", invalid="ignore"):
        est, intercept, _, _ = mr_egger(insts)

    assert est.applicable is False
    assert "not identified" in est.reason
    assert math.isnan(est.estimate) and math.isnan(intercept)


# ---------------------------------------------------------------------------
# Steiger directionality: the verdict must not depend on the units of the traits
# ---------------------------------------------------------------------------

def _steiger_input(unit=1.0, n_exp=None, n_out=None, seed=7):
    rng = np.random.default_rng(seed)
    out = []
    for k in range(10):
        bx = float(rng.uniform(0.05, 0.4))
        out.append(Instrument(
            snp=f"rs{k}", effect_allele="A", other_allele="G", eaf=0.3,
            beta_exposure=bx, se_exposure=bx / 8, pval_exposure=1e-10,
            beta_outcome=bx * 0.5 * unit, se_outcome=0.02 * unit, pval_outcome=0.01,
            f_statistic=64.0, n_exposure=n_exp, n_outcome=n_out))
    return out


def test_the_steiger_verdict_is_invariant_to_the_units_of_the_outcome():
    """Rescaling the outcome changes no science, so it must change no verdict.

    `by -> k*by, sy -> k*sy` is a change of units -- mmol/L to mg/dL. The old
    variance-explained term `2*eaf*(1-eaf)*beta^2` has no trait variance in it, so it
    scaled with k and flipped `correct` from True to False between k=1 and k=3, with p
    astronomically small on BOTH sides. z-statistics are unit-free, so the verdict is
    now constant across five orders of magnitude.
    """
    verdicts = {u: steiger_test(_steiger_input(unit=u))[0]
                for u in (0.1, 1.0, 3.0, 10.0, 100.0, 1000.0)}
    assert set(verdicts.values()) == {True}, verdicts


def test_steiger_reports_a_p_value_only_when_it_has_the_sample_sizes():
    """Without n there is no basis for one, and the old 0.01 had no derivation."""
    correct, p, note = steiger_test(_steiger_input())
    assert correct is True and p is None and "no sample sizes" in note

    correct, p, note = steiger_test(_steiger_input(n_exp=100_000, n_out=100_000))
    assert correct is True
    assert p is not None and 0.0 <= p <= 1.0
    assert note == ""


def test_steiger_detects_a_genuinely_reversed_direction():
    """The false-positive direction: it must still be able to say REVERSED."""
    insts = _steiger_input()
    reversed_insts = [
        Instrument(snp=i.snp, effect_allele=i.effect_allele, other_allele=i.other_allele,
                   eaf=i.eaf, beta_exposure=i.beta_outcome, se_exposure=i.se_outcome,
                   pval_exposure=i.pval_outcome, beta_outcome=i.beta_exposure,
                   se_outcome=i.se_exposure, pval_outcome=i.pval_exposure,
                   f_statistic=i.f_statistic)
        for i in insts
    ]
    assert steiger_test(insts)[0] is True
    assert steiger_test(reversed_insts)[0] is False


# ---------------------------------------------------------------------------
# Weighted mode: the SE must react to the dispersion of the ratios
# ---------------------------------------------------------------------------

def _mode_input(spread, n=15, seed=3):
    rng = np.random.default_rng(seed)
    out = []
    for k in range(n):
        bx = float(rng.uniform(0.05, 0.4))
        ratio = 0.5 + float(rng.normal(0, spread))
        out.append(Instrument(f"rs{k}", "A", "G", 0.3, bx, bx / 8, 1e-10,
                              bx * ratio, 0.02, 0.01, 64.0))
    return out


def test_the_weighted_mode_se_reflects_the_ratios_and_the_closed_form_did_not():
    """The old SE was `2 / sum(|bx|/sy)` -- the instrument count and the outcome standard
    errors, and nothing about the ratios whose mode it reports.

    So it returns THE SAME NUMBER for instruments that all agree and for instruments
    that split into two opposed clusters: 0.01333 for both datasets below, whose ratios
    are one tight cluster around 0.50 and an 8-versus-7 split between 0.20 and 0.80. The
    parametric bootstrap of Hartwig et al. 2017 reacts to the ratios, and it is several
    times wider than the closed form in both cases, so the closed form was also
    over-confident.

    What is NOT asserted is which of the two gets the wider SE. Measured with the
    reference bandwidth rule over three seeds, the tight cluster came out at 0.057 to
    0.063 and the split at 0.040 to 0.045: a two-camp split with a decisive majority pins
    the mode on the larger camp, whose ratios nearly coincide, while the ratios of a
    tight cluster scatter under resampling by more than their spread. An earlier version
    of this test asserted the opposite ordering; it held only for a bandwidth rule the
    reference method does not use, which is exactly the kind of invariant a test should
    not carry.

    The ratios carry a small jitter so no two coincide: exact ties make the median
    absolute deviation zero, and the reference rule then floors the bandwidth at 1e-8,
    which is the reference behaviour but not a realistic dataset.
    """
    rng = np.random.default_rng(2026)

    def _jittered(ratios):
        out = []
        for k, r in enumerate(ratios):
            bx = float(rng.uniform(0.15, 0.25))
            out.append(Instrument(f"rs{k}", "A", "G", 0.3, bx, bx / 8, 1e-10,
                                  bx * r, 0.02, 0.01, 64.0))
        return out

    tight = weighted_mode(_jittered(0.5 + rng.normal(0, 0.02, 15)), n_boot=300)
    split = weighted_mode(_jittered(np.r_[0.2 + rng.normal(0, 0.02, 8),
                                          0.8 + rng.normal(0, 0.02, 7)]), n_boot=300)

    # the OLD formula, evaluated on the same inputs, is blind to the ratios entirely
    old_se = 1.0 / (sum(1.0 / (0.02 / 0.2) for _ in range(15)) * 0.5)
    assert tight.se > 2 * old_se and split.se > 2 * old_se, (old_se, tight.se, split.se)
    assert abs(tight.se - split.se) > 0.2 * min(tight.se, split.se), (tight.se, split.se)
    assert tight.estimate == pytest.approx(0.5, abs=0.06)
    assert split.estimate == pytest.approx(0.2, abs=0.03)     # the larger camp


def test_the_weighted_mode_bandwidth_scales_with_the_data_not_a_constant():
    """A fixed 0.5 stops the mode being a mode.

    Fifteen instruments split 9 to 6 between ratios near 0.20 and 0.80 (jittered so no
    two coincide). The mode is 0.20, by construction -- that is the larger camp, and
    separating it from the weighted mean is the entire reason to run this estimator.
    The inverse-variance weighted mean here is 0.4463.

    With the old absolute bandwidth of 0.5, larger than the 0.6 gap between the camps,
    the kernels merge and the "mode" comes out at 0.3139, most of the way to the mean.
    With the bandwidth rule of Hartwig et al. 2017 eq. 7 it recovers the larger camp at
    phi = 0.3, 0.5 and at the reference default of 1 (0.1962 at each; the rule gives a
    bandwidth of 0.0092 at phi = 1 on these ratios).
    """
    from mendelian_randomisation import _weighted_mode_point, mbe_bandwidth

    rng = np.random.default_rng(7)
    ratios_in = np.r_[0.2 + rng.normal(0, 0.01, 9), 0.8 + rng.normal(0, 0.01, 6)]
    insts = []
    rng2 = np.random.default_rng(3)
    for k, r in enumerate(ratios_in):
        bx = float(rng2.uniform(0.15, 0.25))
        insts.append(Instrument(f"rs{k}", "A", "G", 0.3, bx, bx / 8, 1e-10,
                                bx * r, 0.02, 0.01, 64.0))
    bx = np.array([i.beta_exposure for i in insts]); by = np.array([i.beta_outcome for i in insts])
    sx = np.array([i.se_exposure for i in insts]); sy = np.array([i.se_outcome for i in insts])
    ratios = by / bx
    se_r = np.sqrt(sy ** 2 / bx ** 2 + by ** 2 * sx ** 2 / bx ** 4)
    weights = (1 / se_r ** 2) / np.sum(1 / se_r ** 2)
    mean = ivw(insts).estimate

    absolute = _weighted_mode_point(ratios, weights, 0.5)
    assert abs(absolute - mean) < 0.15, (absolute, mean)          # pulled onto the mean
    for phi in (0.3, 0.5, 1.0):
        proportional = _weighted_mode_point(ratios, weights, mbe_bandwidth(ratios, phi))
        assert proportional == pytest.approx(0.20, abs=0.02), (phi, proportional)
        assert abs(proportional - mean) > 0.2

    # and the shipped estimator is scale-free, which the absolute bandwidth was not
    scaled = [Instrument(i.snp, i.effect_allele, i.other_allele, i.eaf, i.beta_exposure,
                         i.se_exposure, i.pval_exposure, i.beta_outcome * 100,
                         i.se_outcome * 100, i.pval_outcome, i.f_statistic)
              for i in insts]
    a = weighted_mode(insts, n_boot=100).estimate
    b = weighted_mode(scaled, n_boot=100).estimate / 100
    assert a == pytest.approx(b, rel=0.05), (a, b)


# ---------------------------------------------------------------------------
# Weighted mode and Steiger: the terms that tie the code to the cited methods.
# Each of these pins ONE term a plausible re-implementation could silently change
# without any of the behavioural tests above noticing.
# ---------------------------------------------------------------------------

def test_the_weighted_mode_weights_are_inverse_variance_not_inverse_se():
    """Hartwig et al. 2017 eq. 5: `w_j = se_j^-2 / sum(se^-2)`.

    Five imprecise instruments at a ratio of 0.20 (se 0.02) against one precise
    instrument at 0.80 (se 0.005). Inverse-variance weights put 40000 on the precise
    one against 5 x 2500 = 12500 on the camp, so the mode is 0.80. Inverse-SE weights,
    which this function used to apply, give 200 against 5 x 50 = 250, and the mode
    flips to 0.20. Same data, opposite answer; only the exponent differs.
    """
    from mendelian_randomisation import _weighted_mode_point

    ratios = np.array([0.2, 0.201, 0.199, 0.2, 0.2005, 0.8])
    se = np.array([0.02] * 5 + [0.005])
    inv_var = 1.0 / se ** 2
    weights = inv_var / inv_var.sum()
    # a wide, fixed bandwidth so the answer is about the weights and nothing else
    assert _weighted_mode_point(ratios, weights, 0.05) == pytest.approx(0.8, abs=0.02)

    # and the shipped estimator reproduces that, so the weights it builds are these
    insts = [Instrument(f"rs{k}", "A", "G", 0.3, 0.2, 1e-6, 1e-10, 0.2 * r, 0.2 * s_, 0.01, 64.0)
             for k, (r, s_) in enumerate(zip(ratios, se))]
    est = weighted_mode(insts, phi=6.0, n_boot=50)   # phi large enough to match the 0.05 above
    assert est.estimate == pytest.approx(0.8, abs=0.03), est.estimate


def test_the_weighted_mode_se_is_the_scaled_mad_of_the_bootstrap_and_p_is_from_t():
    """The reference implementation reports `1.4826 * mad(draws)` and a t-test on L - 1
    degrees of freedom (TwoSampleMR `mr_mode`: `se_Mode <- apply(beta_Mode.boot, 2,
    stats::mad)`, `P_Mode <- pt(..., df = length(b_exp) - 1) * 2`)."""
    from mendelian_randomisation import _mad

    assert _mad(np.array([1.0, 2.0, 3.0, 4.0, 100.0])) == pytest.approx(1.4826 * 1.0)

    # The SE goes THROUGH that helper. Pinned by substitution rather than by
    # re-deriving the bootstrap here: swap `_mad` for a sentinel and the reported SE
    # must be the sentinel's value. Without this, `np.std(draws)` in its place passed
    # every other test in this file (mutation-checked 2026-09-12).
    import mendelian_randomisation as mr_mod
    real_mad = mr_mod._mad
    mr_mod._mad = lambda draws: 0.123456
    try:
        assert weighted_mode(_mode_input(0.15, n=12), n_boot=20).se == pytest.approx(0.123456)
    finally:
        mr_mod._mad = real_mad

    insts = _mode_input(0.15, n=12)
    est = weighted_mode(insts, n_boot=200)
    expected_p = 2 * stats.t.sf(abs(est.estimate / est.se), df=len(insts) - 1)
    assert est.pvalue == pytest.approx(expected_p, rel=1e-9)
    # not the normal approximation: at df = 11 the two differ by more than rounding
    assert est.pvalue != pytest.approx(2 * stats.norm.sf(abs(est.estimate / est.se)), rel=1e-3)


def test_steiger_p_value_follows_the_reference_conversion_and_aggregation():
    """Pins the two terms that make this the TwoSampleMR Steiger test rather than a
    near neighbour: `r2 = z^2 / (z^2 + n - 2)` per SNP (`get_r_from_bsen`) and the
    per-SNP sample sizes aggregated by their MEAN (`mr_steiger`), then Fisher's z on
    two independent correlations with `1/(n-3)` variances.

    Small sample sizes on purpose: at n = 20 the `-2` moves r2 by 10% for z = 1, so a
    drift to `z^2/(z^2+n)` is visible here and invisible at 100,000. Kept small enough
    that the summed r2 on each side stays below 1 (the production code clamps at 1 for
    atanh, and the expectation below would otherwise need the same clamp)."""
    insts = [
        Instrument("rs1", "A", "G", 0.3, 0.2, 0.2, 1e-10, 0.05, 0.1, 0.01, 64.0,
                   n_exposure=20, n_outcome=40),
        Instrument("rs2", "A", "G", 0.3, 0.3, 0.2, 1e-10, 0.05, 0.1, 0.01, 64.0,
                   n_exposure=30, n_outcome=50),
    ]
    correct, p, note = steiger_test(insts)
    assert correct is True and note == ""

    z_exp = np.array([0.2 / 0.2, 0.3 / 0.2]); z_out = np.array([0.5, 0.5])
    n_exp = np.array([20.0, 30.0]); n_out = np.array([40.0, 50.0])
    r_exp = math.sqrt(float(np.sum(z_exp ** 2 / (z_exp ** 2 + n_exp - 2))))
    r_out = math.sqrt(float(np.sum(z_out ** 2 / (z_out ** 2 + n_out - 2))))
    se = math.sqrt(1 / (n_exp.mean() - 3) + 1 / (n_out.mean() - 3))
    expected = 2 * stats.norm.sf(abs((math.atanh(r_exp) - math.atanh(r_out)) / se))
    assert p == pytest.approx(expected, rel=1e-9)


def test_the_weighted_mode_ratio_se_carries_the_exposure_uncertainty():
    """The reference weighted mode uses the second-order delta-method ratio SE,
    `sqrt(sy^2/bx^2 + by^2*sx^2/bx^4)` (TwoSampleMR `mr_mode`, the column "not assuming
    NOME"). The first-order `sy/|bx|` ignores the exposure SE entirely.

    Five precise instruments at a ratio of 0.20 against one at 0.80 whose outcome SE is
    tiny but whose exposure SE is enormous. Second order: the exposure term makes its
    ratio SE large, it is down-weighted, the mode is 0.20. First order: its ratio SE is
    tiny, it dominates, the mode is 0.80.
    """
    camp = [Instrument(f"rs{k}", "A", "G", 0.3, 0.2, 0.01, 1e-10, 0.2 * r, 0.02, 0.01, 64.0)
            for k, r in enumerate([0.2, 0.201, 0.199, 0.2, 0.2005])]
    loud = Instrument("rs9", "A", "G", 0.3, 0.2, 0.5, 1e-10, 0.16, 0.0005, 0.01, 64.0)
    est = weighted_mode(camp + [loud], phi=6.0, n_boot=50)
    assert est.estimate == pytest.approx(0.2, abs=0.03), est.estimate



# ---------------------------------------------------------------------------
# The empty instrument set, and the two boundaries the guards above leave unpinned
#
# Same class of defect as the rest of this file -- an estimator producing a number
# where the method does not apply -- at the one instrument count the guards above do
# not cover. `MIN_SENSITIVITY_INSTRUMENTS` gates n = 1 and n = 2; nothing gates n = 0,
# which is reachable from the CLI with a well-formed input file whose `instruments`
# array is empty, and from any upstream step (p-value threshold, LD clumping,
# harmonisation) that filters every SNP out.
#
# Measured on the code before this section: `ivw([])` returned
# `applicable=True, estimate=nan, se=inf`, and the pipeline printed all four
# "Running ..." lines before dying with an unhandled ZeroDivisionError from
# `compute_i_squared_gx`, having written no report and no JSON.
# ---------------------------------------------------------------------------


def test_ivw_is_not_applicable_without_instruments():
    """IVW is defined down to a single instrument, where it is the Wald ratio, and no
    further. It previously returned `applicable=True` with `estimate=nan, se=inf`,
    which is the one shape `MREstimate` exists to make impossible."""
    est = ivw([])

    assert est.applicable is False
    assert est.n_snps == 0
    assert "at least one instrument" in est.reason
    assert math.isnan(est.estimate) and math.isnan(est.se)


def test_i_squared_gx_returns_zero_rather_than_dividing_by_zero():
    """`(q - df) / q` with q = 0. At n = 0 the existing `q <= df or df == 0` test
    misses it, because df is -1: `0.0 <= -1` is False and `-1 == 0` is False, so both
    ways past the guard are open and the division runs. This was the unhandled
    ZeroDivisionError that ended the pipeline."""
    assert compute_i_squared_gx([]) == 0.0
    # n = 1 (df = 0) already returned 0.0 and must continue to
    assert compute_i_squared_gx(_egger_input([0.3])) == 0.0


def test_steiger_does_not_assume_sample_sizes_it_was_never_given():
    """`all(v is not None and v > 3 for v in n_exp + n_out)` over an empty list is
    vacuously True, so the no-instrument case took the branch that requires sample
    sizes, reached `np.mean([])`, and returned a NaN p-value -- which the JSON writer's
    `allow_nan=False` then refuses, turning a bad input into a crash at the last step."""
    correct, p, note = steiger_test([])

    assert p is None, f"a p-value was reported from no instruments: {p!r}"
    assert correct is False
    assert "no sample sizes supplied" in note


def test_the_pipeline_refuses_an_empty_instrument_set_before_writing_anything(tmp_path):
    """Zero instruments is not an analysis whose estimators are unavailable; it is the
    absence of the analysis. `applicable=False` is the right answer for an estimator
    that cannot run on a real instrument set, but there is no report to write when
    there are no instruments at all, so the pipeline declines instead of writing one.

    The message names the likely cause, because the realistic route here is not an
    empty file but upstream filtering that removed every SNP.
    """
    out = tmp_path / "run"
    with pytest.raises(NoInstrumentsError) as excinfo:
        run_pipeline([], exposure="X", outcome="Y", output_dir=out, demo=True)

    assert "at least one" in str(excinfo.value)
    assert "harmonis" in str(excinfo.value)
    # It refuses before creating the output tree, so a failed run leaves no artefacts
    # to be mistaken for a completed one.
    assert not out.exists(), sorted(p.name for p in out.iterdir()) if out.exists() else None
    # ValueError, so callers already catching that keep working.
    assert isinstance(excinfo.value, ValueError)


def test_the_cli_refuses_an_empty_instrument_file_without_a_traceback(tmp_path):
    """The artefact-level statement. A well-formed file with an empty `instruments`
    array used to exit on an unhandled ZeroDivisionError with a stack trace, after
    reporting `0 instruments loaded` and four estimators "Running". A user reading
    that cannot tell a bad input from a broken skill."""
    import subprocess

    bad = tmp_path / "empty.json"
    bad.write_text(json.dumps({"exposure": "X", "outcome": "Y", "instruments": []}))

    proc = subprocess.run(
        [sys.executable, str(SKILL_DIR / "mendelian_randomisation.py"),
         "--instruments", str(bad), "--output", str(tmp_path / "out")],
        capture_output=True, text=True,
    )

    assert proc.returncode != 0
    assert "Traceback" not in proc.stderr, proc.stderr
    assert "ZeroDivisionError" not in proc.stderr
    assert "at least one" in proc.stderr
    assert not (tmp_path / "out" / "report.md").exists()


@pytest.mark.parametrize("estimator", [weighted_median, weighted_mode],
                         ids=["weighted_median", "weighted_mode"])
@pytest.mark.parametrize("n", [1, 2])
def test_the_order_statistic_estimators_are_not_applicable_below_three(estimator, n):
    """Their n < 3 guards were reachable only through the single two-instrument
    end-to-end test above, so n = 1 was never exercised on either of them directly."""
    est = estimator(_egger_input([0.2, 0.5][:n]))

    assert est.applicable is False
    assert est.n_snps == n
    assert f"at least {MIN_SENSITIVITY_INSTRUMENTS} instruments" in est.reason
    assert math.isnan(est.estimate) and math.isnan(est.se)


@pytest.mark.parametrize("estimator", [weighted_median, weighted_mode],
                         ids=["weighted_median", "weighted_mode"])
def test_the_order_statistic_estimators_apply_at_exactly_three(estimator):
    """The boundary in the other direction, which nothing currently pins: every test
    that expects a number from these two uses 12, 15 or 30 instruments, so widening
    the guard to `< 4` would not fail a single one of them. Three is the count the
    reference admits, so three must produce an estimate."""
    est = estimator(_egger_input([0.2, 0.35, 0.5]))

    assert est.applicable is True, est.reason
    assert est.n_snps == MIN_SENSITIVITY_INSTRUMENTS
    assert math.isfinite(est.estimate) and math.isfinite(est.se)
    assert est.estimate == pytest.approx(0.5, abs=0.05)


def test_a_report_that_skips_one_estimator_still_corroborates_across_the_others(tmp_path):
    """The mixed case, which no existing test reaches: one estimator declines while
    the others stand.

    Every end-to-end case above is all-or-nothing -- the demo, where all four apply,
    and the two-instrument run, where only IVW does. The report's comparison
    partitions the estimators into `comparable` and `skipped`, and the branch where
    both are non-empty is never executed. It has to name the estimators that actually
    ran rather than the fixed list of four, and it must not fall through to "the IVW
    estimate stands alone", which is true only when nothing else applied.

    Five instruments with near-identical exposure effects: Egger's slope is
    unidentified, the median and the mode are ordinary order statistics of five
    well-defined ratios.
    """
    insts = _egger_input([0.4 + k * 1e-11 for k in range(5)])
    out = tmp_path / "run"
    out.mkdir()
    run_pipeline(insts, exposure="X", outcome="Y", output_dir=out, demo=True)

    result = json.loads((out / "result.json").read_text())
    applicable = {e["method"]: e.get("applicable", True) for e in result["estimates"]}
    assert applicable == {"IVW": True, "MR-Egger": False,
                          "Weighted Median": True, "Weighted Mode": True}, applicable

    report = (out / "report.md").read_text()
    assert "MR-Egger was not computed" in report
    assert "| MR-Egger | not computed |" in report
    # Names the three that ran, and does not claim agreement with the one that did not.
    assert ("consistent estimates across IVW, Weighted Median, Weighted Mode"
            in report), [l for l in report.splitlines() if "consistent" in l]
    assert "stands alone" not in report
