"""Tests for the coloc.abf arithmetic.

Four layers:

1. **An independent re-derivation.** The posterior combine is computed a second
   time directly, without logarithms, on inputs small enough not to overflow, and
   the two must agree. This catches an algebra slip in the log-sum-exp and
   log-difference path, which no amount of self-consistency testing would.
2. **Parity with R coloc.** A frozen fixture holds the five posteriors R coloc
   5.2.3 produced for one fixed input. The R package is GPL and ClawBio is MIT,
   so it is never installed, imported or vendored here; the test reads the stored
   numbers (see `capture_r_reference.R` for how to re-derive them).
3. **Properties, over generated inputs.** That the posteriors are a distribution,
   that they do not depend on which allele each trait reported its effect for,
   that they do not depend on the order of the rows, and that a degenerate
   standard error raises rather than returning something confident.
4. **Constructed controls.** A region built with one shared causal variant must
   come out as H4, and one built with two distinct causal variants as H3.
"""

from __future__ import annotations

import json
import math
import sys
from pathlib import Path

import numpy as np
import pytest
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from coloc_abf_math import (  # noqa: E402
    CC_PRIOR_SD,
    QUANT_PRIOR_SD,
    ColocDataset,
    ColocError,
    coloc_sensitivity,
    combine_abf,
    estimate_sdy,
    per_snp_log_abf,
    prior_variance,
    run_coloc_abf,
)

R_REFERENCE = Path(__file__).parent / "fixtures" / "coloc_abf_r_reference.json"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _combine_directly(l1, l2, p1, p2, p12):
    """The posterior combine again, in linear space, with no shared code.

    Valid only while the Bayes factors are small enough to exponentiate, which
    the callers below ensure. H3 is the off-diagonal mass, hence the subtraction.
    """
    bf1 = sum(math.exp(v) for v in l1)
    bf2 = sum(math.exp(v) for v in l2)
    bf_both = sum(math.exp(a + b) for a, b in zip(l1, l2))
    unnormalised = [
        1.0,
        p1 * bf1,
        p2 * bf2,
        p1 * p2 * (bf1 * bf2 - bf_both),
        p12 * bf_both,
    ]
    total = sum(unnormalised)
    return {f"H{i}": v / total for i, v in enumerate(unnormalised)}


def _region(n=120, causal1=60, causal2=60, peak1=8.0, peak2=7.0,
            se1=0.03, se2=0.02, decay=6.0, seed=11):
    """Two traits over one window, each peaking at the index it is given."""
    rng = np.random.default_rng(seed)
    idx = np.arange(n)
    snp = [f"1_{100_000 + 500 * i}_A_G" for i in idx]
    z1 = peak1 * np.exp(-np.abs(idx - causal1) / decay) + rng.normal(0, 1, n)
    z2 = peak2 * np.exp(-np.abs(idx - causal2) / decay) + rng.normal(0, 1, n)
    d1 = ColocDataset(snp=snp, beta=z1 * se1, se=np.full(n, se1),
                      trait_type="quant", sdY=1.0, label="trait one")
    d2 = ColocDataset(snp=snp, beta=z2 * se2, se=np.full(n, se2),
                      trait_type="quant", sdY=1.0, label="trait two")
    return d1, d2


# ---------------------------------------------------------------------------
# 1. The independent re-derivation
# ---------------------------------------------------------------------------

class TestAgainstADirectRederivation:

    def test_combine_matches_the_linear_space_computation(self):
        l1 = [0.4, 2.1, 5.0, 0.2, -0.3, 1.1]
        l2 = [0.1, 1.7, 4.4, 0.3, 0.0, 0.9]
        got = combine_abf(np.array(l1), np.array(l2), p1=1e-4, p2=1e-4, p12=1e-5)
        want = _combine_directly(l1, l2, 1e-4, 1e-4, 1e-5)
        for key in ("H0", "H1", "H2", "H3", "H4"):
            assert got[key] == pytest.approx(want[key], rel=1e-12, abs=1e-18)

    def test_combine_matches_when_there_is_no_signal(self):
        l1 = [-0.2, -0.1, 0.0, -0.4]
        l2 = [0.1, -0.3, -0.2, 0.0]
        got = combine_abf(np.array(l1), np.array(l2))
        want = _combine_directly(l1, l2, 1e-4, 1e-4, 1e-5)
        for key in ("H0", "H1", "H2", "H3", "H4"):
            assert got[key] == pytest.approx(want[key], rel=1e-12, abs=1e-18)
        assert got["H0"] > 0.99

    def test_per_variant_bayes_factor_matches_the_formula(self):
        beta, se, w = 0.18, 0.035, (0.15 * 1.0) ** 2
        v = se ** 2
        r = w / (w + v)
        want = 0.5 * (math.log(1 - r) + r * (beta / se) ** 2)
        got = float(per_snp_log_abf(np.array([beta]), np.array([se]), w)[0])
        assert got == pytest.approx(want, rel=1e-14)


# ---------------------------------------------------------------------------
# 2. Parity with R coloc
# ---------------------------------------------------------------------------

class TestParityWithRColoc:

    @staticmethod
    def _run_fixture(fixture):
        inp = fixture["inputs"]
        d1 = ColocDataset(snp=inp["snp"], beta=np.array(inp["beta1"]),
                          se=np.array(inp["se1"]), trait_type=inp["type1"],
                          sdY=inp["sdY1"])
        d2 = ColocDataset(snp=inp["snp"], beta=np.array(inp["beta2"]),
                          se=np.array(inp["se2"]), trait_type=inp["type2"],
                          sdY=inp["sdY2"])
        r = run_coloc_abf(d1, d2, p1=inp["p1"], p2=inp["p2"], p12=inp["p12"])
        return r.pp

    def test_matches_the_captured_r_posteriors(self):
        fixture = json.loads(R_REFERENCE.read_text())
        reference = fixture["r_reference"]
        assert reference["present"], (
            "the R reference is marked absent; re-capture it with "
            "tests/capture_r_reference.R"
        )
        got = self._run_fixture(fixture)
        for key in ("PP.H0", "PP.H1", "PP.H2", "PP.H3", "PP.H4"):
            assert got[key] == pytest.approx(reference[key], rel=1e-4, abs=1e-6), (
                f"{key}: this implementation gives {got[key]}, R coloc "
                f"{reference['coloc_version']} gave {reference[key]}"
            )

    def test_output_has_not_drifted_from_the_frozen_pin(self):
        fixture = json.loads(R_REFERENCE.read_text())
        got = self._run_fixture(fixture)
        for key, value in fixture["our_reference"].items():
            assert got[key] == pytest.approx(value, rel=1e-9, abs=1e-15), (
                f"{key}: {got[key]} has drifted from the frozen {value}. The code "
                "changed, so explain why before refreshing the fixture."
            )


# ---------------------------------------------------------------------------
# 3. Properties
# ---------------------------------------------------------------------------

_BETA = st.floats(min_value=-1.0, max_value=1.0, allow_nan=False, allow_infinity=False)
_SE = st.floats(min_value=0.005, max_value=0.5, allow_nan=False, allow_infinity=False)
_SIGN = st.sampled_from([1.0, -1.0])


def _datasets(betas1, ses1, betas2, ses2):
    n = len(betas1)
    snp = [f"1_{1000 + i}_A_G" for i in range(n)]
    return (
        ColocDataset(snp=snp, beta=np.array(betas1), se=np.array(ses1),
                     trait_type="quant", sdY=1.0, label="one"),
        ColocDataset(snp=snp, beta=np.array(betas2), se=np.array(ses2),
                     trait_type="quant", sdY=1.0, label="two"),
    )


@st.composite
def _two_traits(draw, min_size=2, max_size=25):
    n = draw(st.integers(min_value=min_size, max_value=max_size))
    betas1 = draw(st.lists(_BETA, min_size=n, max_size=n))
    ses1 = draw(st.lists(_SE, min_size=n, max_size=n))
    betas2 = draw(st.lists(_BETA, min_size=n, max_size=n))
    ses2 = draw(st.lists(_SE, min_size=n, max_size=n))
    return betas1, ses1, betas2, ses2


class TestProperties:

    @given(_two_traits())
    @settings(deadline=None, max_examples=60,
              suppress_health_check=[HealthCheck.too_slow])
    def test_posteriors_are_a_probability_distribution(self, traits):
        d1, d2 = _datasets(*traits)
        pp = run_coloc_abf(d1, d2).pp
        for key, value in pp.items():
            assert 0.0 <= value <= 1.0, f"{key} outside [0, 1]: {value}"
            assert math.isfinite(value), f"{key} is not finite: {value}"
        assert sum(pp.values()) == pytest.approx(1.0, rel=1e-9)

    @given(_two_traits(), st.lists(_SIGN, min_size=25, max_size=25),
           st.lists(_SIGN, min_size=25, max_size=25))
    @settings(deadline=None, max_examples=60,
              suppress_health_check=[HealthCheck.too_slow])
    def test_posteriors_ignore_which_allele_each_trait_reported(
        self, traits, flips1, flips2
    ):
        """The Bayes factor reads the effect only as z squared.

        So flipping the effect allele of any variant, in either trait,
        independently, cannot change the posteriors. Anyone harmonising alleles
        purely to satisfy this computation is doing work it does not need.
        """
        betas1, ses1, betas2, ses2 = traits
        n = len(betas1)
        d1, d2 = _datasets(betas1, ses1, betas2, ses2)
        flipped1 = [b * s for b, s in zip(betas1, flips1[:n])]
        flipped2 = [b * s for b, s in zip(betas2, flips2[:n])]
        f1, f2 = _datasets(flipped1, ses1, flipped2, ses2)
        assert run_coloc_abf(d1, d2).pp == run_coloc_abf(f1, f2).pp

    @given(_two_traits(min_size=3), st.integers(min_value=0, max_value=10_000))
    @settings(deadline=None, max_examples=60,
              suppress_health_check=[HealthCheck.too_slow])
    def test_posteriors_ignore_the_order_of_the_rows(self, traits, seed):
        betas1, ses1, betas2, ses2 = traits
        d1, d2 = _datasets(betas1, ses1, betas2, ses2)
        order = np.random.default_rng(seed).permutation(len(betas1))
        s1, s2 = _datasets(
            [betas1[i] for i in order], [ses1[i] for i in order],
            [betas2[i] for i in order], [ses2[i] for i in order],
        )
        # The shuffled datasets carry shuffled identifiers too, so the two runs
        # colocalise the same variants in a different order.
        shuffled = ColocDataset(
            snp=[f"1_{1000 + i}_A_G" for i in order],
            beta=s1.beta, se=s1.se, trait_type="quant", sdY=1.0, label="one")
        shuffled2 = ColocDataset(
            snp=[f"1_{1000 + i}_A_G" for i in order],
            beta=s2.beta, se=s2.se, trait_type="quant", sdY=1.0, label="two")
        a = run_coloc_abf(d1, d2).pp
        b = run_coloc_abf(shuffled, shuffled2).pp
        for key in a:
            assert a[key] == pytest.approx(b[key], rel=1e-9, abs=1e-15)

    @pytest.mark.parametrize("bad_se", [0.0, -0.01, float("nan"), float("inf")])
    def test_a_degenerate_standard_error_raises(self, bad_se):
        """Never a confident answer from a standard error that cannot be used."""
        snp = ["1_1_A_G", "1_2_A_G", "1_3_A_G"]
        d1 = ColocDataset(snp=snp, beta=np.array([0.1, 0.2, 0.3]),
                          se=np.array([0.02, bad_se, 0.03]),
                          trait_type="quant", sdY=1.0, label="one")
        d2 = ColocDataset(snp=snp, beta=np.array([0.1, 0.2, 0.3]),
                          se=np.array([0.02, 0.02, 0.03]),
                          trait_type="quant", sdY=1.0, label="two")
        with pytest.raises(ColocError):
            run_coloc_abf(d1, d2)

    @pytest.mark.parametrize("bad_beta", [float("nan"), float("inf")])
    def test_a_non_finite_effect_raises(self, bad_beta):
        snp = ["1_1_A_G", "1_2_A_G"]
        d1 = ColocDataset(snp=snp, beta=np.array([0.1, bad_beta]),
                          se=np.array([0.02, 0.02]),
                          trait_type="quant", sdY=1.0, label="one")
        d2 = ColocDataset(snp=snp, beta=np.array([0.1, 0.2]),
                          se=np.array([0.02, 0.02]),
                          trait_type="quant", sdY=1.0, label="two")
        with pytest.raises(ColocError):
            run_coloc_abf(d1, d2)


# ---------------------------------------------------------------------------
# 4. Constructed controls
# ---------------------------------------------------------------------------

class TestConstructedControls:

    def test_a_shared_causal_variant_gives_h4(self):
        d1, d2 = _region(causal1=60, causal2=60)
        result = run_coloc_abf(d1, d2)
        assert result.h4 > 0.9, result.pp
        assert result.h4 > result.h3
        assert result.n_snps == 120
        assert result.best_snp == d1.snp[60]

    def test_distinct_causal_variants_give_h3(self):
        d1, d2 = _region(causal1=30, causal2=90)
        result = run_coloc_abf(d1, d2)
        assert result.h3 > 0.9, result.pp
        assert result.h3 > result.h4

    def test_only_one_trait_associated_gives_h1_or_h2(self):
        d1, d2 = _region(causal1=60, causal2=60, peak2=0.0)
        result = run_coloc_abf(d1, d2)
        assert result.h1 > 0.5, result.pp


# ---------------------------------------------------------------------------
# Priors, sdY and the intersection
# ---------------------------------------------------------------------------

class TestPriorsAndSdY:

    def test_prior_variance_follows_the_trait_type(self):
        assert prior_variance("cc", 99.0) == pytest.approx(CC_PRIOR_SD ** 2)
        assert prior_variance("quant", 2.0) == pytest.approx((QUANT_PRIOR_SD * 2.0) ** 2)

    def test_sdy_estimate_recovers_a_known_trait_sd(self):
        rng = np.random.default_rng(3)
        n_samples = 40_000
        maf = rng.uniform(0.05, 0.5, 400)
        true_sdy = 1.7
        se = true_sdy / np.sqrt(n_samples * 2.0 * maf * (1.0 - maf))
        assert estimate_sdy(se, maf, n_samples) == pytest.approx(true_sdy, rel=1e-6)

    def test_a_quantitative_trait_without_sdy_or_frequencies_raises(self):
        snp = ["1_1_A_G", "1_2_A_G"]
        d = ColocDataset(snp=snp, beta=np.array([0.1, 0.2]), se=np.array([0.02, 0.02]),
                         trait_type="quant", label="no sdY")
        with pytest.raises(ColocError, match="needs sdY"):
            run_coloc_abf(d, d)

    def test_a_case_control_trait_needs_no_sdy(self):
        snp = ["1_1_A_G", "1_2_A_G", "1_3_A_G"]
        d = ColocDataset(snp=snp, beta=np.array([0.1, 0.4, 0.05]),
                         se=np.array([0.02, 0.02, 0.02]),
                         trait_type="cc", label="case control")
        result = run_coloc_abf(d, d)
        assert result.sdy_used["trait1"] is None
        assert result.prior_sd["trait1"] == pytest.approx(CC_PRIOR_SD)

    def test_a_repeated_identifier_is_refused(self):
        with pytest.raises(ColocError, match="more than once"):
            ColocDataset(snp=["1_1_A_G", "1_1_A_G"], beta=np.array([0.1, 0.1]),
                         se=np.array([0.02, 0.02]), trait_type="cc", label="dupes")

    def test_no_shared_variants_is_refused(self):
        d1 = ColocDataset(snp=["1_1_A_G"], beta=np.array([0.1]), se=np.array([0.02]),
                          trait_type="cc", label="one")
        d2 = ColocDataset(snp=["2_1_A_G"], beta=np.array([0.1]), se=np.array([0.02]),
                          trait_type="cc", label="two")
        with pytest.raises(ColocError, match="no shared variants"):
            run_coloc_abf(d1, d2)

    def test_only_the_intersection_is_used_and_said_so(self):
        snp1 = ["1_1_A_G", "1_2_A_G", "1_3_A_G", "1_4_A_G"]
        snp2 = ["1_2_A_G", "1_3_A_G"]
        d1 = ColocDataset(snp=snp1, beta=np.array([0.1, 0.4, 0.3, 0.05]),
                          se=np.full(4, 0.02), trait_type="cc", label="one")
        d2 = ColocDataset(snp=snp2, beta=np.array([0.4, 0.3]),
                          se=np.full(2, 0.02), trait_type="cc", label="two")
        result = run_coloc_abf(d1, d2)
        assert result.n_snps == 2
        assert any("shared variants" in note for note in result.notes)

    @pytest.mark.parametrize("bad", [0.0, 1.0, -1e-5, 2.0])
    def test_a_prior_outside_the_unit_interval_is_refused(self, bad):
        l = np.array([1.0, 2.0])
        with pytest.raises(ColocError, match="between 0 and 1"):
            combine_abf(l, l, p12=bad)


class TestSensitivity:

    def test_a_strong_shared_signal_is_prior_robust(self):
        d1, d2 = _region(causal1=60, causal2=60, peak1=10.0, peak2=9.0)
        sens = coloc_sensitivity(d1, d2)
        assert sens.passes_at_default
        assert sens.robust
        assert not sens.prior_driven
        assert all(pt.rule_passes for pt in sens.points)

    def test_no_signal_fails_the_rule_at_every_prior(self):
        d1, d2 = _region(causal1=30, causal2=90)
        sens = coloc_sensitivity(d1, d2)
        assert not sens.passes_at_default
        assert sens.passing_p12_range is None
        assert "No colocalisation" in " ".join(sens.notes)

    def test_a_borderline_call_is_reported_as_prior_driven(self):
        """A verdict that changes inside the swept range must say so.

        The parameters are chosen to land PP.H4 near the rule at the default
        p12, which is exactly the case a single posterior would hide.
        """
        found = None
        for peak2 in (4.0, 4.2, 4.4, 4.6, 4.8, 5.0, 5.2, 5.4):
            d1, d2 = _region(causal1=60, causal2=60, peak1=6.0, peak2=peak2, seed=5)
            sens = coloc_sensitivity(d1, d2)
            if sens.prior_driven:
                found = sens
                break
        assert found is not None, (
            "no borderline case found in the swept peak heights; if the "
            "arithmetic changed, re-pick the parameters"
        )
        assert found.prior_driven
        assert not found.robust
        if found.passes_at_default:
            assert "prior-driven" in found.verdict

    def test_the_default_prior_is_evaluated_exactly(self):
        d1, d2 = _region()
        sens = coloc_sensitivity(d1, d2, p12_default=3e-6)
        assert any(pt.p12 == pytest.approx(3e-6) for pt in sens.points)

    def test_the_sweep_agrees_with_running_the_combine_at_that_prior(self):
        d1, d2 = _region()
        sens = coloc_sensitivity(d1, d2)
        for point in sens.points:
            direct = run_coloc_abf(d1, d2, p12=point.p12)
            assert direct.h4 == pytest.approx(point.h4, rel=1e-12)
            assert direct.h3 == pytest.approx(point.h3, rel=1e-12)
