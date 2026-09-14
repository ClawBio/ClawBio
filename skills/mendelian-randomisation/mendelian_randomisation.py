#!/usr/bin/env python3
"""Mendelian Randomisation — two-sample MR with full sensitivity analysis.

Implements IVW, MR-Egger, weighted median, and weighted mode estimators with
Cochran's Q, Egger intercept, Steiger directionality, F-statistic diagnostics,
leave-one-out analysis, and publication-ready visualisation.

References:
    Burgess et al. (2013) Genet Epidemiol 37:658-665 (IVW)
    Bowden et al. (2015) Int J Epidemiol 44:512-525 (MR-Egger)
    Bowden et al. (2016) Genet Epidemiol 40:304-314 (Weighted median)
    Verbanck et al. (2018) Nature Genetics 50:693-698 (MR-PRESSO)

Usage:
    python mendelian_randomisation.py --demo --output /tmp/mr_demo
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path

import numpy as np
from scipy import stats

SCRIPT_DIR = Path(__file__).resolve().parent

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

MIN_F_STAT = 10
EAF_PALINDROME_THRESHOLD = 0.42


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class Instrument:
    snp: str
    effect_allele: str
    other_allele: str
    eaf: float
    beta_exposure: float
    se_exposure: float
    pval_exposure: float
    beta_outcome: float
    se_outcome: float
    pval_outcome: float
    f_statistic: float = 0.0
    # Sample sizes, optional and read from the input when present. The Steiger
    # directionality test needs them to express each association as a variance
    # explained; without them it can still order the two sides, but only under an
    # assumption it then has to state. See `steiger_test`.
    n_exposure: int | None = None
    n_outcome: int | None = None

    @property
    def is_palindromic(self) -> bool:
        pairs = {frozenset({"A", "T"}), frozenset({"C", "G"})}
        return frozenset({self.effect_allele, self.other_allele}) in pairs

    @property
    def palindromic_ambiguous(self) -> bool:
        return self.is_palindromic and EAF_PALINDROME_THRESHOLD < self.eaf < (1 - EAF_PALINDROME_THRESHOLD)

    @property
    def weak_instrument(self) -> bool:
        return self.f_statistic < MIN_F_STAT


# The three sensitivity estimators need at least this many instruments. MR-Egger fits
# a slope AND an intercept, so below 3 there is no residual degree of freedom; the
# weighted median and weighted mode are order statistics of the per-SNP ratios, and
# TwoSampleMR (Hemani et al. 2018, the reference implementation of all three) returns
# NA for each of them below 3 SNPs. IVW is the one estimator defined at n=1, where it
# reduces to the single Wald ratio.
MIN_SENSITIVITY_INSTRUMENTS = 3
MIN_EGGER_INSTRUMENTS = MIN_SENSITIVITY_INSTRUMENTS  # kept for callers that import it

# IVW is the only estimator here defined at n = 1. At n = 0 every weighted sum is zero,
# so the estimate is 0/0 and the standard error divides by zero; nothing above gates it,
# because `MIN_SENSITIVITY_INSTRUMENTS` starts at 1.
MIN_IVW_INSTRUMENTS = 1

# Dimensionless conditioning floor for the Egger slope, replacing an absolute one.
# rel = Var_w(bx) / E_w[bx^2] lies in [0, 1] and is invariant to the units of the data.
# Below sqrt(machine epsilon) the Egger SE is inflated by more than ~8,000x, so the
# estimate is uninformative rather than merely imprecise.
EGGER_MIN_RELATIVE_VARIANCE = math.sqrt(sys.float_info.epsilon)  # ~1.49e-8


class NoInstrumentsError(ValueError):
    """Raised when an analysis is requested with no instruments at all.

    Distinct from `applicable=False`, which describes an estimator that cannot be
    computed on a real instrument set and still leaves a report worth writing. Zero
    instruments is the absence of the analysis itself, so there is nothing to report
    and the pipeline declines before it writes anything. Subclasses `ValueError` so
    callers already catching that keep working.
    """


@dataclass
class MREstimate:
    method: str
    estimate: float
    se: float
    ci_lower: float
    ci_upper: float
    pvalue: float
    n_snps: int
    # An estimator can be UNDEFINED on a given instrument set rather than merely
    # imprecise. Without somewhere to say that, the only options are to raise (which
    # discards the estimates already computed correctly) or to emit a number that reads
    # exactly like a result. `applicable=False` carries `reason` instead, and every
    # consumer below skips the row and prints the reason in its place.
    applicable: bool = True
    reason: str = ""

    @classmethod
    def not_applicable(cls, method: str, n_snps: int, reason: str) -> "MREstimate":
        """An estimator that does not apply to this instrument set.

        The numeric fields are NaN on purpose: any consumer that ignores `applicable`
        and formats them anyway produces a visible `nan` rather than a plausible number,
        and the JSON writer refuses to serialise them at all.
        """
        return cls(method=method, estimate=float("nan"), se=float("nan"),
                   ci_lower=float("nan"), ci_upper=float("nan"),
                   pvalue=float("nan"), n_snps=n_snps,
                   applicable=False, reason=reason)


@dataclass
class SensitivityResults:
    cochran_q: float = 0.0
    cochran_q_pvalue: float = 1.0
    cochran_q_df: int = 0
    egger_intercept: float = 0.0
    egger_intercept_se: float = 0.0
    egger_intercept_pvalue: float = 1.0
    mean_f_statistic: float = 0.0
    min_f_statistic: float = 0.0
    n_weak_instruments: int = 0
    i_squared_gx: float = 0.0
    steiger_correct_direction: bool = True
    # None when the test could order the two sides but had no basis for a p-value; the
    # reason is in `steiger_note`, and every writer prints that instead of a number.
    steiger_pvalue: float | None = None
    steiger_note: str = ""


# ---------------------------------------------------------------------------
# MR estimators
# ---------------------------------------------------------------------------
def ivw(instruments: list[Instrument]) -> MREstimate:
    """Inverse-Variance Weighted estimator (multiplicative random effects).

    Defined down to a single instrument, where it is the Wald ratio, and no further.
    With none, the weighted sums are all zero and the estimate is 0/0 with an infinite
    standard error, so it is declared rather than returned.
    """
    if len(instruments) < MIN_IVW_INSTRUMENTS:
        return MREstimate.not_applicable(
            "IVW", len(instruments),
            f"IVW requires at least one instrument; this analysis has {len(instruments)}")

    bx = np.array([i.beta_exposure for i in instruments])
    by = np.array([i.beta_outcome for i in instruments])
    sy = np.array([i.se_outcome for i in instruments])

    w = 1.0 / (sy ** 2)
    beta_ivw = np.sum(w * bx * by) / np.sum(w * bx ** 2)

    residuals = by - beta_ivw * bx
    if len(instruments) == 1:
        phi = 1.0
    else:
        phi = max(1.0, np.sum(w * residuals ** 2) / (len(instruments) - 1))
    se_ivw = math.sqrt(phi / np.sum(w * bx ** 2))

    z = beta_ivw / se_ivw
    pval = 2 * stats.norm.sf(abs(z))

    return MREstimate(
        method="IVW", estimate=float(beta_ivw), se=float(se_ivw),
        ci_lower=float(beta_ivw - 1.96 * se_ivw),
        ci_upper=float(beta_ivw + 1.96 * se_ivw),
        pvalue=float(pval), n_snps=len(instruments),
    )


def egger_weighted_dispersion(w: np.ndarray, bx: np.ndarray) -> float:
    """`sum_w * sum_i w_i (bx_i - xbar_w)^2` -- the MR-Egger slope denominator.

    Its own function so the property it exists for can be tested directly. That
    property is STRUCTURAL: every term is a product of non-negative numbers, so the
    result cannot be negative in IEEE-754, and it is exactly zero precisely when every
    `bx` is equal.

    The algebraically identical expanded form, `sum_w*sum_wbx2 - sum_wbx**2`, has
    neither guarantee. It is a difference of two large nearly-equal quantities, so as
    the exposure effects converge it collapses onto a cancellation remainder whose sign
    is arbitrary -- and a negative one reaches `math.sqrt` as a domain error. That is
    not visible through `mr_egger` once the conditioning check is in place, because the
    check refuses those inputs anyway; it is visible here.
    """
    sum_w = np.sum(w)
    xbar_w = np.sum(w * bx) / sum_w
    return float(sum_w * np.sum(w * (bx - xbar_w) ** 2))


def mr_egger(instruments: list[Instrument]) -> tuple[MREstimate, float, float, float]:
    """MR-Egger regression. Returns (estimate, intercept, intercept_se, intercept_p).

    Bowden J, Davey Smith G, Burgess S 2015, Int J Epidemiol 44(2):512-525
    (doi:10.1093/ije/dyv080; PMID 26050253).

    Two INDEPENDENT conditions have to hold, and neither implies the other:

    STATISTICAL -- at least `MIN_EGGER_INSTRUMENTS`. Egger fits two parameters, so below
    three there is no residual degree of freedom, no matter how clean the data is.

    NUMERICAL -- at least two distinct `beta_exposure` values, checked as a dimensionless
    ratio. The slope is unidentified when every `bx` coincides, and that is not an n < 3
    problem: with a relative spread around 1e-10 the shipped code raised on roughly a
    third of well-sized inputs at n = 5 and n = 10.

    Below either, the estimator does not apply, and it says so instead of returning a
    number. Returning one was the failure this replaces: at n = 1, over 5,000 draws on
    each of two effect-size distributions, the old code returned a finite slope with a
    median SE of 2.1e+07 in ~27% of cases, raised `math domain error` in ~27%, and hit
    its own guard in ~46% -- a one-ulp sign coin flip rather than a property of the data.
    """
    bx = np.array([i.beta_exposure for i in instruments])
    by = np.array([i.beta_outcome for i in instruments])
    sy = np.array([i.se_outcome for i in instruments])

    # Orient every instrument so its exposure effect is positive, flipping the
    # outcome effect with it (the same variant, other allele). The slope is
    # invariant to this; the INTERCEPT is not -- it is the mean outcome effect at
    # zero exposure effect, so on un-oriented inputs it depends on which allele
    # each GWAS happened to report, i.e. on the allele coding rather than on the
    # data. This is what TwoSampleMR's mr_egger_regression does before its fit
    # (b_out <- b_out * sign(b_exp); b_exp <- abs(b_exp)). R's sign(0) is 0, so an
    # instrument with an exposure effect of exactly zero has its outcome effect
    # zeroed as well there; np.sign does the same, which keeps the two fits equal on
    # that input too (an earlier version kept such an instrument's outcome effect,
    # with a comment claiming the reference did likewise; it does not).
    orient = np.sign(bx)
    bx = np.abs(bx)
    by = by * orient

    w = 1.0 / (sy ** 2)
    n = len(instruments)

    if n < MIN_EGGER_INSTRUMENTS:
        return (
            MREstimate.not_applicable(
                "MR-Egger", n,
                f"MR-Egger requires at least {MIN_EGGER_INSTRUMENTS} instruments "
                f"(it fits a slope and an intercept); this analysis has {n}"),
            float("nan"), float("nan"), float("nan"),
        )

    sum_w = np.sum(w)
    sum_wbx = np.sum(w * bx)
    sum_wbx2 = np.sum(w * bx ** 2)
    sum_wby = np.sum(w * by)

    # Centered (Lagrange) form: denom = sum_w * sum_i w_i (bx_i - xbar_w)^2, a sum of
    # non-negative terms that is zero exactly when every bx is equal. The algebraically
    # identical expanded form, sum_w*sum_wbx2 - sum_wbx**2, is a difference of two large
    # nearly-equal quantities, so in floating point it lands on a cancellation remainder
    # whose SIGN is arbitrary -- which is where the negative values under the square root
    # came from. The NUMERATOR is centered for the same reason: centering the denominator
    # alone fixes the sign and leaves the slope's accuracy degrading as the spread
    # narrows.
    xbar_w = sum_wbx / sum_w
    ybar_w = sum_wby / sum_w
    dx = bx - xbar_w
    denom = egger_weighted_dispersion(w, bx)
    numer = sum_w * np.sum(w * dx * (by - ybar_w))

    # Dimensionless conditioning test. `denom` carries the data's units: rescaling the
    # outcome alone, which changes the slope but not the conditioning at all, moves it
    # across hundreds of orders of magnitude, so no absolute threshold can be a criterion.
    # This ratio is invariant to that rescaling.
    scale = sum_w * sum_wbx2
    rel = denom / scale if scale > 0 else 0.0
    # Written as a negated `>` rather than `rel < threshold` so NaN also fails it. A
    # comparison against NaN is False either way, and only this direction turns that
    # into a refusal rather than into passing the check.
    if not (rel > EGGER_MIN_RELATIVE_VARIANCE):
        return (
            MREstimate.not_applicable(
                "MR-Egger", n,
                "MR-Egger is not identified on these instruments: their exposure "
                f"effects are too close to identical (relative variance {rel:.2e}, "
                f"below {EGGER_MIN_RELATIVE_VARIANCE:.2e}), so the slope has no "
                "informative standard error"),
            float("nan"), float("nan"), float("nan"),
        )

    slope = numer / denom
    intercept = ybar_w - slope * xbar_w

    fitted = intercept + slope * bx
    residuals = by - fitted
    # df cannot be <= 0 here, since n >= MIN_EGGER_INSTRUMENTS above. Guarded anyway:
    # this function is public and directly callable, and that is the difference between
    # "unreachable from our pipeline" and "cannot happen".
    df = n - 2
    phi = 1.0 if df <= 0 else max(1.0, np.sum(w * residuals ** 2) / df)

    se_slope = math.sqrt(phi * sum_w / denom)
    se_intercept = math.sqrt(phi * sum_wbx2 / denom)

    # t reference on n - 2 df, not a normal: both standard errors come from the
    # fit's residual variance, itself an estimate on n - 2 degrees of freedom
    # (two parameters fitted), and the t carries that uncertainty. At n = 3 there
    # is one residual degree of freedom and the p-value is correspondingly wide;
    # a normal here reported p = 0.006 on three demo instruments where the t
    # gives 0.22. TwoSampleMR's mr_egger_regression uses pt(., n - 2) for both.
    # IVW keeps its normal reference (its SE is the analytic one, no residual
    # variance estimated), as in that implementation.
    t_slope = slope / se_slope
    p_slope = 2 * stats.t.sf(abs(t_slope), df=df)
    t_int = intercept / se_intercept
    p_int = 2 * stats.t.sf(abs(t_int), df=df)

    estimate = MREstimate(
        method="MR-Egger", estimate=float(slope), se=float(se_slope),
        ci_lower=float(slope - 1.96 * se_slope),
        ci_upper=float(slope + 1.96 * se_slope),
        pvalue=float(p_slope), n_snps=n,
    )
    return estimate, float(intercept), float(se_intercept), float(p_int)


def weighted_median(instruments: list[Instrument], n_boot: int = 1000) -> MREstimate:
    """Weighted median estimator (Bowden et al., 2016).

    Not applicable below `MIN_SENSITIVITY_INSTRUMENTS`: the median of one or two ratios
    is not a robust estimator of anything, and reporting it alongside IVW at n=1 let the
    report certify a single Wald ratio as "consistent across methods".
    """
    if len(instruments) < MIN_SENSITIVITY_INSTRUMENTS:
        return MREstimate.not_applicable(
            "Weighted Median", len(instruments),
            f"Weighted median requires at least {MIN_SENSITIVITY_INSTRUMENTS} instruments; "
            f"this analysis has {len(instruments)}")
    bx = np.array([i.beta_exposure for i in instruments])
    by = np.array([i.beta_outcome for i in instruments])
    sy = np.array([i.se_outcome for i in instruments])

    ratios = by / bx
    weights = 1.0 / (sy ** 2 / bx ** 2)
    weights = weights / np.sum(weights)

    order = np.argsort(ratios)
    ratios_sorted = ratios[order]
    weights_sorted = weights[order]
    cum_weights = np.cumsum(weights_sorted)
    idx = np.searchsorted(cum_weights, 0.5)
    beta_wm = float(ratios_sorted[min(idx, len(ratios_sorted) - 1)])

    rng = np.random.default_rng(42)
    boot_estimates = []
    for _ in range(n_boot):
        by_boot = by + rng.normal(0, sy)
        bx_boot = bx
        r_boot = by_boot / bx_boot
        o = np.argsort(r_boot)
        cw = np.cumsum(weights[o])
        j = np.searchsorted(cw, 0.5)
        boot_estimates.append(float(r_boot[o[min(j, len(r_boot) - 1)]]))
    se_wm = float(np.std(boot_estimates))
    z = beta_wm / se_wm if se_wm > 0 else 0
    pval = 2 * stats.norm.sf(abs(z))

    return MREstimate(
        method="Weighted Median", estimate=beta_wm, se=se_wm,
        ci_lower=beta_wm - 1.96 * se_wm, ci_upper=beta_wm + 1.96 * se_wm,
        pvalue=float(pval), n_snps=len(instruments),
    )


def _mad(x: np.ndarray) -> float:
    """Median absolute deviation, scaled by 1.4826 for consistency with the SD under
    normality (the constant R's `mad()` applies by default)."""
    med = float(np.median(x))
    return 1.4826 * float(np.median(np.abs(x - med)))


def mbe_bandwidth(ratios: np.ndarray, phi: float) -> float:
    """`h = phi * s`, with `s` the modified Silverman rule Hartwig et al. 2017 eq. 7 use:
    `s = 0.9 * min(sd, 1.4826 * mad) / L^(1/5)`, floored at 1e-8 as in the reference
    implementation so a set of identical ratios still has a width."""
    sd = float(np.std(ratios, ddof=1)) if len(ratios) > 1 else 0.0
    s = 0.9 * min(sd, _mad(ratios)) / len(ratios) ** 0.2
    return max(1e-8, s * phi)


def _weighted_mode_point(ratios: np.ndarray, weights: np.ndarray,
                         bandwidth: float) -> float:
    """The mode of the weighted normal-kernel density over the ratios (Hartwig et al.
    2017 eq. 6), with `weights` already standardised to sum to 1."""
    span = float(np.max(ratios) - np.min(ratios))
    pad = max(span, 3.0 * bandwidth) if span or bandwidth else 1.0
    x_grid = np.linspace(float(np.min(ratios)) - pad, float(np.max(ratios)) + pad, 2000)
    density = np.zeros_like(x_grid)
    for r, w in zip(ratios, weights):
        density += w * stats.norm.pdf(x_grid, loc=r, scale=bandwidth)
    return float(x_grid[int(np.argmax(density))])


def weighted_mode(instruments: list[Instrument], phi: float = 1.0,
                  n_boot: int = 1000, seed: int = 0) -> MREstimate:
    """Weighted mode estimator.

    Hartwig FP, Davey Smith G, Bowden J 2017, Int J Epidemiol 46(6):1985-1998
    (doi:10.1093/ije/dyx102; PMID 29040600), which specifies both a bandwidth
    proportional to the spread of the ratios and a bootstrapped standard error.

    This follows the paper's weighted MBE and its reference implementation
    (TwoSampleMR `mr_weighted_mode`, Hemani et al. 2018) term by term:

    - ratio SEs by the delta method, `sqrt(sy^2/bx^2 + by^2*sx^2/bx^4)` (the
      "not assuming NOME" column the reference uses for the weighted mode);
    - standardised inverse-variance weights, eq. 5: `w_j = se_j^-2 / sum(se^-2)`;
    - bandwidth `h = phi * s` with the modified Silverman rule, eq. 7:
      `s = 0.9 * min(sd, 1.4826*mad) / L^(1/5)`, and the reference's `phi = 1`;
    - standard error = 1.4826 x the median absolute deviation of a parametric
      bootstrap of the ratios (each ratio redrawn from N(ratio, se_ratio)), the
      bandwidth recomputed on every draw;
    - p-value from a t distribution on L - 1 degrees of freedom, as the reference does.

    Why this replaces what was here. THE BANDWIDTH WAS AN ABSOLUTE 0.5 whatever the data
    looked like. On instruments whose ratios have a standard deviation of 0.28 -- an
    ordinary MR scale -- 0.5 is nearly twice the entire spread, so the kernels merge into
    one blob and the "mode" slides onto the weighted mean: measured 0.5469 against an
    IVW estimate of 0.5281, while the same data at a bandwidth of 0.1 gives 0.6679. An
    estimator whose whole purpose is to disagree with the mean when most instruments
    agree with each other cannot have a smoothing width that swamps the disagreement.
    A data-scaled bandwidth also makes the estimator invariant to the units of the data,
    as the mean and the median already were.

    THE STANDARD ERROR IS A BOOTSTRAP. The previous `2 / sum(|bx_i|/sy_i)` is a function
    of the instrument count and the outcome standard errors and NOTHING ELSE -- it does
    not look at the dispersion of the ratios whose mode it is reporting. Measured over
    ratio standard deviations from 0.0000 to 0.8265, an 800-fold change in how spread
    out the estimates are, the reported SE did not move at all: 0.02079 at n=10 and
    0.00429 at n=40 in every case. It was most confident exactly where the point
    estimate was least stable -- at the widest dispersion the mode itself moved from
    0.9673 to 0.3958 between two samples whose SEs were 0.0208 and 0.0043.

    The bootstrap redraws each ratio from its own reported uncertainty and re-derives
    the mode, so the SE reflects how much the mode actually moves.

    The weights were also `1/se`; the paper's eq. 5 and the reference implementation use
    `1/se^2`. Corrected here, since the function now claims to be that estimator.
    """
    n = len(instruments)
    if n < MIN_SENSITIVITY_INSTRUMENTS:
        return MREstimate.not_applicable(
            "Weighted Mode", n,
            f"Weighted mode requires at least {MIN_SENSITIVITY_INSTRUMENTS} instruments; "
            f"this analysis has {n}")
    bx = np.array([i.beta_exposure for i in instruments])
    by = np.array([i.beta_outcome for i in instruments])
    sx = np.array([i.se_exposure for i in instruments])
    sy = np.array([i.se_outcome for i in instruments])

    ratios = by / bx
    # Delta-method ratio SE, second order ("not assuming NOME"), as the reference uses
    # for the weighted mode; the NOME variant drops the second term.
    se_ratios = np.sqrt(sy ** 2 / bx ** 2 + by ** 2 * sx ** 2 / bx ** 4)
    inv_var = 1.0 / se_ratios ** 2
    weights = inv_var / np.sum(inv_var)

    beta_mode = _weighted_mode_point(ratios, weights, mbe_bandwidth(ratios, phi))

    rng = np.random.default_rng(seed)
    draws = np.empty(n_boot)
    for b in range(n_boot):
        ratios_b = rng.normal(ratios, se_ratios)
        draws[b] = _weighted_mode_point(ratios_b, weights, mbe_bandwidth(ratios_b, phi))
    se_mode = _mad(draws)

    t_stat = beta_mode / se_mode if se_mode > 0 else 0.0
    pval = float(2 * stats.t.sf(abs(t_stat), df=n - 1))

    return MREstimate(
        method="Weighted Mode", estimate=beta_mode, se=se_mode,
        ci_lower=beta_mode - 1.96 * se_mode, ci_upper=beta_mode + 1.96 * se_mode,
        pvalue=float(pval), n_snps=len(instruments),
    )


# ---------------------------------------------------------------------------
# Sensitivity tests
# ---------------------------------------------------------------------------
def cochran_q(instruments: list[Instrument], ivw_est: MREstimate) -> tuple[float, float, int]:
    """Cochran's Q test for heterogeneity."""
    bx = np.array([i.beta_exposure for i in instruments])
    by = np.array([i.beta_outcome for i in instruments])
    sy = np.array([i.se_outcome for i in instruments])
    w = 1.0 / (sy ** 2)
    residuals = by - ivw_est.estimate * bx
    q = float(np.sum(w * residuals ** 2))
    df = len(instruments) - 1
    p = float(stats.chi2.sf(q, df)) if df > 0 else 1.0
    return q, p, df


def steiger_test(instruments: list[Instrument]) -> tuple[bool, float | None, str]:
    """Steiger directionality test. Returns (correct_direction, p_value_or_None, note).

    Hemani G, Tilling K, Davey Smith G 2017, PLoS Genet 13(11):e1007081
    (doi:10.1371/journal.pgen.1007081; PMID 29149188).

    Compares how much variance the instruments explain in the exposure against how much
    they explain in the outcome; more in the exposure supports exposure -> outcome.

    THE VARIANCE EXPLAINED IS COMPUTED FROM THE Z-STATISTIC, WHICH IS UNIT-FREE.
    The previous form, `r2 = 2*eaf*(1-eaf)*beta^2`, is a variance explained only if the
    trait happens to have variance 1: there is no division by the trait's variance and
    no sample size anywhere in it. So the verdict moved when the outcome was expressed
    in different units, which is not a scientific property of anything. Measured on ten
    instruments with a true ratio of 0.5, rescaling the outcome alone -- mmol/L to
    mg/dL, say -- flipped `correct` from True to False between factors of 1 and 3, with
    p astronomically small on BOTH sides, so it reported near-certainty in opposite
    directions depending on a unit choice.

    `r2 = z^2 / (z^2 + n - 2)` is the conversion for a continuous trait (the F statistic
    of a one-predictor regression on n - 2 residual degrees of freedom; TwoSampleMR
    `get_r_from_bsen`), and z is invariant to the units of beta because the standard
    error carries the same units. It is applied to BOTH sides. For a binary outcome
    reported in log odds, TwoSampleMR's `mr_steiger` instead converts the outcome side
    with `get_r_from_lor`, from the effect allele frequency, the case and control
    counts and the population prevalence; this input carries none of those, so the
    outcome's variance explained is on the continuous-trait scale here and the note
    says so whenever a p-value is reported.

    Sample sizes are OPTIONAL, and what is reported depends on what is available:

    - both present: a variance explained per side, and a p-value from the Fisher
      z-transform difference of two INDEPENDENT correlations -- which is what two-sample
      MR has by construction, the exposure and outcome coming from different studies.
      (The paper states the one-sample form, Steiger's Z for correlated correlations
      within one population; its two-sample implementation, TwoSampleMR `mr_steiger`,
      uses the independent-samples test, and so does this.)
    - absent: with equal sample sizes the comparison reduces to |z_exposure| >
      |z_outcome|, so the DIRECTION is still well defined and still unit-free. The
      p-value is not: it is returned as None with the assumption stated in `note`,
      rather than as a number from the previous `sqrt(r2_exp + r2_out) * 0.01`, whose
      0.01 has no derivation.
    """
    z_exp = np.array([i.beta_exposure / i.se_exposure if i.se_exposure else 0.0
                      for i in instruments])
    z_out = np.array([i.beta_outcome / i.se_outcome if i.se_outcome else 0.0
                      for i in instruments])

    n_exp = [i.n_exposure for i in instruments]
    n_out = [i.n_outcome for i in instruments]
    # `bool(instruments)` first: `all()` over an empty list is vacuously True, which
    # sent the no-instrument case down the branch that needs sample sizes, through
    # `np.mean([])`, and out with a NaN p-value.
    have_n = bool(instruments) and all(v is not None and v > 3 for v in n_exp + n_out)

    if not have_n:
        correct = bool(np.sum(z_exp ** 2) > np.sum(z_out ** 2))
        return correct, None, (
            "no sample sizes supplied, so the direction is read from the z-statistics "
            "under the assumption that the exposure and outcome studies are of "
            "comparable size")

    r2_exp = float(np.sum(z_exp ** 2 / (z_exp ** 2 + np.array(n_exp, dtype=float) - 2.0)))
    r2_out = float(np.sum(z_out ** 2 / (z_out ** 2 + np.array(n_out, dtype=float) - 2.0)))
    correct = r2_exp > r2_out

    # Fisher z on each side, then the difference of two independent correlations.
    # Clamped below 1 because atanh is undefined at exactly 1, which a very strong
    # instrument set can reach after summing.
    r_exp = min(math.sqrt(max(r2_exp, 0.0)), 1.0 - 1e-12)
    r_out = min(math.sqrt(max(r2_out, 0.0)), 1.0 - 1e-12)
    # Per-instrument sample sizes are aggregated by their mean, as TwoSampleMR
    # `mr_steiger` does (`n = mean(n_exp), n2 = mean(n_out)`).
    n1, n2 = float(np.mean(n_exp)), float(np.mean(n_out))
    se = math.sqrt(1.0 / (n1 - 3.0) + 1.0 / (n2 - 3.0))
    z_stat = (math.atanh(r_exp) - math.atanh(r_out)) / se
    p = float(2 * stats.norm.sf(abs(z_stat)))
    return correct, p, (
        "variance explained on both sides is converted from the z-statistic with the "
        "continuous-trait formula; a binary outcome in log odds would need case and "
        "control counts and the prevalence (TwoSampleMR get_r_from_lor), which this "
        "input does not carry")


def compute_i_squared_gx(instruments: list[Instrument]) -> float:
    """I² for instrument-exposure associations (Bowden et al., 2016)."""
    # `df <= 0` rather than `df == 0`, and BEFORE the weighted mean: with no instruments
    # df is -1, so neither `q <= df` (0.0 <= -1 is False) nor `df == 0` held and
    # `(q - df) / q` divided by a zero q. Returning here also keeps `bx_bar` from
    # evaluating 0/0, which is meaningless and emits a RuntimeWarning.
    df = len(instruments) - 1
    if df <= 0:
        return 0.0

    bx = np.array([i.beta_exposure for i in instruments])
    sx = np.array([i.se_exposure for i in instruments])
    w = 1.0 / (sx ** 2)
    bx_bar = np.sum(w * bx) / np.sum(w)
    q = float(np.sum(w * (bx - bx_bar) ** 2))
    if q <= df:
        return 0.0
    return max(0.0, float((q - df) / q))


def leave_one_out(instruments: list[Instrument]) -> list[tuple[str, MREstimate]]:
    """Leave-one-out IVW analysis."""
    results = []
    for i, inst in enumerate(instruments):
        subset = instruments[:i] + instruments[i + 1:]
        if len(subset) < 2:
            continue
        est = ivw(subset)
        results.append((inst.snp, est))
    return results


def run_sensitivity(instruments: list[Instrument], ivw_est: MREstimate) -> SensitivityResults:
    """Run full sensitivity analysis battery."""
    q, q_p, q_df = cochran_q(instruments, ivw_est)
    f_stats = [i.f_statistic for i in instruments]
    steiger_dir, steiger_p, steiger_note = steiger_test(instruments)
    i2_gx = compute_i_squared_gx(instruments)

    return SensitivityResults(
        cochran_q=q, cochran_q_pvalue=q_p, cochran_q_df=q_df,
        mean_f_statistic=float(np.mean(f_stats)),
        min_f_statistic=float(np.min(f_stats)),
        n_weak_instruments=sum(1 for f in f_stats if f < MIN_F_STAT),
        i_squared_gx=i2_gx,
        steiger_correct_direction=steiger_dir,
        steiger_pvalue=steiger_p,
        steiger_note=steiger_note,
    )


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------
def scatter_plot(instruments: list[Instrument], estimates: list[MREstimate], path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    bx = [i.beta_exposure for i in instruments]
    by = [i.beta_outcome for i in instruments]
    sx = [i.se_exposure for i in instruments]
    sy = [i.se_outcome for i in instruments]

    fig, ax = plt.subplots(figsize=(8, 7))
    ax.errorbar(bx, by, xerr=sx, yerr=sy, fmt="o", markersize=5, color="#2166ac",
                ecolor="#bdbdbd", alpha=0.7, capsize=0, label="Instruments")

    x_range = np.linspace(min(bx) - 0.01, max(bx) + 0.01, 100)
    colours = {"IVW": "#d32f2f", "MR-Egger": "#ff9800", "Weighted Median": "#4caf50", "Weighted Mode": "#9c27b0"}
    for est in estimates:
        # A not-applicable estimator has no slope to draw; plotting NaN silently omits
        # the line but still consumes a legend entry, which reads as "drawn at zero".
        if est.method in colours and est.applicable:
            ax.plot(x_range, est.estimate * x_range, color=colours[est.method],
                    linewidth=1.5, label=f"{est.method} ({est.estimate:.3f})")

    ax.axhline(0, color="grey", linewidth=0.5)
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_xlabel("SNP effect on exposure")
    ax.set_ylabel("SNP effect on outcome")
    ax.set_title("MR Scatter Plot")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def forest_plot(instruments: list[Instrument], ivw_est: MREstimate, path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    ratios = [i.beta_outcome / i.beta_exposure for i in instruments]
    se_ratios = [i.se_outcome / abs(i.beta_exposure) for i in instruments]
    labels = [i.snp for i in instruments]

    fig, ax = plt.subplots(figsize=(8, max(4, len(instruments) * 0.3)))
    y_pos = list(range(len(instruments)))

    ax.errorbar(ratios, y_pos, xerr=[1.96 * s for s in se_ratios], fmt="o", color="#2166ac",
                markersize=4, ecolor="#bdbdbd", capsize=0)
    ax.axvline(ivw_est.estimate, color="#d32f2f", linewidth=1.5, linestyle="--", label=f"IVW: {ivw_est.estimate:.3f}")
    ax.axvline(0, color="grey", linewidth=0.5)
    ax.set_yticks(y_pos)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("Causal estimate (Wald ratio)")
    ax.set_title("MR Forest Plot")
    ax.legend(fontsize=8)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def funnel_plot(instruments: list[Instrument], ivw_est: MREstimate, path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    ratios = [i.beta_outcome / i.beta_exposure for i in instruments]
    precision = [abs(i.beta_exposure) / i.se_outcome for i in instruments]

    fig, ax = plt.subplots(figsize=(7, 6))
    ax.scatter(ratios, precision, s=25, color="#2166ac", alpha=0.7)
    ax.axvline(ivw_est.estimate, color="#d32f2f", linewidth=1.5, linestyle="--", label=f"IVW: {ivw_est.estimate:.3f}")
    ax.set_xlabel("Causal estimate (Wald ratio)")
    ax.set_ylabel("Precision (|bx| / se_outcome)")
    ax.set_title("MR Funnel Plot")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


def leave_one_out_plot(loo_results: list[tuple[str, MREstimate]], ivw_all: MREstimate, path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    labels = [snp for snp, _ in loo_results] + ["All"]
    estimates = [est.estimate for _, est in loo_results] + [ivw_all.estimate]
    ci_lo = [est.ci_lower for _, est in loo_results] + [ivw_all.ci_lower]
    ci_hi = [est.ci_upper for _, est in loo_results] + [ivw_all.ci_upper]

    fig, ax = plt.subplots(figsize=(8, max(4, len(labels) * 0.3)))
    y = list(range(len(labels)))
    xerr_lo = [e - lo for e, lo in zip(estimates, ci_lo)]
    xerr_hi = [hi - e for e, hi in zip(estimates, ci_hi)]

    colours = ["#2166ac"] * len(loo_results) + ["#d32f2f"]
    for i, (est, lo, hi, c) in enumerate(zip(estimates, xerr_lo, xerr_hi, colours)):
        ax.errorbar(est, i, xerr=[[lo], [hi]], fmt="o", color=c, markersize=4, capsize=2)

    ax.axvline(ivw_all.estimate, color="#d32f2f", linewidth=0.8, linestyle=":")
    ax.set_yticks(y)
    ax.set_yticklabels(labels, fontsize=6)
    ax.set_xlabel("IVW estimate (leave-one-out)")
    ax.set_title("Leave-One-Out Analysis")
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
def generate_report(
    instruments: list[Instrument],
    estimates: list[MREstimate],
    sensitivity: SensitivityResults,
    egger_intercept: float,
    egger_intercept_p: float,
    loo: list[tuple[str, MREstimate]],
    exposure: str,
    outcome: str,
    output_dir: Path,
    demo: bool,
) -> None:
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    output_dir.mkdir(parents=True, exist_ok=True)
    for sub in ("tables", "figures", "reproducibility"):
        (output_dir / sub).mkdir(exist_ok=True)

    _write_mr_table(estimates, output_dir / "tables" / "mr_results.tsv")
    _write_sensitivity_table(sensitivity, egger_intercept, egger_intercept_p, output_dir / "tables" / "sensitivity.tsv")
    _write_instruments_table(instruments, output_dir / "tables" / "harmonised_instruments.tsv")
    _write_report_md(instruments, estimates, sensitivity, egger_intercept, egger_intercept_p, exposure, outcome, output_dir, ts, demo)
    _write_result_json(estimates, sensitivity, egger_intercept, egger_intercept_p, exposure, outcome, output_dir, ts, demo)
    _write_repro(output_dir, ts, demo)


def _write_mr_table(estimates, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["method", "estimate", "se", "ci_lower", "ci_upper", "pvalue", "n_snps", "note"])
        for e in estimates:
            if not e.applicable:
                # "not_applicable" in every numeric cell, never a formatted NaN: a
                # spreadsheet renders `nan` in an estimate column as a value someone
                # will try to read.
                w.writerow([e.method, "not_applicable", "not_applicable", "not_applicable",
                            "not_applicable", "not_applicable", e.n_snps, e.reason])
                continue
            w.writerow([e.method, f"{e.estimate:.6f}", f"{e.se:.6f}", f"{e.ci_lower:.6f}", f"{e.ci_upper:.6f}", f"{e.pvalue:.2e}", e.n_snps, ""])


def _steiger_interpretation(correct: bool, pvalue: float | None, *, markdown: bool) -> str:
    """The Steiger row's interpretation. With no sample sizes the direction is a
    comparison of summed z-squared, which stands in for variance explained only
    if the two studies are of comparable size: z-squared grows with n, so when
    the outcome GWAS is the larger study (the usual case for a disease outcome)
    its share is overstated and the comparison leans toward "reversed". There is
    no test behind it either way, so it is reported as consistent or reversed,
    NOT as confirmed: the comparison was made and gave a direction, and what is
    missing is its significance, which needs the sample sizes. "Confirmed" is
    reserved for a computed p-value."""
    if pvalue is None:
        if correct:
            return ("Direction consistent with exposure → outcome; significance not "
                    "assessable without sample sizes")
        return ("**WARNING: direction consistent with reverse causation; significance "
                "not assessable without sample sizes**" if markdown else
                "WARNING: direction consistent with reverse causation; significance "
                "not assessable without sample sizes")
    if correct:
        return "Exposure → Outcome confirmed"
    return "**WARNING: reverse causation**" if markdown else "WARNING: reversed causal direction"


def _write_sensitivity_table(s, egger_int, egger_p, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["test", "statistic", "pvalue", "interpretation"])
        w.writerow(["Cochran_Q", f"{s.cochran_q:.2f}", f"{s.cochran_q_pvalue:.4f}", "Significant = heterogeneity" if s.cochran_q_pvalue < 0.05 else "No significant heterogeneity"])
        if math.isnan(egger_int) or math.isnan(egger_p):
            note = "MR-Egger did not apply to this instrument set, so there is no directional-pleiotropy test"
            w.writerow(["Egger_intercept", "not_applicable", "not_applicable", note])
        else:
            w.writerow(["Egger_intercept", f"{egger_int:.6f}", f"{egger_p:.4f}", "Significant = directional pleiotropy" if egger_p < 0.05 else "No evidence of directional pleiotropy"])
        w.writerow(["Mean_F_statistic", f"{s.mean_f_statistic:.1f}", "N/A", f"{'WEAK' if s.mean_f_statistic < MIN_F_STAT else 'Strong'} instruments"])
        w.writerow(["Min_F_statistic", f"{s.min_f_statistic:.1f}", "N/A", f"{s.n_weak_instruments} weak instruments (F<{MIN_F_STAT})"])
        w.writerow(["I_squared_GX", f"{s.i_squared_gx:.4f}", "N/A", "SIMEX recommended" if s.i_squared_gx < 0.9 else "No SIMEX needed"])
        _sp = f"{s.steiger_pvalue:.4f}" if s.steiger_pvalue is not None else "not_applicable"
        _si = _steiger_interpretation(s.steiger_correct_direction, s.steiger_pvalue, markdown=False)
        w.writerow(["Steiger_direction", "Correct" if s.steiger_correct_direction else "REVERSED", _sp, f"{_si}{'; ' + s.steiger_note if s.steiger_note else ''}"])


def _write_instruments_table(instruments, path):
    with open(path, "w", newline="") as f:
        w = csv.writer(f, delimiter="\t")
        w.writerow(["SNP", "effect_allele", "other_allele", "eaf", "beta_exp", "se_exp", "pval_exp", "beta_out", "se_out", "pval_out", "f_stat", "palindromic", "weak"])
        for i in instruments:
            w.writerow([i.snp, i.effect_allele, i.other_allele, i.eaf, i.beta_exposure, i.se_exposure, f"{i.pval_exposure:.2e}", i.beta_outcome, i.se_outcome, f"{i.pval_outcome:.2e}", f"{i.f_statistic:.1f}", i.is_palindromic, i.weak_instrument])


def _write_report_md(instruments, estimates, sens, egger_int, egger_p, exposure, outcome, output_dir, ts, demo):
    n_palindromic = sum(1 for i in instruments if i.palindromic_ambiguous)
    lines = [
        "# Mendelian Randomisation Report", "",
        f"**Generated**: {ts}",
        f"**Exposure**: {exposure}",
        f"**Outcome**: {outcome}",
        f"**Instruments**: {len(instruments)} SNPs",
        f"**Mode**: {'Demo (cached data, offline)' if demo else 'Live'}", "",
        "## MR Estimates", "",
        "| Method | Estimate | SE | 95% CI | P-value |",
        "|--------|----------|----|--------|---------|",
    ]
    for e in estimates:
        if not e.applicable:
            # Never a formatted NaN: the TSV and the JSON already say "not
            # applicable" for this row, and the reason follows the table.
            lines.append(f"| {e.method} | not computed | not computed | not computed | not computed |")
            continue
        lines.append(f"| {e.method} | {e.estimate:.4f} | {e.se:.4f} | [{e.ci_lower:.4f}, {e.ci_upper:.4f}] | {e.pvalue:.2e} |")
    lines.append("")

    lines.extend([
        "## Sensitivity Analysis", "",
        "| Test | Result | P-value | Interpretation |",
        "|------|--------|---------|----------------|",
        f"| Cochran's Q | {sens.cochran_q:.2f} (df={sens.cochran_q_df}) | {sens.cochran_q_pvalue:.4f} | {'Heterogeneity detected' if sens.cochran_q_pvalue < 0.05 else 'No significant heterogeneity'} |",
        (f"| Egger intercept | {egger_int:.4f} | {egger_p:.4f} | "
         f"{'Directional pleiotropy' if egger_p < 0.05 else 'No directional pleiotropy'} |"
         if not (math.isnan(egger_int) or math.isnan(egger_p)) else
         "| Egger intercept | not computed | not computed | MR-Egger did not apply |"),
        f"| Mean F-statistic | {sens.mean_f_statistic:.1f} | — | {'**WARNING: weak instruments**' if sens.mean_f_statistic < MIN_F_STAT else 'Strong instruments'} |",
        f"| Weak instruments (F<{MIN_F_STAT}) | {sens.n_weak_instruments}/{len(instruments)} | — | {'**WARNING**' if sens.n_weak_instruments > 0 else 'None'} |",
        f"| I²_GX | {sens.i_squared_gx:.4f} | — | {'SIMEX correction recommended' if sens.i_squared_gx < 0.9 else 'Adequate'} |",
        (f"| Steiger direction | {'Correct' if sens.steiger_correct_direction else '**REVERSED**'} | "
         f"{f'{sens.steiger_pvalue:.4f}' if sens.steiger_pvalue is not None else 'not computed'} | "
         f"{_steiger_interpretation(sens.steiger_correct_direction, sens.steiger_pvalue, markdown=True)}"
         f"{'; ' + sens.steiger_note if sens.steiger_note else ''} |"),
        "",
    ])

    if n_palindromic > 0:
        lines.extend([f"**WARNING**: {n_palindromic} palindromic SNP(s) with ambiguous EAF (0.42–0.58) — these were retained but may introduce bias. Manual review recommended.", ""])

    lines.extend([
        "## Interpretation", "",
        f"The IVW estimate suggests a {'positive' if estimates[0].estimate > 0 else 'negative'} causal effect of {exposure} on {outcome} ",
        f"(beta = {estimates[0].estimate:.4f}, 95% CI [{estimates[0].ci_lower:.4f}, {estimates[0].ci_upper:.4f}], P = {estimates[0].pvalue:.2e}). ",
        "",
    ])
    # Compare only estimators that PRODUCED an estimate, and name the ones that did not.
    # A not-applicable row previously entered this comparison as a number, so a two-
    # instrument run whose Egger SE was infinite still certified the result "robust".
    comparable = [e for e in estimates[1:] if e.applicable]
    skipped = [e for e in estimates if not e.applicable]
    consistent = all(abs(e.estimate - estimates[0].estimate) < 2 * estimates[0].se
                     for e in comparable)
    if not comparable:
        lines.append("No sensitivity estimator applies to this instrument set, so the "
                     "IVW estimate stands alone and is not corroborated.")
    elif consistent:
        names = ", ".join(["IVW"] + [e.method for e in comparable])
        lines.append(f"Sensitivity analyses show consistent estimates across {names}, "
                     "supporting a robust causal inference.")
    else:
        lines.append("**Caution**: Estimates differ across methods, suggesting potential violations of MR assumptions. Interpret with care.")
    for e in skipped:
        lines.append(f"\n{e.method} was not computed: {e.reason}.")
    lines.extend(["", "---", "", f"*{DISCLAIMER}*", ""])
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def _write_result_json(estimates, sens, egger_int, egger_p, exposure, outcome, output_dir, ts, demo):
    result = {
        "tool": "ClawBio Mendelian Randomisation",
        "version": "0.1.0",
        "timestamp": ts,
        "mode": "demo" if demo else "live",
        "exposure": exposure,
        "outcome": outcome,
        "estimates": [
            {"method": e.method, "estimate": round(e.estimate, 6), "se": round(e.se, 6),
             "pvalue": f"{e.pvalue:.2e}", "n_snps": e.n_snps}
            if e.applicable else
            {"method": e.method, "applicable": False, "reason": e.reason,
             "n_snps": e.n_snps}
            for e in estimates
        ],
        "sensitivity": {
            "cochran_q": round(sens.cochran_q, 2), "cochran_q_p": round(sens.cochran_q_pvalue, 4),
            # None, not NaN: `allow_nan=False` below would otherwise refuse to write the
            # file at all on exactly the runs this change exists to make well-behaved.
            # JSON null is the honest encoding of "there is no such test here".
            "egger_intercept": None if math.isnan(egger_int) else round(egger_int, 6),
            "egger_intercept_p": None if math.isnan(egger_p) else round(egger_p, 4),
            "mean_f_stat": round(sens.mean_f_statistic, 1),
            "n_weak": sens.n_weak_instruments,
            "i_squared_gx": round(sens.i_squared_gx, 4),
            "steiger_correct": sens.steiger_correct_direction,
            "steiger_p": sens.steiger_pvalue,
            "steiger_note": sens.steiger_note or None,
        },
        "disclaimer": DISCLAIMER,
    }
    # allow_nan=False: Python emits `Infinity` and `NaN` bare by default, and RFC 8259
    # has neither, so a strict parser rejects the WHOLE file rather than one field. A
    # two-instrument run used to write `"se": Infinity` and exit 0. Now any non-finite
    # anywhere in the document raises here instead of shipping unparseable JSON -- for
    # every future one, not only this one.
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, allow_nan=False), encoding="utf-8")


def _write_repro(output_dir, ts, demo):
    d = output_dir / "reproducibility"
    d.mkdir(exist_ok=True)
    (d / "commands.sh").write_text(f"#!/usr/bin/env bash\n# Generated: {ts}\npython mendelian_randomisation.py {'--demo' if demo else ''} --output {output_dir}\n", encoding="utf-8")
    (d / "software_versions.json").write_text(json.dumps({"python": sys.version, "numpy": np.__version__, "scipy": stats.scipy.__version__ if hasattr(stats, 'scipy') else "unknown", "generated": ts}, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Pipeline
# ---------------------------------------------------------------------------
def load_demo_instruments() -> tuple[list[Instrument], str, str]:
    cache = SCRIPT_DIR / "example_data" / "demo_instruments.json"
    with open(cache) as f:
        data = json.load(f)
    instruments = [Instrument(
        snp=s["SNP"], effect_allele=s["effect_allele"], other_allele=s["other_allele"],
        eaf=s["eaf"], beta_exposure=s["beta_exposure"], se_exposure=s["se_exposure"],
        pval_exposure=s["pval_exposure"], beta_outcome=s["beta_outcome"],
        se_outcome=s["se_outcome"], pval_outcome=s["pval_outcome"],
        f_statistic=s["f_statistic"],
        n_exposure=s.get("n_exposure"), n_outcome=s.get("n_outcome"),
    ) for s in data["instruments"]]
    return instruments, data["exposure"], data["outcome"]


def run_pipeline(instruments: list[Instrument], exposure: str, outcome: str, output_dir: Path, demo: bool = False) -> dict:
    # Before the output tree is created, so a refused run leaves nothing behind that
    # could be mistaken for a completed one.
    if not instruments:
        raise NoInstrumentsError(
            "No instruments to analyse. Mendelian randomisation needs at least one "
            "harmonised instrument for the exposure. Check that the input file's "
            "`instruments` array is populated, and that upstream filtering "
            "(p-value threshold, LD clumping, harmonisation) has not removed every SNP."
        )

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    weak = [i for i in instruments if i.weak_instrument]
    if weak:
        print(f"  WARNING: {len(weak)} instrument(s) with F-statistic < {MIN_F_STAT} — weak instrument bias possible", file=sys.stderr)

    palindromic = [i for i in instruments if i.palindromic_ambiguous]
    if palindromic:
        print(f"  WARNING: {len(palindromic)} palindromic SNP(s) with ambiguous EAF — flagged for manual review", file=sys.stderr)

    print(f"[MR] Running IVW ({len(instruments)} instruments)...")
    ivw_est = ivw(instruments)

    print("[MR] Running MR-Egger...")
    egger_est, egger_int, egger_int_se, egger_int_p = mr_egger(instruments)

    print("[MR] Running Weighted Median...")
    wm_est = weighted_median(instruments)

    print("[MR] Running Weighted Mode...")
    wmode_est = weighted_mode(instruments)

    estimates = [ivw_est, egger_est, wm_est, wmode_est]

    print("[MR] Running sensitivity analysis...")
    sens = run_sensitivity(instruments, ivw_est)
    sens.egger_intercept = egger_int
    sens.egger_intercept_se = egger_int_se
    sens.egger_intercept_pvalue = egger_int_p

    print("[MR] Leave-one-out analysis...")
    loo = leave_one_out(instruments)

    print("[MR] Generating plots...")
    scatter_plot(instruments, estimates, output_dir / "figures" / "scatter.png")
    forest_plot(instruments, ivw_est, output_dir / "figures" / "forest.png")
    funnel_plot(instruments, ivw_est, output_dir / "figures" / "funnel.png")
    leave_one_out_plot(loo, ivw_est, output_dir / "figures" / "leave_one_out.png")

    print("[MR] Generating report...")
    generate_report(instruments, estimates, sens, egger_int, egger_int_p, loo, exposure, outcome, output_dir, demo)

    print(f"[MR] IVW estimate: {ivw_est.estimate:.4f} (P={ivw_est.pvalue:.2e})")
    print(f"[MR] Report written to: {output_dir / 'report.md'}")

    return {
        "ivw_estimate": ivw_est.estimate,
        "ivw_pvalue": ivw_est.pvalue,
        "n_instruments": len(instruments),
        "cochran_q_p": sens.cochran_q_pvalue,
        "egger_intercept_p": egger_int_p,
        "n_weak": sens.n_weak_instruments,
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Mendelian Randomisation — two-sample MR")
    parser.add_argument("--demo", action="store_true", help="Run with cached BMI->T2D demo data (offline)")
    parser.add_argument("--instruments", type=str, help="JSON file with harmonised instruments")
    parser.add_argument("--output", type=str, required=True, help="Output directory")

    args = parser.parse_args()

    if args.demo:
        instruments, exposure, outcome = load_demo_instruments()
    elif args.instruments:
        with open(args.instruments) as f:
            data = json.load(f)
        instruments = [Instrument(
            snp=s["SNP"], effect_allele=s["effect_allele"], other_allele=s["other_allele"],
            eaf=s["eaf"], beta_exposure=s["beta_exposure"], se_exposure=s["se_exposure"],
            pval_exposure=s["pval_exposure"], beta_outcome=s["beta_outcome"],
            se_outcome=s["se_outcome"], pval_outcome=s["pval_outcome"],
            f_statistic=s["f_statistic"],
            n_exposure=s.get("n_exposure"), n_outcome=s.get("n_outcome"),
        ) for s in data["instruments"]]
        exposure = data.get("exposure", "Exposure")
        outcome = data.get("outcome", "Outcome")
    else:
        parser.error("Provide --demo or --instruments <json>")

    output_dir = Path(args.output)
    print(f"[MR] Starting MR pipeline: {exposure} -> {outcome}")
    print(f"[MR] {len(instruments)} instruments loaded ({'demo/cached' if args.demo else 'user-provided'})")

    try:
        run_pipeline(instruments, exposure, outcome, output_dir, demo=args.demo)
    except NoInstrumentsError as exc:
        # An empty instrument set is a problem with the input, not a fault in the
        # skill, so it exits the way argparse exits on any other bad argument rather
        # than on a traceback the user has to interpret.
        parser.error(str(exc))


if __name__ == "__main__":
    main()
