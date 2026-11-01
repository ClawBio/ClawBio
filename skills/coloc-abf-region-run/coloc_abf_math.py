#!/usr/bin/env python3
"""coloc.abf: the single-causal-variant colocalisation computation, in Python.

A clean-room implementation of the coloc.abf model. It imports nothing from the R
`coloc` package and does not shell out to R; the arithmetic is a Wakefield
approximate Bayes factor per variant plus the Giambartolomei posterior combine,
both small enough to write directly from the papers.

Method
------
**Per-variant approximate Bayes factor (Wakefield 2009).** For each variant, with
`z = beta / se`, `V = se^2`, and prior effect variance `W`:

    r    = W / (W + V)
    lABF = 0.5 * ( log(1 - r) + r * z^2 )        # natural log

`W` is `(0.15 * sdY)^2` on the trait-SD-standardised scale for a quantitative
trait and `0.2^2` on the log-odds scale for a case-control trait. A quantitative
trait needs `sdY`; when it is not supplied it is estimated from MAF, se and N by
regressing `1 / Var(beta)` on `2 N MAF (1 - MAF)` through the origin.

**Posterior combine (Giambartolomei 2014).** Under the assumption of at most one
causal variant per trait in the region, the five hypotheses are combined from the
two traits' per-variant lABFs (`l1`, `l2`, aligned to the SAME variants) with
per-variant priors `p1`, `p2`, `p12`:

    lH0 = 0
    lH1 = log(p1)           + logsum(l1)
    lH2 = log(p2)           + logsum(l2)
    lH3 = log(p1) + log(p2) + logdiff( logsum(l1) + logsum(l2), logsum(l1 + l2) )
    lH4 = log(p12)          + logsum(l1 + l2)
    PP  = softmax([lH0, lH1, lH2, lH3, lH4])

`logsum` is log-sum-exp and `logdiff(a, b) = log(exp(a) - exp(b))`. H3 is the
off-diagonal mass (each trait causal at a different variant), H4 the diagonal
(one variant causal for both).

**Sensitivity (Wallace 2020).** A single PP.H4 is sensitive to `p12`, the prior
that one variant is causal for both traits. The posteriors are therefore
recomputed across a range of `p12` and the conclusion is reported as
prior-driven when it does not hold across that range. Because the per-variant
Bayes factors do not depend on the priors, the sweep reuses them and only repeats
the cheap combine.

This module is pure: no file I/O, no network, no plotting, and no threshold
beyond the sensitivity rule it is handed.

Citations
---------
- Giambartolomei C, Vukcevic D, Schadt EE, et al. Bayesian test for
  colocalisation between pairs of genetic association studies using summary
  statistics. PLoS Genet. 2014;10(5):e1004383. doi:10.1371/journal.pgen.1004383.
  PMID 24830394.
- Wakefield J. Bayes factors for genome-wide association studies: comparison
  with P-values. Genet Epidemiol. 2009;33(1):79-86. doi:10.1002/gepi.20359.
  PMID 18642345.
- Wallace C. Eliciting priors and relaxing the single causal variant assumption
  in colocalisation analyses. PLoS Genet. 2020;16(4):e1008720.
  doi:10.1371/journal.pgen.1008720. PMID 32310995.
"""

from __future__ import annotations

import math
from collections import Counter
from dataclasses import dataclass, field
from typing import Any, Literal, Sequence

import numpy as np

TraitType = Literal["quant", "cc"]

# Per-variant priors. p12 is the one worth arguing about, and the one the
# sensitivity sweep is built around.
DEFAULT_P1 = 1e-4
DEFAULT_P2 = 1e-4
DEFAULT_P12 = 1e-5

# Prior effect-size standard deviations.
QUANT_PRIOR_SD = 0.15  # standardised scale; W = (0.15 * sdY)^2
CC_PRIOR_SD = 0.2      # log-odds scale; W = 0.2^2 = 0.04

# The rule the sensitivity sweep is evaluated against, and the p12 range it is
# swept over: one order of magnitude either side of the default p12.
DEFAULT_RULE_H4_MIN = 0.5
DEFAULT_P12_GRID_LOW = 1e-6
DEFAULT_P12_GRID_HIGH = 1e-4
DEFAULT_P12_GRID_POINTS = 9


class ColocError(ValueError):
    """Raised when two datasets cannot support a coloc.abf computation."""


@dataclass(frozen=True)
class ColocDataset:
    """One trait's regional summary statistics.

    `snp` holds one identifier per variant, the same identifiers the partner
    dataset uses for the same variants. `beta` and `se` are the per-variant
    effect estimate and its standard error.

    With `trait_type="quant"`, supply `sdY`, or supply `maf` and `n` so it can be
    estimated. With `trait_type="cc"`, the prior is fixed on the log-odds scale
    and `sdY`, `maf` and `n` are not read.
    """

    snp: Sequence[str]
    beta: np.ndarray
    se: np.ndarray
    trait_type: TraitType
    label: str = ""
    sdY: float | None = None
    maf: np.ndarray | None = None
    n: float | np.ndarray | None = None

    def __post_init__(self) -> None:
        beta = np.asarray(self.beta, dtype=float)
        se = np.asarray(self.se, dtype=float)
        object.__setattr__(self, "beta", beta)
        object.__setattr__(self, "se", se)
        if self.maf is not None:
            object.__setattr__(self, "maf", np.asarray(self.maf, dtype=float))
        n_snp = len(self.snp)
        if not (len(beta) == len(se) == n_snp):
            raise ColocError(
                f"dataset {self.label!r}: snp/beta/se length mismatch "
                f"({n_snp}/{len(beta)}/{len(se)})"
            )
        # One row per variant. The combine sums one Bayes factor per variant, so
        # a repeated identifier would be counted twice.
        dupes = [s for s, k in Counter(self.snp).items() if k > 1]
        if dupes:
            raise ColocError(
                f"dataset {self.label!r}: {len(dupes)} variant identifier(s) appear "
                f"more than once (e.g. {dupes[:3]}). coloc.abf sums one Bayes factor "
                f"per variant, so a repeat would be counted twice. Three usual causes: "
                f"the source wrote the same variant twice; one rsID names both alleles "
                f"of a multi-allelic site, so identify rows by chromosome, position and "
                f"alleles instead; or the rows span more than one molecular trait, such "
                f"as several transcripts of one gene, in which case pick one."
            )
        if self.trait_type not in ("quant", "cc"):
            raise ColocError(
                f"dataset {self.label!r}: trait_type must be 'quant' or 'cc', "
                f"got {self.trait_type!r}"
            )


@dataclass(frozen=True)
class ColocAbfResult:
    """The five posteriors plus the provenance of how they were computed.

    `pp` carries H0 to H4 by name. `n_snps` is the number of SHARED variants the
    posteriors were computed over; coloc.abf is defined on the intersection
    only. `best_snp` carries the largest share of the H4 mass, which makes it the
    candidate shared causal variant. `n_colocalising_variants` is the size of the
    smallest set of variants holding 95% of that mass: when it is very small the
    result rests on few variants.
    """

    h0: float
    h1: float
    h2: float
    h3: float
    h4: float
    n_snps: int
    n_colocalising_variants: int
    best_snp: str | None
    best_snp_pp_h4: float | None
    priors: dict[str, float]
    prior_sd: dict[str, float]
    sdy_used: dict[str, float | None]
    sdy_estimated: dict[str, bool]
    notes: list[str] = field(default_factory=list)
    provenance: dict[str, Any] = field(default_factory=dict)

    @property
    def pp(self) -> dict[str, float]:
        return {
            "PP.H0": self.h0,
            "PP.H1": self.h1,
            "PP.H2": self.h2,
            "PP.H3": self.h3,
            "PP.H4": self.h4,
        }


@dataclass(frozen=True)
class AlignedLogAbf:
    """Two traits' per-variant lABFs over their shared variants.

    Everything in coloc.abf that does NOT depend on p1, p2 or p12. Computed once
    and reused by the posterior combine and by the whole sensitivity sweep.
    """

    shared_snps: list[str]
    l1: np.ndarray
    l2: np.ndarray
    sdy_used: dict[str, float | None]
    sdy_estimated: dict[str, bool]
    n_full1: int
    n_full2: int


@dataclass(frozen=True)
class SensitivityPoint:
    """One p12 on the grid: its posteriors, and whether the rule holds there."""

    p12: float
    h3: float
    h4: float
    rule_passes: bool


@dataclass(frozen=True)
class SensitivityResult:
    """Whether the H4 conclusion survives the swept range of p12.

    `prior_driven` is True when the rule's verdict CHANGES across the range, so
    the conclusion depends on the prior rather than on the data.
    `passing_p12_range` is the sub-interval of p12 over which the rule holds, or
    None when it never holds.
    """

    rule: str
    rule_h4_min: float
    p12_default: float
    h4_at_default: float
    passes_at_default: bool
    prior_driven: bool
    robust: bool
    passing_p12_range: tuple[float, float] | None
    points: list[SensitivityPoint] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)

    @property
    def verdict(self) -> str:
        """A short disposition for a log line or a report heading."""
        if not self.passes_at_default:
            return "no shared signal (H4 below the rule at the default p12)"
        if self.robust:
            lo = self.points[0].p12
            hi = self.points[-1].p12
            return f"robust (the H4 call holds across p12 {lo:.0e} to {hi:.0e})"
        return "prior-driven (the H4 call flips inside the swept p12 range)"


# ---------------------------------------------------------------------------
# Per-variant Bayes factors
# ---------------------------------------------------------------------------

def per_snp_log_abf(beta: np.ndarray, se: np.ndarray, w: float) -> np.ndarray:
    """Per-variant natural-log Wakefield approximate Bayes factor.

    `lABF = 0.5 * ( log(1 - r) + r * z^2 )`, with `r = W / (W + V)`,
    `z = beta / se` and `V = se^2`. Vectorised over variants. Standard errors
    must be finite and positive: a zero se gives `V = 0`, hence `r = 1` and
    `log(0)`.

    The expression reads beta only as `z^2`, so the result does not depend on
    which allele each trait reported its effect for.
    """
    beta = np.asarray(beta, dtype=float)
    se = np.asarray(se, dtype=float)
    if se.size and (np.any(se <= 0) or not np.all(np.isfinite(se))):
        raise ColocError("all standard errors must be finite and positive")
    if beta.size and not np.all(np.isfinite(beta)):
        raise ColocError("all effect estimates must be finite")
    v = se ** 2
    r = w / (w + v)
    z = beta / se
    return 0.5 * (np.log1p(-r) + r * z ** 2)


def estimate_sdy(se: np.ndarray, maf: np.ndarray, n: float | np.ndarray) -> float:
    """Estimate a quantitative trait's standard deviation from MAF, se and N.

    For a quantitative trait, `Var(beta_j) ~= sdY^2 / (N_j Var(X_j))` with
    `Var(X_j) = 2 MAF_j (1 - MAF_j)`. So `2 N_j MAF_j (1 - MAF_j)` is linear in
    `1 / Var(beta_j)` through the origin with slope `sdY^2`, and the estimate is
    the square root of the ordinary-least-squares slope through the origin,
    `sum(x y) / sum(x^2)`.
    """
    se = np.asarray(se, dtype=float)
    maf = np.asarray(maf, dtype=float)
    n_arr = np.asarray(n, dtype=float)
    if not se.size:
        raise ColocError("cannot estimate sdY from an empty dataset")
    if np.any(se <= 0) or not np.all(np.isfinite(se)):
        raise ColocError("estimating sdY needs finite, positive standard errors")
    if not np.all(np.isfinite(maf)) or np.any((maf <= 0) | (maf >= 1)):
        raise ColocError("estimating sdY needs MAF strictly between 0 and 1")
    if not np.all(np.isfinite(n_arr)) or np.any(n_arr <= 0):
        raise ColocError("estimating sdY needs a positive sample size")
    oneover = 1.0 / se ** 2
    nvx = 2.0 * n_arr * maf * (1.0 - maf)
    slope = float(np.sum(oneover * nvx) / np.sum(oneover * oneover))
    if slope <= 0 or not math.isfinite(slope):
        raise ColocError(
            "estimating sdY produced a non-positive slope; supply sdY explicitly"
        )
    return math.sqrt(slope)


def prior_variance(trait_type: TraitType, sdy: float) -> float:
    """The Wakefield prior effect variance `W`.

    `(0.15 * sdY)^2` for a quantitative trait; `0.2^2` for a case-control trait,
    which ignores sdY.
    """
    if trait_type == "cc":
        return CC_PRIOR_SD ** 2
    return (QUANT_PRIOR_SD * sdy) ** 2


def _logsumexp(x: np.ndarray) -> float:
    x = np.asarray(x, dtype=float)
    m = float(np.max(x))
    if not math.isfinite(m):
        return m
    return m + math.log(float(np.sum(np.exp(x - m))))


def _logdiff(a: float, b: float) -> float:
    """`log(exp(a) - exp(b))`, which needs `a >= b`."""
    if a < b:
        raise ColocError(f"logdiff requires a >= b, got a={a}, b={b}")
    if a == b:
        return -math.inf
    return a + math.log1p(-math.exp(b - a))


def combine_abf(
    l1: np.ndarray,
    l2: np.ndarray,
    *,
    p1: float = DEFAULT_P1,
    p2: float = DEFAULT_P2,
    p12: float = DEFAULT_P12,
) -> dict[str, float]:
    """Combine two traits' per-variant lABFs into PP.H0 to PP.H4.

    `l1` and `l2` are the per-variant natural-log Bayes factors of the two
    traits, aligned to the same variants in the same order. Returns the five
    posteriors keyed H0 to H4.
    """
    l1 = np.asarray(l1, dtype=float)
    l2 = np.asarray(l2, dtype=float)
    if l1.shape != l2.shape:
        raise ColocError(f"l1 and l2 must be aligned; got {l1.shape} vs {l2.shape}")
    if l1.size == 0:
        raise ColocError("no shared variants to colocalise")
    for name, p in (("p1", p1), ("p2", p2), ("p12", p12)):
        if not (0.0 < p < 1.0) or not math.isfinite(p):
            raise ColocError(f"{name} must be strictly between 0 and 1, got {p}")

    lsum = l1 + l2
    ls1 = _logsumexp(l1)
    ls2 = _logsumexp(l2)
    lsboth = _logsumexp(lsum)

    lh0 = 0.0
    lh1 = math.log(p1) + ls1
    lh2 = math.log(p2) + ls2
    lh3 = math.log(p1) + math.log(p2) + _logdiff(ls1 + ls2, lsboth)
    lh4 = math.log(p12) + lsboth

    all_l = np.array([lh0, lh1, lh2, lh3, lh4], dtype=float)
    pp = np.exp(all_l - _logsumexp(all_l))
    return {
        "H0": float(pp[0]),
        "H1": float(pp[1]),
        "H2": float(pp[2]),
        "H3": float(pp[3]),
        "H4": float(pp[4]),
    }


def _resolve_w(d: ColocDataset) -> tuple[float, float | None, bool]:
    """Resolve `(W, sdY_used, sdY_was_estimated)` for one dataset."""
    if d.trait_type == "cc":
        return prior_variance("cc", 1.0), None, False
    if d.sdY is not None:
        if not math.isfinite(d.sdY) or d.sdY <= 0:
            raise ColocError(f"dataset {d.label!r}: sdY must be positive")
        return prior_variance("quant", d.sdY), d.sdY, False
    if d.maf is not None and d.n is not None:
        sdy = estimate_sdy(d.se, d.maf, d.n)
        return prior_variance("quant", sdy), sdy, True
    raise ColocError(
        f"dataset {d.label!r}: a quantitative trait needs sdY, or maf and n so "
        "that sdY can be estimated"
    )


def aligned_log_abf(dataset1: ColocDataset, dataset2: ColocDataset) -> AlignedLogAbf:
    """Intersect two datasets by variant and compute each trait's per-variant lABF.

    coloc.abf is defined on variants measured in both traits, so this takes the
    intersection, preserving the order of `dataset2`, resolves each trait's prior
    variance, and returns the two aligned lABF vectors.
    """
    idx1 = {s: i for i, s in enumerate(dataset1.snp)}
    shared = [s for s in dataset2.snp if s in idx1]
    if not shared:
        raise ColocError(
            f"no shared variants between {dataset1.label!r} and {dataset2.label!r}. "
            "coloc.abf needs both traits measured at the same variants, so give both "
            "files identifiers on one genome build."
        )
    take1 = np.array([idx1[s] for s in shared])
    idx2 = {s: i for i, s in enumerate(dataset2.snp)}
    take2 = np.array([idx2[s] for s in shared])

    w1, sdy1, est1 = _resolve_w(dataset1)
    w2, sdy2, est2 = _resolve_w(dataset2)

    l1 = per_snp_log_abf(dataset1.beta[take1], dataset1.se[take1], w1)
    l2 = per_snp_log_abf(dataset2.beta[take2], dataset2.se[take2], w2)
    return AlignedLogAbf(
        shared_snps=shared,
        l1=l1,
        l2=l2,
        sdy_used={"trait1": sdy1, "trait2": sdy2},
        sdy_estimated={"trait1": est1, "trait2": est2},
        n_full1=len(dataset1.snp),
        n_full2=len(dataset2.snp),
    )


def run_coloc_abf(
    dataset1: ColocDataset,
    dataset2: ColocDataset,
    *,
    p1: float = DEFAULT_P1,
    p2: float = DEFAULT_P2,
    p12: float = DEFAULT_P12,
    aligned: AlignedLogAbf | None = None,
) -> ColocAbfResult:
    """Run coloc.abf on two datasets of regional summary statistics.

    Intersects them by variant, computes each trait's per-variant lABF with its
    resolved prior variance, combines into PP.H0 to PP.H4, and reports the
    candidate shared causal variant. Pass `aligned` to reuse per-variant Bayes
    factors already computed, as the sensitivity sweep does.
    """
    if aligned is None:
        aligned = aligned_log_abf(dataset1, dataset2)
    shared, l1, l2 = aligned.shared_snps, aligned.l1, aligned.l2

    pp = combine_abf(l1, l2, p1=p1, p2=p2, p12=p12)

    # How the H4 mass is spread over variants.
    lsum = l1 + l2
    snp_h4 = np.exp(lsum - _logsumexp(lsum))
    best = int(np.argmax(snp_h4))
    n_coloc_variants = int(np.searchsorted(np.cumsum(np.sort(snp_h4)[::-1]), 0.95) + 1)
    n_coloc_variants = min(n_coloc_variants, len(shared))

    notes: list[str] = []
    n_full1, n_full2 = aligned.n_full1, aligned.n_full2
    if len(shared) < max(n_full1, n_full2):
        notes.append(
            f"computed over {len(shared)} shared variants; the first trait carried "
            f"{n_full1} and the second {n_full2}. Variants present in only one trait "
            "are not used, because coloc.abf needs both traits at the same variant."
        )

    return ColocAbfResult(
        h0=pp["H0"], h1=pp["H1"], h2=pp["H2"], h3=pp["H3"], h4=pp["H4"],
        n_snps=len(shared),
        n_colocalising_variants=n_coloc_variants,
        best_snp=shared[best],
        best_snp_pp_h4=float(snp_h4[best]),
        priors={"p1": p1, "p2": p2, "p12": p12},
        prior_sd={
            "trait1": QUANT_PRIOR_SD if dataset1.trait_type == "quant" else CC_PRIOR_SD,
            "trait2": QUANT_PRIOR_SD if dataset2.trait_type == "quant" else CC_PRIOR_SD,
        },
        sdy_used=dict(aligned.sdy_used),
        sdy_estimated=dict(aligned.sdy_estimated),
        notes=notes,
        provenance={
            "method": "coloc.abf (Wakefield approximate Bayes factor, Giambartolomei combine)",
            "trait1_type": dataset1.trait_type,
            "trait2_type": dataset2.trait_type,
        },
    )


def coloc_sensitivity(
    dataset1: ColocDataset,
    dataset2: ColocDataset,
    *,
    p1: float = DEFAULT_P1,
    p2: float = DEFAULT_P2,
    p12_default: float = DEFAULT_P12,
    rule_h4_min: float = DEFAULT_RULE_H4_MIN,
    p12_grid_low: float = DEFAULT_P12_GRID_LOW,
    p12_grid_high: float = DEFAULT_P12_GRID_HIGH,
    n_points: int = DEFAULT_P12_GRID_POINTS,
    aligned: AlignedLogAbf | None = None,
) -> SensitivityResult:
    """Sweep p12 and report whether the H4 conclusion survives the range.

    The per-variant Bayes factors do not depend on the priors, so the sweep
    reuses them and repeats only the combine. The default grid brackets the
    default p12 by an order of magnitude either way.
    """
    if aligned is None:
        aligned = aligned_log_abf(dataset1, dataset2)
    l1, l2 = aligned.l1, aligned.l2

    grid = np.logspace(np.log10(p12_grid_low), np.log10(p12_grid_high), n_points)
    # Evaluate the default p12 exactly, so the headline number is one of the
    # swept points rather than an interpolation of them.
    grid = np.unique(np.concatenate([grid, [p12_default]]))

    points: list[SensitivityPoint] = []
    for p12 in grid:
        pp = combine_abf(l1, l2, p1=p1, p2=p2, p12=float(p12))
        points.append(
            SensitivityPoint(
                p12=float(p12),
                h3=pp["H3"],
                h4=pp["H4"],
                rule_passes=pp["H4"] >= rule_h4_min,
            )
        )

    default_pt = min(points, key=lambda pt: abs(pt.p12 - p12_default))
    passes_at_default = default_pt.rule_passes

    passing = [pt.p12 for pt in points if pt.rule_passes]
    passing_range = (min(passing), max(passing)) if passing else None

    prior_driven = len({pt.rule_passes for pt in points}) > 1
    robust = passes_at_default and not prior_driven

    notes: list[str] = []
    if prior_driven and passes_at_default:
        notes.append(
            f"The H4 call is prior-driven: the rule (PP.H4 >= {rule_h4_min}) does not "
            f"hold across the whole swept p12 range ({p12_grid_low:g} to "
            f"{p12_grid_high:g}), so treat the colocalisation as suggestive rather "
            "than established."
        )
    elif not passes_at_default:
        notes.append(
            f"No colocalisation at the default p12={p12_default:g}: "
            f"PP.H4={default_pt.h4:.4f} is below the rule {rule_h4_min}."
        )

    return SensitivityResult(
        rule=f"PP.H4 >= {rule_h4_min}",
        rule_h4_min=rule_h4_min,
        p12_default=p12_default,
        h4_at_default=default_pt.h4,
        passes_at_default=passes_at_default,
        prior_driven=prior_driven,
        robust=robust,
        passing_p12_range=passing_range,
        points=points,
        notes=notes,
    )


__all__ = [
    "AlignedLogAbf",
    "CC_PRIOR_SD",
    "ColocAbfResult",
    "ColocDataset",
    "ColocError",
    "DEFAULT_P1",
    "DEFAULT_P2",
    "DEFAULT_P12",
    "DEFAULT_P12_GRID_HIGH",
    "DEFAULT_P12_GRID_LOW",
    "DEFAULT_P12_GRID_POINTS",
    "DEFAULT_RULE_H4_MIN",
    "QUANT_PRIOR_SD",
    "SensitivityPoint",
    "SensitivityResult",
    "TraitType",
    "aligned_log_abf",
    "coloc_sensitivity",
    "combine_abf",
    "estimate_sdy",
    "per_snp_log_abf",
    "prior_variance",
    "run_coloc_abf",
]
