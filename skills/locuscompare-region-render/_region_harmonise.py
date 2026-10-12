"""Odds ratio -> log-odds effect size, and the tag that says which one a row carries.

A case-control GWAS outcome often publishes an odds ratio and no beta. Any consumer that
needs a beta (a Wald ratio, an effect-size scatter) derives beta = ln(OR) and a standard
error here, so every consumer applies the same rule and gets the same number.

Rows are typed on the small `OddsRatioRow` Protocol rather than on any one skill's variant
class: anything that carries the seven fields below and is a dataclass (the conversion
returns `dataclasses.replace` copies) works. The module imports only the standard library,
so it can be copied unchanged into any skill that needs it; keep it that way.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, replace
from statistics import NormalDist
from typing import Literal, Protocol, TypeVar


class OddsRatioRow(Protocol):
    """The fields the conversion reads and writes. A dataclass instance is required at
    runtime, because filled rows are returned as `dataclasses.replace` copies."""

    variant_id: str
    beta: float | None
    se: float | None
    p_value: float | None
    odds_ratio: float | None
    ci_lower: float | None
    ci_upper: float | None


RowT = TypeVar("RowT", bound=OddsRatioRow)


# ---------------------------------------------------------------------------
# log(OR) -> beta derivation
# ---------------------------------------------------------------------------
# A binary/case-control GWAS outcome reports an ODDS RATIO, not a beta. The conversion
# lives on the consumer side, not in the region fetchers, whose job is a faithful fetch
# and parse (they pass the raw odds_ratio / CI columns through). Apply it to the whole
# fetched region, after the fetch and before any lead or proxy selection, so every later
# step sees the derived betas. OR-vs-beta is ASSUMED to be an
# outcome-wide (study-level) property, so that a study is all-native or all-OR-derived
# (unmeasured across studies). A consumer that can see both kinds in one window must not
# rely on it: the LocusCompare figure labels such a window BETA_SOURCE_MIXED.
_NORM = NormalDist()
_OR_BETA_CI95_Z = 1.959963984540054  # Phi^-1(0.975): the two-sided 95% normal quantile
_P_MAX_USABLE = 1.0 - 1e-12
_OR_Q_MAX = 1.0 - 1e-15  # largest 1-p/2 quantile inv_cdf resolves (z~8); caps tiny-p underflow
_OR_CI_P_DISAGREE_REL = 0.10  # flag CI-vs-p SE relative disagreement above this (annotate)

BETA_SOURCE_NATIVE = "native"
BETA_SOURCE_OR_DERIVED = "or_derived"
# A window in which some outcome rows carry a reported beta and others a beta derived
# from an odds ratio. Not a scale: the reported betas' scale is unknown, so no scale
# label is given for it (see `outcome_effect_scale_label`).
BETA_SOURCE_MIXED = "mixed"

# Human labels for the beta-source codes. The code stays in structured output for
# tooling; the label rides alongside it for a reader.
BETA_SOURCE_LABELS: dict[str, str] = {
    BETA_SOURCE_NATIVE: "effect size as reported by the study",
    BETA_SOURCE_OR_DERIVED: "log odds ratio, derived from the reported odds ratio",
    BETA_SOURCE_MIXED: (
        "mixed: some effect sizes as reported by the study, others log odds ratios "
        "derived from reported odds ratios"
    ),
}

# Why `_or_beta_se` returned None for a row that carries an odds ratio. Written once so
# every surface that states it (panel, notes, SKILL.md) says the same true thing: the SE
# needs a CI, or a usable p-value AND a non-zero log odds ratio (|ln OR| / z is 0 at
# OR = 1, which is not a standard error).
OR_NO_SE_REASON = (
    "no standard error could be derived (no confidence interval, and either no usable "
    "p-value or an odds ratio of exactly 1)"
)

# The axis / banner wording for an outcome whose effect sizes are log odds.
LOG_ODDS_DERIVED_LABEL = BETA_SOURCE_LABELS[BETA_SOURCE_OR_DERIVED]
LOG_ODDS_FINNGEN_LABEL = "log odds ratio, as reported by FinnGen"

# FinnGen publishes disease endpoints as log-odds betas with no odds-ratio column, so a
# FinnGen effect size is on the log-odds scale although its beta source is `native`.
# Measured 2026-10-08 on the FinnGen R12 manifest (finngen_R12_manifest.tsv): 2,466 of
# 2,469 endpoints carry both cases and controls; the 3 without controls are the
# inverse-rank-normalised quantitative traits BMI_IRN, HEIGHT_IRN and WEIGHT_IRN, which
# is what the suffix below excludes.
_FINNGEN_STUDY_PREFIX = "FINNGEN_"
_FINNGEN_QUANTITATIVE_SUFFIX = "_IRN"


def beta_source_label(beta_source: str | None) -> str | None:
    """Human label for a beta-source code; the raw code for an unknown one; None for None."""
    if beta_source is None:
        return None
    return BETA_SOURCE_LABELS.get(beta_source, beta_source)


def outcome_effect_scale_label(
    beta_source: str | None, outcome_study_id: str | None = None
) -> str | None:
    """The scale an outcome's effect sizes are on, when we can state it; else None.

    - `or_derived`: log odds ratio, derived from the reported odds ratio.
    - `native` on a FinnGen case-control endpoint: log odds ratio, as FinnGen reports it.
    - anything else: None. A native beta from any other study is on whatever scale the
      study used, and we do not know it, so we say nothing rather than guess.
    """
    if beta_source == BETA_SOURCE_OR_DERIVED:
        return LOG_ODDS_DERIVED_LABEL
    if beta_source == BETA_SOURCE_NATIVE and outcome_study_id:
        sid = outcome_study_id.upper()
        if sid.startswith(_FINNGEN_STUDY_PREFIX) and not sid.endswith(
            _FINNGEN_QUANTITATIVE_SUFFIX
        ):
            return LOG_ODDS_FINNGEN_LABEL
    return None


@dataclass(frozen=True)
class BetaSourceEvent:
    """Provenance for one outcome variant whose beta was derived from an odds ratio.

    Emitted so an OR-derived beta is never silently indistinguishable from a native one;
    keyed by the row's `variant_id` so each event maps back to its source row.
    `se_source` records which SE-precedence branch fired; `ci_p_flag` marks a CI-vs-p SE
    disagreement > 10% (which CI rounding can cause; annotate, never drop).
    """

    variant_id: str
    se_source: Literal["ci", "p"]
    ci_p_reldiff: float | None
    ci_p_flag: bool


def _or_beta_se(
    odds_ratio: float | None,
    ci_lower: float | None,
    ci_upper: float | None,
    p_value: float | None,
) -> tuple[float, float, str, float | None, bool] | None:
    """(beta, se, se_source, ci_p_reldiff, ci_p_flag) from an odds ratio, or None.

    beta = ln(OR) on the row's own effect allele; the caller's allele harmonisation flips
    it afterwards like any reported beta, and it is never flipped here. SE precedence:
    CI-derived (ln(hi)-ln(lo))/(2*z95) first, else |ln(OR)|/z with z=Phi^-1(1-p/2). A
    reported SE is not used, since its scale (log odds or not) cannot be confirmed.
    Returns None when the OR is unusable or no SE can be derived (no usable CI, and either
    no usable p or OR exactly 1, where |ln(OR)|/z is 0).
    """
    # Finite only: an infinite odds ratio (a separated case-control fit, or a parse of
    # "inf") would give an infinite beta and SE that a consumer then plots or divides by.
    if odds_ratio is None or not (odds_ratio > 0.0) or not math.isfinite(odds_ratio):
        return None
    beta = math.log(odds_ratio)

    se_ci: float | None = None
    if (
        ci_lower is not None and ci_upper is not None
        and ci_lower > 0.0 and ci_upper > 0.0 and ci_upper > ci_lower
        and math.isfinite(ci_upper)
    ):
        cand = (math.log(ci_upper) - math.log(ci_lower)) / (2.0 * _OR_BETA_CI95_Z)
        if cand > 0.0:
            se_ci = cand

    se_p: float | None = None
    if p_value is not None and 0.0 < p_value < _P_MAX_USABLE:
        # A genome-wide-significant binary hit has p far below float resolution, so
        # 1 - p/2 rounds to exactly 1.0 and NormalDist.inv_cdf(1.0) raises. Cap the
        # quantile at the largest value strictly < 1 that inv_cdf resolves (z ~ 8),
        # the honest float64 ceiling on |z|. It makes the SE of such a hit an upper
        # bound, a conservative error.
        q = min(1.0 - p_value / 2.0, _OR_Q_MAX)
        z = _NORM.inv_cdf(q)  # |z| from two-sided p
        if z > 0.0:
            cand = abs(beta) / z
            if cand > 0.0:
                se_p = cand

    if se_ci is not None:
        se, se_source = se_ci, "ci"
    elif se_p is not None:
        se, se_source = se_p, "p"
    else:
        return None

    ci_p_reldiff: float | None = None
    ci_p_flag = False
    if se_ci is not None and se_p is not None:
        # Relative to the kept (primary) CI-derived SE; since z = |beta|/SE, this equals
        # the relative difference of the two implied z-statistics.
        ci_p_reldiff = abs(se_ci - se_p) / se_ci
        ci_p_flag = ci_p_reldiff > _OR_CI_P_DISAGREE_REL
    return beta, se, se_source, ci_p_reldiff, ci_p_flag


def derive_beta_from_or(
    outcome_rows: list[RowT],
) -> tuple[list[RowT], list[BetaSourceEvent]]:
    """Fill beta=ln(OR)+SE on OUTCOME rows that report an odds ratio but no beta.

    Region-level, idempotent (only fills a row whose `beta is None` and which carries an
    `odds_ratio`), applied before any lead or proxy selection so later steps see the
    derived betas. Returns (new_rows, events); `events` is non-empty iff any row was
    converted. Rows are otherwise returned unchanged (native-beta studies produce zero
    events).
    """
    out: list[RowT] = []
    events: list[BetaSourceEvent] = []
    for r in outcome_rows:
        if r.beta is None and r.odds_ratio is not None:
            conv = _or_beta_se(r.odds_ratio, r.ci_lower, r.ci_upper, r.p_value)
            if conv is not None:
                beta, se, se_source, reldiff, flag = conv
                out.append(replace(r, beta=beta, se=se))
                events.append(BetaSourceEvent(
                    variant_id=r.variant_id, se_source=se_source,
                    ci_p_reldiff=reldiff, ci_p_flag=flag,
                ))
                continue
        out.append(r)
    return out, events


__all__ = [
    "BETA_SOURCE_LABELS",
    "BETA_SOURCE_MIXED",
    "BETA_SOURCE_NATIVE",
    "BETA_SOURCE_OR_DERIVED",
    "BetaSourceEvent",
    "LOG_ODDS_DERIVED_LABEL",
    "LOG_ODDS_FINNGEN_LABEL",
    "OR_NO_SE_REASON",
    "OddsRatioRow",
    "beta_source_label",
    "derive_beta_from_or",
    "outcome_effect_scale_label",
]
