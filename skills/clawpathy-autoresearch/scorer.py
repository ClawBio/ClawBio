"""Reproduction error scorer for clawpathy-autoresearch.

Score = mean absolute error across concrete numerical targets.
Lower is better. Zero = perfect reproduction.

Metrics:
- p_value_error: normalised |target - reproduced| / target for -log10(p)
- or_error: absolute |target - reproduced| for odds ratios
- freq_error: absolute |target - reproduced| for effect allele frequencies
- locus_count_error: normalised |target - reproduced| / target for total loci
- variant_missing_penalty: 1.0 per missing variant (not found by rsID)
- direction_error: 1.0 per incorrect effect direction call
"""
from __future__ import annotations

import importlib.util
from dataclasses import dataclass
from pathlib import Path
from typing import Any


# Penalty added per missing variant (not found in agent output)
MISSING_VARIANT_PENALTY = 1.0

# Penalty per wrong direction call
DIRECTION_PENALTY = 1.0

# Weights for combining error components into total
WEIGHTS = {
    "p_value": 0.20,
    "or": 0.25,
    "freq": 0.10,
    "locus_count": 0.15,
    "missing": 0.20,
    "direction": 0.10,
}


@dataclass
class ErrorBreakdown:
    """Breakdown of reproduction error components. All >= 0, lower is better."""

    p_value_error: float
    or_error: float
    freq_error: float
    locus_count_error: float
    variant_missing_penalty: float
    direction_error: float
    total: float

    def to_dict(self) -> dict[str, float]:
        return {
            "p_value_error": self.p_value_error,
            "or_error": self.or_error,
            "freq_error": self.freq_error,
            "locus_count_error": self.locus_count_error,
            "variant_missing_penalty": self.variant_missing_penalty,
            "direction_error": self.direction_error,
            "total": self.total,
        }


def reproduction_error(
    ground_truth: dict[str, Any],
    agent_output: dict[str, Any],
) -> ErrorBreakdown:
    """Compute reproduction error between ground truth and agent output.

    All error components are >= 0. Total is a weighted sum.
    Perfect reproduction returns total == 0.
    """
    gt_variants = ground_truth.get("lead_variants", [])
    agent_variants = agent_output.get("variants_found", [])
    gt_loci = ground_truth.get("total_loci", 0)
    agent_loci = agent_output.get("total_loci_reported", 0)

    if not gt_variants:
        lce = abs(gt_loci - agent_loci) / max(gt_loci, 1)
        return ErrorBreakdown(0.0, 0.0, 0.0, lce, 0.0, 0.0, lce * WEIGHTS["locus_count"])

    agent_by_rsid = {v["rsid"]: v for v in agent_variants}

    p_errors = []
    or_errors = []
    freq_errors = []
    direction_errors = []
    missing_count = 0

    for gt_v in gt_variants:
        rsid = gt_v["rsid"]
        if rsid not in agent_by_rsid:
            missing_count += 1
            continue

        av = agent_by_rsid[rsid]

        # P-value error: normalised by target
        gt_p = gt_v.get("neg_log10_p", 0)
        ag_p = av.get("neg_log10_p", 0)
        if gt_p > 0:
            p_errors.append(abs(gt_p - ag_p) / gt_p)
        else:
            p_errors.append(abs(gt_p - ag_p))

        # OR error: absolute difference
        gt_or = gt_v.get("odds_ratio", 1.0)
        ag_or = av.get("odds_ratio", 1.0)
        or_errors.append(abs(gt_or - ag_or))

        # Frequency error: absolute difference
        gt_freq = gt_v.get("effect_allele_freq", 0.0)
        ag_freq = av.get("effect_allele_freq", 0.0)
        freq_errors.append(abs(gt_freq - ag_freq))

        # Direction error: binary
        gt_dir = gt_v.get("effect_direction", "")
        ag_dir = av.get("effect_direction", "")
        if gt_dir and ag_dir and gt_dir != ag_dir:
            direction_errors.append(DIRECTION_PENALTY)
        else:
            direction_errors.append(0.0)

    n = len(gt_variants)

    p_value_error = sum(p_errors) / len(p_errors) if p_errors else 0.0
    or_error = sum(or_errors) / len(or_errors) if or_errors else 0.0
    freq_error = sum(freq_errors) / len(freq_errors) if freq_errors else 0.0
    direction_error = sum(direction_errors) / n
    variant_missing_penalty = missing_count * MISSING_VARIANT_PENALTY / n

    # Locus count error: normalised
    if gt_loci > 0:
        locus_count_error = abs(gt_loci - agent_loci) / gt_loci
    else:
        locus_count_error = float(agent_loci)

    # Weighted total
    total = (
        WEIGHTS["p_value"] * p_value_error
        + WEIGHTS["or"] * or_error
        + WEIGHTS["freq"] * freq_error
        + WEIGHTS["locus_count"] * locus_count_error
        + WEIGHTS["missing"] * variant_missing_penalty
        + WEIGHTS["direction"] * direction_error
    )

    return ErrorBreakdown(
        p_value_error=p_value_error,
        or_error=or_error,
        freq_error=freq_error,
        locus_count_error=locus_count_error,
        variant_missing_penalty=variant_missing_penalty,
        direction_error=direction_error,
        total=total,
    )


def load_scorer(scorer_path: Path) -> callable:
    """Dynamically import a scorer.py and return its score() function.

    The scorer module must define: score(skill_output: dict, ground_truth: dict) -> float
    """
    scorer_path = Path(scorer_path)
    if not scorer_path.exists():
        raise FileNotFoundError(f"Scorer not found: {scorer_path}")

    spec = importlib.util.spec_from_file_location("workspace_scorer", scorer_path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)

    if not hasattr(module, "score"):
        raise AttributeError(
            f"Scorer at {scorer_path} must define a score() function"
        )
    return module.score
