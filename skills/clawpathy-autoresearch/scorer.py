"""Hybrid scoring engine for clawpathy-autoresearch.

Combines automated numerical checks (variant recovery, effect sizes, locus counts)
with LLM-as-judge qualitative scoring.
"""
from __future__ import annotations

import math
from dataclasses import dataclass
from typing import Any


@dataclass
class ScoreResult:
    """Breakdown of automated scoring components."""

    variant_recovery: float  # 0-1: fraction of ground truth variants found
    direction_accuracy: float  # 0-1: correct risk/protective calls
    effect_size_accuracy: float  # 0-1: OR within expected range
    locus_count_accuracy: float  # 0-1: reported loci count vs expected
    automated_total: float  # 0-10: weighted combination

    @property
    def breakdown(self) -> dict[str, float]:
        return {
            "variant_recovery": self.variant_recovery,
            "direction_accuracy": self.direction_accuracy,
            "effect_size_accuracy": self.effect_size_accuracy,
            "locus_count_accuracy": self.locus_count_accuracy,
            "automated_total": self.automated_total,
        }


def automated_score(
    ground_truth: dict[str, Any],
    agent_output: dict[str, Any],
) -> ScoreResult:
    """Score agent output against ground truth using automated checks.

    Components (each 0-1, then combined to 0-10):
    - variant_recovery: fraction of expected lead variants found by rsID
    - direction_accuracy: fraction with correct effect direction
    - effect_size_accuracy: fraction with OR inside expected range
    - locus_count_accuracy: 1 - |expected - reported| / expected (floored at 0)
    """
    gt_variants = ground_truth.get("lead_variants", [])
    agent_variants = agent_output.get("variants_found", [])
    gt_loci = ground_truth.get("total_loci", 0)
    agent_loci = agent_output.get("total_loci_reported", 0)

    if not gt_variants:
        return ScoreResult(0.0, 0.0, 0.0, 0.0, 0.0)

    # Build lookup of agent variants by rsID
    agent_by_rsid = {v["rsid"]: v for v in agent_variants}

    found = 0
    direction_correct = 0
    effect_correct = 0

    for gt_v in gt_variants:
        rsid = gt_v["rsid"]
        if rsid in agent_by_rsid:
            found += 1
            av = agent_by_rsid[rsid]

            # Direction check
            if av.get("effect_direction") == gt_v.get("effect_direction"):
                direction_correct += 1

            # Effect size check
            or_range = gt_v.get("or_range", [0, 0])
            agent_or = av.get("odds_ratio", 0)
            if or_range[0] <= agent_or <= or_range[1]:
                effect_correct += 1

    n = len(gt_variants)
    variant_recovery = found / n
    direction_accuracy = direction_correct / n if found > 0 else 0.0
    effect_size_accuracy = effect_correct / n if found > 0 else 0.0

    # Locus count accuracy
    if gt_loci > 0:
        locus_count_accuracy = max(0.0, 1.0 - abs(gt_loci - agent_loci) / gt_loci)
    else:
        locus_count_accuracy = 1.0 if agent_loci == 0 else 0.0

    # Weighted combination to 0-10 scale
    # Weights: variant recovery 35%, direction 25%, effect size 25%, locus count 15%
    automated_total = (
        0.35 * variant_recovery
        + 0.25 * direction_accuracy
        + 0.25 * effect_size_accuracy
        + 0.15 * locus_count_accuracy
    ) * 10.0

    return ScoreResult(
        variant_recovery=variant_recovery,
        direction_accuracy=direction_accuracy,
        effect_size_accuracy=effect_size_accuracy,
        locus_count_accuracy=locus_count_accuracy,
        automated_total=automated_total,
    )


def llm_judge_score(
    paper_results_text: str,
    agent_output_text: str,
    model: str = "claude-sonnet-4-20250514",
) -> tuple[float, str]:
    """Score agent output using an LLM judge.

    Returns (score 0-10, reasoning text).
    Falls back to 5.0 if the API call fails.
    """
    try:
        import anthropic
    except ImportError:
        return 5.0, "LLM judge unavailable (anthropic SDK not installed)"

    rubric = (
        "Score the agent's reproduction of this GWAS paper's findings on a 0-10 scale.\n\n"
        "Rubric:\n"
        "- Correct lead variants identified (0-2)\n"
        "- Effect directions and sizes accurate (0-2)\n"
        "- Biological pathways / enrichments identified (0-2)\n"
        "- Methodology appropriate for the paper (0-2)\n"
        "- Overall coherence and completeness (0-2)\n\n"
        "Respond with ONLY a JSON object: {\"score\": <number>, \"reasoning\": \"<text>\"}"
    )

    prompt = (
        f"## Paper's actual results\n\n{paper_results_text}\n\n"
        f"## Agent's reproduction output\n\n{agent_output_text}\n\n"
        f"## Scoring rubric\n\n{rubric}"
    )

    try:
        client = anthropic.Anthropic()
        response = client.messages.create(
            model=model,
            max_tokens=500,
            messages=[{"role": "user", "content": prompt}],
        )
        import json

        result = json.loads(response.content[0].text)
        return float(result["score"]), result["reasoning"]
    except Exception as e:
        return 5.0, f"LLM judge error: {e}"


def combine_scores(
    automated: float,
    llm: float,
    auto_weight: float = 0.6,
    llm_weight: float = 0.4,
) -> float:
    """Combine automated and LLM scores into a final 0-10 score."""
    combined = auto_weight * automated + llm_weight * llm
    return max(0.0, min(10.0, combined))
