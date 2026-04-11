"""Tests for hybrid scorer."""
from __future__ import annotations

import pytest

from skills.clawpathy_autoresearch.scorer import (
    automated_score,
    combine_scores,
    ScoreResult,
)


@pytest.fixture
def ground_truth() -> dict:
    return {
        "lead_variants": [
            {
                "rsid": "rs356182",
                "gene": "SNCA",
                "p_value_order": -45,
                "effect_direction": "risk",
                "or_range": [1.25, 1.40],
            },
            {
                "rsid": "rs34311866",
                "gene": "LRRK2",
                "p_value_order": -30,
                "effect_direction": "risk",
                "or_range": [1.15, 1.30],
            },
        ],
        "qualitative_findings": [
            "Identified immune system enrichment",
            "Lysosomal pathway implicated",
            "90 independent risk loci total",
        ],
        "total_loci": 90,
        "ancestry": "European",
    }


@pytest.fixture
def perfect_agent_output() -> dict:
    return {
        "variants_found": [
            {
                "rsid": "rs356182",
                "gene": "SNCA",
                "p_value": 1e-50,
                "effect_direction": "risk",
                "odds_ratio": 1.32,
            },
            {
                "rsid": "rs34311866",
                "gene": "LRRK2",
                "p_value": 1e-35,
                "effect_direction": "risk",
                "odds_ratio": 1.22,
            },
        ],
        "total_loci_reported": 90,
        "qualitative_summary": (
            "Found immune system enrichment and lysosomal pathway involvement. "
            "Identified 90 independent risk loci in European ancestry."
        ),
    }


@pytest.fixture
def partial_agent_output() -> dict:
    return {
        "variants_found": [
            {
                "rsid": "rs356182",
                "gene": "SNCA",
                "p_value": 1e-50,
                "effect_direction": "risk",
                "odds_ratio": 1.32,
            },
        ],
        "total_loci_reported": 45,
        "qualitative_summary": "Found some immune enrichment.",
    }


def test_automated_score_perfect(ground_truth, perfect_agent_output):
    result = automated_score(ground_truth, perfect_agent_output)
    assert isinstance(result, ScoreResult)
    assert result.variant_recovery == 1.0
    assert result.direction_accuracy == 1.0
    assert result.effect_size_accuracy == 1.0
    assert result.locus_count_accuracy == 1.0
    assert result.automated_total == pytest.approx(10.0, abs=0.01)


def test_automated_score_partial(ground_truth, partial_agent_output):
    result = automated_score(ground_truth, partial_agent_output)
    assert result.variant_recovery == 0.5
    assert result.locus_count_accuracy == 0.5
    assert result.automated_total < 10.0


def test_automated_score_empty_output(ground_truth):
    empty = {"variants_found": [], "total_loci_reported": 0, "qualitative_summary": ""}
    result = automated_score(ground_truth, empty)
    assert result.variant_recovery == 0.0
    assert result.automated_total == pytest.approx(0.0, abs=0.01)


def test_combine_scores():
    auto = 7.5
    llm = 8.0
    combined = combine_scores(auto, llm, auto_weight=0.6, llm_weight=0.4)
    expected = 0.6 * 7.5 + 0.4 * 8.0  # 7.7
    assert combined == pytest.approx(expected, abs=0.01)


def test_combine_scores_clamps_to_scale():
    combined = combine_scores(10.0, 10.0)
    assert combined <= 10.0
    combined = combine_scores(0.0, 0.0)
    assert combined >= 0.0
