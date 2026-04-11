"""Tests for reproduction error scorer.

Score = mean absolute error across concrete numerical targets.
Lower is better. Zero = perfect reproduction.
"""
from __future__ import annotations

import pytest
from pathlib import Path

from skills.clawpathy_autoresearch.scorer import (
    reproduction_error,
    ErrorBreakdown,
    load_scorer,
)


@pytest.fixture
def ground_truth() -> dict:
    return {
        "lead_variants": [
            {
                "rsid": "rs356182",
                "gene": "SNCA",
                "neg_log10_p": 45.0,
                "odds_ratio": 1.33,
                "effect_allele_freq": 0.37,
                "effect_direction": "risk",
            },
            {
                "rsid": "rs34311866",
                "gene": "LRRK2",
                "neg_log10_p": 30.0,
                "odds_ratio": 1.22,
                "effect_allele_freq": 0.02,
                "effect_direction": "risk",
            },
        ],
        "total_loci": 90,
    }


@pytest.fixture
def perfect_output() -> dict:
    return {
        "variants_found": [
            {
                "rsid": "rs356182",
                "neg_log10_p": 45.0,
                "odds_ratio": 1.33,
                "effect_allele_freq": 0.37,
                "effect_direction": "risk",
            },
            {
                "rsid": "rs34311866",
                "neg_log10_p": 30.0,
                "odds_ratio": 1.22,
                "effect_allele_freq": 0.02,
                "effect_direction": "risk",
            },
        ],
        "total_loci_reported": 90,
    }


@pytest.fixture
def partial_output() -> dict:
    return {
        "variants_found": [
            {
                "rsid": "rs356182",
                "neg_log10_p": 40.0,
                "odds_ratio": 1.28,
                "effect_allele_freq": 0.35,
                "effect_direction": "risk",
            },
        ],
        "total_loci_reported": 60,
    }


def test_perfect_reproduction_is_zero(ground_truth, perfect_output):
    result = reproduction_error(ground_truth, perfect_output)
    assert isinstance(result, ErrorBreakdown)
    assert result.total == pytest.approx(0.0, abs=0.001)
    assert result.p_value_error == pytest.approx(0.0, abs=0.001)
    assert result.or_error == pytest.approx(0.0, abs=0.001)
    assert result.freq_error == pytest.approx(0.0, abs=0.001)
    assert result.locus_count_error == pytest.approx(0.0, abs=0.001)
    assert result.variant_missing_penalty == pytest.approx(0.0, abs=0.001)
    assert result.direction_error == pytest.approx(0.0, abs=0.001)


def test_partial_output_has_positive_error(ground_truth, partial_output):
    result = reproduction_error(ground_truth, partial_output)
    assert result.total > 0.0
    assert result.variant_missing_penalty > 0.0  # missed rs34311866
    assert result.locus_count_error > 0.0  # 60 vs 90


def test_empty_output_has_maximum_penalty(ground_truth):
    empty = {"variants_found": [], "total_loci_reported": 0}
    result = reproduction_error(ground_truth, empty)
    assert result.total > 0.0
    assert result.variant_missing_penalty > 0.0


def test_wrong_direction_penalised(ground_truth):
    output = {
        "variants_found": [
            {
                "rsid": "rs356182",
                "neg_log10_p": 45.0,
                "odds_ratio": 1.33,
                "effect_allele_freq": 0.37,
                "effect_direction": "protective",  # wrong
            },
            {
                "rsid": "rs34311866",
                "neg_log10_p": 30.0,
                "odds_ratio": 1.22,
                "effect_allele_freq": 0.02,
                "effect_direction": "risk",
            },
        ],
        "total_loci_reported": 90,
    }
    result = reproduction_error(ground_truth, output)
    assert result.direction_error > 0.0
    assert result.total > 0.0


def test_p_value_error_is_normalised(ground_truth):
    """P-value error should be normalised by target so huge p-values don't dominate."""
    output = {
        "variants_found": [
            {
                "rsid": "rs356182",
                "neg_log10_p": 50.0,  # off by 5 from target 45
                "odds_ratio": 1.33,
                "effect_allele_freq": 0.37,
                "effect_direction": "risk",
            },
            {
                "rsid": "rs34311866",
                "neg_log10_p": 35.0,  # off by 5 from target 30
                "odds_ratio": 1.22,
                "effect_allele_freq": 0.02,
                "effect_direction": "risk",
            },
        ],
        "total_loci_reported": 90,
    }
    result = reproduction_error(ground_truth, output)
    # Normalised: |50-45|/45 = 0.111, |35-30|/30 = 0.167, mean = 0.139
    assert result.p_value_error == pytest.approx(0.139, abs=0.01)


def test_error_breakdown_dict(ground_truth, partial_output):
    result = reproduction_error(ground_truth, partial_output)
    d = result.to_dict()
    assert "total" in d
    assert "p_value_error" in d
    assert "or_error" in d
    assert "freq_error" in d
    assert "locus_count_error" in d
    assert "variant_missing_penalty" in d
    assert "direction_error" in d


def test_load_scorer_from_file(tmp_path: Path):
    scorer_code = '''
def score(skill_output: dict, ground_truth: dict) -> float:
    return abs(skill_output.get("value", 0) - ground_truth.get("value", 0))
'''
    scorer_path = tmp_path / "scorer.py"
    scorer_path.write_text(scorer_code)

    score_fn = load_scorer(scorer_path)
    result = score_fn({"value": 10}, {"value": 7})
    assert result == pytest.approx(3.0)


def test_load_scorer_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_scorer(Path("/nonexistent/scorer.py"))


def test_load_scorer_missing_score_function(tmp_path: Path):
    scorer_path = tmp_path / "scorer.py"
    scorer_path.write_text("x = 1\n")

    with pytest.raises(AttributeError, match="score"):
        load_scorer(scorer_path)
