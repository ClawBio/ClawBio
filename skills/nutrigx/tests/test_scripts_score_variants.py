"""Tests for the 2.0 pipeline step skills/nutrigx/scripts/score_variants.py."""
import importlib.util
import json
from pathlib import Path

import pytest
from hypothesis import given, strategies as st

SKILL = Path(__file__).resolve().parent.parent
_spec = importlib.util.spec_from_file_location(
    "nutrigx_scripts_score_variants", SKILL / "scripts" / "score_variants.py"
)
sv = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(sv)

PANEL = json.loads((SKILL / "data" / "snp_panel.json").read_text())


def _call(risk_count):
    return {"status": "found", "risk_count": risk_count, "normalised": "AG"}


def test_het_rs4988235_lactose_is_low():
    res = sv.compute_nutrient_risk_scores({"rs4988235": _call(1)}, PANEL)
    assert res["lactose"]["score"] == 0.0
    assert res["lactose"]["category"] == "Low"


def test_hom_risk_rs4988235_lactose_is_elevated():
    res = sv.compute_nutrient_risk_scores({"rs4988235": _call(2)}, PANEL)
    assert res["lactose"]["category"] == "Elevated"


@given(st.sampled_from(["additive", "dominant_protective"]))
def test_monotone_in_risk_count(mode):
    scores = [sv.snp_raw_score(c, inheritance=mode) for c in (0, 1, 2)]
    assert scores == sorted(scores)
    assert all(0.0 <= s <= 1.0 for s in scores)


@given(st.sampled_from([0, 1]))
def test_dominant_protective_zero_below_two(rc):
    assert sv.snp_raw_score(rc, inheritance="dominant_protective") == 0.0
