"""ClawBio-only integration: DnaSP's INTENTS.json descriptor for the intent planner.

This file exists only in the ClawBio copy of the skill; the standalone package
and the canonical mirror do not ship INTENTS.json.
"""

import json
from pathlib import Path

import pytest

skill_intents = pytest.importorskip("clawbio.skill_intents")

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]


def _registry():
    return skill_intents.load_default_skill_registry(PROJECT_ROOT)


def test_descriptor_is_valid_and_discovered():
    raw = json.loads((SKILL_DIR / "INTENTS.json").read_text(encoding="utf-8"))
    assert raw["schema"] == skill_intents.SCHEMA
    assert raw["skill"] == "dnasp"
    descriptors = skill_intents.load_skill_intent_descriptors(_registry(), PROJECT_ROOT)
    dnasp = [d for d in descriptors if d.get("skill") == "dnasp"]
    assert len(dnasp) == 1
    assert "dnasp" in dnasp[0]["aliases"]
    assert all(r["demo_policy"] == "only_when_explicit" for r in dnasp[0]["routes"])


def test_explicit_demo_request_plans_the_dnasp_demo():
    plan = skill_intents.plan_skill_intent(
        "run the dnasp demo", None, None, [], _registry(), PROJECT_ROOT
    )
    assert plan.skill == "dnasp"
    assert plan.intent_id == "demo_report"
    assert plan.status == "planned"
    assert any("--demo" in execution.argv for execution in plan.executions)


def test_analysis_request_without_demo_wording_never_plans_the_demo():
    plan = skill_intents.plan_skill_intent(
        "compute Tajima's D on my alignment", None, None, [], _registry(), PROJECT_ROOT
    )
    assert not any("--demo" in execution.argv for execution in plan.executions)
