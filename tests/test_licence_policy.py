"""Licence policy for skills, enforced against the generated catalogue.

Why this exists: a PROPRIETARY skill and a GPL-3.0 skill inside an MIT wheel is
the first thing a procurement reader hits, and a blank licence is worse because
nobody can tell what they are allowed to do. The policy gates the CODE licence
only. Data licences are reported, not gated, because they are declared by 97
contributors and a gate none of them can pass today is a gate that gets
disabled tomorrow.
"""
from __future__ import annotations

import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
CATALOG = ROOT / "skills" / "catalog.json"

# Code licences that may ship inside the MIT wheel.
WHEEL_LICENCES = {"MIT", "Apache-2.0"}

# Skills whose code licence is incompatible with the wheel. They remain in the
# repository as folders labelled by licence and are excluded from the wheel.
# Every entry needs a reason so the exception cannot grow silently.
FOLDER_ONLY = {
    "fastreer": "GPL-3.0: copyleft cannot be redistributed inside an MIT wheel",
    "wes-clinical-report-en": "PROPRIETARY: declared by the contributor; not redistributable",
}


def _skills() -> list[dict]:
    return json.loads(CATALOG.read_text())["skills"]


def test_every_skill_declares_a_code_licence():
    blank = [s["name"] for s in _skills() if not (s.get("license") or "").strip()]
    assert blank == [], f"skills with no code licence: {blank}"


def test_wheel_skills_carry_a_wheel_compatible_licence():
    offenders = [
        (s["name"], s.get("license"))
        for s in _skills()
        if s["name"] not in FOLDER_ONLY and (s.get("license") or "") not in WHEEL_LICENCES
    ]
    assert offenders == [], f"licences that cannot ship in the MIT wheel: {offenders}"


def test_folder_only_exceptions_are_real_and_still_needed():
    by_name = {s["name"]: s for s in _skills()}
    for name, reason in FOLDER_ONLY.items():
        assert name in by_name, f"{name} is listed as an exception but is not in the catalogue"
        lic = by_name[name].get("license") or ""
        assert lic not in WHEEL_LICENCES, (
            f"{name} is now {lic!r}; remove it from FOLDER_ONLY, the exception is stale")
        assert lic.split(":")[0].strip() in reason, f"{name}: reason must name the licence"


def test_wheel_exclusion_matches_the_policy():
    """The build hook and the policy must name the same skills, or one of them
    is lying about what ships."""
    # Read the constant from the AST rather than importing the module: hatch_build
    # imports hatchling, a build-time dependency that need not be installed
    # wherever the tests run.
    import ast
    tree = ast.parse((ROOT / "hatch_build.py").read_text())
    values = [
        ast.literal_eval(node.value)
        for node in ast.walk(tree)
        if isinstance(node, ast.Assign)
        and any(isinstance(t, ast.Name) and t.id == "WHEEL_EXCLUDED_SKILLS" for t in node.targets)
    ]
    assert values, "WHEEL_EXCLUDED_SKILLS not defined in hatch_build.py"
    assert set(values[0]) == set(FOLDER_ONLY)


def test_data_licence_is_reported_not_gated(capsys):
    """Published count, so the gap is visible on every run without blocking anyone."""
    undeclared = [s["name"] for s in _skills()
                  if (s.get("data_license") or "none").lower() in {"none", ""}]
    print(f"\n[licence policy] data_license undeclared on {len(undeclared)} of {len(_skills())} skills")
    assert isinstance(undeclared, list)
