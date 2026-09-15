"""Tests for scripts/generate_skills_summary.py.

The generator's value is that it fails loudly. A silently-passing validator
would let a new skill drop out of SKILLS_SUMMARY.md unnoticed, which is the
exact failure the script exists to prevent, so the negative cases below matter
more than the happy path.
"""

from __future__ import annotations

import copy
import importlib.util
import subprocess
import sys
from pathlib import Path

import pytest

REPO_ROOT = Path(__file__).resolve().parents[1]
SCRIPT_PATH = REPO_ROOT / "scripts" / "generate_skills_summary.py"


def _load_generator():
    spec = importlib.util.spec_from_file_location("generate_skills_summary", SCRIPT_PATH)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def gen():
    return _load_generator()


@pytest.fixture(scope="module")
def skills(gen):
    return gen.collect_skills()


@pytest.fixture(scope="module")
def taxonomy(gen):
    return gen.load_taxonomy()


# --------------------------------------------------------------------------- #
# Coverage validation — the negative cases
# --------------------------------------------------------------------------- #


def test_current_taxonomy_covers_every_skill(gen, taxonomy, skills):
    """The committed taxonomy must classify the live registry exactly."""
    gen.validate_coverage(taxonomy, skills)


def test_unclassified_skill_is_rejected(gen, taxonomy, skills):
    """Adding a skill without classifying it must break, not silently drop it."""
    augmented = dict(skills)
    augmented["brand-new-skill"] = {
        "name": "brand-new-skill",
        "tier": gen.TIER_AGENT,
        "description": "",
        "script_name": None,
        "script_exists": None,
        "folder": "brand-new-skill",
        "has_script": True,
        "spec_only": False,
    }
    with pytest.raises(SystemExit) as excinfo:
        gen.validate_coverage(taxonomy, augmented)
    assert excinfo.value.code == 1


def test_stale_taxonomy_entry_is_rejected(gen, taxonomy, skills):
    """A domain entry naming a deleted skill must break."""
    stale = copy.deepcopy(taxonomy)
    stale["domains"][0]["members"].append("skill-that-was-deleted")
    with pytest.raises(SystemExit) as excinfo:
        gen.validate_coverage(stale, skills)
    assert excinfo.value.code == 1


def test_double_assignment_is_rejected(gen, taxonomy, skills):
    """A skill listed under two domains would inflate the totals."""
    duplicated = copy.deepcopy(taxonomy)
    already_assigned = duplicated["domains"][0]["members"][0]
    duplicated["domains"][1]["members"].append(already_assigned)
    with pytest.raises(SystemExit) as excinfo:
        gen.validate_coverage(duplicated, skills)
    assert excinfo.value.code == 1


# --------------------------------------------------------------------------- #
# Registry parity
# --------------------------------------------------------------------------- #


def test_tiers_match_the_cli_registry(gen, skills):
    """The report's two tiers must be the same ones `clawbio.py list` prints."""
    from clawbio.cli import SKILLS, _skill_md_only_directories

    registered = {name for name, s in skills.items() if s["tier"] == gen.TIER_REGISTERED}
    agent = {name for name, s in skills.items() if s["tier"] == gen.TIER_AGENT}

    assert registered == set(SKILLS)
    assert agent == {path.name for path in _skill_md_only_directories()}


def test_alias_and_folder_namespaces_do_not_collide(gen, skills):
    """The taxonomy keys on one flat namespace; a collision would mislabel a skill."""
    gen.assert_flat_namespace(skills)


def test_spec_only_classification_matches_cli_indicator(gen, skills):
    """`spec_only` must mirror the CLI's `[spec only]` / `[has script]` badge."""
    from clawbio.cli import _skill_md_only_directories

    for skill_dir in _skill_md_only_directories():
        entry = skills[skill_dir.name]
        assert entry["spec_only"] is not entry["has_script"]
        assert entry["has_script"] == gen._has_python_script(skill_dir)


# --------------------------------------------------------------------------- #
# Rendering
# --------------------------------------------------------------------------- #


def test_render_is_deterministic(gen, taxonomy, skills):
    """Reruns must not churn the committed file."""
    assert gen.render(taxonomy, skills) == gen.render(taxonomy, skills)


def test_every_skill_appears_in_the_report(gen, taxonomy, skills):
    content = gen.render(taxonomy, skills)
    missing = [name for name in skills if f"`{name}`" not in content]
    assert not missing, f"absent from rendered report: {missing}"


def test_missing_scripts_are_flagged(gen, taxonomy, skills):
    """A registered skill whose script is gone must be called out, not hidden."""
    broken = [
        s["name"]
        for s in skills.values()
        if s["tier"] == gen.TIER_REGISTERED and not s["script_exists"]
    ]
    content = gen.render(taxonomy, skills)
    if broken:
        assert "script missing" in content
        for name in broken:
            assert name in content
    else:
        assert "Every registered skill resolves to a script on disk." in content


def test_totals_are_internally_consistent(gen, taxonomy, skills):
    """The stated total must equal the sum of the per-domain counts."""
    content = gen.render(taxonomy, skills)
    per_domain = sum(len(d["members"]) for d in taxonomy["domains"])
    assert per_domain == len(skills)
    assert f"**Total: {len(skills)} skills" in content


# --------------------------------------------------------------------------- #
# CLI surface
# --------------------------------------------------------------------------- #


def test_check_mode_reports_committed_file_is_current():
    """`--check` guards against a stale SKILLS_SUMMARY.md landing on main."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--check"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, (
        "SKILLS_SUMMARY.md is out of date; "
        f"run python scripts/generate_skills_summary.py\n{result.stderr}"
    )


def test_output_to_stdout_does_not_touch_the_committed_file(tmp_path):
    before = (REPO_ROOT / "SKILLS_SUMMARY.md").read_text(encoding="utf-8")
    result = subprocess.run(
        [sys.executable, str(SCRIPT_PATH), "--output", "-"],
        cwd=REPO_ROOT,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    assert result.stdout.startswith("# ClawBio Skill Library")
    assert (REPO_ROOT / "SKILLS_SUMMARY.md").read_text(encoding="utf-8") == before
