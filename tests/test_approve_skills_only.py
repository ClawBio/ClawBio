"""Tests for the skills-only CI auto-approval gate.

Issue #360. Fork PRs from first-time contributors park in `action_required`
and surface as "no checks reported", which is indistinguishable from a pass in
`gh pr checks`. Five PRs sat in that state in one week, and #329 went five
remediation rounds with CI having never run at all.

GitHub has no path-scoped approval policy, so this is automation rather than a
setting. `is_skills_only` is the security boundary: everything it returns True
for gets its workflow run approved with no human in the loop. It must fail
closed on anything it does not positively recognise.

What the boundary does and does not buy, stated plainly so nobody mistakes it
for more: `skills/` is exactly where contributor code lives, and CI already
executes it (pytest imports fork test files, and the skill-harness job runs a
fork's shell). So this is not "no untrusted execution". It is "a PR that runs
untrusted code cannot also rewrite the pipeline that runs it, the core package,
or dependency resolution".
"""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest


def _load():
    p = Path(__file__).resolve().parents[1] / "scripts" / "approve_skills_only_runs.py"
    spec = importlib.util.spec_from_file_location("approve_skills_only_runs", p)
    mod = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(mod)
    return mod


MOD = _load()


class TestSkillsOnlyBoundary:
    """Everything this admits runs untrusted code without human review."""

    def test_a_plain_skill_change_is_admitted(self):
        assert MOD.is_skills_only(["skills/my-skill/my_skill.py"]) is True

    def test_several_files_in_one_skill(self):
        assert MOD.is_skills_only([
            "skills/a/a.py", "skills/a/SKILL.md", "skills/a/tests/test_a.py",
        ]) is True

    def test_two_skills_at_once(self):
        assert MOD.is_skills_only(["skills/a/a.py", "skills/b/b.R"]) is True

    def test_the_generated_catalog_is_admitted(self):
        """Contributors regenerate it; it grants no execution."""
        assert MOD.is_skills_only(["skills/catalog.json", "skills/a/a.py"]) is True

    @pytest.mark.parametrize("path", [
        ".github/workflows/ci.yml",
        ".github/workflows/approve-skills-only.yml",
        "pyproject.toml",
        "uv.lock",
        "clawbio/cli.py",
        "clawbio/common/checksums.py",
        "scripts/generate_catalog.py",
        "tests/test_generate_catalog.py",
        "conftest.py",
        "Makefile",
        ".github/dependabot.yml",
    ])
    def test_anything_outside_skills_is_refused(self, path):
        assert MOD.is_skills_only([path]) is False

    def test_one_bad_file_poisons_the_whole_set(self):
        """The dangerous shape: hide a workflow edit behind real skill work."""
        assert MOD.is_skills_only([
            "skills/a/a.py",
            "skills/a/SKILL.md",
            ".github/workflows/ci.yml",
        ]) is False

    def test_an_empty_file_list_is_refused(self):
        """Fail closed. An empty list means the diff could not be read, not
        that the PR is harmless."""
        assert MOD.is_skills_only([]) is False

    @pytest.mark.parametrize("path", [
        "skills/../pyproject.toml",
        "skills/a/../../clawbio/cli.py",
        "../skills/a/a.py",
        "skills/./../../.github/workflows/ci.yml",
    ])
    def test_path_traversal_is_refused(self, path):
        assert MOD.is_skills_only([path]) is False

    @pytest.mark.parametrize("path", [
        "skills-evil/a.py",
        "skillsfoo/a.py",
        "myskills/a.py",
        "skills.py",
        "skills",
    ])
    def test_prefix_confusion_is_refused(self, path):
        """`startswith("skills")` without the separator admits `skills-evil/`."""
        assert MOD.is_skills_only([path]) is False

    @pytest.mark.parametrize("path", [
        "SKILLS/a/a.py",
        "Skills/a/a.py",
    ])
    def test_case_variants_are_refused(self, path):
        """macOS and Windows checkouts are case-insensitive; the boundary is
        not. Refuse rather than guess which filesystem resolves it."""
        assert MOD.is_skills_only([path]) is False

    def test_an_absolute_path_is_refused(self):
        assert MOD.is_skills_only(["/etc/passwd"]) is False

    def test_a_bare_file_at_repo_root_is_refused(self):
        assert MOD.is_skills_only(["README.md"]) is False


class TestRunPairing:
    """Only the run matching an open PR's CURRENT head may be approved.

    All four runs pending on this repo when the gate was written were stale
    commits on one branch, superseded by later pushes. Approving those spends
    runner minutes on code nobody is reviewing.
    """

    def test_a_run_at_an_open_prs_head_is_paired(self):
        runs = [{"id": 1, "head_sha": "aaa"}]
        prs = [{"number": 42, "headRefOid": "aaa"}]
        assert MOD.pair_runs_to_prs(runs, prs) == [(1, 42)]

    def test_a_stale_run_is_not_paired(self):
        runs = [{"id": 1, "head_sha": "old"}]
        prs = [{"number": 42, "headRefOid": "new"}]
        assert MOD.pair_runs_to_prs(runs, prs) == []

    def test_a_run_with_no_open_pr_is_not_paired(self):
        runs = [{"id": 1, "head_sha": "aaa"}]
        assert MOD.pair_runs_to_prs(runs, []) == []

    def test_pairs_are_stable_and_deduplicated(self):
        runs = [{"id": 2, "head_sha": "a"}, {"id": 1, "head_sha": "a"}]
        prs = [{"number": 7, "headRefOid": "a"}]
        assert MOD.pair_runs_to_prs(runs, prs) == [(1, 7), (2, 7)]
