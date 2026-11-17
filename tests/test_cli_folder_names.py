"""`clawbio run` accepts a skill's folder name as well as its CLI alias.

Agents (Hermes, OpenClaw, pi, Claude Code) and ClawHub show the folder name
(`gwas-prs`); the CLI registry is keyed by alias (`prs`). Both must work.
"""

import sys
import types

import pytest

from clawbio import cli


def _runnable_folders() -> dict[str, set[str]]:
    """Folder name -> aliases whose script lives in it, read from the filesystem."""
    folders: dict[str, set[str]] = {}
    for alias, info in cli.SKILLS.items():
        folders.setdefault(info["script"].parent.name, set()).add(alias)
    return {f: a for f, a in folders.items() if (cli.SKILLS_DIR / f / "SKILL.md").is_file()}


@pytest.fixture
def captured_run(monkeypatch):
    """Stub the skill subprocess; record the command and timeout it was given."""
    seen = {}

    def fake_run(cmd, *args, **kwargs):
        seen["cmd"] = [str(part) for part in cmd]
        seen["timeout"] = kwargs.get("timeout")
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    return seen


@pytest.mark.parametrize("folder", sorted(_runnable_folders()))
def test_every_runnable_folder_resolves_to_an_entry_in_that_folder(folder):
    alias = cli.resolve_skill_name(folder)
    assert alias in cli.SKILLS
    assert cli.SKILLS[alias]["script"].parent.name == folder


@pytest.mark.parametrize("alias", sorted(cli.SKILLS))
def test_aliases_resolve_to_themselves(alias):
    assert cli.resolve_skill_name(alias) == alias


@pytest.mark.parametrize(
    "folder", sorted(p.parent.name for p in cli.SKILLS_DIR.glob("*/SKILL.md"))
)
def test_no_skill_folder_is_shadowed_by_another_skills_alias(folder):
    # A folder named like another skill's alias (e.g. a new skills/prs/) would
    # resolve to that other skill. Pick a different folder name, or rename the alias.
    if folder in cli.SKILLS:
        owner = cli.SKILLS[folder]["script"].parent.name
        assert owner == folder, f"skills/{folder}/ is shadowed by alias '{folder}' of skills/{owner}/"


def test_folder_with_two_entries_resolves_to_the_first():
    # pharmgx-reporter backs both pharmgx and drugphoto; the folder means the main mode.
    assert cli.resolve_skill_name("pharmgx-reporter") == "pharmgx"


def test_unknown_name_passes_through_unchanged():
    assert cli.resolve_skill_name("definitely-not-a-real-skill") == "definitely-not-a-real-skill"


def test_run_skill_accepts_folder_name(captured_run, tmp_path):
    result = cli.run_skill("gwas-prs", demo=True, output_dir=str(tmp_path))

    assert result["success"] is True
    assert any(part.endswith("gwas_prs.py") for part in captured_run["cmd"])


def test_cli_run_accepts_folder_name_and_applies_its_timeout(captured_run, monkeypatch, tmp_path):
    # just-prs-mcp is the folder of alias just-prs, which declares a 3900 s timeout.
    monkeypatch.setattr(sys, "argv", ["clawbio", "run", "just-prs-mcp", "--demo", "--output", str(tmp_path)])
    with pytest.raises(SystemExit) as exit_info:
        cli.main()
    if exit_info.value.code not in (None, 0):
        pytest.fail(f"clawbio run just-prs-mcp exited {exit_info.value.code}")

    assert captured_run["timeout"] == cli.SKILLS["just-prs"]["default_timeout_seconds"]


def test_nfcore_help_shortcut_accepts_folder_name(captured_run, monkeypatch):
    monkeypatch.setattr(sys, "argv", ["clawbio", "run", "nfcore-sarek-wrapper", "--help"])
    cli.main()

    assert captured_run["cmd"][-1] == "--help"
    assert captured_run["cmd"][1] == str(cli.SKILLS["sarek-pipeline"]["script"])


def test_orchestrator_routes_a_folder_name_it_does_not_map(captured_run, tmp_path):
    # route_to_clawbio passes names missing from its SKILL_REGISTRY_MAP straight to run_skill.
    import importlib.util

    path = cli.SKILLS_DIR / "bio-orchestrator" / "orchestrator.py"
    spec = importlib.util.spec_from_file_location("orchestrator_under_test", path)
    orchestrator = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(orchestrator)
    assert "nfcore-sarek-wrapper" not in orchestrator.SKILL_REGISTRY_MAP

    results = orchestrator.route_to_clawbio(["nfcore-sarek-wrapper"], output_dir=str(tmp_path))

    assert results["nfcore-sarek-wrapper"]["success"] is True
    assert str(cli.SKILLS["sarek-pipeline"]["script"]) in captured_run["cmd"]


def test_folder_without_a_run_entry_says_so():
    folder = cli._skill_md_only_directories()[0].name
    result = cli.run_skill(folder, demo=True)

    assert result["success"] is False
    assert "no `clawbio run` entry" in result["stderr"]
    assert "SKILL.md" in result["stderr"]


def test_typo_suggests_close_names():
    result = cli.run_skill("gwas-pr", demo=True)

    assert result["success"] is False
    assert "Unknown skill 'gwas-pr'" in result["stderr"]
    assert "gwas-prs" in result["stderr"]
