"""Characterisation tests for the INT-001 extra-args allowlist in ``run_skill``.

The filter silently *drops* anything outside a skill's ``allowed_extra_flags``
(it does not raise or report). These tests pin that behaviour at the argv handed
to ``subprocess.run`` so the filter can be moved out of ``run_skill`` without
changing what reaches a skill.
"""
from __future__ import annotations

import types

import pytest

from clawbio import cli


def _forwarded_cmd(monkeypatch, tmp_path, skill, extra_args):
    captured: dict[str, list[str]] = {}

    def fake_run(cmd, *args, **kwargs):
        captured["cmd"] = list(cmd)
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    cli.run_skill(skill, demo=True, output_dir=str(tmp_path / "out"), extra_args=list(extra_args))
    assert "cmd" in captured, "subprocess.run was never reached"
    return captured["cmd"]


# pharmgx has an empty allowlist, prs a value-taking one, bigquery has valueless flags.
SKILLS_UNDER_TEST = ["pharmgx", "prs", "bigquery", "clinpgx"]


@pytest.mark.parametrize("skill", SKILLS_UNDER_TEST)
def test_unlisted_flag_and_its_value_are_dropped(monkeypatch, tmp_path, skill):
    cmd = _forwarded_cmd(monkeypatch, tmp_path, skill, ["--not-a-real-flag", "evil-value"])
    assert "--not-a-real-flag" not in cmd
    assert "evil-value" not in cmd


@pytest.mark.parametrize("skill", SKILLS_UNDER_TEST)
def test_runner_owned_flags_cannot_be_overridden(monkeypatch, tmp_path, skill):
    cmd = _forwarded_cmd(
        monkeypatch,
        tmp_path,
        skill,
        ["--input", "/etc/passwd", "--output=/tmp/evil", "--demo"],
    )
    assert "/etc/passwd" not in cmd
    assert "/tmp/evil" not in " ".join(cmd)
    assert cmd.count("--output") == 1 and cmd[cmd.index("--output") + 1] == str(tmp_path / "out")
    assert cmd.count("--demo") <= 1


def test_allowed_value_flag_forwarded_with_its_value(monkeypatch, tmp_path):
    cmd = _forwarded_cmd(monkeypatch, tmp_path, "clinpgx", ["--gene", "CYP2D6", "--bogus", "x"])
    assert cmd[cmd.index("--gene") + 1] == "CYP2D6"
    assert "--bogus" not in cmd and "x" not in cmd


def test_allowed_equals_form_forwarded(monkeypatch, tmp_path):
    cmd = _forwarded_cmd(monkeypatch, tmp_path, "clinpgx", ["--gene=CYP2D6"])
    assert "--gene=CYP2D6" in cmd


def test_blocked_demo_does_not_swallow_next_flag(monkeypatch, tmp_path):
    cmd = _forwarded_cmd(monkeypatch, tmp_path, "clinpgx", ["--demo", "--gene", "CYP2D6"])
    assert cmd[cmd.index("--gene") + 1] == "CYP2D6"


def test_blocked_flag_after_first_position_skips_only_itself(monkeypatch, tmp_path):
    cmd = _forwarded_cmd(monkeypatch, tmp_path, "clinpgx", ["--no-cache", "--output=/tmp/evil", "--gene", "CYP2D6"])
    assert cmd[cmd.index("--gene") + 1] == "CYP2D6"


def test_valueless_flag_does_not_swallow_next_token(monkeypatch, tmp_path):
    cmd = _forwarded_cmd(monkeypatch, tmp_path, "bigquery", ["--dry-run", "swallowed"])
    assert "--dry-run" in cmd
    assert "swallowed" not in cmd
