"""Regression tests: nf-core-native (snake_case) flag spellings reach the wrappers.

Repository bug being guarded:

The launcher's INT-001 extra-args allowlist filter (``run_skill``) matches flag
tokens *exactly*, and its vocabulary is fully hyphenated (``--skip-tools``,
``--fasta-fai``, ``--known-indels`` …). nf-core's own pipeline parameters are
snake_case (``--skip_tools``, ``--fasta_fai`` …), so a user copying an upstream
nf-core command had the token silently dropped by the filter before it ever
reached the wrapper.

After the fix, for the three nf-core pipeline skills the filter treats the two
spellings as equivalent and forwards the wrapper's canonical hyphen spelling
(which every wrapper parser already accepts). Non-pipeline skills are untouched.
"""
from __future__ import annotations

from importlib import import_module
import types

import pytest

from clawbio import cli


PIPELINES = ["scrnaseq-pipeline", "rnaseq-pipeline", "sarek-pipeline"]


def _forwarded_cmd(monkeypatch, tmp_path, skill, extra_args):
    """Return the argv ``run_skill`` would hand to its skill-process adapter.

    ``run_skill`` wraps the subprocess in a broad ``except`` that folds any error
    into its result dict, so the stub captures the command and returns a benign
    completed-process stand-in instead of raising.
    """
    captured: dict[str, list[str]] = {}

    def fake_run(cmd, *args, **kwargs):
        captured["cmd"] = list(cmd)
        return types.SimpleNamespace(returncode=0, stdout="", stderr="")

    launcher = import_module("clawbio.common.nfcore.launcher")

    monkeypatch.setattr(cli.subprocess, "run", fake_run)
    monkeypatch.setattr(launcher, "run_wrapper", fake_run)
    monkeypatch.setattr(
        cli.subprocess, "Popen",
        lambda *args, **kwargs: pytest.fail("A routing test must not start a real skill"),
    )
    cli.run_skill(
        skill,
        demo=True,
        output_dir=str(tmp_path / "out"),
        extra_args=list(extra_args),
    )
    assert "cmd" in captured, "Skill-process adapter was never reached"
    return captured["cmd"]


def _value_after(cmd, flag):
    idx = cmd.index(flag)
    return cmd[idx + 1]


def test_sarek_underscore_value_flag_forwarded_as_hyphen(monkeypatch, tmp_path):
    """``--skip_tools baserecalibrator`` (nf-core native) must reach the wrapper as
    the canonical ``--skip-tools baserecalibrator`` instead of being dropped."""
    cmd = _forwarded_cmd(
        monkeypatch, tmp_path, "sarek-pipeline",
        ["--skip_tools", "baserecalibrator"],
    )
    assert "--skip-tools" in cmd, f"underscore spelling dropped: {cmd!r}"
    assert "--skip_tools" not in cmd
    assert _value_after(cmd, "--skip-tools") == "baserecalibrator"


def test_sarek_underscore_equals_form_forwarded_as_hyphen(monkeypatch, tmp_path):
    """The ``--flag=value`` form must be canonicalised on the flag part only."""
    cmd = _forwarded_cmd(
        monkeypatch, tmp_path, "sarek-pipeline",
        ["--skip_tools=baserecalibrator"],
    )
    assert "--skip-tools=baserecalibrator" in cmd, cmd


def test_sarek_underscore_bool_flag_does_not_swallow_next_token(monkeypatch, tmp_path):
    """A value-less nf-core bool (``--joint_germline``) must forward as
    ``--joint-germline`` without consuming the following flag as its value."""
    cmd = _forwarded_cmd(
        monkeypatch, tmp_path, "sarek-pipeline",
        ["--joint_germline", "--tools", "haplotypecaller"],
    )
    assert "--joint-germline" in cmd, cmd
    assert "--tools" in cmd and _value_after(cmd, "--tools") == "haplotypecaller"


def test_sarek_hyphen_spelling_still_works(monkeypatch, tmp_path):
    """The existing hyphen spelling must remain unaffected."""
    cmd = _forwarded_cmd(
        monkeypatch, tmp_path, "sarek-pipeline",
        ["--skip-tools", "baserecalibrator"],
    )
    assert "--skip-tools" in cmd
    assert _value_after(cmd, "--skip-tools") == "baserecalibrator"


def test_unknown_underscore_flag_still_dropped(monkeypatch, tmp_path):
    """Canonicalisation must not blanket-accept unknown flags — an invented
    ``--totally_made_up`` is still filtered out."""
    cmd = _forwarded_cmd(
        monkeypatch, tmp_path, "sarek-pipeline",
        ["--totally_made_up", "value"],
    )
    assert "--totally-made-up" not in cmd
    assert "--totally_made_up" not in cmd


@pytest.mark.parametrize("skill", PIPELINES)
def test_every_pipeline_accepts_native_underscore_spelling(monkeypatch, tmp_path, skill):
    """Parity: for each pipeline pick a value-taking hyphenated flag from its own
    allowlist, and prove its nf-core underscore spelling is forwarded as the
    hyphen form. Data-driven so it stays correct as allowlists evolve."""
    info = cli.SKILLS[skill]
    allowed = info.get("allowed_extra_flags", set())
    without_values = info.get("allowed_extra_flags_without_values", set())
    candidates = sorted(
        f for f in allowed
        if f.startswith("--") and "-" in f[2:] and f not in without_values
    )
    assert candidates, f"{skill} has no multi-word value flag to exercise"
    hyphen = candidates[0]
    underscore = "--" + hyphen[2:].replace("-", "_")

    cmd = _forwarded_cmd(monkeypatch, tmp_path, skill, [underscore, "VALUE"])

    assert hyphen in cmd, f"{skill}: {underscore} not forwarded as {hyphen}: {cmd!r}"
    assert _value_after(cmd, hyphen) == "VALUE"


def test_non_pipeline_skill_not_canonicalised(monkeypatch, tmp_path):
    """Scope guard: the canonicalisation is limited to nf-core pipeline skills, so
    a non-pipeline skill keeps its exact-match filter (clinpgx's ``--no-cache`` is
    not reachable via ``--no_cache``)."""
    cmd = _forwarded_cmd(
        monkeypatch, tmp_path, "clinpgx",
        ["--no_cache", "x"],
    )
    assert "--no-cache" not in cmd
    assert "--no_cache" not in cmd
