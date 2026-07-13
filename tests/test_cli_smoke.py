"""Smoke tests for the clawbio.py CLI entry point."""

import subprocess
import sys
from pathlib import Path

import pytest

CLAWBIO = Path(__file__).resolve().parents[1] / "clawbio.py"


def _run(*args):
    return subprocess.run(
        [sys.executable, str(CLAWBIO), *args],
        capture_output=True,
        text=True,
    )


def test_parser_constructs_without_error():
    """Parser must build cleanly — no duplicate-flag argparse.ArgumentError."""
    result = _run("list")
    assert result.returncode == 0, result.stderr


def test_run_help_exits_zero():
    result = _run("run", "--help")
    assert result.returncode == 0, result.stderr


def test_upload_help_exits_zero():
    result = _run("upload", "--help")
    assert result.returncode == 0, result.stderr


def test_version_flag_prints_version_and_exits_zero():
    from clawbio import __version__

    result = _run("--version")
    assert result.returncode == 0, result.stderr
    assert __version__ in result.stdout


def test_skill_flags_pass_through():
    """Unknown flags must reach the skill, not be rejected by clawbio.py."""
    result = _run("run", "pharmgx", "--unknown-future-flag", "value", "--help")
    # argparse in the skill will handle it; clawbio.py must not crash first
    assert "argparse.ArgumentError" not in result.stderr
    assert "conflicting option string" not in result.stderr


NFCORE_PIPELINES = ("sarek-pipeline", "rnaseq-pipeline", "scrnaseq-pipeline")


@pytest.mark.parametrize("skill", NFCORE_PIPELINES)
def test_pipeline_help_is_delegated_to_the_wrapper_parser(skill):
    """`clawbio.py run <nf-core pipeline> --help` must show the WRAPPER's own help.

    The pipeline wrappers own large, schema-derived CLIs. The launcher only declares a
    subset of their flags as argparse options (the rest ride an allowlist), so rendering
    the launcher's own help hides real, working flags. `--allow-remote-inputs` is the
    clearest case: it works on all three pipelines, but was only visible in --help for
    sarek, because help delegation was hardcoded to sarek-pipeline alone. A user reading
    `clawbio.py run rnaseq-pipeline --help` had no way to discover it.

    Delegating for every nf-core pipeline keeps the launcher's help from drifting away
    from the wrapper parser -- which is exactly what the delegation was introduced to do.
    """
    result = _run("run", skill, "--help")
    assert result.returncode == 0
    assert "--allow-remote-inputs" in result.stdout, (
        f"{skill} --help does not document --allow-remote-inputs, but the flag works"
    )
    # A launcher-rendered help would advertise the launcher, not the wrapper.
    assert "ClawBio — Bioinformatics Skills Runner" not in result.stdout
