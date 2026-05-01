"""Smoke tests for the clawbio.py CLI entry point."""

import subprocess
import sys
from pathlib import Path

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


def test_skill_flags_pass_through():
    """Unknown flags must reach the skill, not be rejected by clawbio.py."""
    result = _run("run", "pharmgx", "--unknown-future-flag", "value", "--help")
    # argparse in the skill will handle it; clawbio.py must not crash first
    assert "argparse.ArgumentError" not in result.stderr
    assert "conflicting option string" not in result.stderr
