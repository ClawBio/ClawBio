"""The wrapper must be runnable directly, as AGENTS.md documents.

Python only puts the *script's* directory on ``sys.path``, not the repo root, so
``from clawbio.common...`` at module top would fail with ``ModuleNotFoundError:
clawbio`` unless the wrapper bootstraps the repo root itself. This guards that the
bootstrap stays in place (parity with nfcore-sarek-wrapper).
"""
from __future__ import annotations

import os
import subprocess
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent.parent
_WRAPPER = _SKILL_DIR / "nfcore_scrnaseq_wrapper.py"


def test_wrapper_help_runs_directly_without_pythonpath(tmp_path):
    """Run the wrapper from an unrelated cwd with PYTHONPATH stripped, so only the
    wrapper's own repo-root bootstrap can satisfy ``import clawbio``. It must print
    usage and exit 0 — never ModuleNotFoundError."""
    env = {k: v for k, v in os.environ.items() if k != "PYTHONPATH"}
    result = subprocess.run(
        [sys.executable, str(_WRAPPER), "--help"],
        cwd=str(tmp_path),
        env=env,
        stdin=subprocess.DEVNULL,
        capture_output=True,
        text=True,
        timeout=180,
    )
    assert "ModuleNotFoundError" not in result.stderr, (
        "Wrapper does not bootstrap the repo root; direct invocation fails:\n"
        + result.stderr
    )
    assert result.returncode == 0, f"--help exited {result.returncode}:\n{result.stderr}"
    assert "usage" in result.stdout.lower()
