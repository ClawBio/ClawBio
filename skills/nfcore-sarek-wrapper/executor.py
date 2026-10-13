# ruff: noqa: E402
from __future__ import annotations

from clawbio.common.nfcore import execution as _shared_execution

from pathlib import Path
import os
import signal
import subprocess
import sys

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import is_under_tmp


_PROCESS_TERMINATION_GRACE_SECONDS = 10


def _macos_tmp_failure_hint(output_dir: Path) -> str:
    """Extra fix hint when a run fails with --output under /tmp on macOS.

    Colima/Docker on macOS does not share /tmp into the VM, a frequent cause of
    '.command.run: No such file or directory' mid-run. Preflight already WARNs about
    this; we repeat the actionable hint on failure so the cause is obvious. Empty
    string on other platforms or when --output is not under /tmp. Mirrors the
    nfcore-scrnaseq / nfcore-rnaseq executor hint so the three wrappers behave the
    same on macOS.
    """
    if sys.platform != "darwin":
        return ""
    if not is_under_tmp(output_dir):
        return ""
    return (
        " On macOS, --output is under /tmp, which Docker/Colima does not share into its VM; "
        "this commonly surfaces as '.command.run: No such file or directory'. "
        "Move --output to a path under your HOME directory and re-run."
    )


_MEMORY_FAILURE_SIGNATURES = _shared_execution._MEMORY_FAILURE_SIGNATURES
_NETWORK_FAILURE_SIGNATURES = _shared_execution._NETWORK_FAILURE_SIGNATURES
_CONFIG_PARSE_FAILURE_SIGNATURES = _shared_execution._CONFIG_PARSE_FAILURE_SIGNATURES


def _read_log_tail(path: Path, limit: int = 65536) -> str:
    return _shared_execution._read_log_tail(path, limit)


def _environment_failure_hints(stdout_path: Path, stderr_path: Path, nextflow_log_path: Path | None = None) -> str:
    return _shared_execution.environment_failure_hints(stdout_path, stderr_path, nextflow_log_path)


def execute_nextflow(command: list[str], *, cwd: Path, output_dir: Path, timeout_seconds: int | None) -> dict[str, object]:
    return _shared_execution.execute_nextflow(
        command, cwd=cwd, output_dir=output_dir, timeout_seconds=timeout_seconds,
        timeout_fix='Increase resources, inspect the pipeline, or retry with a smaller dataset.', terminate_process=_terminate_process_tree,
        terminate_on_interrupt=False,
        failure_hint=lambda stdout, stderr, log: _macos_tmp_failure_hint(output_dir) + _environment_failure_hints(stdout, stderr, log),
    )


def _terminate_process_tree(proc: subprocess.Popen) -> None:
    _shared_execution.terminate_process_tree(proc, grace_seconds=_PROCESS_TERMINATION_GRACE_SECONDS)
