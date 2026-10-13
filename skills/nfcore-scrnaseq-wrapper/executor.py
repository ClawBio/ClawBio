from __future__ import annotations

from clawbio.common.nfcore import execution as _shared_execution

import os
import signal
import subprocess
import sys
from pathlib import Path

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


_MEMORY_FAILURE_SIGNATURES = _shared_execution._MEMORY_FAILURE_SIGNATURES
_NETWORK_FAILURE_SIGNATURES = _shared_execution._NETWORK_FAILURE_SIGNATURES
_CONFIG_PARSE_FAILURE_SIGNATURES = _shared_execution._CONFIG_PARSE_FAILURE_SIGNATURES


def _read_log_tail(path: Path, limit: int = 65536) -> str:
    return _shared_execution._read_log_tail(path, limit)


def _environment_failure_hints(stdout_path: Path, stderr_path: Path, nextflow_log_path: Path | None = None) -> str:
    return _shared_execution.environment_failure_hints(stdout_path, stderr_path, nextflow_log_path)


def _cellbender_failure_hint(stdout_path: Path, stderr_path: Path) -> str:
    """Hint when CellBender background removal is the failing process.

    CELLBENDER_REMOVEBACKGROUND estimates ambient RNA from the droplet-count
    distribution and errors on very small or test datasets (a common symptom is
    ``IndexError: index -100 is out of bounds``). CellBender is optional, so a
    re-run with ``--skip-cellbender`` lets the rest of the pipeline produce its
    count matrices. Best-effort: empty string when CellBender is not the failing
    process (or the logs are unreadable).
    """
    blob = (_read_log_tail(stdout_path) + "\n" + _read_log_tail(stderr_path)).lower()
    if "cellbender_removebackground" in blob and "error executing process" in blob:
        return (
            " The failing process is CellBender background removal, which estimates "
            "ambient RNA from the droplet distribution and does not work on very small "
            "or test datasets (a common symptom is 'IndexError: index -100 is out of "
            "bounds'). CellBender is optional — re-run with --skip-cellbender to let the "
            "rest of the pipeline produce its count matrices."
        )
    return ""


def execute_nextflow(command: list[str], *, cwd: Path, output_dir: Path, timeout_seconds: int | None) -> dict[str, object]:
    return _shared_execution.execute_nextflow(
        command, cwd=cwd, output_dir=output_dir, timeout_seconds=timeout_seconds,
        timeout_fix='Increase resources, inspect the pipeline, or retry with a smaller dataset.', terminate_process=_terminate_process_tree,
        terminate_on_interrupt=True,
        failure_hint=lambda stdout, stderr, log: _macos_tmp_failure_hint(output_dir) + _environment_failure_hints(stdout, stderr, log) + _cellbender_failure_hint(stdout, stderr),
    )


def _macos_tmp_failure_hint(output_dir: Path) -> str:
    """Extra fix hint when a run fails with --output under /tmp on macOS.

    Colima/Docker on macOS does not share /tmp into the VM, a frequent cause of
    '.command.run: No such file or directory' mid-run. Preflight already WARNs about
    this; we repeat the actionable hint on failure so the cause is obvious. Empty
    string on other platforms or when --output is not under /tmp.
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


def _terminate_process_tree(proc: subprocess.Popen) -> None:
    if sys.platform == "win32":
        # Windows has no os.killpg; ``taskkill /T`` terminates the whole child tree
        # (Nextflow plus its spawned task processes) instead of orphaning them when
        # only the parent is killed. Best-effort and guarded: fall back to
        # proc.kill() if taskkill is unavailable. Native Windows is not officially
        # supported (preflight warns); WSL2 reports "linux" and takes the POSIX path
        # below (audit F-7).
        try:
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                capture_output=True,
                timeout=_PROCESS_TERMINATION_GRACE_SECONDS,
            )
        except (OSError, subprocess.SubprocessError):
            pass
        proc.kill()
        proc.wait()
        return
    try:
        pgid = os.getpgid(proc.pid)
        os.killpg(pgid, signal.SIGTERM)
        proc.wait(timeout=_PROCESS_TERMINATION_GRACE_SECONDS)
        return
    except (OSError, ProcessLookupError, subprocess.TimeoutExpired):
        pass
    try:
        pgid = os.getpgid(proc.pid)
        os.killpg(pgid, signal.SIGKILL)
        proc.wait()
        return
    except (OSError, ProcessLookupError):
        pass
    proc.kill()
    proc.wait()
