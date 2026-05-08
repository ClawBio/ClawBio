from __future__ import annotations

import os
import signal
import subprocess
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from errors import ErrorCode, SkillError

_PROCESS_TERMINATION_GRACE_SECONDS = 10


def execute_nextflow(
    command: list[str],
    *,
    cwd: Path,
    output_dir: Path,
    timeout_seconds: int,
) -> dict[str, object]:
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / "stdout.txt"
    stderr_path = logs_dir / "stderr.txt"

    popen_kwargs: dict = {}
    if sys.platform != "win32":
        popen_kwargs["start_new_session"] = True

    with stdout_path.open("w", encoding="utf-8") as stdout_fh, \
         stderr_path.open("w", encoding="utf-8") as stderr_fh:
        try:
            proc = subprocess.Popen(
                command,
                cwd=str(cwd),
                stdout=stdout_fh,
                stderr=stderr_fh,
                **popen_kwargs,
            )
        except OSError as exc:
            for path in (stdout_path, stderr_path):
                path.unlink(missing_ok=True)
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.EXECUTION_FAILED,
                message=f"Could not launch Nextflow: {exc}",
                fix="Ensure nextflow is installed and available on PATH.",
                details={"command": command[0] if command else "", "error": str(exc)},
            ) from exc

        try:
            exit_code = proc.wait(timeout=timeout_seconds)
        except subprocess.TimeoutExpired:
            _terminate_process_tree(proc)
            raise SkillError(
                stage="execution",
                error_code=ErrorCode.EXECUTION_FAILED,
                message="Nextflow execution timed out.",
                fix="Increase resources, inspect the pipeline, or retry with a smaller dataset.",
                details={"timeout_seconds": timeout_seconds},
            )

    if exit_code != 0:
        raise SkillError(
            stage="execution",
            error_code=ErrorCode.EXECUTION_FAILED,
            message="Nextflow execution failed.",
            fix="Inspect logs/stdout.txt and logs/stderr.txt, then correct the failing input or environment.",
            details={
                "exit_code": exit_code,
                "stdout": str(stdout_path),
                "stderr": str(stderr_path),
            },
        )

    return {
        "exit_code": exit_code,
        "stdout_path": str(stdout_path),
        "stderr_path": str(stderr_path),
    }


def _terminate_process_tree(proc: subprocess.Popen) -> None:
    if sys.platform != "win32":
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
