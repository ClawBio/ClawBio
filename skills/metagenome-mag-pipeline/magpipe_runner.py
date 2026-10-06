"""Build and run the snakemake command.

The command is modelled on upstream's `test/run_test.sh`, which is the author's
own smoke-test invocation, so the wrapper runs the workflow the way its author
verified it runs.

`run_snakemake` starts the pipeline in its own process group: snakemake spawns
conda, and conda spawns the tools. Killing only the parent on a timeout leaves
those running, holding the output directory and the disk.
"""

from __future__ import annotations

import os
import shlex
import signal
import subprocess
import sys
import time
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_errors import ErrorCode, SkillError

STAGE = "execution"

_TERMINATION_GRACE_SECONDS = 10
_LOG_TAIL_LINES = 40

# A conda environment build can die with EDQUOT on a shared filesystem that
# reports terabytes free; the identical rerun then succeeds. Retry only when the
# log shows BOTH the failed environment build and a quota marker, never on a
# bare non-zero exit.
_ENV_QUOTA_MARKERS = ("Errno 122", "Disk quota exceeded")
_LOCK_MARKERS = ("LockException", "Directory cannot be locked")


def _log_mentions(tail: list[str], *markers: str) -> bool:
    text = "\n".join(tail)
    return any(marker in text for marker in markers)


def _looks_like_env_quota(tail: list[str]) -> bool:
    return "CreateCondaEnvironmentException" in "\n".join(
        tail
    ) and _log_mentions(tail, *_ENV_QUOTA_MARKERS)


def build_snakemake_command(
    *,
    pipeline_dir: Path,
    config_path: Path,
    work_dir: Path,
    conda_prefix: Path,
    cores: int,
    mem_mb: int,
    dry_run: bool,
) -> tuple[list[str], str]:
    """The argv to run, and the same command as a copy-pasteable string."""
    argv = [
        "snakemake",
        "-s",
        str(Path(pipeline_dir) / "workflow" / "Snakefile"),
        "--configfile",
        str(Path(config_path)),
        "--directory",
        str(Path(work_dir)),
        "--use-conda",
        "--conda-frontend",
        "conda",
        "--conda-prefix",
        str(Path(conda_prefix)),
        "--cores",
        str(int(cores)),
        "--resources",
        f"mem_mb={int(mem_mb)}",
        "--rerun-incomplete",
        "--printshellcmds",
    ]
    if dry_run:
        argv.append("--dry-run")
    return argv, shlex.join(argv)


def _terminate_process_group(process: subprocess.Popen) -> None:
    """SIGTERM the whole group, then SIGKILL whatever is left."""
    try:
        group = os.getpgid(process.pid)
    except (OSError, ProcessLookupError):
        process.kill()
        return
    for number in (signal.SIGTERM, signal.SIGKILL):
        try:
            os.killpg(group, number)
        except (OSError, ProcessLookupError):
            break
        try:
            process.wait(timeout=_TERMINATION_GRACE_SECONDS)
            return
        except subprocess.TimeoutExpired:
            continue
    try:
        process.kill()
        process.wait(timeout=_TERMINATION_GRACE_SECONDS)
    except (OSError, subprocess.SubprocessError):
        pass


def _log_tail(path: Path, lines: int = _LOG_TAIL_LINES) -> list[str]:
    try:
        content = path.read_text(encoding="utf-8", errors="replace").splitlines()
    except OSError:
        return []
    return content[-lines:]


def run_snakemake(
    argv: list[str], *, log_path: Path, timeout_seconds: float, max_attempts: int = 2
) -> int:
    """Run snakemake, streaming combined output to `log_path`. Raises on failure.

    A failure whose log shows a conda environment build dying on disk quota is
    retried (up to `max_attempts` total attempts): that error is transient
    filesystem pressure, not a broken configuration. Every other non-zero exit
    raises on the first attempt.
    """
    log_path = Path(log_path)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    popen_kwargs: dict = {}
    if sys.platform != "win32":
        popen_kwargs["start_new_session"] = True

    attempts = max(1, int(max_attempts))
    started = time.monotonic()
    for attempt in range(1, attempts + 1):
        with log_path.open("w" if attempt == 1 else "a", encoding="utf-8") as handle:
            if attempt == 1:
                handle.write(f"$ {shlex.join(argv)}\n")
            else:
                handle.write(
                    f"\n$ retry attempt {attempt}/{attempts} after a transient "
                    f"environment failure\n$ {shlex.join(argv)}\n"
                )
            handle.flush()
            try:
                process = subprocess.Popen(
                    argv,
                    stdout=handle,
                    stderr=subprocess.STDOUT,
                    **popen_kwargs,
                )
            except OSError as error:
                raise SkillError(
                    stage=STAGE,
                    error_code=ErrorCode.EXECUTION_FAILED,
                    message=f"Could not launch snakemake: {error}",
                    fix="Install snakemake 9.x and conda, then re-run.",
                    details={"command": argv[0] if argv else "", "error": str(error)},
                ) from error
            try:
                returncode = process.wait(timeout=timeout_seconds)
            except subprocess.TimeoutExpired:
                _terminate_process_group(process)
                raise SkillError(
                    stage=STAGE,
                    error_code=ErrorCode.EXECUTION_TIMEOUT,
                    message=(
                        f"snakemake did not finish within "
                        f"{timeout_seconds / 3600:.2f} hours and was terminated."
                    ),
                    fix="Raise --timeout-hours, or run on a machine with more cores and "
                    "memory. Snakemake resumes from <output>/pipeline/workdir with "
                    "--rerun-incomplete, so re-running continues rather than restarting.",
                    details={
                        "timeout_seconds": timeout_seconds,
                        "log": str(log_path),
                        "log_tail": _log_tail(log_path),
                    },
                ) from None

        if returncode == 0:
            return returncode
        tail = _log_tail(log_path)
        if attempt < attempts and _looks_like_env_quota(tail):
            continue
        if _log_mentions(tail, *_LOCK_MARKERS):
            raise SkillError(
                stage=STAGE,
                error_code=ErrorCode.EXECUTION_FAILED,
                message=(
                    f"snakemake exited with status {returncode}: the working "
                    "directory is locked, likely by another run still using it."
                ),
                fix=(
                    "Wait for the other run to finish. If no snakemake process is "
                    "using that directory, clear the stale lock by re-running "
                    "snakemake with --unlock on the same directory, then re-run "
                    "this command; it resumes with --rerun-incomplete."
                ),
                details={
                    "exit_code": returncode,
                    "duration_seconds": round(time.monotonic() - started, 2),
                    "log": str(log_path),
                    "log_tail": tail,
                },
            )

        raise SkillError(
            stage=STAGE,
            error_code=ErrorCode.EXECUTION_FAILED,
            message=f"snakemake exited with status {returncode}.",
            fix=(
                f"Read {log_path.name} for the failing rule. Snakemake keeps its "
                "state under <output>/pipeline/workdir, so re-running the same "
                "command resumes rather than starting again."
            ),
            details={
                "exit_code": returncode,
                "duration_seconds": round(time.monotonic() - started, 2),
                "log": str(log_path),
                "log_tail": tail,
            },
        )
    return returncode
