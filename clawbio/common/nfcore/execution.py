"""Nextflow process mechanics with explicit caller-specific failure policy."""

from __future__ import annotations

import os
import signal
import subprocess
import sys
from collections.abc import Callable
from pathlib import Path

from .errors import SkillError

_MEMORY_FAILURE_SIGNATURES = (
    "process requirement exceeds available memory",
    "process requirement exceeds available cpus",
)
_NETWORK_FAILURE_SIGNATURES = (
    "network is unreachable",
    "java.net.connectexception",
    "java.net.unknownhostexception",
    "no route to host",
    "connection timed out",
    "temporary failure in name resolution",
)
_CONFIG_PARSE_FAILURE_SIGNATURES = (
    "unable to parse config file",
    "configparseexception",
)


def _read_log_tail(path: Path, limit: int = 65536) -> str:
    """Best-effort tail read of a captured log; never raises."""
    try:
        with path.open("r", encoding="utf-8", errors="replace") as fh:
            return fh.read()[-limit:]
    except OSError:
        return ""


def environment_failure_hints(
    stdout_path: Path, stderr_path: Path, nextflow_log_path: Path | None = None
) -> str:
    """Append actionable hints when the captured logs show a known environment
    failure — the host being smaller than a process request, or an unreachable
    network. Diagnosed from the actual Nextflow error text (never predicting
    pipeline resource requirements), so no thresholds are fabricated. Empty string
    when no known signature is present. Shared verbatim across the three wrappers.

    `.nextflow.log` is scanned in addition to stdout/stderr: when Nextflow fails while
    parsing the config (before the pipeline starts) the top-level stderr may show only
    "Unable to parse config file", while the underlying cause ("Network is unreachable")
    is recorded only in `.nextflow.log`. Without it the IPv6/NAT64 hint never shows.
    """
    parts = [_read_log_tail(stdout_path), _read_log_tail(stderr_path)]
    if nextflow_log_path is not None:
        parts.append(_read_log_tail(nextflow_log_path))
    blob = "\n".join(parts).lower()
    hints: list[str] = []
    if any(sig in blob for sig in _MEMORY_FAILURE_SIGNATURES):
        hints.append(
            " A process requested more memory or CPUs than this machine provides "
            "('Process requirement exceeds available memory' in the logs). Cap "
            "requests to your host with a Nextflow config passed via -c, e.g. a file "
            "with `process { resourceLimits = [ memory: '12.GB', cpus: 4 ] }`, or run "
            "on a larger machine or an HPC/cloud executor."
        )
    if any(sig in blob for sig in _NETWORK_FAILURE_SIGNATURES):
        hints.append(
            " Nextflow could not reach the network ('Network is unreachable' or a Java "
            "connection exception in the logs). Confirm outbound HTTPS and DNS to "
            "github.com are reachable; on IPv6-only / NAT64 hosts the JVM prefers IPv4 "
            "by default, so export NXF_OPTS='-Djava.net.preferIPv6Addresses=true' and "
            "re-run."
        )
    if any(sig in blob for sig in _CONFIG_PARSE_FAILURE_SIGNATURES):
        hints.append(
            " Nextflow could not parse the pipeline config; on nf-core this usually "
            "means it could not fetch the remote nf-core/configs 'nfcore_custom.config' "
            "(the `includeConfig ... ? <url> : '/dev/null'` line). For a fully local run "
            "(local --input, references and an already-pulled pipeline) set "
            "NXF_OFFLINE=true so Nextflow skips the remote include; otherwise ensure "
            "outbound HTTPS/DNS to raw.githubusercontent.com is reachable."
        )
    return "".join(hints)


def execute_nextflow(
    command: list[str],
    *,
    cwd: Path,
    output_dir: Path,
    timeout_seconds: int | None,
    timeout_fix: str,
    terminate_process: Callable[[subprocess.Popen], None],
    terminate_on_interrupt: bool,
    failure_hint: Callable[[Path, Path, Path], str],
) -> dict[str, object]:
    # timeout_seconds=None disables the wall-clock cap entirely (long HPC/cloud
    # runs whose walltime is enforced by the scheduler); proc.wait(timeout=None)
    # blocks until completion.
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / "stdout.txt"
    stderr_path = logs_dir / "stderr.txt"

    popen_kwargs: dict = {}
    if sys.platform != "win32":
        popen_kwargs["start_new_session"] = True

    with (
        stdout_path.open("w", encoding="utf-8") as stdout_fh,
        stderr_path.open("w", encoding="utf-8") as stderr_fh,
    ):
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
                error_code="EXECUTION_FAILED",
                message=f"Could not launch Nextflow: {exc}",
                fix="Ensure nextflow is installed and available on PATH.",
                details={"command": command[0] if command else "", "error": str(exc)},
            ) from exc

        try:
            exit_code = proc.wait(timeout=timeout_seconds)
        except KeyboardInterrupt:
            if terminate_on_interrupt:
                terminate_process(proc)
            raise
        except subprocess.TimeoutExpired:
            terminate_process(proc)
            raise SkillError(
                stage="execution",
                error_code="EXECUTION_FAILED",
                message="Nextflow execution timed out.",
                fix=timeout_fix,
                details={"timeout_seconds": timeout_seconds},
            )

    if exit_code != 0:
        raise SkillError(
            stage="execution",
            error_code="EXECUTION_FAILED",
            message="Nextflow execution failed.",
            fix=(
                "Inspect logs/stdout.txt and logs/stderr.txt, then correct the failing input or environment."
                + failure_hint(stdout_path, stderr_path, cwd / ".nextflow.log")
            ),
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


def terminate_process_tree(proc: subprocess.Popen, *, grace_seconds: int = 10) -> None:
    if sys.platform != "win32":
        try:
            pgid = os.getpgid(proc.pid)
            os.killpg(pgid, signal.SIGTERM)
            proc.wait(timeout=grace_seconds)
            return
        except (OSError, ProcessLookupError, subprocess.TimeoutExpired):
            pass
        try:
            pgid = os.getpgid(proc.pid)
            os.killpg(pgid, signal.SIGKILL)
            proc.wait()
            return
        except (OSError, ProcessLookupError, subprocess.TimeoutExpired):
            pass
    else:
        try:
            subprocess.run(
                ["taskkill", "/T", "/F", "/PID", str(proc.pid)],
                stdout=subprocess.DEVNULL,
                stderr=subprocess.DEVNULL,
                check=False,
            )
            return
        except (OSError, subprocess.SubprocessError):
            pass
    proc.kill()
    try:
        proc.wait()
    except subprocess.TimeoutExpired:
        pass
