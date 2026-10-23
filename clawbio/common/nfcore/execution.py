"""Nextflow process mechanics; callers supply only pipeline-specific advice."""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import threading
from collections.abc import Callable, Mapping
from contextlib import contextmanager
from pathlib import Path

from .errors import SkillError
from .paths import is_under_tmp

_TERMINATION_GRACE_SECONDS = 10

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
# nf-schema's header above the list of refused parameters, one "* --name (value): reason" each.
_NF_SCHEMA_REFUSAL_HEADER = "the following invalid input values have been detected"


def read_log_tail(path: Path, limit: int = 65536) -> str:
    """The last ``limit`` bytes of a captured log, decoded; never raises.

    Seeks rather than reading the whole file: a long run's Nextflow logs can reach
    gigabytes, and the failure hints only scan their end.
    """
    try:
        with path.open("rb") as fh:
            fh.seek(0, os.SEEK_END)
            fh.seek(max(0, fh.tell() - limit))
            return fh.read(limit).decode("utf-8", errors="replace")
    except OSError:
        return ""


def macos_tmp_failure_hint(output_dir: Path) -> str:
    """Extra fix hint when a run fails with --output under /tmp on macOS.

    Colima/Docker on macOS does not share /tmp into the VM, a frequent cause of
    '.command.run: No such file or directory' mid-run. Preflight already WARNs about
    this; the actionable hint is repeated on failure so the cause is obvious. Empty
    string on other platforms or when --output is not under /tmp.
    """
    if sys.platform != "darwin" or not is_under_tmp(output_dir):
        return ""
    return (
        " On macOS, --output is under /tmp, which Docker/Colima does not share into its VM; "
        "this commonly surfaces as '.command.run: No such file or directory'. "
        "Move --output to a path under your HOME directory and re-run."
    )


def nf_schema_refusals(text: str) -> list[str]:
    """Each "--name (value): reason" nf-schema listed under its refusal header, once."""
    refused: list[str] = []
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if _NF_SCHEMA_REFUSAL_HEADER not in line.lower():
            continue
        listing = False
        for entry in lines[index + 1:]:
            entry = entry.strip()
            if entry.startswith("* --"):
                listing = True
                if entry[2:] not in refused:
                    refused.append(entry[2:])
            elif entry or listing:
                break
    return refused


def environment_failure_hints(
    stdout_path: Path, stderr_path: Path, nextflow_log_path: Path | None = None
) -> str:
    """Append actionable hints when the captured logs show a known environment
    failure — the host being smaller than a process request, or an unreachable
    network. Diagnosed from the actual Nextflow error text (never predicting
    pipeline resource requirements), so no thresholds are fabricated. Empty string
    when no known signature is present.

    `.nextflow.log` is scanned in addition to stdout/stderr: when Nextflow fails while
    parsing the config (before the pipeline starts) the top-level stderr may show only
    "Unable to parse config file", while the underlying cause ("Network is unreachable")
    is recorded only in `.nextflow.log`. Without it the IPv6/NAT64 hint never shows.
    """
    parts = [read_log_tail(stdout_path), read_log_tail(stderr_path)]
    if nextflow_log_path is not None:
        parts.append(read_log_tail(nextflow_log_path))
    blob = "\n".join(parts).lower()
    hints: list[str] = []
    refused = nf_schema_refusals("\n".join(parts))
    if refused:
        hints.append(
            " nf-schema refused these parameters against the pinned pipeline's "
            "nextflow_schema.json: " + "; ".join(refused) + ". Correct them in your flags "
            "or --params-file and re-run."
        )
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


class _Terminated(BaseException):
    """SIGTERM arrived while Nextflow ran (BaseException: no broad handler eats it)."""

    def __init__(self, signum: int) -> None:
        super().__init__(signum)
        self.signum = signum


def _raise_terminated(signum, frame) -> None:
    raise _Terminated(signum)


@contextmanager
def _sigterm_raises():
    """While Nextflow runs, make SIGTERM raise _Terminated instead of ending Python.

    SIGTERM (a plain `kill`, an agent or a scheduler stopping the wrapper) would
    otherwise leave Nextflow running in its own session. The handler only raises:
    stopping Nextflow must happen after the stack unwinds out of proc.wait(), which
    holds Popen's non-reentrant wait lock. Only the main thread can install it.
    """
    if threading.current_thread() is not threading.main_thread():
        yield
        return
    previous = signal.signal(signal.SIGTERM, _raise_terminated)
    try:
        yield
    finally:
        signal.signal(signal.SIGTERM, previous)


def execute_nextflow(
    command: list[str],
    *,
    cwd: Path,
    output_dir: Path,
    timeout_seconds: int | None,
    timeout_fix: str,
    pipeline_failure_hint: Callable[[Path, Path], str] | None = None,
    env: Mapping[str, str] | None = None,
) -> dict[str, object]:
    """Run Nextflow with logs on disk and its whole process group under our control.

    Nextflow starts in its own session so a timeout can signal the group, which also
    means neither the terminal's Ctrl+C nor a SIGTERM to the wrapper reaches it: both
    must terminate the group explicitly, or Nextflow outlives the wrapper still
    holding the .nextflow lock. timeout_seconds=None disables the wall-clock cap
    (long HPC/cloud runs whose walltime the scheduler enforces). ``env`` is what the
    pipeline release needs on top of the caller's environment, and wins over it.
    """
    logs_dir = output_dir / "logs"
    logs_dir.mkdir(parents=True, exist_ok=True)
    stdout_path = logs_dir / "stdout.txt"
    stderr_path = logs_dir / "stderr.txt"

    popen_kwargs: dict = {}
    if env:
        popen_kwargs["env"] = {**os.environ, **env}
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
            with _sigterm_raises():
                exit_code = proc.wait(timeout=timeout_seconds)
        except KeyboardInterrupt:
            terminate_process_tree(proc)
            raise
        except _Terminated as stop:
            terminate_process_tree(proc)
            raise SystemExit(128 + stop.signum) from None
        except subprocess.TimeoutExpired:
            terminate_process_tree(proc)
            raise SkillError(
                stage="execution",
                error_code="EXECUTION_FAILED",
                message="Nextflow execution timed out.",
                fix=timeout_fix,
                details={"timeout_seconds": timeout_seconds},
            )

    if exit_code != 0:
        fix = (
            "Inspect logs/stdout.txt and logs/stderr.txt, then correct the failing input or environment."
            + macos_tmp_failure_hint(output_dir)
            + environment_failure_hints(stdout_path, stderr_path, cwd / ".nextflow.log")
        )
        if pipeline_failure_hint is not None:
            fix += pipeline_failure_hint(stdout_path, stderr_path)
        raise SkillError(
            stage="execution",
            error_code="EXECUTION_FAILED",
            message="Nextflow execution failed.",
            fix=fix,
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


def terminate_process_tree(
    proc: subprocess.Popen, *, grace_seconds: int = _TERMINATION_GRACE_SECONDS
) -> None:
    """Stop Nextflow and the task processes it spawned, never only the parent."""
    if sys.platform == "win32":
        # Windows has no os.killpg; ``taskkill /T`` terminates the whole child tree.
        # Bounded so a hung taskkill cannot block the wrapper, and followed by
        # proc.kill() in case taskkill is unavailable. Native Windows is not
        # officially supported (preflight warns); WSL2 reports "linux".
        try:
            subprocess.run(
                ["taskkill", "/F", "/T", "/PID", str(proc.pid)],
                capture_output=True,
                timeout=grace_seconds,
            )
        except (OSError, subprocess.SubprocessError):
            # taskkill missing, failing or hung: proc.kill() below still stops the child.
            pass
        proc.kill()
        proc.wait()
        return
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGTERM)
        proc.wait(timeout=grace_seconds)
        return
    except (OSError, subprocess.TimeoutExpired):
        pass
    try:
        os.killpg(os.getpgid(proc.pid), signal.SIGKILL)
        proc.wait()
        return
    except OSError:
        pass
    proc.kill()
    proc.wait()
