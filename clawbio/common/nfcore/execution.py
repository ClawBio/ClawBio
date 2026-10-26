"""Nextflow process mechanics; callers supply only pipeline-specific advice."""

from __future__ import annotations

import os
import signal
import subprocess
import sys
import threading
import time
from collections.abc import Callable, Mapping
from contextlib import contextmanager
from pathlib import Path

from .errors import SkillError
from .paths import is_under_tmp

_TERMINATION_GRACE_SECONDS = 10
# How often a wait checks for a recorded SIGINT/SIGTERM.
_STOP_POLL_SECONDS = 0.2

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
    """Each "--name (value): reason" nf-schema listed under its refusal header, once.

    A refused samplesheet's detail is on the "-> Entry N: ..." lines nf-schema
    writes under its "* --input" item; they are kept with that item.
    """
    refused: list[str] = []
    lines = text.splitlines()
    for index, line in enumerate(lines):
        if _NF_SCHEMA_REFUSAL_HEADER not in line.lower():
            continue
        item: list[str] | None = None
        items: list[str] = []
        for entry in lines[index + 1:]:
            entry = entry.strip()
            if entry.startswith("* --"):
                item = [entry[2:]]
                items.append(item)
            elif entry.startswith("->") and item is not None:
                item.append(entry[2:].strip())
            elif entry or item is not None:
                break
        for parts in items:
            joined = " ".join(parts[:1]) + "".join(f" {part};" for part in parts[1:])
            joined = joined.rstrip(";")
            if joined not in refused:
                refused.append(joined)
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
            " nf-schema refused these against the pinned pipeline's schemas "
            "(nextflow_schema.json, and assets/schema_input.json for the samplesheet): "
            + "; ".join(refused) + ". Correct them in your flags, --params-file or "
            "samplesheet and re-run."
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


@contextmanager
def _stops_recorded():
    """While Nextflow runs, record SIGTERM and SIGINT instead of raising on them.

    Either would otherwise end Python and leave Nextflow running in its own session.
    The handlers must not raise: an exception raised between Popen.wait() taking and
    releasing its non-reentrant wait lock leaves the lock held, and every later wait
    on Nextflow, including the cleanup's, then blocks forever. SIGINT is recorded only
    while Python's default handler is installed. Only the main thread can install
    handlers; elsewhere the list stays empty.
    """
    received: list[int] = []
    if threading.current_thread() is not threading.main_thread():
        yield received
        return

    def record(signum, frame) -> None:
        received.append(signum)

    stops = [signal.SIGTERM]
    if signal.getsignal(signal.SIGINT) is signal.default_int_handler:
        stops.append(signal.SIGINT)
    previous = {signum: signal.signal(signum, record) for signum in stops}
    try:
        yield received
    finally:
        for signum, handler in previous.items():
            signal.signal(signum, handler)


def _wait_unless_stopped(
    proc: subprocess.Popen, timeout_seconds: int | None, received: list[int]
) -> int | None:
    """Nextflow's exit code, or None once a stop is recorded; waits in short slices."""
    deadline = None if timeout_seconds is None else time.monotonic() + timeout_seconds
    while not received:
        remaining = None if deadline is None else deadline - time.monotonic()
        if remaining is not None and remaining <= 0:
            raise subprocess.TimeoutExpired(proc.args, timeout_seconds)
        slice_seconds = _STOP_POLL_SECONDS if remaining is None else min(_STOP_POLL_SECONDS, remaining)
        try:
            return proc.wait(timeout=slice_seconds)
        except subprocess.TimeoutExpired:
            continue
    return None


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
            with _stops_recorded() as received:
                exit_code = _wait_unless_stopped(proc, timeout_seconds, received)
        except KeyboardInterrupt:
            # SIGINT under a caller's own handler, which was left in place.
            terminate_process_tree(proc)
            raise
        except subprocess.TimeoutExpired:
            terminate_process_tree(proc)
            raise SkillError(
                stage="execution",
                error_code="EXECUTION_FAILED",
                message="Nextflow execution timed out.",
                fix=timeout_fix,
                details={"timeout_seconds": timeout_seconds},
            )
        if received:
            terminate_process_tree(proc)
            if received[0] == signal.SIGINT:
                raise KeyboardInterrupt
            raise SystemExit(128 + received[0])

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
    proc: subprocess.Popen, *, grace_seconds: float = _TERMINATION_GRACE_SECONDS
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
    # Nextflow starts its own session, so its pid is the group id. getpgid would
    # fail once Nextflow has exited, while the tasks it left in the group run on.
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        deadline = time.monotonic() + grace_seconds
        while True:
            proc.poll()  # Reap the leader, but keep waiting for its task group.
            try:
                os.killpg(proc.pid, 0)
            except ProcessLookupError:
                proc.wait()
                return
            remaining = deadline - time.monotonic()
            if remaining <= 0:
                break
            time.sleep(min(0.05, remaining))
    except (OSError, subprocess.TimeoutExpired):
        # A failed graceful stop still requires the forceful group cleanup below.
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        return
    except OSError:
        # The group may be gone; fall back to stopping and reaping its leader.
        pass
    proc.kill()
    proc.wait()
