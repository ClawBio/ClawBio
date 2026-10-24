"""Nextflow process mechanics; callers supply only pipeline-specific advice."""

from __future__ import annotations

import os
import re
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
# Why Nextflow stopped a run that never got going: .nextflow.log records an error(...) from
# the pipeline's own input checks as the session's abort cause; Nextflow and its launcher
# print their refusals (manifest.nextflowVersion, the Java check) as ERROR lines.
_SESSION_ABORT_CAUSE = "Session aborted -- Cause:"
# An abort before any session exists (e.g. a plugin requiring a newer Nextflow) is logged
# only as the exception Nextflow uses for refusals it explains to the user. After a
# session abort the same exception also reports that the run's HTML reports could not be
# written, which is not why it stopped, so it is read only when no session abort is.
_ABORT_EXCEPTION = "nextflow.exception.AbortOperationException:"
_ERROR_LINE_PREFIXES = ("ERROR ~ ", "ERROR: ")
# Lines that say a run failed without saying why, or that nf_schema_refusals reports already.
_UNINFORMATIVE_STOP_REASONS = (
    "pipeline failed. please refer to troubleshooting docs",
    "validation of pipeline parameters failed",
    "validation of samplesheet failed",
)
# stderr lines that are never the pipeline's account of what it refused.
_STDERR_NOISE = re.compile(r"^(?:WARN|NOTE|Picked up|Nextflow \S+ is available)", re.IGNORECASE)
_ANSI_ESCAPE = re.compile(r"\x1b\[[0-9;]*[A-Za-z]")
_STDERR_DETAIL_LINES = 3
# Nextflow's observers fail to write their HTML reports after an abort; never the reason.
_REPORT_NOT_WRITTEN = re.compile(r"file already exists: .* -- enable the '\w+\.overwrite' option", re.IGNORECASE)
# An input refusal is logged early, so the start of .nextflow.log is read as well as its
# tail (a long DEBUG log can push the abort cause out of the tail).
_LOG_HEAD_BYTES = 1 << 20


def read_log_head(path: Path, limit: int = _LOG_HEAD_BYTES) -> str:
    """The first ``limit`` bytes of a captured log, decoded; never raises."""
    try:
        with path.open("rb") as fh:
            return fh.read(limit).decode("utf-8", errors="replace")
    except OSError:
        return ""


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


def nextflow_stop_reasons(stdout: str, stderr: str, nextflow_log: str = "") -> list[str]:
    """Why Nextflow stopped the run, once each, followed by what the pipeline printed to
    stderr just before (e.g. sarek names the offending samplesheet entry there).

    Empty when nothing says why: stderr alone is never quoted.
    """
    reasons: list[str] = []

    def add(reason: str) -> None:
        reason = _ANSI_ESCAPE.sub("", reason).strip()
        if _REPORT_NOT_WRITTEN.search(reason):
            return
        if reason and reason not in reasons and not reason.lower().startswith(_UNINFORMATIVE_STOP_REASONS):
            reasons.append(reason)

    log_lines = nextflow_log.splitlines()
    for line in log_lines:
        if _SESSION_ABORT_CAUSE in line:
            add(line.split(_SESSION_ABORT_CAUSE, 1)[1])
    if not reasons:
        for line in log_lines:
            if _ABORT_EXCEPTION in line:
                add(line.split(_ABORT_EXCEPTION, 1)[1])
    for line in (stdout + "\n" + stderr).splitlines():
        line = _ANSI_ESCAPE.sub("", line).strip()
        for prefix in _ERROR_LINE_PREFIXES:
            if line.startswith(prefix):
                add(line[len(prefix):])
    if not reasons:
        return []
    detail = [
        line for line in (_ANSI_ESCAPE.sub("", raw).strip() for raw in stderr.splitlines())
        if line and not _STDERR_NOISE.match(line) and not line.startswith(_ERROR_LINE_PREFIXES)
    ]
    for line in detail[-_STDERR_DETAIL_LINES:]:
        add(line)
    return reasons


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
        parts.append(read_log_head(nextflow_log_path) + "\n" + read_log_tail(nextflow_log_path))
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
    # nf-schema's refusal is the reason Nextflow stopped; its entries need no second listing.
    stopped = [] if refused else nextflow_stop_reasons(parts[0], parts[1], parts[2] if len(parts) > 2 else "")
    if stopped:
        hints.append(
            " Nextflow stopped the run: " + "; ".join(reason.rstrip(".") for reason in stopped)
            + ". Correct what it names and re-run."
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
    # Nextflow starts its own session, so its pid is the group id. getpgid would
    # fail once Nextflow has exited, while the tasks it left in the group run on.
    try:
        os.killpg(proc.pid, signal.SIGTERM)
        proc.wait(timeout=grace_seconds)
        return
    except (OSError, subprocess.TimeoutExpired):
        pass
    try:
        os.killpg(proc.pid, signal.SIGKILL)
        proc.wait()
        return
    except OSError:
        pass
    proc.kill()
    proc.wait()
