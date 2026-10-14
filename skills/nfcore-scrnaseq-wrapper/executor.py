from __future__ import annotations

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
    from clawbio.common.nfcore import read_log_tail

    return read_log_tail(path, limit)


def _environment_failure_hints(
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


def execute_nextflow(
    command: list[str],
    *,
    cwd: Path,
    output_dir: Path,
    timeout_seconds: int | None,
) -> dict[str, object]:
    # timeout_seconds=None disables the wall-clock cap entirely (long HPC/cloud
    # runs whose walltime is enforced by the scheduler). Any positive int caps the
    # run and kills the process tree on expiry.
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
                error_code=ErrorCode.EXECUTION_FAILED,
                message=f"Could not launch Nextflow: {exc}",
                fix="Ensure nextflow is installed and available on PATH.",
                details={"command": command[0] if command else "", "error": str(exc)},
            ) from exc

        try:
            exit_code = proc.wait(timeout=timeout_seconds)
        except KeyboardInterrupt:
            _terminate_process_tree(proc)
            raise
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
            fix=(
                "Inspect logs/stdout.txt and logs/stderr.txt, then correct the failing input or environment."
                + _macos_tmp_failure_hint(output_dir)
                + _environment_failure_hints(stdout_path, stderr_path, cwd / ".nextflow.log")
                + _cellbender_failure_hint(stdout_path, stderr_path)
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
