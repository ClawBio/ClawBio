"""Tests for executor.py — log placement and process handling."""
from __future__ import annotations

from pathlib import Path
import os
import shlex
import signal
import sys
import threading
import time

import pytest

from clawbio.common.nfcore import execution as shared_execution
from errors import SkillError
from executor import execute_nextflow


def test_logs_written_under_root_logs_dir(tmp_path: Path):
    """Run logs must live at ``<output_dir>/logs/`` — the same location as the
    nfcore-rnaseq and nfcore-scrnaseq wrappers, so all three share one layout.
    """
    output_dir = tmp_path / "out"
    output_dir.mkdir()

    result = execute_nextflow(
        ["sh", "-c", "echo hello; echo oops 1>&2"],
        cwd=output_dir,
        output_dir=output_dir,
        timeout_seconds=30,
    )

    logs = output_dir / "logs"
    assert (logs / "stdout.txt").exists()
    assert (logs / "stderr.txt").exists()
    # Logs no longer live inside the reproducibility bundle.
    assert not (output_dir / "reproducibility" / "logs").exists()
    assert result["exit_code"] == 0
    assert (logs / "stdout.txt").read_text().strip() == "hello"
    assert (logs / "stderr.txt").read_text().strip() == "oops"
    # Reported paths point at the root-level logs/ directory.
    assert "logs/stdout.txt" in Path(str(result["stdout_path"])).as_posix()
    assert "reproducibility/logs" not in Path(str(result["stdout_path"])).as_posix()


def test_macos_tmp_failure_hint_appended_on_failure(tmp_path, monkeypatch):
    """On macOS, a failed run with --output under /tmp must append the actionable
    Colima / 'No such file or directory' hint to the fix (parity with
    nfcore-scrnaseq / nfcore-rnaseq)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr(shared_execution, "is_under_tmp", lambda _p: True)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["sh", "-c", "exit 1"], cwd=output_dir, output_dir=output_dir, timeout_seconds=30)
    assert "No such file or directory" in exc.value.fix
    assert "home" in exc.value.fix.lower()


def test_macos_tmp_failure_hint_absent_when_not_under_tmp(tmp_path, monkeypatch):
    """The macOS /tmp hint must not appear for a non-/tmp output directory."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr(shared_execution, "is_under_tmp", lambda _p: False)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["sh", "-c", "exit 1"], cwd=output_dir, output_dir=output_dir, timeout_seconds=30)
    assert "No such file or directory" not in exc.value.fix


def test_memory_limit_failure_hint(tmp_path):
    """A run that fails because a process exceeds host memory must append an
    actionable hint pointing at Nextflow resourceLimits (Issue: nf-core defaults
    request more memory than a small host provides)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    sig = "Process requirement exceeds available memory -- req: 72 GB; avail: 62.8 GB"
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", f"echo {shlex.quote(sig)} 1>&2; exit 1"],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "resourceLimits" in exc.value.fix


def test_network_unreachable_failure_hint(tmp_path):
    """A run that fails with an unreachable network must append the IPv6/NAT64
    NXF_OPTS hint (Issue: IPv6-only / NAT64 hosts where the JVM prefers IPv4)."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    sig = "Network is unreachable (connect failed)"
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", f"echo {shlex.quote(sig)} 1>&2; exit 1"],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "preferIPv6Addresses" in exc.value.fix


def test_network_hint_detected_from_nextflow_log(tmp_path):
    """When Nextflow fails during config parsing, the real cause ('Network is
    unreachable') is frequently only in `.nextflow.log` (written to the launch cwd),
    not in stdout/stderr. The environment hint must scan `.nextflow.log` too, or the
    IPv6/NAT64 `NXF_OPTS='-Djava.net.preferIPv6Addresses=true'` hint never shows on the
    hosts that need it."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    cmd = (
        "printf 'Caused by: java.net.SocketException: Network is unreachable' > .nextflow.log; "
        "echo 'Unable to parse config file' 1>&2; exit 1"
    )
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", cmd],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "preferIPv6Addresses" in exc.value.fix


def test_no_environment_hint_without_signature(tmp_path):
    """A generic failure with no known signature must not append the env hints."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", "echo boom 1>&2; exit 1"],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "resourceLimits" not in exc.value.fix
    assert "preferIPv6Addresses" not in exc.value.fix


def test_config_parse_failure_hint(tmp_path):
    """A run that fails because Nextflow could not fetch/parse the remote nf-core
    config (the `includeConfig … ? <url> : '/dev/null'` line) must point at
    NXF_OFFLINE for a fully-local run."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    sig = "Unable to parse config file: '/root/.nextflow/assets/nf-core/rnaseq/nextflow.config'"
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", f"echo {shlex.quote(sig)} 1>&2; exit 1"],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "NXF_OFFLINE" in exc.value.fix


def test_no_config_parse_hint_without_signal(tmp_path):
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", "echo plain-failure 1>&2; exit 1"],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "NXF_OFFLINE" not in exc.value.fix


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX process groups and SIGINT")
def test_ctrl_c_terminates_nextflow_instead_of_orphaning_it(tmp_path):
    """Nextflow runs in its own session, so the terminal's SIGINT never reaches it.
    Ctrl+C must therefore stop the child explicitly, or Nextflow keeps running (and
    holds the .nextflow lock) after the wrapper has exited 130."""
    pid_file = tmp_path / "child.pid"

    def interrupt_once_child_runs() -> None:
        deadline = time.monotonic() + 10
        while not pid_file.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        if pid_file.exists():
            os.kill(os.getpid(), signal.SIGINT)

    threading.Thread(target=interrupt_once_child_runs, daemon=True).start()
    with pytest.raises(KeyboardInterrupt):
        execute_nextflow(
            ["sh", "-c", f"echo $$ > {shlex.quote(str(pid_file))}; exec sleep 60"],
            cwd=tmp_path, output_dir=tmp_path, timeout_seconds=None,
        )
    with pytest.raises(ProcessLookupError):
        os.kill(int(pid_file.read_text()), 0)
