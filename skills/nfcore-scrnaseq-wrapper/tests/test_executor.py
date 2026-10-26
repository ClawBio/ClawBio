from __future__ import annotations

import shlex
import subprocess
import sys
import time
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from clawbio.common.nfcore import execution as shared_execution
from errors import SkillError
from executor import execute_nextflow


class _MockProc:
    """Minimal Popen replacement that avoids spawning a real process."""

    def __init__(
        self,
        returncode: int = 0,
        runs_until_killed: bool = False,
        keyboard_interrupt_on_first_wait: bool = False,
    ):
        self._rc = returncode
        self._runs_until_killed = runs_until_killed
        self._keyboard_interrupt_on_first = keyboard_interrupt_on_first_wait
        self._wait_calls = 0
        self.killed = False
        self.pid = 1234
        self.args = ["nextflow", "run", "test"]

    def wait(self, timeout=None):
        self._wait_calls += 1
        if self._runs_until_killed and not self.killed:
            time.sleep(timeout)  # A real wait blocks this long before giving up.
            raise subprocess.TimeoutExpired(cmd=["nextflow"], timeout=timeout)
        if self._keyboard_interrupt_on_first and self._wait_calls == 1:
            raise KeyboardInterrupt
        return self._rc

    def kill(self):
        self.killed = True


def test_execute_nextflow_creates_log_files_on_success(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    result = execute_nextflow(
        ["nextflow", "run", "test"],
        cwd=tmp_path,
        output_dir=tmp_path,
        timeout_seconds=60,
    )
    assert result["exit_code"] == 0
    assert (tmp_path / "logs" / "stdout.txt").exists()
    assert (tmp_path / "logs" / "stderr.txt").exists()


def test_execute_nextflow_raises_on_nonzero_exit(tmp_path, monkeypatch):
    proc = _MockProc(returncode=1)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["nextflow", "run", "test"],
            cwd=tmp_path,
            output_dir=tmp_path,
            timeout_seconds=60,
        )
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert exc.value.details["exit_code"] == 1


def _record_taskkill(monkeypatch) -> list[list[str]]:
    """Stub subprocess.run so the win32 ``taskkill /T`` tree-kill is recorded
    instead of spawning a real process (audit F-7)."""
    calls: list[list[str]] = []

    def fake_run(args, *a, **kw):
        calls.append(list(args))
        return subprocess.CompletedProcess(args, 0)

    monkeypatch.setattr(subprocess, "run", fake_run)
    return calls


def test_execute_nextflow_kills_process_on_timeout(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0, runs_until_killed=True)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(sys, "platform", "win32")
    taskkill_calls = _record_taskkill(monkeypatch)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["nextflow", "run", "test"],
            cwd=tmp_path,
            output_dir=tmp_path,
            timeout_seconds=1,
        )
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert proc.killed, "Process must be killed on timeout"
    # The whole child tree is terminated via taskkill /T, not just the parent.
    assert ["taskkill", "/F", "/T", "/PID", "1234"] in taskkill_calls


def test_execute_nextflow_kills_process_on_keyboard_interrupt(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0, keyboard_interrupt_on_first_wait=True)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(sys, "platform", "win32")
    taskkill_calls = _record_taskkill(monkeypatch)

    with pytest.raises(KeyboardInterrupt):
        execute_nextflow(
            ["nextflow", "run", "test"],
            cwd=tmp_path,
            output_dir=tmp_path,
            timeout_seconds=60,
        )

    assert proc.killed, "Process must be killed when the wrapper is interrupted"
    assert ["taskkill", "/F", "/T", "/PID", "1234"] in taskkill_calls


def test_execute_nextflow_result_contains_log_paths(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    result = execute_nextflow(
        ["nextflow", "run", "test"],
        cwd=tmp_path,
        output_dir=tmp_path,
        timeout_seconds=60,
    )
    assert "stdout_path" in result
    assert "stderr_path" in result
    assert Path(result["stdout_path"]).name == "stdout.txt"
    assert Path(result["stderr_path"]).name == "stderr.txt"


def test_execute_nextflow_cleans_up_empty_logs_on_popen_failure(tmp_path, monkeypatch):
    """Empty log files must be removed when Popen itself fails."""

    def raise_os_error(*args, **kwargs):
        raise OSError("nextflow: command not found")

    monkeypatch.setattr(subprocess, "Popen", raise_os_error)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["nextflow", "run", "test"],
            cwd=tmp_path,
            output_dir=tmp_path,
            timeout_seconds=60,
        )
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert not (tmp_path / "logs" / "stdout.txt").exists()
    assert not (tmp_path / "logs" / "stderr.txt").exists()


def test_macos_tmp_failure_hint_added_under_tmp(monkeypatch):
    """On macOS, an EXECUTION_FAILED for an --output under /tmp must include the
    actionable Colima/Docker hint (audit Hallazgo 2)."""
    monkeypatch.setattr(sys, "platform", "darwin")
    hint = shared_execution.macos_tmp_failure_hint(Path("/tmp/scrnaseq_run"))
    assert "/tmp" in hint and "HOME" in hint


def test_macos_tmp_failure_hint_empty_under_home(monkeypatch):
    monkeypatch.setattr(sys, "platform", "darwin")
    assert shared_execution.macos_tmp_failure_hint(Path.home() / "scrnaseq_run") == ""


def test_macos_tmp_failure_hint_empty_on_linux(monkeypatch):
    monkeypatch.setattr(sys, "platform", "linux")
    assert shared_execution.macos_tmp_failure_hint(Path("/tmp/scrnaseq_run")) == ""


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


def test_cellbender_failure_hint(tmp_path):
    """When CellBender background removal is the failing process (it errors on
    very small/test datasets), the fix must point at --skip-cellbender."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    err = (
        "Error executing process > "
        "'NFCORE_SCRNASEQ:SCRNASEQ:H5AD_REMOVEBACKGROUND_BARCODES_CELLBENDER_ANNDATA:"
        "CELLBENDER_REMOVEBACKGROUND (Sample_X)'\n"
        "IndexError: index -100 is out of bounds for axis 0 with size 88"
    )
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", f"echo {shlex.quote(err)} 1>&2; exit 1"],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "--skip-cellbender" in exc.value.fix


def test_no_cellbender_hint_without_signal(tmp_path):
    """A generic failure must not mention CellBender."""
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    with pytest.raises(SkillError) as exc:
        execute_nextflow(
            ["sh", "-c", "echo boom 1>&2; exit 1"],
            cwd=output_dir, output_dir=output_dir, timeout_seconds=30,
        )
    assert "--skip-cellbender" not in exc.value.fix


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
