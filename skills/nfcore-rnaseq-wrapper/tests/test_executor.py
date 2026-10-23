from __future__ import annotations

from pathlib import Path
import os
import shlex
import signal
import subprocess
import sys
import threading
import time

import pytest

from clawbio.common.nfcore import execution as shared_execution

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))


def _purge_foreign_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and _SKILL_DIR not in module_file.parents and module_file != _SKILL_DIR / f"{name}.py":
            sys.modules.pop(name, None)


def _purge_local_modules(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if module is not None and (module_file == _SKILL_DIR / f"{name}.py" or _SKILL_DIR in module_file.parents):
            sys.modules.pop(name, None)


def _remove_skill_dir_from_sys_path() -> None:
    while str(_SKILL_DIR) in sys.path:
        sys.path.remove(str(_SKILL_DIR))


_purge_foreign_modules("errors", "executor")

from errors import SkillError
from executor import execute_nextflow

_purge_local_modules("errors", "executor")
_remove_skill_dir_from_sys_path()


class _MockProc:
    def __init__(self, returncode: int = 0, timeout_on_first_wait: bool = False):
        self._rc = returncode
        self._timeout_on_first = timeout_on_first_wait
        self._wait_calls = 0
        self.killed = False
        self.pid = 1234

    def wait(self, timeout=None):
        self._wait_calls += 1
        if self._timeout_on_first and self._wait_calls == 1:
            raise subprocess.TimeoutExpired(cmd=["nextflow"], timeout=timeout or 1)
        return self._rc

    def kill(self):
        self.killed = True


def test_execute_nextflow_creates_log_files_on_success(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    result = execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert result["exit_code"] == 0
    assert Path(result["stdout_path"]).name == "stdout.txt"
    assert Path(result["stderr_path"]).name == "stderr.txt"
    assert (tmp_path / "logs" / "stdout.txt").exists()
    assert (tmp_path / "logs" / "stderr.txt").exists()


def test_execute_nextflow_raises_on_nonzero_exit_and_keeps_logs(tmp_path, monkeypatch):
    proc = _MockProc(returncode=1)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert exc.value.details["exit_code"] == 1
    assert (tmp_path / "logs" / "stdout.txt").exists()
    assert (tmp_path / "logs" / "stderr.txt").exists()


def test_execute_nextflow_failure_appends_macos_tmp_hint(tmp_path, monkeypatch):
    """On macOS, a failed run with --output under /tmp must append the actionable
    Colima / 'No such file or directory' hint to the fix (parity with
    nfcore-scrnaseq), so the silent cause of the failure is surfaced."""
    proc = _MockProc(returncode=1)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr(shared_execution, "is_under_tmp", lambda _p: True)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert "No such file or directory" in exc.value.fix
    assert "home" in exc.value.fix.lower()


def test_execute_nextflow_failure_no_tmp_hint_when_not_under_tmp(tmp_path, monkeypatch):
    """The macOS /tmp hint must not appear for a non-/tmp output directory."""
    proc = _MockProc(returncode=1)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(sys, "platform", "darwin")
    monkeypatch.setattr(shared_execution, "is_under_tmp", lambda _p: False)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert "No such file or directory" not in exc.value.fix


def test_execute_nextflow_starts_new_session_on_posix(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0)
    popen_kwargs = {}

    def fake_popen(*args, **kwargs):
        popen_kwargs.update(kwargs)
        return proc

    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(subprocess, "Popen", fake_popen)
    execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert popen_kwargs["start_new_session"] is True


def test_execute_nextflow_terminates_process_group_on_posix_timeout(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0, timeout_on_first_wait=True)
    # Nextflow starts its own session, so its pid is the process-group id.
    pgid = proc.pid
    signals_sent: list[tuple[int, int]] = []

    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(os, "killpg", lambda sent_pgid, sig: signals_sent.append((sent_pgid, sig)))

    with pytest.raises(SkillError):
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=1)

    assert signals_sent
    assert all(sent_pgid == pgid for sent_pgid, _ in signals_sent)
    assert proc.killed is False


@pytest.mark.skipif(sys.platform == "win32", reason="POSIX process groups and SIGINT")
def test_ctrl_c_terminates_nextflow_instead_of_orphaning_it(tmp_path):
    """Nextflow runs in its own session, so the terminal's SIGINT never reaches it.
    Ctrl+C must therefore stop the child explicitly, or Nextflow keeps running (and
    holds the .nextflow lock) after the wrapper has exited 130."""
    pid_file = tmp_path / "child.pid"
    main_thread = threading.main_thread().ident

    def interrupt_once_child_runs() -> None:
        deadline = time.monotonic() + 30
        while not pid_file.exists() and time.monotonic() < deadline:
            time.sleep(0.01)
        if pid_file.exists():
            # What the terminal does on Ctrl+C, aimed at the thread blocked in
            # waitpid so the kernel cannot pick this helper thread instead.
            signal.pthread_kill(main_thread, signal.SIGINT)

    # Hermetic: other tests in a full run may leave a different SIGINT handler.
    previous_handler = signal.signal(signal.SIGINT, signal.default_int_handler)
    try:
        threading.Thread(target=interrupt_once_child_runs, daemon=True).start()
        with pytest.raises(KeyboardInterrupt):
            execute_nextflow(
                ["sh", "-c", f"echo $$ > {shlex.quote(str(pid_file))}; exec sleep 60"],
                cwd=tmp_path, output_dir=tmp_path, timeout_seconds=None,
            )
    finally:
        signal.signal(signal.SIGINT, previous_handler)
    with pytest.raises(ProcessLookupError):
        os.kill(int(pid_file.read_text()), 0)


def test_execute_nextflow_cleans_up_empty_logs_on_popen_failure(tmp_path, monkeypatch):
    def raise_os_error(*args, **kwargs):
        raise OSError("nextflow: command not found")

    monkeypatch.setattr(subprocess, "Popen", raise_os_error)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert not (tmp_path / "logs" / "stdout.txt").exists()
    assert not (tmp_path / "logs" / "stderr.txt").exists()


def test_execute_nextflow_timeout_error_reports_timeout_seconds(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0, timeout_on_first_wait=True)
    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: subprocess.CompletedProcess(a[0], 0))
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=7)
    assert exc.value.details["timeout_seconds"] == 7


def test_timeout_error_fix_mentions_leftover_containers(tmp_path, monkeypatch):
    """Killing Nextflow does not stop containers spawned by the Docker daemon (they
    are not in Nextflow's process group). The timeout error must tell the user to
    check for and clean up any leftover containers (audit F7)."""
    proc = _MockProc(returncode=0, timeout_on_first_wait=True)
    monkeypatch.setattr(sys, "platform", "win32")
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(subprocess, "run", lambda *a, **kw: subprocess.CompletedProcess(a[0], 0))
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=1)
    assert "container" in exc.value.fix.lower()


def test_terminate_process_tree_falls_back_to_sigkill_when_sigterm_times_out(monkeypatch):
    """terminate_process_tree must send SIGKILL when SIGTERM's wait() times out."""
    import signal as signal_module

    signals_sent: list[int] = []
    wait_calls = [0]

    class _SlowProc:
        def __init__(self):
            self.pid = 9999
            self.killed = False

        def wait(self, timeout=None):
            wait_calls[0] += 1
            # First wait (grace period after SIGTERM) times out.
            # Second wait (after SIGKILL) completes immediately.
            if wait_calls[0] == 1:
                raise subprocess.TimeoutExpired(cmd=["nextflow"], timeout=timeout or 1)
            return -signal_module.SIGKILL

        def kill(self):
            self.killed = True

    proc = _SlowProc()

    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(
        os,
        "killpg",
        lambda sent_pgid, sig: signals_sent.append(sig),
    )

    shared_execution.terminate_process_tree(proc)

    # SIGTERM must be sent first, then SIGKILL as the fallback — order matters.
    assert signals_sent == [signal_module.SIGTERM, signal_module.SIGKILL], (
        f"Expected [SIGTERM, SIGKILL], got {signals_sent}"
    )
    # proc.kill() (last-resort fallback) must NOT have been called —
    # killpg(SIGKILL) should have succeeded.
    assert proc.killed is False


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
    # stdout/stderr carry only the top-level config-parse message; the network cause is
    # in .nextflow.log (which Nextflow writes to its cwd).
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
