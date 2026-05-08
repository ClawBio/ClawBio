from __future__ import annotations

import subprocess
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from errors import SkillError
from executor import execute_nextflow


class _MockProc:
    """Minimal Popen replacement that avoids spawning a real process."""

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
    assert (tmp_path / "logs" / "stdout.txt").exists()
    assert (tmp_path / "logs" / "stderr.txt").exists()


def test_execute_nextflow_raises_on_nonzero_exit(tmp_path, monkeypatch):
    proc = _MockProc(returncode=1)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert exc.value.details["exit_code"] == 1


def test_execute_nextflow_log_files_created_on_failure(tmp_path, monkeypatch):
    proc = _MockProc(returncode=1)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    with pytest.raises(SkillError):
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert (tmp_path / "logs" / "stdout.txt").exists()
    assert (tmp_path / "logs" / "stderr.txt").exists()


def test_execute_nextflow_kills_process_on_timeout(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0, timeout_on_first_wait=True)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr(sys, "platform", "win32")
    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=1)
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert proc.killed, "Process must be killed on timeout"


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


def test_execute_nextflow_terminates_process_group_on_timeout(tmp_path, monkeypatch):
    """killpg must receive the process GROUP id (pgid), not the process id (pid)."""
    proc = _MockProc(returncode=0, timeout_on_first_wait=True)
    proc.pid = 1234
    _PGID = 5678  # intentionally different from pid

    signals_sent: list[tuple[int, int]] = []

    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr("executor.os.getpgid", lambda pid: _PGID)
    monkeypatch.setattr("executor.os.killpg", lambda pgid, sig: signals_sent.append((pgid, sig)))

    with pytest.raises(SkillError) as exc:
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=1)

    assert exc.value.error_code == "EXECUTION_FAILED"
    assert signals_sent, "Timeout should signal the process group on POSIX"
    assert all(pgid == _PGID for pgid, _ in signals_sent), \
        f"killpg must receive pgid={_PGID}, not pid={proc.pid}"
    assert not proc.killed


def test_execute_nextflow_timeout_falls_back_to_process_kill(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0, timeout_on_first_wait=True)
    proc.pid = 1234
    monkeypatch.setattr(sys, "platform", "linux")
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    monkeypatch.setattr("executor.os.killpg", lambda *a, **kw: (_ for _ in ()).throw(OSError("no group")))
    with pytest.raises(SkillError):
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=1)
    assert proc.killed


def test_execute_nextflow_result_contains_log_paths(tmp_path, monkeypatch):
    proc = _MockProc(returncode=0)
    monkeypatch.setattr(subprocess, "Popen", lambda *a, **kw: proc)
    result = execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
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
        execute_nextflow(["nextflow", "run", "test"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=60)
    assert exc.value.error_code == "EXECUTION_FAILED"
    assert not (tmp_path / "logs" / "stdout.txt").exists()
    assert not (tmp_path / "logs" / "stderr.txt").exists()
