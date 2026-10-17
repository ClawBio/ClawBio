"""The root launcher must let nf-core wrappers stop their detached workflows."""
from __future__ import annotations

import json
import os
import signal
import subprocess
import sys
import textwrap
import time
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[2]
ALIASES = ("rnaseq-pipeline", "sarek-pipeline", "scrnaseq-pipeline")
pytestmark = pytest.mark.skipif(sys.platform == "win32", reason="POSIX process signals")


def _alive(pid: int) -> bool:
    try:
        os.kill(pid, 0)
    except ProcessLookupError:
        return False
    return True


def _stop_fixture(directory: Path, process: subprocess.Popen) -> None:
    """Never leave a failing regression's own processes behind."""
    for filename in ("wrapper.pid", "native.pid"):
        path = directory / filename
        if path.exists():
            try:
                os.kill(int(path.read_text()), signal.SIGTERM)
            except ProcessLookupError:
                pass  # A successful launcher has already reaped this process.
    if process.poll() is None:
        process.kill()
    process.wait(timeout=10)


@pytest.mark.parametrize("alias", ALIASES)
@pytest.mark.parametrize("stop", ("SIGTERM", "SIGINT", "SIGINT-twice", "SIGINT-caught", "timeout"))
def test_root_pipeline_stop_cleans_detached_workflow(alias, stop, tmp_path):
    """Use the real shared executor with a small scientific-process stand-in.

    Its delayed SIGTERM makes a premature parent SIGKILL observable. No Nextflow
    installation is needed for CI; the server reproduction uses real Nextflow.
    """
    native = tmp_path / "native.py"
    native.write_text(textwrap.dedent("""
        import signal, sys, time
        from pathlib import Path
        directory = Path(__file__).parent
        def finish(signum, frame):
            time.sleep(0.35)
            raise SystemExit(0)
        signal.signal(signal.SIGTERM, finish)
        (directory / 'ready').write_text('ready')
        time.sleep(120)
    """))
    wrapper = tmp_path / "wrapper.py"
    wrapper.write_text(textwrap.dedent("""
        import argparse, os, sys
        from pathlib import Path
        from clawbio.common.nfcore import execution
        directory = Path(__file__).parent
        parser = argparse.ArgumentParser()
        parser.add_argument('--output', required=True)
        parser.add_argument('--demo', action='store_true')
        args = parser.parse_args()
        (directory / 'wrapper.pid').write_text(str(os.getpid()))
        original = execution.subprocess.Popen
        def record(*args, **kwargs):
            process = original(*args, **kwargs)
            (directory / 'native.pid').write_text(str(process.pid))
            return process
        execution.subprocess.Popen = record
        try:
            execution.execute_nextflow(
                [sys.executable, str(directory / 'native.py')], cwd=directory,
                output_dir=Path(args.output), timeout_seconds=60, timeout_fix='fixture',
            )
        finally:
            (directory / 'cleanup-finished').write_text('finished')
    """))
    driver = tmp_path / "driver.py"
    driver.write_text(textwrap.dedent(f"""
        import json, sys
        from pathlib import Path
        sys.path.insert(0, {str(ROOT)!r})
        from clawbio import cli
        directory = Path(__file__).parent
        alias = {alias!r}
        cli.SKILLS[alias] = {{**cli.SKILLS[alias], 'script': directory / 'wrapper.py',
                             'demo_args': ['--demo']}}
        try:
            result = cli.run_skill(alias, demo=True, output_dir=str(directory / 'out'),
                                   timeout={3 if stop == 'timeout' else 30})
        except KeyboardInterrupt:
            if {stop!r} != 'SIGINT-caught':
                raise
            (directory / 'interrupt-caught').write_text('caught')
        else:
            (directory / 'runner-result.json').write_text(json.dumps(result))
    """))
    env = {**os.environ, "PYTHONPATH": str(ROOT),
           "CLAWBIO_AUDIT_LOG": str(tmp_path / "audit.jsonl")}
    with (tmp_path / "driver.log").open("w") as log:
        process = subprocess.Popen(
            [sys.executable, str(driver)], env=env, stdout=log,
            stderr=subprocess.STDOUT, start_new_session=True,
        )
        try:
            deadline = time.monotonic() + 15
            while not (tmp_path / "ready").exists():
                assert process.poll() is None, (tmp_path / "driver.log").read_text()
                assert time.monotonic() < deadline, "Fixture did not start"
                time.sleep(0.02)
            if stop != "timeout":
                process.send_signal(getattr(signal, stop.split("-")[0]))
                if stop.endswith("-twice"):
                    time.sleep(0.05)
                    if process.poll() is None:
                        process.send_signal(signal.SIGINT)
            returncode = process.wait(timeout=15)
            assert (tmp_path / "cleanup-finished").exists(), "Wrapper was killed before cleanup"
            assert not _alive(int((tmp_path / "native.pid").read_text()))
            if stop == "timeout":
                result = json.loads((tmp_path / "runner-result.json").read_text())
                assert returncode == 0
                assert result["success"] is False and result["exit_code"] == -1
                assert "Timed out after 3 seconds" in result["stderr"]
            elif stop == "SIGINT-caught":
                assert returncode == 0 and (tmp_path / "interrupt-caught").exists()
            elif stop.startswith("SIGINT"):
                assert returncode in (128 + signal.SIGINT, -signal.SIGINT)
            else:
                assert returncode == 128 + signal.SIGTERM
        finally:
            _stop_fixture(tmp_path, process)


@pytest.mark.parametrize("exit_code", (0, 7))
def test_managed_launcher_preserves_capture_status_and_signal_handlers(exit_code, tmp_path):
    from clawbio.common.nfcore.launcher import run_wrapper

    previous = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
    result = run_wrapper(
        [sys.executable, "-c", f"import sys; print('out'); print('err',file=sys.stderr); sys.exit({exit_code})"],
        cwd=str(tmp_path), timeout=5,
    )
    assert (result.returncode, result.stdout, result.stderr) == (exit_code, "out\n", "err\n")
    assert {s: signal.getsignal(s) for s in previous} == previous


def test_managed_launcher_restores_handlers_after_spawn_failure(tmp_path):
    from clawbio.common.nfcore.launcher import run_wrapper

    previous = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
    with pytest.raises(FileNotFoundError):
        run_wrapper([str(tmp_path / "missing-executable")], cwd=str(tmp_path), timeout=5)
    assert {s: signal.getsignal(s) for s in previous} == previous


def test_managed_launcher_can_run_from_a_worker_thread(tmp_path):
    from concurrent.futures import ThreadPoolExecutor
    from clawbio.common.nfcore.launcher import run_wrapper

    previous = {s: signal.getsignal(s) for s in (signal.SIGINT, signal.SIGTERM)}
    with ThreadPoolExecutor(max_workers=1) as pool:
        result = pool.submit(
            run_wrapper, [sys.executable, "-c", "print('thread')"], cwd=str(tmp_path), timeout=5,
        ).result(timeout=10)
    assert result.returncode == 0 and result.stdout == "thread\n"
    assert {s: signal.getsignal(s) for s in previous} == previous
