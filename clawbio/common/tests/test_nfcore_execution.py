"""Real subprocess execution and injectable failure-policy contracts."""
import importlib
from pathlib import Path
import subprocess
import sys
import pytest


def _run(api, command, output, **overrides):
    return api.execute_nextflow(command, cwd=output, output_dir=output, timeout_seconds=overrides.pop('timeout_seconds', 10), timeout_fix='site-specific timeout advice', terminate_process=api.terminate_process_tree, terminate_on_interrupt=False, failure_hint=lambda stdout, stderr, log: api.environment_failure_hints(stdout, stderr, log), **overrides)


def test_execution_captures_real_stdout_stderr(tmp_path):
    api = importlib.import_module('clawbio.common.nfcore.execution')
    result = _run(api, [sys.executable, '-c', 'import sys; print("hello"); print("oops",file=sys.stderr)'], tmp_path, timeout_seconds=None)
    assert result == {'exit_code': 0, 'stdout_path': str(tmp_path/'logs/stdout.txt'), 'stderr_path': str(tmp_path/'logs/stderr.txt')}
    assert (tmp_path/'logs/stdout.txt').read_text() == 'hello\n'
    assert (tmp_path/'logs/stderr.txt').read_text() == 'oops\n'


def test_launch_failure_removes_empty_logs(tmp_path):
    api = importlib.import_module('clawbio.common.nfcore.execution')
    errors = importlib.import_module('clawbio.common.nfcore.errors')
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [str(tmp_path/'missing-executable')], tmp_path)
    assert failure.value.stage == 'preflight'
    assert failure.value.error_code == 'EXECUTION_FAILED'
    assert not list((tmp_path/'logs').iterdir())


def test_nonzero_execution_reads_nextflow_log_hint(tmp_path):
    api = importlib.import_module('clawbio.common.nfcore.execution')
    errors = importlib.import_module('clawbio.common.nfcore.errors')
    (tmp_path/'.nextflow.log').write_text('java.net.UnknownHostException')
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [sys.executable, '-c', 'raise SystemExit(7)'], tmp_path)
    assert failure.value.details['exit_code'] == 7
    assert 'IPv6' in failure.value.fix


def test_timeout_uses_supplied_advice_and_cleans_up(tmp_path):
    api = importlib.import_module('clawbio.common.nfcore.execution')
    errors = importlib.import_module('clawbio.common.nfcore.errors')
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [sys.executable, '-c', 'import time; time.sleep(60)'], tmp_path, timeout_seconds=0.02)
    assert failure.value.fix == 'site-specific timeout advice'
    assert failure.value.details == {'timeout_seconds': 0.02}
