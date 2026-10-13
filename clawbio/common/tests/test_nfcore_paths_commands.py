"""URI preservation, command argv and display quoting contracts."""
import importlib
from pathlib import Path
import pytest
from hypothesis import given, strategies as st


def test_command_preserves_order_and_remote_work_dir(tmp_path):
    api = importlib.import_module('clawbio.common.nfcore.commands')
    argv, display = api.build_nextflow_command(pipeline_source={'source_kind': 'remote_repo', 'source_ref': 'nf-core/example', 'resolved_version': '1.2.3'}, profile='test,docker', params_path=Path('params with space.yaml'), resume=True, work_dir='s3://bucket/work', extra_configs=[Path('a.config'), Path('a.config')])
    assert argv == ['nextflow', 'run', 'nf-core/example', '-r', '1.2.3', '-profile', 'test,docker', '-params-file', 'params with space.yaml', '-work-dir', 's3://bucket/work', '-c', 'a.config', '-c', 'a.config', '-resume']
    assert "'params with space.yaml'" in display


def test_local_checkout_omits_revision():
    api = importlib.import_module('clawbio.common.nfcore.commands')
    argv, _ = api.build_nextflow_command(pipeline_source={'source_kind': 'local_checkout', 'source_ref': '/local/pipeline', 'resolved_version': 'commit'}, profile='docker', params_path=Path('params.yaml'), resume=False)
    assert argv == ['nextflow', 'run', '/local/pipeline', '-profile', 'docker', '-params-file', 'params.yaml']


def test_windows_display_does_not_change_argv(monkeypatch):
    api = importlib.import_module('clawbio.common.nfcore.commands')
    monkeypatch.setattr(api.sys, 'platform', 'win32')
    argv, display = api.build_nextflow_command(pipeline_source={'source_kind': 'remote_repo', 'source_ref': 'nf-core/example', 'resolved_version': '1'}, profile='docker', params_path=Path('params with space.yaml'), resume=False)
    assert argv[-1] == 'params with space.yaml'
    assert display.endswith('"params with space.yaml"')


@given(st.sampled_from(['s3', 'gs', 'https', 'ftp']), st.text(alphabet='abcXYZ012 -_é', min_size=1, max_size=20))
def test_uri_paths_remain_verbatim(scheme, suffix):
    api = importlib.import_module('clawbio.common.nfcore.paths')
    value = f'{scheme}://bucket/{suffix}'
    assert api.is_remote_uri(value)
    assert api.resolve_path_or_uri(value, base_dir=Path('/unused')) == value


def test_local_resolution_uses_explicit_base_and_segment_containment(tmp_path):
    api = importlib.import_module('clawbio.common.nfcore.paths')
    assert api.resolve_path_or_uri('reads.fq', base_dir=tmp_path) == tmp_path / 'reads.fq'
    assert api.is_within(tmp_path / 'data' / 'reads.fq', tmp_path / 'data')
    assert not api.is_within(tmp_path / 'data2', tmp_path / 'data')
    assert api.is_under_tmp(Path('/tmp/child'))
