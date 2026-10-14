"""Real subprocess execution and injectable failure-policy contracts."""

import importlib
import sys
from pathlib import Path

import pytest


def _run(api, command, output, **overrides):
    return api.execute_nextflow(
        command,
        cwd=output,
        output_dir=output,
        timeout_seconds=overrides.pop("timeout_seconds", 10),
        timeout_fix="site-specific timeout advice",
        **overrides,
    )


def test_execution_captures_real_stdout_stderr(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    result = _run(
        api,
        [
            sys.executable,
            "-c",
            'import sys; print("hello"); print("oops",file=sys.stderr)',
        ],
        tmp_path,
        timeout_seconds=None,
    )
    assert result == {
        "exit_code": 0,
        "stdout_path": str(tmp_path / "logs/stdout.txt"),
        "stderr_path": str(tmp_path / "logs/stderr.txt"),
    }
    assert (tmp_path / "logs/stdout.txt").read_text() == "hello\n"
    assert (tmp_path / "logs/stderr.txt").read_text() == "oops\n"


def test_launch_failure_removes_empty_logs(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [str(tmp_path / "missing-executable")], tmp_path)
    assert failure.value.stage == "preflight"
    assert failure.value.error_code == "EXECUTION_FAILED"
    assert not list((tmp_path / "logs").iterdir())


def test_nonzero_execution_reads_nextflow_log_hint(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    (tmp_path / ".nextflow.log").write_text("java.net.UnknownHostException")
    with pytest.raises(errors.SkillError) as failure:
        _run(api, [sys.executable, "-c", "raise SystemExit(7)"], tmp_path)
    assert failure.value.details["exit_code"] == 7
    assert "IPv6" in failure.value.fix


def test_timeout_uses_supplied_advice_and_cleans_up(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    with pytest.raises(errors.SkillError) as failure:
        _run(
            api,
            [sys.executable, "-c", "import time; time.sleep(60)"],
            tmp_path,
            timeout_seconds=0.02,
        )
    assert failure.value.fix == "site-specific timeout advice"
    assert failure.value.details == {"timeout_seconds": 0.02}


def test_pipeline_hint_follows_the_shared_environment_hints(tmp_path):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    errors = importlib.import_module("clawbio.common.nfcore.errors")
    (tmp_path / ".nextflow.log").write_text("Network is unreachable")
    with pytest.raises(errors.SkillError) as failure:
        _run(
            api,
            [sys.executable, "-c", "raise SystemExit(1)"],
            tmp_path,
            pipeline_failure_hint=lambda stdout, stderr: " PIPELINE-HINT",
        )
    fix = failure.value.fix
    assert fix.startswith("Inspect logs/stdout.txt")
    assert fix.index("preferIPv6Addresses") < fix.index(" PIPELINE-HINT")
    assert fix.endswith(" PIPELINE-HINT")


@pytest.mark.parametrize(
    ("platform", "under_tmp", "expected"),
    [("darwin", True, True), ("darwin", False, False), ("linux", True, False)],
)
def test_macos_tmp_hint_only_for_macos_output_under_tmp(
    platform, under_tmp, expected, monkeypatch
):
    api = importlib.import_module("clawbio.common.nfcore.execution")
    monkeypatch.setattr(sys, "platform", platform)
    monkeypatch.setattr(api, "is_under_tmp", lambda _path: under_tmp)
    hint = api.macos_tmp_failure_hint(Path("/tmp/run"))
    assert bool(hint) is expected
    if expected:
        assert "/tmp" in hint and "HOME" in hint


def test_pipeline_environment_overrides_the_callers(tmp_path, monkeypatch):
    """A pipeline that needs a specific Nextflow setting (e.g. the v1 syntax parser
    for releases predating strict syntax) must get it even if the user's shell
    exports the opposite; everything else is inherited."""
    api = importlib.import_module("clawbio.common.nfcore.execution")
    monkeypatch.setenv("NXF_SYNTAX_PARSER", "v2")
    monkeypatch.setenv("CLAWBIO_INHERITED", "yes")
    _run(
        api,
        ["sh", "-c", 'printf "%s %s" "$NXF_SYNTAX_PARSER" "$CLAWBIO_INHERITED"'],
        tmp_path,
        env={"NXF_SYNTAX_PARSER": "v1"},
    )
    assert (tmp_path / "logs/stdout.txt").read_text() == "v1 yes"
