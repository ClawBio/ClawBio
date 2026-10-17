"""How every nf-core wrapper reports a failure, and the shared preflight checks."""

import json

import pytest

from clawbio.common.nfcore import failures, runtime
from clawbio.common.nfcore.errors import SkillError


def _error(**overrides):
    fields = dict(
        stage="preflight",
        error_code="MISSING_INPUT",
        message="No samplesheet.",
        fix="Pass --input.",
        details={"input": "/x/sheet.csv"},
    )
    fields.update(overrides)
    return SkillError(**fields)


def test_skill_error_is_reported_on_all_three_channels(tmp_path, capsys):
    assert failures.report_skill_error(tmp_path, _error(), verbose=False) == 1
    out, err = capsys.readouterr()
    assert "  code:    MISSING_INPUT" in out and "  fix:     Pass --input." in out
    assert "/x/sheet.csv" not in out  # details only with --verbose
    assert json.loads(err) == _error().to_dict()
    assert json.loads((tmp_path / "result.json").read_text()) == _error().to_dict()


def test_verbose_skill_error_shows_details(tmp_path, capsys):
    failures.report_skill_error(tmp_path, _error(), verbose=True)
    assert '"input": "/x/sheet.csv"' in capsys.readouterr().out


def test_refused_resume_reports_without_touching_the_earlier_result(tmp_path, capsys):
    (tmp_path / "result.json").write_text('{"status": "success"}')
    failures.report_skill_error(tmp_path, _error(error_code="RESUME_DRIFT"), verbose=False)
    assert (tmp_path / "result.json").read_text() == '{"status": "success"}'
    assert json.loads(capsys.readouterr().err)["error_code"] == "RESUME_DRIFT"


def test_unexpected_error_keeps_the_traceback_for_the_bug_report(tmp_path, capsys):
    try:
        raise RuntimeError("boom")
    except RuntimeError as exc:
        assert failures.report_unexpected_error(tmp_path, exc, verbose=False) == 1
    out, err = capsys.readouterr()
    payload = json.loads((tmp_path / "result.json").read_text())
    assert payload["error_code"] == "UNEXPECTED_ERROR"
    assert payload["stage"] == "internal"
    assert payload["details"]["exception_type"] == "RuntimeError"
    assert "RuntimeError: boom" in payload["details"]["traceback"]
    assert json.loads(err) == payload
    assert "Re-run with --verbose" in out and "Traceback" not in out


@pytest.mark.parametrize("value", ["true", "1", "YES", " on "])
def test_demo_refuses_to_start_offline(value, monkeypatch):
    monkeypatch.setenv("NXF_OFFLINE", value)
    with pytest.raises(SkillError) as failure:
        runtime.require_network_for_demo()
    assert failure.value.error_code == "DEMO_REQUIRES_NETWORK"
    assert failure.value.details == {"NXF_OFFLINE": value}


@pytest.mark.parametrize("value", [None, "", "false", "0"])
def test_demo_runs_online(value, monkeypatch):
    if value is None:
        monkeypatch.delenv("NXF_OFFLINE", raising=False)
    else:
        monkeypatch.setenv("NXF_OFFLINE", value)
    runtime.require_network_for_demo()


@pytest.mark.parametrize(
    ("text", "expected"),
    [
        ("  N E X T F L O W\n  version 26.04.3 build 5928", "26.04.3"),
        ('openjdk version "17.0.10" 2024-01-16', "17.0.10"),
        ("version 25.10", "25.10"),
        ("no version here", ""),
    ],
)
def test_reported_version_keeps_the_raw_token(text, expected):
    # 26.04.3 must not become 26.4.3: that is not a release NXF_VER can fetch.
    assert runtime.reported_version(text) == expected
