"""How every nf-core wrapper reports a failure."""

import json

from clawbio.common.nfcore import failures
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
