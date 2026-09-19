"""Tests for clawbio.common.audit."""

import json
import stat
import sys
from pathlib import Path

import pytest

from clawbio.common.audit import write, skill_run, tool_call


def test_write_appends_jsonl_record(tmp_path):
    log = tmp_path / "audit.jsonl"
    write("user_event", skill="pharmgx", version="0.2.0", log_path=log)
    records = [json.loads(l) for l in log.read_text().splitlines()]
    assert len(records) == 1
    assert records[0]["event"] == "user_event"
    assert records[0]["skill"] == "pharmgx"
    assert "timestamp" in records[0]


def test_write_appends_multiple_records(tmp_path):
    log = tmp_path / "audit.jsonl"
    write("user_event", skill="pharmgx", log_path=log)
    write("user_event", skill="gwas-prs", log_path=log)
    records = [json.loads(l) for l in log.read_text().splitlines()]
    assert len(records) == 2
    assert records[1]["skill"] == "gwas-prs"


def test_write_creates_parent_dirs(tmp_path):
    log = tmp_path / "nested" / "deep" / "audit.jsonl"
    write("user_event", skill="pharmgx", log_path=log)
    assert log.exists()


def test_write_silently_ignores_oserror(tmp_path):
    log = tmp_path / "audit.jsonl"
    log.parent.chmod(stat.S_IREAD | stat.S_IEXEC)
    try:
        write("user_event", skill="pharmgx", log_path=log)
    finally:
        log.parent.chmod(stat.S_IRWXU)


def test_skill_run_writes_single_otel_span(tmp_path):
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", input_checksum="abc", log_path=log):
        pass
    records = [json.loads(l) for l in log.read_text().splitlines()]
    assert len(records) == 1
    assert records[0]["event"] == "skill_run"
    assert records[0]["status"] == "OK"


def test_skill_run_record_has_duration(tmp_path):
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", input_checksum="abc", log_path=log):
        pass
    record = json.loads(log.read_text().strip())
    assert "duration_ms" in record
    assert record["duration_ms"] >= 0


def test_skill_run_record_has_trace_and_span_ids(tmp_path):
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", input_checksum="abc", log_path=log) as span_id:
        pass
    record = json.loads(log.read_text().strip())
    assert len(record["span_id"]) == 16
    assert len(record["trace_id"]) == 32
    assert record["span_id"] == span_id


def test_skill_run_failed_sets_error_status(tmp_path):
    log = tmp_path / "audit.jsonl"
    with pytest.raises(ValueError):
        with skill_run("pharmgx", "0.2.0", input_checksum="abc", log_path=log):
            raise ValueError("boom")
    record = json.loads(log.read_text().strip())
    assert record["status"] == "ERROR"
    assert record["error"] == "boom"


def test_skill_run_record_has_required_attributes(tmp_path):
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        pass
    record = json.loads(log.read_text().strip())
    assert record["gen_ai.agent.id"] == "pharmgx"
    assert record["gen_ai.agent.version"] == "0.2.0"
    assert "timestamp" in record


def test_tool_call_follows_genai_spec(tmp_path):
    log = tmp_path / "audit.jsonl"
    with skill_run("gwas-lookup", "0.1.0", input_checksum="abc", log_path=log):
        with tool_call("opengwas_api", rsid="rs3798220", log_path=log):
            pass
    records = [json.loads(l) for l in log.read_text().splitlines()]
    skill = next(r for r in records if r["event"] == "skill_run")
    tool = next(r for r in records if r["event"] == "execute_tool opengwas_api")
    assert tool["trace_id"] == skill["trace_id"]
    assert tool["parent_span_id"] == skill["span_id"]
    assert "duration_ms" in tool


def test_tool_call_failed_sets_error_status(tmp_path):
    log = tmp_path / "audit.jsonl"
    with pytest.raises(RuntimeError):
        with skill_run("gwas-lookup", "0.1.0", input_checksum="abc", log_path=log):
            with tool_call("opengwas_api", log_path=log):
                raise RuntimeError("timeout")
    records = [json.loads(l) for l in log.read_text().splitlines()]
    tool = next(r for r in records if r["event"] == "execute_tool opengwas_api")
    assert tool["status"] == "ERROR"
    assert tool["error.type"] == "RuntimeError"


def test_tool_call_runs_subprocess_via_cmd(tmp_path):
    log = tmp_path / "audit.jsonl"
    with skill_run("seq-wrangler", "0.1.0", input_checksum="abc", log_path=log):
        with tool_call("echo", cmd=["echo", "hello"], log_path=log):
            pass
    records = [json.loads(l) for l in log.read_text().splitlines()]
    cli = next(r for r in records if r["event"] == "execute_tool echo")
    assert "duration_ms" in cli
    assert cli["exit_code"] == 0


def test_tool_call_captures_exit_code_on_failure(tmp_path):
    log = tmp_path / "audit.jsonl"
    with pytest.raises(Exception):
        with skill_run("seq-wrangler", "0.1.0", input_checksum="abc", log_path=log):
            with tool_call("false", cmd=["false"], log_path=log):
                pass
    records = [json.loads(l) for l in log.read_text().splitlines()]
    cli = next(r for r in records if r["event"] == "execute_tool false")
    assert cli["status"] == "ERROR"
    assert "exit_code" in cli


@pytest.mark.skipif(sys.platform != "darwin", reason="chflags is macOS-only")
def test_uappend_flag_prevents_truncation(tmp_path):
    log = tmp_path / "audit.jsonl"
    write("user_event", skill="pharmgx", log_path=log)
    with pytest.raises(OSError):
        log.write_text("wiped")


def test_otlp_exporter_is_opt_in(tmp_path, monkeypatch):
    """Off by default; CLAWBIO_OTLP_ENDPOINT constructs an OTLP exporter."""
    from opentelemetry.exporter.otlp.proto.http import trace_exporter
    from opentelemetry.sdk.trace.export import SpanExportResult

    endpoints = []

    class _Recorder:
        def __init__(self, endpoint=None, **kwargs):
            endpoints.append(endpoint)

        def export(self, spans):
            return SpanExportResult.SUCCESS

        def shutdown(self):
            pass

        def force_flush(self, timeout_millis=None):
            return True

    monkeypatch.setattr(trace_exporter, "OTLPSpanExporter", _Recorder)

    monkeypatch.delenv("CLAWBIO_OTLP_ENDPOINT", raising=False)
    with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "a.jsonl"):
        pass
    assert endpoints == []

    monkeypatch.setenv("CLAWBIO_OTLP_ENDPOINT", "http://localhost:4318")
    with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "b.jsonl"):
        pass
    assert endpoints == ["http://localhost:4318/v1/traces"]


def test_skill_run_records_provenance_attributes(tmp_path):
    """input_checksum, input_file and output_dir reach the log, as documented."""
    log = tmp_path / "audit.jsonl"
    with skill_run(
        "pharmgx",
        "0.2.0",
        input_checksum="abc123",
        input_file="demo_patient.txt",
        output_dir=str(tmp_path / "pharmgx_demo"),
        log_path=log,
    ):
        pass
    record = json.loads(log.read_text().strip())
    assert record["clawbio.input.checksum"] == "abc123"
    assert record["clawbio.input.file"] == "demo_patient.txt"
    assert record["clawbio.output.dir"] == str(tmp_path / "pharmgx_demo")


def test_skill_run_omits_provenance_keys_when_not_given(tmp_path):
    """Callers that pass nothing keep the historical record shape."""
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        pass
    record = json.loads(log.read_text().strip())
    for key in ("clawbio.input.checksum", "clawbio.input.file", "clawbio.output.dir"):
        assert key not in record


def test_span_kinds_are_set_for_trace_viewers(tmp_path):
    """Phoenix renders spans as "unknown" without openinference.span.kind."""
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        with tool_call("bcftools_view", log_path=log):
            pass
    records = [json.loads(line) for line in log.read_text().splitlines()]
    tool = next(r for r in records if r["event"] == "execute_tool bcftools_view")
    root = next(r for r in records if r["event"] == "skill_run")
    assert root["openinference.span.kind"] == "AGENT"
    assert tool["openinference.span.kind"] == "TOOL"
    assert tool["tool.name"] == "bcftools_view"
