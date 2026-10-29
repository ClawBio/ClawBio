"""Tests for clawbio.common.audit."""

import json
import stat
from pathlib import Path

import subprocess
import sys
import pytest

from clawbio.common import audit
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


def test_otlp_exporter_is_opt_in(tmp_path, monkeypatch):
    """Off by default; CLAWBIO_OTLP_ENDPOINT constructs an OTLP exporter."""
    trace_exporter = pytest.importorskip(
        "opentelemetry.exporter.otlp.proto.http.trace_exporter", reason="needs the otlp extra (uv sync --extra otlp)"
    )
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


def test_missing_otlp_extra_warns_and_keeps_the_analysis_running(tmp_path, monkeypatch, capsys):
    """Without clawbio[otlp], a set endpoint must cost the traces, not the run."""
    monkeypatch.setitem(sys.modules, "opentelemetry.exporter.otlp.proto.http.trace_exporter", None)
    monkeypatch.setenv("CLAWBIO_OTLP_ENDPOINT", "http://localhost:4318")
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        pass
    assert "clawbio[otlp]" in capsys.readouterr().err
    assert log.read_text().strip(), "local audit log should still be written"


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
    # The file and the output dir live in input.value / output.value: one key
    # each, so the masking flags cover them.
    assert record["input.value"] == "demo_patient.txt"
    assert record["output.value"] == str(tmp_path / "pharmgx_demo")
    assert "clawbio.input.file" not in record
    assert "clawbio.output.dir" not in record


def test_skill_run_omits_provenance_keys_when_not_given(tmp_path):
    """Callers that pass nothing keep the historical record shape."""
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        pass
    record = json.loads(log.read_text().strip())
    for key in ("clawbio.input.checksum", "input.value", "output.value"):
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
    # A skill run is a fixed pipeline, OpenInference's CHAIN; AGENT means a
    # reasoning loop that chooses its own steps.
    assert root["openinference.span.kind"] == "CHAIN"
    assert tool["openinference.span.kind"] == "TOOL"
    assert tool["tool.name"] == "bcftools_view"


def test_input_output_values_alias_existing_attributes(tmp_path):
    """Phoenix fills its Input/Output panels from input.value / output.value.
    Both are aliases, so nothing new reaches the log or the collector."""
    log = tmp_path / "audit.jsonl"
    out = tmp_path / "pharmgx_demo"
    with skill_run(
        "pharmgx", "0.2.0", input_file="demo_patient.txt", output_dir=str(out), log_path=log
    ):
        with tool_call("bcftools_view", cmd=["echo", "hi"], log_path=log):
            pass
    records = [json.loads(line) for line in log.read_text().splitlines()]
    root = next(r for r in records if r["event"] == "skill_run")
    tool = next(r for r in records if r["event"] == "execute_tool bcftools_view")
    assert root["input.value"] == "demo_patient.txt"
    assert root["output.value"] == str(out)
    assert tool["input.value"] == "echo hi"
    assert tool["output.value"] == "exit_code=0"


def test_input_output_values_omitted_when_there_is_nothing_to_alias(tmp_path):
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        with tool_call("no_subprocess", log_path=log):
            pass
    for record in (json.loads(line) for line in log.read_text().splitlines()):
        assert "input.value" not in record
        assert "output.value" not in record


def test_openinference_hide_flags_redact_the_alias_values(tmp_path, monkeypatch):
    """OPENINFERENCE_HIDE_INPUTS/OUTPUTS replace input.value / output.value with
    the spec's __REDACTED__ marker. Scoped to the OpenInference attributes, as
    the spec defines them: the clawbio.* keys are a separate decision."""
    monkeypatch.setenv("OPENINFERENCE_HIDE_INPUTS", "true")
    monkeypatch.setenv("OPENINFERENCE_HIDE_OUTPUTS", "true")
    log = tmp_path / "audit.jsonl"
    with skill_run(
        "pharmgx", "0.2.0", input_file="demo_patient.txt",
        output_dir=str(tmp_path / "out"), log_path=log,
    ):
        with tool_call("bcftools_view", cmd=["echo", "hi"], log_path=log):
            pass
    records = [json.loads(line) for line in log.read_text().splitlines()]
    root = next(r for r in records if r["event"] == "skill_run")
    tool = next(r for r in records if r["event"] == "execute_tool bcftools_view")
    assert root["input.value"] == "__REDACTED__"
    assert root["output.value"] == "__REDACTED__"
    assert tool["input.value"] == "__REDACTED__"
    assert tool["output.value"] == "__REDACTED__"
    # Nothing else carries the command either.
    assert "gen_ai.tool.call.arguments" not in tool
    # Still recorded: the span exists, only the value is hidden.
    assert root["gen_ai.agent.id"] == "pharmgx"
    assert tool["exit_code"] == 0
    # Nothing else carries the path, so hiding it hides it.
    assert "clawbio.input.file" not in root


def test_failed_tool_call_does_not_leak_hidden_values(tmp_path, monkeypatch):
    """The error path carried the command back: stderr echoes what a tool was
    given, and CalledProcessError's text embeds the whole argv."""
    monkeypatch.setenv("OPENINFERENCE_HIDE_INPUTS", "true")
    log = tmp_path / "audit.jsonl"
    secret = "Doe_J_sample.vcf"
    with pytest.raises(subprocess.CalledProcessError):
        with skill_run("pharmgx", "0.2.0", log_path=log):
            with tool_call("reader", cmd=["cat", secret], log_path=log):
                pass
    assert secret not in log.read_text()


def test_failure_does_not_leak_through_otel_exception_recording(tmp_path, monkeypatch):
    """The SDK records an exception event and overwrites the status description
    with the raw exception text. Neither reaches the JSONL log, so both have to
    be read off the spans the exporter is handed."""
    captured = []

    class _Capture(audit._JsonlExporter):
        def export(self, spans):
            captured.extend(spans)
            return super().export(spans)

    monkeypatch.setattr(audit, "_JsonlExporter", _Capture)
    monkeypatch.setenv("OPENINFERENCE_HIDE_INPUTS", "true")
    secret = "Doe_J_sample.vcf"
    with pytest.raises(subprocess.CalledProcessError):
        with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "audit.jsonl"):
            with tool_call("reader", cmd=["cat", secret], log_path=tmp_path / "audit.jsonl"):
                pass

    assert captured, "no spans captured"
    for span in captured:
        assert not span.events, f"{span.name} recorded {[e.name for e in span.events]}"
        assert secret not in (span.status.description or "")


def test_audit_log_path_comes_from_the_environment(tmp_path, monkeypatch):
    """A skill calls skill_run without a log_path, so the environment is the
    only way to point the log somewhere other than $HOME."""
    log = tmp_path / "from_env.jsonl"
    monkeypatch.setenv("CLAWBIO_AUDIT_LOG", str(log))
    with skill_run("fine-mapping", "0.3.0"):
        with tool_call("abf_fit"):
            pass
    events = [json.loads(line)["event"] for line in log.read_text().splitlines()]
    assert events == ["execute_tool abf_fit", "skill_run"]


def test_in_process_tool_call_reports_attrs_as_its_input(tmp_path):
    """A phase with no subprocess has no command to show, so its keyword
    arguments are the input. One key, not a copy in each form."""
    log = tmp_path / "audit.jsonl"
    with skill_run("fine-mapping", "0.3.0", log_path=log):
        with tool_call("susie_fit", n_variants=200, max_signals=10, log_path=log):
            pass
    tool = next(json.loads(l) for l in log.read_text().splitlines()
                if json.loads(l)["event"] == "execute_tool susie_fit")
    assert json.loads(tool["input.value"]) == {"n_variants": 200, "max_signals": 10}
    assert tool["input.mime_type"] == "application/json"
    assert "n_variants" not in tool


def test_hide_inputs_redacts_attrs_passed_alongside_cmd(tmp_path, monkeypatch):
    """With a cmd, every **attrs kwarg was written raw past both hide flags,
    so `tool_call("vep", cmd=[...], input_vcf=path)` shipped the path."""
    monkeypatch.setenv("OPENINFERENCE_HIDE_INPUTS", "true")
    log = tmp_path / "audit.jsonl"
    secret = "Doe_J_sample.vcf"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        with tool_call("vep", cmd=["echo", "hi"], input_vcf=secret, log_path=log):
            pass
    tool = next(json.loads(l) for l in log.read_text().splitlines()
                if json.loads(l)["event"] == "execute_tool vep")
    assert tool["input_vcf"] == "__REDACTED__"


def test_otlp_endpoint_may_be_given_with_or_without_the_signal_path(tmp_path, monkeypatch):
    """OTel habit is to set the full /v1/traces URL; appending it again 404s."""
    trace_exporter = pytest.importorskip(
        "opentelemetry.exporter.otlp.proto.http.trace_exporter", reason="needs the otlp extra (uv sync --extra otlp)"
    )
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
    for value in ("http://localhost:4318", "http://localhost:4318/v1/traces/"):
        monkeypatch.setenv("CLAWBIO_OTLP_ENDPOINT", value)
        with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "a.jsonl"):
            pass
    assert endpoints == ["http://localhost:4318/v1/traces"] * 2


@pytest.fixture
def exported(monkeypatch):
    """Spans as the exporter saw them; the JSONL record leaves out the resource."""
    spans = []

    class _Recorder(audit._JsonlExporter):
        def export(self, batch):
            spans.extend(batch)
            return super().export(batch)

    monkeypatch.setattr(audit, "_JsonlExporter", _Recorder)
    monkeypatch.delenv("OTEL_RESOURCE_ATTRIBUTES", raising=False)
    monkeypatch.delenv("OTEL_SERVICE_NAME", raising=False)
    return spans


def test_resource_defaults_to_clawbio_when_nobody_names_it(tmp_path, exported):
    with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "a.jsonl"):
        pass
    attrs = exported[0].resource.attributes
    assert attrs["service.name"] == "clawbio"
    assert attrs["openinference.project.name"] == "clawbio"


def test_caller_names_the_project_and_service(tmp_path, exported, monkeypatch):
    # An agent that runs clawbio files its spans under its own project; the
    # standard OTel variables are how a parent process says so.
    monkeypatch.setenv("OTEL_RESOURCE_ATTRIBUTES", "openinference.project.name=hermes-agent,deployment=mini")
    monkeypatch.setenv("OTEL_SERVICE_NAME", "hermes")
    with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "a.jsonl"):
        pass
    attrs = exported[0].resource.attributes
    assert attrs["openinference.project.name"] == "hermes-agent"
    assert attrs["service.name"] == "hermes"
    assert attrs["deployment"] == "mini"


def test_spans_name_the_emitting_module_as_their_scope(tmp_path, exported):
    # The scope says which library emitted a span, whoever owns the project.
    with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "a.jsonl"):
        with tool_call("fetch"):
            pass
    assert {s.instrumentation_scope.name for s in exported} == {"clawbio.common.audit"}



@pytest.mark.parametrize("env, project, service", [
    ({"OTEL_RESOURCE_ATTRIBUTES": "openinference.project.name=hermes-agent"}, "hermes-agent", "clawbio"),
    ({"OTEL_SERVICE_NAME": "hermes"}, "clawbio", "hermes"),
])
def test_environment_overrides_only_what_it_sets(tmp_path, exported, monkeypatch, env, project, service):
    # Setting one variable must leave the other field at clawbio's fallback,
    # never "unknown_service".
    for var, value in env.items():
        monkeypatch.setenv(var, value)
    with skill_run("pharmgx", "0.2.0", log_path=tmp_path / "a.jsonl"):
        pass
    attrs = exported[0].resource.attributes
    assert attrs["openinference.project.name"] == project
    assert attrs["service.name"] == service


_PARENT_TRACE = "4bf92f3577b34da6a3ce929d0e0e4736"
_PARENT_SPAN = "00f067aa0ba902b7"


def test_a_nested_skill_run_stays_under_its_outer_run_despite_traceparent(tmp_path, monkeypatch):
    # full-profile calls run_skill once per skill inside its own skill_run.
    monkeypatch.setenv("TRACEPARENT", f"00-{_PARENT_TRACE}-{_PARENT_SPAN}-01")
    log = tmp_path / "audit.jsonl"
    with skill_run("full-profile", "0.1.0", log_path=log) as outer:
        with skill_run("pharmgx", "0.2.0", log_path=log):
            pass
    inner, profile = (json.loads(line) for line in log.read_text().splitlines())
    assert inner["parent_span_id"] == outer
    assert profile["parent_span_id"] == _PARENT_SPAN


_UNSAMPLED = "00-00000000000000000000000000000001-0000000000000001-00"


@pytest.mark.parametrize("value, parent", [
    # An agent that runs clawbio as a subprocess can only pass its span this way.
    (f"00-{_PARENT_TRACE}-{_PARENT_SPAN}-01", (_PARENT_TRACE, _PARENT_SPAN)),
    ("", None),
    ("not-a-traceparent", None),
    (f"00-{'0' * 32}-{_PARENT_SPAN}-01", None),  # all-zero trace id is invalid
    (f"00-{_PARENT_TRACE}-{_PARENT_SPAN}", None),  # flags missing
    (f"ff-{_PARENT_TRACE}-{_PARENT_SPAN}-01", None),  # version ff is forbidden
    # Unsampled still nests, and ALWAYS_ON keeps the audit record.
    (_UNSAMPLED, ("00000000000000000000000000000001", "0000000000000001")),
])
def test_a_traceparent_gives_a_root_or_the_named_parent(tmp_path, monkeypatch, value, parent):
    # The value is set by whatever launched us; a bad one must not break the
    # run, and must not produce a parent that the header did not name.
    monkeypatch.setenv("TRACEPARENT", value)
    log = tmp_path / "audit.jsonl"
    with skill_run("pharmgx", "0.2.0", log_path=log):
        pass
    record = json.loads(log.read_text().strip())
    assert record["status"] == "OK"
    if parent:
        assert (record["trace_id"], record["parent_span_id"]) == parent
    else:
        assert "parent_span_id" not in record
