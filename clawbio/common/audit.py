"""Central audit log for ClawBio.

Aligns with OpenTelemetry GenAI semantic conventions (pre-stable):
  https://github.com/open-telemetry/semantic-conventions-genai
Re-verify the gen_ai.* attribute names once that repo tags a release.
Also emits OpenInference openinference.span.kind/tool.name, which is what
Phoenix classifies on.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Sequence

from opentelemetry import context as _otel_context
from opentelemetry.sdk.trace import ReadableSpan, TracerProvider
from opentelemetry.sdk.resources import Resource
from opentelemetry.sdk.trace.export import (
    BatchSpanProcessor,
    SimpleSpanProcessor,
    SpanExporter,
    SpanExportResult,
)
from opentelemetry.trace import StatusCode

_DEFAULT_LOG = Path.home() / ".clawbio" / "audit.jsonl"
_REDACTED = "__REDACTED__"


def _hide_either(value: str) -> str:
    """For text that can quote both sides, such as a failure message."""
    return _hide(_hide(value, "OPENINFERENCE_HIDE_INPUTS"), "OPENINFERENCE_HIDE_OUTPUTS")


def _hide(value: str, env_var: str) -> str:
    """OpenInference's masking flags, per its configuration spec.

    The marker is kept rather than the attribute dropped, so a reader can tell
    the value was hidden on purpose. Scoped to the input.value/output.value
    aliases the spec covers; clawbio.* and gen_ai.* keys are untouched.
    """
    if value and os.environ.get(env_var, "").strip().lower() == "true":
        return _REDACTED
    return value
_TRACER_KEY = _otel_context.create_key("clawbio.tracer")


def _set_append_only(path: Path) -> None:
    if sys.platform != "darwin":
        return
    try:
        subprocess.run(["chflags", "uappend", str(path)], check=False, capture_output=True)
    except OSError:
        pass


def _ns_to_iso(ns: int) -> str:
    return datetime.fromtimestamp(ns / 1e9, tz=timezone.utc).isoformat()


class _JsonlExporter(SpanExporter):
    """Exports OTEL spans as JSONL to a local file."""

    def __init__(self, log_path: Path) -> None:
        self._log_path = Path(log_path)

    def export(self, spans: Sequence[ReadableSpan]) -> SpanExportResult:
        self._log_path.parent.mkdir(parents=True, exist_ok=True)
        try:
            with self._log_path.open("a", encoding="utf-8") as f:
                for span in spans:
                    record = {
                        "timestamp": _ns_to_iso(span.start_time),
                        "event": span.name,
                        "trace_id": f"{span.context.trace_id:032x}",
                        "span_id": f"{span.context.span_id:016x}",
                        "duration_ms": round((span.end_time - span.start_time) / 1e6, 3),
                        "status": span.status.status_code.name,
                    }
                    if span.parent:
                        record["parent_span_id"] = f"{span.parent.span_id:016x}"
                    record.update(dict(span.attributes or {}))
                    f.write(json.dumps(record, default=str) + "\n")
            _set_append_only(self._log_path)
        except OSError:
            pass
        return SpanExportResult.SUCCESS

    def shutdown(self) -> None:
        pass


def write(event: str, *, log_path: Path | str = _DEFAULT_LOG, **kwargs) -> None:
    """Write a simple point-in-time audit record as JSONL."""
    log_path = Path(log_path)
    entry = {"timestamp": datetime.now(timezone.utc).isoformat(), "event": event, **kwargs}
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(entry, default=str) + "\n")
        _set_append_only(log_path)
    except OSError:
        pass


@contextmanager
def skill_run(
    skill: str,
    version: str,
    input_checksum: str = "",
    input_file: str = "",
    output_dir: str = "",
    log_path: Path | str = _DEFAULT_LOG,
):
    """Root trace for a skill invocation. Yields the span_id (16-char hex).

    PII warning: ``input_file``, ``output_dir``, and any future kwargs are written
    to ``~/.clawbio/audit.jsonl``, and, when ``CLAWBIO_OTLP_ENDPOINT`` is set, are
    also sent to that collector along with every child span. A local collector
    keeps them on this machine; a remote endpoint does not. Callers must scrub
    patient identifiers (VCF paths, sample IDs, free-text fields) before passing
    them here.

    OPENINFERENCE_HIDE_INPUTS / OPENINFERENCE_HIDE_OUTPUTS redact input.value
    and output.value, and either flag also redacts a failure's error text and
    captured stderr, both of which quote the command. That is every value these
    spans carry apart from the input checksum. Still not a substitute for
    scrubbing: a caller that puts an identifier in a skill name or an attrs
    key is outside their reach.
    """
    provider = TracerProvider(resource=Resource.create({
        "service.name": "clawbio",
        # Phoenix groups traces by this; without it everything lands in "default".
        "openinference.project.name": "clawbio",
    }))
    provider.add_span_processor(SimpleSpanProcessor(_JsonlExporter(Path(log_path))))
    endpoint = os.environ.get("CLAWBIO_OTLP_ENDPOINT")
    if endpoint:
        # Opt-in, and deliberately not OTEL_EXPORTER_OTLP_ENDPOINT: an org-wide
        # setting would otherwise point spans carrying output paths at whatever
        # collector the org runs, which may not be on this machine.
        from opentelemetry.exporter.otlp.proto.http.trace_exporter import OTLPSpanExporter

        provider.add_span_processor(
            BatchSpanProcessor(
                # Short timeout: a dead collector must not make shutdown() look
                # like a hung analysis.
                OTLPSpanExporter(
                    endpoint=f"{endpoint.rstrip('/')}/v1/traces", timeout=5
                )
            )
        )
    tracer = provider.get_tracer("clawbio")

    ctx = _otel_context.set_value(_TRACER_KEY, tracer)
    token = _otel_context.attach(ctx)
    try:
        # The SDK would otherwise record an exception event and overwrite the
        # status description, both carrying the raw text past the hide flags.
        with tracer.start_as_current_span(
            "skill_run", record_exception=False, set_status_on_exception=False
        ) as span:
            span.set_attribute("gen_ai.agent.id", skill)
            span.set_attribute("gen_ai.agent.version", version)
            # Phoenix classifies on this key alone; without it: "unknown".
            span.set_attribute("openinference.span.kind", "AGENT")
            # Each omitted when empty. The checksum keeps a clawbio.* namespace
            # because no convention has an equivalent; the file and the output
            # dir are the viewers' Input and Output, one key each, so the
            # masking flags actually hide them.
            for key, value in (
                ("clawbio.input.checksum", input_checksum),
                ("input.value", _hide(input_file, "OPENINFERENCE_HIDE_INPUTS")),
                ("output.value", _hide(output_dir, "OPENINFERENCE_HIDE_OUTPUTS")),
            ):
                if value:
                    span.set_attribute(key, value)
            try:
                yield f"{span.context.span_id:016x}"
                span.set_status(StatusCode.OK)
            except Exception as exc:
                # The text quotes the failing command, so an input flag hides it.
                span.set_attribute("error", _hide_either(str(exc)))
                span.set_status(StatusCode.ERROR, _hide_either(str(exc)))
                raise
    finally:
        _otel_context.detach(token)
        provider.shutdown()  # flush; a short run exits before the next batch tick


@contextmanager
def tool_call(
    name: str,
    *,
    cmd: List[str] | None = None,
    log_path: Path | str = _DEFAULT_LOG,
    **attrs,
):
    """Child span for a tool or CLI call.

    Span name: ``execute_tool {name}``
    Pass cmd to run a subprocess and capture its exit code automatically.

    PII warning: ``cmd`` tokens, captured ``stderr`` and ``**attrs`` are written
    verbatim to the audit log, and go to the OTLP collector too when the
    enclosing ``skill_run`` has ``CLAWBIO_OTLP_ENDPOINT`` set. Callers must scrub
    file paths, sample IDs, and any patient-identifiable values before passing
    them here.
    """
    tracer = _otel_context.get_value(_TRACER_KEY)
    if tracer is None:
        yield None
        return

    with tracer.start_as_current_span(
        f"execute_tool {name}", record_exception=False, set_status_on_exception=False
    ) as span:
        span.set_attribute("openinference.span.kind", "TOOL")
        span.set_attribute("tool.name", name)
        if cmd is not None:
            span.set_attribute(
                "input.value", _hide(" ".join(cmd), "OPENINFERENCE_HIDE_INPUTS")
            )
        for k, v in attrs.items():
            span.set_attribute(k, str(v))
        try:
            if cmd is not None:
                result = subprocess.run(cmd, capture_output=True, text=True)
                span.set_attribute("exit_code", result.returncode)
                span.set_attribute("output.value", _hide(
                    f"exit_code={result.returncode}", "OPENINFERENCE_HIDE_OUTPUTS"
                ))
                if result.returncode != 0:
                    span.set_attribute("error.type", "NonZeroExit")
                    # stderr quotes the arguments it was given as often as it
                    # reports a result, so either flag hides it.
                    span.set_attribute("stderr", _hide_either(result.stderr[:500]))
                    span.set_status(StatusCode.ERROR, f"exit {result.returncode}")
                    raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
            yield f"{span.context.span_id:016x}"
            span.set_status(StatusCode.OK)
        except subprocess.CalledProcessError:
            raise
        except Exception as exc:
            span.set_attribute("error.type", type(exc).__name__)
            span.set_attribute("error", _hide_either(str(exc)))
            span.set_status(StatusCode.ERROR, _hide_either(str(exc)))
            raise
