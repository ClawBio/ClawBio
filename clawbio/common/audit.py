"""Central audit log for ClawBio.

Aligns with OpenTelemetry GenAI semantic conventions (pre-stable):
  https://github.com/open-telemetry/semantic-conventions-genai
Re-verify the gen_ai.* attribute names once that repo tags a release.
"""

from __future__ import annotations

import json
import subprocess
import sys
from contextlib import contextmanager
from datetime import datetime, timezone
from pathlib import Path
from typing import List, Sequence

from opentelemetry import context as _otel_context
from opentelemetry.sdk.trace import ReadableSpan, TracerProvider
from opentelemetry.sdk.trace.export import SimpleSpanProcessor, SpanExporter, SpanExportResult
from opentelemetry.trace import StatusCode

_DEFAULT_LOG = Path.home() / ".clawbio" / "audit.jsonl"
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
    to ``~/.clawbio/audit.jsonl``. Callers must scrub patient identifiers (VCF
    paths, sample IDs, free-text fields) before passing them here.
    """
    provider = TracerProvider()
    provider.add_span_processor(SimpleSpanProcessor(_JsonlExporter(Path(log_path))))
    tracer = provider.get_tracer("clawbio")

    ctx = _otel_context.set_value(_TRACER_KEY, tracer)
    token = _otel_context.attach(ctx)
    try:
        with tracer.start_as_current_span("skill_run") as span:
            span.set_attribute("gen_ai.agent.id", skill)
            span.set_attribute("gen_ai.agent.version", version)
            try:
                yield f"{span.context.span_id:016x}"
                span.set_status(StatusCode.OK)
            except Exception as exc:
                span.set_attribute("error", str(exc))
                span.set_status(StatusCode.ERROR, str(exc))
                raise
    finally:
        _otel_context.detach(token)


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

    PII warning: ``cmd`` tokens and ``**attrs`` are written verbatim to the audit
    log. Callers must scrub file paths, sample IDs, and any patient-identifiable
    values before passing them here.
    """
    tracer = _otel_context.get_value(_TRACER_KEY)
    if tracer is None:
        yield None
        return

    with tracer.start_as_current_span(f"execute_tool {name}") as span:
        if cmd is not None:
            span.set_attribute("gen_ai.tool.call.arguments", " ".join(cmd))
        for k, v in attrs.items():
            span.set_attribute(k, str(v))
        try:
            if cmd is not None:
                result = subprocess.run(cmd, capture_output=True, text=True)
                span.set_attribute("exit_code", result.returncode)
                if result.returncode != 0:
                    span.set_attribute("error.type", "NonZeroExit")
                    span.set_attribute("stderr", result.stderr[:500])
                    span.set_status(StatusCode.ERROR, f"exit {result.returncode}")
                    raise subprocess.CalledProcessError(result.returncode, cmd, result.stderr)
            yield f"{span.context.span_id:016x}"
            span.set_status(StatusCode.OK)
        except subprocess.CalledProcessError:
            raise
        except Exception as exc:
            span.set_attribute("error.type", type(exc).__name__)
            span.set_attribute("error", str(exc))
            span.set_status(StatusCode.ERROR, str(exc))
            raise
