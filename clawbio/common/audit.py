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
from contextvars import ContextVar
from datetime import datetime, timezone
from pathlib import Path
from typing import TYPE_CHECKING, List, Sequence

# The SDK is imported inside skill_run and tool_call: skills reach write()
# through report.py, and must not need opentelemetry installed to start.
if TYPE_CHECKING:
    from opentelemetry.sdk.trace import ReadableSpan
    from opentelemetry.sdk.trace.export import SpanExportResult

_DEFAULT_LOG = Path.home() / ".clawbio" / "audit.jsonl"
_REDACTED = "__REDACTED__"


def _resolve_log(log_path: Path | str | None) -> Path:
    """Explicit argument, else CLAWBIO_AUDIT_LOG, else ~/.clawbio/audit.jsonl.

    Read per call rather than at import: a skill calls skill_run without a
    path, so the environment is the only way to redirect the log.
    """
    return Path(log_path or os.environ.get("CLAWBIO_AUDIT_LOG") or _DEFAULT_LOG)


def _traces_url(endpoint: str) -> str:
    """OTel habit is to set the full signal URL, so accept it either way."""
    endpoint = endpoint.rstrip("/")
    return endpoint if endpoint.endswith("/v1/traces") else f"{endpoint}/v1/traces"


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
_TRACER = ContextVar("clawbio.audit.tracer", default=None)


def _ns_to_iso(ns: int) -> str:
    return datetime.fromtimestamp(ns / 1e9, tz=timezone.utc).isoformat()


class _JsonlExporter:
    """Exports OTEL spans as JSONL to a local file.

    Duck-typed rather than subclassing SpanExporter, so the class can be
    defined without importing the SDK; SimpleSpanProcessor only calls
    export() and shutdown().
    """

    def __init__(self, log_path: Path) -> None:
        self._log_path = Path(log_path)

    def export(self, spans: Sequence[ReadableSpan]) -> SpanExportResult:
        from opentelemetry.sdk.trace.export import SpanExportResult

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
        except OSError:
            pass
        return SpanExportResult.SUCCESS

    def shutdown(self) -> None:
        pass


def write(event: str, *, log_path: Path | str | None = None, **kwargs) -> None:
    """Write a simple point-in-time audit record as JSONL."""
    log_path = _resolve_log(log_path)
    entry = {"timestamp": datetime.now(timezone.utc).isoformat(), "event": event, **kwargs}
    try:
        log_path.parent.mkdir(parents=True, exist_ok=True)
        with log_path.open("a", encoding="utf-8") as f:
            f.write(json.dumps(entry, default=str) + "\n")
    except OSError:
        pass


@contextmanager
def skill_run(
    skill: str,
    version: str,
    input_checksum: str = "",
    input_file: str = "",
    output_dir: str = "",
    log_path: Path | str | None = None,
):
    """Span for a skill invocation: a root, or a child of the caller's span
    from TRACEPARENT. Yields the span_id (16-char hex).

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
    from opentelemetry.sdk.resources import OTELResourceDetector, Resource
    from opentelemetry.sdk.trace import TracerProvider
    from opentelemetry.sdk.trace.export import BatchSpanProcessor, SimpleSpanProcessor
    from opentelemetry.sdk.trace.sampling import ALWAYS_ON
    from opentelemetry.trace import StatusCode, get_current_span

    # Fallbacks for a plain CLI run. A caller names its own project and service
    # through OTEL_RESOURCE_ATTRIBUTES / OTEL_SERVICE_NAME. Resource.create() lets
    # our values beat those, so the environment is merged back on top; only what
    # it actually sets, or service.name would become "unknown_service".
    # ALWAYS_ON: a caller's TRACEPARENT can say "not sampled", but the audit
    # record is kept for every run regardless of what the caller exports.
    provider = TracerProvider(sampler=ALWAYS_ON, resource=Resource.create({
        "service.name": "clawbio",
        # Phoenix groups traces by this; without it everything lands in "default".
        "openinference.project.name": "clawbio",
    }).merge(OTELResourceDetector().detect()))
    provider.add_span_processor(SimpleSpanProcessor(_JsonlExporter(_resolve_log(log_path))))
    endpoint = os.environ.get("CLAWBIO_OTLP_ENDPOINT")
    if endpoint:
        # Opt-in, and deliberately not OTEL_EXPORTER_OTLP_ENDPOINT: an org-wide
        # setting would otherwise point spans carrying output paths at whatever
        # collector the org runs, which may not be on this machine.
        try:
            from opentelemetry.exporter.otlp.proto.http.trace_exporter import OTLPSpanExporter
        except ImportError:
            # Tracing is optional; the analysis it describes is not.
            print(
                "clawbio: CLAWBIO_OTLP_ENDPOINT is set but the exporter is not "
                "installed (pip install 'clawbio[otlp]'); spans stay in the "
                "local audit log only.",
                file=sys.stderr,
            )
        else:
            provider.add_span_processor(
                BatchSpanProcessor(
                    # Short timeout: a dead collector must not make shutdown() look
                    # like a hung analysis.
                    OTLPSpanExporter(endpoint=_traces_url(endpoint), timeout=5)
                )
            )
    tracer = provider.get_tracer(__name__)

    # An agent that launched us as a subprocess can only hand over its span
    # through the environment. An open span wins (full-profile nests skill_run in
    # skill_run); the propagator drops a malformed value, leaving a root span.
    parent = None
    if not get_current_span().get_span_context().is_valid and os.environ.get("TRACEPARENT"):
        from opentelemetry.trace.propagation.tracecontext import TraceContextTextMapPropagator

        parent = TraceContextTextMapPropagator().extract({
            "traceparent": os.environ["TRACEPARENT"],
            "tracestate": os.environ.get("TRACESTATE", ""),
        })

    token = _TRACER.set(tracer)
    try:
        # The SDK would otherwise record an exception event and overwrite the
        # status description, both carrying the raw text past the hide flags.
        with tracer.start_as_current_span(
            "skill_run", context=parent, record_exception=False, set_status_on_exception=False
        ) as span:
            span.set_attribute("gen_ai.agent.id", skill)
            span.set_attribute("gen_ai.agent.version", version)
            # Phoenix classifies on this key alone; without it: "unknown".
            span.set_attribute("openinference.span.kind", "CHAIN")
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
        _TRACER.reset(token)
        provider.shutdown()  # flush; a short run exits before the next batch tick


@contextmanager
def tool_call(
    name: str,
    *,
    cmd: List[str] | None = None,
    log_path: Path | str | None = None,
    **attrs,
):
    """Child span for a tool or CLI call.

    Span name: ``execute_tool {name}``
    Pass cmd to run a subprocess and capture its exit code automatically.

    PII warning: ``cmd`` tokens, captured ``stderr`` and ``**attrs`` are written
    to the audit log, and go to the OTLP collector too when the enclosing
    ``skill_run`` has ``CLAWBIO_OTLP_ENDPOINT`` set. OPENINFERENCE_HIDE_INPUTS
    redacts the command and the attrs; either hide flag redacts ``stderr``.
    The keys themselves are never hidden, so callers must still scrub file
    paths, sample IDs, and any patient-identifiable values before passing
    them here.
    """
    tracer = _TRACER.get()
    if tracer is None:
        yield None
        return
    from opentelemetry.trace import StatusCode

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
                # These are arguments to the command, so the input flag covers
                # them: without this they were the one value written raw.
                span.set_attribute(k, _hide(str(v), "OPENINFERENCE_HIDE_INPUTS"))
        elif attrs:
            # No subprocess, so no command to show: the arguments are the input.
            span.set_attribute("input.value", _hide(
                json.dumps(attrs, default=str), "OPENINFERENCE_HIDE_INPUTS"
            ))
            span.set_attribute("input.mime_type", "application/json")
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
