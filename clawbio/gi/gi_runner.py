"""Shared CLI runner for the six gi-* skills.

Each skill's ``gi_<task>.py`` is a ~20-line config that calls
``run_skill(task=..., demo_path=..., async_mode=...)``. The runner handles
arg parsing, FASTA → predict → ``{data, meta}``, and writes
``report.md`` + ``result.json`` + ``reproducibility/`` per ClawBio
convention.
"""

from __future__ import annotations

import argparse
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, Optional

import requests

from clawbio.gi.gi_client import Client, FastaError, GIError, read_fasta

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical "
    "device and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# --- Sequence-length contract -------------------------------------------
#
# Every task has its own admission floor, published as ``minLength`` on that
# task's request schema and enforced before any model loads. The values below
# are a LOCAL MIRROR so a too-short FASTA costs no API key, no request and no
# 422; the authority is the served schema at
# ``GET https://api.genomicintelligence.ai/v1/openapi.json``. Re-read it if a
# rejection here disagrees with the server.
#
# Floors are admission control, not biology: they are the strictest thing any
# model for that task needs, and they say nothing about whether the model saw
# real sequence. See TASK_CONTEXT_WINDOW_BP below.
TASK_MIN_BP: Dict[str, int] = {
    "promoter": 300,
    "splice": 100,
    "enhancer": 50,
    "chromatin": 200,
    "annotation": 1000,
    "expression": 9198,
}
REQUEST_MAX_BP = 500_000

# ``sequence_name`` is a display-only field, ``maxLength: 128`` on every
# per-task request schema.
SEQUENCE_NAME_MAX_CHARS = 128

# Extra operator-facing context appended to a too-short rejection.
TASK_MIN_BP_HINT: Dict[str, str] = {
    "expression": (
        "the model scores exactly one TSS-centred window of that size "
        "(TSS ± 4599). Submit at least a full window, or submit a longer "
        "locus and pass --tss-index."
    ),
    "annotation": "the gene finder needs a genomic region, not a single exon.",
}

# The default model's own sliding window, from ``bio_spec.context_window_bp``
# on ``GET /v1/tasks/{task}/models``. A request above the floor but shorter
# than this is ACCEPTED and scored — against a window padded out to the
# context window. So the floor admits, and this is what tells you whether the
# model saw real sequence. ``None`` = no sliding window (annotation and
# expression report ``context_window_bp: null``).
TASK_CONTEXT_WINDOW_BP: Dict[str, Optional[int]] = {
    "promoter": 2000,  # the default model; 300 bp-context models exist too
    "splice": 15000,
    "enhancer": 249,
    "chromatin": 1000,
    "annotation": None,
    "expression": None,
}

# Expression scores exactly one TSS-centred 9,198 bp window (radius 4,599),
# cut server-side; a longer locus is allowed up to REQUEST_MAX_BP provided
# ``tss_index`` says where to cut. Unlike the other tasks, expression does not
# pad — its floor and its window are the same number.
EXPRESSION_WINDOW_BP = TASK_MIN_BP["expression"]
EXPRESSION_TSS_RADIUS = EXPRESSION_WINDOW_BP // 2  # 4599
EXPRESSION_MAX_BP = REQUEST_MAX_BP  # kept as an alias; the cap is task-wide


def validate_sequence_length(task: str, sequence: str) -> Optional[str]:
    """Check a sequence against the task's published length bounds.

    Returns ``None`` when the request would be accepted, otherwise the
    operator-facing reason. Both bounds are ``422 validation_failed`` at the
    API (over-max is *not* a 413 — 413 is the 16 MiB raw-body cap), and both
    are counted on the whitespace-stripped nucleotide string.
    """
    n = len(sequence)
    floor = TASK_MIN_BP.get(task)
    if floor is not None and n < floor:
        hint = TASK_MIN_BP_HINT.get(task)
        msg = f"sequence is {n:,} bp; {task} needs at least {floor} bp"
        return f"{msg} — {hint}" if hint else msg
    if n > REQUEST_MAX_BP:
        return f"sequence is {n:,} bp; the maximum is {REQUEST_MAX_BP} bp"
    return None


def regime_note(task: str, sequence_length: int) -> Optional[str]:
    """Warn when the submission is in range but shorter than the model's window.

    Not an error: the API accepts it and returns a score. The model simply
    sees the sequence padded out to its context window, so the result is a
    padded-window score rather than a full-context one.
    """
    context = TASK_CONTEXT_WINDOW_BP.get(task)
    if context is None or sequence_length >= context:
        return None
    return (
        f"sequence is {sequence_length:,} bp, shorter than the default {task} "
        f"model's {context:,} bp context window — the API accepts and scores it, "
        f"but against a window padded out to {context:,} bp. Treat the result as "
        f"a padded-window score. Check bio_spec.context_window_bp on "
        f"GET /v1/tasks/{task}/models if you passed --model."
    )


def validate_expression_input(sequence: str, tss_index: Optional[int]) -> Optional[str]:
    """Check an expression submission against the API contract locally.

    Returns ``None`` when the request would be accepted, otherwise the
    operator-facing reason. Called before the client is built so an
    off-window FASTA costs no API key, no request, and no 422.

    ``sequence`` must already be the whitespace-stripped nucleotide string
    (``read_fasta`` returns exactly that) — the API counts length and
    interprets ``tss_index`` on the stripped string, not on file characters.

    The length bounds come from :func:`validate_sequence_length`; what is
    expression-specific is the ``tss_index`` rule.
    """
    n = len(sequence)
    problem = validate_sequence_length("expression", sequence)
    if problem:
        return problem
    if tss_index is None:
        if n != EXPRESSION_WINDOW_BP:
            return (
                f"--tss-index is required unless the sequence is exactly "
                f"{EXPRESSION_WINDOW_BP} bp (got {n:,} bp); it is the 0-based TSS offset "
                f"into the sequence, counted after whitespace is stripped"
            )
        return None
    low, high = EXPRESSION_TSS_RADIUS, n - EXPRESSION_TSS_RADIUS
    if not low <= tss_index <= high:
        return (
            f"--tss-index {tss_index} is outside the allowed range [{low}, {high}] for a "
            f"{n:,} bp sequence — the model needs a full ±{EXPRESSION_TSS_RADIUS} bp window "
            f"around the TSS; submit more flanking sequence"
        )
    return None


def _parse_args(task: str) -> argparse.Namespace:
    p = argparse.ArgumentParser(description=f"ClawBio gi-{task}: {task} prediction via Genomic Intelligence API.")
    p.add_argument("--input", type=Path, dest="input_file", help="Input FASTA (single record)")
    p.add_argument("--output", type=Path, default=Path(f"/tmp/gi-{task}"), help="Output directory")
    p.add_argument("--demo", action="store_true", help="Run with the bundled example FASTA")
    p.add_argument("--model", type=str, default=None, help="Override default model name")
    p.add_argument("--description", type=str, default=None, help="Cell type / assay context (required by gi-expression; ignored by other tasks)")
    p.add_argument("--tss-index", type=int, default=None, dest="tss_index", help="0-based TSS offset into the sequence, counted after whitespace is stripped (gi-expression only; required unless the sequence is exactly 9198 bp)")
    p.add_argument("--api-key", type=str, default=None, help="Override GI_API_KEY env (otherwise uses env; raises if unset — see each SKILL.md Authentication section)")
    p.add_argument("--base-url", type=str, default=None, help="Override GI_BASE_URL (default: https://api.genomicintelligence.ai)")
    return p.parse_args()


def _resolve_input(args: argparse.Namespace, demo_path: Path, task: str) -> Path:
    if args.demo or args.input_file is None:
        if not demo_path.exists():
            print(f"Error: bundled demo fixture missing at {demo_path}", file=sys.stderr)
            sys.exit(1)
        return demo_path
    if not args.input_file.exists():
        print(f"Error: --input file not found: {args.input_file}", file=sys.stderr)
        sys.exit(1)
    return args.input_file


class ResponseShapeError(RuntimeError):
    """A 2xx body whose nested fields contradict their documented types.

    Distinct from ``GIError``, which covers what the API itself reported. This
    is a well-formed envelope carrying a field the contract says is an object
    or an array and that arrived as something else.
    """


def _as_obj(v: Any, field: str) -> Dict[str, Any]:
    """Read a response field documented as an object.

    ``Client._require_envelope`` guarantees ``data`` is a non-empty object; it
    deliberately does not police per-task fields nested inside it, because it
    is shared by six tasks and must not encode any one task's schema. So the
    checking happens here.

    Absent or null is legitimate — a task that has no ``prediction`` omits it —
    and becomes ``{}``. A field that is *present with the wrong type* is a
    malformed response and is reported as one.

    Two failure modes pull in opposite directions here. The ``x or {}`` idiom
    this replaces handled null and absent but not a truthy wrong type: a
    ``summary`` arriving as a string passed ``or {}`` untouched and then raised
    AttributeError on ``.get``, in the report writer, which runs after
    ``run_skill``'s try/except has closed — a traceback that reads as a client
    bug. Substituting ``{}`` for it instead fixes the traceback and creates
    something worse: a zero-valued report written to disk and an OK line on
    stderr, so a bad response is indistinguishable from a real prediction of
    nothing. Raising a typed error that ``run_skill`` turns into the same
    diagnostic it gives any other malformed response is neither.
    """
    if v is None:
        return {}
    if not isinstance(v, dict):
        raise ResponseShapeError(f"{field} should be an object, got {type(v).__name__}")
    return v


def _as_objs(v: Any, field: str) -> list:
    """Same, for a field documented as an array of objects.

    A truthy non-list (a bare string) is iterable, so ``or []`` let it through
    and the row loop iterated its characters; non-object elements fail the same
    way. Neither is silently dropped — an array whose elements are the wrong
    type is a malformed response, and a report missing rows it should have had
    is exactly the silent wrong answer this is here to prevent.
    """
    if v is None:
        return []
    if not isinstance(v, list):
        raise ResponseShapeError(f"{field} should be an array, got {type(v).__name__}")
    for i, x in enumerate(v):
        if not isinstance(x, dict):
            raise ResponseShapeError(
                f"{field}[{i}] should be an object, got {type(x).__name__}"
            )
    return v


def _as_bounds(v: Any, field: str) -> Optional[list]:
    """Same, for a field documented as a two-element [start, end) pair.

    ``scored_window`` was read positionally with no guard, inside the very
    block that exists to refuse wrong-typed response fields: a dict is truthy,
    so ``window[0]`` raised an uncaught ``KeyError(0)`` past ``run_skill``'s
    handler rather than the named diagnostic every other field gets. A short
    list is the same hazard one index along.
    """
    if v is None:
        return None
    if not isinstance(v, (list, tuple)):
        raise ResponseShapeError(f"{field} should be an array, got {type(v).__name__}")
    if len(v) != 2:
        raise ResponseShapeError(f"{field} should have 2 elements, got {len(v)}")
    for i, x in enumerate(v):
        if not isinstance(x, int) or isinstance(x, bool):
            raise ResponseShapeError(
                f"{field}[{i}] should be an integer, got {type(x).__name__}"
            )
    return list(v)


def _summarize(task: str, body: Dict[str, Any]) -> Dict[str, Any]:
    """Pick the most useful headline numbers per task from `data`."""
    data = _as_obj(body.get("data"), "data")
    summary = _as_obj(data.get("summary"), "data.summary")
    out: Dict[str, Any] = {"task": task, "model": data.get("model")}
    if task == "promoter":
        out["promoter_windows"] = summary.get("promoter_windows")
        out["total_windows"] = summary.get("total_windows")
        out["regions"] = _as_objs(data.get("regions"), "data.regions")
    elif task == "splice":
        out["sites_found"] = summary.get("total_sites", summary.get("sites_found"))
        out["donor_sites"] = summary.get("donor_sites")
        out["acceptor_sites"] = summary.get("acceptor_sites")
        out["sites"] = _as_objs(data.get("sites"), "data.sites")
    elif task == "enhancer":
        out["windows_processed"] = summary.get("total_windows", summary.get("windows_processed"))
        out["dev_score_max"] = summary.get("dev_score_max")
        out["hk_score_max"] = summary.get("hk_score_max")
    elif task == "chromatin":
        out["windows_processed"] = summary.get("total_windows", summary.get("windows_processed"))
        out["total_annotations"] = summary.get("total_annotations")
    elif task == "expression":
        pred = _as_obj(data.get("prediction"), "data.prediction")
        out["log_tpm"] = pred.get("expression_log_tpm")
        out["tpm"] = pred.get("expression_tpm")
        # Windowing provenance — the API cuts the scored 9,198 bp window
        # itself, so this is the only way to confirm it cut where you meant.
        inp = _as_obj(data.get("input"), "data.input")
        out["tss_index"] = inp.get("tss_index")
        out["scored_window"] = inp.get("scored_window")
        # The submitted length comes from meta.sequence_length, which every
        # response carries and which the revision and the consumer pins cover.
        # data.input.submitted_sequence_length held the same number in an
        # untyped echo covered by nothing, and is being removed; reading it
        # would have dropped the ", of N bp submitted" clause silently, since
        # the report renders it under an isinstance guard.
        meta = _as_obj(body.get("meta"), "meta")
        submitted = meta.get("sequence_length")
        if submitted is None:
            # Older responses carried the number only in the echo. Falling
            # back keeps the clause through either deploy order.
            submitted = inp.get("submitted_sequence_length")
        out["submitted_sequence_length"] = submitted
    elif task == "annotation":
        out["transcripts_found"] = summary.get("total_transcripts", summary.get("transcripts_found"))
        out["transcripts"] = _as_objs(data.get("transcripts"), "data.transcripts")
    out["raw_summary"] = summary
    return out


def _write_report(task: str, summary: Dict[str, Any], body: Dict[str, Any], output_dir: Path, input_path: Path, sequence_name: str, sequence_length: int, elapsed_ms: float, rate_limit: Optional[Dict[str, str]] = None) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    # ``rate_limit`` is the response's own RateLimit-* headers, which every
    # gi SKILL.md tells the reader to consult for the live per-key allowance
    # and which nothing here used to keep: the parsed body does not carry
    # them, so writing only ``full_response`` discarded them.
    (output_dir / "result.json").write_text(json.dumps(
        {"summary": summary, "full_response": body, "rate_limit": rate_limit or {}},
        indent=2,
    ))

    meta = _as_obj(body.get("meta"), "meta")
    model = summary.get("model") or "—"
    lines = [
        f"# gi-{task} report",
        "",
        f"- **Sequence**: `{sequence_name}` ({sequence_length:,} bp)",
        f"- **Input file**: `{input_path}`",
        f"- **Model**: `{model}`",
        f"- **Inference time**: {meta.get('inference_time_ms', elapsed_ms):.0f} ms",
        f"- **Request ID**: `{meta.get('request_id', '—')}`",
        f"- **Generated**: {datetime.now(timezone.utc).isoformat(timespec='seconds')}",
        "",
        "## Headline result",
        "",
    ]
    if task == "promoter":
        lines.append(f"- Promoter windows: **{summary.get('promoter_windows', 0)}** / {summary.get('total_windows', 0)} total")
        regions = _as_objs(summary.get("regions"), "data.regions")
        if regions:
            lines.append("")
            lines.append("| Window | Start | End | Probability |")
            lines.append("|---|---|---|---|")
            for r in regions[:20]:
                lines.append(f"| {r.get('window_index','-')} | {r.get('start','-')} | {r.get('end','-')} | {r.get('probability','-'):.3f} |" if isinstance(r.get('probability'), (int, float)) else f"| {r.get('window_index','-')} | {r.get('start','-')} | {r.get('end','-')} | {r.get('probability','-')} |")
    elif task == "splice":
        lines.append(f"- Splice sites found: **{summary.get('sites_found') or 0}** ({summary.get('donor_sites') or 0} donor + {summary.get('acceptor_sites') or 0} acceptor)")
        sites = _as_objs(summary.get("sites"), "data.sites")[:20]
        if sites:
            lines.append("")
            # Field names are the API's: name / start / end / site_type / score.
            # This table read position/kind/probability, which the API has never
            # returned, so three of four columns rendered as "-" in every report.
            # Span rather than a single position: start/end bound one
            # variable-width tokenizer token and the junction lies inside it, so
            # printing one endpoint states a boundary the response does not give.
            lines.append("| Site | Span (bp) | Type | Strand | Score |")
            lines.append("|---|---|---|---|---|")
            for s in sites:
                start, end = s.get("start"), s.get("end")
                span = f"{start}-{end}" if start is not None and end is not None else "-"
                score = s.get("score")
                score = f"{score:.4f}" if isinstance(score, (int, float)) else (score or "-")
                lines.append(
                    f"| {s.get('name','-')} | {span} | {s.get('site_type','-')} | "
                    f"{s.get('strand','-')} | {score} |"
                )
    elif task == "enhancer":
        lines.append(f"- Windows processed: **{summary.get('windows_processed') or 0}**")
        dev = summary.get("dev_score_max"); hk = summary.get("hk_score_max")
        if dev is not None:
            lines.append(f"- Max developmental-enhancer score: **{dev:.3f}**" if isinstance(dev, (int, float)) else f"- Max developmental-enhancer score: **{dev}**")
        if hk is not None:
            lines.append(f"- Max housekeeping-enhancer score: **{hk:.3f}**" if isinstance(hk, (int, float)) else f"- Max housekeeping-enhancer score: **{hk}**")
    elif task == "chromatin":
        lines.append(f"- Windows processed: **{summary.get('windows_processed') or 0}**")
        lines.append(f"- Total annotations across all tracks: **{summary.get('total_annotations') or 0}**")
    elif task == "expression":
        log_tpm = summary.get("log_tpm")
        tpm = summary.get("tpm")
        if log_tpm is not None:
            lines.append(f"- Predicted expression: **{log_tpm:.4f} log(TPM+1)**" + (f" ≈ {tpm:.2f} TPM" if isinstance(tpm, (int, float)) else ""))
        else:
            lines.append("- See `result.json` for the full prediction payload.")
        window = _as_bounds(summary.get("scored_window"), "data.input.scored_window")
        if window:
            submitted = summary.get("submitted_sequence_length")
            lines.append(
                f"- Scored window: **[{window[0]}, {window[1]})** (TSS index {summary.get('tss_index')}"
                + (f", of {submitted:,} bp submitted" if isinstance(submitted, int) else "")
                + ") — check this is the window you meant; an offset that is merely wrong still returns a confident result."
            )
    elif task == "annotation":
        lines.append(f"- Transcripts found: **{summary.get('transcripts_found') or 0}**")
        tx = _as_objs(summary.get("transcripts"), "data.transcripts")[:20]
        if tx:
            lines.append("")
            lines.append("| Transcript | Start | End | Strand |")
            lines.append("|---|---|---|---|")
            for t in tx:
                lines.append(f"| {t.get('transcript_id','-')} | {t.get('start','-')} | {t.get('end','-')} | {t.get('strand','-')} |")

    lines += [
        "",
        "## Reproducibility",
        "",
        f"- `reproducibility/command.sh` — exact invocation",
        f"- `result.json` — the full `{{data, meta}}` response, the summary this report was built from, and the `RateLimit-*` headers that response carried",
        "",
        "## API",
        "",
        f"`POST /v1/tasks/{task}/predict` on `https://api.genomicintelligence.ai` — see <https://docs.genomicintelligence.ai>.",
        "",
        "---",
        "",
        f"_{DISCLAIMER}_",
        "",
    ]
    (output_dir / "report.md").write_text("\n".join(lines))

    repro = output_dir / "reproducibility"
    repro.mkdir(exist_ok=True)
    cmd = f"python skills/gi-{task}/gi_{task.replace('-', '_')}.py --input {input_path} --output {output_dir}\n"
    (repro / "command.sh").write_text("#!/usr/bin/env bash\nset -euo pipefail\n" + cmd)
    (repro / "command.sh").chmod(0o755)
    (repro / "environment.json").write_text(json.dumps({
        "skill": f"gi-{task}",
        "skill_version": "0.1.0",
        "api_base_url": os.environ.get("GI_BASE_URL", "https://api.genomicintelligence.ai"),
        "model": summary.get("model"),
        "request_id": meta.get("request_id"),
        "timestamp": datetime.now(timezone.utc).isoformat(timespec="seconds"),
    }, indent=2))


def run_skill(*, task: str, demo_path: Path, async_mode: bool = False, default_model: Optional[str] = None, default_options: Optional[Dict[str, Any]] = None) -> int:
    """Skill entry-point. Returns process exit code."""
    args = _parse_args(task)
    input_path = _resolve_input(args, demo_path, task)
    try:
        sequence_name, sequence = read_fasta(input_path)
    except FastaError as e:
        print(f"[gi-{task}] invalid input — {e}", file=sys.stderr)
        return 1
    if not sequence:
        print(f"Error: parsed an empty sequence from {input_path}", file=sys.stderr)
        return 1
    # ``sequence_name`` is capped at 128 chars in every request schema; a long
    # FASTA header would otherwise be a 422 over a display-only field.
    if len(sequence_name) > SEQUENCE_NAME_MAX_CHARS:
        sequence_name = sequence_name[:SEQUENCE_NAME_MAX_CHARS]
        print(f"Warning: FASTA name truncated to {SEQUENCE_NAME_MAX_CHARS} characters (API limit)", file=sys.stderr)
    # Length gate first, for every task: each has its own published floor, and
    # rejecting locally costs no API key, no request and no 422.
    if task == "expression":
        problem = validate_expression_input(sequence, args.tss_index)
    else:
        problem = validate_sequence_length(task, sequence)
        if args.tss_index is not None:
            print(f"Warning: --tss-index applies only to gi-expression; ignoring it for {task}", file=sys.stderr)
    if problem:
        print(f"Error: {problem}", file=sys.stderr)
        return 1
    note = regime_note(task, len(sequence))
    if note:
        print(f"Warning: {note}", file=sys.stderr)

    model = args.model or default_model
    options: Dict[str, Any] = dict(default_options or {})
    if args.description is not None:
        # ``options`` is typed and closed per task (additionalProperties:
        # false). Only expression declares ``description``; forwarding it
        # anywhere else is a hard 422 extra_forbidden, not a silent ignore.
        if task == "expression":
            options["description"] = args.description
        else:
            print(f"Warning: --description applies only to gi-expression; ignoring it for {task}", file=sys.stderr)

    tss_index = args.tss_index if task == "expression" else None

    # Built last: everything above is local validation, so a bad request never
    # needs a key.
    client = Client(api_key=args.api_key, base_url=args.base_url)

    print(f"[gi-{task}] sequence_name={sequence_name} length={len(sequence):,} bp model={model or 'default'} mode={'async' if async_mode else 'sync'}", file=sys.stderr)
    started = time.monotonic()
    try:
        if async_mode:
            job_id = client.submit_async(task, sequence=sequence, sequence_name=sequence_name, model=model, options=options or None, tss_index=tss_index)
            print(f"[gi-{task}] submitted job_id={job_id}", file=sys.stderr)
            def _progress(p: Dict[str, Any]) -> None:
                pct = p.get("percent")
                msg = p.get("message", "")
                if pct is not None:
                    print(f"  {pct:>3}% {msg}", file=sys.stderr)
            body = client.wait_for_job(job_id, on_progress=_progress)
        else:
            body = client.predict(task, sequence=sequence, sequence_name=sequence_name, model=model, options=options or None, tss_index=tss_index)
    except GIError as e:
        print(f"[gi-{task}] API error: {e}", file=sys.stderr)
        return 2
    except requests.RequestException as e:
        # Connection refused, DNS failure, TLS error, read timeout — routine
        # for a hosted service and not the caller's bug. Surface a diagnostic
        # rather than a traceback that reads like a client defect.
        print(f"[gi-{task}] network error reaching the API: {type(e).__name__}: {e}", file=sys.stderr)
        return 2
    except TimeoutError as e:
        # Raised by wait_for_job when a job outlives its deadline.
        print(f"[gi-{task}] timed out waiting for the job: {e}", file=sys.stderr)
        return 2
    except KeyError as e:
        # A 2xx whose body is missing a field we index (e.g. data.job_id on an
        # async submit). Malformed upstream response, not a usage error.
        print(f"[gi-{task}] unexpected API response shape: missing {e}", file=sys.stderr)
        return 2
    elapsed_ms = (time.monotonic() - started) * 1000.0
    # The summary and the report writer run outside the block above, so their
    # own view of a malformed response needs its own handler. Without one a
    # wrong-typed nested field either raised a traceback or — once the helpers
    # coerced it — wrote a zero-valued report and printed the OK line. Both are
    # worse than exiting 2 with the offending field named.
    try:
        summary = _summarize(task, body)
        # getattr because gi_runner.Client is the monkeypatch seam the suites
        # replace, and a test double is not obliged to carry the attribute.
        _write_report(task, summary, body, args.output, input_path, sequence_name, len(sequence), elapsed_ms, rate_limit=getattr(client, "last_rate_limit", None))
    except ResponseShapeError as e:
        print(f"[gi-{task}] unexpected API response shape: {e}", file=sys.stderr)
        return 2
    print(f"[gi-{task}] OK — wrote {args.output}/report.md ({elapsed_ms:.0f} ms wall)", file=sys.stderr)
    return 0
