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

from clawbio.gi.gi_client import Client, GIError, read_fasta

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical "
    "device and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# /v1/tasks/expression/predict scores exactly one TSS-centred window of
# 9,198 bp (radius 4,599). Longer submissions are allowed up to 500 kbp
# provided ``tss_index`` says where to cut; everything else is a 422.
# Mirrors EXPRESSION_MIN_BP / EXPRESSION_MAX_BP / EXPRESSION_TSS_RADIUS
# in the API. Contract: https://docs.genomicintelligence.ai
EXPRESSION_WINDOW_BP = 9198
EXPRESSION_TSS_RADIUS = EXPRESSION_WINDOW_BP // 2  # 4599
EXPRESSION_MAX_BP = 500_000


def validate_expression_input(sequence: str, tss_index: Optional[int]) -> Optional[str]:
    """Check an expression submission against the API contract locally.

    Returns ``None`` when the request would be accepted, otherwise the
    operator-facing reason. Called before the client is built so an
    off-window FASTA costs no API key, no request, and no 422.

    ``sequence`` must already be the whitespace-stripped nucleotide string
    (``read_fasta`` returns exactly that) — the API counts length and
    interprets ``tss_index`` on the stripped string, not on file characters.
    """
    n = len(sequence)
    if n < EXPRESSION_WINDOW_BP:
        return (
            f"sequence is {n:,} bp; expression needs at least {EXPRESSION_WINDOW_BP} bp — "
            f"the model scores exactly one TSS-centred window of that size "
            f"(TSS ± {EXPRESSION_TSS_RADIUS}). Submit at least a full window, or submit a "
            f"longer locus and pass --tss-index."
        )
    if n > EXPRESSION_MAX_BP:
        return f"sequence is {n:,} bp; the maximum is {EXPRESSION_MAX_BP} bp"
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


def _summarize(task: str, body: Dict[str, Any]) -> Dict[str, Any]:
    """Pick the most useful headline numbers per task from `data`."""
    data = body.get("data") or {}
    summary = data.get("summary") or {}
    out: Dict[str, Any] = {"task": task, "model": data.get("model")}
    if task == "promoter":
        out["promoter_windows"] = summary.get("promoter_windows")
        out["total_windows"] = summary.get("total_windows")
        out["regions"] = data.get("regions") or []
    elif task == "splice":
        out["sites_found"] = summary.get("total_sites", summary.get("sites_found"))
        out["donor_sites"] = summary.get("donor_sites")
        out["acceptor_sites"] = summary.get("acceptor_sites")
        out["sites"] = data.get("sites") or []
    elif task == "enhancer":
        out["windows_processed"] = summary.get("total_windows", summary.get("windows_processed"))
        out["dev_score_max"] = summary.get("dev_score_max")
        out["hk_score_max"] = summary.get("hk_score_max")
    elif task == "chromatin":
        out["windows_processed"] = summary.get("total_windows", summary.get("windows_processed"))
        out["total_annotations"] = summary.get("total_annotations")
    elif task == "expression":
        pred = data.get("prediction") or {}
        out["log_tpm"] = pred.get("expression_log_tpm")
        out["tpm"] = pred.get("expression_tpm")
        # Windowing provenance — the API cuts the scored 9,198 bp window
        # itself, so this is the only way to confirm it cut where you meant.
        inp = data.get("input") or {}
        out["tss_index"] = inp.get("tss_index")
        out["scored_window"] = inp.get("scored_window")
        out["submitted_sequence_length"] = inp.get("submitted_sequence_length")
    elif task == "annotation":
        out["transcripts_found"] = summary.get("total_transcripts", summary.get("transcripts_found"))
        out["transcripts"] = data.get("transcripts") or []
    out["raw_summary"] = summary
    return out


def _write_report(task: str, summary: Dict[str, Any], body: Dict[str, Any], output_dir: Path, input_path: Path, sequence_name: str, sequence_length: int, elapsed_ms: float) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "result.json").write_text(json.dumps({"summary": summary, "full_response": body}, indent=2))

    meta = body.get("meta") or {}
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
        regions = summary.get("regions") or []
        if regions:
            lines.append("")
            lines.append("| Window | Start | End | Probability |")
            lines.append("|---|---|---|---|")
            for r in regions[:20]:
                lines.append(f"| {r.get('window_index','-')} | {r.get('start','-')} | {r.get('end','-')} | {r.get('probability','-'):.3f} |" if isinstance(r.get('probability'), (int, float)) else f"| {r.get('window_index','-')} | {r.get('start','-')} | {r.get('end','-')} | {r.get('probability','-')} |")
    elif task == "splice":
        lines.append(f"- Splice sites found: **{summary.get('sites_found') or 0}** ({summary.get('donor_sites') or 0} donor + {summary.get('acceptor_sites') or 0} acceptor)")
        sites = (summary.get("sites") or [])[:20]
        if sites:
            lines.append("")
            lines.append("| Position | Kind | Strand | Probability |")
            lines.append("|---|---|---|---|")
            for s in sites:
                lines.append(f"| {s.get('position','-')} | {s.get('kind','-')} | {s.get('strand','-')} | {s.get('probability','-')} |")
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
        window = summary.get("scored_window")
        if window:
            submitted = summary.get("submitted_sequence_length")
            lines.append(
                f"- Scored window: **[{window[0]}, {window[1]})** (TSS index {summary.get('tss_index')}"
                + (f", of {submitted:,} bp submitted" if isinstance(submitted, int) else "")
                + ") — check this is the window you meant; an offset that is merely wrong still returns a confident result."
            )
    elif task == "annotation":
        lines.append(f"- Transcripts found: **{summary.get('transcripts_found') or 0}**")
        tx = (summary.get("transcripts") or [])[:20]
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
        f"- `result.json` — full `{{data, meta}}` response from the API",
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
    sequence_name, sequence = read_fasta(input_path)
    if not sequence:
        print(f"Error: parsed an empty sequence from {input_path}", file=sys.stderr)
        return 1
    if task == "expression":
        problem = validate_expression_input(sequence, args.tss_index)
        if problem:
            print(f"Error: {problem}", file=sys.stderr)
            return 1
    elif args.tss_index is not None:
        print(f"Warning: --tss-index applies only to gi-expression; ignoring it for {task}", file=sys.stderr)

    client = Client(api_key=args.api_key, base_url=args.base_url)
    model = args.model or default_model
    options: Dict[str, Any] = dict(default_options or {})
    if args.description is not None:
        options["description"] = args.description

    tss_index = args.tss_index if task == "expression" else None

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
    elapsed_ms = (time.monotonic() - started) * 1000.0
    summary = _summarize(task, body)
    _write_report(task, summary, body, args.output, input_path, sequence_name, len(sequence), elapsed_ms)
    print(f"[gi-{task}] OK — wrote {args.output}/report.md ({elapsed_ms:.0f} ms wall)", file=sys.stderr)
    return 0
