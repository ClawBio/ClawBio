#!/usr/bin/env python3
"""clawbio step: preflight.

Runs the ClawBio 1.x preflight validation suite (``preflight.run_preflight``,
ported unchanged) against the curated inputs this DAG exposes: Java and
Nextflow version checks, the execution-backend check for `--profile`
(Docker/conda/Singularity/...), reference/index presence and conflict checks,
protocol-vs-preset compatibility, and the local-first remote-input policy. It
also eagerly creates and validates `output_dir` — the directory a real run
would eventually write results into — exactly as ClawBio 1.x's preflight
does.

**Deliberate departure from ClawBio 1.x**: the 1.x wrapper aborts the whole
process the instant any preflight check fails, specifically so Nextflow is
never invoked against an invalid setup. This step keeps that per-check logic
byte-for-byte (nothing about *what* is validated changes), but never lets a
failure abort the step itself — it always writes preflight_report.json and
exits 0. That is what lets `build_params`/`build_command` still run and
produce the exact command a real run needs, even in this sandbox where
Nextflow is never installed (so the Nextflow version check is *expected* to
fail here — see SKILL.md). ``ok: false`` in the report means what it always
meant: do not run the composed command until every failing check is fixed.
"""
from __future__ import annotations

import argparse
import json
import sys
import traceback
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from errors import SkillError
from preflight import run_preflight
from wrapper_args import build_args, parse_bool, parse_optional_str


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Run ClawBio 1.x's ported preflight checks and report pass/fail without aborting the DAG."
    )
    parser.add_argument("--pipeline-source", required=True, help="pipeline_source.json from resolve_pipeline_source")
    parser.add_argument("--samplesheet-summary", required=True, help="samplesheet_summary.json from validate_samplesheet")
    parser.add_argument("--preset", required=True)
    parser.add_argument("--protocol", default="none")
    parser.add_argument("--genome", default="none")
    parser.add_argument("--fasta", default="none")
    parser.add_argument("--gtf", default="none")
    parser.add_argument("--profile", required=True)
    parser.add_argument("--allow-remote-inputs", default="false")
    parser.add_argument("--skip-cellbender", default="false")
    parser.add_argument("--output-dir", required=True, help="Where a real Nextflow run would write results")
    parser.add_argument("--out", required=True)
    args = parser.parse_args()

    pipeline_source = json.loads(Path(args.pipeline_source).read_text(encoding="utf-8"))
    samplesheet_summary = json.loads(Path(args.samplesheet_summary).read_text(encoding="utf-8"))

    wrapper_args = build_args(
        preset=args.preset,
        protocol=parse_optional_str(args.protocol),
        genome=parse_optional_str(args.genome),
        fasta=parse_optional_str(args.fasta),
        gtf=parse_optional_str(args.gtf),
        profile=args.profile,
        allow_remote_inputs=parse_bool(args.allow_remote_inputs),
        skip_cellbender=parse_bool(args.skip_cellbender),
        output=args.output_dir,
        resume=False,
        demo=False,
    )

    report: dict[str, object]
    try:
        result = run_preflight(
            wrapper_args,
            pipeline_source=pipeline_source,
            samplesheet_summary=samplesheet_summary,
        )
        report = {"ok": True, **result}
        print("[preflight] PASSED — all checks succeeded; the composed command is safe to run.")
    except SkillError as exc:
        report = {"ok": False, "failure": exc.to_dict()}
        print(
            f"[preflight] FAILED at check {exc.error_code} ({exc.stage}): {exc.message}\n"
            f"  fix: {exc.fix}",
            file=sys.stderr,
        )
    except Exception as exc:  # pragma: no cover - defensive; never crash the step
        report = {
            "ok": False,
            "failure": {
                "error_code": "UNEXPECTED_ERROR",
                "message": str(exc),
                "traceback": traceback.format_exc(),
            },
        }
        print(f"[preflight] FAILED unexpectedly: {exc}", file=sys.stderr)

    Path(args.out).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out).write_text(json.dumps(report, indent=2, default=str), encoding="utf-8")
    print(json.dumps(report, indent=2, default=str))
    return 0  # always 0 — see module docstring


if __name__ == "__main__":
    raise SystemExit(main())
