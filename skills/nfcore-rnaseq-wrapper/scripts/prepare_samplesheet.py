#!/usr/bin/env python3
"""DAG step: prepare_samplesheet.

Ports nfcore_rnaseq_wrapper.py's ``_prepare_samplesheet`` dispatch (demo /
self-contained-test-profile / user samplesheet) onto the real, unmodified
``samplesheet_builder.validate_and_normalize_samplesheet`` for the "user"
path, and small verbatim-ported stub writers for the demo/noinput paths.
Writes the normalized samplesheet CSV at
``<output>/reproducibility/samplesheet.{valid,demo,noinput}.csv`` -- the same
path convention as the real wrapper -- plus a JSON summary consumed by the
next DAG step (``preflight``).
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from _isolated_imports import purge_foreign_bare_modules  # noqa: E402

purge_foreign_bare_modules("errors", "schemas", "samplesheet_builder", "args_common")

from args_common import build_args, none_or  # noqa: E402
from errors import SkillError  # noqa: E402
from samplesheet_builder import validate_and_normalize_samplesheet  # noqa: E402


def _write_demo_samplesheet(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"sample,fastq_1,fastq_2,strandedness\n")


def _write_noinput_samplesheet_stub(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"# Input provided by self-contained nf-core test profile (no --input required)\n")


def _final_samplesheet_path(output_dir: Path, *, demo: bool, noinput: bool) -> Path:
    if demo:
        name = "samplesheet.demo.csv"
    elif noinput:
        name = "samplesheet.noinput.csv"
    else:
        name = "samplesheet.valid.csv"
    return output_dir / "reproducibility" / name


def main() -> int:
    parser = argparse.ArgumentParser(description="Validate/normalize the nf-core/rnaseq samplesheet, or write the demo stub.")
    parser.add_argument("--output", required=True)
    parser.add_argument("--demo", required=True)
    parser.add_argument("--input", required=True, help="Samplesheet CSV path, or the sentinel 'none'.")
    parser.add_argument("--skip-alignment", required=True)
    parser.add_argument("--profile", required=True)
    parser.add_argument("--out-summary", required=True, help="Path to write samplesheet_summary.json")
    parsed = parser.parse_args()

    args = build_args(
        output=parsed.output,
        demo=parsed.demo,
        input=parsed.input,
        skip_alignment=parsed.skip_alignment,
        profile=parsed.profile,
    )
    output_dir = Path(args.output).expanduser().resolve()

    try:
        if args.demo:
            normalized_samplesheet = _final_samplesheet_path(output_dir, demo=True, noinput=False)
            _write_demo_samplesheet(normalized_samplesheet)
            summary = {
                "normalized_path": str(normalized_samplesheet),
                "sample_count": 0,
                "sample_names": [],
                "fastq_paths": [],
                "bam_paths": [],
                "unknown_columns": [],
                "strandedness_counts": {},
                "mode": "demo",
            }
        elif getattr(args, "_noinput", False):
            normalized_samplesheet = _final_samplesheet_path(output_dir, demo=False, noinput=True)
            _write_noinput_samplesheet_stub(normalized_samplesheet)
            summary = {
                "normalized_path": str(normalized_samplesheet),
                "sample_count": 0,
                "sample_names": [],
                "fastq_paths": [],
                "bam_paths": [],
                "unknown_columns": [],
                "strandedness_counts": {},
                "mode": "noinput",
            }
        else:
            if not none_or(parsed.input):
                raise SkillError(
                    stage="validation",
                    error_code="MISSING_INPUT",
                    message="An input samplesheet is required unless demo=true or a self-contained test profile is used.",
                    fix="Pass input=<samplesheet.csv>, or demo=true.",
                    details={},
                )
            normalized_samplesheet = _final_samplesheet_path(output_dir, demo=False, noinput=False)
            summary = validate_and_normalize_samplesheet(
                Path(args.input).expanduser().resolve(),
                normalized_samplesheet,
                skip_alignment=args.skip_alignment,
            )
            summary["normalized_path"] = str(summary["normalized_path"])
            summary["fastq_paths"] = [str(p) for p in summary["fastq_paths"]]
            summary["bam_paths"] = [str(p) for p in summary["bam_paths"]]
            summary["mode"] = "user"
    except SkillError as exc:
        out_path = Path(parsed.out_summary)
        out_path.parent.mkdir(parents=True, exist_ok=True)
        out_path.write_text(json.dumps(exc.to_dict(), indent=2) + "\n", encoding="utf-8")
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    out_summary_path = Path(parsed.out_summary)
    out_summary_path.parent.mkdir(parents=True, exist_ok=True)
    out_summary_path.write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(f"[prepare_samplesheet] mode={summary['mode']} samples={summary['sample_count']} -> {normalized_samplesheet}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
