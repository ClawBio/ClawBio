#!/usr/bin/env python3
"""clawbio step: validate_samplesheet.

Parses and normalizes the input samplesheet CSV exactly as a real run would:
resolves FASTQ paths relative to the CSV's parent directory, normalizes
sample-name whitespace, validates required/preset-conditional columns and the
`sample_type`/`feature_type` enums, rejects duplicate rows and inconsistent
repeated-sample metadata, and writes a normalized copy with the header order
nf-core/scrnaseq 4.1.0's `assets/schema_input.json` expects. Ported unchanged
from ClawBio 1.x's ``samplesheet_builder.py``.

Unlike `preflight`, this is a genuine hard-fail boundary: a malformed
samplesheet is a real user error and should stop the DAG here, matching
ClawBio 1.x's behavior (and requiring no Nextflow/Docker to catch).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPTS_DIR = Path(__file__).resolve().parent
if str(_SCRIPTS_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPTS_DIR))

from errors import SkillError
from samplesheet_builder import validate_and_normalize_samplesheet
from wrapper_args import parse_optional_int


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Validate and normalize an nf-core/scrnaseq samplesheet CSV."
    )
    parser.add_argument("--input", required=True)
    parser.add_argument("--preset", required=True)
    parser.add_argument("--expected-cells", default="none")
    parser.add_argument("--out", required=True, help="Path for the normalized samplesheet CSV")
    parser.add_argument("--out-summary", required=True, help="Path for the samplesheet_summary.json")
    args = parser.parse_args()

    try:
        expected_cells_override = parse_optional_int(args.expected_cells)
    except ValueError:
        print(
            json.dumps(
                {
                    "ok": False,
                    "error_code": "INVALID_SAMPLESHEET",
                    "message": f"--expected-cells must be an integer or 'none', got {args.expected_cells!r}",
                },
                indent=2,
            ),
            file=sys.stderr,
        )
        return 1

    input_path = Path(args.input).expanduser().resolve()
    output_path = Path(args.out)

    try:
        summary = validate_and_normalize_samplesheet(
            input_path,
            output_path,
            expected_cells_override=expected_cells_override,
            preset=args.preset,
        )
    except SkillError as exc:
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    Path(args.out_summary).parent.mkdir(parents=True, exist_ok=True)
    Path(args.out_summary).write_text(json.dumps(summary, indent=2, default=str), encoding="utf-8")

    print(
        f"[validate_samplesheet] {summary['sample_count']} sample(s), "
        f"{len(summary['unknown_columns'])} unknown column(s) omitted"
    )
    print(json.dumps(summary, indent=2, default=str))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
