#!/usr/bin/env python3
# nfcore_sarek_wrapper / scripts/samplesheet_validate.py
"""Stage 1 of the pipeline.yaml DAG: validate + normalize a Sarek samplesheet.

Thin CLI wrapper around the faithfully-ported
``samplesheet_builder.validate_and_normalize_samplesheet`` (same directory).
All samplesheet-shape rules (required columns per --step, FASTQ/Spring/BAM/CRAM
mutual exclusivity, tumor/normal pairing, contamination requirements for
varlociraptor, ...) run for real here — this stage needs no Nextflow, Java, or
container runtime, only Python + the samplesheet + any local FASTQ/BAM/CRAM/VCF
paths it references (their *existence* is checked; their content is not).

Demo mode (``--demo true``): skips samplesheet validation entirely and writes
a placeholder normalized-samplesheet/report pair, mirroring 1.x's ``--demo``
behavior of ignoring ``--input`` (the upstream `-profile test` dataset ships
its own remote samplesheet that this port does not fetch or execute — see
SKILL.md).
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

_SCRIPT_DIR = Path(__file__).resolve().parent
if str(_SCRIPT_DIR) not in sys.path:
    sys.path.insert(0, str(_SCRIPT_DIR))

from errors import SkillError
from samplesheet_builder import validate_and_normalize_samplesheet


def _tool_list(raw: str) -> list[str]:
    if not raw or raw.strip().lower() in ("none", ""):
        return []
    return [t.strip().lower() for t in raw.split(",") if t.strip()]


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Validate + normalize a Sarek samplesheet")
    parser.add_argument("--input", required=True, help="Samplesheet path, or 'none' when --demo true")
    parser.add_argument("--step", required=True, help="Sarek --step (mapping/markduplicates/.../annotate)")
    parser.add_argument("--tools", required=True, help="Comma-separated --tools, or 'none'")
    parser.add_argument("--demo", required=True, help="true/false: upstream -profile test mode")
    parser.add_argument("--out-samplesheet", required=True, help="Output path for the normalized CSV")
    parser.add_argument("--out-report", required=True, help="Output path for the samplesheet report JSON")
    args = parser.parse_args(argv)

    demo = args.demo.strip().lower() in ("1", "true", "yes", "on")
    out_samplesheet = Path(args.out_samplesheet)
    out_report = Path(args.out_report)
    out_samplesheet.parent.mkdir(parents=True, exist_ok=True)
    out_report.parent.mkdir(parents=True, exist_ok=True)

    if demo:
        # 1.x --demo clears --input unconditionally: the upstream `-profile test`
        # dataset ships its own samplesheet + references (fetched by Nextflow
        # itself over the network). This port does not execute Nextflow (see
        # SKILL.md), so there is nothing to validate or normalize here; record
        # that fact instead of fabricating sample rows.
        out_samplesheet.write_text(
            "# --demo true: the real wrapper ignores --input and lets the upstream\n"
            "# nf-core/sarek `-profile test` supply its own samplesheet + references\n"
            "# at Nextflow launch time. This port does not execute Nextflow, so no\n"
            "# samplesheet is normalized in demo mode.\n",
            encoding="utf-8",
        )
        report = {
            "mode": "demo",
            "normalized_path": None,
            "step": args.step,
            "sample_count": 0,
            "patient_count": 0,
            "analysis_mode": "germline",
            "pairings": [],
            "rows_by_patient": {},
            "note": (
                "--demo true: samplesheet validation skipped; upstream -profile "
                "test supplies its own dataset at real Nextflow launch time."
            ),
        }
        out_report.write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
        print(f"[samplesheet_validate] --demo true: skipped validation, wrote placeholder report to {out_report}")
        return 0

    try:
        report = validate_and_normalize_samplesheet(
            Path(args.input),
            out_samplesheet,
            step=args.step,
            tools=_tool_list(args.tools),
        )
    except SkillError as exc:
        print(json.dumps(exc.to_dict(), indent=2), file=sys.stderr)
        return 1

    # JSON-serialize: normalized_path is a Path, rows_by_patient values may
    # embed non-JSON-native types (all are already str/int/float/None in
    # practice, but Path objects are not) — cast defensively.
    serializable = dict(report)
    serializable["mode"] = "real"
    serializable["normalized_path"] = str(report["normalized_path"])
    for key in ("fastq_paths", "bam_paths", "cram_paths", "vcf_paths", "spring_paths", "tables"):
        serializable[key] = [str(p) for p in report.get(key, [])]
    out_report.write_text(json.dumps(serializable, indent=2, default=str) + "\n", encoding="utf-8")

    print(
        f"[samplesheet_validate] OK: {report['sample_count']} sample(s), "
        f"{report['patient_count']} patient(s), analysis_mode={report['analysis_mode']!r}"
    )
    for w in serializable.get("unknown_columns") or []:
        print(f"  · unknown column (passed through, not validated): {w}", file=sys.stderr)
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
