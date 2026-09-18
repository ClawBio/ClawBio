#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import statistics
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

SKILL_DIR = Path(__file__).resolve().parent
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)
REQUIRED_COLUMNS = {
    "sample_id",
    "batch",
    "total_reads",
    "mapped_pct",
    "duplicate_pct",
    "mitochondrial_pct",
    "contamination_pct",
    "complexity_score",
}
NUMERIC_COLUMNS = REQUIRED_COLUMNS - {"sample_id", "batch"}
OPTIONAL_NUMERIC_COLUMNS = {"fingerprint_match_pct"}
OPTIONAL_TEXT_COLUMNS = {"expected_sex", "observed_sex"}


def load_metrics(path: Path) -> list[dict[str, str | float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input is missing required columns: {sorted(missing)}")
        rows: list[dict[str, str | float]] = []
        fieldnames = set(reader.fieldnames or [])
        for raw in reader:
            row: dict[str, str | float] = {"sample_id": raw["sample_id"], "batch": raw["batch"]}
            for column in sorted(NUMERIC_COLUMNS):
                try:
                    row[column] = float(raw[column])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Column {column} must be numeric") from exc
            for column in sorted(OPTIONAL_NUMERIC_COLUMNS & fieldnames):
                value = raw.get(column, "")
                if value != "":
                    try:
                        row[column] = float(value)
                    except (TypeError, ValueError) as exc:
                        raise ValueError(f"Column {column} must be numeric") from exc
            for column in sorted(OPTIONAL_TEXT_COLUMNS & fieldnames):
                row[column] = raw.get(column, "")
            rows.append(row)
    if not rows:
        raise ValueError("Input contains no samples")
    return rows


def _median_abs_deviation(values: list[float]) -> float:
    median = statistics.median(values)
    return statistics.median([abs(value - median) for value in values]) or 1.0


def analyze_records(rows: list[dict[str, str | float]]) -> dict:
    read_values = [float(row["total_reads"]) for row in rows]
    read_median = statistics.median(read_values)
    read_mad = _median_abs_deviation(read_values)
    samples = []
    issue_counts: dict[str, int] = {}
    for row in rows:
        issues: list[str] = []
        expected_sex = str(row.get("expected_sex", "")).strip().lower()
        observed_sex = str(row.get("observed_sex", "")).strip().lower()
        if expected_sex and observed_sex and expected_sex not in {"unknown", "na"} and observed_sex not in {"unknown", "na"} and expected_sex != observed_sex:
            issues.append("sex_mismatch")
        if "fingerprint_match_pct" in row and float(row["fingerprint_match_pct"]) < 95.0:
            issues.append("identity_mismatch")
        if float(row["complexity_score"]) < 0.60 or float(row["duplicate_pct"]) > 35.0:
            issues.append("low_complexity")
        if float(row["contamination_pct"]) > 5.0:
            issues.append("contamination")
        z_like = (read_median - float(row["total_reads"])) / read_mad
        if z_like > 2.5 or float(row["mapped_pct"]) < 80.0 or float(row["mitochondrial_pct"]) > 15.0:
            issues.append("batch_shift")
        if not issues and float(row["mapped_pct"]) < 90.0:
            issues.append("mapping_drop")
        dominant = issues[0] if issues else "none"
        for issue in issues:
            issue_counts[issue] = issue_counts.get(issue, 0) + 1
        samples.append(
            {
                **row,
                "read_depth_delta_pct": round((float(row["total_reads"]) - read_median) / read_median * 100.0, 2),
                "issues": issues,
                "dominant_issue": dominant,
                "status": "flagged" if issues else "pass",
            }
        )
    flagged = [sample for sample in samples if sample["status"] == "flagged"]
    return {
        "skill": "sample-qc-triage",
        "summary": {
            "sample_count": len(samples),
            "flagged_count": len(flagged),
            "median_reads": read_median,
            "issue_counts": issue_counts,
        },
        "samples": samples,
        "disclaimer": DISCLAIMER,
    }


def _write_csv(path: Path, rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    fieldnames = [
        "sample_id", "batch", "status", "dominant_issue", "expected_sex", "observed_sex",
        "fingerprint_match_pct", "total_reads", "mapped_pct", "duplicate_pct", "mitochondrial_pct",
        "contamination_pct", "complexity_score", "read_depth_delta_pct",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({field: row.get(field, "") for field in fieldnames})


def repro_command(input_path: Path, output_dir: Path, demo: bool) -> ReproCommand:
    args: list[str | ReproPath] = ["--demo"] if demo else [
        "--input",
        ReproPath(input_path, "repo_root" if input_path.is_relative_to(_PROJECT_ROOT) else "auto"),
    ]
    args += ["--output", ReproPath(output_dir, "output_dir")]
    return ReproCommand(
        script_path=Path("skills/sample-qc-triage/sample_qc_triage.py"),
        args=args,
        comment="Reproduce this sample-qc-triage run",
    )


def write_repro_bundle(input_path: Path, output_dir: Path, demo: bool) -> None:
    write_environment_yml(
        output_dir,
        env_name="clawbio-sample-qc-triage",
        pip_deps=[],
        python_version="3.11",
    )
    write_portable_commands_sh(
        output_dir,
        repro_command(input_path, output_dir, demo),
        repo_root=_PROJECT_ROOT,
    )
    write_checksums(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            output_dir / "tables" / "sample_flags.csv",
        ],
        output_dir,
        anchor=output_dir,
    )


def write_outputs(result: dict, input_path: Path, output_dir: Path, demo: bool) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: output directory already exists and files may be overwritten: {output_dir}", file=sys.stderr)
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)
    _write_csv(output_dir / "tables" / "sample_flags.csv", result["samples"])
    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    rows = [
        "# Sample QC Triage Report",
        "",
        f"**Input**: `{input_path}`",
        f"**Mode**: {'Synthetic demo data' if demo else 'User-provided local data'}",
        f"**Samples**: {result['summary']['sample_count']}",
        f"**Flagged samples**: {result['summary']['flagged_count']}",
        "",
        "| Sample | Batch | Status | Dominant issue | Sex check | Fingerprint % | Mapped % | Duplicate % | Contamination % |",
        "|---|---|---|---|---|---:|---:|---:|---:|",
    ]
    for sample in result["samples"]:
        sex_check = ""
        if sample.get("expected_sex") or sample.get("observed_sex"):
            sex_check = f"{sample.get('expected_sex', '')}->{sample.get('observed_sex', '')}"
        rows.append(
            f"| {sample['sample_id']} | {sample['batch']} | {sample['status']} | {sample['dominant_issue']} | "
            f"{sex_check} | {sample.get('fingerprint_match_pct', '')} | "
            f"{sample['mapped_pct']} | {sample['duplicate_pct']} | {sample['contamination_pct']} |"
        )
    rows.extend([
        "",
        "## Interpretation",
        "",
        "Flagged samples are candidates for local review. Identity and sex mismatch flags are operational QC signals, not forensic proof.",
        "",
        DISCLAIMER,
        "",
    ])
    (output_dir / "report.md").write_text("\n".join(rows), encoding="utf-8")
    write_repro_bundle(input_path, output_dir, demo)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Sample QC Triage")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("sample_qc_triage_out"))
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args(argv)
    input_path = SKILL_DIR / "demo_qc_metrics.csv" if args.demo else args.input
    if input_path is None:
        parser.error("--input is required unless --demo is used")
    try:
        result = analyze_records(load_metrics(input_path))
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    write_outputs(result, input_path, args.output, args.demo)
    print(f"Sample QC Triage wrote {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
