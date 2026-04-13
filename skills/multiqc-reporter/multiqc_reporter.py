#!/usr/bin/env python3
"""
multiqc_reporter.py — ClawBio MultiQC skill
============================================
Wraps the multiqc CLI to aggregate QC reports from any bioinformatics tool
outputs, then generates a ClawBio-style markdown summary from the JSON data.

Usage
-----
    python skills/multiqc-reporter/multiqc_reporter.py \
        --input <dir> [<dir2> ...] --output <report_dir>

    python skills/multiqc-reporter/multiqc_reporter.py --demo --output /tmp/multiqc_demo
"""

from __future__ import annotations

import argparse
import json
import shutil
import subprocess
import sys
import tempfile
from datetime import datetime, timezone
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)


def check_multiqc() -> str:
    path = shutil.which("multiqc")
    if not path:
        print("ERROR: multiqc not found on PATH.\nInstall with:  pip install multiqc", file=sys.stderr)
        sys.exit(1)
    return path


def validate_input_dirs(dirs: list[str]) -> list[Path]:
    paths: list[Path] = []
    for d in dirs:
        p = Path(d)
        if not p.is_dir():
            print(f"ERROR: Input directory not found: {d}", file=sys.stderr)
            sys.exit(1)
        paths.append(p)
    return paths


def run_multiqc(input_dirs: list[Path], output_dir: Path) -> None:
    cmd = ["multiqc"] + [str(d) for d in input_dirs] + ["--outdir", str(output_dir)]
    result = subprocess.run(cmd)
    if result.returncode != 0:
        print("ERROR: multiqc exited with a non-zero status.", file=sys.stderr)
        sys.exit(result.returncode)


def parse_general_stats(output_dir: Path) -> dict[str, dict[str, object]]:
    """Flatten multiqc_data.json report_general_stats_data into {sample: {metric: value}}."""
    stats_file = output_dir / "multiqc_data" / "multiqc_data.json"
    if not stats_file.exists():
        return {}
    raw: dict = json.loads(stats_file.read_text())
    samples: dict[str, dict[str, object]] = {}
    for _tool, tool_samples in raw.get("report_general_stats_data", {}).items():
        if not isinstance(tool_samples, dict):
            continue
        for sample, metrics in tool_samples.items():
            if isinstance(metrics, dict):
                samples.setdefault(sample, {}).update(metrics)
    return samples


def write_report(
    output_dir: Path,
    input_dirs: list[Path],
    samples: dict[str, dict[str, object]],
) -> None:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    input_str = ", ".join(str(d) for d in input_dirs)

    lines: list[str] = [
        "# MultiQC Report",
        "",
        f"**Date**: {now}  ",
        f"**Input directories**: {input_str}  ",
        "",
        "## Per-Sample QC",
        "",
    ]

    if samples:
        all_metrics = sorted({k for s in samples.values() for k in s})
        lines.append("| Sample | " + " | ".join(all_metrics) + " |")
        lines.append("|--------|" + "|".join(["-----"] * len(all_metrics)) + "|")
        for sample in sorted(samples):
            vals = [str(samples[sample].get(m, "—")) for m in all_metrics]
            lines.append(f"| {sample} | " + " | ".join(vals) + " |")
    else:
        lines.append(
            "_No general stats data found in `multiqc_data/multiqc_data.json`. "
            "Check MultiQC logs and that modules produced parseable output._"
        )

    lines += [
        "",
        "## Outputs",
        "",
        "- `multiqc_report.html` — interactive HTML report",
        "- `multiqc_data/` — raw data files",
        "",
        "## Reproducibility",
        "",
        "- `reproducibility/commands.sh` — replay this ClawBio MultiQC run",
        "- `reproducibility/environment.yml` — suggested conda environment",
        "- `reproducibility/checksums.sha256` — key outputs",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
    ]

    (output_dir / "report.md").write_text("\n".join(lines) + "\n")


def make_demo_fastqc_files(dest: Path) -> None:
    """Write minimal FastQC stubs for 3 samples so MultiQC populates general stats."""
    samples = [
        ("SAMPLE_01", 1_000_000, 49, 5.5),
        ("SAMPLE_02", 920_000, 50, 15.0),
        ("SAMPLE_03", 880_000, 48, 7.5),
    ]
    for name, total, pct_gc, pct_dup in samples:
        sample_dir = dest / f"{name}_fastqc"
        sample_dir.mkdir(parents=True, exist_ok=True)
        pct_dedup = 100.0 - pct_dup
        content = (
            f"##FastQC\t0.11.9\n"
            f">>Basic Statistics\tpass\n"
            f"#Measure\tValue\n"
            f"Filename\t{name}.fastq.gz\n"
            f"File type\tConventional base calls\n"
            f"Encoding\tSanger / Illumina 1.9\n"
            f"Total Sequences\t{total}\n"
            f"Sequences flagged as poor quality\t0\n"
            f"Sequence length\t150\n"
            f"%GC\t{pct_gc}\n"
            f">>END_MODULE\n"
            f">>Per base sequence quality\tpass\n"
            f">>END_MODULE\n"
            f">>Per sequence quality scores\tpass\n"
            f">>END_MODULE\n"
            f">>Per base sequence content\tpass\n"
            f">>END_MODULE\n"
            f">>Per sequence GC content\tpass\n"
            f">>END_MODULE\n"
            f">>Per base N content\tpass\n"
            f">>END_MODULE\n"
            f">>Sequence Length Distribution\tpass\n"
            f">>END_MODULE\n"
            f">>Sequence Duplication Levels\tpass\n"
            f"#Total Deduplicated Percentage\t{pct_dedup:.1f}\n"
            f">>END_MODULE\n"
            f">>Overrepresented sequences\tpass\n"
            f">>END_MODULE\n"
            f">>Adapter Content\tpass\n"
            f">>END_MODULE\n"
        )
        (sample_dir / "fastqc_data.txt").write_text(content)


def _multiqc_pip_spec() -> str:
    try:
        from importlib.metadata import version
        return f"multiqc=={version('multiqc')}"
    except Exception:
        return "multiqc"


def _output_paths_to_checksum(output_dir: Path) -> list[Path]:
    paths: list[Path] = []
    for rel in ("report.md", "multiqc_report.html"):
        p = output_dir / rel
        if p.is_file():
            paths.append(p)
    mdir = output_dir / "multiqc_data"
    if mdir.is_dir():
        paths.extend(p for p in sorted(mdir.rglob("*")) if p.is_file())
    return paths


def cli_command_for_repro(output_dir: Path, *, demo: bool, input_dirs: list[Path]) -> str:
    script = str((_PROJECT_ROOT / "skills" / "multiqc-reporter" / "multiqc_reporter.py").resolve())
    parts = [sys.executable, script]
    if demo:
        parts.append("--demo")
    else:
        parts += ["--input"] + [str(p.resolve()) for p in input_dirs]
    parts += ["--output", str(output_dir.resolve())]
    return " ".join(parts)


def write_reproducibility_bundle(output_dir: Path, *, demo: bool, input_dirs: list[Path]) -> None:
    write_commands_sh(output_dir, cli_command_for_repro(output_dir, demo=demo, input_dirs=input_dirs))
    write_environment_yml(output_dir, "clawbio-multiqc-reporter", pip_deps=[_multiqc_pip_spec()], python_version="3.11")
    write_checksums(_output_paths_to_checksum(output_dir), output_dir, anchor=output_dir)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Aggregate QC reports with MultiQC and generate a ClawBio markdown summary"
    )
    parser.add_argument("--input", nargs="+", metavar="DIR",
        help="One or more directories to scan (MultiQC auto-detects tool outputs)")
    parser.add_argument("--output", required=True, metavar="DIR",
        help="Output directory for report.md, HTML report, and data files")
    parser.add_argument("--demo", action="store_true",
        help="Run with synthetic FastQC demo data (no --input required)")
    args = parser.parse_args()

    if not args.demo and not args.input:
        parser.error("--input is required unless --demo is used")

    check_multiqc()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    if args.demo:
        with tempfile.TemporaryDirectory() as tmp:
            make_demo_fastqc_files(Path(tmp))
            run_multiqc([Path(tmp)], output_dir)
        input_dirs_display: list[Path] = [Path("(demo synthetic FastQC data)")]
    else:
        input_dirs_display = validate_input_dirs(args.input)
        run_multiqc(input_dirs_display, output_dir)

    samples = parse_general_stats(output_dir)
    write_report(output_dir, input_dirs_display, samples)
    write_reproducibility_bundle(output_dir, demo=args.demo, input_dirs=[] if args.demo else input_dirs_display)
    print(f"Report written to: {output_dir / 'report.md'}")


if __name__ == "__main__":
    main()
