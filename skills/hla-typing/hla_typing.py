#!/usr/bin/env python3
"""Hla Typing - HLA allele typing from WGS/WES VCF data."""

import argparse
import json
import shlex
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)
DISCLAIMER = ("ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.")


def parse_args():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", type=Path, dest="input_file", help="Input file path")
    parser.add_argument("--output", type=Path, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
    return parser.parse_args()


def validate_input(input_path: Path) -> dict:
    """Validate and parse the input file. Returns parsed data dict."""
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)
    lines = []
    with open(input_path) as f:
        for line in f:
            line = line.strip()
            if line and not line.startswith("#"):
                lines.append(line)
    return {"lines": lines, "source": str(input_path)}


def run_analysis(data: dict) -> dict:
    """Core analysis logic. Returns result dict."""
    # TODO: implement core hla-typing logic
    return {
        "skill": "hla-typing",
        "version": "0.1.0",
        "source": data.get("source", "unknown"),
        "variants_processed": len(data.get("lines", [])),
        "findings": [],
        "status": "skeleton"
    }


def write_report(result: dict, output_dir: Path) -> None:
    """Write report.md and result.json to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # result.json
    with open(output_dir / "result.json", "w") as f:
        json.dump(result, f, indent=2)

    # report.md
    n = result.get("variants_processed", 0)
    findings = result.get("findings", [])
    report = [
        "# Hla Typing Report",
        "",
        f"**Input**: {result.get('source', 'unknown')}",
        f"**Variants processed**: {n}",
        f"**Findings**: {len(findings)}",
        "",
        "## Results",
        "",
        "| Locus | Finding | Confidence |",
        "|-------|---------|------------|",
    ]
    for f_ in findings:
        report.append(f"| {f_.get('locus', '-')} | {f_.get('finding', '-')} | {f_.get('confidence', '-')} |")
    if not findings:
        report.append("| - | No findings (skeleton implementation) | - |")
    report.extend([
        "",
        "## Summary",
        "",
        f"Analysis completed on {n} variants. {len(findings)} findings reported.",
        "",
        f"*{DISCLAIMER}*",
        "",
    ])
    with open(output_dir / "report.md", "w") as f:
        f.write("\n".join(report))

    print(f"Report written to {output_dir / 'report.md'}")
    print(f"Results written to {output_dir / 'result.json'}")


def run_demo(output_dir: Path) -> None:
    """Run with built-in synthetic demo data."""
    demo_input = SKILL_DIR / "demo_input.txt"
    if not demo_input.exists():
        print("Error: demo data not found", file=sys.stderr)
        sys.exit(1)
    data = validate_input(demo_input)
    result = run_analysis(data)
    write_report(result, output_dir)


def write_repro_bundle(output_dir: Path, run_args: list[str]) -> None:
    """Write reproducibility/ (commands.sh, environment.yml, checksums.sha256)."""
    write_commands_sh(
        output_dir,
        shlex.join(["python", "skills/hla-typing/hla_typing.py", *run_args]),
    )
    # The script uses only the standard library.
    env = write_environment_yml(output_dir, "clawbio-hla-typing", pip_deps=[])
    write_checksums(
        [output_dir / "report.md", output_dir / "result.json", env],
        output_dir,
        anchor=output_dir,
    )


def main():
    args = parse_args()
    if args.demo:
        output = (args.output or Path("/tmp") / "hla_typing" / "demo").resolve()
        run_demo(output)
        run_args = ["--demo", "--output", str(output)]
    elif args.input_file:
        data = validate_input(args.input_file)
        result = run_analysis(data)
        output = (args.output or args.input_file.parent / "output").resolve()
        write_report(result, output)
        run_args = ["--input", str(args.input_file.resolve()), "--output", str(output)]
    else:
        print("Error: provide --input <file> or --demo", file=sys.stderr)
        sys.exit(1)
    write_repro_bundle(output, run_args)


if __name__ == "__main__":
    main()
