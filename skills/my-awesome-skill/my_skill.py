#!/usr/bin/env python3
"""my-awesome-skill: one-line description."""

import argparse
import json
from pathlib import Path
from datetime import datetime

import pandas as pd


def parse_input(input_path: Path) -> pd.DataFrame:
    """Read and validate the input file."""
    df = pd.read_csv(input_path, sep=" ", comment="#")
    # df = pd.read_csv(input_path, sep="\t", comment="#")
    required_cols = {"rsid", "chromosome", "position", "genotype"}
    if not required_cols.issubset(df.columns):
        missing = required_cols - set(df.columns)
        raise ValueError(f"Missing columns: {missing}")
    return df


def analyse(df: pd.DataFrame) -> dict:
    """Core analysis logic. This is where your domain expertise goes."""
    results = {
        "total_variants": len(df),
        "chromosomes": sorted(df["chromosome"].unique().tolist()),
        "timestamp": datetime.now().isoformat(),
    }
    # Add your analysis here
    return results


def generate_report(results: dict, output_dir: Path) -> Path:
    """Write the markdown report."""
    output_dir.mkdir(parents=True, exist_ok=True)
    report_path = output_dir / "report.md"

    lines = [
        "# My Awesome Skill Report",
        "",
        f"**Generated**: {results['timestamp']}",
        "",
        "## Summary",
        "",
        f"- Total variants analysed: {results['total_variants']}",
        f"- Chromosomes covered: {', '.join(str(c) for c in results['chromosomes'])}",
        "",
        "---",
        "",
        '*ClawBio is a research and educational tool. It is not a medical device '
        'and does not provide clinical diagnoses. Consult a healthcare professional '
        'before making any medical decisions.*',
    ]

    report_path.write_text("\n".join(lines))
    return report_path


import matplotlib
matplotlib.use("Agg")  # Non-interactive backend
import matplotlib.pyplot as plt

def plot_summary(results: dict, output_dir: Path):
    """Generate a summary figure."""
    fig, ax = plt.subplots(figsize=(8, 5))
    # Your plot logic here
    ax.set_title("My Awesome Skill: Summary")
    fig_path = output_dir / "figures"
    fig_path.mkdir(exist_ok=True)
    fig.savefig(fig_path / "summary.png", dpi=150, bbox_inches="tight")
    plt.close(fig)
    print(f"Figure saved: {fig_path / 'summary.png'}")


def main():
    parser = argparse.ArgumentParser(description="My Awesome Skill")
    parser.add_argument("--input", type=Path, required=False, help="Input file")
    parser.add_argument("--output", type=Path, default=Path("/tmp/my_skill_output"), help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with built-in demo data")
    args = parser.parse_args()

    if args.demo:
        # Use the demo data bundled with the skill
        input_path = Path(__file__).parent / "demo_input.txt"
    elif args.input:
        input_path = args.input
    else:
        parser.error("Provide --input <file> or use --demo")

    print(f"Reading: {input_path}")
    df = parse_input(input_path)

    print(f"Analysing {len(df)} variants...")
    results = analyse(df)

    report = generate_report(results, args.output)
    print(f"Report written to: {report}")

    plot_summary(results, args.output)


if __name__ == "__main__":
    main()