#!/usr/bin/env python3
"""Drug Repurposing Screen — pooled viability screen to prioritised candidates.

Format-agnostic, objective-driven analysis of compound × sample panels.
Supports --demo (bundled toy dataset) or --bundle with schema/objective YAML.

Usage:
    python drug_repurposing_screen.py --demo --output /tmp/drs_demo
    python drug_repurposing_screen.py --bundle ./demo --schema ./demo/schema.yaml \\
        --objective ./demo/objective.yaml --output ./out
"""

from __future__ import annotations

import argparse
import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent
REPO_ROOT = SCRIPT_DIR.parents[1]
sys.path.insert(0, str(SCRIPT_DIR))
if str(REPO_ROOT) not in sys.path:
    sys.path.insert(0, str(REPO_ROOT))

from clawbio.common.reproducibility import (
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)
from screen_engine import DISCLAIMER, PipelineResult, load_yaml, run_pipeline


def write_report(result: PipelineResult, output_dir: Path, objective: dict) -> None:
    lines = [
        "# Drug Repurposing Screen Report",
        "",
        f"**Objective:** {objective.get('name', 'unspecified')}",
        f"**Generated:** {datetime.now(timezone.utc).strftime('%Y-%m-%d %H:%M UTC')}",
        "",
        "## Summary",
        "",
        f"- Samples screened: {result.summary['n_samples']}",
        f"- Compounds tested: {result.summary['n_compounds']}",
        f"- Primary hits: {result.summary['n_hits']}",
        f"- Context-selective compounds: {result.summary['n_context_selective']}",
        f"- Top candidate: `{result.summary.get('top_compound', 'n/a')}`",
        "",
        "## Top prioritised candidates",
        "",
    ]
    top = result.priority.head(10)
    if len(top):
        show_cols = [
            c for c in ("rank", "compound_id", "compound_name", "selectivity_class",
                        "priority", "feature", "feature_type", "clinical_phase")
            if c in top.columns
        ]
        sub = top[show_cols]
        lines.append("| " + " | ".join(show_cols) + " |")
        lines.append("| " + " | ".join(["---"] * len(show_cols)) + " |")
        for _, row in sub.iterrows():
            lines.append("| " + " | ".join(str(row[c]) for c in show_cols) + " |")
    else:
        lines.append("_No candidates ranked._")
    lines.extend(["", "## Disclaimer", "", DISCLAIMER])
    (output_dir / "report.md").write_text("\n".join(lines))


def write_reproducibility(output_dir: Path, args: argparse.Namespace) -> None:
    repro = output_dir / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    if args.schema and Path(args.schema).exists():
        shutil.copy(args.schema, repro / "schema.yaml")
    if args.objective and Path(args.objective).exists():
        shutil.copy(args.objective, repro / "objective.yaml")

    cmd_args: list[str | ReproPath] = ["--output", ReproPath(output_dir, anchor="output_dir")]
    if args.demo:
        cmd_args.append("--demo")
    else:
        for flag in ("--bundle", "--schema", "--objective"):
            cmd_args += [flag, ReproPath(Path(getattr(args, flag[2:])).resolve(), anchor="repo_root")]
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=SCRIPT_DIR.relative_to(REPO_ROOT) / "drug_repurposing_screen.py",
            args=cmd_args,
            comment="Reproduce this drug-repurposing-screen run",
        ),
        repo_root=REPO_ROOT,
    )
    write_environment_yml(
        output_dir,
        env_name="clawbio-drug-repurposing-screen",
        pip_deps=["numpy>=1.24", "pandas>=2.0", "scipy>=1.10", "pyyaml>=6.0", "pyarrow>=14.0"],
    )
    # Written last so it covers everything above. report.md/report.html carry a
    # wall-clock timestamp, so they stay out to keep replays byte-identical.
    files = [output_dir / "result.json"]
    files += [output_dir / "tables" / n for n in
              ("selectivity.csv", "biomarker_univariate_all_matrices.csv", "priority_table.csv")]
    files += [output_dir / "cache" / f"{n}.parquet" for n in
              ("qc_primary", "primary_hits", "selectivity", "biomarkers", "priority")]
    files += [repro / n for n in ("schema.yaml", "objective.yaml", "commands.sh", "environment.yml")]
    write_checksums(files, output_dir, anchor=output_dir)


def write_outputs(result: PipelineResult, output_dir: Path, objective: dict, args: argparse.Namespace) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    tables = output_dir / "tables"
    cache = output_dir / "cache"
    figures = output_dir / "figures"
    tables.mkdir(exist_ok=True)
    cache.mkdir(exist_ok=True)
    figures.mkdir(exist_ok=True)

    result.qc_primary.to_parquet(cache / "qc_primary.parquet", index=False)
    result.primary_hits.to_parquet(cache / "primary_hits.parquet", index=False)
    result.selectivity.to_parquet(cache / "selectivity.parquet", index=False)
    result.biomarkers.to_parquet(cache / "biomarkers.parquet", index=False)
    result.priority.to_parquet(cache / "priority.parquet", index=False)

    result.selectivity.to_csv(tables / "selectivity.csv", index=False)
    result.biomarkers.to_csv(tables / "biomarker_univariate_all_matrices.csv", index=False)
    result.priority.to_csv(tables / "priority_table.csv", index=False)

    write_report(result, output_dir, objective)
    (output_dir / "report.html").write_text(
        "<html><body><pre>" + (output_dir / "report.md").read_text() + "</pre></body></html>"
    )
    (output_dir / "result.json").write_text(json.dumps({
        "summary": result.summary,
        "priority_table": result.priority.head(20).to_dict(orient="records"),
        "objective": objective.get("name"),
    }, indent=2))
    write_reproducibility(output_dir, args)


def resolve_demo_paths() -> tuple[Path, Path, Path]:
    demo = SCRIPT_DIR / "demo"
    return demo, demo / "schema.yaml", demo / "objective.yaml"


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="Drug repurposing screen analysis")
    parser.add_argument("--demo", action="store_true", help="Run bundled toy dataset")
    parser.add_argument("--bundle", type=Path, help="Path to screen bundle directory")
    parser.add_argument("--schema", type=Path, help="schema.yaml path")
    parser.add_argument("--objective", type=Path, help="objective.yaml path (required)")
    parser.add_argument("--output", type=Path, required=True, help="Output directory")
    parser.add_argument("--resume", action="store_true", help="Use cached parquet if present")
    args = parser.parse_args(argv)

    if args.demo:
        bundle, schema, objective = resolve_demo_paths()
    else:
        if not args.bundle or not args.schema or not args.objective:
            parser.error("--bundle, --schema, and --objective are required unless --demo is set")
        bundle, schema, objective = args.bundle, args.schema, args.objective

    if not schema.exists():
        print(f"Schema not found: {schema}", file=sys.stderr)
        return 1
    if not objective.exists():
        print(f"Objective not found: {objective}", file=sys.stderr)
        return 1

    objective_data = load_yaml(objective)
    cache_dir = args.output / "cache"
    if args.resume and (cache_dir / "priority.parquet").exists():
        print("Resume: cache hit — skipping rebuild")
        return 0

    print(f"Running drug-repurposing-screen on {bundle} …")
    result = run_pipeline(bundle, schema, objective, cache_dir=cache_dir)
    write_outputs(result, args.output, objective_data, args)
    print(f"Done. {result.summary['n_hits']} hits; top: {result.summary.get('top_compound')}")
    print(f"Outputs: {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
