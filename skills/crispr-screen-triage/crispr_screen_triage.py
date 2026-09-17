#!/usr/bin/env python3

from __future__ import annotations

import argparse
import csv
import json
import math
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
REQUIRED_COLUMNS = {"guide_id", "gene", "control_count", "treatment_count", "essentiality", "druggability"}


def load_counts(path: Path) -> list[dict[str, str | float]]:
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        missing = REQUIRED_COLUMNS - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"Input is missing required columns: {sorted(missing)}")
        rows: list[dict[str, str | float]] = []
        for raw in reader:
            row: dict[str, str | float] = {"guide_id": raw["guide_id"], "gene": raw["gene"]}
            for column in sorted(REQUIRED_COLUMNS - {"guide_id", "gene"}):
                try:
                    row[column] = float(raw[column])
                except (TypeError, ValueError) as exc:
                    raise ValueError(f"Column {column} must be numeric") from exc
            rows.append(row)
    if not rows:
        raise ValueError("Input contains no guides")
    return rows


def _median(values: list[float]) -> float:
    ordered = sorted(values)
    middle = len(ordered) // 2
    if len(ordered) % 2:
        return ordered[middle]
    return (ordered[middle - 1] + ordered[middle]) / 2.0


def triage_genes(rows: list[dict[str, str | float]]) -> dict:
    guides = []
    grouped: dict[str, list[dict[str, str | float | bool]]] = {}

    for row in rows:
        control = float(row["control_count"])
        treated = float(row["treatment_count"])
        log2fc = math.log2((treated + 1.0) / (control + 1.0))
        guide = {**row, "log2_fold_change": round(log2fc, 3), "low_count": control + treated < 100}
        guides.append(guide)
        grouped.setdefault(str(row["gene"]), []).append(guide)

    genes = []
    for gene, gene_guides in grouped.items():
        median_log2fc = _median([float(guide["log2_fold_change"]) for guide in gene_guides])
        essentiality = sum(float(guide["essentiality"]) for guide in gene_guides) / len(gene_guides)
        druggability = sum(float(guide["druggability"]) for guide in gene_guides) / len(gene_guides)
        depletion = max(0.0, -median_log2fc)
        score = (0.55 * depletion) + (0.25 * druggability) + (0.20 * essentiality)
        priority = "high" if score >= 1.35 and median_log2fc <= -1.0 else "medium" if score >= 0.75 else "watch"
        genes.append(
            {
                "gene": gene,
                "guide_count": len(gene_guides),
                "median_log2_fold_change": round(median_log2fc, 3),
                "essentiality": round(essentiality, 3),
                "druggability": round(druggability, 3),
                "low_count_guides": sum(1 for guide in gene_guides if guide["low_count"]),
                "priority_score": round(score, 3),
                "priority": priority,
            }
        )

    genes.sort(key=lambda item: (-float(item["priority_score"]), str(item["gene"])))
    top_hits = [gene for gene in genes if gene["priority"] in {"high", "medium"}]
    return {
        "skill": "crispr-screen-triage",
        "summary": {"gene_count": len(genes), "guide_count": len(guides), "top_hit_count": len(top_hits)},
        "guides": guides,
        "genes": genes,
        "top_hits": top_hits,
        "disclaimer": DISCLAIMER,
    }


def _write_csv(path: Path, rows: list[dict], fieldnames: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
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
        script_path=Path("skills/crispr-screen-triage/crispr_screen_triage.py"),
        args=args,
        comment="Reproduce this crispr-screen-triage run",
    )


def write_repro_bundle(input_path: Path, output_dir: Path, demo: bool) -> None:
    write_environment_yml(
        output_dir,
        env_name="clawbio-crispr-screen-triage",
        pip_deps=[],
        python_version="3.10",
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
            output_dir / "tables" / "triaged_genes.csv",
            output_dir / "tables" / "guide_metrics.csv",
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
    _write_csv(
        output_dir / "tables" / "triaged_genes.csv",
        result["genes"],
        [
            "gene",
            "guide_count",
            "median_log2_fold_change",
            "essentiality",
            "druggability",
            "low_count_guides",
            "priority_score",
            "priority",
        ],
    )
    _write_csv(
        output_dir / "tables" / "guide_metrics.csv",
        result["guides"],
        [
            "guide_id",
            "gene",
            "control_count",
            "treatment_count",
            "log2_fold_change",
            "essentiality",
            "druggability",
            "low_count",
        ],
    )
    (output_dir / "result.json").write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    rows = [
        "# CRISPR Screen Triage Report",
        "",
        f"**Input**: `{input_path}`",
        f"**Mode**: {'Synthetic demo data' if demo else 'User-provided local data'}",
        f"**Guides scored**: {result['summary']['guide_count']}",
        f"**Genes scored**: {result['summary']['gene_count']}",
        "",
        "| Rank | Gene | Guides | Median log2FC | Essentiality | Druggability | Score | Priority |",
        "|---:|---|---:|---:|---:|---:|---:|---|",
    ]
    for index, gene in enumerate(result["genes"], 1):
        rows.append(
            f"| {index} | {gene['gene']} | {gene['guide_count']} | {gene['median_log2_fold_change']} | "
            f"{gene['essentiality']} | {gene['druggability']} | {gene['priority_score']} | {gene['priority']} |"
        )
    rows.extend(
        [
            "",
            "## Interpretation",
            "",
            "Top hits show replicated guide-level depletion plus user-supplied essentiality and druggability annotations.",
            "",
            DISCLAIMER,
            "",
        ]
    )
    (output_dir / "report.md").write_text("\n".join(rows), encoding="utf-8")
    write_repro_bundle(input_path, output_dir, demo)


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description="CRISPR Screen Triage")
    parser.add_argument("--input", type=Path)
    parser.add_argument("--output", type=Path, default=Path("crispr_screen_triage_out"))
    parser.add_argument("--demo", action="store_true")
    args = parser.parse_args(argv)
    input_path = SKILL_DIR / "demo_screen_counts.csv" if args.demo else args.input
    if input_path is None:
        parser.error("--input is required unless --demo is used")
    try:
        result = triage_genes(load_counts(input_path))
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 2
    write_outputs(result, input_path, args.output, args.demo)
    print(f"CRISPR Screen Triage wrote {args.output / 'report.md'}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
