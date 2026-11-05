#!/usr/bin/env python3
"""Stage 7: fan-in — render report.md and result.json.

Reads every stage's `*_meta.json` sidecar (written by align/trim/
select_model/build_tree/root) plus render.py's table/figure/tree outputs
and assembles them into ClawBio's contract outputs, mirroring 1.x's
`report.md` + `result.json` writer at the end of `main()`. 1.x's own
reproducibility-bundle writer (`write_reproducibility_bundle`/
`write_checksums`) is intentionally NOT reproduced — ClawBio 2.0's
`clawbio bundle` verb now owns that job for the whole run.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

import common

BOOTSTRAP_LABELS = {
    "ufboot": "UFBoot (1 000 replicates)",
    "standard": "Standard bootstrap (100 replicates)",
    "all": "Triple support: UFBoot + aLRT + aBayes",
}


def load_json(path: Path) -> dict:
    return json.loads(path.read_text())


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path, help="Original input FASTA (for filename display)")
    p.add_argument("--align-meta", required=True, type=Path)
    p.add_argument("--trim-meta", required=True, type=Path)
    p.add_argument("--model-meta", required=True, type=Path)
    p.add_argument("--tree-meta", required=True, type=Path)
    p.add_argument("--root-meta", required=True, type=Path)
    p.add_argument("--tree", required=True, type=Path, help="Final Newick tree (from render)")
    p.add_argument("--table", required=True, type=Path, help="Branch-support CSV (from render)")
    p.add_argument("--figure", required=True, type=Path, help="Phylogram PNG (from render)")
    p.add_argument("--out-report", required=True, type=Path, help="Output report.md path")
    p.add_argument("--out-result", required=True, type=Path, help="Output result.json path")
    args = p.parse_args()

    args.out_report.parent.mkdir(parents=True, exist_ok=True)
    args.out_result.parent.mkdir(parents=True, exist_ok=True)
    output_dir = args.out_report.parent

    print("[report] assembling report.md + result.json")

    align_meta = load_json(args.align_meta)
    trim_meta = load_json(args.trim_meta)
    model_meta = load_json(args.model_meta)
    tree_meta = load_json(args.tree_meta)
    root_meta = load_json(args.root_meta)

    num_taxa = align_meta.get("num_taxa", 0)
    aligner_used = align_meta.get("aligner")
    trimmed = bool(trim_meta.get("trimmed", False))
    model = tree_meta.get("model") or model_meta.get("model", "GTR+G")
    engine_used = tree_meta.get("engine_used", "iqtree2")
    bootstrap = tree_meta.get("bootstrap", "ufboot")
    run_mode = tree_meta.get("run_mode", "live")
    root_mode = root_meta.get("root_mode")
    outgroup = root_meta.get("outgroup") or tree_meta.get("outgroup")

    pipeline_steps = [
        align_meta.get("stage", ""),
        trim_meta.get("stage", ""),
        model_meta.get("stage", ""),
        tree_meta.get("stage", ""),
        root_meta.get("stage", ""),
    ]

    newick = args.tree.read_text().strip()
    table = pd.read_csv(args.table).to_dict("records")

    def rel(path: Path) -> str:
        try:
            return str(path.relative_to(output_dir))
        except ValueError:
            return path.name

    # ── result.json (ClawBio contract) ──────────────────────────────────
    chat_summary = [
        f"Phylogenetic tree inferred for {num_taxa} taxa.",
        f"Substitution model: **{model}** (selected via ModelFinder).",
        f"Tree engine: {engine_used} | Bootstrap: {bootstrap} | "
        f"Trimming: {'trimAl -automated1' if trimmed else 'skipped'}.",
        f"Internal node support values in `{rel(args.table)}`.",
        f"Proportional phylogram at `{rel(args.figure)}`.",
    ]
    preferred_artifacts = [
        {"type": "figure", "path": rel(args.figure)},
        {"type": "report", "path": rel(args.out_report)},
        {"type": "table", "path": rel(args.table)},
        {"type": "tree", "path": rel(args.tree)},
    ]
    suggested_actions = [
        f"Open {rel(args.figure)} to inspect the tree topology.",
        "Use --root midpoint to apply midpoint rooting (requires ETE3).",
        "Try --bootstrap all for triple support (UFBoot + aLRT + aBayes).",
        "Annotate the tree in FigTree or ggtree (R) for publication quality.",
    ]
    if engine_used == "precomputed":
        suggested_actions.insert(0, "Install IQ-TREE2 or RAxML-NG for live inference.")

    result_data = {
        "engine": engine_used,
        "aligner": aligner_used,
        "aligned": align_meta.get("aligned", False),
        "trimmed": trimmed,
        "model": model,
        "bootstrap_mode": bootstrap,
        "num_taxa": num_taxa,
        "tree": newick,
        "input_file": str(args.input),
        "output_dir": str(output_dir),
        "pipeline_steps": pipeline_steps,
        "run_mode": run_mode,
        "status": "success",
        "chat_summary_lines": chat_summary,
        "preferred_artifacts": preferred_artifacts,
        "workflow_state": "completed",
        "suggested_actions": suggested_actions,
        "contract_alerts": (
            ["Tree is pre-computed fallback — install IQ-TREE2 for live inference"]
            if engine_used == "precomputed"
            else []
        ),
    }
    args.out_result.write_text(json.dumps(result_data, indent=2))

    # ── report.md ────────────────────────────────────────────────────────
    bs_label = BOOTSTRAP_LABELS.get(bootstrap, bootstrap)
    rooting_label = (
        "midpoint" if root_mode == "midpoint"
        else (f"outgroup: {outgroup}" if outgroup else "unrooted")
    )
    report_lines = [
        "# Phylogenetics Builder Report",
        "",
        "## Pipeline Summary",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        f"| Input | `{args.input.name}` |",
        f"| Taxa | {num_taxa} |",
        f"| Aligner | {aligner_used if aligner_used else 'pre-aligned'} |",
        f"| Trimming | {'trimAl -automated1' if trimmed else 'skipped'} |",
        f"| Substitution model | `{model}` |",
        f"| Tree engine | {engine_used} |",
        f"| Bootstrap | {bs_label} |",
        f"| Rooting | {rooting_label} |",
        "",
        "## Pipeline Steps",
        "",
    ]
    for step in pipeline_steps:
        report_lines.append(f"- `{step}`")
    report_lines += [
        "",
        "## Branch Lengths & Support Values",
        "",
        "| Node / Taxon | Branch Length | Support |",
        "|:-------------|:-------------:|:-------:|",
    ]
    for row in table:
        report_lines.append(f"| {row['node']} | {row['length']:.5f} | {row['support']} |")
    report_lines += [
        "",
        "## Outputs",
        "",
        "| Artifact | Path |",
        "|----------|------|",
        f"| Newick tree | `{rel(args.tree)}` |",
        f"| Phylogram figure | `{rel(args.figure)}` |",
        f"| Branch support table | `{rel(args.table)}` |",
        "| Reproducibility bundle | `.clawbio-run/` (see `clawbio bundle`) |",
        "",
        "---",
        "",
        f"*{common.DISCLAIMER}*",
        "",
    ]
    args.out_report.write_text("\n".join(report_lines))

    print(f"  Saved: {args.out_report}, {args.out_result}")


if __name__ == "__main__":
    main()
