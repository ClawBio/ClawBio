#!/usr/bin/env python3
"""Step 3/7 of scrna_orchestrator: normalize, log1p, and select HVGs.

Faithful port of `run_preprocess()` from the 1.x skill.

**This is the DAG's answer to the optional/skip design question the task
description flags**: does clustering *strictly need* doublet-filtered data,
or does it merely benefit from it? Doublet removal is a legitimate, commonly
skipped step (sorted single-cell protocols, negligible doublet rates) — the
data is still scientifically valid without it, just includes some predicted
doublets. So `pipeline.yaml` gives this step `prereqs: [qc_filter,
detect_doublets]` (ordering only for detect_doublets — wait until it's done
*or skipped* before proceeding) and `requires: {qc_filtered_h5ad: present}`
(the one truly mandatory data dependency, from the non-optional qc_filter
step). It deliberately does NOT declare `requires: {doublet_filtered_h5ad:
present}` — that would block this step (and everything downstream) forever
the moment `detect_doublets` is skipped, defeating the whole point of the
`optional: true` feature.

Instead, this script does the branch itself (the sanctioned place for
conditional logic per CLAUDE.md: "no Turing-complete YAML ... write a Python
script and make that script the step body"): it always receives
`qc_filtered.h5ad` via `--qc-h5ad`, and opportunistically looks for a sibling
`doublet_filtered.h5ad` in the same output directory — present only if
`detect_doublets` actually ran (not skipped). All steps in this run share one
flat `{out}` directory (see `qc_filter.py`, `detect_doublets.py`
`command_template`s), so this sibling-file lookup is deterministic: by the
time this step is unblocked, `detect_doublets` has definitively resolved to
`done` or `skipped`, so there's no race.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import import_scanpy  # noqa: E402


def run_preprocess(adata, n_top_hvg: int):
    """Normalise, log-transform, and prepare full-gene + HVG branches. Verbatim from 1.x skill."""
    sc = import_scanpy()

    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    sc.pp.highly_variable_genes(adata_norm, n_top_genes=n_top_hvg, flavor="seurat")

    n_hvg = int(adata_norm.var["highly_variable"].sum())
    if n_hvg == 0:
        raise ValueError("No highly variable genes found.")

    adata_hvg = adata_norm[:, adata_norm.var["highly_variable"]].copy()
    return adata_norm, adata_hvg, n_hvg


def main() -> None:
    parser = argparse.ArgumentParser(description="scrna_orchestrator: preprocess step")
    parser.add_argument("--qc-h5ad", required=True, help="qc_filtered.h5ad from the qc_filter step (mandatory)")
    parser.add_argument("--n-top-hvg", type=int, default=2000)
    parser.add_argument("--norm-out", required=True, help="output path for normalized.h5ad (full-gene)")
    parser.add_argument("--hvg-out", required=True, help="output path for hvg.h5ad (HVG subset)")
    parser.add_argument("--summary-out", required=True, help="output path for preprocess_summary.json")
    args = parser.parse_args()

    sc = import_scanpy()
    qc_path = Path(args.qc_h5ad)

    # Opportunistic sibling lookup: use doublet-filtered cells if detect_doublets ran,
    # otherwise fall back to the mandatory qc_filtered.h5ad. See module docstring.
    doublet_path = qc_path.parent / "doublet_filtered.h5ad"
    if doublet_path.is_file():
        adata = sc.read_h5ad(doublet_path)
        input_source = "doublet_filtered"
    else:
        adata = sc.read_h5ad(qc_path)
        input_source = "qc_filtered"

    try:
        adata_norm, adata_hvg, n_hvg = run_preprocess(adata, n_top_hvg=args.n_top_hvg)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    norm_path = Path(args.norm_out)
    norm_path.parent.mkdir(parents=True, exist_ok=True)
    adata_norm.write_h5ad(norm_path)

    hvg_path = Path(args.hvg_out)
    hvg_path.parent.mkdir(parents=True, exist_ok=True)
    adata_hvg.write_h5ad(hvg_path)

    summary = {
        "input_source": input_source,
        "doublet_detection_ran": input_source == "doublet_filtered",
        "n_cells": int(adata_norm.n_obs),
        "n_genes": int(adata_norm.n_vars),
        "n_hvg": n_hvg,
        "n_top_hvg_requested": args.n_top_hvg,
    }
    summary_path = Path(args.summary_out)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2))

    print(f"preprocess complete: input_source={input_source}, n_hvg={n_hvg}")
    print(f"  wrote {norm_path}")
    print(f"  wrote {hvg_path}")
    print(f"  wrote {summary_path}")


if __name__ == "__main__":
    main()
