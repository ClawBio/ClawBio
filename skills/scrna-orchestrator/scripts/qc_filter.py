#!/usr/bin/env python3
"""Step 1/7 of scrna_orchestrator: load raw counts and apply QC filtering.

Faithful port of `load_data()` + `qc_filter()` from the ClawBio 1.x
`scrna-orchestrator/scrna_orchestrator.py` (h5ad + 10x Matrix Market input;
mitochondrial-percentage + min-genes/min-cells filtering). This is the DAG's
only mandatory upstream step — every later step (optional or not) ultimately
traces its `requires` gate back to this step's `qc_filtered_h5ad` emit.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import detect_processed_input_reason, import_scanpy, resolve_input_source, load_10x_mtx_data  # noqa: E402


def qc_filter(adata, min_genes: int, min_cells: int, max_mt_pct: float):
    """Compute QC metrics and apply filtering. Verbatim logic from the 1.x skill."""
    sc = import_scanpy()

    adata = adata.copy()
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, log1p=False, inplace=True)

    stats = {"n_cells_before": int(adata.n_obs), "n_genes_before": int(adata.n_vars)}

    adata = adata[adata.obs["n_genes_by_counts"] >= min_genes, :].copy()
    adata = adata[adata.obs["pct_counts_mt"] <= max_mt_pct, :].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)

    stats["n_cells_after"] = int(adata.n_obs)
    stats["n_genes_after"] = int(adata.n_vars)

    if adata.n_obs == 0:
        raise ValueError("Filtering removed all cells. Adjust QC thresholds.")
    if adata.n_vars == 0:
        raise ValueError("Filtering removed all genes. Adjust QC thresholds.")

    return adata, stats


def main() -> None:
    parser = argparse.ArgumentParser(description="scrna_orchestrator: qc_filter step")
    parser.add_argument("--input", required=True, help="raw-count .h5ad or 10x Matrix Market input")
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--min-cells", type=int, default=3)
    parser.add_argument("--max-mt-pct", type=float, default=20.0)
    parser.add_argument("--out", required=True, help="output path for qc_filtered.h5ad")
    parser.add_argument("--stats-out", required=True, help="output path for qc_stats.json")
    args = parser.parse_args()

    sc = import_scanpy()
    source_info = resolve_input_source(Path(args.input))
    if source_info["format"] == "h5ad":
        adata = sc.read_h5ad(Path(source_info["input_path"]))
    else:
        adata = load_10x_mtx_data(source_info)

    reason = detect_processed_input_reason(adata, expected_input="raw-count .h5ad or 10x single-cell input")
    if reason:
        print(f"ERROR: {reason}", file=sys.stderr)
        sys.exit(1)

    try:
        filtered, stats = qc_filter(
            adata, min_genes=args.min_genes, min_cells=args.min_cells, max_mt_pct=args.max_mt_pct
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_h5ad(out_path)

    stats["input_format"] = source_info["format"]
    stats["input_files"] = [str(p) for p in source_info["files"]]
    stats_path = Path(args.stats_out)
    stats_path.parent.mkdir(parents=True, exist_ok=True)
    stats_path.write_text(json.dumps(stats, indent=2))

    print(f"qc_filter complete: {stats['n_cells_after']} cells x {stats['n_genes_after']} genes retained")
    print(f"  wrote {out_path}")
    print(f"  wrote {stats_path}")


if __name__ == "__main__":
    main()
