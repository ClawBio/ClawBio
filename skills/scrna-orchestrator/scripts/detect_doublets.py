#!/usr/bin/env python3
"""Step 2/7 of scrna_orchestrator: OPTIONAL Scrublet doublet detection.

Faithful port of `run_doublet_detection()` from the 1.x skill (method is
always "scrublet" here — the original CLI's `--doublet-method none` choice
*is* this step's `clawbio skip <id> --reason ...`). Runs on the QC-filtered
raw counts and removes predicted doublets before normalization/clustering.

This step is declared `optional: true` in pipeline.yaml: sorted single-cell
protocols or datasets with negligible doublet rates can legitimately skip it.
`preprocess.py` (the next step) discovers whether this step ran by checking
for `doublet_filtered.h5ad` at the conventional path next to its declared
`qc_filtered_h5ad` input, rather than pipeline.yaml gating on this step's
emit — see preprocess.py's module docstring for why.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import import_scanpy  # noqa: E402


def run_doublet_detection(adata, random_state: int):
    """Run Scrublet on QC-filtered raw counts. Verbatim logic from the 1.x skill."""
    import scrublet  # noqa: F401  (import-check, matches source's _import_scrublet)

    sc = import_scanpy()

    adata = adata.copy()
    obs_before = adata.obs.copy()
    n_prin_comps = min(30, adata.n_obs - 1, adata.n_vars - 1)
    if n_prin_comps < 2:
        raise ValueError(
            "Scrublet requires at least 3 cells and 3 genes after QC filtering. "
            f"Got n_obs={adata.n_obs}, n_vars={adata.n_vars}."
        )

    sc.pp.scrublet(adata, log_transform=False, n_prin_comps=n_prin_comps, random_state=random_state, verbose=False)

    for column in obs_before.columns:
        if column not in adata.obs.columns:
            adata.obs[column] = obs_before[column].reindex(adata.obs_names)

    predicted = adata.obs["predicted_doublet"].fillna(False).astype(bool)
    n_cells_scored = int(adata.n_obs)
    n_predicted_doublets = int(predicted.sum())
    filtered = adata[~predicted.to_numpy(), :].copy()
    if filtered.n_obs == 0:
        raise ValueError("Doublet detection removed all cells. Re-run with `clawbio skip detect_doublets`.")

    summary = {
        "method": "scrublet",
        "n_cells_scored": n_cells_scored,
        "n_predicted_doublets": n_predicted_doublets,
        "n_cells_retained": int(filtered.n_obs),
        "predicted_doublet_rate": round(n_predicted_doublets / max(1, n_cells_scored), 4),
    }
    return filtered, summary


def main() -> None:
    parser = argparse.ArgumentParser(description="scrna_orchestrator: detect_doublets step (optional)")
    parser.add_argument("--qc-h5ad", required=True, help="qc_filtered.h5ad from the qc_filter step")
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--out", required=True, help="output path for doublet_filtered.h5ad")
    parser.add_argument("--summary-out", required=True, help="output path for doublet_summary.json")
    args = parser.parse_args()

    sc = import_scanpy()
    adata = sc.read_h5ad(Path(args.qc_h5ad))

    try:
        filtered, summary = run_doublet_detection(adata, random_state=args.random_state)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    filtered.write_h5ad(out_path)

    summary_path = Path(args.summary_out)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2))

    print(
        f"detect_doublets complete: {summary['n_predicted_doublets']} predicted doublets "
        f"({summary['predicted_doublet_rate']:.1%}) of {summary['n_cells_scored']} cells; "
        f"{summary['n_cells_retained']} retained"
    )
    print(f"  wrote {out_path}")
    print(f"  wrote {summary_path}")


if __name__ == "__main__":
    main()
