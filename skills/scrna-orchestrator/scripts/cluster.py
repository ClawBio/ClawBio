#!/usr/bin/env python3
"""Step 4/7 of scrna_orchestrator: embed (PCA or a latent rep) and Leiden-cluster.

Faithful port of `run_embedding_cluster()` from the 1.x skill. Graph basis is
parameter-driven (`--use-rep auto|none|<obsm key>`), not a DAG branch: this is
ordinary `{var}` substitution into one step's script, which is exactly the
kind of conditional logic CLAUDE.md says belongs in Python, not YAML. Default
`auto` looks for `X_scvi` in `.obsm` (the source skill's optional
latent-downstream mode, e.g. from an `integrated.h5ad` produced by scVI
upstream) and falls back to the standard scale+PCA branch when absent — which
is what the bundled demo fixture exercises.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import import_scanpy  # noqa: E402

DEFAULT_LATENT_REP = "X_scvi"


def _resolve_use_rep(adata, use_rep: str) -> str:
    normalized = (use_rep or "auto").strip()
    if normalized == "none":
        return ""
    if normalized == "auto":
        return DEFAULT_LATENT_REP if DEFAULT_LATENT_REP in adata.obsm else ""
    if normalized not in adata.obsm:
        raise ValueError(f"Requested latent representation not found in adata.obsm: {normalized}")
    return normalized


def run_embedding_cluster(adata, n_pcs, n_neighbors, leiden_resolution, random_state, *, use_rep=""):
    """Compute graph neighbors, UMAP, and Leiden clusters. Verbatim from the 1.x skill."""
    sc = import_scanpy()

    adata = adata.copy()
    if use_rep:
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, use_rep=use_rep)
        n_pcs_eff = 0
    else:
        sc.pp.scale(adata, max_value=10)
        n_pcs_eff = min(n_pcs, adata.n_obs - 1, adata.n_vars - 1)
        if n_pcs_eff < 1:
            raise ValueError(
                "PCA requires at least 2 cells and 2 genes after QC/HVG selection. "
                f"Got n_obs={adata.n_obs}, n_vars={adata.n_vars}."
            )
        sc.tl.pca(adata, n_comps=n_pcs_eff, random_state=random_state, svd_solver="arpack")
        sc.pp.neighbors(adata, n_neighbors=n_neighbors, n_pcs=n_pcs_eff)
    sc.tl.umap(adata, random_state=random_state)
    sc.tl.leiden(adata, resolution=leiden_resolution, random_state=random_state)

    return adata, n_pcs_eff


def main() -> None:
    parser = argparse.ArgumentParser(description="scrna_orchestrator: cluster step")
    parser.add_argument("--hvg-h5ad", required=True, help="hvg.h5ad from the preprocess step")
    parser.add_argument("--norm-h5ad", required=True, help="normalized.h5ad from the preprocess step")
    parser.add_argument("--use-rep", default="auto", help="`auto`, `none`, or an `.obsm` key such as `X_scvi`")
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--n-neighbors", type=int, default=15)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--random-state", type=int, default=0)
    parser.add_argument("--out", required=True, help="output path for clustered.h5ad")
    parser.add_argument("--summary-out", required=True, help="output path for cluster_summary.json")
    args = parser.parse_args()

    sc = import_scanpy()
    norm_adata = sc.read_h5ad(Path(args.norm_h5ad))

    try:
        resolved_use_rep = _resolve_use_rep(norm_adata, args.use_rep)
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    graph_input = norm_adata if resolved_use_rep else sc.read_h5ad(Path(args.hvg_h5ad))

    try:
        adata_graph, n_pcs_eff = run_embedding_cluster(
            graph_input,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            leiden_resolution=args.leiden_resolution,
            random_state=args.random_state,
            use_rep=resolved_use_rep,
        )
    except ValueError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata_graph.write_h5ad(out_path)

    n_clusters = int(adata_graph.obs["leiden"].nunique())
    summary = {
        "graph_basis": resolved_use_rep or "pca",
        "n_pcs_eff": n_pcs_eff,
        "n_neighbors": args.n_neighbors,
        "leiden_resolution": args.leiden_resolution,
        "n_clusters": n_clusters,
        "n_cells": int(adata_graph.n_obs),
    }
    summary_path = Path(args.summary_out)
    summary_path.parent.mkdir(parents=True, exist_ok=True)
    summary_path.write_text(json.dumps(summary, indent=2))

    print(f"cluster complete: graph_basis={summary['graph_basis']}, n_clusters={n_clusters}")
    print(f"  wrote {out_path}")
    print(f"  wrote {summary_path}")


if __name__ == "__main__":
    main()
