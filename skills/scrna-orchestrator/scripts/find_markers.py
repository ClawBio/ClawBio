#!/usr/bin/env python3
"""Step 5/7 of scrna_orchestrator: attach Leiden labels and run marker discovery.

Faithful port of `attach_leiden_labels()` + `run_markers()` + the
`cluster_summary`/`markers_top` table-writing pieces of `write_tables()` from
the 1.x skill. Always runs (mandatory step): cluster-vs-rest Wilcoxon marker
detection is core to every scRNA run, not an optional add-on like doublet
detection or annotation.
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import import_scanpy  # noqa: E402


def _cluster_sort_key(cluster: str):
    try:
        return (0, int(cluster))
    except ValueError:
        return (1, cluster)


def attach_leiden_labels(graph_adata, expression_adata):
    """Attach graph-derived Leiden labels onto another AnnData with matching cells. Verbatim from 1.x skill."""
    attached = expression_adata.copy()
    leiden_labels = graph_adata.obs["leiden"].astype(str).reindex(attached.obs_names)
    if bool(leiden_labels.isna().any()):
        raise ValueError("Internal error: missing Leiden labels for normalized cells.")
    category_order = sorted(pd.Index(leiden_labels.dropna().unique()), key=_cluster_sort_key)
    attached.obs["leiden"] = pd.Categorical(leiden_labels.tolist(), categories=category_order, ordered=False)
    return attached


def run_markers(adata, top_markers: int = 10):
    """Run rank_genes_groups and return full + top marker tables. Verbatim from 1.x skill."""
    sc = import_scanpy()

    adata = adata.copy()
    sc.tl.rank_genes_groups(adata, groupby="leiden", method="wilcoxon", pts=True)

    clusters = list(adata.obs["leiden"].cat.categories)
    dfs = []
    for cluster in clusters:
        df = sc.get.rank_genes_groups_df(adata, group=cluster)
        df.insert(0, "cluster", cluster)
        dfs.append(df)

    markers_all = pd.concat(dfs, axis=0, ignore_index=True)
    markers_top = (
        markers_all.sort_values(["cluster", "scores"], ascending=[True, False])
        .groupby("cluster", as_index=False, group_keys=False)
        .head(top_markers)
        .reset_index(drop=True)
    )
    return adata, markers_all, markers_top


def main() -> None:
    parser = argparse.ArgumentParser(description="scrna_orchestrator: find_markers step")
    parser.add_argument("--clustered-h5ad", required=True, help="clustered.h5ad from the cluster step")
    parser.add_argument("--norm-h5ad", required=True, help="normalized.h5ad from the preprocess step")
    parser.add_argument("--top-markers", type=int, default=10)
    parser.add_argument("--out", required=True, help="output path for markers.h5ad")
    parser.add_argument("--top-csv-out", required=True, help="output path for markers_top.csv")
    parser.add_argument("--all-csv-out", required=True, help="output path for markers_all.csv")
    parser.add_argument("--cluster-summary-out", required=True, help="output path for cluster_summary.csv")
    args = parser.parse_args()

    sc = import_scanpy()
    graph_adata = sc.read_h5ad(Path(args.clustered_h5ad))
    norm_adata = sc.read_h5ad(Path(args.norm_h5ad))

    expression_adata = attach_leiden_labels(graph_adata, norm_adata)
    adata_markers, markers_all, markers_top = run_markers(expression_adata, top_markers=args.top_markers)

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    adata_markers.write_h5ad(out_path)

    top_csv_path = Path(args.top_csv_out)
    top_csv_path.parent.mkdir(parents=True, exist_ok=True)
    markers_top.to_csv(top_csv_path, index=False)

    all_csv_path = Path(args.all_csv_out)
    all_csv_path.parent.mkdir(parents=True, exist_ok=True)
    markers_all.to_csv(all_csv_path, index=False)

    cluster_counts = graph_adata.obs["leiden"].value_counts().sort_index()
    cluster_summary = pd.DataFrame(
        {
            "cluster": cluster_counts.index.astype(str),
            "n_cells": cluster_counts.values.astype(int),
            "proportion": (cluster_counts.values / max(1, int(graph_adata.n_obs))).round(4),
        }
    )
    cluster_summary_path = Path(args.cluster_summary_out)
    cluster_summary_path.parent.mkdir(parents=True, exist_ok=True)
    cluster_summary.to_csv(cluster_summary_path, index=False)

    print(f"find_markers complete: {len(markers_top)} top marker rows across {cluster_counts.shape[0]} clusters")
    print(f"  wrote {out_path}")
    print(f"  wrote {top_csv_path}")
    print(f"  wrote {all_csv_path}")
    print(f"  wrote {cluster_summary_path}")


if __name__ == "__main__":
    main()
