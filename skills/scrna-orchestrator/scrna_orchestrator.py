#!/usr/bin/env python3
"""ClawBio scRNA Orchestrator (MVP).

Scanpy-based single-cell RNA-seq pipeline:
QC/filtering -> normalisation/log1p -> HVG -> PCA/neighbors/UMAP ->
Leiden clustering -> marker detection.

Usage:
    python scrna_orchestrator.py --input sample.h5ad --output report_dir
    python scrna_orchestrator.py --demo --output demo_report
"""

from __future__ import annotations

import argparse
import json
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file
from clawbio.common.report import (
    generate_report_footer,
    generate_report_header,
    write_result_json,
)

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions."
)


def _import_scanpy():
    """Import scanpy lazily with a clear user-facing error."""
    try:
        import scanpy as sc  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "scanpy is required for scrna-orchestrator. "
            "Install it with: pip install scanpy anndata"
        ) from exc
    return sc


def build_demo_adata(random_state: int):
    """Create deterministic synthetic AnnData demo data."""
    from anndata import AnnData  # type: ignore

    rng = np.random.default_rng(random_state)
    n_cells = 240
    n_genes = 480
    n_clusters = 3
    cells_per_cluster = n_cells // n_clusters

    base_profiles = []
    for i in range(n_clusters):
        base = rng.gamma(shape=2.0, scale=1.2, size=n_genes)
        marker_start = 40 + i * 20
        marker_end = marker_start + 15
        base[marker_start:marker_end] += 6.0
        base_profiles.append(base)

    expr_blocks = []
    labels = []
    for i, base in enumerate(base_profiles):
        lam = np.clip(base, 0.05, None)
        counts = rng.poisson(lam=lam, size=(cells_per_cluster, n_genes))
        libsize_scale = rng.lognormal(mean=0.0, sigma=0.35, size=(cells_per_cluster, 1))
        counts = np.round(counts * libsize_scale).astype(np.int32)
        expr_blocks.append(counts)
        labels.extend([f"cluster_{i}"] * cells_per_cluster)

    x = np.vstack(expr_blocks)
    gene_names = [f"Gene{i:03d}" for i in range(n_genes)]
    for i in range(20):
        gene_names[i] = f"MT-GENE{i:02d}"

    obs = pd.DataFrame(
        {
            "sample_id": [f"cell_{i:03d}" for i in range(n_cells)],
            "demo_truth": labels,
        },
        index=[f"cell_{i:03d}" for i in range(n_cells)],
    )
    var = pd.DataFrame(index=gene_names)
    return AnnData(X=x, obs=obs, var=var)


def load_data(input_path: str | None, demo: bool, random_state: int):
    """Load AnnData from .h5ad or build demo data."""
    sc = _import_scanpy()

    if demo:
        adata = build_demo_adata(random_state)
        return adata, None, True

    if not input_path:
        raise ValueError("Provide --input <file.h5ad> or --demo.")

    path = Path(input_path)
    if not path.exists():
        raise FileNotFoundError(f"Input file not found: {path}")
    if path.suffix.lower() != ".h5ad":
        raise ValueError(
            f"Only .h5ad is supported in MVP. Received: {path.name}"
        )

    adata = sc.read_h5ad(path)
    return adata, path, False


def qc_filter(
    adata,
    min_genes: int,
    min_cells: int,
    max_mt_pct: float,
) -> tuple[Any, dict[str, int]]:
    """Compute QC metrics and apply filtering."""
    sc = _import_scanpy()

    adata = adata.copy()
    adata.var_names_make_unique()
    adata.var["mt"] = adata.var_names.str.upper().str.startswith("MT-")
    sc.pp.calculate_qc_metrics(
        adata,
        qc_vars=["mt"],
        percent_top=None,
        log1p=False,
        inplace=True,
    )

    stats = {
        "n_cells_before": int(adata.n_obs),
        "n_genes_before": int(adata.n_vars),
    }

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


def run_preprocess(adata, n_top_hvg: int):
    """Normalise, log-transform, and select HVGs."""
    sc = _import_scanpy()

    adata = adata.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    sc.pp.highly_variable_genes(adata, n_top_genes=n_top_hvg, flavor="seurat")

    n_hvg = int(adata.var["highly_variable"].sum())
    if n_hvg == 0:
        raise ValueError("No highly variable genes found.")

    adata = adata[:, adata.var["highly_variable"]].copy()
    return adata, n_hvg


def run_embedding_cluster(
    adata,
    n_pcs: int,
    n_neighbors: int,
    leiden_resolution: float,
    random_state: int,
):
    """Compute PCA, graph neighbors, UMAP, and Leiden clusters."""
    sc = _import_scanpy()

    adata = adata.copy()
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


def run_markers(adata, top_markers: int = 10):
    """Run rank_genes_groups and return full + top marker tables."""
    sc = _import_scanpy()

    adata = adata.copy()
    sc.tl.rank_genes_groups(
        adata,
        groupby="leiden",
        method="wilcoxon",
        pts=True,
    )

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


def plot_core_figures(adata, markers_top: pd.DataFrame, figures_dir: Path) -> list[str]:
    """Create QC/UMAP/marker plots."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc = _import_scanpy()
    figures_dir.mkdir(parents=True, exist_ok=True)

    created: list[str] = []

    sc.pl.violin(
        adata,
        ["n_genes_by_counts", "total_counts", "pct_counts_mt"],
        jitter=0.4,
        multi_panel=True,
        show=False,
    )
    plt.tight_layout()
    qc_path = figures_dir / "qc_violin.png"
    plt.savefig(qc_path, dpi=180, bbox_inches="tight")
    plt.close("all")
    created.append(qc_path.name)

    sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False)
    plt.tight_layout()
    umap_path = figures_dir / "umap_leiden.png"
    plt.savefig(umap_path, dpi=180, bbox_inches="tight")
    plt.close("all")
    created.append(umap_path.name)

    marker_genes = (
        markers_top.groupby("cluster", as_index=False)
        .head(3)["names"]
        .dropna()
        .astype(str)
        .tolist()
    )
    marker_genes = list(dict.fromkeys(marker_genes))[:20]
    if marker_genes:
        dot = sc.pl.dotplot(
            adata,
            var_names=marker_genes,
            groupby="leiden",
            show=False,
            return_fig=True,
        )
        marker_path = figures_dir / "marker_dotplot.png"
        dot.savefig(marker_path, dpi=180)
        plt.close("all")
        created.append(marker_path.name)

    return created


def write_tables(
    adata,
    markers_top: pd.DataFrame,
    tables_dir: Path,
) -> tuple[Path, Path, Path]:
    """Write cluster summary and marker tables."""
    tables_dir.mkdir(parents=True, exist_ok=True)

    cluster_counts = adata.obs["leiden"].value_counts().sort_index()
    cluster_summary = pd.DataFrame(
        {
            "cluster": cluster_counts.index.astype(str),
            "n_cells": cluster_counts.values.astype(int),
            "proportion": (cluster_counts.values / max(1, int(adata.n_obs))).round(4),
        }
    )
    cluster_path = tables_dir / "cluster_summary.csv"
    cluster_summary.to_csv(cluster_path, index=False)

    csv_path = tables_dir / "markers_top.csv"
    tsv_path = tables_dir / "markers_top.tsv"
    markers_top.to_csv(csv_path, index=False)
    markers_top.to_csv(tsv_path, sep="\t", index=False)

    return cluster_path, csv_path, tsv_path


def render_report(
    output_dir: Path,
    input_path: Path | None,
    is_demo: bool,
    qc_stats: dict[str, int],
    n_hvg: int,
    n_clusters: int,
    n_pcs_eff: int,
    params: dict[str, Any],
) -> Path:
    """Create markdown report.md."""
    input_files = [input_path] if input_path else []
    header = generate_report_header(
        title="scRNA Orchestrator Report",
        skill_name="scrna-orchestrator",
        input_files=input_files,
        extra_metadata={
            "Mode": "demo" if is_demo else "input",
            "Cells (before QC)": str(qc_stats["n_cells_before"]),
            "Cells (after QC)": str(qc_stats["n_cells_after"]),
            "Genes (after QC)": str(qc_stats["n_genes_after"]),
            "Leiden clusters": str(n_clusters),
            "HVG selected": str(n_hvg),
        },
    )

    body = f"""## Summary

- Cells before QC: **{qc_stats["n_cells_before"]}**
- Cells after QC: **{qc_stats["n_cells_after"]}**
- Genes before QC: **{qc_stats["n_genes_before"]}**
- Genes after QC: **{qc_stats["n_genes_after"]}**
- HVGs selected: **{n_hvg}**
- Leiden clusters: **{n_clusters}**

## Core Figures

![QC Violin](figures/qc_violin.png)
![UMAP Leiden](figures/umap_leiden.png)
![Marker Dotplot](figures/marker_dotplot.png)

## Tables

- `tables/cluster_summary.csv`
- `tables/markers_top.csv`
- `tables/markers_top.tsv`

## Methods

- QC/filtering: `min_genes={params["min_genes"]}`, `min_cells={params["min_cells"]}`, `max_mt_pct={params["max_mt_pct"]}`
- Normalisation: total-count normalisation (`target_sum=1e4`) + `log1p`
- Feature selection: `n_top_hvg={params["n_top_hvg"]}`
- Embedding: `n_pcs={n_pcs_eff}`, `n_neighbors={params["n_neighbors"]}`, UMAP
- Clustering: Leiden `resolution={params["leiden_resolution"]}`
- Marker analysis: `scanpy.tl.rank_genes_groups` (Wilcoxon, cluster-vs-rest)

## Reproducibility

See:
- `reproducibility/commands.sh`
- `reproducibility/environment.yml`
- `reproducibility/checksums.sha256`
"""

    report_path = output_dir / "report.md"
    report_path.write_text(header + body + generate_report_footer(), encoding="utf-8")
    return report_path


def write_reproducibility(
    output_dir: Path,
    input_path: Path | None,
    is_demo: bool,
    args: argparse.Namespace,
) -> None:
    """Write commands.sh, environment.yml, and checksums.sha256."""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    quoted_output = shlex.quote(str(output_dir))

    if is_demo:
        cmd_line = (
            "python skills/scrna-orchestrator/scrna_orchestrator.py "
            f"--demo --output {quoted_output}"
        )
    else:
        if input_path is None:
            raise ValueError("input_path is required when --demo is not used.")
        quoted_input = shlex.quote(str(input_path))
        cmd_line = (
            "python skills/scrna-orchestrator/scrna_orchestrator.py "
            f"--input {quoted_input} --output {quoted_output}"
        )

    commands = f"""#!/usr/bin/env bash
# Reproducibility script — ClawBio scRNA Orchestrator
# Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

set -euo pipefail

{cmd_line}
"""
    (repro_dir / "commands.sh").write_text(commands, encoding="utf-8")

    env_yml = """name: clawbio-scrna
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - scanpy=1.10.2
  - anndata=0.10.8
  - numpy=1.26.4
  - pandas=2.2.2
  - matplotlib=3.8.4
  - seaborn=0.13.2
  - leidenalg=0.10.2
  - python-igraph=0.11.6
"""
    (repro_dir / "environment.yml").write_text(env_yml, encoding="utf-8")

    checksum_targets: list[Path] = []
    if input_path and input_path.exists():
        checksum_targets.append(input_path)
    checksum_targets.extend(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            output_dir / "tables" / "cluster_summary.csv",
            output_dir / "tables" / "markers_top.csv",
            output_dir / "tables" / "markers_top.tsv",
            output_dir / "figures" / "qc_violin.png",
            output_dir / "figures" / "umap_leiden.png",
            output_dir / "figures" / "marker_dotplot.png",
        ]
    )

    lines: list[str] = []
    for path in checksum_targets:
        if not path.exists():
            continue
        rel = path.relative_to(output_dir) if path.is_relative_to(output_dir) else path.name
        lines.append(f"{sha256_file(path)}  {rel}")
    (repro_dir / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(args: argparse.Namespace) -> dict[str, Any]:
    """Run the full scRNA MVP pipeline."""
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    figures_dir.mkdir(exist_ok=True)
    tables_dir.mkdir(exist_ok=True)

    adata, input_path, is_demo = load_data(args.input, args.demo, args.random_state)
    adata_qc, qc_stats = qc_filter(
        adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mt_pct=args.max_mt_pct,
    )
    adata_pp, n_hvg = run_preprocess(adata_qc, n_top_hvg=args.n_top_hvg)
    adata_emb, n_pcs_eff = run_embedding_cluster(
        adata_pp,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        leiden_resolution=args.leiden_resolution,
        random_state=args.random_state,
    )
    adata_markers, markers_all, markers_top = run_markers(
        adata_emb,
        top_markers=args.top_markers,
    )

    cluster_path, markers_csv, markers_tsv = write_tables(adata_markers, markers_top, tables_dir)
    created_figures = plot_core_figures(adata_markers, markers_top, figures_dir)

    n_clusters = int(adata_markers.obs["leiden"].nunique())
    params = {
        "min_genes": args.min_genes,
        "min_cells": args.min_cells,
        "max_mt_pct": args.max_mt_pct,
        "n_top_hvg": args.n_top_hvg,
        "n_pcs": args.n_pcs,
        "n_neighbors": args.n_neighbors,
        "leiden_resolution": args.leiden_resolution,
        "random_state": args.random_state,
    }
    report_path = render_report(
        output_dir=output_dir,
        input_path=input_path,
        is_demo=is_demo,
        qc_stats=qc_stats,
        n_hvg=n_hvg,
        n_clusters=n_clusters,
        n_pcs_eff=n_pcs_eff,
        params=params,
    )

    write_result_json(
        output_dir=output_dir,
        skill="scrna",
        version="0.1.0",
        summary={
            "n_cells_before": qc_stats["n_cells_before"],
            "n_cells_after": qc_stats["n_cells_after"],
            "n_genes_before": qc_stats["n_genes_before"],
            "n_genes_after": qc_stats["n_genes_after"],
            "n_hvg": n_hvg,
            "n_clusters": n_clusters,
        },
        data={
            "cluster_labels": sorted(adata_markers.obs["leiden"].astype(str).unique().tolist()),
            "tables": [cluster_path.name, markers_csv.name, markers_tsv.name],
            "figures": created_figures,
            "disclaimer": DISCLAIMER,
        },
        input_checksum=sha256_file(input_path) if input_path else "",
    )

    write_reproducibility(output_dir, input_path, is_demo, args)

    return {
        "report_path": report_path,
        "output_dir": output_dir,
        "n_clusters": n_clusters,
        "n_cells_after": qc_stats["n_cells_after"],
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description="ClawBio scRNA Orchestrator — Scanpy QC/clustering/markers MVP",
    )
    parser.add_argument("--input", "-i", help="Input AnnData file (.h5ad)")
    parser.add_argument("--output", "-o", default="scrna_report", help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data")
    parser.add_argument("--min-genes", type=int, default=200, help="Minimum genes per cell")
    parser.add_argument("--min-cells", type=int, default=3, help="Minimum cells per gene")
    parser.add_argument("--max-mt-pct", type=float, default=20.0, help="Maximum mitochondrial percentage")
    parser.add_argument("--n-top-hvg", type=int, default=2000, help="Number of highly variable genes")
    parser.add_argument("--n-pcs", type=int, default=50, help="Number of principal components")
    parser.add_argument("--n-neighbors", type=int, default=15, help="Number of neighbors for graph construction")
    parser.add_argument("--leiden-resolution", type=float, default=1.0, help="Leiden resolution")
    parser.add_argument("--random-state", type=int, default=0, help="Random seed")
    parser.add_argument("--top-markers", type=int, default=10, help="Top markers per cluster")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.demo and not args.input:
        print("ERROR: Provide --input <file.h5ad> or --demo", file=sys.stderr)
        sys.exit(1)

    try:
        result = run_pipeline(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    print("\nscRNA Orchestrator complete")
    print(f"  Report: {result['report_path']}")
    print(f"  Output: {result['output_dir']}")
    print(f"  Cells after QC: {result['n_cells_after']}")
    print(f"  Leiden clusters: {result['n_clusters']}")


if __name__ == "__main__":
    main()
