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
import os
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
DEMO_SOURCE_ENV = "CLAWBIO_SCRNA_DEMO_SOURCE"


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


def load_demo_adata(random_state: int, demo_source_policy: str | None = None):
    """Load real PBMC3k demo data, falling back to synthetic data when needed."""
    sc = _import_scanpy()
    policy = (demo_source_policy or os.getenv(DEMO_SOURCE_ENV, "auto")).strip().lower()
    if policy not in {"auto", "pbmc3k", "synthetic"}:
        policy = "auto"

    if policy == "synthetic":
        return build_demo_adata(random_state), "synthetic_forced"

    try:
        adata = sc.datasets.pbmc3k()
        adata.var_names_make_unique()
        sc.pp.filter_cells(adata, min_counts=1)
        if adata.n_obs == 0:
            raise ValueError("PBMC3k demo had no cells after filtering min_counts=1.")
        return adata, "pbmc3k_raw"
    except Exception as exc:
        print(
            f"WARNING: Failed to load PBMC3k demo ({exc}); falling back to synthetic demo data.",
            file=sys.stderr,
        )
        return build_demo_adata(random_state), "synthetic_fallback"


def _sample_expression_values(x, max_values: int = 200_000) -> np.ndarray:
    """Sample expression values from dense/sparse matrices without densifying sparse input."""
    try:
        from scipy import sparse  # type: ignore
    except Exception:
        sparse = None

    if sparse is not None and sparse.issparse(x):
        values = np.asarray(x.data).ravel()
    else:
        values = np.asarray(x).ravel()

    if values.size > max_values:
        step = max(1, values.size // max_values)
        values = values[::step][:max_values]

    return values.astype(np.float64, copy=False)


def detect_processed_input_reason(adata) -> str | None:
    """Detect whether input looks preprocessed (log-normalized/scaled) instead of raw counts."""
    values = _sample_expression_values(adata.X)
    if values.size == 0:
        return None

    finite = values[np.isfinite(values)]
    if finite.size == 0:
        return None

    uns_markers = {
        "neighbors",
        "pca",
        "umap",
        "rank_genes_groups",
        "draw_graph",
        "louvain",
    }
    uns_hits = sorted(key for key in uns_markers if key in adata.uns)
    has_negative = bool(np.any(finite < -1e-8))
    frac_non_integer = float(np.mean(np.abs(finite - np.rint(finite)) > 1e-6))
    max_val = float(np.max(finite))

    reason: str | None = None
    if has_negative:
        reason = "Detected negative expression values, indicating scaled/transformed input."
    elif frac_non_integer > 0.20 and (max_val <= 50.0 or bool(uns_hits)):
        reason = (
            "Detected mostly non-integer expression values that look like normalized/log-transformed input."
        )

    if reason is None:
        return None

    if uns_hits:
        reason += f" Found processed-analysis metadata in adata.uns: {', '.join(uns_hits)}."
    reason += (
        " This skill expects raw-count .h5ad input. `pbmc3k_processed` is not supported; "
        "use raw counts (e.g., `scanpy.datasets.pbmc3k()`)."
    )
    return reason


def load_data(input_path: str | None, demo: bool, random_state: int):
    """Load AnnData from .h5ad or build demo data."""
    sc = _import_scanpy()

    if demo:
        adata, demo_source = load_demo_adata(random_state)
        return adata, None, True, demo_source

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
    processed_reason = detect_processed_input_reason(adata)
    if processed_reason:
        raise ValueError(processed_reason)
    return adata, path, False, None


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


def resolve_de_request(args: argparse.Namespace) -> dict[str, str] | None:
    """Validate DE arguments and return the request when enabled."""
    provided = {
        "--de-groupby": args.de_groupby,
        "--de-group1": args.de_group1,
        "--de-group2": args.de_group2,
    }
    provided_count = sum(1 for value in provided.values() if value)
    if provided_count == 0:
        return None
    if provided_count != 3:
        missing = [flag for flag, value in provided.items() if not value]
        raise ValueError(
            "DE requires --de-groupby, --de-group1, and --de-group2 together. "
            f"Missing: {', '.join(missing)}."
        )

    return {
        "groupby": str(args.de_groupby),
        "group1": str(args.de_group1),
        "group2": str(args.de_group2),
    }


def run_two_group_de(
    adata,
    *,
    groupby: str,
    group1: str,
    group2: str,
    top_genes: int,
) -> tuple[pd.DataFrame, pd.DataFrame, dict[str, Any]]:
    """Run two-group DE with Wilcoxon and return full/top tables plus summary."""
    sc = _import_scanpy()

    if group1 == group2:
        raise ValueError("--de-group1 and --de-group2 must be different values.")

    if groupby not in adata.obs.columns:
        available_cols = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
        raise ValueError(
            f"DE groupby column not found in adata.obs: {groupby}. "
            f"Available columns: {available_cols}."
        )

    groups = adata.obs[groupby].astype(str)
    available_groups = sorted(groups.dropna().unique().tolist())
    missing_groups = [g for g in (group1, group2) if g not in available_groups]
    if missing_groups:
        raise ValueError(
            f"DE group value(s) not found in {groupby}: {', '.join(missing_groups)}. "
            f"Available groups: {', '.join(available_groups)}."
        )

    mask = groups.isin([group1, group2]).to_numpy()
    if int(mask.sum()) < 2:
        raise ValueError(
            f"DE requires at least 2 cells across {group1} and {group2} in {groupby}."
        )

    adata_de = adata[mask].copy()
    de_groups = groups[mask].tolist()
    adata_de.obs["_de_group"] = pd.Categorical(
        de_groups,
        categories=[group1, group2],
        ordered=True,
    )

    n_group1 = int(sum(1 for g in de_groups if g == group1))
    n_group2 = int(sum(1 for g in de_groups if g == group2))
    if n_group1 == 0 or n_group2 == 0:
        raise ValueError(
            f"DE comparison requires both groups to have cells. Got {group1}={n_group1}, {group2}={n_group2}."
        )

    sc.tl.rank_genes_groups(
        adata_de,
        groupby="_de_group",
        groups=[group1],
        reference=group2,
        method="wilcoxon",
        pts=True,
    )

    de_full = sc.get.rank_genes_groups_df(adata_de, group=group1).reset_index(drop=True)
    if de_full.empty:
        raise ValueError(
            f"DE did not return any genes for {group1} vs {group2} in {groupby}."
        )

    de_top = (
        de_full.sort_values("scores", ascending=False)
        .head(top_genes)
        .reset_index(drop=True)
    )

    summary = {
        "enabled": True,
        "groupby": groupby,
        "group1": group1,
        "group2": group2,
        "n_cells_group1": n_group1,
        "n_cells_group2": n_group2,
        "n_genes_full": int(len(de_full)),
        "top_table": "de_top.csv",
        "full_table": "de_full.csv",
        "top_gene_names": de_top["names"].dropna().astype(str).tolist(),
        "volcano_plot": "",
    }
    return de_full, de_top, summary


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


def write_de_tables(
    de_full: pd.DataFrame,
    de_top: pd.DataFrame,
    tables_dir: Path,
) -> tuple[Path, Path]:
    """Write DE full/top tables."""
    tables_dir.mkdir(parents=True, exist_ok=True)
    de_full_path = tables_dir / "de_full.csv"
    de_top_path = tables_dir / "de_top.csv"
    de_full.to_csv(de_full_path, index=False)
    de_top.to_csv(de_top_path, index=False)
    return de_full_path, de_top_path


def plot_de_volcano(
    de_full: pd.DataFrame,
    figures_dir: Path,
    *,
    group1: str,
    group2: str,
) -> Path:
    """Create DE volcano plot from full two-group DE results."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figures_dir.mkdir(parents=True, exist_ok=True)

    if "logfoldchanges" not in de_full.columns:
        raise ValueError("DE table missing required column: logfoldchanges")

    p_col = "pvals_adj" if "pvals_adj" in de_full.columns else "pvals"
    if p_col not in de_full.columns:
        raise ValueError("DE table missing required p-value column (pvals_adj or pvals).")

    log_fc = pd.to_numeric(de_full["logfoldchanges"], errors="coerce").to_numpy(dtype=np.float64)
    pvals = pd.to_numeric(de_full[p_col], errors="coerce").to_numpy(dtype=np.float64)
    pvals = np.clip(pvals, 1e-300, 1.0)
    neg_log10 = -np.log10(pvals)
    finite_mask = np.isfinite(log_fc) & np.isfinite(neg_log10)
    if int(finite_mask.sum()) == 0:
        raise ValueError("No finite DE points available for volcano plot.")

    sig_mask = finite_mask & (pvals < 0.05) & (np.abs(log_fc) >= 1.0)
    up_mask = sig_mask & (log_fc > 0)
    down_mask = sig_mask & (log_fc < 0)
    nonsig_mask = finite_mask & (~sig_mask)

    fig, ax = plt.subplots(figsize=(7.2, 5.6))
    ax.scatter(
        log_fc[nonsig_mask],
        neg_log10[nonsig_mask],
        s=10,
        c="#94a3b8",
        alpha=0.65,
        edgecolors="none",
        label="Not significant",
    )
    if int(up_mask.sum()) > 0:
        ax.scatter(
            log_fc[up_mask],
            neg_log10[up_mask],
            s=16,
            c="#dc2626",
            alpha=0.85,
            edgecolors="none",
            label=f"Up in {group1}",
        )
    if int(down_mask.sum()) > 0:
        ax.scatter(
            log_fc[down_mask],
            neg_log10[down_mask],
            s=16,
            c="#2563eb",
            alpha=0.85,
            edgecolors="none",
            label=f"Up in {group2}",
        )

    ax.axvline(-1.0, color="#64748b", linewidth=0.9, linestyle="--")
    ax.axvline(1.0, color="#64748b", linewidth=0.9, linestyle="--")
    ax.axhline(-np.log10(0.05), color="#64748b", linewidth=0.9, linestyle="--")
    ax.set_xlabel("log2 fold change")
    y_label = "-log10(adjusted p-value)" if p_col == "pvals_adj" else "-log10(p-value)"
    ax.set_ylabel(y_label)
    ax.set_title(f"DE Volcano: {group1} vs {group2}")
    ax.legend(loc="upper right", frameon=False)
    ax.grid(alpha=0.18, linewidth=0.5)

    plot_path = figures_dir / "de_volcano.png"
    fig.tight_layout()
    fig.savefig(plot_path, dpi=180, bbox_inches="tight")
    plt.close(fig)
    return plot_path


def render_report(
    output_dir: Path,
    input_path: Path | None,
    is_demo: bool,
    demo_source: str | None,
    qc_stats: dict[str, int],
    n_hvg: int,
    n_clusters: int,
    n_pcs_eff: int,
    params: dict[str, Any],
    de_summary: dict[str, Any],
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
            "Demo source": demo_source if is_demo and demo_source else "n/a",
        },
    )

    top_genes = de_summary.get("top_gene_names", [])
    if de_summary.get("enabled"):
        volcano_plot_name = str(de_summary.get("volcano_plot", "")).strip()
        if volcano_plot_name:
            volcano_line = f"- Volcano plot: `figures/{volcano_plot_name}`"
            volcano_image = f"\n![DE Volcano](figures/{volcano_plot_name})\n"
        else:
            volcano_line = "- Volcano plot: not generated (use `--de-volcano`)"
            volcano_image = ""

        de_top_preview = (
            "\n".join([f"- `{gene}`" for gene in top_genes[:10]])
            if top_genes
            else "- None"
        )
        de_section = f"""## Differential Expression (Two-Group)

- Grouping column: `{de_summary["groupby"]}`
- Comparison: `{de_summary["group1"]}` vs `{de_summary["group2"]}`
- Cells in groups: `{de_summary["group1"]}={de_summary["n_cells_group1"]}`, `{de_summary["group2"]}={de_summary["n_cells_group2"]}`
- Genes in full DE table: **{de_summary["n_genes_full"]}**
- Full DE table: `tables/{de_summary["full_table"]}`
- Top DE table: `tables/{de_summary["top_table"]}`
{volcano_line}

Top DE genes by score:
{de_top_preview}
{volcano_image}
"""
        de_methods = (
            "- Differential expression: `scanpy.tl.rank_genes_groups` "
            f"(Wilcoxon, `{de_summary['group1']}` vs `{de_summary['group2']}`, "
            f"`groupby={de_summary['groupby']}`)"
        )
        if volcano_plot_name:
            de_methods += "; volcano plot with thresholds `p<0.05`, `|log2FC|>=1`"
    else:
        de_section = """## Differential Expression (Two-Group)

- Not enabled for this run (use `--de-groupby --de-group1 --de-group2`).
"""
        de_methods = "- Differential expression: not enabled"

    de_tables = ""
    if de_summary.get("enabled"):
        de_tables = (
            "- `tables/de_full.csv`\n"
            "- `tables/de_top.csv`\n"
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
{de_tables}

{de_section}

## Methods

- QC/filtering: `min_genes={params["min_genes"]}`, `min_cells={params["min_cells"]}`, `max_mt_pct={params["max_mt_pct"]}`
- Normalisation: total-count normalisation (`target_sum=1e4`) + `log1p`
- Feature selection: `n_top_hvg={params["n_top_hvg"]}`
- Embedding: `n_pcs={n_pcs_eff}`, `n_neighbors={params["n_neighbors"]}`, UMAP
- Clustering: Leiden `resolution={params["leiden_resolution"]}`
- Marker analysis: `scanpy.tl.rank_genes_groups` (Wilcoxon, cluster-vs-rest)
{de_methods}

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

    if args.de_groupby and args.de_group1 and args.de_group2:
        de_flags = [
            "--de-groupby",
            str(args.de_groupby),
            "--de-group1",
            str(args.de_group1),
            "--de-group2",
            str(args.de_group2),
            "--de-top-genes",
            str(args.de_top_genes),
        ]
        if args.de_volcano:
            de_flags.append("--de-volcano")
        cmd_line = f"{cmd_line} {' '.join(shlex.quote(flag) for flag in de_flags)}"

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
            output_dir / "tables" / "de_full.csv",
            output_dir / "tables" / "de_top.csv",
            output_dir / "figures" / "qc_violin.png",
            output_dir / "figures" / "umap_leiden.png",
            output_dir / "figures" / "marker_dotplot.png",
            output_dir / "figures" / "de_volcano.png",
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

    de_request = resolve_de_request(args)
    if args.de_top_genes < 1:
        raise ValueError("--de-top-genes must be >= 1.")
    if args.de_volcano and de_request is None:
        raise ValueError("--de-volcano requires --de-groupby, --de-group1, and --de-group2.")

    de_summary: dict[str, Any] = {
        "enabled": False,
        "groupby": "",
        "group1": "",
        "group2": "",
        "n_cells_group1": 0,
        "n_cells_group2": 0,
        "n_genes_full": 0,
        "top_table": "",
        "full_table": "",
        "top_gene_names": [],
        "volcano_plot": "",
    }

    adata, input_path, is_demo, demo_source = load_data(args.input, args.demo, args.random_state)
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

    de_table_paths: list[Path] = []
    if de_request:
        # Run DE on log1p-normalized (unscaled) expression values while reusing
        # clustering labels from the embedding step.
        adata_de_input = adata_pp.copy()
        leiden_labels = adata_emb.obs["leiden"].reindex(adata_de_input.obs_names)
        if bool(leiden_labels.isna().any()):
            raise ValueError("Internal error: missing Leiden labels for DE input cells.")
        adata_de_input.obs["leiden"] = leiden_labels.astype(str)
        de_full, de_top, de_summary = run_two_group_de(
            adata_de_input,
            groupby=de_request["groupby"],
            group1=de_request["group1"],
            group2=de_request["group2"],
            top_genes=args.de_top_genes,
        )
        de_full_path, de_top_path = write_de_tables(de_full, de_top, tables_dir)
        de_table_paths = [de_full_path, de_top_path]
        if args.de_volcano:
            volcano_path = plot_de_volcano(
                de_full,
                figures_dir,
                group1=de_summary["group1"],
                group2=de_summary["group2"],
            )
            de_summary["volcano_plot"] = volcano_path.name

    cluster_path, markers_csv, markers_tsv = write_tables(adata_markers, markers_top, tables_dir)
    created_figures = plot_core_figures(adata_markers, markers_top, figures_dir)
    if de_summary.get("volcano_plot"):
        created_figures.append(str(de_summary["volcano_plot"]))

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
        demo_source=demo_source,
        qc_stats=qc_stats,
        n_hvg=n_hvg,
        n_clusters=n_clusters,
        n_pcs_eff=n_pcs_eff,
        params=params,
        de_summary=de_summary,
    )

    tables_written = [cluster_path.name, markers_csv.name, markers_tsv.name]
    if de_table_paths:
        tables_written.extend(path.name for path in de_table_paths)

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
            "tables": tables_written,
            "figures": created_figures,
            "demo_source": demo_source if is_demo else "not_demo",
            "de": {
                "enabled": bool(de_summary["enabled"]),
                "groupby": de_summary["groupby"] if de_summary["enabled"] else "",
                "group1": de_summary["group1"] if de_summary["enabled"] else "",
                "group2": de_summary["group2"] if de_summary["enabled"] else "",
                "n_genes_full": int(de_summary["n_genes_full"]) if de_summary["enabled"] else 0,
                "top_table": de_summary["top_table"] if de_summary["enabled"] else "",
                "volcano_plot": de_summary["volcano_plot"] if de_summary["enabled"] else "",
            },
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
        description="ClawBio scRNA Orchestrator — Scanpy QC, clustering, markers, and optional two-group DE",
    )
    parser.add_argument("--input", "-i", help="Input AnnData file (.h5ad)")
    parser.add_argument("--output", "-o", default="scrna_report", help="Output directory")
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run demo data (PBMC3k raw preferred, fallback to synthetic)",
    )
    parser.add_argument("--min-genes", type=int, default=200, help="Minimum genes per cell")
    parser.add_argument("--min-cells", type=int, default=3, help="Minimum cells per gene")
    parser.add_argument("--max-mt-pct", type=float, default=20.0, help="Maximum mitochondrial percentage")
    parser.add_argument("--n-top-hvg", type=int, default=2000, help="Number of highly variable genes")
    parser.add_argument("--n-pcs", type=int, default=50, help="Number of principal components")
    parser.add_argument("--n-neighbors", type=int, default=15, help="Number of neighbors for graph construction")
    parser.add_argument("--leiden-resolution", type=float, default=1.0, help="Leiden resolution")
    parser.add_argument("--random-state", type=int, default=0, help="Random seed")
    parser.add_argument("--top-markers", type=int, default=10, help="Top markers per cluster")
    parser.add_argument("--de-groupby", default=None, help="obs column for two-group DE")
    parser.add_argument("--de-group1", default=None, help="Group 1 value for DE")
    parser.add_argument("--de-group2", default=None, help="Group 2 reference value for DE")
    parser.add_argument("--de-top-genes", type=int, default=50, help="Top DE genes to include in summary table")
    parser.add_argument("--de-volcano", action="store_true", help="Generate optional DE volcano plot")
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
