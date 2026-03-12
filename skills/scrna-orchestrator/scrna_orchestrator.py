#!/usr/bin/env python3
"""ClawBio scRNA Orchestrator.

Scanpy-based single-cell RNA-seq pipeline:
QC/filtering -> optional doublet detection -> normalisation/log1p ->
optional CellTypist annotation -> HVG -> PCA/neighbors/UMAP ->
Leiden clustering -> marker detection.

Usage:
    python scrna_orchestrator.py --input sample.h5ad --output report_dir
    python scrna_orchestrator.py --input filtered_feature_bc_matrix --output report_dir
    python scrna_orchestrator.py --demo --output demo_report
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
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
DEFAULT_CELLTYPIST_MODEL = "Immune_All_Low"


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


def _import_scrublet():
    """Import scrublet for optional doublet detection."""
    try:
        import scrublet  # noqa: F401
    except ImportError as exc:
        raise RuntimeError(
            "scrublet is required for --doublet-method scrublet. "
            "Install it with: pip install scrublet"
        ) from exc


def _import_celltypist():
    """Import celltypist for optional annotation."""
    try:
        import celltypist  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "celltypist is required for --annotate celltypist. "
            "Install it with: pip install celltypist"
        ) from exc
    return celltypist


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
        " This skill expects raw-count single-cell input. `pbmc3k_processed` is not supported; "
        "use raw counts (e.g., `scanpy.datasets.pbmc3k()`)."
    )
    return reason


def _cluster_sort_key(cluster: str) -> tuple[int, Any]:
    """Sort cluster labels numerically when possible."""
    try:
        return (0, int(cluster))
    except ValueError:
        return (1, cluster)


def _normalize_celltypist_model_name(model_name: str) -> str:
    """Normalize a CellTypist model identifier to a .pkl filename when needed."""
    normalized = model_name.strip()
    if not normalized:
        return f"{DEFAULT_CELLTYPIST_MODEL}.pkl"
    if "/" not in normalized and "\\" not in normalized and not normalized.endswith(".pkl"):
        normalized = f"{normalized}.pkl"
    return normalized


def resolve_celltypist_model_path(celltypist, model_name: str) -> Path:
    """Resolve a CellTypist model to a local path without triggering downloads."""
    normalized = _normalize_celltypist_model_name(model_name)
    if "/" in normalized or "\\" in normalized:
        path = Path(normalized).expanduser()
    else:
        path = Path(celltypist.models.models_path) / normalized

    if path.exists():
        return path

    raise RuntimeError(
        "CellTypist model "
        f"'{normalized}' was not found locally at {path}. "
        "Install it first with: "
        f"python -c \"import celltypist; celltypist.models.download_models(model='{normalized}')\". "
        "Runtime downloads are disabled for this skill."
    )


def _is_h5ad_input_path(path: Path) -> bool:
    """Return True when the input points to an AnnData file."""
    return path.is_file() and path.suffix.lower() == ".h5ad"


def _is_matrix_market_input_path(path: Path) -> bool:
    """Return True when the input is a matrix.mtx or matrix.mtx.gz file."""
    name = path.name.lower()
    return path.is_file() and (name.endswith(".mtx") or name.endswith(".mtx.gz"))


def resolve_10x_mtx_source(path: Path) -> dict[str, Any]:
    """Resolve a 10x Matrix Market input from a directory or matrix file."""
    source_dir = path if path.is_dir() else path.parent
    if not source_dir.exists():
        raise FileNotFoundError(f"Input path not found: {path}")

    if path.is_dir():
        matrix_candidates = sorted(
            [
                candidate
                for candidate in source_dir.iterdir()
                if candidate.is_file() and candidate.name.endswith(("matrix.mtx", "matrix.mtx.gz"))
            ]
        )
        if not matrix_candidates:
            raise ValueError(
                "10x Matrix Market input requires a directory containing "
                "`matrix.mtx` or `matrix.mtx.gz`."
            )
        if len(matrix_candidates) > 1:
            candidate_names = ", ".join(candidate.name for candidate in matrix_candidates)
            raise ValueError(
                "Multiple 10x matrix files were found in the input directory. "
                "Point `--input` to a specific `matrix.mtx` or `matrix.mtx.gz` file. "
                f"Found: {candidate_names}"
            )
        matrix_path = matrix_candidates[0]
    else:
        if not _is_matrix_market_input_path(path):
            raise ValueError(
                "10x Matrix Market input must be a `matrix.mtx(.gz)` file or a directory "
                "containing one."
            )
        matrix_path = path

    if matrix_path.name.endswith("matrix.mtx.gz"):
        prefix = matrix_path.name[: -len("matrix.mtx.gz")]
        compressed = True
        features_path = source_dir / f"{prefix}features.tsv.gz"
        barcodes_path = source_dir / f"{prefix}barcodes.tsv.gz"
        missing = [candidate.name for candidate in (features_path, barcodes_path) if not candidate.exists()]
        if missing:
            raise ValueError(
                "Compressed 10x Matrix Market input requires matching `features.tsv.gz` "
                f"and `barcodes.tsv.gz` files. Missing: {', '.join(missing)}"
            )
        input_files = [matrix_path, barcodes_path, features_path]
        feature_table_path = features_path
    else:
        prefix = matrix_path.name[: -len("matrix.mtx")]
        compressed = False
        barcodes_path = source_dir / f"{prefix}barcodes.tsv"
        features_path = source_dir / f"{prefix}features.tsv"
        genes_path = source_dir / f"{prefix}genes.tsv"
        if not barcodes_path.exists():
            raise ValueError(
                "Uncompressed 10x Matrix Market input requires a matching `barcodes.tsv` file."
            )
        if features_path.exists():
            feature_table_path = features_path
        elif genes_path.exists():
            feature_table_path = genes_path
        else:
            raise ValueError(
                "Uncompressed 10x Matrix Market input requires either `features.tsv` "
                "or legacy `genes.tsv` alongside `matrix.mtx`."
            )
        input_files = [matrix_path, barcodes_path, feature_table_path]

    return {
        "format": "10x_mtx",
        "reader_path": source_dir,
        "files": input_files,
        "compressed": compressed,
        "prefix": prefix,
        "matrix_path": matrix_path,
        "barcodes_path": barcodes_path,
        "feature_table_path": feature_table_path,
    }


def resolve_input_source(path: Path) -> dict[str, Any]:
    """Resolve supported input types into a normalized metadata bundle."""
    if not path.exists():
        raise FileNotFoundError(f"Input path not found: {path}")

    if _is_h5ad_input_path(path):
        return {
            "format": "h5ad",
            "reader_path": path,
            "files": [path],
            "compressed": False,
            "prefix": "",
        }

    if path.is_dir() or _is_matrix_market_input_path(path):
        return resolve_10x_mtx_source(path)

    raise ValueError(
        "Supported inputs are raw-count `.h5ad` and 10x Matrix Market inputs "
        "(`matrix.mtx`, `matrix.mtx.gz`, or a directory containing them). "
        f"Received: {path.name}"
    )


def compute_input_checksum(input_source: dict[str, Any] | None) -> str:
    """Compute a stable checksum for one or more input files."""
    if input_source is None:
        return ""

    input_files = sorted((Path(path) for path in input_source["files"]), key=lambda path: path.name)
    if len(input_files) == 1:
        return sha256_file(input_files[0])

    digest = hashlib.sha256()
    for path in input_files:
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()


def load_10x_mtx_data(source_info: dict[str, Any]):
    """Load 10x Matrix Market input without depending on Scanpy's reader signature."""
    from anndata import AnnData  # type: ignore
    from scipy import sparse  # type: ignore
    from scipy.io import mmread  # type: ignore

    matrix_path = Path(source_info["matrix_path"])
    barcodes_path = Path(source_info["barcodes_path"])
    feature_table_path = Path(source_info["feature_table_path"])

    if matrix_path.suffix == ".gz":
        with gzip.open(matrix_path, "rb") as handle:
            matrix = mmread(handle)
    else:
        matrix = mmread(str(matrix_path))

    if sparse.issparse(matrix):
        x = matrix.T.tocsr()
    else:
        x = np.asarray(matrix).T

    barcodes = pd.read_csv(
        barcodes_path,
        header=None,
        sep="\t",
        compression="infer",
    )
    features = pd.read_csv(
        feature_table_path,
        header=None,
        sep="\t",
        compression="infer",
    )

    if features.shape[1] < 2:
        raise ValueError(
            "10x feature table must contain at least gene ID and gene symbol columns."
        )

    obs_names = pd.Index(barcodes.iloc[:, 0].astype(str), dtype="object")
    obs_names.name = None
    var_names = pd.Index(features.iloc[:, 1].astype(str), dtype="object")
    var_names.name = None
    obs = pd.DataFrame(index=obs_names)
    var = pd.DataFrame(index=var_names)
    var["gene_ids"] = features.iloc[:, 0].astype(str).to_numpy()
    if features.shape[1] >= 3:
        var["feature_types"] = features.iloc[:, 2].astype(str).to_numpy()

    adata = AnnData(X=x, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def load_data(input_path: str | None, demo: bool, random_state: int):
    """Load AnnData from supported inputs or build demo data."""
    sc = _import_scanpy()

    if demo:
        adata, demo_source = load_demo_adata(random_state)
        return adata, None, True, demo_source, None

    if not input_path:
        raise ValueError("Provide --input <input.h5ad|matrix.mtx|10x_dir> or --demo.")

    path = Path(input_path)
    source_info = resolve_input_source(path)
    if source_info["format"] == "h5ad":
        adata = sc.read_h5ad(path)
    else:
        adata = load_10x_mtx_data(source_info)
    processed_reason = detect_processed_input_reason(adata)
    if processed_reason:
        raise ValueError(processed_reason)
    return adata, path, False, None, source_info


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


def run_doublet_detection(
    adata,
    method: str,
    random_state: int,
) -> tuple[Any, dict[str, Any] | None]:
    """Run optional doublet detection on QC-filtered raw counts."""
    if method == "none":
        return adata.copy(), None

    _import_scrublet()
    sc = _import_scanpy()

    adata = adata.copy()
    obs_before = adata.obs.copy()
    n_prin_comps = min(30, adata.n_obs - 1, adata.n_vars - 1)
    if n_prin_comps < 2:
        raise ValueError(
            "Scrublet requires at least 3 cells and 3 genes after QC filtering. "
            f"Got n_obs={adata.n_obs}, n_vars={adata.n_vars}."
        )

    sc.pp.scrublet(
        adata,
        log_transform=False,
        n_prin_comps=n_prin_comps,
        random_state=random_state,
        verbose=False,
    )

    for column in obs_before.columns:
        if column not in adata.obs.columns:
            adata.obs[column] = obs_before[column].reindex(adata.obs_names)

    predicted = adata.obs["predicted_doublet"].fillna(False).astype(bool)
    n_cells_scored = int(adata.n_obs)
    n_predicted_doublets = int(predicted.sum())
    filtered = adata[~predicted.to_numpy(), :].copy()
    if filtered.n_obs == 0:
        raise ValueError("Doublet detection removed all cells. Re-run with --doublet-method none.")

    summary = {
        "method": method,
        "n_cells_scored": n_cells_scored,
        "n_predicted_doublets": n_predicted_doublets,
        "n_cells_retained": int(filtered.n_obs),
        "predicted_doublet_rate": round(n_predicted_doublets / max(1, n_cells_scored), 4),
    }
    return filtered, summary


def run_preprocess(adata, n_top_hvg: int):
    """Normalise, log-transform, and prepare full-gene + HVG branches."""
    sc = _import_scanpy()

    adata_norm = adata.copy()
    sc.pp.normalize_total(adata_norm, target_sum=1e4)
    sc.pp.log1p(adata_norm)
    sc.pp.highly_variable_genes(adata_norm, n_top_genes=n_top_hvg, flavor="seurat")

    n_hvg = int(adata_norm.var["highly_variable"].sum())
    if n_hvg == 0:
        raise ValueError("No highly variable genes found.")

    adata_hvg = adata_norm[:, adata_norm.var["highly_variable"]].copy()
    return adata_norm, adata_hvg, n_hvg


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


def aggregate_cluster_annotations(
    clusters: pd.Series,
    predicted_labels: pd.Series,
    conf_scores: pd.Series,
    model_name: str,
) -> pd.DataFrame:
    """Aggregate per-cell CellTypist predictions to cluster-level summaries."""
    frame = pd.DataFrame(
        {
            "cluster": clusters.astype(str),
            "predicted_cell_type": predicted_labels.astype(str),
            "conf_score": conf_scores.astype(float),
        },
        index=clusters.index,
    )

    rows = []
    cluster_order = sorted(frame["cluster"].unique().tolist(), key=_cluster_sort_key)
    for cluster in cluster_order:
        group = frame.loc[frame["cluster"] == cluster]
        counts = group["predicted_cell_type"].value_counts()
        max_count = int(counts.max())
        winning_labels = sorted(counts[counts == max_count].index.astype(str).tolist())
        majority_label = winning_labels[0]
        majority_mask = group["predicted_cell_type"] == majority_label
        rows.append(
            {
                "cluster": cluster,
                "n_cells": int(group.shape[0]),
                "predicted_cell_type": majority_label,
                "support_fraction": round(float(majority_mask.mean()), 4),
                "mean_confidence": round(
                    float(group.loc[majority_mask, "conf_score"].mean()),
                    4,
                ),
                "annotation_model": model_name,
            }
        )

    return pd.DataFrame(rows)


def run_celltypist_annotation(adata, model_name: str) -> tuple[pd.DataFrame, dict[str, Any]]:
    """Run optional local CellTypist annotation and return cluster-level summaries."""
    celltypist = _import_celltypist()
    model_path = resolve_celltypist_model_path(celltypist, model_name)
    model = celltypist.models.Model.load(str(model_path))

    overlap = int(np.isin(np.asarray(adata.var_names, dtype=object), model.features).sum())
    if overlap == 0:
        raise RuntimeError(
            "CellTypist annotation requires human gene symbols overlapping the local model. "
            f"Found 0 overlapping genes with {model_path.name}."
        )

    adata = adata.copy()
    adata.var_names_make_unique()
    result = celltypist.annotate(
        adata,
        model=model,
        majority_voting=False,
    )

    predicted = result.predicted_labels["predicted_labels"].reindex(adata.obs_names)
    conf_scores = result.probability_matrix.max(axis=1).reindex(adata.obs_names)
    annotations = aggregate_cluster_annotations(
        clusters=adata.obs["leiden"],
        predicted_labels=predicted,
        conf_scores=conf_scores,
        model_name=model_path.name,
    )
    metadata = {
        "backend": "celltypist",
        "model": model_path.name,
        "model_path": str(model_path),
        "overlap_genes": overlap,
        "putative": True,
        "n_clusters_annotated": int(annotations.shape[0]),
    }
    return annotations, metadata


def plot_core_figures(adata, markers_top: pd.DataFrame, figures_dir: Path) -> list[Path]:
    """Create QC/UMAP/marker plots."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc = _import_scanpy()
    figures_dir.mkdir(parents=True, exist_ok=True)

    created: list[Path] = []

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
    created.append(qc_path)

    sc.pl.umap(adata, color="leiden", legend_loc="on data", show=False)
    plt.tight_layout()
    umap_path = figures_dir / "umap_leiden.png"
    plt.savefig(umap_path, dpi=180, bbox_inches="tight")
    plt.close("all")
    created.append(umap_path)

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
        created.append(marker_path)

    return created


def write_tables(
    adata,
    markers_top: pd.DataFrame,
    tables_dir: Path,
    doublet_summary: dict[str, Any] | None = None,
    annotation_table: pd.DataFrame | None = None,
) -> dict[str, Path]:
    """Write cluster, marker, and optional feature tables."""
    tables_dir.mkdir(parents=True, exist_ok=True)

    table_paths: dict[str, Path] = {}
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
    table_paths["cluster_summary"] = cluster_path

    csv_path = tables_dir / "markers_top.csv"
    tsv_path = tables_dir / "markers_top.tsv"
    markers_top.to_csv(csv_path, index=False)
    markers_top.to_csv(tsv_path, sep="\t", index=False)
    table_paths["markers_top_csv"] = csv_path
    table_paths["markers_top_tsv"] = tsv_path

    if doublet_summary is not None:
        doublet_path = tables_dir / "doublet_summary.csv"
        pd.DataFrame([doublet_summary]).to_csv(doublet_path, index=False)
        table_paths["doublet_summary"] = doublet_path

    if annotation_table is not None:
        annotation_path = tables_dir / "cluster_annotations.csv"
        annotation_table.to_csv(annotation_path, index=False)
        table_paths["cluster_annotations"] = annotation_path

    return table_paths


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
    input_source: dict[str, Any] | None,
    is_demo: bool,
    demo_source: str | None,
    qc_stats: dict[str, int],
    n_hvg: int,
    n_clusters: int,
    n_pcs_eff: int,
    params: dict[str, Any],
    table_paths: dict[str, Path],
    figure_paths: list[Path],
    n_cells_analyzed: int,
    de_summary: dict[str, Any],
    doublet_summary: dict[str, Any] | None = None,
    annotation_table: pd.DataFrame | None = None,
    annotation_info: dict[str, Any] | None = None,
) -> Path:
    """Create markdown report.md."""
    input_files = [Path(path) for path in input_source["files"]] if input_source else []
    header = generate_report_header(
        title="scRNA Orchestrator Report",
        skill_name="scrna-orchestrator",
        input_files=input_files,
        extra_metadata={
            "Mode": "demo" if is_demo else "input",
            "Input format": "demo" if is_demo else input_source["format"] if input_source else "unknown",
            "Cells (before QC)": str(qc_stats["n_cells_before"]),
            "Cells (after QC)": str(n_cells_analyzed),
            "Genes (after QC)": str(qc_stats["n_genes_after"]),
            "Leiden clusters": str(n_clusters),
            "HVG selected": str(n_hvg),
            "Doublet method": params["doublet_method"],
            "Annotation backend": params["annotate"],
            "Demo source": demo_source if is_demo and demo_source else "n/a",
        },
    )

    figure_labels = {
        "qc_violin.png": "QC Violin",
        "umap_leiden.png": "UMAP Leiden",
        "marker_dotplot.png": "Marker Dotplot",
    }
    lines = ["## Summary", ""]
    lines.append(f"- Cells before QC: **{qc_stats['n_cells_before']}**")
    lines.append(f"- Cells after QC: **{n_cells_analyzed}**")
    if doublet_summary is not None:
        lines.append(f"- Cells after filter-only QC: **{qc_stats['n_cells_after']}**")
        lines.append(f"- Predicted doublets: **{doublet_summary['n_predicted_doublets']}**")
        lines.append(
            "- Doublet detection: "
            f"**{doublet_summary['method']}** ({doublet_summary['predicted_doublet_rate']:.1%} predicted)"
        )
    lines.append(f"- Genes before QC: **{qc_stats['n_genes_before']}**")
    lines.append(f"- Genes after QC: **{qc_stats['n_genes_after']}**")
    lines.append(f"- HVGs selected: **{n_hvg}**")
    lines.append(f"- Leiden clusters: **{n_clusters}**")
    if annotation_info is not None:
        lines.append(f"- Clusters annotated: **{annotation_info['n_clusters_annotated']}**")

    lines.extend(["", "## Core Figures", ""])
    for figure_path in figure_paths:
        if figure_path.name not in figure_labels:
            continue
        lines.append(f"![{figure_labels[figure_path.name]}](figures/{figure_path.name})")

    lines.extend(["", "## Tables", ""])
    for table_path in table_paths.values():
        lines.append(f"- `tables/{table_path.name}`")

    if doublet_summary is not None:
        lines.extend(["", "## Doublet Detection", ""])
        lines.append("- Method: `scanpy.pp.scrublet`")
        lines.append(f"- Cells scored: **{doublet_summary['n_cells_scored']}**")
        lines.append(f"- Predicted doublets: **{doublet_summary['n_predicted_doublets']}**")
        lines.append(f"- Cells retained: **{doublet_summary['n_cells_retained']}**")
        lines.append("- Summary table: `tables/doublet_summary.csv`")

    if annotation_table is not None and annotation_info is not None:
        lines.extend(["", "## Cell Type Annotation", ""])
        lines.append("- Backend: `CellTypist`")
        lines.append(f"- Model: `{annotation_info['model']}`")
        lines.append(f"- Overlapping genes with model: **{annotation_info['overlap_genes']}**")
        lines.append("- Labels are **putative**, model-based assignments and should be manually reviewed.")
        lines.append("- Summary table: `tables/cluster_annotations.csv`")
        lines.append("")
        for row in annotation_table.itertuples(index=False):
            lines.append(
                f"- Cluster `{row.cluster}` -> **{row.predicted_cell_type}** "
                f"(support={row.support_fraction:.2f}, mean_confidence={row.mean_confidence:.2f})"
            )

    lines.extend(["", "## Differential Expression (Two-Group)", ""])
    top_genes = de_summary.get("top_gene_names", [])
    if de_summary.get("enabled"):
        lines.append(f"- Grouping column: `{de_summary['groupby']}`")
        lines.append(f"- Comparison: `{de_summary['group1']}` vs `{de_summary['group2']}`")
        lines.append(
            f"- Cells in groups: `{de_summary['group1']}={de_summary['n_cells_group1']}`, "
            f"`{de_summary['group2']}={de_summary['n_cells_group2']}`"
        )
        lines.append(f"- Genes in full DE table: **{de_summary['n_genes_full']}**")
        lines.append(f"- Full DE table: `tables/{de_summary['full_table']}`")
        lines.append(f"- Top DE table: `tables/{de_summary['top_table']}`")
        volcano_plot_name = str(de_summary.get("volcano_plot", "")).strip()
        if volcano_plot_name:
            lines.append(f"- Volcano plot: `figures/{volcano_plot_name}`")
        else:
            lines.append("- Volcano plot: not generated (use `--de-volcano`)")
        lines.append("")
        lines.append("Top DE genes by score:")
        if top_genes:
            lines.extend([f"- `{gene}`" for gene in top_genes[:10]])
        else:
            lines.append("- None")
        if volcano_plot_name:
            lines.extend(["", f"![DE Volcano](figures/{volcano_plot_name})"])
        de_methods = (
            "- Differential expression: `scanpy.tl.rank_genes_groups` "
            f"(Wilcoxon, `{de_summary['group1']}` vs `{de_summary['group2']}`, "
            f"`groupby={de_summary['groupby']}`)"
        )
        if volcano_plot_name:
            de_methods += "; volcano plot with thresholds `p<0.05`, `|log2FC|>=1`"
    else:
        lines.append("- Not enabled for this run (use `--de-groupby --de-group1 --de-group2`).")
        de_methods = "- Differential expression: not enabled"

    lines.extend(["", "## Methods", ""])
    lines.append(
        "- QC/filtering: "
        f"`min_genes={params['min_genes']}`, `min_cells={params['min_cells']}`, `max_mt_pct={params['max_mt_pct']}`"
    )
    if doublet_summary is not None:
        lines.append(
            "- Doublet detection: "
            "`scanpy.pp.scrublet` on QC-filtered raw counts before normalization/clustering"
        )
    lines.append("- Normalisation: total-count normalisation (`target_sum=1e4`) + `log1p`")
    lines.append(f"- Feature selection: `n_top_hvg={params['n_top_hvg']}`")
    lines.append(
        "- Embedding: "
        f"`n_pcs={n_pcs_eff}`, `n_neighbors={params['n_neighbors']}`, UMAP"
    )
    lines.append(f"- Clustering: Leiden `resolution={params['leiden_resolution']}`")
    lines.append("- Marker analysis: `scanpy.tl.rank_genes_groups` (Wilcoxon, cluster-vs-rest)")
    if annotation_info is not None:
        lines.append(
            "- Annotation: "
            f"`CellTypist` model `{annotation_info['model']}` on normalized/log1p full-gene expression"
        )
    lines.append(de_methods)

    lines.extend(["", "## Reproducibility", "", "See:"])
    lines.append("- `reproducibility/commands.sh`")
    lines.append("- `reproducibility/environment.yml`")
    lines.append("- `reproducibility/checksums.sha256`")

    report_path = output_dir / "report.md"
    report_path.write_text(header + "\n".join(lines) + generate_report_footer(), encoding="utf-8")
    return report_path


def build_repro_command(
    output_dir: Path,
    input_path: Path | None,
    is_demo: bool,
    args: argparse.Namespace,
) -> str:
    """Build a reproducible CLI command for commands.sh."""
    parts = ["python", "skills/scrna-orchestrator/scrna_orchestrator.py"]
    if is_demo:
        parts.append("--demo")
    else:
        if input_path is None:
            raise ValueError("input_path is required when --demo is not used.")
        parts.extend(["--input", str(input_path)])

    parts.extend(["--output", str(output_dir)])

    tunable_defaults = [
        ("--min-genes", args.min_genes, 200),
        ("--min-cells", args.min_cells, 3),
        ("--max-mt-pct", args.max_mt_pct, 20.0),
        ("--n-top-hvg", args.n_top_hvg, 2000),
        ("--n-pcs", args.n_pcs, 50),
        ("--n-neighbors", args.n_neighbors, 15),
        ("--leiden-resolution", args.leiden_resolution, 1.0),
        ("--random-state", args.random_state, 0),
        ("--top-markers", args.top_markers, 10),
    ]
    for flag, value, default in tunable_defaults:
        if value != default:
            parts.extend([flag, str(value)])

    if args.doublet_method != "none":
        parts.extend(["--doublet-method", args.doublet_method])
    if args.annotate != "none":
        parts.extend(["--annotate", args.annotate])
        parts.extend(["--annotation-model", args.annotation_model])
    if args.de_groupby and args.de_group1 and args.de_group2:
        parts.extend(["--de-groupby", str(args.de_groupby)])
        parts.extend(["--de-group1", str(args.de_group1)])
        parts.extend(["--de-group2", str(args.de_group2)])
        parts.extend(["--de-top-genes", str(args.de_top_genes)])
        if args.de_volcano:
            parts.append("--de-volcano")

    return " ".join(shlex.quote(part) for part in parts)


def write_reproducibility(
    output_dir: Path,
    input_path: Path | None,
    input_source: dict[str, Any] | None,
    is_demo: bool,
    args: argparse.Namespace,
    table_paths: dict[str, Path],
    figure_paths: list[Path],
) -> None:
    """Write commands.sh, environment.yml, and checksums.sha256."""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    cmd_line = build_repro_command(output_dir, input_path, is_demo, args)

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
  - pip
  - pip:
      - scrublet==0.2.3
      - celltypist==1.7.1
"""
    (repro_dir / "environment.yml").write_text(env_yml, encoding="utf-8")

    checksum_targets: list[Path] = []
    if input_source is not None:
        checksum_targets.extend(Path(path) for path in input_source["files"] if Path(path).exists())
    checksum_targets.extend(
        [
            output_dir / "report.md",
            output_dir / "result.json",
            *table_paths.values(),
            *figure_paths,
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
    """Run the full scRNA pipeline."""
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

    adata, input_path, is_demo, demo_source, input_source = load_data(
        args.input,
        args.demo,
        args.random_state,
    )
    adata_qc, qc_stats = qc_filter(
        adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mt_pct=args.max_mt_pct,
    )
    adata_clean, doublet_summary = run_doublet_detection(
        adata_qc,
        method=args.doublet_method,
        random_state=args.random_state,
    )
    adata_norm, adata_hvg, n_hvg = run_preprocess(adata_clean, n_top_hvg=args.n_top_hvg)
    adata_emb, n_pcs_eff = run_embedding_cluster(
        adata_hvg,
        n_pcs=args.n_pcs,
        n_neighbors=args.n_neighbors,
        leiden_resolution=args.leiden_resolution,
        random_state=args.random_state,
    )
    adata_markers, _, markers_top = run_markers(
        adata_emb,
        top_markers=args.top_markers,
    )

    leiden_labels = adata_markers.obs["leiden"].astype(str).reindex(adata_norm.obs_names)
    if bool(leiden_labels.isna().any()):
        raise ValueError("Internal error: missing Leiden labels for normalized cells.")
    adata_norm.obs["leiden"] = leiden_labels

    annotation_table = None
    annotation_info = None
    if args.annotate != "none":
        annotation_table, annotation_info = run_celltypist_annotation(
            adata_norm.copy(),
            model_name=args.annotation_model,
        )

    de_full = None
    de_top = None
    if de_request:
        de_full, de_top, de_summary = run_two_group_de(
            adata_norm.copy(),
            groupby=de_request["groupby"],
            group1=de_request["group1"],
            group2=de_request["group2"],
            top_genes=args.de_top_genes,
        )

    table_paths = write_tables(
        adata_markers,
        markers_top,
        tables_dir,
        doublet_summary=doublet_summary,
        annotation_table=annotation_table,
    )
    if de_request and de_full is not None and de_top is not None:
        de_full_path, de_top_path = write_de_tables(de_full, de_top, tables_dir)
        table_paths["de_full"] = de_full_path
        table_paths["de_top"] = de_top_path

    figure_paths = plot_core_figures(adata_markers, markers_top, figures_dir)
    if de_request and de_full is not None and args.de_volcano:
        volcano_path = plot_de_volcano(
            de_full,
            figures_dir,
            group1=de_summary["group1"],
            group2=de_summary["group2"],
        )
        de_summary["volcano_plot"] = volcano_path.name
        figure_paths.append(volcano_path)

    n_clusters = int(adata_markers.obs["leiden"].nunique())
    n_cells_analyzed = int(adata_markers.n_obs)
    params = {
        "min_genes": args.min_genes,
        "min_cells": args.min_cells,
        "max_mt_pct": args.max_mt_pct,
        "n_top_hvg": args.n_top_hvg,
        "n_pcs": args.n_pcs,
        "n_neighbors": args.n_neighbors,
        "leiden_resolution": args.leiden_resolution,
        "random_state": args.random_state,
        "doublet_method": args.doublet_method,
        "annotate": args.annotate,
        "annotation_model": args.annotation_model,
    }
    report_path = render_report(
        output_dir=output_dir,
        input_path=input_path,
        input_source=input_source,
        is_demo=is_demo,
        demo_source=demo_source,
        qc_stats=qc_stats,
        n_hvg=n_hvg,
        n_clusters=n_clusters,
        n_pcs_eff=n_pcs_eff,
        params=params,
        table_paths=table_paths,
        figure_paths=figure_paths,
        n_cells_analyzed=n_cells_analyzed,
        de_summary=de_summary,
        doublet_summary=doublet_summary,
        annotation_table=annotation_table,
        annotation_info=annotation_info,
    )

    summary = {
        "n_cells_before": qc_stats["n_cells_before"],
        "n_cells_after": n_cells_analyzed,
        "n_genes_before": qc_stats["n_genes_before"],
        "n_genes_after": qc_stats["n_genes_after"],
        "n_hvg": n_hvg,
        "n_clusters": n_clusters,
    }
    if doublet_summary is not None:
        summary["n_predicted_doublets"] = doublet_summary["n_predicted_doublets"]
    if annotation_info is not None:
        summary["n_clusters_annotated"] = annotation_info["n_clusters_annotated"]

    data: dict[str, Any] = {
        "cluster_labels": sorted(
            adata_markers.obs["leiden"].astype(str).unique().tolist(),
            key=_cluster_sort_key,
        ),
        "input": {
            "format": "demo" if is_demo else input_source["format"] if input_source else "unknown",
            "files": [] if input_source is None else [Path(path).name for path in input_source["files"]],
        },
        "tables": [path.name for path in table_paths.values()],
        "figures": [path.name for path in figure_paths],
        "demo_source": demo_source if is_demo else "not_demo",
        "de": {
            "enabled": bool(de_summary["enabled"]),
            "groupby": de_summary["groupby"] if de_summary["enabled"] else "",
            "group1": de_summary["group1"] if de_summary["enabled"] else "",
            "group2": de_summary["group2"] if de_summary["enabled"] else "",
            "n_genes_full": int(de_summary["n_genes_full"]) if de_summary["enabled"] else 0,
            "full_table": de_summary["full_table"] if de_summary["enabled"] else "",
            "top_table": de_summary["top_table"] if de_summary["enabled"] else "",
            "volcano_plot": de_summary["volcano_plot"] if de_summary["enabled"] else "",
        },
        "disclaimer": DISCLAIMER,
    }
    if doublet_summary is not None:
        doublet_table = table_paths.get("doublet_summary")
        data["doublet"] = {
            **doublet_summary,
            "table": doublet_table.name if doublet_table else "",
        }
    if annotation_info is not None:
        annotation_path = table_paths.get("cluster_annotations")
        data["annotation"] = {
            **annotation_info,
            "table": annotation_path.name if annotation_path else "",
        }

    write_result_json(
        output_dir=output_dir,
        skill="scrna",
        version="0.1.0",
        summary=summary,
        data=data,
        input_checksum=compute_input_checksum(input_source),
    )

    write_reproducibility(
        output_dir,
        input_path,
        input_source,
        is_demo,
        args,
        table_paths=table_paths,
        figure_paths=figure_paths,
    )

    return {
        "report_path": report_path,
        "output_dir": output_dir,
        "n_clusters": n_clusters,
        "n_cells_after": n_cells_analyzed,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "ClawBio scRNA Orchestrator — Scanpy QC/clustering/markers with optional "
            "doublet detection, CellTypist annotation, and two-group DE"
        ),
    )
    parser.add_argument(
        "--input",
        "-i",
        help="Input raw-count .h5ad, matrix.mtx(.gz), or 10x Matrix Market directory",
    )
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
    parser.add_argument(
        "--doublet-method",
        choices=("none", "scrublet"),
        default="none",
        help="Optional doublet detection method",
    )
    parser.add_argument(
        "--annotate",
        choices=("none", "celltypist"),
        default="none",
        help="Optional cell type annotation backend",
    )
    parser.add_argument(
        "--annotation-model",
        default=DEFAULT_CELLTYPIST_MODEL,
        help="Local CellTypist model name or path (used with --annotate celltypist)",
    )
    parser.add_argument("--de-groupby", default=None, help="obs column for two-group DE")
    parser.add_argument("--de-group1", default=None, help="Group 1 value for DE")
    parser.add_argument("--de-group2", default=None, help="Group 2 reference value for DE")
    parser.add_argument("--de-top-genes", type=int, default=50, help="Top DE genes to include in summary table")
    parser.add_argument("--de-volcano", action="store_true", help="Generate optional DE volcano plot")
    return parser.parse_args()


def main() -> None:
    args = parse_args()
    if not args.demo and not args.input:
        print("ERROR: Provide --input <input.h5ad|matrix.mtx|10x_dir> or --demo", file=sys.stderr)
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
