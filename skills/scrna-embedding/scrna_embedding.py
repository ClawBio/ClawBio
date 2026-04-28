#!/usr/bin/env python3
"""ClawBio scRNA Embedding — scVI/scANVI latent embedding and integration."""

from __future__ import annotations

import argparse
import shlex
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import numpy as np
import pandas as pd
from sklearn.metrics import silhouette_score

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file
from clawbio.common.report import (
    DISCLAIMER,
    generate_report_footer,
    generate_report_header,
    write_result_json,
)
from clawbio.common.reproducibility import (
    ReproCommand,
    ReproPath,
    write_portable_commands_sh,
)
from clawbio.common.scrna_io import compute_input_checksum, load_count_adata

EMBEDDING_ARTIFACT_KEY = "clawbio_scrna_embedding"
DEFAULT_DOWNSTREAM_REP = "X_scvi"
DEFAULT_COUNTS_LAYER = "counts"


def model_display_name(method: str) -> str:
    """Return the user-facing model name for a backend identifier."""
    if method == "scanvi":
        return "scANVI"
    return "scVI"


def _import_scanpy():
    """Import scanpy lazily with a clear user-facing error."""
    try:
        import scanpy as sc  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "scanpy is required for scrna-embedding. "
            "Install it with: pip install scanpy anndata"
        ) from exc
    return sc


def _import_scvi():
    """Import scvi-tools lazily with a clear user-facing error."""
    try:
        import scvi  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "scvi-tools is required for scrna-embedding. "
            "Install it with: pip install scvi-tools torch"
        ) from exc
    return scvi


def build_demo_adata(random_state: int):
    """Create deterministic synthetic AnnData with cluster structure and batch effect."""
    from anndata import AnnData  # type: ignore

    rng = np.random.default_rng(random_state)
    n_batches = 2
    n_clusters = 3
    cells_per_combo = 36
    n_genes = 640

    base_profiles = []
    for cluster_idx in range(n_clusters):
        base = rng.gamma(shape=2.3, scale=1.0, size=n_genes)
        marker_start = 45 + cluster_idx * 30
        marker_end = marker_start + 18
        base[marker_start:marker_end] += 7.0
        base_profiles.append(base)

    batch_effects = []
    for batch_idx in range(n_batches):
        effect = np.ones(n_genes, dtype=np.float64)
        shift_start = 220 + batch_idx * 35
        shift_end = shift_start + 24
        effect[shift_start:shift_end] *= 1.8 + 0.15 * batch_idx
        effect[shift_start + 70 : shift_end + 70] *= 0.72
        batch_effects.append(effect)

    rows = []
    truth = []
    batches = []
    for batch_idx, batch_effect in enumerate(batch_effects):
        for cluster_idx, base in enumerate(base_profiles):
            lam = np.clip(base * batch_effect, 0.05, None)
            counts = rng.poisson(lam=lam, size=(cells_per_combo, n_genes))
            libsize_scale = rng.lognormal(mean=0.0, sigma=0.28, size=(cells_per_combo, 1))
            counts = np.round(counts * libsize_scale).astype(np.int32)
            rows.append(counts)
            truth.extend([f"cluster_{cluster_idx}"] * cells_per_combo)
            batches.extend([f"batch_{batch_idx}"] * cells_per_combo)

    x = np.vstack(rows)
    obs_names = [f"cell_{idx:03d}" for idx in range(x.shape[0])]
    gene_names = [f"Gene{idx:03d}" for idx in range(n_genes)]
    for idx in range(20):
        gene_names[idx] = f"MT-GENE{idx:02d}"

    obs = pd.DataFrame(
        {
            "sample_id": obs_names,
            "demo_truth": truth,
            "demo_batch": batches,
        },
        index=obs_names,
    )
    var = pd.DataFrame(index=pd.Index(gene_names, dtype="object"))
    return AnnData(X=x, obs=obs, var=var)


def load_demo_adata(random_state: int):
    """Load deterministic local synthetic demo data."""
    return build_demo_adata(random_state), "synthetic_scrna_embedding_demo"


def load_data(input_path: str | None, demo: bool, random_state: int, layer: str | None):
    """Load AnnData from supported raw-count input or synthetic demo data."""
    sc = _import_scanpy()

    if demo:
        adata, demo_source = load_demo_adata(random_state)
        return adata, None, True, demo_source, None

    if not input_path:
        raise ValueError("Provide --input <input.h5ad|matrix.mtx|10x_dir> or --demo.")

    adata, input_source = load_count_adata(
        input_path,
        h5ad_loader=sc.read_h5ad,
        expected_input="raw-count .h5ad or 10x single-cell input",
        layer=layer,
    )
    return adata, Path(input_path), False, None, input_source


def ensure_output_dir(output_dir: Path) -> None:
    """Create an output directory, refusing to overwrite a non-empty one."""
    if output_dir.exists() and any(output_dir.iterdir()):
        raise ValueError(
            f"Output directory already exists and is not empty: {output_dir}. "
            "Choose a new directory."
        )
    output_dir.mkdir(parents=True, exist_ok=True)


def validate_args(args: argparse.Namespace) -> None:
    """Validate basic CLI invariants."""
    if args.method not in {"scvi", "scanvi"}:
        raise ValueError(f"Unsupported embedding method: {args.method}")
    if args.method == "scanvi":
        if not args.labels_key:
            raise ValueError("--labels-key is required when --method scanvi.")
        if not args.unlabeled_category:
            raise ValueError("--unlabeled-category is required when --method scanvi.")
    else:
        if args.labels_key:
            raise ValueError("--labels-key is only supported when --method scanvi.")
        if args.unlabeled_category:
            raise ValueError("--unlabeled-category is only supported when --method scanvi.")
    if args.latent_dim < 2:
        raise ValueError("--latent-dim must be >= 2.")
    if args.max_epochs < 1:
        raise ValueError("--max-epochs must be >= 1.")
    if args.n_neighbors < 2:
        raise ValueError("--n-neighbors must be >= 2.")


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


def prepare_training_data(adata_counts, n_top_hvg: int) -> tuple[Any, Any, int]:
    """Build full-gene log-normalized data and an HVG raw-count training view."""
    sc = _import_scanpy()

    adata_norm_full = adata_counts.copy()
    sc.pp.normalize_total(adata_norm_full, target_sum=1e4)
    sc.pp.log1p(adata_norm_full)
    sc.pp.highly_variable_genes(adata_norm_full, n_top_genes=n_top_hvg, flavor="seurat")

    hvg_mask = adata_norm_full.var["highly_variable"].astype(bool).to_numpy()
    n_hvg = int(hvg_mask.sum())
    if n_hvg == 0:
        raise ValueError("No highly variable genes found.")

    adata_model = adata_counts[:, hvg_mask].copy()
    return adata_norm_full, adata_model, n_hvg


def prepare_batch_obs(adata, batch_key: str | None) -> tuple[Any, str]:
    """Validate and standardize batch metadata when provided."""
    adata = adata.copy()
    if not batch_key:
        return adata, ""

    if batch_key not in adata.obs.columns:
        available_cols = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
        raise ValueError(
            f"Batch key not found in adata.obs: {batch_key}. Available columns: {available_cols}."
        )

    adata.obs[batch_key] = adata.obs[batch_key].astype(str)
    return adata, batch_key


def prepare_label_obs(
    adata,
    *,
    method: str,
    labels_key: str | None,
    unlabeled_category: str | None,
) -> tuple[Any, dict[str, Any]]:
    """Validate and standardize label metadata for scANVI."""
    adata = adata.copy()
    empty_info = {
        "labels_key": "",
        "unlabeled_category": "",
        "n_unlabeled_cells": 0,
        "unlabeled_category_present": False,
        "observed_labels": [],
    }
    if method != "scanvi":
        return adata, empty_info

    assert labels_key is not None
    assert unlabeled_category is not None
    if labels_key not in adata.obs.columns:
        available_cols = ", ".join(sorted(map(str, adata.obs.columns.tolist())))
        raise ValueError(
            f"Labels key not found in adata.obs: {labels_key}. Available columns: {available_cols}."
        )

    raw_labels = adata.obs[labels_key]
    if raw_labels.isna().any():
        raise ValueError(
            f"Labels column '{labels_key}' contains missing values. "
            "Encode unlabeled cells explicitly and set --unlabeled-category to that value."
        )

    normalized_labels = raw_labels.astype(str)
    observed_labels = list(pd.Index(pd.unique(normalized_labels), dtype="object"))
    categories = observed_labels.copy()
    if unlabeled_category not in categories:
        categories.append(unlabeled_category)

    adata.obs[labels_key] = pd.Categorical(normalized_labels, categories=categories)
    n_unlabeled_cells = int((normalized_labels == unlabeled_category).sum())
    return adata, {
        "labels_key": labels_key,
        "unlabeled_category": unlabeled_category,
        "n_unlabeled_cells": n_unlabeled_cells,
        "unlabeled_category_present": n_unlabeled_cells > 0,
        "observed_labels": observed_labels,
    }


def infer_model_device(model) -> str:
    """Infer the device that a trained scVI model used."""
    device = getattr(model, "device", None)
    if device is None and getattr(model, "module", None) is not None:
        try:
            device = next(model.module.parameters()).device
        except Exception:
            device = None

    device_text = str(device) if device is not None else "unknown"
    if device_text.startswith("cuda"):
        return "gpu"
    return device_text


def train_scvi_model(
    adata_model,
    *,
    batch_key: str | None,
    latent_dim: int,
    max_epochs: int,
    accelerator: str,
    random_state: int,
):
    """Train an scVI model on the HVG raw-count view."""
    scvi = _import_scvi()

    scvi.settings.seed = random_state
    setup_kwargs: dict[str, Any] = {}
    if batch_key:
        setup_kwargs["batch_key"] = batch_key

    scvi.model.SCVI.setup_anndata(adata_model, **setup_kwargs)
    model = scvi.model.SCVI(adata_model, n_latent=latent_dim)
    model.train(
        max_epochs=max_epochs,
        accelerator=accelerator,
        enable_progress_bar=False,
    )
    return model, infer_model_device(model)


def train_scanvi_model(
    scvi_model,
    *,
    adata_model,
    labels_key: str,
    unlabeled_category: str,
    max_epochs: int,
    accelerator: str,
):
    """Refine a trained SCVI model with scANVI supervision."""
    scvi = _import_scvi()

    model = scvi.model.SCANVI.from_scvi_model(
        scvi_model,
        adata=adata_model,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
    )
    model.train(
        max_epochs=max_epochs,
        accelerator=accelerator,
        enable_progress_bar=False,
    )
    return model, infer_model_device(model)


def train_embedding_model(
    adata_model,
    *,
    method: str,
    batch_key: str | None,
    labels_key: str | None,
    unlabeled_category: str | None,
    latent_dim: int,
    max_epochs: int,
    accelerator: str,
    random_state: int,
):
    """Train the requested latent model and return the fitted model/device."""
    scvi_model, accelerator_used = train_scvi_model(
        adata_model,
        batch_key=batch_key,
        latent_dim=latent_dim,
        max_epochs=max_epochs,
        accelerator=accelerator,
        random_state=random_state,
    )
    if method == "scvi":
        return scvi_model, accelerator_used

    assert labels_key is not None
    assert unlabeled_category is not None
    scanvi_model, accelerator_used = train_scanvi_model(
        scvi_model,
        adata_model=adata_model,
        labels_key=labels_key,
        unlabeled_category=unlabeled_category,
        max_epochs=max_epochs,
        accelerator=accelerator,
    )
    return scanvi_model, accelerator_used


def run_latent_embedding(
    adata_norm_full,
    latent: np.ndarray,
    *,
    n_neighbors: int,
    random_state: int,
):
    """Attach latent embedding and run neighbors/UMAP on top of it."""
    sc = _import_scanpy()

    adata_latent = adata_norm_full.copy()
    adata_latent.obsm["X_scvi"] = latent
    sc.pp.neighbors(adata_latent, n_neighbors=n_neighbors, use_rep="X_scvi")
    sc.tl.umap(adata_latent, random_state=random_state)
    return adata_latent


def prepare_latent_plot_labels(
    adata,
    *,
    batch_key: str,
    labels_key: str = "",
) -> tuple[Any, str]:
    """Pick or derive a categorical obs column for the latent UMAP legend."""
    sc = _import_scanpy()

    if labels_key and labels_key in adata.obs.columns:
        adata.obs[labels_key] = adata.obs[labels_key].astype(str)
        return adata, labels_key
    if "demo_truth" in adata.obs.columns:
        adata.obs["demo_truth"] = adata.obs["demo_truth"].astype(str)
        return adata, "demo_truth"

    if "scvi_latent_group" not in adata.obs.columns:
        sc.tl.leiden(adata, key_added="scvi_latent_group")
    adata.obs["scvi_latent_group"] = adata.obs["scvi_latent_group"].astype(str)
    return adata, "scvi_latent_group"


def write_tables(
    adata,
    tables_dir: Path,
    *,
    batch_key: str,
    latent_color_key: str,
) -> dict[str, Path]:
    """Write scVI-specific latent embedding export."""
    tables_dir.mkdir(parents=True, exist_ok=True)

    latent = np.asarray(adata.obsm["X_scvi"])
    latent_df = pd.DataFrame(
        latent,
        index=adata.obs_names,
        columns=[f"latent_{idx + 1}" for idx in range(latent.shape[1])],
    )
    latent_df.insert(0, "cell_id", adata.obs_names.astype(str))
    latent_df["umap_1"] = np.asarray(adata.obsm["X_umap"])[:, 0]
    latent_df["umap_2"] = np.asarray(adata.obsm["X_umap"])[:, 1]
    if batch_key:
        latent_df[batch_key] = adata.obs[batch_key].astype(str).to_numpy()
    if latent_color_key:
        latent_df[latent_color_key] = adata.obs[latent_color_key].astype(str).to_numpy()

    latent_path = tables_dir / "latent_embeddings.csv"
    latent_df.to_csv(latent_path, index=False)
    return {"latent_embeddings": latent_path}


def compute_batch_mixing_metrics(adata, *, batch_key: str) -> pd.DataFrame:
    """Compute lightweight batch mixing diagnostics from latent neighbors."""
    if not batch_key:
        return pd.DataFrame(
            [
                {"metric": "batch_key", "value": "none", "interpretation": "Batch integration not requested."},
            ]
        )

    if "connectivities" not in adata.obsp:
        raise ValueError("Neighbors graph missing; cannot compute batch mixing metrics.")

    batches = adata.obs[batch_key].astype(str)
    connectivity = adata.obsp["connectivities"].tocsr()
    diff_batch_fractions: list[float] = []
    entropies: list[float] = []

    for idx in range(connectivity.shape[0]):
        row = connectivity[idx]
        neighbor_idx = row.indices[row.indices != idx]
        if neighbor_idx.size == 0:
            continue
        neighbor_batches = batches.iloc[neighbor_idx].to_numpy()
        same_fraction = float(np.mean(neighbor_batches == batches.iloc[idx]))
        diff_batch_fractions.append(1.0 - same_fraction)

        counts = pd.Series(neighbor_batches).value_counts(normalize=True)
        entropy = float(-(counts * np.log2(counts)).sum())
        max_entropy = np.log2(max(2, int(batches.nunique())))
        entropies.append(entropy / max_entropy if max_entropy > 0 else 0.0)

    silhouette_value = np.nan
    if int(batches.nunique()) > 1 and adata.n_obs > int(batches.nunique()):
        try:
            silhouette_value = float(silhouette_score(np.asarray(adata.obsm["X_scvi"]), batches.to_numpy()))
        except Exception:
            silhouette_value = np.nan

    metrics = [
        {
            "metric": "n_batches",
            "value": int(batches.nunique()),
            "interpretation": "Number of observed batch labels.",
        },
        {
            "metric": "mean_cross_batch_neighbor_fraction",
            "value": round(float(np.mean(diff_batch_fractions)) if diff_batch_fractions else 0.0, 4),
            "interpretation": "Higher is better for integration; measures how often neighbors come from other batches.",
        },
        {
            "metric": "mean_neighbor_batch_entropy",
            "value": round(float(np.mean(entropies)) if entropies else 0.0, 4),
            "interpretation": "Higher is better for mixing; normalized entropy of batch labels across neighbors.",
        },
        {
            "metric": "batch_silhouette",
            "value": "" if np.isnan(silhouette_value) else round(silhouette_value, 4),
            "interpretation": "Lower is better for integration; silhouette of batch labels in `X_scvi`.",
        },
    ]
    return pd.DataFrame(metrics)


def plot_core_figures(
    adata,
    figures_dir: Path,
    *,
    method: str,
    batch_key: str,
    latent_color_key: str,
) -> list[Path]:
    """Create latent-space plots for the selected embedding backend."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc = _import_scanpy()
    figures_dir.mkdir(parents=True, exist_ok=True)
    method_name = model_display_name(method)

    created: list[Path] = []

    sc.pl.umap(
        adata,
        color=latent_color_key,
        legend_loc="right margin",
        title=f"{method_name} latent UMAP ({latent_color_key})",
        show=False,
    )
    plt.tight_layout()
    latent_path = figures_dir / "umap_scvi_latent.png"
    plt.savefig(latent_path, dpi=180, bbox_inches="tight")
    plt.close("all")
    created.append(latent_path)

    if batch_key:
        sc.pl.umap(
            adata,
            color=batch_key,
            legend_loc="right margin",
            title=f"{method_name} latent UMAP ({batch_key})",
            show=False,
        )
        plt.tight_layout()
        batch_path = figures_dir / "umap_scvi_batch.png"
        plt.savefig(batch_path, dpi=180, bbox_inches="tight")
        plt.close("all")
        created.append(batch_path)

    return created


def build_downstream_scrna_command(path: str = "integrated.h5ad") -> str:
    """Return the recommended downstream scRNA command."""
    return f"clawbio.py run scrna --input {path} --use-rep {DEFAULT_DOWNSTREAM_REP}"


def write_integrated_h5ad(
    adata,
    output_dir: Path,
    *,
    raw_counts_adata,
    params: dict[str, Any],
) -> Path:
    """Write integrated AnnData with latent embedding and downstream metadata."""
    path = output_dir / "integrated.h5ad"
    artifact = adata.copy()
    counts_matrix = raw_counts_adata.X
    artifact.layers[DEFAULT_COUNTS_LAYER] = (
        counts_matrix.copy() if hasattr(counts_matrix, "copy") else counts_matrix
    )
    artifact.uns[EMBEDDING_ARTIFACT_KEY] = {
        "source_skill": "scrna-embedding",
        "version": "0.1.0",
        "preferred_rep": DEFAULT_DOWNSTREAM_REP,
        "counts_layer": DEFAULT_COUNTS_LAYER,
        "x_matrix_kind": "log1p_normalized",
        "neighbors_rep": DEFAULT_DOWNSTREAM_REP,
        "downstream_command": build_downstream_scrna_command(),
        "params": {
            "method": params["method"],
            "batch_key": params["batch_key"],
            "labels_key": params["labels_key"],
            "unlabeled_category": params["unlabeled_category"],
            "latent_dim": params["latent_dim"],
            "n_neighbors": params["n_neighbors"],
        },
    }
    artifact.write_h5ad(path)
    return path


def render_report(
    *,
    output_dir: Path,
    input_source: dict[str, Any] | None,
    is_demo: bool,
    demo_source: str | None,
    qc_stats: dict[str, int],
    n_hvg: int,
    params: dict[str, Any],
    label_info: dict[str, Any],
    accelerator_used: str,
    batch_key: str,
    latent_color_key: str,
    batch_metrics_path: Path | None,
    downstream_command: str,
) -> Path:
    """Create markdown report.md."""
    method_name = model_display_name(params["method"])
    header = generate_report_header(
        title="scRNA Embedding Report",
        skill_name="scrna-embedding",
        input_files=[Path(path) for path in input_source["files"]] if input_source else [],
        extra_metadata={
            "Mode": "demo" if is_demo else "input",
            "Input format": input_source["format"] if input_source else "demo",
            "Method": params["method"],
            "Batch key": batch_key or "none",
            "Labels key": label_info["labels_key"] or "none",
            "Unlabeled category": label_info["unlabeled_category"] or "none",
            "Unlabeled cells (after QC)": str(label_info["n_unlabeled_cells"]),
            "Cells (before QC)": str(qc_stats["n_cells_before"]),
            "Cells (after QC)": str(qc_stats["n_cells_after"]),
            "Genes (after QC)": str(qc_stats["n_genes_after"]),
            "Latent dim": str(params["latent_dim"]),
            "HVG selected": str(n_hvg),
            "Accelerator used": accelerator_used,
            "Demo source": demo_source if is_demo and demo_source else "n/a",
        },
    )

    batch_section = ""
    if batch_key:
        batch_section = (
            f"![{method_name} Batch UMAP](figures/umap_scvi_batch.png)\n"
            f"- Figure file: `figures/umap_scvi_batch.png`\n"
            f"- Colored by: `{batch_key}`\n\n"
        )

    label_section = ""
    if label_info["labels_key"]:
        unlabeled_note = (
            f"Observed **{label_info['n_unlabeled_cells']}** cells matching "
            f"`{label_info['unlabeled_category']}` after QC."
            if label_info["unlabeled_category_present"]
            else f"No cells matched `{label_info['unlabeled_category']}` after QC; "
            "the run proceeded with fully labeled cells."
        )
        label_section = (
            "## Label Metadata\n\n"
            f"- Labels key: `{label_info['labels_key']}`\n"
            f"- Unlabeled category: `{label_info['unlabeled_category']}`\n"
            f"- {unlabeled_note}\n\n"
        )

    body = f"""## Summary

- Method: **{params["method"]}**
- Cells before QC: **{qc_stats["n_cells_before"]}**
- Cells after QC: **{qc_stats["n_cells_after"]}**
- Genes before QC: **{qc_stats["n_genes_before"]}**
- Genes after QC: **{qc_stats["n_genes_after"]}**
- HVGs selected: **{n_hvg}**
- Latent dimensions: **{params["latent_dim"]}**
- Accelerator used: **{accelerator_used}**

## Core Figures

![{method_name} Latent UMAP](figures/umap_scvi_latent.png)
- Figure file: `figures/umap_scvi_latent.png`
- Colored by: `{latent_color_key}`

{batch_section}
{label_section}## Tables

- `tables/latent_embeddings.csv`
{f"- `tables/{batch_metrics_path.name}`" if batch_metrics_path is not None else ""}

## Key Outputs

- `integrated.h5ad` with `obsm["X_scvi"]`, log-normalized `X`, and raw counts in `layers["counts"]`

## Downstream Workflow

- Recommended next step for clustering, annotation, and contrastive markers:
  `{downstream_command}`
- `integrated.h5ad` carries ClawBio metadata in `uns["{EMBEDDING_ARTIFACT_KEY}"]` so `scrna-orchestrator` can auto-detect latent downstream mode

## Methods

- Input validation: raw-count `.h5ad` or 10x Matrix Market input only
- QC/filtering: `min_genes={params["min_genes"]}`, `min_cells={params["min_cells"]}`, `max_mt_pct={params["max_mt_pct"]}`
- HVG selection: `n_top_hvg={params["n_top_hvg"]}`
- Embedding method: `{params["model_class"]}`
- Training: `latent_dim={params["latent_dim"]}`, `max_epochs={params["max_epochs"]}`, `accelerator={params["accelerator"]}`
- Labels: `labels_key={label_info["labels_key"] or "none"}`, `unlabeled_category={label_info["unlabeled_category"] or "none"}`
- Neighbors graph: `use_rep="X_scvi"`, `n_neighbors={params["n_neighbors"]}`
- Visualization labels: latent UMAP uses `{latent_color_key}`; if no annotation exists, `scvi_latent_group` is derived from the latent graph for display only
- Output focus: latent embedding export, stable integrated artifact generation, and optional batch-view diagnostics only

## Reproducibility

See:
- `reproducibility/commands.sh`
- `reproducibility/environment.yml`
- `reproducibility/checksums.sha256`
"""

    report_path = output_dir / "report.md"
    report_path.write_text(header + body + generate_report_footer(), encoding="utf-8")
    return report_path


def repro_command_for_bundle(
    output_dir: Path,
    *,
    input_source: "dict[str, Any] | None",
    is_demo: bool,
    method: str,
    layer: "str | None",
    batch_key: "str | None",
    labels_key: "str | None",
    unlabeled_category: "str | None",
    min_genes: int,
    min_cells: int,
    max_mt_pct: float,
    n_top_hvg: int,
    latent_dim: int,
    max_epochs: int,
    n_neighbors: int,
    random_state: int,
    accelerator: str,
) -> ReproCommand:
    args: list = []
    if is_demo:
        args.append("--demo")
    elif input_source is not None:
        args += ["--input", ReproPath(input_source["input_path"], anchor="auto")]
    args += ["--output", ReproPath(output_dir, anchor="output_dir")]
    if method != "scvi":
        args += ["--method", method]
    if layer:
        args += ["--layer", layer]
    if batch_key:
        args += ["--batch-key", batch_key]
    if labels_key:
        args += ["--labels-key", labels_key]
    if unlabeled_category:
        args += ["--unlabeled-category", unlabeled_category]
    defaults = {
        "--min-genes": 200, "--min-cells": 3, "--max-mt-pct": 20.0,
        "--n-top-hvg": 2000, "--latent-dim": 10, "--max-epochs": 20,
        "--n-neighbors": 15, "--random-state": 0, "--accelerator": "auto",
    }
    values = {
        "--min-genes": min_genes, "--min-cells": min_cells, "--max-mt-pct": max_mt_pct,
        "--n-top-hvg": n_top_hvg, "--latent-dim": latent_dim, "--max-epochs": max_epochs,
        "--n-neighbors": n_neighbors, "--random-state": random_state, "--accelerator": accelerator,
    }
    for flag, value in values.items():
        if value != defaults[flag]:
            args += [flag, str(value)]
    return ReproCommand(
        script_path=Path("skills/scrna-embedding/scrna_embedding.py"),
        args=args,
        comment="Replay this ClawBio scrna-embedding run",
    )


def write_reproducibility(
    *,
    output_dir: Path,
    input_source: "dict[str, Any] | None",
    is_demo: bool,
    args: argparse.Namespace,
    table_paths: "dict[str, Path]",
    figure_paths: "list[Path]",
    integrated_path: Path,
) -> None:
    """Write commands.sh, environment.yml, and checksums.sha256."""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    env_yml = """name: clawbio-scrna-embedding
channels:
  - conda-forge
  - bioconda
  - defaults
dependencies:
  - python=3.11
  - scanpy=1.10.2
  - anndata=0.12.10
  - numpy=1.26.4
  - pandas=2.2.2
  - matplotlib=3.8.4
  - seaborn=0.13.2
  - leidenalg=0.10.2
  - python-igraph=0.11.6
  - pytorch
  - pip
  - pip:
      - scvi-tools
"""
    (repro_dir / "environment.yml").write_text(env_yml, encoding="utf-8")

    write_portable_commands_sh(
        output_dir,
        repro_command_for_bundle(
            output_dir,
            input_source=input_source,
            is_demo=is_demo,
            method=args.method,
            layer=args.layer,
            batch_key=args.batch_key,
            labels_key=args.labels_key,
            unlabeled_category=args.unlabeled_category,
            min_genes=args.min_genes,
            min_cells=args.min_cells,
            max_mt_pct=args.max_mt_pct,
            n_top_hvg=args.n_top_hvg,
            latent_dim=args.latent_dim,
            max_epochs=args.max_epochs,
            n_neighbors=args.n_neighbors,
            random_state=args.random_state,
            accelerator=args.accelerator,
        ),
        repo_root=_PROJECT_ROOT,
    )

    # report.md and result.json contain timestamps — exclude from checksums.
    # Only fingerprint deterministic science outputs: tables, figures, and the
    # integrated AnnData file.
    checksum_targets: list[Path] = [
        integrated_path,
        *table_paths.values(),
        *figure_paths,
    ]

    lines: list[str] = []
    for path in checksum_targets:
        if not path.exists():
            continue
        rel = path.relative_to(output_dir) if path.is_relative_to(output_dir) else path.name
        lines.append(f"{sha256_file(path)}  {rel}")
    (repro_dir / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def run_pipeline(args: argparse.Namespace) -> dict[str, Any]:
    """Run the full scVI/scANVI embedding pipeline."""
    validate_args(args)

    output_dir = Path(args.output)
    ensure_output_dir(output_dir)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"
    figures_dir.mkdir(exist_ok=True)
    tables_dir.mkdir(exist_ok=True)

    adata, input_path, is_demo, demo_source, input_source = load_data(
        args.input,
        args.demo,
        args.random_state,
        args.layer,
    )
    adata, batch_key = prepare_batch_obs(adata, args.batch_key)
    adata_qc, qc_stats = qc_filter(
        adata,
        min_genes=args.min_genes,
        min_cells=args.min_cells,
        max_mt_pct=args.max_mt_pct,
    )
    adata_qc, batch_key = prepare_batch_obs(adata_qc, batch_key)
    adata_qc, label_info = prepare_label_obs(
        adata_qc,
        method=args.method,
        labels_key=args.labels_key,
        unlabeled_category=args.unlabeled_category,
    )
    adata_norm_full, adata_model, n_hvg = prepare_training_data(adata_qc, args.n_top_hvg)
    adata_model, _ = prepare_batch_obs(adata_model, batch_key)

    model, accelerator_used = train_embedding_model(
        adata_model,
        method=args.method,
        batch_key=batch_key or None,
        labels_key=label_info["labels_key"] or None,
        unlabeled_category=label_info["unlabeled_category"] or None,
        latent_dim=args.latent_dim,
        max_epochs=args.max_epochs,
        accelerator=args.accelerator,
        random_state=args.random_state,
    )
    latent = np.asarray(model.get_latent_representation(), dtype=np.float32)
    adata_latent = run_latent_embedding(
        adata_norm_full,
        latent,
        n_neighbors=args.n_neighbors,
        random_state=args.random_state,
    )
    adata_latent, latent_color_key = prepare_latent_plot_labels(
        adata_latent,
        batch_key=batch_key,
        labels_key=label_info["labels_key"],
    )

    table_paths = write_tables(
        adata_latent,
        tables_dir,
        batch_key=batch_key,
        latent_color_key=latent_color_key,
    )
    batch_metrics_path: Path | None = None
    batch_metrics_summary: dict[str, Any] = {}
    if batch_key:
        batch_metrics = compute_batch_mixing_metrics(adata_latent, batch_key=batch_key)
        batch_metrics_path = tables_dir / "batch_mixing_metrics.csv"
        batch_metrics.to_csv(batch_metrics_path, index=False)
        table_paths["batch_mixing_metrics"] = batch_metrics_path
        metric_map = dict(zip(batch_metrics["metric"], batch_metrics["value"]))
        batch_metrics_summary = {
            "n_batches": int(metric_map.get("n_batches", 0)),
            "mean_cross_batch_neighbor_fraction": metric_map.get("mean_cross_batch_neighbor_fraction", ""),
            "mean_neighbor_batch_entropy": metric_map.get("mean_neighbor_batch_entropy", ""),
            "batch_silhouette": metric_map.get("batch_silhouette", ""),
        }
    figure_paths = plot_core_figures(
        adata_latent,
        figures_dir,
        method=args.method,
        batch_key=batch_key,
        latent_color_key=latent_color_key,
    )
    params = {
        "method": args.method,
        "model_class": "scvi.model.SCANVI" if args.method == "scanvi" else "scvi.model.SCVI",
        "layer": args.layer or "",
        "batch_key": batch_key,
        "labels_key": label_info["labels_key"],
        "unlabeled_category": label_info["unlabeled_category"],
        "min_genes": args.min_genes,
        "min_cells": args.min_cells,
        "max_mt_pct": args.max_mt_pct,
        "n_top_hvg": args.n_top_hvg,
        "latent_dim": args.latent_dim,
        "max_epochs": args.max_epochs,
        "n_neighbors": args.n_neighbors,
        "accelerator": args.accelerator,
        "random_state": args.random_state,
    }
    integrated_path = write_integrated_h5ad(
        adata_latent,
        output_dir,
        raw_counts_adata=adata_qc,
        params=params,
    )
    downstream_command = build_downstream_scrna_command(integrated_path.name)
    report_path = render_report(
        output_dir=output_dir,
        input_source=input_source,
        is_demo=is_demo,
        demo_source=demo_source,
        qc_stats=qc_stats,
        n_hvg=n_hvg,
        params=params,
        label_info=label_info,
        accelerator_used=accelerator_used,
        batch_key=batch_key,
        latent_color_key=latent_color_key,
        batch_metrics_path=batch_metrics_path,
        downstream_command=downstream_command,
    )

    tables_written = [path.name for path in table_paths.values()]
    figures_written = [path.name for path in figure_paths]
    write_result_json(
        output_dir=output_dir,
        skill="scrna-embedding",
        version="0.1.0",
        summary={
            "method": args.method,
            "input_format": input_source["format"] if input_source else "demo",
            "batch_key": batch_key,
            "labels_key": label_info["labels_key"],
            "unlabeled_category": label_info["unlabeled_category"],
            "n_unlabeled_cells": label_info["n_unlabeled_cells"],
            "unlabeled_category_present": label_info["unlabeled_category_present"],
            "cells_before_qc": qc_stats["n_cells_before"],
            "cells_after_qc": qc_stats["n_cells_after"],
            "genes_after_qc": qc_stats["n_genes_after"],
            "n_hvg": n_hvg,
            "latent_dim": args.latent_dim,
            "accelerator_used": accelerator_used,
            "latent_plot_color_by": latent_color_key,
            "batch_mixing": batch_metrics_summary,
            "downstream_scrna_command": downstream_command,
        },
        data={
            "method": args.method,
            "layer": args.layer or "",
            "input_format": input_source["format"] if input_source else "demo",
            "batch_key": batch_key,
            "labels_key": label_info["labels_key"],
            "unlabeled_category": label_info["unlabeled_category"],
            "n_unlabeled_cells": label_info["n_unlabeled_cells"],
            "unlabeled_category_present": label_info["unlabeled_category_present"],
            "observed_labels": label_info["observed_labels"],
            "cells_before_qc": qc_stats["n_cells_before"],
            "cells_after_qc": qc_stats["n_cells_after"],
            "genes_after_qc": qc_stats["n_genes_after"],
            "n_hvg": n_hvg,
            "latent_dim": args.latent_dim,
            "accelerator_used": accelerator_used,
            "latent_plot_color_by": latent_color_key,
            "batch_mixing": batch_metrics_summary,
            "tables": tables_written,
            "figures": figures_written,
            "integrated_h5ad": integrated_path.name,
            "counts_layer": DEFAULT_COUNTS_LAYER,
            "artifact_metadata_key": EMBEDDING_ARTIFACT_KEY,
            "preferred_downstream_rep": DEFAULT_DOWNSTREAM_REP,
            "downstream_scrna_command": downstream_command,
            "demo_source": demo_source if is_demo else "not_demo",
            "disclaimer": DISCLAIMER,
        },
        input_checksum=compute_input_checksum(input_source),
    )

    write_reproducibility(
        output_dir=output_dir,
        input_source=input_source,
        is_demo=is_demo,
        args=args,
        table_paths=table_paths,
        figure_paths=figure_paths,
        integrated_path=integrated_path,
    )

    return {
        "report_path": report_path,
        "output_dir": output_dir,
        "n_cells_after": qc_stats["n_cells_after"],
        "input_path": input_path,
    }


def parse_args() -> argparse.Namespace:
    parser = argparse.ArgumentParser(
        description=(
            "ClawBio scRNA Embedding — scVI/scANVI latent embedding and batch-aware diagnostics"
        ),
    )
    parser.add_argument(
        "--input",
        "-i",
        help="Input raw-count `.h5ad`, `matrix.mtx(.gz)`, or 10x Matrix Market directory",
    )
    parser.add_argument("--output", "-o", default="scrna_embedding_report", help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run local synthetic demo data")
    parser.add_argument("--method", choices=("scvi", "scanvi"), default="scvi", help="Embedding backend")
    parser.add_argument("--layer", default=None, help="Raw-count layer to use for `.h5ad` input")
    parser.add_argument("--batch-key", default=None, help="obs column for batch-aware integration")
    parser.add_argument("--labels-key", default=None, help="obs column with labels for scANVI")
    parser.add_argument(
        "--unlabeled-category",
        default=None,
        help="Category value representing unlabeled cells for scANVI",
    )
    parser.add_argument("--min-genes", type=int, default=200, help="Minimum genes per cell")
    parser.add_argument("--min-cells", type=int, default=3, help="Minimum cells per gene")
    parser.add_argument("--max-mt-pct", type=float, default=20.0, help="Maximum mitochondrial percentage")
    parser.add_argument("--n-top-hvg", type=int, default=2000, help="Number of highly variable genes")
    parser.add_argument("--latent-dim", type=int, default=10, help="Latent dimensionality for scVI/scANVI")
    parser.add_argument("--max-epochs", type=int, default=20, help="Maximum scVI/scANVI training epochs")
    parser.add_argument("--n-neighbors", type=int, default=15, help="Neighbors for graph construction")
    parser.add_argument("--random-state", type=int, default=0, help="Random seed")
    parser.add_argument(
        "--accelerator",
        choices=("auto", "cpu", "gpu", "mps"),
        default="auto",
        help="Training accelerator passed to scvi-tools",
    )
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

    print("\nscRNA Embedding complete")
    print(f"  Report: {result['report_path']}")
    print(f"  Output: {result['output_dir']}")
    print(f"  Cells after QC: {result['n_cells_after']}")


if __name__ == "__main__":
    main()
