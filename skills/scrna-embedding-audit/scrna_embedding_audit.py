#!/usr/bin/env python3
"""scrna-embedding-audit: score cell embeddings against an in-run PCA baseline.

Usage:
    python scrna_embedding_audit.py --input cells.h5ad --output <dir> \
        --embeddings X_scvi,X_geneformer --labels-key cell_type [--batch-key donor]
    python scrna_embedding_audit.py --demo --output /tmp/embedding_audit_demo
"""

from __future__ import annotations

import argparse
import importlib.util
import json
import platform
import shlex
import sys
import time
from importlib import metadata
from pathlib import Path

import numpy as np

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file
from clawbio.common.reproducibility import (
    ReproCommand,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)
from clawbio.common.scrna_io import (
    compute_input_checksum,
    load_count_adata,
    resolve_input_source,
)

SKILL_DIR = Path(__file__).resolve().parent
SKILL_NAME = "scrna-embedding-audit"
SKILL_VERSION = "0.1.0"
DEMO_SPEC = SKILL_DIR / "examples" / "demo_spec.json"

# Guards below are this skill's own choices, not upstream figures.
MIN_CELLS = 50
MIN_MAX_CELLS = 50
MIN_LISI_NEIGHBORS = 3
MIN_RESAMPLES = 10
RESAMPLE_FRACTION = 0.5
DUPLICATE_WARN_FRACTION = 0.01
EXCLUDED_OBSM_KEYS = ("X_umap", "X_tsne", "X_spatial", "X_pca")
EXCLUDED_OBSM_PREFIXES = ("X_draw_graph_",)
MIN_EMBEDDING_COLUMNS = 3
BASELINE_KEY = "baseline:pca"
CONTROL_KEY = "control:random"
PROVENANCE_USER = "user obsm"
PROVENANCE_BASELINE = "computed from counts in this run"
PROVENANCE_CONTROL = (
    "seeded random control (i.i.d. normal), excluded from verdict counts"
)

REPLAY_PACKAGES = (
    "scanpy",
    "anndata",
    "numpy",
    "pandas",
    "scipy",
    "scikit-learn",
    "matplotlib",
    "h5py",
    "pyarrow",
)

REFERENCES = {
    "scib_metrics_fixture_version": "0.6.1",
    "luecken_2022": "10.1038/s41592-021-01336-8",
    "korsunsky_2019_harmony_lisi": "10.1038/s41592-019-0619-0",
    "kedzierska_2025_zero_shot_fm_benchmark": "10.1186/s13059-025-03574-x",
}


def _load_sibling(name: str):
    spec = importlib.util.spec_from_file_location(
        f"clawbio_{name}", SKILL_DIR / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


metrics_mod = _load_sibling("embedding_audit_metrics")
report_mod = _load_sibling("embedding_audit_report")
METRIC_NAMES = metrics_mod.METRIC_NAMES


class InputError(ValueError):
    """Input that the skill refuses with a named reason (exit code 2)."""


# --- demo --------------------------------------------------------------------


def load_demo_spec() -> dict:
    return json.loads(DEMO_SPEC.read_text(encoding="utf-8"))


def _pca_numpy(matrix: np.ndarray, n_comps: int, random_state: int) -> np.ndarray:
    """Deterministic PCA by SVD on a centred dense matrix (demo only)."""
    centred = matrix - matrix.mean(axis=0)
    u, s, _ = np.linalg.svd(centred, full_matrices=False)
    comps = u[:, :n_comps] * s[:n_comps]
    # fix the sign convention so the demo does not depend on LAPACK's choice
    signs = np.sign(comps[np.argmax(np.abs(comps), axis=0), np.arange(comps.shape[1])])
    signs[signs == 0] = 1.0
    return comps * signs


def generate_demo_adata(seed: int | None = None, spec: dict | None = None):
    """Synthetic counts with a known cell-type and batch structure plus three embeddings.

    All three embeddings are computed from the simulated counts (recipes in
    ``examples/demo_spec.json``); none is drawn from the labels directly.
    """
    import anndata as ad
    import pandas as pd

    spec = spec or load_demo_spec()
    seed = spec["seed"] if seed is None else seed
    rng = np.random.default_rng(seed)
    shared = list(spec["shared_cell_types"])
    batches = list(spec["batches"])
    n_per = int(spec["cells_per_type_per_batch"])
    single = spec["single_batch_cell_type"]
    types = shared + [single["name"]]
    n_genes = int(spec["n_genes"])

    cell_type: list[str] = []
    donor: list[str] = []
    for t in shared:
        for b in batches:
            cell_type += [t] * n_per
            donor += [b] * n_per
    cell_type += [single["name"]] * int(single["n_cells"])
    donor += [single["batch"]] * int(single["n_cells"])
    labels = np.array(cell_type, dtype=object)
    batch = np.array(donor, dtype=object)
    n = labels.shape[0]

    lognorm = spec["gene_mean_log_normal"]
    base_mean = rng.lognormal(
        mean=lognorm["mean"], sigma=lognorm["sigma"], size=n_genes
    )
    k_markers = int(spec["markers_per_type"])
    marker_sets = {
        t: np.arange(i * k_markers, (i + 1) * k_markers) for i, t in enumerate(types)
    }
    n_batch_genes = int(spec["n_batch_genes"])
    batch_genes = np.arange(n_genes - n_batch_genes, n_genes)

    mean = np.tile(base_mean, (n, 1))
    for t, genes in marker_sets.items():
        mean[np.ix_(labels == t, genes)] *= float(spec["marker_fold_change"])
    shifted = batch == spec["batch_shifted_in"]
    mean[np.ix_(shifted, batch_genes)] *= float(spec["batch_fold_change"])
    dispersion = float(spec["negative_binomial_dispersion"])
    lam = rng.gamma(shape=dispersion, scale=mean / dispersion)
    counts = rng.poisson(lam).astype(np.float32)

    totals = counts.sum(axis=1, keepdims=True)
    totals[totals == 0] = 1.0
    lognormed = np.log1p(counts / totals * float(spec["normalize_target_sum"]))
    n_dims = int(spec["n_embedding_dims"])

    def per_batch_centre(m: np.ndarray) -> np.ndarray:
        out = m.copy()
        for b in batches:
            mask = batch == b
            out[mask] -= out[mask].mean(axis=0)
        return out

    x_corrected = _pca_numpy(per_batch_centre(lognormed), n_dims, seed)
    weighted = lognormed.copy()
    weighted[:, batch_genes] *= float(spec["batch_gene_weight"])
    x_batchy = _pca_numpy(weighted, n_dims, seed)
    all_markers = np.concatenate(list(marker_sets.values()))
    keep = np.setdiff1d(np.arange(n_genes), all_markers)
    x_over = _pca_numpy(per_batch_centre(lognormed[:, keep]), n_dims, seed)

    obs = pd.DataFrame(
        {
            spec["obs_columns"]["labels"]: labels,
            spec["obs_columns"]["batch"]: batch,
        },
        index=pd.Index([f"cell_{i:04d}" for i in range(n)], dtype=object),
    )
    var = pd.DataFrame(
        index=pd.Index([f"gene_{i:04d}" for i in range(n_genes)], dtype=object)
    )
    adata = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.obsm["X_corrected"] = x_corrected
    adata.obsm["X_batchy"] = x_batchy
    adata.obsm["X_overcorrected"] = x_over
    return adata


# --- input handling ----------------------------------------------------------


def _obs_values(adata, key: str, what: str) -> np.ndarray:
    import pandas as pd

    if key not in adata.obs.columns:
        available = ", ".join(map(str, adata.obs.columns[:40]))
        raise InputError(
            f"{what} column {key!r} not found in obs; available columns: {available}"
        )
    series = adata.obs[key]
    missing = pd.isna(series).to_numpy()
    values = np.array(
        [None if m else str(v) for v, m in zip(series.to_numpy(), missing)],
        dtype=object,
    )
    return values


def _is_scorable_embedding(value) -> bool:
    arr = np.asarray(value)
    return arr.ndim == 2 and np.issubdtype(arr.dtype, np.number)


def select_embeddings(
    adata, requested: list[str] | None, warnings: list[str]
) -> list[str]:
    available = list(adata.obsm.keys())
    if requested:
        missing = [k for k in requested if k not in adata.obsm]
        if missing:
            raise InputError(
                f"embedding key(s) not found in obsm: {', '.join(missing)}; available: {', '.join(available) or 'none'}"
            )
        for key in requested:
            if not _is_scorable_embedding(adata.obsm[key]):
                raise InputError(f"obsm[{key!r}] is not a 2-D numeric array")
        if "X_pca" in requested:
            warnings.append(
                "obsm['X_pca'] is scored on request; it is the same method as the baseline and is "
                "expected to tie or differ only by preprocessing choices."
            )
        return list(requested)
    chosen: list[str] = []
    excluded: list[str] = []
    for key in available:
        if not key.startswith("X_") or not _is_scorable_embedding(adata.obsm[key]):
            continue
        if key in EXCLUDED_OBSM_KEYS or key.startswith(EXCLUDED_OBSM_PREFIXES):
            excluded.append(key)
            continue
        if np.asarray(adata.obsm[key]).shape[1] < MIN_EMBEDDING_COLUMNS:
            excluded.append(key)
            continue
        chosen.append(key)
    if excluded:
        warnings.append(
            "obsm keys skipped by default (layouts or fewer than "
            f"{MIN_EMBEDDING_COLUMNS} columns; pass --embeddings to score them): {', '.join(excluded)}"
        )
    if not chosen:
        raise InputError(
            "no scorable embedding found in obsm (keys starting with 'X_' with at least "
            f"{MIN_EMBEDDING_COLUMNS} columns); available: {', '.join(available) or 'none'}"
        )
    return chosen


def load_input(path: Path, *, counts_layer: str | None):
    """Load an h5ad through the shared raw-count loader and resolve the counts source."""
    import anndata as ad
    import h5py

    if path.suffix != ".h5ad" or not path.is_file():
        raise InputError(f"--input must be an existing .h5ad file, got {path}")
    try:
        # list the layer names without loading the file (anndata's backed mode is
        # not used: it fails on some files under pandas 3)
        with h5py.File(path, "r") as handle:
            layers = list(handle["layers"].keys()) if "layers" in handle else []
    except OSError as exc:
        raise InputError(f"cannot read h5ad: {exc}") from exc
    layer = counts_layer
    if layer is None and "counts" in layers:
        layer = "counts"
    try:
        adata, source = load_count_adata(
            path,
            h5ad_loader=ad.read_h5ad,
            expected_input="raw UMI counts in layers['counts'] or X",
            layer=layer,
        )
    except ValueError as exc:
        raise InputError(str(exc)) from exc
    count_source = f"layers[{layer!r}]" if layer else "X"
    return adata, source, count_source


def _counts_dense_or_sparse(x):
    from scipy import sparse

    if sparse.issparse(x):
        return sparse.csr_matrix(x, dtype=np.float64)
    return np.asarray(x, dtype=np.float64)


def _check_counts(counts) -> tuple[float, np.ndarray]:
    """Reject negative or non-finite counts; return (non_integer_fraction, per-cell totals)."""
    from scipy import sparse

    values = counts.data if sparse.issparse(counts) else np.asarray(counts).ravel()
    if values.size and not np.all(np.isfinite(values)):
        raise InputError("counts contain non-finite values")
    if values.size and np.any(values < 0):
        raise InputError(
            "counts contain negative values; pass --counts-layer pointing at raw counts"
        )
    non_integer = (
        float(np.mean(np.abs(values - np.rint(values)) > 1e-6)) if values.size else 0.0
    )
    totals = np.asarray(counts.sum(axis=1)).ravel()
    return non_integer, totals


# --- baseline ----------------------------------------------------------------


def compute_pca_baseline(
    counts,
    var_names,
    *,
    batch: np.ndarray | None,
    n_top_hvg: int,
    n_pcs: int,
    random_state: int,
) -> tuple[np.ndarray, dict]:
    """Standard scanpy/Seurat recipe on the scored cells, from a fresh AnnData."""
    import anndata as ad
    import pandas as pd
    import scanpy as sc
    from scipy import sparse

    matrix = sparse.csr_matrix(counts, dtype=np.float64)
    var = pd.DataFrame(index=pd.Index([str(v) for v in var_names], dtype=object))
    obs = pd.DataFrame(
        index=pd.Index([str(i) for i in range(matrix.shape[0])], dtype=object)
    )
    adata = ad.AnnData(X=matrix, obs=obs, var=var)
    if batch is not None:
        adata.obs["batch"] = pd.Categorical([str(b) for b in batch])
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    n_genes = adata.n_vars
    hvg_used = n_genes
    batch_aware = False
    if n_genes > n_top_hvg:
        sc.pp.highly_variable_genes(
            adata,
            n_top_genes=n_top_hvg,
            flavor="seurat",
            batch_key="batch" if batch is not None else None,
        )
        adata = adata[:, adata.var["highly_variable"].to_numpy()].copy()
        hvg_used = adata.n_vars
        batch_aware = batch is not None
    dense = adata.X.toarray() if sparse.issparse(adata.X) else np.asarray(adata.X)
    adata.X = np.asarray(dense, dtype=np.float64)
    sc.pp.scale(adata, max_value=10)
    n_comps = int(min(n_pcs, adata.n_obs - 1, adata.n_vars - 1))
    if n_comps < 2:
        raise InputError(
            "too few cells or genes for a PCA baseline with at least 2 components"
        )
    sc.pp.pca(
        adata,
        n_comps=n_comps,
        svd_solver="arpack",
        random_state=random_state,
        dtype="float64",
        mask_var=None,
    )
    embedding = np.ascontiguousarray(adata.obsm["X_pca"], dtype=np.float64)
    recipe = {
        "name": "pca",
        "provenance": PROVENANCE_BASELINE,
        "recipe": (
            "normalize_total(1e4) -> log1p -> highly_variable_genes(seurat, n_top_genes"
            + (", batch_key)" if batch_aware else ")")
            + " when n_genes > n_top_hvg -> scale(max_value=10) -> pca(arpack, float64); "
            "scanpy/Seurat standard recipe, not scib-metrics' X_pre"
        ),
        "n_genes_input": int(n_genes),
        "n_hvg_used": int(hvg_used),
        "hvg_batch_aware": batch_aware,
        "n_pcs_used": n_comps,
    }
    return embedding, recipe


# --- audit -------------------------------------------------------------------


def run_audit(
    adata,
    *,
    labels_key: str,
    batch_key: str | None,
    embedding_keys: list[str] | None,
    exclude_labels: list[str],
    count_source: str,
    n_top_hvg: int,
    n_pcs: int,
    lisi_neighbors: int,
    graph_neighbors: int,
    transfer_k: int,
    max_cells: int,
    n_resamples: int,
    random_state: int,
) -> dict:
    timings: dict[str, float] = {}
    warnings: list[str] = []
    t0 = time.perf_counter()

    labels_all = _obs_values(adata, labels_key, "labels")
    batch_all = _obs_values(adata, batch_key, "batch") if batch_key else None
    keys = select_embeddings(adata, embedding_keys, warnings)

    counts = _counts_dense_or_sparse(adata.X)
    non_integer_fraction, totals = _check_counts(counts)
    if non_integer_fraction > 0:
        warnings.append(
            f"{non_integer_fraction:.1%} of count values are not integers (EM-based quantification?); "
            "the baseline treats them as counts"
        )

    keep = np.ones(adata.n_obs, dtype=bool)
    dropped: dict[str, int] = {}

    def drop(mask: np.ndarray, reason: str) -> None:
        mask = mask & keep
        if mask.any():
            dropped[reason] = int(mask.sum())
            keep[mask] = False

    drop(np.array([v is None for v in labels_all]), "missing_label")
    if exclude_labels:
        drop(np.isin(labels_all.astype(str), exclude_labels), "excluded_label")
    if batch_all is not None:
        drop(np.array([v is None for v in batch_all]), "missing_batch")
    drop(totals <= 0, "zero_total_counts")

    n_kept = int(keep.sum())
    if n_kept < MIN_CELLS:
        raise InputError(
            f"only {n_kept} scorable cells after dropping {dict(dropped)}; this skill needs at least {MIN_CELLS}"
        )

    labels_kept = labels_all[keep].astype(str)
    batch_kept = batch_all[keep].astype(str) if batch_all is not None else None
    subsampled = n_kept > max_cells
    if subsampled:
        local = metrics_mod.stratified_subsample(
            labels_kept, batch_kept, n_cells=max_cells, random_state=random_state
        )
        warnings.append(
            f"{n_kept} cells exceed --max-cells {max_cells}; one stratified subsample of {len(local)} cells is scored"
        )
    else:
        local = np.arange(n_kept)
    scored_idx = np.flatnonzero(keep)[local]
    labels = labels_kept[local]
    batch = batch_kept[local] if batch_kept is not None else None
    n_scored = int(scored_idx.size)
    timings["load_and_select"] = time.perf_counter() - t0

    label_counts = {
        str(k): int(v) for k, v in zip(*np.unique(labels, return_counts=True))
    }
    batch_counts = (
        {str(k): int(v) for k, v in zip(*np.unique(batch, return_counts=True))}
        if batch is not None
        else {}
    )
    single_batch_fraction = None
    if batch is not None:
        single = [
            lab for lab in label_counts if len(np.unique(batch[labels == lab])) == 1
        ]
        single_batch_fraction = (
            float(np.mean(np.isin(labels, single))) if single else 0.0
        )
        if single:
            warnings.append(
                f"labels present in a single batch: {', '.join(single)} ({single_batch_fraction:.1%} of scored cells); "
                "iLISI penalises keeping them separate"
            )

    strata_pairs = (
        list(zip(labels, batch))
        if batch is not None
        else [(lab, None) for lab in labels]
    )
    strata_all = (
        list(zip(labels_kept, batch_kept))
        if batch_kept is not None
        else [(lab, None) for lab in labels_kept]
    )
    strata: list[dict] = []
    for lab, bat in sorted(set(strata_all), key=lambda p: (p[0], p[1] or "")):
        strata.append(
            {
                "label": lab,
                "batch": bat,
                "n_total": int(sum(1 for p in strata_all if p == (lab, bat))),
                "n_scored": int(sum(1 for p in strata_pairs if p == (lab, bat))),
            }
        )

    t1 = time.perf_counter()
    baseline, baseline_info = compute_pca_baseline(
        counts[scored_idx],
        adata.var_names,
        batch=batch,
        n_top_hvg=n_top_hvg,
        n_pcs=n_pcs,
        random_state=random_state,
    )
    timings["baseline_pca"] = time.perf_counter() - t1
    control = np.random.default_rng(random_state).standard_normal(
        (n_scored, baseline.shape[1])
    )

    embeddings: dict[str, np.ndarray] = {}
    provenance: dict[str, str] = {}
    for key in keys:
        embeddings[key] = np.ascontiguousarray(
            np.asarray(adata.obsm[key])[scored_idx], dtype=np.float64
        )
        provenance[key] = PROVENANCE_USER
    embeddings[BASELINE_KEY] = baseline
    provenance[BASELINE_KEY] = PROVENANCE_BASELINE
    embeddings[CONTROL_KEY] = control
    provenance[CONTROL_KEY] = PROVENANCE_CONTROL

    half_size = len(
        metrics_mod.stratified_subsample(
            labels, batch, fraction=RESAMPLE_FRACTION, random_state=random_state + 1
        )
    )
    k_lisi = int(min(lisi_neighbors, half_size))
    k_graph = int(min(graph_neighbors, half_size))
    if k_lisi < lisi_neighbors or k_graph < graph_neighbors:
        warnings.append(
            f"neighbour counts reduced to the half-sample size {half_size}: lisi {k_lisi}, graph {k_graph}"
        )
    k_transfer = int(transfer_k)

    def score(
        x: np.ndarray, lab: np.ndarray, bat: np.ndarray | None, seed: int
    ) -> dict:
        return metrics_mod.score_embedding(
            x,
            lab,
            bat,
            lisi_neighbors=k_lisi,
            graph_neighbors=k_graph,
            transfer_k=k_transfer,
            random_state=seed,
        )

    t2 = time.perf_counter()
    full: dict[str, dict] = {
        key: score(x, labels, batch, random_state) for key, x in embeddings.items()
    }
    timings["full_scoring"] = time.perf_counter() - t2

    for key in keys:
        details = full[key]["details"]
        dup = details.get("duplicate_row_fraction", 0.0)
        if dup > DUPLICATE_WARN_FRACTION:
            warnings.append(
                f"{key}: {dup:.1%} of cells share identical coordinates; k-NN ties make metrics order-dependent"
            )
        if details.get("lisi_degenerate_cells", 0):
            warnings.append(
                f"{key}: {details['lisi_degenerate_cells']} cells had a degenerate LISI kernel (equidistant neighbours)"
            )

    t3 = time.perf_counter()
    resample_values: dict[str, dict[str, list[float | None]]] = {
        key: {m: [] for m in METRIC_NAMES} for key in embeddings
    }
    resample_sizes: list[int] = []
    for r in range(n_resamples):
        seed_r = random_state + 1 + r
        idx = metrics_mod.stratified_subsample(
            labels, batch, fraction=RESAMPLE_FRACTION, random_state=seed_r
        )
        resample_sizes.append(int(idx.size))
        lab_r = labels[idx]
        bat_r = batch[idx] if batch is not None else None
        for key, x in embeddings.items():
            out = score(x[idx], lab_r, bat_r, seed_r)
            for m in METRIC_NAMES:
                resample_values[key][m].append(out["metrics"][m])
    timings["resampling"] = time.perf_counter() - t3

    verdicts: dict[str, dict] = {}
    for key in keys:
        per_metric: dict[str, dict] = {}
        counts_by_group = {
            g: {
                "higher": 0,
                "lower": 0,
                "not separable": 0,
                "tied": 0,
                "not evaluated": 0,
                "n_metrics": 0,
            }
            for g in ("bio", "batch", "transfer")
        }
        for m in METRIC_NAMES:
            emb_vals = resample_values[key][m]
            base_vals = resample_values[BASELINE_KEY][m]
            diffs = [
                (e - b) if (e is not None and b is not None) else None
                for e, b in zip(emb_vals, base_vals)
            ]
            v = metrics_mod.verdict(
                diffs, full[key]["metrics"][m], full[BASELINE_KEY]["metrics"][m]
            )
            if v["verdict"] == "not evaluated":
                reason = (
                    full[key]["reasons"].get(m)
                    or full[BASELINE_KEY]["reasons"].get(m)
                    or v["reason"]
                )
                v["reason"] = reason
            v["full_value"] = full[key]["metrics"][m]
            v["baseline_value"] = full[BASELINE_KEY]["metrics"][m]
            v["group"] = metrics_mod.METRIC_GROUPS[m]
            per_metric[m] = v
            counts_by_group[v["group"]][v["verdict"]] += 1
            counts_by_group[v["group"]]["n_metrics"] += 1
        verdicts[key] = {"per_metric": per_metric, "counts_by_group": counts_by_group}

    timings["total"] = time.perf_counter() - t0
    return {
        "n_cells": int(adata.n_obs),
        "cells_dropped": dropped,
        "cells_scored": n_scored,
        "subsampled": bool(subsampled),
        "strata": strata,
        "count_source": count_source,
        "non_integer_fraction": non_integer_fraction,
        "labels": {
            "key": labels_key,
            "n_unique": len(label_counts),
            "counts": label_counts,
        },
        "batch": {
            "key": batch_key,
            "n_unique": len(batch_counts),
            "counts": batch_counts,
            "fraction_cells_in_single_batch_labels": single_batch_fraction,
        },
        "baseline": baseline_info,
        "effective": {
            "lisi_neighbors": k_lisi,
            "graph_neighbors": k_graph,
            "transfer_k": k_transfer,
            "perplexity": int(k_lisi // 3),
            "n_pcs": baseline_info["n_pcs_used"],
            "n_hvg_used": baseline_info["n_hvg_used"],
            "resample_fraction": RESAMPLE_FRACTION,
            "resample_sizes": sorted(set(resample_sizes)),
            "diff_tolerance": metrics_mod.DIFF_TOLERANCE,
            "sd_factor": metrics_mod.SD_FACTOR,
            "ceiling": metrics_mod.CEILING,
        },
        "embedding_keys": keys,
        "embeddings": {
            key: {
                "shape": [int(s) for s in embeddings[key].shape],
                "provenance": provenance[key],
                "metrics": full[key]["metrics"],
                "reasons": full[key]["reasons"],
                "details": full[key]["details"],
            }
            for key in embeddings
        },
        "resamples": resample_values,
        "verdicts": verdicts,
        "metric_definitions": metrics_mod.METRIC_DEFINITIONS,
        "metric_definitions_omitted": metrics_mod.OMITTED_METRICS,
        "metric_groups": metrics_mod.METRIC_GROUPS,
        "references": REFERENCES,
        "warnings": warnings,
        "timings_seconds": {k: round(v, 3) for k, v in timings.items()},
    }


# --- reproducibility bundle --------------------------------------------------


def _analysis_environment() -> dict[str, str]:
    """Pin the installed analysis dependency closure, excluding unrelated extras."""
    from packaging.requirements import Requirement
    from packaging.utils import canonicalize_name

    pending = list(REPLAY_PACKAGES)
    versions: dict[str, str] = {}
    while pending:
        name = canonicalize_name(pending.pop())
        if name in versions:
            continue
        try:
            versions[name] = metadata.version(name)
        except metadata.PackageNotFoundError:
            continue
        for raw in metadata.requires(name) or []:
            requirement = Requirement(raw)
            if requirement.marker is None or requirement.marker.evaluate({"extra": ""}):
                pending.append(requirement.name)
    return dict(sorted(versions.items()))


def _threadpool_info() -> list[dict]:
    try:
        from threadpoolctl import threadpool_info

        return [
            {
                k: v
                for k, v in entry.items()
                if k in ("user_api", "internal_api", "num_threads", "version")
            }
            for entry in threadpool_info()
        ]
    except (
        ImportError,
        RuntimeError,
        ValueError,
    ):  # pragma: no cover - optional diagnostic
        return []


def _input_identity(input_path: Path | None, *, seed: int) -> dict:
    if input_path is None:
        return {
            "kind": "synthetic_demo",
            "seed": seed,
            "recipe_sha256": sha256_file(DEMO_SPEC),
        }
    if input_path.suffix != ".h5ad" or not input_path.is_file():
        raise InputError(f"--input must be an existing .h5ad file, got {input_path}")
    try:
        source = resolve_input_source(input_path)
    except ValueError as exc:
        raise InputError(str(exc)) from exc
    return {
        "kind": "measured_input",
        "sha256": compute_input_checksum(source),
        "files": [
            {"path": p.name, "size_bytes": p.stat().st_size, "sha256": sha256_file(p)}
            for p in source["files"]
        ],
    }


def _write_repro_bundle(
    output_dir: Path,
    *,
    demo: bool,
    input_path: Path | None,
    checksum_paths: list[Path],
    parameters: dict,
    input_identity: dict,
) -> None:
    preflight: list[str] = []
    args: list = []
    if demo:
        args.append("--demo")
    else:
        assert input_path is not None
        preflight.append(
            ': "${INPUT_PATH:?Set INPUT_PATH to the h5ad used for this run}"'
        )
        args += [
            "--input",
            '"${INPUT_PATH}"',
            "--expected-input-sha256",
            input_identity["sha256"],
        ]
    for name, value in parameters.items():
        if value is None or value is False:
            continue
        flag = "--" + name.replace("_", "-")
        if value is True:
            args.append(flag)
        elif isinstance(value, list):
            args += [flag, shlex.quote(",".join(map(str, value)))]
        else:
            args += [flag, shlex.quote(str(value))]
    args += ["--output", '"${REPLAY_OUTPUT:-$OUTPUT_DIR/replay}"']
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills/scrna-embedding-audit/scrna_embedding_audit.py"),
            args=args,
            comment="Replay this ClawBio scrna-embedding-audit run",
            preflight=preflight,
        ),
        repo_root=_PROJECT_ROOT,
    )
    packages = _analysis_environment()
    write_environment_yml(
        output_dir,
        env_name="clawbio-scrna-embedding-audit",
        conda_deps=[],
        pip_deps=[f"{name}=={version}" for name, version in packages.items()],
        python_version=platform.python_version(),
    )
    source_paths = [
        Path(__file__).resolve(),
        SKILL_DIR / "embedding_audit_metrics.py",
        SKILL_DIR / "embedding_audit_report.py",
        SKILL_DIR / "fixtures" / "scib_metrics_0.6.1_expected.json",
        _PROJECT_ROOT / "clawbio/common/reproducibility.py",
        _PROJECT_ROOT / "clawbio/common/scrna_io.py",
        _PROJECT_ROOT / "uv.lock",
    ]
    manifest = {
        "skill": SKILL_NAME,
        "skill_version": SKILL_VERSION,
        "parameters": parameters,
        "input": input_identity,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": packages,
            "threadpool_info": _threadpool_info(),
        },
        "source_files": [
            {"path": str(p.relative_to(_PROJECT_ROOT)), "sha256": sha256_file(p)}
            for p in source_paths
            if p.is_file()
        ],
        "replay_note": "Use the same source and pinned environment; replay writes a new directory and verifies the input hash.",
    }
    manifest_path = output_dir / "reproducibility/run_manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    checksum_paths = checksum_paths + [
        manifest_path,
        output_dir / "reproducibility/commands.sh",
        output_dir / "reproducibility/environment.yml",
    ]
    missing = [str(p) for p in checksum_paths if not p.is_file()]
    if missing:
        raise RuntimeError(
            "output files promised but not written: " + ", ".join(missing)
        )
    write_checksums(checksum_paths, output_dir, anchor=output_dir)


# --- CLI ---------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Score cell embeddings in an h5ad against an in-run PCA baseline (scIB-style metrics)."
    )
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument(
        "--input",
        help="h5ad with embeddings in obsm and raw counts in layers['counts'] or X",
    )
    group.add_argument(
        "--demo", action="store_true", help="run on synthetic data (offline)"
    )
    parser.add_argument(
        "--output", default="./embedding_audit_output", help="output directory"
    )
    parser.add_argument(
        "--embeddings",
        help="comma-separated obsm keys to score (default: all X_* embeddings except layouts)",
    )
    parser.add_argument(
        "--labels-key", help="obs column with cell-type labels (required with --input)"
    )
    parser.add_argument(
        "--batch-key", help="obs column with batch; without it batch metrics are null"
    )
    parser.add_argument(
        "--exclude-labels", help="comma-separated label values to drop before scoring"
    )
    parser.add_argument(
        "--counts-layer",
        help="layer holding raw counts (default: 'counts' if present, else X)",
    )
    parser.add_argument(
        "--n-top-hvg", type=int, default=2000, help="HVGs for the PCA baseline"
    )
    parser.add_argument(
        "--n-pcs", type=int, default=50, help="PCA components for the baseline"
    )
    parser.add_argument(
        "--lisi-neighbors", type=int, default=90, help="k (incl. self) for cLISI/iLISI"
    )
    parser.add_argument(
        "--graph-neighbors",
        type=int,
        default=15,
        help="k (incl. self) for graph connectivity",
    )
    parser.add_argument(
        "--transfer-k",
        type=int,
        default=15,
        help="k for leave-one-batch-out label transfer",
    )
    parser.add_argument(
        "--max-cells",
        type=int,
        default=5000,
        help="cells scored; above this a stratified subsample is used",
    )
    parser.add_argument(
        "--n-resamples", type=int, default=50, help="half-samples for the verdict rule"
    )
    parser.add_argument(
        "--random-state",
        type=int,
        default=0,
        help="seed for subsampling, KMeans and PCA",
    )
    parser.add_argument(
        "--expected-input-sha256",
        help="refuse to run if the input hash differs (replay guard)",
    )
    parser.add_argument(
        "--overwrite",
        action="store_true",
        help="allow writing into a non-empty output directory",
    )
    return parser


def _validate_bounds(args, parser: argparse.ArgumentParser) -> None:
    if args.lisi_neighbors < MIN_LISI_NEIGHBORS:
        parser.error(f"--lisi-neighbors must be at least {MIN_LISI_NEIGHBORS}")
    if args.graph_neighbors < 2:
        parser.error("--graph-neighbors must be at least 2")
    if args.transfer_k < 1:
        parser.error("--transfer-k must be at least 1")
    if args.max_cells < MIN_MAX_CELLS:
        parser.error(f"--max-cells must be at least {MIN_MAX_CELLS}")
    if args.n_resamples < MIN_RESAMPLES:
        parser.error(f"--n-resamples must be at least {MIN_RESAMPLES}")
    if args.n_top_hvg < 2 or args.n_pcs < 2:
        parser.error("--n-top-hvg and --n-pcs must be at least 2")


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    _validate_bounds(args, parser)
    output_dir = Path(args.output).expanduser().resolve()
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        if not output_dir.is_dir() or not args.overwrite:
            parser.error(
                "Output directory is not empty; choose a new directory or use --overwrite"
            )
        print(
            f"WARNING: --overwrite will replace report files in {output_dir}",
            file=sys.stderr,
        )
    if args.demo:
        if args.counts_layer or args.expected_input_sha256:
            parser.error("--counts-layer and --expected-input-sha256 require --input")
        input_path = None
        source_label = "synthetic demo (examples/demo_spec.json)"
    else:
        if not args.labels_key:
            parser.error("--labels-key is required with --input")
        input_path = Path(args.input).expanduser().resolve()
        source_label = str(input_path)

    embedding_keys = (
        [k.strip() for k in args.embeddings.split(",") if k.strip()]
        if args.embeddings
        else None
    )
    exclude_labels = (
        [k.strip() for k in args.exclude_labels.split(",") if k.strip()]
        if args.exclude_labels
        else []
    )

    try:
        identity = _input_identity(input_path, seed=args.random_state)
        if (
            args.expected_input_sha256
            and identity["sha256"] != args.expected_input_sha256
        ):
            raise InputError(
                "Input SHA-256 differs from the recorded run; refusing replay"
            )
        if args.demo:
            spec = load_demo_spec()
            adata = generate_demo_adata(seed=args.random_state, spec=spec)
            labels_key = spec["obs_columns"]["labels"]
            batch_key = spec["obs_columns"]["batch"]
            count_source = "layers['counts']"
        else:
            adata, _source, count_source = load_input(
                input_path, counts_layer=args.counts_layer
            )
            labels_key = args.labels_key
            batch_key = args.batch_key
        result = run_audit(
            adata,
            labels_key=labels_key,
            batch_key=batch_key,
            embedding_keys=embedding_keys,
            exclude_labels=exclude_labels,
            count_source=count_source,
            n_top_hvg=args.n_top_hvg,
            n_pcs=args.n_pcs,
            lisi_neighbors=args.lisi_neighbors,
            graph_neighbors=args.graph_neighbors,
            transfer_k=args.transfer_k,
            max_cells=args.max_cells,
            n_resamples=args.n_resamples,
            random_state=args.random_state,
        )
    except InputError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc

    parameters = {
        "embeddings": embedding_keys,
        "labels_key": None if args.demo else args.labels_key,
        "batch_key": None if args.demo else args.batch_key,
        "exclude_labels": exclude_labels or None,
        "counts_layer": args.counts_layer,
        "n_top_hvg": args.n_top_hvg,
        "n_pcs": args.n_pcs,
        "lisi_neighbors": args.lisi_neighbors,
        "graph_neighbors": args.graph_neighbors,
        "transfer_k": args.transfer_k,
        "max_cells": args.max_cells,
        "n_resamples": args.n_resamples,
        "random_state": args.random_state,
    }
    result["parameters"] = parameters
    result["input_provenance"] = identity
    checksum_paths, _report, _result = report_mod.write_outputs(
        result,
        output_dir,
        source_label=source_label,
        demo=args.demo,
        skill_name=SKILL_NAME,
        skill_version=SKILL_VERSION,
        baseline_key=BASELINE_KEY,
        control_key=CONTROL_KEY,
    )
    _write_repro_bundle(
        output_dir,
        demo=args.demo,
        input_path=input_path,
        checksum_paths=checksum_paths,
        parameters=parameters,
        input_identity=identity,
    )
    print(f"[{SKILL_NAME}] Done → {output_dir}/report.md")
    for key in result["embedding_keys"]:
        counts = result["verdicts"][key]["counts_by_group"]
        summary = "; ".join(
            f"{g}: {c['higher']} higher / {c['lower']} lower / {c['not separable']} not separable / {c['tied']} tied of {c['n_metrics']}"
            for g, c in counts.items()
        )
        print(f"  {key}: {summary}")


if __name__ == "__main__":
    main()
