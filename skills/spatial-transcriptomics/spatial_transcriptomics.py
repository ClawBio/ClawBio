#!/usr/bin/env python3
"""Visium spatial transcriptomics: QC, clustering, markers, spatial stats.

Usage:
    python spatial_transcriptomics.py --input <outs_or_h5ad> --output <dir>
    python spatial_transcriptomics.py --demo --output /tmp/spatial_demo
"""

from __future__ import annotations

import argparse
import csv
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

import numpy as np

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.report import DISCLAIMER
from clawbio.common.reproducibility import (
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

SKILL_DIR = Path(__file__).resolve().parent
SKILL_VERSION = "0.1.0"
MIN_SPOTS = 10
PUBLIC_VISIUM_ID = "V1_Human_Lymph_Node"
PUBLIC_VISIUM_MATRIX_URL = (
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/"
    "V1_Human_Lymph_Node/V1_Human_Lymph_Node_filtered_feature_bc_matrix.tar.gz"
)
PUBLIC_VISIUM_SPATIAL_URL = (
    "https://cf.10xgenomics.com/samples/spatial-exp/1.1.0/"
    "V1_Human_Lymph_Node/V1_Human_Lymph_Node_spatial.tar.gz"
)
DEMO_N_ROW = 8
DEMO_N_COL = 8
DEMO_N_GENES = 40
DEMO_SEED = 7
MARKER_LEFT = ("EPCAM", "KRT8", "KRT18", "CDH1")
MARKER_RIGHT = ("COL1A1", "VIM", "DCN", "LUM")
SPATIAL_NEIGHBORS = 6
NHOOD_PERMS = 50
CO_OCCURRENCE_BINS = 6

REPLAY_PIP_DEPENDENCIES = (
    "scanpy>=1.10",
    "leidenalg>=0.10",
    "numpy>=1.24",
    "pandas>=2.0",
    "matplotlib>=3.7",
    "scikit-learn>=1.3",
    "scipy>=1.10",
    "opentelemetry-sdk>=1.20,<2",
)


class InsufficientSpotsError(ValueError):
    """Raised when too few in-tissue spots remain after load/QC."""


def generate_demo_adata(*, seed: int = DEMO_SEED):
    """Two-domain Visium-like grid. Offline, no download."""
    import anndata as ad
    from scipy import sparse

    rng = np.random.default_rng(seed)
    n_row, n_col = DEMO_N_ROW, DEMO_N_COL
    n_spots = n_row * n_col
    genes = list(MARKER_LEFT) + list(MARKER_RIGHT)
    genes += [f"GENE{i:02d}" for i in range(len(genes), DEMO_N_GENES)]
    n_genes = len(genes)

    rows, cols, barcodes, domains = [], [], [], []
    for i in range(n_row):
        for j in range(n_col):
            rows.append(i)
            cols.append(j)
            barcodes.append(f"SPOT-{i}-{j}")
            domains.append("epithelium" if j < n_col // 2 else "stroma")

    X = rng.poisson(2.0, size=(n_spots, n_genes)).astype(np.int32)
    left_idx = [genes.index(g) for g in MARKER_LEFT]
    right_idx = [genes.index(g) for g in MARKER_RIGHT]
    for s, domain in enumerate(domains):
        if domain == "epithelium":
            X[s, left_idx] = rng.poisson(28.0, size=len(left_idx))
        else:
            X[s, right_idx] = rng.poisson(28.0, size=len(right_idx))

    adata = ad.AnnData(sparse.csr_matrix(X))
    adata.obs_names = np.array(barcodes, dtype=object)
    adata.var_names = np.array(genes, dtype=object)
    adata.obs["array_row"] = np.asarray(rows, dtype=int)
    adata.obs["array_col"] = np.asarray(cols, dtype=int)
    adata.obs["in_tissue"] = np.ones(n_spots, dtype=int)
    adata.obs["domain"] = np.array(domains, dtype=object)
    adata.obsm["spatial"] = np.column_stack(
        [np.asarray(cols, dtype=float) * 100.0, np.asarray(rows, dtype=float) * 100.0]
    )
    adata.uns["spatial_transcriptomics_demo"] = True
    return adata


def write_visium_outs(adata, dest: Path) -> Path:
    """Write a minimal SpaceRanger `outs/` tree from an AnnData object."""
    import gzip
    from scipy import sparse
    from scipy.io import mmwrite

    dest = Path(dest)
    mtx_dir = dest / "filtered_feature_bc_matrix"
    spatial_dir = dest / "spatial"
    mtx_dir.mkdir(parents=True, exist_ok=True)
    spatial_dir.mkdir(parents=True, exist_ok=True)

    X = adata.X
    if sparse.issparse(X):
        gene_by_cell = X.T.tocsc()
    else:
        gene_by_cell = sparse.csc_matrix(np.asarray(X).T)
    with gzip.open(mtx_dir / "matrix.mtx.gz", "wb") as handle:
        mmwrite(handle, gene_by_cell)

    with gzip.open(mtx_dir / "barcodes.tsv.gz", "wt", encoding="utf-8") as handle:
        handle.write("\n".join(map(str, adata.obs_names)) + "\n")
    with gzip.open(mtx_dir / "features.tsv.gz", "wt", encoding="utf-8", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t")
        for name in adata.var_names:
            writer.writerow([f"{name}.1", name, "Gene Expression"])

    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    with (spatial_dir / "tissue_positions.csv").open(
        "w", encoding="utf-8", newline=""
    ) as fh:
        writer = csv.writer(fh)
        writer.writerow(
            [
                "barcode",
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ]
        )
        for i, barcode in enumerate(adata.obs_names):
            writer.writerow(
                [
                    barcode,
                    int(adata.obs["in_tissue"].iloc[i]) if "in_tissue" in adata.obs else 1,
                    int(adata.obs["array_row"].iloc[i]) if "array_row" in adata.obs else i,
                    int(adata.obs["array_col"].iloc[i]) if "array_col" in adata.obs else i,
                    int(coords[i, 1]),
                    int(coords[i, 0]),
                ]
            )
    (spatial_dir / "scalefactors_json.json").write_text(
        json.dumps(
            {
                "spot_diameter_fullres": 50.0,
                "tissue_hires_scalef": 1.0,
                "fiducial_diameter_fullres": 100.0,
                "tissue_lowres_scalef": 1.0,
            }
        )
        + "\n",
        encoding="utf-8",
    )
    return dest


def _resolve_visium_root(path: Path) -> Path:
    path = Path(path)
    if (path / "spatial").is_dir() and (
        (path / "filtered_feature_bc_matrix").is_dir()
        or (path / "filtered_feature_bc_matrix.h5").exists()
    ):
        return path
    if (path / "outs").is_dir():
        return _resolve_visium_root(path / "outs")
    raise ValueError(
        f"{path} is not a SpaceRanger outs directory: expected spatial/ plus "
        "filtered_feature_bc_matrix/ (or the .h5 matrix)."
    )


def _load_positions(spatial_dir: Path, barcodes) -> tuple[np.ndarray, np.ndarray]:
    positions_csv = spatial_dir / "tissue_positions.csv"
    positions_list = spatial_dir / "tissue_positions_list.csv"
    if positions_csv.exists():
        import pandas as pd

        table = pd.read_csv(positions_csv)
    elif positions_list.exists():
        import pandas as pd

        table = pd.read_csv(
            positions_list,
            header=None,
            names=[
                "barcode",
                "in_tissue",
                "array_row",
                "array_col",
                "pxl_row_in_fullres",
                "pxl_col_in_fullres",
            ],
        )
    else:
        raise ValueError(f"No tissue_positions.csv under {spatial_dir}")

    table["barcode"] = table["barcode"].astype(str)
    table = table.set_index("barcode")
    missing = [b for b in barcodes if b not in table.index]
    if missing:
        raise ValueError(
            f"{len(missing)} matrix barcodes have no spatial position "
            f"(example {missing[0]})."
        )
    table = table.loc[list(barcodes)]
    coords = np.column_stack(
        [
            table["pxl_col_in_fullres"].to_numpy(dtype=float),
            table["pxl_row_in_fullres"].to_numpy(dtype=float),
        ]
    )
    in_tissue = table["in_tissue"].to_numpy(dtype=int)
    return coords, in_tissue


def public_visium_cache_dir() -> Path:
    """Durable cache for the public 10x lymph-node Visium outs tree."""
    return Path.home() / ".cache" / "clawbio" / "visium" / PUBLIC_VISIUM_ID / "outs"


_VISIUM_TAR_PREFIXES = ("filtered_feature_bc_matrix/", "spatial/")


def _safe_extract_visium_tar(tar, dest: Path) -> None:
    """Extract a 10x Visium tarball after rejecting path traversal and links.

    CodeQL flags raw ``extractall`` of a downloaded archive. Members must stay
    under ``dest`` and start with the two SpaceRanger prefixes this skill reads.
    """
    import tarfile

    dest = dest.resolve()
    for member in tar.getmembers():
        name = member.name.replace("\\", "/")
        while name.startswith("./"):
            name = name[2:]
        if not name or name.startswith("/") or ".." in name.split("/"):
            raise OSError(f"Refusing unsafe tar member {member.name!r}")
        if member.issym() or member.islnk():
            raise OSError(f"Refusing link member {member.name!r}")
        if not any(name == prefix.rstrip("/") or name.startswith(prefix) for prefix in _VISIUM_TAR_PREFIXES):
            raise OSError(f"Unexpected tar member {member.name!r}")
        target = (dest / name).resolve()
        if target != dest and dest not in target.parents:
            raise OSError(f"Tar member escapes destination: {member.name!r}")
        member.name = name
    extract_kw = {}
    if hasattr(tarfile, "data_filter"):
        extract_kw["filter"] = "data"
    tar.extractall(dest, **extract_kw)


def _visium_outs_ready(dest: Path) -> bool:
    matrix = dest / "filtered_feature_bc_matrix" / "matrix.mtx.gz"
    spatial = dest / "spatial"
    positions = spatial / "tissue_positions.csv"
    positions_list = spatial / "tissue_positions_list.csv"
    return matrix.is_file() and spatial.is_dir() and (positions.is_file() or positions_list.is_file())


def ensure_public_visium_outs(dest: Path | None = None) -> Path:
    """Download 10x Human Lymph Node Visium filtered matrix + spatial tarballs.

    Extracts so ``filtered_feature_bc_matrix/`` and ``spatial/`` sit together.
    Reuses the cache when those files already exist. Network failure raises
    ``OSError`` rather than synthesizing an outs tree.
    """
    import tarfile
    import urllib.request

    dest = Path(dest) if dest is not None else public_visium_cache_dir()
    dest = dest.resolve()
    if _visium_outs_ready(dest):
        return dest

    dest.mkdir(parents=True, exist_ok=True)
    tarball_dir = dest.parent / "tarballs"
    tarball_dir.mkdir(parents=True, exist_ok=True)
    downloads = (
        (PUBLIC_VISIUM_MATRIX_URL, tarball_dir / "filtered_feature_bc_matrix.tar.gz"),
        (PUBLIC_VISIUM_SPATIAL_URL, tarball_dir / "spatial.tar.gz"),
    )
    try:
        for url, path in downloads:
            if path.is_file() and path.stat().st_size > 0:
                continue
            urllib.request.urlretrieve(url, path)
            if path.stat().st_size == 0:
                raise OSError(f"Downloaded empty file from {url}")
        for _url, path in downloads:
            with tarfile.open(path, "r:gz") as tar:
                _safe_extract_visium_tar(tar, dest)
    except OSError:
        raise
    except Exception as exc:  # noqa: BLE001 — surface as download failure
        raise OSError(f"Failed to download public Visium outs from 10x Genomics: {exc}") from exc
    if not _visium_outs_ready(dest):
        raise OSError(
            f"Extracted {dest} but did not find filtered_feature_bc_matrix/matrix.mtx.gz "
            "plus spatial/tissue_positions*.csv"
        )
    return dest


def load_spatial(path: Path):
    """Load Visium `outs/` or an h5ad that already carries `obsm['spatial']`."""
    import anndata as ad
    import scanpy as sc

    path = Path(path)
    if path.is_file() and path.suffix == ".h5ad":
        adata = ad.read_h5ad(path)
        if "spatial" not in adata.obsm:
            raise ValueError(
                f"{path.name} has no obsm['spatial']. This skill analyses measured "
                "spot coordinates (Visium). For H&E-to-expression prediction use "
                "deepspot-m; for dissociated scRNA-seq use scrna-orchestrator."
            )
        return adata

    if not path.exists():
        raise FileNotFoundError(path)

    root = _resolve_visium_root(path)
    mtx_dir = root / "filtered_feature_bc_matrix"
    if not mtx_dir.is_dir():
        raise ValueError(
            f"{root} has spatial/ but no filtered_feature_bc_matrix/ directory. "
            "The HDF5 matrix is not read in this version; pass the mtx folder."
        )
    adata = sc.read_10x_mtx(mtx_dir, var_names="gene_symbols")
    coords, in_tissue = _load_positions(root / "spatial", adata.obs_names)
    adata.obsm["spatial"] = coords
    adata.obs["in_tissue"] = in_tissue
    adata = adata[adata.obs["in_tissue"] == 1].copy()
    return adata


def knn_indices(coords: np.ndarray, n_neighbors: int) -> np.ndarray:
    """kNN index matrix of shape (n, k), self excluded."""
    from sklearn.neighbors import NearestNeighbors

    n = coords.shape[0]
    k = min(n_neighbors, max(n - 1, 1))
    nn = NearestNeighbors(n_neighbors=k + 1, algorithm="auto")
    nn.fit(coords)
    idx = nn.kneighbors(coords, return_distance=False)
    return idx[:, 1:]


def moran_i(values: np.ndarray, knn_idx: np.ndarray) -> float:
    """Row-standardised kNN Moran's I.

    I = (zᵀ W z) / (zᵀ z) with W row-standardised, so each neighbour weight
    is 1/k. Matches the kNN path documented by Squidpy's spatial_autocorr
    (Moran 1950).
    """
    x = np.asarray(values, dtype=float)
    z = x - x.mean()
    denom = float(np.dot(z, z))
    if denom == 0.0 or knn_idx.size == 0:
        return 0.0
    lag = z[knn_idx].mean(axis=1)
    return float(np.dot(z, lag) / denom)


def neighbor_pair_counts(
    knn_idx: np.ndarray, labels: np.ndarray, clusters: list[str]
) -> np.ndarray:
    """Count directed (cluster_i, cluster_j) neighbour pairs."""
    index = {c: i for i, c in enumerate(clusters)}
    n = len(clusters)
    counts = np.zeros((n, n), dtype=float)
    for i, neighbours in enumerate(knn_idx):
        src = index[str(labels[i])]
        for j in neighbours:
            counts[src, index[str(labels[j])]] += 1.0
    return counts


def nhood_enrichment(
    knn_idx: np.ndarray,
    labels: np.ndarray,
    n_perms: int = NHOOD_PERMS,
    seed: int = DEMO_SEED,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Permutation z-score of cluster–cluster spatial neighbour counts.

    Same estimator as Squidpy `gr.nhood_enrichment`: observed neighbour-pair
    counts versus labels shuffled on the fixed spatial graph.
    """
    labels = np.asarray(labels).astype(str)
    clusters = sorted(set(labels.tolist()))
    observed = neighbor_pair_counts(knn_idx, labels, clusters)
    rng = np.random.default_rng(seed)
    null = np.stack(
        [
            neighbor_pair_counts(knn_idx, rng.permutation(labels), clusters)
            for _ in range(n_perms)
        ]
    )
    mean = null.mean(axis=0)
    std = null.std(axis=0)
    z = (observed - mean) / np.where(std == 0.0, 1.0, std)
    z = np.where(std == 0.0, 0.0, z)
    return clusters, z, observed


def co_occurrence(
    coords: np.ndarray,
    labels: np.ndarray,
    n_bins: int = CO_OCCURRENCE_BINS,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Cluster co-occurrence versus distance, normalised by cluster frequency.

    For distance bin d, score[i, j, d] = P(cluster j at distance d from i) /
    P(cluster j). Values > 1 mean j is enriched around i at that distance.
    Follows the ratio used by Squidpy `gr.co_occurrence`.
    """
    labels = np.asarray(labels).astype(str)
    clusters = sorted(set(labels.tolist()))
    n_c = len(clusters)
    n = coords.shape[0]
    if n < 3:
        empty = np.zeros((n_c, n_c, n_bins))
        return empty, np.zeros(n_bins), clusters

    diff = coords[:, None, :] - coords[None, :, :]
    dist = np.sqrt((diff ** 2).sum(axis=2))
    np.fill_diagonal(dist, np.nan)
    finite = dist[np.isfinite(dist)]
    edges = np.quantile(finite, np.linspace(0.0, 1.0, n_bins + 1))
    edges[0] = 0.0
    freq = np.array([(labels == c).mean() for c in clusters], dtype=float)

    scores = np.zeros((n_c, n_c, n_bins))
    for b in range(n_bins):
        lo, hi = edges[b], edges[b + 1]
        mask = (dist >= lo) & (dist < hi) if b < n_bins - 1 else (dist >= lo) & (dist <= hi)
        for i, ci in enumerate(clusters):
            src = labels == ci
            if not src.any():
                continue
            for j, cj in enumerate(clusters):
                if freq[j] == 0.0:
                    continue
                hits = mask[src] & (labels[None, :] == cj)
                denom = mask[src].sum()
                if denom == 0:
                    continue
                scores[i, j, b] = (hits.sum() / denom) / freq[j]
    centres = 0.5 * (edges[:-1] + edges[1:])
    return scores, centres, clusters


def _dense_column(adata, gene: str) -> np.ndarray:
    matrix = adata[:, gene].X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix).reshape(-1)


def run_pipeline(
    adata,
    *,
    min_genes: int = 5,
    min_cells: int = 1,
    n_top_hvg: int = 2000,
    n_pcs: int = 8,
    n_neighbors: int = 8,
    leiden_resolution: float = 0.5,
    random_state: int = DEMO_SEED,
    top_markers: int = 5,
):
    """QC → Leiden → Wilcoxon markers → Moran I / nhood / co-occurrence."""
    import scanpy as sc

    adata = adata.copy()
    adata.var_names_make_unique()
    n_in = adata.n_obs
    if n_in < MIN_SPOTS:
        raise InsufficientSpotsError(
            f"Only {n_in} spots loaded (minimum {MIN_SPOTS}). "
            "This is not enough for a spatial neighbourhood graph."
        )

    # percent_top defaults include 50/100/200/500; a 40-gene demo is smaller
    # than those cutoffs and Scanpy raises IndexError.
    sc.pp.calculate_qc_metrics(adata, percent_top=None, inplace=True)
    sc.pp.filter_cells(adata, min_genes=min_genes)
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if adata.n_obs < MIN_SPOTS:
        raise InsufficientSpotsError(
            f"Only {adata.n_obs} spots remain after QC (minimum {MIN_SPOTS})."
        )

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    hvg_n = min(n_top_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg_n, flavor="seurat")
    sc.pp.pca(
        adata,
        n_comps=min(n_pcs, adata.n_obs - 1, adata.n_vars - 1),
        random_state=random_state,
    )
    sc.pp.neighbors(
        adata,
        n_neighbors=min(n_neighbors, adata.n_obs - 1),
        random_state=random_state,
    )
    sc.tl.umap(adata, random_state=random_state)
    try:
        sc.tl.leiden(
            adata,
            resolution=leiden_resolution,
            random_state=random_state,
            flavor="igraph",
            n_iterations=2,
            directed=False,
        )
    except TypeError:
        sc.tl.leiden(
            adata, resolution=leiden_resolution, random_state=random_state
        )

    marker_adata = adata
    if adata.n_vars > 80 and "highly_variable" in adata.var:
        marker_adata = adata[:, adata.var["highly_variable"]].copy()
    sc.tl.rank_genes_groups(
        marker_adata, groupby="leiden", method="wilcoxon", pts=True, use_raw=False
    )

    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    knn_idx = knn_indices(coords, SPATIAL_NEIGHBORS)
    if adata.n_vars > 80 and "highly_variable" in adata.var:
        moran_adata = adata[:, adata.var["highly_variable"]]
    else:
        moran_adata = adata
    moran_X = moran_adata.X
    if hasattr(moran_X, "toarray"):
        moran_X = moran_X.toarray()
    moran_X = np.asarray(moran_X, dtype=float)
    moran_rows = [
        {"gene": str(gene), "moran_i": moran_i(moran_X[:, i], knn_idx)}
        for i, gene in enumerate(moran_adata.var_names)
    ]
    moran_rows.sort(key=lambda row: row["moran_i"], reverse=True)

    clusters, zscore, observed = nhood_enrichment(
        knn_idx, adata.obs["leiden"].to_numpy(), seed=random_state
    )
    occ, distances, occ_clusters = co_occurrence(coords, adata.obs["leiden"].to_numpy())

    names = marker_adata.uns["rank_genes_groups"]["names"]
    scores = marker_adata.uns["rank_genes_groups"]["scores"]
    pvals = marker_adata.uns["rank_genes_groups"]["pvals_adj"]
    marker_rows = []
    groups = list(names.dtype.names)
    for group in groups:
        for rank in range(min(top_markers, len(names[group]))):
            marker_rows.append(
                {
                    "cluster": str(group),
                    "gene": str(names[group][rank]),
                    "score": float(scores[group][rank]),
                    "pvals_adj": float(pvals[group][rank]),
                    "rank": rank + 1,
                }
            )

    return {
        "adata": adata,
        "knn_idx": knn_idx,
        "moran": moran_rows,
        "nhood_clusters": clusters,
        "nhood_zscore": zscore,
        "nhood_counts": observed,
        "co_occurrence": occ,
        "co_occurrence_distance": distances,
        "co_occurrence_clusters": occ_clusters,
        "markers": marker_rows,
        "n_spots_in": n_in,
    }


def _write_figures(adata, output_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    figs = output_dir / "figures"
    figs.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    leiden = adata.obs["leiden"].astype(str)
    clusters = sorted(leiden.unique())
    colour_map = {c: plt.cm.tab10(i % 10) for i, c in enumerate(clusters)}
    colours = [colour_map[c] for c in leiden]

    if "X_umap" in adata.obsm:
        umap = adata.obsm["X_umap"]
        fig, ax = plt.subplots(figsize=(5, 4))
        ax.scatter(umap[:, 0], umap[:, 1], c=colours, s=28, edgecolors="none")
        ax.set_xlabel("UMAP1")
        ax.set_ylabel("UMAP2")
        ax.set_title("Leiden clusters (UMAP)")
        fig.tight_layout()
        path = figs / "umap_leiden.png"
        fig.savefig(path, dpi=120)
        plt.close(fig)
        written.append(path)

    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    fig, ax = plt.subplots(figsize=(5, 4))
    ax.scatter(coords[:, 0], coords[:, 1], c=colours, s=36, edgecolors="none")
    ax.set_xlabel("spatial x")
    ax.set_ylabel("spatial y")
    ax.set_title("Leiden clusters (spatial)")
    ax.set_aspect("equal", adjustable="box")
    fig.tight_layout()
    path = figs / "spatial_leiden.png"
    fig.savefig(path, dpi=120)
    plt.close(fig)
    written.append(path)
    return written


def _write_tables(result: dict, output_dir: Path) -> list[Path]:
    import pandas as pd

    tables = output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []

    markers_path = tables / "markers_top.csv"
    pd.DataFrame(result["markers"]).to_csv(markers_path, index=False)
    written.append(markers_path)

    moran_path = tables / "moran_i.csv"
    pd.DataFrame(result["moran"]).to_csv(moran_path, index=False)
    written.append(moran_path)

    clusters = result["nhood_clusters"]
    z = result["nhood_zscore"]
    nhood_path = tables / "nhood_enrichment.csv"
    frame = pd.DataFrame(z, index=clusters, columns=clusters)
    frame.to_csv(nhood_path)
    written.append(nhood_path)
    return written


def _write_repro_bundle(
    output_dir: Path,
    *,
    demo: bool,
    input_path: Path | None,
    checksum_paths: list[Path],
) -> None:
    preflight: list[str] = []
    args: list = []
    if demo:
        args.append("--demo")
    else:
        assert input_path is not None
        try:
            input_path.resolve().relative_to(_PROJECT_ROOT)
            args += ["--input", ReproPath(input_path, anchor="repo_root")]
        except ValueError:
            preflight.append(
                ': "${INPUT_PATH:?Set INPUT_PATH to the Visium outs directory or h5ad used for this run}"'
            )
            args += ["--input", '"${INPUT_PATH}"']
    args += ["--output", ReproPath(output_dir, anchor="output_dir")]
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills/spatial-transcriptomics/spatial_transcriptomics.py"),
            args=args,
            comment="Replay this ClawBio spatial-transcriptomics run",
            preflight=preflight,
        ),
        repo_root=_PROJECT_ROOT,
    )
    write_environment_yml(
        output_dir,
        env_name="clawbio-spatial-transcriptomics",
        conda_deps=[],
        pip_deps=list(REPLAY_PIP_DEPENDENCIES),
        python_version="3.11",
    )
    write_checksums(checksum_paths, output_dir, anchor=output_dir)


def generate_report(
    result: dict,
    output_dir: Path,
    *,
    source_label: str,
    demo: bool,
) -> tuple[list[Path], Path, Path]:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = result["adata"]
    figure_paths = _write_figures(adata, output_dir)
    table_paths = _write_tables(result, output_dir)
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    n_clusters = adata.obs["leiden"].nunique()
    top_moran = result["moran"][:8]
    demo_note = (
        "> Demo mode. Spots are a synthetic 8×8 grid with an epithelium-like "
        "left domain (EPCAM/KRT) and a stroma-like right domain (COL1A1/VIM). "
        "No Visium dataset was downloaded.\n"
        if demo
        else ""
    )
    lines = [
        "# Spatial Transcriptomics Report" + (" (demo)" if demo else ""),
        "",
        f"**Generated**: {timestamp}  ",
        f"**Skill**: spatial-transcriptomics v{SKILL_VERSION}  ",
        f"**Input**: {source_label}  ",
        f"**Spots**: {adata.n_obs} (from {result['n_spots_in']} loaded)  ",
        f"**Genes**: {adata.n_vars}  ",
        f"**Leiden clusters**: {n_clusters}  ",
        "",
        demo_note,
        "## 1. Summary",
        "",
        (
            "Scanpy ran QC, normalisation, PCA, UMAP and Leiden clustering, then "
            + "Wilcoxon cluster markers. Spatial statistics use a k-nearest-neighbour "
            + "graph on `obsm['spatial']`: Moran's I per gene, neighbourhood "
            + "enrichment z-scores, and distance-binned cluster co-occurrence. These "
            + "are the estimators Squidpy documents (`spatial_autocorr`, "
            + "`nhood_enrichment`, `co_occurrence`); they are computed here without "
            + "importing Squidpy so the skill does not pull spatialdata/dask into the "
            + "ClawBio lockfile."
        ),
        "",
        "## 2. Spatially variable genes (Moran's I)",
        "",
        "| Gene | Moran's I |",
        "|---|---|",
    ]
    for row in top_moran:
        lines.append(f"| {row['gene']} | {row['moran_i']:.3f} |")
    lines += [
        "",
        "## 3. Cluster markers (Wilcoxon, top per Leiden cluster)",
        "",
        "| Cluster | Gene | Score | adj. p |",
        "|---|---|---|---|",
    ]
    for row in result["markers"]:
        lines.append(
            f"| {row['cluster']} | {row['gene']} | {row['score']:.2f} | {row['pvals_adj']:.3g} |"
        )
    clusters = result["nhood_clusters"]
    z = result["nhood_zscore"]
    lines += [
        "",
        "## 4. Neighbourhood enrichment (z-score vs label permutation)",
        "",
        "| From \\ To | " + " | ".join(clusters) + " |",
        "|" + "|".join(["---"] * (len(clusters) + 1)) + "|",
    ]
    for i, src in enumerate(clusters):
        cells = " | ".join(f"{z[i, j]:.2f}" for j in range(len(clusters)))
        lines.append(f"| {src} | {cells} |")
    lines += [
        "",
        (
            "Positive z on the diagonal means a cluster's spots sit next to each "
            + "other more than a random labelling of the same spatial graph would."
        ),
        "",
        "## Output files",
        "",
        "| File | Description |",
        "|---|---|",
        "| `result.json` | Machine-readable summary |",
        "| `figures/umap_leiden.png` | Leiden on UMAP |",
        "| `figures/spatial_leiden.png` | Leiden on spot coordinates |",
        "| `tables/markers_top.csv` | Wilcoxon markers |",
        "| `tables/moran_i.csv` | Moran's I for every gene |",
        "| `tables/nhood_enrichment.csv` | Neighbourhood enrichment z-scores |",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(line for line in lines if line is not None), encoding="utf-8")

    payload = {
        "skill": "spatial-transcriptomics",
        "skill_version": SKILL_VERSION,
        "demo": demo,
        "input": source_label,
        "n_spots": int(adata.n_obs),
        "n_genes": int(adata.n_vars),
        "n_clusters": int(n_clusters),
        "leiden_clusters": [str(c) for c in sorted(adata.obs["leiden"].astype(str).unique())],
        "top_moran": top_moran,
        "markers": result["markers"],
        "nhood_clusters": result["nhood_clusters"],
        "nhood_zscore": result["nhood_zscore"].tolist(),
        "co_occurrence_distance": result["co_occurrence_distance"].tolist(),
        "timestamp": timestamp,
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")

    checksum_paths = figure_paths + table_paths
    checksum_paths += [
        output_dir / "reproducibility" / "commands.sh",
        output_dir / "reproducibility" / "environment.yml",
    ]
    return checksum_paths, report_path, result_path


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Visium spatial transcriptomics (Scanpy + spatial stats)",
    )
    p.add_argument("--input", help="SpaceRanger outs/ directory or spatial h5ad")
    p.add_argument("--output", default="./spatial_output", help="Output directory")
    p.add_argument("--demo", action="store_true", help="Run on the bundled synthetic grid")
    p.add_argument("--min-genes", type=int, default=5)
    p.add_argument("--min-cells", type=int, default=1)
    p.add_argument("--leiden-resolution", type=float, default=0.5)
    p.add_argument("--n-top-hvg", type=int, default=2000)
    p.add_argument("--random-state", type=int, default=DEMO_SEED)
    return p


def main(argv: list[str] | None = None) -> None:
    args = _build_parser().parse_args(argv)
    output_dir = Path(args.output)
    if args.demo:
        adata = generate_demo_adata(seed=args.random_state)
        source_label = "synthetic 8x8 Visium-like grid (demo)"
        input_path = None
    else:
        if not args.input:
            print("ERROR: Provide --input <outs_or_h5ad> or --demo", file=sys.stderr)
            sys.exit(1)
        input_path = Path(args.input)
        adata = load_spatial(input_path)
        source_label = str(input_path.resolve())

    try:
        result = run_pipeline(
            adata,
            min_genes=args.min_genes,
            min_cells=args.min_cells,
            n_top_hvg=args.n_top_hvg,
            leiden_resolution=args.leiden_resolution,
            random_state=args.random_state,
        )
    except InsufficientSpotsError as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)

    checksum_paths, _report, _result = generate_report(
        result, output_dir, source_label=source_label, demo=args.demo
    )
    _write_repro_bundle(
        output_dir,
        demo=args.demo,
        input_path=input_path,
        checksum_paths=checksum_paths,
    )
    print(f"[spatial-transcriptomics] Done → {output_dir}/report.md")
    print(
        f"  {result['adata'].n_obs} spots, "
        f"{result['adata'].obs['leiden'].nunique()} Leiden clusters"
    )


if __name__ == "__main__":
    main()
