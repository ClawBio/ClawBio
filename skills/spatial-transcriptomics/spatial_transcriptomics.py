#!/usr/bin/env python3
"""Visium spatial transcriptomics: QC, clustering, markers, spatial stats.

Usage:
    python spatial_transcriptomics.py --input <outs_or_h5ad> --output <dir>
    python spatial_transcriptomics.py --demo --output /tmp/spatial_demo
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import platform
import shlex
import sys
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

REPLAY_PACKAGES = (
    "scanpy",
    "leidenalg",
    "igraph",
    "numpy",
    "pandas",
    "matplotlib",
    "scikit-learn",
    "scipy",
    "anndata",
    "umap-learn",
    "pynndescent",
    "numba",
    "llvmlite",
    "h5py",
    "opentelemetry-sdk",
)


def _load_sibling(name: str):
    spec = importlib.util.spec_from_file_location(
        f"clawbio_{name}", SKILL_DIR / f"{name}.py"
    )
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


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
    with gzip.open(
        mtx_dir / "features.tsv.gz", "wt", encoding="utf-8", newline=""
    ) as handle:
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
                    int(adata.obs["in_tissue"].iloc[i])
                    if "in_tissue" in adata.obs
                    else 1,
                    int(adata.obs["array_row"].iloc[i])
                    if "array_row" in adata.obs
                    else i,
                    int(adata.obs["array_col"].iloc[i])
                    if "array_col" in adata.obs
                    else i,
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
_VISIUM_FILES = (
    "filtered_feature_bc_matrix/matrix.mtx.gz",
    "filtered_feature_bc_matrix/barcodes.tsv.gz",
    "filtered_feature_bc_matrix/features.tsv.gz",
    "spatial/tissue_positions.csv",
    "spatial/tissue_positions_list.csv",
    "spatial/scalefactors_json.json",
)
_PUBLIC_TARBALL_HASHES = {
    "filtered_feature_bc_matrix.tar.gz": "93f7424de945eb886db17e5184d2112c77855bd54c330b2208a846870595e4e8",
    "spatial.tar.gz": "812808883366ff9623dc8354847a7211b0d922b2bfc4c9359d6e12e993ea6a73",
}


def _safe_extract_visium_tar(tar, dest: Path) -> None:
    """Copy only known regular data files, never archive-supplied paths/metadata."""
    import shutil

    dest = dest.resolve()
    targets = {name: dest / name for name in _VISIUM_FILES}
    selected = []
    seen = set()
    for member in tar.getmembers():
        name = member.name.replace("\\", "/")
        while name.startswith("./"):
            name = name[2:]
        if not name or name.startswith("/") or ".." in name.split("/"):
            raise OSError(f"Refusing unsafe tar member {member.name!r}")
        if not (member.isfile() or member.isdir()):
            raise OSError(f"Refusing non-regular tar member {member.name!r}")
        if not any(
            name == prefix.rstrip("/") or name.startswith(prefix)
            for prefix in _VISIUM_TAR_PREFIXES
        ):
            raise OSError(f"Unexpected tar member {member.name!r}")
        if member.isdir() or name not in targets:
            continue  # Tissue images are not input to this skill.
        target = targets[name]
        if name in seen or target.is_symlink() or dest not in target.resolve().parents:
            raise OSError(f"Refusing unsafe or duplicate tar member {member.name!r}")
        if member.size > 1_000_000_000:
            raise OSError(f"Refusing oversized tar member {member.name!r}")
        seen.add(name)
        selected.append((member, target))
    for member, target in selected:
        target.parent.mkdir(parents=True, exist_ok=True)
        if dest not in target.resolve().parents:
            raise OSError("Tar member escapes destination")
        source = tar.extractfile(member)
        if source is None:
            raise OSError(f"Cannot read tar member {member.name!r}")
        with source, target.open("xb") as handle:
            shutil.copyfileobj(source, handle)


def _visium_outs_ready(dest: Path) -> bool:
    matrix = dest / "filtered_feature_bc_matrix" / "matrix.mtx.gz"
    spatial = dest / "spatial"
    positions = spatial / "tissue_positions.csv"
    positions_list = spatial / "tissue_positions_list.csv"
    required = [
        matrix,
        matrix.parent / "barcodes.tsv.gz",
        matrix.parent / "features.tsv.gz",
    ]
    return all(p.is_file() and p.stat().st_size > 0 for p in required) and any(
        p.is_file() and p.stat().st_size > 0 for p in (positions, positions_list)
    )


def ensure_public_visium_outs(dest: Path | None = None) -> Path:
    """Download 10x Human Lymph Node Visium filtered matrix + spatial tarballs.

    Extracts so ``filtered_feature_bc_matrix/`` and ``spatial/`` sit together.
    Reuses the cache when those files already exist. Network failure raises
    ``OSError`` rather than synthesizing an outs tree.
    """
    import tarfile
    import tempfile
    import urllib.request

    dest = Path(dest) if dest is not None else public_visium_cache_dir()
    dest = dest.resolve()
    if _visium_outs_ready(dest):
        return dest

    dest.parent.mkdir(parents=True, exist_ok=True)
    tarball_dir = dest.parent / "tarballs"
    tarball_dir.mkdir(parents=True, exist_ok=True)
    downloads = (
        (PUBLIC_VISIUM_MATRIX_URL, tarball_dir / "filtered_feature_bc_matrix.tar.gz"),
        (PUBLIC_VISIUM_SPATIAL_URL, tarball_dir / "spatial.tar.gz"),
    )
    try:
        for url, path in downloads:
            expected = _PUBLIC_TARBALL_HASHES[path.name]
            if path.is_file() and sha256_file(path) == expected:
                continue
            with tempfile.TemporaryDirectory(
                prefix="visium-download-", dir=tarball_dir
            ) as temporary:
                pending = Path(temporary) / path.name
                urllib.request.urlretrieve(url, pending)
                if sha256_file(pending) != expected:
                    raise OSError(
                        f"SHA-256 mismatch for public Visium archive {path.name}"
                    )
                pending.replace(path)
        with tempfile.TemporaryDirectory(
            prefix="visium-extract-", dir=dest.parent
        ) as temporary:
            staging = Path(temporary) / "outs"
            for _url, path in downloads:
                with tarfile.open(path, "r:gz") as tar:
                    _safe_extract_visium_tar(tar, staging)
            if not _visium_outs_ready(staging):
                raise OSError(
                    "Public Visium archive is missing matrix, barcodes, features or positions"
                )
            # Preserve a partial cache for inspection rather than deleting user files.
            if dest.exists():
                backup = (
                    Path(tempfile.mkdtemp(prefix="visium-incomplete-", dir=dest.parent))
                    / "outs"
                )
                dest.rename(backup)
            staging.rename(dest)
    except OSError:
        raise
    except Exception as exc:
        raise OSError(
            f"Failed to download public Visium outs from 10x Genomics: {exc}"
        ) from exc
    if not _visium_outs_ready(dest):
        raise OSError(
            f"Extracted {dest} but did not find filtered_feature_bc_matrix/matrix.mtx.gz "
            "plus spatial/tissue_positions*.csv"
        )
    return dest


def _validate_spatial_counts(adata) -> None:
    """Enforce the raw-count and two-dimensional coordinate input contract."""
    from scipy import sparse

    if "spatial" not in adata.obsm:
        raise ValueError(
            "Input has no obsm['spatial']; use scrna-orchestrator for dissociated scRNA-seq."
        )
    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    if coords.shape != (adata.n_obs, 2) or not np.isfinite(coords).all():
        raise ValueError(
            "obsm['spatial'] must contain finite x,y coordinates with shape (n_spots, 2)"
        )
    if adata.n_vars < 2 or adata.X is None:
        raise ValueError("Input must contain raw counts for at least two genes")
    values = (
        adata.X.data if sparse.issparse(adata.X) else np.asarray(adata.X).reshape(-1)
    )
    if (
        not np.isfinite(values).all()
        or np.any(values < 0)
        or not np.allclose(values, np.rint(values), rtol=0, atol=1e-6)
    ):
        raise ValueError(
            "Expected finite, nonnegative integer raw counts; use --counts-layer for processed h5ad"
        )
    if "log1p" in adata.uns:
        raise ValueError(
            "Input is marked log-transformed; select raw counts with --counts-layer"
        )
    if not adata.obs_names.is_unique:
        raise ValueError("Spot barcodes must be unique")


def load_spatial(path: Path, *, counts_layer: str | None = None):
    """Load Visium `outs/` or an h5ad that already carries `obsm['spatial']`."""
    import anndata as ad
    import scanpy as sc

    path = Path(path)
    if path.is_file() and path.suffix == ".h5ad":
        adata, _source = load_count_adata(
            path,
            h5ad_loader=ad.read_h5ad,
            expected_input="raw measured spot counts",
            layer=counts_layer,
        )
        if "spatial" not in adata.obsm:
            raise ValueError(
                f"{path.name} has no obsm['spatial']. This skill analyses measured "
                "spot coordinates (Visium). For H&E-to-expression prediction use "
                "deepspot-m; for dissociated scRNA-seq use scrna-orchestrator."
            )
        if counts_layer:
            adata.uns.pop("log1p", None)
        adata.uns["spatial_transcriptomics_count_source"] = (
            f"layers[{counts_layer}]" if counts_layer else "X"
        )
        _validate_spatial_counts(adata)
        return adata

    if not path.exists():
        raise FileNotFoundError(path)
    if counts_layer:
        raise ValueError("--counts-layer is only supported for h5ad input")

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
    adata.uns["spatial_transcriptomics_count_source"] = "SpaceRanger matrix"
    _validate_spatial_counts(adata)
    return adata


_stats = _load_sibling("spatial_stats")
knn_indices = _stats.knn_indices
moran_i = _stats.moran_i
neighbor_pair_counts = _stats.neighbor_pair_counts
nhood_enrichment = _stats.nhood_enrichment
co_occurrence = _stats.co_occurrence


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
    nhood_perms: int = NHOOD_PERMS,
    max_pct_mt: float | None = None,
):
    """QC → Leiden → Wilcoxon markers → Moran I / nhood / co-occurrence."""
    import scanpy as sc

    _validate_spatial_counts(adata)
    for name, value, minimum in (
        ("min_genes", min_genes, 1),
        ("min_cells", min_cells, 1),
        ("n_top_hvg", n_top_hvg, 2),
        ("n_pcs", n_pcs, 1),
        ("n_neighbors", n_neighbors, 2),
        ("nhood_perms", nhood_perms, 2),
        ("top_markers", top_markers, 1),
        ("random_state", random_state, 0),
    ):
        if value < minimum:
            raise ValueError(f"{name} must be at least {minimum}")
    if not np.isfinite(leiden_resolution) or leiden_resolution <= 0:
        raise ValueError("leiden_resolution must be finite and positive")
    if max_pct_mt is not None and not 0 <= max_pct_mt <= 100:
        raise ValueError("max_pct_mt must be between 0 and 100")
    adata = adata.copy()
    adata.var_names_make_unique()
    n_in = adata.n_obs
    n_genes_in = adata.n_vars
    count_source = adata.uns.get("spatial_transcriptomics_count_source", "X")
    parameters = {
        "min_genes": min_genes,
        "min_cells": min_cells,
        "n_top_hvg": n_top_hvg,
        "n_pcs": n_pcs,
        "n_neighbors": n_neighbors,
        "leiden_resolution": leiden_resolution,
        "random_state": random_state,
        "top_markers": top_markers,
        "nhood_perms": nhood_perms,
        "max_pct_mt": max_pct_mt,
    }
    warnings = []
    if n_in < MIN_SPOTS:
        raise InsufficientSpotsError(
            f"Only {n_in} spots loaded (minimum {MIN_SPOTS}). "
            "This is not enough for a spatial neighbourhood graph."
        )

    # percent_top defaults include 50/100/200/500; a 40-gene demo is smaller
    # than those cutoffs and Scanpy raises IndexError.
    adata.var["mt"] = np.asarray(
        adata.var_names.str.upper().str.startswith("MT-"), dtype=bool
    )
    mt_gene_count = int(adata.var["mt"].sum())
    sc.pp.calculate_qc_metrics(adata, qc_vars=["mt"], percent_top=None, inplace=True)
    qc_table = adata.obs[["n_genes_by_counts", "total_counts", "pct_counts_mt"]].copy()
    qc_table.index.name = "barcode"
    qc_table["pct_counts_mt"] = qc_table["pct_counts_mt"].fillna(0.0)
    keep = (qc_table["n_genes_by_counts"] >= min_genes) & (qc_table["total_counts"] > 0)
    if max_pct_mt is not None and mt_gene_count:
        keep &= qc_table["pct_counts_mt"] <= max_pct_mt
    qc_table["pass_qc"] = keep
    adata = adata[keep.to_numpy()].copy()
    sc.pp.filter_genes(adata, min_cells=min_cells)
    if adata.n_obs < MIN_SPOTS:
        raise InsufficientSpotsError(
            f"Only {adata.n_obs} spots remain after QC (minimum {MIN_SPOTS})."
        )
    if adata.n_vars < 2:
        raise ValueError("Fewer than two genes remain after counts QC")
    if not mt_gene_count:
        warnings.append(
            "No MT- gene symbols were found; mitochondrial QC is unavailable for these gene identifiers."
        )
        qc_table["pct_counts_mt"] = np.nan
        if max_pct_mt is not None:
            warnings.append(
                "The requested mitochondrial ceiling could not be applied because no MT- genes were identified."
            )
    elif max_pct_mt is None:
        warnings.append(
            "Mitochondrial percentages are reported but not filtered; inspect QC before biological interpretation."
        )
    retained_qc = qc_table.loc[keep]
    qc_summary = {
        "n_spots_loaded": n_in,
        "n_spots_after_qc": int(adata.n_obs),
        "n_genes_loaded": n_genes_in,
        "n_genes_after_qc": int(adata.n_vars),
        "mt_gene_count": mt_gene_count,
        "median_total_counts": float(retained_qc["total_counts"].median()),
        "median_n_genes_by_counts": float(retained_qc["n_genes_by_counts"].median()),
        "median_pct_counts_mt": float(retained_qc["pct_counts_mt"].median())
        if mt_gene_count
        else None,
        "max_pct_counts_mt": float(retained_qc["pct_counts_mt"].max())
        if mt_gene_count
        else None,
    }

    adata.layers["counts"] = adata.X.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    hvg_n = min(n_top_hvg, adata.n_vars)
    sc.pp.highly_variable_genes(adata, n_top_genes=hvg_n, flavor="seurat")
    n_hvg = int(adata.var["highly_variable"].sum())
    if n_hvg < 2:
        raise ValueError("Fewer than two variable genes; cannot compute expression PCA")
    effective_pcs = min(n_pcs, adata.n_obs - 1, n_hvg - 1)
    sc.pp.pca(
        adata,
        n_comps=effective_pcs,
        random_state=random_state,
    )
    sc.pp.neighbors(
        adata,
        n_neighbors=min(n_neighbors, adata.n_obs - 1),
        use_rep="X_pca",
        n_pcs=effective_pcs,
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
        sc.tl.leiden(adata, resolution=leiden_resolution, random_state=random_state)

    marker_adata = adata
    if adata.n_vars > 80 and "highly_variable" in adata.var:
        marker_adata = adata[:, adata.var["highly_variable"]].copy()
    group_counts = adata.obs["leiden"].value_counts()
    marker_groups = group_counts[
        (group_counts >= 2) & (group_counts < adata.n_obs)
    ].index.tolist()
    if len(marker_groups) < len(group_counts):
        warnings.append(
            "Wilcoxon markers require at least two spots per cluster and a nonempty rest group; unsupported groups are omitted."
        )
    if marker_groups:
        sc.tl.rank_genes_groups(
            marker_adata,
            groupby="leiden",
            groups=marker_groups,
            method="wilcoxon",
            pts=True,
            use_raw=False,
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
    moran_rows.sort(
        key=lambda row: row["moran_i"] if np.isfinite(row["moran_i"]) else -np.inf,
        reverse=True,
    )

    clusters, zscore, observed = nhood_enrichment(
        knn_idx, adata.obs["leiden"].to_numpy(), n_perms=nhood_perms, seed=random_state
    )
    occ, distances, occ_clusters = co_occurrence(coords, adata.obs["leiden"].to_numpy())

    marker_rows = []
    ranked = marker_adata.uns.get("rank_genes_groups") if marker_groups else None
    if ranked is not None:
        names, scores, pvals = ranked["names"], ranked["scores"], ranked["pvals_adj"]
    for group in names.dtype.names if ranked is not None else []:
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
        "parameters": parameters,
        "effective_embedding": {
            "n_pcs": effective_pcs,
            "n_neighbors": min(n_neighbors, adata.n_obs - 1),
        },
        "analysis_scope": {
            "gene_selection": "hvg" if adata.n_vars > 80 else "all",
            "n_genes_tested": len(moran_rows),
            "genes_tested": list(map(str, moran_adata.var_names)),
            "moran_inference": "descriptive_no_p_values",
        },
        "qc_summary": qc_summary,
        "qc_table": qc_table,
        "warnings": warnings,
        "count_source": count_source,
    }


generate_report = _load_sibling("spatial_report").generate_report


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
            ': "${INPUT_PATH:?Set INPUT_PATH to the Visium outs directory or h5ad used for this run}"'
        )
        args += [
            "--input",
            '"${INPUT_PATH}"',
            "--expected-input-sha256",
            input_identity["sha256"],
        ]
    for name, value in parameters.items():
        if value is not None:
            args += ["--" + name.replace("_", "-"), shlex.quote(str(value))]
    args += ["--output", '"${REPLAY_OUTPUT:-$OUTPUT_DIR/replay}"']
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path(
                "skills/spatial-transcriptomics/spatial_transcriptomics.py"
            ),
            args=args,
            comment="Replay this ClawBio spatial-transcriptomics run",
            preflight=preflight,
        ),
        repo_root=_PROJECT_ROOT,
    )
    packages = _analysis_environment()
    write_environment_yml(
        output_dir,
        env_name="clawbio-spatial-transcriptomics",
        conda_deps=[],
        pip_deps=[f"{name}=={version}" for name, version in packages.items()],
        python_version=platform.python_version(),
    )
    source_paths = [
        Path(__file__).resolve(),
        SKILL_DIR / "spatial_stats.py",
        SKILL_DIR / "spatial_report.py",
        _PROJECT_ROOT / "clawbio/common/reproducibility.py",
        _PROJECT_ROOT / "clawbio/common/scrna_io.py",
        _PROJECT_ROOT / "uv.lock",
    ]
    manifest = {
        "skill": "spatial-transcriptomics",
        "skill_version": SKILL_VERSION,
        "parameters": parameters,
        "input": input_identity,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "packages": packages,
        },
        "source_files": [
            {"path": str(p.relative_to(_PROJECT_ROOT)), "sha256": sha256_file(p)}
            for p in source_paths
            if p.is_file()
        ],
        "replay_note": "Use the same source and pinned environment; replay writes a new directory and verifies input hashes.",
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
    write_checksums(checksum_paths, output_dir, anchor=output_dir)


def _analysis_environment() -> dict[str, str]:
    """Pin the installed analysis dependency closure, excluding unrelated extras."""
    from packaging.requirements import Requirement
    from packaging.utils import canonicalize_name

    pending = list(REPLAY_PACKAGES)
    versions = {}
    while pending:
        name = canonicalize_name(pending.pop())
        if name in versions:
            continue
        versions[name] = metadata.version(name)
        for raw in metadata.requires(name) or []:
            requirement = Requirement(raw)
            if requirement.marker is None or requirement.marker.evaluate({"extra": ""}):
                pending.append(requirement.name)
    return dict(sorted(versions.items()))


def _input_identity(input_path: Path | None, *, seed: int) -> dict:
    if input_path is None:
        recipe = SKILL_DIR / "examples/demo_spec.json"
        return {
            "kind": "synthetic_demo",
            "seed": seed,
            "recipe_sha256": sha256_file(recipe),
        }
    input_path = input_path.expanduser().resolve()
    if input_path.is_file() and input_path.suffix == ".h5ad":
        source = resolve_input_source(input_path)
        anchor = input_path.parent
    else:
        anchor = _resolve_visium_root(input_path)
        source = resolve_input_source(anchor / "filtered_feature_bc_matrix")
        spatial = anchor / "spatial"
        positions = spatial / "tissue_positions.csv"
        if not positions.is_file():
            positions = spatial / "tissue_positions_list.csv"
        source["files"].append(positions)
    return {
        "kind": "measured_input",
        "sha256": compute_input_checksum(source),
        "files": [
            {
                "path": str(p.relative_to(anchor)),
                "size_bytes": p.stat().st_size,
                "sha256": sha256_file(p),
            }
            for p in source["files"]
        ],
    }


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="Visium spatial transcriptomics (Scanpy + spatial stats)",
    )
    source = p.add_mutually_exclusive_group(required=True)
    source.add_argument(
        "--input", help="SpaceRanger outs/ directory or raw-count spatial h5ad"
    )
    p.add_argument("--output", default="./spatial_output", help="Output directory")
    source.add_argument(
        "--demo", action="store_true", help="Run on the bundled synthetic grid"
    )
    p.add_argument("--min-genes", type=int, default=5)
    p.add_argument("--min-cells", type=int, default=1)
    p.add_argument("--leiden-resolution", type=float, default=0.5)
    p.add_argument("--n-top-hvg", type=int, default=2000)
    p.add_argument("--random-state", type=int, default=DEMO_SEED)
    p.add_argument("--n-pcs", type=int, default=8)
    p.add_argument(
        "--n-neighbors",
        type=int,
        default=8,
        help="Expression PCA graph neighbours; spatial graph remains k=6",
    )
    p.add_argument("--nhood-perms", type=int, default=NHOOD_PERMS)
    p.add_argument("--top-markers", type=int, default=5)
    p.add_argument(
        "--max-pct-mt",
        type=float,
        help="Optional maximum mitochondrial count percentage (MT- gene symbols)",
    )
    p.add_argument(
        "--counts-layer", help="Explicit raw-count layer in a processed h5ad"
    )
    p.add_argument(
        "--expected-input-sha256", help="Verify the recorded input digest before replay"
    )
    p.add_argument(
        "--overwrite",
        action="store_true",
        help="Explicitly allow replacing report files in a nonempty output directory",
    )
    return p


def main(argv: list[str] | None = None) -> None:
    parser = _build_parser()
    args = parser.parse_args(argv)
    output_dir = Path(args.output).expanduser().resolve()
    if output_dir.exists() and (not output_dir.is_dir() or any(output_dir.iterdir())):
        if not output_dir.is_dir() or not args.overwrite:
            parser.error(
                "Output directory is not empty; choose a new directory or use --overwrite"
            )
        print(
            f"WARNING: --overwrite will replace spatial report files in {output_dir}",
            file=sys.stderr,
        )
    if args.demo:
        if args.counts_layer or args.expected_input_sha256:
            parser.error("--counts-layer and --expected-input-sha256 require --input")
        source_label = "synthetic 8x8 Visium-like grid (demo)"
        input_path = None
    else:
        if not args.input:
            print("ERROR: Provide --input <outs_or_h5ad> or --demo", file=sys.stderr)
            sys.exit(1)
        input_path = Path(args.input).expanduser().resolve()
        source_label = str(input_path.resolve())

    try:
        identity = _input_identity(input_path, seed=args.random_state)
        if (
            args.expected_input_sha256
            and identity["sha256"] != args.expected_input_sha256
        ):
            raise ValueError(
                "Input SHA-256 differs from the recorded run; refusing replay"
            )
        adata = (
            generate_demo_adata(seed=args.random_state)
            if args.demo
            else load_spatial(input_path, counts_layer=args.counts_layer)
        )
        result = run_pipeline(
            adata,
            min_genes=args.min_genes,
            min_cells=args.min_cells,
            n_top_hvg=args.n_top_hvg,
            leiden_resolution=args.leiden_resolution,
            random_state=args.random_state,
            n_pcs=args.n_pcs,
            n_neighbors=args.n_neighbors,
            nhood_perms=args.nhood_perms,
            top_markers=args.top_markers,
            max_pct_mt=args.max_pct_mt,
        )
    except (ValueError, OSError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        raise SystemExit(2) from exc

    result["parameters"]["counts_layer"] = args.counts_layer
    result["input_provenance"] = identity
    checksum_paths, _report, _result = generate_report(
        result, output_dir, source_label=source_label, demo=args.demo
    )
    _write_repro_bundle(
        output_dir,
        demo=args.demo,
        input_path=input_path,
        checksum_paths=checksum_paths,
        parameters=result["parameters"],
        input_identity=identity,
    )
    print(f"[spatial-transcriptomics] Done → {output_dir}/report.md")
    print(
        f"  {result['adata'].n_obs} spots, "
        f"{result['adata'].obs['leiden'].nunique()} Leiden clusters"
    )


if __name__ == "__main__":
    main()
