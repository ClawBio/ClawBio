"""Shared single-cell input parsing and validation helpers."""

from __future__ import annotations

import gzip
import hashlib
from pathlib import Path
from typing import Any, Callable

import numpy as np
import pandas as pd

from clawbio.common.checksums import sha256_file


def _sample_expression_values(x, max_values: int = 200_000) -> np.ndarray:
    """Sample expression values from dense or sparse matrices without densifying."""
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


def detect_processed_input_reason(
    adata,
    *,
    expected_input: str,
    layer: str | None = None,
) -> str | None:
    """Detect whether input looks preprocessed rather than raw counts."""
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

    if layer:
        reason += f" Checked layer `{layer}` as the requested raw-count source."
    if uns_hits:
        reason += f" Found processed-analysis metadata in adata.uns: {', '.join(uns_hits)}."
    reason += (
        f" This skill expects {expected_input}. `pbmc3k_processed` is not supported; "
        "use raw counts (e.g., `scanpy.datasets.pbmc3k()`)."
    )
    return reason


def _split_10x_prefix(filename: str) -> str | None:
    if filename.endswith("matrix.mtx.gz"):
        return filename[: -len("matrix.mtx.gz")]
    if filename.endswith("matrix.mtx"):
        return filename[: -len("matrix.mtx")]
    return None


def _pick_existing(paths: list[Path]) -> Path | None:
    for path in paths:
        if path.exists():
            return path
    return None


def resolve_input_source(path: Path) -> dict[str, Any]:
    """Resolve supported single-cell input paths into a structured source descriptor."""
    path = path.expanduser().resolve()
    if not path.exists():
        raise FileNotFoundError(f"Input not found: {path}")

    suffixes = "".join(path.suffixes).lower()
    if suffixes == ".h5ad" or path.suffix.lower() == ".h5ad":
        return {
            "format": "h5ad",
            "input_path": path,
            "files": [path],
        }

    if path.is_dir():
        matrix_path = _pick_existing([path / "matrix.mtx", path / "matrix.mtx.gz"])
        if matrix_path is None:
            raise ValueError(
                "10x Matrix Market directory must contain `matrix.mtx` or `matrix.mtx.gz`."
            )
        prefix = _split_10x_prefix(matrix_path.name) or ""
    else:
        prefix = _split_10x_prefix(path.name)
        matrix_path = path
        if prefix is None:
            raise ValueError(
                "Unsupported input. Provide a raw-count `.h5ad`, a 10x directory, or `matrix.mtx(.gz)`."
            )

    parent = matrix_path.parent
    barcodes_path = _pick_existing(
        [parent / f"{prefix}barcodes.tsv", parent / f"{prefix}barcodes.tsv.gz"]
    )
    features_path = _pick_existing(
        [
            parent / f"{prefix}features.tsv",
            parent / f"{prefix}features.tsv.gz",
            parent / f"{prefix}genes.tsv",
            parent / f"{prefix}genes.tsv.gz",
        ]
    )
    missing = []
    if barcodes_path is None:
        missing.append(f"{prefix}barcodes.tsv(.gz)")
    if features_path is None:
        missing.append(f"{prefix}features.tsv(.gz) or {prefix}genes.tsv(.gz)")
    if missing:
        raise ValueError(f"Missing required 10x sidecar file(s): {', '.join(missing)}.")

    return {
        "format": "10x_mtx",
        "input_path": path,
        "matrix_path": matrix_path,
        "barcodes_path": barcodes_path,
        "features_path": features_path,
        "files": [matrix_path, barcodes_path, features_path],
    }


def load_10x_mtx_data(source_info: dict[str, Any]):
    """Load a 10x Matrix Market dataset into AnnData."""
    from anndata import AnnData  # type: ignore
    from scipy import io as scipy_io  # type: ignore
    from scipy import sparse  # type: ignore

    matrix_path = Path(source_info["matrix_path"])
    if matrix_path.suffix == ".gz":
        with gzip.open(matrix_path, "rb") as handle:
            matrix = scipy_io.mmread(handle)
    else:
        matrix = scipy_io.mmread(str(matrix_path))

    matrix = sparse.csr_matrix(matrix).transpose().tocsr()

    barcodes = pd.read_csv(
        source_info["barcodes_path"],
        sep="\t",
        header=None,
        compression="infer",
    )
    features = pd.read_csv(
        source_info["features_path"],
        sep="\t",
        header=None,
        compression="infer",
    )

    if matrix.shape[0] != len(barcodes):
        raise ValueError(
            "10x input mismatch: matrix columns do not match number of barcodes."
        )
    if matrix.shape[1] != len(features):
        raise ValueError(
            "10x input mismatch: matrix rows do not match number of features."
        )

    obs_names = pd.Index(barcodes.iloc[:, 0].astype(str), dtype="object")
    obs_names.name = None

    feature_name_col = 1 if features.shape[1] >= 2 else 0
    var_names = pd.Index(features.iloc[:, feature_name_col].astype(str), dtype="object")
    var_names.name = None

    obs = pd.DataFrame(index=obs_names)
    var = pd.DataFrame(index=var_names)
    var["gene_ids"] = features.iloc[:, 0].astype(str).to_numpy()
    if features.shape[1] >= 3:
        var["feature_types"] = features.iloc[:, 2].astype(str).to_numpy()

    adata = AnnData(X=matrix, obs=obs, var=var)
    adata.var_names_make_unique()
    return adata


def load_count_adata(
    input_path: str | Path,
    *,
    h5ad_loader: Callable[[Path], Any],
    expected_input: str,
    layer: str | None = None,
) -> tuple[Any, dict[str, Any]]:
    """Load supported count input and validate it looks like raw counts."""
    source_info = resolve_input_source(Path(input_path))
    if source_info["format"] == "h5ad":
        adata = h5ad_loader(Path(source_info["input_path"]))
        if layer:
            if layer not in adata.layers:
                raise ValueError(f"Requested layer not found in `.h5ad`: {layer}")
            adata = adata.copy()
            layer_matrix = adata.layers[layer]
            adata.X = layer_matrix.copy() if hasattr(layer_matrix, "copy") else layer_matrix
        processed_reason = detect_processed_input_reason(
            adata,
            expected_input=expected_input,
            layer=layer,
        )
        if processed_reason:
            raise ValueError(processed_reason)
        source_info["selected_layer"] = layer or ""
        return adata, source_info

    if layer:
        raise ValueError("--layer is only supported for `.h5ad` input.")

    adata = load_10x_mtx_data(source_info)
    processed_reason = detect_processed_input_reason(
        adata,
        expected_input=expected_input,
    )
    if processed_reason:
        raise ValueError(processed_reason)
    source_info["selected_layer"] = ""
    return adata, source_info


def compute_input_checksum(input_source: dict[str, Any] | None) -> str:
    """Compute a stable checksum across one or more input files."""
    if not input_source:
        return ""

    digest = hashlib.sha256()
    for path in sorted((Path(p) for p in input_source.get("files", [])), key=lambda p: str(p)):
        digest.update(path.name.encode("utf-8"))
        digest.update(b"\0")
        digest.update(sha256_file(path).encode("utf-8"))
        digest.update(b"\0")
    return digest.hexdigest()
