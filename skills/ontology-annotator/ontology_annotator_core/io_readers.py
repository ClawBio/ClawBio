"""io_readers.py — Load a metadata table (CSV, TSV, or .h5ad `obs`) for
ontology annotation.
"""

from __future__ import annotations

from pathlib import Path

import pandas as pd

DEFAULT_COLUMN_ONTOLOGY_MAP = {
    "tissue": "uberon",
    "cell_type": "cl",
    "cell type": "cl",
    "celltype": "cl",
    "disease": "mondo",
    "trait": "efo",
    "phenotype": "efo",
}


def load_table(input_path: Path) -> pd.DataFrame:
    """Load a metadata table from `.csv`, `.tsv`, or `.h5ad`."""
    suffix = input_path.suffix.lower()
    if suffix == ".h5ad":
        return load_h5ad_obs(input_path)
    if suffix == ".tsv":
        return pd.read_csv(input_path, sep="\t")
    if suffix == ".csv":
        return pd.read_csv(input_path)
    # Unknown extension: sniff comma first, then tab.
    try:
        return pd.read_csv(input_path)
    except Exception:
        return pd.read_csv(input_path, sep="\t")


def load_h5ad_obs(input_path: Path) -> pd.DataFrame:
    """Load the `.obs` metadata table out of an AnnData `.h5ad` file.

    `anndata` is an optional dependency (only .h5ad input needs it) — raises
    a RuntimeError with an install hint rather than a raw ImportError.
    """
    try:
        import anndata as ad  # type: ignore
    except ImportError as exc:
        raise RuntimeError(
            "anndata is required to read .h5ad input. Install it with: pip install anndata"
        ) from exc

    adata = ad.read_h5ad(input_path)
    obs = adata.obs.copy()
    obs.insert(0, "cell_id", obs.index.astype(str))
    return obs.reset_index(drop=True)


def parse_columns_arg(columns_arg: str) -> dict[str, str]:
    """Parse `--columns col:ontology,col2:ontology2` into `{col: ontology}`."""
    parsed: dict[str, str] = {}
    for chunk in columns_arg.split(","):
        chunk = chunk.strip()
        if not chunk:
            continue
        if ":" not in chunk:
            raise ValueError(
                f"--columns entry {chunk!r} is missing ':<ontology>' (expected col:ontology)"
            )
        col, ontology = chunk.split(":", 1)
        col, ontology = col.strip(), ontology.strip().lower()
        if not col or not ontology:
            raise ValueError(f"--columns entry {chunk!r} has an empty column or ontology")
        parsed[col] = ontology
    return parsed


def infer_columns(df: pd.DataFrame) -> dict[str, str]:
    """Auto-detect columns to annotate when `--columns` was not given.

    Matches df column names (case-insensitively) against
    DEFAULT_COLUMN_ONTOLOGY_MAP. Returns {} if nothing matches, so the
    caller can raise a clear "specify --columns" error rather than silently
    annotating the wrong thing.
    """
    inferred: dict[str, str] = {}
    lower_lookup = {c.lower(): c for c in df.columns}
    for candidate_name, ontology in DEFAULT_COLUMN_ONTOLOGY_MAP.items():
        actual = lower_lookup.get(candidate_name)
        if actual and actual not in inferred:
            inferred[actual] = ontology
    return inferred
