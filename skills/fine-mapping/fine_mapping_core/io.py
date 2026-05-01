"""
io.py — Input parsing for fine-mapping.

Handles sumstats TSV/CSV/TXT and LD matrix loading (.npy or .tsv).
"""

from __future__ import annotations

import numpy as np
import pandas as pd
from pathlib import Path


# ---------------------------------------------------------------------------
# Sumstats
# ---------------------------------------------------------------------------

REQUIRED_COLS_Z = {"rsid", "z"}
REQUIRED_COLS_BETA = {"rsid", "beta", "se"}
OPTIONAL_COLS = {"chr", "pos", "p", "maf", "n", "a1", "a2"}


def load_sumstats(path: Path, chr_: str | None = None, start: int | None = None, end: int | None = None) -> pd.DataFrame:
    """Load GWAS summary statistics from a delimited file.

    Accepts TSV, CSV, or whitespace-delimited text. Detects whether the file
    provides z-scores directly (column ``z``) or beta + se (both required).
    Computes z = beta / se when only beta+se are present.

    After loading, optionally filters to a genomic window if chr/start/end
    are supplied.

    Returns a DataFrame with columns:
        rsid, chr, pos, z, beta (if available), se (if available),
        p (if available), maf (if available), n (if available)
    """
    raw = _read_delimited(path)
    raw.columns = [c.strip().lower() for c in raw.columns]

    # Rename common aliases
    _alias(raw, {"snp": "rsid", "variant_id": "rsid", "id": "rsid",
                 "effect_allele": "a1", "other_allele": "a2",
                 "effect_size": "beta", "standard_error": "se",
                 "zscore": "z", "z_score": "z",
                 "chromosome": "chr", "position": "pos", "bp": "pos",
                 "pval": "p", "p_value": "p", "pvalue": "p", "p.value": "p"})

    has_z = "z" in raw.columns
    has_beta_se = ("beta" in raw.columns and "se" in raw.columns)

    if not has_z and not has_beta_se:
        raise ValueError(
            "Sumstats must contain either a 'z' column or both 'beta' and 'se' columns. "
            f"Found columns: {list(raw.columns)}"
        )

    if not has_z:
        # Validate SE before computing z = beta / se
        invalid_se = (raw["se"] <= 0) | raw["se"].isna()
        if invalid_se.any():
            n_bad = int(invalid_se.sum())
            print(f"  [io] WARNING: {n_bad} rows have zero, negative, or missing SE. "
                  "These will produce invalid z-scores and be dropped.")
        raw["z"] = raw["beta"] / raw["se"]

    if "rsid" not in raw.columns:
        raise ValueError("Sumstats must contain an 'rsid' (or 'snp'/'variant_id') column.")

    # Coerce numeric columns
    for col in ["z", "beta", "se", "pos", "p", "maf", "n"]:
        if col in raw.columns:
            raw[col] = pd.to_numeric(raw[col], errors="coerce")

    # Drop rows with missing or non-finite z
    n_before = len(raw)
    raw = raw[raw["z"].notna() & np.isfinite(raw["z"])]
    if len(raw) < n_before:
        print(f"  [io] Dropped {n_before - len(raw)} rows with missing/infinite z-scores")

    # Locus window filter
    if chr_ is not None and start is not None and end is not None:
        if "chr" not in raw.columns or "pos" not in raw.columns:
            raise ValueError("Locus window filter requires 'chr' and 'pos' columns in sumstats.")
        chr_clean = str(chr_).lstrip("chr")
        raw["_chr_clean"] = raw["chr"].astype(str).str.lstrip("chr")
        raw = raw[(raw["_chr_clean"] == chr_clean) & (raw["pos"] >= start) & (raw["pos"] <= end)]
        raw = raw.drop(columns=["_chr_clean"])
        print(f"  [io] Locus filter chr{chr_clean}:{start:,}-{end:,} → {len(raw)} variants")

    raw = raw.reset_index(drop=True)
    return raw


def _read_delimited(path: Path) -> pd.DataFrame:
    """Try TSV then CSV then whitespace-delimited."""
    path = Path(path)
    for sep in ["\t", ",", r"\s+"]:
        try:
            df = pd.read_csv(path, sep=sep, engine="python", comment="#")
            if df.shape[1] > 1:
                return df
        except Exception:
            continue
    raise ValueError(f"Could not parse {path} as TSV, CSV, or whitespace-delimited text.")


def _alias(df: pd.DataFrame, mapping: dict[str, str]) -> None:
    """Rename columns in-place using alias mapping (only if target not already present)."""
    for old, new in mapping.items():
        if old in df.columns and new not in df.columns:
            df.rename(columns={old: new}, inplace=True)


# ---------------------------------------------------------------------------
# LD matrix
# ---------------------------------------------------------------------------


def load_ld(path: Path, n_variants: int) -> np.ndarray:
    """Load a pre-computed LD (correlation) matrix.

    Accepts:
      - ``.npy`` files: loaded directly via ``np.load``
      - ``.tsv`` / ``.csv`` / ``.txt`` files: loaded via pandas

    The matrix must be square with dimension == n_variants.
    """
    path = Path(path)
    if path.suffix == ".npy":
        R = np.load(path)
    else:
        # Detect whether the file has a string header row (e.g. rsIDs) and skip it
        peek = pd.read_csv(path, sep=None, header=None, engine="python", nrows=1)
        first_val = str(peek.iloc[0, 0])
        try:
            float(first_val)
            hdr = None
        except ValueError:
            hdr = 0  # first row is a header — skip it
        R = pd.read_csv(path, sep=None, header=hdr, index_col=None, engine="python").values

    R = np.asarray(R, dtype=float)

    if R.ndim != 2 or R.shape[0] != R.shape[1]:
        raise ValueError(f"LD matrix must be square; got shape {R.shape}")
    if R.shape[0] != n_variants:
        raise ValueError(
            f"LD matrix dimension ({R.shape[0]}) != number of variants in sumstats ({n_variants}). "
            "Ensure the LD matrix rows/columns correspond to the variants in the same order."
        )

    # Clip to [-1, 1] to guard against floating-point drift
    R = np.clip(R, -1.0, 1.0)
    # Ensure exact symmetry
    R = (R + R.T) / 2.0
    np.fill_diagonal(R, 1.0)
    return R
