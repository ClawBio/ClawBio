"""
abf.py — Approximate Bayes Factors for single-causal-variant fine-mapping.

Reference: Wakefield (2009) Am J Hum Genet, doi:10.1016/j.ajhg.2008.12.010
"""

from __future__ import annotations

import numpy as np
import pandas as pd


DEFAULT_W = 0.04  # prior variance on effect size (σ = 0.2 on log-OR scale)


def compute_abf(df: pd.DataFrame, w: float = DEFAULT_W) -> np.ndarray:
    """Compute Approximate Bayes Factors and posterior inclusion probabilities.

    Parameters
    ----------
    df:
        DataFrame with column ``z`` (z-scores). Optionally ``se`` (standard
        error) and ``n`` (sample size). If ``se`` is present, V_i = se_i^2.
        Otherwise V_i = 1/n_eff where n_eff is the median sample size or 10000.
    w:
        Prior variance W (Wakefield 2009 default 0.04).

    Returns
    -------
    pip : np.ndarray, shape (n_variants,)
        Posterior inclusion probabilities assuming at most one causal variant.
    """
    z = df["z"].values.astype(float)
    n_variants = len(z)

    if "se" in df.columns and df["se"].notna().all():
        se_vals = df["se"].values.astype(float)
        if np.any(se_vals <= 0):
            raise ValueError(
                "Standard errors must be positive. Got zero or negative values: "
                "se=0 produces V=0 which makes log ABF undefined."
            )
        V = se_vals ** 2
    elif "n" in df.columns and df["n"].notna().any():
        n_vals = df["n"].values.astype(float)
        if np.any(n_vals <= 0):
            raise ValueError(
                "Sample sizes must be positive. Got zero or negative values."
            )
        # Per-variant V_i = 1/n_i (not median-collapsed)
        V = 1.0 / n_vals
    else:
        # Default: assume n_eff = 10 000
        V = np.full(n_variants, 1.0 / 10_000.0)

    log_abf = _log_abf(z, V, w)

    # Normalise to PIPs (uniform prior across variants)
    log_abf_shifted = log_abf - np.max(log_abf)
    abf = np.exp(log_abf_shifted)
    pip = abf / abf.sum()
    return pip


def _log_abf(z: np.ndarray, V: np.ndarray, w: float) -> np.ndarray:
    """Compute log ABF for each variant (Wakefield 2009 eq. 3).

    log ABF_i = 0.5 * log(V_i / (V_i + W))
                + 0.5 * z_i^2 * W / (V_i + W)
    """
    r = w / (V + w)
    log_abf = 0.5 * np.log(1.0 - r) + 0.5 * z**2 * r
    return log_abf
