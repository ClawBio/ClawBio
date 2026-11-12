#!/usr/bin/env python
"""celltype-specificity-profiler — a ClawBio skill.

Given a gene, report how cell-type-specific its expression is across a
single-cell RNA-seq atlas: the tau specificity index, the expression
bimodality coefficient, and the cell types that drive the signal.

Implements two single-cell features from *The Virtual Biotech*
(Zhang et al. 2026, bioRxiv 10.64898/2026.02.23.707551):
  - tau cell-type specificity index (Yanai et al. 2005)
  - expression bimodality coefficient (skewness/kurtosis; psychometrics transfer)

Core math is factored into pure functions (compute_tau, bimodality_coefficient)
so it is importable and testable without I/O.
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
from pathlib import Path
from typing import Sequence

import numpy as np

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

SKILL_NAME = "celltype-specificity-profiler"

# tau cut for the cell-type-specific vs broadly-expressed call. Zhang et al. 2026
# derived this as the midpoint of a K-means (k=2) split of their *trial-level* tau
# distribution (tau = 0.69, Extended Methods) — a cohort-specific boundary, NOT a
# universal threshold. Used here only as an interpretive default; tau itself is the
# primary, source-independent output.
TAU_THRESHOLD = 0.69
TAU_THRESHOLD_NOTE = (
    "cohort-specific K-means(k=2) midpoint of the trial-level tau distribution in "
    "Zhang et al. 2026 (tau=0.69); an interpretive default, not a universal cutoff"
)
LOW_EXPRESSION_FRACTION = 0.01
# Zhang et al. 2026 (Extended Methods) excludes cell types with <20 cells from tau,
# since per-cell-type means from tiny populations are unreliable and, through the
# max-normalization, can distort tau.
MIN_CELLS_FOR_TAU = 20

# Trial-success priors from Zhang et al. 2026 (correlational, not causal). Odds
# ratios verified verbatim against the paper's Results: Phase I→II OR 1.27
# (95% CI 1.22–1.33); primary-endpoint OR 1.11 (95% CI 1.09–1.14).
TRIAL_PRIOR = {
    "note": "Odds ratios from Zhang et al. 2026 (bioRxiv 10.64898/2026.02.23.707551)",
    "phase_I_to_II_OR": 1.27,
    "primary_endpoint_OR": 1.11,
    "lower_AE_rate": True,
}


# --------------------------------------------------------------------------- #
# Pure functions (the algorithm). No I/O — importable for tests.
# --------------------------------------------------------------------------- #
def compute_tau(means: Sequence[float]) -> float:
    """Tau cell-type specificity index (Yanai et al. 2005).

        tau = sum_i (1 - x_i / x_max) / (n - 1)

    over n cell types, where x_i is the (non-negative) mean expression in cell
    type i. Returns a value in [0, 1]: 0 = ubiquitous, 1 = restricted to one
    cell type.
    """
    x = np.asarray(means, dtype=float)
    x = x[~np.isnan(x)]
    n = x.size
    if n < 2:
        return float("nan")
    x = np.clip(x, 0.0, None)
    x_max = x.max()
    if x_max <= 0:
        return 0.0
    return float(np.sum(1.0 - x / x_max) / (n - 1))


def bimodality_coefficient(values: Sequence[float]) -> float:
    """Sarle's bimodality coefficient over a set of values.

        BC = (skew^2 + 1) / (kurtosis + 3*(n-1)^2 / ((n-2)*(n-3)))

    where skew and kurtosis are the sample (bias-corrected, excess) statistics.
    Higher BC suggests a more bimodal ("on/off") distribution. Returns NaN when
    n < 4 (the correction term is undefined).
    """
    v = np.asarray(values, dtype=float)
    v = v[~np.isnan(v)]
    n = v.size
    if n < 4:
        return float("nan")
    mean = v.mean()
    sd = v.std(ddof=1)
    if sd == 0:
        return float("nan")
    z = (v - mean) / sd
    # bias-corrected sample skewness (G1)
    g1 = (n / ((n - 1) * (n - 2))) * np.sum(z**3)
    # bias-corrected sample excess kurtosis (G2)
    g2 = ((n * (n + 1)) / ((n - 1) * (n - 2) * (n - 3))) * np.sum(z**4) - (
        3 * (n - 1) ** 2
    ) / ((n - 2) * (n - 3))
    denom = g2 + 3 * (n - 1) ** 2 / ((n - 2) * (n - 3))
    if denom == 0:
        return float("nan")
    return float((g1**2 + 1) / denom)


# --------------------------------------------------------------------------- #
# Atlas handling
# --------------------------------------------------------------------------- #
DEMO_GENE = "MS4A1"  # real B-cell marker — cell-type-specific in PBMCs


def _load_demo():
    """Load a small **real** public dataset as a stand-in for the expression
    matrix an upstream ClawBio skill (`scrna-embedding`) hands over in the chain.

    Uses scanpy's bundled, real pbmc3k 10x dataset (2,638 cells, 8 annotated
    cell types). No synthetic data. Returns (adata, gene, atlas_name).
    """
    import scanpy as sc

    proc = sc.datasets.pbmc3k_processed()  # real, cached after first download
    # .raw holds log-normalized (non-negative) expression for all genes;
    # .X in the processed object is z-scored (has negatives) and unsuitable for tau.
    adata = proc.raw.to_adata() if proc.raw is not None else proc
    adata.obs["cell_type"] = proc.obs["louvain"].astype(str).values
    return adata, DEMO_GENE, "pbmc3k (10x, real; scanpy bundled)"


def _resolve_cell_type_key(adata, requested: str | None) -> str:
    """Find the obs column holding cell-type labels."""
    if requested and requested in adata.obs:
        return requested
    for candidate in ("cell_type", "celltype", "cell_type_name", "louvain", "leiden", "CellType"):
        if candidate in adata.obs:
            return candidate
    raise SystemExit(
        f"no cell-type column found in atlas obs ({list(adata.obs.columns)}); "
        "pass --cell-type-key"
    )


def _resolve_gene(var_names, gene: str) -> str:
    """Resolve a gene symbol against the atlas var index, with a small alias map.
    Fails loudly on a genuinely missing symbol rather than returning zeros."""
    names = list(var_names)
    if gene in names:
        return gene
    aliases = {"B7-H3": "CD276", "B7H3": "CD276", "CD276": "B7-H3"}
    alt = aliases.get(gene)
    if alt and alt in names:
        return alt
    raise SystemExit(
        f"gene {gene!r} not found in atlas var index "
        f"(checked alias {alt!r}). Available example symbols: {names[:5]}..."
    )


def _per_celltype_stats(adata, gene: str, cell_type_key: str = "cell_type"):
    """Return per-cell-type stats list, sorted by mean_expr descending."""
    import numpy as _np

    col = adata[:, gene].X
    if hasattr(col, "toarray"):
        col = col.toarray()
    col = _np.asarray(col).ravel().astype(float)

    cell_types = adata.obs[cell_type_key].astype(str).values
    stats = []
    for ct in sorted(set(cell_types)):
        mask = cell_types == ct
        vals = col[mask]
        n_cells = int(mask.sum())
        expressing = vals > 0
        pct = float(expressing.mean()) if n_cells else 0.0
        stats.append(
            {
                "cell_type": ct,
                "mean_expr": round(float(vals.mean()), 4) if n_cells else 0.0,
                "median_expr": round(float(_np.median(vals)), 4) if n_cells else 0.0,
                "pct_expressing": round(pct, 4),
                "n_cells": n_cells,
            }
        )
    stats.sort(key=lambda d: d["mean_expr"], reverse=True)
    return stats, col


def profile_gene(adata, gene: str, atlas_name: str, cell_type_key: str = "cell_type"):
    """Compute the full profile dict (without the trial_prior block)."""
    stats, col = _per_celltype_stats(adata, gene, cell_type_key)

    # tau is computed only over cell types with >= MIN_CELLS_FOR_TAU cells
    # (Zhang et al. 2026). Smaller cell types are still reported in
    # per_celltype_stats for transparency, just excluded from the tau means.
    tau_stats = [s for s in stats if s["n_cells"] >= MIN_CELLS_FOR_TAU]
    means = [s["mean_expr"] for s in tau_stats]
    tau = compute_tau(means)
    n_excluded = len(stats) - len(tau_stats)

    expressing_vals = col[col > 0]
    bc = bimodality_coefficient(expressing_vals)

    frac_expressing = float((col > 0).mean()) if col.size else 0.0
    low_expression = frac_expressing < LOW_EXPRESSION_FRACTION

    interpretation = (
        f"cell-type-specific (tau > {TAU_THRESHOLD})"
        if (not np.isnan(tau) and tau > TAU_THRESHOLD)
        else "broadly expressed"
    )

    top = [
        {
            "cell_type": s["cell_type"],
            "mean_expr": s["mean_expr"],
            "pct_expressing": s["pct_expressing"],
        }
        for s in stats[:3]
    ]

    return {
        "skill": SKILL_NAME,
        "gene": gene,
        "atlas": atlas_name,
        "tau": round(tau, 4) if not np.isnan(tau) else None,
        "tau_threshold": TAU_THRESHOLD,
        "tau_threshold_note": TAU_THRESHOLD_NOTE,
        "n_cell_types_used_for_tau": len(tau_stats),
        "n_cell_types_excluded_small": n_excluded,
        "bimodality_coefficient": round(bc, 4) if not np.isnan(bc) else None,
        "interpretation": interpretation,
        "low_expression": bool(low_expression),
        "top_cell_types": top,
        "per_celltype_stats": stats,
    }


# --------------------------------------------------------------------------- #
# Reproducibility bundle
# --------------------------------------------------------------------------- #
def _write_reproducibility(out_dir: str, argv: Sequence[str], output_files):
    write_commands_sh(
        out_dir, shlex.join(["python", f"skills/{SKILL_NAME}/profiler.py", *argv])
    )
    env_path = write_environment_yml(
        out_dir,
        f"clawbio-{SKILL_NAME}",
        pip_deps=["scanpy", "anndata", "numpy>=1.23", "scipy>=1.9", "pandas>=2.0"],
    )
    write_checksums([*output_files, env_path], out_dir, anchor=out_dir)


def _write_csv(path: str, stats):
    import csv

    fields = ["cell_type", "mean_expr", "median_expr", "pct_expressing", "n_cells"]
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fields)
        writer.writeheader()
        for row in stats:
            writer.writerow(row)


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #
def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="profiler.py",
        description="Per-gene cell-type specificity profiler (tau + bimodality).",
    )
    p.add_argument("--gene", type=str, default=None, help="HGNC symbol, e.g. MS4A1")
    p.add_argument(
        "--atlas",
        type=str,
        default=None,
        help="Path to an .h5ad expression matrix (in the chain: upstream scrna-embedding output)",
    )
    p.add_argument(
        "--cell-type-key",
        type=str,
        default="cell_type",
        help="obs column holding cell-type labels (default: cell_type)",
    )
    p.add_argument("--tissue", type=str, default=None, help="Restrict to a tissue label")
    p.add_argument(
        "--trial-prior",
        action="store_true",
        help="Attach Zhang et al. 2026 trial-success odds ratios",
    )
    p.add_argument(
        "--demo",
        action="store_true",
        help="Run on a small REAL public dataset (scanpy pbmc3k) — no synthetic data",
    )
    p.add_argument(
        "--output",
        "--out",
        dest="out",
        type=str,
        default="./output",
        help="Output directory (--output is the ClawBio runner convention; --out is accepted too)",
    )
    return p


def main(argv: Sequence[str] | None = None) -> int:
    argv = list(sys.argv[1:] if argv is None else argv)
    args = build_parser().parse_args(argv)

    import anndata as ad

    if args.demo:
        adata, demo_gene, atlas_name = _load_demo()
        gene = args.gene or demo_gene
    else:
        if not args.atlas:
            raise SystemExit("--atlas is required unless --demo is given")
        if not args.gene:
            raise SystemExit("--gene is required unless --demo is given")
        if not os.path.exists(args.atlas):
            raise SystemExit(f"atlas not found: {args.atlas}")
        adata = ad.read_h5ad(args.atlas)
        atlas_name = str(adata.uns.get("atlas_name", os.path.basename(args.atlas)))
        gene = args.gene

    cell_type_key = _resolve_cell_type_key(adata, args.cell_type_key)

    if args.tissue:
        if "tissue" not in adata.obs:
            raise SystemExit("atlas has no 'tissue' column but --tissue was given")
        adata = adata[adata.obs["tissue"].astype(str) == args.tissue].copy()
        if adata.n_obs == 0:
            raise SystemExit(f"no cells with tissue == {args.tissue!r}")

    resolved = _resolve_gene(adata.var_names, gene)
    profile = profile_gene(adata, resolved, atlas_name, cell_type_key)
    # report under the symbol the user asked for
    profile["gene"] = gene if resolved == gene else resolved

    if args.trial_prior:
        tau = profile["tau"]
        specific = tau is not None and tau > TAU_THRESHOLD
        block = dict(TRIAL_PRIOR)
        block["lower_AE_rate"] = bool(specific)
        profile["trial_prior"] = block

    out_dir = os.path.abspath(args.out)
    os.makedirs(out_dir, exist_ok=True)
    profile_path = os.path.join(out_dir, "profile.json")
    csv_path = os.path.join(out_dir, "per_celltype.csv")

    with open(profile_path, "w") as fh:
        json.dump(profile, fh, indent=2)
    _write_csv(csv_path, profile["per_celltype_stats"])
    _write_reproducibility(
        out_dir,
        argv,
        [profile_path, csv_path],
    )

    print(json.dumps(profile, indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
