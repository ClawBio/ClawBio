#!/usr/bin/env python3
"""
finemapping_benchmark.py — Swappable fine-mapping pipeline benchmark.

Runs multiple fine-mapping methods (ABF, SuSiE, and future methods) on
the same synthetic locus data, scores each against injected causal signals,
and picks the best. This is LOI Milestone 1: the first swappable pipeline.

The pattern: agent swaps methods, runs cascade, benchmarks, keeps best.

Usage:
    python tests/benchmark/finemapping_benchmark.py
    python tests/benchmark/finemapping_benchmark.py --output /tmp/fm_bench
    python tests/benchmark/finemapping_benchmark.py --methods abf,susie
"""

from __future__ import annotations

import json
import sys
import time
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# Add project root and fine-mapping skill to path
PROJECT_ROOT = Path(__file__).resolve().parents[2]
SKILL_DIR = PROJECT_ROOT / "skills" / "fine-mapping"
sys.path.insert(0, str(PROJECT_ROOT))
sys.path.insert(0, str(SKILL_DIR))

from fine_mapping_core.abf import compute_abf
from fine_mapping_core.susie import run_susie
from fine_mapping_core.susie_inf import run_susie_inf, cred_inf
from fine_mapping_core.credible_sets import build_credible_set_abf, build_credible_sets_susie

# ---------------------------------------------------------------------------
# Synthetic benchmark locus generator
# ---------------------------------------------------------------------------

# Known causal variant indices (ground truth for this benchmark)
CAUSAL_INDICES = [60, 140]
CAUSAL_EFFECTS = [0.25, 0.20]


def make_benchmark_locus(
    n_variants: int = 200,
    n_samples: int = 5000,
    causal_indices: list[int] | None = None,
    causal_effects: list[float] | None = None,
    seed: int = 42,
) -> tuple[pd.DataFrame, np.ndarray, list[int]]:
    """Generate a synthetic GWAS locus with known causal signals.

    Returns (sumstats_df, ld_matrix, causal_indices).
    """
    rng = np.random.default_rng(seed)
    causal_idx = causal_indices or CAUSAL_INDICES
    effects = causal_effects or CAUSAL_EFFECTS

    # LD structure: two blocks with AR(1) correlation
    rho = 0.8
    block_size = n_variants // 2
    R = np.eye(n_variants)

    for i in range(n_variants):
        for j in range(i + 1, n_variants):
            same_block = (i < block_size) == (j < block_size)
            dist = abs(i - j)
            if same_block:
                r = rho ** dist
            else:
                r = 0.05 * rho ** dist
            R[i, j] = r
            R[j, i] = r

    # Simulate genotypes via Cholesky
    L = np.linalg.cholesky(R + 1e-6 * np.eye(n_variants))
    Z_raw = rng.standard_normal((n_samples, n_variants))
    X = Z_raw @ L.T
    # Binarize to approximate genotypes
    X = (X > 0).astype(float)

    # Simulate phenotype with causal effects
    y = np.zeros(n_samples)
    for idx, eff in zip(causal_idx, effects):
        y += eff * X[:, idx]
    y += rng.standard_normal(n_samples) * 0.5

    # Compute summary statistics
    betas = np.zeros(n_variants)
    ses = np.zeros(n_variants)
    pvals = np.zeros(n_variants)

    for i in range(n_variants):
        xi = X[:, i]
        xi_centered = xi - xi.mean()
        denom = (xi_centered ** 2).sum()
        if denom == 0:
            betas[i] = 0
            ses[i] = 1
            pvals[i] = 1
            continue
        b = (xi_centered * (y - y.mean())).sum() / denom
        resid = y - y.mean() - b * xi_centered
        se = np.sqrt((resid ** 2).sum() / ((n_samples - 2) * denom))
        z = b / se if se > 0 else 0
        from scipy import stats
        p = 2 * stats.norm.sf(abs(z))
        betas[i] = b
        ses[i] = se
        pvals[i] = p

    df = pd.DataFrame({
        "rsid": [f"rs_bench_{i:04d}" for i in range(n_variants)],
        "chr": "2",
        "pos": [127890000 + i * 500 for i in range(n_variants)],
        "ref": "G",
        "alt": "A",
        "beta": betas,
        "se": ses,
        "pvalue": pvals,
        "z": betas / np.where(ses > 0, ses, 1),
        "n": n_samples,
    })

    return df, R, causal_idx


# ---------------------------------------------------------------------------
# Method runners
# ---------------------------------------------------------------------------


def run_method_abf(
    sumstats: pd.DataFrame, ld: np.ndarray, **kwargs
) -> dict:
    """Run ABF fine-mapping and return results."""
    t0 = time.time()
    pips = compute_abf(sumstats, w=0.04)
    credsets = build_credible_set_abf(pips, sumstats, coverage=0.95)
    elapsed = time.time() - t0

    # Extract indices from credible set variants
    rsid_to_idx = {r: i for i, r in enumerate(sumstats["rsid"])}
    cs_indices = []
    for cs in credsets:
        for v in cs.get("variants", []):
            idx = rsid_to_idx.get(v.get("rsid"))
            if idx is not None:
                cs_indices.append(idx)
    cs_indices = sorted(set(cs_indices))

    # Total coverage from credible sets
    coverage = credsets[0].get("coverage", 0.95) if credsets else 0.0

    return {
        "method": "ABF",
        "pips": pips.tolist(),
        "credible_set_indices": cs_indices,
        "credible_set_size": len(cs_indices),
        "coverage": coverage,
        "elapsed": round(elapsed, 3),
    }


def run_method_susie(
    sumstats: pd.DataFrame, ld: np.ndarray, **kwargs
) -> dict:
    """Run SuSiE fine-mapping and return results."""
    t0 = time.time()
    z = sumstats["z"].values
    n = int(sumstats["n"].iloc[0])

    susie_result = run_susie(z=z, R=ld, n=n, L=10, max_iter=100, tol=1e-3)
    credsets = build_credible_sets_susie(
        susie_result["alpha"], sumstats, R=ld, coverage=0.95, min_purity=0.5
    )
    elapsed = time.time() - t0

    pips = susie_result["pip"].tolist()
    rsid_to_idx = {r: i for i, r in enumerate(sumstats["rsid"])}
    cs_indices = []
    for cs in credsets:
        for v in cs.get("variants", []):
            idx = rsid_to_idx.get(v.get("rsid"))
            if idx is not None:
                cs_indices.append(idx)
    cs_indices = sorted(set(cs_indices))

    # Extract indices per credible set for reporting
    cs_details = []
    for cs in credsets:
        indices = [rsid_to_idx[v["rsid"]] for v in cs.get("variants", []) if v.get("rsid") in rsid_to_idx]
        cs_details.append({"indices": indices, "coverage": cs.get("coverage", 0),
                           "purity": cs.get("purity", 0)})

    return {
        "method": "SuSiE",
        "pips": pips,
        "credible_sets": cs_details,
        "n_credible_sets": len(credsets),
        "credible_set_indices": cs_indices,
        "credible_set_size": len(cs_indices),
        "elapsed": round(elapsed, 3),
    }


def run_method_susieinf(
    sumstats: pd.DataFrame, ld: np.ndarray, **kwargs
) -> dict:
    """Run SuSiE-inf fine-mapping and return results."""
    t0 = time.time()
    z = sumstats["z"].values
    n = int(sumstats["n"].iloc[0])

    result = run_susie_inf(z=z, R=ld, n=n, L=10, max_iter=100, tol=1e-3)
    credsets = cred_inf(result["alpha"], R=ld, coverage=0.95, purity=0.5)
    elapsed = time.time() - t0

    rsid_to_idx = {r: i for i, r in enumerate(sumstats["rsid"])}
    cs_indices = sorted({
        i
        for cs in credsets
        for i in cs
    })

    cs_details = [
        {"indices": cs, "size": len(cs)}
        for cs in credsets
    ]

    return {
        "method": "SuSiE-inf",
        "pips": result["pip"].tolist(),
        "credible_sets": cs_details,
        "n_credible_sets": len(credsets),
        "credible_set_indices": cs_indices,
        "credible_set_size": len(cs_indices),
        "elapsed": round(elapsed, 3),
    }


# Registry of available methods
METHODS = {
    "abf": run_method_abf,
    "susie": run_method_susie,
    "susieinf": run_method_susieinf,
    # Future methods: "finemap", "polyfun" will be added here
}


# ---------------------------------------------------------------------------
# Scoring
# ---------------------------------------------------------------------------


def score_method(result: dict, causal_indices: list[int], n_variants: int) -> dict:
    """Score a method's results against known causal signals.

    Metrics:
    - causal_in_credset: how many true causal variants are in credible sets
    - causal_pip_sum: sum of PIPs at causal positions (higher = better)
    - credset_size: smaller is better (more precise)
    - precision: causal_in_credset / credset_size
    - recall: causal_in_credset / n_causal
    - pip_rank: average rank of causal variants by PIP (lower = better)
    """
    pips = np.array(result["pips"])
    cs_indices = set(result.get("credible_set_indices", []))
    n_causal = len(causal_indices)

    # How many causal variants captured in credible set
    causal_captured = sum(1 for c in causal_indices if c in cs_indices)

    # PIP at causal positions
    causal_pips = [pips[c] for c in causal_indices]
    causal_pip_sum = sum(causal_pips)

    # Rank of causal variants (1-indexed, lower = better)
    pip_order = np.argsort(-pips)
    ranks = np.empty_like(pip_order)
    ranks[pip_order] = np.arange(1, len(pips) + 1)
    causal_ranks = [int(ranks[c]) for c in causal_indices]
    avg_rank = sum(causal_ranks) / len(causal_ranks) if causal_ranks else n_variants

    cs_size = result.get("credible_set_size", 0)
    precision = causal_captured / cs_size if cs_size > 0 else 0.0
    recall = causal_captured / n_causal if n_causal > 0 else 0.0
    f1 = (2 * precision * recall / (precision + recall)) if (precision + recall) > 0 else 0.0

    # Composite score: weighted combination
    # Higher is better. Reward recall and PIP concentration, penalise large credible sets.
    composite = (
        0.4 * recall
        + 0.3 * (causal_pip_sum / n_causal)
        + 0.2 * precision
        + 0.1 * (1.0 - avg_rank / n_variants)
    )

    return {
        "method": result["method"],
        "causal_captured": causal_captured,
        "causal_total": n_causal,
        "causal_pips": [round(p, 4) for p in causal_pips],
        "causal_pip_sum": round(causal_pip_sum, 4),
        "causal_ranks": causal_ranks,
        "avg_causal_rank": round(avg_rank, 1),
        "credible_set_size": cs_size,
        "precision": round(precision, 4),
        "recall": round(recall, 4),
        "f1": round(f1, 4),
        "composite_score": round(composite, 4),
        "elapsed": result.get("elapsed", 0),
    }


# ---------------------------------------------------------------------------
# Benchmark runner
# ---------------------------------------------------------------------------


def run_benchmark(
    methods: list[str] | None = None,
    seed: int = 42,
    output_dir: Path | None = None,
) -> dict:
    """Run all specified methods on the same synthetic locus and compare.

    Returns a dict with per-method scores and the winning method.
    """
    method_names = methods or list(METHODS.keys())
    print(f"Generating benchmark locus (seed={seed})...")
    sumstats, ld, causal_idx = make_benchmark_locus(seed=seed)
    n_variants = len(sumstats)

    print(f"Causal variants at indices {causal_idx}")
    print(f"Running {len(method_names)} methods: {', '.join(method_names)}")
    print()

    results = []
    for name in method_names:
        if name not in METHODS:
            print(f"  Unknown method: {name}, skipping")
            continue
        print(f"  [{name}] running... ", end="", flush=True)
        try:
            raw = METHODS[name](sumstats, ld)
            scored = score_method(raw, causal_idx, n_variants)
            print(f"done ({scored['elapsed']}s) "
                  f"recall={scored['recall']:.2f} "
                  f"precision={scored['precision']:.2f} "
                  f"composite={scored['composite_score']:.4f}")
            results.append(scored)
        except Exception as e:
            print(f"FAILED: {e}")
            results.append({
                "method": name,
                "error": str(e),
                "composite_score": 0,
            })

    # Pick winner
    valid = [r for r in results if "error" not in r]
    if valid:
        winner = max(valid, key=lambda r: r["composite_score"])
    else:
        winner = None

    benchmark = {
        "seed": seed,
        "n_variants": n_variants,
        "causal_indices": causal_idx,
        "methods": results,
        "winner": winner["method"] if winner else None,
        "winner_score": winner["composite_score"] if winner else 0,
    }

    # Summary
    print()
    print("=" * 60)
    print("BENCHMARK RESULTS")
    print("=" * 60)
    print(f"{'Method':<12} {'Recall':<8} {'Prec':<8} {'F1':<8} {'CS Size':<10} {'Composite':<10}")
    print("-" * 60)
    for r in sorted(results, key=lambda x: -x.get("composite_score", 0)):
        if "error" in r:
            print(f"{r['method']:<12} ERROR: {r['error'][:40]}")
        else:
            print(f"{r['method']:<12} {r['recall']:<8.4f} {r['precision']:<8.4f} "
                  f"{r['f1']:<8.4f} {r['credible_set_size']:<10} {r['composite_score']:<10.4f}")
    print("-" * 60)
    if winner:
        print(f"WINNER: {winner['method']} (composite={winner['composite_score']:.4f})")
    print()

    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "finemapping_benchmark.json").write_text(
            json.dumps(benchmark, indent=2)
        )
        print(f"Results written to {output_dir}/finemapping_benchmark.json")

    return benchmark


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    import argparse

    parser = argparse.ArgumentParser(description="ClawBio fine-mapping benchmark")
    parser.add_argument("--methods", type=str, default=None,
                        help="Comma-separated methods (default: all)")
    parser.add_argument("--seed", type=int, default=42)
    parser.add_argument("--output", type=str, default=None)
    args = parser.parse_args()

    methods = args.methods.split(",") if args.methods else None
    output_dir = Path(args.output) if args.output else None

    benchmark = run_benchmark(methods=methods, seed=args.seed, output_dir=output_dir)
    sys.exit(0 if benchmark.get("winner") else 1)


if __name__ == "__main__":
    main()
