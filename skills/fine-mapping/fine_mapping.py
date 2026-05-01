#!/usr/bin/env python3
"""
SuSiE Fine-Mapper — Statistical fine-mapping of GWAS loci.

Implements Approximate Bayes Factors (ABF, Wakefield 2009) and SuSiE
(Wang et al. 2020) in pure Python/numpy. No R or external SuSiE package
required.

Usage:
    # ABF (no LD needed)
    python fine_mapping.py --sumstats locus.tsv --output /tmp/finemapping

    # SuSiE with LD matrix
    python fine_mapping.py --sumstats locus.tsv --ld ld_matrix.npy --output /tmp/finemapping

    # Locus window filter
    python fine_mapping.py --sumstats gwas_full.tsv \\
        --chr 1 --start 109000000 --end 110000000 \\
        --ld ld_matrix.npy --output /tmp/finemapping

    # Demo (synthetic 200-variant locus, two causal signals)
    python fine_mapping.py --demo --output /tmp/finemapping_demo
"""

from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Path setup
# ---------------------------------------------------------------------------

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from fine_mapping_core.io import load_sumstats, load_ld
from fine_mapping_core.abf import compute_abf
from fine_mapping_core.susie import run_susie
from fine_mapping_core.credible_sets import build_credible_sets_susie, build_credible_set_abf
from fine_mapping_core.report import (
    generate_markdown, write_tables, generate_figures,
    write_reproducibility, DISCLAIMER,
)

try:
    from clawbio.common.report import write_result_json
    _HAS_COMMON = True
except ImportError:
    _HAS_COMMON = False


# ---------------------------------------------------------------------------
# Demo data generator
# ---------------------------------------------------------------------------


def make_demo_data(seed: int = 42) -> tuple[pd.DataFrame, np.ndarray]:
    """Generate a synthetic 200-variant locus with two injected causal signals.

    Signal 1: variant index 60 (effect size ~0.25)
    Signal 2: variant index 140 (effect size ~0.20)

    LD is simulated as an AR(1) structure with ρ = 0.8 in two LD blocks,
    with weaker inter-block correlation.
    """
    rng = np.random.default_rng(seed)
    n_variants = 200
    n_samples = 5000

    # Positions along chromosome 1 (simulated region 109.0–109.2 Mb)
    positions = np.linspace(109_000_000, 109_200_000, n_variants, dtype=int)
    rsids = [f"rs_demo_{i+1:04d}" for i in range(n_variants)]

    # LD matrix: two LD blocks with AR(1) structure
    R = _make_block_ld(n_variants, block_size=100, rho=0.8, inter_rho=0.05, rng=rng)

    # True causal variants
    causal = {60: 0.25, 140: 0.20}

    # Simulate z-scores: z_i = sqrt(n) * sum_j R_ij * beta_j + noise
    true_betas = np.zeros(n_variants)
    for idx, beta in causal.items():
        true_betas[idx] = beta

    signal = np.sqrt(n_samples) * (R @ true_betas)
    noise = rng.normal(0, 1, n_variants)
    z = signal + noise

    # Compute p-values from z using normal CDF approximation
    # Abramowitz & Stegun rational approximation for erfc
    def _norm_sf(x):
        """Survival function of standard normal via erfc."""
        return 0.5 * _erfc(x / np.sqrt(2))

    def _erfc(x):
        """Complementary error function approximation (A&S 7.1.26)."""
        t = 1.0 / (1.0 + 0.3275911 * np.abs(x))
        poly = t * (0.254829592 + t * (-0.284496736 + t * (1.421413741
               + t * (-1.453152027 + t * 1.061405429))))
        return poly * np.exp(-(x**2))

    p = 2 * _norm_sf(np.abs(z))

    df = pd.DataFrame({
        "rsid": rsids,
        "chr": "1",
        "pos": positions,
        "z": z,
        "se": np.full(n_variants, 1.0 / np.sqrt(n_samples)),
        "n": n_samples,
        "p": p,
    })

    return df, R


def _make_block_ld(n: int, block_size: int, rho: float, inter_rho: float, rng) -> np.ndarray:
    """Block AR(1) LD matrix."""
    R = np.full((n, n), inter_rho)
    for start in range(0, n, block_size):
        end = min(start + block_size, n)
        size = end - start
        for i in range(size):
            for j in range(size):
                R[start + i, start + j] = rho ** abs(i - j)
    np.fill_diagonal(R, 1.0)
    # Ensure PSD via nearest correlation
    R = _nearPD(R)
    return R


def _nearPD(A: np.ndarray, eps: float = 1e-6) -> np.ndarray:
    """Project matrix to nearest positive semi-definite correlation matrix."""
    eigvals, eigvecs = np.linalg.eigh(A)
    eigvals = np.maximum(eigvals, eps)
    A_psd = eigvecs @ np.diag(eigvals) @ eigvecs.T
    # Re-scale to correlation matrix
    d = np.sqrt(np.diag(A_psd))
    A_psd = A_psd / np.outer(d, d)
    np.fill_diagonal(A_psd, 1.0)
    return A_psd


# ---------------------------------------------------------------------------
# Core pipeline
# ---------------------------------------------------------------------------


def run_finemapping(
    sumstats_path: Optional[Path],
    ld_path: Optional[Path],
    output_dir: Path,
    chr_: Optional[str] = None,
    start: Optional[int] = None,
    end: Optional[int] = None,
    max_signals: int = 10,
    coverage: float = 0.95,
    min_purity: float = 0.5,
    w: float = 0.04,
    make_figures: bool = True,
    gene_track: bool = False,
    demo: bool = False,
) -> dict:
    """Run the full fine-mapping pipeline.

    Returns a results dict with keys: method, pip, credible_sets, params.
    """
    if coverage <= 0 or coverage > 1:
        raise ValueError("Coverage must be in (0, 1], got %s" % coverage)
    if min_purity < 0 or min_purity > 1:
        raise ValueError("min_purity must be in [0, 1], got %s" % min_purity)

    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Load data ---
    if demo:
        print("  Generating synthetic demo locus (200 variants, 2 causal signals)...")
        df, R = make_demo_data()
        input_label = "synthetic_demo"
    else:
        print(f"  Loading sumstats: {sumstats_path}")
        df = load_sumstats(sumstats_path, chr_=chr_, start=start, end=end)
        input_label = str(sumstats_path)

        R = None
        if ld_path:
            print(f"  Loading LD matrix: {ld_path}")
            R = load_ld(ld_path, len(df))

    if len(df) == 0:
        raise ValueError("No variants remain after filtering. Check locus window parameters.")

    print(f"  Variants: {len(df):,}")

    # --- Step 2: Fine-mapping ---
    if R is not None:
        print(f"  Method: SuSiE (L={max_signals}, coverage={coverage:.0%})")
        n_eff = int(df["n"].median()) if "n" in df.columns and df["n"].notna().any() else 10_000
        # null_weight: prior probability of no effect for each single-effect
        # regression. Prevents phantom PIPs on null loci. Value 1/(L+1) gives
        # equal prior to "no effect" as to each of L possible effects.
        null_wt = 1.0 / (max_signals + 1)
        result = run_susie(
            z=df["z"].values,
            R=R,
            n=n_eff,
            L=max_signals,
            w=w,
            min_purity=min_purity,
            null_weight=null_wt,
        )
        pip = result["pip"]
        method = "SuSiE"
        params = {
            "L": max_signals,
            "coverage": coverage,
            "min_purity": min_purity,
            "w": w,
            "converged": result["converged"],
            "n_iter": result["n_iter"],
            "n_eff": n_eff,
        }
        if result["converged"]:
            print(f"  SuSiE converged in {result['n_iter']} iterations")
        else:
            print(f"  SuSiE did not converge in {result['n_iter']} iterations (ELBO tolerance {1e-3})")

        df["pip"] = pip
        credible_sets = build_credible_sets_susie(
            alpha=result["alpha"], df=df, R=R,
            coverage=coverage, min_purity=min_purity,
        )
    else:
        print(f"  Method: ABF (no LD matrix; W={w})")
        pip = compute_abf(df, w=w)
        method = "ABF"
        params = {"coverage": coverage, "w": w}
        df["pip"] = pip
        credible_sets = build_credible_set_abf(pip=pip, df=df, coverage=coverage)

    # --- Step 3: Annotate CS membership ---
    df["cs_membership"] = ""
    for cs in credible_sets:
        for v in cs["variants"]:
            mask = df["rsid"] == v["rsid"]
            df.loc[mask, "cs_membership"] = cs["cs_id"]

    # Print summary
    print(f"\n  Results:")
    print(f"    Lead variant: {df.loc[df['pip'].idxmax(), 'rsid']} (PIP = {df['pip'].max():.4f})")
    print(f"    Credible sets: {len(credible_sets)}")
    for cs in credible_sets:
        purity_str = f", purity={cs['purity']:.3f}" if cs.get("purity") is not None else ""
        flag = " [LOW PURITY]" if not cs.get("pure", True) else ""
        print(f"      {cs['cs_id']}: {cs['size']} variants, "
              f"coverage={cs['coverage']*100:.1f}%{purity_str}{flag}")
    print(f"    Variants with PIP ≥ 0.5: {(df['pip'] >= 0.5).sum()}")
    print(f"    Variants with PIP ≥ 0.1: {(df['pip'] >= 0.1).sum()}")
    print()

    # --- Step 4: Write outputs ---
    print("  Writing report...")
    report_md = generate_markdown(df, credible_sets, method, params, input_path=input_label)
    (output_dir / "report.md").write_text(report_md)

    print("  Writing tables...")
    write_tables(output_dir, df, credible_sets)

    if make_figures:
        print("  Generating figures...")
        generate_figures(output_dir, df, credible_sets, R=R, gene_track=gene_track)

    print("  Writing fine_mapping.json...")
    results = {
        "method": method,
        "params": params,
        "n_variants": len(df),
        "n_credible_sets": len(credible_sets),
        "lead_rsid": str(df.loc[df["pip"].idxmax(), "rsid"]),
        "lead_pip": float(df["pip"].max()),
        "credible_sets": credible_sets,
        "pip_above_0_1": int((df["pip"] >= 0.1).sum()),
        "pip_above_0_5": int((df["pip"] >= 0.5).sum()),
    }
    (output_dir / "fine_mapping.json").write_text(json.dumps(results, indent=2, default=str))

    if _HAS_COMMON:
        write_result_json(
            output_dir=output_dir,
            skill="fine-mapping",
            version="0.1.0",
            summary={
                "method": method,
                "n_variants": len(df),
                "n_credible_sets": len(credible_sets),
                "lead_rsid": results["lead_rsid"],
                "lead_pip": results["lead_pip"],
            },
            data=results,
        )

    # Reproducibility bundle
    cmd = _build_cmd(sumstats_path, ld_path, output_dir, chr_, start, end,
                     max_signals, coverage, min_purity, w, demo)
    write_reproducibility(output_dir, cmd, params)

    print(f"  Report: {output_dir / 'report.md'}")
    print(f"  Full output: {output_dir}/")
    print(f"\n  {DISCLAIMER}")

    return results


def _build_cmd(sumstats_path, ld_path, output_dir, chr_, start, end,
               max_signals, coverage, min_purity, w, demo) -> str:
    parts = ["python skills/fine-mapping/fine_mapping.py"]
    if demo:
        parts.append("--demo")
    else:
        parts.append(f"--sumstats {sumstats_path}")
        if ld_path:
            parts.append(f"--ld {ld_path}")
        if chr_:
            parts += [f"--chr {chr_}", f"--start {start}", f"--end {end}"]
    parts += [
        f"--max-signals {max_signals}",
        f"--coverage {coverage}",
        f"--min-purity {min_purity}",
        f"--prior-variance {w}",
        f"--output {output_dir}",
    ]
    return " \\\n  ".join(parts)


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="SuSiE Fine-Mapper — statistical fine-mapping of GWAS loci"
    )
    parser.add_argument("--sumstats", help="GWAS summary statistics file (TSV/CSV/TXT)")
    parser.add_argument("--ld", help="Pre-computed LD matrix (.npy or .tsv)")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true",
                        help="Run with synthetic demo data (200 variants, 2 causal signals)")
    parser.add_argument("--chr", dest="chr_", help="Chromosome for locus filter (e.g. 1)")
    parser.add_argument("--start", type=int, help="Start position for locus filter")
    parser.add_argument("--end", type=int, help="End position for locus filter")
    parser.add_argument("--max-signals", type=int, default=10,
                        help="Max causal signals L for SuSiE (default: 10)")
    parser.add_argument("--coverage", type=float, default=0.95,
                        help="Credible set coverage threshold (default: 0.95)")
    parser.add_argument("--min-purity", type=float, default=0.5,
                        help="Min average pairwise |r| within CS (default: 0.5)")
    parser.add_argument("--prior-variance", type=float, default=0.04,
                        help="Prior variance W for ABF/SuSiE (default: 0.04, Wakefield 2009)")
    parser.add_argument("--no-figures", action="store_true", help="Skip figure generation")
    parser.add_argument("--gene-track", action="store_true",
                        help="Fetch gene annotations from Ensembl and add a gene track below the regional association plot (requires internet)")

    args = parser.parse_args()

    if not args.sumstats and not args.demo:
        parser.print_help()
        print("\nError: provide --sumstats or --demo")
        sys.exit(1)

    print("SuSiE Fine-Mapper")
    print("=" * 60)
    print()

    run_finemapping(
        sumstats_path=Path(args.sumstats) if args.sumstats else None,
        ld_path=Path(args.ld) if args.ld else None,
        output_dir=Path(args.output),
        chr_=args.chr_,
        start=args.start,
        end=args.end,
        max_signals=args.max_signals,
        coverage=args.coverage,
        min_purity=args.min_purity,
        w=args.prior_variance,
        make_figures=not args.no_figures,
        gene_track=args.gene_track,
        demo=args.demo,
    )


if __name__ == "__main__":
    main()
