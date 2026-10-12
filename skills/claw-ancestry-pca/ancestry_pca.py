#!/usr/bin/env python3
"""Ancestry Decomposition PCA: place cohorts in global genetic context.

Usage:
    python ancestry_pca.py --input cohort.vcf --pop-map pops.csv --output results/
    python ancestry_pca.py --demo --output results/
    python ancestry_pca.py --demo          # text summary to stdout

Computes PCA on VCF genotype data, coloured by population labels, and
produces a multi-panel figure (PC1v2, PC2v3, PC1v3, scree plot) plus a
markdown report with tables and reproducibility metadata.
"""

from __future__ import annotations

import argparse
import json
import shlex
import sys
from collections import Counter
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Shared library imports
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.parsers import parse_vcf_matrix
from clawbio.common.checksums import sha256_file as _sha256_file
from clawbio.common.report import (
    generate_report_header,
    generate_report_footer,
    write_result_json,
)
from clawbio.common.reproducibility import (
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------
SKILL_DIR = Path(__file__).resolve().parent
DEMO_VCF = _PROJECT_ROOT / "examples" / "demo_populations.vcf"
DEMO_POP_MAP = _PROJECT_ROOT / "examples" / "demo_population_map.csv"

# Colourblind-friendly palette (Wong 2011 + extras)
POP_COLOURS = {
    "AFR": "#E69F00",
    "AMR": "#56B4E9",
    "EAS": "#009E73",
    "EUR": "#F0E442",
    "SAS": "#0072B2",
    "OCE": "#D55E00",
    "MID": "#CC79A7",
}
_FALLBACK_COLOURS = [
    "#1f77b4", "#ff7f0e", "#2ca02c", "#d62728", "#9467bd",
    "#8c564b", "#e377c2", "#7f7f7f", "#bcbd22", "#17becf",
]


def _get_colour(pop: str, idx: int = 0) -> str:
    """Return a colour for a population label."""
    upper = pop.upper()
    if upper in POP_COLOURS:
        return POP_COLOURS[upper]
    return _FALLBACK_COLOURS[idx % len(_FALLBACK_COLOURS)]


# ===================================================================
# Population map
# ===================================================================


def load_population_map(
    filepath: Optional[Path], samples: List[str]
) -> Dict[str, str]:
    """Parse a CSV/TSV population map or infer from sample ID prefixes.

    Accepts files with columns like sample_id/population (flexible naming).
    """
    if filepath and Path(filepath).exists():
        filepath = Path(filepath)
        sep = "\t" if filepath.suffix in (".tsv", ".txt") else ","
        df = pd.read_csv(filepath, sep=sep)
        col_map = {}
        pop_found = sid_found = False
        for col in df.columns:
            lower = col.lower().strip()
            if not pop_found and lower in (
                "population", "ancestry", "pop", "superpopulation",
            ):
                col_map[col] = "population"
                pop_found = True
            elif not sid_found and lower in (
                "sample_id", "sample", "id", "iid",
            ):
                col_map[col] = "sample_id"
                sid_found = True
        df = df.rename(columns=col_map)
        if "population" not in df.columns:
            raise ValueError(
                "No population column found. Columns: %s" % list(df.columns)
            )
        if "sample_id" not in df.columns:
            raise ValueError(
                "No sample_id column found. Columns: %s" % list(df.columns)
            )
        return dict(zip(df["sample_id"].astype(str), df["population"].astype(str)))

    # Infer from sample ID prefix (e.g. AFR_001 -> AFR)
    pop_map = {}
    for s in samples:
        prefix = s.split("_")[0].upper()
        pop_map[s] = prefix if prefix else "UNKNOWN"
    return pop_map


# ===================================================================
# PCA computation
# ===================================================================


def compute_pca(
    geno_matrix: np.ndarray, n_components: int = 10
) -> Tuple[np.ndarray, np.ndarray]:
    """Run PCA on a genotype matrix (0/1/2/-1 encoding).

    Missing values (-1) are mean-imputed per variant.

    Returns:
        coords: (n_samples, n_components) PC coordinates
        explained_var: explained variance ratio per component
    """
    from sklearn.decomposition import PCA

    mat = geno_matrix.astype(np.float64).copy()

    # Mean-impute missing values per variant
    for j in range(mat.shape[1]):
        col = mat[:, j]
        missing = col == -1
        if missing.any():
            valid_mean = col[~missing].mean() if (~missing).any() else 0.0
            col[missing] = valid_mean

    n_components = min(n_components, mat.shape[0] - 1, mat.shape[1])
    n_components = max(1, n_components)
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(mat)
    return coords, pca.explained_variance_ratio_


# ===================================================================
# Visualisation
# ===================================================================


def plot_pca_composite(
    coords: np.ndarray,
    explained_var: np.ndarray,
    sample_pops: List[str],
    output_path: Path,
) -> None:
    """4-panel PCA composite figure.

    Panel A: PC1 vs PC2 (main structure)
    Panel B: PC2 vs PC3 (finer structure)
    Panel C: PC1 vs PC3 (third axis)
    Panel D: Scree plot (variance explained)
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pops_unique = sorted(set(sample_pops))
    colour_map = {
        pop: _get_colour(pop, i) for i, pop in enumerate(pops_unique)
    }

    n_pcs = coords.shape[1]

    fig, axes = plt.subplots(2, 2, figsize=(14, 12))
    fig.suptitle("Ancestry Decomposition PCA", fontsize=16, fontweight="bold")

    # Helper for scatter panels
    def _scatter(ax, pc_x, pc_y, title):
        if pc_x >= n_pcs or pc_y >= n_pcs:
            ax.set_visible(False)
            return
        for pop in pops_unique:
            mask = [p == pop for p in sample_pops]
            ax.scatter(
                coords[mask, pc_x], coords[mask, pc_y],
                c=colour_map[pop], label=pop,
                s=50, alpha=0.8, edgecolors="white", linewidths=0.3,
            )
        var_x = explained_var[pc_x] * 100 if pc_x < len(explained_var) else 0
        var_y = explained_var[pc_y] * 100 if pc_y < len(explained_var) else 0
        ax.set_xlabel("PC%d (%.1f%%)" % (pc_x + 1, var_x))
        ax.set_ylabel("PC%d (%.1f%%)" % (pc_y + 1, var_y))
        ax.set_title(title)
        ax.grid(True, alpha=0.3)
        ax.legend(title="Population", fontsize=7, loc="best", framealpha=0.9)

    # Panel A: PC1 vs PC2
    _scatter(axes[0, 0], 0, 1, "A. PC1 vs PC2")

    # Panel B: PC2 vs PC3
    _scatter(axes[0, 1], 1, 2, "B. PC2 vs PC3")

    # Panel C: PC1 vs PC3
    _scatter(axes[1, 0], 0, 2, "C. PC1 vs PC3")

    # Panel D: Scree plot
    ax_scree = axes[1, 1]
    n_show = min(len(explained_var), 10)
    pcs = range(1, n_show + 1)
    var_pcts = [v * 100 for v in explained_var[:n_show]]
    cum_var = np.cumsum(var_pcts)

    ax_scree.bar(pcs, var_pcts, color="#2196F3", alpha=0.8, label="Individual")
    ax_scree.plot(pcs, cum_var, "o-", color="#F44336", label="Cumulative")
    ax_scree.set_xlabel("Principal Component")
    ax_scree.set_ylabel("Variance Explained (%)")
    ax_scree.set_title("D. Scree Plot")
    ax_scree.set_xticks(list(pcs))
    ax_scree.legend(loc="center right")
    ax_scree.grid(True, axis="y", alpha=0.3)

    plt.tight_layout(rect=[0, 0, 1, 0.96])
    fig.savefig(output_path, dpi=200)
    plt.close(fig)


# ===================================================================
# Report generation
# ===================================================================


def generate_report(
    input_path: Path,
    pop_map_path: Optional[Path],
    n_samples: int,
    n_variants: int,
    n_components: int,
    explained_var: np.ndarray,
    pop_counts: Dict[str, int],
    output_dir: Path,
    figures_generated: bool,
) -> str:
    """Generate the ancestry PCA markdown report."""
    header = generate_report_header(
        title="Ancestry Decomposition PCA",
        skill_name="claw-ancestry-pca",
        input_files=[input_path] + ([pop_map_path] if pop_map_path else []),
        extra_metadata={
            "Samples": str(n_samples),
            "Populations": str(len(pop_counts)),
            "Variants analysed": str(n_variants),
            "Components computed": str(n_components),
        },
    )

    # Variance table
    n_show = min(len(explained_var), 10)
    var_rows = []
    cum = 0.0
    for i in range(n_show):
        pct = explained_var[i] * 100
        cum += pct
        var_rows.append("| PC%d | %.2f%% | %.2f%% |" % (i + 1, pct, cum))
    var_table = "\n".join(var_rows)

    # Population table
    total = sum(pop_counts.values())
    pop_rows = []
    for pop in sorted(pop_counts.keys()):
        count = pop_counts[pop]
        pct = count / total * 100 if total > 0 else 0
        pop_rows.append("| %s | %d | %.1f%% |" % (pop, count, pct))
    pop_table = "\n".join(pop_rows)

    # Figure reference
    fig_section = ""
    if figures_generated:
        fig_section = """## PCA Plots

![PCA Composite](figures/pca_composite.png)

**Panel A** (PC1 vs PC2): Main axis of population structure.
**Panel B** (PC2 vs PC3): Finer population differentiation.
**Panel C** (PC1 vs PC3): Third axis of variation.
**Panel D** (Scree plot): Variance explained per component.
"""

    body = """## Variance Explained

| Component | Variance | Cumulative |
|-----------|----------|------------|
%(var_table)s

## Population Composition

| Population | Count | Proportion |
|------------|-------|------------|
%(pop_table)s

%(fig_section)s
## Interpretation

- **PC1** typically captures the largest axis of global differentiation
- **PC2** separates major continental groups
- **PC3** often reveals finer substructure within continental groups
- Examine the scree plot to assess how many PCs carry meaningful signal

## Methods

- **Tool**: ClawBio Ancestry PCA v0.1.0
- **PCA**: scikit-learn PCA on mean-imputed genotype matrix (0/1/2 encoding)
- **Missing data**: Per-variant mean imputation
- **Input format**: VCF with biallelic SNPs

## References

- Mallick, S. et al. (2016). The Simons Genome Diversity Project. Nature, 538, 201-206.
- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio

## Reproducibility

```bash
python skills/claw-ancestry-pca/ancestry_pca.py \\
    --input %(input_name)s %(pop_map_flag)s--output %(output_name)s
```
""" % {
        "var_table": var_table,
        "pop_table": pop_table,
        "fig_section": fig_section,
        "input_name": input_path.name,
        "pop_map_flag": ("--pop-map %s " % pop_map_path.name) if pop_map_path else "",
        "output_name": output_dir.name,
    }

    footer = generate_report_footer()
    return header + body + footer


# ===================================================================
# Main pipeline
# ===================================================================


def run_analysis(
    input_path: Path,
    pop_map_path: Optional[Path],
    output_dir: Path,
    n_components: int = 10,
    no_figures: bool = False,
    demo: bool = False,
) -> dict:
    """Full ancestry PCA pipeline with report output."""
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"

    if not no_figures:
        figures_dir.mkdir(exist_ok=True)
    tables_dir.mkdir(exist_ok=True)

    # Parse VCF
    print("Parsing VCF...")
    samples, variant_ids, geno_matrix = parse_vcf_matrix(input_path)
    n_samples, n_variants = geno_matrix.shape
    print("  %d samples, %d variants" % (n_samples, n_variants))

    # Population assignments
    pop_map = load_population_map(pop_map_path, samples)
    sample_pops = [pop_map.get(s, "UNKNOWN") for s in samples]
    pop_counts = dict(Counter(sample_pops))
    pops = sorted(pop_counts.keys())
    print("  Populations: %s" % ", ".join(
        "%s (n=%d)" % (p, pop_counts[p]) for p in pops
    ))

    # PCA
    print("Computing PCA (%d components)..." % n_components)
    coords, explained_var = compute_pca(geno_matrix, n_components)
    actual_components = coords.shape[1]
    print("  PC1: %.1f%%  PC2: %.1f%%" % (
        explained_var[0] * 100,
        explained_var[1] * 100 if len(explained_var) > 1 else 0,
    ))

    # Save tables
    coord_df = pd.DataFrame(
        coords,
        columns=["PC%d" % (i + 1) for i in range(actual_components)],
    )
    coord_df.insert(0, "sample_id", samples)
    coord_df.insert(1, "population", sample_pops)
    coord_df.to_csv(tables_dir / "pc_coordinates.csv", index=False)

    var_df = pd.DataFrame({
        "component": ["PC%d" % (i + 1) for i in range(len(explained_var))],
        "variance_explained": explained_var,
        "cumulative": np.cumsum(explained_var),
    })
    var_df.to_csv(tables_dir / "variance_explained.csv", index=False)

    # Figures
    figures_generated = False
    if not no_figures:
        print("Generating figures...")
        try:
            plot_pca_composite(
                coords, explained_var, sample_pops,
                figures_dir / "pca_composite.png",
            )
            figures_generated = True
        except ImportError as e:
            print("  Warning: %s — figures skipped." % e, file=sys.stderr)

    # Report
    print("Generating report...")
    report_text = generate_report(
        input_path=input_path,
        pop_map_path=pop_map_path,
        n_samples=n_samples,
        n_variants=n_variants,
        n_components=actual_components,
        explained_var=explained_var,
        pop_counts=pop_counts,
        output_dir=output_dir,
        figures_generated=figures_generated,
    )
    report_path = output_dir / "report.md"
    report_path.write_text(report_text)

    # result.json
    write_result_json(
        output_dir=output_dir,
        skill="claw-ancestry-pca",
        version="0.1.0",
        summary={
            "n_samples": n_samples,
            "n_variants": n_variants,
            "n_populations": len(pop_counts),
            "n_components": actual_components,
            "pc1_variance": round(float(explained_var[0]) * 100, 2),
            "pc2_variance": round(float(explained_var[1]) * 100, 2) if len(explained_var) > 1 else 0,
        },
        data={
            "population_counts": pop_counts,
            "variance_explained": [round(float(v), 6) for v in explained_var],
        },
        input_checksum=_sha256_file(input_path) if input_path.exists() else "",
    )

    # Reproducibility bundle
    cmd = ["python", "skills/claw-ancestry-pca/ancestry_pca.py"]
    if demo:
        cmd.append("--demo")
    else:
        cmd += ["--input", str(input_path)]
        if pop_map_path:
            cmd += ["--pop-map", str(pop_map_path)]
    cmd += ["--output", str(output_dir), "--n-components", str(n_components)]
    if no_figures:
        cmd.append("--no-figures")
    write_commands_sh(output_dir, shlex.join(cmd))
    env_path = write_environment_yml(
        output_dir, "clawbio-claw-ancestry-pca",
        pip_deps=["numpy", "pandas", "scikit-learn", "matplotlib"],
    )
    write_checksums(
        [report_path, output_dir / "result.json",
         tables_dir / "pc_coordinates.csv", tables_dir / "variance_explained.csv",
         figures_dir / "pca_composite.png", env_path],
        output_dir, anchor=output_dir,
    )

    print("\nDone.")
    print("  Report: %s" % report_path)
    if figures_generated:
        print("  Figures: %s" % figures_dir)

    return {
        "n_samples": n_samples,
        "n_variants": n_variants,
        "n_components": actual_components,
        "coords": coords,
        "explained_var": explained_var,
        "sample_pops": sample_pops,
        "pop_counts": pop_counts,
    }


def run_summary(
    input_path: Path,
    pop_map_path: Optional[Path],
    n_components: int = 10,
) -> str:
    """Quick text summary to stdout (no files written)."""
    samples, variant_ids, geno_matrix = parse_vcf_matrix(input_path)
    n_samples, n_variants = geno_matrix.shape

    pop_map = load_population_map(pop_map_path, samples)
    sample_pops = [pop_map.get(s, "UNKNOWN") for s in samples]
    pop_counts = dict(Counter(sample_pops))

    coords, explained_var = compute_pca(geno_matrix, n_components)

    lines = []
    lines.append("ANCESTRY DECOMPOSITION PCA")
    lines.append("")
    lines.append("Input: %s" % input_path.name)
    lines.append("Samples: %d" % n_samples)
    lines.append("Variants: %d" % n_variants)
    lines.append("Populations: %d" % len(pop_counts))
    lines.append("")
    lines.append("== VARIANCE EXPLAINED ==")
    cum = 0.0
    for i in range(min(len(explained_var), 5)):
        pct = explained_var[i] * 100
        cum += pct
        lines.append("  PC%d: %.2f%% (cumulative: %.2f%%)" % (i + 1, pct, cum))
    lines.append("")
    lines.append("== POPULATION COUNTS ==")
    for pop in sorted(pop_counts.keys()):
        lines.append("  %s: %d" % (pop, pop_counts[pop]))

    return "\n".join(lines)


# ===================================================================
# CLI
# ===================================================================


def main() -> None:
    parser = argparse.ArgumentParser(
        description="ClawBio Ancestry PCA: population structure from VCF genotype data"
    )
    parser.add_argument("--input", "-i", help="Input VCF file (.vcf or .vcf.gz)")
    parser.add_argument("--pop-map", "-p", default=None,
                        help="Population map CSV/TSV (columns: sample_id, population)")
    parser.add_argument("--output", "-o", default=None,
                        help="Output directory (enables full report + figures)")
    parser.add_argument("--demo", action="store_true",
                        help="Run with bundled demo data (50 samples, 5 pops)")
    parser.add_argument("--n-components", "-n", type=int, default=10,
                        help="Number of PCs to compute (default: 10)")
    parser.add_argument("--no-figures", action="store_true",
                        help="Skip figure generation")
    args = parser.parse_args()

    if args.demo:
        input_path = DEMO_VCF
        pop_map_path = DEMO_POP_MAP
    elif args.input:
        input_path = Path(args.input)
        pop_map_path = Path(args.pop_map) if args.pop_map else None
    else:
        parser.error("Provide --input or --demo")
        return

    if not input_path.exists():
        print("Error: input file not found: %s" % input_path, file=sys.stderr)
        sys.exit(1)

    # Stdout summary mode (no --output)
    if not args.output:
        text = run_summary(input_path, pop_map_path, args.n_components)
        print(text)
        sys.exit(0)

    # Full report mode
    output_dir = Path(args.output)
    run_analysis(
        input_path=input_path,
        pop_map_path=pop_map_path,
        output_dir=output_dir,
        n_components=args.n_components,
        no_figures=args.no_figures,
        demo=args.demo,
    )


if __name__ == "__main__":
    main()
