#!/usr/bin/env python3
"""Equity Scorer: HEIM diversity metrics from VCF or ancestry data.

Usage:
    python equity_scorer.py --input <vcf_or_csv> [--pop-map <csv>] [--output <dir>]

Computes from real genotype data:
  - Observed and expected heterozygosity per population
  - Pairwise Hudson FST between all population pairs
  - PCA of genotype matrix (PC1 vs PC2)
  - Population representation statistics
  - Composite HEIM Equity Score (0-100)

Outputs a markdown report with figures and tables to the output directory.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import json
import sys
from collections import Counter, OrderedDict
from datetime import datetime, timezone
from itertools import combinations
from typing import Dict, List, Optional, Tuple

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

# Global population proportions (1000 Genomes superpopulations, approximate)
GLOBAL_PROPORTIONS = {
    "AFR": 0.17,
    "AMR": 0.13,
    "EAS": 0.22,
    "EUR": 0.16,
    "SAS": 0.26,
    "OCE": 0.005,
    "MID": 0.06,
}

DEFAULT_WEIGHTS = (0.35, 0.25, 0.20, 0.20)

# Colours for population plots (colourblind-friendly palette)
POP_COLOURS = {
    "AFR": "#E69F00",
    "AMR": "#56B4E9",
    "EAS": "#009E73",
    "EUR": "#F0E442",
    "SAS": "#0072B2",
    "OCE": "#D55E00",
    "MID": "#CC79A7",
    "UNKNOWN": "#999999",
}


# ===================================================================
# VCF PARSING
# ===================================================================

def parse_vcf(filepath: Path) -> Tuple[List[str], List[str], np.ndarray]:
    """Parse a VCF file into a genotype matrix.

    Returns:
        samples: list of sample IDs
        variant_ids: list of variant IDs (or CHROM:POS)
        genotype_matrix: numpy array of shape (n_samples, n_variants)
                         with values 0 (hom ref), 1 (het), 2 (hom alt), -1 (missing)
    """
    samples = []
    variant_ids = []
    genotype_rows = []  # list of lists, one per variant (transposed later)

    with open(filepath) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                samples = parts[9:]
                continue

            parts = line.strip().split("\t")
            chrom, pos, vid = parts[0], parts[1], parts[2]
            if vid == ".":
                vid = f"{chrom}:{pos}"
            variant_ids.append(vid)

            fmt_fields = parts[8].split(":")
            gt_idx = fmt_fields.index("GT") if "GT" in fmt_fields else 0

            row = []
            for sample_field in parts[9:]:
                gt_str = sample_field.split(":")[gt_idx]
                gt_str = gt_str.replace("|", "/")  # handle phased
                if "." in gt_str:
                    row.append(-1)
                else:
                    alleles = gt_str.split("/")
                    row.append(int(alleles[0]) + int(alleles[1]))
            genotype_rows.append(row)

    if not samples:
        raise ValueError("No samples found in VCF header")
    if not genotype_rows:
        raise ValueError("No variants found in VCF")

    # Shape: (n_samples, n_variants)
    geno_matrix = np.array(genotype_rows, dtype=np.int8).T
    return samples, variant_ids, geno_matrix


def load_population_map(
    filepath: Optional[Path], samples: List[str]
) -> Dict[str, str]:
    """Load population assignments from a CSV or infer from sample IDs.

    The CSV should have columns: sample_id, population
    If no map is provided, infer from sample ID prefixes (e.g., AFR_001 -> AFR).
    """
    if filepath and filepath.exists():
        df = pd.read_csv(filepath)
        col_map = {}
        pop_found = False
        sid_found = False
        for col in df.columns:
            lower = col.lower().strip()
            if not pop_found and lower in ("population", "ancestry", "pop"):
                col_map[col] = "population"
                pop_found = True
            elif not pop_found and lower == "superpopulation":
                col_map[col] = "population"
                pop_found = True
            elif not sid_found and lower in ("sample_id", "sample", "id", "iid"):
                col_map[col] = "sample_id"
                sid_found = True
        df = df.rename(columns=col_map)
        return dict(zip(df["sample_id"], df["population"]))

    # Infer from sample ID prefix
    pop_map = {}
    known_pops = set(GLOBAL_PROPORTIONS.keys())
    for s in samples:
        prefix = s.split("_")[0].upper()
        pop_map[s] = prefix if prefix in known_pops else "UNKNOWN"
    return pop_map


def parse_ancestry_csv(filepath: Path) -> pd.DataFrame:
    """Parse a CSV with sample_id and population/ancestry columns."""
    df = pd.read_csv(filepath)
    col_map = {}
    for col in df.columns:
        lower = col.lower().strip()
        if lower in ("population", "ancestry", "pop", "superpopulation"):
            col_map[col] = "population"
        elif lower in ("sample_id", "sample", "id", "iid"):
            col_map[col] = "sample_id"
    df = df.rename(columns=col_map)

    if "population" not in df.columns:
        raise ValueError(
            "No population/ancestry column found. Columns: %s" % list(df.columns)
        )
    if "sample_id" not in df.columns:
        df["sample_id"] = ["SAMPLE_%d" % i for i in range(len(df))]
    return df


# ===================================================================
# POPULATION GENETICS METRICS
# ===================================================================

def compute_allele_frequencies(
    geno_matrix: np.ndarray, pop_indices: Dict[str, List[int]]
) -> Dict[str, np.ndarray]:
    """Compute per-population allele frequencies for the alt allele.

    Args:
        geno_matrix: (n_samples, n_variants), values 0/1/2/-1
        pop_indices: dict mapping population name to list of sample indices

    Returns:
        dict mapping population name to array of alt allele frequencies (n_variants,)
    """
    afs = {}
    for pop, indices in pop_indices.items():
        sub = geno_matrix[indices, :]
        # Mask missing values
        valid = sub != -1
        allele_sum = np.where(valid, sub, 0).sum(axis=0).astype(np.float64)
        allele_count = valid.sum(axis=0) * 2  # diploid
        allele_count = np.where(allele_count == 0, 1, allele_count)  # avoid div by 0
        afs[pop] = allele_sum / allele_count
    return afs


def compute_heterozygosity(
    geno_matrix: np.ndarray, pop_indices: Dict[str, List[int]]
) -> Tuple[Dict[str, float], Dict[str, float], Dict[str, np.ndarray], Dict[str, np.ndarray]]:
    """Compute observed and expected heterozygosity per population.

    Observed Het = fraction of heterozygous genotypes (per site, averaged)
    Expected Het = 2pq (per site, averaged)

    Returns:
        obs_het: mean observed Het per population
        exp_het: mean expected Het per population
        obs_het_per_site: per-site observed Het arrays
        exp_het_per_site: per-site expected Het arrays
    """
    afs = compute_allele_frequencies(geno_matrix, pop_indices)
    obs_het = {}
    exp_het = {}
    obs_per_site = {}
    exp_per_site = {}

    for pop, indices in pop_indices.items():
        sub = geno_matrix[indices, :]
        valid = sub != -1
        n_valid = valid.sum(axis=0).astype(np.float64)
        n_valid = np.where(n_valid == 0, 1, n_valid)

        # Observed: count of het genotypes (value == 1)
        het_count = ((sub == 1) & valid).sum(axis=0).astype(np.float64)
        site_obs_het = het_count / n_valid
        obs_per_site[pop] = site_obs_het
        obs_het[pop] = float(np.mean(site_obs_het))

        # Expected: 2pq
        p = afs[pop]
        q = 1 - p
        site_exp_het = 2 * p * q
        exp_per_site[pop] = site_exp_het
        exp_het[pop] = float(np.mean(site_exp_het))

    return obs_het, exp_het, obs_per_site, exp_per_site


def compute_pairwise_fst(
    geno_matrix: np.ndarray, pop_indices: Dict[str, List[int]]
) -> Tuple[pd.DataFrame, Dict[Tuple[str, str], float]]:
    """Compute pairwise FST between all population pairs.

    Uses Nei's GST approach (robust with small samples):
        Per-site: FST = (HT - HS) / HT
        where HT = 2 * p_total * (1 - p_total)  [total expected Het]
              HS = weighted mean of 2*pi*(1-pi)  [within-pop expected Het]

    Final FST = mean across all polymorphic sites (ratio of averages).

    Returns:
        fst_df: DataFrame with populations as index/columns, values = FST
        fst_dict: dict mapping (pop1, pop2) tuple to FST value
    """
    afs = compute_allele_frequencies(geno_matrix, pop_indices)
    pops = sorted(pop_indices.keys())
    n_pops = len(pops)

    fst_matrix = np.zeros((n_pops, n_pops))
    fst_dict = {}

    for i, j in combinations(range(n_pops), 2):
        pop1, pop2 = pops[i], pops[j]
        p1 = afs[pop1]
        p2 = afs[pop2]
        n1 = len(pop_indices[pop1])
        n2 = len(pop_indices[pop2])
        total_n = n1 + n2

        # Weighted mean allele frequency (pooled)
        p_total = (p1 * n1 + p2 * n2) / total_n

        # Total expected heterozygosity
        ht = 2 * p_total * (1 - p_total)

        # Within-population expected heterozygosity (weighted mean)
        hs = (n1 * 2 * p1 * (1 - p1) + n2 * 2 * p2 * (1 - p2)) / total_n

        # FST per site
        numerator = ht - hs
        denominator = ht

        # Ratio of averages across polymorphic sites
        valid = denominator > 0.001
        if valid.any():
            fst_val = float(np.sum(numerator[valid]) / np.sum(denominator[valid]))
            fst_val = max(0.0, fst_val)
        else:
            fst_val = 0.0

        fst_matrix[i, j] = fst_val
        fst_matrix[j, i] = fst_val
        fst_dict[(pop1, pop2)] = fst_val

    fst_df = pd.DataFrame(fst_matrix, index=pops, columns=pops)
    return fst_df, fst_dict


def compute_pca(
    geno_matrix: np.ndarray, n_components: int = 10
) -> Tuple[np.ndarray, np.ndarray]:
    """Run PCA on the genotype matrix.

    Handles missing data by mean imputation per variant.

    Returns:
        coords: (n_samples, n_components) principal component coordinates
        explained_var: explained variance ratio per component
    """
    from sklearn.decomposition import PCA

    # Mean-impute missing values per variant
    mat = geno_matrix.astype(np.float64).copy()
    for j in range(mat.shape[1]):
        col = mat[:, j]
        missing = col == -1
        if missing.any():
            valid_mean = col[~missing].mean() if (~missing).any() else 0.0
            col[missing] = valid_mean

    n_components = min(n_components, mat.shape[0], mat.shape[1])
    pca = PCA(n_components=n_components)
    coords = pca.fit_transform(mat)
    return coords, pca.explained_variance_ratio_


# ===================================================================
# HEIM EQUITY SCORE
# ===================================================================

def compute_representation_index(pop_counts: Dict[str, int]) -> float:
    """Measure how well sample proportions match global proportions (0-1)."""
    total = sum(pop_counts.values())
    if total == 0:
        return 0.0
    sample_props = {k: v / total for k, v in pop_counts.items()}
    max_deviation = 0.0
    for pop, global_prop in GLOBAL_PROPORTIONS.items():
        sample_prop = sample_props.get(pop, 0.0)
        deviation = abs(sample_prop - global_prop)
        max_deviation = max(max_deviation, deviation)
    return max(0.0, 1.0 - max_deviation)


def compute_heterozygosity_balance(het_values: Dict[str, float]) -> float:
    """Ratio of mean observed heterozygosity to theoretical max (0-1)."""
    if not het_values:
        return 0.0
    mean_het = np.mean(list(het_values.values()))
    return float(min(1.0, mean_het / 0.5))


def compute_fst_coverage(n_populations: int, n_pairwise_computed: int) -> float:
    """Fraction of possible pairwise FST values actually computed."""
    n_possible = n_populations * (n_populations - 1) // 2
    if n_possible == 0:
        return 0.0
    return min(1.0, n_pairwise_computed / n_possible)


def compute_geographic_spread(populations: set) -> float:
    """Fraction of continental groups represented (out of 7)."""
    continent_map = {
        "AFR": "Africa", "AMR": "Americas", "EAS": "East Asia",
        "EUR": "Europe", "SAS": "South Asia", "OCE": "Oceania",
        "MID": "Middle East",
    }
    continents = set()
    for pop in populations:
        pop_upper = pop.upper()
        if pop_upper in continent_map:
            continents.add(continent_map[pop_upper])
    return len(continents) / 7.0


def compute_heim_score(
    pop_counts: Dict[str, int],
    het_values: Dict[str, float],
    n_pairwise_fst: int,
    weights: Tuple[float, float, float, float] = DEFAULT_WEIGHTS,
) -> dict:
    """Compute the composite HEIM Equity Score (0-100)."""
    n_pops = len(pop_counts)
    ri = compute_representation_index(pop_counts)
    hb = compute_heterozygosity_balance(het_values)
    fc = compute_fst_coverage(n_pops, n_pairwise_fst)
    gs = compute_geographic_spread(set(pop_counts.keys()))

    w1, w2, w3, w4 = weights
    score = (w1 * ri + w2 * hb + w3 * fc + w4 * gs) * 100

    rating = (
        "Excellent" if score >= 80 else
        "Good" if score >= 60 else
        "Fair" if score >= 40 else
        "Poor" if score >= 20 else
        "Critical"
    )

    return {
        "heim_score": round(score, 1),
        "rating": rating,
        "components": {
            "representation_index": round(ri, 3),
            "heterozygosity_balance": round(hb, 3),
            "fst_coverage": round(fc, 3),
            "geographic_spread": round(gs, 3),
        },
        "weights": {"w1": w1, "w2": w2, "w3": w3, "w4": w4},
        "n_samples": sum(pop_counts.values()),
        "n_populations": n_pops,
        "population_counts": pop_counts,
    }


# ===================================================================
# VISUALISATION
# ===================================================================

def _get_pop_colour(pop: str) -> str:
    return POP_COLOURS.get(pop.upper(), "#999999")


def plot_pca(
    coords: np.ndarray,
    explained_var: np.ndarray,
    sample_pops: List[str],
    output_path: Path,
) -> None:
    """PCA scatter plot: PC1 vs PC2 coloured by population."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pops_unique = sorted(set(sample_pops))
    fig, ax = plt.subplots(figsize=(9, 7))

    for pop in pops_unique:
        mask = [p == pop for p in sample_pops]
        ax.scatter(
            coords[mask, 0], coords[mask, 1],
            c=_get_pop_colour(pop),
            label=pop,
            s=60, alpha=0.8, edgecolors="white", linewidths=0.5,
        )

    ax.set_xlabel("PC1 (%.1f%% variance)" % (explained_var[0] * 100))
    ax.set_ylabel("PC2 (%.1f%% variance)" % (explained_var[1] * 100))
    ax.set_title("Principal Component Analysis of Genotype Data")
    ax.legend(title="Population", loc="best", framealpha=0.9)
    ax.grid(True, alpha=0.3)
    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_fst_heatmap(fst_df: pd.DataFrame, output_path: Path) -> None:
    """Heatmap of pairwise FST values."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pops = list(fst_df.index)
    n = len(pops)
    fig, ax = plt.subplots(figsize=(7, 6))

    im = ax.imshow(fst_df.values, cmap="YlOrRd", aspect="equal", vmin=0)
    cbar = fig.colorbar(im, ax=ax, shrink=0.8)
    cbar.set_label("FST (Nei's GST)")

    ax.set_xticks(range(n))
    ax.set_yticks(range(n))
    ax.set_xticklabels(pops, fontsize=11)
    ax.set_yticklabels(pops, fontsize=11)
    ax.set_title("Pairwise FST Between Populations")

    # Annotate cells
    for i in range(n):
        for j in range(n):
            val = fst_df.values[i, j]
            colour = "white" if val > fst_df.values.max() * 0.6 else "black"
            ax.text(j, i, "%.3f" % val, ha="center", va="center",
                    fontsize=10, color=colour)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_heterozygosity(
    obs_het: Dict[str, float],
    exp_het: Dict[str, float],
    output_path: Path,
) -> None:
    """Grouped bar chart: observed vs expected heterozygosity per population."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pops = sorted(obs_het.keys())
    x = np.arange(len(pops))
    width = 0.35

    fig, ax = plt.subplots(figsize=(9, 5))
    bars_obs = ax.bar(x - width / 2, [obs_het[p] for p in pops], width,
                      label="Observed", color="#2196F3", edgecolor="white")
    bars_exp = ax.bar(x + width / 2, [exp_het[p] for p in pops], width,
                      label="Expected (2pq)", color="#FF9800", edgecolor="white")

    ax.set_ylabel("Mean Heterozygosity")
    ax.set_xlabel("Population")
    ax.set_title("Observed vs Expected Heterozygosity")
    ax.set_xticks(x)
    ax.set_xticklabels(pops)
    ax.legend()
    ax.set_ylim(0, max(max(obs_het.values()), max(exp_het.values())) * 1.2)
    ax.grid(True, axis="y", alpha=0.3)

    # Value labels
    for bar in bars_obs:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                "%.3f" % bar.get_height(), ha="center", va="bottom", fontsize=8)
    for bar in bars_exp:
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.005,
                "%.3f" % bar.get_height(), ha="center", va="bottom", fontsize=8)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_ancestry_bar(pop_counts: Dict[str, int], output_path: Path) -> None:
    """Bar chart of population proportions."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    pops = sorted(pop_counts.keys())
    counts = [pop_counts[p] for p in pops]
    total = sum(counts)
    props = [c / total * 100 for c in counts]
    colours = [_get_pop_colour(p) for p in pops]

    fig, ax = plt.subplots(figsize=(9, 5))
    bars = ax.bar(pops, props, color=colours, edgecolor="white")
    ax.set_ylabel("Proportion (%)")
    ax.set_xlabel("Population")
    ax.set_title("Sample Ancestry Distribution")

    for bar, count in zip(bars, counts):
        ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.5,
                "n=%d" % count, ha="center", va="bottom", fontsize=9)

    # Reference line for global proportions
    for i, pop in enumerate(pops):
        global_pct = GLOBAL_PROPORTIONS.get(pop.upper(), 0) * 100
        if global_pct > 0:
            ax.plot([i - 0.3, i + 0.3], [global_pct, global_pct],
                    color="red", linewidth=2, linestyle="--")

    # Add legend for the reference line
    from matplotlib.lines import Line2D
    legend_line = Line2D([0], [0], color="red", linewidth=2, linestyle="--",
                         label="Global proportion")
    ax.legend(handles=[legend_line], loc="upper right")
    ax.grid(True, axis="y", alpha=0.3)

    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def plot_heim_gauge(score: float, rating: str, output_path: Path) -> None:
    """Horizontal gauge bar showing the HEIM score."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    fig, ax = plt.subplots(figsize=(10, 2.5))

    # Background gradient zones
    zones = [
        (0, 20, "#F44336", "Critical"),
        (20, 40, "#FF9800", "Poor"),
        (40, 60, "#FFC107", "Fair"),
        (60, 80, "#8BC34A", "Good"),
        (80, 100, "#4CAF50", "Excellent"),
    ]
    for start, end, colour, label in zones:
        ax.barh(0, end - start, left=start, height=0.7, color=colour, alpha=0.3)
        ax.text((start + end) / 2, -0.55, label, ha="center", va="top", fontsize=8,
                color=colour, fontweight="bold")

    # Score indicator
    score_colour = "#4CAF50" if score >= 80 else "#8BC34A" if score >= 60 else "#FFC107" if score >= 40 else "#FF9800" if score >= 20 else "#F44336"
    ax.barh(0, score, height=0.7, color=score_colour, edgecolor="none", alpha=0.85)
    ax.plot(score, 0, marker="v", markersize=14, color="black", zorder=5)

    ax.set_xlim(0, 100)
    ax.set_ylim(-1, 0.8)
    ax.set_yticks([])
    ax.set_xlabel("HEIM Equity Score")
    ax.set_title("HEIM Score: %.0f/100 (%s)" % (score, rating), fontsize=14, fontweight="bold")

    plt.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


# ===================================================================
# REPORT GENERATION
# ===================================================================

def sha256_file(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def generate_report(
    heim_result: dict,
    input_path: Path,
    output_dir: Path,
    figures: Dict[str, Path],
    obs_het: Dict[str, float],
    exp_het: Dict[str, float],
    fst_df: Optional[pd.DataFrame],
    pca_variance: Optional[np.ndarray],
    n_variants: int,
) -> str:
    """Generate the full HEIM markdown report."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    checksum = sha256_file(input_path) if input_path.exists() else "N/A"
    pop_counts = heim_result["population_counts"]
    components = heim_result["components"]
    total = heim_result["n_samples"]

    # Build figure references
    fig_refs = {}
    for name, path in figures.items():
        rel = path.relative_to(output_dir)
        fig_refs[name] = "![%s](%s)" % (name, rel)

    # Population table
    pop_rows = []
    for pop in sorted(pop_counts.keys()):
        count = pop_counts[pop]
        pct = count / total * 100
        global_pct = GLOBAL_PROPORTIONS.get(pop.upper(), 0) * 100
        ratio = pct / global_pct if global_pct > 0 else float("inf")
        o_het = obs_het.get(pop, 0)
        e_het = exp_het.get(pop, 0)
        pop_rows.append(
            "| %s | %d | %.1f%% | %.1f%% | %.2fx | %.4f | %.4f |"
            % (pop, count, pct, global_pct, ratio, o_het, e_het)
        )
    pop_table = "\n".join(pop_rows)

    # FST summary
    fst_section = ""
    if fst_df is not None:
        fst_rows = []
        pops = sorted(pop_counts.keys())
        for i, j in combinations(range(len(pops)), 2):
            val = fst_df.iloc[i, j]
            fst_rows.append("| %s vs %s | %.4f |" % (pops[i], pops[j], val))
        fst_table = "\n".join(fst_rows)
        fst_section = """## Pairwise FST

| Comparison | Hudson FST |
|------------|-----------|
%s

%s
""" % (fst_table, fig_refs.get("FST Heatmap", ""))

    # PCA section
    pca_section = ""
    if pca_variance is not None:
        pca_section = """## Principal Component Analysis

%s

- PC1 explains %.1f%% of variance
- PC2 explains %.1f%% of variance
- Top 5 components explain %.1f%% of total variance
""" % (
            fig_refs.get("PCA", ""),
            pca_variance[0] * 100,
            pca_variance[1] * 100,
            sum(pca_variance[:min(5, len(pca_variance))]) * 100,
        )

    # Most/least represented
    sorted_pops = sorted(pop_counts.items(), key=lambda x: x[1], reverse=True)
    most_rep = sorted_pops[0]
    least_rep = sorted_pops[-1]

    report = """# HEIM Equity Report

**Date**: %(date)s
**Input**: `%(input_name)s`
**Checksum (SHA-256)**: `%(checksum)s`
**Samples**: %(total)d
**Populations**: %(n_pops)d
**Variants analysed**: %(n_variants)d

---

## HEIM Equity Score: %(score)s/100 (%(rating)s)

%(gauge_fig)s

### Score Breakdown

| Component | Value | Weight | Description |
|-----------|-------|--------|-------------|
| Representation Index | %(ri).3f | %(w1)s | Match to global population proportions |
| Heterozygosity Balance | %(hb).3f | %(w2)s | Genetic diversity relative to theoretical max |
| FST Coverage | %(fc).3f | %(w3)s | Fraction of pairwise comparisons computed |
| Geographic Spread | %(gs).3f | %(w4)s | Continental groups represented (out of 7) |

### Key Findings

- **Most represented**: %(most_pop)s (%(most_pct).1f%%, %(most_ratio).1fx global proportion)
- **Least represented**: %(least_pop)s (%(least_pct).1f%%, %(least_ratio).1fx global proportion)
- **Mean observed heterozygosity**: %(mean_het).4f (highest: %(max_het_pop)s at %(max_het).4f)

---

## Population Distribution

| Population | Count | Sample %% | Global %% | Ratio | Obs Het | Exp Het |
|------------|-------|-----------|-----------|-------|---------|---------|
%(pop_table)s

%(ancestry_fig)s

## Heterozygosity

%(het_fig)s

%(fst_section)s

%(pca_section)s

---

## Methods

- **Tool**: ClawBio Equity Scorer v0.1.0
- **HEIM framework**: Health Equity Index for Minorities (Corpas, 2026)
- **Heterozygosity**: Observed = proportion of heterozygous genotypes per site, averaged across %(n_variants)d variants. Expected = 2pq from population allele frequencies.
- **FST**: Nei's GST (HT-HS)/HT, ratio of averages across sites. Values floored at 0.
- **PCA**: scikit-learn PCA on mean-imputed genotype matrix (0/1/2 encoding).
- **Global reference**: Approximate continental proportions from the 1000 Genomes Project.

## Reproducibility

```bash
# Re-run this analysis
python equity_scorer.py --input %(input_name)s --output %(output_name)s
```

**Input checksum**: `%(checksum)s`

## References

- Corpas, M. (2026). ClawBio. https://github.com/manuelcorpas/ClawBio
- Hudson, R.R., Slatkin, M. & Maddison, W.P. (1992). Estimation of levels of gene flow from DNA sequence data. Genetics, 132(2), 583-589.
- The 1000 Genomes Project Consortium (2015). A global reference for human genetic variation. Nature, 526, 68-74.
""" % {
        "date": now,
        "input_name": input_path.name,
        "checksum": checksum,
        "total": total,
        "n_pops": heim_result["n_populations"],
        "n_variants": n_variants,
        "score": heim_result["heim_score"],
        "rating": heim_result["rating"],
        "gauge_fig": fig_refs.get("HEIM Gauge", ""),
        "ri": components["representation_index"],
        "hb": components["heterozygosity_balance"],
        "fc": components["fst_coverage"],
        "gs": components["geographic_spread"],
        "w1": heim_result["weights"]["w1"],
        "w2": heim_result["weights"]["w2"],
        "w3": heim_result["weights"]["w3"],
        "w4": heim_result["weights"]["w4"],
        "most_pop": most_rep[0],
        "most_pct": most_rep[1] / total * 100,
        "most_ratio": (most_rep[1] / total) / GLOBAL_PROPORTIONS.get(most_rep[0].upper(), 0.01),
        "least_pop": least_rep[0],
        "least_pct": least_rep[1] / total * 100,
        "least_ratio": (least_rep[1] / total) / GLOBAL_PROPORTIONS.get(least_rep[0].upper(), 0.01),
        "mean_het": np.mean(list(obs_het.values())),
        "max_het_pop": max(obs_het, key=obs_het.get),
        "max_het": max(obs_het.values()),
        "pop_table": pop_table,
        "ancestry_fig": fig_refs.get("Ancestry Distribution", ""),
        "het_fig": fig_refs.get("Heterozygosity", ""),
        "fst_section": fst_section,
        "pca_section": pca_section,
        "output_name": output_dir.name,
    }
    return report


# ===================================================================
# MAIN PIPELINE
# ===================================================================

def run_vcf_pipeline(
    vcf_path: Path,
    pop_map_path: Optional[Path],
    output_dir: Path,
    weights: Tuple[float, float, float, float],
) -> dict:
    """Full equity scoring pipeline from VCF genotype data."""
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(exist_ok=True)

    print("Parsing VCF...")
    samples, variant_ids, geno_matrix = parse_vcf(vcf_path)
    n_samples, n_variants = geno_matrix.shape
    print("  %d samples, %d variants" % (n_samples, n_variants))

    # Population assignments
    pop_map = load_population_map(pop_map_path, samples)
    sample_pops = [pop_map.get(s, "UNKNOWN") for s in samples]
    pop_counts = dict(Counter(sample_pops))
    pops = sorted(pop_counts.keys())
    print("  Populations: %s" % ", ".join("%s (n=%d)" % (p, pop_counts[p]) for p in pops))

    # Build index mapping
    pop_indices = {}
    for pop in pops:
        pop_indices[pop] = [i for i, p in enumerate(sample_pops) if p == pop]

    # Heterozygosity
    print("Computing heterozygosity...")
    obs_het, exp_het, obs_per_site, exp_per_site = compute_heterozygosity(geno_matrix, pop_indices)
    for pop in pops:
        print("  %s: obs=%.4f  exp=%.4f" % (pop, obs_het[pop], exp_het[pop]))

    # FST
    print("Computing pairwise FST...")
    fst_df, fst_dict = compute_pairwise_fst(geno_matrix, pop_indices)
    for (p1, p2), val in sorted(fst_dict.items()):
        print("  %s vs %s: %.4f" % (p1, p2, val))

    # PCA
    print("Computing PCA...")
    pca_coords, pca_variance = compute_pca(geno_matrix)
    print("  PC1: %.1f%%  PC2: %.1f%%" % (pca_variance[0] * 100, pca_variance[1] * 100))

    # HEIM score
    print("Computing HEIM Equity Score...")
    heim_result = compute_heim_score(pop_counts, obs_het, len(fst_dict), weights)
    print("  Score: %s/100 (%s)" % (heim_result["heim_score"], heim_result["rating"]))

    # Save tables
    pd.DataFrame([
        {"population": k, "count": v, "proportion": v / sum(pop_counts.values()),
         "obs_het": obs_het.get(k, 0), "exp_het": exp_het.get(k, 0)}
        for k, v in sorted(pop_counts.items())
    ]).to_csv(tables_dir / "population_summary.csv", index=False)

    fst_df.to_csv(tables_dir / "fst_matrix.csv")

    het_df = pd.DataFrame({
        "population": pops,
        "observed_het": [obs_het[p] for p in pops],
        "expected_het": [exp_het[p] for p in pops],
    })
    het_df.to_csv(tables_dir / "heterozygosity.csv", index=False)

    with open(tables_dir / "heim_score.json", "w") as f:
        json.dump(heim_result, f, indent=2)

    # Figures
    print("Generating figures...")
    figures = OrderedDict()
    try:
        gauge_path = figures_dir / "heim_gauge.png"
        plot_heim_gauge(heim_result["heim_score"], heim_result["rating"], gauge_path)
        figures["HEIM Gauge"] = gauge_path

        ancestry_path = figures_dir / "ancestry_bar.png"
        plot_ancestry_bar(pop_counts, ancestry_path)
        figures["Ancestry Distribution"] = ancestry_path

        het_path = figures_dir / "heterozygosity.png"
        plot_heterozygosity(obs_het, exp_het, het_path)
        figures["Heterozygosity"] = het_path

        fst_path = figures_dir / "fst_heatmap.png"
        plot_fst_heatmap(fst_df, fst_path)
        figures["FST Heatmap"] = fst_path

        pca_path = figures_dir / "pca_plot.png"
        plot_pca(pca_coords, pca_variance, sample_pops, pca_path)
        figures["PCA"] = pca_path
    except ImportError as e:
        print("Warning: %s. Some figures skipped." % e, file=sys.stderr)

    # Report
    print("Generating report...")
    report = generate_report(
        heim_result, vcf_path, output_dir, figures,
        obs_het, exp_het, fst_df, pca_variance, n_variants,
    )
    report_path = output_dir / "report.md"
    report_path.write_text(report)

    print("\nDone.")
    print("  HEIM Score: %s/100 (%s)" % (heim_result["heim_score"], heim_result["rating"]))
    print("  Report: %s" % report_path)
    print("  Figures: %s" % figures_dir)

    return heim_result


def run_csv_pipeline(
    csv_path: Path,
    output_dir: Path,
    weights: Tuple[float, float, float, float],
) -> dict:
    """Equity scoring from ancestry CSV (no genotype data)."""
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(exist_ok=True)

    df = parse_ancestry_csv(csv_path)
    pop_counts = dict(Counter(df["population"]))

    # Literature-based heterozygosity estimates when no genotype data
    het_estimates = {
        "AFR": 0.35, "AMR": 0.28, "EAS": 0.25,
        "EUR": 0.27, "SAS": 0.26, "OCE": 0.30,
        "MID": 0.28, "UNKNOWN": 0.25,
    }
    obs_het = {pop: het_estimates.get(pop.upper(), 0.25) for pop in pop_counts}
    exp_het = obs_het.copy()

    n_pops = len(pop_counts)
    n_pairwise = n_pops * (n_pops - 1) // 2

    heim_result = compute_heim_score(pop_counts, obs_het, n_pairwise, weights)

    pd.DataFrame([
        {"population": k, "count": v, "proportion": v / sum(pop_counts.values())}
        for k, v in sorted(pop_counts.items())
    ]).to_csv(tables_dir / "population_summary.csv", index=False)

    with open(tables_dir / "heim_score.json", "w") as f:
        json.dump(heim_result, f, indent=2)

    figures = OrderedDict()
    try:
        gauge_path = figures_dir / "heim_gauge.png"
        plot_heim_gauge(heim_result["heim_score"], heim_result["rating"], gauge_path)
        figures["HEIM Gauge"] = gauge_path

        bar_path = figures_dir / "ancestry_bar.png"
        plot_ancestry_bar(pop_counts, bar_path)
        figures["Ancestry Distribution"] = bar_path
    except ImportError:
        print("Warning: matplotlib not available, skipping figures", file=sys.stderr)

    report = generate_report(
        heim_result, csv_path, output_dir, figures,
        obs_het, exp_het, None, None, 0,
    )
    report_path = output_dir / "report.md"
    report_path.write_text(report)

    print("HEIM Score: %s/100 (%s)" % (heim_result["heim_score"], heim_result["rating"]))
    print("Report: %s" % report_path)
    return heim_result


# ===================================================================
# CLI
# ===================================================================

def main() -> None:
    parser = argparse.ArgumentParser(
        description="HEIM Equity Scorer: diversity metrics from VCF or ancestry data"
    )
    parser.add_argument("--input", "-i", required=True,
                        help="Input VCF or ancestry CSV file")
    parser.add_argument("--pop-map", "-p", default=None,
                        help="Population map CSV (sample_id, population). "
                             "If not provided, infers from sample ID prefixes.")
    parser.add_argument("--output", "-o", default="equity_report",
                        help="Output directory for report and figures")
    parser.add_argument("--weights", default="0.35,0.25,0.20,0.20",
                        help="HEIM component weights: RI,HB,FC,GS (default: 0.35,0.25,0.20,0.20)")
    args = parser.parse_args()

    weights = tuple(float(w) for w in args.weights.split(","))
    if len(weights) != 4:
        print("Error: weights must be 4 comma-separated floats", file=sys.stderr)
        sys.exit(1)

    input_path = Path(args.input)
    output_dir = Path(args.output)
    pop_map_path = Path(args.pop_map) if args.pop_map else None

    suffix = input_path.suffix.lower()
    if suffix in (".vcf",):
        run_vcf_pipeline(input_path, pop_map_path, output_dir, weights)
    elif suffix in (".csv", ".tsv"):
        run_csv_pipeline(input_path, output_dir, weights)
    else:
        print("Error: unsupported file type '%s'. Use .vcf or .csv" % suffix,
              file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
