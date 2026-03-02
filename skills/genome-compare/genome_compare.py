#!/usr/bin/env python3
"""
ClawBio Genome Comparator
=========================
Compare your genome to George Church (PGP-1) and estimate ancestry composition.

Usage:
    python genome_compare.py --input your_23andme.txt --output results/
    python genome_compare.py --demo --output results/

The reference genome is George Church's public 23andMe data from the
Personal Genome Project (hu43860C, CC0 licensed).
"""

from __future__ import annotations

import argparse
import gzip
import hashlib
import json
import sys
from collections import defaultdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

SKILL_DIR = Path(__file__).resolve().parent
DATA_DIR = SKILL_DIR / "data"
REFERENCE_FILE = DATA_DIR / "george_church_23andme.txt.gz"
MANUEL_ANCESTRY_FILE = DATA_DIR / "manuel_ancestry.json"
AIMS_PANEL_FILE = DATA_DIR / "aims_panel.json"
DEMO_PATIENT_FILE = DATA_DIR / "manuel_corpas_23andme.txt.gz"

# --------------------------------------------------------------------------- #
# IBS reference values for context
# --------------------------------------------------------------------------- #

IBS_REFERENCE = [
    {"relationship": "Identical twins / same person", "ibs": 1.000},
    {"relationship": "Parent-child", "ibs": 0.850},
    {"relationship": "Full siblings", "ibs": 0.820},
    {"relationship": "Half-siblings", "ibs": 0.780},
    {"relationship": "First cousins", "ibs": 0.760},
    {"relationship": "Unrelated, same population", "ibs": 0.735},
    {"relationship": "European vs South Asian", "ibs": 0.720},
    {"relationship": "European vs East Asian", "ibs": 0.710},
    {"relationship": "European vs African", "ibs": 0.695},
    {"relationship": "African vs East Asian", "ibs": 0.680},
]

CHROMOSOMES = [str(c) for c in range(1, 23)] + ["X", "Y", "MT"]

# --------------------------------------------------------------------------- #
# Parsing
# --------------------------------------------------------------------------- #


def _open_file(filepath: str | Path):
    """Open a file, handling .gz transparently. Retries on iCloud Errno 11."""
    import time as _time
    filepath = str(filepath)
    for attempt in range(3):
        try:
            if filepath.endswith(".gz"):
                return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
            return open(filepath, encoding="utf-8", errors="replace")
        except OSError as e:
            if e.errno == 11 and attempt < 2:
                print(f"  iCloud lock (Errno 11), retrying ({attempt + 1}/3)...")
                _time.sleep(1 + attempt)
            else:
                raise


def parse_23andme_extended(filepath: str | Path) -> Tuple[dict, dict]:
    """Parse a 23andMe file returning genotypes AND chromosome/position metadata.

    Returns:
        genotypes: {rsid: genotype_str}
        positions: {rsid: {"chrom": str, "pos": int}}
    """
    genotypes = {}
    positions = {}
    with _open_file(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rsid, chrom, pos, geno = parts[0], parts[1], parts[2], parts[3]
            if not (rsid.startswith("rs") or rsid.startswith("i")):
                continue
            geno = geno.replace("-", "").replace("--", "")
            if len(geno) < 1:
                continue
            genotypes[rsid] = geno
            try:
                positions[rsid] = {"chrom": chrom, "pos": int(pos)}
            except ValueError:
                positions[rsid] = {"chrom": chrom, "pos": 0}
    return genotypes, positions


# --------------------------------------------------------------------------- #
# IBS computation
# --------------------------------------------------------------------------- #


def _ibs_at_site(geno_a: str, geno_b: str) -> int:
    """Compute IBS (0, 1, or 2) between two diploid genotype calls.

    Handles haploid calls on X/Y/MT (single character).
    """
    if len(geno_a) == 1 and len(geno_b) == 1:
        return 2 if geno_a == geno_b else 0
    if len(geno_a) < 2 or len(geno_b) < 2:
        return 0

    a_alleles = list(geno_a[:2])
    b_alleles = list(geno_b[:2])

    # Multiset matching: count shared alleles
    remaining = list(a_alleles)
    match_count = 0
    for b in b_alleles:
        if b in remaining:
            remaining.remove(b)
            match_count += 1
    return match_count


def compute_ibs(
    geno_a: dict,
    geno_b: dict,
    positions_a: Optional[dict] = None,
) -> Tuple[float, int, int, dict]:
    """Compute Identity By State between two genotype dictionaries.

    Returns:
        ibs_score: float in [0, 1]
        n_overlap: number of shared rsIDs
        n_concordant: number of sites with IBS=2
        per_chrom: {chrom: {"ibs": float, "n_overlap": int, "n_concordant": int}}
    """
    # Find overlapping rsIDs
    shared = set(geno_a.keys()) & set(geno_b.keys())

    if not shared:
        return 0.0, 0, 0, {}

    # Per-chromosome accumulators
    chrom_ibs_sum: Dict[str, int] = defaultdict(int)
    chrom_count: Dict[str, int] = defaultdict(int)
    chrom_concordant: Dict[str, int] = defaultdict(int)

    total_ibs = 0
    n_concordant = 0

    for rsid in shared:
        ibs = _ibs_at_site(geno_a[rsid], geno_b[rsid])
        total_ibs += ibs
        if ibs == 2:
            n_concordant += 1

        # Assign to chromosome
        chrom = "?"
        if positions_a and rsid in positions_a:
            chrom = positions_a[rsid]["chrom"]
        chrom_ibs_sum[chrom] += ibs
        chrom_count[chrom] += 1
        if ibs == 2:
            chrom_concordant[chrom] += 1

    n_overlap = len(shared)
    ibs_score = total_ibs / (2 * n_overlap)

    # Per-chromosome breakdown
    per_chrom = {}
    for chrom in CHROMOSOMES:
        if chrom in chrom_count and chrom_count[chrom] > 0:
            per_chrom[chrom] = {
                "ibs": chrom_ibs_sum[chrom] / (2 * chrom_count[chrom]),
                "n_overlap": chrom_count[chrom],
                "n_concordant": chrom_concordant[chrom],
            }

    return ibs_score, n_overlap, n_concordant, per_chrom


# --------------------------------------------------------------------------- #
# Ancestry estimation
# --------------------------------------------------------------------------- #


def load_aims_panel(panel_path: Path) -> Tuple[list, list]:
    """Load AIMs panel. Returns (markers, population_names)."""
    import time as _time
    for attempt in range(3):
        try:
            with open(panel_path) as f:
                data = json.load(f)
            populations = data["meta"]["populations"]
            markers = data["markers"]
            return markers, populations
        except OSError as e:
            if e.errno == 11 and attempt < 2:
                print(f"  iCloud lock reading AIMs panel (Errno 11), retrying ({attempt + 1}/3)...")
                _time.sleep(1 + attempt)
            else:
                raise


def _count_alt_alleles(geno: str, ref: str, alt: str) -> Optional[int]:
    """Count alt alleles in a diploid genotype given ref/alt annotation.

    Returns 0 (hom ref), 1 (het), 2 (hom alt), or None if ambiguous.
    Skips ambiguous complement-pair SNPs (A/T, C/G) to avoid strand errors.
    """
    if len(geno) < 2:
        return None

    complement = {"A": "T", "T": "A", "C": "G", "G": "C"}
    alleles = [geno[0].upper(), geno[1].upper()]
    ref_u = ref.upper()
    alt_u = alt.upper()

    # Skip complement-pair SNPs (A/T or C/G) — strand can't be resolved
    if complement.get(ref_u) == alt_u:
        return None

    # Try direct match
    count = sum(1 for a in alleles if a == alt_u)
    ref_count = sum(1 for a in alleles if a == ref_u)
    if count + ref_count == 2:
        return count

    # Try complement (23andMe may report on minus strand)
    alt_c = complement.get(alt_u, alt_u)
    ref_c = complement.get(ref_u, ref_u)
    count = sum(1 for a in alleles if a == alt_c)
    ref_count = sum(1 for a in alleles if a == ref_c)
    if count + ref_count == 2:
        return count

    return None  # Can't resolve


def estimate_ancestry(
    genotypes: dict,
    aims_panel: list,
    populations: list,
) -> dict:
    """Estimate ancestry composition using maximum likelihood + EM.

    For each AIM, computes P(genotype | population freq) using binomial
    genotype model, then runs EM to estimate admixture proportions.

    Returns:
        {
            "continental": {"EUR": 0.977, "AFR": 0.003, ...},
            "n_aims_used": int,
            "n_aims_total": int,
            "method": "EM admixture (binomial genotype likelihood)",
        }
    """
    pop_keys = ["afr", "eur", "eas", "sas", "amr"]
    pop_labels = populations  # ["AFR", "EUR", "EAS", "SAS", "AMR"]
    K = len(pop_keys)

    # Build genotype and frequency matrices using ref/alt allele annotations
    geno_counts = []  # alt allele counts (0, 1, 2)
    freq_matrix = []  # (n_aims, K) alt allele frequencies

    for marker in aims_panel:
        rsid = marker["rsid"]
        if rsid not in genotypes:
            continue
        geno = genotypes[rsid]
        ref = marker.get("ref", "")
        alt = marker.get("alt", "")
        if not ref or not alt:
            continue

        count = _count_alt_alleles(geno, ref, alt)
        if count is None:
            continue

        freqs = [marker.get(k, 0.5) for k in pop_keys]
        geno_counts.append(count)
        freq_matrix.append(freqs)

    n_aims = len(geno_counts)
    if n_aims == 0:
        return {
            "continental": {p: 1.0 / K for p in pop_labels},
            "n_aims_used": 0,
            "n_aims_total": len(aims_panel),
            "method": "uniform (no AIMs found in genotype data)",
        }

    geno_arr = np.array(geno_counts)
    freq_arr = np.clip(np.array(freq_matrix), 0.001, 0.999)

    # EM admixture
    Q = np.ones(K) / K  # initial mixing proportions

    for _ in range(200):
        # E-step: compute P(genotype | pop) for each AIM and pop
        # P(g=0|p) = (1-p)^2, P(g=1|p) = 2p(1-p), P(g=2|p) = p^2
        log_lik = np.zeros((n_aims, K))
        for k in range(K):
            p = freq_arr[:, k]
            for g_val in [0, 1, 2]:
                mask = geno_arr == g_val
                if g_val == 0:
                    log_lik[mask, k] = 2 * np.log(1 - p[mask])
                elif g_val == 1:
                    log_lik[mask, k] = np.log(2) + np.log(p[mask]) + np.log(1 - p[mask])
                else:
                    log_lik[mask, k] = 2 * np.log(p[mask])

        # Weighted by Q
        log_weighted = log_lik + np.log(Q)[np.newaxis, :]
        # Normalise per-AIM (log-sum-exp for stability)
        max_lw = log_weighted.max(axis=1, keepdims=True)
        exp_lw = np.exp(log_weighted - max_lw)
        responsibilities = exp_lw / exp_lw.sum(axis=1, keepdims=True)

        # M-step
        Q_new = responsibilities.mean(axis=0)
        Q_new = np.clip(Q_new, 1e-6, None)
        Q_new /= Q_new.sum()

        if np.allclose(Q, Q_new, atol=1e-7):
            break
        Q = Q_new

    continental = {pop_labels[k]: round(float(Q[k]), 4) for k in range(K)}

    return {
        "continental": continental,
        "n_aims_used": n_aims,
        "n_aims_total": len(aims_panel),
        "method": "EM admixture (binomial genotype likelihood)",
    }


# --------------------------------------------------------------------------- #
# Visualisation
# --------------------------------------------------------------------------- #


def plot_chromosome_ibs(per_chrom: dict, output_path: Path) -> None:
    """Bar chart of IBS score per chromosome."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    chroms = [c for c in CHROMOSOMES if c in per_chrom]
    scores = [per_chrom[c]["ibs"] for c in chroms]
    overlaps = [per_chrom[c]["n_overlap"] for c in chroms]

    fig, ax1 = plt.subplots(figsize=(14, 5))
    bars = ax1.bar(chroms, scores, color="#2196F3", alpha=0.8, label="IBS score")
    ax1.set_xlabel("Chromosome")
    ax1.set_ylabel("IBS Score", color="#2196F3")
    ax1.set_ylim(0.5, 1.0)
    ax1.tick_params(axis="y", labelcolor="#2196F3")

    ax2 = ax1.twinx()
    ax2.plot(chroms, overlaps, "o-", color="#FF5722", alpha=0.7, label="Overlap count")
    ax2.set_ylabel("SNPs overlapping", color="#FF5722")
    ax2.tick_params(axis="y", labelcolor="#FF5722")

    ax1.set_title("Identity By State — Per Chromosome")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_ancestry_composition(
    user_ancestry: dict,
    reference_name: str,
    reference_ancestry: dict,
    output_path: Path,
) -> None:
    """Side-by-side stacked horizontal bar: user vs reference ancestry."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    pops = list(user_ancestry.keys())
    user_vals = [user_ancestry[p] for p in pops]
    ref_vals = [reference_ancestry.get(p, 0) for p in pops]

    colors = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336", "#607D8B"]
    while len(colors) < len(pops):
        colors.append("#888888")

    fig, ax = plt.subplots(figsize=(10, 4))

    y_positions = [1, 0]
    labels = ["You", reference_name]

    for i, (vals, y) in enumerate(zip([user_vals, ref_vals], y_positions)):
        left = 0
        for j, (pop, val) in enumerate(zip(pops, vals)):
            bar = ax.barh(y, val, left=left, color=colors[j], edgecolor="white", height=0.5)
            if val > 0.05:
                ax.text(
                    left + val / 2, y, f"{pop}\n{val:.1%}",
                    ha="center", va="center", fontsize=7, fontweight="bold", color="white",
                )
            left += val

    ax.set_yticks(y_positions)
    ax.set_yticklabels(labels, fontsize=12)
    ax.set_xlim(0, 1.0)
    ax.set_xlabel("Proportion")
    ax.set_title("Ancestry Composition — You vs Reference")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_ancestry_pie(ancestry: dict, output_path: Path) -> None:
    """Pie chart of ancestry composition."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    pops = [p for p, v in ancestry.items() if v > 0.005]
    vals = [ancestry[p] for p in pops]
    colors = ["#2196F3", "#4CAF50", "#FF9800", "#9C27B0", "#F44336", "#607D8B"]
    while len(colors) < len(pops):
        colors.append("#888888")

    fig, ax = plt.subplots(figsize=(7, 7))
    wedges, texts, autotexts = ax.pie(
        vals, labels=pops, autopct="%1.1f%%", colors=colors[:len(pops)],
        pctdistance=0.8, startangle=90,
    )
    for t in autotexts:
        t.set_fontsize(9)
    ax.set_title("Your Estimated Ancestry Composition")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


def plot_ibs_context(ibs_score: float, output_path: Path) -> None:
    """Horizontal gauge showing where IBS falls on relationship spectrum."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig, ax = plt.subplots(figsize=(12, 3))

    # Draw reference markers
    for ref in IBS_REFERENCE:
        ax.axvline(ref["ibs"], color="#BDBDBD", linestyle="--", alpha=0.5)
        ax.text(
            ref["ibs"], 1.15, ref["relationship"],
            ha="center", va="bottom", fontsize=6, rotation=45,
        )

    # Draw the user's score
    ax.axvline(ibs_score, color="#F44336", linewidth=3, label=f"Your IBS: {ibs_score:.4f}")
    ax.scatter([ibs_score], [0.5], color="#F44336", s=200, zorder=5)

    ax.set_xlim(0.65, 1.02)
    ax.set_ylim(0, 1.5)
    ax.set_yticks([])
    ax.set_xlabel("IBS Score")
    ax.set_title("Where Your Genetic Similarity Falls")
    ax.legend(loc="lower right", fontsize=10)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150, bbox_inches="tight")
    plt.close(fig)


# --------------------------------------------------------------------------- #
# Report generation
# --------------------------------------------------------------------------- #


def _sha256(filepath: Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()[:16]


def _interpret_ibs(score: float) -> str:
    """Return a human-readable interpretation of the IBS score."""
    for i, ref in enumerate(IBS_REFERENCE):
        if score >= ref["ibs"]:
            if i == 0:
                return f"Your IBS score ({score:.4f}) is consistent with **{ref['relationship'].lower()}**."
            prev = IBS_REFERENCE[i - 1]
            return (
                f"Your IBS score ({score:.4f}) falls between "
                f"**{ref['relationship'].lower()}** ({ref['ibs']:.3f}) and "
                f"**{prev['relationship'].lower()}** ({prev['ibs']:.3f})."
            )
    last = IBS_REFERENCE[-1]
    return (
        f"Your IBS score ({score:.4f}) is below typical "
        f"**{last['relationship'].lower()}** ({last['ibs']:.3f}), "
        f"suggesting distant or cross-continental ancestry."
    )


def generate_report(
    input_path: Path,
    reference_path: Path,
    ibs_score: float,
    n_overlap: int,
    n_concordant: int,
    per_chrom: dict,
    user_ancestry: dict,
    manuel_ancestry: dict,
    figures: dict,
    output_dir: Path,
    is_demo: bool,
) -> str:
    """Generate the genome comparison markdown report."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    input_hash = _sha256(input_path) if input_path.exists() else "n/a"
    ref_hash = _sha256(reference_path) if reference_path.exists() else "n/a"

    # Ancestry flat for display
    continental = user_ancestry.get("continental", {})
    top_pop = max(continental, key=continental.get) if continental else "Unknown"
    top_pct = continental.get(top_pop, 0) * 100

    lines = []
    lines.append("# Genome Comparison Report")
    lines.append("")
    lines.append(f"**Date**: {now}  ")
    lines.append(f"**Reference**: George Church (PGP-1, hu43860C)  ")
    if is_demo:
        lines.append("**Input**: Manuel Corpas ([PGP-UK uk6D0CFA](https://my.personalgenomes.org.uk/profile/uk6D0CFA), Corpasome, demo mode)  ")
    else:
        lines.append(f"**Input**: {input_path.name}  ")
    lines.append("")

    # Summary box
    lines.append("## Summary")
    lines.append("")
    lines.append(f"| Metric | Value |")
    lines.append(f"|--------|-------|")
    lines.append(f"| **IBS Score** | {ibs_score:.4f} |")
    lines.append(f"| **SNPs Overlapping** | {n_overlap:,} |")
    lines.append(f"| **Genotypes Identical (IBS=2)** | {n_concordant:,} ({n_concordant/n_overlap*100:.1f}%) |")
    lines.append(f"| **Top Ancestry** | {top_pop} ({top_pct:.1f}%) |")
    lines.append(f"| **AIMs Used** | {user_ancestry.get('n_aims_used', 0)} / {user_ancestry.get('n_aims_total', 0)} |")
    lines.append("")

    # IBS interpretation
    lines.append("## Identity By State Analysis")
    lines.append("")
    lines.append(_interpret_ibs(ibs_score))
    lines.append("")
    lines.append("### Population Context")
    lines.append("")
    lines.append("| Relationship | Typical IBS |")
    lines.append("|---|---|")
    for ref in IBS_REFERENCE:
        marker = " **<--**" if abs(ref["ibs"] - ibs_score) < 0.015 else ""
        lines.append(f"| {ref['relationship']} | {ref['ibs']:.3f}{marker} |")
    lines.append("")

    if "ibs_context" in figures:
        lines.append(f"![IBS Context](figures/{figures['ibs_context']})")
        lines.append("")

    # Chromosome breakdown
    lines.append("## Chromosome Breakdown")
    lines.append("")
    lines.append("| Chr | SNPs | IBS Score | Identical (IBS=2) |")
    lines.append("|-----|------|-----------|-------------------|")
    for chrom in CHROMOSOMES:
        if chrom in per_chrom:
            c = per_chrom[chrom]
            lines.append(
                f"| {chrom} | {c['n_overlap']:,} | {c['ibs']:.4f} | {c['n_concordant']:,} |"
            )
    lines.append("")

    if "chromosome_ibs" in figures:
        lines.append(f"![Chromosome IBS](figures/{figures['chromosome_ibs']})")
        lines.append("")

    # Ancestry composition
    lines.append("## Ancestry Composition")
    lines.append("")
    lines.append("### Your Estimated Ancestry")
    lines.append("")
    lines.append("| Population | Proportion |")
    lines.append("|------------|-----------|")
    for pop in sorted(continental, key=continental.get, reverse=True):
        pct = continental[pop] * 100
        lines.append(f"| {pop} | {pct:.1f}% |")
    lines.append("")
    lines.append(f"*Method: {user_ancestry.get('method', 'EM admixture')}*")
    lines.append("")

    if "ancestry_pie" in figures:
        lines.append(f"![Ancestry Pie](figures/{figures['ancestry_pie']})")
        lines.append("")

    # In demo mode, show known ground truth for validation
    if is_demo and manuel_ancestry:
        mc = manuel_ancestry.get("ancestry_composition", {})
        lines.append("### Ground Truth — Manuel Corpas (23andMe Ancestry Report)")
        lines.append("")
        lines.append(
            "*Since this demo uses Manuel Corpas as the input, we can compare "
            "the estimation above against his actual 23andMe ancestry report:*"
        )
        lines.append("")
        lines.append("| Region | Proportion |")
        lines.append("|--------|-----------|")
        for region, data in mc.items():
            if isinstance(data, dict):
                total = data.get("total", 0)
                lines.append(f"| **{region}** | **{total*100:.1f}%** |")
                for sub, val in data.items():
                    if sub == "total":
                        continue
                    if isinstance(val, dict):
                        lines.append(f"| &emsp;{sub} | {val.get('total', 0)*100:.1f}% |")
                        for subsub, vv in val.items():
                            if subsub == "total":
                                continue
                            lines.append(f"| &emsp;&emsp;{subsub} | {vv*100:.1f}% |")
                    else:
                        lines.append(f"| &emsp;{sub} | {val*100:.1f}% |")
            else:
                lines.append(f"| **{region}** | **{data*100:.1f}%** |")
        lines.append("")

        hg = manuel_ancestry.get("haplogroups", {})
        if hg:
            lines.append(f"**Maternal haplogroup**: {hg.get('maternal', '?')}  ")
            lines.append(f"**Paternal haplogroup**: {hg.get('paternal', '?')}  ")
            lines.append("")

        neanderthal = manuel_ancestry.get("neanderthal", {})
        if neanderthal:
            lines.append(f"**Neanderthal**: {neanderthal.get('description', '')}  ")
            lines.append("")

    if "ancestry_composition" in figures:
        lines.append(f"![Ancestry Comparison](figures/{figures['ancestry_composition']})")
        lines.append("")

    # George Church info
    lines.append("## About the Reference Genome")
    lines.append("")
    lines.append(
        "**George Church** (hu43860C) is PGP-1 — the first participant in the "
        "[Personal Genome Project](https://pgp.med.harvard.edu/). He is a professor "
        "of genetics at Harvard Medical School and a pioneer of genomics. His 23andMe "
        "data (569,226 SNPs, build 36) is CC0 public domain."
    )
    lines.append("")

    # Methods
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "**IBS (Identity By State)**: For each SNP present in both files, alleles are "
        "compared. IBS per site = number of shared alleles (0, 1, or 2). Overall IBS = "
        "sum(per-site IBS) / (2 * N). Range: 0.0 (no alleles shared) to 1.0 (identical)."
    )
    lines.append("")
    lines.append(
        "**Ancestry estimation**: EM admixture algorithm using ancestry-informative "
        "markers (AIMs) with population allele frequencies from 1000 Genomes Phase 3. "
        "Binomial genotype likelihood model with iterative expectation-maximisation. "
        "This is an approximation — for clinical-grade ancestry, use ADMIXTURE or RFMix."
    )
    lines.append("")
    lines.append(
        "**Reference IBS values** are approximate ranges from population genetics "
        "literature (1000 Genomes, HGDP). Actual values depend on SNP panel and "
        "population composition."
    )
    lines.append("")

    # Disclaimer
    lines.append("## Disclaimer")
    lines.append("")
    lines.append(
        "*ClawBio is a research and educational tool. It is not a medical device and "
        "does not provide clinical diagnoses. Consult a healthcare professional before "
        "making any medical decisions.*"
    )
    lines.append("")

    # Reproducibility
    lines.append("## Reproducibility")
    lines.append("")
    lines.append(f"- Input checksum: `{input_hash}`")
    lines.append(f"- Reference checksum: `{ref_hash}`")
    lines.append(f"- Generated: {now}")
    lines.append("")

    report_text = "\n".join(lines)

    # Write report
    report_path = output_dir / "report.md"
    report_path.write_text(report_text, encoding="utf-8")

    return report_text


# --------------------------------------------------------------------------- #
# Main pipeline
# --------------------------------------------------------------------------- #


def run_comparison(
    input_path: Path,
    output_dir: Path,
    reference_path: Optional[Path] = None,
    aims_path: Optional[Path] = None,
    no_figures: bool = False,
    is_demo: bool = False,
) -> dict:
    """Run full genome comparison pipeline."""
    if reference_path is None:
        reference_path = REFERENCE_FILE
    if aims_path is None:
        aims_path = AIMS_PANEL_FILE

    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    if not no_figures:
        fig_dir.mkdir(parents=True, exist_ok=True)

    # Parse both genomes
    print(f"  Parsing input: {input_path.name}")
    user_geno, user_pos = parse_23andme_extended(input_path)
    print(f"    {len(user_geno):,} SNPs loaded")

    print(f"  Parsing reference: {reference_path.name}")
    ref_geno, ref_pos = parse_23andme_extended(reference_path)
    print(f"    {len(ref_geno):,} SNPs loaded")

    # Compute IBS
    print("  Computing IBS...")
    ibs_score, n_overlap, n_concordant, per_chrom = compute_ibs(
        user_geno, ref_geno, user_pos
    )
    print(f"    IBS = {ibs_score:.4f} across {n_overlap:,} shared SNPs")

    # Estimate ancestry
    print("  Estimating ancestry...")
    aims_panel, populations = load_aims_panel(aims_path)
    user_ancestry = estimate_ancestry(user_geno, aims_panel, populations)
    print(f"    Used {user_ancestry['n_aims_used']}/{user_ancestry['n_aims_total']} AIMs")
    for pop, pct in sorted(
        user_ancestry["continental"].items(), key=lambda x: -x[1]
    ):
        print(f"    {pop}: {pct*100:.1f}%")

    # Load Manuel's known ancestry
    manuel_ancestry = {}
    if MANUEL_ANCESTRY_FILE.exists():
        import time as _time
        for _attempt in range(3):
            try:
                with open(MANUEL_ANCESTRY_FILE) as f:
                    manuel_ancestry = json.load(f)
                break
            except OSError as e:
                if e.errno == 11 and _attempt < 2:
                    print(f"  iCloud lock reading ancestry file (Errno 11), retrying ({_attempt + 1}/3)...")
                    _time.sleep(1 + _attempt)
                else:
                    raise

    # Generate figures
    figures = {}
    if not no_figures:
        try:
            plot_chromosome_ibs(per_chrom, fig_dir / "chromosome_ibs.png")
            figures["chromosome_ibs"] = "chromosome_ibs.png"

            plot_ancestry_pie(user_ancestry["continental"], fig_dir / "ancestry_pie.png")
            figures["ancestry_pie"] = "ancestry_pie.png"

            plot_ibs_context(ibs_score, fig_dir / "ibs_context.png")
            figures["ibs_context"] = "ibs_context.png"

            # Side-by-side: user vs Manuel's continental
            manuel_continental = {}
            mc = manuel_ancestry.get("ancestry_composition", {})
            pop_map = {
                "European": "EUR",
                "Western Asian & North African": "AMR",
                "Sub-Saharan African": "AFR",
                "East Asian & Native American": "EAS",
                "Unassigned": "SAS",
            }
            for region, data in mc.items():
                pop_key = pop_map.get(region, region)
                if isinstance(data, dict):
                    manuel_continental[pop_key] = data.get("total", 0)
                else:
                    manuel_continental[pop_key] = data

            plot_ancestry_composition(
                user_ancestry["continental"],
                "Manuel Corpas",
                manuel_continental,
                fig_dir / "ancestry_comparison.png",
            )
            figures["ancestry_composition"] = "ancestry_comparison.png"
        except Exception as e:
            print(f"  Warning: figure generation failed: {e}")

    # Generate report
    print("  Generating report...")
    generate_report(
        input_path=input_path,
        reference_path=reference_path,
        ibs_score=ibs_score,
        n_overlap=n_overlap,
        n_concordant=n_concordant,
        per_chrom=per_chrom,
        user_ancestry=user_ancestry,
        manuel_ancestry=manuel_ancestry,
        figures=figures,
        output_dir=output_dir,
        is_demo=is_demo,
    )

    return {
        "ibs_score": ibs_score,
        "n_overlap": n_overlap,
        "n_concordant": n_concordant,
        "per_chrom": per_chrom,
        "ancestry": user_ancestry,
        "figures": figures,
    }


# --------------------------------------------------------------------------- #
# CLI
# --------------------------------------------------------------------------- #


def main():
    parser = argparse.ArgumentParser(
        description="ClawBio Genome Comparator: compare your genome to George Church (PGP-1)"
    )
    parser.add_argument("--input", "-i", help="Path to your 23andMe/AncestryDNA file")
    parser.add_argument(
        "--reference", "-r", default=None,
        help="Reference genome file (default: George Church PGP-1)"
    )
    parser.add_argument("--output", "-o", default="genome_compare_results", help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with Manuel Corpas as input (demo)")
    parser.add_argument("--aims-panel", default=None, help="Path to AIMs panel JSON")
    parser.add_argument("--no-figures", action="store_true", help="Skip figure generation")
    args = parser.parse_args()

    if args.demo:
        input_path = DEMO_PATIENT_FILE
    elif args.input:
        input_path = Path(args.input)
    else:
        parser.error("Provide --input or --demo")
        return

    if not Path(input_path).exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    output_dir = Path(args.output)
    reference = Path(args.reference) if args.reference else None
    aims = Path(args.aims_panel) if args.aims_panel else None

    result = run_comparison(
        input_path=Path(input_path),
        output_dir=output_dir,
        reference_path=reference,
        aims_path=aims,
        no_figures=args.no_figures,
        is_demo=args.demo,
    )

    print(f"\n  Done. Report at: {output_dir / 'report.md'}")
    sys.exit(0)


if __name__ == "__main__":
    main()
