#!/usr/bin/env python3
"""
ClawBio Genome Comparator (ClawBio 2.0 port)
=============================================
Compare your genome to George Church (PGP-1) and estimate ancestry composition.

Usage:
    python genome_compare.py --input your_23andme.txt --output results/

The reference genome is a bundled subsample of George Church's public 23andMe
data from the Personal Genome Project (hu43860C, CC0 licensed). See
../SKILL.md for how this subsample was derived from the full 1.x file and
why (repo-size reasons, not a fidelity concession to the demo's numbers).

Port of ClawBio 1.x `genome-compare` (`genome_compare.py`) for ClawBio 2.0.
The core analysis logic (IBS computation, EM admixture ancestry estimation,
plotting, report + result.json generation) is unchanged from the 1.x source.

Adapted for the CLI surface only:
  - Dropped the `from clawbio.common.parsers import ...` / `.checksums` /
    `.report` imports and the `sys.path` PROJECT_ROOT hack (that module
    belongs to the 1.x repo layout). Below is a local, dependency-free
    reimplementation of just the pieces this script actually uses: a 23andMe
    parser (the only format `_parse_genotype_file` ever requested — format
    was hardcoded to "23andme" in the original, so the auto-detect /
    AncestryDNA / MyHeritage / VCF branches of the 1.x shared parser were
    never reachable from this script and are not ported), sha256 hashing,
    and the result.json envelope writer. Same pattern as
    skills/rnaseq_de/scripts/rnaseq_de.py and skills/claw_ancestry_pca.
  - Dropped `stage_from_icloud` (a macOS iCloud Drive file-locking
    workaround irrelevant to CLI-driven runs where the pipeline supplies
    ordinary local paths) — same precedent as claw_ancestry_pca, which
    never carried it forward either.
  - `--demo` is no longer a CLI flag. `pipeline.yaml` has no
    conditional/boolean-flag mechanism (see CLAUDE.md's "no Turing-complete
    YAML" constraint), so `command_template` cannot conditionally pass
    `--demo` depending on which file `{input}` resolves to. Per CLAUDE.md's
    guidance ("if a step needs branching logic, write a Python script and
    make that script the step body"), the branch moved into the script:
    `_is_demo_input()` below content-hashes the input file against the
    bundled `data/manuel_corpas_23andme.txt.gz` fixture and only then adds
    the "known ground truth" report section. Demo mode is now just
    `--input data/manuel_corpas_23andme.txt.gz` (see SKILL.md's Run
    section) — no separate flag needed, and real users' own files never get
    a spuriously-labeled ground-truth section.
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
from typing import Dict, Optional, Tuple

import numpy as np

# --------------------------------------------------------------------------- #
# Local dependency-free reimplementations of clawbio.common.* (1.x only)
# --------------------------------------------------------------------------- #


def open_genetic_file(filepath: str | Path):
    """Open a file, handling .gz transparently."""
    filepath = str(filepath)
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    return open(filepath, encoding="utf-8", errors="replace")


def parse_23andme(filepath: str | Path) -> Tuple[dict, dict]:
    """Parse a 23andMe raw-data file (.txt or .txt.gz).

    Returns:
        genotypes: {rsid: genotype_str}
        positions: {rsid: {"chrom": str, "pos": int}}
    """
    genotypes: Dict[str, str] = {}
    positions: Dict[str, dict] = {}
    with open_genetic_file(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rsid, chrom, pos_str, geno = parts[0], parts[1], parts[2], parts[3]
            if not (rsid.startswith("rs") or rsid.startswith("i")):
                continue
            geno = geno.replace("-", "")
            if len(geno) < 1:
                continue
            try:
                pos = int(pos_str)
            except ValueError:
                pos = 0
            genotypes[rsid] = geno
            positions[rsid] = {"chrom": chrom, "pos": pos}
    return genotypes, positions


def sha256_file(filepath: str | Path) -> str:
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_hex(filepath: str | Path, length: int = 16) -> str:
    return sha256_file(filepath)[:length]


def write_result_json(
    output_dir: str | Path,
    skill: str,
    version: str,
    summary: dict,
    data: dict,
    input_checksum: str = "",
) -> Path:
    """Write the standardized result.json envelope alongside report.md.

    Local re-implementation of ClawBio 1.x's
    `clawbio.common.report.write_result_json` (envelope-shape compatible;
    this skill only ever passes summary/data/input_checksum, so `datasets`,
    `status`, and `ok` — which 1.x's version supports as opt-in extras —
    are omitted here).
    """
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    envelope = {
        "skill": skill,
        "version": version,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
        "summary": summary,
        "data": data,
    }

    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(envelope, indent=2, default=str))
    return result_path


# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

SKILL_DIR = Path(__file__).resolve().parent.parent
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


def _parse_genotype_file(filepath: str | Path) -> Tuple[dict, dict]:
    return parse_23andme(filepath)


def _is_demo_input(input_path: Path) -> bool:
    """Detect whether the input file is the bundled Manuel Corpas demo fixture.

    Compared by content hash (not path equality) so it still matches if the
    CLI resolves/copies the path differently. See module docstring for why
    this replaces the 1.x `--demo` flag.
    """
    if not DEMO_PATIENT_FILE.exists() or not input_path.exists():
        return False
    try:
        return sha256_file(input_path) == sha256_file(DEMO_PATIENT_FILE)
    except OSError:
        return False


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
    shared = set(geno_a.keys()) & set(geno_b.keys())

    if not shared:
        return 0.0, 0, 0, {}

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

        chrom = "?"
        if positions_a and rsid in positions_a:
            chrom = positions_a[rsid]["chrom"]
        chrom_ibs_sum[chrom] += ibs
        chrom_count[chrom] += 1
        if ibs == 2:
            chrom_concordant[chrom] += 1

    n_overlap = len(shared)
    ibs_score = total_ibs / (2 * n_overlap)

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
    with open(panel_path) as f:
        data = json.load(f)
    populations = data["meta"]["populations"] if "meta" in data else data["populations"]
    markers = data["markers"]
    return markers, populations


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

    if complement.get(ref_u) == alt_u:
        return None

    count = sum(1 for a in alleles if a == alt_u)
    ref_count = sum(1 for a in alleles if a == ref_u)
    if count + ref_count == 2:
        return count

    alt_c = complement.get(alt_u, alt_u)
    ref_c = complement.get(ref_u, ref_u)
    count = sum(1 for a in alleles if a == alt_c)
    ref_count = sum(1 for a in alleles if a == ref_c)
    if count + ref_count == 2:
        return count

    return None


def estimate_ancestry(
    genotypes: dict,
    aims_panel: list,
    populations: list,
) -> dict:
    """Estimate ancestry composition using maximum likelihood + EM.

    For each AIM, computes P(genotype | population freq) using binomial
    genotype model, then runs EM to estimate admixture proportions.
    """
    pop_keys = ["afr", "eur", "eas", "sas", "amr"]
    pop_labels = populations
    K = len(pop_keys)

    geno_counts = []
    freq_matrix = []

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

    Q = np.ones(K) / K
    MAX_ITER = 200
    _em_converged = False

    for _em_iter in range(MAX_ITER):
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

        log_weighted = log_lik + np.log(Q)[np.newaxis, :]
        max_lw = log_weighted.max(axis=1, keepdims=True)
        exp_lw = np.exp(log_weighted - max_lw)
        responsibilities = exp_lw / exp_lw.sum(axis=1, keepdims=True)

        Q_new = responsibilities.mean(axis=0)
        Q_new = np.clip(Q_new, 1e-6, None)
        Q_new /= Q_new.sum()

        if np.allclose(Q, Q_new, atol=1e-7):
            _em_converged = True
            break
        Q = Q_new

    if not _em_converged:
        print(
            f"WARNING: EM admixture did not converge in {MAX_ITER} iterations; "
            "ancestry proportions may be unreliable.",
            file=sys.stderr,
        )

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
    ax1.bar(chroms, scores, color="#2196F3", alpha=0.8, label="IBS score")
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

    for vals, y in zip([user_vals, ref_vals], y_positions):
        left = 0
        for j, (pop, val) in enumerate(zip(pops, vals)):
            ax.barh(y, val, left=left, color=colors[j], edgecolor="white", height=0.5)
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
    _wedges, _texts, autotexts = ax.pie(
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

    for ref in IBS_REFERENCE:
        ax.axvline(ref["ibs"], color="#BDBDBD", linestyle="--", alpha=0.5)
        ax.text(
            ref["ibs"], 1.15, ref["relationship"],
            ha="center", va="bottom", fontsize=6, rotation=45,
        )

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
    input_hash = sha256_hex(str(input_path)) if input_path.exists() else "n/a"
    ref_hash = sha256_hex(str(reference_path)) if reference_path.exists() else "n/a"

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

    lines.append("## About the Reference Genome")
    lines.append("")
    lines.append(
        "**George Church** (hu43860C) is PGP-1 — the first participant in the "
        "[Personal Genome Project](https://pgp.med.harvard.edu/). He is a professor "
        "of genetics at Harvard Medical School and a pioneer of genomics. His full "
        "23andMe data (569,226 SNPs, build 36) is CC0 public domain; this ClawBio 2.0 "
        "port bundles a ~6,000-SNP subsample for repo-size reasons (see SKILL.md)."
    )
    lines.append("")

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

    lines.append("## Disclaimer")
    lines.append("")
    lines.append(
        "*ClawBio is a research and educational tool. It is not a medical device and "
        "does not provide clinical diagnoses. Consult a healthcare professional before "
        "making any medical decisions.*"
    )
    lines.append("")

    lines.append("## Reproducibility")
    lines.append("")
    lines.append(f"- Input checksum: `{input_hash}`")
    lines.append(f"- Reference checksum: `{ref_hash}`")
    lines.append(f"- Generated: {now}")
    lines.append("")

    report_text = "\n".join(lines)

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

    print(f"  [diag] SKILL_DIR = {SKILL_DIR}", file=sys.stderr)
    print(f"  [diag] DATA_DIR  = {DATA_DIR}", file=sys.stderr)
    print(f"  [diag] input     = {input_path} (exists={input_path.exists()})", file=sys.stderr)
    print(f"  [diag] reference = {reference_path} (exists={reference_path.exists()})", file=sys.stderr)
    print(f"  [diag] aims      = {aims_path} (exists={aims_path.exists()})", file=sys.stderr)
    print(f"  [diag] output    = {output_dir}", file=sys.stderr)
    print(f"  [diag] numpy     = {np.__version__}", file=sys.stderr)
    try:
        import matplotlib
        print(f"  [diag] matplotlib = {matplotlib.__version__}", file=sys.stderr)
    except ImportError:
        print(f"  [diag] matplotlib = NOT INSTALLED", file=sys.stderr)

    output_dir.mkdir(parents=True, exist_ok=True)
    fig_dir = output_dir / "figures"
    if not no_figures:
        fig_dir.mkdir(parents=True, exist_ok=True)

    print(f"  Parsing input: {input_path.name}")
    user_geno, user_pos = _parse_genotype_file(input_path)
    print(f"    {len(user_geno):,} SNPs loaded")

    print(f"  Parsing reference: {reference_path.name}")
    ref_geno, ref_pos = _parse_genotype_file(reference_path)
    print(f"    {len(ref_geno):,} SNPs loaded")

    print("  Computing IBS...")
    ibs_score, n_overlap, n_concordant, per_chrom = compute_ibs(
        user_geno, ref_geno, user_pos
    )
    print(f"    IBS = {ibs_score:.4f} across {n_overlap:,} shared SNPs")

    print("  Estimating ancestry...")
    aims_panel, populations = load_aims_panel(aims_path)
    user_ancestry = estimate_ancestry(user_geno, aims_panel, populations)
    print(f"    Used {user_ancestry['n_aims_used']}/{user_ancestry['n_aims_total']} AIMs")
    for pop, pct in sorted(
        user_ancestry["continental"].items(), key=lambda x: -x[1]
    ):
        print(f"    {pop}: {pct*100:.1f}%")

    manuel_ancestry = {}
    if MANUEL_ANCESTRY_FILE.exists():
        with open(MANUEL_ANCESTRY_FILE) as f:
            manuel_ancestry = json.load(f)

    figures = {}
    if not no_figures:
        try:
            plot_chromosome_ibs(per_chrom, fig_dir / "chromosome_ibs.png")
            figures["chromosome_ibs"] = "chromosome_ibs.png"

            plot_ancestry_pie(user_ancestry["continental"], fig_dir / "ancestry_pie.png")
            figures["ancestry_pie"] = "ancestry_pie.png"

            plot_ibs_context(ibs_score, fig_dir / "ibs_context.png")
            figures["ibs_context"] = "ibs_context.png"

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

    continental = user_ancestry.get("continental", {})
    top_pop = max(continental, key=continental.get) if continental else "Unknown"
    write_result_json(
        output_dir=output_dir,
        skill="genome_compare",
        version="0.2.0",
        summary={
            "ibs_score": ibs_score,
            "n_overlap": n_overlap,
            "n_concordant": n_concordant,
            "top_ancestry": top_pop,
        },
        data={
            "ibs_score": ibs_score,
            "per_chrom": per_chrom,
            "ancestry": user_ancestry,
        },
        input_checksum=sha256_hex(str(input_path)),
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
    parser.add_argument("--input", "-i", required=True, help="Path to your 23andMe raw-data file (.txt or .txt.gz)")
    parser.add_argument(
        "--reference", "-r", default=None,
        help="Reference genome file (default: bundled George Church PGP-1 subsample)"
    )
    parser.add_argument("--output", "-o", required=True, help="Output directory for report.md, result.json, figures/")
    parser.add_argument("--aims-panel", default=None, help="Path to AIMs panel JSON")
    parser.add_argument("--no-figures", action="store_true", help="Skip figure generation")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    output_dir = Path(args.output)
    reference = Path(args.reference) if args.reference else None
    aims = Path(args.aims_panel) if args.aims_panel else None

    run_comparison(
        input_path=input_path,
        output_dir=output_dir,
        reference_path=reference,
        aims_path=aims,
        no_figures=args.no_figures,
        is_demo=_is_demo_input(input_path),
    )

    print(f"\n  Done. Report at: {output_dir / 'report.md'}")
    sys.exit(0)


if __name__ == "__main__":
    import traceback
    try:
        main()
    except Exception:
        traceback.print_exc(file=sys.stderr)
        sys.exit(1)
