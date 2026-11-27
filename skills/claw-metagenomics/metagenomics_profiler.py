#!/usr/bin/env python3
"""
ClawBio Metagenomics Profiler
Shotgun metagenomics pipeline: Kraken2/Bracken taxonomy + RGI/CARD resistome
+ HUMAnN3 functional profiling.

Usage:
    python metagenomics_profiler.py --r1 R1.fastq.gz --r2 R2.fastq.gz --output report/
    python metagenomics_profiler.py --demo --output demo_report/
"""

from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)
from clawbio.common.textio import write_text_lf_atomic  # noqa: E402

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

SCRIPT_DIR = Path(__file__).resolve().parent

KRAKEN2_CONFIDENCE = 0.2
BRACKEN_THRESHOLD = 10       # minimum reads for Bracken re-estimation
BRACKEN_LEVEL = "S"          # species
DEFAULT_READ_LENGTH = 150
RGI_CRITERIA = ("Perfect", "Strict")

SAFETY_DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions. Antimicrobial resistance findings "
    "must be confirmed by culture-based susceptibility testing."
)

# WHO Critical Priority Pathogens and their associated ARG families / drug classes
WHO_PRIORITY_ARGS: Dict[str, Dict[str, Any]] = {
    "Critical": {
        "pathogens": [
            "Acinetobacter baumannii",
            "Pseudomonas aeruginosa",
            "Enterobacteriaceae",
            "Escherichia coli",
            "Klebsiella pneumoniae",
        ],
        "drug_classes": ["carbapenem", "cephalosporin"],
        "arg_families": [
            "NDM", "OXA-48", "KPC", "VIM", "IMP",
            "CTX-M", "SHV", "TEM",
        ],
    },
    "High": {
        "pathogens": [
            "Enterococcus faecium",
            "Staphylococcus aureus",
            "Helicobacter pylori",
            "Campylobacter",
            "Salmonella",
            "Neisseria gonorrhoeae",
        ],
        "drug_classes": [
            "vancomycin", "methicillin", "clarithromycin",
            "fluoroquinolone",
        ],
        "arg_families": [
            "VanA", "VanB", "mecA", "mecC",
            "QnrA", "QnrB", "QnrS", "GyrA", "ParC",
        ],
    },
    "Medium": {
        "pathogens": [
            "Streptococcus pneumoniae",
            "Haemophilus influenzae",
            "Shigella",
        ],
        "drug_classes": ["penicillin", "ampicillin"],
        "arg_families": ["PBP", "TEM", "ROB"],
    },
}

# Colourblind-friendly palette for taxonomy plots
TAXON_COLOURS = [
    "#4C72B0", "#DD8452", "#55A868", "#C44E52", "#8172B3",
    "#937860", "#DA8BC3", "#8C8C8C", "#CCB974", "#64B5CD",
    "#E69F00", "#56B4E9", "#009E73", "#F0E442", "#0072B2",
    "#D55E00", "#CC79A7", "#999999", "#BFBFBF", "#A6761D",
]

# ---------------------------------------------------------------------------
# Demo data: pre-computed Peru sewage metagenomics (6 samples, 3 sites)
# ---------------------------------------------------------------------------

DEMO_TAXONOMY = {
    "species": [
        "Escherichia coli", "Klebsiella pneumoniae", "Pseudomonas aeruginosa",
        "Acinetobacter baumannii", "Enterococcus faecium", "Bacteroides fragilis",
        "Salmonella enterica", "Staphylococcus aureus", "Clostridioides difficile",
        "Campylobacter jejuni", "Prevotella copri", "Faecalibacterium prausnitzii",
        "Bifidobacterium longum", "Lactobacillus rhamnosus", "Helicobacter pylori",
        "Vibrio cholerae", "Shigella flexneri", "Proteus mirabilis",
        "Serratia marcescens", "Citrobacter freundii",
    ],
    "Lima": [12.3, 8.7, 5.1, 3.9, 2.8, 4.5, 1.9, 1.2, 0.8, 0.6,
             3.2, 5.1, 2.1, 1.0, 0.4, 0.3, 1.1, 0.9, 0.7, 0.5],
    "Cusco": [8.1, 5.2, 3.8, 2.1, 1.5, 6.2, 0.8, 0.5, 1.2, 0.9,
              5.8, 7.3, 3.4, 2.1, 0.7, 0.1, 0.4, 1.3, 0.6, 0.8],
    "Iquitos": [15.6, 11.3, 7.2, 5.8, 3.9, 3.1, 3.2, 2.1, 0.5, 1.4,
                2.1, 3.8, 1.5, 0.8, 1.1, 1.8, 2.3, 1.1, 1.2, 0.9],
}

DEMO_RESISTOME = {
    "gene": [
        "NDM-1", "OXA-48", "KPC-3", "CTX-M-15", "CTX-M-27",
        "mecA", "VanA", "VanB", "QnrB", "QnrS",
        "TEM-1", "SHV-12", "GyrA_S83L", "ParC_S80I", "aac(6')-Ib",
        "aph(3'')-Ib", "sul1", "sul2", "dfrA12", "tet(A)",
        "erm(B)", "mph(A)", "catA1", "floR",
    ],
    "drug_class": [
        "carbapenem", "carbapenem", "carbapenem", "cephalosporin", "cephalosporin",
        "methicillin", "vancomycin", "vancomycin", "fluoroquinolone", "fluoroquinolone",
        "penicillin", "cephalosporin", "fluoroquinolone", "fluoroquinolone", "aminoglycoside",
        "aminoglycoside", "sulfonamide", "sulfonamide", "trimethoprim", "tetracycline",
        "macrolide", "macrolide", "phenicol", "phenicol",
    ],
    "criteria": [
        "Perfect", "Strict", "Strict", "Perfect", "Strict",
        "Strict", "Perfect", "Strict", "Strict", "Strict",
        "Perfect", "Strict", "Strict", "Strict", "Strict",
        "Strict", "Perfect", "Perfect", "Strict", "Perfect",
        "Strict", "Strict", "Perfect", "Strict",
    ],
    "Lima": [
        18, 12, 8, 45, 22,
        5, 3, 1, 15, 9,
        52, 28, 11, 7, 19,
        14, 38, 29, 16, 33,
        12, 8, 6, 4,
    ],
    "Cusco": [
        4, 2, 1, 18, 8,
        2, 1, 0, 6, 3,
        31, 12, 4, 2, 8,
        6, 22, 15, 9, 19,
        7, 3, 2, 1,
    ],
    "Iquitos": [
        25, 18, 14, 58, 31,
        8, 5, 3, 22, 14,
        67, 35, 18, 12, 27,
        21, 51, 38, 23, 45,
        18, 12, 9, 7,
    ],
}

DEMO_PATHWAYS = {
    "pathway": [
        "PWY-7219: adenosine ribonucleotides de novo biosynthesis",
        "PWY-7221: guanosine ribonucleotides de novo biosynthesis",
        "PWY-6123: inosine-5'-phosphate biosynthesis I",
        "GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)",
        "PWY-5100: pyruvate fermentation to acetate and lactate II",
        "TCA: TCA cycle I (prokaryotic)",
        "PWY-6163: chorismate biosynthesis from 3-dehydroquinate",
        "PWY-6122: 5-aminoimidazole ribonucleotide biosynthesis II",
        "PANTOSYN-PWY: pantothenate and coenzyme A biosynthesis I",
        "PWY-7234: inosine-5'-phosphate biosynthesis III",
    ],
    "mean_abundance": [
        0.0182, 0.0165, 0.0148, 0.0137, 0.0121,
        0.0118, 0.0095, 0.0089, 0.0082, 0.0076,
    ],
}


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def run_command(
    cmd: List[str],
    description: str,
    cwd: Optional[Path] = None,
    env: Optional[dict] = None,
    critical: bool = True,
) -> subprocess.CompletedProcess:
    """Run a subprocess command with proper error handling.

    Args:
        cmd: Command and arguments to run.
        description: Human-readable description for log output.
        cwd: Working directory for the subprocess.
        env: Environment variables for the subprocess.
        critical: If True (default), raise RuntimeError on non-zero exit code.
                  If False, log a warning and continue.
    """
    print(f"  [{description}] {' '.join(cmd[:4])}...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            cwd=cwd,
            env=env,
            timeout=7200,  # 2-hour timeout for large metagenomes
        )
        if result.returncode != 0:
            stderr_tail = ""
            if result.stderr:
                stderr_lines = result.stderr.strip().split("\n")
                stderr_tail = "\n".join(
                    f"    STDERR: {line}" for line in stderr_lines[-5:]
                )
            if critical:
                msg = (
                    f"{description} failed with exit code {result.returncode}."
                )
                if stderr_tail:
                    msg += f"\n{stderr_tail}"
                raise RuntimeError(msg)
            else:
                print(
                    f"  WARNING: {description} returned non-zero exit code "
                    f"{result.returncode}",
                    file=sys.stderr,
                )
                if stderr_tail:
                    print(stderr_tail, file=sys.stderr)
        return result
    except FileNotFoundError:
        print(f"  ERROR: '{cmd[0]}' not found. Is it installed and on PATH?",
              file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print(f"  ERROR: {description} timed out after 2 hours.", file=sys.stderr)
        sys.exit(1)


def detect_threads() -> int:
    """Auto-detect thread count from available CPUs."""
    cpu_count = os.cpu_count()
    return max(1, cpu_count - 1) if cpu_count and cpu_count > 1 else 1


def check_tool(name: str) -> bool:
    """Check if an external tool is available on PATH."""
    return shutil.which(name) is not None


# ---------------------------------------------------------------------------
# Pipeline step 1: Kraken2 + Bracken taxonomy
# ---------------------------------------------------------------------------

def run_kraken2(
    r1: Path,
    r2: Optional[Path],
    output_dir: Path,
    db_path: Path,
    threads: int,
    confidence: float = KRAKEN2_CONFIDENCE,
) -> Path:
    """Run Kraken2 taxonomic classification."""
    report_path = output_dir / "kraken2_report.txt"
    output_path = output_dir / "kraken2_output.txt"

    cmd = [
        "kraken2",
        "--db", str(db_path),
        "--threads", str(threads),
        "--confidence", str(confidence),
        "--report", str(report_path),
        "--output", str(output_path),
    ]

    if r2 is not None:
        cmd.extend(["--paired", str(r1), str(r2)])
    else:
        cmd.append(str(r1))

    run_command(cmd, "Kraken2 classification", critical=True)

    if not report_path.exists():
        raise FileNotFoundError(
            f"Kraken2 did not produce the expected report file: {report_path}"
        )
    return report_path


def run_bracken(
    kraken_report: Path,
    output_dir: Path,
    db_path: Path,
    read_length: int = DEFAULT_READ_LENGTH,
    level: str = BRACKEN_LEVEL,
    threshold: int = BRACKEN_THRESHOLD,
) -> Path:
    """Run Bracken abundance re-estimation."""
    bracken_output = output_dir / f"bracken_{level.lower()}.tsv"

    cmd = [
        "bracken",
        "-d", str(db_path),
        "-i", str(kraken_report),
        "-o", str(bracken_output),
        "-r", str(read_length),
        "-l", level,
        "-t", str(threshold),
    ]

    run_command(cmd, "Bracken re-estimation")

    if not bracken_output.exists():
        print(f"  ERROR: Bracken did not produce the expected output file: "
              f"{bracken_output}", file=sys.stderr)
    return bracken_output


def parse_bracken_output(bracken_path: Path) -> pd.DataFrame:
    """Parse Bracken species-level output into a DataFrame."""
    df = pd.read_csv(bracken_path, sep="\t")
    # Bracken output columns: name, taxonomy_id, taxonomy_lvl,
    # kraken_assigned_reads, added_reads, new_est_reads, fraction_total_reads
    if "fraction_total_reads" in df.columns:
        df = df.sort_values("fraction_total_reads", ascending=False)
        df["relative_abundance"] = df["fraction_total_reads"] * 100
    elif "new_est_reads" in df.columns:
        total = df["new_est_reads"].sum()
        df["relative_abundance"] = (df["new_est_reads"] / total * 100) if total > 0 else 0
        df = df.sort_values("relative_abundance", ascending=False)
    return df


# ---------------------------------------------------------------------------
# Pipeline step 2: RGI resistome profiling
# ---------------------------------------------------------------------------

def run_rgi(
    r1: Path,
    r2: Optional[Path],
    output_dir: Path,
    threads: int,
) -> Path:
    """Run RGI MAIN for antimicrobial resistance gene detection."""
    # RGI operates on a single concatenated input or an assembly
    # For read-level: we use rgi bwt (read-based)
    rgi_output = output_dir / "rgi_results"

    if r2 is not None:
        cmd = [
            "rgi", "bwt",
            "--read_one", str(r1),
            "--read_two", str(r2),
            "--output_file", str(rgi_output),
            "--threads", str(threads),
            "--aligner", "bowtie2",
            "--include_wildcard",
        ]
    else:
        cmd = [
            "rgi", "bwt",
            "--read_one", str(r1),
            "--output_file", str(rgi_output),
            "--threads", str(threads),
            "--aligner", "bowtie2",
            "--include_wildcard",
        ]

    run_command(cmd, "RGI resistome profiling", critical=True)

    # RGI bwt outputs: <prefix>.allele_mapping_data.txt
    allele_file = Path(str(rgi_output) + ".allele_mapping_data.txt")
    gene_file = Path(str(rgi_output) + ".gene_mapping_data.txt")

    # Prefer allele-level results, fall back to gene-level
    if allele_file.exists():
        return allele_file
    elif gene_file.exists():
        return gene_file
    else:
        raise FileNotFoundError(
            f"RGI did not produce expected output files: "
            f"{allele_file} or {gene_file}"
        )


def parse_rgi_output(rgi_path: Path) -> pd.DataFrame:
    """Parse RGI results into a DataFrame.

    If the file does not exist, returns an empty DataFrame with a
    ``_analysis_status`` column set to ``"FAILED"`` so that downstream
    code can distinguish "no file produced" from "zero ARGs detected".
    """
    if not rgi_path.exists():
        df = pd.DataFrame(columns=["_analysis_status"])
        df.attrs["_analysis_status"] = "FAILED"
        return df

    df = pd.read_csv(rgi_path, sep="\t")
    # Filter to Perfect and Strict hits only
    if "Cut_Off" in df.columns:
        df = df[df["Cut_Off"].isin(RGI_CRITERIA)]
    return df


def classify_who_priority(resistome_df: pd.DataFrame) -> pd.DataFrame:
    """Classify detected ARGs by WHO critical priority tier."""
    if resistome_df.empty:
        return pd.DataFrame(columns=["gene", "drug_class", "who_priority"])

    results = []
    gene_col = "Best_Hit_ARO" if "Best_Hit_ARO" in resistome_df.columns else "gene"
    drug_col = "Drug Class" if "Drug Class" in resistome_df.columns else "drug_class"

    for _, row in resistome_df.iterrows():
        gene_name = str(row.get(gene_col, ""))
        drug_class = str(row.get(drug_col, "")).lower()
        priority = "Not classified"

        for tier, info in WHO_PRIORITY_ARGS.items():
            for family in info["arg_families"]:
                if family.lower() in gene_name.lower():
                    priority = tier
                    break
            if priority != "Not classified":
                break
            for dc in info["drug_classes"]:
                if dc in drug_class:
                    priority = tier
                    break
            if priority != "Not classified":
                break

        results.append({
            "gene": gene_name,
            "drug_class": drug_class,
            "who_priority": priority,
        })

    return pd.DataFrame(results)


# ---------------------------------------------------------------------------
# Pipeline step 3: HUMAnN3 functional profiling
# ---------------------------------------------------------------------------

def run_humann3(
    r1: Path,
    r2: Optional[Path],
    output_dir: Path,
    threads: int,
    db_path: Optional[Path] = None,
) -> Optional[Path]:
    """Run HUMAnN3 for functional pathway profiling.

    Returns the path to the pathway abundance file, or None if no output
    was produced.
    """
    humann_dir = output_dir / "humann3"
    humann_dir.mkdir(parents=True, exist_ok=True)

    # HUMAnN3 expects a single concatenated input
    if r2 is not None:
        concat_path = output_dir / "concat_reads.fastq.gz"
        if not concat_path.exists():
            print("  [Concatenating R1+R2 for HUMAnN3]...")
            with open(concat_path, "wb") as out_f:
                for fq in (r1, r2):
                    with open(fq, "rb") as in_f:
                        for chunk in iter(lambda: in_f.read(65536), b""):
                            out_f.write(chunk)
        input_path = concat_path
    else:
        input_path = r1

    cmd = [
        "humann",
        "--input", str(input_path),
        "--output", str(humann_dir),
        "--threads", str(threads),
    ]

    if db_path is not None:
        cmd.extend(["--nucleotide-database", str(db_path / "chocophlan")])
        cmd.extend(["--protein-database", str(db_path / "uniref")])

    run_command(cmd, "HUMAnN3 functional profiling")

    # Find pathway abundance file
    pathabundance_files = list(humann_dir.glob("*pathabundance.tsv"))
    if pathabundance_files:
        return pathabundance_files[0]

    print("  WARNING: HUMAnN3 pathway abundance file not found.", file=sys.stderr)
    return None


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------

def plot_taxonomy_barchart(
    taxonomy_df: pd.DataFrame,
    output_path: Path,
    top_n: int = 20,
    site_columns: Optional[List[str]] = None,
) -> None:
    """Stacked bar chart of top N species by relative abundance."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if site_columns and all(c in taxonomy_df.columns for c in site_columns):
        # Multi-site: stacked grouped bar chart
        top_species = taxonomy_df.head(top_n)
        species = top_species["species"].tolist()
        n_species = len(species)
        n_sites = len(site_columns)
        x = np.arange(n_sites)
        bar_width = 0.7

        fig, ax = plt.subplots(figsize=(10, 7))

        bottom = np.zeros(n_sites)
        for i, sp in enumerate(species):
            values = [top_species[col].iloc[i] for col in site_columns]
            colour = TAXON_COLOURS[i % len(TAXON_COLOURS)]
            ax.bar(x, values, bar_width, bottom=bottom, label=sp, color=colour)
            bottom += np.array(values)

        ax.set_xlabel("Collection Site", fontsize=12)
        ax.set_ylabel("Relative Abundance (%)", fontsize=12)
        ax.set_title("Taxonomic Composition (Top %d Species)" % top_n, fontsize=14)
        ax.set_xticks(x)
        ax.set_xticklabels(site_columns, fontsize=11)
        ax.legend(
            bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=8,
            title="Species", title_fontsize=9,
        )
        ax.grid(True, axis="y", alpha=0.3)

    elif "relative_abundance" in taxonomy_df.columns:
        # Single sample: horizontal bar chart
        top_species = taxonomy_df.head(top_n)
        fig, ax = plt.subplots(figsize=(10, 7))
        species = top_species["name"].tolist() if "name" in top_species.columns else top_species.index.tolist()
        abundances = top_species["relative_abundance"].tolist()
        colours = [TAXON_COLOURS[i % len(TAXON_COLOURS)] for i in range(len(species))]

        y_pos = np.arange(len(species))
        ax.barh(y_pos, abundances, color=colours, edgecolor="white", linewidth=0.5)
        ax.set_yticks(y_pos)
        ax.set_yticklabels(species, fontsize=9)
        ax.invert_yaxis()
        ax.set_xlabel("Relative Abundance (%)", fontsize=12)
        ax.set_title("Taxonomic Composition (Top %d Species)" % top_n, fontsize=14)
        ax.grid(True, axis="x", alpha=0.3)
    else:
        return

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_resistome_heatmap(
    resistome_df: pd.DataFrame,
    output_path: Path,
    site_columns: Optional[List[str]] = None,
) -> None:
    """Heatmap of ARG families by drug class."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import seaborn as sns

    if site_columns and all(c in resistome_df.columns for c in site_columns):
        # Multi-site heatmap: rows = genes, columns = sites
        genes = resistome_df["gene"].tolist()
        data = resistome_df[site_columns].values
        heat_df = pd.DataFrame(data, index=genes, columns=site_columns)
    elif "RPKM" in resistome_df.columns or "Reads" in resistome_df.columns:
        # Single sample: pivot by drug class vs gene
        val_col = "RPKM" if "RPKM" in resistome_df.columns else "Reads"
        gene_col = "Best_Hit_ARO" if "Best_Hit_ARO" in resistome_df.columns else "gene"
        drug_col = "Drug Class" if "Drug Class" in resistome_df.columns else "drug_class"
        heat_df = resistome_df.pivot_table(
            values=val_col, index=gene_col, columns=drug_col,
            aggfunc="sum", fill_value=0,
        )
    else:
        return

    fig, ax = plt.subplots(figsize=(12, max(6, len(heat_df) * 0.35)))
    sns.heatmap(
        heat_df,
        ax=ax,
        cmap="YlOrRd",
        linewidths=0.5,
        linecolor="white",
        annot=True,
        fmt=".0f",
        cbar_kws={"label": "Read Count"},
    )
    ax.set_title("Antimicrobial Resistance Gene Profile", fontsize=14)
    ax.set_ylabel("Resistance Gene", fontsize=12)
    ax.set_xlabel("Collection Site" if site_columns else "Drug Class", fontsize=12)
    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


def plot_who_critical_args(
    who_df: pd.DataFrame,
    output_path: Path,
) -> None:
    """Summary figure of WHO-critical priority ARG detections."""
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if who_df.empty:
        return

    # Count genes per priority tier
    tier_counts = who_df["who_priority"].value_counts()
    tier_order = ["Critical", "High", "Medium", "Not classified"]
    tier_colours = {
        "Critical": "#C44E52",
        "High": "#DD8452",
        "Medium": "#CCB974",
        "Not classified": "#8C8C8C",
    }

    counts = [tier_counts.get(t, 0) for t in tier_order]
    colours = [tier_colours[t] for t in tier_order]

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.5]})

    # Panel A: Pie chart of priority tiers
    ax_pie = axes[0]
    nonzero = [(t, c, col) for t, c, col in zip(tier_order, counts, colours) if c > 0]
    if nonzero:
        labels, vals, cols = zip(*nonzero)
        wedges, texts, autotexts = ax_pie.pie(
            vals, labels=labels, colors=cols, autopct="%1.0f%%",
            startangle=90, textprops={"fontsize": 10},
        )
        for at in autotexts:
            at.set_fontweight("bold")
    ax_pie.set_title("WHO Priority Classification", fontsize=13, fontweight="bold")

    # Panel B: Breakdown of Critical-tier genes by drug class
    ax_bar = axes[1]
    critical_df = who_df[who_df["who_priority"] == "Critical"]
    if not critical_df.empty:
        drug_counts = critical_df["drug_class"].value_counts()
        drug_classes = drug_counts.index.tolist()
        drug_vals = drug_counts.values.tolist()
        bar_colours = [tier_colours["Critical"]] * len(drug_classes)

        y_pos = np.arange(len(drug_classes))
        ax_bar.barh(y_pos, drug_vals, color=bar_colours, edgecolor="white")
        ax_bar.set_yticks(y_pos)
        ax_bar.set_yticklabels(
            [dc.title() for dc in drug_classes], fontsize=10,
        )
        ax_bar.invert_yaxis()
        ax_bar.set_xlabel("Number of ARGs Detected", fontsize=11)
        ax_bar.set_title(
            "WHO-Critical ARGs by Drug Class", fontsize=13, fontweight="bold",
        )
        ax_bar.grid(True, axis="x", alpha=0.3)

        # Annotate bars with gene names
        for i, dc in enumerate(drug_classes):
            genes = critical_df[critical_df["drug_class"] == dc]["gene"].tolist()
            gene_str = ", ".join(genes[:5])
            if len(genes) > 5:
                gene_str += f" (+{len(genes) - 5})"
            ax_bar.text(
                drug_vals[i] + 0.1, i, gene_str,
                va="center", fontsize=8, style="italic",
            )
    else:
        ax_bar.text(
            0.5, 0.5, "No WHO-Critical ARGs detected",
            ha="center", va="center", fontsize=14, transform=ax_bar.transAxes,
        )
        ax_bar.set_title(
            "WHO-Critical ARGs by Drug Class", fontsize=13, fontweight="bold",
        )

    plt.tight_layout()
    fig.savefig(output_path, dpi=300, bbox_inches="tight")
    plt.close(fig)
    print(f"  Saved: {output_path.name}")


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------

def generate_report(
    output_dir: Path,
    taxonomy_df: pd.DataFrame,
    resistome_df: pd.DataFrame,
    who_df: pd.DataFrame,
    pathways_df: Optional[pd.DataFrame],
    figures: Dict[str, Path],
    input_files: List[Path],
    is_demo: bool,
    rgi_failed: bool = False,
    humann_skipped: bool = False,
    humann_failed: bool = False,
) -> Path:
    """Generate the full markdown report.

    Args:
        rgi_failed: True if RGI did not produce output (report says FAILED).
        humann_skipped: True if HUMAnN3 was intentionally skipped by user.
        humann_failed: True if HUMAnN3 ran but did not produce output.
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    report_path = output_dir / "report.md"

    # Figure references
    fig_refs = {}
    for name, path in figures.items():
        try:
            rel = path.relative_to(output_dir)
        except ValueError:
            rel = path.name
        fig_refs[name] = f"![{name}]({rel})"

    # Taxonomy summary
    if "species" in taxonomy_df.columns:
        top5 = taxonomy_df.head(5)
        tax_lines = []
        for _, row in top5.iterrows():
            tax_lines.append(f"  - *{row['species']}*")
    elif "name" in taxonomy_df.columns:
        top5 = taxonomy_df.head(5)
        tax_lines = []
        for _, row in top5.iterrows():
            tax_lines.append(
                f"  - *{row['name']}* ({row.get('relative_abundance', 0):.1f}%)"
            )
    else:
        tax_lines = ["  - No taxonomy data available"]

    # Resistome summary
    if rgi_failed:
        resistome_section = """## Antimicrobial Resistance Profile

> **ANALYSIS FAILED**: RGI did not produce output. Resistome could not be profiled.
> Re-run with valid inputs or check that the CARD database is correctly installed.

"""
    else:
        n_args = len(resistome_df) if not resistome_df.empty else 0
        n_critical = len(who_df[who_df["who_priority"] == "Critical"]) if not who_df.empty else 0
        n_high = len(who_df[who_df["who_priority"] == "High"]) if not who_df.empty else 0
        resistome_section = f"""## Antimicrobial Resistance Profile

{fig_refs.get("resistome", "")}

Total ARG detections: **{n_args}**
- WHO-Critical priority: **{n_critical}**
- WHO-High priority: **{n_high}**

### WHO-Critical ARG Summary

{fig_refs.get("who_critical", "")}

"""

    # Pathway summary
    pathway_section = ""
    if humann_failed:
        pathway_section = """## Functional Pathways (HUMAnN3)

> **FUNCTIONAL PROFILING FAILED**: HUMAnN3 did not produce pathway output.
> Check database installation and input file compatibility.

"""
    elif humann_skipped:
        pathway_section = """## Functional Pathways (HUMAnN3)

*Functional profiling was skipped (no database specified or --skip-functional used).*

"""
    elif pathways_df is not None and not pathways_df.empty:
        n_pathways = len(pathways_df)
        top_pw = pathways_df.iloc[0]
        pw_name = top_pw.get("pathway", top_pw.get("# Pathway", "Unknown"))
        pathway_section = f"""## Functional Pathways (HUMAnN3)

Total pathways detected: **{n_pathways}**

Top pathway: {pw_name}

"""

    # Input checksums
    checksum_lines = []
    for fp in input_files:
        if fp.exists():
            cs = sha256_file(fp)
            checksum_lines.append(f"  - `{fp.name}`: `{cs}`")
        else:
            checksum_lines.append(f"  - `{fp.name}`: demo data (no file)")

    report = f"""# Metagenomics Profiling Report

**Date**: {now}
**Tool**: ClawBio Metagenomics Profiler v0.1.0
**Mode**: {"Demo (pre-computed Peru sewage data)" if is_demo else "Full pipeline"}

---

## Taxonomic Composition

{fig_refs.get("taxonomy", "")}

### Top Species

{chr(10).join(tax_lines)}

{resistome_section}{pathway_section}---

## Methods

- **Taxonomic classification**: Kraken2 (confidence threshold {KRAKEN2_CONFIDENCE}) with Bracken species-level re-estimation (minimum {BRACKEN_THRESHOLD} reads)
- **Resistome profiling**: RGI against the CARD database (Perfect + Strict criteria only)
- **Functional profiling**: HUMAnN3 with MetaCyc pathway stratification
- **WHO Priority mapping**: ARGs classified against the 2024 WHO Bacterial Priority Pathogens List
- **Figures**: 300 dpi, publication-quality

## Input Checksums

{chr(10).join(checksum_lines)}

## Disclaimer

> {SAFETY_DISCLAIMER}

## Citations

- Wood, D.E., Lu, J. & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257.
- Lu, J. et al. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- Alcock, B.P. et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51(D1), D419-D430.
- Beghini, F. et al. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife*, 10, e65088.
- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio
"""

    report_path.write_text(report, encoding="utf-8")
    print(f"  Saved: report.md")
    return report_path


# ---------------------------------------------------------------------------
# Reproducibility bundle
# ---------------------------------------------------------------------------

# The pipeline's external tools and Python libraries, for environment.yml.
CONDA_DEPS = [
    "kraken2>=2.1.3",
    "bracken>=2.9",
    "rgi>=6.0",
    "humann>=3.8",
    "pandas>=2.0",
    "numpy>=1.24",
    "matplotlib>=3.7",
    "seaborn>=0.12",
    "scipy>=1.10",
    "biopython>=1.82",
]

# A real run checks the tools before it starts; the demo needs none of them.
TOOL_PREFLIGHT = [
    "# Verify tool versions",
    'kraken2 --version 2>/dev/null || echo "WARNING: kraken2 not found"',
    'bracken -v 2>/dev/null || echo "WARNING: bracken not found"',
    'rgi main --version 2>/dev/null || echo "WARNING: rgi not found"',
    'humann --version 2>/dev/null || echo "WARNING: humann not found"',
]


def _quoted_path(value: str) -> str:
    """A path for commands.sh: under $CLAWBIO_ROOT when inside the checkout, else absolute."""
    path = Path(value).expanduser().resolve()
    if path.is_relative_to(_PROJECT_ROOT):
        return f'"$CLAWBIO_ROOT"/{shlex.quote(path.relative_to(_PROJECT_ROOT).as_posix())}'
    return shlex.quote(str(path))


def _replay_args(args: argparse.Namespace, output_dir: Path, is_demo: bool) -> list[str | ReproPath]:
    if is_demo:
        return ["--demo", "--output", ReproPath(output_dir, anchor="output_dir")]
    replay: list[str | ReproPath] = []
    for flag, value in (("--r1", args.r1), ("--r2", args.r2), ("--input", args.input)):
        if value:
            replay += [flag, _quoted_path(value)]
    if args.skip_functional:
        replay.append("--skip-functional")
    for flag, value in (("--kraken2-db", args.kraken2_db),
                        ("--humann-db", getattr(args, "humann_db", None))):
        if value:
            replay += [flag, _quoted_path(value)]
    if args.read_length != DEFAULT_READ_LENGTH:
        replay += ["--read-length", str(args.read_length)]
    return replay + ["--output", ReproPath(output_dir, anchor="output_dir")]


def write_reproducibility_bundle(
    output_dir: Path,
    input_files: List[Path],
    args: argparse.Namespace,
    is_demo: bool,
) -> None:
    """Write reproducibility/ with the shared helpers: commands.sh, environment.yml,
    checksums.sha256 and, for a real run, inputs.json with each input's SHA-256."""
    output_dir = Path(output_dir).resolve()
    repro_dir = output_dir / "reproducibility"
    outputs = sorted(p for p in output_dir.rglob("*") if p.is_file() and repro_dir not in p.parents)
    repro_dir.mkdir(parents=True, exist_ok=True)

    bundle = []
    if not is_demo:
        # Inputs live outside the output directory, so their digests are recorded
        # here by flag rather than in checksums.sha256.
        inputs = {
            flag: {"file": Path(value).name, "sha256": sha256_file(Path(value))}
            for flag, value in (("r1", args.r1), ("r2", args.r2), ("input", args.input))
            if value and Path(value).is_file()
        }
        inputs_path = repro_dir / "inputs.json"
        write_text_lf_atomic(inputs_path, json.dumps(inputs, indent=2) + "\n")
        bundle.append(inputs_path)

    bundle.append(write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills/claw-metagenomics/metagenomics_profiler.py"),
            args=_replay_args(args, output_dir, is_demo),
            comment="Reproduce this ClawBio Metagenomics Profiler run",
            preflight=[] if is_demo else TOOL_PREFLIGHT,
        ),
        repo_root=_PROJECT_ROOT,
    ))
    bundle.append(write_environment_yml(
        output_dir,
        env_name="clawbio-metagenomics",
        pip_deps=[],
        conda_deps=CONDA_DEPS,
        channels=["conda-forge", "bioconda"],
    ))
    write_checksums([*outputs, *bundle], output_dir, anchor=output_dir)
    print("  Saved: reproducibility/ (commands.sh, environment.yml, checksums.sha256)")


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------

def run_demo(output_dir: Path) -> None:
    """Generate full report from pre-computed Peru sewage metagenomics data."""
    print("Metagenomics Profiler -- ClawBio")
    print("=" * 40)
    print("Mode: demo (pre-computed Peru sewage data)")
    print("Samples: 6 (3 sites: Lima, Cusco, Iquitos)\n")

    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(exist_ok=True)

    # --- Taxonomy ---
    print("Generating taxonomy data...")
    tax_df = pd.DataFrame(DEMO_TAXONOMY)
    tax_df.to_csv(tables_dir / "taxonomy_species.tsv", sep="\t", index=False)

    classified_pct = 94.2
    top_sp = tax_df.iloc[0]
    print(f"  Total classified: {classified_pct}%")
    print(f"  Top species: {top_sp['species']} "
          f"(Lima: {top_sp['Lima']}%, Cusco: {top_sp['Cusco']}%, "
          f"Iquitos: {top_sp['Iquitos']}%)")

    # --- Resistome ---
    print("Generating resistome data...")
    res_df = pd.DataFrame(DEMO_RESISTOME)
    res_df.to_csv(tables_dir / "resistome_profile.tsv", sep="\t", index=False)

    total_perfect = sum(1 for c in DEMO_RESISTOME["criteria"] if c == "Perfect")
    total_strict = sum(1 for c in DEMO_RESISTOME["criteria"] if c == "Strict")
    n_drug_classes = len(set(DEMO_RESISTOME["drug_class"]))
    print(f"  Total ARG hits: {len(res_df)} "
          f"(Perfect: {total_perfect}, Strict: {total_strict})")
    print(f"  Drug classes: {n_drug_classes}")

    # WHO classification for demo data
    who_records = []
    for i, gene in enumerate(DEMO_RESISTOME["gene"]):
        drug_class = DEMO_RESISTOME["drug_class"][i]
        priority = "Not classified"
        for tier, info in WHO_PRIORITY_ARGS.items():
            for family in info["arg_families"]:
                if family.lower() in gene.lower():
                    priority = tier
                    break
            if priority != "Not classified":
                break
            for dc in info["drug_classes"]:
                if dc in drug_class:
                    priority = tier
                    break
            if priority != "Not classified":
                break
        who_records.append({
            "gene": gene,
            "drug_class": drug_class,
            "who_priority": priority,
        })

    who_df = pd.DataFrame(who_records)
    who_df.to_csv(tables_dir / "who_priority_args.tsv", sep="\t", index=False)

    n_critical = len(who_df[who_df["who_priority"] == "Critical"])
    critical_genes = who_df[who_df["who_priority"] == "Critical"]["gene"].tolist()
    print(f"  WHO-Critical ARGs detected: {n_critical}")
    print(f"    - {', '.join(critical_genes)}")

    # --- Functional pathways ---
    print("Generating pathway data...")
    pw_df = pd.DataFrame(DEMO_PATHWAYS)
    pw_df.to_csv(tables_dir / "pathway_abundance.tsv", sep="\t", index=False)
    print(f"  Total pathways: {len(pw_df)}")
    print(f"  Top: {pw_df.iloc[0]['pathway']}")

    # --- Figures ---
    print("\nGenerating figures...")
    figures: Dict[str, Path] = {}

    try:
        # Figure 1: Taxonomy bar chart
        tax_fig_path = figures_dir / "taxonomy_barplot.png"
        site_cols = ["Lima", "Cusco", "Iquitos"]
        plot_taxonomy_barchart(tax_df, tax_fig_path, top_n=20, site_columns=site_cols)
        figures["taxonomy"] = tax_fig_path

        # Figure 2: Resistome heatmap
        res_fig_path = figures_dir / "resistome_heatmap.png"
        plot_resistome_heatmap(res_df, res_fig_path, site_columns=site_cols)
        figures["resistome"] = res_fig_path

        # Figure 3: WHO-critical ARG summary
        who_fig_path = figures_dir / "who_critical_args.png"
        plot_who_critical_args(who_df, who_fig_path)
        figures["who_critical"] = who_fig_path

    except ImportError as e:
        print(f"  WARNING: {e}. Some figures skipped.", file=sys.stderr)

    # --- Report ---
    print("\nGenerating report...")
    generate_report(
        output_dir=output_dir,
        taxonomy_df=tax_df,
        resistome_df=res_df,
        who_df=who_df,
        pathways_df=pw_df,
        figures=figures,
        input_files=[Path("demo_peru_sewage.fastq.gz")],
        is_demo=True,
    )

    # --- Reproducibility ---
    demo_args = argparse.Namespace(
        r1=None, r2=None, input=None, output=str(output_dir),
        skip_functional=False, kraken2_db=None, read_length=DEFAULT_READ_LENGTH,
    )
    write_reproducibility_bundle(
        output_dir=output_dir,
        input_files=[Path("demo_peru_sewage.fastq.gz")],
        args=demo_args,
        is_demo=True,
    )

    print(f"\nDone. Output: {output_dir}")
    print(f"  Report:  {output_dir / 'report.md'}")
    print(f"  Figures: {figures_dir}/")
    print(f"  Tables:  {tables_dir}/")


# ---------------------------------------------------------------------------
# Full pipeline
# ---------------------------------------------------------------------------

def run_pipeline(args: argparse.Namespace) -> None:
    """Run the full metagenomics profiling pipeline."""
    print("Metagenomics Profiler -- ClawBio")
    print("=" * 40)

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    figures_dir.mkdir(exist_ok=True)
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(exist_ok=True)

    threads = detect_threads()
    print(f"Threads: {threads}")

    # Resolve input files
    if args.r1:
        r1 = Path(args.r1)
        r2 = Path(args.r2) if args.r2 else None
        input_files = [r1] + ([r2] if r2 else [])
    elif args.input:
        r1 = Path(args.input)
        r2 = None
        input_files = [r1]
    else:
        print("ERROR: Provide --r1/--r2 or --input.", file=sys.stderr)
        sys.exit(1)

    # Validate inputs exist
    for fp in input_files:
        if not fp.exists():
            print(f"ERROR: Input file not found: {fp}", file=sys.stderr)
            sys.exit(1)

    # Check tools
    required_tools = ["kraken2", "bracken", "rgi"]
    if not args.skip_functional:
        required_tools.append("humann")

    missing = [t for t in required_tools if not check_tool(t)]
    if missing:
        print(f"ERROR: Required tools not found on PATH: {', '.join(missing)}",
              file=sys.stderr)
        print("Install with: conda install -c bioconda " + " ".join(missing),
              file=sys.stderr)
        sys.exit(1)

    # Resolve Kraken2 database
    kraken2_db = Path(args.kraken2_db) if args.kraken2_db else None
    if kraken2_db is None:
        env_db = os.environ.get("KRAKEN2_DB")
        if env_db:
            kraken2_db = Path(env_db)

    if kraken2_db is None or not kraken2_db.exists():
        print("ERROR: Kraken2 database not found. Set --kraken2-db or $KRAKEN2_DB.",
              file=sys.stderr)
        sys.exit(1)

    print(f"Kraken2 DB: {kraken2_db}")
    print(f"Input: {', '.join(str(f) for f in input_files)}\n")

    # --- Step 1: Taxonomy ---
    print("Step 1/3: Taxonomic classification")
    kraken_report = run_kraken2(r1, r2, output_dir, kraken2_db, threads)
    bracken_output = run_bracken(
        kraken_report, output_dir, kraken2_db,
        read_length=args.read_length,
    )
    tax_df = parse_bracken_output(bracken_output)
    tax_df.to_csv(tables_dir / "taxonomy_species.tsv", sep="\t", index=False)
    print(f"  Species detected: {len(tax_df)}")
    if not tax_df.empty:
        top = tax_df.iloc[0]
        name_col = "name" if "name" in tax_df.columns else tax_df.columns[0]
        print(f"  Top species: {top[name_col]} "
              f"({top.get('relative_abundance', 0):.1f}%)\n")

    # --- Step 2: Resistome ---
    print("Step 2/3: Resistome profiling")
    rgi_failed = False
    try:
        rgi_output = run_rgi(r1, r2, output_dir, threads)
        res_df = parse_rgi_output(rgi_output)
    except (RuntimeError, FileNotFoundError) as exc:
        print(f"  ERROR: RGI failed: {exc}", file=sys.stderr)
        rgi_failed = True
        res_df = pd.DataFrame()

    if rgi_failed:
        who_df = pd.DataFrame(columns=["gene", "drug_class", "who_priority"])
        print("  RESISTOME ANALYSIS FAILED -- see errors above.\n")
    else:
        who_df = classify_who_priority(res_df)
        res_df.to_csv(tables_dir / "resistome_profile.tsv", sep="\t", index=False)
        who_df.to_csv(tables_dir / "who_priority_args.tsv", sep="\t", index=False)
        n_critical = len(who_df[who_df["who_priority"] == "Critical"]) if not who_df.empty else 0
        print(f"  ARG hits: {len(res_df)} (WHO-Critical: {n_critical})\n")

    # --- Step 3: Functional (optional) ---
    pw_df = None
    humann_skipped = False
    humann_failed = False
    if not args.skip_functional:
        print("Step 3/3: Functional profiling")
        humann_db = Path(args.humann_db) if args.humann_db else None
        if humann_db is None:
            env_db = os.environ.get("HUMANN_DB")
            if env_db:
                humann_db = Path(env_db)
            else:
                print("  NOTE: No HUMAnN3 database specified (--humann-db or "
                      "$HUMANN_DB). HUMAnN3 will be skipped.", file=sys.stderr)
                humann_skipped = True
        if humann_db is not None and not humann_db.exists():
            print(f"  WARNING: HUMAnN3 database path does not exist: {humann_db}",
                  file=sys.stderr)
            humann_skipped = True

        if not humann_skipped:
            pw_output = run_humann3(r1, r2, output_dir, threads, humann_db)
            if pw_output is not None and pw_output.exists():
                pw_df = pd.read_csv(pw_output, sep="\t", comment="#")
                pw_df.to_csv(tables_dir / "pathway_abundance.tsv", sep="\t", index=False)
                print(f"  Pathways: {len(pw_df)}\n")
            else:
                humann_failed = True
                print("  ERROR: HUMAnN3 did not produce pathway output.\n",
                      file=sys.stderr)
        else:
            print("  Step 3/3: Functional profiling (SKIPPED -- no database)\n")
    else:
        humann_skipped = True
        print("Step 3/3: Functional profiling (SKIPPED)\n")

    # --- Figures ---
    print("Generating figures...")
    figures: Dict[str, Path] = {}
    try:
        tax_fig_path = figures_dir / "taxonomy_barplot.png"
        plot_taxonomy_barchart(tax_df, tax_fig_path)
        figures["taxonomy"] = tax_fig_path

        if not rgi_failed and not res_df.empty:
            res_fig_path = figures_dir / "resistome_heatmap.png"
            plot_resistome_heatmap(res_df, res_fig_path)
            figures["resistome"] = res_fig_path

        if not rgi_failed and not who_df.empty:
            who_fig_path = figures_dir / "who_critical_args.png"
            plot_who_critical_args(who_df, who_fig_path)
            figures["who_critical"] = who_fig_path
    except ImportError as e:
        print(f"  WARNING: {e}. Some figures skipped.", file=sys.stderr)

    # --- Report ---
    print("\nGenerating report...")
    generate_report(
        output_dir=output_dir,
        taxonomy_df=tax_df,
        resistome_df=res_df,
        who_df=who_df,
        pathways_df=pw_df,
        figures=figures,
        input_files=input_files,
        is_demo=False,
        rgi_failed=rgi_failed,
        humann_skipped=humann_skipped,
        humann_failed=humann_failed,
    )

    # --- Reproducibility ---
    write_reproducibility_bundle(
        output_dir=output_dir,
        input_files=input_files,
        args=args,
        is_demo=False,
    )

    print(f"\nDone. Output: {output_dir}")
    print(f"  Report:  {output_dir / 'report.md'}")
    print(f"  Figures: {figures_dir}/")
    print(f"  Tables:  {tables_dir}/")
    print(f"\n{SAFETY_DISCLAIMER}")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(
        description=(
            "ClawBio Metagenomics Profiler: Kraken2/Bracken taxonomy + "
            "RGI/CARD resistome + HUMAnN3 functional profiling"
        ),
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Full pipeline\n"
            "  python metagenomics_profiler.py --r1 R1.fq.gz --r2 R2.fq.gz -o report/\n"
            "\n"
            "  # Demo mode (no external tools needed)\n"
            "  python metagenomics_profiler.py --demo -o demo_report/\n"
            "\n"
            "  # Skip HUMAnN3 (faster)\n"
            "  python metagenomics_profiler.py --r1 R1.fq.gz --r2 R2.fq.gz "
            "--skip-functional -o report/\n"
        ),
    )

    input_group = parser.add_argument_group("Input")
    input_group.add_argument(
        "--r1", default=None,
        help="Forward reads FASTQ (R1)",
    )
    input_group.add_argument(
        "--r2", default=None,
        help="Reverse reads FASTQ (R2)",
    )
    input_group.add_argument(
        "--input", "-i", default=None,
        help="Single concatenated/interleaved FASTQ (alternative to --r1/--r2)",
    )

    output_group = parser.add_argument_group("Output")
    output_group.add_argument(
        "--output", "-o", default="metagenomics_report",
        help="Output directory (default: metagenomics_report)",
    )

    pipeline_group = parser.add_argument_group("Pipeline options")
    pipeline_group.add_argument(
        "--demo", action="store_true",
        help="Run demo with pre-computed Peru sewage data (no tools needed)",
    )
    pipeline_group.add_argument(
        "--skip-functional", action="store_true",
        help="Skip HUMAnN3 functional profiling (faster)",
    )
    pipeline_group.add_argument(
        "--kraken2-db", default=None,
        help="Path to Kraken2 database (or set $KRAKEN2_DB)",
    )
    pipeline_group.add_argument(
        "--humann-db", default=None,
        help="Path to HUMAnN3 databases (or set $HUMANN_DB)",
    )
    pipeline_group.add_argument(
        "--read-length", type=int, default=DEFAULT_READ_LENGTH,
        help=f"Read length for Bracken (default: {DEFAULT_READ_LENGTH})",
    )

    args = parser.parse_args()

    if args.demo:
        run_demo(Path(args.output))
    elif args.r1 or args.input:
        run_pipeline(args)
    else:
        parser.print_help()
        print("\nERROR: Provide --r1/--r2, --input, or --demo.", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()
