"""Shared constants and pipeline logic for the claw_metagenomics stage scripts.

Ported faithfully from ClawBio 1.x's metagenomics_profiler.py
(https://github.com/ClawBio/ClawBio/tree/main/skills/claw-metagenomics). The
science (Kraken2/Bracken taxonomy, RGI/CARD resistome, WHO priority
classification, HUMAnN3 functional profiling, figure generation) is
unchanged. What changed is the CLI surface: the 1.x script was one
monolithic `run_pipeline()`/`run_demo()` entry point; here the same logic is
split across five stage scripts (taxonomy.py, resistome.py, functional.py,
figures.py, report.py) so each maps onto one clawbio pipeline.yaml step, with
intermediate results handed off as files instead of in-memory DataFrames.
The internal reproducibility-bundle writer from 1.x is dropped — ClawBio 2.0's
own `clawbio bundle` verb now owns that job.
"""
from __future__ import annotations

import os
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any, Dict, List, Optional

import numpy as np
import pandas as pd

# ---------------------------------------------------------------------------
# Constants (faithful to 1.x)
# ---------------------------------------------------------------------------

KRAKEN2_CONFIDENCE = 0.2
BRACKEN_THRESHOLD = 10       # minimum reads for Bracken re-estimation
BRACKEN_LEVEL = "S"          # species
DEFAULT_READ_LENGTH = 150
RGI_CRITERIA = ("Perfect", "Strict")

# Sentinel filenames: passing the bundled demo fixture as --r1/--r2 selects
# demo mode for every stage, mirroring 1.x's --demo flag but decided in
# Python (per-stage, from the input paths) rather than branching in YAML.
DEMO_R1_NAME = "demo_R1.fastq.gz"
DEMO_R2_NAME = "demo_R2.fastq.gz"

# Sentinel values for optional external databases.
SENTINEL_DEMO = "demo"
SENTINEL_SKIP = "skip"

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

def is_demo_input(r1: Path, r2: Optional[Path]) -> bool:
    """True if r1/r2 are the bundled demo fixtures (data/demo_R{1,2}.fastq.gz).

    All three profiling stages call this on their own --r1/--r2 args, so
    demo-vs-real is decided once, consistently, in Python — not via a
    separate flag threaded through pipeline.yaml.
    """
    if r1.name != DEMO_R1_NAME:
        return False
    if r2 is not None and r2.name != DEMO_R2_NAME:
        return False
    return True


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


def write_status_table(path: Path, status: str) -> None:
    """Write a placeholder TSV carrying only a `_status` sentinel column.

    Used when a stage cannot produce real data (tool missing, run failed,
    functional profiling explicitly skipped) so the step still emits the
    file clawbio's `done --output` contract requires — downstream stages
    (figures, report) check for `_status` to render a graceful section
    instead of crashing on missing columns.
    """
    df = pd.DataFrame({"_status": [status]})
    df.to_csv(path, sep="\t", index=False)


def read_status(df: pd.DataFrame) -> Optional[str]:
    """Return the `_status` sentinel value if this table is a placeholder."""
    if "_status" in df.columns and len(df.columns) == 1:
        return str(df["_status"].iloc[0])
    return None


# ---------------------------------------------------------------------------
# Stage: Kraken2 + Bracken taxonomy
# ---------------------------------------------------------------------------

def run_kraken2(
    r1: Path,
    r2: Optional[Path],
    work_dir: Path,
    db_path: Path,
    threads: int,
    confidence: float = KRAKEN2_CONFIDENCE,
) -> Path:
    """Run Kraken2 taxonomic classification."""
    report_path = work_dir / "kraken2_report.txt"
    output_path = work_dir / "kraken2_output.txt"

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
    work_dir: Path,
    db_path: Path,
    read_length: int = DEFAULT_READ_LENGTH,
    level: str = BRACKEN_LEVEL,
    threshold: int = BRACKEN_THRESHOLD,
) -> Path:
    """Run Bracken abundance re-estimation."""
    bracken_output = work_dir / f"bracken_{level.lower()}.tsv"

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
# Stage: RGI resistome profiling
# ---------------------------------------------------------------------------

def run_rgi(
    r1: Path,
    r2: Optional[Path],
    work_dir: Path,
    threads: int,
) -> Path:
    """Run RGI MAIN for antimicrobial resistance gene detection."""
    # RGI operates on a single concatenated input or an assembly
    # For read-level: we use rgi bwt (read-based)
    rgi_output = work_dir / "rgi_results"

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
    """Parse RGI results into a DataFrame, filtered to Perfect/Strict hits."""
    df = pd.read_csv(rgi_path, sep="\t")
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


def demo_who_classification() -> pd.DataFrame:
    """WHO classification for the embedded demo resistome table."""
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
    return pd.DataFrame(who_records)


# ---------------------------------------------------------------------------
# Stage: HUMAnN3 functional profiling
# ---------------------------------------------------------------------------

def run_humann3(
    r1: Path,
    r2: Optional[Path],
    work_dir: Path,
    threads: int,
    db_path: Optional[Path] = None,
) -> Optional[Path]:
    """Run HUMAnN3 for functional pathway profiling.

    Returns the path to the pathway abundance file, or None if no output
    was produced.
    """
    humann_dir = work_dir / "humann3"
    humann_dir.mkdir(parents=True, exist_ok=True)

    # HUMAnN3 expects a single concatenated input
    if r2 is not None:
        concat_path = work_dir / "concat_reads.fastq.gz"
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
        # Nothing plottable — still leave a valid (empty) figure file so the
        # step's declared emit exists.
        fig, ax = plt.subplots(figsize=(6, 2))
        ax.text(0.5, 0.5, "No taxonomy data available", ha="center", va="center")
        ax.axis("off")

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

    heat_df = None
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

    fig, ax = plt.subplots(
        figsize=(12, max(6, len(heat_df) * 0.35)) if heat_df is not None and len(heat_df) else (8, 4)
    )
    if heat_df is not None and not heat_df.empty:
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
    else:
        ax.text(0.5, 0.5, "No resistome data available", ha="center", va="center")
        ax.axis("off")

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

    fig, axes = plt.subplots(1, 2, figsize=(14, 6), gridspec_kw={"width_ratios": [1, 1.5]})

    if who_df.empty:
        for ax in axes:
            ax.text(0.5, 0.5, "No WHO-priority data available", ha="center", va="center")
            ax.axis("off")
        plt.tight_layout()
        fig.savefig(output_path, dpi=300, bbox_inches="tight")
        plt.close(fig)
        print(f"  Saved: {output_path.name}")
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
