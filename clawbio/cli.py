#!/usr/bin/env python3
"""
clawbio.py — ClawBio Bioinformatics Skills Runner
Standalone CLI and importable module for running ClawBio skills.

Usage:
    python clawbio.py list
    python clawbio.py run pharmgx --demo
    python clawbio.py run equity --input data.vcf
    python clawbio.py run pharmgx --input patient.txt --output ./results
    python clawbio.py upload --input patient.txt --patient-id PT001
    python clawbio.py run pharmgx --profile profiles/PT001.json --output ./results
    python clawbio.py run full-profile --profile profiles/PT001.json --output ./results

Importable:
    # With the repository checkout on sys.path:
    from clawbio import run_skill, list_skills, upload_profile
    result = run_skill("pharmgx", demo=True)
"""

import argparse
import json
import os
import re
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

from clawbio import __version__
from clawbio.contract_alerts import append_contract_alert_log, normalise_contract_alerts

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

# Resolve the ClawBio content root in two layouts:
#   - dev checkout:    skills/ lives at the repository root (one level above this
#                      package), so CLAWBIO_DIR points there.
#   - installed wheel: skills/ and examples/ are bundled inside the package, so
#                      CLAWBIO_DIR points at the package directory itself.
_PKG_DIR = Path(__file__).resolve().parent
_INSTALLED = (_PKG_DIR / "skills").is_dir()
CLAWBIO_DIR = _PKG_DIR if _INSTALLED else _PKG_DIR.parent
SKILLS_DIR = CLAWBIO_DIR / "skills"
EXAMPLES_DIR = CLAWBIO_DIR / "examples"
# Read-only package data must never be written to once installed; route output
# and patient profiles to the current working directory in that case.
_WRITE_ROOT = Path.cwd() if _INSTALLED else CLAWBIO_DIR
DEFAULT_OUTPUT_ROOT = _WRITE_ROOT / "output"
PROFILES_DIR = _WRITE_ROOT / "profiles"

# Python binary — use the same interpreter that launched the ClawBio CLI
PYTHON = sys.executable

SCRNASEQ_BACKEND_PROFILES = {
    "docker",
    "conda",
    "mamba",
    "singularity",
    "apptainer",
    "podman",
    "shifter",
    "charliecloud",
    "wave",
    "gpu",
    "debug",
    "arm64",
    "emulate_amd64",
    "test",
    "test_full",
    "test_cellrangermulti",
    "test_multiome",
}


def _is_scrnaseq_backend_profile(value: str | None) -> bool:
    if not value:
        return False
    components = [part.strip() for part in str(value).split(",") if part.strip()]
    return bool(components) and all(
        component in SCRNASEQ_BACKEND_PROFILES or re.fullmatch(r"[A-Za-z0-9_.-]+", component)
        for component in components
    )

# --------------------------------------------------------------------------- #
# ANSI color support
# --------------------------------------------------------------------------- #

def _use_color() -> bool:
    return sys.stdout.isatty() and os.environ.get("NO_COLOR") is None

_COLOR = _use_color()

BOLD    = "\033[1m"  if _COLOR else ""
DIM     = "\033[2m"  if _COLOR else ""
RED     = "\033[31m" if _COLOR else ""
GREEN   = "\033[32m" if _COLOR else ""
YELLOW  = "\033[33m" if _COLOR else ""
CYAN    = "\033[36m" if _COLOR else ""
WHITE   = "\033[37m" if _COLOR else ""
BG_RED  = "\033[41m" if _COLOR else ""
RESET   = "\033[0m"  if _COLOR else ""


def colorize_report_line(line: str) -> str:
    """Apply ANSI color to a report line based on clinical significance."""
    stripped = line.strip()
    if not stripped:
        return line
    if stripped.startswith("#"):
        return f"{CYAN}{BOLD}{line}{RESET}"
    upper = stripped.upper()
    # Special: warfarin + avoid → red background
    if "WARFARIN" in upper and "AVOID" in upper:
        return f"{BG_RED}{WHITE}{BOLD}{line}{RESET}"
    if "AVOID" in upper:
        return f"{RED}{BOLD}{line}{RESET}"
    if "CAUTION" in upper:
        return f"{YELLOW}{line}{RESET}"
    if "STANDARD" in upper or "| OK" in upper or "NORMAL" in upper:
        return f"{GREEN}{line}{RESET}"
    if stripped.startswith("---") or stripped.startswith("===") or stripped.startswith("| ---"):
        return f"{DIM}{line}{RESET}"
    return line


def print_boxed_header(title: str):
    """Print a Unicode rounded-box header."""
    w = len(title) + 4
    print(f"{CYAN}╭{'─' * w}╮{RESET}")
    print(f"{CYAN}│  {BOLD}{title}{RESET}{CYAN}  │{RESET}")
    print(f"{CYAN}╰{'─' * w}╯{RESET}")


def _parse_md_table(text: str, header_start: str) -> list[list[str]]:
    """Extract data rows from a markdown table identified by its header."""
    rows = []
    found = False
    for line in text.splitlines():
        if header_start in line:
            found = True
            continue
        if found:
            if line.strip().startswith("| ---") or line.strip().startswith("|---"):
                continue
            if line.strip().startswith("|") and line.count("|") >= 3:
                rows.append([c.strip() for c in line.split("|")[1:-1]])
            elif rows:
                break
    return rows


def format_pharmgx_preview(report_text: str, report_path: str):
    """Render a rich, biologically insightful pharmgx report for the terminal."""
    lines = report_text.splitlines()

    # --- Extract metadata ---
    meta = {}
    for line in lines:
        for key in ("Pharmacogenomic SNPs found", "Genes profiled",
                     "Drugs assessed", "Input", "Format detected"):
            if f"**{key}**" in line:
                meta[key] = line.split(":", 1)[-1].strip().strip("`* ")

    # --- Extract gene profile rows ---
    gene_rows = _parse_md_table(report_text, "| Gene | Full Name |")

    # --- Extract drug summary rows ---
    summary = {}
    for row in _parse_md_table(report_text, "| Category | Count |"):
        if len(row) >= 2:
            summary[row[0]] = row[1]

    # --- Extract actionable alerts ---
    avoid_drugs, caution_drugs = [], []
    section = None
    for line in lines:
        if "AVOID / USE ALTERNATIVE:" in line:
            section = "avoid"
        elif "USE WITH CAUTION:" in line:
            section = "caution"
        elif line.startswith("---") or (line.startswith("##") and "Actionable" not in line):
            section = None
        elif section and line.strip().startswith("- **"):
            m = re.match(r'- \*\*(.+?)\*\* \((.+?)\) \[(.+?)]: (.+)', line.strip())
            if m:
                entry = {"drug": m[1], "brand": m[2], "genes": m[3], "rec": m[4]}
                (avoid_drugs if section == "avoid" else caution_drugs).append(entry)

    # === RENDER ===
    W = 60
    snps = meta.get("Pharmacogenomic SNPs found", "?")
    n_genes = meta.get("Genes profiled", "?")
    n_drugs = meta.get("Drugs assessed", "?")
    fmt = meta.get("Format detected", "unknown")

    # ── Header ──
    print(f"\n{CYAN}╭{'─' * W}╮{RESET}")
    print(f"{CYAN}│{RESET}  {BOLD}{CYAN}ClawBio PharmGx Report{RESET}"
          f"{' ' * (W - 24)}{CYAN}│{RESET}")
    print(f"{CYAN}│{RESET}  {DIM}Corpasome (CC0) · doi:10.6084/m9.figshare.693052{RESET}"
          f"{' ' * (W - 51)}{CYAN}│{RESET}")
    print(f"{CYAN}╰{'─' * W}╯{RESET}")
    print()
    print(f"  {BOLD}{n_genes}{RESET} genes  {DIM}·{RESET}  "
          f"{BOLD}{snps}{RESET} SNPs  {DIM}·{RESET}  "
          f"{BOLD}{n_drugs}{RESET} drugs  {DIM}·{RESET}  "
          f"{DIM}{fmt} format{RESET}")

    # ── Critical findings ──
    if avoid_drugs:
        print(f"\n  {BG_RED}{WHITE}{BOLD} {'▲ CRITICAL FINDING':^{W - 4}} {RESET}")
        print(f"  {RED}{'─' * W}{RESET}")
        for d in avoid_drugs:
            print(f"    {RED}{BOLD}{d['drug']}{RESET} ({d['brand']})  "
                  f"{DIM}[{d['genes']}]{RESET}")
            if d["drug"].lower() == "warfarin":
                print()
                print(f"    {YELLOW}{BOLD}VKORC1{RESET}{YELLOW} rs9923231 {BOLD}TT{RESET}"
                      f"  {DIM}→{RESET}  Both copies carry the sensitivity allele.")
                print(f"    {DIM}This patient produces less vitamin K epoxide reductase,{RESET}")
                print(f"    {DIM}making them hyper-responsive to warfarin's mechanism.{RESET}")
                print()
                print(f"    {YELLOW}{BOLD}CYP2C9{RESET}{YELLOW} *1/*2 {DIM}(rs1799853 CT){RESET}"
                      f"  {DIM}→{RESET}  Intermediate Metabolizer.")
                print(f"    {DIM}Warfarin is cleared ~40% more slowly than in *1/*1 carriers,{RESET}")
                print(f"    {DIM}causing the drug to accumulate at standard doses.{RESET}")
                print()
                print(f"    {RED}{BOLD}Combined effect:{RESET}  "
                      f"Standard doses risk {RED}{BOLD}life-threatening bleeding{RESET}.")
                print(f"    CPIC guidelines recommend {BOLD}50–80% dose reduction{RESET} or")
                print(f"    switching to a DOAC (apixaban, rivaroxaban).")
            else:
                print(f"    {d['rec']}")
        print(f"  {RED}{'─' * W}{RESET}")

    # ── Gene profile ──
    print(f"\n  {CYAN}{BOLD}Gene Profile{RESET}")
    print(f"  {DIM}{'─' * (W - 2)}{RESET}")
    for row in gene_rows:
        if len(row) < 4:
            continue
        gene, _, diplotype, phenotype = row[:4]
        # Split off "(X/Y SNPs tested)" qualifier from diplotype for cleaner display
        dip_match = re.match(r'^(.+?)\s*(\(\d/\d SNPs tested\))?$', diplotype)
        dip_core = dip_match[1] if dip_match else diplotype
        dip_note = f" {DIM}{dip_match[2]}{RESET}" if dip_match and dip_match[2] else ""
        # Choose color by phenotype category
        if "Unknown" in phenotype or "unmapped" in phenotype:
            pc = YELLOW
            phenotype_short = "Unknown"
            extra = f"  {DIM}(needs clinical testing){RESET}"
        elif "High" in phenotype:
            pc, phenotype_short, extra = RED, phenotype, ""
        elif "Poor" in phenotype:
            pc, phenotype_short, extra = RED, phenotype, ""
        elif "Intermediate" in phenotype:
            pc, phenotype_short, extra = YELLOW, "Intermediate", ""
        elif "Non-expressor" in phenotype:
            pc, phenotype_short, extra = DIM, "Non-expressor", ""
        else:
            pc, phenotype_short, extra = GREEN, "Normal", ""
        wmark = f"  {RED}← warfarin{RESET}" if gene in ("CYP2C9", "VKORC1") else ""
        print(f"  {BOLD}{gene:<10}{RESET} {DIM}{dip_core:<12}{RESET}"
              f" {pc}{phenotype_short}{RESET}{extra}{dip_note}{wmark}")

    # ── Drug summary ──
    print(f"\n  {CYAN}{BOLD}Drug Summary{RESET}")
    print(f"  {DIM}{'─' * (W - 2)}{RESET}")
    buckets = [
        ("Avoid / use alternative", RED,    BOLD),
        ("Use with caution",       YELLOW, ""),
        ("Standard dosing",        GREEN,  ""),
        ("Insufficient data",      DIM,    ""),
    ]
    for cat, color, bld in buckets:
        count = summary.get(cat, "0")
        b = BOLD if bld else ""
        print(f"  {color}{b}■{RESET}  {color}{count:>2} {cat}{RESET}")

    # ── Caution list ──
    if caution_drugs:
        print()
        names = [f"{YELLOW}{BOLD}{d['drug']}{RESET}" for d in caution_drugs]
        print(f"  {YELLOW}Caution:{RESET} {f'{DIM}, {RESET}'.join(names)}")

    # ── Footer ──
    print(f"\n  {DIM}Full report → {report_path}{RESET}")
    print(f"  {DIM}Disclaimer: research/educational use only — not a medical device{RESET}")
    print(f"{BOLD}{'━' * W}{RESET}")


# --------------------------------------------------------------------------- #
# Skills registry
# --------------------------------------------------------------------------- #

SKILLS = {
    "pharmgx": {
        "script": SKILLS_DIR / "pharmgx-reporter" / "pharmgx_reporter.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "pharmgx-reporter" / "demo_patient.txt"),
        ],
        "description": "Pharmacogenomics reporter (12 genes, 31 SNPs, 51 drugs)",
        "allowed_extra_flags": {"--weights"},
        "api_module": "skills.pharmgx-reporter.api",
        "accepts_genotypes": True,
    },
    "equity": {
        "script": SKILLS_DIR / "equity-scorer" / "equity_scorer.py",
        "demo_args": [
            "--input",
            str(EXAMPLES_DIR / "demo_populations.vcf"),
            "--pop-map",
            str(EXAMPLES_DIR / "demo_population_map.csv"),
        ],
        "description": "HEIM equity scorer (FST, heterozygosity, population representation)",
        "allowed_extra_flags": {"--weights", "--pop-map"},
        "accepts_genotypes": False,  # needs VCF/CSV file, not genotype dict
    },
    "nutrigx": {
        "script": SKILLS_DIR / "nutrigx" / "nutrigx.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "nutrigx" / "tests" / "synthetic_patient.csv"),
        ],
        "description": "Nutrigenomics advisor (diet, vitamins, caffeine, lactose)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": True,
    },
    "dnasp": {
        "script": SKILLS_DIR / "dnasp" / "dnasp.py",
        "demo_args": ["--demo"],
        "description": "DnaSP 6 population genetics (Pi, Tajima's D, Fu & Li, Fay & Wu, MK, Ka/Ks, Fst, and more)",
        "allowed_extra_flags": {
            "--fasta", "--outgroup", "--pop-map", "--window", "--step",
            "--all", "--pi", "--theta", "--tajima", "--fuliD", "--fuliF",
            "--hka", "--mk", "--kaks", "--r2", "--fufs", "--sfs",
            "--tstv", "--codon", "--faywu", "--fst",
            "--n-sim", "--sim-seed",
        },
        "accepts_genotypes": False,
    },
    "metagenomics": {
        "script": SKILLS_DIR / "claw-metagenomics" / "metagenomics_profiler.py",
        "demo_args": ["--demo"],
        "description": "Metagenomics profiler (Kraken2, RGI/CARD, HUMAnN3)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": False,
    },
    "analyze-fasta": {
        "script": SKILLS_DIR / "analyze-fasta" / "analyze_fasta.py",
        "demo_args": ["--demo"],
        "description": "Single FASTA analyzer (auto-detect nucleotide/protein, GC, ORFs, MW, pI, GRAVY)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": False,
    },
    "phylo": {
        "script": SKILLS_DIR / "phylogenetics-builder" / "phylogenetics_builder.py",
        "demo_args": ["--demo"],
        "description": "Build maximum-likelihood phylogenetic trees from aligned FASTA data using IQ-TREE 2",
        "allowed_extra_flags": set(),
        "accepts_genotypes": False,
    },
    "scrnaQ3�M4T4 =w��r�����ם��h_8��m�֥                ("umi_dedup_tool", "--umi-dedup-tool"),
                ("umitools_extract_method", "--umitools-extract-method"),
                ("umitools_bc_pattern", "--umitools-bc-pattern"),
                ("umitools_bc_pattern2", "--umitools-bc-pattern2"),
                ("umitools_umi_separator", "--umitools-umi-separator"),
                ("umitools_grouping_method", "--umitools-grouping-method"),
                ("salmon_quant_libtype", "--salmon-quant-libtype"),
                ("extra_star_align_args", "--extra-star-align-args"),
                ("extra_bowtie2_align_args", "--extra-bowtie2-align-args"),
                ("extra_salmon_quant_args", "--extra-salmon-quant-args"),
                ("extra_kallisto_quant_args", "--extra-kallisto-quant-args"),
                ("rsem_extra_args", "--rsem-extra-args"),
                ("contaminant_screening", "--contaminant-screening"),
                ("contaminant_screening_input", "--contaminant-screening-input"),
                ("kraken_db", "--kraken-db"),
                ("bracken_precision", "--bracken-precision"),
                ("sylph_db", "--sylph-db"),
                ("sylph_taxonomy", "--sylph-taxonomy"),
                ("bbsplit_fasta_list", "--bbsplit-fasta-list"),
                ("bbsplit_index", "--bbsplit-index"),
                ("rseqc_modules", "--rseqc-modules"),
                ("gtf_extra_attributes", "--gtf-extra-attributes"),
                ("gtf_group_features", "--gtf-group-features"),
                ("featurecounts_group_type", "--featurecounts-group-type"),
                ("featurecounts_feature_type", "--featurecounts-feature-type"),
                ("seq_platform", "--seq-platform"),
                ("multiqc_config", "--multiqc-config"),
                ("multiqc_logo", "--multiqc-logo"),
                ("multiqc_methods_description", "--multiqc-methods-description"),
                ("email_on_fail", "--email-on-fail"),
                ("igenomes_base", "--igenomes-base"),
                ("sortmerna_index", "--sortmerna-index"),
                ("ribo_database_manifest", "--ribo-database-manifest"),
                ("hisat2_build_memory", "--hisat2-build-memory"),
                ("gpu_container_options", "--gpu-container-options"),
                ("extra_fqlint_args", "--extra-fqlint-args"),
                ("publish_dir_mode", "--publish-dir-mode"),
                ("metadata", "--metadata"),
                ("formula", "--formula"),
                ("contrast", "--contrast"),
                ("downstream_output", "--downstream-output"),
            ):
                attr, flag = value_flag
                value = getattr(args, attr, None)
                if value:
                    extra.extend([flag, value])
            for int_flag in (
                ("pseudo_aligner_kmer_size", "--pseudo-aligner-kmer-size"),
                ("min_trimmed_reads", "--min-trimmed-reads"),
                ("umi_discard_read", "--umi-discard-read"),
                ("kallisto_quant_fraglen", "--kallisto-quant-fraglen"),
                ("kallisto_quant_fraglen_sd", "--kallisto-quant-fraglen-sd"),
            ):
                attr, flag = int_flag
                value = getattr(args, attr, None)
                if value is not None:
                    extra.extend([flag, str(value)])
            for float_flag in (
                ("min_mapped_reads", "--min-mapped-reads"),
                ("stranded_threshold", "--stranded-threshold"),
                ("unstranded_threshold", "--unstranded-threshold"),
            ):
                attr, flag = float_flag
                value = getattr(args, attr, None)
                if value is not None:
                    extra.extend([flag, str(value)])
            for bool_flag in (
                ("prokaryotic", "--prokaryotic"),
                ("rapid_quant", "--rapid-quant"),
                ("arm", "--arm"),
                ("gencode", "--gencode"),
                ("remove_ribo_rna", "--remove-ribo-rna"),
                ("with_umi", "--with-umi"),
                ("skip_umi_extract", "--skip-umi-extract"),
                ("umitools_dedup_stats", "--umitools-dedup-stats"),
                ("umitools_dedup_primary_only", "--umitools-dedup-primary-only"),
                ("bam_csi_index", "--bam-csi-index"),
                ("stringtie_ignore_gtf", "--stringtie-ignore-gtf"),
                ("gffread_transcript_fasta", "--gffread-transcript-fasta"),
                ("use_rustqc", "--use-rustqc"),
                ("use_parabricks_star", "--use-parabricks-star"),
                ("use_sentieon_star", "--use-sentieon-star"),
                ("use_gpu_ribodetector", "--use-gpu-ribodetector"),
                ("skip_trimming", "--skip-trimming"),
                ("skip_alignment", "--skip-alignment"),
                ("skip_pseudo_alignment", "--skip-pseudo-alignment"),
                ("skip_quantification_merge", "--skip-quantification-merge"),
                ("skip_markduplicates", "--skip-markduplicates"),
                ("skip_bigwig", "--skip-bigwig"),
                ("skip_stringtie", "--skip-stringtie"),
                ("skip_dupradar", "--skip-dupradar"),
                ("skip_qualimap", "--skip-qualimap"),
                ("skip_rseqc", "--skip-rseqc"),
                ("skip_biotype_qc", "--skip-biotype-qc"),
                ("skip_deseq2_qc", "--skip-deseq2-qc"),
                ("enable_preseq", "--enable-preseq"),
                ("skip_qc", "--skip-qc"),
                ("save_trimmed", "--save-trimmed"),
                ("save_unaligned", "--save-unaligned"),
                ("save_merged_fastq", "--save-merged-fastq"),
                ("save_non_ribo_reads", "--save-non-ribo-reads"),
                ("save_umi_intermeds", "--save-umi-intermeds"),
                ("save_kraken_assignments", "--save-kraken-assignments"),
                ("save_kraken_unassigned", "--save-kraken-unassigned"),
                ("skip_bbsplit", "--skip-bbsplit"),
                ("save_bbsplit_reads", "--save-bbsplit-reads"),
                ("skip_linting", "--skip-linting"),
                ("skip_gtf_filter", "--skip-gtf-filter"),
                ("skip_gtf_transcript_filter", "--skip-gtf-transcript-filter"),
            ):
                attr, flag = bool_flag
                if getattr(args, attr, False):
                    extra.append(flag)
            # --deseq2-vst is tri-state: True (force on), False (--no-deseq2-vst),
            # None (default — don't forward either direction).
            if getattr(args, "deseq2_vst", None) is False:
                extra.append("--no-deseq2-vst")
            elif getattr(args, "deseq2_vst", None) is True:
                extra.append("--deseq2-vst")
        if args.skill == "sarek-pipeline":
            for value_flag in (
                ("pipeline_local", "--pipeline-local"),
                ("aligner", "--aligner"),
                ("seq_platform", "--seq-platform"),
                ("email_on_fail", "--email-on-fail"),
                ("multiqc_config", "--multiqc-config"),
                ("multiqc_logo", "--multiqc-logo"),
                ("multiqc_methods_description", "--multiqc-methods-description"),
                ("igenomes_base", "--igenomes-base"),
                ("bbsplit_fasta_list", "--bbsplit-fasta-list"),
                ("bbsplit_index", "--bbsplit-index"),
                ("publish_dir_mode", "--publish-dir-mode"),
            ):
                attr, flag = value_flag
                value = getattr(args, attr, None)
                if value:
                    extra.extend([flag, value])
            for attr, flag in (
                ("arm", "--arm"),
                ("save_trimmed", "--save-trimmed"),
                ("save_bbsplit_reads", "--save-bbsplit-reads"),
            ):
                if getattr(args, attr, False):
                    extra.append(flag)
        if getattr(args, "drug", None):
            extra.extend(["--drug", args.drug])
        if getattr(args, "dose", None):
            extra.extend(["--dose", args.dose])
        if getattr(args, "trait", None):
            extra.extend(["--trait", args.trait])
        if getattr(args, "pgs_id", None):
            extra.extend(["--pgs-id", args.pgs_id])
        if getattr(args, "gene", None):
            extra.extend(["--gene", args.gene])
        if getattr(args, "genes", None):
            extra.extend(["--genes", args.genes])
        if getattr(args, "rsid", None):
            extra.extend(["--rsid", args.rsid])
        if getattr(args, "skip", None):
            extra.extend(["--skip", args.skip])
        if getattr(args, "query", None):
            extra.extend(["--query", args.query])
        if getattr(args, "location", None):
            extra.extend(["--location", args.location])
        if getattr(args, "max_rows", None) is not None:
            extra.extend(["--max-rows", str(args.max_rows)])
        if getattr(args, "max_bytes_billed", None) is not None:
            extra.extend(["--max-bytes-billed", str(args.max_bytes_billed)])
        if getattr(args, "param", None):
            for param in args.param:
                extra.extend(["--param", param])
        if getattr(args, "dry_run", False):
            extra.append("--dry-run")
        if getattr(args, "list_datasets", None):
            extra.extend(["--list-datasets", args.list_datasets])
        if getattr(args, "list_tables", None):
            extra.extend(["--list-tables", args.list_tables])
        if getattr(args, "describe", None):
            extra.extend(["--describe", args.describe])
        if getattr(args, "preview", None) is not None:
            extra.extend(["--preview", str(args.preview)])
        if getattr(args, "count_only", False):
            extra.append("--count-only")
        if getattr(args, "paper", None):
            extra.extend(["--paper", args.paper])
        if getattr(args, "note", None):
            for note in args.note:
                extra.extend(["--note", note])
        if getattr(args, "geo_id", None):
            extra.extend(["--geo-id", args.geo_id])
        if getattr(args, "clocks", None):
            extra.extend(["--clocks", args.clocks])
        if getattr(args, "metadata_cols", None):
            extra.extend(["--metadata-cols", args.metadata_cols])
        if getattr(args, "imputer_strategy", None):
            extra.extend(["--imputer-strategy", args.imputer_strategy])
        if getattr(args, "skip_epicv2_aggregation", False):
            extra.append("--skip-epicv2-aggregation")
        if getattr(args, "verbose", False):
            extra.append("--verbose")
        if getattr(args, "vcf", None):
            extra.extend(["--vcf", args.vcf])
        if getattr(args, "qc", None):
            extra.extend(["--qc", args.qc])
        if getattr(args, "sample_sheet", None):
            extra.extend(["--sample-sheet", args.sample_sheet])
        if getattr(args, "metadata_provider", None):
            extra.extend(["--metadata-provider", args.metadata_provider])
        if getattr(args, "ica_project_id", None):
            extra.extend(["--ica-project-id", args.ica_project_id])
        if getattr(args, "ica_run_id", None):
            extra.extend(["--ica-run-id", args.ica_run_id])
        if getattr(args, "counts", None):
            extra.extend(["--counts", args.counts])
        if args.skill != "rnaseq-pipeline" and getattr(args, "metadata", None):
            extra.extend(["--metadata", args.metadata])
        if args.skill != "rnaseq-pipeline" and getattr(args, "formula", None):
            extra.extend(["--formula", args.formula])
        if args.skill != "rnaseq-pipeline" and getattr(args, "contrast", None):
            extra.extend(["--contrast", args.contrast])
        if getattr(args, "backend", None):
            extra.extend(["--backend", args.backend])
        if getattr(args, "min_count", None) is not None:
            extra.extend(["--min-count", str(args.min_count)])
        if getattr(args, "min_samples", None) is not None:
            extra.extend(["--min-samples", str(args.min_samples)])
        if getattr(args, "mode", None):
            extra.extend(["--mode", args.mode])
        if getattr(args, "adata", None):
            extra.extend(["--adata", args.adata])
        if getattr(args, "top_genes", None) is not None:
            extra.extend(["--top-genes", str(args.top_genes)])
        if getattr(args, "label_top", None) is not None:
            extra.extend(["--label-top", str(args.label_top)])
        if getattr(args, "padj_threshold", None) is not None:
            extra.extend(["--padj-threshold", str(args.padj_threshold)])
        if getattr(args, "lfc_threshold", None) is not None:
            extra.extend(["--lfc-threshold", str(args.lfc_threshold)])
        if getattr(args, "min_basemean", None) is not None:
            extra.extend(["--min-basemean", str(args.min_basemean)])
        if getattr(args, "method", None):
            extra.extend(["--method", args.method])
        if getattr(args, "layer", None):
            extra.extend(["--layer", args.layer])
        if getattr(args, "batch_key", None):
            extra.extend(["--batch-key", args.batch_key])
        if getattr(args, "labels_key", None):
            extra.extend(["--labels-key", args.labels_key])
        if getattr(args, "unlabeled_category", None):
            extra.extend(["--unlabeled-category", args.unlabeled_category])
        if getattr(args, "min_genes", None) is not None:
            extra.extend(["--min-genes", str(args.min_genes)])
        if getattr(args, "min_cells", None) is not None:
            extra.extend(["--min-cells", str(args.min_cells)])
        if getattr(args, "max_mt_pct", None) is not None:
            extra.extend(["--max-mt-pct", str(args.max_mt_pct)])
        if getattr(args, "n_top_hvg", None) is not None:
            extra.extend(["--n-top-hvg", str(args.n_top_hvg)])
        if getattr(args, "n_pcs", None) is not None:
            extra.extend(["--n-pcs", str(args.n_pcs)])
        if getattr(args, "latent_dim", None) is not None:
            extra.extend(["--latent-dim", str(args.latent_dim)])
        if getattr(args, "max_epochs", None) is not None:
            extra.extend(["--max-epochs", str(args.max_epochs)])
        if getattr(args, "n_neighbors", None) is not None:
            extra.extend(["--n-neighbors", str(args.n_neighbors)])
        if getattr(args, "use_rep", None):
            extra.extend(["--use-rep", args.use_rep])
        if getattr(args, "leiden_resolution", None) is not None:
            extra.extend(["--leiden-resolution", str(args.leiden_resolution)])
        if getattr(args, "random_state", None) is not None:
            extra.extend(["--random-state", str(args.random_state)])
        if getattr(args, "top_markers", None) is not None:
            extra.extend(["--top-markers", str(args.top_markers)])
        if getattr(args, "accelerator", None):
            extra.extend(["--accelerator", args.accelerator])
        if getattr(args, "contrast_groupby", None):
            extra.extend(["--contrast-groupby", args.contrast_groupby])
        if getattr(args, "contrast_scope", None):
            extra.extend(["--contrast-scope", args.contrast_scope])
        if getattr(args, "contrast_clusterby", None):
            extra.extend(["--contrast-clusterby", args.contrast_clusterby])
        if getattr(args, "contrast_top_genes", None) is not None:
            extra.extend(["--contrast-top-genes", str(args.contrast_top_genes)])
        if getattr(args, "doublet_method", None):
            extra.extend(["--doublet-method", args.doublet_method])
        if getattr(args, "annotate", None):
            extra.extend(["--annotate", args.annotate])
        if getattr(args, "annotation_model", None):
            extra.extend(["--annotation-model", args.annotation_model])
        if getattr(args, "search", None):
            extra.extend(["--search", args.search])
        if getattr(args, "recommend", None):
            extra.extend(["--recommend", args.recommend])
        if getattr(args, "workflow", None):
            extra.extend(["--workflow", args.workflow])
        if getattr(args, "package_details", None):
            extra.extend(["--package-details", args.package_details])
        if getattr(args, "docs_search", None):
            extra.extend(["--docs-search", args.docs_search])
        if getattr(args, "package_docs", None):
            extra.extend(["--package-docs", args.package_docs])
        if getattr(args, "list_domains", False):
            extra.append("--list-domains")
        if getattr(args, "setup", False):
            extra.append("--setup")
        if getattr(args, "install", None):
            extra.extend(["--install", args.install])
        if getattr(args, "skill_format", None):
            extra.extend(["--format", args.skill_format])
        if getattr(args, "container", None):
            extra.extend(["--container", args.container])
        if getattr(args, "modality", None):
            extra.extend(["--modality", args.modality])
        if getattr(args, "max_results", None) is not None:
            extra.extend(["--max-results", str(args.max_results)])
        # flow-bio skill flags
        if getattr(args, "flow_search", None):
            extra.extend(["--search", args.flow_search])
        if getattr(args, "pipelines", False):
            extra.append("--pipelines")
        if getattr(args, "samples", False):
            extra.append("--samples")
        if getattr(args, "projects", False):
            extra.append("--projects")
        if getattr(args, "executions", False):
            extra.append("--executions")
        if getattr(args, "organisms", False):
            extra.append("--organisms")
        if getattr(args, "sample_types", False):
            extra.append("--sample-types")
        if getattr(args, "data", False):
            extra.append("--data")
        if getattr(args, "metadata_attributes", False):
            extra.append("--metadata-attributes")
        if getattr(args, "search_samples", None):
            extra.append("--search-samples")
            extra.extend(args.search_samples)
        if getattr(args, "upload_sample", False):
            extra.append("--upload-sample")
        if getattr(args, "name", None):
            extra.extend(["--name", args.name])
        if getattr(args, "reads1", None):
            extra.extend(["--reads1", args.reads1])
        if getattr(args, "reads2", None):
            extra.extend(["--reads2", args.reads2])
        if getattr(args, "organism", None):
            extra.extend(["--organism", args.organism])
        if getattr(args, "project", None):
            extra.extend(["--project", args.project])
        if getattr(args, "run_pipeline", None):
            extra.extend(["--run-pipeline", args.run_pipeline])
        if getattr(args, "run_samples", None):
            extra.extend(["--run-samples", args.run_samples])
        if getattr(args, "run_data", None):
            extra.extend(["--run-data", args.run_data])
        if getattr(args, "run_params", None):
            extra.extend(["--run-params", args.run_params])
        if getattr(args, "genome", None):
            extra.extend(["--genome", args.genome])
        if getattr(args, "pipeline_detail", None):
            extra.extend(["--pipeline", args.pipeline_detail])
        if getattr(args, "sample_detail", None):
            extra.extend(["--sample", args.sample_detail])
        if getattr(args, "execution_detail", None):
            extra.extend(["--execution", args.execution_detail])
        if getattr(args, "json", False):
            extra.append("--json")

        run_timeout = args.timeout
        if args.timeout == 300:
            run_timeout = SKILLS.get(args.skill, {}).get("default_timeout_seconds", args.timeout)

        result = run_skill(
            skill_name=args.skill,
            input_path=args.input_path,
            output_dir=args.output_dir,
            demo=args.demo,
            extra_args=extra or None,
            timeout=run_timeout,
            profile_path=getattr(args, "profile_path", None),
        )

        # Summary mode: skill printed text to stdout — relay it directly
        if result["output_dir"] is None and result["success"] and result["stdout"]:
            print(result["stdout"], end="")
            sys.exit(0)

        print()
        if result["success"]:
            print(f"  {GREEN}{BOLD}Status:   OK{RESET} {DIM}(exit {result['exit_code']}){RESET}")
        else:
            print(f"  {RED}{BOLD}Status:   FAILED{RESET} {DIM}(exit {result['exit_code']}){RESET}")
        print(f"  {DIM}Duration: {result['duration_seconds']}s{RESET}")
        if result["output_dir"]:
            print(f"  Output:   {result['output_dir']}")
        if result["files"]:
            print(f"  Files:    {', '.join(result['files'])}")
        # Show a preview of the report if one was generated
        if result["success"] and result["output_dir"]:
            report = Path(result["output_dir"]) / "report.md"
            if report.exists():
                text = report.read_text()
                if args.skill == "pharmgx":
                    format_pharmgx_preview(text, str(report))
                else:
                    lines = text.splitlines()
                    print()
                    print_boxed_header("Report Preview")
                    for ln in lines[:40]:
                        print(colorize_report_line(ln))
                    remaining = max(0, len(lines) - 40)
                    if remaining:
                        print(f"\n  {DIM}... ({remaining} more lines in {report}){RESET}")
                    print(f"{BOLD}{'━' * 60}{RESET}")
        if not result["success"] and result["stderr"]:
            print(f"\n  {RED}Error:{RESET}\n{result['stderr'][-800:]}")
        sys.exit(0 if result["success"] else 1)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
