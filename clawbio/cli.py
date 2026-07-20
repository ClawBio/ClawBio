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
    "scrna": {
        "script": SKILLS_DIR / "scrna-orchestrator" / "scrna_orchestrator.py",
        "demo_args": ["--demo"],
        "description": "scRNA Orchestrator (Scanpy QC, doublet detection, clustering, annotation, optional latent downstream mode, dataset-level + within-cluster contrastive markers)",
        "allowed_extra_flags": {
            "--min-genes",
            "--min-cells",
            "--max-mt-pct",
            "--n-top-hvg",
            "--n-pcs",
            "--n-neighbors",
            "--use-rep",
            "--leiden-resolution",
            "--random-state",
            "--top-markers",
            "--contrast-groupby",
            "--contrast-scope",
            "--contrast-clusterby",
            "--contrast-top-genes",
            "--doublet-method",
            "--annotate",
            "--annotation-model",
        },
        "accepts_genotypes": False,
    },
    "scrna-embedding": {
        "script": SKILLS_DIR / "scrna-embedding" / "scrna_embedding.py",
        "demo_args": ["--demo"],
        "description": "scRNA Embedding (scVI/scANVI latent embedding, optional batch integration, stable integrated h5ad export)",
        "allowed_extra_flags": {
            "--method",
            "--layer",
            "--batch-key",
            "--labels-key",
            "--unlabeled-category",
            "--min-genes",
            "--min-cells",
            "--max-mt-pct",
            "--n-top-hvg",
            "--latent-dim",
            "--max-epochs",
            "--n-neighbors",
            "--random-state",
            "--accelerator",
        },
        "accepts_genotypes": False,
    },
    "celltype-specificity-profiler": {
        "script": SKILLS_DIR / "celltype-specificity-profiler" / "profiler.py",
        "demo_args": ["--demo"],
        "description": "Cell-Type Specificity Profiler (per-gene tau specificity index + expression bimodality coefficient from an annotated single-cell atlas)",
        "allowed_extra_flags": {
            "--gene",
            "--atlas",
            "--cell-type-key",
            "--tissue",
            "--trial-prior",
        },
        "accepts_genotypes": False,
    },
    "compare": {
        "script": SKILLS_DIR / "genome-compare" / "genome_compare.py",
        "demo_args": ["--demo"],
        "description": "Genome comparator (IBS vs George Church + ancestry estimation)",
        "allowed_extra_flags": {"--no-figures", "--aims-panel", "--reference"},
        "summary_default": True,
        "accepts_genotypes": True,
    },
    "drugphoto": {
        "script": SKILLS_DIR / "pharmgx-reporter" / "pharmgx_reporter.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "genome-compare" / "data" / "manuel_corpas_23andme.txt.gz"),
        ],
        "description": "Drug photo analysis (single-drug PGx lookup from photo identification)",
        "allowed_extra_flags": {"--drug", "--dose"},
        "summary_default": True,
        "accepts_genotypes": True,
    },
    "prs": {
        "script": SKILLS_DIR / "gwas-prs" / "gwas_prs.py",
        "demo_args": ["--demo"],
        "description": "GWAS Polygenic Risk Score calculator (PGS Catalog, 3000+ scores)",
        "allowed_extra_flags": {"--trait", "--pgs-id", "--min-overlap", "--max-variants", "--build"},
        "accepts_genotypes": True,
    },
    "clinpgx": {
        "script": SKILLS_DIR / "clinpgx" / "clinpgx.py",
        "demo_args": ["--demo"],
        "description": "ClinPGx API query (gene-drug interactions, CPIC guidelines, drug labels)",
        "allowed_extra_flags": {"--gene", "--genes", "--drug", "--drugs", "--no-cache"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "gwas": {
        "script": SKILLS_DIR / "gwas-lookup" / "gwas_lookup.py",
        "demo_args": ["--demo"],
        "description": "GWAS Lookup — federated variant query across 9 genomic databases",
        "allowed_extra_flags": {"--rsid", "--skip", "--no-figures", "--no-cache", "--max-hits"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "bigquery": {
        "script": SKILLS_DIR / "bigquery-public" / "bigquery_public.py",
        "demo_args": ["--demo"],
        "description": "BigQuery Public — read-only SQL bridge for public datasets with local outputs",
        "allowed_extra_flags": {
            "--query",
            "--location",
            "--max-rows",
            "--max-bytes-billed",
            "--param",
            "--dry-run",
            "--list-datasets",
            "--list-tables",
            "--describe",
            "--preview",
            "--count-only",
            "--paper",
            "--note",
        },
        "allowed_extra_flags_without_values": {"--dry-run", "--count-only"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "profile": {
        "script": SKILLS_DIR / "profile-report" / "profile_report.py",
        "demo_args": ["--demo"],
        "description": "Unified personal genomic profile report",
        "allowed_extra_flags": {"--profile"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "galaxy": {
        "script": SKILLS_DIR / "galaxy-bridge" / "galaxy_bridge.py",
        "demo_args": ["--demo"],
        "description": "Galaxy tool discovery and execution (8,000+ bioinformatics tools)",
        "allowed_extra_flags": {"--search", "--list-categories", "--tool-details", "--run", "--max-results"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "bioc": {
        "script": SKILLS_DIR / "bioconductor-bridge" / "bioconductor_bridge.py",
        "demo_args": ["--demo"],
        "description": "Bioconductor package discovery, workflow recommendation, setup, and starter code generation",
        "allowed_extra_flags": {
            "--search",
            "--recommend",
            "--workflow",
            "--package-details",
            "--docs-search",
            "--package-docs",
            "--list-domains",
            "--setup",
            "--install",
            "--format",
            "--modality",
            "--container",
            "--max-results",
        },
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "illumina": {
        "script": SKILLS_DIR / "illumina-bridge" / "illumina_bridge.py",
        "demo_args": ["--demo"],
        "description": "Illumina / DRAGEN bundle import and metadata normalization",
        "allowed_extra_flags": {
            "--vcf",
            "--qc",
            "--sample-sheet",
            "--metadata-provider",
            "--ica-project-id",
            "--ica-run-id",
        },
        "accepts_genotypes": False,
    },
    "data-extract": {
        "script": SKILLS_DIR / "data-extractor" / "data_extractor.py",
        "demo_args": ["--demo"],
        "description": "Extract numerical data from scientific figure images (Claude vision + OpenCV)",
        "allowed_extra_flags": {"--web", "--port", "--plot-type"},
        "api_module": "skills.data-extractor.data_extractor_api",
        "accepts_genotypes": False,
    },
    "rnaseq": {
        "script": SKILLS_DIR / "rnaseq-de" / "rnaseq_de.py",
        "demo_args": ["--demo"],
        "description": "Bulk/pseudo-bulk RNA-seq differential expression (QC + PCA + DE)",
        "allowed_extra_flags": {
            "--counts",
            "--metadata",
            "--formula",
            "--contrast",
            "--backend",
            "--min-count",
            "--min-samples",
        },
    },
    "scrnaseq-pipeline": {
        "script": SKILLS_DIR / "nfcore-scrnaseq-wrapper" / "nfcore_scrnaseq_wrapper.py",
        "demo_args": ["--demo"],
        "description": "Wrapper de preprocessing scRNA FASTQ-to-h5ad vía scrnaseq/Nextflow",
        # Keep the ClawBio runner timeout above the wrapper's internal Nextflow
        # timeout so the wrapper can terminate the process group cleanly first.
        "default_timeout_seconds": 60 * 60 * 12 + 10 * 60,
        "max_output_files_listed": 50,
        "allowed_extra_flags": {
            "--check",
            "-v",
            "--verbose",
            "--no-banner",
            "-c",
            "--config",
            "--nextflow-config",
            "--profile",
            "--pipeline-version",
            "--allow-dirty-pipeline",
            "--require-local-pipeline",
            "--allow-pipeline-version-override",
            "--allow-remote-inputs",
            "--trust-config-params",
            "--preset",
            "--aligner",
            "--protocol",
            "--email",
            "--multiqc-title",
            "--expected-cells",
            "--timeout-hours",
            "--work-dir",
            "--allow-conda-cellranger",
            "--resume",
            "--genome",
            "--igenomes-base",
            "--igenomes-ignore",
            "--email-on-fail",
            "--multiqc-config",
            "--multiqc-logo",
            "--multiqc-methods-description",
            "--publish-dir-mode",
            "--trace-report-suffix",
            "--monochrome-logs",
            "--save-reference",
            "--save-align-intermeds",
            "--no-save-align-intermeds",
            "--skip-cellbender",
            "--skip-fastqc",
            "--skip-emptydrops",
            "--skip-multiqc",
            "--skip-cellranger-renaming",
            "--skip-cellrangermulti-vdjref",
            "--run-downstream",
            "--skip-downstream",
            "--fasta",
            "--gtf",
            "--transcript-fasta",
            "--txp2gene",
            "--simpleaf-index",
            "--simpleaf-umi-resolution",
            "--kallisto-index",
            "--kb-workflow",
            "--kb-t1c",
            "--kb-t2c",
            "--star-index",
            "--star-feature",
            "--star-ignore-sjdbgtf",
            "--seq-center",
            "--cellranger-index",
            "--cellranger-vdj-index",
            "--cellrangerarc-config",
            "--cellrangerarc-reference",
            "--barcode-whitelist",
            "--motifs",
            "--gex-frna-probe-set",
            "--gex-target-panel",
            "--gex-cmo-set",
            "--fb-reference",
            "--vdj-inner-enrichment-primers",
            "--gex-barcode-sample-assignment",
            "--cellranger-multi-barcodes",
        },
        "allowed_extra_flags_without_values": {
            "--check",
            "-v",
            "--verbose",
            "--no-banner",
            "--allow-dirty-pipeline",
            "--require-local-pipeline",
            "--allow-pipeline-version-override",
            "--allow-remote-inputs",
            "--trust-config-params",
            "--resume",
            "--allow-conda-cellranger",
            "--skip-cellbender",
            "--skip-fastqc",
            "--skip-emptydrops",
            "--skip-multiqc",
            "--skip-cellranger-renaming",
            "--skip-cellrangermulti-vdjref",
            "--run-downstream",
            "--skip-downstream",
            "--save-reference",
            "--save-align-intermeds",
            "--no-save-align-intermeds",
            "--igenomes-ignore",
            "--star-ignore-sjdbgtf",
            "--monochrome-logs",
        },
        "accepts_genotypes": False,
    },
    "rnaseq-pipeline": {
        "script": SKILLS_DIR / "nfcore-rnaseq-wrapper" / "nfcore_rnaseq_wrapper.py",
        "demo_args": ["--demo"],
        "description": "Wrapper de preprocessing bulk RNA-seq FASTQ-to-counts vía nf-core/rnaseq",
        "default_timeout_seconds": 60 * 60 * 12 + 10 * 60,
        "max_output_files_listed": 50,
        "allowed_extra_flags": {
            "--check",
            "-v",
            "--verbose",
            "--no-banner",
            "--profile",
            "--pipeline-version",
            "--allow-pipeline-version-override",
            "--pipeline-local",
            "--timeout-hours",
            "--resume",
            "--aligner",
            "--pseudo-aligner",
            "--pseudo-aligner-kmer-size",
            "--prokaryotic",
            "--genome",
            "--fasta",
            "--gtf",
            "--gff",
            "--gene-bed",
            "--transcript-fasta",
            "--additional-fasta",
            "--splicesites",
            "--star-index",
            "--hisat2-index",
            "--rsem-index",
            "--salmon-index",
            "--kallisto-index",
            "--bowtie2-index",
            "--gencode",
            "--trimmer",
            "--extra-trimgalore-args",
            "--extra-fastp-args",
            "--min-trimmed-reads",
            "--remove-ribo-rna",
            "--ribo-removal-tool",
            "--with-umi",
            "--umi-dedup-tool",
            "--umitools-extract-method",
            "--umitools-bc-pattern",
            "--umitools-bc-pattern2",
            "--umitools-umi-separator",
            "--umi-discard-read",
            "--umitools-grouping-method",
            "--skip-umi-extract",
            "--seq-center",
            "--seq-platform",
            "--min-mapped-reads",
            "--star-ignore-sjdbgtf",
            "--salmon-quant-libtype",
            "--extra-star-align-args",
            "--extra-bowtie2-align-args",
            "--extra-salmon-quant-args",
            "--extra-kallisto-quant-args",
            "--kallisto-quant-fraglen",
            "--kallisto-quant-fraglen-sd",
            "--bam-csi-index",
            "--stringtie-ignore-gtf",
            "--stranded-threshold",
            "--unstranded-threshold",
            "--contaminant-screening",
            "--contaminant-screening-input",
            "--kraken-db",
            "--bracken-precision",
            "--sylph-db",
            "--sylph-taxonomy",
            "--bbsplit-fasta-list",
            "--bbsplit-index",
            "--skip-bbsplit",
            "--save-kraken-assignments",
            "--save-kraken-unassigned",
            "--featurecounts-group-type",
            "--featurecounts-feature-type",
            "--gtf-extra-attributes",
            "--gtf-group-features",
            "--deseq2-vst",
            "--no-deseq2-vst",
            "--rseqc-modules",
            "--skip-trimming",
            "--skip-alignment",
            "--skip-pseudo-alignment",
            "--skip-quantification-merge",
            "--skip-markduplicates",
            "--skip-bigwig",
            "--skip-stringtie",
            "--skip-fastqc",
            "--skip-dupradar",
            "--skip-qualimap",
            "--skip-rseqc",
            "--skip-biotype-qc",
            "--skip-deseq2-qc",
            "--skip-multiqc",
            "--enable-preseq",
            "--skip-qc",
            "--save-reference",
            "--save-trimmed",
            "--save-align-intermeds",
            "--save-unaligned",
            "--save-merged-fastq",
            "--save-non-ribo-reads",
            "--save-umi-intermeds",
            "--email",
            "--email-on-fail",
            "--multiqc-title",
            "--multiqc-config",
            "--multiqc-logo",
            "--multiqc-methods-description",
            "--rsem-extra-args",
            "--run-downstream",
            "--skip-downstream",
            "--metadata",
            "--formula",
            "--contrast",
            "--downstream-output",
            "--igenomes-base",
            "--sortmerna-index",
            "--ribo-database-manifest",
            "--hisat2-build-memory",
            "--gpu-container-options",
            "--extra-fqlint-args",
            "--publish-dir-mode",
            "--rapid-quant",
            "--arm",
            "--save-bbsplit-reads",
            "--gffread-transcript-fasta",
            "--use-rustqc",
            "--use-parabricks-star",
            "--use-sentieon-star",
            "--use-gpu-ribodetector",
            "--skip-linting",
            "--skip-gtf-filter",
            "--skip-gtf-transcript-filter",
            "--umitools-dedup-stats",
            "--umitools-dedup-primary-only",
            "--nextflow-config",
            "-c",
            "--config",
            "--work-dir",
            "--allow-remote-inputs",
        },
        "allowed_extra_flags_without_values": {
            "--check",
            "-v",
            "--verbose",
            "--no-banner",
            "--allow-pipeline-version-override",
            "--allow-remote-inputs",
            "--resume",
            "--prokaryotic",
            "--rapid-quant",
            "--arm",
            "--gencode",
            "--remove-ribo-rna",
            "--with-umi",
            "--skip-umi-extract",
            "--umitools-dedup-stats",
            "--umitools-dedup-primary-only",
            "--star-ignore-sjdbgtf",
            "--bam-csi-index",
            "--stringtie-ignore-gtf",
            "--gffread-transcript-fasta",
            "--deseq2-vst",
            "--no-deseq2-vst",
            "--skip-trimming",
            "--skip-alignment",
            "--skip-pseudo-alignment",
            "--skip-quantification-merge",
            "--skip-markduplicates",
            "--skip-bigwig",
            "--skip-stringtie",
            "--skip-fastqc",
            "--skip-dupradar",
            "--skip-qualimap",
            "--skip-rseqc",
            "--skip-biotype-qc",
            "--skip-deseq2-qc",
            "--skip-multiqc",
            "--enable-preseq",
            "--skip-qc",
            "--skip-linting",
            "--skip-gtf-filter",
            "--skip-gtf-transcript-filter",
            "--save-reference",
            "--save-trimmed",
            "--save-align-intermeds",
            "--save-unaligned",
            "--save-merged-fastq",
            "--save-non-ribo-reads",
            "--save-umi-intermeds",
            "--save-bbsplit-reads",
            "--skip-bbsplit",
            "--save-kraken-assignments",
            "--save-kraken-unassigned",
            "--run-downstream",
            "--skip-downstream",
            "--use-rustqc",
            "--use-parabricks-star",
            "--use-sentieon-star",
            "--use-gpu-ribodetector",
        },
        "accepts_genotypes": False,
    },
    "sarek-pipeline": {
        "script": SKILLS_DIR / "nfcore-sarek-wrapper" / "nfcore_sarek_wrapper.py",
        "demo_args": ["--demo"],
        "description": "Wrapper de variant calling germinal y somático vía nf-core/sarek 3.8.1",
        "default_timeout_seconds": 60 * 60 * 12 + 10 * 60,
        "max_output_files_listed": 50,
        # Keep this allowlist aligned with nfcore_sarek_wrapper.build_parser().
        "allowed_extra_flags": set("""
            --check --resume --arm --gpu --spark-profile --mutect-profile
            --run-downstream --downstream-skill --profile --nextflow-config -c --config
            --pipeline-version --allow-pipeline-version-override --pipeline-local --params-file --no-banner
            --verbose --timeout-hours --work-dir --allow-remote-inputs --extra-param --step --tools --skip-tools --aligner
            --no-intervals --wes --joint-germline --joint-mutect2
            --only-paired-variant-calling --ignore-soft-clipped-bases
            --filter-vcfs --normalize-vcfs --snv-consensus-calling
            --concatenate-vcfs --build-only-index --download-cache
            --use-gatk-spark --seq-platform --seq-center --email
            --email-on-fail --publish-dir-mode --outdir-cache --multiqc-title
            --multiqc-config --multiqc-logo --multiqc-methods-description
            --hook-url --trace-report-suffix --max-multiqc-email-size
            --input-restart --trim-fastq --trim-nextseq --clip-r1 --clip-r2
            --three-prime-clip-r1 --three-prime-clip-r2 --length-required
            --split-fastq --save-split-fastqs --save-trimmed
            --umi-read-structure --group-by-umi-strategy --umi-location
            --umi-tag --umi-length --umi-base-skip --umi-in-read-header
            --sentieon-consensus --save-mapped --save-output-as-bam
            --markduplicates-pixel-distance --nucleotides-per-second
            --ascat-min-base-qual --ascat-min-counts --ascat-min-map-qual
            --ascat-ploidy --ascat-purity --ascat-genome --cf-coeff
            --cf-contamination --cf-contamination-adjustment --cf-minqual
            --cf-mincov --cf-ploidy --cf-window --cf-chrom-len
            --cnvkit-reference --freebayes-filter
            --sentieon-haplotyper-emit-mode --sentieon-dnascope-emit-mode
            --sentieon-dnascope-pcr-indel-model --gatk-pcr-indel-model
            --varlociraptor-chunk-size --varlociraptor-scenario-tumor-only
            --varlociraptor-scenario-somatic --varlociraptor-scenario-germline
            --consensus-min-count --bcftools-filter-criteria --bcftools-columns
            --bcftools-header-lines --snpeff-db --snpeff-cache --vep-cache
            --vep-cache-version --vep-genome --vep-species --vep-version
            --vep-out-format --vep-custom-args --vep-include-fasta
            --vep-condel --vep-dbnsfp --vep-loftee --vep-mastermind
            --vep-phenotypes --vep-spliceai --vep-spliceregion --dbnsfp
            --dbnsfp-tbi --dbnsfp-consequence --dbnsfp-fields
            --mastermind-file --mastermind-mutations --mastermind-var-iden
            --mastermind-url --phenotypes-file --phenotypes-file-tbi
            --phenotypes-include-types --spliceai-snv --spliceai-snv-tbi
            --spliceai-indel --spliceai-indel-tbi --bcftools-annotations
            --bcftools-annotations-tbi --condel-config --snpsift-databases
            --genome --igenomes-base --igenomes-ignore --fasta --fasta-fai
            --dict --bwa --bwamem2 --dragmap --dbsnp --dbsnp-tbi
            --dbsnp-vqsr --known-indels --known-indels-tbi
            --known-indels-vqsr --known-snps --known-snps-tbi
            --known-snps-vqsr --germline-resource --germline-resource-tbi
            --pon --pon-tbi --intervals --ascat-alleles --ascat-loci
            --ascat-loci-gc --ascat-loci-rt --chr-dir --mappability
            --msisensor2-models --msisensorpro-scan --ngscheckmate-bed
            --sentieon-dnascope-model --bbsplit-fasta-list --bbsplit-index
            --save-reference --save-bbsplit-reads
        """.split()),
        "allowed_extra_flags_without_values": set("""
            --check --resume --arm --gpu --spark-profile --mutect-profile
            --allow-pipeline-version-override --allow-remote-inputs
            --run-downstream --no-banner --verbose --no-intervals --wes
            --joint-germline --joint-mutect2 --only-paired-variant-calling
            --ignore-soft-clipped-bases --filter-vcfs --normalize-vcfs
            --snv-consensus-calling --concatenate-vcfs --build-only-index
            --download-cache --trim-fastq --trim-nextseq --save-split-fastqs
            --save-trimmed --umi-in-read-header --sentieon-consensus
            --save-mapped --save-output-as-bam --cf-contamination-adjustment
            --vep-include-fasta --vep-condel --vep-dbnsfp --vep-loftee
            --vep-mastermind --vep-phenotypes --vep-spliceai
            --vep-spliceregion --mastermind-mutations --mastermind-var-iden
            --mastermind-url --igenomes-ignore --save-reference
            --save-bbsplit-reads
        """.split()),
        # The wrapper itself decides which native no-input modes are legal.
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "rdoutlier": {
        "script": SKILLS_DIR / "rare-disease-rnaseq" / "rare_disease_rnaseq.py",
        "demo_args": ["--demo"],
        "description": "Rare-disease blood RNA-seq outlier detection (NGRL-style: cases vs control panel + disease-gene filter)",
        "allowed_extra_flags": {
            "--counts",
            "--cases",
            "--controls",
            "--panel",
            "--z-threshold",
            "--output",
            "--seed",
        },
    },
    "methylation": {
        "script": SKILLS_DIR / "methylation-clock" / "methylation_clock.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "methylation-clock" / "data" / "GSE139307_small.csv.gz"),
        ],
        "description": "Epigenetic age from methylation clocks (PyAging)",
        "no_input_required": True,
        "allowed_extra_flags": {
            "--geo-id",
            "--clocks",
            "--metadata-cols",
            "--imputer-strategy",
            "--skip-epicv2-aggregation",
            "--verbose",
        },
    },
    "diffviz": {
        "script": SKILLS_DIR / "diff-visualizer" / "diff_visualizer.py",
        "demo_args": ["--demo"],
        "description": "Differential expression visualizer (bulk RNA-seq + scRNA downstream figure/report pack)",
        "allowed_extra_flags": {
            "--mode",
            "--counts",
            "--metadata",
            "--adata",
            "--top-genes",
            "--label-top",
            "--padj-threshold",
            "--lfc-threshold",
            "--min-basemean",
        },
        "accepts_genotypes": False,
    },
    "protocols-io": {
        "script": SKILLS_DIR / "protocols-io" / "protocols_io.py",
        "demo_args": ["--demo"],
        "description": "protocols.io bridge — search, browse, and retrieve scientific protocols via REST API",
        "allowed_extra_flags": {
            "--login",
            "--search",
            "--protocol",
            "--steps",
            "--dump",
            "--page-size",
            "--page",
            "--filter",
        },
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "acmg": {
        "script": SKILLS_DIR / "clinical-variant-reporter" / "clinical_variant_reporter.py",
        "demo_args": ["--demo"],
        "description": "ACMG/AMP clinical variant classifier (28-criteria, SF v3.2 screening)",
        "allowed_extra_flags": {"--genes", "--assembly"},
        "accepts_genotypes": False,
    },
    "llm-bench": {
        "script": SKILLS_DIR / "llm-biobank-bench" / "llm_biobank_bench.py",
        "demo_args": ["--demo"],
        "description": "Benchmark LLMs on UK Biobank knowledge retrieval (4 tasks, 6 models)",
        "allowed_extra_flags": {
            "--task",
            "--models",
            "--schema19",
            "--schema27",
        },
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "mr": {
        "script": SKILLS_DIR / "mendelian-randomisation" / "mendelian_randomisation.py",
        "demo_args": ["--demo"],
        "description": "Mendelian Randomisation — two-sample MR with IVW, Egger, weighted median/mode + full sensitivity",
        "allowed_extra_flags": {"--instruments"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "eqtl-region": {
        "script": SKILLS_DIR / "eqtl-catalogue-region-fetch" / "eqtl_catalogue_region_fetch.py",
        "demo_args": ["--demo"],
        "description": "eQTL Catalogue region fetch — tabix-on-FTP cis-QTL summary stats per genomic window",
        "allowed_extra_flags": {"--list-demos", "--no-cache"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "gwas-region": {
        "script": SKILLS_DIR / "gwas-catalog-region-fetch" / "gwas_catalog_region_fetch.py",
        "demo_args": ["--demo"],
        "description": "GWAS Catalog region fetch — tabix-on-FTP harmonised summary stats per genomic window",
        "allowed_extra_flags": {"--list-demos", "--no-cache"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "ld-region": {
        "script": SKILLS_DIR / "ld-1000g-region-compute" / "ld_1000g_region_compute.py",
        "demo_args": ["--demo"],
        "description": "1000G LD region compute — plink 1.9 r² between a lead and partners in a region for one super-population",
        "allowed_extra_flags": {"--list-demos", "--no-cache", "--super-pop", "--panel"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "ukb-ppp-region": {
        "script": SKILLS_DIR / "ukb-ppp-region-fetch" / "ukb_ppp_region_fetch.py",
        "demo_args": ["--demo"],
        "description": "UKB-PPP region fetch: per-variant plasma cis-pQTL summary stats per genomic window (Sun 2023, Synapse-backed)",
        "allowed_extra_flags": {"--list-demos", "--no-cache"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "locuscompare-region": {
        "script": SKILLS_DIR / "locuscompare-region-render" / "cli.py",
        "demo_args": ["--demo"],
        "description": "LocusCompare regional diagnostic — 4-panel coloc visualisation (Liu 2019) for one (eqtl × gwas) pair",
        "allowed_extra_flags": {"--list-demos", "--no-cache", "--super-pop"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "affprot": {
        "script": SKILLS_DIR / "affinity-proteomics" / "affinity_proteomics.py",
        "demo_args": ["--demo", "--platform", "olink"],
        "description": "Affinity proteomics — Olink NPX + SomaLogic SomaScan differential abundance",
        "allowed_extra_flags": {
            "--platform", "--meta", "--group-col", "--contrast",
            "--fdr", "--fc", "--top-n", "--test",
        },
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "gwas-pipe": {
        "script": SKILLS_DIR / "gwas-pipeline" / "gwas_pipeline.py",
        "demo_args": ["--demo"],
        "description": "GWAS pipeline — PLINK2 QC + REGENIE two-step association (Manhattan, QQ, lead variants)",
        "allowed_extra_flags": {
            "--bed", "--bgen", "--pheno", "--covar",
            "--trait-type", "--trait",
            "--geno", "--mind", "--maf", "--hwe",
        },
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "flow": {
        "script": SKILLS_DIR / "flow-bio" / "flow_bio.py",
        "demo_args": ["--demo"],
        "description": "Flow.bio API bridge (pipelines, samples, projects, executions)",
        "allowed_extra_flags": {
            "--login", "--username", "--password", "--token", "--url",
            "--pipelines", "--samples", "--projects", "--executions",
            "--organisms", "--sample-types", "--data",
            "--pipeline", "--sample", "--execution",
            "--metadata-attributes",
            "--pipeline-detail", "--sample-detail", "--execution-detail",
            "--search", "--search-samples", "--upload-sample", "--name", "--sample-type",
            "--reads1", "--reads2", "--organism", "--project",
            "--run-pipeline", "--run-samples", "--run-data", "--run-params",
            "--genome", "--json",
        },
        "no_input_required": True,
        "accepts_genotypes": False,
    },
    "sample-qc": {
        "script": SKILLS_DIR / "sample-qc-triage" / "sample_qc_triage.py",
        "demo_args": ["--demo"],
        "description": "Sample QC triage (identity, sex, contamination, batch-shift outlier triage)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": False,
    },
    "crispr-triage": {
        "script": SKILLS_DIR / "crispr-screen-triage" / "crispr_screen_triage.py",
        "demo_args": ["--demo"],
        "description": "CRISPR screen triage (deterministic guide-level hit ranking)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": False,
    },
    "marker-map": {
        "script": SKILLS_DIR / "marker-dominance-mapper" / "marker_dominance_mapper.py",
        "demo_args": ["--demo"],
        "description": "Marker dominance mapper (marker-based spot regions + SVG map)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": False,
    },
    "cnv-acmg": {
        "script": SKILLS_DIR / "cnv-acmg-classifier" / "cnv_acmg_classifier.py",
        "demo_args": ["--demo"],
        "description": "CNV/SV ACMG classifier (ClinGen 2019 / Riggs 2020 dosage point framework)",
        "allowed_extra_flags": {"--dosage-map", "--gene-model"},
        "accepts_genotypes": False,
    },
    "fastreer": {
        "script": SKILLS_DIR / "fastreer" / "fastreer.py",
        "demo_args": ["--demo"],
        "description": "fastreer: phylogenetic trees and distance matrices from VCF/FASTA",
        "allowed_extra_flags": {
            "--command", "--threads", "--mem", "--bootstrap",
            "--kmer", "--window-bp", "--window-variants", "--timeout", "--verbose",
        },
        "no_input_required": False,
        "accepts_genotypes": False,
    },
    "bioqc": {
        "script": SKILLS_DIR / "bioqc-mcp" / "bioqc_mcp.py",
        "demo_args": ["--demo"],
        "description": "BioQC quality control & custom visualizer (wraps FastQC + MultiQC, 20+ chart types, dual CLI/MCP server)",
        "allowed_extra_flags": {
            "--input", "--output", "--threads", "--mode", "--chart-type",
            "--chart-data", "--title", "--x-label", "--y-label", "--style",
            "--width", "--height",
        },
        "accepts_genotypes": False,
    },
    "pathway-enricher": {
        "script": SKILLS_DIR / "pathway-enricher" / "pathway_enricher.py",
        "demo_args": ["--demo"],
        "description": "Gene-set pathway enrichment via Enrichr (KEGG, GO BP/MF/CC, Reactome, WikiPathways) — bubble charts & ranked tables",
        "allowed_extra_flags": {
            "--input", "--output", "--databases", "--top-n",
        },
        "accepts_genotypes": False,
    },
    "repurposing": {
        "script": SKILLS_DIR / "drug-repurposing-screen" / "drug_repurposing_screen.py",
        "demo_args": ["--demo"],
        "description": "Drug repurposing screen — pooled viability QC, hits, selectivity, biomarkers, priority table",
        "allowed_extra_flags": {"--bundle", "--schema", "--objective", "--resume"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
}

try:
    from clawbio.skill_intents import DescriptorError, augment_skill_registry_with_descriptors

    SKILLS = augment_skill_registry_with_descriptors(SKILLS, CLAWBIO_DIR)
except DescriptorError as exc:
    # Descriptor routing is optional; keep the static registry usable if a
    # descriptor is malformed or violates descriptor security constraints.
    print(f"Warning: ignored invalid skill intent descriptor: {exc}", file=sys.stderr)

# Skills that run in the full-profile pipeline (order matters)
FULL_PROFILE_PIPELINE = ["pharmgx", "nutrigx", "prs", "compare"]

# --------------------------------------------------------------------------- #
# list_skills
# --------------------------------------------------------------------------- #



def _discover_skills() -> dict[str, dict]:
    """Auto-discover skill directories that aren't in the static SKILLS dict."""
    discovered = {}
    if not SKILLS_DIR.exists():
        return discovered
    
    for skill_dir in sorted(SKILLS_DIR.iterdir()):
        if not skill_dir.is_dir() or skill_dir.name.startswith("."):
            continue
        if skill_dir.name in SKILLS:
            continue
        
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.exists():
            continue
        
        has_script = any(
            f.suffix == ".py" and f.name != "__init__.py"
            for f in sorted(skill_dir.rglob("*.py"))
            if "tests" not in str(f.relative_to(skill_dir)).split("/")[0:1]
            and "__pycache__" not in str(f)
        )
        
        description = ""
        try:
            text = skill_md.read_text(encoding="utf-8")
            import re
            desc_match = re.search(r"^description:\s*(.+)", text, re.MULTILINE)
            if desc_match:
                description = desc_match.group(1).strip().strip('"').strip("'")
        except Exception:
            pass
        
        if not description:
            description = f"(auto-discovered from {skill_dir.name}/)"
        
        discovered[skill_dir.name] = {
            "has_script": has_script,
            "description": description,
        }
    
    return discovered


def list_skills() -> dict:
    """Print available skills and return the registry dict."""
    print(f"{BOLD}ClawBio Skills{RESET}")
    print(f"{'═' * 60}")
    
    print(f"{BOLD}Registered Skills{RESET} ({len(SKILLS)}):")
    for name, info in SKILLS.items():
        script_exists = info["script"].exists()
        status = f"{GREEN}OK{RESET}" if script_exists else f"{RED}MISSING{RESET}"
        print(f"  {BOLD}{name:<20}{RESET} {info['description']}")
        print(f"  {'':20} {DIM}script: {info['script'].name}{RESET} [{status}]")
        print()
    
    discovered = _discover_skills()
    if discovered:
        print(f"{BOLD}Auto-discovered Skills{RESET} ({len(discovered)}):")
        for name, info in sorted(discovered.items()):
            indicator = f"{GREEN}has script{RESET}" if info["has_script"] else f"{YELLOW}spec only{RESET}"
            print(f"  {DIM}{name:<20}{RESET} {info['description']}")
            print(f"  {'':20} {DIM}[{indicator}]{RESET}")
            print()
    
    total = len(SKILLS) + len(discovered)
    print(f"{DIM}Total: {total} skills ({len(SKILLS)} registered, {len(discovered)} auto-discovered){RESET}")
    print()
    print(f"{DIM}Run a skill:  clawbio run <skill> --demo{RESET}")
    print(f"{DIM}With input:   clawbio run <skill> --input <file>{RESET}")
    print(f"{DIM}Upload once:  clawbio upload --input <file> --patient-id PT001{RESET}")
    print(f"{DIM}Full profile: clawbio run full-profile --profile profiles/PT001.json{RESET}")
    return SKILLS


def upload_profile(
    input_path: str,
    patient_id: str = "",
    fmt: str = "auto",
) -> dict:
    """Parse a genetic file and save a PatientProfile.

    Returns a dict with profile path and metadata.
    """
    # Lazy import to avoid requiring clawbio package for basic subprocess usage
    if str(CLAWBIO_DIR) not in sys.path:
        sys.path.insert(0, str(CLAWBIO_DIR))
    from clawbio.common.profile import PatientProfile

    profile = PatientProfile.from_genetic_file(input_path, patient_id=patient_id, fmt=fmt)
    PROFILES_DIR.mkdir(parents=True, exist_ok=True)
    pid = profile.metadata["patient_id"]
    profile_path = PROFILES_DIR / f"{pid}.json"
    profile.save(profile_path)

    return {
        "success": True,
        "profile_path": str(profile_path),
        "patient_id": pid,
        "genotype_count": profile.genotype_count,
        "checksum": profile.metadata["checksum"],
    }


# --------------------------------------------------------------------------- #
# run_skill
# --------------------------------------------------------------------------- #


def _load_structured_skill_result(out_dir: Path | None) -> tuple[dict | None, Path | None]:
    """Load a skill's result.json envelope when present and valid."""
    if out_dir is None:
        return None, None
    result_json_path = out_dir / "result.json"
    if not result_json_path.exists():
        return None, None
    try:
        payload = json.loads(result_json_path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError):
        return None, result_json_path
    if not isinstance(payload, dict):
        return None, result_json_path
    return payload, result_json_path


def _load_report_markdown(out_dir: Path | None) -> str | None:
    """Read the primary markdown report from an output bundle when present."""
    if out_dir is None or not out_dir.exists():
        return None
    for pattern in ("report.md", "*_report.md", "*.md"):
        for md_file in sorted(out_dir.glob(pattern)):
            if md_file.name.startswith("."):
                continue
            try:
                return md_file.read_text(encoding="utf-8")
            except OSError:
                continue
    return None


def _promote_structured_result_fields(result: dict, out_dir: Path | None) -> None:
    """Attach parsed result.json fields to the top-level run result."""
    payload, result_json_path = _load_structured_skill_result(out_dir)
    if payload is not None:
        result["skill_result_json"] = payload
    if result_json_path is not None:
        result["result_json_path"] = str(result_json_path)

    if isinstance(payload, dict):
        # Structured result fields form the small skill-to-ClawBio display and
        # action contract:
        # - chat_summary_lines: concise, skill-authored text for chat UIs
        # - preferred_artifacts: generated files the UI should surface first
        # - suggested_actions: deterministic next-step requests to offer later
        # - workflow_state: skill-emitted state identity/lifecycle metadata
        # - contract_alerts: structured contract/path discrepancy alerts
        # - report_md: full markdown report text embedded in result.json
        for field in (
            "chat_summary_lines",
            "preferred_artifacts",
            "suggested_actions",
            "workflow_state",
            "contract_alerts",
            "report_md",
        ):
            if field in payload:
                result[field] = (
                    normalise_contract_alerts(payload[field])
                    if field == "contract_alerts"
                    else payload[field]
                )

        if out_dir is not None and result.get("contract_alerts"):
            append_contract_alert_log(
                out_dir / "contract_alerts.jsonl",
                result["contract_alerts"],
                run_id=out_dir.name,
                skill=result.get("skill") or None,
            )

    if "report_md" not in result:
        report_md = _load_report_markdown(out_dir)
        if report_md is not None:
            result["report_md"] = report_md


_NFCORE_PIPELINE_SKILLS = {"scrnaseq-pipeline", "rnaseq-pipeline", "sarek-pipeline"}


def _canonical_extra_flag(flag: str) -> str:
    """Normalise a long CLI flag to its hyphenated spelling.

    nf-core parameters are snake_case (``--skip_tools``) while the pipeline
    wrappers expose them as hyphenated flags (``--skip-tools``). Mapping the two
    to a single canonical form lets the INT-001 allowlist recognise an
    nf-core-native command copied from the upstream docs without enumerating
    every underscore variant. Only the long-flag name is touched; short flags and
    values are returned unchanged.
    """
    if flag.startswith("--"):
        return "--" + flag[2:].replace("_", "-")
    return flag


def run_skill(
    skill_name: str,
    input_path: str | None = None,
    output_dir: str | None = None,
    demo: bool = False,
    extra_args: list[str] | None = None,
    timeout: int = 300,
    profile_path: str | None = None,
) -> dict:
    """
    Run a ClawBio skill as a subprocess.

    Returns a structured dict with success status, output paths, and logs.
    Importable by any agent (RoboTerri, RoboIsaac, Claude Code).
    """
    # Handle full-profile virtual skill
    if skill_name == "full-profile":
        return _run_full_profile(
            profile_path=profile_path,
            input_path=input_path,
            output_dir=output_dir,
            timeout=timeout,
        )

    # Validate skill
    skill_info = SKILLS.get(skill_name)
    if not skill_info:
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": None,
            "files": [],
            "stdout": "",
            "stderr": f"Unknown skill '{skill_name}'. Available: {list(SKILLS.keys())}",
            "duration_seconds": 0,
        }

    script_path = skill_info["script"]
    if not script_path.exists():
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": None,
            "files": [],
            "stdout": "",
            "stderr": f"Script not found: {script_path}",
            "duration_seconds": 0,
        }

    # If --profile is given, resolve the input file from the profile
    resolved_input = input_path
    if resolved_input:
        input_candidate = Path(resolved_input).expanduser()
        if not input_candidate.is_absolute():
            input_candidate = Path.cwd() / input_candidate
        resolved_input = str(input_candidate.resolve())
    if profile_path and not input_path and not demo:
        if str(CLAWBIO_DIR) not in sys.path:
            sys.path.insert(0, str(CLAWBIO_DIR))
        from clawbio.common.profile import PatientProfile
        profile = PatientProfile.load(profile_path)
        stored_input = profile.metadata.get("input_file", "")
        if stored_input:
            # Resolve relative paths against CLAWBIO_DIR
            p = Path(stored_input)
            if not p.is_absolute():
                p = CLAWBIO_DIR / p
            if p.exists():
                resolved_input = str(p.resolve())

    # Build output directory
    summary_mode = skill_info.get("summary_default", False) and not output_dir
    if summary_mode:
        out_dir = None
    elif output_dir:
        out_dir = Path(output_dir).expanduser().resolve()
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = DEFAULT_OUTPUT_ROOT / f"{skill_name}_{ts}"
    if out_dir:
        output_error = _ensure_output_directory(out_dir)
        if output_error:
            return {
                "skill": skill_name,
                "success": False,
                "exit_code": -1,
                "output_dir": str(out_dir),
                "files": [],
                "stdout": "",
                "stderr": json.dumps(output_error, indent=2),
                "duration_seconds": 0,
            }

    # Build command
    cmd = [PYTHON, str(script_path)]

    if demo:
        cmd.extend(skill_info["demo_args"])
    elif resolved_input:
        cmd.extend(["--input", str(resolved_input)])
    elif not skill_info.get("no_input_required"):
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": str(out_dir) if out_dir else None,
            "files": [],
            "stdout": "",
            "stderr": "No input provided. Use --demo, --input <file>, or --profile <path>.",
            "duration_seconds": 0,
        }

    if out_dir:
        cmd.extend(["--output", str(out_dir)])

    # SEC INT-001: filter extra_args against per-skill allowlist
    if extra_args:
        allowed = skill_info.get("allowed_extra_flags", set())
        flags_without_values = skill_info.get("allowed_extra_flags_without_values", set())
        blocked = {"--input", "--output", "--demo"}
        # nf-core parameters are snake_case; the pipeline wrappers expose them as
        # hyphenated flags. For those skills, treat the two spellings as
        # equivalent so an nf-core-native command is not silently dropped, and
        # forward the wrapper's canonical hyphen spelling (which its parser
        # already accepts). Other skills keep exact-match behaviour unchanged.
        canonicalize = skill_name in _NFCORE_PIPELINE_SKILLS

        def _key(token: str) -> str:
            name = token.split("=", 1)[0]
            return _canonical_extra_flag(name) if canonicalize else name

        filtered = []
        i = 0
        while i < len(extra_args):
            token = extra_args[i]
            flag = _key(token)
            if flag in blocked:
                i += 2 if "=" not in token and i + 1 < len(extra_args) else i + 1
                continue
            if flag in allowed:
                _name, sep, value = token.partition("=")
                filtered.append(f"{flag}={value}" if sep else flag)
                if (
                    not sep
                    and flag not in flags_without_values
                    and i + 1 < len(extra_args)
                    and _key(extra_args[i + 1]) not in allowed
                    and _key(extra_args[i + 1]) not in blocked
                ):
                    filtered.append(extra_args[i + 1])
                    i += 1
            i += 1
        cmd.extend(filtered)

    # Run subprocess
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(script_path.parent),
        )
        duration = round(time.time() - t0, 2)
    except subprocess.TimeoutExpired:
        duration = round(time.time() - t0, 2)
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": str(out_dir) if out_dir else None,
            "files": [],
            "stdout": "",
            "stderr": f"Timed out after {timeout} seconds.",
            "duration_seconds": duration,
        }
    except Exception as e:
        duration = round(time.time() - t0, 2)
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": str(out_dir) if out_dir else None,
            "files": [],
            "stdout": "",
            "stderr": str(e),
            "duration_seconds": duration,
        }

    # Collect output files
    if out_dir and out_dir.exists():
        max_files = int(skill_info.get("max_output_files_listed", 200))
        all_output_files = sorted(f.name for f in out_dir.rglob("*") if f.is_file())
        output_files = all_output_files[:max_files]
        if len(all_output_files) > max_files:
            output_files.append(f"... {len(all_output_files) - max_files} more files")
    else:
        output_files = []

    result = {
        "skill": skill_name,
        "success": proc.returncode == 0,
        "exit_code": proc.returncode,
        "output_dir": str(out_dir) if out_dir else None,
        "files": output_files,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "duration_seconds": duration,
    }

    if result["success"]:
        _promote_structured_result_fields(result, out_dir)

    # If profile was used, store the result back into it
    if profile_path and result["success"] and out_dir:
        _store_result_in_profile(profile_path, skill_name, out_dir)

    return result


def _ensure_output_directory(out_dir: Path) -> dict[str, object] | None:
    if out_dir.exists() and not out_dir.is_dir():
        return {
            "ok": False,
            "stage": "preflight",
            "error_code": "OUTPUT_DIR_NOT_WRITABLE",
            "message": "Output path exists but is not a directory.",
            "fix": "Choose a directory path for --output, or remove/rename the existing file.",
            "details": {"output": str(out_dir)},
        }
    try:
        out_dir.mkdir(parents=True, exist_ok=True)
    except OSError as exc:
        return {
            "ok": False,
            "stage": "preflight",
            "error_code": "OUTPUT_DIR_NOT_WRITABLE",
            "message": "Output directory could not be created.",
            "fix": "Choose a writable output location.",
            "details": {"output": str(out_dir), "error": str(exc)},
        }
    return None


# --------------------------------------------------------------------------- #
# Full-profile pipeline
# --------------------------------------------------------------------------- #


def _run_full_profile(
    profile_path: str | None,
    input_path: str | None,
    output_dir: str | None,
    timeout: int = 300,
) -> dict:
    """Run all genotype-consuming skills sequentially, accumulating results."""
    if not profile_path and not input_path:
        return {
            "skill": "full-profile",
            "success": False,
            "exit_code": -1,
            "output_dir": None,
            "files": [],
            "stdout": "",
            "stderr": "full-profile requires --profile or --input.",
            "duration_seconds": 0,
        }

    # Create profile if only input was given
    if not profile_path and input_path:
        upload_result = upload_profile(input_path)
        if not upload_result["success"]:
            return {
                "skill": "full-profile",
                "success": False,
                "exit_code": -1,
                "output_dir": None,
                "files": [],
                "stdout": "",
                "stderr": "Failed to create profile from input file.",
                "duration_seconds": 0,
            }
        profile_path = upload_result["profile_path"]

    # Setup output
    if output_dir:
        out_dir = Path(output_dir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = DEFAULT_OUTPUT_ROOT / f"full_profile_{ts}"
    out_dir.mkdir(parents=True, exist_ok=True)

    t0 = time.time()
    all_results = {}
    all_files = []
    combined_stdout = []
    combined_stderr = []
    any_failure = False

    for skill_name in FULL_PROFILE_PIPELINE:
        skill_out = out_dir / skill_name
        print(f"  Running {skill_name}...")
        result = run_skill(
            skill_name=skill_name,
            profile_path=profile_path,
            output_dir=str(skill_out),
            timeout=timeout,
        )
        all_results[skill_name] = {
            "success": result["success"],
            "exit_code": result["exit_code"],
            "files": result["files"],
        }
        if result["stdout"]:
            combined_stdout.append(f"=== {skill_name} ===\n{result['stdout']}")
        if result["stderr"]:
            combined_stderr.append(f"=== {skill_name} ===\n{result['stderr']}")
        all_files.extend(result["files"])
        if not result["success"]:
            any_failure = True
            print(f"    WARNING: {skill_name} failed (exit {result['exit_code']})")

    duration = round(time.time() - t0, 2)

    # Write aggregate summary
    summary = {
        "pipeline": FULL_PROFILE_PIPELINE,
        "profile": profile_path,
        "results": all_results,
        "completed_at": datetime.now().isoformat(),
    }
    summary_path = out_dir / "pipeline_summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, default=str))

    return {
        "skill": "full-profile",
        "success": not any_failure,
        "exit_code": 0 if not any_failure else 1,
        "output_dir": str(out_dir),
        "files": all_files + ["pipeline_summary.json"],
        "stdout": "\n\n".join(combined_stdout),
        "stderr": "\n\n".join(combined_stderr),
        "duration_seconds": duration,
    }


def _store_result_in_profile(profile_path: str, skill_name: str, out_dir: Path) -> None:
    """Load result.json from a skill's output and store it in the profile."""
    try:
        if str(CLAWBIO_DIR) not in sys.path:
            sys.path.insert(0, str(CLAWBIO_DIR))
        from clawbio.common.profile import PatientProfile

        result_json = out_dir / "result.json"
        if not result_json.exists():
            return

        profile = PatientProfile.load(profile_path)
        result_data = json.loads(result_json.read_text())
        profile.add_skill_result(skill_name, result_data)
        profile.save(profile_path)
    except Exception:
        pass  # Don't fail the main pipeline for profile storage issues


# --------------------------------------------------------------------------- #
# CLI entry point
# --------------------------------------------------------------------------- #


def main():
    # Pipeline wrappers own large, schema-derived CLIs. Delegate their help so
    # `clawbio.py run <pipeline> --help` cannot drift from the wrapper parser.
    #
    # This applies to EVERY nf-core pipeline, not just sarek. The launcher declares only a
    # subset of each wrapper's flags as argparse options; the rest are forwarded through the
    # per-skill allowlist. Rendering the launcher's own help therefore hides flags that
    # genuinely work — `--allow-remote-inputs` was the case in point: supported by all three
    # wrappers, but discoverable in --help only for sarek, because the delegation below was
    # hardcoded to sarek-pipeline.
    if (
        len(sys.argv) >= 4
        and sys.argv[1] == "run"
        and sys.argv[2] in _NFCORE_PIPELINE_SKILLS
        and any(arg in {"-h", "--help"} for arg in sys.argv[3:])
    ):
        subprocess.run([PYTHON, str(SKILLS[sys.argv[2]]["script"]), "--help"], check=False)
        return

    parser = argparse.ArgumentParser(
        description="ClawBio — Bioinformatics Skills Runner",
    )
    parser.add_argument(
        "--version", action="version", version=f"%(prog)s {__version__}"
    )
    sub = parser.add_subparsers(dest="command")

    # list
    sub.add_parser("list", help="List available skills")

    # upload
    upload_parser = sub.add_parser("upload", help="Upload genetic data and create a patient profile")
    upload_parser.add_argument("--input", required=True, dest="input_path", help="Path to genetic data file")
    upload_parser.add_argument("--patient-id", default="", help="Patient identifier (default: derived from filename)")
    upload_parser.add_argument("--format", default="auto", help="File format: auto, 23andme, ancestry, vcf")

    # run
    run_parser = sub.add_parser("run", help="Run a skill")
    run_parser.add_argument("skill", help="Skill name (e.g. pharmgx, equity, full-profile)")
    run_parser.add_argument("--demo", action="store_true", help="Run with demo data")
    run_parser.add_argument("--input", dest="input_path", help="Path to input file")
    run_parser.add_argument("--output", dest="output_dir", help="Output directory")
    run_parser.add_argument("--profile", dest="profile_path", metavar="PROFILE", help="Execution backend profile or comma-separated nf-core profile list for pipeline skills, or path to patient profile JSON for genomics skills")
    run_parser.add_argument(
        "--timeout", type=int, default=300, help="Timeout in seconds (default: 300)"
    )
    run_parser.add_argument("--drug", default=None, help="Drug name for single-drug lookup (drugphoto skill)")
    run_parser.add_argument("--dose", default=None, help="Visible dose from packaging (e.g. '50mg')")
    run_parser.add_argument("--trait", default=None, help="Trait search term for PRS skill")
    run_parser.add_argument("--pgs-id", default=None, help="PGS Catalog score ID for PRS skill")
    run_parser.add_argument("--gene", default=None, help="Gene symbol for ClinPGx skill")
    run_parser.add_argument("--genes", default=None, help="Comma-separated gene symbols for ClinPGx")
    run_parser.add_argument("--rsid", default=None, help="rsID for GWAS lookup skill (e.g. rs3798220)")
    run_parser.add_argument("--skip", default=None, help="Comma-separated API names to skip (gwas-lookup skill)")
    run_parser.add_argument("--query", default=None, help="Inline SQL query for bigquery skill")
    run_parser.add_argument("--location", default=None, help="BigQuery location (e.g. US, EU)")
    run_parser.add_argument("--max-rows", type=int, default=None, help="Maximum number of query rows for bigquery skill")
    run_parser.add_argument(
        "--max-bytes-billed",
        type=int,
        default=None,
        help="Maximum billed bytes safeguard for bigquery skill",
    )
    run_parser.add_argument(
        "--param",
        action="append",
        default=None,
        help="Repeatable bigquery parameter in name=type:value format",
    )
    run_parser.add_argument("--dry-run", action="store_true", help="BigQuery dry-run (estimate bytes only)")
    run_parser.add_argument("--list-datasets", default=None, help="List BigQuery datasets for a project")
    run_parser.add_argument("--list-tables", default=None, help="List BigQuery tables for a dataset (project.dataset)")
    run_parser.add_argument("--describe", default=None, help="Describe a BigQuery table schema (project.dataset.table)")
    run_parser.add_argument("--preview", type=int, default=None, help="Preview wrapper row limit for bigquery skill")
    run_parser.add_argument("--count-only", action="store_true", help="Return only row count for bigquery skill")
    run_parser.add_argument("--paper", default=None, help="Paper reference/DOI/URL/path for bigquery provenance")
    run_parser.add_argument("--note", action="append", default=None, help="Repeatable provenance note for bigquery skill")
    run_parser.add_argument("--geo-id", default=None, help="GEO accession for methylation clock skill")
    run_parser.add_argument("--clocks", default=None, help="Comma-separated clock names for methylation skill")
    run_parser.add_argument("--metadata-cols", default=None, help="Comma-separated metadata columns for methylation skill")
    run_parser.add_argument("--imputer-strategy", default=None, help="Imputer strategy for methylation skill")
    run_parser.add_argument("--skip-epicv2-aggregation", action="store_true", help="Skip EPICv2 probe aggregation")
    run_parser.add_argument("--verbose", action="store_true", help="Verbose output for skill backends")
    run_parser.add_argument("--vcf", default=None, help="Explicit VCF override for illumina skill")
    run_parser.add_argument("--qc", default=None, help="Explicit QC metrics override for illumina skill")
    run_parser.add_argument("--sample-sheet", default=None, help="Explicit SampleSheet override for illumina skill")
    run_parser.add_argument(
        "--metadata-provider",
        default=None,
        help="Optional metadata provider for illumina skill (none or ica)",
    )
    run_parser.add_argument("--ica-project-id", default=None, help="ICA project ID for illumina skill")
    run_parser.add_argument("--ica-run-id", default=None, help="ICA analysis/run ID for illumina skill")
    run_parser.add_argument("--counts", default=None, help="Counts matrix for rnaseq/diffviz bulk workflows")
    run_parser.add_argument("--metadata", default=None, help="Sample metadata for rnaseq/diffviz bulk workflows")
    run_parser.add_argument("--formula", default=None, help="Design formula for rnaseq skill")
    run_parser.add_argument("--contrast", default=None, help="Contrast for rnaseq skill: factor,numerator,denominator")
    run_parser.add_argument("--backend", default=None, help="Backend for rnaseq skill (auto|pydeseq2|simple)")
    run_parser.add_argument("--min-count", type=int, default=None, help="Minimum count threshold for rnaseq skill")
    run_parser.add_argument("--min-samples", type=int, default=None, help="Minimum samples threshold for rnaseq skill")
    run_parser.add_argument("--check", action="store_true", help="Preflight-only mode for pipeline skills (scrnaseq-pipeline, rnaseq-pipeline)")
    # Nextflow config passthrough for the nf-core pipeline skills. The three
    # spellings (-c / --config / --nextflow-config) are aliases of one repeatable
    # option — matching every wrapper's own parser — so a config supplied with any
    # spelling reaches scrnaseq-pipeline, rnaseq-pipeline, and sarek-pipeline alike.
    run_parser.add_argument(
        "-c",
        "--config",
        "--nextflow-config",
        dest="nextflow_config",
        action="append",
        metavar="CONFIG",
        default=None,
        help="Additional Nextflow config file(s) forwarded to nf-core pipeline skills as -c "
        "(aliases: -c/--config/--nextflow-config; repeatable).",
    )
    run_parser.add_argument("--pipeline-version", default=None, help="Pinned pipeline version/tag for pipeline skills")
    run_parser.add_argument("--allow-dirty-pipeline", action="store_true", help="Allow an intentionally modified local scrnaseq checkout")
    run_parser.add_argument("--require-local-pipeline", action="store_true", help="Require a verifiable local scrnaseq checkout")
    run_parser.add_argument("--preset", default=None, help="Curated preset for pipeline skills")
    run_parser.add_argument("--protocol", default=None, help="Protocol value for pipeline skills")
    run_parser.add_argument("--email", default=None, help="Email address for pipeline completion notification")
    run_parser.add_argument("--email-on-fail", default=None, help="Email address for failure notification only")
    run_parser.add_argument("--multiqc-title", default=None, help="Custom MultiQC title for pipeline skills")
    run_parser.add_argument("--expected-cells", type=int, default=None, help="expected_cells override for scrnaseq-pipeline")
    run_parser.add_argument("--work-dir", default=None, help="Nextflow work directory for scrnaseq-pipeline")
    run_parser.add_argument("--allow-conda-cellranger", action="store_true", help="Allow Cell Ranger presets with conda/mamba when a trusted site config supplies Cell Ranger")
    run_parser.add_argument("--resume", action="store_true", help="Enable resume policy for pipeline skills")
    run_parser.add_argument("--save-reference", action="store_true", help="Save built reference indexes for pipeline skills")
    run_parser.add_argument(
        "--save-align-intermeds",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Save alignment intermediates for pipeline skills; use --no-save-align-intermeds to force off",
    )
    run_parser.add_argument("--skip-cellbender", action="store_true", help="Disable cellbender for scrnaseq-pipeline")
    run_parser.add_argument("--skip-fastqc", action="store_true", help="Skip FastQC for pipeline skills")
    run_parser.add_argument(
        "--skip-emptydrops",
        action="store_true",
        help="Deprecated alias for --skip-cellbender in scrnaseq-pipeline",
    )
    run_parser.add_argument("--skip-multiqc", action="store_true", help="Skip MultiQC for pipeline skills")
    run_parser.add_argument("--skip-cellranger-renaming", action="store_true", help="Skip CellRanger sample renaming")
    run_parser.add_argument("--skip-cellrangermulti-vdjref", action="store_true", help="Skip CellRanger Multi VDJ reference build")
    run_parser.add_argument("--run-downstream", action="store_true", help="Opt in to scrna_orchestrator handoff after scrnaseq-pipeline")
    run_parser.add_argument("--skip-downstream", action="store_true", help="Compatibility flag for pipeline downstream handoff")
    run_parser.add_argument("--fasta", default=None, help="Genome FASTA for pipeline skills")
    run_parser.add_argument("--gtf", default=None, help="Annotation GTF for scrnaseq-pipeline")
    run_parser.add_argument("--transcript-fasta", default=None, help="Transcript FASTA for scrnaseq-pipeline")
    run_parser.add_argument("--txp2gene", default=None, help="Transcript-to-gene map for scrnaseq-pipeline")
    run_parser.add_argument("--simpleaf-index", default=None, help="Prebuilt simpleaf index for scrnaseq-pipeline")
    run_parser.add_argument("--simpleaf-umi-resolution", default=None, help="simpleaf UMI resolution strategy")
    run_parser.add_argument("--kallisto-index", default=None, help="Prebuilt kallisto index for scrnaseq-pipeline")
    run_parser.add_argument("--kb-workflow", default=None, help="Kallisto workflow for scrnaseq-pipeline")
    run_parser.add_argument("--kb-t1c", default=None, help="Kallisto cDNA transcripts-to-capture file")
    run_parser.add_argument("--kb-t2c", default=None, help="Kallisto intron transcripts-to-capture file")
    run_parser.add_argument("--star-index", default=None, help="Prebuilt STAR index for scrnaseq-pipeline")
    run_parser.add_argument("--star-feature", default=None, help="STARsolo feature type for scrnaseq-pipeline")
    run_parser.add_argument("--star-ignore-sjdbgtf", action="store_true", help="Disable STAR SJDB GTF usage")
    run_parser.add_argument("--seq-center", default=None, help="Sequencing center for scrnaseq-pipeline")
    run_parser.add_argument("--cellranger-index", default=None, help="Prebuilt cellranger index for scrnaseq-pipeline")
    run_parser.add_argument("--cellranger-vdj-index", default=None, help="Prebuilt CellRanger VDJ reference index")
    run_parser.add_argument("--cellrangerarc-config", default=None, help="CellRanger ARC config file")
    run_parser.add_argument("--cellrangerarc-reference", default=None, help="CellRanger ARC reference name")
    run_parser.add_argument("--barcode-whitelist", default=None, help="Barcode whitelist override for scrnaseq-pipeline")
    run_parser.add_argument("--motifs", default=None, help="Motif file for CellRanger ARC")
    run_parser.add_argument("--gex-frna-probe-set", default=None, help="CellRanger Multi fixed RNA probe set")
    run_parser.add_argument("--gex-target-panel", default=None, help="CellRanger Multi target panel")
    run_parser.add_argument("--gex-cmo-set", default=None, help="CellRanger Multi CMO set")
    run_parser.add_argument("--fb-reference", default=None, help="Feature barcoding reference CSV")
    run_parser.add_argument("--vdj-inner-enrichment-primers", default=None, help="VDJ inner enrichment primers file")
    run_parser.add_argument("--gex-barcode-sample-assignment", default=None, help="GEX barcode sample assignment CSV")
    run_parser.add_argument("--cellranger-multi-barcodes", default=None, help="CellRanger Multi barcodes samplesheet")
    run_parser.add_argument("--downstream-output", default=None, help="Output directory for rnaseq-pipeline automatic rnaseq handoff")

    # rnaseq-pipeline (nfcore-rnaseq-wrapper) specific flags. Declared here so they
    # appear in `python clawbio.py run rnaseq-pipeline --help` and parse with the
    # correct types instead of surviving only via parse_known_args + the registry
    # allowlist.
    run_parser.add_argument("--pipeline-local", default=None, help="Local nf-core pipeline checkout for rnaseq-pipeline or sarek-pipeline")
    run_parser.add_argument("--aligner", default=None, help="Aligner for pipeline skills: scrnaseq-pipeline uses nf-core/scrnaseq aligners; rnaseq-pipeline uses star_salmon|star_rsem|hisat2|bowtie2_salmon; sarek-pipeline uses its aligner override")
    run_parser.add_argument("--pseudo-aligner", default=None, help="Pseudo-aligner for rnaseq-pipeline (salmon|kallisto)")
    run_parser.add_argument("--pseudo-aligner-kmer-size", type=int, default=None, help="K-mer size for pseudo-aligner index")
    run_parser.add_argument("--prokaryotic", action="store_true", help="Compose -profile prokaryotic for rnaseq-pipeline")
    run_parser.add_argument("--gff", default=None, help="GFF3 annotation for rnaseq-pipeline")
    run_parser.add_argument("--gene-bed", default=None, help="Gene BED file for rnaseq-pipeline")
    run_parser.add_argument("--additional-fasta", default=None, help="Additional FASTA (spike-ins) for rnaseq-pipeline")
    run_parser.add_argument("--splicesites", default=None, help="HISAT2 splice sites file")
    run_parser.add_argument("--hisat2-index", default=None, help="Prebuilt HISAT2 index for rnaseq-pipeline")
    run_parser.add_argument("--rsem-index", default=None, help="Prebuilt RSEM index for rnaseq-pipeline")
    run_parser.add_argument("--bowtie2-index", default=None, help="Prebuilt Bowtie2 index for rnaseq-pipeline")
    run_parser.add_argument("--salmon-index", default=None, help="Prebuilt Salmon index for rnaseq-pipeline")
    run_parser.add_argument("--gencode", action="store_true", help="GENCODE GTF flag for rnaseq-pipeline")
    run_parser.add_argument("--trimmer", default=None, help="Trimmer for rnaseq-pipeline (trimgalore|fastp)")
    run_parser.add_argument("--extra-trimgalore-args", default=None, help="Extra args passed to TrimGalore")
    run_parser.add_argument("--extra-fastp-args", default=None, help="Extra args passed to fastp")
    run_parser.add_argument("--min-trimmed-reads", type=int, default=None, help="Minimum trimmed reads per sample")
    run_parser.add_argument("--remove-ribo-rna", action="store_true", help="Enable rRNA removal for rnaseq-pipeline")
    run_parser.add_argument("--ribo-removal-tool", default=None, help="rRNA removal tool (sortmerna|ribodetector|bowtie2)")
    run_parser.add_argument("--with-umi", action="store_true", help="Enable UMI handling for rnaseq-pipeline")
    run_parser.add_argument("--umi-dedup-tool", default=None, help="UMI dedup tool (umitools|umicollapse)")
    run_parser.add_argument("--umitools-extract-method", default=None, choices=["string", "regex"], help="umi_tools extract method (nf-core schema: string|regex)")
    run_parser.add_argument("--umitools-bc-pattern", default=None, help="umi_tools barcode pattern")
    run_parser.add_argument("--umitools-bc-pattern2", default=None, help="umi_tools barcode pattern (read 2)")
    run_parser.add_argument("--umi-discard-read", type=int, default=None, help="UMI: discard read 1 or 2")
    run_parser.add_argument("--umitools-umi-separator", default=None, help="umi_tools UMI separator")
    run_parser.add_argument("--umitools-grouping-method", default=None, choices=["unique", "percentile", "cluster", "adjacency", "directional"], help="umi_tools grouping method (nf-core schema enum)")
    run_parser.add_argument("--skip-umi-extract", action="store_true", help="Skip UMI extraction (UMIs already in headers)")
    run_parser.add_argument("--salmon-quant-libtype", default=None, choices=["A", "IS", "ISF", "ISR", "IU", "MS", "MSF", "MSR", "MU", "OS", "OSF", "OSR", "OU", "SF", "SR", "U"], help="Force Salmon library type (nf-core schema enum)")
    run_parser.add_argument("--extra-star-align-args", default=None, help="Extra args passed to STAR (star_salmon only)")
    run_parser.add_argument("--extra-bowtie2-align-args", default=None, help="Extra args passed to Bowtie2 (bowtie2_salmon only)")
    run_parser.add_argument("--extra-salmon-quant-args", default=None, help="Extra args passed to salmon quant")
    run_parser.add_argument("--extra-kallisto-quant-args", default=None, help="Extra args passed to kallisto quant")
    run_parser.add_argument("--rsem-extra-args", default=None, help="Extra args forwarded to rsem-calculate-expression (star_rsem only)")
    run_parser.add_argument("--kallisto-quant-fraglen", type=int, default=None, help="Kallisto fragment length (single-end)")
    run_parser.add_argument("--kallisto-quant-fraglen-sd", type=int, default=None, help="Kallisto fragment length SD (single-end)")
    run_parser.add_argument("--min-mapped-reads", type=float, default=None, help="Minimum %% unique mapped reads (default 5)")
    run_parser.add_argument("--stranded-threshold", type=float, default=None, help="Salmon strand inference threshold (default 0.8)")
    run_parser.add_argument("--unstranded-threshold", type=float, default=None, help="Salmon strand inference unstranded threshold (default 0.1)")
    run_parser.add_argument("--contaminant-screening", default=None, choices=["kraken2", "kraken2_bracken", "sylph"], help="Contaminant screening engine for rnaseq-pipeline")
    run_parser.add_argument("--contaminant-screening-input", default=None, choices=["trimmed", "unmapped"], help="Reads used for contaminant screening")
    run_parser.add_argument("--kraken-db", default=None, help="Kraken2 database for contaminant screening")
    run_parser.add_argument("--bracken-precision", default=None, choices=["D", "P", "C", "O", "F", "G", "S"], help="Bracken precision level (D P C O F G S)")
    run_parser.add_argument("--sylph-db", default=None, help="Sylph database for contaminant screening")
    run_parser.add_argument("--sylph-taxonomy", default=None, help="Sylph taxonomy TSV")
    run_parser.add_argument("--bbsplit-fasta-list", default=None, help="BBSplit FASTA list")
    run_parser.add_argument("--bbsplit-index", default=None, help="Prebuilt BBSplit index")
    run_parser.add_argument("--skip-bbsplit", action="store_true", help="Skip BBSplit even when --bbsplit-fasta-list or --bbsplit-index is provided")
    run_parser.add_argument("--save-kraken-assignments", action="store_true", help="Save Kraken assignment outputs")
    run_parser.add_argument("--save-kraken-unassigned", action="store_true", help="Save Kraken unassigned reads")
    run_parser.add_argument("--bam-csi-index", action="store_true", help="Use CSI BAM index instead of BAI")
    run_parser.add_argument("--stringtie-ignore-gtf", action="store_true", help="StringTie de novo assembly without GTF")
    run_parser.add_argument("--gtf-extra-attributes", default=None, help="GTF extra attributes for Salmon")
    run_parser.add_argument("--gtf-group-features", default=None, help="GTF grouping attribute for Salmon")
    run_parser.add_argument("--featurecounts-group-type", default=None, help="featureCounts biotype group attribute")
    run_parser.add_argument("--featurecounts-feature-type", default=None, help="featureCounts feature type")
    run_parser.add_argument("--deseq2-vst", action="store_true", default=None, help="Force DESeq2 VST on")
    run_parser.add_argument("--no-deseq2-vst", dest="deseq2_vst", action="store_false", help="Disable DESeq2 VST")
    run_parser.add_argument("--rseqc-modules", default=None, help="Comma-separated RSeQC modules to run")
    run_parser.add_argument("--skip-trimming", action="store_true", help="Skip adapter trimming")
    run_parser.add_argument("--skip-alignment", action="store_true", help="Skip all alignment-based processes")
    run_parser.add_argument("--skip-pseudo-alignment", action="store_true", help="Skip pseudo-alignment processes")
    run_parser.add_argument("--skip-quantification-merge", action="store_true", help="Per-sample TSVs instead of merged matrix")
    run_parser.add_argument("--skip-markduplicates", action="store_true", help="Skip Picard MarkDuplicates")
    run_parser.add_argument("--skip-bigwig", action="store_true", help="Skip bigWig creation")
    run_parser.add_argument("--skip-stringtie", action="store_true", help="Skip StringTie")
    run_parser.add_argument("--skip-dupradar", action="store_true", help="Skip dupRadar")
    run_parser.add_argument("--skip-qualimap", action="store_true", help="Skip Qualimap")
    run_parser.add_argument("--skip-rseqc", action="store_true", help="Skip RSeQC")
    run_parser.add_argument("--skip-biotype-qc", action="store_true", help="Skip featureCounts biotype QC")
    run_parser.add_argument("--skip-deseq2-qc", action="store_true", help="Skip DESeq2 PCA/heatmap")
    run_parser.add_argument("--skip-qc", action="store_true", help="Skip all QC except MultiQC")
    run_parser.add_argument("--enable-preseq", action="store_true", help="Enable Preseq complexity estimation")
    run_parser.add_argument("--save-trimmed", action="store_true", help="Save trimmed FASTQs")
    run_parser.add_argument("--save-unaligned", action="store_true", help="Save unaligned reads")
    run_parser.add_argument("--save-merged-fastq", action="store_true", help="Save merged technical-replicate FASTQs")
    run_parser.add_argument("--save-non-ribo-reads", action="store_true", help="Save non-rRNA reads")
    run_parser.add_argument("--save-umi-intermeds", action="store_true", help="Save UMI intermediate files")
    run_parser.add_argument("--seq-platform", default=None, help="Sequencing platform for BAM read group")
    run_parser.add_argument("--multiqc-config", default=None, help="MultiQC config file for rnaseq-pipeline")
    run_parser.add_argument("--multiqc-logo", default=None, help="MultiQC logo image for rnaseq-pipeline")
    run_parser.add_argument("--multiqc-methods-description", default=None, help="MultiQC methods description YAML for rnaseq-pipeline")
    run_parser.add_argument("--trace-report-suffix", default=None, help="Suffix appended to Nextflow trace report filenames for pipeline skills")
    run_parser.add_argument("--monochrome-logs", action="store_true", help="Disable ANSI colors in Nextflow logs for pipeline skills")
    # rnaseq-pipeline — profile modifiers and architecture / GPU acceleration
    run_parser.add_argument("--arm", action="store_true", help="Compose -profile arm64 for ARM64-native containers (Apple Silicon, AWS Graviton)")
    run_parser.add_argument("--rapid-quant", action="store_true", help="Compose -profile rapid_quant for rnaseq-pipeline (per-sample pseudo-alignment only)")
    run_parser.add_argument("--gpu-container-options", default=None, help="GPU container runtime flags for rnaseq-pipeline (e.g. '--gpus all' for Docker, '--nv' for Singularity)")
    run_parser.add_argument("--use-rustqc", action="store_true", help="Enable experimental RustQC consolidated QC for rnaseq-pipeline")
    run_parser.add_argument("--use-parabricks-star", action="store_true", help="Accelerate STAR + MarkDuplicates with Parabricks (GPU) for rnaseq-pipeline")
    run_parser.add_argument("--use-sentieon-star", action="store_true", help="Accelerate STAR with Sentieon for rnaseq-pipeline (requires SENTIEON_LICENSE_BASE64 secret)")
    run_parser.add_argument("--use-gpu-ribodetector", action="store_true", help="Use GPU acceleration for Ribodetector rRNA removal")
    # rnaseq-pipeline — references and index overrides not already covered above
    run_parser.add_argument("--igenomes-base", default=None, help="iGenomes base path override (S3 URI or local mirror) for pipeline skills")
    run_parser.add_argument("--sortmerna-index", default=None, help="Prebuilt SortMeRNA index for rnaseq-pipeline")
    run_parser.add_argument("--ribo-database-manifest", default=None, help="SortMeRNA rRNA database manifest for rnaseq-pipeline")
    run_parser.add_argument("--hisat2-build-memory", default=None, help="HISAT2 index build memory (e.g. '200.GB') for rnaseq-pipeline")
    run_parser.add_argument("--gffread-transcript-fasta", action="store_true", help="Use gffread instead of RSEM to derive transcript FASTA")
    # rnaseq-pipeline — UMI extras
    run_parser.add_argument("--umitools-dedup-stats", action="store_true", help="Generate UMI-tools dedup statistics for rnaseq-pipeline")
    run_parser.add_argument("--umitools-dedup-primary-only", action="store_true", help="Restrict UMI-tools dedup to primary alignments")
    # rnaseq-pipeline — linting / GTF filtering / extras
    run_parser.add_argument("--extra-fqlint-args", default=None, help="Extra fq-lint arguments for rnaseq-pipeline")
    run_parser.add_argument("--skip-linting", action="store_true", help="Skip fq-lint pre-trimming FASTQ checks")
    run_parser.add_argument("--skip-gtf-filter", action="store_true", help="Skip GTF scaffold/transcript ID filtering")
    run_parser.add_argument("--skip-gtf-transcript-filter", action="store_true", help="Skip the transcript_id check inside the GTF filter")
    # rnaseq-pipeline — output and publishing controls
    run_parser.add_argument("--save-bbsplit-reads", action="store_true", help="Save BBSplit-separated FASTQ files for rnaseq-pipeline")
    run_parser.add_argument("--publish-dir-mode", default=None, choices=["symlink", "rellink", "link", "copy", "copyNoFollow", "move"], help="Nextflow publishDir mode for pipeline skills (pipeline default: copy)")
    # Nextflow config passthrough (-c/--config/--nextflow-config) is declared once,
    # above, as a single repeatable option shared by every nf-core pipeline skill.
    run_parser.add_argument("--mode", default=None, help="Mode for diffviz skill (auto|bulk|scrna)")
    run_parser.add_argument("--adata", default=None, help="AnnData input for enhanced diffviz scRNA plots")
    run_parser.add_argument("--top-genes", type=int, default=None, help="Top genes/markers to display in diffviz")
    run_parser.add_argument("--label-top", type=int, default=None, help="Label top hits in diffviz plots")
    run_parser.add_argument(
        "--padj-threshold",
        type=float,
        default=None,
        help="Adjusted p-value threshold for diffviz significance highlighting",
    )
    run_parser.add_argument(
        "--lfc-threshold",
        type=float,
        default=None,
        help="Absolute log fold-change threshold for diffviz significance highlighting",
    )
    run_parser.add_argument(
        "--min-basemean",
        type=float,
        default=None,
        help="Minimum baseMean retained in diffviz bulk display plots/tables",
    )
    run_parser.add_argument("--method", default=None, help="Embedding backend (scrna-embedding skill)")
    run_parser.add_argument("--layer", default=None, help="Raw-count layer for `.h5ad` input (scrna-embedding skill)")
    run_parser.add_argument("--batch-key", default=None, help="obs batch column for integration (scrna-embedding skill)")
    run_parser.add_argument("--labels-key", default=None, help="obs label column for scANVI (scrna-embedding skill)")
    run_parser.add_argument(
        "--unlabeled-category",
        default=None,
        help="Category value representing unlabeled cells for scANVI (scrna-embedding skill)",
    )
    run_parser.add_argument("--min-genes", type=int, default=None, help="Minimum genes per cell (scrna/scrna-embedding skill)")
    run_parser.add_argument("--min-cells", type=int, default=None, help="Minimum cells per gene (scrna/scrna-embedding skill)")
    run_parser.add_argument(
        "--max-mt-pct",
        type=float,
        default=None,
        help="Maximum mitochondrial percentage (scrna/scrna-embedding skill)",
    )
    run_parser.add_argument(
        "--n-top-hvg",
        type=int,
        default=None,
        help="Number of highly variable genes to keep (scrna/scrna-embedding skill)",
    )
    run_parser.add_argument("--n-pcs", type=int, default=None, help="Number of PCA components (scrna skill)")
    run_parser.add_argument("--latent-dim", type=int, default=None, help="Latent dimensionality (scrna-embedding skill)")
    run_parser.add_argument("--max-epochs", type=int, default=None, help="Max training epochs (scrna-embedding skill)")
    run_parser.add_argument(
        "--n-neighbors",
        type=int,
        default=None,
        help="Neighbors for graph construction (scrna/scrna-embedding skill)",
    )
    run_parser.add_argument(
        "--use-rep",
        default=None,
        help="Graph representation key or mode such as `auto`, `none`, or `X_scvi` (scrna skill)",
    )
    run_parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=None,
        help="Leiden resolution (scrna skill)",
    )
    run_parser.add_argument("--random-state", type=int, default=None, help="Random seed (scrna/scrna-embedding skill)")
    run_parser.add_argument(
        "--top-markers",
        type=int,
        default=None,
        help="Top markers per cluster (scrna skill)",
    )
    run_parser.add_argument(
        "--accelerator",
        default=None,
        help="Training accelerator (scrna-embedding skill)",
    )
    run_parser.add_argument(
        "--contrast-groupby",
        default=None,
        help="obs column for contrastive marker analysis (scrna skill)",
    )
    run_parser.add_argument(
        "--contrast-scope",
        default=None,
        help="Contrast scope: dataset, within-cluster, or both (scrna skill)",
    )
    run_parser.add_argument(
        "--contrast-clusterby",
        default=None,
        help="Cluster/partition column for within-cluster contrasts (scrna skill)",
    )
    run_parser.add_argument(
        "--contrast-top-genes",
        type=int,
        default=None,
        help="Top contrastive marker genes in summary table (scrna skill)",
    )
    run_parser.add_argument(
        "--doublet-method",
        default=None,
        help="Optional doublet detection method for scrna skill",
    )
    run_parser.add_argument(
        "--annotate",
        default=None,
        help="Optional annotation backend for scrna skill",
    )
    run_parser.add_argument(
        "--annotation-model",
        default=None,
        help="Local CellTypist model name or path for scrna skill",
    )
    run_parser.add_argument("--search", default=None, help="Search query (bioc / galaxy skills)")
    run_parser.add_argument("--recommend", default=None, help="Recommendation query for bioc skill")
    run_parser.add_argument("--workflow", default=None, help="Workflow query for bioc skill")
    run_parser.add_argument("--package-details", default=None, help="Bioconductor package name for bioc skill")
    run_parser.add_argument("--docs-search", default=None, help="Documentation search query for bioc skill")
    run_parser.add_argument("--package-docs", default=None, help="Fetch package documentation for bioc skill")
    run_parser.add_argument("--list-domains", action="store_true", help="List supported Bioconductor domains")
    run_parser.add_argument("--setup", action="store_true", help="Inspect local Bioconductor setup")
    run_parser.add_argument("--install", default=None, help="Comma-separated Bioconductor packages to install")
    run_parser.add_argument("--format", dest="skill_format", default=None, help="Input format hint for bioc skill")
    run_parser.add_argument("--container", default=None, help="Canonical object/container hint for bioc skill")
    run_parser.add_argument("--modality", default=None, help="Modality hint for bioc skill")
    run_parser.add_argument("--max-results", type=int, default=None, help="Maximum bioc search/recommendation results")
    # flow-bio skill flags
    run_parser.add_argument("--flow-search", dest="flow_search", default=None, help="Search query (flow skill)")
    run_parser.add_argument("--pipelines", action="store_true", help="List pipelines (flow skill)")
    run_parser.add_argument("--samples", action="store_true", help="List samples (flow skill)")
    run_parser.add_argument("--projects", action="store_true", help="List projects (flow skill)")
    run_parser.add_argument("--executions", action="store_true", help="List executions (flow skill)")
    run_parser.add_argument("--organisms", action="store_true", help="List organisms (flow skill)")
    run_parser.add_argument("--sample-types", action="store_true", help="List sample types (flow skill)")
    run_parser.add_argument("--data", action="store_true", help="List data (flow skill)")
    run_parser.add_argument("--metadata-attributes", action="store_true", help="List metadata attributes (flow skill)")
    run_parser.add_argument("--search-samples", nargs="+", default=None, help="Search samples by metadata key=value pairs (flow skill)")
    run_parser.add_argument("--upload-sample", action="store_true", help="Upload a sample (flow skill)")
    run_parser.add_argument("--name", default=None, help="Sample name for upload (flow skill)")
    run_parser.add_argument("--reads1", default=None, help="First reads file (flow skill)")
    run_parser.add_argument("--reads2", default=None, help="Second reads file (flow skill)")
    run_parser.add_argument("--organism", default=None, help="Organism name or ID (flow skill)")
    run_parser.add_argument("--project", default=None, help="Project ID (flow skill)")
    run_parser.add_argument("--run-pipeline", default=None, help="Pipeline version ID to run (flow skill)")
    run_parser.add_argument("--run-samples", default=None, help="Comma-separated sample IDs for pipeline (flow skill)")
    run_parser.add_argument("--run-data", default=None, help="Comma-separated data IDs for pipeline (flow skill)")
    run_parser.add_argument("--run-params", default=None, help="Pipeline parameters as JSON string (flow skill)")
    run_parser.add_argument("--genome", default=None, help="Genome ID for pipeline run (flow/scrnaseq skill)")
    run_parser.add_argument("--pipeline-detail", default=None, dest="pipeline_detail", help="Get pipeline details by ID (flow skill)")
    run_parser.add_argument("--sample-detail", default=None, dest="sample_detail", help="Get sample details by ID (flow skill)")
    run_parser.add_argument("--execution-detail", default=None, dest="execution_detail", help="Get execution details by ID (flow skill)")
    run_parser.add_argument("--json", action="store_true", help="Output raw JSON (flow skill)")

    args, extra = parser.parse_known_args()

    if args.command == "list":
        list_skills()

    elif args.command == "upload":
        result = upload_profile(
            input_path=args.input_path,
            patient_id=args.patient_id,
            fmt=args.format,
        )
        if result["success"]:
            print(f"  Profile created: {result['profile_path']}")
            print(f"  Patient ID:      {result['patient_id']}")
            print(f"  Genotypes:       {result['genotype_count']}")
            print(f"  Checksum:        {result['checksum'][:16]}")
        else:
            print("  Upload failed.")
            sys.exit(1)

    elif args.command == "run":
        skill_backend_profile = None
        if args.skill == "scrnaseq-pipeline" and _is_scrnaseq_backend_profile(getattr(args, "profile_path", None)):
            skill_backend_profile = args.profile_path
            args.profile_path = None
        elif args.skill in {"rnaseq-pipeline", "sarek-pipeline"} and getattr(args, "profile_path", None) is not None:
            skill_backend_profile = args.profile_path
            args.profile_path = None

        # Build extra_args from skill-specific flags
        extra = list(extra)
        if getattr(args, "check", False):
            extra.append("--check")
        # Nextflow config files (-c/--config/--nextflow-config) apply to every
        # nf-core pipeline wrapper, which all accept those spellings as aliases of
        # one repeatable option. Forward each entry as --nextflow-config so a config
        # supplied with any spelling reaches the wrapper regardless of pipeline.
        if args.skill in {"scrnaseq-pipeline", "rnaseq-pipeline", "sarek-pipeline"}:
            for cfg in getattr(args, "nextflow_config", None) or []:
                extra.extend(["--nextflow-config", cfg])
        if skill_backend_profile:
            extra.extend(["--profile", skill_backend_profile])
        if getattr(args, "pipeline_version", None):
            extra.extend(["--pipeline-version", args.pipeline_version])
        if getattr(args, "preset", None):
            extra.extend(["--preset", args.preset])
        if getattr(args, "protocol", None):
            extra.extend(["--protocol", args.protocol])
        if getattr(args, "email", None):
            extra.extend(["--email", args.email])
        if getattr(args, "multiqc_title", None):
            extra.extend(["--multiqc-title", args.multiqc_title])
        if args.skill == "scrnaseq-pipeline":
            if getattr(args, "allow_dirty_pipeline", False):
                extra.append("--allow-dirty-pipeline")
            if getattr(args, "require_local_pipeline", False):
                extra.append("--require-local-pipeline")
            if getattr(args, "allow_conda_cellranger", False):
                extra.append("--allow-conda-cellranger")
            if getattr(args, "aligner", None):
                extra.extend(["--aligner", args.aligner])
            for value_flag in (
                ("email_on_fail", "--email-on-fail"),
                ("multiqc_config", "--multiqc-config"),
                ("multiqc_logo", "--multiqc-logo"),
                ("multiqc_methods_description", "--multiqc-methods-description"),
                ("igenomes_base", "--igenomes-base"),
                ("publish_dir_mode", "--publish-dir-mode"),
                ("trace_report_suffix", "--trace-report-suffix"),
                ("work_dir", "--work-dir"),
            ):
                attr, flag = value_flag
                value = getattr(args, attr, None)
                if value:
                    extra.extend([flag, value])
            if getattr(args, "monochrome_logs", False):
                extra.append("--monochrome-logs")
        if getattr(args, "expected_cells", None) is not None:
            extra.extend(["--expected-cells", str(args.expected_cells)])
        if getattr(args, "resume", False):
            extra.append("--resume")
        if getattr(args, "save_reference", False):
            extra.append("--save-reference")
        if getattr(args, "save_align_intermeds", None) is True:
            extra.append("--save-align-intermeds")
        elif getattr(args, "save_align_intermeds", None) is False:
            extra.append("--no-save-align-intermeds")
        if getattr(args, "skip_cellbender", False):
            extra.append("--skip-cellbender")
        if getattr(args, "skip_fastqc", False):
            extra.append("--skip-fastqc")
        if getattr(args, "skip_emptydrops", False):
            extra.append("--skip-emptydrops")
        if getattr(args, "skip_multiqc", False):
            extra.append("--skip-multiqc")
        if getattr(args, "skip_cellranger_renaming", False):
            extra.append("--skip-cellranger-renaming")
        if getattr(args, "skip_cellrangermulti_vdjref", False):
            extra.append("--skip-cellrangermulti-vdjref")
        if getattr(args, "run_downstream", False):
            extra.append("--run-downstream")
        if getattr(args, "skip_downstream", False):
            extra.append("--skip-downstream")
        if getattr(args, "fasta", None):
            extra.extend(["--fasta", args.fasta])
        if getattr(args, "gtf", None):
            extra.extend(["--gtf", args.gtf])
        if getattr(args, "transcript_fasta", None):
            extra.extend(["--transcript-fasta", args.transcript_fasta])
        if getattr(args, "txp2gene", None):
            extra.extend(["--txp2gene", args.txp2gene])
        if getattr(args, "simpleaf_index", None):
            extra.extend(["--simpleaf-index", args.simpleaf_index])
        if getattr(args, "simpleaf_umi_resolution", None):
            extra.extend(["--simpleaf-umi-resolution", args.simpleaf_umi_resolution])
        if getattr(args, "kallisto_index", None):
            extra.extend(["--kallisto-index", args.kallisto_index])
        if getattr(args, "kb_workflow", None):
            extra.extend(["--kb-workflow", args.kb_workflow])
        if getattr(args, "kb_t1c", None):
            extra.extend(["--kb-t1c", args.kb_t1c])
        if getattr(args, "kb_t2c", None):
            extra.extend(["--kb-t2c", args.kb_t2c])
        if getattr(args, "star_index", None):
            extra.extend(["--star-index", args.star_index])
        if getattr(args, "star_feature", None):
            extra.extend(["--star-feature", args.star_feature])
        if getattr(args, "star_ignore_sjdbgtf", False):
            extra.append("--star-ignore-sjdbgtf")
        if getattr(args, "seq_center", None):
            extra.extend(["--seq-center", args.seq_center])
        if getattr(args, "cellranger_index", None):
            extra.extend(["--cellranger-index", args.cellranger_index])
        if getattr(args, "cellranger_vdj_index", None):
            extra.extend(["--cellranger-vdj-index", args.cellranger_vdj_index])
        if getattr(args, "cellrangerarc_config", None):
            extra.extend(["--cellrangerarc-config", args.cellrangerarc_config])
        if getattr(args, "cellrangerarc_reference", None):
            extra.extend(["--cellrangerarc-reference", args.cellrangerarc_reference])
        if getattr(args, "barcode_whitelist", None):
            extra.extend(["--barcode-whitelist", args.barcode_whitelist])
        if getattr(args, "motifs", None):
            extra.extend(["--motifs", args.motifs])
        if getattr(args, "gex_frna_probe_set", None):
            extra.extend(["--gex-frna-probe-set", args.gex_frna_probe_set])
        if getattr(args, "gex_target_panel", None):
            extra.extend(["--gex-target-panel", args.gex_target_panel])
        if getattr(args, "gex_cmo_set", None):
            extra.extend(["--gex-cmo-set", args.gex_cmo_set])
        if getattr(args, "fb_reference", None):
            extra.extend(["--fb-reference", args.fb_reference])
        if getattr(args, "vdj_inner_enrichment_primers", None):
            extra.extend(["--vdj-inner-enrichment-primers", args.vdj_inner_enrichment_primers])
        if getattr(args, "gex_barcode_sample_assignment", None):
            extra.extend(["--gex-barcode-sample-assignment", args.gex_barcode_sample_assignment])
        if getattr(args, "cellranger_multi_barcodes", None):
            extra.extend(["--cellranger-multi-barcodes", args.cellranger_multi_barcodes])

        # rnaseq-pipeline (nfcore-rnaseq-wrapper) flag forwarding. Only forward when
        # the active skill is rnaseq-pipeline so we don't pollute other skills'
        # extra args with unrelated values.
        if args.skill == "rnaseq-pipeline":
            if getattr(args, "pipeline_local", None):
                extra.extend(["--pipeline-local", args.pipeline_local])
            for value_flag in (
                ("aligner", "--aligner"),
                ("pseudo_aligner", "--pseudo-aligner"),
                ("gff", "--gff"),
                ("gene_bed", "--gene-bed"),
                ("additional_fasta", "--additional-fasta"),
                ("splicesites", "--splicesites"),
                ("hisat2_index", "--hisat2-index"),
                ("rsem_index", "--rsem-index"),
                ("bowtie2_index", "--bowtie2-index"),
                ("salmon_index", "--salmon-index"),
                ("trimmer", "--trimmer"),
                ("extra_trimgalore_args", "--extra-trimgalore-args"),
                ("extra_fastp_args", "--extra-fastp-args"),
                ("ribo_removal_tool", "--ribo-removal-tool"),
                ("umi_dedup_tool", "--umi-dedup-tool"),
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
