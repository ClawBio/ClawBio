#!/usr/bin/env python3
"""
clawbio.py — ClawBio Bioinformatics Skills Runner
==================================================
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

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

CLAWBIO_DIR = Path(__file__).resolve().parent
SKILLS_DIR = CLAWBIO_DIR / "skills"
EXAMPLES_DIR = CLAWBIO_DIR / "examples"
DEFAULT_OUTPUT_ROOT = CLAWBIO_DIR / "output"
PROFILES_DIR = CLAWBIO_DIR / "profiles"

# Python binary — use the same interpreter that launched clawbio.py
PYTHON = sys.executable

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
        "script": SKILLS_DIR / "nutrigx_advisor" / "nutrigx_advisor.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "nutrigx_advisor" / "tests" / "synthetic_patient.csv"),
        ],
        "description": "Nutrigenomics advisor (diet, vitamins, caffeine, lactose)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": True,
    },
    "metagenomics": {
        "script": SKILLS_DIR / "claw-metagenomics" / "metagenomics_profiler.py",
        "demo_args": ["--demo"],
        "description": "Metagenomics profiler (Kraken2, RGI/CARD, HUMAnN3)",
        "allowed_extra_flags": set(),
        "accepts_genotypes": False,
    },
    "scrna": {
        "script": SKILLS_DIR / "scrna-orchestrator" / "scrna_orchestrator.py",
        "demo_args": ["--demo"],
        "description": "scRNA Orchestrator (Scanpy QC, clustering, marker detection, optional two-group DE + volcano)",
        "allowed_extra_flags": {
            "--min-genes",
            "--min-cells",
            "--max-mt-pct",
            "--n-top-hvg",
            "--n-pcs",
            "--n-neighbors",
            "--leiden-resolution",
            "--random-state",
            "--top-markers",
            "--de-groupby",
            "--de-group1",
            "--de-group2",
            "--de-top-genes",
            "--de-volcano",
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
    "profile": {
        "script": SKILLS_DIR / "profile-report" / "profile_report.py",
        "demo_args": ["--demo"],
        "description": "Unified personal genomic profile report",
        "allowed_extra_flags": {"--profile"},
        "no_input_required": True,
        "accepts_genotypes": False,
    },
}

# Skills that run in the full-profile pipeline (order matters)
FULL_PROFILE_PIPELINE = ["pharmgx", "nutrigx", "prs", "compare"]

# --------------------------------------------------------------------------- #
# list_skills
# --------------------------------------------------------------------------- #


def list_skills() -> dict:
    """Print available skills and return the registry dict."""
    print(f"{BOLD}ClawBio Skills{RESET}")
    print(f"{'═' * 55}")
    for name, info in SKILLS.items():
        script_exists = info["script"].exists()
        status = f"{GREEN}OK{RESET}" if script_exists else f"{RED}MISSING{RESET}"
        print(f"  {BOLD}{name:<15}{RESET} {info['description']}")
        print(f"  {'':15} {DIM}script: {info['script'].name}{RESET} [{status}]")
        print()
    print(f"{DIM}Run a skill:  python clawbio.py run <skill> --demo{RESET}")
    print(f"{DIM}With input:   python clawbio.py run <skill> --input <file>{RESET}")
    print(f"{DIM}Upload once:  python clawbio.py upload --input <file> --patient-id PT001{RESET}")
    print(f"{DIM}Full profile: python clawbio.py run full-profile --profile profiles/PT001.json{RESET}")
    return SKILLS


# --------------------------------------------------------------------------- #
# upload_profile
# --------------------------------------------------------------------------- #


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
        out_dir = Path(output_dir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = DEFAULT_OUTPUT_ROOT / f"{skill_name}_{ts}"
    if out_dir:
        out_dir.mkdir(parents=True, exist_ok=True)

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
        blocked = {"--input", "--output", "--demo"}
        filtered = []
        i = 0
        while i < len(extra_args):
            flag = extra_args[i].split("=")[0]
            if flag in blocked:
                i += 2 if "=" not in extra_args[i] and i + 1 < len(extra_args) else i + 1
                continue
            if flag in allowed:
                filtered.append(extra_args[i])
                if "=" not in extra_args[i] and i + 1 < len(extra_args) and not extra_args[i + 1].startswith("-"):
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
        output_files = sorted(
            [f.name for f in out_dir.rglob("*") if f.is_file()],
        )
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

    # If profile was used, store the result back into it
    if profile_path and result["success"] and out_dir:
        _store_result_in_profile(profile_path, skill_name, out_dir)

    return result


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
    parser = argparse.ArgumentParser(
        description="ClawBio — Bioinformatics Skills Runner",
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
    run_parser.add_argument("--profile", dest="profile_path", help="Path to patient profile JSON")
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
    run_parser.add_argument("--min-genes", type=int, default=None, help="Minimum genes per cell (scrna skill)")
    run_parser.add_argument("--min-cells", type=int, default=None, help="Minimum cells per gene (scrna skill)")
    run_parser.add_argument(
        "--max-mt-pct",
        type=float,
        default=None,
        help="Maximum mitochondrial percentage (scrna skill)",
    )
    run_parser.add_argument(
        "--n-top-hvg",
        type=int,
        default=None,
        help="Number of highly variable genes to keep (scrna skill)",
    )
    run_parser.add_argument("--n-pcs", type=int, default=None, help="Number of PCA components (scrna skill)")
    run_parser.add_argument(
        "--n-neighbors",
        type=int,
        default=None,
        help="Neighbors for graph construction (scrna skill)",
    )
    run_parser.add_argument(
        "--leiden-resolution",
        type=float,
        default=None,
        help="Leiden resolution (scrna skill)",
    )
    run_parser.add_argument("--random-state", type=int, default=None, help="Random seed (scrna skill)")
    run_parser.add_argument(
        "--top-markers",
        type=int,
        default=None,
        help="Top markers per cluster (scrna skill)",
    )
    run_parser.add_argument("--de-groupby", default=None, help="obs column for DE (scrna skill)")
    run_parser.add_argument("--de-group1", default=None, help="Group 1 value for DE (scrna skill)")
    run_parser.add_argument("--de-group2", default=None, help="Group 2 reference value for DE (scrna skill)")
    run_parser.add_argument(
        "--de-top-genes",
        type=int,
        default=None,
        help="Top DE genes in summary table (scrna skill)",
    )
    run_parser.add_argument(
        "--de-volcano",
        action="store_true",
        help="Generate DE volcano plot (scrna skill)",
    )

    args = parser.parse_args()

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
        # Build extra_args from skill-specific flags
        extra = []
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
        if getattr(args, "n_neighbors", None) is not None:
            extra.extend(["--n-neighbors", str(args.n_neighbors)])
        if getattr(args, "leiden_resolution", None) is not None:
            extra.extend(["--leiden-resolution", str(args.leiden_resolution)])
        if getattr(args, "random_state", None) is not None:
            extra.extend(["--random-state", str(args.random_state)])
        if getattr(args, "top_markers", None) is not None:
            extra.extend(["--top-markers", str(args.top_markers)])
        if getattr(args, "de_groupby", None):
            extra.extend(["--de-groupby", args.de_groupby])
        if getattr(args, "de_group1", None):
            extra.extend(["--de-group1", args.de_group1])
        if getattr(args, "de_group2", None):
            extra.extend(["--de-group2", args.de_group2])
        if getattr(args, "de_top_genes", None) is not None:
            extra.extend(["--de-top-genes", str(args.de_top_genes)])
        if getattr(args, "de_volcano", False):
            extra.append("--de-volcano")

        result = run_skill(
            skill_name=args.skill,
            input_path=args.input_path,
            output_dir=args.output_dir,
            demo=args.demo,
            extra_args=extra or None,
            timeout=args.timeout,
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
