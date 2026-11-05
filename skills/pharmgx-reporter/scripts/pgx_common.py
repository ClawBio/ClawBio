#!/usr/bin/env python3
"""Shared PGx data tables and pure logic for the pharmgx_reporter skill.

Ported from ClawBio 1.x `skills/pharmgx-reporter/pharmgx_reporter.py`
(https://github.com/ClawBio/ClawBio/tree/main/skills/pharmgx-reporter).
The star-allele calling, phenotype mapping, CPIC drug guideline tables, and
build-detection logic are unchanged from 1.x — only the module boundary is
new. 1.x shipped this as one ~2300-line script; this port splits it across
DAG steps (`parse_and_profile.py` -> `lookup_drugs.py` -> `enrich_clinpgx.py`
-> `report.py`, plus the standalone `dosage_card.py`), and those step
scripts all import this module rather than duplicating the ~450 lines of
gene/drug tables four times. This is a within-skill convenience import (all
files live under this one skill's `scripts/`), not a cross-skill dependency
— each DAG step's `command_template` still invokes a single,
independently-runnable script.

Two 1.x shared-library imports are inlined here (dependency-free) because
`clawbio.common.*` is a 1.x-only package tree that does not exist in the
2.0 repo — the 2.0 `clawbio` package is the CLI, not a skill-library:

  - `clawbio.common.parsers` -> `detect_format`/`parse_23andme`/
    `parse_ancestry`/`genotypes_to_simple`/`genotypes_to_positions` below.
    Trimmed to the two formats this skill's SKILL.md documents (23andMe,
    AncestryDNA); the 1.x shared parser also handles MyHeritage/VCF and
    iCloud-Drive staging, neither used by this skill.
  - `clawbio.common.checksums` / `clawbio.common.report` -> `sha256_file`,
    `sha256_hex`, `DISCLAIMER`, `write_result_json` below (same envelope
    shape as 1.x).

NOT ported: `clawbio.common.html_report` (HtmlReportBuilder / self-contained
HTML report) and `clawbio.common.rec_shortener` (guideline-table HTML
parsing used only to build per-drug `source_recs` text during ClinPGx
enrichment). Both are cosmetic/presentational additions on top of the core
CPIC lookup this port preserves faithfully; omitting them keeps this port's
scope to the deterministic gene-calling and drug-recommendation logic, and
matches the report_md + result_json (no HTML) convention already used by
`skills/gwas_lookup` in this repo.
"""

from __future__ import annotations

import csv
import gzip
import hashlib
import json
from dataclasses import dataclass, asdict
from datetime import datetime, timezone
from pathlib import Path
from typing import Any, Optional

# ---------------------------------------------------------------------------
# Disclaimer (clawbio.common.report.DISCLAIMER, inlined)
# ---------------------------------------------------------------------------

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# ---------------------------------------------------------------------------
# Checksums (clawbio.common.checksums, inlined)
# ---------------------------------------------------------------------------


def sha256_file(filepath: str | Path) -> str:
    """Compute SHA-256 checksum of a file (full 64-char hex digest)."""
    h = hashlib.sha256()
    with open(filepath, "rb") as f:
        for chunk in iter(lambda: f.read(8192), b""):
            h.update(chunk)
    return h.hexdigest()


def sha256_hex(filepath: str | Path, length: int = 16) -> str:
    """Truncated SHA-256 hex digest (default 16 chars) — report convention."""
    return sha256_file(filepath)[:length]


# ---------------------------------------------------------------------------
# result.json (clawbio.common.report.write_result_json, inlined)
# ---------------------------------------------------------------------------


def write_result_json(
    output_dir: str | Path,
    skill: str,
    version: str,
    summary: dict[str, Any],
    data: dict[str, Any],
    input_checksum: str = "",
) -> Path:
    """Write the standard ClawBio result.json envelope (same shape as 1.x)."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    envelope = {
        "skill": skill,
        "version": version,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
        "datasets": {},
        "summary": summary,
        "data": data,
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(envelope, indent=2, default=str))
    return result_path


# ---------------------------------------------------------------------------
# Genetic file parsing (clawbio.common.parsers, inlined + trimmed)
# ---------------------------------------------------------------------------


@dataclass
class GenotypeRecord:
    chrom: str = ""
    pos: int = 0
    genotype: str = ""
    allele1: str = ""
    allele2: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


def open_genetic_file(filepath: str | Path):
    """Open a file, handling .gz transparently."""
    filepath = str(filepath)
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    return open(filepath, encoding="utf-8", errors="replace")


def detect_format(filepath: str | Path) -> str:
    """Auto-detect genetic file format from header. Returns "23andme",
    "ancestry", or "unknown" (this skill supports DTC array formats only;
    see SKILL.md's Input Formats table)."""
    filepath = Path(filepath)
    with open_genetic_file(filepath) as f:
        for line in f:
            lower = line.lower()
            if "rsid" in lower and "chromosome" in lower and "genotype" in lower:
                return "23andme"
            if "rsid" in lower and "allele1" in lower:
                return "ancestry"
            if not line.startswith("#"):
                break
    return "unknown"


def parse_23andme(filepath: str | Path) -> dict[str, GenotypeRecord]:
    """Parse 23andMe raw data file. Returns {rsid: GenotypeRecord}."""
    genotypes: dict[str, GenotypeRecord] = {}
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
            allele1 = geno[0] if len(geno) >= 1 else ""
            allele2 = geno[1] if len(geno) >= 2 else ""
            genotypes[rsid] = GenotypeRecord(
                chrom=chrom, pos=pos, genotype=geno, allele1=allele1, allele2=allele2
            )
    return genotypes


def parse_ancestry(filepath: str | Path) -> dict[str, GenotypeRecord]:
    """Parse AncestryDNA raw data file. Returns {rsid: GenotypeRecord}."""
    genotypes: dict[str, GenotypeRecord] = {}
    lines = []
    with open_genetic_file(filepath) as f:
        for line in f:
            if not line.startswith("#"):
                lines.append(line)
    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        rsid = (row.get("rsid") or "").strip()
        allele1 = (row.get("allele1") or "").strip()
        allele2 = (row.get("allele2") or "").strip()
        chrom = (row.get("chromosome") or row.get("chr") or "").strip()
        pos_str = (row.get("position") or row.get("pos") or "0").strip()
        if not rsid.startswith("rs"):
            continue
        geno = allele1 + allele2
        try:
            pos = int(pos_str)
        except ValueError:
            pos = 0
        genotypes[rsid] = GenotypeRecord(
            chrom=chrom, pos=pos, genotype=geno, allele1=allele1, allele2=allele2
        )
    return genotypes


def genotypes_to_simple(records: dict[str, GenotypeRecord]) -> dict[str, str]:
    return {rsid: rec.genotype for rsid, rec in records.items()}


def genotypes_to_positions(records: dict[str, GenotypeRecord]) -> dict[str, dict]:
    return {rsid: {"chrom": rec.chrom, "pos": rec.pos} for rsid, rec in records.items()}


def parse_genetic_file(filepath: str | Path, fmt: str = "auto") -> dict[str, GenotypeRecord]:
    filepath = Path(filepath)
    if fmt == "auto":
        fmt = detect_format(filepath)
    if fmt == "23andme":
        return parse_23andme(filepath)
    if fmt == "ancestry":
        return parse_ancestry(filepath)
    raise ValueError(f"Unknown or undetectable format: {fmt}")


# ---------------------------------------------------------------------------
# Strand utilities
# ---------------------------------------------------------------------------

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def flip_genotype(genotype: str) -> str:
    """Return the complement strand genotype (e.g. 'GA' -> 'CT')."""
    return "".join(COMPLEMENT.get(b, b) for b in genotype)


# ---------------------------------------------------------------------------
# 1. PGx SNP definitions
# ---------------------------------------------------------------------------

PGX_SNPS = {
    # CYP2C19
    "rs4244285":  {"gene": "CYP2C19", "allele": "*2",  "effect": "no_function"},
    "rs4986893":  {"gene": "CYP2C19", "allele": "*3",  "effect": "no_function"},
    "rs12248560": {"gene": "CYP2C19", "allele": "*17", "effect": "increased_function"},
    "rs28399504": {"gene": "CYP2C19", "allele": "*4",  "effect": "no_function"},
    # CYP2D6
    "rs3892097":  {"gene": "CYP2D6", "allele": "*4",  "effect": "no_function"},
    "rs5030655":  {"gene": "CYP2D6", "allele": "*6",  "effect": "no_function"},
    "rs16947":    {"gene": "CYP2D6", "allele": "*2",  "effect": "normal_function"},
    "rs1065852":  {"gene": "CYP2D6", "allele": "*10", "effect": "decreased_function"},
    "rs28371725": {"gene": "CYP2D6", "allele": "*41", "effect": "decreased_function"},
    # CYP2C9
    "rs1799853":  {"gene": "CYP2C9", "allele": "*2", "effect": "decreased_function"},
    "rs1057910":  {"gene": "CYP2C9", "allele": "*3", "effect": "decreased_function"},
    # VKORC1
    "rs9923231":  {"gene": "VKORC1", "allele": "-1639G>A", "effect": "decreased_expression"},
    # SLCO1B1
    "rs4149056":  {"gene": "SLCO1B1", "allele": "*5", "effect": "decreased_function"},
    # DPYD
    "rs3918290":  {"gene": "DPYD", "allele": "*2A",  "effect": "no_function"},
    "rs55886062": {"gene": "DPYD", "allele": "*13",  "effect": "no_function"},
    "rs67376798": {"gene": "DPYD", "allele": "D949V", "effect": "decreased_function"},
    # TPMT
    "rs1800460":  {"gene": "TPMT", "allele": "*3B", "effect": "no_function"},
    "rs1142345":  {"gene": "TPMT", "allele": "*3C", "effect": "no_function"},
    "rs1800462":  {"gene": "TPMT", "allele": "*2",  "effect": "no_function"},
    # UGT1A1
    "rs8175347":  {"gene": "UGT1A1", "allele": "*28", "effect": "decreased_function", "note": "TA-repeat proxy"},
    "rs4148323":  {"gene": "UGT1A1", "allele": "*6",  "effect": "decreased_function"},
    # CYP3A5
    "rs776746":    {"gene": "CYP3A5", "allele": "*3", "effect": "no_function"},
    "rs10264272":  {"gene": "CYP3A5", "allele": "*6", "effect": "no_function"},
    "rs41303343":  {"gene": "CYP3A5", "allele": "*7", "effect": "no_function"},
    # CYP2B6
    "rs3745274":  {"gene": "CYP2B6", "allele": "*9",  "effect": "decreased_function"},
    "rs28399499": {"gene": "CYP2B6", "allele": "*18", "effect": "no_function"},
    # NUDT15
    "rs116855232": {"gene": "NUDT15", "allele": "*3", "effect": "no_function"},
    "rs147390019": {"gene": "NUDT15", "allele": "*2", "effect": "decreased_function"},
    # CYP1A2
    "rs762551":   {"gene": "CYP1A2", "allele": "*1F", "effect": "increased_function"},
    "rs2069514":  {"gene": "CYP1A2", "allele": "*1C", "effect": "decreased_function"},
    # MTHFR
    "rs1801133":  {"gene": "MTHFR", "allele": "677T",  "effect": "decreased_function"},
    "rs1801131":  {"gene": "MTHFR", "allele": "1298C", "effect": "decreased_function"},
}

# ---------------------------------------------------------------------------
# 2. Gene definitions with phenotype rules
# ---------------------------------------------------------------------------

GENE_DEFS = {
    "CYP2C19": {
        "name": "Cytochrome P450 2C19",
        "function": "Drug metabolism",
        "ref": "*1",
        "variants": {
            "rs4244285":  {"allele": "*2",  "alt": "A", "effect": "no_function"},
            "rs4986893":  {"allele": "*3",  "alt": "A", "effect": "no_function"},
            "rs12248560": {"allele": "*17", "alt": "T", "effect": "increased_function"},
            "rs28399504": {"allele": "*4",  "alt": "G", "effect": "no_function"},
        },
        "phenotypes": {
            "Ultrarapid Metabolizer":  ["*17/*17"],
            "Rapid Metabolizer":       ["*1/*17"],
            "Normal Metabolizer":      ["*1/*1"],
            "Intermediate Metabolizer": ["*1/*2", "*1/*3", "*2/*17", "*1/*4"],
            "Poor Metabolizer":        ["*2/*2", "*2/*3", "*3/*3", "*2/*4", "*3/*4", "*4/*4"],
        },
    },
    "CYP2D6": {
        "name": "Cytochrome P450 2D6",
        "function": "Drug metabolism (25% of all drugs)",
        "ref": "*1",
        "variants": {
            "rs3892097":  {"allele": "*4",  "alt": "T", "effect": "no_function"},
            "rs5030655":  {"allele": "*6",  "alt": "C", "effect": "no_function"},
            "rs16947":    {"allele": "*2",  "alt": "A", "effect": "normal_function"},
            "rs1065852":  {"allele": "*10", "alt": "T", "effect": "decreased_function"},
            "rs28371725": {"allele": "*41", "alt": "T", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["*1/*1", "*1/*2", "*2/*2", "*1/*10"],
            "Intermediate Metabolizer": ["*1/*4", "*1/*41", "*2/*41", "*10/*10", "*4/*10", "*10/*41", "*41/*41"],
            "Poor Metabolizer":         ["*4/*4", "*4/*6", "*6/*6", "*4/*41"],
        },
    },
    "CYP2C9": {
        "name": "Cytochrome P450 2C9",
        "function": "Warfarin and NSAID metabolism",
        "ref": "*1",
        "variants": {
            "rs1799853": {"allele": "*2", "alt": "T", "effect": "decreased_function"},
            "rs1057910": {"allele": "*3", "alt": "C", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["*1/*1"],
            "Intermediate Metabolizer": ["*1/*2", "*1/*3", "*2/*2"],
            "Poor Metabolizer":         ["*2/*3", "*3/*3"],
        },
    },
    "VKORC1": {
        "name": "Vitamin K Epoxide Reductase",
        "function": "Warfarin target enzyme",
        "ref": "G",
        "type": "genotype",
        "rsid": "rs9923231",
        "variants": {
            "rs9923231": {"allele": "A", "alt": "A", "effect": "decreased_expression"},
        },
        "phenotypes": {
            "Normal Warfarin Sensitivity":       ["GG", "CC"],
            "Intermediate Warfarin Sensitivity":  ["GA", "AG", "CT", "TC"],
            "High Warfarin Sensitivity":          ["AA", "TT"],
        },
    },
    "SLCO1B1": {
        "name": "Solute Carrier Organic Anion Transporter 1B1",
        "function": "Hepatic statin uptake",
        "ref": "T",
        "type": "genotype",
        "rsid": "rs4149056",
        "variants": {
            "rs4149056": {"allele": "*5", "alt": "C", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal Function":       ["TT"],
            "Decreased Function":    ["TC", "CT"],
            "Poor Function":         ["CC"],
        },
    },
    "DPYD": {
        "name": "Dihydropyrimidine Dehydrogenase",
        "function": "Fluoropyrimidine metabolism",
        "ref": "Normal",
        "type": "dpyd",
        "variants": {
            "rs3918290":  {"allele": "*2A",  "alt": "T", "effect": "no_function"},
            "rs55886062": {"allele": "*13",  "alt": "C", "effect": "no_function"},
            "rs67376798": {"allele": "D949V", "alt": "A", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["Normal/Normal"],
            "Intermediate Metabolizer": ["Normal/*2A", "Normal/*13", "Normal/D949V"],
            "Poor Metabolizer":         ["*2A/*2A", "*2A/*13", "*13/*13"],
        },
    },
    "TPMT": {
        "name": "Thiopurine S-Methyltransferase",
        "function": "Thiopurine metabolism",
        "ref": "*1",
        "variants": {
            "rs1800460": {"allele": "*3B", "alt": "T", "effect": "no_function"},
            "rs1142345": {"allele": "*3C", "alt": "C", "effect": "no_function"},
            "rs1800462": {"allele": "*2",  "alt": "G", "effect": "no_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["*1/*1"],
            "Intermediate Metabolizer": ["*1/*2", "*1/*3A", "*1/*3B", "*1/*3C"],
            "Poor Metabolizer":         ["*2/*2", "*2/*3A", "*3A/*3A", "*3B/*3B", "*3C/*3C",
                                         "*3B/*3C", "*2/*3B", "*2/*3C"],
        },
    },
    "UGT1A1": {
        "name": "UDP-Glucuronosyltransferase 1A1",
        "function": "Irinotecan and bilirubin metabolism",
        "ref": "*1",
        "variants": {
            "rs8175347": {"allele": "*28", "alt": "T", "effect": "decreased_function"},
            "rs4148323": {"allele": "*6",  "alt": "A", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["*1/*1"],
            "Intermediate Metabolizer": ["*1/*28", "*1/*6", "*6/*28"],
            "Poor Metabolizer":         ["*28/*28", "*6/*6"],
        },
    },
    "CYP3A5": {
        "name": "Cytochrome P450 3A5",
        "function": "Tacrolimus metabolism",
        "ref": "*1",
        "variants": {
            "rs776746":   {"allele": "*3", "alt": "G", "effect": "no_function"},
            "rs10264272": {"allele": "*6", "alt": "A", "effect": "no_function"},
            "rs41303343": {"allele": "*7", "alt": "T", "effect": "no_function"},
        },
        "phenotypes": {
            "CYP3A5 Expressor":          ["*1/*1"],
            "Intermediate Expressor":     ["*1/*3", "*1/*6", "*1/*7"],
            "CYP3A5 Non-expressor":       ["*3/*3", "*3/*6", "*6/*6", "*3/*7"],
        },
    },
    "CYP2B6": {
        "name": "Cytochrome P450 2B6",
        "function": "Efavirenz and methadone metabolism",
        "ref": "*1",
        "variants": {
            "rs3745274":  {"allele": "*9",  "alt": "T", "effect": "decreased_function"},
            "rs28399499": {"allele": "*18", "alt": "C", "effect": "no_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["*1/*1"],
            "Intermediate Metabolizer": ["*1/*9", "*1/*18", "*9/*18"],
            "Poor Metabolizer":         ["*9/*9", "*18/*18"],
        },
    },
    "NUDT15": {
        "name": "Nudix Hydrolase 15",
        "function": "Thiopurine metabolism",
        "ref": "*1",
        "variants": {
            "rs116855232": {"allele": "*3", "alt": "T", "effect": "no_function"},
            "rs147390019": {"allele": "*2", "alt": "A", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["*1/*1"],
            "Intermediate Metabolizer": ["*1/*2", "*1/*3"],
            "Poor Metabolizer":         ["*2/*2", "*2/*3", "*3/*3"],
        },
    },
    "CYP1A2": {
        "name": "Cytochrome P450 1A2",
        "function": "Caffeine and clozapine metabolism",
        "ref": "*1",
        "variants": {
            "rs762551":  {"allele": "*1F", "alt": "A", "effect": "increased_function"},
            "rs2069514": {"allele": "*1C", "alt": "A", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Ultrarapid Metabolizer":   ["*1F/*1F"],
            "Normal Metabolizer":       ["*1/*1", "*1/*1F"],
            "Intermediate Metabolizer": ["*1/*1C", "*1C/*1F"],
            "Poor Metabolizer":         ["*1C/*1C"],
        },
    },
    "MTHFR": {
        "name": "Methylenetetrahydrofolate Reductase",
        "function": "Folate metabolism; affects methotrexate toxicity",
        "type": "mthfr",
        "rsid_677": "rs1801133",
        "rsid_1298": "rs1801131",
        "variants": {
            "rs1801133": {"allele": "677T",  "alt": "T", "effect": "decreased_function"},
            "rs1801131": {"allele": "1298C", "alt": "C", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal MTHFR enzyme activity (677CC)":          ["677CC/1298AA", "677CC/1298NOT_TESTED",
                                                              "677NOT_TESTED/1298AA"],
            "Reduced MTHFR enzyme activity (677CT)":         ["677CT/1298AA", "677CC/1298AC"],
            "Strongly reduced MTHFR enzyme activity (677TT)": ["677TT/1298AA", "677CT/1298AC"],
        },
    },
}

# ---------------------------------------------------------------------------
# 3. CPIC drug guidelines
# ---------------------------------------------------------------------------

GUIDELINES = {
    # --- CYP2C19 drugs ---
    "Clopidogrel": {
        "brand": "Plavix", "class": "Antiplatelet Agent", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "standard", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Omeprazole": {
        "brand": "Prilosec", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Pantoprazole": {
        "brand": "Protonix", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Lansoprazole": {
        "brand": "Prevacid", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Esomeprazole": {
        "brand": "Nexium", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Dexlansoprazole": {
        "brand": "Dexilant", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Citalopram": {
        "brand": "Celexa", "class": "SSRI Antidepressant", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Escitalopram": {
        "brand": "Lexapro", "class": "SSRI Antidepressant", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Sertraline": {
        "brand": "Zoloft", "class": "SSRI Antidepressant", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "standard", "rapid_metabolizer": "standard",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "standard",
        },
    },
    "Voriconazole": {
        "brand": "Vfend", "class": "Antifungal", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": "caution", "rapid_metabolizer": "caution",
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    # --- CYP2D6 drugs ---
    "Codeine": {
        "brand": "Tylenol w/ Codeine", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Tramadol": {
        "brand": "Ultram", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Hydrocodone": {
        "brand": "Vicodin", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    "Oxycodone": {
        "brand": "OxyContin", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Tamoxifen": {
        "brand": "Nolvadex", "class": "SERM (Oncology)", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "standard", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Amitriptyline": {
        "brand": "Elavil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Nortriptyline": {
        "brand": "Pamelor", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Desipramine": {
        "brand": "Norpramin", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Imipramine": {
        "brand": "Tofranil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Doxepin": {
        "brand": "Sinequan", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Trimipramine": {
        "brand": "Surmontil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Clomipramine": {
        "brand": "Anafranil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "avoid", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Paroxetine": {
        "brand": "Paxil", "class": "SSRI Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Fluoxetine": {
        "brand": "Prozac", "class": "SSRI Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Venlafaxine": {
        "brand": "Effexor", "class": "SNRI Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "standard", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Metoprolol": {
        "brand": "Lopressor", "class": "Beta-Blocker", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Ondansetron": {
        "brand": "Zofran", "class": "Antiemetic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Risperidone": {
        "brand": "Risperdal", "class": "Antipsychotic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    "Aripiprazole": {
        "brand": "Abilify", "class": "Antipsychotic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "standard", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Haloperidol": {
        "brand": "Haldol", "class": "Antipsychotic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    "Atomoxetine": {
        "brand": "Strattera", "class": "ADHD Medication", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "standard", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    # --- CYP2C9 drugs ---
    "Phenytoin": {
        "brand": "Dilantin", "class": "Antiepileptic", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Celecoxib": {
        "brand": "Celebrex", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Flurbiprofen": {
        "brand": "Ansaid", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Piroxicam": {
        "brand": "Feldene", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Meloxicam": {
        "brand": "Mobic", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    # --- Warfarin (multi-gene) ---
    "Warfarin": {
        "brand": "Coumadin", "class": "Anticoagulant", "genes": ["CYP2C9", "VKORC1"],
        "special": "warfarin",
    },
    # --- SLCO1B1 drugs ---
    "Simvastatin": {
        "brand": "Zocor", "class": "Statin", "gene": "SLCO1B1",
        "recs": {
            "normal_function": "standard",
            "decreased_function": "caution", "poor_function": "avoid",
        },
    },
    "Atorvastatin": {
        "brand": "Lipitor", "class": "Statin", "gene": "SLCO1B1",
        "recs": {
            "normal_function": "standard",
            "decreased_function": "caution", "poor_function": "caution",
        },
    },
    "Rosuvastatin": {
        "brand": "Crestor", "class": "Statin", "gene": "SLCO1B1",
        "recs": {
            "normal_function": "standard",
            "decreased_function": "standard", "poor_function": "standard",
        },
    },
    "Pravastatin": {
        "brand": "Pravachol", "class": "Statin", "gene": "SLCO1B1",
        "recs": {
            "normal_function": "standard",
            "decreased_function": "standard", "poor_function": "standard",
        },
    },
    # --- DPYD drugs ---
    "Fluorouracil": {
        "brand": "5-FU", "class": "Antineoplastic", "gene": "DPYD",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Capecitabine": {
        "brand": "Xeloda", "class": "Antineoplastic", "gene": "DPYD",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    # --- TPMT / NUDT15 drugs ---
    "Azathioprine": {
        "brand": "Imuran", "class": "Immunosuppressant", "gene": "TPMT",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Mercaptopurine": {
        "brand": "Purinethol", "class": "Immunosuppressant", "gene": "TPMT",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Thioguanine": {
        "brand": "Tabloid", "class": "Immunosuppressant", "gene": "TPMT",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    # --- UGT1A1 drug ---
    "Irinotecan": {
        "brand": "Camptosar", "class": "Antineoplastic", "gene": "UGT1A1",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    "Atazanavir": {
        "brand": "Reyataz", "class": "Antiretroviral", "gene": "UGT1A1",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    # --- CYP3A5 drug ---
    "Tacrolimus": {
        "brand": "Prograf", "class": "Immunosuppressant", "gene": "CYP3A5",
        "recs": {
            "extensive_metabolizer": "caution",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "standard",
        },
    },
    # --- CYP2B6 drug ---
    "Efavirenz": {
        "brand": "Sustiva", "class": "Antiretroviral", "gene": "CYP2B6",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    # --- CYP1A2 drugs ---
    "Clozapine": {
        "brand": "Clozaril", "class": "Antipsychotic", "gene": "CYP1A2",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    # --- MTHFR drug ---
    "Methotrexate": {
        "brand": "Rheumatrex / Trexall", "class": "DMARD / Antineoplastic", "gene": "MTHFR",
        "recs": {
            "normal_activity": "standard",
            "intermediate_activity": "caution",
            "reduced_activity": "caution",
        },
    },
    # --- Additional CYP2C9 NSAIDs ---
    "Diclofenac": {
        "brand": "Voltaren", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "standard",
        },
    },
    "Ibuprofen": {
        "brand": "Advil / Motrin", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    "Naproxen": {
        "brand": "Aleve / Naprosyn", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    # --- Additional CYP2D6 drugs ---
    "Dextromethorphan": {
        "brand": "Robitussin DM", "class": "Antitussive", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
    "Propafenone": {
        "brand": "Rythmol", "class": "Antiarrhythmic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": "caution", "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "avoid",
        },
    },
    # --- Additional CYP2B6 drugs ---
    "Methadone": {
        "brand": "Dolophine", "class": "Opioid Analgesic", "gene": "CYP2B6",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "caution", "poor_metabolizer": "caution",
        },
    },
    "Bupropion": {
        "brand": "Wellbutrin / Zyban", "class": "Antidepressant / Smoking Cessation", "gene": "CYP2B6",
        "recs": {
            "normal_metabolizer": "standard",
            "intermediate_metabolizer": "standard", "poor_metabolizer": "caution",
        },
    },
}

# ---------------------------------------------------------------------------
# 3b. Single-drug lookup helpers (dosage_card.py / drug_photo chaining)
# ---------------------------------------------------------------------------


def _levenshtein(s1, s2):
    """Minimal Levenshtein distance (no external deps)."""
    if len(s1) < len(s2):
        return _levenshtein(s2, s1)
    if len(s2) == 0:
        return len(s1)
    prev = range(len(s2) + 1)
    for i, c1 in enumerate(s1):
        curr = [i + 1]
        for j, c2 in enumerate(s2):
            curr.append(min(prev[j + 1] + 1, curr[j] + 1,
                            prev[j] + (0 if c1 == c2 else 1)))
        prev = curr
    return prev[-1]


def resolve_drug_name(query):
    """Resolve a drug query (brand or generic, fuzzy) to a GUIDELINES key.

    Tiers: 1) exact generic, 2) exact brand, 3) substring, 4) Levenshtein <= 2.
    Returns the canonical generic name or None.
    """
    q = query.strip().lower()

    for name in GUIDELINES:
        if name.lower() == q:
            return name

    for name, info in GUIDELINES.items():
        if info.get("brand", "").lower() == q:
            return name

    for name, info in GUIDELINES.items():
        if q in name.lower() or q in info.get("brand", "").lower():
            return name
    for name, info in GUIDELINES.items():
        if name.lower() in q or info.get("brand", "").lower() in q:
            return name

    best, best_dist = None, 3
    for name, info in GUIDELINES.items():
        d = _levenshtein(q, name.lower())
        if d < best_dist:
            best, best_dist = name, d
        d2 = _levenshtein(q, info.get("brand", "").lower())
        if d2 < best_dist:
            best, best_dist = name, d2
    return best


def lookup_single_drug(drug_name, profiles):
    """Look up a single resolved drug against gene profiles."""
    info = GUIDELINES.get(drug_name)
    if not info:
        return None

    if info.get("special") == "warfarin":
        classification, warfarin_note = get_warfarin_rec(profiles)
        cyp2c9 = profiles.get("CYP2C9", {})
        vkorc1 = profiles.get("VKORC1", {})
        result = {
            "drug": drug_name, "brand": info["brand"], "class": info["class"],
            "gene": "CYP2C9 + VKORC1",
            "diplotype": f"CYP2C9 {cyp2c9.get('diplotype', '?')} / VKORC1 {vkorc1.get('diplotype', '?')}",
            "phenotype": f"CYP2C9 {cyp2c9.get('phenotype', '?')} / VKORC1 {vkorc1.get('phenotype', '?')}",
            "classification": classification,
        }
        if warfarin_note:
            result["note"] = warfarin_note
        return result

    gene = info["gene"]
    if gene not in profiles:
        return {
            "drug": drug_name, "brand": info["brand"], "class": info["class"],
            "gene": gene, "diplotype": "NOT_TESTED", "phenotype": "Indeterminate",
            "classification": "indeterminate",
        }

    prof = profiles[gene]
    pheno_key = phenotype_to_key(prof["phenotype"])
    recs = info.get("recs", {})
    classification = recs.get(pheno_key, "indeterminate")

    return {
        "drug": drug_name, "brand": info["brand"], "class": info["class"],
        "gene": gene, "diplotype": prof["diplotype"], "phenotype": prof["phenotype"],
        "classification": classification,
    }


_CLASS_LABELS = {
    "standard": "STANDARD DOSING",
    "caution": "USE WITH CAUTION",
    "avoid": "AVOID — DO NOT USE",
    "indeterminate": "INSUFFICIENT DATA",
}


def format_dosage_card(result, visible_dose=None):
    """Format a single-drug lookup result as a visual text card."""
    cls_label = _CLASS_LABELS.get(result["classification"], result["classification"].upper())
    bar = "━" * 35

    cl = result["classification"]
    _CLS_TEXT = {
        "standard": "Standard dosing expected to be effective.",
        "caution": "Dose adjustment or monitoring may be needed.",
        "avoid": "Consider alternative medication.",
        "indeterminate": "Insufficient data for recommendation.",
    }
    rec_text = result.get("note") or _CLS_TEXT.get(cl, "")
    if visible_dose:
        if cl == "standard":
            rec_text = f"Your genotype supports {result['drug']} {visible_dose} as prescribed."
        elif cl == "caution":
            rec_text = f"{visible_dose} may need adjustment."
        elif cl == "avoid":
            rec_text = f"Your genotype contraindicates {result['drug']} {visible_dose}."

    words = rec_text.split()
    rec_lines = []
    current = "  "
    for w in words:
        if len(current) + len(w) + 1 > 44:
            rec_lines.append(current)
            current = "  " + w
        else:
            current += (" " if len(current) > 2 else "") + w
    if current.strip():
        rec_lines.append(current)
    rec_block = "\n".join(rec_lines)

    card = f"""{bar}
  DRUG PHOTO ANALYSIS
{bar}

  Identified: {result['drug']} ({result['brand']})
  Class: {result['class']}

  YOUR GENETIC PROFILE
  Gene: {result['gene']}
  Diplotype: {result['diplotype']}
  Phenotype: {result['phenotype']}

  RECOMMENDATION: {cls_label}
{rec_block}

  Source: FDA Table of Pharmacogenomic
  Biomarkers in Drug Labeling & CPIC
  Guidelines (cpicpgx.org)

  DISCLAIMER: Research/educational use only.
  Consult a healthcare professional.
{bar}"""
    return card


# ---------------------------------------------------------------------------
# 4. Reference genome build detection
# ---------------------------------------------------------------------------

# Source: Ensembl release 112 / dbSNP b156. Verified 2026-07-30. Both builds
# are stored so a GRCh37 file is positively identified rather than inferred
# from "did not match GRCh38", which cannot distinguish an older build from
# a corrupt one.
_PGX_BUILD_POSITIONS = {
    "rs4244285":  {"chrom": "10", "GRCh37": 96541616, "GRCh38": 94781859},  # CYP2C19*2
    "rs3892097":  {"chrom": "22", "GRCh37": 42524947, "GRCh38": 42128945},  # CYP2D6*4
    "rs1799853":  {"chrom": "10", "GRCh37": 96702047, "GRCh38": 94942290},  # CYP2C9*2
    "rs9923231":  {"chrom": "16", "GRCh37": 31107689, "GRCh38": 31096368},  # VKORC1 -1639G>A
    "rs1801133":  {"chrom": "1",  "GRCh37": 11856378, "GRCh38": 11796321},  # MTHFR C677T
}

_POS_TOLERANCE = 0
_BUILDS = ("GRCh38", "GRCh37")


def _normalise_chrom(value) -> Optional[str]:
    """Normalise 'chr10', '10', 10 -> '10'. Returns None if unusable."""
    if value is None:
        return None
    text = str(value).strip()
    if text.lower().startswith("chr"):
        text = text[3:]
    return text or None


def detect_reference_genome(positions: dict) -> Optional[str]:
    """Identify the genome build the input coordinates are annotated against.

    Returns "GRCh38", "GRCh37", "unknown_build", or None (no usable panel
    position). Calls themselves are rsID-based, not coordinate-based, so
    GRCh37 is fully callable — only `unknown_build` (coordinates match
    neither build) withholds gene calls.
    """
    votes = {build: 0 for build in _BUILDS}
    checked = 0

    for rsid, entry in _PGX_BUILD_POSITIONS.items():
        if rsid not in positions:
            continue
        record = positions[rsid] or {}
        actual_pos = record.get("pos")
        if actual_pos is None or actual_pos <= 0:
            continue

        actual_chrom = _normalise_chrom(record.get("chrom"))
        if actual_chrom is not None and actual_chrom != entry["chrom"]:
            checked += 1
            continue

        checked += 1
        for build in _BUILDS:
            if abs(actual_pos - entry[build]) <= _POS_TOLERANCE:
                votes[build] += 1

    if checked == 0:
        return None

    best = max(votes, key=lambda build: votes[build])
    if votes[best] == 0:
        return "unknown_build"
    runner_up = max(votes[b] for b in _BUILDS if b != best)
    if votes[best] == runner_up:
        return "unknown_build"
    return best


def parse_file(path):
    """Parse a genetic data file and extract PGx-relevant SNPs.

    Returns (fmt, total_snps, pgx_dict, ref_genome).
    """
    fmt = detect_format(path)
    records = parse_genetic_file(str(path), fmt=fmt) if fmt != "unknown" else {}
    snps = genotypes_to_simple(records)
    positions = genotypes_to_positions(records)
    snps = {rsid: gt.upper() for rsid, gt in snps.items() if gt and gt not in ("--", "00")}

    ref_genome = detect_reference_genome(positions)

    pgx = {}
    for rsid, info in PGX_SNPS.items():
        if rsid in snps:
            pgx[rsid] = {"genotype": snps[rsid], **info}

    return fmt, len(snps), pgx, ref_genome


# ---------------------------------------------------------------------------
# 5. Star allele caller
# ---------------------------------------------------------------------------

def call_diplotype(gene, pgx_snps):
    gdef = GENE_DEFS[gene]

    if gdef.get("type") == "mthfr":
        rsid_677 = gdef["rsid_677"]
        rsid_1298 = gdef["rsid_1298"]
        gt_677 = pgx_snps[rsid_677]["genotype"] if rsid_677 in pgx_snps else "NOT_TESTED"
        gt_1298 = pgx_snps[rsid_1298]["genotype"] if rsid_1298 in pgx_snps else "NOT_TESTED"
        if gt_677 == "NOT_TESTED" and gt_1298 == "NOT_TESTED":
            return "NOT_TESTED"
        if gt_677 != "NOT_TESTED":
            alleles_677 = set(gt_677.upper())
            on_plus_strand = alleles_677 <= {"G", "A"}
            if on_plus_strand:
                gt_677 = "".join(sorted(flip_genotype(gt_677)))
            else:
                gt_677 = "".join(sorted(gt_677.upper()))
        if gt_1298 != "NOT_TESTED":
            alleles_1298 = set(gt_1298.upper())
            on_plus_strand = alleles_1298 <= {"T", "G"}
            if on_plus_strand:
                gt_1298 = "".join(sorted(flip_genotype(gt_1298)))
            else:
                gt_1298 = "".join(sorted(gt_1298.upper()))
        return f"677{gt_677}/1298{gt_1298}"

    if gdef.get("type") == "genotype":
        rsid = gdef["rsid"]
        if rsid in pgx_snps:
            return pgx_snps[rsid]["genotype"]
        return "NOT_TESTED"

    sv_rsids = {r for r, v in gdef["variants"].items() if v["alt"].upper() in ("DEL", "INS", "TA7")}
    gene_rsids = [r for r in gdef["variants"].keys() if r not in sv_rsids]
    tested = [r for r in gene_rsids if r in pgx_snps]

    if not tested:
        sv_in_data = [r for r in sv_rsids if r in pgx_snps]
        if sv_in_data:
            return "NOT_TESTED"
        return "NOT_TESTED"

    detected = []
    sv_untestable = []
    for rsid, vdef in gdef["variants"].items():
        if rsid in pgx_snps:
            gt = pgx_snps[rsid]["genotype"]
            alt = vdef["alt"].upper()
            if alt in ("DEL", "INS", "TA7"):
                is_het = len(set(gt)) > 1
                if is_het:
                    sv_untestable.append({"rsid": rsid, "allele": vdef["allele"], "alt": alt})
                continue
            alt_count = gt.count(alt)
            if alt_count > 0:
                detected.append({"rsid": rsid, "allele": vdef["allele"],
                                 "copies": alt_count, "effect": vdef["effect"]})

    if sv_untestable:
        sv_desc = ", ".join(f"{s['allele']}({s['rsid']})" for s in sv_untestable)
        return f"Indeterminate (structural variant not assessed: {sv_desc})"

    if gdef.get("type") == "dpyd":
        if not detected:
            if len(tested) == len(gene_rsids):
                return "Normal/Normal"
            return f"Normal/Normal ({len(tested)}/{len(gene_rsids)} SNPs tested)"
        v = detected[0]
        if v["copies"] == 2:
            return f"{v['allele']}/{v['allele']}"
        return f"Normal/{v['allele']}"

    if not detected:
        if len(tested) == len(gene_rsids):
            return f"{gdef['ref']}/{gdef['ref']}"
        return f"{gdef['ref']}/{gdef['ref']} ({len(tested)}/{len(gene_rsids)} SNPs tested)"

    detected.sort(key=lambda v: (0 if v["effect"] == "no_function" else 1))

    het_variants = [v for v in detected if v["copies"] == 1]
    if len(het_variants) >= 2:
        has_no_function = any(v["effect"] == "no_function" for v in het_variants)
        if has_no_function:
            allele_desc = " + ".join(f"{v['allele']}({v['rsid']})" for v in het_variants)
            return f"Indeterminate (phase ambiguity: {allele_desc})"

    a1_parts, a2_parts = [], []
    for v in detected:
        if v["copies"] == 2:
            a1_parts.append(v["allele"])
            a2_parts.append(v["allele"])
        elif v["copies"] == 1:
            if not a1_parts:
                a1_parts.append(v["allele"])
            elif not a2_parts:
                a2_parts.append(v["allele"])

    a1 = a1_parts[0] if a1_parts else gdef["ref"]
    a2 = a2_parts[0] if a2_parts else gdef["ref"]
    alleles = sorted([a1, a2])
    return "/".join(alleles)


def call_phenotype(gene, diplotype):
    import re as _re

    if diplotype == "NOT_TESTED":
        return "Indeterminate (not genotyped)"

    if diplotype.startswith("Indeterminate"):
        return diplotype

    gdef = GENE_DEFS[gene]
    norm = diplotype.upper()

    has_partial = "(" in norm and "SNPS TESTED" in norm
    match_str = norm.split("(")[0].strip()

    if has_partial:
        cov_match = _re.search(r"\((\d+)/(\d+)\s+SNPS", norm)
        if cov_match:
            tested, total = int(cov_match.group(1)), int(cov_match.group(2))
            if tested < total:
                return f"Indeterminate (incomplete coverage: {tested}/{total} SNPs tested)"

    for desc, conditions in gdef["phenotypes"].items():
        for cond in conditions:
            if match_str == cond.upper():
                return desc
            parts = cond.split("/")
            if len(parts) == 2 and match_str == f"{parts[1]}/{parts[0]}".upper():
                return desc
    return f"Unknown (unmapped diplotype: {diplotype})"


# ---------------------------------------------------------------------------
# 6. Drug recommendation lookup
# ---------------------------------------------------------------------------

def phenotype_to_key(phenotype_desc):
    """Map phenotype description to GUIDELINES rec key."""
    mapping = {
        "Normal Metabolizer": "normal_metabolizer",
        "Intermediate Metabolizer": "intermediate_metabolizer",
        "Poor Metabolizer": "poor_metabolizer",
        "Rapid Metabolizer": "rapid_metabolizer",
        "Ultrarapid Metabolizer": "ultrarapid_metabolizer",
        "Normal Warfarin Sensitivity": "normal_warfarin_sensitivity",
        "Intermediate Warfarin Sensitivity": "intermediate_warfarin_sensitivity",
        "High Warfarin Sensitivity": "high_warfarin_sensitivity",
        "Normal Function": "normal_function",
        "Decreased Function": "decreased_function",
        "Intermediate Function": "intermediate_function",
        "Poor Function": "poor_function",
        "CYP3A5 Expressor": "extensive_metabolizer",
        "Intermediate Expressor": "intermediate_metabolizer",
        "CYP3A5 Non-expressor": "poor_metabolizer",
        "Normal MTHFR enzyme activity (677CC)": "normal_activity",
        "Reduced MTHFR enzyme activity (677CT)": "intermediate_activity",
        "Strongly reduced MTHFR enzyme activity (677TT)": "reduced_activity",
    }
    key = mapping.get(phenotype_desc)
    if key:
        return key
    stripped = phenotype_desc.split("(")[0].strip() if "(" in phenotype_desc else phenotype_desc
    key = mapping.get(stripped)
    if key:
        return key
    for label, val in mapping.items():
        if label.startswith(stripped):
            return val
    return "indeterminate"


def get_warfarin_rec(profiles):
    cyp2c9_data = profiles.get("CYP2C9", {})
    vkorc1_data = profiles.get("VKORC1", {})
    cyp2c9 = cyp2c9_data.get("phenotype", "")
    vkorc1 = vkorc1_data.get("phenotype", "")

    if "indeterminate" in cyp2c9.lower() or "not genotyped" in cyp2c9.lower() or not cyp2c9:
        return "indeterminate", "CYP2C9 not genotyped. Cannot provide genotype-guided warfarin dosing. Clinical testing recommended."
    if "indeterminate" in vkorc1.lower() or "not genotyped" in vkorc1.lower() or not vkorc1:
        return "indeterminate", "VKORC1 not genotyped. Cannot provide genotype-guided warfarin dosing. Clinical testing recommended."
    if "unknown" in cyp2c9.lower() or "unknown" in vkorc1.lower():
        return "indeterminate", "CYP2C9/VKORC1 phenotype could not be determined. Clinical testing recommended."

    cyp2c9_normal = "normal" in cyp2c9.lower()
    vkorc1_normal = "normal" in vkorc1.lower()

    if cyp2c9_normal and vkorc1_normal:
        return "standard", None
    elif "poor" in cyp2c9.lower() or "high" in vkorc1.lower():
        return "avoid", None
    else:
        return "caution", None


def lookup_drugs(profiles):
    results = {"standard": [], "caution": [], "avoid": [], "indeterminate": []}

    for drug_name, drug in GUIDELINES.items():
        if drug.get("special") == "warfarin":
            classification, warfarin_note = get_warfarin_rec(profiles)
            entry = {
                "drug": drug_name, "brand": drug["brand"],
                "class": drug["class"], "gene": "CYP2C9+VKORC1",
                "classification": classification,
            }
            if warfarin_note:
                entry["note"] = warfarin_note
            results.setdefault(classification, []).append(entry)
            continue

        gene = drug["gene"]
        if gene not in profiles:
            results["indeterminate"].append({
                "drug": drug_name, "brand": drug["brand"],
                "class": drug["class"], "gene": gene,
                "classification": "indeterminate",
            })
            continue

        pheno_key = phenotype_to_key(profiles[gene]["phenotype"])

        if pheno_key == "indeterminate":
            results["indeterminate"].append({
                "drug": drug_name, "brand": drug["brand"],
                "class": drug["class"], "gene": gene,
                "classification": "indeterminate",
            })
            continue

        recs = drug.get("recs", {})
        classification = recs.get(pheno_key, "indeterminate")

        results.setdefault(classification, []).append({
            "drug": drug_name, "brand": drug["brand"],
            "class": drug["class"], "gene": gene,
            "classification": classification,
        })

    return results
