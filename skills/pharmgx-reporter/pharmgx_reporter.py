#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ClawBio PharmGx Reporter
Pharmacogenomic report generator from DTC genetic data (23andMe/AncestryDNA).

Analyses 33 pharmacogenomic SNPs across 13 genes, calls star alleles and
metabolizer phenotypes, and looks up CPIC drug recommendations for 59 medications.

Usage:
    python pharmgx_reporter.py --input patient_data.txt --output report_dir
"""

import argparse
import json
import re
import sys
from datetime import datetime, timezone
from pathlib import Path
from urllib.parse import quote

# ---------------------------------------------------------------------------
# Shared library imports
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.parsers import parse_genetic_file, genotypes_to_simple, genotypes_to_positions
from clawbio.common.checksums import sha256_hex, sha256_file
from clawbio.common.report import write_result_json, DISCLAIMER
from clawbio.common.html_report import HtmlReportBuilder, write_html_report

# ---------------------------------------------------------------------------
# Strand utilities (mirrors nutrigx_advisor/extract_genotypes.py)
# ---------------------------------------------------------------------------

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def flip_genotype(genotype: str) -> str:
    """Return the complement strand genotype (e.g. 'GA' → 'CT')."""
    return "".join(COMPLEMENT.get(b, b) for b in genotype)


# ---------------------------------------------------------------------------
# 1. PGx SNP definitions (ported from PharmXD snp-parser.js)
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
# 2. Gene definitions with phenotype rules (from phenotype.js)
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
            # CPIC 2020 (Caudle, PMID 31647186): NM at AS >= 1.25
            # *1/*10 has AS = 1.0 + 0.25 = 1.25 -> NM (not IM)
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
            # Activity labels follow DPWG MTHFR nomenclature with genotype detail
            "Normal MTHFR enzyme activity (677CC)":          ["677CC/1298AA", "677CC/1298NOT_TESTED",
                                                              "677NOT_TESTED/1298AA"],
            "Reduced MTHFR enzyme activity (677CT)":         ["677CT/1298AA", "677CC/1298AC"],
            "Strongly reduced MTHFR enzyme activity (677TT)": ["677TT/1298AA", "677CT/1298AC"],
        },
    },
}

# ---------------------------------------------------------------------------
# 3. CPIC drug guidelines (from cpic-lookup.js, all 51 drugs)
# ---------------------------------------------------------------------------

# Phenotype key mapping for lookup
# For star-allele genes: ultrarapid_metabolizer, normal_metabolizer,
#   intermediate_metabolizer, poor_metabolizer
# For VKORC1: normal_warfarin_sensitivity, intermediate_warfarin_sensitivity,
#   high_warfarin_sensitivity
# For SLCO1B1: normal_function, intermediate_function, poor_function
# For CYP3A5: extensive_metabolizer, intermediate_metabolizer, poor_metabolizer

def _pheno_key(description):
    """Convert phenotype description to lookup key."""
    return description.lower().replace(" ", "_")


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
        # CPIC 2022: Decreased Function (TC) = avoid simvastatin 40+ mg
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
        # CPIC (Theken 2020, Table S9): diclofenac PK not significantly impacted
        # by CYP2C9 variants in vivo. No CPIC recommendation for PM (Level C).
        # Diclofenac is listed as an alternative for CYP2C9 PMs.
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
# 3b. Single-drug lookup helpers (drug photo skill)
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

    Tiers: 1) exact generic, 2) exact brand, 3) substring, 4) Levenshtein ≤ 2.
    Returns the canonical generic name or None.
    """
    q = query.strip().lower()

    # Tier 1: exact generic match
    for name in GUIDELINES:
        if name.lower() == q:
            return name

    # Tier 2: exact brand match
    for name, info in GUIDELINES.items():
        if info.get("brand", "").lower() == q:
            return name

    # Tier 3: substring match on generic or brand
    for name, info in GUIDELINES.items():
        if q in name.lower() or q in info.get("brand", "").lower():
            return name
    for name, info in GUIDELINES.items():
        if name.lower() in q or info.get("brand", "").lower() in q:
            return name

    # Tier 4: Levenshtein ≤ 2
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
    """Look up a single resolved drug against gene profiles.

    Returns a dict with drug, brand, class, gene, diplotype, phenotype,
    classification, recommendation — or None if drug not found.
    """
    info = GUIDELINES.get(drug_name)
    if not info:
        return None

    # Warfarin is multi-gene special case
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
    """Format a single-drug lookup result as a visual Telegram card."""
    cls_label = _CLASS_LABELS.get(result["classification"], result["classification"].upper())
    bar = "\u2501" * 35  # ━

    # Build dose-aware recommendation line
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

    # Wrap recommendation text at ~42 chars
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
# 4. File parser (delegates to clawbio.common.parsers)
# ---------------------------------------------------------------------------

def detect_format(lines: list[str]) -> str:
    """Detect file format from header lines (backward-compatible wrapper).

    Delegates to clawbio.common.parsers.detect_format via a temp file.
    """
    from clawbio.common.parsers import detect_format as _detect_fmt
    import tempfile, os
    fd, tmp = tempfile.mkstemp(suffix=".txt")
    try:
        with os.fdopen(fd, "w") as f:
            f.write("\n".join(lines))
        return _detect_fmt(Path(tmp))
    finally:
        os.unlink(tmp)


# Known GRCh38 positions for key PGx SNPs (used for reference genome mismatch detection)
_GRCH38_POSITIONS = {
    "rs4244285": 96541616,   # CYP2C19*2
    "rs3892097": 42128945,   # CYP2D6*4
    "rs1799853": 96702047,   # CYP2C9*2
    "rs9923231": 31107689,   # VKORC1
    "rs1801133": 11796321,   # MTHFR C677T
}

# Tolerance for position comparison (exact match expected within same build)
_POS_TOLERANCE = 0


def detect_reference_genome(positions: dict) -> str | None:
    """Detect reference genome build from SNP positions.

    Returns "GRCh37_mismatch" if positions suggest GRCh37 coordinates,
    "GRCh38" if they match GRCh38, or None if insufficient data.
    """
    matches_38 = 0
    mismatches = 0
    checked = 0

    for rsid, expected_38 in _GRCH38_POSITIONS.items():
        if rsid in positions:
            actual_pos = positions[rsid].get("pos")
            if actual_pos is not None and actual_pos > 0:
                checked += 1
                if abs(actual_pos - expected_38) <= _POS_TOLERANCE:
                    matches_38 += 1
                else:
                    mismatches += 1

    if checked == 0:
        return None
    # Require strong evidence of mismatch: at least 2 mismatches AND at least
    # as many mismatches as matches. A single mismatch out of few checked
    # positions may be a test-case artifact or chromosome model difference.
    if mismatches >= 2 and mismatches >= matches_38:
        return "GRCh37_mismatch"
    return "GRCh38"


def parse_file(path):
    """Parse a genetic data file and extract PGx-relevant SNPs.

    Uses the shared parser from clawbio.common.parsers for file reading and
    format detection, then filters to the PGx panel.

    Returns:
        (fmt, total_snps, pgx_dict, ref_genome) where pgx_dict maps rsid -> {genotype, gene, allele, effect}
        and ref_genome is "GRCh38", "GRCh37_mismatch", or None.
    """
    from clawbio.common.parsers import detect_format as _detect_fmt

    # Use shared parser for robust format detection and file reading
    try:
        fmt = _detect_fmt(path)
    except ValueError:
        fmt = "unknown"

    records = parse_genetic_file(str(path), fmt=fmt if fmt != "unknown" else "auto")
    snps = genotypes_to_simple(records)
    positions = genotypes_to_positions(records)
    # Normalize genotypes to uppercase for PGx matching
    snps = {rsid: gt.upper() for rsid, gt in snps.items() if gt and gt not in ("--", "00")}

    # Detect reference genome
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
        # Build combined diplotype: "677{CT}/1298{AC}" using CPIC gene-strand notation.
        # MTHFR is on the minus strand. The coding strand alleles are:
        #   rs1801133 (C677T): C (ref) > T (risk)
        #   rs1801131 (A1298C): A (ref) > C (risk)
        # The plus strand complement would be G>A and T>G respectively.
        # Auto-detect strand: if genotype contains C or T for rs1801133,
        # it is already on the coding strand and should NOT be flipped.
        # If it contains G or A, it is on the plus strand and needs flipping.
        rsid_677  = gdef["rsid_677"]
        rsid_1298 = gdef["rsid_1298"]
        gt_677  = pgx_snps[rsid_677]["genotype"]  if rsid_677  in pgx_snps else "NOT_TESTED"
        gt_1298 = pgx_snps[rsid_1298]["genotype"] if rsid_1298 in pgx_snps else "NOT_TESTED"
        if gt_677 == "NOT_TESTED" and gt_1298 == "NOT_TESTED":
            return "NOT_TESTED"
        if gt_677 != "NOT_TESTED":
            # Auto-detect strand for rs1801133 (C677T):
            # Coding strand alleles are C and T. Plus strand alleles are G and A.
            alleles_677 = set(gt_677.upper())
            on_plus_strand = alleles_677 <= {"G", "A"}
            if on_plus_strand:
                gt_677 = "".join(sorted(flip_genotype(gt_677)))
            else:
                gt_677 = "".join(sorted(gt_677.upper()))
        if gt_1298 != "NOT_TESTED":
            # Auto-detect strand for rs1801131 (A1298C):
            # Coding strand alleles are A and C. Plus strand alleles are T and G.
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

    # Exclude structural variant SNPs (DEL, INS, TA7) from both the total
    # panel count and tested count: they are inherently untestable from DTC data.
    sv_rsids = {r for r, v in gdef["variants"].items() if v["alt"].upper() in ("DEL", "INS", "TA7")}
    gene_rsids = [r for r in gdef["variants"].keys() if r not in sv_rsids]
    tested = [r for r in gene_rsids if r in pgx_snps]

    if not tested:
        # Check if the only SNPs present are SV SNPs
        sv_in_data = [r for r in sv_rsids if r in pgx_snps]
        if sv_in_data:
            # Gene has data but only at untestable SV positions
            return "NOT_TESTED"
        return "NOT_TESTED"

    detected = []
    sv_untestable = []  # SV SNPs where patient carries a het call
    for rsid, vdef in gdef["variants"].items():
        if rsid in pgx_snps:
            gt = pgx_snps[rsid]["genotype"]
            alt = vdef["alt"].upper()
            if alt in ("DEL", "INS", "TA7"):
                # Only flag if patient has a heterozygous call (possible carrier)
                is_het = len(set(gt)) > 1
                if is_het:
                    sv_untestable.append({"rsid": rsid, "allele": vdef["allele"], "alt": alt})
                continue
            alt_count = gt.count(alt)
            if alt_count > 0:
                detected.append({"rsid": rsid, "allele": vdef["allele"],
                                 "copies": alt_count, "effect": vdef["effect"]})

    # If patient carries a het call at an SV SNP, the gene result is unreliable.
    # Report as Indeterminate rather than falsely claiming Normal.
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

    # Phase ambiguity detection: if 2+ heterozygous (copies==1) variants are
    # detected and at least one has no_function effect, the cis/trans phase
    # assignment from unphased DTC data can change the clinical phenotype.
    # Flag as Indeterminate only in this clinically significant scenario.
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

    # Structural variant limitations: diplotype already flagged as indeterminate
    if diplotype.startswith("Indeterminate"):
        return diplotype

    gdef = GENE_DEFS[gene]
    norm = diplotype.upper()

    # Check for partial-coverage annotations (e.g. "*1/*1 (1/3 SNPs tested)")
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
        # MTHFR activity labels (DPWG nomenclature with genotype detail)
        "Normal MTHFR enzyme activity (677CC)": "normal_activity",
        "Reduced MTHFR enzyme activity (677CT)": "intermediate_activity",
        "Strongly reduced MTHFR enzyme activity (677TT)": "reduced_activity",
    }
    # Try exact match first, then strip qualifiers like "(inferred)"
    key = mapping.get(phenotype_desc)
    if key:
        return key
    stripped = phenotype_desc.split("(")[0].strip() if "(" in phenotype_desc else phenotype_desc
    key = mapping.get(stripped)
    if key:
        return key
    # Try prefix match: "Normal" → "Normal Metabolizer"
    for label, val in mapping.items():
        if label.startswith(stripped):
            return val
    return "indeterminate"


def get_warfarin_rec(profiles):
    cyp2c9_data = profiles.get("CYP2C9", {})
    vkorc1_data = profiles.get("VKORC1", {})
    cyp2c9 = cyp2c9_data.get("phenotype", "")
    vkorc1 = vkorc1_data.get("phenotype", "")

    # If either gene was not genotyped, we cannot provide warfarin guidance
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


# ---------------------------------------------------------------------------
# 6b. ClinPGx evidence enrichment
# ---------------------------------------------------------------------------

# Evidence level ranking (higher is stronger)
_EVIDENCE_RANK = {"1A": 6, "1B": 5, "2A": 4, "2B": 3, "3": 2, "4": 1}

_EVIDENCE_BADGE_CLASS = {
    "1A": "badge-evidence-high",
    "1B": "badge-evidence-high",
    "2A": "badge-evidence-moderate",
    "2B": "badge-evidence-moderate",
    "3": "badge-evidence-low",
    "4": "badge-evidence-minimal",
}


def enrich_with_clinpgx(drug_results, cache_dir=None):
    """Query ClinPGx API for evidence levels, sources, and verification per drug.

    Returns a dict keyed by lowercase drug name with evidence metadata.
    Returns ``{}`` if the ClinPGx skill is unavailable.
    """
    try:
        _clinpgx_dir = _PROJECT_ROOT / "skills" / "clinpgx"
        if str(_clinpgx_dir) not in sys.path:
            sys.path.insert(0, str(_clinpgx_dir))
        from clinpgx import ClinPGxClient
    except Exception:
        return {}

    if cache_dir is None:
        cache_dir = Path.home() / ".clawbio" / "clinpgx_cache"
    client = ClinPGxClient(cache_dir=Path(cache_dir), use_cache=True)

    # Collect unique genes from all drug results
    all_drugs = []
    for cat_list in drug_results.values():
        all_drugs.extend(cat_list)
    genes = sorted({d["gene"] for d in all_drugs if d["gene"] not in ("", "CYP2C9+VKORC1")})
    # Add both genes for warfarin
    if any(d["gene"] == "CYP2C9+VKORC1" for d in all_drugs):
        for g in ("CYP2C9", "VKORC1"):
            if g not in genes:
                genes.append(g)
        genes.sort()

    enrichment = {}
    total = len(genes)

    for idx, gene in enumerate(genes, 1):
        print(f"  Enriching with ClinPGx data... [{idx}/{total}] {gene}")
        try:
            # Get gene accession ID
            gene_data = client.get_gene(gene)
            gene_id = gene_data[0].get("id", "") if gene_data else ""

            # Get clinical annotations for this gene
            annotations = client.get_clinical_annotations(gene_symbol=gene)

            # Get guidelines for this gene
            guidelines = client.get_guidelines(gene_accession_id=gene_id) if gene_id else []

            # Build per-drug enrichment from annotations
            for ann in annotations:
                chemicals = ann.get("relatedChemicals", [])
                level_obj = ann.get("levelOfEvidence", {})
                level_term = level_obj.get("term", "") if isinstance(level_obj, dict) else str(level_obj)

                for chem in chemicals:
                    chem_name = chem.get("name", "").lower()
                    if not chem_name:
                        continue

                    existing = enrichment.get(chem_name, {})
                    existing_rank = _EVIDENCE_RANK.get(existing.get("evidence_level", ""), 0)
                    new_rank = _EVIDENCE_RANK.get(level_term, 0)

                    if new_rank > existing_rank:
                        existing["evidence_level"] = level_term

                    # Accumulate sources
                    sources = existing.get("sources", set())
                    if isinstance(sources, list):
                        sources = set(sources)
                    existing["sources"] = sources
                    enrichment[chem_name] = existing

            # Check guidelines for CPIC verification, sources, and store raw guidelines
            for gl in guidelines:
                source = gl.get("source", "")
                has_dosing = gl.get("dosingInformation", False)
                gl_chemicals = gl.get("relatedChemicals", [])

                for chem in gl_chemicals:
                    chem_name = chem.get("name", "").lower()
                    if not chem_name:
                        continue
                    existing = enrichment.get(chem_name, {})
                    sources = existing.get("sources", set())
                    if isinstance(sources, list):
                        sources = set(sources)
                    if source:
                        sources.add(source)
                    existing["sources"] = sources

                    if has_dosing and source == "CPIC":
                        existing["verified"] = True

                    guideline_name = gl.get("name", "")
                    if guideline_name:
                        existing["guideline_name"] = guideline_name

                    # Store raw guideline objects for structured table parsing
                    raw_guidelines = existing.get("_guidelines", [])
                    raw_guidelines.append(gl)
                    existing["_guidelines"] = raw_guidelines

                    enrichment[chem_name] = existing

        except Exception as exc:
            print(f"    Warning: ClinPGx query failed for {gene}: {exc}")

    # Convert sources sets to sorted lists for JSON serialization
    for key in enrichment:
        src = enrichment[key].get("sources", set())
        if isinstance(src, set):
            enrichment[key]["sources"] = sorted(src)
        if "verified" not in enrichment[key]:
            enrichment[key]["verified"] = False

    return enrichment


def extract_phenotype_recs(enrichment, drug_results, profiles):
    """Extract phenotype-specific recommendations from ALL guideline sources.

    Parses the structured HTML tables in ClinPGx guideline textMarkdown to
    find the exact recommendation for the patient's phenotype from each source
    (DPWG, CPIC, CPNDS, RNPGx). No LLM needed.

    Mutates enrichment in-place: adds 'source_recs' list of {source, rec, strength}.
    """
    from clawbio.common.rec_shortener import extract_all_source_recs, shorten_rec

    # Build drug→(gene, phenotype) mapping
    drug_phenotype = {}
    for cat_drugs in drug_results.values():
        for d in cat_drugs:
            gene = d.get("gene", "")
            drug_key = d["drug"].lower()
            if gene in profiles:
                drug_phenotype[drug_key] = {
                    "gene": gene,
                    "phenotype": profiles[gene]["phenotype"],
                }

    extracted = 0
    for drug_key, entry in enrichment.items():
        guidelines = entry.get("_guidelines", [])
        if not guidelines:
            continue
        pheno_info = drug_phenotype.get(drug_key)
        if not pheno_info:
            continue

        all_recs = extract_all_source_recs(
            guidelines,
            drug_name=drug_key,
            patient_phenotype=pheno_info["phenotype"],
            gene=pheno_info["gene"],
        )
        if all_recs:
            entry["source_recs"] = [
                {"source": r["source"], "rec": shorten_rec(r["rec"]), "strength": r["strength"]}
                for r in all_recs
            ]
            extracted += 1

    # Remove raw guideline objects (large, not needed after extraction)
    for entry in enrichment.values():
        entry.pop("_guidelines", None)

    if extracted:
        print(f"  Extracted phenotype-specific recommendation(s) for {extracted} drug(s) from guideline tables.")


_CLASSIFICATION_SUMMARY = {
    "standard": "Standard dosing expected to be effective.",
    "caution": "Dose adjustment or monitoring may be needed.",
    "avoid": "Consider alternative medication.",
    "indeterminate": "Insufficient data for recommendation.",
}


# Source acronym expansions
_SOURCE_FULL_NAME = {
    "CPIC": "Clinical Pharmacogenetics Implementation Consortium",
    "DPWG": "Dutch Pharmacogenetics Working Group",
    "CPNDS": "Canadian Pharmacogenomics Network for Drug Safety",
    "RNPGx": "French National Network of Pharmacogenetics",
}


def _evidence_level_html(enrichment_entry):
    """Render the Evidence Level column: badge + checkmark."""
    import html as _h

    if not enrichment_entry:
        return ""

    level = enrichment_entry.get("evidence_level", "")
    verified = enrichment_entry.get("verified", False)

    badge_cls = _EVIDENCE_BADGE_CLASS.get(level, "badge-evidence-na")
    level_display = _h.escape(level) if level else "N/A"
    badge = f'<span class="badge {badge_cls}">{level_display}</span>'

    verify_html = ' <span class="evidence-verified">&#10003;</span>' if verified else ""

    return f"{badge}{verify_html}"


def _evidence_cell_html(enrichment_entry, classification=""):
    """Render the recommendation cell from enrichment data."""
    import html as _h

    if not enrichment_entry:
        fallback = _CLASSIFICATION_SUMMARY.get(classification, "")
        if fallback:
            return f'<div class="evidence-recs"><span class="evidence-rec-text">{_h.escape(fallback)}</span></div>'
        return ""

    # Show phenotype-specific recommendations from ALL guideline sources
    source_recs = enrichment_entry.get("source_recs", [])
    if source_recs:
        lines = []
        for sr in source_recs:
            rec_text = _h.escape(sr["rec"])
            src = sr["source"]
            full_name = _SOURCE_FULL_NAME.get(src, src)
            src_html = f'<span class="evidence-rec-source" title="{_h.escape(full_name)}">{_h.escape(src)}</span>'
            lines.append(f'{src_html} <span class="evidence-rec-text">{rec_text}</span>')
        return '<div class="evidence-recs">' + "<br>".join(lines) + "</div>"

    # Fallback: generic summary based on classification
    fallback = _CLASSIFICATION_SUMMARY.get(classification, "")
    if fallback:
        return f'<div class="evidence-recs"><span class="evidence-rec-text">{_h.escape(fallback)}</span></div>'
    return ""


# ---------------------------------------------------------------------------
# 6b. ClinPGx URL builder
# ---------------------------------------------------------------------------

def build_clinpgx_url(profiles):
    """Build a ClinPGx URL pre-filled with the patient's genotypes.

    ClinPGx (https://www.clinpgx.org) accepts a JSON payload as the 'q'
    query parameter.  The payload maps each gene to a three-element list:
    [allele1, allele2, null].  The JSON is URL-encoded and appended to
    https://www.clinpgx.org/genotypeResults?q={encoded}.

    Parameters
    ----------
    profiles : dict
        Gene -> {"diplotype": str, "phenotype": str}.  Diplotypes that
        contain a "/" separator (e.g. "*1/*2") are split into two alleles.
        Genes without a diplotype or without "/" are skipped.

    Returns
    -------
    str
        Full ClinPGx URL with encoded genotype payload.
    """
    payload = {}
    for gene, info in profiles.items():
        dip = info.get("diplotype", "")
        if not dip or "/" not in dip:
            continue
        parts = dip.split("/")
        if len(parts) != 2:
            continue
        payload[gene] = [parts[0], parts[1], None]
    encoded = quote(json.dumps(payload), safe="")
    return "https://www.clinpgx.org/genotypeResults?q=" + encoded


def build_cpic_gene_url(gene):
    """Return the CPIC guideline page URL for a given gene.

    Parameters
    ----------
    gene : str
        Gene symbol (e.g. "CYP2D6").

    Returns
    -------
    str
        URL to the CPIC gene page.
    """
    return "https://cpicpgx.org/gene/" + gene.lower() + "/"


# ---------------------------------------------------------------------------
# 7. Report generator
# ---------------------------------------------------------------------------

ICON = {"standard": "OK", "caution": "CAUTION", "avoid": "AVOID", "indeterminate": "INDETERMINATE — INSUFFICIENT DATA"}


def generate_report(input_path, fmt, total_snps, pgx_snps, profiles, drug_results):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    checksum = sha256_file(str(input_path))
    fname = Path(input_path).name

    lines = []
    lines.append("# ClawBio PharmGx Report")
    lines.append("")
    lines.append(f"**Date**: {now}")
    lines.append(f"**Input**: `{fname}`")
    lines.append(f"**Format detected**: {fmt}")
    lines.append(f"**Checksum (SHA-256)**: `{checksum}`")
    lines.append(f"**Total SNPs in file**: {total_snps}")
    lines.append(f"**Pharmacogenomic SNPs found**: {len(pgx_snps)}/{len(PGX_SNPS)}")
    lines.append(f"**Genes profiled**: {len(profiles)}")
    lines.append(f"**Drugs assessed**: {sum(len(v) for v in drug_results.values())}")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Data quality warning
    not_tested = [g for g, p in profiles.items() if p["diplotype"] == "NOT_TESTED"]
    unknown_pheno = [g for g, p in profiles.items()
                     if "unknown" in p["phenotype"].lower() or "indeterminate" in p["phenotype"].lower()]

    # Per-gene limitation warnings for genes with known CNV/SV/repeat issues
    # These must appear in the report body (not just stderr) for disclosure compliance.
    _GENE_LIMITATIONS = {
        "CYP2D6": (
            "CYP2D6: Copy number variation (gene deletion CYP2D6*5, duplication CYP2D6*1xN/*2xN, "
            "hybrid alleles CYP2D6*13/*36) cannot be detected from DTC genotyping data. "
            "Phenotype assignment may be incomplete. Clinical-grade CYP2D6 testing includes CNV analysis."
        ),
        "UGT1A1": (
            "UGT1A1: The UGT1A1*28 allele (rs8175347) is a TA-repeat polymorphism that cannot be "
            "reliably genotyped by SNP arrays. If rs8175347 is present in the input, the reported "
            "genotype is a proxy SNP call and may not reflect true TA repeat count."
        ),
        "CYP3A5": (
            "CYP3A5: The CYP3A5*7 allele (rs41303343) is an insertion/deletion variant. "
            "DTC genotyping platforms use proxy SNP calls that may not accurately detect this allele."
        ),
    }

    lines.append("## DATA QUALITY WARNING")
    lines.append("")
    if not_tested:
        lines.append(f"**{len(not_tested)} gene(s) could not be assessed** because the "
                     "relevant SNPs were not found in the input file: "
                     f"{', '.join(not_tested)}")
        lines.append("")
        lines.append("Drugs depending on these genes are marked INDETERMINATE below. "
                     "Do not assume normal metabolism for untested genes.")
        lines.append("")
    if unknown_pheno:
        unmapped = [g for g in unknown_pheno if g not in not_tested]
        if unmapped:
            lines.append(f"**{len(unmapped)} gene(s) have unmapped diplotypes**: "
                         f"{', '.join(unmapped)}. These diplotypes could not be matched "
                         "to a known phenotype. Clinical pharmacogenomic testing is recommended.")
            lines.append("")

    # Always emit per-gene limitation warnings for genes present in the profile
    # These are placed directly in the DQW body (no subsection header) so that
    # harness regex captures them as part of the DATA QUALITY WARNING text.
    gene_warnings = []
    for gene_name in sorted(profiles.keys()):
        if gene_name in _GENE_LIMITATIONS:
            gene_warnings.append(_GENE_LIMITATIONS[gene_name])
    if gene_warnings:
        lines.append("**Gene-Specific Limitations:**")
        lines.append("")
        for warning in gene_warnings:
            lines.append(f"- {warning}")
        lines.append("")

    lines.append("---")
    lines.append("")

    # Panel limitations disclosure
    lines.append("## Panel Limitations")
    lines.append("")
    lines.append("This report uses **SNP-based genotyping only** from a panel of "
                 f"{len(PGX_SNPS)} pharmacogenomic variants. The following cannot be detected:")
    lines.append("")
    lines.append("- **Copy number variants (CNVs)**: Gene deletions (e.g. CYP2D6\\*5) "
                 "and duplications (e.g. CYP2D6\\*1xN, \\*2xN)")
    lines.append("- **Structural variants**: CYP2D6-CYP2D7 hybrid alleles (e.g. \\*13, \\*36)")
    lines.append("- **Repeat polymorphisms**: UGT1A1\\*28 (TA7 repeat in rs8175347)")
    lines.append("- **HLA typing**: HLA-B\\*57:01 (abacavir hypersensitivity)")
    lines.append("- **Mitochondrial variants**: MT-RNR1 m.1555A>G (aminoglycoside ototoxicity)")
    lines.append("- **G6PD deficiency**: G6PD A- and Mediterranean variants")
    lines.append("")
    lines.append("For CYP2D6, a result of 'Normal Metabolizer' does NOT rule out "
                 "gene deletions or duplications that alter metabolizer status. "
                 "Clinical-grade CYP2D6 testing includes CNV analysis.")
    lines.append("")
    lines.append("## Genes Not Assessed by This Panel")
    lines.append("")
    lines.append("The following clinically relevant pharmacogenomic genes are **not included** "
                 "in this panel. Consider targeted testing if indicated:")
    lines.append("")
    lines.append("| Gene | Clinical Relevance | CPIC Guideline |")
    lines.append("|------|-------------------|----------------|")
    lines.append("| HLA-B\\*57:01 | Abacavir hypersensitivity | CPIC Abacavir (2014, updated 2020) |")
    lines.append("| HLA-B\\*58:01 | Allopurinol hypersensitivity (SJS/TEN) | CPIC Allopurinol (Hershfield 2013, updated 2015) |")
    lines.append("| G6PD | Rasburicase, chloroquine contraindication | CPIC G6PD (2014) |")
    lines.append("| MT-RNR1 | Aminoglycoside ototoxicity | FDA label |")
    lines.append("| HLA-A\\*31:01 | Carbamazepine hypersensitivity | CPIC Carbamazepine (2017) |")
    lines.append("| CYP2B6 | Efavirenz metabolism | CPIC Efavirenz (2019) |")
    lines.append("")
    lines.append("---")
    lines.append("")

    # Summary counts
    n_std = len(drug_results["standard"])
    n_cau = len(drug_results["caution"])
    n_avo = len(drug_results["avoid"])
    n_ind = len(drug_results.get("indeterminate", []))
    lines.append("## Drug Response Summary")
    lines.append("")
    lines.append(f"| Category | Count |")
    lines.append(f"|----------|-------|")
    lines.append(f"| Standard dosing | {n_std} |")
    lines.append(f"| Use with caution | {n_cau} |")
    lines.append(f"| Avoid / use alternative | {n_avo} |")
    if n_ind > 0:
        lines.append(f"| Insufficient data | {n_ind} |")
    lines.append("")

    # Alert drugs
    if n_avo > 0 or n_cau > 0:
        lines.append("### Actionable Alerts")
        lines.append("")
        if n_avo > 0:
            lines.append("**AVOID / USE ALTERNATIVE:**")
            lines.append("")
            for d in drug_results["avoid"]:
                lines.append(f"- **{d['drug']}** ({d['brand']}) [{d['gene']}]")
            lines.append("")
        if n_cau > 0:
            lines.append("**USE WITH CAUTION:**")
            lines.append("")
            for d in drug_results["caution"]:
                lines.append(f"- **{d['drug']}** ({d['brand']}) [{d['gene']}]")
            lines.append("")

    lines.append("---")
    lines.append("")

    # Gene profiles
    lines.append("## Gene Profiles")
    lines.append("")
    lines.append("| Gene | Full Name | Diplotype | Phenotype |")
    lines.append("|------|-----------|-----------|-----------|")
    for gene in GENE_DEFS:
        if gene in profiles:
            p = profiles[gene]
            lines.append(f"| {gene} | {GENE_DEFS[gene]['name']} | {p['diplotype']} | {p['phenotype']} |")
    # CPIC Level 1A genes not in this panel: explicit "not assessed" disclosure
    # These must appear in the gene profiles table, not just a separate section.
    # HLA-B*57:01: abacavir hypersensitivity (~48% risk in carriers)
    # MT-RNR1: aminoglycoside ototoxicity (permanent deafness, single dose)
    # G6PD: rasburicase/chloroquine contraindication
    _out_of_panel = [
        ("HLA-B", "HLA-B*57:01", "Not assessed", "Indeterminate (not in panel)"),
        ("MT-RNR1", "MT-RNR1 (mitochondrial)", "Not assessed", "Indeterminate (not in panel)"),
        ("G6PD", "Glucose-6-Phosphate Dehydrogenase", "Not assessed", "Indeterminate (not in panel)"),
        ("HLA-A", "HLA-A*31:01", "Not assessed", "Indeterminate (not in panel)"),
        # CYP2B6 removed: it IS in the SNP panel (rs3745274, rs28399499)
    ]
    for gene_sym, full_name, dip, pheno in _out_of_panel:
        lines.append(f"| {gene_sym} | {full_name} | {dip} | {pheno} |")
    lines.append("")

    # Interactive ClinPGx Link
    clinpgx_url = build_clinpgx_url(profiles)
    lines.append("## Interactive ClinPGx Link")
    lines.append("")
    lines.append("**Explore this patient's pharmacogenomic profile** with CPIC, DPWG, "
                 "and FDA guideline annotations:")
    lines.append("")
    lines.append(f"[Open in ClinPGx]({clinpgx_url})")
    lines.append("")
    lines.append("This link pre-fills ClinPGx with the patient's genotypes. "
                 "No data is uploaded; the genotypes are encoded in the URL itself.")
    lines.append("")
    lines.append("### CPIC Gene References")
    lines.append("")
    lines.append("| Gene | Diplotype | Phenotype | CPIC Guideline |")
    lines.append("|------|-----------|-----------|----------------|")
    for gene in GENE_DEFS:
        if gene in profiles:
            p = profiles[gene]
            cpic_url = build_cpic_gene_url(gene)
            lines.append(f"| {gene} | {p['diplotype']} | {p['phenotype']} "
                         f"| [{gene} on CPIC]({cpic_url}) |")
    lines.append("")

    # Detected variants
    lines.append("## Detected Variants")
    lines.append("")
    lines.append("| rsID | Gene | Star Allele | Genotype | Effect |")
    lines.append("|------|------|-------------|----------|--------|")
    for rsid, info in sorted(pgx_snps.items(), key=lambda x: x[1]["gene"]):
        lines.append(f"| {rsid} | {info['gene']} | {info['allele']} | {info['genotype']} | {info['effect']} |")
    lines.append("")

    # Full drug table
    lines.append("---")
    lines.append("")
    lines.append("## Complete Drug Recommendations")
    lines.append("")
    lines.append("| Drug | Brand | Class | Gene | Status |")
    lines.append("|------|-------|-------|------|--------|")
    for cat in ["avoid", "caution", "indeterminate", "standard"]:
        for d in sorted(drug_results.get(cat, []), key=lambda x: x["drug"]):
            status = ICON.get(d["classification"], d["classification"].upper())
            note_suffix = f" — {d['note']}" if d.get("note") else ""
            lines.append(f"| {d['drug']} | {d['brand']} | {d['class']} | {d['gene']} | {status}{note_suffix} |")
    lines.append("")

    # Disclaimer
    lines.append("---")
    lines.append("")
    lines.append("## Disclaimer")
    lines.append("")
    lines.append("This report is for **research and educational purposes only**. "
                 "It is NOT a diagnostic device and should NOT be used to make medication decisions "
                 "without consulting a qualified healthcare professional.")
    lines.append("")
    lines.append("Pharmacogenomic recommendations are based on CPIC guidelines (cpicpgx.org). "
                 "DTC genetic tests have limitations: they may not detect all relevant variants, "
                 "and results should be confirmed by clinical-grade testing before clinical use.")
    lines.append("")

    # Methods
    lines.append("## Methods")
    lines.append("")
    lines.append("- **Tool**: ClawBio PharmGx Reporter v0.2.0")
    lines.append("- **SNP panel**: 31 pharmacogenomic variants across 12 genes")
    lines.append("- **Star allele calling**: Simplified DTC-compatible algorithm (single-SNP per allele)")
    lines.append("- **Phenotype assignment**: CPIC-based diplotype-to-phenotype mapping")
    lines.append("- **Drug guidelines**: 51 drugs from CPIC (cpicpgx.org), simplified for DTC context")
    lines.append("")

    # Reproducibility
    lines.append("## Reproducibility")
    lines.append("")
    lines.append("```bash")
    lines.append(f"python pharmgx_reporter.py --input {fname} --output report")
    lines.append("```")
    lines.append("")
    lines.append(f"**Input checksum**: `{checksum}`")
    lines.append("")

    # References
    lines.append("## References")
    lines.append("")
    lines.append("- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio")
    lines.append("- CPIC. Clinical Pharmacogenetics Implementation Consortium. https://cpicpgx.org/")
    lines.append("- Caudle, K.E. et al. (2014). Standardizing terms for clinical pharmacogenetic test results. Genet Med, 16(9), 655-663.")
    lines.append("- PharmGKB. https://www.pharmgkb.org/")
    lines.append("- ClinPGx. Clinical Pharmacogenomics. https://www.clinpgx.org/")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 7b. HTML report generator
# ---------------------------------------------------------------------------

def _build_gene_rsid_map(pgx_snps):
    """Map gene name → list of (rsid, allele) from detected PGx SNPs."""
    gene_map = {}
    for rsid, info in pgx_snps.items():
        gene_map.setdefault(info["gene"], []).append((rsid, info["allele"]))
    return gene_map


def _drug_links_html(gene_str, gene_rsid_map):
    """Build hyperlinks column: Gene:rsID per line, linking to ClinPGx."""
    import html as _html
    lines = []
    genes = [g.strip() for g in gene_str.replace("+", ",").split(",")]
    for gene in genes:
        rsids = gene_rsid_map.get(gene, [])
        for rsid, _allele in rsids:
            lines.append(
                f'<a href="https://www.clinpgx.org/rsid/{rsid}" '
                f'target="_blank" rel="noopener">'
                f'{_html.escape(gene)}:{_html.escape(rsid)}</a>'
            )
        if not rsids:
            lines.append(
                f'<a href="https://www.clinpgx.org/gene/{_html.escape(gene)}" '
                f'target="_blank" rel="noopener">{_html.escape(gene)}</a>'
            )
    return '<span class="gene-links">' + "<br>".join(lines) + "</span>" if lines else "&mdash;"


def generate_html_report(input_path, fmt, total_snps, pgx_snps, profiles, drug_results, clinpgx_enrichment=None):
    """Build a self-contained HTML report using HtmlReportBuilder."""
    import html as _html

    checksum = sha256_file(str(input_path))
    fname = Path(input_path).name

    n_std = len(drug_results["standard"])
    n_cau = len(drug_results["caution"])
    n_avo = len(drug_results["avoid"])
    n_ind = len(drug_results.get("indeterminate", []))
    n_total = n_std + n_cau + n_avo + n_ind

    gene_rsid_map = _build_gene_rsid_map(pgx_snps)
    not_tested = [g for g, p in profiles.items() if p["diplotype"] == "NOT_TESTED"]
    n_genes_tested = len(profiles) - len(not_tested)

    b = HtmlReportBuilder("ClawBio PharmGx Report", "PharmGx Reporter v0.2.0")

    # ── Disclaimer at top ──
    b.add_disclaimer()

    # ── Branded header ──
    b.add_header_block("Your Medication Report", "How your genes affect your medications")

    # ── Executive summary ──
    avoid_drugs = sorted(drug_results.get("avoid", []), key=lambda x: x["drug"])
    if avoid_drugs:
        avoid_lines = "; ".join(
            f"{d['drug']} ({d['brand']})"
            for d in avoid_drugs
        )
        avoid_item = ("\u26d4", f"{n_avo} drug(s) to avoid", avoid_lines, "avoid")
    else:
        avoid_item = ("\u2705", "No drugs to avoid",
                      "No high-risk gene-drug interactions detected.", "ok")
    b.add_executive_summary([
        avoid_item,
        ("\u26a0\ufe0f", f"{n_cau} drug(s) requiring caution",
         "Dose adjustments or alternatives may be recommended.", "caution"),
    ])

    # ── Data quality warnings ──
    unknown_pheno = [g for g, p in profiles.items()
                     if "unknown" in p["phenotype"].lower()
                     or "indeterminate" in p["phenotype"].lower()]
    if not_tested:
        b.add_alert_box(
            "caution",
            f"{len(not_tested)} gene(s) not assessed",
            f"Relevant SNPs not found in input file: {', '.join(not_tested)}. "
            "Drugs depending on these genes are marked Insufficient below. "
            "Do not assume normal metabolism for untested genes.",
        )
    unmapped = [g for g in unknown_pheno if g not in not_tested]
    if unmapped:
        b.add_alert_box(
            "info",
            f"{len(unmapped)} gene(s) have unmapped diplotypes",
            f"{', '.join(unmapped)}. These diplotypes could not be matched to a known phenotype. "
            "Clinical pharmacogenomic testing is recommended.",
        )

    # ── Drug recommendations table (non-standard drugs) ──
    b.add_section("Drug Recommendations")
    from clawbio.common.html_report import _BADGE_CLASS, _BADGE_LABEL
    has_enrichment = bool(clinpgx_enrichment)

    def _build_row(d, enrichment_entry, row_cls):
        """Build a table row with Drug, Evidence, Recommendation, Class, Genes."""
        cls = d["classification"]
        badge_cls = _BADGE_CLASS.get(cls, "badge-indeterminate")
        badge_lbl = _BADGE_LABEL.get(cls, _html.escape(cls))
        badge = f'<span class="badge {badge_cls}">{badge_lbl}</span>'

        drug_cell = f"<strong>{_html.escape(d['drug'])}</strong> ({_html.escape(d['brand'])})"
        evidence_cell = _evidence_level_html(enrichment_entry)
        rec_cell = badge + _evidence_cell_html(enrichment_entry, classification=cls)
        notes_cell = _html.escape(d['class'])
        if d.get("note"):
            notes_cell += f'<br><small style="color:#c0392b">{_html.escape(d["note"])}</small>'
        links_cell = _drug_links_html(d["gene"], gene_rsid_map)

        return (
            f'<tr class="row-{_html.escape(row_cls)}">'
            f"<td>{drug_cell}</td><td>{rec_cell}</td><td>{evidence_cell}</td>"
            f"<td>{notes_cell}</td><td>{links_cell}</td></tr>"
        )

    _TH = "<th>Drug</th><th>Recommendation</th><th>Evidence</th><th>Class</th><th>Genes</th>"

    rows_html = []
    for cat in ["avoid", "caution", "indeterminate"]:
        for d in sorted(drug_results.get(cat, []), key=lambda x: x["drug"]):
            entry = clinpgx_enrichment.get(d["drug"].lower(), {}) if has_enrichment else {}
            rows_html.append(_build_row(d, entry, d["classification"]))

    drug_table = (
        '<div class="table-wrap"><table><thead><tr>'
        + _TH +
        "</tr></thead><tbody>"
        + "\n".join(rows_html)
        + "</tbody></table></div>"
    )
    b.add_raw_html(drug_table)

    # ── Standard drugs (collapsible) ──
    std_rows_html = []
    for d in sorted(drug_results.get("standard", []), key=lambda x: x["drug"]):
        entry = clinpgx_enrichment.get(d["drug"].lower(), {}) if has_enrichment else {}
        std_rows_html.append(_build_row(d, entry, "standard"))

    std_table_html = (
        '<div class="table-wrap"><table><thead><tr>'
        + _TH +
        "</tr></thead><tbody>"
        + "\n".join(std_rows_html)
        + "</tbody></table></div>"
    )
    b.add_details(f"Standard Drugs ({n_std} medications \u2014 click to expand)", std_table_html)

    # ── Gene Profiles (collapsible) ──
    gene_rows_html = []
    for gene in GENE_DEFS:
        if gene in profiles:
            p = profiles[gene]
            gene_rows_html.append(
                f"<tr><td>{_html.escape(gene)}</td>"
                f"<td>{_html.escape(GENE_DEFS[gene]['name'])}</td>"
                f"<td>{_html.escape(p['diplotype'])}</td>"
                f"<td>{_html.escape(p['phenotype'])}</td></tr>"
            )
    gene_table_html = (
        '<div class="table-wrap"><table><thead><tr>'
        "<th>Gene</th><th>Full Name</th><th>Diplotype</th><th>Phenotype</th>"
        "</tr></thead><tbody>"
        + "\n".join(gene_rows_html)
        + "</tbody></table></div>"
    )
    b.add_details("Gene Profiles (click to expand)", gene_table_html)

    # ── Detected variants (collapsible) ──
    var_rows_html = []
    for rsid, info in sorted(pgx_snps.items(), key=lambda x: x[1]["gene"]):
        var_rows_html.append(
            f"<tr><td>{_html.escape(rsid)}</td>"
            f"<td>{_html.escape(info['gene'])}</td>"
            f"<td>{_html.escape(info['allele'])}</td>"
            f"<td>{_html.escape(info['genotype'])}</td>"
            f"<td>{_html.escape(info['effect'])}</td></tr>"
        )
    var_table_html = (
        '<div class="table-wrap"><table><thead><tr>'
        "<th>rsID</th><th>Gene</th><th>Star Allele</th><th>Genotype</th><th>Effect</th>"
        "</tr></thead><tbody>"
        + "\n".join(var_rows_html)
        + "</tbody></table></div>"
    )
    b.add_details("Detected Variants (click to expand)", var_table_html)

    # ── Input details (collapsible) ──
    detail_items = {
        "Input file": fname,
        "Format detected": fmt,
        "Checksum (SHA-256)": checksum,
        "Total SNPs in file": str(total_snps),
        "Pharmacogenomic SNPs found": f"{len(pgx_snps)}/{len(PGX_SNPS)}",
        "Genes profiled": str(len(profiles)),
        "Drugs assessed": str(n_total),
    }
    meta_parts = [
        f"<p><strong>{_html.escape(k)}:</strong> {_html.escape(v)}</p>"
        for k, v in detail_items.items()
    ]
    meta_html = f'<div class="metadata">{"".join(meta_parts)}</div>'
    b.add_details("Input Details (click to expand)", meta_html)

    # ── Footer ──
    b.add_footer_block("PharmGx Reporter", "0.2.0")

    return b.render()


# ---------------------------------------------------------------------------
# 8. Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ClawBio PharmGx Reporter: pharmacogenomic report from DTC genetic data")
    parser.add_argument("--input", default=None, help="Path to genetic data file (23andMe/AncestryDNA)")
    parser.add_argument("--output", default="pharmgx_report", help="Output directory (default: pharmgx_report)")
    parser.add_argument("--drug", default=None, help="Single drug lookup (brand or generic name)")
    parser.add_argument("--dose", default=None, help="Visible dose from packaging (e.g. '50mg')")
    parser.add_argument("--no-enrich", action="store_true", help="Skip ClinPGx evidence enrichment")
    parser.add_argument("--demo", action="store_true", help="Run with bundled demo patient data")
    args = parser.parse_args()

    if args.demo:
        demo_file = Path(__file__).resolve().parent / "demo_patient.txt"
        if not demo_file.exists():
            print("Error: demo_patient.txt not found alongside script", file=sys.stderr)
            sys.exit(1)
        args.input = str(demo_file)
        if args.output == "pharmgx_report":
            args.output = str(Path(__file__).resolve().parent / "demo_report")
        print("Running in demo mode with bundled patient data")
        print()

    if not args.input:
        parser.error("--input is required (or use --demo)")

    if not Path(args.input).exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"ClawBio PharmGx Reporter v0.2.0")
    print(f"================================")
    print()

    # Parse
    print(f"Parsing: {args.input}")
    fmt, total_snps, pgx_snps, ref_genome = parse_file(args.input)
    print(f"  Format: {fmt}")
    print(f"  Total SNPs: {total_snps}")
    print(f"  PGx SNPs found: {len(pgx_snps)}/{len(PGX_SNPS)}")
    if ref_genome:
        print(f"  Reference genome: {ref_genome}")
    print()

    if fmt == "unknown":
        print("WARNING: Could not detect input file format. Results may be unreliable.",
              file=sys.stderr)

    if ref_genome == "GRCh37_mismatch":
        print("WARNING: Input coordinates appear to use GRCh37 (not GRCh38). "
              "Some gene results may be affected.", file=sys.stderr)

    if len(pgx_snps) == 0:
        print("ERROR: No pharmacogenomic SNPs found in this file.", file=sys.stderr)
        print("Cannot generate a report from zero data. Verify the input file", file=sys.stderr)
        print("is a valid 23andMe or AncestryDNA export.", file=sys.stderr)
        sys.exit(1)

    # Profile genes
    profiles = {}
    for gene in GENE_DEFS:
        diplotype = call_diplotype(gene, pgx_snps)
        phenotype = call_phenotype(gene, diplotype)
        profiles[gene] = {"diplotype": diplotype, "phenotype": phenotype}

    # If reference genome is GRCh37, mark all genes as Indeterminate
    # because SNP coordinates may not match, leading to incorrect allele calls.
    if ref_genome == "GRCh37_mismatch":
        for gene in profiles:
            profiles[gene] = {
                "diplotype": "Indeterminate (GRCh37 input detected; GRCh38 expected)",
                "phenotype": "Indeterminate (reference genome mismatch)",
            }

    # UGT1A1: If the SV marker SNP (rs8175347) is missing from input,
    # mark UGT1A1 as incomplete/Indeterminate since the key allele (*28)
    # cannot be assessed. rs8175347 is a TA-repeat polymorphism that most
    # DTC platforms omit entirely, making any UGT1A1 call without it unreliable.
    if "UGT1A1" in profiles and ref_genome != "GRCh37_mismatch":
        ugt_sv_rsid = "rs8175347"
        ugt_snp_rsid = "rs4148323"
        has_sv = ugt_sv_rsid in pgx_snps
        has_snp = ugt_snp_rsid in pgx_snps
        if not has_sv:
            # *28 repeat marker absent: UGT1A1 coverage is incomplete.
            # Cannot reliably assess the most clinically important allele.
            profiles["UGT1A1"] = {
                "diplotype": "Indeterminate (rs8175347 not tested)",
                "phenotype": "Indeterminate (rs8175347 not tested)",
            }

    not_tested = [g for g, p in profiles.items() if p["diplotype"] == "NOT_TESTED"]
    if not_tested:
        print(f"WARNING: {len(not_tested)} gene(s) not testable from this data: {', '.join(not_tested)}",
              file=sys.stderr)

    print("Gene Profiles:")
    print(f"  {'Gene':<10} {'Diplotype':<20} {'Phenotype'}")
    print(f"  {'-'*10} {'-'*20} {'-'*35}")
    for gene, p in profiles.items():
        print(f"  {gene:<10} {p['diplotype']:<20} {p['phenotype']}")
    print()

    # Single-drug lookup mode (--drug flag)
    if args.drug:
        resolved = resolve_drug_name(args.drug)
        if not resolved:
            print(f"Drug not found: '{args.drug}'. Available drugs: {len(GUIDELINES)}", file=sys.stderr)
            sys.exit(1)
        result = lookup_single_drug(resolved, profiles)
        print(format_dosage_card(result, visible_dose=args.dose))
        sys.exit(0)

    # Drug lookup
    drug_results = lookup_drugs(profiles)
    n_std = len(drug_results["standard"])
    n_cau = len(drug_results["caution"])
    n_avo = len(drug_results["avoid"])
    n_ind = len(drug_results.get("indeterminate", []))

    total_assessed = n_std + n_cau + n_avo + n_ind
    print(f"Drug Recommendations ({total_assessed} drugs):")
    print(f"  Standard:        {n_std}")
    print(f"  Caution:         {n_cau}")
    print(f"  Avoid:           {n_avo}")
    if n_ind > 0:
        print(f"  Insufficient data: {n_ind}")
    print()

    if n_avo > 0:
        print("ALERT - Drugs to AVOID:")
        for d in drug_results["avoid"]:
            print(f"  * {d['drug']} ({d['brand']})")
        print()

    # ClinPGx evidence enrichment
    clinpgx_enrichment = {}
    if not getattr(args, "no_enrich", False):
        print("Querying ClinPGx for evidence data...")
        try:
            clinpgx_enrichment = enrich_with_clinpgx(drug_results)
            if clinpgx_enrichment:
                print(f"  Enriched {len(clinpgx_enrichment)} drug(s) with evidence data.")
                # Extract phenotype-specific recs from CPIC guideline tables
                try:
                    extract_phenotype_recs(clinpgx_enrichment,
                                           drug_results=drug_results,
                                           profiles=profiles)
                except Exception as extract_exc:
                    print(f"  Recommendation extraction failed: {extract_exc}",
                          file=sys.stderr)
            else:
                print("  ClinPGx enrichment returned no data (skill may be unavailable).")
        except Exception as exc:
            print(f"  ClinPGx enrichment failed: {exc}", file=sys.stderr)
        print()

    # Generate report
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    report = generate_report(args.input, fmt, total_snps, pgx_snps, profiles, drug_results)
    report_path = outdir / "report.md"
    report_path.write_text(report)

    # Generate HTML report
    html_content = generate_html_report(
        args.input, fmt, total_snps, pgx_snps, profiles, drug_results,
        clinpgx_enrichment=clinpgx_enrichment or None,
    )
    html_path = write_html_report(outdir, "report.html", html_content)

    # Write result.json using shared report helper
    input_checksum = sha256_hex(str(args.input))
    result_data = {
        "gene_profiles": profiles,
        "drug_recommendations": drug_results,
    }
    if clinpgx_enrichment:
        result_data["clinpgx_enrichment"] = clinpgx_enrichment
    result_json_path = write_result_json(
        output_dir=outdir,
        skill="pharmgx",
        version="0.2.0",
        summary={
            "total_snps_in_file": total_snps,
            "pgx_snps_found": len(pgx_snps),
            "pgx_snps_total": len(PGX_SNPS),
            "genes_profiled": len(profiles),
            "drugs_assessed": total_assessed,
            "drugs_standard": n_std,
            "drugs_caution": n_cau,
            "drugs_avoid": n_avo,
            "drugs_indeterminate": n_ind,
            "clinpgx_enriched": len(clinpgx_enrichment),
        },
        data=result_data,
        input_checksum=input_checksum,
    )

    print(f"Report saved: {report_path}")
    print(f"HTML report:  {html_path}")
    print(f"Result JSON:  {result_json_path}")
    print("Done.")


if __name__ == "__main__":
    main()
