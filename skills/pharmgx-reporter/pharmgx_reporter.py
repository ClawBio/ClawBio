#!/usr/bin/env python3
"""
ClawBio PharmGx Reporter
Pharmacogenomic report generator from DTC genetic data (23andMe/AncestryDNA).

Analyses 31 pharmacogenomic SNPs across 12 genes, calls star alleles and
metabolizer phenotypes, and looks up CPIC drug recommendations for 51 medications.

Usage:
    python pharmgx_reporter.py --input patient_data.txt --output report_dir
"""

import argparse
import gzip
import hashlib
import os
import re
import sys
from datetime import datetime, timezone
from pathlib import Path

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
    "rs8175347":  {"gene": "UGT1A1", "allele": "*28", "effect": "decreased_function"},
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
            "Ultrarapid Metabolizer":  ["*17/*17", "*1/*17"],
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
            "rs5030655":  {"allele": "*6",  "alt": "DEL", "effect": "no_function"},
            "rs16947":    {"allele": "*2",  "alt": "A", "effect": "normal_function"},
            "rs1065852":  {"allele": "*10", "alt": "T", "effect": "decreased_function"},
            "rs28371725": {"allele": "*41", "alt": "T", "effect": "decreased_function"},
        },
        "phenotypes": {
            "Normal Metabolizer":       ["*1/*1", "*1/*2", "*2/*2"],
            "Intermediate Metabolizer": ["*1/*4", "*1/*10", "*1/*41", "*10/*10", "*4/*10", "*10/*41", "*41/*41"],
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
            "Intermediate Function": ["TC", "CT"],
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
            "Poor Metabolizer":         ["*2/*2", "*2/*3A", "*3A/*3A", "*3B/*3B", "*3C/*3C"],
        },
    },
    "UGT1A1": {
        "name": "UDP-Glucuronosyltransferase 1A1",
        "function": "Irinotecan and bilirubin metabolism",
        "ref": "*1",
        "variants": {
            "rs8175347": {"allele": "*28", "alt": "TA7", "effect": "decreased_function"},
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
            "rs41303343": {"allele": "*7", "alt": "INS", "effect": "no_function"},
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
            "ultrarapid_metabolizer": ("standard", "Use recommended dose."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Consider alternative antiplatelet therapy (prasugrel, ticagrelor)."),
            "poor_metabolizer":       ("avoid", "Use alternative antiplatelet therapy (prasugrel or ticagrelor)."),
        },
    },
    "Omeprazole": {
        "brand": "Prilosec", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "For H. pylori: increase dose 50-100%. For GERD: standard dose."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "For chronic therapy >12 wk: consider 50% dose reduction."),
        },
    },
    "Pantoprazole": {
        "brand": "Protonix", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "For H. pylori: increase dose 50-100%."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "For chronic therapy >12 wk: consider 50% dose reduction."),
        },
    },
    "Lansoprazole": {
        "brand": "Prevacid", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "For H. pylori: increase dose 50-100%."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "For chronic therapy: consider 50% dose reduction."),
        },
    },
    "Esomeprazole": {
        "brand": "Nexium", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "For H. pylori: increase dose 50-100%."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "For chronic therapy: consider dose reduction."),
        },
    },
    "Dexlansoprazole": {
        "brand": "Dexilant", "class": "Proton Pump Inhibitor", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "For H. pylori: increase dose 50-100%."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "For chronic therapy: consider dose reduction."),
        },
    },
    "Citalopram": {
        "brand": "Celexa", "class": "SSRI Antidepressant", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "Select alternative SSRI not dependent on CYP2C19."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "Consider 50% reduction. Max 20 mg/day."),
        },
    },
    "Escitalopram": {
        "brand": "Lexapro", "class": "SSRI Antidepressant", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "Select alternative SSRI or titrate to max dose."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "Consider 50% dose reduction."),
        },
    },
    "Sertraline": {
        "brand": "Zoloft", "class": "SSRI Antidepressant", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("standard", "Use recommended starting dose."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("standard", "Use recommended starting dose."),
        },
    },
    "Voriconazole": {
        "brand": "Vfend", "class": "Antifungal", "gene": "CYP2C19",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "Use alternative antifungal; voriconazole may be ineffective."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("standard", "Use recommended dose."),
            "poor_metabolizer":       ("caution", "Increased exposure; consider dose reduction or TDM."),
        },
    },
    # --- CYP2D6 drugs ---
    "Codeine": {
        "brand": "Tylenol w/ Codeine", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid codeine. Risk of morphine toxicity."),
            "normal_metabolizer":     ("standard", "Use label-recommended dosing."),
            "intermediate_metabolizer": ("caution", "Use with caution; may have reduced analgesia."),
            "poor_metabolizer":       ("avoid", "Avoid codeine. Select alternative analgesic."),
        },
    },
    "Tramadol": {
        "brand": "Ultram", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid tramadol. Increased toxicity risk."),
            "normal_metabolizer":     ("standard", "Use label-recommended dosing."),
            "intermediate_metabolizer": ("caution", "Use with caution; possible reduced analgesia."),
            "poor_metabolizer":       ("avoid", "Avoid tramadol. Select alternative analgesic."),
        },
    },
    "Hydrocodone": {
        "brand": "Vicodin", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "Use with caution; increased active metabolite."),
            "normal_metabolizer":     ("standard", "Use label-recommended dosing."),
            "intermediate_metabolizer": ("caution", "May have reduced analgesia; monitor response."),
            "poor_metabolizer":       ("caution", "Reduced analgesia likely; consider alternative."),
        },
    },
    "Oxycodone": {
        "brand": "OxyContin", "class": "Opioid Analgesic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "Use with caution; increased active metabolite."),
            "normal_metabolizer":     ("standard", "Use label-recommended dosing."),
            "intermediate_metabolizer": ("standard", "Use label-recommended dosing."),
            "poor_metabolizer":       ("caution", "Reduced analgesia possible; consider alternative."),
        },
    },
    "Tamoxifen": {
        "brand": "Nolvadex", "class": "SERM (Oncology)", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("standard", "Use recommended dose."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Reduced efficacy possible. Avoid strong CYP2D6 inhibitors."),
            "poor_metabolizer":       ("avoid", "Consider aromatase inhibitor or higher tamoxifen dose with TDM."),
        },
    },
    "Amitriptyline": {
        "brand": "Elavil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid TCA; likely ineffective."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("caution", "Consider 25% dose reduction."),
            "poor_metabolizer":       ("avoid", "Avoid TCA. If necessary, reduce dose 50%."),
        },
    },
    "Nortriptyline": {
        "brand": "Pamelor", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid TCA; reduced efficacy."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("caution", "Consider 25% dose reduction."),
            "poor_metabolizer":       ("avoid", "Avoid TCA or reduce dose 50%."),
        },
    },
    "Desipramine": {
        "brand": "Norpramin", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid TCA; likely ineffective."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("caution", "Consider 25% dose reduction."),
            "poor_metabolizer":       ("avoid", "Avoid TCA or reduce dose 50%."),
        },
    },
    "Imipramine": {
        "brand": "Tofranil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid TCA; likely ineffective."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("caution", "Consider 25% dose reduction."),
            "poor_metabolizer":       ("avoid", "Avoid TCA or reduce dose 50%."),
        },
    },
    "Doxepin": {
        "brand": "Sinequan", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid TCA; likely ineffective."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("caution", "Consider 25% dose reduction."),
            "poor_metabolizer":       ("avoid", "Avoid TCA or reduce dose 50%."),
        },
    },
    "Trimipramine": {
        "brand": "Surmontil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid TCA; likely ineffective."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("caution", "Consider 25% dose reduction."),
            "poor_metabolizer":       ("avoid", "Avoid TCA or reduce dose 50%."),
        },
    },
    "Clomipramine": {
        "brand": "Anafranil", "class": "Tricyclic Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("avoid", "Avoid TCA; likely ineffective."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("caution", "Consider 25% dose reduction."),
            "poor_metabolizer":       ("avoid", "Avoid TCA or reduce dose 50%."),
        },
    },
    "Paroxetine": {
        "brand": "Paxil", "class": "SSRI Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "Select alternative or titrate to response."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "Consider 50% dose reduction or alternative SSRI."),
        },
    },
    "Fluoxetine": {
        "brand": "Prozac", "class": "SSRI Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "Select alternative or titrate to response."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "Consider 50% dose reduction or alternative."),
        },
    },
    "Venlafaxine": {
        "brand": "Effexor", "class": "SNRI Antidepressant", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("standard", "Use recommended dose."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("standard", "Use recommended dose."),
            "poor_metabolizer":       ("caution", "Consider 50% dose reduction or switch to desvenlafaxine."),
        },
    },
    "Metoprolol": {
        "brand": "Lopressor", "class": "Beta-Blocker", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "May need higher dose or alternative beta-blocker."),
            "normal_metabolizer":     ("standard", "Use recommended starting dose."),
            "intermediate_metabolizer": ("standard", "Use recommended starting dose."),
            "poor_metabolizer":       ("caution", "Consider 50% dose reduction or alternative beta-blocker."),
        },
    },
    "Ondansetron": {
        "brand": "Zofran", "class": "Antiemetic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "May have reduced antiemetic effect."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("standard", "Use recommended dose."),
            "poor_metabolizer":       ("caution", "May have reduced antiemetic effect; consider alternative."),
        },
    },
    "Risperidone": {
        "brand": "Risperdal", "class": "Antipsychotic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "May need higher dose."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Consider dose reduction."),
            "poor_metabolizer":       ("caution", "Consider 50% dose reduction."),
        },
    },
    "Aripiprazole": {
        "brand": "Abilify", "class": "Antipsychotic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("standard", "Use recommended dose."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("standard", "Use recommended dose."),
            "poor_metabolizer":       ("caution", "Reduce dose to 75% of usual."),
        },
    },
    "Haloperidol": {
        "brand": "Haldol", "class": "Antipsychotic", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("caution", "May need higher dose."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Consider dose reduction."),
            "poor_metabolizer":       ("caution", "Reduce dose; monitor for side effects."),
        },
    },
    "Atomoxetine": {
        "brand": "Strattera", "class": "ADHD Medication", "gene": "CYP2D6",
        "recs": {
            "ultrarapid_metabolizer": ("standard", "Use recommended dose."),
            "normal_metabolizer":     ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("standard", "Use recommended dose."),
            "poor_metabolizer":       ("caution", "Start at lower dose; 2-fold higher exposure expected."),
        },
    },
    # --- CYP2C9 drugs ---
    "Phenytoin": {
        "brand": "Dilantin", "class": "Antiepileptic", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Reduce dose 25%. Monitor levels closely."),
            "poor_metabolizer":         ("avoid", "Reduce dose 50%. Consider alternative anticonvulsant."),
        },
    },
    "Celecoxib": {
        "brand": "Celebrex", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Start at lowest dose."),
            "poor_metabolizer":         ("avoid", "Use lowest dose or avoid; consider alternative NSAID."),
        },
    },
    "Flurbiprofen": {
        "brand": "Ansaid", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Start at lowest dose."),
            "poor_metabolizer":         ("avoid", "Use lowest dose or consider alternative."),
        },
    },
    "Piroxicam": {
        "brand": "Feldene", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Start at lowest dose; monitor."),
            "poor_metabolizer":         ("avoid", "Avoid or use lowest dose."),
        },
    },
    "Meloxicam": {
        "brand": "Mobic", "class": "NSAID", "gene": "CYP2C9",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Start at lowest dose."),
            "poor_metabolizer":         ("caution", "Use lowest dose; monitor GI and renal effects."),
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
            "normal_function":       ("standard", "Use desired starting dose."),
            "intermediate_function": ("caution", "Lower dose or alternative statin. 4.5x myopathy risk."),
            "poor_function":         ("avoid", "Use alternative statin (pravastatin, rosuvastatin) or max 20 mg/day."),
        },
    },
    "Atorvastatin": {
        "brand": "Lipitor", "class": "Statin", "gene": "SLCO1B1",
        "recs": {
            "normal_function":       ("standard", "Use desired starting dose."),
            "intermediate_function": ("caution", "Consider CK surveillance."),
            "poor_function":         ("caution", "Lower starting dose or alternative statin."),
        },
    },
    "Rosuvastatin": {
        "brand": "Crestor", "class": "Statin", "gene": "SLCO1B1",
        "recs": {
            "normal_function":       ("standard", "Use desired starting dose."),
            "intermediate_function": ("standard", "Use desired starting dose."),
            "poor_function":         ("standard", "Preferred alternative to simvastatin."),
        },
    },
    "Pravastatin": {
        "brand": "Pravachol", "class": "Statin", "gene": "SLCO1B1",
        "recs": {
            "normal_function":       ("standard", "Use desired starting dose."),
            "intermediate_function": ("standard", "Use desired starting dose."),
            "poor_function":         ("standard", "Preferred alternative to simvastatin."),
        },
    },
    # --- DPYD drugs ---
    "Fluorouracil": {
        "brand": "5-FU", "class": "Antineoplastic", "gene": "DPYD",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Reduce dose 50%. Monitor toxicity closely."),
            "poor_metabolizer":         ("avoid", "Avoid fluorouracil. Select alternative agent."),
        },
    },
    "Capecitabine": {
        "brand": "Xeloda", "class": "Antineoplastic", "gene": "DPYD",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Reduce dose 50%. Monitor toxicity."),
            "poor_metabolizer":         ("avoid", "Avoid capecitabine. Select alternative."),
        },
    },
    # --- TPMT / NUDT15 drugs ---
    "Azathioprine": {
        "brand": "Imuran", "class": "Immunosuppressant", "gene": "TPMT",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Reduce dose 30-70%. Monitor blood counts weekly."),
            "poor_metabolizer":         ("avoid", "Avoid or reduce to 10% dose. Consider alternative."),
        },
    },
    "Mercaptopurine": {
        "brand": "Purinethol", "class": "Immunosuppressant", "gene": "TPMT",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Reduce dose 30-70%. Monitor blood counts."),
            "poor_metabolizer":         ("avoid", "Avoid or reduce to 10% dose."),
        },
    },
    "Thioguanine": {
        "brand": "Tabloid", "class": "Immunosuppressant", "gene": "TPMT",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Reduce dose 30-70%."),
            "poor_metabolizer":         ("avoid", "Avoid or drastically reduce dose."),
        },
    },
    # --- UGT1A1 drug ---
    "Irinotecan": {
        "brand": "Camptosar", "class": "Antineoplastic", "gene": "UGT1A1",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Consider dose reduction if prior toxicity."),
            "poor_metabolizer":         ("avoid", "Reduce initial dose by at least one level."),
        },
    },
    "Atazanavir": {
        "brand": "Reyataz", "class": "Antiretroviral", "gene": "UGT1A1",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("standard", "Use recommended dose. May develop jaundice."),
            "poor_metabolizer":         ("caution", "Higher risk of jaundice. Monitor bilirubin."),
        },
    },
    # --- CYP3A5 drug ---
    "Tacrolimus": {
        "brand": "Prograf", "class": "Immunosuppressant", "gene": "CYP3A5",
        "recs": {
            "extensive_metabolizer":    ("caution", "Increase dose 1.5-2x. Titrate to target trough."),
            "intermediate_metabolizer": ("caution", "Increase dose 1.5x. Titrate to target trough."),
            "poor_metabolizer":         ("standard", "Use recommended dose (most patients)."),
        },
    },
    # --- CYP2B6 drug ---
    "Efavirenz": {
        "brand": "Sustiva", "class": "Antiretroviral", "gene": "CYP2B6",
        "recs": {
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("caution", "Consider dose reduction to 400 mg if toxicity occurs."),
            "poor_metabolizer":         ("caution", "Reduce dose to 400 mg or consider alternative."),
        },
    },
    # --- CYP1A2 drugs ---
    "Clozapine": {
        "brand": "Clozaril", "class": "Antipsychotic", "gene": "CYP1A2",
        "recs": {
            "ultrarapid_metabolizer":   ("caution", "May need higher dose; monitor levels."),
            "normal_metabolizer":       ("standard", "Use recommended dose."),
            "intermediate_metabolizer": ("standard", "Use recommended dose."),
            "poor_metabolizer":         ("caution", "Consider dose reduction; monitor for toxicity."),
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
        classification, rec = get_warfarin_rec(profiles)
        cyp2c9 = profiles.get("CYP2C9", {})
        vkorc1 = profiles.get("VKORC1", {})
        return {
            "drug": drug_name, "brand": info["brand"], "class": info["class"],
            "gene": "CYP2C9 + VKORC1",
            "diplotype": f"CYP2C9 {cyp2c9.get('diplotype', '?')} / VKORC1 {vkorc1.get('diplotype', '?')}",
            "phenotype": f"CYP2C9 {cyp2c9.get('phenotype', '?')} / VKORC1 {vkorc1.get('phenotype', '?')}",
            "classification": classification, "recommendation": rec,
        }

    gene = info["gene"]
    if gene not in profiles:
        return {
            "drug": drug_name, "brand": info["brand"], "class": info["class"],
            "gene": gene, "diplotype": "NOT_TESTED", "phenotype": "Indeterminate",
            "classification": "indeterminate",
            "recommendation": "Gene not profiled. No recommendation available.",
        }

    prof = profiles[gene]
    pheno_key = phenotype_to_key(prof["phenotype"])
    recs = info.get("recs", {})

    if pheno_key in recs:
        classification, rec = recs[pheno_key]
    elif pheno_key == "indeterminate":
        classification = "indeterminate"
        rec = f"Phenotype indeterminate ({prof['phenotype']}). Cannot assess."
    else:
        classification = "indeterminate"
        rec = f"Phenotype '{prof['phenotype']}' not in guidelines."

    return {
        "drug": drug_name, "brand": info["brand"], "class": info["class"],
        "gene": gene, "diplotype": prof["diplotype"], "phenotype": prof["phenotype"],
        "classification": classification, "recommendation": rec,
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
    rec_text = result["recommendation"]
    if visible_dose:
        cl = result["classification"]
        if cl == "standard":
            rec_text = f"Your genotype supports {result['drug']} {visible_dose} as prescribed."
        elif cl == "caution":
            rec_text = f"{visible_dose} may need adjustment. {result['recommendation']}"
        elif cl == "avoid":
            rec_text = f"Your genotype contraindicates {result['drug']} {visible_dose}. {result['recommendation']}"

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
# 4. File parser
# ---------------------------------------------------------------------------

def detect_format(lines):
    for line in lines[:20]:
        if line.startswith("# rsid"):
            return "23andme"
        if "RSID" in line and "CHROMOSOME" in line:
            return "ancestrydna"
        if "rsid" in line and "chromosome" in line:
            return "generic"
    for line in lines[:20]:
        if line.startswith("#") or line.strip() == "":
            continue
        parts = line.split("\t")
        if len(parts) >= 4 and parts[0].startswith("rs"):
            return "23andme"
        parts = line.split(",")
        if len(parts) >= 4 and parts[0].startswith("rs"):
            return "ancestrydna"
    return "unknown"


def parse_file(path):
    if str(path).endswith(".gz"):
        with gzip.open(path, "rt", errors="replace") as fh:
            content = fh.read()
    else:
        content = Path(path).read_text()
    lines = content.split("\n")
    fmt = detect_format(lines)

    snps = {}
    for line in lines:
        if line.startswith("#") or line.strip() == "":
            continue
        if "rsid" in line.lower() and "chromosome" in line.lower():
            continue
        parts = line.split("\t") if "\t" in line else line.split(",")
        if len(parts) >= 4:
            rsid = parts[0].strip()
            if not rsid.startswith("rs"):
                continue
            if len(parts) == 5:
                genotype = parts[3].strip() + parts[4].strip()
            else:
                genotype = parts[3].strip()
            if genotype and genotype not in ("--", "00"):
                snps[rsid] = genotype.upper()

    pgx = {}
    for rsid, info in PGX_SNPS.items():
        if rsid in snps:
            pgx[rsid] = {"genotype": snps[rsid], **info}

    return fmt, len(snps), pgx


# ---------------------------------------------------------------------------
# 5. Star allele caller
# ---------------------------------------------------------------------------

def call_diplotype(gene, pgx_snps):
    gdef = GENE_DEFS[gene]

    if gdef.get("type") == "genotype":
        rsid = gdef["rsid"]
        if rsid in pgx_snps:
            return pgx_snps[rsid]["genotype"]
        return "NOT_TESTED"

    # Count how many of this gene's SNPs were actually present in the file
    gene_rsids = list(gdef["variants"].keys())
    tested = [r for r in gene_rsids if r in pgx_snps]

    if not tested:
        return "NOT_TESTED"

    detected = []
    for rsid, vdef in gdef["variants"].items():
        if rsid in pgx_snps:
            gt = pgx_snps[rsid]["genotype"]
            alt = vdef["alt"].upper()
            if alt in ("DEL", "INS", "TA7"):
                print(f"  WARNING: {gene} {rsid} has structural variant "
                      f"alt={alt}, cannot interpret from DTC data",
                      file=sys.stderr)
                continue
            alt_count = gt.count(alt)
            if alt_count > 0:
                detected.append({"rsid": rsid, "allele": vdef["allele"],
                                 "copies": alt_count, "effect": vdef["effect"]})

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
    if diplotype == "NOT_TESTED":
        return "Indeterminate (not genotyped)"

    gdef = GENE_DEFS[gene]
    norm = diplotype.upper()

    # Strip partial-coverage annotations for matching (e.g. "*1/*1 (2/4 SNPs tested)")
    match_str = norm.split("(")[0].strip()

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
        "Ultrarapid Metabolizer": "ultrarapid_metabolizer",
        "Normal Warfarin Sensitivity": "normal_warfarin_sensitivity",
        "Intermediate Warfarin Sensitivity": "intermediate_warfarin_sensitivity",
        "High Warfarin Sensitivity": "high_warfarin_sensitivity",
        "Normal Function": "normal_function",
        "Intermediate Function": "intermediate_function",
        "Poor Function": "poor_function",
        "CYP3A5 Expressor": "extensive_metabolizer",
        "Intermediate Expressor": "intermediate_metabolizer",
        "CYP3A5 Non-expressor": "poor_metabolizer",
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
        return "standard", "Use warfarin dosing algorithm. Standard dose range expected."
    elif "poor" in cyp2c9.lower() or "high" in vkorc1.lower():
        return "avoid", "Significantly reduce dose (50-80% reduction). Consider DOAC alternative."
    else:
        return "caution", "Reduce initial dose. Use genotype-guided dosing algorithm."


def lookup_drugs(profiles):
    results = {"standard": [], "caution": [], "avoid": [], "indeterminate": []}

    for drug_name, drug in GUIDELINES.items():
        if drug.get("special") == "warfarin":
            classification, rec = get_warfarin_rec(profiles)
            results.setdefault(classification, []).append({
                "drug": drug_name, "brand": drug["brand"],
                "class": drug["class"], "gene": "CYP2C9+VKORC1",
                "recommendation": rec, "classification": classification,
            })
            continue

        gene = drug["gene"]
        if gene not in profiles:
            results["indeterminate"].append({
                "drug": drug_name, "brand": drug["brand"],
                "class": drug["class"], "gene": gene,
                "recommendation": "Gene not profiled. No recommendation available.",
                "classification": "indeterminate",
            })
            continue

        pheno_key = phenotype_to_key(profiles[gene]["phenotype"])

        if pheno_key == "indeterminate":
            results["indeterminate"].append({
                "drug": drug_name, "brand": drug["brand"],
                "class": drug["class"], "gene": gene,
                "recommendation": f"Gene phenotype indeterminate ({profiles[gene]['phenotype']}). Cannot assess.",
                "classification": "indeterminate",
            })
            continue

        recs = drug.get("recs", {})

        if pheno_key in recs:
            classification, rec = recs[pheno_key]
        else:
            classification = "indeterminate"
            rec = f"Phenotype '{profiles[gene]['phenotype']}' not covered by available guidelines. Consult clinical pharmacogenomics."

        results.setdefault(classification, []).append({
            "drug": drug_name, "brand": drug["brand"],
            "class": drug["class"], "gene": gene,
            "recommendation": rec, "classification": classification,
        })

    return results


# ---------------------------------------------------------------------------
# 7. Report generator
# ---------------------------------------------------------------------------

ICON = {"standard": "OK", "caution": "CAUTION", "avoid": "AVOID", "indeterminate": "INSUFFICIENT DATA"}


def generate_report(input_path, fmt, total_snps, pgx_snps, profiles, drug_results):
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    checksum = hashlib.sha256(Path(input_path).read_bytes()).hexdigest()
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
    if not_tested or unknown_pheno:
        lines.append("## DATA QUALITY WARNING")
        lines.append("")
        if not_tested:
            lines.append(f"**{len(not_tested)} gene(s) could not be assessed** because the "
                         "relevant SNPs were not found in the input file: "
                         f"{', '.join(not_tested)}")
            lines.append("")
            lines.append("Drugs depending on these genes are marked INSUFFICIENT DATA below. "
                         "Do not assume normal metabolism for untested genes.")
            lines.append("")
        if unknown_pheno:
            unmapped = [g for g in unknown_pheno if g not in not_tested]
            if unmapped:
                lines.append(f"**{len(unmapped)} gene(s) have unmapped diplotypes**: "
                             f"{', '.join(unmapped)}. These diplotypes could not be matched "
                             "to a known phenotype. Clinical pharmacogenomic testing is recommended.")
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
                lines.append(f"- **{d['drug']}** ({d['brand']}) [{d['gene']}]: {d['recommendation']}")
            lines.append("")
        if n_cau > 0:
            lines.append("**USE WITH CAUTION:**")
            lines.append("")
            for d in drug_results["caution"]:
                lines.append(f"- **{d['drug']}** ({d['brand']}) [{d['gene']}]: {d['recommendation']}")
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
    lines.append("| Drug | Brand | Class | Gene | Status | Recommendation |")
    lines.append("|------|-------|-------|------|--------|----------------|")
    for cat in ["avoid", "caution", "indeterminate", "standard"]:
        for d in sorted(drug_results.get(cat, []), key=lambda x: x["drug"]):
            status = ICON.get(d["classification"], d["classification"].upper())
            lines.append(f"| {d['drug']} | {d['brand']} | {d['class']} | {d['gene']} | {status} | {d['recommendation']} |")
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
    lines.append("- **Tool**: ClawBio PharmGx Reporter v0.1.0")
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
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# 8. Main
# ---------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="ClawBio PharmGx Reporter: pharmacogenomic report from DTC genetic data")
    parser.add_argument("--input", required=True, help="Path to genetic data file (23andMe/AncestryDNA)")
    parser.add_argument("--output", default="pharmgx_report", help="Output directory (default: pharmgx_report)")
    parser.add_argument("--drug", default=None, help="Single drug lookup (brand or generic name)")
    parser.add_argument("--dose", default=None, help="Visible dose from packaging (e.g. '50mg')")
    parser.add_argument("--html", action="store_true", help="Also generate a rich-text HTML report")
    args = parser.parse_args()

    if not Path(args.input).exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    print(f"ClawBio PharmGx Reporter v0.1.0")
    print(f"================================")
    print()

    # Parse
    print(f"Parsing: {args.input}")
    fmt, total_snps, pgx_snps = parse_file(args.input)
    print(f"  Format: {fmt}")
    print(f"  Total SNPs: {total_snps}")
    print(f"  PGx SNPs found: {len(pgx_snps)}/{len(PGX_SNPS)}")
    print()

    if fmt == "unknown":
        print("WARNING: Could not detect input file format. Results may be unreliable.",
              file=sys.stderr)

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
            print(f"  * {d['drug']} ({d['brand']}): {d['recommendation']}")
        print()

    # Generate report
    outdir = Path(args.output)
    outdir.mkdir(parents=True, exist_ok=True)
    report = generate_report(args.input, fmt, total_snps, pgx_snps, profiles, drug_results)
    report_path = outdir / "report.md"
    report_path.write_text(report)
    print(f"Report saved: {report_path}")

    # HTML report (always generated alongside markdown)
    from html_report import generate_html_report
    html = generate_html_report(
        input_path=args.input,
        fmt=fmt,
        total_snps=total_snps,
        pgx_snps=pgx_snps,
        profiles=profiles,
        drug_results=drug_results,
        gene_defs=GENE_DEFS,
        pgx_snp_defs=PGX_SNPS,
    )
    html_path = outdir / "report.html"
    html_path.write_text(html)
    print(f"HTML report: {html_path}")

    print("Done.")


if __name__ == "__main__":
    main()
