"""
test_nutrigx.py — Automated test suite for NutriGx Advisor
Run with: pytest tests/test_nutrigx.py -v

Uses a FIXED synthetic patient (synthetic_patient.csv) with known genotypes
so that all assertions are deterministic and reproducible. This file is NOT
meant to showcase the skill — use examples/generate_patient.py for varied demos.
"""

import json
import sys
from pathlib import Path

import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

sys.path.insert(0, str(Path(__file__).parent.parent))

from parse_input import parse_genetic_file
from extract_genotypes import extract_snp_genotypes
from score_variants import compute_nutrient_risk_scores


SKILL_DIR = Path(__file__).parent.parent
SYNTHETIC = Path(__file__).parent / "synthetic_patient.csv"
PANEL     = SKILL_DIR / "data" / "snp_panel.json"


def load_panel():
    with open(PANEL) as f:
        return json.load(f)


# ── Parsing ────────────────────────────────────────────────────────────────────

def test_parse_23andme():
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    assert len(table) >= 20
    assert "rs1801133" in table
    assert table["rs1801133"] in ("CT", "TC")


def test_all_panel_snps_present():
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    genotyped = sum(1 for v in calls.values() if v["status"] in ("found", "allele_mismatch"))
    assert genotyped == len(panel), f"Expected all {len(panel)} SNPs genotyped, got {genotyped}"


# ── Extraction ─────────────────────────────────────────────────────────────────

def test_mthfr_heterozygous():
    """Fixed patient has MTHFR C677T = CT (1 risk allele)."""
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    mthfr = calls["rs1801133"]
    assert mthfr["status"] == "found"
    assert mthfr["risk_count"] == 1


def test_vdr_homozygous_risk():
    """Fixed patient has VDR TaqI = CC (2 risk alleles) → drives Elevated vitamin D score."""
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    vdr = calls["rs731236"]
    assert vdr["status"] == "found"
    assert vdr["risk_count"] == 2


def test_aldh2_ref_homozygous():
    """Fixed patient has ALDH2 = GG (0 risk alleles). Flagged as allele_mismatch
    because GG doesn't contain risk allele A even after strand flip."""
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    aldh2 = calls["rs671"]
    assert aldh2["status"] in ("found", "allele_mismatch")
    # allele_mismatch sets risk_count to None; found sets it to 0
    assert aldh2["risk_count"] in (0, None)


# ── Scoring ────────────────────────────────────────────────────────────────────

def test_scores_structure():
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    scores = compute_nutrient_risk_scores(calls, panel)
    assert "folate" in scores
    assert "vitamin_d" in scores
    assert "omega3" in scores
    assert "alcohol" in scores
    for domain, data in scores.items():
        if data["score"] is not None:
            assert 0.0 <= data["score"] <= 10.0
        assert data["category"] in ("Low", "Moderate", "Elevated", "Unknown")


def test_vitamin_d_elevated():
    """VDR TaqI hom risk → Vitamin D expected Elevated."""
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    scores = compute_nutrient_risk_scores(calls, panel)
    assert scores["vitamin_d"]["category"] == "Elevated"


def test_alcohol_low_or_moderate():
    """ALDH2 GG ref hom → Alcohol expected Low or Moderate (allele_mismatch
    may inflate score slightly since the SNP contribution is uncertain)."""
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    scores = compute_nutrient_risk_scores(calls, panel)
    assert scores["alcohol"]["category"] in ("Low", "Moderate")


def test_folate_not_low():
    """MTHFR C677T het → Folate should be Moderate or Elevated, not Low."""
    panel = load_panel()
    table = parse_genetic_file(str(SYNTHETIC), fmt="23andme")
    calls = extract_snp_genotypes(table, panel)
    scores = compute_nutrient_risk_scores(calls, panel)
    assert scores["folate"]["category"] in ("Moderate", "Elevated")


# ── Hom-ref and empty input (clawbio_bench ng_07, ng_09) ─────────────────────

def test_homref_scores_zero():
    """ng_09: Hom-ref genotype (all ref alleles) should score 0, not allele_mismatch."""
    panel = [
        {"rsid": "rs1801133", "gene": "MTHFR", "risk_allele": "T",
         "ref_allele": "C", "nutrient_domain": "folate", "weight": 1.0},
    ]
    # Genotype CC = homozygous reference (no risk allele T)
    genotype_table = {"rs1801133": "CC"}
    calls = extract_snp_genotypes(genotype_table, panel)
    assert calls["rs1801133"]["status"] == "found", (
        f"Hom-ref should be 'found', got '{calls['rs1801133']['status']}'"
    )
    assert calls["rs1801133"]["risk_count"] == 0, (
        f"Hom-ref should have risk_count=0, got {calls['rs1801133']['risk_count']}"
    )


def test_empty_input_exits_cleanly(tmp_path):
    """ng_07: Empty input file should produce a clean error, not a traceback."""
    import subprocess
    empty_file = tmp_path / "empty.txt"
    empty_file.write_text("")
    result = subprocess.run(
        [sys.executable, str(SKILL_DIR / "nutrigx.py"),
         "--input", str(empty_file), "--output", str(tmp_path / "out")],
        capture_output=True, text=True,
    )
    assert result.returncode != 0, "Empty input should fail"
    assert "ERROR" in result.stderr, (
        f"Should print ERROR message, got stderr: {result.stderr}"
    )


def _lactose_score(genotype: str) -> dict:
    panel = load_panel()
    calls = extract_snp_genotypes({"rs4988235": genotype}, panel)
    scores = compute_nutrient_risk_scores(calls, panel)
    return {
        "call": calls["rs4988235"],
        "score": scores["lactose"],
    }


def test_rs4988235_panel_polarity_and_citation():
    entry = next(s for s in load_panel() if s["rsid"] == "rs4988235")
    assert entry["ref_allele"] == "A"
    assert entry["risk_allele"] == "G"
    assert entry["pmid"] == "11788828"
    assert entry["inheritance"] == "dominant_protective"


def test_rs4988235_aa_and_tt_are_persistent():
    for genotype in ("AA", "TT"):
        result = _lactose_score(genotype)
        assert result["call"]["status"] == "found"
        assert result["call"]["risk_count"] == 0
        assert result["score"]["category"] == "Low"
        assert result["score"]["score"] == 0.0


def test_rs4988235_ag_and_ct_are_persistent_dominant():
    for genotype in ("AG", "GA", "CT", "TC"):
        result = _lactose_score(genotype)
        assert result["call"]["status"] == "found"
        assert result["call"]["risk_count"] == 1
        assert result["score"]["category"] == "Low"
        assert result["score"]["score"] == 0.0


def test_rs4988235_gg_and_cc_are_non_persistent():
    for genotype in ("GG", "CC"):
        result = _lactose_score(genotype)
        assert result["call"]["status"] == "found"
        assert result["call"]["risk_count"] == 2
        assert result["score"]["category"] == "Elevated"
        assert result["score"]["score"] == 10.0


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/nutrigx/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("nutrigx", tmp_path) == snapshot({'nutrigx_report.md': """\
# NutriGx Personalised Nutrition Report

**Generated**: <TIMESTAMP> UTC  \n\
**Tool**: ClawBio NutriGx Advisor v0.2.0  \n\
**Input**: `synthetic_patient.csv`  \n\

> **Disclaimer**: This report is for research and educational purposes only. It does not constitute medical advice. Consult a registered dietitian or clinical geneticist before making significant dietary changes or starting supplements.

## Executive Summary

Analysed **28** of **28** panel SNPs from your genetic data.

**1 nutrient domain(s) at elevated genetic risk:**
- Vitamin D (score: 6.72/10)

**7 domain(s) at moderate genetic risk:**
- Folate / B-Vitamins (score: 5.0/10)
- Omega-3 / LC-PUFA (score: 4.04/10)
- Vitamin A (Beta-carotene) (score: 5.0/10)
- Vitamin B6 (score: 5.0/10)
- Carbohydrate Metabolism (score: 5.0/10)
- Caffeine Metabolism (score: 5.0/10)
- Antioxidant / Detox (score: 5.0/10)

---

## Nutrient Risk Score Overview

| Nutrient Domain | Score (0–10) | Risk Category |
|-----------------|:------------:|:-------------:|
| Vitamin D | 6.72 | 🔴 Elevated (3/3 SNPs tested) |
| Folate / B-Vitamins | 5.0 | 🟡 Moderate (3/3 SNPs tested) |
| Vitamin A (Beta-carotene) | 5.0 | 🟡 Moderate (2/2 SNPs tested) |
| Vitamin B6 | 5.0 | 🟡 Moderate (1/1 SNPs tested) |
| Carbohydrate Metabolism | 5.0 | 🟡 Moderate (2/2 SNPs tested) |
| Caffeine Metabolism | 5.0 | 🟡 Moderate (2/2 SNPs tested) |
| Antioxidant / Detox | 5.0 | 🟡 Moderate (4/4 SNPs tested) |
| Omega-3 / LC-PUFA | 4.04 | 🟡 Moderate (4/4 SNPs tested) |
| Alcohol Metabolism | 2.35 | 🟢 Low (2/2 SNPs tested) |
| Fat Metabolism | 1.57 | 🟢 Low (3/3 SNPs tested) |
| Vitamin C | 0.0 | 🟢 Low (1/1 SNPs tested) |
| Lactose Tolerance | 0.0 | 🟢 Low (1/1 SNPs tested) |

---

## Detailed Findings by Nutrient Domain

### Vitamin D

**Risk Score**: 6.72/10 — **Elevated**  \n\
**SNPs tested**: 3 | **Not on chip**: 0 | **Coverage**: 3/3 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| VDR | rs2228570 | `TC` | 1/2 | Decreased Vitamin D Response |
| VDR | rs731236 | `CC` | 2/2 | Decreased Vitamin D Response |
| GC | rs4588 | `CA` | 1/2 | Decreased Vitamin D Binding |

**Recommendation**

> Genetic predisposition to low vitamin D binding efficiency. Test and maintain 25(OH)D >100 nmol/L. D3+K2 co-supplementation advisable.

---

### Folate / B-Vitamins

**Risk Score**: 5.0/10 — **Moderate**  \n\
**SNPs tested**: 3 | **Not on chip**: 0 | **Coverage**: 3/3 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| MTHFR | rs1801133 | `CT` | 1/2 | Decreased Folate Conversion |
| MTHFR | rs1801131 | `AC` | 1/2 | Decreased Folate Conversion |
| MTR | rs1805087 | `AG` | 1/2 | Increased Homocysteine Risk |

**Recommendation**

> Consider increasing dietary folate. If homocysteine is elevated, 5-MTHF (methylfolate) supplement preferred over folic acid.

---

### Vitamin A (Beta-carotene)

**Risk Score**: 5.0/10 — **Moderate**  \n\
**SNPs tested**: 2 | **Not on chip**: 0 | **Coverage**: 2/2 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| BCMO1 | rs7501331 | `CT` | 1/2 | Decreased Beta Carotene Conversion |
| BCMO1 | rs12934922 | `AT` | 1/2 | Decreased Beta Carotene Conversion |

**Recommendation**

> Reduced BCMO1 activity. Include preformed vitamin A (liver, eggs, dairy) alongside carotenoid-rich vegetables.

---

### Vitamin B6

**Risk Score**: 5.0/10 — **Moderate**  \n\
**SNPs tested**: 1 | **Not on chip**: 0 | **Coverage**: 1/1 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| ALPL | rs1256335 | `CT` | 1/2 | Decreased B6 Utilisation |

**Recommendation**

> No specific recommendation available.

---

### Carbohydrate Metabolism

**Risk Score**: 5.0/10 — **Moderate**  \n\
**SNPs tested**: 2 | **Not on chip**: 0 | **Coverage**: 2/2 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| FTO | rs9939609 | `AT` | 1/2 | Increased Fat Mass Energy Dysregulation |
| TCF7L2 | rs7903146 | `CT` | 1/2 | Impaired Carbohydrate Metabolism |

**Recommendation**

> Glycaemic index awareness recommended. Favour whole-grain, low-GI carbohydrates.

---

### Caffeine Metabolism

**Risk Score**: 5.0/10 — **Moderate**  \n\
**SNPs tested**: 2 | **Not on chip**: 0 | **Coverage**: 2/2 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| CYP1A2 | rs762551 | `CA` | 1/2 | Slow Caffeine Metabolism |
| AHR | rs4410790 | `CT` | 1/2 | Reduced Cyp1A2 Induction |

**Recommendation**

> Intermediate metaboliser. Caffeine ≤200 mg/day recommended.

---

### Antioxidant / Detox

**Risk Score**: 5.0/10 — **Moderate**  \n\
**SNPs tested**: 4 | **Not on chip**: 0 | **Coverage**: 4/4 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| SOD2 | rs4880 | `CT` | 1/2 | Decreased Mitochondrial Antioxidant |
| GPX1 | rs1050450 | `CT` | 1/2 | Decreased Glutathione Peroxidase |
| NQO1 | rs1800566 | `CT` | 1/2 | Decreased Coq10 Regeneration |
| COMT | rs4680 | `GA` | 1/2 | Altered Methylation Catechol Metabolism |

**Recommendation**

> Some oxidative stress pathway variants. Increase dietary antioxidants; selenium and CoQ10 food sources beneficial.

---

### Omega-3 / LC-PUFA

**Risk Score**: 4.04/10 — **Moderate**  \n\
**SNPs tested**: 4 | **Not on chip**: 0 | **Coverage**: 4/4 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| FADS1 | rs174546 | `TC` | 1/2 | Decreased Lc Pufa Synthesis |
| FADS2 | rs1535 | `AG` | 1/2 | Altered Omega6 Omega3 Ratio |
| ELOVL2 | rs953413 | `AG` | 1/2 | Lower Epa To Dha Conversion |
| APOE | rs429358 | `TT` | 0/2 | Increased Ldl On Saturated Fat |

**Recommendation**

> Variants associated with somewhat lower conversion of ALA and EPA to DHA. Favour direct EPA/DHA sources (oily fish, algae oil) over relying on ALA from plant foods.

---

### Alcohol Metabolism

**Risk Score**: 2.35/10 — **Low**  \n\
**SNPs tested**: 2 | **Not on chip**: 0 | **Coverage**: 2/2 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| ADH1B | rs1229984 | `AG` | 1/2 | Altered Acetaldehyde Accumulation |
| ALDH2 | rs671 | `GG` | 0/2 | Impaired Acetaldehyde Clearance |

**Recommendation**

> Standard alcohol metabolism. General guidelines apply.

---

### Fat Metabolism

**Risk Score**: 1.57/10 — **Low**  \n\
**SNPs tested**: 3 | **Not on chip**: 0 | **Coverage**: 3/3 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| PPARG | rs1801282 | `CG` | 1/2 | Reduced Fat Oxidation |
| APOA5 | rs662799 | `TT` | 0/2 | Elevated Triglycerides On Fat |
| APOE | rs7412 | `CC` | 0/2 | Apoe Epsilon2 Lower Ldl |

**Recommendation**

> No specific fat quality concerns from genetics alone.

---

### Vitamin C

**Risk Score**: 0.0/10 — **Low**  \n\
**SNPs tested**: 1 | **Not on chip**: 0 | **Coverage**: 1/1 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| SLC23A1 | rs33972313 | `CC` | 0/2 | Decreased Vitamin C Reabsorption |

**Recommendation**

> Standard dietary vitamin C (citrus, peppers, kiwi) sufficient.

---

### Lactose Tolerance

**Risk Score**: 0.0/10 — **Low**  \n\
**SNPs tested**: 1 | **Not on chip**: 0 | **Coverage**: 1/1 SNPs tested

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| MCM6 | rs4988235 | `GA` | 1/2 | Lactase Non Persistence |

**Recommendation**

> Likely lactase-persistent. Can tolerate dairy normally.

---

## Supplement Interaction Notes

| Supplement | Relevant Genes | Caution |
|------------|---------------|---------|
| Folic acid (synthetic) | MTHFR | Prefer 5-MTHF if MTHFR variant present |
| Vitamin D3 | VDR, GC | K2 co-supplementation improves utilisation |
| Fish oil / EPA+DHA | FADS1, FADS2, ELOVL2 | Higher dose needed if FADS variants present |
| Vitamin A (retinol) | BCMO1 | Do not exceed UL (3000 μg RAE/day) |
| Iron | HFE (not in panel) | Not assessed; request separate iron panel |
| CoQ10 | NQO1 | Ubiquinol form preferred if NQO1 P187S homozygous |

---

## Reproducibility

This report was generated deterministically. See `commands.sh` and `environment.yml`
in the `reproducibility/` subdirectory to reproduce this analysis on any machine.

## References

SNP-nutrient associations sourced from GWAS Catalog, ClinVar, and CPIC guidelines.
Full citations available in `skills/nutrigx/SKILL.md`.
""", 'result.json': """\
{
  "skill": "nutrigx",
  "version": "0.2.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "sha256:353c96a2f52acfc0",
  "datasets": {},
  "summary": {
    "total_variants_loaded": 28,
    "panel_snps_tested": 28,
    "panel_snps_missing": 0,
    "panel_size": 28,
    "elevated_domains": [
      "vitamin_d"
    ],
    "moderate_domains": [
      "folate",
      "omega3",
      "vitamin_a",
      "vitamin_b6",
      "carbohydrate",
      "caffeine",
      "antioxidant"
    ],
    "domains_assessed": 12
  },
  "data": {
    "risk_scores": {
      "folate": {
        "score": 5.0,
        "category": "Moderate",
        "contributing_snps": [
          {
            "rsid": "rs1801133",
            "gene": "MTHFR",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.85,
            "effect_direction": "decreased_folate_conversion"
          },
          {
            "rsid": "rs1801131",
            "gene": "MTHFR",
            "genotype": "AC",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.5,
            "effect_direction": "decreased_folate_conversion"
          },
          {
            "rsid": "rs1805087",
            "gene": "MTR",
            "genotype": "AG",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.4,
            "effect_direction": "increased_homocysteine_risk"
          }
        ],
        "tested_snps": 3,
        "missing_snps": 0,
        "coverage": "3/3 SNPs tested"
      },
      "vitamin_d": {
        "score": 6.72,
        "category": "Elevated",
        "contributing_snps": [
          {
            "rsid": "rs2228570",
            "gene": "VDR",
            "genotype": "TC",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.6,
            "effect_direction": "decreased_vitamin_d_response"
          },
          {
            "rsid": "rs731236",
            "gene": "VDR",
            "genotype": "CC",
            "risk_count": 2,
            "raw_score": 1.0,
            "weight": 0.55,
            "effect_direction": "decreased_vitamin_d_response"
          },
          {
            "rsid": "rs4588",
            "gene": "GC",
            "genotype": "CA",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.45,
            "effect_direction": "decreased_vitamin_d_binding"
          }
        ],
        "tested_snps": 3,
        "missing_snps": 0,
        "coverage": "3/3 SNPs tested"
      },
      "omega3": {
        "score": 4.04,
        "category": "Moderate",
        "contributing_snps": [
          {
            "rsid": "rs174546",
            "gene": "FADS1",
            "genotype": "TC",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.8,
            "effect_direction": "decreased_lc_pufa_synthesis"
          },
          {
            "rsid": "rs1535",
            "gene": "FADS2",
            "genotype": "AG",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.7,
            "effect_direction": "altered_omega6_omega3_ratio"
          },
          {
            "rsid": "rs953413",
            "gene": "ELOVL2",
            "genotype": "AG",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.6,
            "effect_direction": "lower_epa_to_dha_conversion"
          },
          {
            "rsid": "rs429358",
            "gene": "APOE",
            "genotype": "TT",
            "risk_count": 0,
            "raw_score": 0.0,
            "weight": 0.5,
            "effect_direction": "increased_ldl_on_saturated_fat"
          }
        ],
        "tested_snps": 4,
        "missing_snps": 0,
        "coverage": "4/4 SNPs tested"
      },
      "vitamin_a": {
        "score": 5.0,
        "category": "Moderate",
        "contributing_snps": [
          {
            "rsid": "rs7501331",
            "gene": "BCMO1",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.7,
            "effect_direction": "decreased_beta_carotene_conversion"
          },
          {
            "rsid": "rs12934922",
            "gene": "BCMO1",
            "genotype": "AT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.6,
            "effect_direction": "decreased_beta_carotene_conversion"
          }
        ],
        "tested_snps": 2,
        "missing_snps": 0,
        "coverage": "2/2 SNPs tested"
      },
      "vitamin_c": {
        "score": 0.0,
        "category": "Low",
        "contributing_snps": [
          {
            "rsid": "rs33972313",
            "gene": "SLC23A1",
            "genotype": "CC",
            "risk_count": 0,
            "raw_score": 0.0,
            "weight": 0.5,
            "effect_direction": "decreased_vitamin_c_reabsorption"
          }
        ],
        "tested_snps": 1,
        "missing_snps": 0,
        "coverage": "1/1 SNPs tested"
      },
      "vitamin_b6": {
        "score": 5.0,
        "category": "Moderate",
        "contributing_snps": [
          {
            "rsid": "rs1256335",
            "gene": "ALPL",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.45,
            "effect_direction": "decreased_b6_utilisation"
          }
        ],
        "tested_snps": 1,
        "missing_snps": 0,
        "coverage": "1/1 SNPs tested"
      },
      "carbohydrate": {
        "score": 5.0,
        "category": "Moderate",
        "contributing_snps": [
          {
            "rsid": "rs9939609",
            "gene": "FTO",
            "genotype": "AT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.7,
            "effect_direction": "increased_fat_mass_energy_dysregulation"
          },
          {
            "rsid": "rs7903146",
            "gene": "TCF7L2",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.75,
            "effect_direction": "impaired_carbohydrate_metabolism"
          }
        ],
        "tested_snps": 2,
        "missing_snps": 0,
        "coverage": "2/2 SNPs tested"
      },
      "fat_metabolism": {
        "score": 1.57,
        "category": "Low",
        "contributing_snps": [
          {
            "rsid": "rs1801282",
            "gene": "PPARG",
            "genotype": "CG",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.55,
            "effect_direction": "reduced_fat_oxidation"
          },
          {
            "rsid": "rs662799",
            "gene": "APOA5",
            "genotype": "TT",
            "risk_count": 0,
            "raw_score": 0.0,
            "weight": 0.65,
            "effect_direction": "elevated_triglycerides_on_fat"
          },
          {
            "rsid": "rs7412",
            "gene": "APOE",
            "genotype": "CC",
            "risk_count": 0,
            "raw_score": 0.0,
            "weight": 0.55,
            "effect_direction": "apoe_epsilon2_lower_ldl"
          }
        ],
        "tested_snps": 3,
        "missing_snps": 0,
        "coverage": "3/3 SNPs tested"
      },
      "caffeine": {
        "score": 5.0,
        "category": "Moderate",
        "contributing_snps": [
          {
            "rsid": "rs762551",
            "gene": "CYP1A2",
            "genotype": "CA",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.85,
            "effect_direction": "slow_caffeine_metabolism"
          },
          {
            "rsid": "rs4410790",
            "gene": "AHR",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.4,
            "effect_direction": "reduced_cyp1a2_induction"
          }
        ],
        "tested_snps": 2,
        "missing_snps": 0,
        "coverage": "2/2 SNPs tested"
      },
      "alcohol": {
        "score": 2.35,
        "category": "Low",
        "contributing_snps": [
          {
            "rsid": "rs1229984",
            "gene": "ADH1B",
            "genotype": "AG",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.8,
            "effect_direction": "altered_acetaldehyde_accumulation"
          },
          {
            "rsid": "rs671",
            "gene": "ALDH2",
            "genotype": "GG",
            "risk_count": 0,
            "raw_score": 0.0,
            "weight": 0.9,
            "effect_direction": "impaired_acetaldehyde_clearance"
          }
        ],
        "tested_snps": 2,
        "missing_snps": 0,
        "coverage": "2/2 SNPs tested"
      },
      "lactose": {
        "score": 0.0,
        "category": "Low",
        "contributing_snps": [
          {
            "rsid": "rs4988235",
            "gene": "MCM6",
            "genotype": "GA",
            "risk_count": 1,
            "raw_score": 0.0,
            "weight": 0.95,
            "effect_direction": "lactase_non_persistence"
          }
        ],
        "tested_snps": 1,
        "missing_snps": 0,
        "coverage": "1/1 SNPs tested"
      },
      "antioxidant": {
        "score": 5.0,
        "category": "Moderate",
        "contributing_snps": [
          {
            "rsid": "rs4880",
            "gene": "SOD2",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.55,
            "effect_direction": "decreased_mitochondrial_antioxidant"
          },
          {
            "rsid": "rs1050450",
            "gene": "GPX1",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.5,
            "effect_direction": "decreased_glutathione_peroxidase"
          },
          {
            "rsid": "rs1800566",
            "gene": "NQO1",
            "genotype": "CT",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.6,
            "effect_direction": "decreased_coq10_regeneration"
          },
          {
            "rsid": "rs4680",
            "gene": "COMT",
            "genotype": "GA",
            "risk_count": 1,
            "raw_score": 0.5,
            "weight": 0.45,
            "effect_direction": "altered_methylation_catechol_metabolism"
          }
        ],
        "tested_snps": 4,
        "missing_snps": 0,
        "coverage": "4/4 SNPs tested"
      }
    },
    "snp_calls": {
      "rs1801133": {
        "rsid": "rs1801133",
        "gene": "MTHFR",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "folate"
      },
      "rs1801131": {
        "rsid": "rs1801131",
        "gene": "MTHFR",
        "status": "found",
        "genotype": "AC",
        "normalised": "AC",
        "risk_allele": "C",
        "risk_count": 1,
        "nutrient_domain": "folate"
      },
      "rs1805087": {
        "rsid": "rs1805087",
        "gene": "MTR",
        "status": "found",
        "genotype": "AG",
        "normalised": "AG",
        "risk_allele": "G",
        "risk_count": 1,
        "nutrient_domain": "folate"
      },
      "rs2228570": {
        "rsid": "rs2228570",
        "gene": "VDR",
        "status": "found",
        "genotype": "TC",
        "normalised": "TC",
        "risk_allele": "C",
        "risk_count": 1,
        "nutrient_domain": "vitamin_d"
      },
      "rs731236": {
        "rsid": "rs731236",
        "gene": "VDR",
        "status": "found",
        "genotype": "CC",
        "normalised": "CC",
        "risk_allele": "C",
        "risk_count": 2,
        "nutrient_domain": "vitamin_d"
      },
      "rs4588": {
        "rsid": "rs4588",
        "gene": "GC",
        "status": "found",
        "genotype": "CA",
        "normalised": "CA",
        "risk_allele": "A",
        "risk_count": 1,
        "nutrient_domain": "vitamin_d"
      },
      "rs174546": {
        "rsid": "rs174546",
        "gene": "FADS1",
        "status": "found",
        "genotype": "TC",
        "normalised": "TC",
        "risk_allele": "C",
        "risk_count": 1,
        "nutrient_domain": "omega3"
      },
      "rs1535": {
        "rsid": "rs1535",
        "gene": "FADS2",
        "status": "found",
        "genotype": "AG",
        "normalised": "AG",
        "risk_allele": "A",
        "risk_count": 1,
        "nutrient_domain": "omega3"
      },
      "rs953413": {
        "rsid": "rs953413",
        "gene": "ELOVL2",
        "status": "found",
        "genotype": "AG",
        "normalised": "AG",
        "risk_allele": "A",
        "risk_count": 1,
        "nutrient_domain": "omega3"
      },
      "rs429358": {
        "rsid": "rs429358",
        "gene": "APOE",
        "status": "found",
        "genotype": "TT",
        "normalised": "TT",
        "risk_allele": "C",
        "risk_count": 0,
        "nutrient_domain": "omega3"
      },
      "rs7501331": {
        "rsid": "rs7501331",
        "gene": "BCMO1",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "vitamin_a"
      },
      "rs12934922": {
        "rsid": "rs12934922",
        "gene": "BCMO1",
        "status": "found",
        "genotype": "AT",
        "normalised": "AT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "vitamin_a"
      },
      "rs33972313": {
        "rsid": "rs33972313",
        "gene": "SLC23A1",
        "status": "found",
        "genotype": "CC",
        "normalised": "CC",
        "risk_allele": "T",
        "risk_count": 0,
        "nutrient_domain": "vitamin_c"
      },
      "rs1256335": {
        "rsid": "rs1256335",
        "gene": "ALPL",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "vitamin_b6"
      },
      "rs9939609": {
        "rsid": "rs9939609",
        "gene": "FTO",
        "status": "found",
        "genotype": "AT",
        "normalised": "AT",
        "risk_allele": "A",
        "risk_count": 1,
        "nutrient_domain": "carbohydrate"
      },
      "rs7903146": {
        "rsid": "rs7903146",
        "gene": "TCF7L2",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "carbohydrate"
      },
      "rs1801282": {
        "rsid": "rs1801282",
        "gene": "PPARG",
        "status": "found",
        "genotype": "CG",
        "normalised": "CG",
        "risk_allele": "G",
        "risk_count": 1,
        "nutrient_domain": "fat_metabolism"
      },
      "rs662799": {
        "rsid": "rs662799",
        "gene": "APOA5",
        "status": "found",
        "genotype": "TT",
        "normalised": "TT",
        "risk_allele": "C",
        "risk_count": 0,
        "nutrient_domain": "fat_metabolism"
      },
      "rs7412": {
        "rsid": "rs7412",
        "gene": "APOE",
        "status": "found",
        "genotype": "CC",
        "normalised": "CC",
        "risk_allele": "T",
        "risk_count": 0,
        "nutrient_domain": "fat_metabolism"
      },
      "rs762551": {
        "rsid": "rs762551",
        "gene": "CYP1A2",
        "status": "found",
        "genotype": "CA",
        "normalised": "CA",
        "risk_allele": "C",
        "risk_count": 1,
        "nutrient_domain": "caffeine"
      },
      "rs4410790": {
        "rsid": "rs4410790",
        "gene": "AHR",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "caffeine"
      },
      "rs1229984": {
        "rsid": "rs1229984",
        "gene": "ADH1B",
        "status": "found",
        "genotype": "AG",
        "normalised": "AG",
        "risk_allele": "A",
        "risk_count": 1,
        "nutrient_domain": "alcohol"
      },
      "rs671": {
        "rsid": "rs671",
        "gene": "ALDH2",
        "status": "found",
        "genotype": "GG",
        "normalised": "GG",
        "risk_allele": "A",
        "risk_count": 0,
        "nutrient_domain": "alcohol"
      },
      "rs4988235": {
        "rsid": "rs4988235",
        "gene": "MCM6",
        "status": "found",
        "genotype": "GA",
        "normalised": "GA",
        "risk_allele": "G",
        "risk_count": 1,
        "nutrient_domain": "lactose"
      },
      "rs4880": {
        "rsid": "rs4880",
        "gene": "SOD2",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "C",
        "risk_count": 1,
        "nutrient_domain": "antioxidant"
      },
      "rs1050450": {
        "rsid": "rs1050450",
        "gene": "GPX1",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "antioxidant"
      },
      "rs1800566": {
        "rsid": "rs1800566",
        "gene": "NQO1",
        "status": "found",
        "genotype": "CT",
        "normalised": "CT",
        "risk_allele": "T",
        "risk_count": 1,
        "nutrient_domain": "antioxidant"
      },
      "rs4680": {
        "rsid": "rs4680",
        "gene": "COMT",
        "status": "found",
        "genotype": "GA",
        "normalised": "GA",
        "risk_allele": "A",
        "risk_count": 1,
        "nutrient_domain": "antioxidant"
      }
    }
  }
}\
"""})
