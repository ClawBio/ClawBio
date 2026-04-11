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
        [sys.executable, str(SKILL_DIR / "nutrigx_advisor.py"),
         "--input", str(empty_file), "--output", str(tmp_path / "out")],
        capture_output=True, text=True,
    )
    assert result.returncode != 0, "Empty input should fail"
    assert "ERROR" in result.stderr, (
        f"Should print ERROR message, got stderr: {result.stderr}"
    )
