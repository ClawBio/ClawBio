"""
test_genome_compare.py — Automated test suite for ClawBio Genome Comparator

Run with: pytest skills/genome-compare/tests/test_genome_compare.py -v

Uses FIXED demo data (Manuel Corpas vs George Church) so all assertions
are deterministic and reproducible.
"""

import json
import sys
from pathlib import Path

import numpy as np
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from genome_compare import (
    REFERENCE_FILE,
    DEMO_PATIENT_FILE,
    AIMS_PANEL_FILE,
    MANUEL_ANCESTRY_FILE,
    IBS_REFERENCE,
    _ibs_at_site,
    _count_alt_alleles,
    _parse_genotype_file as parse_23andme_extended,
    compute_ibs,
    load_aims_panel,
    estimate_ancestry,
    generate_report,
    run_comparison,
)


# ------------------------------------------------------------------ #
# Shared parses
# ------------------------------------------------------------------ #
#
# Both reference files are ~5 MB gzipped and were re-parsed by every test that
# needed them -- 14 parses across 28 tests, which was the whole of this file's
# ~42s runtime (`--durations` showed the cost spread evenly, not one slow test).
#
# Module scope is safe here because neither compute_ibs nor estimate_ancestry
# mutates the structures it is handed, and no test writes to them. If that ever
# changes, these must go back to per-test parses or hand out copies.


@pytest.fixture(scope="module")
def reference_genome():
    """Parsed George Church reference: (genotypes, positions)."""
    return parse_23andme_extended(REFERENCE_FILE)


@pytest.fixture(scope="module")
def demo_patient():
    """Parsed Manuel Corpas demo patient: (genotypes, positions)."""
    return parse_23andme_extended(DEMO_PATIENT_FILE)


@pytest.fixture(scope="module")
def demo_comparison(tmp_path_factory):
    """One full demo run, shared by the report tests: (result, output_dir).

    Four tests ran this identical pipeline and then asserted on different parts
    of it, at ~4.5s each. tmp_path_factory rather than tmp_path because the
    latter is function-scoped and cannot be used from a module-scoped fixture.
    """
    out = tmp_path_factory.mktemp("genome_compare_demo") / "report"
    result = run_comparison(
        input_path=DEMO_PATIENT_FILE,
        output_dir=out,
        no_figures=True,
        is_demo=True,
    )
    return result, out


# ------------------------------------------------------------------ #
# Parsing tests
# ------------------------------------------------------------------ #


def test_parse_reference_loads(reference_genome):
    """George Church's file parses and has >500k SNPs."""
    geno, pos = reference_genome
    assert len(geno) > 500_000, f"Expected >500k SNPs, got {len(geno)}"


def test_parse_demo_patient(demo_patient):
    """Manuel's file (demo patient) parses and has >500k SNPs."""
    geno, pos = demo_patient
    assert len(geno) > 500_000, f"Expected >500k SNPs, got {len(geno)}"


def test_parse_returns_chromosome_metadata(demo_patient):
    """Positions dict has chrom/pos for parsed SNPs."""
    geno, pos = demo_patient
    # Check a few random entries
    for rsid in list(geno.keys())[:10]:
        assert rsid in pos, f"{rsid} missing from positions"
        assert "chrom" in pos[rsid]
        assert "pos" in pos[rsid]


def test_parse_skips_comments(reference_genome):
    """Lines starting with # are skipped."""
    geno, pos = reference_genome
    for rsid in geno:
        assert not rsid.startswith("#")


# ------------------------------------------------------------------ #
# IBS unit tests
# ------------------------------------------------------------------ #


def test_ibs_identical_homozygous():
    """Same homozygous genotype → IBS = 2."""
    assert _ibs_at_site("AA", "AA") == 2


def test_ibs_identical_heterozygous():
    """Same heterozygous genotype → IBS = 2."""
    assert _ibs_at_site("AG", "AG") == 2


def test_ibs_het_order_independent():
    """AT vs TA → IBS = 2 (order doesn't matter)."""
    assert _ibs_at_site("AT", "TA") == 2


def test_ibs_one_shared():
    """AT vs AA → IBS = 1."""
    assert _ibs_at_site("AT", "AA") == 1


def test_ibs_none_shared():
    """AA vs TT → IBS = 0."""
    assert _ibs_at_site("AA", "TT") == 0


def test_ibs_haploid_match():
    """Haploid A vs A → IBS = 2."""
    assert _ibs_at_site("A", "A") == 2


def test_ibs_haploid_mismatch():
    """Haploid A vs T → IBS = 0."""
    assert _ibs_at_site("A", "T") == 0


def test_compute_ibs_score_range(demo_patient, reference_genome):
    """IBS score is in [0, 1]."""
    geno_a, pos_a = demo_patient
    geno_b, _ = reference_genome
    score, n_overlap, n_concordant, per_chrom = compute_ibs(geno_a, geno_b, pos_a)
    assert 0.0 <= score <= 1.0
    assert n_overlap > 0
    assert n_concordant >= 0


def test_compute_ibs_self_gives_one(demo_patient):
    """Comparing a genome against itself gives IBS = 1.0."""
    geno, pos = demo_patient
    score, n_overlap, n_concordant, _ = compute_ibs(geno, geno, pos)
    assert score == 1.0
    assert n_concordant == n_overlap


# ------------------------------------------------------------------ #
# Alt allele counting
# ------------------------------------------------------------------ #


def test_count_alt_direct():
    """Direct allele match: AG with ref=A, alt=G → 1."""
    assert _count_alt_alleles("AG", "A", "G") == 1


def test_count_alt_hom_ref():
    """Homozygous ref: AA with ref=A, alt=G → 0."""
    assert _count_alt_alleles("AA", "A", "G") == 0


def test_count_alt_hom_alt():
    """Homozygous alt: GG with ref=A, alt=G → 2."""
    assert _count_alt_alleles("GG", "A", "G") == 2


def test_count_alt_skips_ambiguous():
    """A/T complement pair is ambiguous → None."""
    assert _count_alt_alleles("AA", "A", "T") is None


def test_count_alt_skips_cg_ambiguous():
    """C/G complement pair is ambiguous → None."""
    assert _count_alt_alleles("CC", "C", "G") is None


# ------------------------------------------------------------------ #
# Per-chromosome breakdown
# ------------------------------------------------------------------ #


def test_per_chromosome_keys(demo_patient, reference_genome):
    """All autosomes with data appear in breakdown."""
    geno_a, pos_a = demo_patient
    geno_b, _ = reference_genome
    _, _, _, per_chrom = compute_ibs(geno_a, geno_b, pos_a)
    # At minimum, chromosomes 1-22 should be present
    for c in [str(i) for i in range(1, 23)]:
        assert c in per_chrom, f"Chromosome {c} missing from breakdown"


def test_per_chromosome_ibs_range(demo_patient, reference_genome):
    """Each chromosome's IBS is in [0, 1]."""
    geno_a, pos_a = demo_patient
    geno_b, _ = reference_genome
    _, _, _, per_chrom = compute_ibs(geno_a, geno_b, pos_a)
    for chrom, data in per_chrom.items():
        assert 0.0 <= data["ibs"] <= 1.0, f"Chr {chrom} IBS out of range"


# ------------------------------------------------------------------ #
# Ancestry estimation
# ------------------------------------------------------------------ #


def test_load_aims_panel():
    """AIMs panel loads with expected fields."""
    markers, pops = load_aims_panel(AIMS_PANEL_FILE)
    assert len(markers) > 30
    assert "AFR" in pops
    assert "EUR" in pops
    for m in markers:
        assert "rsid" in m
        assert "ref" in m
        assert "alt" in m


def test_ancestry_proportions_sum_to_one(demo_patient):
    """Ancestry proportions sum to approximately 1.0."""
    geno, _ = demo_patient
    markers, pops = load_aims_panel(AIMS_PANEL_FILE)
    result = estimate_ancestry(geno, markers, pops)
    total = sum(result["continental"].values())
    assert abs(total - 1.0) < 0.01, f"Sum = {total}"


def test_ancestry_returns_all_populations(demo_patient):
    """All 5 superpopulations present in output."""
    geno, _ = demo_patient
    markers, pops = load_aims_panel(AIMS_PANEL_FILE)
    result = estimate_ancestry(geno, markers, pops)
    for pop in ["AFR", "EUR", "EAS", "SAS", "AMR"]:
        assert pop in result["continental"]


def test_demo_patient_eur_is_top(demo_patient):
    """For Manuel Corpas (European), EUR should be the top ancestry."""
    geno, _ = demo_patient
    markers, pops = load_aims_panel(AIMS_PANEL_FILE)
    result = estimate_ancestry(geno, markers, pops)
    top = max(result["continental"], key=result["continental"].get)
    assert top == "EUR", f"Expected EUR as top, got {top}"


# ------------------------------------------------------------------ #
# Report tests
# ------------------------------------------------------------------ #


def test_report_contains_key_sections(demo_comparison):
    """Report has expected markdown headers."""
    _, out = demo_comparison
    report = (out / "report.md").read_text()
    assert "# Genome Comparison Report" in report
    assert "## Summary" in report
    assert "## Identity By State Analysis" in report
    assert "## Ancestry Composition" in report


def test_report_contains_disclaimer(demo_comparison):
    """Standard ClawBio disclaimer is present."""
    _, out = demo_comparison
    report = (out / "report.md").read_text()
    assert "not a medical device" in report


def test_report_contains_methods(demo_comparison):
    """Methods section is present."""
    _, out = demo_comparison
    report = (out / "report.md").read_text()
    assert "## Methods" in report
    assert "Identity By State" in report


# ------------------------------------------------------------------ #
# End-to-end
# ------------------------------------------------------------------ #


def test_end_to_end_demo(demo_comparison):
    """Full pipeline runs and produces expected outputs."""
    result, out = demo_comparison
    assert result["ibs_score"] > 0.6
    assert result["ibs_score"] < 0.9
    assert result["n_overlap"] > 400_000
    assert result["n_concordant"] > 0
    assert "EUR" in result["ancestry"]["continental"]
    assert (out / "report.md").exists()
