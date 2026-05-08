"""
tests/test_methylation_cycle.py
ClawBio · claw-methylation-cycle v0.1.1

Test suite — required per ClawBio TDD mandate.
Covers: parsing, compute_* functions with known I/O pairs, Safety Rule 6,
compound heterozygosity detection, NMC/BH4 scores, and missing SNP handling.
"""

import sys
import tempfile
from pathlib import Path

import pytest

# Allow running from repo root or tests/ directory
sys.path.insert(0, str(Path(__file__).parent.parent))
from methylation_cycle import (
    PANEL,
    analyse,
    compute_mthfr_combined,
    count_risk_alleles,
    estimate_activity,
    parse_genotype_file,
    _dopamine_impact,
    _serotonin_impact,
)


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DEMO_INPUT = """\
# ClawBio demo genotype file - synthetic data
# rsid\tchromosome\tposition\tgenotype
rs1801133\t1\t11856378\tTT
rs1801131\t1\t11854476\tCC
rs1801394\t5\t7870860\tAG
rs1805087\t1\t236894897\tAG
rs234706\t21\t43351116\tCT
rs3733890\t5\t78396840\tAG
rs1979277\t17\t18025410\tCT
rs4680\t22\t19951271\tAA
rs819147\t20\t33764554\tTT
"""

PARTIAL_INPUT = """\
# Only MTHFR SNPs present
rs1801133\t1\t11856378\tCT
rs1801131\t1\t11854476\tAC
"""

WILDTYPE_INPUT = """\
rs1801133\t1\t11856378\tCC
rs1801131\t1\t11854476\tAA
rs1801394\t5\t7870860\tAA
rs1805087\t1\t236894897\tAA
rs234706\t21\t43351116\tCC
rs3733890\t5\t78396840\tGG
rs1979277\t17\t18025410\tCC
rs4680\t22\t19951271\tGG
rs819147\t20\t33764554\tCC
"""


def write_temp(content: str) -> Path:
    tmp = tempfile.NamedTemporaryFile(
        mode="w", suffix=".txt", delete=False, encoding="utf-8"
    )
    tmp.write(content)
    tmp.close()
    return Path(tmp.name)


# ---------------------------------------------------------------------------
# parse_genotype_file
# ---------------------------------------------------------------------------

class TestParseGenotypeFile:
    def test_parses_demo_input(self):
        path = write_temp(DEMO_INPUT)
        result = parse_genotype_file(path)
        assert "rs1801133" in result
        assert result["rs1801133"] == "TT"

    def test_skips_comment_lines(self):
        path = write_temp(DEMO_INPUT)
        result = parse_genotype_file(path)
        # Comment lines begin with '#'; none should appear as keys
        for key in result:
            assert not key.startswith("#")

    def test_returns_uppercase_genotype(self):
        content = "rs4680\t22\t19951271\taa\n"
        path = write_temp(content)
        result = parse_genotype_file(path)
        assert result["rs4680"] == "AA"

    def test_empty_file_returns_empty_dict(self):
        path = write_temp("# only comments\n")
        result = parse_genotype_file(path)
        assert result == {}


# ---------------------------------------------------------------------------
# count_risk_alleles
# ---------------------------------------------------------------------------

class TestCountRiskAlleles:
    def test_homozygous_risk(self):
        assert count_risk_alleles("TT", "T") == 2

    def test_heterozygous(self):
        assert count_risk_alleles("CT", "T") == 1

    def test_wildtype(self):
        assert count_risk_alleles("CC", "T") == 0

    def test_different_allele(self):
        assert count_risk_alleles("AG", "G") == 1


# ---------------------------------------------------------------------------
# estimate_activity
# ---------------------------------------------------------------------------

class TestEstimateActivity:
    def test_wildtype_returns_100(self):
        snp = PANEL["rs1801133"]
        assert estimate_activity(snp, 0) == 100

    def test_heterozygous_mthfr_c677t(self):
        # MTHFR C677T het: ~65% (doi:10.1038/ng0595-111)
        snp = PANEL["rs1801133"]
        assert estimate_activity(snp, 1) == 65

    def test_homozygous_mthfr_c677t(self):
        # MTHFR C677T hom: ~30% (doi:10.1038/ng0595-111)
        snp = PANEL["rs1801133"]
        assert estimate_activity(snp, 2) == 30

    def test_heterozygous_comt(self):
        # COMT Val158Met het: ~65% (doi:10.1097/00008571-199606000-00007)
        snp = PANEL["rs4680"]
        assert estimate_activity(snp, 1) == 65

    def test_homozygous_comt(self):
        # COMT Val158Met hom: ~25% (doi:10.1523/JNEUROSCI.4106-03.2004)
        snp = PANEL["rs4680"]
        assert estimate_activity(snp, 2) == 25


# ---------------------------------------------------------------------------
# compute_mthfr_combined
# ---------------------------------------------------------------------------

class TestComputeMthfrCombined:
    def test_wildtype_both(self):
        # No risk alleles: 100% activity
        assert compute_mthfr_combined(0, 0) == 100

    def test_c677t_heterozygous_only(self):
        # C677T het only: ~65% (doi:10.1038/ng0595-111)
        assert compute_mthfr_combined(1, 0) == 65

    def test_c677t_homozygous_only(self):
        # C677T hom only: ~30% (doi:10.1038/ng0595-111)
        assert compute_mthfr_combined(2, 0) == 30

    def test_a1298c_heterozygous_only(self):
        # A1298C het only: ~80% (doi:10.1086/301927)
        assert compute_mthfr_combined(0, 1) == 80

    def test_a1298c_homozygous_only(self):
        # A1298C hom only: ~60% (doi:10.1086/301927)
        assert compute_mthfr_combined(0, 2) == 60

    def test_compound_het(self):
        # C677T het + A1298C het: ~15% synergistic reduction
        # source: doi:10.1086/301927 (van der Put et al. 1998)
        assert compute_mthfr_combined(1, 1) == 15

    def test_compound_hom_c677t_het_a1298c(self):
        # C677T hom + A1298C het: ~12% (most severe configuration)
        # source: doi:10.1086/301927
        assert compute_mthfr_combined(2, 1) == 12

    def test_compound_het_c677t_hom_a1298c(self):
        # C677T het + A1298C hom: ~20%
        # source: doi:10.1086/301927
        assert compute_mthfr_combined(1, 2) == 20


# ---------------------------------------------------------------------------
# analyse - Safety Rule 6: missing SNPs
# ---------------------------------------------------------------------------

class TestSafetyRule6:
    def test_missing_snp_activity_is_none(self):
        """Missing SNPs must have activity_pct=None, not 100."""
        result = analyse({})
        for rsid, r in result["gene_results"].items():
            assert r["activity_pct"] is None, (
                f"{rsid} has activity_pct={r['activity_pct']}, expected None"
            )

    def test_missing_snp_status_is_not_assessed(self):
        result = analyse({})
        for rsid, r in result["gene_results"].items():
            assert r["status"] == "not_assessed"

    def test_nmc_is_none_when_no_snps_assessed(self):
        """NMC must be None (not inflated by wildtype assumptions) when panel is empty."""
        result = analyse({})
        assert result["summary"]["net_methylation_capacity"] is None

    def test_coverage_pct_is_zero_when_no_snps(self):
        result = analyse({})
        assert result["summary"]["nmc_coverage_pct"] == 0

    def test_partial_panel_coverage_less_than_100(self):
        path = write_temp(PARTIAL_INPUT)
        genotypes = parse_genotype_file(path)
        result = analyse(genotypes)
        assert result["summary"]["nmc_coverage_pct"] < 100

    def test_full_panel_coverage_is_100(self):
        path = write_temp(WILDTYPE_INPUT)
        genotypes = parse_genotype_file(path)
        result = analyse(genotypes)
        assert result["summary"]["nmc_coverage_pct"] == 100


# ---------------------------------------------------------------------------
# analyse - NMC / BH4 scores with known inputs
# ---------------------------------------------------------------------------

class TestAnalyseScores:
    def test_wildtype_nmc_is_100(self):
        path = write_temp(WILDTYPE_INPUT)
        genotypes = parse_genotype_file(path)
        result = analyse(genotypes)
        assert result["summary"]["net_methylation_capacity"] == 100

    def test_wildtype_bh4_is_100(self):
        path = write_temp(WILDTYPE_INPUT)
        genotypes = parse_genotype_file(path)
        result = analyse(genotypes)
        assert result["summary"]["bh4_axis_capacity"] == 100

    def test_compound_het_detected_in_demo(self):
        path = write_temp(DEMO_INPUT)
        genotypes = parse_genotype_file(path)
        result = analyse(genotypes)
        # Demo has MTHFR TT (hom C677T) + CC (hom A1298C) -> compound het
        assert result["summary"]["mthfr_compound_heterozygous"] is True

    def test_demo_nmc_below_50(self):
        """Demo data has severe compound het; NMC should be substantially below 50."""
        path = write_temp(DEMO_INPUT)
        genotypes = parse_genotype_file(path)
        result = analyse(genotypes)
        nmc = result["summary"]["net_methylation_capacity"]
        assert nmc is not None
        assert nmc < 50

    def test_analyse_output_structure(self):
        path = write_temp(WILDTYPE_INPUT)
        genotypes = parse_genotype_file(path)
        result = analyse(genotypes)
        assert "metadata" in result
        assert "summary" in result
        assert "gene_results" in result
        assert "missing_rsids" in result
        summary_keys = {
            "net_methylation_capacity", "nmc_coverage_pct",
            "bh4_axis_capacity", "mthfr_combined_activity",
            "mthfr_compound_heterozygous",
            "dopamine_synthesis_impact", "serotonin_synthesis_impact",
        }
        assert summary_keys.issubset(result["summary"].keys())


# ---------------------------------------------------------------------------
# Neurotransmitter impact - dopamine and serotonin differ
# ---------------------------------------------------------------------------

class TestNeurotransmitterImpact:
    def test_dopamine_normal_at_high_bh4(self):
        assert _dopamine_impact(90) == "Within Normal Range"

    def test_serotonin_normal_at_high_bh4(self):
        assert _serotonin_impact(90) == "Within Normal Range"

    def test_serotonin_more_sensitive_than_dopamine(self):
        """
        At intermediate BH4 (e.g. 50%), serotonin should show greater impact.
        TPH2 has higher Km for BH4 than TH.
        source: doi:10.1042/BJ20031542 (Fitzpatrick, Biochem J 2004)
        """
        # At BH4=50: dopamine -> Moderately Reduced; serotonin -> Moderately Reduced
        # At BH4=62: dopamine -> Within Normal Range; serotonin -> Moderately Reduced
        assert _dopamine_impact(62) == "Within Normal Range"
        assert _serotonin_impact(62) == "Moderately Reduced"

    def test_severely_reduced_both_at_very_low_bh4(self):
        assert _dopamine_impact(20) == "Severely Reduced"
        assert _serotonin_impact(20) == "Severely Reduced"
