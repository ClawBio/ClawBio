"""
tests/test_methylation_cycle.py
ClawBio · claw-methylation-cycle v0.1.2

Test suite — required per ClawBio TDD mandate.
Covers: parsing, compute_* functions with known I/O pairs, Safety Rule 6,
compound heterozygosity detection, NMC/BH4 scores, and missing SNP handling.
"""

import sys
import tempfile
from pathlib import Path

import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

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


# ---------------------------------------------------------------------------
# Forward-strand (23andMe / AncestryDNA) genotypes for minus-strand genes
# ---------------------------------------------------------------------------

# Same person as a compound heterozygote (C677T + A1298C), written the way a real
# 23andMe file writes it: GRCh37 forward strand. MTHFR, CBS and SHMT1 sit on the
# minus strand, so C677T reads A/G and A1298C reads T/G.
FORWARD_STRAND_INPUT = """\
# rsid\tchromosome\tposition\tgenotype
rs1801133\t1\t11856378\tAG
rs1801131\t1\t11854476\tGT
rs1801394\t5\t7870973\tAG
rs1805087\t1\t237048500\tAA
rs234706\t21\t44485350\tAG
rs3733890\t5\t78421959\tGG
rs1979277\t17\t18232096\tAA
rs4680\t22\t19951271\tAG
"""


class TestForwardStrandGenotypes:
    def _result(self, text):
        return analyse(parse_genotype_file(write_temp(text)))

    def test_mthfr_c677t_forward_het_counted(self):
        assert self._result(FORWARD_STRAND_INPUT)["gene_results"]["rs1801133"]["n_risk_alleles"] == 1

    def test_mthfr_a1298c_forward_het_counted(self):
        assert self._result(FORWARD_STRAND_INPUT)["gene_results"]["rs1801131"]["n_risk_alleles"] == 1

    def test_compound_het_detected_on_forward_strand(self):
        assert self._result(FORWARD_STRAND_INPUT)["summary"]["mthfr_compound_heterozygous"] is True

    def test_forward_strand_compound_het_not_reported_as_normal(self):
        assert self._result(FORWARD_STRAND_INPUT)["summary"]["mthfr_combined_activity"] < 100

    def test_cbs_and_shmt1_forward_strand(self):
        g = self._result(FORWARD_STRAND_INPUT)["gene_results"]
        assert g["rs234706"]["n_risk_alleles"] == 1
        assert g["rs1979277"]["n_risk_alleles"] == 2

    def test_forward_strand_wildtype_is_zero(self):
        text = "rs1801133\t1\t11856378\tGG\nrs1801131\t1\t11854476\tTT\n"
        g = self._result(text)["gene_results"]
        assert g["rs1801133"]["n_risk_alleles"] == 0
        assert g["rs1801131"]["n_risk_alleles"] == 0

    def test_forward_strand_homozygous_c677t(self):
        g = self._result("rs1801133\t1\t11856378\tAA\n")["gene_results"]
        assert g["rs1801133"]["n_risk_alleles"] == 2

    def test_coding_strand_input_still_works(self):
        g = self._result("rs1801133\t1\t11856378\tCT\nrs1801131\t1\t11854476\tAC\n")["gene_results"]
        assert g["rs1801133"]["n_risk_alleles"] == 1
        assert g["rs1801131"]["n_risk_alleles"] == 1

    def test_plus_strand_genes_unchanged(self):
        # MTRR, MTR, BHMT, COMT are on the plus strand: no strand logic applies
        g = self._result(FORWARD_STRAND_INPUT)["gene_results"]
        assert g["rs1801394"]["n_risk_alleles"] == 1
        assert g["rs1805087"]["n_risk_alleles"] == 0
        assert g["rs4680"]["n_risk_alleles"] == 1


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/claw-methylation-cycle/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("claw-methylation-cycle", tmp_path) == snapshot({'report.md': """\
# ClawBio - Methylation Cycle Clinical Report

**Date**: <TIMESTAMP>
**Tool**: claw-methylation-cycle v0.1.3
**Author**: Samuel Carmona Aguirre <samuel@unimed-consulting.es>
**SNPs assessed**: 9/9

> **Research and educational use only (RUO). Not a diagnostic device.**
Consult a qualified clinician before modifying supplementation or treatment.
Enzymatic activity estimates are population-derived approximations, not direct assays.

---
## Executive Summary

| Metric | Value | Status |
|--------|-------|--------|
| Net Methylation Capacity | 58/100 | LOW - Moderately Reduced |
| NMC Panel Coverage | 100% | FULL |
| BH4 Axis Capacity | 13/100 | CRITICAL - Severely Reduced |
| MTHFR Combined Activity | 15% | COMPOUND HET DETECTED |
| Dopamine Synthesis Impact | Severely Reduced | |
| Serotonin Synthesis Impact | Severely Reduced | |

---
## Enzymatic Activity Profile

| Gene | Enzyme | rsID | Genotype | Activity | Status |
|------|--------|------|----------|----------|--------|
| **MTHFR** | C677T | rs1801133 | `CT` | 65% | mildly reduced (65%) |
| **MTRR** | A66G | rs1801394 | `AG` | 80% | mildly reduced (80%) |
| **MTR** | A2756G | rs1805087 | `AG` | 85% | normal (85%) |
| **CBS** | C699T | rs234706 | `CC` | 100% | normal (100%) |
| **BHMT** | R239Q | rs3733890 | `AG` | 70% | mildly reduced (70%) |
| **SHMT1** | C1420T | rs1979277 | `CC` | 100% | normal (100%) |
| **COMT** | Val158Met | rs4680 | `GA` | 65% | mildly reduced (65%) |
| **AHCY** | AHCY | rs819147 | `CC` | 100% | normal (100%) |
| **MTHFR** | Combined (677+1298) | -- | -- | 15% [COMPOUND HET] | severely reduced (15%) |

---
## BH4 / Neurotransmitter Axis

The **MTHFR -> 5-MTHF -> BH4** axis is the critical connection between
methylation cycle variants and neurotransmitter synthesis.
BH4 (tetrahydrobiopterin) is the essential cofactor for tyrosine hydroxylase
(dopamine synthesis) and tryptophan hydroxylase 2 (serotonin synthesis).
TPH2 has a higher Km for BH4 than TH, making serotonin synthesis more
sensitive to BH4 reduction at moderate deficiencies.
source: doi:10.1042/BJ20031542 (Fitzpatrick, Biochem J 2004)

**Estimated BH4 production capacity**: 13% of normal
- Dopamine synthesis pathway: **Severely Reduced**
- Serotonin synthesis pathway: **Severely Reduced**

Some literature reports an association between BH4 deficiency and
ADHD, depression, and anxiety phenotypes
(Ledford et al., 2021, Nutrients 13(3):768;
 Spuch & Agis-Balboa, 2014, SEBBM 179:18-21).
This genotype indicates reduced BH4 production capacity.
A clinician should contextualise this finding with the patient's clinical history.

---
## Compound Heterozygosity

[COMPOUND HET] **MTHFR Compound Heterozygous detected (C677T + A1298C)**

Both MTHFR variants are present simultaneously. This combination reduces
total MTHFR enzymatic activity more than either variant alone.
It is the most clinically significant single-gene methylation finding.
Literature reports active folate (5-MTHF) as preferred over synthetic folic acid
for this genotype; clinician assessment required.
source: doi:10.1086/301927 (van der Put et al., Am J Hum Genet 1998)

---
## Clinical Recommendations

> **FOR CLINICIAN REVIEW ONLY - do not self-administer.**
> The following nutrients are reported in the peer-reviewed literature for the
> pathways indicated. Dosing and indication require individualised clinical
> assessment by a qualified clinician.

- **Genotype finding** -- MTHFR enzymatic reduction impairs folic acid conversion. Literature supports use of 5-MTHF (methylfolate) as the bioavailable form. Ref: Lamers Y et al. (2004) Am J Clin Nutr 80(5):1234-41.
- **Genotype finding** -- MTHFR compound heterozygous (C677T + A1298C): combined activity significantly reduced (15%). Literature reports 5-MTHF + methylcobalamin (B12) for this genotype profile. Ref: Ledford AW et al. (2021) Nutrients 13(3):768.
- **Genotype finding** -- BH4 capacity estimated at 13% of normal. Riboflavin (B2) is reported in the literature as an MTHFR cofactor supporting BH4 regeneration. Vitamin C is reported to maintain BH4 in its reduced (active) form. Dosing and indication require clinician assessment. Ref: McNulty H et al. (2017) Am J Clin Nutr 106(1):128-36.
- **Genotype finding** -- Dopamine and serotonin synthesis pathways may be affected via estimated BH4 reduction. Some literature reports an association between BH4 deficiency and ADHD, depression, and anxiety phenotypes (Ledford et al., 2021). Where clinically relevant, a clinician may evaluate whether neurodevelopmental symptoms correlate with BH4 capacity for potential non-pharmacological support.
- **Genotype finding** -- MTRR A66G variant present. Literature reports methylcobalamin (B12 active form) as preferred over cyanocobalamin for this genotype. Hydroxocobalamin is an alternative. Ref: Olteanu H et al. (2002) Biochemistry 41(45):13378-85.

---
## Missing Variants

All panel SNPs were found in the input file.

---
## References

- Frosst P et al. (1995). A candidate genetic risk factor: common mutation in MTHFR. Nat Genet 10:111-3. doi:10.1038/ng0595-111
- van der Put NM et al. (1998). Second common mutation in MTHFR. Am J Hum Genet 62:1044-51. doi:10.1086/301927
- Gaughan DJ et al. (2001). MTRR A66G polymorphism and neural tube defects. QJM 94:609-17. doi:10.1093/qjmed/94.11.609
- Harmon DL et al. (1999). MTR A2756G polymorphism. Atherosclerosis 146:295-300. doi:10.1016/S0021-9150(99)00113-0
- Lachman HM et al. (1996). COMT pharmacogenetics. Pharmacogenetics 6:243-50. doi:10.1097/00008571-199606000-00007
- Fitzpatrick PF (2004). Tetrahydropterin-dependent amino acid hydroxylases. Biochem J 380:299-310. doi:10.1042/BJ20031542
- Ledford AW et al. (2021). MTHFR and BH4 in neuropsychiatric disorders. Nutrients 13:768. doi:10.3390/nu13030768
- Lamers Y et al. (2004). Supplementation with [6S]-5-methyltetrahydrofolate or folic acid equally reduces plasma total homocysteine. Am J Clin Nutr 80(5):1234-41.
- McNulty H et al. (2017). Riboflavin lowers homocysteine in children and adults with common MTHFR polymorphism. Am J Clin Nutr 106(1):128-36.
- Olteanu H et al. (2002). Differences in the efficiency of reductive methylation of cob(II)alamin. Biochemistry 41(45):13378-85.
- Slow S et al. (2004). Plasma betaine and homocysteine. Clin Chim Acta 340(1-2):57-67.
- Esteller M. (2014). Introduccion a la epigenetica. SEBBM 179:4-6.
- Spuch C, Agis-Balboa RC. (2014). Epigenetica en neurociencias. SEBBM 179:18-21.
- Carmona Aguirre S. (2014/UNESCO 2016). Holomedicina. UNIMED Consulting.
- ClawBio (2026). https://github.com/ClawBio/ClawBio\
""", 'result.json': """\
{
  "metadata": {
    "tool": "claw-methylation-cycle v0.1.3",
    "author": "Samuel Carmona Aguirre <samuel@unimed-consulting.es>",
    "generated_utc": "<TIMESTAMP>",
    "snps_in_panel": 9,
    "snps_found": 9,
    "snps_missing": 0,
    "missing_rsids": []
  },
  "summary": {
    "net_methylation_capacity": 58,
    "nmc_coverage_pct": 100,
    "bh4_axis_capacity": 13,
    "mthfr_combined_activity": 15,
    "mthfr_compound_heterozygous": true,
    "dopamine_synthesis_impact": "Severely Reduced",
    "serotonin_synthesis_impact": "Severely Reduced"
  },
  "gene_results": {
    "rs1801133": {
      "gene": "MTHFR",
      "variant": "C677T",
      "genotype": "CT",
      "n_risk_alleles": 1,
      "status": "heterozygous",
      "activity_pct": 65,
      "effect": "Decreased Folate to 5-MTHF Conversion",
      "weight": 0.35,
      "mthfr_combined_activity": 15,
      "compound_heterozygous": true
    },
    "rs1801131": {
      "gene": "MTHFR",
      "variant": "A1298C",
      "genotype": "GT",
      "n_risk_alleles": 1,
      "status": "heterozygous",
      "activity_pct": 80,
      "effect": "Decreased MTHFR Activity (modifier)",
      "weight": 0.0
    },
    "rs1801394": {
      "gene": "MTRR",
      "variant": "A66G",
      "genotype": "AG",
      "n_risk_alleles": 1,
      "status": "heterozygous",
      "activity_pct": 80,
      "effect": "Decreased Methionine Synthase Reductase Activity",
      "weight": 0.15
    },
    "rs1805087": {
      "gene": "MTR",
      "variant": "A2756G",
      "genotype": "AG",
      "n_risk_alleles": 1,
      "status": "heterozygous",
      "activity_pct": 85,
      "effect": "Decreased Methionine Synthase Activity",
      "weight": 0.1
    },
    "rs234706": {
      "gene": "CBS",
      "variant": "C699T",
      "genotype": "CC",
      "n_risk_alleles": 0,
      "status": "normal",
      "activity_pct": 100,
      "effect": "Increased CBS Activity (diverts homocysteine to transsulfuration)",
      "weight": 0.05
    },
    "rs3733890": {
      "gene": "BHMT",
      "variant": "R239Q",
      "genotype": "AG",
      "n_risk_alleles": 1,
      "status": "heterozygous",
      "activity_pct": 70,
      "effect": "Decreased Betaine-Homocysteine Methyltransferase Activity",
      "weight": 0.15
    },
    "rs1979277": {
      "gene": "SHMT1",
      "variant": "C1420T",
      "genotype": "CC",
      "n_risk_alleles": 0,
      "status": "normal",
      "activity_pct": 100,
      "effect": "Decreased Serine Hydroxymethyltransferase Activity",
      "weight": 0.05
    },
    "rs4680": {
      "gene": "COMT",
      "variant": "Val158Met",
      "genotype": "GA",
      "n_risk_alleles": 1,
      "status": "heterozygous",
      "activity_pct": 65,
      "effect": "Decreased Catechol-O-Methyltransferase Activity",
      "weight": 0.1
    },
    "rs819147": {
      "gene": "AHCY",
      "variant": "AHCY",
      "genotype": "CC",
      "n_risk_alleles": 0,
      "status": "normal",
      "activity_pct": 100,
      "effect": "Decreased Adenosylhomocysteinase Activity",
      "weight": 0.05
    }
  },
  "missing_rsids": []
}\
"""})
