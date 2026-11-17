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
from inline_snapshot import snapshot

from clawbio.common import demo_runner

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


# ------------------------------------------------------------------ #
# Reproducibility bundle
# ------------------------------------------------------------------ #


def test_demo_writes_reproducibility_bundle(tmp_path):
    """--demo writes commands.sh, environment.yml and checksums that verify."""
    import hashlib
    import subprocess

    out = tmp_path / "out"
    script = Path(__file__).resolve().parent.parent / "genome_compare.py"
    subprocess.run(
        [sys.executable, str(script), "--demo", "--no-figures", "--output", str(out)],
        check=True, capture_output=True,
    )
    repro = out / "reproducibility"
    commands = (repro / "commands.sh").read_text()
    assert "python skills/genome-compare/genome_compare.py --demo" in commands
    assert "numpy" in (repro / "environment.yml").read_text()

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = [line.split("  ", 1)[1] for line in lines]
    assert {"report.md", "result.json", "reproducibility/environment.yml"} <= set(labels)
    for line in lines:
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((out / label).read_bytes()).hexdigest() == digest


# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `pytest skills/genome-compare/tests --inline-snapshot=fix`; it rewrites the
# snapshot() values below, and the PR diff shows which numbers moved.

@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("genome-compare", tmp_path) == snapshot({'report.md': """\
# Genome Comparison Report

**Date**: <TIMESTAMP> UTC  \n\
**Reference**: George Church (PGP-1, hu43860C)  \n\
**Input**: Manuel Corpas ([PGP-UK uk6D0CFA](https://my.personalgenomes.org.uk/profile/uk6D0CFA), Corpasome, demo mode)  \n\

## Summary

| Metric | Value |
|--------|-------|
| **IBS Score** | 0.7418 |
| **SNPs Overlapping** | 544,515 |
| **Genotypes Identical (IBS=2)** | 300,474 (55.2%) |
| **Top Ancestry** | EUR (54.9%) |
| **AIMs Used** | 47 / 65 |

## Identity By State Analysis

Your IBS score (0.7418) falls between **unrelated, same population** (0.735) and **first cousins** (0.760).

### Population Context

| Relationship | Typical IBS |
|---|---|
| Identical twins / same person | 1.000 |
| Parent-child | 0.850 |
| Full siblings | 0.820 |
| Half-siblings | 0.780 |
| First cousins | 0.760 |
| Unrelated, same population | 0.735 **<--** |
| European vs South Asian | 0.720 |
| European vs East Asian | 0.710 |
| European vs African | 0.695 |
| African vs East Asian | 0.680 |

![IBS Context](figures/ibs_context.png)

## Chromosome Breakdown

| Chr | SNPs | IBS Score | Identical (IBS=2) |
|-----|------|-----------|-------------------|
| 1 | 41,030 | 0.7454 | 22,578 |
| 2 | 43,734 | 0.7362 | 23,662 |
| 3 | 36,549 | 0.7371 | 19,639 |
| 4 | 32,263 | 0.7376 | 17,414 |
| 5 | 33,119 | 0.7454 | 18,289 |
| 6 | 35,792 | 0.7489 | 19,941 |
| 7 | 29,142 | 0.7417 | 15,819 |
| 8 | 30,277 | 0.7457 | 16,790 |
| 9 | 25,669 | 0.7411 | 14,095 |
| 10 | 28,313 | 0.7447 | 15,482 |
| 11 | 26,444 | 0.7408 | 14,497 |
| 12 | 26,306 | 0.7435 | 14,260 |
| 13 | 20,151 | 0.7336 | 10,796 |
| 14 | 17,753 | 0.7461 | 9,825 |
| 15 | 16,114 | 0.7318 | 8,654 |
| 16 | 16,421 | 0.7473 | 9,106 |
| 17 | 14,329 | 0.7495 | 8,003 |
| 18 | 16,208 | 0.7412 | 8,854 |
| 19 | 9,454 | 0.7406 | 5,089 |
| 20 | 13,726 | 0.7398 | 7,440 |
| 21 | 7,998 | 0.7336 | 4,256 |
| 22 | 8,207 | 0.7336 | 4,406 |
| X | 13,017 | 0.7000 | 9,108 |
| Y | 1,136 | 0.9974 | 1,133 |
| MT | 1,363 | 0.9817 | 1,338 |

![Chromosome IBS](figures/chromosome_ibs.png)

## Ancestry Composition

### Your Estimated Ancestry

| Population | Proportion |
|------------|-----------|
| EUR | 54.9% |
| AFR | 28.8% |
| EAS | 14.4% |
| SAS | 1.7% |
| AMR | 0.1% |

*Method: EM admixture (binomial genotype likelihood)*

![Ancestry Pie](figures/ancestry_pie.png)

### Ground Truth — Manuel Corpas (23andMe Ancestry Report)

*Since this demo uses Manuel Corpas as the input, we can compare the estimation above against his actual 23andMe ancestry report:*

| Region | Proportion |
|--------|-----------|
| **European** | **97.7%** |
| &emsp;Southern European | 92.0% |
| &emsp;&emsp;Spanish & Portuguese | 79.8% |
| &emsp;&emsp;Italian | 6.6% |
| &emsp;&emsp;Broadly Southern European | 5.6% |
| &emsp;Northwestern European | 3.0% |
| &emsp;&emsp;French & German | 0.6% |
| &emsp;&emsp;Broadly Northwestern European | 2.4% |
| &emsp;Ashkenazi Jewish | 0.5% |
| &emsp;Broadly European | 2.4% |
| **Western Asian & North African** | **0.8%** |
| &emsp;North African & Arabian | 0.3% |
| &emsp;Broadly Western Asian & North African | 0.5% |
| **Sub-Saharan African** | **0.3%** |
| &emsp;West African | 0.2% |
| &emsp;&emsp;Senegambian & Guinean | 0.2% |
| &emsp;Northern East African | 0.1% |
| &emsp;&emsp;Sudanese | 0.1% |
| **East Asian & Native American** | **0.1%** |
| **Unassigned** | **1.0%** |

**Maternal haplogroup**: W  \n\
**Paternal haplogroup**: R-P312  \n\

**Neanderthal**: More Neanderthal variants than 61% of customers  \n\

![Ancestry Comparison](figures/ancestry_comparison.png)

## About the Reference Genome

**George Church** (hu43860C) is PGP-1 — the first participant in the [Personal Genome Project](https://pgp.med.harvard.edu/). He is a professor of genetics at Harvard Medical School and a pioneer of genomics. His 23andMe data (569,226 SNPs, build 36) is CC0 public domain.

## Methods

**IBS (Identity By State)**: For each SNP present in both files, alleles are compared. IBS per site = number of shared alleles (0, 1, or 2). Overall IBS = sum(per-site IBS) / (2 * N). Range: 0.0 (no alleles shared) to 1.0 (identical).

**Ancestry estimation**: EM admixture algorithm using ancestry-informative markers (AIMs) with population allele frequencies from 1000 Genomes Phase 3. Binomial genotype likelihood model with iterative expectation-maximisation. This is an approximation — for clinical-grade ancestry, use ADMIXTURE or RFMix.

**Reference IBS values** are approximate ranges from population genetics literature (1000 Genomes, HGDP). Actual values depend on SNP panel and population composition.

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*

## Reproducibility

- Input checksum: `4dca2c9863d9b9cb`
- Reference checksum: `538d32b5dc9de76d`
- Generated: <TIMESTAMP> UTC
""", 'result.json': """\
{
  "skill": "genome-compare",
  "version": "0.2.0",
  "completed_at": "<TIMESTAMP>",
  "input_checksum": "sha256:4dca2c9863d9b9cb",
  "datasets": {},
  "summary": {
    "ibs_score": 0.7418060108536955,
    "n_overlap": 544515,
    "n_concordant": 300474,
    "top_ancestry": "EUR"
  },
  "data": {
    "ibs_score": 0.7418060108536955,
    "per_chrom": {
      "1": {
        "ibs": 0.7454179868388984,
        "n_overlap": 41030,
        "n_concordant": 22578
      },
      "2": {
        "ibs": 0.7361892349201994,
        "n_overlap": 43734,
        "n_concordant": 23662
      },
      "3": {
        "ibs": 0.7370653095843935,
        "n_overlap": 36549,
        "n_concordant": 19639
      },
      "4": {
        "ibs": 0.7375631528376159,
        "n_overlap": 32263,
        "n_concordant": 17414
      },
      "5": {
        "ibs": 0.7454029409100517,
        "n_overlap": 33119,
        "n_concordant": 18289
      },
      "6": {
        "ibs": 0.7489383102369245,
        "n_overlap": 35792,
        "n_concordant": 19941
      },
      "7": {
        "ibs": 0.7416786768238282,
        "n_overlap": 29142,
        "n_concordant": 15819
      },
      "8": {
        "ibs": 0.7456650262575553,
        "n_overlap": 30277,
        "n_concordant": 16790
      },
      "9": {
        "ibs": 0.741068993727843,
        "n_overlap": 25669,
        "n_concordant": 14095
      },
      "10": {
        "ibs": 0.7447109101826016,
        "n_overlap": 28313,
        "n_concordant": 15482
      },
      "11": {
        "ibs": 0.7407540462864922,
        "n_overlap": 26444,
        "n_concordant": 14497
      },
      "12": {
        "ibs": 0.7434995818444461,
        "n_overlap": 26306,
        "n_concordant": 14260
      },
      "13": {
        "ibs": 0.7335616098456652,
        "n_overlap": 20151,
        "n_concordant": 10796
      },
      "14": {
        "ibs": 0.7460992508308455,
        "n_overlap": 17753,
        "n_concordant": 9825
      },
      "15": {
        "ibs": 0.7317860245749038,
        "n_overlap": 16114,
        "n_concordant": 8654
      },
      "16": {
        "ibs": 0.7472748310090738,
        "n_overlap": 16421,
        "n_concordant": 9106
      },
      "17": {
        "ibs": 0.7495289273501291,
        "n_overlap": 14329,
        "n_concordant": 8003
      },
      "18": {
        "ibs": 0.7412080454096742,
        "n_overlap": 16208,
        "n_concordant": 8854
      },
      "19": {
        "ibs": 0.7405859953458853,
        "n_overlap": 9454,
        "n_concordant": 5089
      },
      "20": {
        "ibs": 0.7398003788430716,
        "n_overlap": 13726,
        "n_concordant": 7440
      },
      "21": {
        "ibs": 0.7335583895973994,
        "n_overlap": 7998,
        "n_concordant": 4256
      },
      "22": {
        "ibs": 0.7336420129158036,
        "n_overlap": 8207,
        "n_concordant": 4406
      },
      "X": {
        "ibs": 0.6999692709533687,
        "n_overlap": 13017,
        "n_concordant": 9108
      },
      "Y": {
        "ibs": 0.9973591549295775,
        "n_overlap": 1136,
        "n_concordant": 1133
      },
      "MT": {
        "ibs": 0.9816581071166545,
        "n_overlap": 1363,
        "n_concordant": 1338
      }
    },
    "ancestry": {
      "continental": {
        "AFR": 0.2883,
        "EUR": 0.5493,
        "EAS": 0.1442,
        "SAS": 0.0173,
        "AMR": 0.0009
      },
      "n_aims_used": 47,
      "n_aims_total": 65,
      "method": "EM admixture (binomial genotype likelihood)"
    }
  }
}\
"""})
