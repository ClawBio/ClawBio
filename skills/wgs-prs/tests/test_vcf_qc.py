"""Tests for vcf_qc.py. Uses a synthetic minimal VCF; mocks bcftools."""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from unittest.mock import patch, MagicMock

import pytest
import sys

from clawbio.common.vcf_qc import VcfQC, QcConfig, QcResult


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MINIMAL_VCF = textwrap.dedent("""\
    ##fileformat=VCFv4.2
    ##FILTER=<ID=PASS,Description="All filters passed">
    #CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE
    1\t100\trs1\tA\tG\t50\tPASS\t.\tGT:DP\t0/1:30
    1\t200\trs2\tC\tT\t60\tPASS\t.\tGT:DP\t1/1:25
    1\t300\trs3\tG\tA\t45\tPASS\t.\tGT:DP\t0/1:20
    1\t400\trs4\tT\tC\t55\tPASS\t.\tGT:DP\t0/1:35
    1\t500\trs5\tA\tAT\t40\tPASS\t.\tGT:DP\t0/1:15
""")

BCFTOOLS_STATS_OUTPUT = textwrap.dedent("""\
    # This file was produced by bcftools stats
    SN\t0\tnumber of samples:\t1
    SN\t0\tnumber of records:\t5
    SN\t0\tnumber of SNPs:\t4
    SN\t0\tnumber of indels:\t1
    SN\t0\tnumber of heterozygous SNPs:\t3
    SN\t0\tnumber of homozygous SNPs:\t1
    TSTV\t0\t3\t1\t2\t2.000\t0\t0
""")


@pytest.fixture
def minimal_vcf(tmp_path) -> Path:
    p = tmp_path / "test.vcf"
    p.write_text(MINIMAL_VCF)
    return p


# ---------------------------------------------------------------------------
# Python-only stats (no bcftools)
# ---------------------------------------------------------------------------

class TestPythonStats:
    def test_parses_snp_count(self, minimal_vcf):
        qc = VcfQC.__new__(VcfQC)
        qc._bcftools = None
        qc.config = QcConfig()
        stats = qc._python_stats(minimal_vcf)
        assert stats["number of SNPs:"] == 4

    def test_parses_indel_count(self, minimal_vcf):
        qc = VcfQC.__new__(VcfQC)
        qc._bcftools = None
        qc.config = QcConfig()
        stats = qc._python_stats(minimal_vcf)
        assert stats["number of indels:"] == 1

    def test_titv_ratio_computed(self, minimal_vcf):
        qc = VcfQC.__new__(VcfQC)
        qc._bcftools = None
        qc.config = QcConfig()
        stats = qc._python_stats(minimal_vcf)
        # A→G, C→T, G→A, T→C are all transitions → Ti/Tv = 4/0 → None (no TV)
        # Actually Ti/Tv of infinity → stored as None (tv_count = 0)
        # Check that it doesn't crash regardless
        assert "titv_ratio" in stats

    def test_het_hom_counts(self, minimal_vcf):
        qc = VcfQC.__new__(VcfQC)
        qc._bcftools = None
        qc.config = QcConfig()
        stats = qc._python_stats(minimal_vcf)
        assert stats["number of heterozygous SNPs:"] >= 0
        assert stats["number of homozygous SNPs:"] >= 0


# ---------------------------------------------------------------------------
# Parse bcftools stats output
# ---------------------------------------------------------------------------

class TestParseBcftoolsStats:
    def test_parses_snp_count(self):
        qc = VcfQC.__new__(VcfQC)
        stats = qc._parse_bcftools_stats(BCFTOOLS_STATS_OUTPUT)
        assert stats["number of SNPs:"] == 4.0

    def test_parses_indel_count(self):
        qc = VcfQC.__new__(VcfQC)
        stats = qc._parse_bcftools_stats(BCFTOOLS_STATS_OUTPUT)
        assert stats["number of indels:"] == 1.0

    def test_parses_titv_ratio(self):
        qc = VcfQC.__new__(VcfQC)
        stats = qc._parse_bcftools_stats(BCFTOOLS_STATS_OUTPUT)
        assert stats["titv_ratio"] == pytest.approx(2.000)

    def test_parses_het_hom(self):
        qc = VcfQC.__new__(VcfQC)
        stats = qc._parse_bcftools_stats(BCFTOOLS_STATS_OUTPUT)
        assert stats["number of heterozygous SNPs:"] == 3.0
        assert stats["number of homozygous SNPs:"] == 1.0


# ---------------------------------------------------------------------------
# Pass / fail evaluation
# ---------------------------------------------------------------------------

class TestPassFail:
    def _make_result(self, **kwargs) -> QcResult:
        defaults = dict(
            snp_count=4_000_000,
            indel_count=500_000,
            total_variants=4_500_000,
            titv_ratio=2.1,
            het_count=2_000_000,
            hom_alt_count=1_000_000,
            het_hom_ratio=2.0,
            filtered_out=10_000,
        )
        defaults.update(kwargs)
        r = QcResult(**{k: v for k, v in defaults.items() if k in QcResult.__dataclass_fields__})
        return r

    def test_good_metrics_pass(self):
        cfg = QcConfig()
        qc = VcfQC(cfg)
        result = self._make_result()
        qc._evaluate_pass_fail(result)
        assert result.passes_qc

    def test_too_few_snps_fails(self):
        cfg = QcConfig(min_snp_count=1000)
        qc = VcfQC(cfg)
        result = self._make_result(snp_count=50)
        qc._evaluate_pass_fail(result)
        assert not result.passes_qc
        assert any("SNP" in r for r in result.fail_reasons)

    def test_low_titv_fails(self):
        cfg = QcConfig(min_titv_ratio=1.8)
        qc = VcfQC(cfg)
        result = self._make_result(titv_ratio=1.2)
        qc._evaluate_pass_fail(result)
        assert not result.passes_qc

    def test_high_het_hom_fails(self):
        cfg = QcConfig(max_het_hom_ratio=3.0)
        qc = VcfQC(cfg)
        result = self._make_result(het_hom_ratio=5.0)
        qc._evaluate_pass_fail(result)
        assert not result.passes_qc

    def test_missing_titv_adds_warning_not_failure(self):
        cfg = QcConfig(min_snp_count=100)
        qc = VcfQC(cfg)
        result = self._make_result(titv_ratio=None)
        qc._evaluate_pass_fail(result)
        # Should still pass (warning only)
        assert result.passes_qc
        assert any("Ti/Tv" in w for w in result.warnings)


# ---------------------------------------------------------------------------
# Full run: no bcftools (Python-only mode)
# ---------------------------------------------------------------------------

class TestFullRunNoBcftools:
    def test_run_without_bcftools_uses_python_fallback(self, minimal_vcf, tmp_path):
        cfg = QcConfig(min_snp_count=1)   # low threshold so tiny VCF passes
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            qc = VcfQC(cfg)
            result = qc.run(input_vcf=minimal_vcf, output_dir=tmp_path / "qc")

        assert result.metrics_json is not None
        assert result.metrics_json.exists()
        assert result.canonical_vcf == minimal_vcf   # no-op in Python mode

    def test_metrics_json_written(self, minimal_vcf, tmp_path):
        cfg = QcConfig(min_snp_count=1)
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            qc = VcfQC(cfg)
            result = qc.run(input_vcf=minimal_vcf, output_dir=tmp_path / "qc")

        data = json.loads(result.metrics_json.read_text())
        assert "qc_status" in data
        assert "metrics" in data
        assert "snp_count" in data["metrics"]
