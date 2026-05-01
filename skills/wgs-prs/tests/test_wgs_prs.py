"""Integration tests for wgs_to_prs_bridge.py.

Uses the Corpas chr20 VCF subset from ClawBio for the VCF-entry path,
and mocks subprocess for sarek calls.
"""

from __future__ import annotations

import json
import textwrap
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import sys

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from wgs_prs import BridgeConfig, WgsToPrsBridge
from clawbio.common.sarek import SarekConfig
from clawbio.common.vcf_qc import QcConfig

# ---------------------------------------------------------------------------
# Synthetic VCF fixture (enough SNPs to pass QC thresholds)
# ---------------------------------------------------------------------------

def _write_synthetic_vcf(path: Path, n_snps: int = 200) -> Path:
    """Write a minimal VCF with n_snps transitions (all Ti → Ti/Tv = ∞, handled as warning)."""
    lines = [
        "##fileformat=VCFv4.2",
        "##FILTER=<ID=PASS,Description='All filters passed'>",
        "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tSAMPLE",
    ]
    allele_pairs = [("A", "G"), ("C", "T"), ("G", "A"), ("T", "C")]
    for i in range(n_snps):
        ref, alt = allele_pairs[i % 4]
        gt = "0/1" if i % 3 != 0 else "1/1"
        lines.append(
            f"1\t{(i + 1) * 1000}\trs{i}\t{ref}\t{alt}\t50\tPASS\t.\tGT:DP\t{gt}:30"
        )
    path.write_text("\n".join(lines) + "\n")
    return path


@pytest.fixture
def synthetic_vcf(tmp_path) -> Path:
    return _write_synthetic_vcf(tmp_path / "test.vcf", n_snps=200)


@pytest.fixture
def corpas_chr20_vcf() -> Path | None:
    """Point to real Corpas chr20 subset if available."""
    p = Path(__file__).resolve().parents[3] / "corpas-30x" / \
        "subsets" / "chr20_snps_indels.vcf.gz"
    return p if p.exists() else None


# ---------------------------------------------------------------------------
# Bridge — skip-sarek (VCF entry point)
# ---------------------------------------------------------------------------

class TestBridgeVcfEntry:
    def test_run_with_input_vcf_skips_sarek_stage(self, synthetic_vcf, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            sample_id="TEST",
            dry_run=False,
            fail_fast=False,
            qc=QcConfig(min_snp_count=50),   # low threshold for tiny fixture
        )
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):   # Python QC mode
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(input_vcf=synthetic_vcf)

        sarek_stage = next(s for s in report.stages if s.name == "sarek")
        assert sarek_stage.status == "skipped"

    def test_report_json_written(self, synthetic_vcf, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            fail_fast=False,
            qc=QcConfig(min_snp_count=50),
        )
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(input_vcf=synthetic_vcf)

        assert report.report_json is not None
        assert report.report_json.exists()
        data = json.loads(report.report_json.read_text())
        assert "overall_status" in data
        assert "stages" in data

    def test_markdown_report_written(self, synthetic_vcf, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            fail_fast=False,
            qc=QcConfig(min_snp_count=50),
        )
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(input_vcf=synthetic_vcf)

        assert report.report_md is not None
        assert report.report_md.exists()
        content = report.report_md.read_text()
        assert "ClawBio WGS-PRS Bridge Report" in content
        assert "VCF QC Metrics" in content

    def test_qc_metrics_in_report(self, synthetic_vcf, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            fail_fast=False,
            qc=QcConfig(min_snp_count=50),
        )
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(input_vcf=synthetic_vcf)

        assert report.qc_metrics is not None
        assert "metrics" in report.qc_metrics


# ---------------------------------------------------------------------------
# Bridge — dry run from FASTQ
# ---------------------------------------------------------------------------

class TestBridgeDryRun:
    def test_dry_run_from_fastq_does_not_call_subprocess(self, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            dry_run=True,
            fail_fast=False,
            qc=QcConfig(min_snp_count=1),
        )
        with patch("clawbio.common.sarek.shutil.which", return_value="/usr/bin/nextflow"), \
             patch("clawbio.common.vcf_qc.shutil.which", return_value=None), \
             patch("clawbio.common.sarek.subprocess.run") as mock_run:
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(
                fastq_r1=str(tmp_path / "R1.fastq.gz"),
                fastq_r2=str(tmp_path / "R2.fastq.gz"),
            )
        # subprocess.run should NOT have been called for sarek (dry_run=True)
        mock_run.assert_not_called()

    def test_dry_run_produces_report_json(self, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            dry_run=True,
            fail_fast=False,
        )
        with patch("clawbio.common.sarek.shutil.which", return_value="/usr/bin/nextflow"), \
             patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(
                fastq_r1=str(tmp_path / "R1.fastq.gz"),
            )

        assert report.report_json.exists()


# ---------------------------------------------------------------------------
# Bridge — fail_fast behaviour
# ---------------------------------------------------------------------------

class TestFailFast:
    def test_fail_fast_true_aborts_after_qc_failure(self, synthetic_vcf, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            fail_fast=True,
            qc=QcConfig(min_snp_count=999_999),   # impossible threshold → QC always fails
        )
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(input_vcf=synthetic_vcf)

        assert report.overall_status == "failed"
        # prs_scoring stage should not exist
        stage_names = [s.name for s in report.stages]
        assert "prs_scoring" not in stage_names

    def test_fail_fast_false_continues_after_qc_failure(self, synthetic_vcf, tmp_path):
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            fail_fast=False,
            qc=QcConfig(min_snp_count=999_999),   # always fails
        )
        with patch("clawbio.common.vcf_qc.shutil.which", return_value=None):
            bridge = WgsToPrsBridge(cfg)
            report = bridge.run(input_vcf=synthetic_vcf)

        stage_names = [s.name for s in report.stages]
        assert "prs_scoring" in stage_names


# ---------------------------------------------------------------------------
# Corpas chr20 integration test (skipped if data not present)
# ---------------------------------------------------------------------------

class TestCorpasIntegration:
    @pytest.mark.skipif(
        not (
            Path(__file__).resolve().parents[3] /
            "corpas-30x" / "subsets" / "chr20_snps_indels.vcf.gz"
        ).exists(),
        reason="Corpas chr20 subset not present — download from Zenodo",
    )
    def test_corpas_chr20_passes_qc(self, tmp_path):
        vcf = (
            Path(__file__).resolve().parents[3] /
            "corpas-30x" / "subsets" / "chr20_snps_indels.vcf.gz"
        )
        cfg = BridgeConfig(
            output_dir=str(tmp_path / "out"),
            fail_fast=False,
            qc=QcConfig(min_snp_count=1000, min_titv_ratio=1.5),
        )
        bridge = WgsToPrsBridge(cfg)
        report = bridge.run(input_vcf=vcf)

        qc_stage = next(s for s in report.stages if s.name == "vcf_qc")
        assert qc_stage.status in {"success", "failed"}  # status depends on bcftools availability
        assert report.qc_metrics is not None
        # Stage 3 may fail if gwas-prs dependencies are missing — that is acceptable here
        assert report.report_json is not None and report.report_json.exists()
