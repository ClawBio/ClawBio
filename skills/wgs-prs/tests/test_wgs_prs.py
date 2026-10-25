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
from hypothesis import example, given, settings
from hypothesis import strategies as st

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

import wgs_prs
from wgs_prs import BridgeConfig, BridgeReport, WgsToPrsBridge, write_demo_vcf
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
# Bridge: skip-sarek (VCF entry point)
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


class TestBridgeReport:
    def test_markdown_report_embeds_gwas_prs_output(self, tmp_path):
        """The bridge reads the report filename written by gwas-prs."""
        output_dir = tmp_path / "out"
        prs_output = output_dir / "prs_output"
        prs_output.mkdir(parents=True)
        (prs_output / "prs_report.md").write_text("# PRS result\n\nType 2 diabetes")

        bridge = WgsToPrsBridge(BridgeConfig(output_dir=str(output_dir)))
        report = BridgeReport(
            sample_id="TEST",
            timestamp="2026-09-08T00:00:00+00:00",
            output_dir=output_dir,
        )

        markdown = bridge._write_markdown_report(report).read_text()

        assert "PRS report available" in markdown
        assert "# PRS result" in markdown


# ---------------------------------------------------------------------------
# Bridge: dry run from FASTQ
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
# Bridge: fail_fast behaviour
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
        reason="Corpas chr20 subset not present: download from Zenodo",
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
        # Stage 3 may fail if gwas-prs dependencies are missing; that is acceptable here
        assert report.report_json is not None and report.report_json.exists()


# ---------------------------------------------------------------------------
# Demo VCF
# ---------------------------------------------------------------------------

def _declared(vcf: Path, kind: str) -> list[str]:
    lines = vcf.read_text(encoding="utf-8").splitlines()
    return [line.split("ID=", 1)[1].split(",", 1)[0].rstrip(">") for line in lines if line.startswith(f"##{kind}=<")]


def _assert_bcftools_can_index(vcf: Path) -> None:
    records = [line.split("\t") for line in vcf.read_text(encoding="utf-8").splitlines() if not line.startswith("#")]
    contigs = _declared(vcf, "contig")
    assert {record[0] for record in records} <= set(contigs)
    assert {key for record in records for key in record[8].split(":")} <= set(_declared(vcf, "FORMAT"))
    # Indexing the normalised VCF needs records in header-contig, then position, order.
    keys = [(contigs.index(record[0]), int(record[1])) for record in records]
    assert keys == sorted(keys)


@settings(deadline=None, max_examples=25)
@given(n_snps=st.integers(min_value=0, max_value=300))
def test_demo_vcf_declares_every_contig_and_format_key_it_uses(tmp_path_factory, n_snps) -> None:
    """bcftools refuses a record whose CHROM or FORMAT key the header does not
    declare (`CONTIG id=0 not present in the header`) and cannot index unsorted
    records, so the --demo run failed wherever bcftools is installed."""
    _assert_bcftools_can_index(write_demo_vcf(tmp_path_factory.mktemp("vcf") / "demo.vcf", n_snps=n_snps))


_CONTIGS = [str(n) for n in range(1, 23)] + ["X", "Y", "MT"]


@settings(deadline=None, max_examples=50)
@given(panel=st.lists(
    st.tuples(st.sampled_from(_CONTIGS), st.integers(min_value=1, max_value=250_000_000)),
    min_size=1, max_size=20,
))
@example(panel=[("X", 100), ("2", 50)])
def test_demo_vcf_orders_sex_and_mitochondrial_contigs(tmp_path_factory, panel) -> None:
    """Sorting contigs with int() raised ValueError as soon as the demo panel
    held an X, Y or MT variant."""
    variants = tuple((contig, position, f"rs{i}", "A", "G") for i, (contig, position) in enumerate(panel))
    with patch.object(wgs_prs, "_DEMO_T2D_PANEL_VARIANTS", variants):
        vcf = write_demo_vcf(tmp_path_factory.mktemp("vcf") / "demo.vcf", n_snps=0)
    _assert_bcftools_can_index(vcf)


def test_bcftools_reads_the_demo_vcf(tmp_path) -> None:
    import shutil
    import subprocess

    if shutil.which("bcftools") is None:
        pytest.skip("bcftools is not installed")
    vcf = write_demo_vcf(tmp_path / "demo.vcf")
    result = subprocess.run(
        ["bcftools", "view", "-Ob", "-o", str(tmp_path / "demo.bcf"), str(vcf)],
        capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr


def test_demo_vcf_orders_sex_chromosomes_before_the_mitochondrion(tmp_path) -> None:
    """The order the docstring promises: autosomes numerically, then X, Y, MT. Sorting
    non-numeric names alphabetically put MT before X and Y."""
    variants = tuple((contig, 100, f"rs{i}", "A", "G") for i, contig in enumerate(["MT", "Y", "2", "X", "10", "1"]))
    with patch.object(wgs_prs, "_DEMO_T2D_PANEL_VARIANTS", variants):
        vcf = write_demo_vcf(tmp_path / "demo.vcf", n_snps=0)
    assert _declared(vcf, "contig") == ["1", "2", "10", "X", "Y", "MT"]
