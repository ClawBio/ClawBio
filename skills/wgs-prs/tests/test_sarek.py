"""Tests for sarek_wrapper.py. Mocks subprocess to avoid needing Nextflow."""

from __future__ import annotations

import csv
import json
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
import sys

from clawbio.common.sarek import SarekConfig, SarekWrapper, build_samplesheet


# ---------------------------------------------------------------------------
# build_samplesheet
# ---------------------------------------------------------------------------

class TestBuildSamplesheet:
    def test_writes_csv_with_correct_columns(self, tmp_path):
        ss = build_samplesheet(
            fastq_r1="/data/sample_R1.fastq.gz",
            fastq_r2="/data/sample_R2.fastq.gz",
            output_path=tmp_path / "samplesheet.csv",
            sample_id="HG001",
            sex="XX",
        )
        assert ss.exists()
        with open(ss) as fh:
            reader = csv.DictReader(fh)
            rows = list(reader)
        assert len(rows) == 1
        row = rows[0]
        assert row["patient"] == "HG001"
        assert row["sample"] == "HG001"
        assert row["sex"] == "XX"
        assert row["fastq_1"].endswith("sample_R1.fastq.gz")
        assert row["fastq_2"].endswith("sample_R2.fastq.gz")

    def test_single_end_fastq_r2_empty(self, tmp_path):
        ss = build_samplesheet(
            fastq_r1="/data/single_R1.fastq.gz",
            fastq_r2=None,
            output_path=tmp_path / "ss.csv",
        )
        with open(ss) as fh:
            reader = csv.DictReader(fh)
            row = list(reader)[0]
        assert row["fastq_2"] == ""

    def test_creates_parent_directories(self, tmp_path):
        nested = tmp_path / "a" / "b" / "ss.csv"
        build_samplesheet("/data/R1.fastq.gz", None, nested)
        assert nested.exists()


# ---------------------------------------------------------------------------
# SarekConfig
# ---------------------------------------------------------------------------

class TestSarekConfig:
    def test_to_nextflow_params_includes_genome_and_tools(self):
        cfg = SarekConfig(genome="GATK.GRCh38", tools=["haplotypecaller"])
        params = cfg.to_nextflow_params()
        assert "--genome" in params
        assert "GATK.GRCh38" in params
        assert "--tools" in params

    def test_skip_bqsr_adds_skip_tools_flag(self):
        cfg = SarekConfig(skip_bqsr=True)
        params = cfg.to_nextflow_params()
        assert "--skip_tools" in params
        assert "baserecalibrator" in params

    def test_joint_germline_flag(self):
        cfg = SarekConfig(joint_germline=True)
        params = cfg.to_nextflow_params()
        assert "--joint_germline" in params


# ---------------------------------------------------------------------------
# SarekWrapper.run: dry run
# ---------------------------------------------------------------------------

class TestSarekWrapperDryRun:
    def test_dry_run_returns_mock_vcf_path_without_subprocess(self, tmp_path):
        cfg = SarekConfig(
            dry_run=True,
            sample_id="HG001",
            output_dir=str(tmp_path / "sarek_out"),
        )
        wrapper = SarekWrapper(cfg)
        # Patch nextflow check so it doesn't fail even if nextflow absent
        with patch("clawbio.common.sarek.shutil.which", return_value="/usr/bin/nextflow"):
            vcf = wrapper.run(
                fastq_r1=str(tmp_path / "R1.fastq.gz"),
                fastq_r2=str(tmp_path / "R2.fastq.gz"),
            )
        assert "HG001" in str(vcf)
        assert vcf.suffix in {".gz", ".vcf"}

    def test_dry_run_writes_samplesheet(self, tmp_path):
        cfg = SarekConfig(
            dry_run=True,
            sample_id="TEST",
            output_dir=str(tmp_path / "out"),
        )
        wrapper = SarekWrapper(cfg)
        with patch("clawbio.common.sarek.shutil.which", return_value="/usr/bin/nextflow"):
            wrapper.run(fastq_r1="/fake/R1.fastq.gz")
        ss = tmp_path / "out" / "samplesheet.csv"
        assert ss.exists()


# ---------------------------------------------------------------------------
# SarekWrapper.run: missing nextflow
# ---------------------------------------------------------------------------

class TestSarekWrapperMissingNextflow:
    def test_raises_runtime_error_if_nextflow_absent(self, tmp_path):
        cfg = SarekConfig(output_dir=str(tmp_path), dry_run=False)
        wrapper = SarekWrapper(cfg)
        with patch("clawbio.common.sarek.shutil.which", return_value=None):
            with pytest.raises(RuntimeError, match="Nextflow not found"):
                wrapper.run(fastq_r1="/fake/R1.fastq.gz")


# ---------------------------------------------------------------------------
# SarekWrapper.check_environment
# ---------------------------------------------------------------------------

class TestCheckEnvironment:
    def test_returns_dict_with_expected_keys(self):
        wrapper = SarekWrapper()
        with patch("clawbio.common.sarek.shutil.which", return_value=None):
            status = wrapper.check_environment()
        assert "nextflow_found" in status
        assert "docker_found" in status
        assert "singularity_found" in status

    def test_nextflow_found_true_when_present(self):
        wrapper = SarekWrapper()
        with patch("clawbio.common.sarek.shutil.which", return_value="/usr/local/bin/nextflow"):
            with patch("clawbio.common.sarek.subprocess.run") as mock_run:
                mock_run.return_value = MagicMock(
                    stdout="Nextflow version 23.10.0", returncode=0
                )
                status = wrapper.check_environment()
        assert status["nextflow_found"] is True


# ---------------------------------------------------------------------------
# SarekWrapper.write_run_manifest
# ---------------------------------------------------------------------------

class TestRunManifest:
    def test_writes_valid_json(self, tmp_path):
        cfg = SarekConfig(genome="GATK.GRCh38", sarek_version="3.4.4")
        wrapper = SarekWrapper(cfg)
        out = tmp_path / "manifest.json"
        wrapper.write_run_manifest(out)
        data = json.loads(out.read_text())
        assert data["genome"] == "GATK.GRCh38"
        assert data["version"] == "3.4.4"
        assert "tools" in data
