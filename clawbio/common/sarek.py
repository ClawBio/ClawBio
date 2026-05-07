#!/usr/bin/env python3
"""
sarek_wrapper.py: nf-core/sarek orchestration layer for ClawBio WGS-PRS bridge
ClawBio WGS-PRS Bridge v0.1.0
Author: David de Lorenzo
License: MIT

Generates a sarek-compatible samplesheet from FASTQ inputs, launches the
nf-core/sarek pipeline via Nextflow, and locates the output VCF(s).

Supports:
  - Paired-end FASTQ → GVCF/VCF via GATK HaplotypeCaller
  - Single-sample and multi-sample modes
  - Dry-run mode (samplesheet generation + command preview only)
  - Custom Nextflow profiles (local, docker, singularity, conda, slurm)

Usage:
    from clawbio.common.sarek import SarekWrapper
    wrapper = SarekWrapper(config)
    vcf_path = wrapper.run(fastq_r1="sample_R1.fastq.gz", fastq_r2="sample_R2.fastq.gz")
"""

from __future__ import annotations

import csv
import json
import logging
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------

@dataclass
class SarekConfig:
    """Runtime configuration for the sarek wrapper."""

    # Nextflow / sarek settings
    nextflow_bin: str = "nextflow"
    sarek_version: str = "3.4.4"
    genome: str = "GATK.GRCh38"          # sarek genome alias
    tools: list[str] = field(default_factory=lambda: ["haplotypecaller"])
    profile: str = "docker"               # docker | singularity | conda | local
    work_dir: str = "/tmp/nf_work"
    output_dir: str = "sarek_output"

    # Sample metadata
    sample_id: str = "SAMPLE"
    sex: str = "XX"                       # XX or XY; used for chrX/Y calling
    status: int = 0                       # 0 = normal, 1 = tumour

    # Resource limits (passed as Nextflow params)
    max_cpus: int = 8
    max_memory: str = "32.GB"
    max_time: str = "24.h"

    # Bridge behaviour
    dry_run: bool = False
    skip_bqsr: bool = False               # Skip base quality score recalibration
    joint_germline: bool = False          # Multi-sample joint genotyping

    def to_nextflow_params(self) -> list[str]:
        """Convert to Nextflow CLI --param value pairs."""
        params = [
            "--genome", self.genome,
            "--tools", ",".join(self.tools),
            "--outdir", self.output_dir,
            "--max_cpus", str(self.max_cpus),
            "--max_memory", self.max_memory,
            "--max_time", self.max_time,
        ]
        if self.skip_bqsr:
            params.append("--skip_tools")
            params.append("baserecalibrator")
        if self.joint_germline:
            params.append("--joint_germline")
        return params


@dataclass
class SarekSample:
    """One row in the nf-core/sarek samplesheet (CSV format)."""

    patient: str
    sex: str          # XX / XY / NA
    status: int       # 0 = normal, 1 = tumour
    sample: str
    lane: str
    fastq_1: str
    fastq_2: str      # empty string for single-end


# ---------------------------------------------------------------------------
# Samplesheet builder
# ---------------------------------------------------------------------------

def build_samplesheet(
    fastq_r1: str | Path,
    fastq_r2: str | Path | None,
    output_path: str | Path,
    sample_id: str = "SAMPLE",
    sex: str = "XX",
    status: int = 0,
    lane: str = "L001",
) -> Path:
    """Generate a nf-core/sarek CSV samplesheet.

    Args:
        fastq_r1: Path to forward reads FASTQ (required).
        fastq_r2: Path to reverse reads FASTQ (None for single-end).
        output_path: Where to write the samplesheet CSV.
        sample_id: Sample identifier used throughout the pipeline.
        sex: Biological sex. "XX" or "XY" (affects variant calling on sex chromosomes).
        status: 0 = normal germline, 1 = tumour.
        lane: Sequencing lane identifier.

    Returns:
        Path to the written samplesheet.
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    sample = SarekSample(
        patient=sample_id,
        sex=sex,
        status=status,
        sample=sample_id,
        lane=lane,
        fastq_1=str(Path(fastq_r1).resolve()),
        fastq_2=str(Path(fastq_r2).resolve()) if fastq_r2 else "",
    )

    with open(output_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            fieldnames=["patient", "sex", "status", "sample", "lane", "fastq_1", "fastq_2"],
        )
        writer.writeheader()
        writer.writerow(asdict(sample))

    log.info("Samplesheet written to %s", output_path)
    return output_path


# ---------------------------------------------------------------------------
# Sarek wrapper
# ---------------------------------------------------------------------------

class SarekWrapper:
    """Orchestrate nf-core/sarek for WGS germline variant calling.

    Example:
        cfg = SarekConfig(profile="docker", output_dir="/results/sarek")
        wrapper = SarekWrapper(cfg)
        vcf = wrapper.run(
            fastq_r1="/data/sample_R1.fastq.gz",
            fastq_r2="/data/sample_R2.fastq.gz",
        )
        print(f"VCF ready at: {vcf}")
    """

    SAREK_NF_PIPELINE = "nf-core/sarek"

    def __init__(self, config: SarekConfig | None = None) -> None:
        self.config = config or SarekConfig()
        self._output_dir = Path(self.config.output_dir)

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        fastq_r1: str | Path,
        fastq_r2: str | Path | None = None,
        samplesheet: str | Path | None = None,
    ) -> Path:
        """Run the full sarek pipeline and return the path to the output VCF.

        Args:
            fastq_r1: Forward reads FASTQ.gz (ignored if samplesheet is provided).
            fastq_r2: Reverse reads FASTQ.gz (optional, ignored if samplesheet provided).
            samplesheet: Pre-built sarek samplesheet CSV (overrides fastq_r1/r2).

        Returns:
            Path to the genotyped VCF (or GVCF if joint_germline mode).

        Raises:
            RuntimeError: If Nextflow is not found or the pipeline fails.
            FileNotFoundError: If the expected VCF is absent after the run.
        """
        self._check_nextflow()

        # Build samplesheet if not provided
        if samplesheet is None:
            ss_path = self._output_dir / "samplesheet.csv"
            samplesheet = build_samplesheet(
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                output_path=ss_path,
                sample_id=self.config.sample_id,
                sex=self.config.sex,
                status=self.config.status,
            )
        else:
            samplesheet = Path(samplesheet)

        cmd = self._build_command(samplesheet)
        self._log_command(cmd)

        if self.config.dry_run:
            log.info("[DRY RUN] Skipping Nextflow execution.")
            return self._mock_vcf_path()

        self._execute(cmd)
        return self._locate_vcf()

    def check_environment(self) -> dict[str, str | bool]:
        """Check whether required tools are available.

        Returns:
            Dict with keys: nextflow_found, nextflow_version, docker_found,
            singularity_found, java_found.
        """
        status: dict[str, str | bool] = {}
        for tool in ["nextflow", "docker", "singularity", "java"]:
            path = shutil.which(tool)
            status[f"{tool}_found"] = path is not None
            if path and tool == "nextflow":
                try:
                    result = subprocess.run(
                        ["nextflow", "-version"],
                        capture_output=True, text=True, timeout=10
                    )
                    status["nextflow_version"] = result.stdout.strip().split("\n")[0]
                except Exception:
                    status["nextflow_version"] = "unknown"
        return status

    # ------------------------------------------------------------------
    # Internal helpers
    # ------------------------------------------------------------------

    def _check_nextflow(self) -> None:
        if not shutil.which(self.config.nextflow_bin):
            raise RuntimeError(
                f"Nextflow not found at '{self.config.nextflow_bin}'. "
                "Install with: curl -s https://get.nextflow.io | bash\n"
                "Or run with dry_run=True to generate the samplesheet without executing."
            )

    def _build_command(self, samplesheet: Path) -> list[str]:
        cmd = [
            self.config.nextflow_bin,
            "run",
            f"{self.SAREK_NF_PIPELINE}",
            "-r", self.config.sarek_version,
            "-profile", self.config.profile,
            "-work-dir", self.config.work_dir,
            "--input", str(samplesheet),
        ]
        cmd += self.config.to_nextflow_params()
        return cmd

    def _log_command(self, cmd: list[str]) -> None:
        pretty = " \\\n    ".join(cmd)
        log.info("Nextflow command:\n    %s", pretty)

    def _execute(self, cmd: list[str]) -> None:
        log.info("Launching nf-core/sarek (this may take several hours for WGS)...")
        start = time.time()
        try:
            result = subprocess.run(
                cmd,
                check=True,
                text=True,
                stdout=subprocess.PIPE,
                stderr=subprocess.STDOUT,
            )
            elapsed = time.time() - start
            log.info("sarek completed in %.0f s", elapsed)
            # Stream last 20 lines of output to log
            lines = (result.stdout or "").splitlines()
            for line in lines[-20:]:
                log.debug("sarek | %s", line)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"nf-core/sarek failed (exit code {exc.returncode}). "
                f"Check {self._output_dir / 'pipeline_info'} for logs."
            ) from exc

    def _locate_vcf(self) -> Path:
        """Find the canonical output VCF from the sarek output directory.

        sarek writes variants to:
          <outdir>/variant_calling/haplotypecaller/<sample>/*.vcf.gz  (single)
          <outdir>/variant_calling/haplotypecaller/<sample>/*.g.vcf.gz (gvcf)
          <outdir>/variant_calling/deepvariant/<sample>/*.vcf.gz
        """
        search_dirs = [
            self._output_dir / "variant_calling" / "haplotypecaller",
            self._output_dir / "variant_calling" / "deepvariant",
            self._output_dir / "variant_calling" / "freebayes",
        ]

        candidates: list[Path] = []
        for d in search_dirs:
            if d.exists():
                candidates += list(d.rglob("*.vcf.gz"))

        # Prefer non-gVCF files
        non_gvcf = [p for p in candidates if ".g.vcf" not in str(p)]
        if non_gvcf:
            vcf = sorted(non_gvcf)[-1]
            log.info("Found output VCF: %s", vcf)
            return vcf

        if candidates:
            vcf = sorted(candidates)[-1]
            log.info("Found output gVCF: %s", vcf)
            return vcf

        raise FileNotFoundError(
            f"No VCF found in {self._output_dir / 'variant_calling'}. "
            "Check that the pipeline completed successfully."
        )

    def _mock_vcf_path(self) -> Path:
        """Return a placeholder VCF path for dry-run mode."""
        return self._output_dir / "variant_calling" / "haplotypecaller" / \
               self.config.sample_id / f"{self.config.sample_id}.vcf.gz"

    # ------------------------------------------------------------------
    # Utility: generate a reproducibility bundle
    # ------------------------------------------------------------------

    def write_run_manifest(self, output_path: str | Path) -> Path:
        """Write a JSON manifest recording pipeline parameters for reproducibility."""
        manifest = {
            "pipeline": self.SAREK_NF_PIPELINE,
            "version": self.config.sarek_version,
            "genome": self.config.genome,
            "tools": self.config.tools,
            "profile": self.config.profile,
            "skip_bqsr": self.config.skip_bqsr,
            "joint_germline": self.config.joint_germline,
            "max_cpus": self.config.max_cpus,
            "max_memory": self.config.max_memory,
        }
        output_path = Path(output_path)
        output_path.write_text(json.dumps(manifest, indent=2))
        log.info("Run manifest written to %s", output_path)
        return output_path


# ---------------------------------------------------------------------------
# CLI (standalone use)
# ---------------------------------------------------------------------------

def _cli() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="Run nf-core/sarek on paired FASTQ files",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--fastq-r1", required=True, help="Forward reads FASTQ.gz")
    parser.add_argument("--fastq-r2", help="Reverse reads FASTQ.gz")
    parser.add_argument("--sample-id", default="SAMPLE", help="Sample identifier")
    parser.add_argument("--sex", default="XX", choices=["XX", "XY"], help="Biological sex")
    parser.add_argument("--genome", default="GATK.GRCh38", help="Reference genome alias")
    parser.add_argument("--profile", default="docker", help="Nextflow execution profile")
    parser.add_argument("--output-dir", default="sarek_output", help="Output directory")
    parser.add_argument("--work-dir", default="/tmp/nf_work", help="Nextflow work directory")
    parser.add_argument("--dry-run", action="store_true", help="Print command only, don't run")
    parser.add_argument("--skip-bqsr", action="store_true", help="Skip base quality recalibration")
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    cfg = SarekConfig(
        sample_id=args.sample_id,
        sex=args.sex,
        genome=args.genome,
        profile=args.profile,
        output_dir=args.output_dir,
        work_dir=args.work_dir,
        dry_run=args.dry_run,
        skip_bqsr=args.skip_bqsr,
    )
    wrapper = SarekWrapper(cfg)
    vcf = wrapper.run(fastq_r1=args.fastq_r1, fastq_r2=args.fastq_r2)
    print(f"\nOutput VCF: {vcf}")


if __name__ == "__main__":
    _cli()
