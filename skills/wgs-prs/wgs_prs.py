#!/usr/bin/env python3
"""
wgs_prs.py — WGS-PRS Skill: End-to-end WGS to Polygenic Risk Scores
ClawBio WGS-PRS Skill v0.1.0
Author: David de Lorenzo
License: MIT

Wires together:
  1. nf-core/sarek  (FASTQ → aligned BAM → genotyped VCF)
  2. vcf_qc         (normalisation → hard filtering → QC metrics)
  3. ClawBio gwas-prs (VCF → polygenic risk scores)
  4. Aggregated report (Markdown + JSON)

The canonical entry-point is the `WgsToPrsBridge` class or the CLI.

Usage:
    # Python API
    from skills.wgs_prs.wgs_prs import WgsToPrsBridge, BridgeConfig
    bridge = WgsToPrsBridge(BridgeConfig(output_dir="results"))
    report = bridge.run(
        fastq_r1="sample_R1.fastq.gz",
        fastq_r2="sample_R2.fastq.gz",
    )

    # Start from an existing VCF (skip sarek)
    report = bridge.run(input_vcf="sample.vcf.gz")

    # CLI
    python skills/wgs-prs/wgs_prs.py --fastq-r1 sample_R1.fastq.gz \\
        --fastq-r2 sample_R2.fastq.gz --output-dir results/
"""

from __future__ import annotations

import json
import logging
import shutil
import subprocess
import sys
import time
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Imports — clawbio.common shared library
# ---------------------------------------------------------------------------

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.sarek import SarekConfig, SarekWrapper
from clawbio.common.vcf_qc import QcConfig, VcfQC, QcResult

log = logging.getLogger(__name__)

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Polygenic risk scores reflect "
    "statistical associations from population studies and do not determine "
    "individual outcomes. Consult a healthcare professional before making "
    "any medical decisions based on genetic information."
)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class BridgeConfig:
    """Master configuration for the WGS-PRS bridge."""

    output_dir: str = "wgs_prs_output"
    sample_id: str = "SAMPLE"
    sex: str = "XX"                        # XX / XY

    # Stage 1 — sarek
    sarek: SarekConfig = field(default_factory=SarekConfig)

    # Stage 2 — VCF QC
    qc: QcConfig = field(default_factory=QcConfig)

    # Stage 3 — PRS
    clawbio_root: str = ""                 # path to ClawBio-0.5.0; auto-detected if empty
    prs_traits: list[str] = field(default_factory=list)   # empty = all curated scores
    pgs_ids: list[str] = field(default_factory=list)      # explicit PGS IDs to compute

    # Bridge behaviour
    skip_sarek: bool = False               # set True when providing a VCF directly
    dry_run: bool = False                  # skip all subprocess calls
    fail_fast: bool = True                 # abort on QC failure (False = warn and continue)
    log_level: str = "INFO"


# ---------------------------------------------------------------------------
# Stage result containers
# ---------------------------------------------------------------------------

@dataclass
class StageResult:
    name: str
    status: str = "pending"   # pending | running | success | skipped | failed
    output: Optional[Path] = None
    error: str = ""
    duration_s: float = 0.0


@dataclass
class BridgeReport:
    sample_id: str
    timestamp: str
    output_dir: Path
    stages: list[StageResult] = field(default_factory=list)
    qc_metrics: Optional[dict] = None
    prs_summary: Optional[dict] = None
    overall_status: str = "pending"
    report_md: Optional[Path] = None
    report_json: Optional[Path] = None

    def succeeded(self) -> bool:
        return self.overall_status == "success"


# ---------------------------------------------------------------------------
# Main orchestrator
# ---------------------------------------------------------------------------

class WgsToPrsBridge:
    """End-to-end WGS → PRS pipeline orchestrator.

    Stages:
      1. Variant calling    nf-core/sarek  (FASTQ → VCF)
      2. VCF QC             bcftools norm + filter + stats
      3. PRS scoring        ClawBio gwas-prs skill
      4. Report generation  Markdown + JSON
    """

    def __init__(self, config: BridgeConfig | None = None) -> None:
        self.config = config or BridgeConfig()
        self._output_dir = Path(self.config.output_dir)
        self._output_dir.mkdir(parents=True, exist_ok=True)
        self._setup_logging()
        self._clawbio_root = self._find_clawbio_root()

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(
        self,
        fastq_r1: str | Path | None = None,
        fastq_r2: str | Path | None = None,
        input_vcf: str | Path | None = None,
        samplesheet: str | Path | None = None,
    ) -> BridgeReport:
        """Execute the full bridge pipeline.

        Call with either FASTQ inputs or a pre-existing VCF:
          bridge.run(fastq_r1=..., fastq_r2=...)    # full pipeline
          bridge.run(input_vcf=...)                  # skip sarek

        Args:
            fastq_r1: Forward reads FASTQ.gz (required unless input_vcf is given).
            fastq_r2: Reverse reads FASTQ.gz (optional, paired-end recommended).
            input_vcf: Pre-existing VCF — skips Stage 1.
            samplesheet: Pre-built sarek samplesheet CSV (overrides fastq_r1/r2).

        Returns:
            BridgeReport with status, metrics, and output file paths.
        """
        report = BridgeReport(
            sample_id=self.config.sample_id,
            timestamp=datetime.now(timezone.utc).isoformat(),
            output_dir=self._output_dir,
        )

        log.info("=" * 60)
        log.info("ClawBio WGS-PRS Bridge  |  sample: %s", self.config.sample_id)
        log.info("Output dir: %s", self._output_dir)
        log.info("=" * 60)

        # ---- Stage 1: Variant calling --------------------------------
        vcf_path: Optional[Path] = None
        if input_vcf:
            log.info("Stage 1: SKIPPED (using provided VCF: %s)", input_vcf)
            vcf_path = Path(input_vcf)
            report.stages.append(StageResult(name="sarek", status="skipped", output=vcf_path))
        else:
            stage1 = self._run_stage1(fastq_r1, fastq_r2, samplesheet)
            report.stages.append(stage1)
            if stage1.status == "failed":
                report.overall_status = "failed"
                self._write_report(report)
                return report
            vcf_path = stage1.output

        # ---- Stage 2: VCF QC -----------------------------------------
        stage2, qc_result = self._run_stage2(vcf_path)
        report.stages.append(stage2)
        report.qc_metrics = self._load_json(qc_result.metrics_json) if qc_result else None

        if stage2.status == "failed" and self.config.fail_fast:
            log.error("QC failed — aborting. Set fail_fast=False to continue despite QC failure.")
            report.overall_status = "failed"
            self._write_report(report)
            return report

        canonical_vcf = qc_result.canonical_vcf if qc_result else vcf_path

        # ---- Stage 3: PRS scoring ------------------------------------
        stage3 = self._run_stage3(canonical_vcf)
        report.stages.append(stage3)
        if stage3.output and stage3.output.exists():
            report.prs_summary = self._load_json(
                stage3.output / "result.json"
            )

        # ---- Stage 4: Report generation ------------------------------
        stage4 = self._run_stage4(report)
        report.stages.append(stage4)
        report.report_md = stage4.output
        report.report_json = self._output_dir / "bridge_report.json"

        # Final status
        failed = [s for s in report.stages if s.status == "failed"]
        report.overall_status = "failed" if failed else "success"

        self._write_report(report)

        log.info("=" * 60)
        log.info("Bridge complete — status: %s", report.overall_status.upper())
        if report.report_md:
            log.info("Report: %s", report.report_md)
        log.info("=" * 60)

        return report

    # ------------------------------------------------------------------
    # Stage 1 — nf-core/sarek
    # ------------------------------------------------------------------

    def _run_stage1(
        self,
        fastq_r1,
        fastq_r2,
        samplesheet,
    ) -> StageResult:
        stage = StageResult(name="sarek")
        stage.status = "running"
        log.info("\n--- Stage 1: Variant calling (nf-core/sarek) ---")

        # Forward bridge config into sarek config
        cfg = self.config.sarek
        cfg.sample_id = self.config.sample_id
        cfg.sex = self.config.sex
        cfg.output_dir = str(self._output_dir / "sarek_output")
        cfg.dry_run = self.config.dry_run

        t0 = time.time()
        try:
            wrapper = SarekWrapper(cfg)
            vcf = wrapper.run(
                fastq_r1=fastq_r1,
                fastq_r2=fastq_r2,
                samplesheet=samplesheet,
            )
            wrapper.write_run_manifest(self._output_dir / "sarek_manifest.json")
            stage.status = "success"
            stage.output = vcf
        except Exception as exc:
            stage.status = "failed"
            stage.error = str(exc)
            log.error("Stage 1 FAILED: %s", exc)

        stage.duration_s = time.time() - t0
        return stage

    # ------------------------------------------------------------------
    # Stage 2 — VCF QC
    # ------------------------------------------------------------------

    def _run_stage2(self, vcf_path: Path | None) -> tuple[StageResult, Optional[QcResult]]:
        stage = StageResult(name="vcf_qc")
        if vcf_path is None:
            stage.status = "failed"
            stage.error = "No VCF available for QC (Stage 1 may have failed)"
            return stage, None

        stage.status = "running"
        log.info("\n--- Stage 2: VCF QC ---")

        qc_dir = self._output_dir / "vcf_qc"
        t0 = time.time()
        try:
            qc = VcfQC(self.config.qc)
            result = qc.run(input_vcf=vcf_path, output_dir=qc_dir)
            stage.output = result.canonical_vcf
            stage.duration_s = time.time() - t0

            if result.passes_qc:
                stage.status = "success"
                log.info("VCF QC PASSED")
            else:
                stage.status = "failed"
                stage.error = "; ".join(result.fail_reasons)
                log.warning("VCF QC FAILED: %s", stage.error)

            return stage, result
        except Exception as exc:
            stage.status = "failed"
            stage.error = str(exc)
            stage.duration_s = time.time() - t0
            log.error("Stage 2 FAILED: %s", exc)
            return stage, None

    # ------------------------------------------------------------------
    # Stage 3 — PRS scoring via ClawBio gwas-prs
    # ------------------------------------------------------------------

    def _run_stage3(self, vcf_path: Path | None) -> StageResult:
        stage = StageResult(name="prs_scoring")
        if vcf_path is None:
            stage.status = "failed"
            stage.error = "No canonical VCF available for PRS scoring"
            return stage

        stage.status = "running"
        log.info("\n--- Stage 3: PRS scoring ---")

        prs_output = self._output_dir / "prs_output"
        t0 = time.time()

        if self.config.dry_run:
            log.info("[DRY RUN] Skipping PRS scoring.")
            stage.status = "skipped"
            stage.output = prs_output
            return stage

        gwas_prs_script = self._find_gwas_prs_script()
        if gwas_prs_script is None:
            stage.status = "failed"
            stage.error = (
                "gwas_prs.py not found. Set clawbio_root in BridgeConfig "
                "or ensure ClawBio-0.5.0 is adjacent to this script."
            )
            log.error(stage.error)
            return stage

        cmd = [
            sys.executable, str(gwas_prs_script),
            "--input", str(vcf_path),
            "--output", str(prs_output),
        ]
        if self.config.pgs_ids:
            cmd += ["--pgs-id", self.config.pgs_ids[0]]
        elif self.config.prs_traits:
            cmd += ["--trait", self.config.prs_traits[0]]

        log.info("Running gwas-prs: %s", " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, text=True)
            stage.status = "success"
            stage.output = prs_output
            log.info("PRS scoring complete → %s", prs_output)
        except subprocess.CalledProcessError as exc:
            stage.status = "failed"
            stage.error = f"gwas_prs.py exited with code {exc.returncode}"
            log.error("Stage 3 FAILED: %s", stage.error)

        stage.duration_s = time.time() - t0
        return stage

    # ------------------------------------------------------------------
    # Stage 4 — Aggregated report
    # ------------------------------------------------------------------

    def _run_stage4(self, report: BridgeReport) -> StageResult:
        stage = StageResult(name="report")
        stage.status = "running"
        log.info("\n--- Stage 4: Report generation ---")

        t0 = time.time()
        try:
            md_path = self._write_markdown_report(report)
            stage.output = md_path
            stage.status = "success"
        except Exception as exc:
            stage.status = "failed"
            stage.error = str(exc)
            log.error("Stage 4 FAILED: %s", exc)

        stage.duration_s = time.time() - t0
        return stage

    def _write_markdown_report(self, report: BridgeReport) -> Path:
        """Compose the final Markdown report from stage outputs."""
        lines = [
            f"# ClawBio WGS-PRS Bridge Report",
            f"",
            f"**Sample:** {report.sample_id}  ",
            f"**Generated:** {report.timestamp}  ",
            f"**Output directory:** `{report.output_dir}`",
            f"",
            f"---",
            f"",
            f"## Pipeline Stages",
            f"",
            f"| Stage | Status | Duration |",
            f"|-------|--------|----------|",
        ]
        for s in report.stages:
            icon = {"success": "✅", "failed": "❌", "skipped": "⏭️", "running": "⏳", "pending": "⏸️"}.get(s.status, "❓")
            dur = f"{s.duration_s:.1f}s" if s.duration_s else "—"
            lines.append(f"| {s.name} | {icon} {s.status} | {dur} |")

        # QC metrics section
        lines += ["", "---", "", "## VCF QC Metrics", ""]
        if report.qc_metrics:
            m = report.qc_metrics.get("metrics", {})
            lines += [
                f"**QC Status:** {report.qc_metrics.get('qc_status', 'N/A')}",
                f"",
                f"| Metric | Value |",
                f"|--------|-------|",
                f"| Total variants | {m.get('total_variants', 'N/A'):,} |" if isinstance(m.get('total_variants'), int) else f"| Total variants | N/A |",
                f"| SNPs | {m.get('snp_count', 'N/A'):,} |" if isinstance(m.get('snp_count'), int) else "| SNPs | N/A |",
                f"| Indels | {m.get('indel_count', 'N/A'):,} |" if isinstance(m.get('indel_count'), int) else "| Indels | N/A |",
                f"| Ti/Tv ratio | {m.get('titv_ratio', 'N/A')} |",
                f"| Het/Hom ratio | {m.get('het_hom_ratio', 'N/A')} |",
                f"| Filtered variants | {m.get('filtered_out', 'N/A')} |",
            ]
            if report.qc_metrics.get("fail_reasons"):
                lines += ["", "**Fail reasons:**"]
                for r in report.qc_metrics["fail_reasons"]:
                    lines.append(f"- {r}")
            if report.qc_metrics.get("warnings"):
                lines += ["", "**Warnings:**"]
                for w in report.qc_metrics["warnings"]:
                    lines.append(f"- ⚠️ {w}")
        else:
            lines.append("_QC metrics not available._")

        # PRS section
        lines += ["", "---", "", "## Polygenic Risk Scores", ""]
        prs_output = self._output_dir / "prs_output"
        prs_report = prs_output / "report.md"
        if prs_report.exists():
            lines.append(f"PRS report available at: `{prs_report}`")
            lines.append("")
            # Embed a short excerpt
            prs_lines = prs_report.read_text().splitlines()
            lines += prs_lines[:40]
            if len(prs_lines) > 40:
                lines.append(f"\n_... ({len(prs_lines) - 40} more lines — see full report)_")
        else:
            lines.append("_PRS scores not available (scoring stage did not complete)._")

        # Errors section
        failed_stages = [s for s in report.stages if s.status == "failed"]
        if failed_stages:
            lines += ["", "---", "", "## Errors", ""]
            for s in failed_stages:
                lines.append(f"**{s.name}:** {s.error}")

        lines += [
            "",
            "---",
            "",
            f"> {DISCLAIMER}",
        ]

        md_path = self._output_dir / "bridge_report.md"
        md_path.write_text("\n".join(lines))
        log.info("Markdown report written to %s", md_path)
        return md_path

    # ------------------------------------------------------------------
    # Internal utilities
    # ------------------------------------------------------------------

    def _write_report(self, report: BridgeReport) -> None:
        """Write the machine-readable JSON report."""
        data = {
            "sample_id": report.sample_id,
            "timestamp": report.timestamp,
            "output_dir": str(report.output_dir),
            "overall_status": report.overall_status,
            "stages": [
                {
                    "name": s.name,
                    "status": s.status,
                    "output": str(s.output) if s.output else None,
                    "error": s.error,
                    "duration_s": round(s.duration_s, 2),
                }
                for s in report.stages
            ],
            "qc_metrics": report.qc_metrics,
        }
        json_path = self._output_dir / "bridge_report.json"
        json_path.write_text(json.dumps(data, indent=2))
        report.report_json = json_path

    def _find_clawbio_root(self) -> Optional[Path]:
        """Return the ClawBio project root (this skill lives inside it)."""
        if self.config.clawbio_root:
            p = Path(self.config.clawbio_root)
            return p if p.exists() else None
        if (_PROJECT_ROOT / "clawbio.py").exists():
            return _PROJECT_ROOT
        return None

    def _find_gwas_prs_script(self) -> Optional[Path]:
        if self._clawbio_root:
            p = self._clawbio_root / "skills" / "gwas-prs" / "gwas_prs.py"
            if p.exists():
                return p
        return None

    def _setup_logging(self) -> None:
        level = getattr(logging, self.config.log_level.upper(), logging.INFO)
        logging.basicConfig(
            level=level,
            format="%(asctime)s  %(levelname)-8s  %(message)s",
            datefmt="%H:%M:%S",
        )

    @staticmethod
    def _load_json(path: Optional[Path]) -> Optional[dict]:
        if path and path.exists():
            try:
                return json.loads(path.read_text())
            except Exception:
                return None
        return None


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------

def _cli() -> None:
    import argparse

    parser = argparse.ArgumentParser(
        description="ClawBio WGS-PRS Bridge — FASTQ/VCF → Polygenic Risk Scores",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
        epilog=(
            "Examples:\n"
            "  # Full pipeline from FASTQ:\n"
            "  %(prog)s --fastq-r1 sample_R1.fq.gz --fastq-r2 sample_R2.fq.gz\n\n"
            "  # Start from a VCF (skip sarek):\n"
            "  %(prog)s --input-vcf sample.vcf.gz\n\n"
            "  # Dry run (no subprocess calls):\n"
            "  %(prog)s --fastq-r1 sample_R1.fq.gz --dry-run\n"
        ),
    )

    # Input options
    inp = parser.add_argument_group("Input (choose one)")
    inp.add_argument("--fastq-r1", help="Forward reads FASTQ.gz")
    inp.add_argument("--fastq-r2", help="Reverse reads FASTQ.gz (omit for single-end)")
    inp.add_argument("--input-vcf", help="Pre-existing VCF (skips sarek)")
    inp.add_argument("--samplesheet", help="Pre-built nf-core/sarek samplesheet CSV")

    # Sample metadata
    meta = parser.add_argument_group("Sample metadata")
    meta.add_argument("--sample-id", default="SAMPLE")
    meta.add_argument("--sex", default="XX", choices=["XX", "XY"])

    # sarek options
    sk = parser.add_argument_group("nf-core/sarek options")
    sk.add_argument("--genome", default="GATK.GRCh38")
    sk.add_argument("--profile", default="docker")
    sk.add_argument("--sarek-version", default="3.4.4")
    sk.add_argument("--skip-bqsr", action="store_true")

    # QC options
    qco = parser.add_argument_group("VCF QC options")
    qco.add_argument("--reference-fasta", default="", help="Reference FASTA for bcftools norm")
    qco.add_argument("--min-qual", type=float, default=30.0)
    qco.add_argument("--min-depth", type=int, default=10)
    qco.add_argument("--no-fail-fast", action="store_true",
                     help="Continue even if QC fails")

    # PRS options
    prs = parser.add_argument_group("PRS options")
    prs.add_argument("--trait", help="Trait to score (e.g. 'type 2 diabetes')")
    prs.add_argument("--pgs-id", help="Specific PGS Catalog ID (e.g. PGS000013)")
    prs.add_argument("--clawbio-root", default="", help="Path to ClawBio-0.5.0 directory")

    # Output / behaviour
    out = parser.add_argument_group("Output")
    out.add_argument("--output-dir", default="wgs_prs_output")
    out.add_argument("--dry-run", action="store_true")
    out.add_argument("--log-level", default="INFO", choices=["DEBUG", "INFO", "WARNING", "ERROR"])

    args = parser.parse_args()

    if not args.fastq_r1 and not args.input_vcf:
        parser.error("Provide either --fastq-r1 (with optional --fastq-r2) or --input-vcf")

    cfg = BridgeConfig(
        output_dir=args.output_dir,
        sample_id=args.sample_id,
        sex=args.sex,
        clawbio_root=args.clawbio_root,
        prs_traits=[args.trait] if args.trait else [],
        pgs_ids=[args.pgs_id] if args.pgs_id else [],
        dry_run=args.dry_run,
        fail_fast=not args.no_fail_fast,
        log_level=args.log_level,
        sarek=SarekConfig(
            genome=args.genome,
            profile=args.profile,
            sarek_version=args.sarek_version,
            skip_bqsr=args.skip_bqsr,
        ),
        qc=QcConfig(
            reference_fasta=args.reference_fasta,
            min_qual=args.min_qual,
            min_depth=args.min_depth,
        ),
    )

    bridge = WgsToPrsBridge(cfg)
    report = bridge.run(
        fastq_r1=args.fastq_r1,
        fastq_r2=args.fastq_r2,
        input_vcf=args.input_vcf,
        samplesheet=args.samplesheet,
    )

    sys.exit(0 if report.succeeded() else 1)


if __name__ == "__main__":
    _cli()
