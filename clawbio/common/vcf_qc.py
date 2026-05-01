#!/usr/bin/env python3
"""
vcf_qc.py — VCF normalisation, hard filtering, and QC metrics for the ClawBio bridge
ClawBio WGS-PRS Bridge v0.1.0
Author: David de Lorenzo
License: MIT

Produces:
  - A normalised, filtered "canonical VCF" ready for PRS scoring
  - A QC metrics JSON with Ti/Tv ratio, het/hom ratio, variant counts, and
    per-score coverage estimates

Requires:
  - bcftools >= 1.17 (for normalisation and stats)

Falls back gracefully if bcftools is unavailable (Python-only mode with
reduced metrics — suitable for quick testing on pre-normalised VCFs).

Usage:
    from clawbio.common.vcf_qc import VcfQC, QcConfig
    qc = VcfQC(QcConfig(reference_fasta="/ref/GRCh38.fa"))
    result = qc.run(input_vcf="sample.vcf.gz", output_dir="qc_output")
    if result.passes_qc:
        print(f"Canonical VCF: {result.canonical_vcf}")
    else:
        print(f"QC FAILED: {result.fail_reasons}")
"""

from __future__ import annotations

import gzip
import json
import logging
import re
import shutil
import subprocess
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Optional

log = logging.getLogger(__name__)


# ---------------------------------------------------------------------------
# Configuration
# ---------------------------------------------------------------------------

@dataclass
class QcConfig:
    """Hard-filter thresholds and QC pass/fail criteria."""

    # bcftools normalisation
    reference_fasta: str = ""                  # required for left-align; empty = skip
    split_multiallelics: bool = True           # decompose multi-allelic sites

    # Hard filters (applied via bcftools filter)
    min_qual: float = 30.0                     # QUAL score threshold
    min_depth: int = 10                        # FORMAT/DP threshold
    max_missing_rate: float = 0.10             # fraction allowed to be ./. calls

    # QC pass/fail thresholds
    min_titv_ratio: float = 1.8                # WGS coding: ~3.0; whole genome: ~2.0
    max_titv_ratio: float = 2.5
    min_het_hom_ratio: float = 1.0             # typical WGS germline range: 1.2–2.0
    max_het_hom_ratio: float = 3.0
    min_snp_count: int = 100                   # sanity check — expect millions for WGS

    # PRS coverage threshold — warn if below this for any score
    min_prs_coverage_pct: float = 50.0


# ---------------------------------------------------------------------------
# Result dataclass
# ---------------------------------------------------------------------------

@dataclass
class QcResult:
    """Output of VcfQC.run()."""

    canonical_vcf: Optional[Path] = None       # filtered, normalised VCF
    metrics_json: Optional[Path] = None        # QC metrics file path
    passes_qc: bool = False
    fail_reasons: list[str] = field(default_factory=list)
    warnings: list[str] = field(default_factory=list)

    # Core metrics
    total_variants: int = 0
    snp_count: int = 0
    indel_count: int = 0
    titv_ratio: Optional[float] = None
    het_count: int = 0
    hom_alt_count: int = 0
    het_hom_ratio: Optional[float] = None
    filtered_out: int = 0

    def summary(self) -> str:
        lines = [
            "=== VCF QC Summary ===",
            f"  Total variants : {self.total_variants:,}",
            f"  SNPs           : {self.snp_count:,}",
            f"  Indels         : {self.indel_count:,}",
            f"  Ti/Tv ratio    : {self.titv_ratio:.3f}" if self.titv_ratio else "  Ti/Tv ratio    : N/A",
            f"  Het/Hom ratio  : {self.het_hom_ratio:.3f}" if self.het_hom_ratio else "  Het/Hom ratio  : N/A",
            f"  Filtered out   : {self.filtered_out:,}",
            f"  QC status      : {'PASS' if self.passes_qc else 'FAIL'}",
        ]
        if self.fail_reasons:
            lines.append(f"  Fail reasons   : {'; '.join(self.fail_reasons)}")
        if self.warnings:
            lines.append(f"  Warnings       : {'; '.join(self.warnings)}")
        return "\n".join(lines)


# ---------------------------------------------------------------------------
# Main class
# ---------------------------------------------------------------------------

class VcfQC:
    """Normalise and QC-check a VCF for downstream PRS scoring.

    Pipeline steps:
      1. bcftools norm  — left-align indels, split multiallelic sites
      2. bcftools filter — apply hard QUAL/DP thresholds
      3. bcftools stats  — compute Ti/Tv, het/hom, variant counts
      4. Pass/fail evaluation against QcConfig thresholds
      5. Write metrics JSON
    """

    def __init__(self, config: QcConfig | None = None) -> None:
        self.config = config or QcConfig()
        self._bcftools = shutil.which("bcftools")
        if not self._bcftools:
            log.warning(
                "bcftools not found — normalisation and hard filtering will be skipped. "
                "Install with: conda install -c bioconda bcftools"
            )

    # ------------------------------------------------------------------
    # Public API
    # ------------------------------------------------------------------

    def run(self, input_vcf: str | Path, output_dir: str | Path) -> QcResult:
        """Run the full QC pipeline on a single-sample VCF.

        Args:
            input_vcf: Path to the input VCF or VCF.gz.
            output_dir: Directory to write canonical VCF and metrics JSON.

        Returns:
            QcResult with pass/fail status and all metrics populated.
        """
        input_vcf = Path(input_vcf)
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)

        result = QcResult()

        if self._bcftools:
            normalised = self._normalise(input_vcf, output_dir)
            filtered, n_filtered = self._hard_filter(normalised, output_dir)
            result.filtered_out = n_filtered
            result.canonical_vcf = filtered
            stats = self._compute_stats(filtered)
        else:
            log.warning("Falling back to Python-only VCF parsing (reduced metrics).")
            result.canonical_vcf = input_vcf
            stats = self._python_stats(input_vcf)

        self._populate_result(result, stats)
        self._evaluate_pass_fail(result)

        metrics_path = output_dir / "qc_metrics.json"
        self._write_metrics(result, metrics_path)
        result.metrics_json = metrics_path

        log.info("\n%s", result.summary())
        return result

    # ------------------------------------------------------------------
    # Step 1: Normalisation
    # ------------------------------------------------------------------

    def _normalise(self, input_vcf: Path, output_dir: Path) -> Path:
        """Left-align indels and optionally split multiallelic sites."""
        out = output_dir / "normalised.vcf.gz"

        cmd = [self._bcftools, "norm"]
        if self.config.reference_fasta:
            cmd += ["-f", self.config.reference_fasta]
        if self.config.split_multiallelics:
            cmd += ["-m", "-any"]           # split multiallelic into biallelic
        cmd += ["-O", "z", "-o", str(out), str(input_vcf)]

        log.info("Normalising VCF: %s", " ".join(cmd))
        self._run_cmd(cmd, step="norm")

        # Index
        self._run_cmd([self._bcftools, "index", str(out)], step="index-norm")
        return out

    # ------------------------------------------------------------------
    # Step 2: Hard filtering
    # ------------------------------------------------------------------

    def _hard_filter(self, normalised: Path, output_dir: Path) -> tuple[Path, int]:
        """Apply QUAL and DP hard filters. Returns (filtered_vcf, n_removed)."""
        out = output_dir / "canonical.vcf.gz"
        cfg = self.config

        # Count before
        n_before = self._count_variants(normalised)

        filter_expr = (
            f"QUAL < {cfg.min_qual} || "
            f"FORMAT/DP < {cfg.min_depth}"
        )
        cmd = [
            self._bcftools, "filter",
            "-e", filter_expr,
            "-s", "CLAWBIO_HARD_FILTER",    # soft-flag rather than remove
            "--IndelGap", "5",
            "-O", "z", "-o", str(out),
            str(normalised),
        ]
        log.info("Applying hard filters: %s", filter_expr)
        self._run_cmd(cmd, step="filter")

        # Remove soft-flagged variants
        passed = output_dir / "canonical_pass.vcf.gz"
        self._run_cmd(
            [self._bcftools, "view", "-f", "PASS,.", "-O", "z", "-o", str(passed), str(out)],
            step="filter-pass",
        )
        self._run_cmd([self._bcftools, "index", str(passed)], step="index-pass")

        n_after = self._count_variants(passed)
        n_filtered = max(0, n_before - n_after)
        log.info("Hard filter removed %d / %d variants", n_filtered, n_before)
        return passed, n_filtered

    # ------------------------------------------------------------------
    # Step 3: Stats via bcftools stats
    # ------------------------------------------------------------------

    def _compute_stats(self, vcf: Path) -> dict:
        """Run bcftools stats and parse the summary numbers table."""
        cmd = [self._bcftools, "stats", str(vcf)]
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
        return self._parse_bcftools_stats(result.stdout)

    def _parse_bcftools_stats(self, stats_output: str) -> dict:
        """Extract key metrics from bcftools stats SN (summary numbers) section."""
        metrics: dict = {}
        for line in stats_output.splitlines():
            if line.startswith("SN"):
                # bcftools stats SN format: SN \t 0 \t key: \t value
                parts = line.split("\t")
                if len(parts) >= 4:
                    key = parts[2].strip()          # e.g. "number of SNPs:"
                    try:
                        val = float(parts[3].strip())
                    except ValueError:
                        continue
                    metrics[key] = val

            # Ti/Tv ratio is in the TSTV section
            if line.startswith("TSTV"):
                parts = line.split("\t")
                if len(parts) >= 5:
                    try:
                        metrics["titv_ratio"] = float(parts[4])
                    except (ValueError, IndexError):
                        pass

        return metrics

    # ------------------------------------------------------------------
    # Step 3 fallback: Python-only stats (no bcftools)
    # ------------------------------------------------------------------

    def _python_stats(self, vcf: Path) -> dict:
        """Parse a VCF with Python and return a best-effort stats dict."""
        metrics: dict = {
            "number of SNPs:": 0,
            "number of indels:": 0,
            "number of heterozygous SNPs:": 0,
            "number of homozygous SNPs:": 0,
            "titv_ratio": None,
        }

        opener = gzip.open if str(vcf).endswith(".gz") else open
        transitions = {"AG", "GA", "CT", "TC"}

        ti_count = 0
        tv_count = 0

        with opener(vcf, "rt") as fh:
            for line in fh:
                if line.startswith("#"):
                    continue
                parts = line.split("\t")
                if len(parts) < 10:
                    continue
                ref, alt = parts[3], parts[4]
                fmt_idx = parts[8].split(":").index("GT") if "GT" in parts[8] else -1
                gt_raw = parts[9].split(":")[fmt_idx] if fmt_idx >= 0 else "."

                # SNP vs indel
                is_snp = len(ref) == 1 and len(alt) == 1 and alt != "."
                if is_snp:
                    metrics["number of SNPs:"] += 1
                    pair = ref.upper() + alt.upper()
                    if pair in transitions:
                        ti_count += 1
                    else:
                        tv_count += 1
                elif alt != ".":
                    metrics["number of indels:"] += 1

                # Het / hom
                gt = re.split(r"[|/]", gt_raw.replace(".", "0"))
                if len(gt) == 2:
                    if gt[0] != gt[1]:
                        metrics["number of heterozygous SNPs:"] += 1
                    elif gt[0] != "0":
                        metrics["number of homozygous SNPs:"] += 1

        if tv_count > 0:
            metrics["titv_ratio"] = round(ti_count / tv_count, 4)

        return metrics

    # ------------------------------------------------------------------
    # Step 4: Populate result and evaluate pass/fail
    # ------------------------------------------------------------------

    def _populate_result(self, result: QcResult, stats: dict) -> None:
        result.snp_count = int(stats.get("number of SNPs:", 0))
        result.indel_count = int(stats.get("number of indels:", 0))
        result.total_variants = result.snp_count + result.indel_count
        result.het_count = int(stats.get("number of heterozygous SNPs:", 0))
        result.hom_alt_count = int(stats.get("number of homozygous SNPs:", 0))
        result.titv_ratio = stats.get("titv_ratio")

        if result.hom_alt_count > 0:
            result.het_hom_ratio = round(result.het_count / result.hom_alt_count, 4)

    def _evaluate_pass_fail(self, result: QcResult) -> None:
        cfg = self.config
        reasons = []
        warnings = []

        if result.snp_count < cfg.min_snp_count:
            reasons.append(
                f"Too few SNPs: {result.snp_count} < {cfg.min_snp_count}"
            )

        if result.titv_ratio is not None:
            if result.titv_ratio < cfg.min_titv_ratio:
                reasons.append(
                    f"Ti/Tv ratio too low: {result.titv_ratio:.3f} < {cfg.min_titv_ratio}"
                )
            elif result.titv_ratio > cfg.max_titv_ratio:
                warnings.append(
                    f"Ti/Tv ratio unusually high: {result.titv_ratio:.3f} > {cfg.max_titv_ratio} "
                    "(possible contamination or amplicon data)"
                )
        else:
            warnings.append("Could not compute Ti/Tv ratio")

        if result.het_hom_ratio is not None:
            if result.het_hom_ratio < cfg.min_het_hom_ratio:
                reasons.append(
                    f"Het/Hom ratio too low: {result.het_hom_ratio:.3f} < {cfg.min_het_hom_ratio} "
                    "(possible inbreeding or contamination)"
                )
            elif result.het_hom_ratio > cfg.max_het_hom_ratio:
                reasons.append(
                    f"Het/Hom ratio too high: {result.het_hom_ratio:.3f} > {cfg.max_het_hom_ratio} "
                    "(possible mixed sample)"
                )
        else:
            warnings.append("Could not compute Het/Hom ratio")

        result.fail_reasons = reasons
        result.warnings = warnings
        result.passes_qc = len(reasons) == 0

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    def _count_variants(self, vcf: Path) -> int:
        try:
            result = subprocess.run(
                [self._bcftools, "view", "-H", str(vcf)],
                capture_output=True, text=True, check=True
            )
            return sum(1 for _ in result.stdout.splitlines())
        except Exception:
            return 0

    def _run_cmd(self, cmd: list[str], step: str = "") -> None:
        log.debug("Running [%s]: %s", step, " ".join(cmd))
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
        except subprocess.CalledProcessError as exc:
            raise RuntimeError(
                f"bcftools step '{step}' failed:\n{exc.stderr}"
            ) from exc

    def _write_metrics(self, result: QcResult, path: Path) -> None:
        metrics = {
            "qc_status": "PASS" if result.passes_qc else "FAIL",
            "fail_reasons": result.fail_reasons,
            "warnings": result.warnings,
            "metrics": {
                "total_variants": result.total_variants,
                "snp_count": result.snp_count,
                "indel_count": result.indel_count,
                "titv_ratio": result.titv_ratio,
                "het_count": result.het_count,
                "hom_alt_count": result.hom_alt_count,
                "het_hom_ratio": result.het_hom_ratio,
                "filtered_out": result.filtered_out,
            },
        }
        path.write_text(json.dumps(metrics, indent=2))
        log.info("QC metrics written to %s", path)


# ---------------------------------------------------------------------------
# CLI (standalone use)
# ---------------------------------------------------------------------------

def _cli() -> None:
    import argparse
    parser = argparse.ArgumentParser(
        description="Normalise and QC-check a VCF for PRS scoring",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )
    parser.add_argument("--input", required=True, help="Input VCF or VCF.gz")
    parser.add_argument("--output-dir", default="vcf_qc_output", help="Output directory")
    parser.add_argument("--reference", default="", help="Reference FASTA for normalisation")
    parser.add_argument("--min-qual", type=float, default=30.0)
    parser.add_argument("--min-depth", type=int, default=10)
    args = parser.parse_args()

    logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")

    cfg = QcConfig(
        reference_fasta=args.reference,
        min_qual=args.min_qual,
        min_depth=args.min_depth,
    )
    qc = VcfQC(cfg)
    result = qc.run(input_vcf=args.input, output_dir=args.output_dir)
    sys.exit(0 if result.passes_qc else 1)


if __name__ == "__main__":
    import sys
    _cli()
