#!/usr/bin/env python3
"""GWAS Pipeline — PLINK2 QC + REGENIE two-step association testing.

Automates end-to-end genome-wide association studies: genotype QC via PLINK2,
whole-genome regression via REGENIE (Step 1 + Step 2), and post-GWAS
visualisation (Manhattan plot, QQ plot, lead variant extraction).

References:
    Mbatchou et al. (2021) Nature Genetics 53:1097-1103 (REGENIE)
    Chang et al. (2015) GigaScience 4:7 (PLINK2)
    Anderson et al. (2010) Nature Protocols 5:1564-1573 (QC thresholds)

Usage:
    python gwas_pipeline.py --demo --output /tmp/gwas_demo
    python gwas_pipeline.py --bed data --pheno pheno.txt --covar covar.txt \
        --trait-type bt --trait Y1 --output results/
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import shutil
import subprocess
import sys
from dataclasses import dataclass
from datetime import datetime, timezone
from pathlib import Path

SCRIPT_DIR = Path(__file__).resolve().parent

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

GWS_THRESHOLD = 5e-8
SUGGESTIVE_THRESHOLD = 1e-5


# ---------------------------------------------------------------------------
# Binary detection
# ---------------------------------------------------------------------------
def _find_binary(name: str) -> str | None:
    """Locate a binary on PATH or in the clawbio-gwas conda env."""
    path = shutil.which(name)
    if path:
        return path
    for conda_root in [Path(sys.prefix).parent, Path("/opt/anaconda3/envs")]:
        candidate = conda_root / "clawbio-gwas" / "bin" / name
        if candidate.exists():
            return str(candidate)
    return None


def _require_binaries() -> tuple[str, str]:
    """Return (plink2_path, regenie_path) or exit with installation instructions."""
    plink2 = _find_binary("plink2")
    regenie = _find_binary("regenie")
    missing = []
    if not plink2:
        missing.append("plink2")
    if not regenie:
        missing.append("regenie")
    if missing:
        print(
            f"ERROR: Required binaries not found: {', '.join(missing)}\n"
            "Install via conda:\n"
            "  CONDA_SUBDIR=osx-64 conda create -n clawbio-gwas "
            "-c conda-forge -c bioconda plink2 regenie -y\n"
            "  conda activate clawbio-gwas",
            file=sys.stderr,
        )
        sys.exit(1)
    return plink2, regenie


def _run_cmd(cmd: list[str], cwd: Path, log: list[str], label: str) -> subprocess.CompletedProcess:
    """Run a subprocess command, log it, and handle errors."""
    cmd_str = " ".join(str(c) for c in cmd)
    log.append(cmd_str)
    print(f"  [{label}] {cmd_str}")
    result = subprocess.run(cmd, capture_output=True, text=True, cwd=str(cwd), timeout=600)
    if result.returncode != 0:
        print(f"  WARNING: {label} exited with code {result.returncode}", file=sys.stderr)
        if result.stderr:
            for line in result.stderr.strip().splitlines()[-5:]:
                print(f"    {line}", file=sys.stderr)
    return result


# ---------------------------------------------------------------------------
# Phase 1: QC via PLINK2
# ---------------------------------------------------------------------------
@dataclass
class QCParams:
    geno: float = 0.02
    mind: float = 0.02
    maf: float = 0.01
    hwe: float = 1e-6
    ld_window: int = 1000
    ld_step: int = 100
    ld_r2: float = 0.9


def run_qc(
    plink2: str,
    bed_prefix: Path,
    work_dir: Path,
    params: QCParams,
    cmd_log: list[str],
) -> tuple[Path, Path]:
    """Run PLINK2 QC cascade and LD pruning. Returns (qc'd bed prefix, prune.in path)."""
    qc_dir = work_dir / "qc"
    qc_dir.mkdir(parents=True, exist_ok=True)

    qc_prefix = qc_dir / "qc_filtered"
    _run_cmd([
        plink2, "--bfile", str(bed_prefix),
        "--geno", str(params.geno),
        "--mind", str(params.mind),
        "--maf", str(params.maf),
        "--hwe", str(params.hwe),
        "--make-bed", "--out", str(qc_prefix),
    ], cwd=work_dir, log=cmd_log, label="QC")

    prune_prefix = qc_dir / "ld_prune"
    _run_cmd([
        plink2, "--bfile", str(qc_prefix),
        "--indep-pairwise", str(params.ld_window), str(params.ld_step), str(params.ld_r2),
        "--out", str(prune_prefix),
    ], cwd=work_dir, log=cmd_log, label="LD-prune")

    prune_in = prune_prefix.with_suffix(".prune.in")
    if not prune_in.exists():
        print("WARNING: LD pruning did not produce .prune.in file", file=sys.stderr)
        prune_in = None

    return qc_prefix, prune_in


def _parse_qc_log(qc_prefix: Path) -> dict:
    """Parse PLINK2 .log file for QC summary statistics."""
    log_file = qc_prefix.with_suffix(".log")
    stats = {"variants_before": 0, "variants_after": 0, "samples_before": 0, "samples_after": 0}
    if not log_file.exists():
        return stats
    text = log_file.read_text()
    for line in text.splitlines():
        if "variants loaded" in line.lower():
            parts = line.strip().split()
            if parts and parts[0].isdigit():
                stats["variants_before"] = int(parts[0])
        if "variants remaining" in line.lower() or ("variants and" in line.lower() and "pass" in line.lower()):
            parts = line.strip().split()
            if parts and parts[0].isdigit():
                stats["variants_after"] = int(parts[0])
        if "samples" in line.lower() and ("loaded" in line.lower() or "remaining" in line.lower()):
            parts = line.strip().split()
            if parts and parts[0].isdigit():
                if stats["samples_before"] == 0:
                    stats["samples_before"] = int(parts[0])
                else:
                    stats["samples_after"] = int(parts[0])
    if stats["samples_after"] == 0:
        stats["samples_after"] = stats["samples_before"]
    if stats["variants_after"] == 0:
        stats["variants_after"] = stats["variants_before"]
    return stats


# ---------------------------------------------------------------------------
# Phase 2: REGENIE Step 1
# ---------------------------------------------------------------------------
def run_step1(
    regenie: str,
    bed_prefix: Path,
    pheno_file: Path,
    covar_file: Path | None,
    trait_type: str,
    trait_name: str,
    prune_in: Path | None,
    work_dir: Path,
    bsize: int,
    cmd_log: list[str],
) -> Path:
    """Run REGENIE Step 1. Returns path to pred.list file."""
    step1_dir = work_dir / "step1"
    step1_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = step1_dir / "fit_out"

    cmd = [
        regenie, "--step", "1",
        "--bed", str(bed_prefix),
        "--phenoFile", str(pheno_file),
        "--phenoCol", trait_name,
        "--bsize", str(bsize),
        "--out", str(out_prefix),
    ]
    if trait_type == "bt":
        cmd.append("--bt")
    if covar_file:
        cmd.extend(["--covarFile", str(covar_file)])
    if prune_in and prune_in.exists():
        cmd.extend(["--extract", str(prune_in)])

    _run_cmd(cmd, cwd=work_dir, log=cmd_log, label="REGENIE-Step1")

    pred_list = out_prefix.parent / f"{out_prefix.name}_pred.list"
    if not pred_list.exists():
        print(f"ERROR: REGENIE Step 1 did not produce {pred_list}", file=sys.stderr)
        sys.exit(1)
    return pred_list


# ---------------------------------------------------------------------------
# Phase 3: REGENIE Step 2
# ---------------------------------------------------------------------------
def run_step2(
    regenie: str,
    bed_prefix: Path,
    bgen_file: Path | None,
    pheno_file: Path,
    covar_file: Path | None,
    trait_type: str,
    trait_name: str,
    pred_list: Path,
    work_dir: Path,
    bsize: int,
    cmd_log: list[str],
) -> Path:
    """Run REGENIE Step 2. Returns path to results .regenie file."""
    step2_dir = work_dir / "step2"
    step2_dir.mkdir(parents=True, exist_ok=True)
    out_prefix = step2_dir / "assoc_out"

    cmd = [
        regenie, "--step", "2",
        "--phenoFile", str(pheno_file),
        "--phenoCol", trait_name,
        "--pred", str(pred_list),
        "--bsize", str(bsize),
        "--out", str(out_prefix),
    ]

    if bgen_file and bgen_file.exists():
        cmd.extend(["--bgen", str(bgen_file)])
    else:
        cmd.extend(["--bed", str(bed_prefix)])

    if trait_type == "bt":
        cmd.extend(["--bt", "--firth", "--approx", "--pThresh", "0.01"])
    if covar_file:
        cmd.extend(["--covarFile", str(covar_file)])

    _run_cmd(cmd, cwd=work_dir, log=cmd_log, label="REGENIE-Step2")

    result_file = step2_dir / f"assoc_out_{trait_name}.regenie"
    if not result_file.exists():
        print(f"ERROR: REGENIE Step 2 did not produce {result_file}", file=sys.stderr)
        sys.exit(1)
    return result_file


# ---------------------------------------------------------------------------
# Phase 4: Post-GWAS processing
# ---------------------------------------------------------------------------
@dataclass
class GWASResult:
    chrom: str
    pos: int
    snp_id: str
    allele0: str
    allele1: str
    a1freq: float
    n: int
    beta: float
    se: float
    chisq: float
    log10p: float

    @property
    def pvalue(self) -> float:
        if self.log10p <= 0:
            return 1.0
        return 10 ** (-self.log10p)


def parse_regenie_output(result_file: Path) -> list[GWASResult]:
    """Parse a .regenie results file into a list of GWASResult."""
    results: list[GWASResult] = []
    with open(result_file) as fh:
        for line in fh:
            if line.startswith("CHROM") or line.startswith("#"):
                continue
            parts = line.strip().split()
            if len(parts) < 13:
                continue
            try:
                log10p = float(parts[12])
            except (ValueError, IndexError):
                continue
            if log10p <= 0:
                continue
            try:
                results.append(GWASResult(
                    chrom=parts[0],
                    pos=int(parts[1]),
                    snp_id=parts[2],
                    allele0=parts[3],
                    allele1=parts[4],
                    a1freq=float(parts[5]),
                    n=int(parts[7]),
                    beta=float(parts[9]),
                    se=float(parts[10]),
                    chisq=float(parts[11]),
                    log10p=log10p,
                ))
            except (ValueError, IndexError):
                continue
    return results


def compute_lambda_gc(results: list[GWASResult]) -> float:
    """Compute genomic inflation factor (lambda GC) from chi-squared statistics."""
    if not results:
        return 1.0
    chisq_values = sorted(r.chisq for r in results if r.chisq > 0)
    if not chisq_values:
        return 1.0
    median_chisq = chisq_values[len(chisq_values) // 2]
    return median_chisq / 0.4549364  # chi2 median for 1 df


def extract_lead_variants(results: list[GWASResult], threshold: float = GWS_THRESHOLD) -> list[GWASResult]:
    """Extract variants reaching genome-wide significance."""
    return sorted(
        [r for r in results if r.pvalue < threshold],
        key=lambda r: r.log10p,
        reverse=True,
    )


def manhattan_plot(results: list[GWASResult], output_path: Path) -> None:
    """Generate a Manhattan plot."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  WARNING: matplotlib/numpy not available — skipping Manhattan plot", file=sys.stderr)
        return

    chrom_order = [str(c) for c in range(1, 23)] + ["X", "Y"]
    chrom_to_idx = {c: i for i, c in enumerate(chrom_order)}

    x_positions = []
    y_values = []
    colours = []
    offset = 0
    chrom_centers = {}
    prev_chrom = None

    sorted_results = sorted(results, key=lambda r: (chrom_to_idx.get(r.chrom, 99), r.pos))

    for r in sorted_results:
        if r.chrom != prev_chrom:
            if prev_chrom is not None:
                offset += 50
            prev_chrom = r.chrom
        x = offset + r.pos
        x_positions.append(x)
        y_values.append(r.log10p)
        idx = chrom_to_idx.get(r.chrom, 0)
        colours.append("#2166ac" if idx % 2 == 0 else "#b2182b")
        if r.chrom not in chrom_centers:
            chrom_centers[r.chrom] = []
        chrom_centers[r.chrom].append(x)

    fig, ax = plt.subplots(figsize=(14, 5))
    ax.scatter(x_positions, y_values, c=colours, s=8, alpha=0.7, edgecolors="none")

    gws_line = -math.log10(GWS_THRESHOLD)
    sug_line = -math.log10(SUGGESTIVE_THRESHOLD)
    ax.axhline(y=gws_line, color="#d32f2f", linestyle="--", linewidth=0.8, label=f"GWS (P=5e-8)")
    ax.axhline(y=sug_line, color="#ff9800", linestyle="--", linewidth=0.8, label=f"Suggestive (P=1e-5)")

    tick_positions = []
    tick_labels = []
    for chrom, positions in chrom_centers.items():
        tick_positions.append(sum(positions) / len(positions))
        tick_labels.append(chrom)
    ax.set_xticks(tick_positions)
    ax.set_xticklabels(tick_labels, fontsize=8)

    ax.set_xlabel("Chromosome")
    ax.set_ylabel("-log10(P)")
    ax.set_title("Manhattan Plot")
    ax.legend(loc="upper right", fontsize=8)
    ax.set_xlim(min(x_positions) - 10, max(x_positions) + 10)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def qq_plot(results: list[GWASResult], output_path: Path, lambda_gc: float) -> None:
    """Generate a QQ plot with lambda GC annotation."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        print("  WARNING: matplotlib/numpy not available — skipping QQ plot", file=sys.stderr)
        return

    observed = sorted([r.log10p for r in results if r.log10p > 0], reverse=True)
    n = len(observed)
    if n == 0:
        return
    expected = [-math.log10((i + 0.5) / n) for i in range(n)]

    fig, ax = plt.subplots(figsize=(6, 6))
    ax.scatter(expected, observed, s=8, alpha=0.6, color="#2166ac", edgecolors="none")

    max_val = max(max(expected), max(observed)) + 0.5
    ax.plot([0, max_val], [0, max_val], "k--", linewidth=0.8)

    ax.set_xlabel("Expected -log10(P)")
    ax.set_ylabel("Observed -log10(P)")
    ax.set_title(f"QQ Plot (λGC = {lambda_gc:.3f})")
    ax.set_xlim(0, max_val)
    ax.set_ylim(0, max_val)
    ax.set_aspect("equal")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
def generate_report(
    results: list[GWASResult],
    lead_variants: list[GWASResult],
    lambda_gc: float,
    qc_stats: dict,
    output_dir: Path,
    cmd_log: list[str],
    trait_name: str,
    trait_type: str,
    demo: bool,
    input_label: str,
) -> None:
    """Generate all output files: report.md, result.json, tables/, reproducibility/."""
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    _write_summary_stats(results, output_dir / "tables" / "gwas_results.tsv")
    if lead_variants:
        _write_lead_variants(lead_variants, output_dir / "tables" / "lead_variants.tsv")

    _write_markdown_report(
        results, lead_variants, lambda_gc, qc_stats,
        output_dir, timestamp, trait_name, trait_type, demo, input_label,
    )
    _write_result_json(
        results, lead_variants, lambda_gc, qc_stats,
        output_dir, timestamp, trait_name, trait_type, demo,
    )
    _write_reproducibility(output_dir, cmd_log, timestamp)


def _write_summary_stats(results: list[GWASResult], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["CHROM", "POS", "ID", "A0", "A1", "A1FREQ", "N", "BETA", "SE", "CHISQ", "LOG10P", "P"])
        for r in results:
            writer.writerow([
                r.chrom, r.pos, r.snp_id, r.allele0, r.allele1,
                f"{r.a1freq:.6f}", r.n, f"{r.beta:.6f}", f"{r.se:.6f}",
                f"{r.chisq:.6f}", f"{r.log10p:.6f}", f"{r.pvalue:.2e}",
            ])


def _write_lead_variants(variants: list[GWASResult], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        writer.writerow(["CHROM", "POS", "ID", "A0", "A1", "BETA", "SE", "P", "LOG10P"])
        for r in variants:
            writer.writerow([
                r.chrom, r.pos, r.snp_id, r.allele0, r.allele1,
                f"{r.beta:.6f}", f"{r.se:.6f}", f"{r.pvalue:.2e}", f"{r.log10p:.4f}",
            ])


def _write_markdown_report(
    results, lead_variants, lambda_gc, qc_stats,
    output_dir, timestamp, trait_name, trait_type, demo, input_label,
) -> None:
    lines: list[str] = []
    lines.append("# GWAS Report — PLINK2 QC + REGENIE Association")
    lines.append("")
    lines.append(f"**Generated**: {timestamp}")
    lines.append(f"**Input**: {input_label}")
    lines.append(f"**Trait**: {trait_name} ({'binary' if trait_type == 'bt' else 'quantitative'})")
    lines.append(f"**Mode**: {'Demo (REGENIE example data)' if demo else 'Live'}")
    lines.append("")

    lines.append("## QC Summary")
    lines.append("")
    lines.append(f"- Variants before QC: {qc_stats.get('variants_before', 'N/A')}")
    lines.append(f"- Variants after QC: {qc_stats.get('variants_after', 'N/A')}")
    lines.append(f"- Samples before QC: {qc_stats.get('samples_before', 'N/A')}")
    lines.append(f"- Samples after QC: {qc_stats.get('samples_after', 'N/A')}")
    lines.append("")

    lines.append("## Association Results Summary")
    lines.append("")
    lines.append(f"- Total variants tested: {len(results)}")
    lines.append(f"- Lambda GC: {lambda_gc:.4f}")
    gws = [r for r in results if r.pvalue < GWS_THRESHOLD]
    sug = [r for r in results if r.pvalue < SUGGESTIVE_THRESHOLD]
    lines.append(f"- Genome-wide significant (P < 5e-8): {len(gws)}")
    lines.append(f"- Suggestive (P < 1e-5): {len(sug)}")
    lines.append("")

    if lead_variants:
        lines.append("## Lead Variants (P < 5e-8)")
        lines.append("")
        lines.append("| CHR | POS | ID | A1 | BETA | SE | P |")
        lines.append("|-----|-----|----|----|------|----|---|")
        for r in lead_variants:
            lines.append(f"| {r.chrom} | {r.pos} | {r.snp_id} | {r.allele1} | {r.beta:.4f} | {r.se:.4f} | {r.pvalue:.2e} |")
        lines.append("")

    if not lead_variants:
        lines.append("## Top Variants (by P-value)")
        lines.append("")
        lines.append("| CHR | POS | ID | A1 | BETA | SE | LOG10P |")
        lines.append("|-----|-----|----|----|------|----|--------|")
        top = sorted(results, key=lambda r: r.log10p, reverse=True)[:10]
        for r in top:
            lines.append(f"| {r.chrom} | {r.pos} | {r.snp_id} | {r.allele1} | {r.beta:.4f} | {r.se:.4f} | {r.log10p:.4f} |")
        lines.append("")

    lines.append("## Methodology")
    lines.append("")
    lines.append("Genotype QC was performed with PLINK2 (Chang et al., GigaScience 2015). ")
    lines.append("Association testing used REGENIE's two-step whole-genome regression ")
    lines.append("(Mbatchou et al., Nature Genetics 2021): Step 1 fits a ridge regression ")
    lines.append("model with LOCO predictions; Step 2 tests each variant conditional on ")
    lines.append("the Step 1 polygenic background.")
    if trait_type == "bt":
        lines.append(" Firth logistic regression with fast approximation was used for binary traits.")
    lines.append("")

    lines.append("---")
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def _write_result_json(
    results, lead_variants, lambda_gc, qc_stats,
    output_dir, timestamp, trait_name, trait_type, demo,
) -> None:
    result = {
        "tool": "ClawBio GWAS Pipeline",
        "version": "0.1.0",
        "method": "PLINK2 QC + REGENIE two-step regression",
        "references": [
            "Mbatchou et al. (2021) PMID 34017140",
            "Chang et al. (2015) PMID 25722852",
        ],
        "timestamp": timestamp,
        "mode": "demo" if demo else "live",
        "trait": trait_name,
        "trait_type": trait_type,
        "qc_summary": qc_stats,
        "total_variants_tested": len(results),
        "lambda_gc": round(lambda_gc, 4),
        "gws_hits": len([r for r in results if r.pvalue < GWS_THRESHOLD]),
        "suggestive_hits": len([r for r in results if r.pvalue < SUGGESTIVE_THRESHOLD]),
        "top_10": [
            {
                "chrom": r.chrom, "pos": r.pos, "id": r.snp_id,
                "beta": round(r.beta, 6), "se": round(r.se, 6),
                "log10p": round(r.log10p, 4), "p": f"{r.pvalue:.2e}",
            }
            for r in sorted(results, key=lambda r: r.log10p, reverse=True)[:10]
        ],
        "disclaimer": DISCLAIMER,
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2, ensure_ascii=False), encoding="utf-8",
    )


def _write_reproducibility(output_dir: Path, cmd_log: list[str], timestamp: str) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(exist_ok=True)
    (repro_dir / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        f"# GWAS Pipeline — reproducibility log\n"
        f"# Generated: {timestamp}\n\n"
        + "\n".join(cmd_log) + "\n",
        encoding="utf-8",
    )

    plink2 = _find_binary("plink2")
    regenie = _find_binary("regenie")
    plink_ver = ""
    regenie_ver = ""
    if plink2:
        try:
            plink_ver = subprocess.run([plink2, "--version"], capture_output=True, text=True, timeout=10).stdout.strip()
        except Exception:
            pass
    if regenie:
        try:
            regenie_ver = subprocess.run([regenie, "--version"], capture_output=True, text=True, timeout=10).stdout.strip()
        except Exception:
            pass

    versions = {
        "plink2": plink_ver,
        "regenie": regenie_ver,
        "python": sys.version,
        "generated": timestamp,
    }
    (repro_dir / "software_versions.json").write_text(
        json.dumps(versions, indent=2), encoding="utf-8",
    )


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def run_pipeline(
    bed_prefix: Path,
    bgen_file: Path | None,
    pheno_file: Path,
    covar_file: Path | None,
    trait_type: str,
    trait_name: str,
    output_dir: Path,
    qc_params: QCParams,
    step1_bsize: int = 100,
    step2_bsize: int = 200,
    demo: bool = False,
) -> dict:
    """Run the full GWAS pipeline: QC -> Step 1 -> Step 2 -> Post-GWAS."""
    plink2, regenie = _require_binaries()
    cmd_log: list[str] = []

    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    print("[GWAS] Phase 1: Genotype QC (PLINK2)")
    qc_prefix, prune_in = run_qc(plink2, bed_prefix, output_dir, qc_params, cmd_log)
    qc_stats = _parse_qc_log(qc_prefix)

    print("[GWAS] Phase 2: REGENIE Step 1 (whole-genome regression)")
    pred_list = run_step1(
        regenie, qc_prefix, pheno_file, covar_file,
        trait_type, trait_name, prune_in, output_dir, step1_bsize, cmd_log,
    )

    print("[GWAS] Phase 3: REGENIE Step 2 (association testing)")
    result_file = run_step2(
        regenie, bed_prefix, bgen_file, pheno_file, covar_file,
        trait_type, trait_name, pred_list, output_dir, step2_bsize, cmd_log,
    )

    print("[GWAS] Phase 4: Post-GWAS processing")
    results = parse_regenie_output(result_file)
    lambda_gc = compute_lambda_gc(results)
    lead_variants = extract_lead_variants(results)

    print(f"  Variants tested: {len(results)}")
    print(f"  Lambda GC: {lambda_gc:.4f}")
    print(f"  GWS hits (P < 5e-8): {len(lead_variants)}")

    manhattan_plot(results, output_dir / "figures" / "manhattan.png")
    qq_plot(results, output_dir / "figures" / "qq_plot.png", lambda_gc)

    input_label = str(bed_prefix)
    generate_report(
        results, lead_variants, lambda_gc, qc_stats, output_dir,
        cmd_log, trait_name, trait_type, demo, input_label,
    )

    print(f"[GWAS] Report written to: {output_dir / 'report.md'}")
    return {
        "total_variants": len(results),
        "lambda_gc": lambda_gc,
        "gws_hits": len(lead_variants),
        "output_dir": str(output_dir),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="GWAS Pipeline — PLINK2 QC + REGENIE association")
    parser.add_argument("--bed", type=str, help="PLINK BED prefix (without .bed extension)")
    parser.add_argument("--bgen", type=str, help="BGEN file for Step 2 (optional)")
    parser.add_argument("--pheno", type=str, help="Phenotype file (FID IID trait columns)")
    parser.add_argument("--covar", type=str, help="Covariate file (FID IID covariate columns)")
    parser.add_argument("--trait-type", choices=["bt", "qt"], default="bt", help="bt=binary, qt=quantitative")
    parser.add_argument("--trait", type=str, default="Y1", help="Phenotype column name")
    parser.add_argument("--output", type=str, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with REGENIE example data")
    parser.add_argument("--geno", type=float, default=0.02, help="Variant missingness threshold")
    parser.add_argument("--mind", type=float, default=0.02, help="Sample missingness threshold")
    parser.add_argument("--maf", type=float, default=0.01, help="Minor allele frequency threshold")
    parser.add_argument("--hwe", type=float, default=1e-6, help="HWE p-value threshold")

    args = parser.parse_args()

    if not args.output:
        parser.error("--output is required")

    if args.demo:
        data_dir = SCRIPT_DIR / "example_data"
        bed_prefix = data_dir / "example"
        bgen_file = data_dir / "example.bgen"
        pheno_file = data_dir / "phenotype_bin.txt"
        covar_file = data_dir / "covariates.txt"
        trait_type = "bt"
        trait_name = "Y1"
    else:
        if not args.bed or not args.pheno:
            parser.error("--bed and --pheno are required (or use --demo)")
        bed_prefix = Path(args.bed)
        bgen_file = Path(args.bgen) if args.bgen else None
        pheno_file = Path(args.pheno)
        covar_file = Path(args.covar) if args.covar else None
        trait_type = args.trait_type
        trait_name = args.trait

    output_dir = Path(args.output)
    if output_dir.exists() and any(output_dir.iterdir()):
        print(f"WARNING: Output directory '{output_dir}' is not empty.", file=sys.stderr)

    qc_params = QCParams(geno=args.geno, mind=args.mind, maf=args.maf, hwe=args.hwe)

    print(f"[GWAS] Starting pipeline ({'demo' if args.demo else 'live'} mode)")
    print(f"[GWAS] Trait: {trait_name} ({'binary' if trait_type == 'bt' else 'quantitative'})")

    summary = run_pipeline(
        bed_prefix=bed_prefix,
        bgen_file=bgen_file,
        pheno_file=pheno_file,
        covar_file=covar_file,
        trait_type=trait_type,
        trait_name=trait_name,
        output_dir=output_dir,
        qc_params=qc_params,
        demo=args.demo,
    )

    print(f"[GWAS] Pipeline complete: {summary['total_variants']} variants tested, "
          f"λGC={summary['lambda_gc']:.4f}, {summary['gws_hits']} GWS hits")


if __name__ == "__main__":
    main()
