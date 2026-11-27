#!/usr/bin/env python3
"""
ClawBio Seq Wrangler
NGS QC, alignment, and BAM processing pipeline.
Wraps FastQC, BWA/Bowtie2/Minimap2, SAMtools, and MultiQC.

Usage:
    # Single sample SE
    python seq_wrangler.py --r1 sample.fastq.gz --index /ref/hg38 --output results/

    # Single sample PE
    python seq_wrangler.py --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz --index /ref/hg38 --output results/

    # Batch mode (samplesheet CSV)
    python seq_wrangler.py --samplesheet samples.csv --index /ref/hg38 --output results/

    # Demo mode (no external tools needed)
    python seq_wrangler.py --demo --output demo_results/
"""
from __future__ import annotations

import argparse
import csv
import json
import os
import re
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Dict, List, Optional, Tuple

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)
from clawbio.common.textio import write_text_lf_atomic  # noqa: E402


# Constants

DEFAULT_MAPQ = 20
DEFAULT_THREADS = os.cpu_count() or 4
SUPPORTED_ALIGNERS = ("bwa", "bowtie2", "minimap2")
SAFETY_DISCLAIMER = (
    "ClawBio Seq Wrangler is a research and educational tool. "
    "Results must be validated before use in clinical or production settings."
)

DEMO_FLAGSTAT = {
    "CTRL_REP1": {
        "total": 10_000_000, "mapped": 9_750_000, "mapped_pct": 97.5,
        "duplicates": 850_000, "dup_pct": 8.7,
        "paired": True, "properly_paired": 9_500_000,
    },
    "TREAT_REP1": {
        "total": 8_500_000, "mapped": 8_200_000, "mapped_pct": 96.5,
        "duplicates": 600_000, "dup_pct": 7.3,
        "paired": False, "properly_paired": 0,
    },
}

DEMO_COVERAGE = {
    "CTRL_REP1": [
        {"chr": "chr1", "mean_depth": 28.4, "breadth": 0.97},
        {"chr": "chr2", "mean_depth": 27.9, "breadth": 0.96},
        {"chr": "chrX", "mean_depth": 14.2, "breadth": 0.91},
        {"chr": "chrM", "mean_depth": 1842.0, "breadth": 1.00},
    ],
    "TREAT_REP1": [
        {"chr": "chr1", "mean_depth": 24.1, "breadth": 0.95},
        {"chr": "chr2", "mean_depth": 23.8, "breadth": 0.94},
        {"chr": "chrX", "mean_depth": 23.5, "breadth": 0.94},
        {"chr": "chrM", "mean_depth": 1560.0, "breadth": 1.00},
    ],
}

DEMO_INSERT_SIZE = {
    "CTRL_REP1": {"mean": 312.4, "std": 48.2, "median": 305.0},
}

DEMO_SAMPLES = [
    {"sample": "CTRL_REP1", "fastq_1": "demo_ctrl_R1.fastq.gz", "fastq_2": "demo_ctrl_R2.fastq.gz"},
    {"sample": "TREAT_REP1", "fastq_1": "demo_treat_R1.fastq.gz", "fastq_2": ""},
]



# Utility

def run_command(
    cmd: List[str],
    description: str,
    critical: bool = True,
    capture_output: bool = True,
) -> subprocess.CompletedProcess:
    print(f"  [{description}] {' '.join(str(c) for c in cmd[:7])}...")
    try:
        result = subprocess.run(
            [str(c) for c in cmd],
            capture_output=capture_output,
            text=True if capture_output else False,
            timeout=7200,
            check=False,
        )
        if result.returncode != 0:
            stderr_tail = ""
            if capture_output and result.stderr:
                lines = result.stderr.strip().split("\n")
                stderr_tail = "\n".join(f"  STDERR: {l}" for l in lines[-5:])
            if critical:
                msg = f"{description} failed (exit {result.returncode})."
                if stderr_tail:
                    msg += f"\n{stderr_tail}"
                raise RuntimeError(msg)
            else:
                print(f"  WARNING: {description} returned non-zero exit {result.returncode}", file=sys.stderr)
                if stderr_tail:
                    print(stderr_tail, file=sys.stderr)
        return result
    except FileNotFoundError:
        print(f"  ERROR: '{cmd[0]}' not found. Is it installed and on PATH?", file=sys.stderr)
        sys.exit(1)
    except subprocess.TimeoutExpired:
        print(f"  ERROR: {description} timed out after 2 hours.", file=sys.stderr)
        sys.exit(1)


def check_tool(name: str) -> bool:
    return shutil.which(name) is not None


def detect_threads() -> int:
    count = os.cpu_count()
    return max(1, count - 1) if count and count > 1 else 1


def load_samplesheet(path: Path) -> List[Dict[str, str]]:
    samples = []
    with open(path, newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        if not {"sample", "fastq_1"}.issubset(set(reader.fieldnames or [])):
            raise ValueError(f"Samplesheet must have columns: sample, fastq_1. Found: {reader.fieldnames}")
        for i, row in enumerate(reader, start=2):
            sample = row["sample"].strip()
            fq1 = row["fastq_1"].strip()
            fq2 = row.get("fastq_2", "").strip()
            if not sample or not fq1:
                print(f"  WARNING: Skipping empty row at line {i}.", file=sys.stderr)
                continue
            samples.append({"sample": sample, "fastq_1": fq1, "fastq_2": fq2})
    if not samples:
        raise ValueError(f"No valid samples found in: {path}")
    return samples


def validate_input_files(samples: List[Dict[str, str]]) -> None:
    for s in samples:
        fq1 = Path(s["fastq_1"])
        if not fq1.exists():
            raise FileNotFoundError(f"FASTQ file not found for sample {s['sample']}: {fq1}")
        if s["fastq_2"]:
            fq2 = Path(s["fastq_2"])
            if not fq2.exists():
                raise FileNotFoundError(f"FASTQ R2 not found for sample {s['sample']}: {fq2}")


def validate_tools(aligner: str, trim: bool, run_multiqc_flag: bool) -> None:
    if not check_tool("samtools"):
        raise EnvironmentError("Required tool missing: samtools")
    if not check_tool(aligner):
        raise EnvironmentError(f"Required aligner missing: {aligner}")
    if trim and not check_tool("fastp"):
        print("  WARNING: --trim requested but fastp is not installed; trimming will be skipped.", file=sys.stderr)
    if run_multiqc_flag and not check_tool("multiqc"):
        print("  WARNING: --run-multiqc requested but multiqc is not installed; MultiQC will be skipped.", file=sys.stderr)



# Step 1 — FastQC

def run_fastqc(fastqs: List[Path], output_dir: Path, threads: int) -> Optional[Path]:
    if not check_tool("fastqc"):
        print("  WARNING: fastqc not found — skipping.", file=sys.stderr)
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    cmd = ["fastqc", "-t", str(threads), "-o", str(output_dir)] + [str(f) for f in fastqs]
    run_command(cmd, "FastQC", critical=False)
    return output_dir



# Step 2 — fastp trimming (optional)

def run_fastp_SE(fastq: Path, output_dir: Path, threads: int) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = fastq.name.replace(".fastq.gz", "").replace(".fq.gz", "")
    trimmed = output_dir / f"{stem}_trimmed.fastq.gz"
    run_command([
        "fastp", "-i", str(fastq), "-o", str(trimmed),
        "-j", str(output_dir / f"{stem}_fastp.json"),
        "-h", str(output_dir / f"{stem}_fastp.html"),
        "-w", str(threads),
    ], f"fastp SE: {fastq.name}")
    return trimmed


def run_fastp_PE(r1: Path, r2: Path, output_dir: Path, threads: int) -> Tuple[Path, Path]:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = r1.name.replace(".fastq.gz", "").replace(".fq.gz", "").replace("_R1", "")
    t1 = output_dir / r1.name.replace(".fastq.gz", "_trimmed.fastq.gz").replace(".fq.gz", "_trimmed.fq.gz")
    t2 = output_dir / r2.name.replace(".fastq.gz", "_trimmed.fastq.gz").replace(".fq.gz", "_trimmed.fq.gz")
    run_command([
        "fastp", "-i", str(r1), "-I", str(r2),
        "-o", str(t1), "-O", str(t2),
        "-j", str(output_dir / f"{stem}_fastp.json"),
        "-h", str(output_dir / f"{stem}_fastp.html"),
        "-w", str(threads),
    ], f"fastp PE: {r1.name}")
    return t1, t2



# Step 3 — Alignment

def align_SE(sample: str, fastq: Path, index: str, aligner: str, output_dir: Path, threads: int) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    sam = output_dir / f"{sample}.sam"
    if aligner == "bwa":
        cmd = ["bwa", "mem", "-t", str(threads), index, str(fastq)]
        with open(sam, "w", encoding="utf-8") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"bwa failed: {result.stderr[-500:]}")
    elif aligner == "bowtie2":
        run_command(["bowtie2", "-p", str(threads), "-x", index, "-U", str(fastq), "-S", str(sam)],
                    f"bowtie2 SE: {sample}")
    elif aligner == "minimap2":
        cmd = ["minimap2", "-ax", "sr", "-t", str(threads), index, str(fastq)]
        with open(sam, "w", encoding="utf-8") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {result.stderr[-500:]}")
    else:
        raise ValueError(f"Unsupported aligner: {aligner}")
    return sam


def align_PE(sample: str, r1: Path, r2: Path, index: str, aligner: str, output_dir: Path, threads: int) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    sam = output_dir / f"{sample}.sam"
    if aligner == "bwa":
        cmd = ["bwa", "mem", "-t", str(threads), index, str(r1), str(r2)]
        with open(sam, "w", encoding="utf-8") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"bwa failed: {result.stderr[-500:]}")
    elif aligner == "bowtie2":
        run_command(["bowtie2", "-p", str(threads), "-x", index, "-1", str(r1), "-2", str(r2), "-S", str(sam)],
                    f"bowtie2 PE: {sample}")
    elif aligner == "minimap2":
        cmd = ["minimap2", "-ax", "sr", "-t", str(threads), index, str(r1), str(r2)]
        with open(sam, "w", encoding="utf-8") as out_f:
            result = subprocess.run(cmd, stdout=out_f, stderr=subprocess.PIPE, text=True, check=False)
        if result.returncode != 0:
            raise RuntimeError(f"minimap2 failed: {result.stderr[-500:]}")
    else:
        raise ValueError(f"Unsupported aligner: {aligner}")
    return sam



# Step 4 — BAM processing

def sam_to_sorted_bam(sam: Path, output_dir: Path, mapq: int, threads: int, keep_sam: bool = False) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    sample = sam.stem
    tmp_bam = output_dir / f"{sample}.tmp.bam"
    sorted_bam = output_dir / f"{sample}_sorted.tmp.bam"

    run_command(["samtools", "view", "-@", str(threads), "-b", "-q", str(mapq), str(sam), "-o", str(tmp_bam)],
                f"samtools view MAPQ≥{mapq}: {sample}")
    run_command(["samtools", "sort", "-@", str(threads), "-m", "2G", "-n", "-o", str(sorted_bam), str(tmp_bam)],
                f"samtools sort -n: {sample}")

    tmp_bam.unlink(missing_ok=True)
    if not keep_sam:
        sam.unlink(missing_ok=True)
    return sorted_bam


def run_fixmate(name_sorted_bam: Path, threads: int) -> Path:
    out = name_sorted_bam.parent / name_sorted_bam.name.replace("_sorted.tmp", "_fixmate.tmp")
    run_command(["samtools", "fixmate", "-@", str(threads), "-m", str(name_sorted_bam), str(out)],
                f"samtools fixmate: {name_sorted_bam.stem}")
    name_sorted_bam.unlink(missing_ok=True)
    return out


def run_coordinate_sort(bam: Path, threads: int, sample: str) -> Path:
    out = bam.parent / f"{sample}_resort.tmp.bam"
    run_command(["samtools", "sort", "-@", str(threads), "-m", "2G", "-o", str(out), str(bam)],
                f"samtools coordinate sort: {sample}")
    bam.unlink(missing_ok=True)
    return out


def run_markdup(bam: Path, output_dir: Path, threads: int, remove_dups: bool = False) -> Path:
    output_dir.mkdir(parents=True, exist_ok=True)
    stem = bam.stem.replace("_fixmate.tmp", "").replace("_sorted.tmp", "").replace("_resort.tmp", "")
    out_bam = output_dir / f"{stem}_sorted.bam"
    cmd = ["samtools", "markdup", "-@", str(threads)]
    if remove_dups:
        cmd.append("-r")
    cmd += [str(bam), str(out_bam)]
    run_command(cmd, f"samtools markdup: {stem}")
    bam.unlink(missing_ok=True)
    return out_bam


def index_bam(bam: Path, threads: int) -> None:
    run_command(["samtools", "index", "-@", str(threads), str(bam)], f"samtools index: {bam.name}")


def get_flagstat(bam: Path) -> Dict:
    result = run_command(["samtools", "flagstat", str(bam)], f"samtools flagstat: {bam.name}")
    stats: Dict = {"raw": result.stdout}
    for line in result.stdout.splitlines():
        if "in total" in line:
            stats["total"] = int(line.split(" + ")[0].strip())
        elif "mapped (" in line and "primary mapped" not in line:
            parts = line.split(" + ")
            stats["mapped"] = int(parts[0].strip())
            pct_str = line.split("(")[-1].split("%")[0].strip()
            try:
                stats["mapped_pct"] = float(pct_str)
            except ValueError:
                stats["mapped_pct"] = 0.0
        elif "duplicates" in line:
            stats["duplicates"] = int(line.split(" + ")[0].strip())
        elif "properly paired" in line:
            stats["properly_paired"] = int(line.split(" + ")[0].strip())
    return stats


def get_coverage(bam: Path) -> List[Dict]:
    result = run_command(["samtools", "coverage", str(bam)], f"samtools coverage: {bam.name}")
    rows: List[Dict] = []
    lines = result.stdout.strip().splitlines()
    if not lines:
        return rows
    headers = lines[0].lstrip("#").split("\t")
    for line in lines[1:]:
        fields = line.split("\t")
        row = dict(zip(headers, fields))
        covbases = float(row.get("covbases", 0) or 0)
        endpos = float(row.get("endpos", 1) or 1)
        rows.append({
            "chr": row.get("rname", ""),
            "mean_depth": float(row.get("meandepth", 0) or 0),
            "breadth": covbases / max(endpos, 1),
        })
    return rows


def get_insert_size(bam: Path) -> Optional[Dict]:
    result = run_command(["samtools", "stats", str(bam)], f"samtools stats: {bam.name}", critical=False)
    if result.returncode != 0:
        return None

    # Weighted moments to avoid huge memory usage
    total_n = 0
    sum_x = 0.0
    sum_x2 = 0.0
    weighted_values: List[Tuple[int, int]] = []

    for line in result.stdout.splitlines():
        if line.startswith("IS\t"):
            parts = line.split("\t")
            if len(parts) >= 3:
                try:
                    size = int(parts[1])
                    count = int(parts[2])
                except ValueError:
                    continue
                if count <= 0:
                    continue
                weighted_values.append((size, count))
                total_n += count
                sum_x += size * count
                sum_x2 += (size ** 2) * count

    if total_n == 0:
        return None

    mean = sum_x / total_n
    variance = max((sum_x2 / total_n) - (mean ** 2), 0.0)
    std = variance ** 0.5

    # Exact weighted median from histogram
    weighted_values.sort(key=lambda x: x[0])
    midpoint = (total_n + 1) / 2
    acc = 0
    median = weighted_values[-1][0]
    for size, count in weighted_values:
        acc += count
        if acc >= midpoint:
            median = size
            break

    return {"mean": round(mean, 2), "std": round(std, 2), "median": float(median)}


def process_sample(
    sample: str,
    fq1: Path,
    fq2: Optional[Path],
    index: str,
    aligner: str,
    output_dir: Path,
    threads: int,
    mapq: int,
    remove_dups: bool,
    trim: bool,
    keep_sam: bool
) -> Dict:
    is_pe = fq2 is not None
    print(f"\n{'='*50}")
    print(f"  Sample: {sample} ({'PE' if is_pe else 'SE'})")
    print(f"{'='*50}")

    align_dir = output_dir / "alignment"
    bam_dir = output_dir / "bam"
    trim_dir = output_dir / "trimmed"

    if trim:
        if not check_tool("fastp"):
            print("  WARNING: fastp not found — skipping trimming.", file=sys.stderr)
        else:
            if is_pe:
                fq1, fq2 = run_fastp_PE(fq1, fq2, trim_dir, threads)
            else:
                fq1 = run_fastp_SE(fq1, trim_dir, threads)

    if is_pe:
        sam = align_PE(sample, fq1, fq2, index, aligner, align_dir, threads)
    else:
        sam = align_SE(sample, fq1, index, aligner, align_dir, threads)

    name_sorted = sam_to_sorted_bam(sam, bam_dir, mapq, threads, keep_sam=keep_sam)
    fixmated = run_fixmate(name_sorted, threads)
    resorted = run_coordinate_sort(fixmated, threads, sample)
    final_bam = run_markdup(resorted, bam_dir, threads, remove_dups=remove_dups)
    index_bam(final_bam, threads)

    flagstat = get_flagstat(final_bam)
    coverage = get_coverage(final_bam)
    insert_size = get_insert_size(final_bam) if is_pe else None

    return {
        "sample": sample,
        "is_pe": is_pe,
        "bam": str(final_bam),
        "flagstat": flagstat,
        "coverage": coverage,
        "insert_size": insert_size,
    }



# Step 5 — MultiQC

def run_multiqc(search_dir: Path, output_dir: Path) -> Optional[Path]:
    if not check_tool("multiqc"):
        print("  WARNING: multiqc not found — skipping.", file=sys.stderr)
        return None
    output_dir.mkdir(parents=True, exist_ok=True)
    run_command(["multiqc", str(search_dir), "-o", str(output_dir), "--force"], "MultiQC", critical=False)
    report = output_dir / "multiqc_report.html"
    return report if report.exists() else None



# Step 6 — Report

def generate_report(
    output_dir: Path,
    sample_stats: List[Dict],
    input_files: List[Path],
    is_demo: bool,
    aligner: str,
    mapq: int,
    remove_dups: bool,
) -> Path:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    report_path = output_dir / "report.md"

    sample_sections = []
    for s in sample_stats:
        name = s["sample"]
        fs = s["flagstat"]
        read_type = "Paired-end" if s["is_pe"] else "Single-end"

        def fmt_int(v: object) -> str:
            return f"{v:,}" if isinstance(v, int) else str(v)

        flagstat_lines = [
            f"| Total reads | {fmt_int(fs.get('total', 'N/A'))} |",
            f"| Mapped reads | {fmt_int(fs.get('mapped', 'N/A'))} ({fs.get('mapped_pct', 0):.1f}%) |",
            f"| Duplicates | {fmt_int(fs.get('duplicates', 'N/A'))} |",
        ]
        if s["is_pe"]:
            flagstat_lines.append(f"| Properly paired | {fmt_int(fs.get('properly_paired', 'N/A'))} |")

        flagstat_table = "\n".join(["| Metric | Value |", "|--------|-------|"] + flagstat_lines)

        cov_rows = s.get("coverage", [])[:5]
        cov_lines = [
            f"| {r['chr']} | {r['mean_depth']:.1f}x | {r['breadth']*100:.1f}% |"
            for r in cov_rows
        ]
        cov_table = "\n".join([
            "| Chromosome | Mean depth | Breadth |",
            "|------------|------------|---------|",
        ] + (cov_lines if cov_lines else ["| — | — | — |"]))

        is_section = ""
        if s["is_pe"] and s.get("insert_size"):
            ins = s["insert_size"]
            is_section = (
                f"\n**Insert size**: mean {ins['mean']} bp | "
                f"median {ins['median']} bp | std {ins['std']} bp\n"
            )

        sample_sections.append(f"""
### {name} ({read_type})

**Alignment statistics**

{flagstat_table}

**Coverage (top chromosomes)**

{cov_table}
{is_section}
**BAM**: `{s['bam']}`
""")

    checksum_lines = []
    for fp in input_files:
        if fp.exists():
            checksum_lines.append(f"- `{fp.name}`: `{sha256_file(fp)}`")
        else:
            checksum_lines.append(f"- `{fp.name}`: demo data (no file on disk)")

    report = f"""# Seq Wrangler Report

**Date**: {now}  
**Tool**: ClawBio Seq Wrangler v0.1.0  
**Mode**: {"Demo (synthetic data)" if is_demo else "Full pipeline"}

---

## Run parameters

| Parameter | Value |
|-----------|-------|
| Aligner | {aligner} |
| MAPQ threshold | ≥{mapq} |
| Duplicate handling | {"Remove (-r)" if remove_dups else "Mark only"} |
| Samples processed | {len(sample_stats)} |

---

## Results per sample
{"".join(sample_sections)}

---

## Methods

- **QC**: FastQC per-sample quality assessment
- **Trimming**: fastp (if --trim flag used)
- **Alignment**: {aligner} with MAPQ ≥ {mapq} filter
- **BAM processing**: samtools sort -n → fixmate → coordinate-sort → markdup{" (-r remove)" if remove_dups else " (mark only)"} → index
- **Coverage**: samtools coverage (per-chromosome mean depth and breadth)
- **Insert size**: samtools stats IS histogram (paired-end only)
- **MultiQC**: aggregated QC report across samples (if --run-multiqc)

## Input checksums

{"".join(checksum_lines) if checksum_lines else "- demo data"}

## Disclaimer

> {SAFETY_DISCLAIMER}
"""
    report_path.write_text(report, encoding="utf-8")
    print(f"  Saved: {report_path}")
    return report_path



# Step 7 — Reproducibility bundle

PIPELINE_TOOLS = ["samtools", "bwa", "bowtie2", "minimap2", "fastqc", "fastp", "multiqc"]

# A real run records the tool versions it replays with; the demo calls no tools.
TOOL_PREFLIGHT = [
    "# Tool versions",
    "samtools --version | head -1 || true",
    "bwa 2>&1 | head -1 || true",
    "bowtie2 --version | head -1 || true",
    "minimap2 --version | head -1 || true",
    "fastqc --version || true",
    "fastp --version || true",
    "multiqc --version || true",
]

# Flags replayed when they differ from the parser default, in this order.
_VALUE_OPTIONS = ["aligner", "genome_build", "mapq"]
_SWITCHES = ["trim", "remove_duplicates", "run_fastqc", "run_multiqc", "keep_sam"]


def _flag(name: str) -> str:
    return "--" + name.replace("_", "-")


def _quoted_path(value: str) -> str:
    """A path for commands.sh: under $CLAWBIO_ROOT when inside the checkout, else absolute."""
    path = Path(value).expanduser().resolve()
    if path.is_relative_to(_PROJECT_ROOT):
        return f'"$CLAWBIO_ROOT"/{shlex.quote(path.relative_to(_PROJECT_ROOT).as_posix())}'
    return shlex.quote(str(path))


def _replay_args(args: argparse.Namespace, output_dir: Path, is_demo: bool) -> list[str | ReproPath]:
    defaults = build_parser()
    replay: list[str | ReproPath] = ["--demo"] if is_demo else []
    if not is_demo:
        if args.samplesheet:
            replay += ["--samplesheet", _quoted_path(args.samplesheet)]
        else:
            replay += ["--r1", _quoted_path(args.r1)]
            if args.r2:
                replay += ["--r2", _quoted_path(args.r2)]
            if getattr(args, "sample_name", None):
                replay += ["--sample-name", shlex.quote(args.sample_name)]
        # An aligner index is a path prefix, not a file, so it is quoted as given.
        replay += ["--index", shlex.quote(str(Path(args.index).expanduser().resolve()))]
    for name in _VALUE_OPTIONS:
        value = getattr(args, name)
        # A real run always names its aligner; anything else only when not the default.
        if (not is_demo and name == "aligner") or value != defaults.get_default(name):
            replay += [_flag(name), shlex.quote(str(value))]
    if not is_demo:  # thread count is per machine, but a real run replays with the same
        replay += ["--threads", str(args.threads)]
    for name in _SWITCHES:
        if getattr(args, name, False):
            replay.append(_flag(name))
    return replay + ["--output", ReproPath(output_dir, anchor="output_dir")]


def write_reproducibility_bundle(output_dir: Path, input_files: List[Path], args: argparse.Namespace, is_demo: bool) -> None:
    """Write reproducibility/ with the shared helpers once every other output exists.

    inputs.json holds the SHA-256 of each input FASTQ (and the samplesheet); the demo
    has no input files. checksums.sha256 covers every output, relative to output_dir.
    """
    output_dir = Path(output_dir).resolve()
    repro_dir = output_dir / "reproducibility"
    outputs = sorted(p for p in output_dir.rglob("*") if p.is_file() and repro_dir not in p.parents)
    repro_dir.mkdir(parents=True, exist_ok=True)

    bundle = []
    if not is_demo:
        inputs: dict[str, object] = {}
        if args.samplesheet:
            sheet = Path(args.samplesheet)
            inputs["samplesheet"] = {"file": sheet.name, "sha256": sha256_file(sheet)}
        inputs["files"] = [{"file": fp.name, "sha256": sha256_file(fp)}
                           for fp in input_files if fp.is_file()]
        inputs_path = repro_dir / "inputs.json"
        write_text_lf_atomic(inputs_path, json.dumps(inputs, indent=2) + "\n")
        bundle.append(inputs_path)

    run_metadata = repro_dir / "run_metadata.json"
    write_text_lf_atomic(run_metadata, json.dumps({
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "is_demo": is_demo,
        "args": vars(args),
        "n_input_files": len(input_files),
    }, indent=2) + "\n")
    bundle.append(run_metadata)

    bundle.append(write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=Path("skills/seq-wrangler/seq_wrangler.py"),
            args=_replay_args(args, output_dir, is_demo),
            comment="Reproduce this Seq Wrangler run",
            preflight=[] if is_demo else TOOL_PREFLIGHT,
        ),
        repo_root=_PROJECT_ROOT,
    ))
    bundle.append(write_environment_yml(
        output_dir, env_name="clawbio-seq-wrangler", pip_deps=[],
        conda_deps=PIPELINE_TOOLS, channels=["conda-forge", "bioconda"],
    ))
    write_checksums([*outputs, *bundle], output_dir, anchor=output_dir)



# Pipeline runners

def resolve_samples_from_args(args: argparse.Namespace) -> List[Dict[str, str]]:
    if args.samplesheet:
        samples = load_samplesheet(Path(args.samplesheet))
    else:
        sample_name = args.sample_name or re.sub(r'[_\.]R[12][_\.].*|[_\.]R[12]$', '', Path(args.r1).stem.replace('.fastq', '').replace('.fq', ''))
        samples = [{
            "sample": sample_name,
            "fastq_1": args.r1,
            "fastq_2": args.r2 or "",
        }]
    return samples


def run_pipeline(args: argparse.Namespace) -> None:
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    samples = resolve_samples_from_args(args)
    validate_input_files(samples)
    validate_tools(args.aligner, args.trim, args.run_multiqc)

    input_fastqs: List[Path] = []
    for s in samples:
        input_fastqs.append(Path(s["fastq_1"]))
        if s["fastq_2"]:
            input_fastqs.append(Path(s["fastq_2"]))

    if args.run_fastqc:
        run_fastqc(input_fastqs, output_dir / "fastqc", args.threads)

    sample_stats: List[Dict] = []
    for s in samples:
        fq1 = Path(s["fastq_1"])
        fq2 = Path(s["fastq_2"]) if s["fastq_2"] else None
        stats = process_sample(
            sample=s["sample"],
            fq1=fq1,
            fq2=fq2,
            index=args.index,
            aligner=args.aligner,
            output_dir=output_dir,
            threads=args.threads,
            mapq=args.mapq,
            remove_dups=args.remove_duplicates,
            trim=args.trim,
            keep_sam=args.keep_sam
        )
        sample_stats.append(stats)

    if args.run_multiqc:
        run_multiqc(output_dir, output_dir / "multiqc")

    report_path = generate_report(
        output_dir=output_dir,
        sample_stats=sample_stats,
        input_files=input_fastqs,
        is_demo=False,
        aligner=args.aligner,
        mapq=args.mapq,
        remove_dups=args.remove_duplicates,
    )

    (output_dir / "summary.json").write_text(json.dumps(sample_stats, indent=2), encoding="utf-8")
    write_reproducibility_bundle(output_dir, input_fastqs, args, is_demo=False)

    print("\nDone.")
    print(f"Report: {report_path}")


def run_demo(args: argparse.Namespace) -> None:
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    sample_stats = []
    for s in DEMO_SAMPLES:
        sample = s["sample"]
        is_pe = bool(s["fastq_2"])
        sample_stats.append({
            "sample": sample,
            "is_pe": is_pe,
            "bam": str(output_dir / "bam" / f"{sample}_sorted.bam"),
            "flagstat": DEMO_FLAGSTAT[sample],
            "coverage": DEMO_COVERAGE[sample],
            "insert_size": DEMO_INSERT_SIZE.get(sample),
        })

    input_files = [Path(x["fastq_1"]) for x in DEMO_SAMPLES]
    report_path = generate_report(
        output_dir=output_dir,
        sample_stats=sample_stats,
        input_files=input_files,
        is_demo=True,
        aligner=args.aligner,
        mapq=args.mapq,
        remove_dups=args.remove_duplicates,
    )

    (output_dir / "summary.json").write_text(json.dumps(sample_stats, indent=2), encoding="utf-8")
    write_reproducibility_bundle(output_dir, input_files, args, is_demo=True)

    print("\nDemo complete.")
    print(f"Report: {report_path}")



# CLI

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Seq Wrangler: FASTQ QC, alignment, and BAM processing for read-to-BAM workflows.",
        formatter_class=argparse.ArgumentDefaultsHelpFormatter,
    )

    io_group = parser.add_argument_group("Input mode")
    io_group.add_argument("--samplesheet", type=str, default=None, help="CSV with columns: sample,fastq_1,fastq_2")
    io_group.add_argument("--r1", type=str, default=None, help="FASTQ R1 (or SE FASTQ)")
    io_group.add_argument("--r2", type=str, default=None, help="FASTQ R2 for paired-end")
    io_group.add_argument("--sample-name", type=str, default=None, help="Sample name for single-sample mode")

    parser.add_argument("--index", type=str, default=None, help="Aligner index prefix (required unless --demo)")
    parser.add_argument("--aligner", choices=SUPPORTED_ALIGNERS, default="bwa", help="Aligner backend")
    parser.add_argument("--genome-build", choices=["GRCh38", "GRCh37"], default="GRCh38",
                        help="Human genome build label for metadata/reporting")
    parser.add_argument("--keep-sam", action="store_true", 
                        help="Keep intermediate SAM files (default: delete after BAM conversion)")

    parser.add_argument("--output", type=str, required=True, help="Output directory")
    parser.add_argument("--threads", type=int, default=detect_threads(), help="Threads")
    parser.add_argument("--mapq", type=int, default=DEFAULT_MAPQ, help="MAPQ filter threshold")

    parser.add_argument("--trim", action="store_true", help="Run fastp trimming before alignment")
    parser.add_argument("--remove-duplicates", action="store_true",
                        help="Use samtools markdup -r (remove duplicates)")
    parser.add_argument("--run-fastqc", action="store_true", help="Run FastQC if available")
    parser.add_argument("--run-multiqc", action="store_true", help="Run MultiQC aggregation if available")

    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data (no external tools)")
    return parser


def validate_args(args: argparse.Namespace) -> None:
    if args.demo:
        return

    if not args.index:
        raise ValueError("--index is required unless --demo is used.")

    using_samplesheet = args.samplesheet is not None
    using_single = args.r1 is not None or args.r2 is not None

    if using_samplesheet and using_single:
        raise ValueError("Use either --samplesheet OR --r1/--r2, not both.")

    if not using_samplesheet and not args.r1:
        raise ValueError("Provide --samplesheet OR --r1 (and optional --r2).")

    if args.r2 and not args.r1:
        raise ValueError("--r2 requires --r1.")


def main() -> None:
    parser = build_parser()
    args = parser.parse_args()

    try:
        validate_args(args)
        print(f"Seq Wrangler starting | Aligner={args.aligner} | Build={args.genome_build} | Threads={args.threads}")
        if args.demo:
            run_demo(args)
        else:
            run_pipeline(args)
    except Exception as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        sys.exit(1)


if __name__ == "__main__":
    main()