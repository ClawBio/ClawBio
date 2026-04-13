#!/usr/bin/env python3
"""Archaic introgression detection skill for ClawBio.

Detects Neanderthal and Denisovan introgression segments from modern
human genomes using IBDmix, Sprime, hmmix, or a pure-Python LOD fallback.
"""

from __future__ import annotations

import argparse
import json
import math
import os
import shutil
import subprocess
import sys
from dataclasses import dataclass, field, asdict
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

JAVA_HOME = os.environ.get(
    "JAVA_HOME",
    "/opt/homebrew/opt/openjdk/libexec/openjdk.jdk/Contents/Home",
)
SPRIME_JAR = os.environ.get("SPRIME_JAR", "/opt/homebrew/bin/sprime.jar")

SKILL_DIR = Path(__file__).resolve().parent
DEMO_MODERN = SKILL_DIR / "examples" / "demo_modern.vcf"
DEMO_ARCHAIC = SKILL_DIR / "examples" / "demo_archaic.vcf"

# ---------------------------------------------------------------------------
# Data classes
# ---------------------------------------------------------------------------


@dataclass
class IntrogressionSegment:
    """A single detected archaic introgression segment."""

    sample: str
    chrom: str
    start: int
    end: int
    archaic_source: str = "Neanderthal"
    method: str = "ibdmix"
    score: float = 0.0
    num_variants: int = 0

    @property
    def length(self) -> int:
        """Length of the segment in base pairs."""
        return self.end - self.start

    def to_dict(self) -> dict:
        """Convert to dictionary including computed properties."""
        d = asdict(self)
        d["length"] = self.length
        return d

    def to_bed(self) -> str:
        """Return a BED-format line (0-based start)."""
        bed_start = self.start - 1 if self.start > 0 else 0
        return f"{self.chrom}\t{bed_start}\t{self.end}\t{self.sample}\t{self.score:.2f}"


# ---------------------------------------------------------------------------
# Tool discovery
# ---------------------------------------------------------------------------


def find_tool(name: str) -> Optional[str]:
    """Find an external tool binary on PATH.

    For sprime, also checks the default SPRIME_JAR location.
    Returns the path string or None.
    """
    name_lower = name.lower()

    if name_lower == "sprime":
        jar = Path(SPRIME_JAR)
        if jar.exists():
            return str(jar)
        # Also try PATH
        found = shutil.which("sprime")
        if found:
            return found
        return None

    # ibdmix, hmmix, or anything else
    found = shutil.which(name_lower)
    if found:
        return found

    # Try capitalised variants
    for variant in [name, name.upper(), name.capitalize()]:
        found = shutil.which(variant)
        if found:
            return found

    return None


def list_available_methods() -> List[str]:
    """Return names of methods whose external tools are available."""
    methods = []
    if find_tool("ibdmix") is not None:
        methods.append("ibdmix")
    if find_tool("sprime") is not None:
        methods.append("sprime")
    if find_tool("hmmix") is not None:
        methods.append("hmmix")
    # Pure-python fallback is always available
    if "ibdmix" not in methods:
        methods.append("ibdmix_fallback")
    return methods


# ---------------------------------------------------------------------------
# VCF parsing
# ---------------------------------------------------------------------------


def parse_vcf_samples(vcf_path: str | Path) -> List[str]:
    """Extract sample names from a VCF header line."""
    vcf_path = Path(vcf_path)
    with open(vcf_path, "r") as fh:
        for line in fh:
            if line.startswith("#CHROM"):
                fields = line.strip().split("\t")
                return fields[9:]
    return []


def parse_vcf_variants(
    vcf_path: str | Path,
) -> List[Dict[str, str]]:
    """Parse variant records from a VCF file.

    Returns a list of dicts with keys: chrom, pos, id, ref, alt, info.
    """
    vcf_path = Path(vcf_path)
    variants: List[Dict[str, str]] = []
    with open(vcf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 8:
                continue
            info_str = fields[7]
            variants.append(
                {
                    "chrom": fields[0],
                    "pos": fields[1],
                    "id": fields[2],
                    "ref": fields[3],
                    "alt": fields[4],
                    "info": info_str,
                }
            )
    return variants


def parse_vcf_genotypes(vcf_path: str | Path) -> np.ndarray:
    """Parse genotypes from a VCF into a numpy array.

    Returns array of shape (n_samples, n_variants) with values:
      0 = hom_ref, 1 = het, 2 = hom_alt, -1 = missing
    """
    vcf_path = Path(vcf_path)
    samples = parse_vcf_samples(vcf_path)
    n_samples = len(samples)

    rows: List[List[int]] = []
    with open(vcf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue
            # Find GT index in FORMAT
            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0

            geno_row: List[int] = []
            for i in range(9, 9 + n_samples):
                if i >= len(fields):
                    geno_row.append(-1)
                    continue
                gt_field = fields[i].split(":")[gt_idx]
                alleles = gt_field.replace("|", "/").split("/")
                if "." in alleles:
                    geno_row.append(-1)
                else:
                    alt_count = sum(int(a) for a in alleles)
                    geno_row.append(min(alt_count, 2))
            rows.append(geno_row)

    if not rows:
        return np.empty((n_samples, 0), dtype=np.int8)

    # Shape: (n_variants, n_samples) then transpose to (n_samples, n_variants)
    arr = np.array(rows, dtype=np.int8)
    return arr.T


def get_shared_positions(
    modern_vcf: str | Path, archaic_vcf: str | Path
) -> List[int]:
    """Return sorted list of positions present in both VCFs."""
    modern_variants = parse_vcf_variants(modern_vcf)
    archaic_variants = parse_vcf_variants(archaic_vcf)

    modern_pos = {v["pos"] for v in modern_variants}
    archaic_pos = {v["pos"] for v in archaic_variants}

    shared = sorted(int(p) for p in modern_pos & archaic_pos)
    return shared


# ---------------------------------------------------------------------------
# Command builders
# ---------------------------------------------------------------------------


def build_ibdmix_command(
    modern_vcf: str,
    archaic_vcf: str,
    output_dir: str,
    lod_threshold: float = 3.0,
) -> List[str]:
    """Build IBDmix command line."""
    ibdmix_bin = find_tool("ibdmix")
    if ibdmix_bin is None:
        return []
    return [
        ibdmix_bin,
        "--modern", modern_vcf,
        "--archaic", archaic_vcf,
        "--output", output_dir,
        "--lod", str(lod_threshold),
    ]


def build_sprime_command(
    modern_vcf: str,
    outgroup_vcf: str,
    output_prefix: str,
) -> List[str]:
    """Build Sprime command line."""
    sprime_path = find_tool("sprime")
    if sprime_path is None:
        return []

    java_bin = str(Path(JAVA_HOME) / "bin" / "java")
    if not Path(java_bin).exists():
        java_bin = "java"

    if sprime_path.endswith(".jar"):
        return [
            java_bin, "-jar", sprime_path,
            "gt=" + modern_vcf,
            "outgroup=" + outgroup_vcf,
            "output=" + output_prefix,
        ]
    return [
        sprime_path,
        "--vcf", modern_vcf,
        "--outgroup", outgroup_vcf,
        "--out", output_prefix,
    ]


def build_hmmix_command(
    modern_vcf: str,
    output_dir: str,
) -> List[str]:
    """Build hmmix command line."""
    hmmix_bin = find_tool("hmmix")
    if hmmix_bin is None:
        return []
    return [
        hmmix_bin,
        "--infile", modern_vcf,
        "--outfolder", output_dir,
    ]


# ---------------------------------------------------------------------------
# IBDmix runner (binary + pure-Python fallback)
# ---------------------------------------------------------------------------


def _lod_score_segment(
    modern_genos: np.ndarray,
    archaic_genos: np.ndarray,
    positions: List[int],
    sample_idx: int,
    lod_threshold: float = 3.0,
) -> List[Tuple[int, int, float, int]]:
    """Pure-Python LOD segment caller for one sample.

    Returns list of (start, end, lod_score, num_variants).
    Uses a simplified LOD model: for each shared variant, the LOD
    contribution is positive when modern and archaic share an allele
    and negative otherwise.
    """
    n_variants = modern_genos.shape[1]
    segments: List[Tuple[int, int, float, int]] = []

    current_lod = 0.0
    seg_start = None
    seg_variants = 0
    max_lod = 0.0

    for j in range(n_variants):
        m_gt = modern_genos[sample_idx, j]
        a_gt = archaic_genos[0, j] if archaic_genos.shape[0] > 0 else -1

        if m_gt < 0 or a_gt < 0:
            continue

        # Simplified LOD contribution
        if m_gt > 0 and a_gt > 0:
            # Both carry alt allele: evidence for shared ancestry
            delta = 0.8 + 0.2 * min(m_gt, a_gt)
        elif m_gt == 0 and a_gt == 0:
            # Both hom-ref: weak positive signal
            delta = 0.1
        else:
            # Discordant: evidence against
            delta = -0.5

        current_lod += delta
        if current_lod > max_lod:
            max_lod = current_lod

        if seg_start is None and current_lod > 0:
            seg_start = positions[j]
            seg_variants = 1
        elif seg_start is not None:
            seg_variants += 1
            if current_lod < 0:
                # End the segment
                if max_lod >= lod_threshold:
                    segments.append(
                        (seg_start, positions[j], max_lod, seg_variants)
                    )
                seg_start = None
                seg_variants = 0
                current_lod = 0.0
                max_lod = 0.0

    # Close any open segment
    if seg_start is not None and max_lod >= lod_threshold:
        segments.append(
            (seg_start, positions[n_variants - 1], max_lod, seg_variants)
        )

    return segments


def run_ibdmix(
    modern_vcf: str | Path,
    archaic_vcf: str | Path,
    output_dir: str | Path,
    lod_threshold: float = 3.0,
    samples: Optional[List[str]] = None,
) -> List[IntrogressionSegment]:
    """Run IBDmix (binary or pure-Python fallback).

    Returns a list of IntrogressionSegment objects.
    """
    modern_vcf = Path(modern_vcf)
    archaic_vcf = Path(archaic_vcf)
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Try external binary first
    ibdmix_bin = find_tool("ibdmix")
    if ibdmix_bin is not None:
        cmd = build_ibdmix_command(
            str(modern_vcf), str(archaic_vcf), str(output_dir), lod_threshold
        )
        try:
            subprocess.run(cmd, check=True, capture_output=True, text=True)
            # Parse IBDmix output (tab-separated: sample, chrom, start, end, lod, n_snps)
            results_file = output_dir / "ibdmix_results.txt"
            if results_file.exists():
                return _parse_ibdmix_output(results_file)
        except (subprocess.CalledProcessError, FileNotFoundError):
            pass  # Fall through to pure-Python

    # Pure-Python fallback
    return _run_ibdmix_fallback(
        modern_vcf, archaic_vcf, output_dir, lod_threshold, samples
    )


def _parse_ibdmix_output(results_file: Path) -> List[IntrogressionSegment]:
    """Parse IBDmix binary output file."""
    segments: List[IntrogressionSegment] = []
    with open(results_file) as fh:
        for line in fh:
            if line.startswith("#") or line.startswith("sample"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 6:
                continue
            segments.append(
                IntrogressionSegment(
                    sample=parts[0],
                    chrom=parts[1],
                    start=int(parts[2]),
                    end=int(parts[3]),
                    archaic_source="Neanderthal",
                    method="ibdmix",
                    score=float(parts[4]),
                    num_variants=int(parts[5]),
                )
            )
    return segments


def _run_ibdmix_fallback(
    modern_vcf: Path,
    archaic_vcf: Path,
    output_dir: Path,
    lod_threshold: float,
    samples: Optional[List[str]] = None,
) -> List[IntrogressionSegment]:
    """Pure-Python LOD-based IBDmix fallback."""
    all_samples = parse_vcf_samples(modern_vcf)
    modern_genos = parse_vcf_genotypes(modern_vcf)
    archaic_genos = parse_vcf_genotypes(archaic_vcf)
    modern_variants = parse_vcf_variants(modern_vcf)
    archaic_variants = parse_vcf_variants(archaic_vcf)

    # Build position arrays
    modern_positions = [int(v["pos"]) for v in modern_variants]
    archaic_positions_set = {int(v["pos"]) for v in archaic_variants}

    # Filter to shared positions
    shared_mask = [i for i, p in enumerate(modern_positions) if p in archaic_positions_set]
    if not shared_mask:
        return []

    positions = [modern_positions[i] for i in shared_mask]
    m_genos = modern_genos[:, shared_mask]

    # Build archaic genotype array at shared positions
    archaic_pos_list = [int(v["pos"]) for v in archaic_variants]
    archaic_pos_to_idx = {p: i for i, p in enumerate(archaic_pos_list)}
    archaic_shared_mask = [archaic_pos_to_idx[p] for p in positions]
    a_genos = archaic_genos[:, archaic_shared_mask]

    # Determine chromosome from first variant
    chrom = modern_variants[shared_mask[0]]["chrom"] if modern_variants else "unknown"

    # Filter samples
    if samples:
        sample_indices = [i for i, s in enumerate(all_samples) if s in samples]
        sample_names = [all_samples[i] for i in sample_indices]
    else:
        sample_indices = list(range(len(all_samples)))
        sample_names = all_samples

    segments: List[IntrogressionSegment] = []
    for idx, sample_idx in enumerate(sample_indices):
        raw_segs = _lod_score_segment(
            m_genos, a_genos, positions, sample_idx, lod_threshold
        )
        for start, end, score, n_var in raw_segs:
            segments.append(
                IntrogressionSegment(
                    sample=sample_names[idx],
                    chrom=chrom,
                    start=start,
                    end=end,
                    archaic_source="Neanderthal",
                    method="ibdmix_fallback",
                    score=round(score, 4),
                    num_variants=n_var,
                )
            )

    return segments


# ---------------------------------------------------------------------------
# Summary and reporting
# ---------------------------------------------------------------------------


def compute_summary(segments: List[IntrogressionSegment]) -> Dict[str, dict]:
    """Compute per-individual summary statistics from segments."""
    by_sample: Dict[str, List[IntrogressionSegment]] = {}
    for seg in segments:
        by_sample.setdefault(seg.sample, []).append(seg)

    summary: Dict[str, dict] = {}
    for sample, segs in by_sample.items():
        lengths = [s.length for s in segs]
        scores = [s.score for s in segs]
        summary[sample] = {
            "total_segments": len(segs),
            "total_length_bp": sum(lengths),
            "mean_segment_length": round(sum(lengths) / len(segs), 1) if segs else 0.0,
            "mean_score": round(sum(scores) / len(scores), 4) if scores else 0.0,
        }
    return summary


def write_report(
    segments: List[IntrogressionSegment],
    output_dir: Path,
    method: str = "ibdmix",
    lod_threshold: float = 3.0,
) -> Path:
    """Write a BED file of introgression segments."""
    output_dir.mkdir(parents=True, exist_ok=True)
    bed_path = output_dir / "segments.bed"
    with open(bed_path, "w") as fh:
        fh.write("#chrom\tstart\tend\tsample\tscore\n")
        for seg in segments:
            fh.write(seg.to_bed() + "\n")
    return bed_path


def write_result_json(
    segments: List[IntrogressionSegment],
    output_dir: Path,
    method: str = "ibdmix",
    lod_threshold: float = 3.0,
) -> Path:
    """Write full JSON results with segments and summary."""
    output_dir.mkdir(parents=True, exist_ok=True)
    samples_in_results = set(s.sample for s in segments)
    summary = compute_summary(segments)

    result = {
        "method": method,
        "lod_threshold": lod_threshold,
        "num_samples": len(samples_in_results),
        "segments": [s.to_dict() for s in segments],
        "summary": summary,
    }

    json_path = output_dir / "introgression_results.json"
    with open(json_path, "w") as fh:
        json.dump(result, fh, indent=2)
    return json_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Build CLI argument parser."""
    parser = argparse.ArgumentParser(
        description="Detect archaic introgression segments in modern human genomes",
    )
    parser.add_argument(
        "--input", "-i",
        type=str,
        help="Path to modern human VCF file",
    )
    parser.add_argument(
        "--archaic", "-a",
        type=str,
        help="Path to archaic reference VCF file",
    )
    parser.add_argument(
        "--output", "-o",
        type=str,
        default="/tmp/introgression_output",
        help="Output directory (default: /tmp/introgression_output)",
    )
    parser.add_argument(
        "--method", "-m",
        type=str,
        default="ibdmix",
        choices=["ibdmix", "sprime", "hmmix"],
        help="Detection method (default: ibdmix)",
    )
    parser.add_argument(
        "--samples", "-s",
        type=str,
        default=None,
        help="Comma-separated list of sample names to analyse",
    )
    parser.add_argument(
        "--lod",
        type=float,
        default=3.0,
        help="LOD score threshold for segment calling (default: 3.0)",
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run with bundled demo data",
    )
    return parser


def main(args: Optional[List[str]] = None) -> int:
    """Main entry point."""
    parser = build_parser()
    opts = parser.parse_args(args)

    if opts.demo:
        opts.input = str(DEMO_MODERN)
        opts.archaic = str(DEMO_ARCHAIC)

    if not opts.input or not opts.archaic:
        parser.error("--input and --archaic are required (or use --demo)")

    if not Path(opts.input).exists():
        print(f"Error: modern VCF not found: {opts.input}", file=sys.stderr)
        return 1
    if not Path(opts.archaic).exists():
        print(f"Error: archaic VCF not found: {opts.archaic}", file=sys.stderr)
        return 1

    output_dir = Path(opts.output)
    sample_list = opts.samples.split(",") if opts.samples else None

    # Run detection
    segments = run_ibdmix(
        opts.input, opts.archaic, output_dir, opts.lod, sample_list
    )

    # Write outputs
    json_path = write_result_json(segments, output_dir, opts.method, opts.lod)
    bed_path = write_report(segments, output_dir, opts.method, opts.lod)

    # Print summary
    summary = compute_summary(segments)
    print(f"Found {len(segments)} introgression segment(s)")
    for sample, stats in summary.items():
        print(
            f"  {sample}: {stats['total_segments']} segments, "
            f"{stats['total_length_bp']} bp total, "
            f"mean LOD {stats['mean_score']}"
        )
    print(f"\nResults written to: {json_path}")
    print(f"BED file written to: {bed_path}")

    return 0


if __name__ == "__main__":
    sys.exit(main())
