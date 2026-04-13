#!/usr/bin/env python3
"""Real-data pipeline for archaic introgression detection.

Supports EIGENSTRAT format (.ind/.snp/.geno) for modern genomes and
VCF for archaic reference panels. Handles both text and binary packed
EIGENSTRAT .geno files.
"""

from __future__ import annotations

import os
import shutil
import struct
import subprocess
import sys
from math import ceil
from pathlib import Path
from typing import Dict, List, Optional, Tuple

import numpy as np

# Sibling import: add parent dir to sys.path so archaic_introgression is importable
_parent = str(Path(__file__).resolve().parent)
if _parent not in sys.path:
    sys.path.insert(0, _parent)

from archaic_introgression import IntrogressionSegment

# ---------------------------------------------------------------------------
# EIGENSTRAT readers
# ---------------------------------------------------------------------------


def read_eigenstrat_ind(ind_path: str | Path) -> List[Dict[str, str]]:
    """Read an EIGENSTRAT .ind file.

    Each line: sample_name sex population
    Returns list of dicts with keys: name, sex, population.
    """
    ind_path = Path(ind_path)
    individuals: List[Dict[str, str]] = []
    with open(ind_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 3:
                individuals.append(
                    {"name": parts[0], "sex": parts[1], "population": parts[2]}
                )
            elif len(parts) == 2:
                individuals.append(
                    {"name": parts[0], "sex": parts[1], "population": "Unknown"}
                )
    return individuals


def read_eigenstrat_snp(
    snp_path: str | Path,
    chrom: Optional[str] = None,
) -> List[Dict[str, str]]:
    """Read an EIGENSTRAT .snp file.

    Each line: snp_id chrom genetic_pos physical_pos ref alt
    If chrom is provided, filter to that chromosome only.
    Returns list of dicts with keys: id, chrom, genetic_pos, physical_pos, ref, alt.
    """
    snp_path = Path(snp_path)
    snps: List[Dict[str, str]] = []
    with open(snp_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 6:
                continue
            record = {
                "id": parts[0],
                "chrom": parts[1],
                "genetic_pos": parts[2],
                "physical_pos": parts[3],
                "ref": parts[4],
                "alt": parts[5],
            }
            if chrom is None or record["chrom"] == chrom:
                snps.append(record)
    return snps


def list_chromosomes(snp_path: str | Path) -> List[str]:
    """List unique chromosomes in a .snp file, in order of first appearance."""
    snp_path = Path(snp_path)
    seen: set = set()
    chroms: List[str] = []
    with open(snp_path, "r") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) >= 2:
                c = parts[1]
                if c not in seen:
                    seen.add(c)
                    chroms.append(c)
    return chroms


def read_eigenstrat_geno(
    geno_path: str | Path,
    n_samples: int,
    n_snps: int,
) -> np.ndarray:
    """Read an EIGENSTRAT .geno file (text or binary packed format).

    Text format: each line is a string of digits (0, 1, 2, 9) with length n_samples.
    Binary packed format: starts with 'GENO' magic header. Each row is
    ceil(n_samples / 4) bytes. Each genotype is 2 bits:
      0 = hom_ref, 1 = het, 2 = hom_alt, 3 = missing (mapped to -1).

    Returns array of shape (n_samples, n_snps) with values 0, 1, 2, -1.
    """
    geno_path = Path(geno_path)

    # Check for binary format by reading first 4 bytes
    with open(geno_path, "rb") as fh:
        magic = fh.read(4)

    if magic == b"GENO":
        return _read_binary_geno(geno_path, n_samples, n_snps)
    else:
        return _read_text_geno(geno_path, n_samples, n_snps)


def _read_text_geno(
    geno_path: Path, n_samples: int, n_snps: int
) -> np.ndarray:
    """Read text-format EIGENSTRAT .geno file."""
    result = np.full((n_samples, n_snps), -1, dtype=np.int8)
    with open(geno_path, "r") as fh:
        snp_idx = 0
        for line in fh:
            line = line.strip()
            if not line:
                continue
            if snp_idx >= n_snps:
                break
            for sample_idx in range(min(len(line), n_samples)):
                ch = line[sample_idx]
                if ch == "9":
                    result[sample_idx, snp_idx] = -1
                elif ch in ("0", "1", "2"):
                    result[sample_idx, snp_idx] = int(ch)
                else:
                    result[sample_idx, snp_idx] = -1
            snp_idx += 1
    return result


def _read_binary_geno(
    geno_path: Path, n_samples: int, n_snps: int
) -> np.ndarray:
    """Read binary packed EIGENSTRAT .geno file.

    Format after 'GENO' header:
    - Each row (SNP) is ceil(n_samples / 4) bytes
    - Each byte encodes 4 genotypes, 2 bits each (LSB first)
    - 0=hom_ref, 1=het, 2=hom_alt, 3=missing
    """
    bytes_per_row = ceil(n_samples / 4)
    result = np.full((n_samples, n_snps), -1, dtype=np.int8)

    with open(geno_path, "rb") as fh:
        # Skip the 4-byte magic header
        fh.seek(4)

        for snp_idx in range(n_snps):
            row_data = fh.read(bytes_per_row)
            if len(row_data) < bytes_per_row:
                break

            sample_idx = 0
            for byte_val in row_data:
                for bit_pos in range(4):
                    if sample_idx >= n_samples:
                        break
                    geno = (byte_val >> (bit_pos * 2)) & 0x03
                    if geno == 3:
                        result[sample_idx, snp_idx] = -1
                    else:
                        result[sample_idx, snp_idx] = geno
                    sample_idx += 1
    return result


# ---------------------------------------------------------------------------
# Archaic genotype extraction
# ---------------------------------------------------------------------------


def extract_archaic_genotypes(
    archaic_vcf: str | Path,
    chrom: Optional[str] = None,
    positions: Optional[List[int]] = None,
) -> Dict[int, int]:
    """Extract archaic genotypes from a VCF file.

    Uses bcftools for .vcf.gz with .tbi index, falls back to pure Python.
    Returns dict mapping position -> genotype (0, 1, 2, or -1 for missing).
    """
    archaic_vcf = Path(archaic_vcf)
    vcf_str = str(archaic_vcf)

    # Try bcftools for indexed VCFs
    if vcf_str.endswith(".vcf.gz"):
        tbi_path = Path(vcf_str + ".tbi")
        if tbi_path.exists() and shutil.which("bcftools"):
            return _extract_with_bcftools(vcf_str, chrom, positions)

    # Pure Python fallback
    return _extract_pure_python(archaic_vcf, chrom, positions)


def _extract_with_bcftools(
    vcf_path: str,
    chrom: Optional[str],
    positions: Optional[List[int]],
) -> Dict[int, int]:
    """Extract genotypes using bcftools query."""
    cmd = ["bcftools", "query", "-f", "%POS\t[%GT]\n"]
    if chrom:
        cmd.extend(["-r", chrom])
    cmd.append(vcf_path)

    try:
        result = subprocess.run(cmd, capture_output=True, text=True, check=True)
    except (subprocess.CalledProcessError, FileNotFoundError):
        return _extract_pure_python(Path(vcf_path), chrom, positions)

    genotypes: Dict[int, int] = {}
    positions_set = set(positions) if positions else None

    for line in result.stdout.strip().split("\n"):
        if not line:
            continue
        parts = line.split("\t")
        if len(parts) < 2:
            continue
        pos = int(parts[0])
        if positions_set and pos not in positions_set:
            continue
        gt_str = parts[1]
        genotypes[pos] = _parse_gt_string(gt_str)

    return genotypes


def _extract_pure_python(
    vcf_path: Path,
    chrom: Optional[str],
    positions: Optional[List[int]],
) -> Dict[int, int]:
    """Extract genotypes using pure Python VCF parsing."""
    genotypes: Dict[int, int] = {}
    positions_set = set(positions) if positions else None

    with open(vcf_path, "r") as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            fields = line.strip().split("\t")
            if len(fields) < 10:
                continue

            record_chrom = fields[0]
            record_pos = int(fields[1])

            if chrom and record_chrom != chrom:
                continue
            if positions_set and record_pos not in positions_set:
                continue

            fmt = fields[8].split(":")
            gt_idx = fmt.index("GT") if "GT" in fmt else 0
            gt_str = fields[9].split(":")[gt_idx]
            genotypes[record_pos] = _parse_gt_string(gt_str)

    return genotypes


def _parse_gt_string(gt_str: str) -> int:
    """Parse a VCF GT string to an integer genotype."""
    alleles = gt_str.replace("|", "/").split("/")
    if "." in alleles:
        return -1
    try:
        return min(sum(int(a) for a in alleles), 2)
    except ValueError:
        return -1


# ---------------------------------------------------------------------------
# Genotype table building
# ---------------------------------------------------------------------------


def build_genotype_table(
    ind_path: str | Path,
    snp_path: str | Path,
    geno_path: str | Path,
    archaic_vcf: str | Path,
    chrom: Optional[str] = None,
) -> Tuple[np.ndarray, List[int], List[str]]:
    """Build an IBDmix-format genotype table from EIGENSTRAT + archaic VCF.

    Returns:
        geno_table: array of shape (n_samples + 1, n_shared_snps)
                    Row 0 is archaic, rows 1..n are modern samples.
        positions: list of shared positions
        sample_names: list of modern sample names
    """
    individuals = read_eigenstrat_ind(ind_path)
    snps = read_eigenstrat_snp(snp_path, chrom=chrom)
    n_samples = len(individuals)

    # Read all SNPs to get total count for geno reading
    all_snps = read_eigenstrat_snp(snp_path)
    n_total_snps = len(all_snps)

    modern_genos = read_eigenstrat_geno(geno_path, n_samples, n_total_snps)

    # Get positions for this chromosome
    if chrom:
        chrom_indices = [i for i, s in enumerate(all_snps) if s["chrom"] == chrom]
        positions = [int(all_snps[i]["physical_pos"]) for i in chrom_indices]
        modern_chrom = modern_genos[:, chrom_indices]
    else:
        positions = [int(s["physical_pos"]) for s in all_snps]
        modern_chrom = modern_genos

    # Extract archaic genotypes
    archaic_genos = extract_archaic_genotypes(archaic_vcf, chrom, positions)

    # Find shared positions
    shared_positions = sorted(set(positions) & set(archaic_genos.keys()))
    if not shared_positions:
        empty = np.empty((n_samples + 1, 0), dtype=np.int8)
        return empty, [], [ind["name"] for ind in individuals]

    # Build position index maps
    pos_to_modern_idx = {p: i for i, p in enumerate(positions)}

    # Build result table
    n_shared = len(shared_positions)
    table = np.full((n_samples + 1, n_shared), -1, dtype=np.int8)

    for j, pos in enumerate(shared_positions):
        # Archaic (row 0)
        table[0, j] = archaic_genos.get(pos, -1)
        # Modern samples (rows 1..n)
        m_idx = pos_to_modern_idx.get(pos)
        if m_idx is not None:
            table[1:, j] = modern_chrom[:, m_idx]

    sample_names = [ind["name"] for ind in individuals]
    return table, shared_positions, sample_names


# ---------------------------------------------------------------------------
# IBDmix on genotype table
# ---------------------------------------------------------------------------


def _write_ibdmix_gt_file(
    geno_table: np.ndarray,
    positions: List[int],
    sample_names: List[str],
    archaic_name: str,
    chrom: str,
    output_path: Path,
) -> Path:
    """Write genotype table in IBDmix's expected format.

    Format: chrom pos ref alt ARCHAIC SAMPLE1 SAMPLE2 ...
    Values: 0=hom_ref, 1=het, 2=hom_alt, 9=missing
    """
    output_path = Path(output_path)
    output_path.parent.mkdir(parents=True, exist_ok=True)

    with open(output_path, "w") as f:
        header = ["chrom", "pos", "ref", "alt", archaic_name] + sample_names
        f.write("\t".join(header) + "\n")

        for j, pos in enumerate(positions):
            row = [chrom, str(pos), "N", "N"]  # ref/alt unknown from EIGENSTRAT
            # Archaic genotype (row 0)
            gt = geno_table[0, j]
            row.append(str(gt if gt >= 0 else 9))
            # Modern genotypes (rows 1..n)
            for i in range(1, geno_table.shape[0]):
                gt = geno_table[i, j]
                row.append(str(gt if gt >= 0 else 9))
            f.write("\t".join(row) + "\n")

    return output_path


def _parse_ibdmix_output(
    output_path: Path, chrom: str
) -> List[IntrogressionSegment]:
    """Parse IBDmix binary output.

    Format: ID chrom start end slod
    """
    segments = []
    with open(output_path) as f:
        for line in f:
            if line.startswith("ID") or line.startswith("#"):
                continue
            parts = line.strip().split("\t")
            if len(parts) < 5:
                continue
            segments.append(
                IntrogressionSegment(
                    sample=parts[0],
                    chrom=parts[1],
                    start=int(parts[2]),
                    end=int(parts[3]),
                    archaic_source="Neandertal",
                    method="IBDmix",
                    score=float(parts[4]),
                    num_variants=0,
                )
            )
    return segments


def run_ibdmix_on_table(
    geno_table: np.ndarray,
    positions: List[int],
    sample_names: List[str],
    chrom: str = "unknown",
    lod_threshold: float = 3.0,
    archaic_name: str = "AltaiNeandertal",
    work_dir: Optional[Path] = None,
) -> List[IntrogressionSegment]:
    """Run IBDmix (binary or pure-Python fallback) on a pre-built genotype table.

    geno_table: shape (n_samples + 1, n_snps), row 0 = archaic
    """
    import tempfile

    ibdmix_bin = shutil.which("ibdmix") or shutil.which("IBDmix")

    if ibdmix_bin is not None and len(positions) > 0:
        # Write genotype table in IBDmix format
        if work_dir is None:
            work_dir = Path(tempfile.mkdtemp(prefix="ibdmix_"))
        else:
            work_dir = Path(work_dir)
            work_dir.mkdir(parents=True, exist_ok=True)

        gt_path = _write_ibdmix_gt_file(
            geno_table, positions, sample_names,
            archaic_name, chrom, work_dir / f"gt_chr{chrom}.txt",
        )
        out_path = work_dir / f"ibdmix_chr{chrom}.txt"

        cmd = [
            ibdmix_bin,
            "--genotype", str(gt_path),
            "--output", str(out_path),
            "--LOD-threshold", str(lod_threshold),
            "--archaic", archaic_name,
        ]

        try:
            result = subprocess.run(
                cmd, capture_output=True, text=True, timeout=3600,
            )
            if result.returncode == 0 and out_path.exists():
                segments = _parse_ibdmix_output(out_path, chrom)
                if segments:
                    return segments
        except (subprocess.TimeoutExpired, subprocess.CalledProcessError):
            pass

    # Pure-Python fallback
    archaic_row = geno_table[0:1, :]
    modern_rows = geno_table[1:, :]

    segments: List[IntrogressionSegment] = []
    for sample_idx in range(len(sample_names)):
        raw_segs = _lod_score_table_segment(
            modern_rows, archaic_row, positions, sample_idx, lod_threshold
        )
        for start, end, score, n_var in raw_segs:
            segments.append(
                IntrogressionSegment(
                    sample=sample_names[sample_idx],
                    chrom=chrom,
                    start=start,
                    end=end,
                    archaic_source="Neandertal",
                    method="ibdmix_fallback",
                    score=round(score, 4),
                    num_variants=n_var,
                )
            )

    return segments


def _lod_score_table_segment(
    modern_genos: np.ndarray,
    archaic_genos: np.ndarray,
    positions: List[int],
    sample_idx: int,
    lod_threshold: float = 3.0,
) -> List[Tuple[int, int, float, int]]:
    """LOD segment caller for genotype table format."""
    n_variants = modern_genos.shape[1]
    segments: List[Tuple[int, int, float, int]] = []

    current_lod = 0.0
    seg_start = None
    seg_variants = 0
    max_lod = 0.0

    for j in range(n_variants):
        m_gt = modern_genos[sample_idx, j]
        a_gt = archaic_genos[0, j]

        if m_gt < 0 or a_gt < 0:
            continue

        if m_gt > 0 and a_gt > 0:
            delta = 0.8 + 0.2 * min(m_gt, a_gt)
        elif m_gt == 0 and a_gt == 0:
            delta = 0.1
        else:
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
                if max_lod >= lod_threshold:
                    segments.append(
                        (seg_start, positions[j], max_lod, seg_variants)
                    )
                seg_start = None
                seg_variants = 0
                current_lod = 0.0
                max_lod = 0.0

    if seg_start is not None and max_lod >= lod_threshold:
        segments.append(
            (seg_start, positions[n_variants - 1], max_lod, seg_variants)
        )

    return segments


# ---------------------------------------------------------------------------
# Full chromosome pipeline
# ---------------------------------------------------------------------------


def run_chromosome(
    ind_path: str | Path,
    snp_path: str | Path,
    geno_path: str | Path,
    archaic_vcf: str | Path,
    chrom: str,
    lod_threshold: float = 3.0,
    archaic_name: str = "AltaiNeandertal",
    work_dir: Optional[Path] = None,
) -> List[IntrogressionSegment]:
    """Run the full introgression detection pipeline for one chromosome.

    Steps:
    1. Read EIGENSTRAT files
    2. Extract archaic genotypes
    3. Build genotype table
    4. Run IBDmix (binary or fallback)
    5. Return segments
    """
    geno_table, positions, sample_names = build_genotype_table(
        ind_path, snp_path, geno_path, archaic_vcf, chrom=chrom
    )

    if len(positions) == 0:
        return []

    segments = run_ibdmix_on_table(
        geno_table, positions, sample_names,
        chrom=chrom, lod_threshold=lod_threshold,
        archaic_name=archaic_name, work_dir=work_dir,
    )

    return segments


# ---------------------------------------------------------------------------
# Parallel multi-chromosome pipeline
# ---------------------------------------------------------------------------


def _run_chrom_worker(args: tuple) -> List[dict]:
    """Worker function for parallel execution. Returns dicts (picklable)."""
    ind_path, snp_path, geno_path, archaic_vcf_pattern, chrom, lod, archaic_name, work_dir = args

    # Resolve archaic VCF path for this chromosome
    archaic_vcf = Path(str(archaic_vcf_pattern).replace("{chrom}", chrom))
    if not archaic_vcf.exists():
        return []

    chrom_work = Path(work_dir) / f"chr{chrom}"
    chrom_work.mkdir(parents=True, exist_ok=True)

    segments = run_chromosome(
        ind_path, snp_path, geno_path, archaic_vcf,
        chrom=chrom, lod_threshold=lod,
        archaic_name=archaic_name, work_dir=chrom_work,
    )
    return [s.to_dict() for s in segments]


def run_full_genome(
    ind_path: str | Path,
    snp_path: str | Path,
    geno_path: str | Path,
    archaic_vcf_pattern: str,
    output_dir: str | Path,
    archaic_name: str = "AltaiNeandertal",
    chroms: Optional[List[str]] = None,
    lod_threshold: float = 3.0,
    n_workers: int = 4,
) -> List[IntrogressionSegment]:
    """Run introgression detection across all chromosomes in parallel.

    Args:
        archaic_vcf_pattern: Path with {chrom} placeholder, e.g.
            '/data/Altai/chr{chrom}_mq25_mapab100.vcf.gz'
        chroms: List of chromosome names to process. If None, auto-detect.
        n_workers: Number of parallel processes.

    Returns:
        List of all IntrogressionSegment across all chromosomes.
    """
    from concurrent.futures import ProcessPoolExecutor, as_completed

    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)

    if chroms is None:
        chroms = list_chromosomes(snp_path)
        # Filter to autosomes (1-22) for introgression analysis
        chroms = [c for c in chroms if c.isdigit() and 1 <= int(c) <= 22]

    print(f"Running {len(chroms)} chromosomes with {n_workers} workers")
    print(f"Archaic: {archaic_name}")
    print(f"VCF pattern: {archaic_vcf_pattern}")

    tasks = [
        (str(ind_path), str(snp_path), str(geno_path),
         archaic_vcf_pattern, chrom, lod_threshold,
         archaic_name, str(output_dir))
        for chrom in chroms
    ]

    all_segment_dicts = []

    with ProcessPoolExecutor(max_workers=n_workers) as executor:
        futures = {
            executor.submit(_run_chrom_worker, task): task[4]
            for task in tasks
        }
        for future in as_completed(futures):
            chrom = futures[future]
            try:
                seg_dicts = future.result()
                all_segment_dicts.extend(seg_dicts)
                print(f"  chr{chrom}: {len(seg_dicts)} segments")
            except Exception as exc:
                print(f"  chr{chrom}: FAILED ({exc})")

    # Convert dicts back to IntrogressionSegment objects
    segments = []
    for d in all_segment_dicts:
        segments.append(IntrogressionSegment(
            sample=d["sample"],
            chrom=d["chrom"],
            start=d["start"],
            end=d["end"],
            archaic_source=d["archaic_source"],
            method=d["method"],
            score=d["score"],
            num_variants=d["num_variants"],
        ))

    # Write combined BED output
    bed_path = output_dir / "all_segments.bed"
    with open(bed_path, "w") as f:
        for s in sorted(segments, key=lambda x: (x.chrom, x.start)):
            f.write(s.to_bed() + "\n")

    # Write summary CSV
    csv_path = output_dir / "segment_summary.csv"
    with open(csv_path, "w") as f:
        f.write("sample,chrom,start,end,length,archaic_source,method,score,num_variants\n")
        for s in segments:
            d = s.to_dict()
            f.write(f"{d['sample']},{d['chrom']},{d['start']},{d['end']},"
                    f"{d['length']},{d['archaic_source']},{d['method']},"
                    f"{d['score']},{d['num_variants']}\n")

    print(f"\nTotal: {len(segments)} segments across {len(chroms)} chromosomes")
    print(f"BED: {bed_path}")
    print(f"CSV: {csv_path}")

    return segments
