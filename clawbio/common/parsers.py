"""Unified genetic file parsing for ClawBio.

Consolidates parsers from:
- skills/nutrigx_advisor/parse_input.py (cleanest modular design)
- skills/genome-compare/genome_compare.py (iCloud staging, gzip, position metadata)
- skills/pharmgx-reporter/pharmgx_reporter.py (gzip handling)
- skills/equity-scorer/equity_scorer.py (VCF matrix parsing)
- skills/gwas-prs/gwas_prs.py (PGS scoring file parser)

Supported formats: 23andMe, AncestryDNA, MyHeritage, VCF.
All parsers return dict[str, GenotypeRecord] for consistency.
"""

from __future__ import annotations

import csv
import gzip
import re
import subprocess
import sys
import tempfile
import time as _time
from dataclasses import dataclass, asdict
from pathlib import Path
from typing import Optional


@dataclass
class GenotypeRecord:
    """Single genotype call at a variant site."""

    chrom: str = ""
    pos: int = 0
    genotype: str = ""
    allele1: str = ""
    allele2: str = ""

    def to_dict(self) -> dict:
        return asdict(self)


# ---------------------------------------------------------------------------
# iCloud staging (from genome-compare)
# ---------------------------------------------------------------------------


def stage_from_icloud(filepath: str | Path) -> Path:
    """If filepath is on iCloud Drive, copy to /tmp to avoid Errno 11 deadlock.

    The macOS ``bird`` daemon holds file locks on iCloud Drive files during
    sync/indexing.  Uses subprocess ``cp`` (which macOS handles differently
    from Python's open/read for iCloud files) with retry.
    """
    filepath = Path(filepath)
    path_str = str(filepath)
    if "Mobile Documents" not in path_str and "com~apple~CloudDocs" not in path_str:
        return filepath  # not on iCloud, use directly

    cache_dir = Path(tempfile.gettempdir()) / "clawbio_cache"
    cache_dir.mkdir(exist_ok=True)
    # Use hash of full path as cache key to avoid collisions when
    # different files share the same name (e.g. multiple input.txt)
    import hashlib
    path_hash = hashlib.md5(str(filepath.resolve()).encode()).hexdigest()[:12]
    cached = cache_dir / f"{path_hash}_{filepath.name}"

    needs_copy = not cached.exists()
    if not needs_copy:
        try:
            needs_copy = filepath.stat().st_mtime > cached.stat().st_mtime
        except OSError:
            needs_copy = True

    if needs_copy:
        print(f"  [stage] copying {filepath.name} to {cached}", file=sys.stderr)
        for attempt in range(4):
            try:
                subprocess.run(
                    ["cp", str(filepath), str(cached)],
                    check=True,
                    capture_output=True,
                    timeout=30,
                )
                break
            except (subprocess.CalledProcessError, OSError) as e:
                if attempt < 3:
                    _time.sleep(2**attempt)
                else:
                    raise OSError(
                        f"Cannot stage {filepath.name} from iCloud after 4 attempts: {e}"
                    )
    return cached


# ---------------------------------------------------------------------------
# File open helper (gzip-transparent)
# ---------------------------------------------------------------------------


def open_genetic_file(filepath: str | Path):
    """Open a file, handling .gz transparently. Stages from iCloud first."""
    filepath = str(stage_from_icloud(Path(filepath)))
    if filepath.endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    return open(filepath, encoding="utf-8", errors="replace")


# ---------------------------------------------------------------------------
# Format detection
# ---------------------------------------------------------------------------


def detect_format(filepath: str | Path) -> str:
    """Auto-detect genetic file format from header.

    Returns one of: "23andme", "ancestry", "vcf", "unknown".
    """
    filepath = Path(filepath)

    # Check compound extensions first
    suffixes = "".join(filepath.suffixes).lower()
    if ".vcf" in suffixes:
        return "vcf"

    with open_genetic_file(filepath) as f:
        for line in f:
            if line.startswith("##fileformat=VCF"):
                return "vcf"
            lower = line.lower()
            if "rsid" in lower and "chromosome" in lower and "genotype" in lower:
                return "23andme"
            if "rsid" in lower and "allele1" in lower:
                return "ancestry"
            # MyHeritage: CSV with RSID,CHROMOSOME,POSITION,RESULT header
            if "rsid" in lower and "chromosome" in lower and "result" in lower:
                return "myheritage"
            if not line.startswith("#"):
                break

    # Fallback: infer from extension
    ext = filepath.suffix.lower()
    if ext == ".vcf":
        return "vcf"

    raise ValueError(
        f"Cannot auto-detect genetic file format for '{filepath}'. "
        f"No recognized header found and extension '{ext}' is ambiguous. "
        f"Please specify --format explicitly (23andme, ancestry, or vcf)."
    )


# ---------------------------------------------------------------------------
# 23andMe parser
# ---------------------------------------------------------------------------


def parse_23andme(filepath: str | Path) -> dict[str, GenotypeRecord]:
    """Parse 23andMe raw data file.

    Returns {rsid: GenotypeRecord} with chromosome/position metadata.
    Handles gzip transparently.
    """
    genotypes: dict[str, GenotypeRecord] = {}
    with open_genetic_file(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rsid, chrom, pos_str, geno = parts[0], parts[1], parts[2], parts[3]
            if not (rsid.startswith("rs") or rsid.startswith("i")):
                continue
            geno = geno.replace("-", "").replace("--", "")
            if len(geno) < 1:
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                pos = 0

            allele1 = geno[0] if len(geno) >= 1 else ""
            allele2 = geno[1] if len(geno) >= 2 else ""

            genotypes[rsid] = GenotypeRecord(
                chrom=chrom,
                pos=pos,
                genotype=geno,
                allele1=allele1,
                allele2=allele2,
            )
    return genotypes


# ---------------------------------------------------------------------------
# AncestryDNA parser
# ---------------------------------------------------------------------------


def parse_ancestry(filepath: str | Path) -> dict[str, GenotypeRecord]:
    """Parse AncestryDNA raw data file.

    Returns {rsid: GenotypeRecord}.
    """
    genotypes: dict[str, GenotypeRecord] = {}

    # Read lines, skipping comment lines
    lines = []
    with open_genetic_file(filepath) as f:
        for line in f:
            if not line.startswith("#"):
                lines.append(line)

    reader = csv.DictReader(lines, delimiter="\t")
    for row in reader:
        rsid = row.get("rsid", "").strip()
        allele1 = row.get("allele1", "").strip()
        allele2 = row.get("allele2", "").strip()
        chrom = row.get("chromosome", row.get("chr", "")).strip()
        pos_str = row.get("position", row.get("pos", "0")).strip()

        if not rsid.startswith("rs"):
            continue

        geno = allele1 + allele2
        try:
            pos = int(pos_str)
        except ValueError:
            pos = 0

        genotypes[rsid] = GenotypeRecord(
            chrom=chrom,
            pos=pos,
            genotype=geno,
            allele1=allele1,
            allele2=allele2,
        )
    return genotypes


# ---------------------------------------------------------------------------
# MyHeritage parser
# ---------------------------------------------------------------------------


def parse_myheritage(filepath: str | Path) -> dict[str, GenotypeRecord]:
    """Parse MyHeritage DNA raw data file.

    Format: CSV with columns RSID,CHROMOSOME,POSITION,RESULT
    Comment lines begin with '#'. RESULT is the concatenated genotype (e.g. 'AA').

    Returns {rsid: GenotypeRecord}.
    """
    genotypes: dict[str, GenotypeRecord] = {}

    lines = []
    with open_genetic_file(filepath) as f:
        for line in f:
            if not line.startswith("#"):
                lines.append(line)

    reader = csv.DictReader(lines)
    for row in reader:
        # Normalise to lowercase keys to handle any header casing
        row = {k.lower(): v for k, v in row.items()}
        rsid = (row.get("rsid") or "").strip()
        chrom = (row.get("chromosome") or "").strip()
        pos_str = (row.get("position") or "0").strip()
        result = (row.get("result") or "").strip()

        if not rsid.startswith("rs"):
            continue
        geno = result.replace("-", "")
        if not geno:
            continue

        try:
            pos = int(pos_str)
        except ValueError:
            pos = 0

        allele1 = geno[0] if len(geno) >= 1 else ""
        allele2 = geno[1] if len(geno) >= 2 else ""

        genotypes[rsid] = GenotypeRecord(
            chrom=chrom,
            pos=pos,
            genotype=geno,
            allele1=allele1,
            allele2=allele2,
        )
    return genotypes


# ---------------------------------------------------------------------------
# VCF parser (single-sample genotype dict)
# ---------------------------------------------------------------------------


def parse_vcf(filepath: str | Path) -> dict[str, GenotypeRecord]:
    """Parse VCF file, extracting GT field for the first sample.

    Returns {rsid: GenotypeRecord} with genotype as concatenated bases.
    """
    genotypes: dict[str, GenotypeRecord] = {}

    with open_genetic_file(filepath) as f:
        for line in f:
            line = line.strip()
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                continue
            parts = line.split("\t")
            if len(parts) < 10:
                continue

            chrom = parts[0]
            pos_str = parts[1]
            rsid = parts[2]
            ref = parts[3]
            alts = parts[4].split(",")
            alleles = [ref] + alts

            if not rsid.startswith("rs"):
                continue

            fmt = parts[8].split(":")
            if "GT" not in fmt:
                continue
            gt_idx = fmt.index("GT")
            sample = parts[9].split(":")[gt_idx]

            # Handle phased (|) or unphased (/)
            indices = re.split(r"[|/]", sample)
            try:
                called = [alleles[int(i)] for i in indices if i != "."]
                geno = "".join(called)
            except (IndexError, ValueError):
                continue

            try:
                pos = int(pos_str)
            except ValueError:
                pos = 0

            allele1 = called[0] if len(called) >= 1 else ""
            allele2 = called[1] if len(called) >= 2 else ""

            genotypes[rsid] = GenotypeRecord(
                chrom=chrom,
                pos=pos,
                genotype=geno,
                allele1=allele1,
                allele2=allele2,
            )
    return genotypes


# ---------------------------------------------------------------------------
# VCF matrix parser (multi-sample, for equity-scorer)
# ---------------------------------------------------------------------------


def parse_vcf_matrix(filepath: str | Path):
    """Parse a VCF file into a genotype matrix for population genetics.

    Returns:
        samples: list of sample IDs
        variant_ids: list of variant IDs (or CHROM:POS)
        genotype_matrix: numpy array of shape (n_samples, n_variants)
                         with values 0 (hom ref), 1 (het), 2 (hom alt), -1 (missing)

    Requires numpy (imported locally to keep base parser dependency-free).
    """
    import numpy as np

    samples: list[str] = []
    variant_ids: list[str] = []
    genotype_rows: list[list[int]] = []

    with open_genetic_file(filepath) as f:
        for line in f:
            if line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                parts = line.strip().split("\t")
                samples = parts[9:]
                continue

            parts = line.strip().split("\t")
            chrom, pos, vid = parts[0], parts[1], parts[2]
            if vid == ".":
                vid = f"{chrom}:{pos}"
            variant_ids.append(vid)

            fmt_fields = parts[8].split(":")
            if "GT" not in fmt_fields:
                raise ValueError(
                    f"GT field not found in FORMAT column at variant {vid} "
                    f"(FORMAT={parts[8]}). Cannot parse genotypes."
                )
            gt_idx = fmt_fields.index("GT")

            row = []
            for sample_field in parts[9:]:
                gt_str = sample_field.split(":")[gt_idx]
                gt_str = gt_str.replace("|", "/")
                if "." in gt_str:
                    row.append(-1)
                else:
                    allele_indices = gt_str.split("/")
                    row.append(int(allele_indices[0]) + int(allele_indices[1]))
            genotype_rows.append(row)

    if not samples:
        raise ValueError("No samples found in VCF header")
    if not genotype_rows:
        raise ValueError("No variants found in VCF")

    geno_matrix = np.array(genotype_rows, dtype=np.int8).T
    return samples, variant_ids, geno_matrix


# ---------------------------------------------------------------------------
# Unified parser entry point
# ---------------------------------------------------------------------------


def parse_genetic_file(
    filepath: str | Path,
    fmt: str = "auto",
) -> dict[str, GenotypeRecord]:
    """Parse genetic data file in any supported format.

    Args:
        filepath: Path to genetic data file (23andMe, AncestryDNA, MyHeritage, VCF).
        fmt: Format hint — "auto", "23andme", "ancestry", "myheritage", or "vcf".

    Returns:
        dict mapping rsid -> GenotypeRecord.
    """
    filepath = Path(filepath)

    if fmt == "auto":
        fmt = detect_format(filepath)

    parsers = {
        "23andme": parse_23andme,
        "ancestry": parse_ancestry,
        "myheritage": parse_myheritage,
        "vcf": parse_vcf,
    }

    if fmt not in parsers:
        raise ValueError(f"Unknown format: {fmt}. Choose from: {list(parsers.keys())}")

    return parsers[fmt](filepath)


# ---------------------------------------------------------------------------
# Convenience: extract simple genotype dict {rsid: genotype_str}
# ---------------------------------------------------------------------------


def genotypes_to_simple(records: dict[str, GenotypeRecord]) -> dict[str, str]:
    """Convert GenotypeRecord dict to simple {rsid: genotype_str} dict.

    This is the backward-compatible format used by existing skill code.
    """
    return {rsid: rec.genotype for rsid, rec in records.items()}


def genotypes_to_positions(
    records: dict[str, GenotypeRecord],
) -> dict[str, dict]:
    """Extract {rsid: {"chrom": str, "pos": int}} from GenotypeRecord dict."""
    return {
        rsid: {"chrom": rec.chrom, "pos": rec.pos}
        for rsid, rec in records.items()
    }
