#!/usr/bin/env python3
"""Stage 2 of the gwas_prs DAG: parse a DTC genotype file.

Port of ClawBio 1.x `gwas-prs`. The 1.x script delegated to the shared
`clawbio.common.parsers` package (23andMe/AncestryDNA/MyHeritage/VCF
parsers consolidated from several 1.x skills); that package does not exist
in this repo, so the two formats gwas-prs's own SKILL.md documents
(23andMe, AncestryDNA) are inlined here directly, dependency-free.

Independent of scores-file resolution (scripts/resolve_scores.py) — this
stage only reads the genotype file — so the two run as DAG siblings.
"""

from __future__ import annotations

import argparse
import csv
import gzip
import json
import sys
from pathlib import Path


def _open_genetic_file(filepath: Path):
    if str(filepath).endswith(".gz"):
        return gzip.open(filepath, "rt", encoding="utf-8", errors="replace")
    return open(filepath, encoding="utf-8", errors="replace")


def detect_format(filepath: Path) -> str:
    """Auto-detect 23andMe vs AncestryDNA from header lines."""
    with _open_genetic_file(filepath) as f:
        for line in f:
            lower = line.lower()
            if "rsid" in lower and "chromosome" in lower and "genotype" in lower:
                return "23andme"
            if "rsid" in lower and "allele1" in lower:
                return "ancestry"
            if not line.startswith("#"):
                break
    raise ValueError(
        f"Cannot auto-detect genotype file format for '{filepath}'. "
        "Expected a 23andMe (rsid/chromosome/.../genotype) or AncestryDNA "
        "(rsid/chromosome/.../allele1/allele2) header."
    )


def parse_23andme(filepath: Path) -> dict[str, str]:
    """Parse a 23andMe raw data file. Returns {rsid: genotype_string}."""
    genotypes: dict[str, str] = {}
    with _open_genetic_file(filepath) as f:
        for line in f:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 4:
                continue
            rsid, geno = parts[0], parts[3]
            if not (rsid.startswith("rs") or rsid.startswith("i")):
                continue
            geno = geno.replace("-", "")
            if not geno:
                continue
            genotypes[rsid] = geno
    return genotypes


def parse_ancestry(filepath: Path) -> dict[str, str]:
    """Parse an AncestryDNA raw data file. Returns {rsid: genotype_string}."""
    genotypes: dict[str, str] = {}
    lines = []
    with _open_genetic_file(filepath) as f:
        for line in f:
            if not line.startswith("#"):
                lines.append(line)

    delimiter = "\t" if "\t" in (lines[0] if lines else "") else ","
    reader = csv.DictReader(lines, delimiter=delimiter)
    for row in reader:
        row = {k.strip().lower(): (v or "").strip() for k, v in row.items()}
        rsid = row.get("rsid", "")
        allele1 = row.get("allele1", "")
        allele2 = row.get("allele2", "")
        if not rsid.startswith("rs"):
            continue
        geno = allele1 + allele2
        if not geno:
            continue
        genotypes[rsid] = geno
    return genotypes


def load_genotypes(filepath: Path) -> tuple[str, dict[str, str]]:
    fmt = detect_format(filepath)
    if fmt == "23andme":
        return fmt, parse_23andme(filepath)
    if fmt == "ancestry":
        return fmt, parse_ancestry(filepath)
    raise ValueError(f"Unsupported format: {fmt}")


def main() -> None:
    parser = argparse.ArgumentParser(description="Parse a 23andMe/AncestryDNA genotype file")
    parser.add_argument("--input", required=True, help="Path to genotype file (.txt or .txt.gz)")
    parser.add_argument("--out", required=True, help="Output path for genotype.json")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}")
        sys.exit(1)

    print(f"Parsing genotype file: {input_path.name}")
    fmt, genotypes = load_genotypes(input_path)
    print(f"  Format: {fmt}")
    print(f"  Total SNPs parsed: {len(genotypes):,}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        "filepath": str(input_path),
        "format": fmt,
        "total_snps": len(genotypes),
        "genotypes": genotypes,
    }, indent=2))
    print(f"Wrote {out_path}")


if __name__ == "__main__":
    main()
