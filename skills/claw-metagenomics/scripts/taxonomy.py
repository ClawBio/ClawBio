#!/usr/bin/env python3
"""Stage 1: Kraken2 + Bracken taxonomic classification.

Demo mode (r1/r2 == the bundled data/demo_R{1,2}.fastq.gz fixtures) skips
Kraken2/Bracken entirely and writes the pre-computed Peru sewage taxonomy
table (6 samples, 3 sites), matching ClawBio 1.x's `--demo` behaviour.

Real mode requires kraken2 + bracken on PATH and a Kraken2 database
(Standard-8 / PlusPF, 8-70 GB — not bundled; see SKILL.md).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import common


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--r1", required=True, type=Path, help="Forward reads FASTQ")
    p.add_argument("--r2", type=Path, default=None, help="Reverse reads FASTQ")
    p.add_argument("--kraken2-db", default=None,
                    help="Path to Kraken2 database, or 'demo' for demo runs")
    p.add_argument("--read-length", type=int, default=common.DEFAULT_READ_LENGTH)
    p.add_argument("--out", required=True, type=Path, help="Output taxonomy TSV path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    work_dir = args.out.parent

    print("[taxonomy] Kraken2 + Bracken species-level classification")

    if common.is_demo_input(args.r1, args.r2):
        print("  [demo mode] using pre-computed Peru sewage taxonomy data")
        tax_df = pd.DataFrame(common.DEMO_TAXONOMY)
        tax_df.to_csv(args.out, sep="\t", index=False)
        print(f"  Species detected: {len(tax_df)}")
        print(f"  Saved: {args.out}")
        return

    for fp in (args.r1, args.r2):
        if fp is not None and not fp.exists():
            print(f"ERROR: input file not found: {fp}", file=sys.stderr)
            sys.exit(1)

    missing = [t for t in ("kraken2", "bracken") if not common.check_tool(t)]
    if missing:
        print(f"ERROR: required tools not found on PATH: {', '.join(missing)}",
              file=sys.stderr)
        print("Install with: conda install -c bioconda " + " ".join(missing),
              file=sys.stderr)
        sys.exit(1)

    kraken2_db = None
    if args.kraken2_db and args.kraken2_db != common.SENTINEL_DEMO:
        kraken2_db = Path(args.kraken2_db)
    if kraken2_db is None or not kraken2_db.exists():
        print("ERROR: Kraken2 database not found. Pass a real --kraken2-db "
              "path (or 'demo' with the bundled demo FASTQ fixtures).",
              file=sys.stderr)
        sys.exit(1)

    threads = common.detect_threads()
    print(f"  Threads: {threads}")
    print(f"  Kraken2 DB: {kraken2_db}")

    kraken_report = common.run_kraken2(args.r1, args.r2, work_dir, kraken2_db, threads)
    bracken_output = common.run_bracken(
        kraken_report, work_dir, kraken2_db, read_length=args.read_length,
    )
    tax_df = common.parse_bracken_output(bracken_output)
    tax_df.to_csv(args.out, sep="\t", index=False)

    print(f"  Species detected: {len(tax_df)}")
    if not tax_df.empty:
        top = tax_df.iloc[0]
        name_col = "name" if "name" in tax_df.columns else tax_df.columns[0]
        print(f"  Top species: {top[name_col]} "
              f"({top.get('relative_abundance', 0):.1f}%)")
    print(f"  Saved: {args.out}")


if __name__ == "__main__":
    main()
