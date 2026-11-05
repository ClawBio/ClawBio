#!/usr/bin/env python3
"""Stage 3: HUMAnN3 functional pathway profiling (MetaCyc/UniRef).

Demo mode (r1/r2 == the bundled data/demo_R{1,2}.fastq.gz fixtures) writes
the pre-computed Peru sewage pathway table instead of running HUMAnN3.

Real mode requires `humann` on PATH and ChocoPhlAn + UniRef90 databases
(~15 GB — not bundled; see SKILL.md). Pass --humann-db skip to skip
functional profiling deliberately (mirrors 1.x's --skip-functional / "no
database configured" graceful degrade) — the step still emits a table, just
one carrying `_status=skipped` instead of pathway rows, so downstream
figures/report steps can render an honest "not run" section instead of
crashing on a missing file.
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
    p.add_argument("--humann-db", default=None,
                    help="Path to HUMAnN3 databases, 'demo', or 'skip'")
    p.add_argument("--out", required=True, type=Path, help="Output pathway TSV path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)
    work_dir = args.out.parent

    print("[functional] HUMAnN3 pathway profiling (MetaCyc stratification)")

    if common.is_demo_input(args.r1, args.r2):
        print("  [demo mode] using pre-computed Peru sewage pathway data")
        pw_df = pd.DataFrame(common.DEMO_PATHWAYS)
        pw_df.to_csv(args.out, sep="\t", index=False)
        print(f"  Total pathways: {len(pw_df)}")
        print(f"  Top: {pw_df.iloc[0]['pathway']}")
        print(f"  Saved: {args.out}")
        return

    for fp in (args.r1, args.r2):
        if fp is not None and not fp.exists():
            print(f"ERROR: input file not found: {fp}", file=sys.stderr)
            sys.exit(1)

    if args.humann_db == common.SENTINEL_SKIP or not args.humann_db:
        print("  Functional profiling SKIPPED (--humann-db skip / not provided)")
        common.write_status_table(args.out, "skipped")
        return

    if not common.check_tool("humann"):
        print("ERROR: required tool not found on PATH: humann", file=sys.stderr)
        print("Install with: conda install -c bioconda humann", file=sys.stderr)
        sys.exit(1)

    humann_db = Path(args.humann_db)
    if not humann_db.exists():
        print(f"  WARNING: HUMAnN3 database path does not exist: {humann_db}",
              file=sys.stderr)
        print("  Functional profiling SKIPPED (database not found)")
        common.write_status_table(args.out, "skipped")
        return

    threads = common.detect_threads()
    print(f"  Threads: {threads}")

    pw_output = common.run_humann3(args.r1, args.r2, work_dir, threads, humann_db)
    if pw_output is not None and pw_output.exists():
        pw_df = pd.read_csv(pw_output, sep="\t", comment="#")
        pw_df.to_csv(args.out, sep="\t", index=False)
        print(f"  Pathways: {len(pw_df)}")
        print(f"  Saved: {args.out}")
    else:
        print("  ERROR: HUMAnN3 did not produce pathway output.", file=sys.stderr)
        common.write_status_table(args.out, "failed")


if __name__ == "__main__":
    main()
