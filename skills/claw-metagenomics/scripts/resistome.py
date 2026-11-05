#!/usr/bin/env python3
"""Stage 2: RGI/CARD resistome profiling + WHO priority classification.

Demo mode (r1/r2 == the bundled data/demo_R{1,2}.fastq.gz fixtures) writes
the pre-computed Peru sewage resistome table (24 ARG hits across 3 sites)
instead of running RGI, matching ClawBio 1.x's `--demo` behaviour.

Real mode requires `rgi` on PATH (CARD database auto-loaded via
`rgi auto_load`, ~500 MB — not bundled; see SKILL.md). If RGI fails, both
output tables are written as `_status=failed` placeholders rather than
crashing the step, mirroring 1.x's graceful degrade.
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
    p.add_argument("--out-resistome", required=True, type=Path,
                    help="Output resistome (raw ARG hits) TSV path")
    p.add_argument("--out-who", required=True, type=Path,
                    help="Output WHO-priority-classified ARG TSV path")
    args = p.parse_args()

    args.out_resistome.parent.mkdir(parents=True, exist_ok=True)
    args.out_who.parent.mkdir(parents=True, exist_ok=True)
    work_dir = args.out_resistome.parent

    print("[resistome] RGI MAIN against CARD (Perfect + Strict hits only)")

    if common.is_demo_input(args.r1, args.r2):
        print("  [demo mode] using pre-computed Peru sewage resistome data")
        res_df = pd.DataFrame(common.DEMO_RESISTOME)
        res_df.to_csv(args.out_resistome, sep="\t", index=False)
        who_df = common.demo_who_classification()
        who_df.to_csv(args.out_who, sep="\t", index=False)

        total_perfect = sum(1 for c in common.DEMO_RESISTOME["criteria"] if c == "Perfect")
        total_strict = sum(1 for c in common.DEMO_RESISTOME["criteria"] if c == "Strict")
        n_critical = len(who_df[who_df["who_priority"] == "Critical"])
        print(f"  Total ARG hits: {len(res_df)} (Perfect: {total_perfect}, Strict: {total_strict})")
        print(f"  WHO-Critical ARGs detected: {n_critical}")
        print(f"  Saved: {args.out_resistome}")
        print(f"  Saved: {args.out_who}")
        return

    for fp in (args.r1, args.r2):
        if fp is not None and not fp.exists():
            print(f"ERROR: input file not found: {fp}", file=sys.stderr)
            sys.exit(1)

    if not common.check_tool("rgi"):
        print("ERROR: required tool not found on PATH: rgi", file=sys.stderr)
        print("Install with: conda install -c bioconda rgi", file=sys.stderr)
        sys.exit(1)

    threads = common.detect_threads()
    print(f"  Threads: {threads}")

    try:
        rgi_output = common.run_rgi(args.r1, args.r2, work_dir, threads)
        res_df = common.parse_rgi_output(rgi_output)
    except (RuntimeError, FileNotFoundError) as exc:
        print(f"  ERROR: RGI failed: {exc}", file=sys.stderr)
        common.write_status_table(args.out_resistome, "failed")
        common.write_status_table(args.out_who, "failed")
        print("  RESISTOME ANALYSIS FAILED -- see errors above.")
        return

    who_df = common.classify_who_priority(res_df)
    res_df.to_csv(args.out_resistome, sep="\t", index=False)
    who_df.to_csv(args.out_who, sep="\t", index=False)

    n_critical = len(who_df[who_df["who_priority"] == "Critical"]) if not who_df.empty else 0
    print(f"  ARG hits: {len(res_df)} (WHO-Critical: {n_critical})")
    print(f"  Saved: {args.out_resistome}")
    print(f"  Saved: {args.out_who}")


if __name__ == "__main__":
    main()
