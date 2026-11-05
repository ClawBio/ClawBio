#!/usr/bin/env python3
"""Stage 5: post-inference rooting.

Outgroup rooting already happened inside build_tree.py (a tree-builder CLI
flag). This stage only applies midpoint rooting (`--root midpoint`), via
ETE3 with a Bio.Phylo fallback, then an identity fallback if both are
unavailable/fail — exactly 1.x's `root_midpoint()`, which never raises.
`--root none` (default) passes the tree through untouched.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import common


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path, help="Newick tree from build_tree")
    p.add_argument("--root", default=common.SENTINEL_NONE, help="'midpoint' or 'none'")
    p.add_argument("--outgroup", default=common.SENTINEL_NONE, help="Outgroup taxa (recorded only; already applied at build_tree)")
    p.add_argument("--out", required=True, type=Path, help="Output (rooted) Newick tree path")
    p.add_argument("--out-meta", required=True, type=Path, help="Output stage-metadata JSON path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    newick = args.input.read_text().strip()
    root_mode = args.root.strip().lower() if args.root else common.SENTINEL_NONE
    outgroup = common.none_if_sentinel(args.outgroup)

    print(f"[root] input={args.input} root={root_mode} outgroup={outgroup}")

    if root_mode == "midpoint" and newick:
        newick = common.root_midpoint(newick)
        stage = "root:midpoint"
    elif outgroup:
        stage = f"root:outgroup({outgroup})"
    else:
        stage = "root:none"

    args.out.write_text(newick + "\n")
    meta = {"stage": stage, "root_mode": root_mode if root_mode != common.SENTINEL_NONE else None, "outgroup": outgroup}
    args.out_meta.write_text(json.dumps(meta, indent=2))
    print(f"  [{stage}] -> {args.out}")


if __name__ == "__main__":
    main()
