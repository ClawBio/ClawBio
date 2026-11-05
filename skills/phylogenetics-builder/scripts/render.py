#!/usr/bin/env python3
"""Stage 6: parse the rooted Newick tree and render outputs.

Ports 1.x's "Parse and write outputs" section: `newick_to_table()` extracts
per-node branch lengths/support into a table (written as CSV here instead
of an in-memory DataFrame), and `draw_tree()` renders a proportional
phylogram via Bio.Phylo + matplotlib (with its own placeholder-image
fallback on render failure — never fails the step).
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd

import common


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--input", required=True, type=Path, help="Rooted Newick tree from root step")
    p.add_argument("--out-tree", required=True, type=Path, help="Output copy of the final Newick tree")
    p.add_argument("--out-table", required=True, type=Path, help="Output branch-support CSV path")
    p.add_argument("--out-figure", required=True, type=Path, help="Output phylogram PNG path")
    args = p.parse_args()

    args.out_tree.parent.mkdir(parents=True, exist_ok=True)
    args.out_table.parent.mkdir(parents=True, exist_ok=True)
    args.out_figure.parent.mkdir(parents=True, exist_ok=True)

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    newick = args.input.read_text().strip()
    print(f"[render] input={args.input}")

    try:
        table = common.newick_to_table(newick)
    except Exception as exc:
        print(f"Error parsing Newick: {exc}", file=sys.stderr)
        sys.exit(1)

    args.out_tree.write_text(newick + "\n")
    pd.DataFrame(table).to_csv(args.out_table, index=False)
    common.draw_tree(newick, args.out_figure)

    print(f"  nodes={len(table)} -> {args.out_tree}, {args.out_table}, {args.out_figure}")


if __name__ == "__main__":
    main()
