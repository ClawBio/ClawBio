#!/usr/bin/env python3
"""Stage 4: publication-quality figures from the taxonomy/resistome tables.

Fan-in step: reads the taxonomy_table + resistome_table + who_table emitted
by the taxonomy and resistome steps, produces the three figures the 1.x
script generates:
  1. Taxonomy bar chart (top 20 species)
  2. Resistome heatmap (ARG families x drug class / site)
  3. WHO-critical ARG summary (priority tier pie + drug-class breakdown)

Multi-site (demo) vs single-sample (real) rendering is auto-detected from
the table's columns (`Lima`/`Cusco`/`Iquitos` present => multi-site), the
same branching the 1.x plotting functions already did internally — no
decision is made in pipeline.yaml.
"""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd

import common

SITE_COLUMNS = ["Lima", "Cusco", "Iquitos"]


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--taxonomy-table", required=True, type=Path)
    p.add_argument("--resistome-table", required=True, type=Path)
    p.add_argument("--who-table", required=True, type=Path)
    p.add_argument("--out-taxonomy", required=True, type=Path)
    p.add_argument("--out-resistome", required=True, type=Path)
    p.add_argument("--out-who", required=True, type=Path)
    args = p.parse_args()

    for out in (args.out_taxonomy, args.out_resistome, args.out_who):
        out.parent.mkdir(parents=True, exist_ok=True)

    print("[figures] generating publication-quality figures (300 dpi)")

    tax_df = pd.read_csv(args.taxonomy_table, sep="\t")
    res_df = pd.read_csv(args.resistome_table, sep="\t")
    who_df = pd.read_csv(args.who_table, sep="\t")

    site_cols = SITE_COLUMNS if all(c in tax_df.columns for c in SITE_COLUMNS) else None
    resistome_failed = common.read_status(res_df) == "failed"

    common.plot_taxonomy_barchart(tax_df, args.out_taxonomy, top_n=20, site_columns=site_cols)

    if resistome_failed:
        print("  [resistome heatmap] resistome step failed — rendering placeholder")
        common.plot_resistome_heatmap(pd.DataFrame(), args.out_resistome, site_columns=site_cols)
        common.plot_who_critical_args(pd.DataFrame(), args.out_who)
    else:
        common.plot_resistome_heatmap(res_df, args.out_resistome, site_columns=site_cols)
        common.plot_who_critical_args(who_df, args.out_who)

    print(f"  Saved: {args.out_taxonomy}")
    print(f"  Saved: {args.out_resistome}")
    print(f"  Saved: {args.out_who}")


if __name__ == "__main__":
    main()
