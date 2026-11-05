#!/usr/bin/env python3
"""Stage 5: render the unified Markdown report.

Fan-in step: reads every table and figure produced by the taxonomy,
resistome, functional, and figures steps and assembles them into a single
report.md, mirroring 1.x's generate_report(). The 1.x script's own
reproducibility-bundle writer (commands.sh/environment.yml/checksums.sha256)
is intentionally NOT reproduced here — ClawBio 2.0's `clawbio bundle` verb
now owns that job for the whole run, derived from the event log rather than
duplicated per skill.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
from pathlib import Path

import pandas as pd

import common


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--taxonomy-table", required=True, type=Path)
    p.add_argument("--resistome-table", required=True, type=Path)
    p.add_argument("--who-table", required=True, type=Path)
    p.add_argument("--pathway-table", required=True, type=Path)
    p.add_argument("--taxonomy-figure", required=True, type=Path)
    p.add_argument("--resistome-figure", required=True, type=Path)
    p.add_argument("--who-figure", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path, help="Output report.md path")
    args = p.parse_args()

    args.out.parent.mkdir(parents=True, exist_ok=True)

    print("[report] assembling report.md")

    tax_df = pd.read_csv(args.taxonomy_table, sep="\t")
    res_df = pd.read_csv(args.resistome_table, sep="\t")
    who_df = pd.read_csv(args.who_table, sep="\t")
    pw_df = pd.read_csv(args.pathway_table, sep="\t")

    # Demo runs carry per-site columns (Lima/Cusco/Iquitos) in the taxonomy
    # table; real runs carry a single-sample relative_abundance column.
    is_demo = "Lima" in tax_df.columns
    rgi_failed = common.read_status(res_df) == "failed"
    pw_status = common.read_status(pw_df)
    humann_skipped = pw_status == "skipped"
    humann_failed = pw_status == "failed"

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")

    def fig_ref(name: str, path: Path) -> str:
        try:
            rel = path.relative_to(args.out.parent)
        except ValueError:
            rel = path.name
        return f"![{name}]({rel})"

    # Taxonomy summary
    if "species" in tax_df.columns:
        tax_lines = [f"  - *{row['species']}*" for _, row in tax_df.head(5).iterrows()]
    elif "name" in tax_df.columns:
        tax_lines = [
            f"  - *{row['name']}* ({row.get('relative_abundance', 0):.1f}%)"
            for _, row in tax_df.head(5).iterrows()
        ]
    else:
        tax_lines = ["  - No taxonomy data available"]

    # Resistome summary
    if rgi_failed:
        resistome_section = """## Antimicrobial Resistance Profile

> **ANALYSIS FAILED**: RGI did not produce output. Resistome could not be profiled.
> Re-run with valid inputs or check that the CARD database is correctly installed.

"""
    else:
        n_args = len(res_df)
        n_critical = len(who_df[who_df["who_priority"] == "Critical"]) if not who_df.empty else 0
        n_high = len(who_df[who_df["who_priority"] == "High"]) if not who_df.empty else 0
        resistome_section = f"""## Antimicrobial Resistance Profile

{fig_ref("resistome", args.resistome_figure)}

Total ARG detections: **{n_args}**
- WHO-Critical priority: **{n_critical}**
- WHO-High priority: **{n_high}**

### WHO-Critical ARG Summary

{fig_ref("who_critical", args.who_figure)}

"""

    # Pathway summary
    if humann_failed:
        pathway_section = """## Functional Pathways (HUMAnN3)

> **FUNCTIONAL PROFILING FAILED**: HUMAnN3 did not produce pathway output.
> Check database installation and input file compatibility.

"""
    elif humann_skipped:
        pathway_section = """## Functional Pathways (HUMAnN3)

*Functional profiling was skipped (no database specified or --humann-db skip used).*

"""
    else:
        n_pathways = len(pw_df)
        top_pw = pw_df.iloc[0]
        pw_name = top_pw.get("pathway", top_pw.get("# Pathway", "Unknown"))
        pathway_section = f"""## Functional Pathways (HUMAnN3)

Total pathways detected: **{n_pathways}**

Top pathway: {pw_name}

"""

    report = f"""# Metagenomics Profiling Report

**Date**: {now}
**Tool**: ClawBio Metagenomics Profiler (2.0 port)
**Mode**: {"Demo (pre-computed Peru sewage data)" if is_demo else "Full pipeline"}

---

## Taxonomic Composition

{fig_ref("taxonomy", args.taxonomy_figure)}

### Top Species

{chr(10).join(tax_lines)}

{resistome_section}{pathway_section}---

## Methods

- **Taxonomic classification**: Kraken2 (confidence threshold {common.KRAKEN2_CONFIDENCE}) with Bracken species-level re-estimation (minimum {common.BRACKEN_THRESHOLD} reads)
- **Resistome profiling**: RGI against the CARD database (Perfect + Strict criteria only)
- **Functional profiling**: HUMAnN3 with MetaCyc pathway stratification
- **WHO Priority mapping**: ARGs classified against the WHO Bacterial Priority Pathogens List
- **Figures**: 300 dpi, publication-quality

## Reproducibility

This report was generated by a `clawbio`-driven run. See `.clawbio-run/commands.sh`,
`checksums.sha256`, and `environment.yml` (written by `clawbio bundle`) for the exact
commands, output hashes, and environment for this run.

## Disclaimer

> {common.SAFETY_DISCLAIMER}

## Citations

- Wood, D.E., Lu, J. & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257.
- Lu, J. et al. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- Alcock, B.P. et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51(D1), D419-D430.
- Beghini, F. et al. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife*, 10, e65088.
- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio
"""

    args.out.write_text(report, encoding="utf-8")
    print(f"  Saved: {args.out}")


if __name__ == "__main__":
    main()
