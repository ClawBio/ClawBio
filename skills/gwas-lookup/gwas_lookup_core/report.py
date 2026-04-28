"""
report.py — Generate markdown report, CSV tables, and figures.
"""

from __future__ import annotations

import csv
import math
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.report import DISCLAIMER


def _fmt_pval(pval: Any) -> str:
    if pval is None:
        return "N/A"
    try:
        p = float(pval)
        if p == 0:
            return "0"
        if p < 1e-300:
            return f"<1e-300"
        return f"{p:.2e}"
    except (TypeError, ValueError):
        return str(pval)


def _fmt_float(val: Any, decimals: int = 3) -> str:
    if val is None:
        return "N/A"
    try:
        return f"{float(val):.{decimals}f}"
    except (TypeError, ValueError):
        return str(val)


# ---------------------------------------------------------------------------
# Markdown report
# ---------------------------------------------------------------------------


def generate_markdown(variant: dict, merged: dict) -> str:
    """Generate a full markdown report from resolved variant + merged results."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    rsid = variant.get("rsid", "unknown")

    lines = [
        f"# GWAS Lookup Report: {rsid}",
        "",
        f"**Date**: {now}",
        f"**Variant**: {rsid}",
        f"**Location (GRCh38)**: chr{variant.get('chr', '?')}:{variant.get('pos_grch38', '?')}",
        f"**Location (GRCh37)**: chr{variant.get('chr', '?')}:{variant.get('pos_grch37', '?')}",
        f"**Alleles**: {variant.get('allele_string', '?')}",
        f"**Class**: {variant.get('var_class', '?')}",
        f"**Consequence**: {variant.get('most_severe_consequence', '?')}",
        f"**Minor allele**: {variant.get('minor_allele', '?')} (MAF: {_fmt_float(variant.get('maf'), 4)})",
        "",
        "---",
        "",
    ]

    # --- VEP consequences ---
    vep = variant.get("vep_info", {})
    if vep.get("status") == "ok" and vep.get("consequences"):
        lines.append("## Variant Effect Prediction (VEP)")
        lines.append("")
        lines.append("| Gene | Consequence | Impact | SIFT | PolyPhen |")
        lines.append("|------|------------|--------|------|----------|")
        seen = set()
        for c in vep["consequences"]:
            gene = c.get("gene_symbol", "")
            terms = ", ".join(c.get("consequence_terms", []))
            key = (gene, terms)
            if key in seen:
                continue
            seen.add(key)
            lines.append(
                f"| {gene} | {terms} | {c.get('impact', '')} "
                f"| {c.get('sift', '')} | {c.get('polyphen', '')} |"
            )
        lines.append("")

    # --- GWAS associations ---
    gwas = merged.get("gwas_associations", [])
    summary = merged.get("summary", {})
    lines.append(f"## GWAS Associations ({summary.get('total_gwas', 0)} total, "
                 f"{summary.get('total_gwas_significant', 0)} genome-wide significant)")
    lines.append("")
    if gwas:
        lines.append("| Trait | P-value | OR/Beta | Risk Allele | Source | Study |")
        lines.append("|-------|---------|---------|-------------|--------|-------|")
        for a in gwas[:50]:
            or_beta = _fmt_float(a.get("or_beta")) if a.get("or_beta") else _fmt_float(a.get("beta"))
            gws = " **" if a.get("genome_wide_significant") else ""
            lines.append(
                f"| {a.get('trait', '')}{gws} | {_fmt_pval(a.get('pval'))} "
                f"| {or_beta} | {a.get('risk_allele', '')} "
                f"| {a.get('source', '')} | {a.get('study', '')} |"
            )
        if len(gwas) > 50:
            lines.append(f"| ... | *{len(gwas) - 50} more in CSV* | | | | |")
    else:
        lines.append("*No GWAS associations found.*")
    lines.append("")

    # --- PheWAS: UKB ---
    phewas = merged.get("phewas", {})
    for biobank, label in [("ukb", "UKB-TOPMed"), ("finngen", "FinnGen"), ("bbj", "Biobank Japan")]:
        phewas_list = phewas.get(biobank, [])
        sig_count = sum(1 for p in phewas_list if p.get("genome_wide_significant"))
        lines.append(f"## PheWAS: {label} ({len(phewas_list)} phenotypes, {sig_count} significant)")
        lines.append("")
        if phewas_list:
            lines.append("| Phenotype | Category | P-value | Beta | MAF |")
            lines.append("|-----------|----------|---------|------|-----|")
            for p in phewas_list[:30]:
                gws = " **" if p.get("genome_wide_significant") else ""
                lines.append(
                    f"| {p.get('phenostring', p.get('phenocode', ''))}{gws} "
                    f"| {p.get('category', '')} | {_fmt_pval(p.get('pval'))} "
                    f"| {_fmt_float(p.get('beta'))} | {_fmt_float(p.get('maf'), 4)} |"
                )
            if len(phewas_list) > 30:
                lines.append(f"| ... | *{len(phewas_list) - 30} more in CSV* | | | |")
        else:
            lines.append(f"*No {label} PheWAS results found.*")
        lines.append("")

    # --- eQTL associations ---
    eqtls = merged.get("eqtl_associations", [])
    lines.append(f"## eQTL Associations ({summary.get('total_eqtls', 0)} total)")
    lines.append("")
    if eqtls:
        lines.append("| Gene | Tissue | P-value | Effect Size | Source |")
        lines.append("|------|--------|---------|-------------|--------|")
        for e in eqtls[:30]:
            lines.append(
                f"| {e.get('gene', '')} | {e.get('tissue_name', e.get('tissue', ''))} "
                f"| {_fmt_pval(e.get('pval'))} | {_fmt_float(e.get('effect_size'))} "
                f"| {e.get('source', '')} |"
            )
        if len(eqtls) > 30:
            lines.append(f"| ... | *{len(eqtls) - 30} more in CSV* | | | |")
    else:
        lines.append("*No eQTL associations found.*")
    lines.append("")

    # --- Credible sets ---
    cred_sets = merged.get("credible_sets", [])
    if cred_sets:
        lines.append(f"## Fine-Mapping / Credible Sets ({len(cred_sets)} sets)")
        lines.append("")
        lines.append("| Trait | Study | Post. Prob. | P-value | 95% CS | 99% CS |")
        lines.append("|-------|-------|-------------|---------|--------|--------|")
        for cs in cred_sets:
            lines.append(
                f"| {cs.get('trait', '')} | {cs.get('study_id', '')} "
                f"| {_fmt_float(cs.get('posterior_probability'))} "
                f"| {_fmt_pval(cs.get('pval'))} "
                f"| {'Yes' if cs.get('is_95_credible') else 'No'} "
                f"| {'Yes' if cs.get('is_99_credible') else 'No'} |"
            )
        lines.append("")

    # --- Data sources ---
    sources = merged.get("data_sources", {})
    lines.append("## Data Sources")
    lines.append("")
    lines.append("| Source | Status |")
    lines.append("|--------|--------|")
    for name, info in sources.items():
        status = info.get("status", "unknown")
        badge = "OK" if status == "ok" else f"WARNING: {status}"
        msg = info.get("message", "")
        if msg and status != "ok":
            badge += f" — {msg[:80]}"
        lines.append(f"| {name} | {badge} |")
    lines.append("")

    # --- Methods + disclaimer ---
    lines.append("## Methods")
    lines.append("")
    lines.append("This report was generated by the ClawBio GWAS Lookup skill, which queries")
    lines.append("9 genomic databases in parallel: Ensembl, GWAS Catalog, Open Targets,")
    lines.append("UKB-TOPMed PheWeb, FinnGen, Biobank Japan PheWeb, GTEx, EBI eQTL Catalogue,")
    lines.append("and LocusZoom PortalDev.")
    lines.append("")
    lines.append("- **Genome-wide significance threshold**: p < 5 x 10^-8")
    lines.append("- **Coordinate system**: GRCh38 (primary), GRCh37 for BBJ")
    lines.append("- **Cache**: 24-hour local file cache")
    lines.append("")
    lines.append("## Disclaimer")
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CSV tables
# ---------------------------------------------------------------------------


def _write_csv(filepath: Path, rows: list[dict]):
    if not rows:
        return
    filepath.parent.mkdir(parents=True, exist_ok=True)
    # Collect all fieldnames across all rows (rows may have different keys)
    fieldnames = list(dict.fromkeys(k for row in rows for k in row.keys()))
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=fieldnames, extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def write_tables(output_dir: Path, merged: dict):
    """Write CSV tables from merged results."""
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    gwas = merged.get("gwas_associations", [])
    if gwas:
        _write_csv(tables_dir / "gwas_associations.csv", gwas)

    phewas = merged.get("phewas", {})
    for key, filename in [("ukb", "phewas_ukb.csv"), ("finngen", "phewas_finngen.csv"), ("bbj", "phewas_bbj.csv")]:
        if phewas.get(key):
            _write_csv(tables_dir / filename, phewas[key])

    eqtls = merged.get("eqtl_associations", [])
    if eqtls:
        _write_csv(tables_dir / "eqtl_associations.csv", eqtls)

    cred_sets = merged.get("credible_sets", [])
    if cred_sets:
        _write_csv(tables_dir / "credible_sets.csv", cred_sets)


# ---------------------------------------------------------------------------
# Figures
# ---------------------------------------------------------------------------


def generate_figures(output_dir: Path, merged: dict, variant: dict):
    """Generate matplotlib figures. Skips gracefully if matplotlib is unavailable."""
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    fig_dir = output_dir / "figures"
    fig_dir.mkdir(parents=True, exist_ok=True)

    # --- GWAS traits dot plot (top 15 by -log10(p)) ---
    gwas = merged.get("gwas_associations", [])
    plotable = []
    for a in gwas:
        pval = a.get("pval")
        if pval is not None and pval > 0:
            plotable.append((a.get("trait", "?")[:40], -math.log10(pval)))

    if plotable:
        plotable.sort(key=lambda x: x[1], reverse=True)
        plotable = plotable[:15]
        traits, logp = zip(*plotable)

        fig, ax = plt.subplots(figsize=(8, max(4, len(traits) * 0.4)))
        colors = ["#d32f2f" if lp > -math.log10(5e-8) else "#1976d2" for lp in logp]
        ax.barh(range(len(traits)), logp, color=colors)
        ax.set_yticks(range(len(traits)))
        ax.set_yticklabels(traits, fontsize=8)
        ax.set_xlabel("-log10(p-value)")
        ax.set_title(f"Top GWAS Associations: {variant.get('rsid', '')}")
        ax.axvline(-math.log10(5e-8), color="red", linestyle="--", alpha=0.5, label="GWS (5e-8)")
        ax.legend(fontsize=7)
        ax.invert_yaxis()
        plt.tight_layout()
        plt.savefig(fig_dir / "gwas_traits_dotplot.png", dpi=150)
        plt.close()

    # --- Population allele frequency bar chart ---
    pops = variant.get("vep_info", {}) if variant.get("vep_info", {}).get("status") == "ok" else {}
    ot_result = merged.get("data_sources", {}).get("open_targets", {})
    # Try to get population frequencies from Open Targets variant result
    # (stored separately in the api_results passed through merged)
    # For now, use Ensembl populations if available
    ens_pops = variant.get("populations", [])
    if ens_pops:
        # Group by superpopulation prefix
        pop_freqs = {}
        minor = variant.get("minor_allele", "")
        for p in ens_pops:
            if p.get("allele") == minor:
                pop_name = p.get("population", "")
                freq = p.get("frequency")
                if freq is not None and pop_name:
                    # Take the first entry per population
                    if pop_name not in pop_freqs:
                        pop_freqs[pop_name] = freq

        if len(pop_freqs) >= 2:
            # Show top 10 by frequency
            sorted_pops = sorted(pop_freqs.items(), key=lambda x: x[1], reverse=True)[:10]
            names, freqs = zip(*sorted_pops)

            fig, ax = plt.subplots(figsize=(8, max(3, len(names) * 0.35)))
            ax.barh(range(len(names)), freqs, color="#7b1fa2")
            ax.set_yticks(range(len(names)))
            ax.set_yticklabels(names, fontsize=8)
            ax.set_xlabel(f"Allele Frequency ({minor})")
            ax.set_title(f"Population Frequencies: {variant.get('rsid', '')}")
            ax.invert_yaxis()
            plt.tight_layout()
            plt.savefig(fig_dir / "allele_freq_populations.png", dpi=150)
            plt.close()


# ---------------------------------------------------------------------------
# Reproducibility bundle
# ---------------------------------------------------------------------------


def write_reproducibility(output_dir: Path, variant: dict, skip_apis: list[str]):
    """Write commands.sh and api_versions.json for reproducibility."""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    rsid = variant.get("rsid", "unknown")
    skip_str = ",".join(skip_apis) if skip_apis else ""
    skip_flag = f" --skip {skip_str}" if skip_str else ""

    commands = f"""#!/bin/bash
# Reproduce this GWAS Lookup report
# Generated: {datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")}

python skills/gwas-lookup/gwas_lookup.py \\
  --rsid {rsid}{skip_flag} \\
  --output {output_dir}
"""
    (repro_dir / "commands.sh").write_text(commands)

    import json
    versions = {
        "tool": "ClawBio GWAS Lookup",
        "version": "0.2.0",
        "apis": {
            "ensembl": "https://rest.ensembl.org",
            "gwas_catalog": "https://www.ebi.ac.uk/gwas/rest/api",
            "open_targets": "https://api.platform.opentargets.org/api/v4",
            "pheweb_ukb": "https://pheweb.org/UKB-TOPMed",
            "finngen": "https://r12.finngen.fi",
            "pheweb_bbj": "https://pheweb.jp",
            "gtex": "https://gtexportal.org/api/v2",
            "eqtl_catalogue": "https://www.ebi.ac.uk/eqtl/api/v3",
            "portaldev": "https://portaldev.sph.umich.edu/api/v1",
        },
        "skipped": skip_apis,
        "generated": datetime.now(timezone.utc).isoformat(),
    }
    (repro_dir / "api_versions.json").write_text(json.dumps(versions, indent=2))
