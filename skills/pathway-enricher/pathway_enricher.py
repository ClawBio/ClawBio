#!/usr/bin/env python3
"""
pathway_enricher.py — ClawBio Pathway Enricher skill
=====================================================
Gene-set pathway enrichment analysis using the Enrichr REST API.
Queries KEGG, GO (BP/MF/CC), Reactome, and WikiPathways in a single run,
then produces publication-quality bubble/bar charts and a Markdown report.

Usage:
    # Demo mode (built-in AD gene set)
    python skills/pathway-enricher/pathway_enricher.py --demo

    # Run with a gene list file
    python skills/pathway-enricher/pathway_enricher.py \\
        --input my_genes.txt --output results/

    # Specify which databases to query
    python skills/pathway-enricher/pathway_enricher.py \\
        --input my_genes.txt --output results/ --databases kegg go_bp reactome
"""

from __future__ import annotations

import argparse
import json
import logging
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

import requests

# ---------------------------------------------------------------------------
# Optional visualisation stack
# ---------------------------------------------------------------------------
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import matplotlib.patches as mpatches
    import numpy as np
    import pandas as pd

    HAS_VIZ = True
except ImportError:
    HAS_VIZ = False

# ---------------------------------------------------------------------------
# Project root → clawbio.common
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

try:
    from clawbio.common.reproducibility import (  # noqa: E402
        ReproCommand,
        ReproPath,
        write_checksums,
        write_environment_yml,
        write_portable_commands_sh,
    )

    HAS_REPRO = True
except (ImportError, Exception):
    HAS_REPRO = False

# ---------------------------------------------------------------------------
# Logging
# ---------------------------------------------------------------------------
logging.basicConfig(level=logging.INFO, format="%(levelname)s  %(message)s")
logger = logging.getLogger("pathway-enricher")

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)

ENRICHR_BASE = "https://maayanlab.cloud/Enrichr"

# Library → friendly display name
LIBRARIES: dict[str, str] = {
    "kegg": ("KEGG_2021_Human", "KEGG 2021 Human"),
    "go_bp": ("GO_Biological_Process_2023", "GO Biological Process 2023"),
    "go_mf": ("GO_Molecular_Function_2023", "GO Molecular Function 2023"),
    "go_cc": ("GO_Cellular_Component_2023", "GO Cellular Component 2023"),
    "reactome": ("Reactome_2022", "Reactome 2022"),
    "wikipathways": ("WikiPathways_2023_Human", "WikiPathways 2023"),
}

# Demo gene list — Alzheimer's Disease hallmark genes
DEMO_GENES: list[str] = [
    "APP", "APOE", "PSEN1", "PSEN2", "BIN1", "CLU",
    "TREM2", "ABCA7", "CR1", "MS4A6A", "PICALM", "SORL1",
    "EPHA1", "CD33", "CD2AP", "FERMT2", "SLC24A4", "RIN3",
    "INPP5D", "MEF2C", "NME8", "ZCWPW1", "CELF1", "HLA-DRB5",
    "PTK2B",
]

ADJ_P_THRESHOLD = 0.05
TOP_N_PER_LIB = 15          # rows shown per library in report
RATE_LIMIT_DELAY = 0.5      # seconds between library queries


# ===========================================================================
# Enrichr API helpers
# ===========================================================================

def _post_gene_list(genes: list[str], description: str = "ClawBio gene set") -> str:
    """Upload a gene list to Enrichr and return the userListId."""
    payload = {
        "list": (None, "\n".join(genes)),
        "description": (None, description),
    }
    resp = requests.post(
        f"{ENRICHR_BASE}/addList",
        files=payload,
        timeout=30,
    )
    resp.raise_for_status()
    data = resp.json()
    user_list_id = str(data["userListId"])
    logger.info("Uploaded gene list to Enrichr → userListId=%s", user_list_id)
    return user_list_id


def _query_library(user_list_id: str, library: str) -> list[dict[str, Any]]:
    """
    Query one Enrichr library and return a list of result dicts.

    Enrichr result row indices:
        0  rank
        1  term name
        2  p-value
        3  z-score
        4  combined score
        5  overlapping genes (list)
        6  adjusted p-value
        7  old p-value
        8  old adjusted p-value
    """
    url = f"{ENRICHR_BASE}/enrich"
    params = {"userListId": user_list_id, "backgroundType": library}
    resp = requests.get(url, params=params, timeout=30)
    resp.raise_for_status()
    raw = resp.json().get(library, [])
    results: list[dict[str, Any]] = []
    for row in raw:
        results.append(
            {
                "rank": row[0],
                "term": row[1],
                "pvalue": row[2],
                "zscore": row[3],
                "combined_score": row[4],
                "genes": row[5],
                "adj_pvalue": row[6],
                "gene_count": len(row[5]),
            }
        )
    return results


def run_enrichr(genes: list[str], database_keys: list[str]) -> dict[str, list[dict]]:
    """Upload genes to Enrichr and query the requested databases."""
    user_list_id = _post_gene_list(genes)
    results: dict[str, list[dict]] = {}
    for key in database_keys:
        lib_id, lib_name = LIBRARIES[key]
        logger.info("Querying %s …", lib_name)
        try:
            rows = _query_library(user_list_id, lib_id)
            results[key] = rows
            logger.info("  → %d terms returned", len(rows))
        except Exception as exc:  # noqa: BLE001
            logger.warning("  ✗ %s failed: %s", lib_name, exc)
            results[key] = []
        time.sleep(RATE_LIMIT_DELAY)
    return results


# ===========================================================================
# Visualisation
# ===========================================================================

_PALETTE = [
    "#6366F1",  # indigo
    "#EC4899",  # pink
    "#10B981",  # emerald
    "#F59E0B",  # amber
    "#3B82F6",  # blue
    "#EF4444",  # red
]


def _make_bubble_chart(
    df: "pd.DataFrame",
    title: str,
    out_path: Path,
    top_n: int = 20,
) -> None:
    """Horizontal bubble chart: combined_score on x, term on y, bubble size = gene_count."""
    if df.empty:
        return

    plot_df = (
        df[df["adj_pvalue"] < ADJ_P_THRESHOLD]
        .sort_values("combined_score", ascending=False)
        .head(top_n)
    )
    if plot_df.empty:
        plot_df = df.sort_values("combined_score", ascending=False).head(top_n)

    # Shorten long term names
    plot_df = plot_df.copy()
    plot_df["term_short"] = plot_df["term"].apply(
        lambda t: t[:55] + "…" if len(t) > 55 else t
    )
    plot_df = plot_df.sort_values("combined_score", ascending=True)

    fig, ax = plt.subplots(figsize=(12, max(4, 0.45 * len(plot_df))))
    scatter = ax.scatter(
        plot_df["combined_score"],
        range(len(plot_df)),
        s=plot_df["gene_count"] * 25,
        c=-np.log10(plot_df["adj_pvalue"].clip(lower=1e-300)),
        cmap="viridis",
        alpha=0.85,
        linewidths=0.5,
        edgecolors="white",
        zorder=3,
    )
    ax.set_yticks(range(len(plot_df)))
    ax.set_yticklabels(plot_df["term_short"], fontsize=9)
    ax.set_xlabel("Combined Score (log p × z-score)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", pad=14)
    ax.grid(axis="x", linestyle="--", alpha=0.4, zorder=0)
    ax.spines[["top", "right"]].set_visible(False)

    cbar = plt.colorbar(scatter, ax=ax, pad=0.01)
    cbar.set_label("−log₁₀(adj. p-value)", fontsize=9)

    # Size legend
    sizes = [1, 5, 10]
    handles = [
        plt.scatter([], [], s=s * 25, c="#6366F1", alpha=0.8, label=f"{s} genes")
        for s in sizes
    ]
    ax.legend(
        handles=handles,
        title="Gene count",
        loc="lower right",
        fontsize=8,
        title_fontsize=8,
        framealpha=0.7,
    )

    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved bubble chart → %s", out_path)


def _make_bar_chart_summary(
    all_results: dict[str, list[dict]],
    database_keys: list[str],
    out_path: Path,
    top_n: int = 8,
) -> None:
    """Multi-panel horizontal bar chart: one panel per database, bars = −log₁₀(adj. p)."""
    panels = [k for k in database_keys if all_results.get(k)]
    if not panels:
        return

    n = len(panels)
    fig, axes = plt.subplots(1, n, figsize=(7 * n, 8), sharey=False)
    if n == 1:
        axes = [axes]

    for ax, key, color in zip(axes, panels, _PALETTE):
        _, lib_name = LIBRARIES[key]
        df = pd.DataFrame(all_results[key])
        if df.empty:
            ax.set_visible(False)
            continue

        sig = df[df["adj_pvalue"] < ADJ_P_THRESHOLD].sort_values("adj_pvalue").head(top_n)
        if sig.empty:
            sig = df.sort_values("adj_pvalue").head(top_n)

        sig = sig.copy()
        sig["term_short"] = sig["term"].apply(lambda t: t[:40] + "…" if len(t) > 40 else t)
        sig["neg_log_p"] = -np.log10(sig["adj_pvalue"].clip(lower=1e-300))
        sig = sig.sort_values("neg_log_p", ascending=True)

        bars = ax.barh(
            range(len(sig)),
            sig["neg_log_p"],
            color=color,
            alpha=0.85,
            edgecolor="white",
            linewidth=0.6,
        )
        ax.set_yticks(range(len(sig)))
        ax.set_yticklabels(sig["term_short"], fontsize=8)
        ax.set_xlabel("−log₁₀(adj. p-value)", fontsize=10)
        ax.set_title(lib_name, fontsize=11, fontweight="bold")
        ax.axvline(x=-np.log10(ADJ_P_THRESHOLD), color="red", linestyle="--", alpha=0.7, linewidth=1.2)
        ax.grid(axis="x", linestyle="--", alpha=0.35, zorder=0)
        ax.spines[["top", "right"]].set_visible(False)

    fig.suptitle(
        "Top Enriched Pathways per Database  (red dashed = adj. p = 0.05)",
        fontsize=13,
        fontweight="bold",
        y=1.01,
    )
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved bar summary → %s", out_path)


def _make_heatmap(
    all_results: dict[str, list[dict]],
    genes: list[str],
    database_keys: list[str],
    out_path: Path,
    top_n: int = 10,
) -> None:
    """Gene × pathway presence heatmap for top terms across all databases."""
    if not HAS_VIZ:
        return

    rows = []
    for key in database_keys:
        df = pd.DataFrame(all_results.get(key, []))
        if df.empty:
            continue
        sig = df[df["adj_pvalue"] < ADJ_P_THRESHOLD].head(top_n)
        if sig.empty:
            sig = df.head(3)
        for _, row in sig.iterrows():
            term_genes = [g.upper() for g in row["genes"]]
            for g in genes:
                rows.append({"pathway": row["term"][:40], "gene": g.upper(), "hit": int(g.upper() in term_genes)})

    if not rows:
        return

    df_heat = pd.DataFrame(rows)
    pivot = df_heat.pivot_table(index="gene", columns="pathway", values="hit", fill_value=0)
    # Keep only genes with at least one hit
    pivot = pivot[pivot.sum(axis=1) > 0]
    if pivot.empty:
        return

    fig, ax = plt.subplots(figsize=(min(24, 0.7 * len(pivot.columns) + 2), max(4, 0.35 * len(pivot))))
    im = ax.imshow(pivot.values, aspect="auto", cmap="YlOrRd", vmin=0, vmax=1)
    ax.set_xticks(range(len(pivot.columns)))
    ax.set_xticklabels(pivot.columns, rotation=45, ha="right", fontsize=7)
    ax.set_yticks(range(len(pivot)))
    ax.set_yticklabels(pivot.index, fontsize=8)
    ax.set_title("Gene × Top Pathway Membership Heatmap", fontsize=12, fontweight="bold")
    plt.colorbar(im, ax=ax, label="In pathway (1=yes, 0=no)")
    plt.tight_layout()
    plt.savefig(out_path, dpi=150, bbox_inches="tight")
    plt.close(fig)
    logger.info("Saved heatmap → %s", out_path)


# ===========================================================================
# Report generation
# ===========================================================================

def _fmt_genes(genes: list[str]) -> str:
    return ", ".join(f"`{g}`" for g in genes[:10]) + (
        f" … (+{len(genes) - 10} more)" if len(genes) > 10 else ""
    )


def _table_md(rows: list[dict], top_n: int = TOP_N_PER_LIB) -> str:
    if not rows:
        return "_No results returned._\n"
    header = "| Rank | Term | adj. p | Combined Score | Genes (overlap) |"
    sep = "|------|------|--------|----------------|-----------------|"
    lines = [header, sep]
    for r in rows[:top_n]:
        g_str = ", ".join(r["genes"][:6])
        if len(r["genes"]) > 6:
            g_str += f" +{len(r['genes'])-6}"
        lines.append(
            f"| {r['rank']} | {r['term'][:60]} | {r['adj_pvalue']:.2e} | {r['combined_score']:.1f} | {g_str} |"
        )
    return "\n".join(lines) + "\n"


def _figure_md(fig_path: Path, output_dir: Path, caption: str) -> str:
    rel = fig_path.relative_to(output_dir)
    return f"![{caption}]({rel})\n*{caption}*\n"


def build_report(
    genes: list[str],
    all_results: dict[str, list[dict]],
    database_keys: list[str],
    figures: dict[str, Path],
    output_dir: Path,
    run_ts: str,
) -> str:
    _, lib_names = zip(*[LIBRARIES[k] for k in database_keys]) if database_keys else ([], [])

    sig_counts = {
        k: sum(1 for r in all_results.get(k, []) if r["adj_pvalue"] < ADJ_P_THRESHOLD)
        for k in database_keys
    }
    total_sig = sum(sig_counts.values())

    sections = [
        f"# 🔬 Pathway Enrichment Report\n",
        f"**Generated:** {run_ts}  \n**Tool:** ClawBio Pathway Enricher v0.1.0  \n**API:** Enrichr (MaayenLab, Icahn School of Medicine at Mount Sinai)\n",
        f"> ⚠️ **Disclaimer:** {DISCLAIMER}\n",
        "---\n",
        "## Summary\n",
        f"| Parameter | Value |\n|-----------|-------|\n"
        f"| Gene list size | **{len(genes)} genes** |\n"
        f"| Databases queried | {', '.join(lib_names)} |\n"
        f"| Total significant pathways (adj. p < 0.05) | **{total_sig}** |\n",
        "\n**Input genes:** " + _fmt_genes(genes) + "\n",
        "---\n",
    ]

    # Summary figure
    if "bar_summary" in figures:
        sections.append("## Overview — Top Pathways per Database\n")
        sections.append(_figure_md(figures["bar_summary"], output_dir, "Top enriched pathways per database"))
        sections.append("\n")

    # Per-database sections
    for key in database_keys:
        _, lib_name = LIBRARIES[key]
        rows = all_results.get(key, [])
        n_sig = sig_counts[key]
        sections.append(f"---\n\n## {lib_name}\n")
        sections.append(f"**Significant terms (adj. p < 0.05):** {n_sig}\n\n")
        if f"bubble_{key}" in figures:
            sections.append(_figure_md(figures[f"bubble_{key}"], output_dir, f"Bubble chart — {lib_name}"))
            sections.append("\n")
        sections.append("### Top Results\n\n")
        sections.append(_table_md(rows))
        sections.append("\n")

    # Heatmap
    if "heatmap" in figures:
        sections.append("---\n\n## Gene × Pathway Membership\n")
        sections.append(_figure_md(figures["heatmap"], output_dir, "Gene × pathway heatmap"))
        sections.append("\n")

    sections.append("---\n\n## Methods\n")
    sections.append(
        "Enrichment analysis performed using the [Enrichr API](https://maayanlab.cloud/Enrichr). "
        "Combined score = −log(p-value) × z-score (Enrichr's default ranking metric). "
        "P-values adjusted using the Benjamini–Hochberg method (performed server-side by Enrichr). "
        f"Significance threshold: adj. p < {ADJ_P_THRESHOLD}.\n"
    )
    sections.append("\n---\n\n> " + DISCLAIMER + "\n")

    return "\n".join(sections)


# ===========================================================================
# CSV writers
# ===========================================================================

def save_tables(
    all_results: dict[str, list[dict]],
    database_keys: list[str],
    tables_dir: Path,
) -> None:
    """Write per-library CSV files."""
    if not HAS_VIZ:
        # Write minimal CSVs without pandas
        for key in database_keys:
            rows = all_results.get(key, [])
            out = tables_dir / f"{key}_enrichment.csv"
            with out.open("w") as fh:
                fh.write("rank,term,pvalue,adj_pvalue,zscore,combined_score,gene_count,genes\n")
                for r in rows:
                    genes_str = "|".join(r["genes"])
                    fh.write(
                        f"{r['rank']},{r['term']!r},{r['pvalue']},{r['adj_pvalue']},"
                        f"{r['zscore']},{r['combined_score']},{r['gene_count']},{genes_str!r}\n"
                    )
        return

    for key in database_keys:
        rows = all_results.get(key, [])
        out = tables_dir / f"{key}_enrichment.csv"
        if not rows:
            pd.DataFrame().to_csv(out, index=False)
            continue
        df = pd.DataFrame(rows)
        df["genes"] = df["genes"].apply(lambda g: "|".join(g))
        df.to_csv(out, index=False)
        logger.info("Saved table → %s", out)


# ===========================================================================
# Result JSON
# ===========================================================================

def build_result_json(
    genes: list[str],
    all_results: dict[str, list[dict]],
    database_keys: list[str],
    run_ts: str,
) -> dict[str, Any]:
    top_per_db: dict[str, Any] = {}
    for key in database_keys:
        rows = all_results.get(key, [])
        sig = [r for r in rows if r["adj_pvalue"] < ADJ_P_THRESHOLD]
        top = (sig or rows)[:5]
        top_per_db[key] = [{"term": r["term"], "adj_pvalue": r["adj_pvalue"], "combined_score": r["combined_score"]} for r in top]

    return {
        "tool": "pathway-enricher",
        "version": "0.1.0",
        "run_timestamp": run_ts,
        "input_genes": genes,
        "databases_queried": database_keys,
        "top_results": top_per_db,
        "significant_counts": {
            k: sum(1 for r in all_results.get(k, []) if r["adj_pvalue"] < ADJ_P_THRESHOLD)
            for k in database_keys
        },
    }


# ===========================================================================
# Main entry point
# ===========================================================================

def parse_gene_file(path: Path) -> list[str]:
    """Read a gene-list file (one per line, or comma-separated). Skip comments."""
    genes: list[str] = []
    with path.open() as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            for part in line.split(","):
                g = part.strip().upper()
                if g:
                    genes.append(g)
    # Deduplicate while preserving order
    seen: set[str] = set()
    unique = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            unique.append(g)
    return unique


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        prog="pathway_enricher",
        description="Gene-set pathway enrichment via Enrichr (KEGG, GO, Reactome, WikiPathways).",
    )
    parser.add_argument("--input", type=Path, help="Path to gene list file (.txt or .csv)")
    parser.add_argument("--output", type=Path, default=Path("output/pathway-enricher"), help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with built-in Alzheimer's gene set")
    parser.add_argument(
        "--databases",
        nargs="+",
        choices=list(LIBRARIES),
        default=list(LIBRARIES),
        help="Databases to query (default: all)",
    )
    parser.add_argument("--top-n", type=int, default=TOP_N_PER_LIB, help="Top terms per library in report")
    args = parser.parse_args(argv)

    # ----- Resolve genes -------------------------------------------------------
    if args.demo:
        genes = DEMO_GENES
        logger.info("Demo mode — using %d Alzheimer's disease hallmark genes", len(genes))
        if args.output == Path("output/pathway-enricher"):
            ts_str = datetime.now(timezone.utc).strftime("%Y%m%dT%H%M%SZ")
            args.output = Path(f"output/pathway-enricher-demo-{ts_str}")
    elif args.input:
        if not args.input.exists():
            logger.error("Input file not found: %s", args.input)
            return 1
        genes = parse_gene_file(args.input)
        logger.info("Parsed %d unique genes from %s", len(genes), args.input)
    else:
        parser.print_help()
        return 1

    if not genes:
        logger.error("No genes found in input.")
        return 1

    # ----- Output dirs ---------------------------------------------------------
    output_dir = args.output
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    run_ts = datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ")
    database_keys = args.databases

    # ----- Enrichr query -------------------------------------------------------
    logger.info("Submitting %d genes to Enrichr …", len(genes))
    try:
        all_results = run_enrichr(genes, database_keys)
    except requests.exceptions.ConnectionError:
        logger.error("Cannot reach Enrichr API. Check your internet connection.")
        return 1
    except requests.exceptions.HTTPError as exc:
        logger.error("Enrichr API error: %s", exc)
        return 1

    # ----- Tables --------------------------------------------------------------
    save_tables(all_results, database_keys, output_dir / "tables")

    # ----- Figures -------------------------------------------------------------
    figures: dict[str, Path] = {}
    if HAS_VIZ:
        import pandas as pd  # noqa: F811 (may not have been imported at top level)
        import numpy as np  # noqa: F811

        for key in database_keys:
            rows = all_results.get(key, [])
            if rows:
                df = pd.DataFrame(rows)
                fig_path = output_dir / "figures" / f"bubble_chart_{key}.png"
                _, lib_name = LIBRARIES[key]
                _make_bubble_chart(df, f"Enriched Pathways — {lib_name}", fig_path, top_n=20)
                figures[f"bubble_{key}"] = fig_path

        bar_path = output_dir / "figures" / "bar_chart_summary.png"
        _make_bar_chart_summary(all_results, database_keys, bar_path, top_n=8)
        figures["bar_summary"] = bar_path

        heatmap_path = output_dir / "figures" / "heatmap_top_pathways.png"
        _make_heatmap(all_results, genes, database_keys, heatmap_path, top_n=6)
        figures["heatmap"] = heatmap_path
    else:
        logger.warning("matplotlib/numpy/pandas not installed — skipping figures. Install with: pip install matplotlib numpy pandas")

    # ----- Markdown report -----------------------------------------------------
    report_md = build_report(genes, all_results, database_keys, figures, output_dir, run_ts)
    (output_dir / "report.md").write_text(report_md)
    logger.info("Saved report → %s/report.md", output_dir)

    # ----- result.json ---------------------------------------------------------
    result = build_result_json(genes, all_results, database_keys, run_ts)
    (output_dir / "result.json").write_text(json.dumps(result, indent=2))

    # ----- Reproducibility pack ------------------------------------------------
    db_arg = " ".join(database_keys)
    input_arg = "--demo" if args.demo else str(args.input)
    if HAS_REPRO:
        try:
            skill_script = Path("skills/pathway-enricher/pathway_enricher.py")
            repro_cmd = ReproCommand(
                script_path=skill_script,
                args=[
                    "--input", input_arg,
                    "--output", str(output_dir),
                    "--databases", *database_keys,
                ],
                comment="Pathway Enricher reproducibility command",
            )
            write_portable_commands_sh(output_dir, repro_cmd, repo_root=_PROJECT_ROOT)
            write_environment_yml(
                output_dir,
                env_name="clawbio-pathway-enricher",
                pip_deps=["requests>=2.28", "matplotlib>=3.5", "numpy>=1.23", "pandas>=1.5"],
            )
            result_json_path = output_dir / "result.json"
            write_checksums(
                [result_json_path],
                output_dir,
                anchor=output_dir,
            )
        except Exception as _repro_exc:  # noqa: BLE001
            logger.warning("Reproducibility pack skipped: %s", _repro_exc)
            # Fallback minimal commands.sh
            _cmd_sh = output_dir / "reproducibility" / "commands.sh"
            _cmd_sh.write_text(
                f"#!/usr/bin/env bash\n"
                f"python skills/pathway-enricher/pathway_enricher.py "
                f"{input_arg} --output {output_dir} --databases {db_arg}\n"
            )
    else:
        # Minimal commands.sh fallback
        (output_dir / "reproducibility" / "commands.sh").write_text(
            f"#!/usr/bin/env bash\n"
            f"python skills/pathway-enricher/pathway_enricher.py "
            f"{input_arg} --output {output_dir} --databases {db_arg}\n"
        )

    # ----- Done ----------------------------------------------------------------
    logger.info("✅  Pathway enrichment complete → %s", output_dir)
    print(f"\n✅  Done! Results written to: {output_dir}")
    print(f"    report.md  |  tables/  |  figures/  |  result.json")

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
