#!/usr/bin/env python3
"""Affinity Proteomics Pipeline — Olink NPX + SomaLogic SomaScan analysis.

Unified analysis for affinity-based proteomics platforms with platform-aware
QC, normalisation, differential abundance testing, and visualisation.

Supported platforms:
  - Olink (Proximity Extension Assay): NPX CSV/Parquet
  - SomaLogic SomaScan (SOMAmer aptamers): ADAT files

References:
    Assarsson et al. (2014) PLOS ONE 9:e95192 (Olink PEA)
    Gold et al. (2010) PLOS ONE 5:e15004 (SOMAmer)

Usage:
    python affinity_proteomics.py --demo --platform olink --output /tmp/olink_demo
    python affinity_proteomics.py --demo --platform somascan --output /tmp/soma_demo
    python affinity_proteomics.py --platform olink --input data.csv --meta meta.csv \
        --group-col Group --contrast "Case,Control" --output results/
"""

from __future__ import annotations

import argparse
import csv
import json
import math
import sys
from dataclasses import dataclass, field
from datetime import datetime, timezone
from pathlib import Path
from typing import Literal

import numpy as np
import pandas as pd

SCRIPT_DIR = Path(__file__).resolve().parent

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)


# ---------------------------------------------------------------------------
# Data structures
# ---------------------------------------------------------------------------
@dataclass
class ProteomicsData:
    """Unified representation of affinity proteomics data."""
    platform: Literal["olink", "somascan"]
    expression: pd.DataFrame   # samples (rows) x proteins (cols), values are NPX or log10(RFU)
    protein_info: pd.DataFrame  # protein metadata (UniProt, Gene, etc.)
    sample_info: pd.DataFrame   # sample metadata (Group, Age, Sex, etc.)
    qc_summary: dict = field(default_factory=dict)


@dataclass
class DiffAbundanceResult:
    protein_id: str
    gene_symbol: str
    uniprot: str
    mean_group1: float
    mean_group2: float
    log2fc: float
    pvalue: float
    padj: float
    significant: bool


# ---------------------------------------------------------------------------
# Olink NPX parser
# ---------------------------------------------------------------------------
def parse_olink_npx(
    npx_path: Path,
    meta_path: Path | None = None,
    group_col: str = "Group",
) -> ProteomicsData:
    """Parse Olink NPX CSV (long format) into ProteomicsData."""
    npx_df = pd.read_csv(npx_path)

    required = {"SampleID", "OlinkID", "NPX"}
    missing = required - set(npx_df.columns)
    if missing:
        raise ValueError(f"NPX file missing columns: {missing}")

    qc_warned = npx_df[npx_df.get("QC_Warning", pd.Series(dtype=str)).eq("WARN")]
    warned_samples = set(qc_warned["SampleID"].unique()) if len(qc_warned) > 0 else set()

    clean = npx_df[~npx_df["SampleID"].isin(warned_samples)].copy()

    pivot = clean.pivot_table(index="SampleID", columns="OlinkID", values="NPX", aggfunc="first")
    pivot = pivot.dropna(axis=1, how="all")

    protein_cols = list(pivot.columns)
    protein_info_rows = []
    for oid in protein_cols:
        row = clean[clean["OlinkID"] == oid].iloc[0] if len(clean[clean["OlinkID"] == oid]) > 0 else {}
        protein_info_rows.append({
            "protein_id": oid,
            "gene_symbol": row.get("Assay", oid) if isinstance(row, pd.Series) else oid,
            "uniprot": row.get("UniProt", "") if isinstance(row, pd.Series) else "",
            "panel": row.get("Panel", "") if isinstance(row, pd.Series) else "",
        })
    protein_info = pd.DataFrame(protein_info_rows).set_index("protein_id")

    if meta_path and meta_path.exists():
        sample_info = pd.read_csv(meta_path).set_index("SampleID")
    else:
        sample_info = pd.DataFrame(index=pivot.index)

    sample_info = sample_info.reindex(pivot.index)

    lod_summary = {}
    if "LOD" in npx_df.columns:
        for oid in protein_cols:
            sub = clean[clean["OlinkID"] == oid]
            if len(sub) > 0:
                lod_val = sub["LOD"].iloc[0]
                above = (sub["NPX"] > lod_val).sum()
                lod_summary[oid] = {"lod": float(lod_val), "pct_above_lod": round(100 * above / len(sub), 1)}

    return ProteomicsData(
        platform="olink",
        expression=pivot,
        protein_info=protein_info,
        sample_info=sample_info,
        qc_summary={
            "total_samples": len(npx_df["SampleID"].unique()),
            "qc_warned_samples": len(warned_samples),
            "samples_after_qc": len(pivot),
            "total_proteins": len(protein_cols),
            "lod_summary": lod_summary,
        },
    )


# ---------------------------------------------------------------------------
# SomaLogic ADAT parser
# ---------------------------------------------------------------------------
def parse_somascan_adat(adat_path: Path, group_col: str = "SampleGroup") -> ProteomicsData:
    """Parse SomaLogic ADAT file into ProteomicsData."""
    try:
        import somadata
    except ImportError:
        print("ERROR: 'somadata' package required. Install: pip install somadata", file=sys.stderr)
        sys.exit(1)

    adat = somadata.read_adat(str(adat_path))

    sample_meta_cols = list(adat.index.names) if adat.index.names[0] is not None else []
    if sample_meta_cols:
        sample_info = adat.index.to_frame(index=False)
        sample_info.index = range(len(sample_info))
    else:
        sample_info = pd.DataFrame(index=range(len(adat)))

    row_check = sample_info["RowCheck"] if "RowCheck" in sample_info.columns else pd.Series(["PASS"] * len(adat))
    sample_type = sample_info["SampleType"] if "SampleType" in sample_info.columns else pd.Series(["Sample"] * len(adat))

    keep_mask = (row_check == "PASS") & (sample_type == "Sample")
    n_removed_rowcheck = (~(row_check == "PASS") & (sample_type == "Sample")).sum()
    n_removed_type = (sample_type != "Sample").sum()

    rfu_matrix = adat.values[keep_mask].astype(float)
    rfu_matrix = np.where(rfu_matrix > 0, rfu_matrix, np.nan)
    log_matrix = np.log10(rfu_matrix)
    expression = pd.DataFrame(log_matrix, columns=range(rfu_matrix.shape[1]))

    col_meta = {}
    if hasattr(adat.columns, 'codes') or isinstance(adat.columns[0], tuple):
        col_names = adat.columns.get_level_values(0) if hasattr(adat.columns, 'get_level_values') else [c[0] for c in adat.columns]
        for i, col_tuple in enumerate(adat.columns):
            if isinstance(col_tuple, tuple) and len(col_tuple) >= 8:
                col_meta[i] = {
                    "seq_id": col_tuple[0],
                    "gene_symbol": col_tuple[7] if len(col_tuple) > 7 else col_tuple[4],
                    "uniprot": col_tuple[5] if len(col_tuple) > 5 else "",
                    "target": col_tuple[3] if len(col_tuple) > 3 else "",
                    "col_check": col_tuple[15] if len(col_tuple) > 15 else "PASS",
                }

    protein_info_rows = []
    col_check_fail = 0
    valid_cols = []
    for i in range(rfu_matrix.shape[1]):
        meta = col_meta.get(i, {})
        cc = meta.get("col_check", "PASS")
        if cc != "PASS":
            col_check_fail += 1
            continue
        valid_cols.append(i)
        protein_info_rows.append({
            "protein_id": meta.get("seq_id", f"analyte_{i}"),
            "gene_symbol": meta.get("gene_symbol", f"gene_{i}"),
            "uniprot": meta.get("uniprot", ""),
            "target": meta.get("target", ""),
        })

    expression = expression[valid_cols].copy()
    expression.columns = [r["protein_id"] for r in protein_info_rows]
    protein_info = pd.DataFrame(protein_info_rows).set_index("protein_id")

    filtered_sample_info = sample_info[keep_mask].reset_index(drop=True)

    return ProteomicsData(
        platform="somascan",
        expression=expression,
        protein_info=protein_info,
        sample_info=filtered_sample_info,
        qc_summary={
            "total_samples": len(adat),
            "removed_non_sample": int(n_removed_type),
            "removed_rowcheck": int(n_removed_rowcheck),
            "samples_after_qc": len(expression),
            "total_analytes": rfu_matrix.shape[1],
            "col_check_fail": col_check_fail,
            "analytes_after_qc": len(valid_cols),
            "transformation": "log10(RFU)",
        },
    )


# ---------------------------------------------------------------------------
# Differential abundance
# ---------------------------------------------------------------------------
def differential_abundance(
    data: ProteomicsData,
    group_col: str,
    contrast: tuple[str, str],
    fdr_threshold: float = 0.05,
    fc_threshold: float = 0.25,
    test: str = "ttest",
) -> list[DiffAbundanceResult]:
    """Run differential abundance analysis between two groups."""
    from scipy import stats
    from statsmodels.stats.multitest import multipletests

    group1_name, group2_name = contrast
    if group_col not in data.sample_info.columns:
        raise ValueError(f"Group column '{group_col}' not found in sample metadata. Available: {list(data.sample_info.columns)}")

    mask1 = data.sample_info[group_col] == group1_name
    mask2 = data.sample_info[group_col] == group2_name

    if mask1.sum() == 0 or mask2.sum() == 0:
        raise ValueError(f"No samples found for contrast {group1_name} vs {group2_name} in column '{group_col}'")

    expr1 = data.expression.loc[mask1.values]
    expr2 = data.expression.loc[mask2.values]

    results_raw: list[dict] = []
    for protein_id in data.expression.columns:
        v1 = expr1[protein_id].dropna().values
        v2 = expr2[protein_id].dropna().values
        if len(v1) < 3 or len(v2) < 3:
            continue

        mean1, mean2 = float(np.mean(v1)), float(np.mean(v2))
        log2fc = (mean1 - mean2) / math.log10(2) if data.platform == "somascan" else mean1 - mean2

        if test == "ttest":
            stat, pval = stats.ttest_ind(v1, v2, equal_var=False)
        else:
            stat, pval = stats.mannwhitneyu(v1, v2, alternative="two-sided")

        if np.isnan(pval):
            pval = 1.0

        info = data.protein_info.loc[protein_id] if protein_id in data.protein_info.index else {}
        results_raw.append({
            "protein_id": protein_id,
            "gene_symbol": info.get("gene_symbol", protein_id) if isinstance(info, pd.Series) else protein_id,
            "uniprot": info.get("uniprot", "") if isinstance(info, pd.Series) else "",
            "mean_group1": mean1,
            "mean_group2": mean2,
            "log2fc": log2fc,
            "pvalue": pval,
        })

    if not results_raw:
        return []

    pvalues = [r["pvalue"] for r in results_raw]
    reject, padj, _, _ = multipletests(pvalues, method="fdr_bh")

    results: list[DiffAbundanceResult] = []
    for r, adj_p, sig in zip(results_raw, padj, reject):
        results.append(DiffAbundanceResult(
            protein_id=r["protein_id"],
            gene_symbol=r["gene_symbol"],
            uniprot=r["uniprot"],
            mean_group1=r["mean_group1"],
            mean_group2=r["mean_group2"],
            log2fc=r["log2fc"],
            pvalue=r["pvalue"],
            padj=float(adj_p),
            significant=bool(sig) and abs(r["log2fc"]) >= fc_threshold,
        ))

    return sorted(results, key=lambda r: r.pvalue)


# ---------------------------------------------------------------------------
# Visualisation
# ---------------------------------------------------------------------------
def volcano_plot(results: list[DiffAbundanceResult], output_path: Path, contrast: tuple[str, str]) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
    except ImportError:
        return

    log2fcs = [r.log2fc for r in results]
    neg_log10p = [-math.log10(max(r.padj, 1e-300)) for r in results]
    colours = ["#d32f2f" if r.significant else "#9e9e9e" for r in results]

    fig, ax = plt.subplots(figsize=(8, 6))
    ax.scatter(log2fcs, neg_log10p, c=colours, s=15, alpha=0.7, edgecolors="none")

    ax.axhline(y=-math.log10(0.05), color="#ff9800", linestyle="--", linewidth=0.8, label="FDR = 0.05")
    ax.axvline(x=0.25, color="#bdbdbd", linestyle=":", linewidth=0.8)
    ax.axvline(x=-0.25, color="#bdbdbd", linestyle=":", linewidth=0.8)

    sig_results = [r for r in results if r.significant]
    top_sig = sorted(sig_results, key=lambda r: r.pvalue)[:10]
    for r in top_sig:
        ax.annotate(r.gene_symbol, (r.log2fc, -math.log10(max(r.padj, 1e-300))),
                     fontsize=7, ha="center", va="bottom")

    ax.set_xlabel("log2 Fold Change")
    ax.set_ylabel("-log10(adjusted P)")
    ax.set_title(f"Volcano Plot: {contrast[0]} vs {contrast[1]}")
    ax.legend(loc="upper right", fontsize=8)
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def heatmap_plot(
    data: ProteomicsData,
    results: list[DiffAbundanceResult],
    group_col: str,
    output_path: Path,
    top_n: int = 30,
) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        import seaborn as sns
    except ImportError:
        return

    sig_sorted = sorted([r for r in results if r.significant], key=lambda r: r.pvalue)[:top_n]
    if len(sig_sorted) < 2:
        sig_sorted = sorted(results, key=lambda r: r.pvalue)[:min(top_n, len(results))]
    if len(sig_sorted) < 2:
        return

    protein_ids = [r.protein_id for r in sig_sorted]
    labels = [r.gene_symbol for r in sig_sorted]

    sub = data.expression[protein_ids].copy()
    sub = (sub - sub.mean()) / sub.std()
    sub.columns = labels

    if group_col in data.sample_info.columns:
        order = data.sample_info[group_col].sort_values().index
        sub = sub.reindex(order)

    fig_height = max(4, len(labels) * 0.3)
    fig, ax = plt.subplots(figsize=(10, fig_height))
    sns.heatmap(sub.T, cmap="RdBu_r", center=0, ax=ax, xticklabels=False, yticklabels=True)
    ax.set_title(f"Top {len(labels)} Differentially Abundant Proteins (z-scored)")
    ax.set_ylabel("")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


def pca_plot(data: ProteomicsData, group_col: str, output_path: Path) -> None:
    try:
        import matplotlib
        matplotlib.use("Agg")
        import matplotlib.pyplot as plt
        from sklearn.decomposition import PCA
    except ImportError:
        return

    expr_clean = data.expression.dropna(axis=1, how="any").fillna(0)
    if expr_clean.shape[1] < 2:
        return

    pca = PCA(n_components=2)
    coords = pca.fit_transform(expr_clean.values)

    fig, ax = plt.subplots(figsize=(7, 6))
    if group_col in data.sample_info.columns:
        groups = data.sample_info[group_col].values
        unique_groups = sorted(set(groups))
        palette = ["#2166ac", "#b2182b", "#4daf4a", "#984ea3", "#ff7f00"]
        for i, g in enumerate(unique_groups):
            mask = groups == g
            ax.scatter(coords[mask, 0], coords[mask, 1], label=g,
                       color=palette[i % len(palette)], s=30, alpha=0.7)
        ax.legend(fontsize=8)
    else:
        ax.scatter(coords[:, 0], coords[:, 1], s=30, alpha=0.7, color="#2166ac")

    ax.set_xlabel(f"PC1 ({pca.explained_variance_ratio_[0]*100:.1f}%)")
    ax.set_ylabel(f"PC2 ({pca.explained_variance_ratio_[1]*100:.1f}%)")
    ax.set_title("PCA of Protein Expression")
    fig.tight_layout()
    fig.savefig(output_path, dpi=150)
    plt.close(fig)


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------
def generate_report(
    data: ProteomicsData,
    results: list[DiffAbundanceResult],
    contrast: tuple[str, str],
    group_col: str,
    output_dir: Path,
    demo: bool,
    input_label: str,
) -> None:
    timestamp = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    _write_diff_table(results, output_dir / "tables" / "diff_abundance.tsv")
    _write_report_md(data, results, contrast, group_col, output_dir, timestamp, demo, input_label)
    _write_result_json(data, results, contrast, output_dir, timestamp, demo)
    _write_reproducibility(output_dir, timestamp, data.platform, demo)


def _write_diff_table(results: list[DiffAbundanceResult], path: Path) -> None:
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["protein_id", "gene_symbol", "uniprot", "mean_group1", "mean_group2",
                     "log2FC", "pvalue", "padj", "significant"])
        for r in results:
            w.writerow([r.protein_id, r.gene_symbol, r.uniprot,
                        f"{r.mean_group1:.4f}", f"{r.mean_group2:.4f}",
                        f"{r.log2fc:.4f}", f"{r.pvalue:.2e}", f"{r.padj:.2e}", r.significant])


def _write_report_md(data, results, contrast, group_col, output_dir, timestamp, demo, input_label):
    sig = [r for r in results if r.significant]
    lines = [
        "# Affinity Proteomics Report", "",
        f"**Generated**: {timestamp}",
        f"**Platform**: {data.platform.upper()}",
        f"**Input**: {input_label}",
        f"**Contrast**: {contrast[0]} vs {contrast[1]} (column: {group_col})",
        f"**Mode**: {'Demo' if demo else 'Live'}", "",
        "## QC Summary", "",
    ]
    for k, v in data.qc_summary.items():
        if k != "lod_summary":
            lines.append(f"- {k}: {v}")
    lines.append("")

    lines.extend([
        "## Differential Abundance Summary", "",
        f"- Total proteins tested: {len(results)}",
        f"- Significant (FDR < 0.05 and |log2FC| >= 0.25): {len(sig)}",
        f"- Up in {contrast[0]}: {sum(1 for r in sig if r.log2fc > 0)}",
        f"- Up in {contrast[1]}: {sum(1 for r in sig if r.log2fc < 0)}", "",
    ])

    if sig:
        lines.extend(["## Top Significant Proteins", "",
                       "| Protein | Gene | log2FC | P-value | Adj P |",
                       "|---------|------|--------|---------|-------|"])
        for r in sig[:20]:
            lines.append(f"| {r.protein_id} | {r.gene_symbol} | {r.log2fc:.3f} | {r.pvalue:.2e} | {r.padj:.2e} |")
        lines.append("")

    lines.extend(["---", "", f"*{DISCLAIMER}*", ""])
    (output_dir / "report.md").write_text("\n".join(lines), encoding="utf-8")


def _write_result_json(data, results, contrast, output_dir, timestamp, demo):
    sig = [r for r in results if r.significant]
    result = {
        "tool": "ClawBio Affinity Proteomics",
        "version": "0.1.0",
        "platform": data.platform,
        "timestamp": timestamp,
        "mode": "demo" if demo else "live",
        "contrast": list(contrast),
        "qc_summary": {k: v for k, v in data.qc_summary.items() if k != "lod_summary"},
        "total_proteins_tested": len(results),
        "significant_proteins": len(sig),
        "top_10": [
            {"protein_id": r.protein_id, "gene": r.gene_symbol, "log2fc": round(r.log2fc, 4),
             "padj": f"{r.padj:.2e}"}
            for r in results[:10]
        ],
        "disclaimer": DISCLAIMER,
    }
    (output_dir / "result.json").write_text(json.dumps(result, indent=2), encoding="utf-8")


def _write_reproducibility(output_dir, timestamp, platform, demo):
    repro = output_dir / "reproducibility"
    repro.mkdir(exist_ok=True)
    cmd = f"python {Path(__file__).name} --demo --platform {platform} --output {output_dir}" if demo else "# See commands below"
    (repro / "commands.sh").write_text(f"#!/usr/bin/env bash\n# Generated: {timestamp}\n{cmd}\n", encoding="utf-8")
    (repro / "software_versions.json").write_text(json.dumps({
        "python": sys.version, "pandas": pd.__version__, "numpy": np.__version__,
        "generated": timestamp,
    }, indent=2), encoding="utf-8")


# ---------------------------------------------------------------------------
# Demo data helpers
# ---------------------------------------------------------------------------
def _get_somascan_demo_path() -> Path:
    """Locate the bundled SomaLogic example ADAT from the somadata package."""
    try:
        import somadata
        pkg_dir = Path(somadata.__file__).parent
        adat_path = pkg_dir / "data" / "example_data.adat"
        if adat_path.exists():
            return adat_path
    except ImportError:
        pass
    raise FileNotFoundError("somadata package not installed or example_data.adat not found")


# ---------------------------------------------------------------------------
# Main pipeline
# ---------------------------------------------------------------------------
def run_pipeline(
    platform: str,
    input_path: Path,
    meta_path: Path | None,
    group_col: str,
    contrast: tuple[str, str],
    output_dir: Path,
    fdr: float = 0.05,
    fc: float = 0.25,
    top_n: int = 30,
    test: str = "ttest",
    demo: bool = False,
) -> dict:
    """Run the full affinity proteomics pipeline."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "figures").mkdir(exist_ok=True)

    print(f"[PROT] Parsing {platform.upper()} data: {input_path}")
    if platform == "olink":
        data = parse_olink_npx(input_path, meta_path, group_col)
    elif platform == "somascan":
        data = parse_somascan_adat(input_path, group_col)
    else:
        raise ValueError(f"Unknown platform: {platform}")

    print(f"[PROT] Samples after QC: {len(data.expression)}, Proteins: {data.expression.shape[1]}")

    print(f"[PROT] Running differential abundance: {contrast[0]} vs {contrast[1]}")
    results = differential_abundance(data, group_col, contrast, fdr, fc, test)
    sig = [r for r in results if r.significant]
    print(f"[PROT] Significant proteins (FDR<{fdr}, |log2FC|>={fc}): {len(sig)}")

    print("[PROT] Generating plots...")
    volcano_plot(results, output_dir / "figures" / "volcano.png", contrast)
    heatmap_plot(data, results, group_col, output_dir / "figures" / "heatmap.png", top_n)
    pca_plot(data, group_col, output_dir / "figures" / "pca.png")

    print("[PROT] Generating report...")
    generate_report(data, results, contrast, group_col, output_dir, demo, str(input_path))

    print(f"[PROT] Report written to: {output_dir / 'report.md'}")
    return {
        "platform": platform,
        "samples": len(data.expression),
        "proteins": data.expression.shape[1],
        "significant": len(sig),
        "output_dir": str(output_dir),
    }


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------
def main() -> None:
    parser = argparse.ArgumentParser(description="Affinity Proteomics — Olink NPX + SomaLogic SomaScan")
    parser.add_argument("--platform", choices=["olink", "somascan"], help="Platform type")
    parser.add_argument("--input", type=str, help="Input file (NPX CSV for Olink, ADAT for SomaLogic)")
    parser.add_argument("--meta", type=str, help="Sample metadata CSV (required for Olink)")
    parser.add_argument("--group-col", type=str, default="Group", help="Column for group comparison")
    parser.add_argument("--contrast", type=str, default="Case,Control", help="Two groups: Group1,Group2")
    parser.add_argument("--output", type=str, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run with built-in demo data")
    parser.add_argument("--fdr", type=float, default=0.05, help="FDR threshold")
    parser.add_argument("--fc", type=float, default=0.25, help="log2 fold-change threshold")
    parser.add_argument("--top-n", type=int, default=30, help="Top N proteins for heatmap")
    parser.add_argument("--test", choices=["ttest", "wilcoxon"], default="ttest")

    args = parser.parse_args()

    if not args.output:
        parser.error("--output is required")

    contrast = tuple(args.contrast.split(","))
    if len(contrast) != 2:
        parser.error("--contrast must be two comma-separated group names")

    if args.demo:
        platform = args.platform or "olink"
        if platform == "olink":
            input_path = SCRIPT_DIR / "example_data" / "olink_demo_npx.csv"
            meta_path = SCRIPT_DIR / "example_data" / "olink_demo_meta.csv"
            group_col = "Group"
            contrast = ("Case", "Control")
        else:
            input_path = _get_somascan_demo_path()
            meta_path = None
            group_col = "Sex"
            contrast = ("F", "M")
    else:
        if not args.input or not args.platform:
            parser.error("--input and --platform required (or use --demo)")
        platform = args.platform
        input_path = Path(args.input)
        meta_path = Path(args.meta) if args.meta else None
        group_col = args.group_col

    output_dir = Path(args.output)
    print(f"[PROT] Starting {platform.upper()} pipeline ({'demo' if args.demo else 'live'} mode)")

    summary = run_pipeline(
        platform=platform, input_path=input_path, meta_path=meta_path,
        group_col=group_col, contrast=contrast, output_dir=output_dir,
        fdr=args.fdr, fc=args.fc, top_n=args.top_n, test=args.test,
        demo=args.demo,
    )
    print(f"[PROT] Complete: {summary['samples']} samples, {summary['proteins']} proteins, "
          f"{summary['significant']} significant")


if __name__ == "__main__":
    main()
