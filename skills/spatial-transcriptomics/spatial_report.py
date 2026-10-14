"""Output rendering for the spatial-transcriptomics analysis.

This module deliberately has no Scanpy dependency at import time.  It accepts
the analysis result produced by :mod:`spatial_transcriptomics` and turns its
numeric outputs into inspectable, strictly valid JSON and report artifacts.
"""

from __future__ import annotations

import json
from datetime import UTC, datetime
from pathlib import Path
from typing import Any

import numpy as np

from clawbio.common.report import DISCLAIMER

SKILL_VERSION = "0.1.0"


def _as_json_value(value: Any) -> Any:
    """Convert NumPy/Pandas values to JSON values, representing non-finite as null."""
    if isinstance(value, np.generic):
        value = value.item()
    if isinstance(value, float):
        return value if np.isfinite(value) else None
    if isinstance(value, np.ndarray):
        return _as_json_value(value.tolist())
    if isinstance(value, dict):
        return {str(key): _as_json_value(item) for key, item in value.items()}
    if isinstance(value, (list, tuple)):
        return [_as_json_value(item) for item in value]
    return value


def _display_number(value: Any, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{number:.{digits}f}" if np.isfinite(number) else "NA"


def _display_pvalue(value: Any) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return "NA"
    return f"{number:.3g}" if np.isfinite(number) else "NA"


def _dense_values(adata, genes: list[str]) -> np.ndarray:
    matrix = adata[:, genes].X
    if hasattr(matrix, "toarray"):
        matrix = matrix.toarray()
    return np.asarray(matrix, dtype=float)


def _write_figures(result: dict, output_dir: Path) -> list[Path]:
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    from matplotlib.lines import Line2D

    adata = result["adata"]
    figures = output_dir / "figures"
    figures.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    leiden = adata.obs["leiden"].astype(str)
    clusters = sorted(leiden.unique())
    colour_map = {cluster: plt.cm.tab20(i % 20) for i, cluster in enumerate(clusters)}
    colours = [colour_map[cluster] for cluster in leiden]
    handles = [
        Line2D(
            [0],
            [0],
            marker="o",
            color="w",
            label=cluster,
            markerfacecolor=colour_map[cluster],
            markersize=7,
        )
        for cluster in clusters
    ]

    if "X_umap" in adata.obsm:
        umap = np.asarray(adata.obsm["X_umap"], dtype=float)
        fig, ax = plt.subplots(figsize=(5.5, 4.2))
        ax.scatter(umap[:, 0], umap[:, 1], c=colours, s=28, edgecolors="none")
        ax.set(xlabel="UMAP1", ylabel="UMAP2", title="Leiden clusters (UMAP)")
        ax.legend(
            handles=handles, title="Cluster", bbox_to_anchor=(1.02, 1), loc="upper left"
        )
        fig.tight_layout()
        path = figures / "umap_leiden.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)

    coords = np.asarray(adata.obsm["spatial"], dtype=float)
    fig, ax = plt.subplots(figsize=(5.5, 4.2))
    ax.scatter(coords[:, 0], coords[:, 1], c=colours, s=36, edgecolors="none")
    ax.set(xlabel="spatial x", ylabel="spatial y", title="Leiden clusters (spatial)")
    ax.set_aspect("equal", adjustable="box")
    ax.invert_yaxis()
    ax.legend(
        handles=handles, title="Cluster", bbox_to_anchor=(1.02, 1), loc="upper left"
    )
    fig.tight_layout()
    path = figures / "spatial_leiden.png"
    fig.savefig(path, dpi=140)
    plt.close(fig)
    written.append(path)

    marker_rows = result.get("markers", [])
    # Include an equal number of leading marker rows from every cluster before
    # adding any remaining markers, so early cluster labels do not monopolise
    # a compact heatmap.
    rows_by_cluster: dict[str, list[dict]] = {}
    for row in marker_rows:
        rows_by_cluster.setdefault(str(row.get("cluster", "")), []).append(row)
    marker_genes: list[str] = []
    for marker_rank in range(2):
        for cluster in clusters:
            for row in rows_by_cluster.get(cluster, []):
                if int(row.get("rank", marker_rank + 1)) != marker_rank + 1:
                    continue
                gene = str(row.get("gene", ""))
                if gene and gene in adata.var_names and gene not in marker_genes:
                    marker_genes.append(gene)
                break
    for row in marker_rows:
        gene = str(row.get("gene", ""))
        if gene and gene in adata.var_names and gene not in marker_genes:
            marker_genes.append(gene)
    marker_genes = marker_genes[:24]
    if marker_genes:
        expression = _dense_values(adata, marker_genes)
        means = np.vstack(
            [
                expression[np.asarray(leiden == cluster), :].mean(axis=0)
                for cluster in clusters
            ]
        )
        fig_width = max(5.5, 0.45 * len(marker_genes) + 2.5)
        fig, ax = plt.subplots(
            figsize=(fig_width, max(3.0, 0.55 * len(clusters) + 1.8))
        )
        image = ax.imshow(means, aspect="auto", cmap="viridis")
        ax.set_xticks(range(len(marker_genes)), marker_genes, rotation=55, ha="right")
        ax.set_yticks(range(len(clusters)), clusters)
        ax.set_xlabel("Marker gene")
        ax.set_ylabel("Leiden cluster")
        ax.set_title("Mean log1p expression of reported markers")
        fig.colorbar(image, ax=ax, label="mean log1p expression")
        fig.tight_layout()
        path = figures / "marker_heatmap.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)
    else:
        # A permitted --overwrite run must not embed a heatmap from a prior
        # analysis when this run has no reportable markers.
        (figures / "marker_heatmap.png").unlink(missing_ok=True)

    qc_table = result.get("qc_table")
    if qc_table is not None and len(qc_table):
        fig, axes = plt.subplots(1, 2, figsize=(8, 3.4))
        total_counts = np.asarray(qc_table["total_counts"], dtype=float)
        total_counts = total_counts[np.isfinite(total_counts)]
        axes[0].hist(total_counts, bins=30, color="#4c78a8")
        axes[0].set(title="Total counts per spot", xlabel="counts", ylabel="spots")
        mitochondrial = np.asarray(qc_table["pct_counts_mt"], dtype=float)
        mitochondrial = mitochondrial[np.isfinite(mitochondrial)]
        if mitochondrial.size:
            axes[1].hist(mitochondrial, bins=30, color="#f58518")
        else:
            axes[1].text(
                0.5,
                0.5,
                "Unavailable:\nno MT- gene symbols",
                ha="center",
                va="center",
                transform=axes[1].transAxes,
            )
        axes[1].set(
            title="Mitochondrial fraction",
            xlabel="percent mitochondrial counts",
            ylabel="spots",
        )
        fig.tight_layout()
        path = figures / "qc_spot_metrics.png"
        fig.savefig(path, dpi=140)
        plt.close(fig)
        written.append(path)
    else:
        (figures / "qc_spot_metrics.png").unlink(missing_ok=True)
    return written


def _write_tables(result: dict, output_dir: Path) -> list[Path]:
    import pandas as pd

    tables = output_dir / "tables"
    tables.mkdir(parents=True, exist_ok=True)
    written: list[Path] = []
    for name, rows in (
        ("markers_top.csv", result.get("markers", [])),
        ("moran_i.csv", result.get("moran", [])),
    ):
        path = tables / name
        pd.DataFrame(rows).replace([np.inf, -np.inf], np.nan).to_csv(path, index=False)
        written.append(path)

    clusters = [str(cluster) for cluster in result.get("nhood_clusters", [])]
    zscore = np.asarray(result.get("nhood_zscore", []), dtype=float)
    path = tables / "nhood_enrichment.csv"
    pd.DataFrame(zscore, index=clusters, columns=clusters).replace(
        [np.inf, -np.inf], np.nan
    ).to_csv(path)
    written.append(path)

    radii = np.asarray(result.get("co_occurrence_distance", []), dtype=float)
    cooc_clusters = [
        str(cluster) for cluster in result.get("co_occurrence_clusters", [])
    ]
    scores = np.asarray(result.get("co_occurrence", []), dtype=float)
    cooc_rows: list[dict[str, Any]] = []
    for src_i, source in enumerate(cooc_clusters):
        for dst_i, target in enumerate(cooc_clusters):
            for radius_i, radius in enumerate(radii):
                score = scores[src_i, dst_i, radius_i]
                cooc_rows.append(
                    {
                        "source_cluster": source,
                        "target_cluster": target,
                        "radius": radius,
                        "ratio": score if np.isfinite(score) else np.nan,
                        "status": "defined"
                        if np.isfinite(score)
                        else "undefined_no_eligible_pairs",
                    }
                )
    path = tables / "co_occurrence.csv"
    pd.DataFrame(cooc_rows).to_csv(path, index=False)
    written.append(path)

    qc_table = result.get("qc_table")
    if qc_table is not None:
        path = tables / "qc_spot_metrics.csv"
        qc_table.replace([np.inf, -np.inf], np.nan).to_csv(
            path, index=True, index_label="barcode"
        )
        written.append(path)
    else:
        (tables / "qc_spot_metrics.csv").unlink(missing_ok=True)
    qc_summary = result.get("qc_summary")
    if qc_summary:
        path = tables / "qc_summary.csv"
        pd.DataFrame([qc_summary]).replace([np.inf, -np.inf], np.nan).to_csv(
            path, index=False
        )
        written.append(path)
    else:
        (tables / "qc_summary.csv").unlink(missing_ok=True)
    return written


def generate_report(
    result: dict, output_dir: Path, *, source_label: str, demo: bool
) -> tuple[list[Path], Path, Path]:
    """Write report artifacts and return every file that belongs in checksums."""
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    adata = result["adata"]
    figures = _write_figures(result, output_dir)
    tables = _write_tables(result, output_dir)
    timestamp = datetime.now(UTC).strftime("%Y-%m-%d %H:%M UTC")
    clusters = [str(cluster) for cluster in result.get("nhood_clusters", [])]
    zscore = np.asarray(result.get("nhood_zscore", []), dtype=float)
    top_moran = result.get("moran", [])[:8]
    scope = result.get("analysis_scope", {})
    qc = result.get("qc_summary", {})
    warnings = [str(item) for item in result.get("warnings", [])]
    marker_path = (
        "figures/marker_heatmap.png"
        if (output_dir / "figures/marker_heatmap.png").exists()
        else None
    )
    demo_note = (
        "> Demo mode uses synthetic spots; it is not measured tissue data.\n"
        if demo
        else ""
    )
    scope_text = (
        f"{scope.get('n_genes_tested', len(result.get('moran', [])))} genes selected as "
        f"`{scope.get('gene_selection', 'unspecified')}` were tested."
    )
    lines = [
        "# Spatial Transcriptomics Report" + (" (demo)" if demo else ""),
        "",
        f"**Generated**: {timestamp}  ",
        f"**Skill**: spatial-transcriptomics v{SKILL_VERSION}  ",
        f"**Input**: {source_label}  ",
        f"**Count source**: {result.get('count_source', 'unspecified')}  ",
        f"**Spots**: {adata.n_obs} (from {result.get('n_spots_in', qc.get('n_spots_loaded', 'NA'))} loaded)  ",
        f"**Genes after QC**: {adata.n_vars}  ",
        f"**Leiden clusters**: {adata.obs['leiden'].nunique()}  ",
        "",
        demo_note,
        "## Summary",
        "",
        "Expression-based clustering uses Scanpy. Coordinate-based spatial statistics use the spot coordinates in `obsm['spatial']`.",
        "",
        "![Leiden UMAP](figures/umap_leiden.png)",
        "",
        "![Spatial Leiden clusters](figures/spatial_leiden.png)",
        "",
        "## QC",
        "",
        f"Loaded {qc.get('n_spots_loaded', result.get('n_spots_in', 'NA'))} spots and retained {qc.get('n_spots_after_qc', adata.n_obs)} after QC.",
        f"Median total counts: {_display_number(qc.get('median_total_counts'))}; median detected genes: {_display_number(qc.get('median_n_genes_by_counts'))}; median mitochondrial fraction: {_display_number(qc.get('median_pct_counts_mt'))}%.",
        "",
        "![QC distributions](figures/qc_spot_metrics.png)"
        if (output_dir / "figures/qc_spot_metrics.png").exists()
        else "",
        "",
        "## Spatially variable genes (descriptive Moran's I)",
        "",
        scope_text,
        "Moran's I is descriptive in this run; no p-values or multiple-testing correction were computed.",
        "",
        "| Gene | Moran's I |",
        "|---|---|",
    ]
    lines += [
        f"| {row.get('gene', 'NA')} | {_display_number(row.get('moran_i'))} |"
        for row in top_moran
    ]
    lines += [
        "",
        "## Cluster markers",
        "",
        "![Marker heatmap](figures/marker_heatmap.png)"
        if marker_path
        else "No valid reported markers were available for a heatmap.",
    ]
    lines += ["", "| Cluster | Gene | Score | adj. p |", "|---|---|---|---|"]
    lines += [
        f"| {row.get('cluster', 'NA')} | {row.get('gene', 'NA')} | {_display_number(row.get('score'), 2)} | {_display_pvalue(row.get('pvals_adj'))} |"
        for row in result.get("markers", [])
    ]
    lines += [
        "",
        "## Neighbourhood enrichment",
        "",
        "| From \\ To | " + " | ".join(clusters) + " |",
        "|" + "|".join(["---"] * (len(clusters) + 1)) + "|",
    ]
    for i, source in enumerate(clusters):
        cells = " | ".join(
            _display_number(zscore[i, j], 2) for j in range(len(clusters))
        )
        lines.append(f"| {source} | {cells} |")
    lines += [
        "",
        f"The null distribution uses {result.get('parameters', {}).get('nhood_perms', 'the configured number of')} label permutations. These z-scores do not provide calibrated extreme-tail p-values. `NA` means the permutation null had zero standard deviation, so a z-score is undefined.",
        "",
        "## Cluster co-occurrence",
        "",
        "[Download the complete co-occurrence table](tables/co_occurrence.csv). For source cluster *i*, target cluster *j*, and radius *r*, the ratio is P(*j* | *i*, eligible pair, 0 < distance ≤ *r*) / P(*j* | eligible pair, 0 < distance ≤ *r*). Values above 1 indicate enrichment. `NA` denotes no eligible pairs.",
        "",
        "## Interpretation limits",
        "",
        "Spots can contain mixed cell populations. These outputs describe spot-level expression and spatial patterns; they do not establish cell types, tissue states, or clinical conclusions.",
    ]
    if warnings:
        lines += ["", "## Analysis warnings", ""] + [
            f"- {warning}" for warning in warnings
        ]
    lines += [
        "",
        "## Output files",
        "",
        "| File | Description |",
        "|---|---|",
        "| `result.json` | Machine-readable results, including full co-occurrence scores |",
    ]
    if marker_path:
        lines.append(
            "| `figures/marker_heatmap.png` | Mean log1p expression for reported markers |"
        )
    lines += [
        "| `tables/co_occurrence.csv` | Cumulative-radius co-occurrence ratios |",
        "| `tables/qc_spot_metrics.csv` | Spot-level QC values |",
        "",
        "---",
        "",
        f"*{DISCLAIMER}*",
        "",
    ]
    report_path = output_dir / "report.md"
    report_path.write_text(
        "\n".join(line for line in lines if line is not None), encoding="utf-8"
    )

    payload = _as_json_value(
        {
            "skill": "spatial-transcriptomics",
            "skill_version": SKILL_VERSION,
            "demo": demo,
            "input": source_label,
            "timestamp": timestamp,
            "n_spots": int(adata.n_obs),
            "n_genes": int(adata.n_vars),
            "n_clusters": int(adata.obs["leiden"].nunique()),
            "leiden_clusters": sorted(adata.obs["leiden"].astype(str).unique()),
            "parameters": result.get("parameters", {}),
            "effective_embedding": result.get("effective_embedding", {}),
            "analysis_scope": scope,
            "qc_summary": qc,
            "count_source": result.get("count_source"),
            "input_provenance": result.get("input_provenance", {}),
            "warnings": warnings,
            "top_moran": top_moran,
            "markers": result.get("markers", []),
            "nhood_clusters": clusters,
            "nhood_zscore": zscore,
            "co_occurrence": {
                "axis_order": ["source_cluster", "target_cluster", "radius"],
                "clusters": result.get("co_occurrence_clusters", []),
                "radii": result.get("co_occurrence_distance", []),
                "radius_lower_exclusive": 0,
                "scores": result.get("co_occurrence", []),
            },
        }
    )
    result_path = output_dir / "result.json"
    result_path.write_text(
        json.dumps(payload, indent=2, allow_nan=False) + "\n", encoding="utf-8"
    )
    return figures + tables + [report_path, result_path], report_path, result_path
