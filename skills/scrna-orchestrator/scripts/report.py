#!/usr/bin/env python3
"""Step 7/7 of scrna_orchestrator: figures, tables, report.md, result.json.

Faithful port of `plot_core_figures()`, `write_tables()`, and `render_report()`
from the 1.x skill, trimmed of the self-written `reproducibility/` bundle
(commands.sh/environment.yml/checksums.sha256) — that's now `clawbio bundle`'s
job at the run level, not a per-skill concern.

`prereqs: [find_markers, annotate_celltypist]` in pipeline.yaml — ordering
only for `annotate_celltypist` (optional). Like `preprocess.py`, this script
discovers optional upstream results (doublet summary, CellTypist annotation)
by checking for their conventional filenames in the shared `{out}` directory
rather than gating on their emits, so the report renders correctly whether
zero, one, or both optional steps ran.
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path
from typing import Any

sys.path.insert(0, str(Path(__file__).resolve().parent))
from _common import (  # noqa: E402
    DISCLAIMER,
    generate_report_footer,
    generate_report_header,
    import_scanpy,
    write_result_json,
)

import pandas as pd  # noqa: E402


def plot_core_figures(graph_adata, marker_adata, markers_top: pd.DataFrame, figures_dir: Path) -> list[Path]:
    """Create QC/UMAP/marker plots. Verbatim from the 1.x skill."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sc = import_scanpy()
    figures_dir.mkdir(parents=True, exist_ok=True)
    created: list[Path] = []

    sc.pl.violin(
        graph_adata, ["n_genes_by_counts", "total_counts", "pct_counts_mt"], jitter=0.4, multi_panel=True, show=False
    )
    plt.tight_layout()
    qc_path = figures_dir / "qc_violin.png"
    plt.savefig(qc_path, dpi=180, bbox_inches="tight")
    plt.close("all")
    created.append(qc_path)

    sc.pl.umap(graph_adata, color="leiden", legend_loc="on data", show=False)
    plt.tight_layout()
    umap_path = figures_dir / "umap_leiden.png"
    plt.savefig(umap_path, dpi=180, bbox_inches="tight")
    plt.close("all")
    created.append(umap_path)

    marker_genes = markers_top.groupby("cluster", as_index=False).head(3)["names"].dropna().astype(str).tolist()
    marker_genes = list(dict.fromkeys(marker_genes))[:20]
    if marker_genes:
        dot = sc.pl.dotplot(marker_adata, var_names=marker_genes, groupby="leiden", show=False, return_fig=True)
        marker_path = figures_dir / "marker_dotplot.png"
        dot.savefig(marker_path, dpi=180)
        plt.close("all")
        created.append(marker_path)

    return created


def render_report(
    *,
    output_dir: Path,
    qc_stats: dict[str, Any],
    n_hvg: int,
    n_clusters: int,
    n_pcs_eff: int,
    params: dict[str, Any],
    table_paths: dict[str, Path],
    figure_paths: list[Path],
    n_cells_analyzed: int,
    doublet_summary: dict[str, Any] | None,
    annotation_table: pd.DataFrame | None,
    annotation_info: dict[str, Any] | None,
) -> Path:
    """Create markdown report.md. Trimmed port of render_report() (contrastive-marker
    sections and the self-written reproducibility bundle dropped; see module docstring)."""
    header = generate_report_header(
        title="scRNA Orchestrator Report",
        skill_name="scrna_orchestrator",
        skill_version="2.0.0",
        extra_metadata={
            "Cells (before QC)": str(qc_stats["n_cells_before"]),
            "Cells (after QC)": str(n_cells_analyzed),
            "Genes (after QC)": str(qc_stats["n_genes_after"]),
            "Leiden clusters": str(n_clusters),
            "HVG selected": str(n_hvg),
            "Graph basis": params["graph_basis"],
            "Doublet detection": "scrublet" if doublet_summary is not None else "skipped",
            "Annotation backend": "celltypist" if annotation_info is not None else "skipped",
        },
    )

    figure_labels = {"qc_violin.png": "QC Violin", "umap_leiden.png": "UMAP Leiden", "marker_dotplot.png": "Marker Dotplot"}
    lines = ["## Summary", ""]
    lines.append(f"- Cells before QC: **{qc_stats['n_cells_before']}**")
    lines.append(f"- Cells after QC: **{n_cells_analyzed}**")
    if doublet_summary is not None:
        lines.append(f"- Predicted doublets: **{doublet_summary['n_predicted_doublets']}**")
        lines.append(
            f"- Doublet detection: **scrublet** ({doublet_summary['predicted_doublet_rate']:.1%} predicted)"
        )
    else:
        lines.append("- Doublet detection: **skipped** (`clawbio skip detect_doublets`)")
    lines.append(f"- Genes before QC: **{qc_stats['n_genes_before']}**")
    lines.append(f"- Genes after QC: **{qc_stats['n_genes_after']}**")
    lines.append(f"- HVGs selected: **{n_hvg}**")
    lines.append(f"- Leiden clusters: **{n_clusters}**")
    if annotation_info is not None:
        lines.append(f"- Clusters annotated: **{annotation_info['n_clusters_annotated']}**")

    lines.extend(["", "## Core Figures", ""])
    for figure_path in figure_paths:
        if figure_path.name not in figure_labels:
            continue
        lines.append(f"![{figure_labels[figure_path.name]}](figures/{figure_path.name})")

    lines.extend(["", "## Tables", ""])
    for table_path in table_paths.values():
        lines.append(f"- `tables/{table_path.name}`")

    if doublet_summary is not None:
        lines.extend(["", "## Doublet Detection", ""])
        lines.append("- Method: `scanpy.pp.scrublet`")
        lines.append(f"- Cells scored: **{doublet_summary['n_cells_scored']}**")
        lines.append(f"- Predicted doublets: **{doublet_summary['n_predicted_doublets']}**")
        lines.append(f"- Cells retained: **{doublet_summary['n_cells_retained']}**")
        lines.append("- Summary table: `tables/doublet_summary.csv`")
    else:
        lines.extend(["", "## Doublet Detection", ""])
        lines.append("- Not enabled for this run (`detect_doublets` was skipped).")

    if annotation_table is not None and annotation_info is not None:
        lines.extend(["", "## Cell Type Annotation", ""])
        lines.append("- Backend: `CellTypist`")
        lines.append(f"- Model: `{annotation_info['model']}`")
        lines.append(f"- Overlapping genes with model: **{annotation_info['overlap_genes']}**")
        lines.append("- Labels are **putative**, model-based assignments and should be manually reviewed.")
        lines.append("- Summary table: `tables/cluster_annotations.csv`")
        lines.append("")
        for row in annotation_table.itertuples(index=False):
            lines.append(
                f"- Cluster `{row.cluster}` -> **{row.predicted_cell_type}** "
                f"(support={row.support_fraction:.2f}, mean_confidence={row.mean_confidence:.2f})"
            )
    else:
        lines.extend(["", "## Cell Type Annotation", ""])
        lines.append("- Not enabled for this run (`annotate_celltypist` was skipped).")

    lines.extend(["", "## Methods", ""])
    lines.append(f"- QC/filtering: `min_genes={params['min_genes']}`, `max_mt_pct={params['max_mt_pct']}`")
    if doublet_summary is not None:
        lines.append("- Doublet detection: `scanpy.pp.scrublet` on QC-filtered raw counts before normalization/clustering")
    lines.append("- Normalisation: total-count normalisation (`target_sum=1e4`) + `log1p`")
    lines.append(f"- Feature selection: `n_top_hvg={params['n_top_hvg']}`")
    if params["graph_basis"] == "pca":
        lines.append(f"- Graph construction: PCA (`n_pcs={n_pcs_eff}`) + neighbors + UMAP")
    else:
        lines.append(f"- Graph construction: latent representation `{params['graph_basis']}` + neighbors + UMAP")
    lines.append(f"- Clustering: Leiden `resolution={params['leiden_resolution']}`")
    lines.append("- Marker analysis: `scanpy.tl.rank_genes_groups` (Wilcoxon, cluster-vs-rest) on normalized full-gene expression")
    if annotation_info is not None:
        lines.append(f"- Annotation: `CellTypist` model `{annotation_info['model']}` on normalized/log1p full-gene expression")

    lines.extend(["", "## Reproducibility", "", "See `clawbio bundle` output: `commands.sh`, `environment.yml`, `checksums.sha256`."])

    report_path = output_dir / "report.md"
    report_path.write_text(header + "\n".join(lines) + generate_report_footer(), encoding="utf-8")
    return report_path


def main() -> None:
    parser = argparse.ArgumentParser(description="scrna_orchestrator: report step")
    parser.add_argument("--clustered-h5ad", required=True)
    parser.add_argument("--markers-h5ad", required=True)
    parser.add_argument("--markers-top-csv", required=True)
    parser.add_argument("--cluster-summary-csv", required=True)
    parser.add_argument("--qc-stats-json", required=True)
    parser.add_argument("--cluster-summary-json", required=True)
    parser.add_argument("--preprocess-summary-json", required=True)
    parser.add_argument("--min-genes", type=int, default=200)
    parser.add_argument("--max-mt-pct", type=float, default=20.0)
    parser.add_argument("--n-top-hvg", type=int, default=2000)
    parser.add_argument("--leiden-resolution", type=float, default=1.0)
    parser.add_argument("--report-out", required=True, help="output path for report.md")
    parser.add_argument("--result-json-out", required=True, help="output path for result.json")
    args = parser.parse_args()

    sc = import_scanpy()
    output_dir = Path(args.report_out).parent
    output_dir.mkdir(parents=True, exist_ok=True)
    figures_dir = output_dir / "figures"
    tables_dir = output_dir / "tables"

    graph_adata = sc.read_h5ad(Path(args.clustered_h5ad))
    marker_adata = sc.read_h5ad(Path(args.markers_h5ad))
    markers_top = pd.read_csv(args.markers_top_csv)
    cluster_summary = pd.read_csv(args.cluster_summary_csv)
    qc_stats = json.loads(Path(args.qc_stats_json).read_text())
    cluster_info = json.loads(Path(args.cluster_summary_json).read_text())
    preprocess_info = json.loads(Path(args.preprocess_summary_json).read_text())

    # Opportunistic sibling lookups for optional upstream steps (see module docstring).
    shared_out_dir = Path(args.clustered_h5ad).parent
    doublet_summary_path = shared_out_dir / "doublet_summary.json"
    doublet_summary = json.loads(doublet_summary_path.read_text()) if doublet_summary_path.is_file() else None

    annotation_csv_path = shared_out_dir / "cluster_annotations.csv"
    annotation_info_path = shared_out_dir / "annotation_info.json"
    if annotation_csv_path.is_file() and annotation_info_path.is_file():
        annotation_table = pd.read_csv(annotation_csv_path)
        annotation_info = json.loads(annotation_info_path.read_text())
    else:
        annotation_table = None
        annotation_info = None

    table_paths: dict[str, Path] = {}
    tables_dir.mkdir(parents=True, exist_ok=True)
    for name, df in (("cluster_summary.csv", cluster_summary), ("markers_top.csv", markers_top)):
        dest = tables_dir / name
        df.to_csv(dest, index=False)
        table_paths[name] = dest
    if doublet_summary is not None:
        dest = tables_dir / "doublet_summary.csv"
        pd.DataFrame([doublet_summary]).to_csv(dest, index=False)
        table_paths["doublet_summary.csv"] = dest
    if annotation_table is not None:
        dest = tables_dir / "cluster_annotations.csv"
        annotation_table.to_csv(dest, index=False)
        table_paths["cluster_annotations.csv"] = dest

    figure_paths = plot_core_figures(graph_adata, marker_adata, markers_top, figures_dir)

    n_clusters = int(graph_adata.obs["leiden"].nunique())
    n_cells_analyzed = int(graph_adata.n_obs)
    params = {
        "min_genes": args.min_genes,
        "max_mt_pct": args.max_mt_pct,
        "n_top_hvg": args.n_top_hvg,
        "leiden_resolution": args.leiden_resolution,
        "graph_basis": cluster_info["graph_basis"],
    }

    report_path = render_report(
        output_dir=output_dir,
        qc_stats=qc_stats,
        n_hvg=preprocess_info["n_hvg"],
        n_clusters=n_clusters,
        n_pcs_eff=cluster_info["n_pcs_eff"],
        params=params,
        table_paths=table_paths,
        figure_paths=figure_paths,
        n_cells_analyzed=n_cells_analyzed,
        doublet_summary=doublet_summary,
        annotation_table=annotation_table,
        annotation_info=annotation_info,
    )

    summary = {
        "n_cells_before": qc_stats["n_cells_before"],
        "n_cells_after": n_cells_analyzed,
        "n_genes_before": qc_stats["n_genes_before"],
        "n_genes_after": qc_stats["n_genes_after"],
        "n_hvg": preprocess_info["n_hvg"],
        "n_clusters": n_clusters,
        "graph_basis": cluster_info["graph_basis"],
        "doublet_detection_ran": doublet_summary is not None,
        "annotation_ran": annotation_info is not None,
    }
    data: dict[str, Any] = {
        "cluster_labels": sorted(graph_adata.obs["leiden"].astype(str).unique().tolist()),
        "tables": [p.name for p in table_paths.values()],
        "figures": [p.name for p in figure_paths],
        "disclaimer": DISCLAIMER,
    }
    if doublet_summary is not None:
        data["doublet"] = doublet_summary
    if annotation_info is not None:
        data["annotation"] = annotation_info

    write_result_json(
        output_dir=output_dir,
        skill="scrna_orchestrator",
        version="2.0.0",
        summary=summary,
        data=data,
    )
    # write_result_json always writes to output_dir/result.json; move/copy to the
    # declared --result-json-out path if the caller asked for a different name.
    result_json_default = output_dir / "result.json"
    result_json_out = Path(args.result_json_out)
    if result_json_default != result_json_out:
        result_json_out.parent.mkdir(parents=True, exist_ok=True)
        result_json_out.write_text(result_json_default.read_text())

    print(f"report complete: {n_clusters} clusters, {n_cells_analyzed} cells")
    print(f"  wrote {report_path}")
    print(f"  wrote {result_json_out}")


if __name__ == "__main__":
    main()
