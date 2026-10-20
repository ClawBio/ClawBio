"""Tests for diff-visualizer skill."""

from __future__ import annotations

import importlib.util
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = SKILL_DIR / "diff_visualizer.py"
EXAMPLES_DIR = SKILL_DIR / "examples"
REPO_ROOT = SKILL_DIR.parent.parent
CLAWBIO_PATH = REPO_ROOT / "clawbio.py"
ORCHESTRATOR_PATH = REPO_ROOT / "skills" / "bio-orchestrator" / "orchestrator.py"
RNASEQ_SCRIPT = REPO_ROOT / "skills" / "rnaseq-de" / "rnaseq_de.py"
RNASEQ_COUNTS = REPO_ROOT / "skills" / "rnaseq-de" / "examples" / "demo_counts.csv"
RNASEQ_META = REPO_ROOT / "skills" / "rnaseq-de" / "examples" / "demo_metadata.csv"


def _run_diffviz(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )


def _run_clawbio(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(CLAWBIO_PATH), "run", "diffviz"] + args,
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )


def _load_orchestrator_module():
    spec = importlib.util.spec_from_file_location("bio_orchestrator_module", ORCHESTRATOR_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _load_diffviz_module():
    spec = importlib.util.spec_from_file_location("diff_visualizer_module", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


def _require_anndata() -> None:
    pytest.importorskip("anndata")


def _write_within_cluster_contrast_table(path: Path) -> None:
    df = pd.DataFrame(
        [
            {
                "cluster": "0",
                "scope": "within-cluster",
                "groupby": "condition",
                "group1": "treated",
                "group2": "control",
                "comparison_id": "treated__vs__control",
                "names": "GeneA",
                "scores": 5.2,
                "logfoldchanges": 1.7,
                "pvals_adj": 0.001,
            },
            {
                "cluster": "0",
                "scope": "within-cluster",
                "groupby": "condition",
                "group1": "treated",
                "group2": "control",
                "comparison_id": "treated__vs__control",
                "names": "GeneB",
                "scores": 4.7,
                "logfoldchanges": 1.2,
                "pvals_adj": 0.004,
            },
            {
                "cluster": "1",
                "scope": "within-cluster",
                "groupby": "condition",
                "group1": "treated",
                "group2": "control",
                "comparison_id": "treated__vs__control",
                "names": "GeneC",
                "scores": 6.1,
                "logfoldchanges": -1.5,
                "pvals_adj": 0.002,
            },
            {
                "cluster": "1",
                "scope": "within-cluster",
                "groupby": "condition",
                "group1": "treated",
                "group2": "control",
                "comparison_id": "treated__vs__control",
                "names": "GeneD",
                "scores": 4.9,
                "logfoldchanges": -1.1,
                "pvals_adj": 0.007,
            },
        ]
    )
    path.parent.mkdir(parents=True, exist_ok=True)
    df.to_csv(path, index=False)


def test_bulk_table_input_outputs(tmp_path: Path):
    output_dir = tmp_path / "bulk_output"
    result = _run_diffviz(
        [
            "--input",
            str(EXAMPLES_DIR / "demo_bulk_de_results.csv"),
            "--output",
            str(output_dir),
        ]
    )
    assert result.returncode == 0, result.stderr

    expected = [
        output_dir / "report.md",
        output_dir / "report.html",
        output_dir / "result.json",
        output_dir / "figures" / "volcano.png",
        output_dir / "figures" / "top_genes_bar.png",
        output_dir / "figures" / "ma_plot.png",
        output_dir / "tables" / "top_genes.csv",
        output_dir / "tables" / "significant_genes.csv",
    ]
    for path in expected:
        assert path.exists(), f"Missing output file: {path}"

    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["mode"] == "bulk"
    assert payload["summary"]["source_kind"] == "bulk-table"
    assert "volcano.png" in payload["summary"]["figures_generated"]


def test_rnaseq_output_dir_autodetected_and_heatmap_skip_is_reported(tmp_path: Path):
    rnaseq_out = tmp_path / "rnaseq_upstream"
    rnaseq_result = subprocess.run(
        [
            sys.executable,
            str(RNASEQ_SCRIPT),
            "--counts",
            str(RNASEQ_COUNTS),
            "--metadata",
            str(RNASEQ_META),
            "--formula",
            "~ batch + condition",
            "--contrast",
            "condition,treated,control",
            "--backend",
            "simple",
            "--output",
            str(rnaseq_out),
        ],
        cwd=str(REPO_ROOT),
        capture_output=True,
        text=True,
    )
    assert rnaseq_result.returncode == 0, rnaseq_result.stderr

    output_dir = tmp_path / "diffviz_from_rnaseq"
    result = _run_diffviz(["--input", str(rnaseq_out), "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr
    assert (output_dir / "figures" / "volcano.png").exists()
    report_text = (output_dir / "report.md").read_text()
    assert "Skipped bulk heatmap because both --counts and --metadata were not provided." in report_text
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["source_kind"] == "rnaseq-output-dir"


def test_scrna_contrast_table_input_outputs(tmp_path: Path):
    output_dir = tmp_path / "scrna_contrast_output"
    result = _run_diffviz(
        [
            "--mode",
            "scrna",
            "--input",
            str(EXAMPLES_DIR / "demo_scrna_contrast.csv"),
            "--output",
            str(output_dir),
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "figures" / "contrast_volcano.png").exists()
    assert (output_dir / "figures" / "top_markers_bar.png").exists()
    assert (output_dir / "tables" / "top_markers.csv").exists()
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["mode"] == "scrna"
    assert payload["summary"]["source_kind"] == "scrna-contrast-table"


def test_scrna_output_dir_with_dataset_contrast_is_autodetected(tmp_path: Path):
    upstream_dir = tmp_path / "scrna_upstream"
    tables_dir = upstream_dir / "tables"
    tables_dir.mkdir(parents=True)
    contrast_path = tables_dir / "contrastive_markers_full.csv"
    contrast_path.write_text((EXAMPLES_DIR / "demo_scrna_contrast.csv").read_text(), encoding="utf-8")

    output_dir = tmp_path / "scrna_from_output_dir"
    result = _run_diffviz(["--mode", "scrna", "--input", str(upstream_dir), "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr
    assert (output_dir / "figures" / "contrast_volcano.png").exists()
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["source_kind"] == "scrna-output-dir"


def test_scrna_output_dir_with_within_cluster_contrast_renders_cluster_panels(tmp_path: Path):
    upstream_dir = tmp_path / "scrna_within_cluster_upstream"
    within_cluster_path = upstream_dir / "tables" / "within_cluster_contrastive_markers_full.csv"
    _write_within_cluster_contrast_table(within_cluster_path)

    output_dir = tmp_path / "scrna_within_cluster_viz"
    result = _run_diffviz(["--mode", "scrna", "--input", str(upstream_dir), "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr
    assert (output_dir / "figures" / "within_cluster_marker_panels.png").exists()
    assert (output_dir / "tables" / "within_cluster_top_markers.csv").exists()
    assert not (output_dir / "figures" / "contrast_volcano.png").exists()
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["source_kind"] == "scrna-output-dir"
    assert payload["summary"]["within_cluster_comparisons"] == 2


@pytest.mark.skipif(sys.version_info < (3, 12), reason="anndata>=0.13 writes pandas 3 strings to .h5ad; it needs Python 3.12+")
def test_scrna_markers_with_adata_generate_enhanced_plots(tmp_path: Path):
    _require_anndata()
    import anndata as ad

    output_dir = tmp_path / "scrna_markers_adata_output"
    adata_path = tmp_path / "tiny_scrna.h5ad"
    X = np.array(
        [
            [5, 4, 1, 0, 0],
            [4, 5, 1, 0, 0],
            [0, 1, 5, 4, 1],
            [0, 0, 4, 5, 1],
            [1, 0, 0, 1, 5],
            [1, 0, 0, 0, 4],
        ],
        dtype=float,
    )
    adata = ad.AnnData(
        X=X,
        obs=pd.DataFrame(
            {"cluster": ["0", "0", "1", "1", "2", "2"]},
            index=[f"cell_{i}" for i in range(X.shape[0])],
        ),
        var=pd.DataFrame(index=["MS4A1", "CD79A", "CD3D", "NKG7", "GNLY"]),
    )
    adata.obsm["X_umap"] = np.array(
        [[0, 0], [0.1, 0.2], [1.0, 0.8], [1.2, 1.1], [2.0, 1.8], [2.1, 2.0]],
        dtype=float,
    )
    adata.write_h5ad(adata_path)

    result = _run_diffviz(
        [
            "--mode",
            "scrna",
            "--input",
            str(EXAMPLES_DIR / "demo_scrna_markers.csv"),
            "--adata",
            str(adata_path),
            "--output",
            str(output_dir),
        ]
    )
    assert result.returncode == 0, result.stderr
    expected = [
        output_dir / "figures" / "marker_rank_bars.png",
        output_dir / "figures" / "marker_dotplot.png",
        output_dir / "figures" / "marker_heatmap.png",
        output_dir / "figures" / "umap_feature_panel.png",
        output_dir / "tables" / "top_markers_by_cluster.csv",
    ]
    for path in expected:
        assert path.exists(), f"Missing output file: {path}"

    payload = json.loads((output_dir / "result.json").read_text())
    assert "adata" in payload["summary"]["enhanced_inputs_used"]


def test_detection_precedence_prefers_finished_tables(tmp_path: Path):
    orchestrator = _load_orchestrator_module()
    diffviz_table = tmp_path / "de_results.csv"
    diffviz_table.write_text(
        "gene,log2FoldChange,padj,baseMean\nGeneA,2.1,0.001,100\nGeneB,-1.2,0.02,90\n",
        encoding="utf-8",
    )
    raw_counts = tmp_path / "counts.csv"
    raw_counts.write_text(
        "gene,s1,s2,s3,s4\nGeneA,10,12,30,32\nGeneB,8,9,18,21\nGeneC,3,4,7,8\n",
        encoding="utf-8",
    )
    raw_meta = tmp_path / "metadata.csv"
    raw_meta.write_text(
        "sample_id,condition\ns1,control\ns2,control\ns3,treated\ns4,treated\n",
        encoding="utf-8",
    )

    assert orchestrator.detect_skill_from_tabular_header(diffviz_table) == "diff-visualizer"
    assert orchestrator.detect_skill_from_tabular_header(raw_counts) == "rnaseq-de"
    assert orchestrator.detect_skill_from_tabular_header(raw_meta) == "rnaseq-de"


def test_non_empty_output_dir_is_rejected(tmp_path: Path):
    output_dir = tmp_path / "occupied"
    output_dir.mkdir()
    (output_dir / "keep.txt").write_text("occupied", encoding="utf-8")

    result = _run_diffviz(
        [
            "--input",
            str(EXAMPLES_DIR / "demo_bulk_de_results.csv"),
            "--output",
            str(output_dir),
        ]
    )
    assert result.returncode != 0
    assert "not empty" in result.stderr


def test_missing_required_columns_fail_with_actionable_error(tmp_path: Path):
    bad_input = tmp_path / "bad.csv"
    bad_input.write_text("gene,score\nGeneA,1.2\n", encoding="utf-8")

    result = _run_diffviz(["--input", str(bad_input), "--output", str(tmp_path / "out")])
    assert result.returncode != 0
    assert "Could not auto-detect input table type" in result.stderr


def test_clawbio_runner_demo_modes_pass(tmp_path: Path):
    bulk_out = tmp_path / "bulk_demo"
    bulk_result = _run_clawbio(["--demo", "--output", str(bulk_out)])
    assert bulk_result.returncode == 0, bulk_result.stderr
    assert (bulk_out / "report.md").exists()
    assert (bulk_out / "report.html").exists()

    scrna_out = tmp_path / "scrna_demo"
    scrna_result = _run_clawbio(["--demo", "--mode", "scrna", "--output", str(scrna_out)])
    assert scrna_result.returncode == 0, scrna_result.stderr
    assert (scrna_out / "figures" / "contrast_volcano.png").exists()
    assert (scrna_out / "figures" / "marker_rank_bars.png").exists()


def test_bulk_volcano_display_cap_handles_extreme_outliers(tmp_path: Path):
    diffviz = _load_diffviz_module()
    df = pd.DataFrame(
        {
            "gene": [f"Gene{i}" for i in range(12)],
            "log2FoldChange": [0.8, -1.2, 1.6, -2.0, 2.4, -3.1, 4.2, -4.8, 5.5, -6.0, 7.1, 43.0],
            "padj": [0.6, 0.04, 0.01, 0.2, 0.03, 0.001, 0.02, 0.15, 0.0005, 0.9, 0.04, 1e-6],
            "baseMean": [120, 90, 150, 140, 180, 210, 240, 80, 300, 60, 260, 22],
        }
    )
    outpath = tmp_path / "volcano.png"
    diffviz.plot_bulk_volcano(
        df,
        outpath,
        padj_threshold=0.05,
        lfc_threshold=1.0,
        label_top=5,
    )
    assert outpath.exists()
    display_cap = diffviz._display_cap(df["log2FoldChange"], minimum=4.0)
    assert display_cap < 43.0


def test_bulk_visualization_filter_respects_min_basemean(tmp_path: Path):
    input_path = tmp_path / "bulk.csv"
    input_path.write_text(
        "\n".join(
            [
                "gene,log2FoldChange,padj,baseMean",
                "LowExprHit,7.5,1e-12,5",
                "HighExprHit,3.2,1e-10,48",
                "HighExprDown,-2.8,1e-6,65",
                "Borderline,1.5,0.03,11",
            ]
        )
        + "\n",
        encoding="utf-8",
    )
    output_dir = tmp_path / "filtered_output"
    result = _run_diffviz(["--input", str(input_path), "--output", str(output_dir), "--min-basemean", "10"])
    assert result.returncode == 0, result.stderr

    top_genes = pd.read_csv(output_dir / "tables" / "top_genes.csv")
    assert "LowExprHit" not in set(top_genes["gene"])
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["bulk_display_filter"]["applied"] is True
    assert payload["summary"]["bulk_display_filter"]["rows_after"] == 3


def test_orchestrator_routes_visualization_queries_without_breaking_analysis_routes():
    orchestrator = _load_orchestrator_module()
    assert orchestrator.detect_skill_from_query("visualize DE results from this table") == "diff-visualizer"
    assert orchestrator.detect_skill_from_query("make a marker heatmap for these markers") == "diff-visualizer"
    assert orchestrator.detect_skill_from_query("run differential expression on this count matrix") == "rnaseq-de"
    assert orchestrator.detect_skill_from_query("find marker genes in this h5ad single-cell dataset") == "scrna-orchestrator"
