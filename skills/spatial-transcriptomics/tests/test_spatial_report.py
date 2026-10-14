"""Regression tests for the delivered spatial-analysis report artifacts."""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("scanpy")

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
sys.path.insert(0, str(SKILL_DIR))
sys.path.insert(0, str(PROJECT_ROOT))

import spatial_transcriptomics as st


def _result() -> dict:
    result = st.run_pipeline(st.generate_demo_adata(), random_state=7)
    adata = result["adata"]
    result.update(
        parameters={"n_top_hvg": 40, "random_state": 7},
        analysis_scope={
            "gene_selection": "all",
            "n_genes_tested": len(result["moran"]),
            "genes_tested": [row["gene"] for row in result["moran"]],
            "moran_inference": "descriptive_no_p_values",
        },
        qc_summary={
            "n_spots_loaded": 64,
            "n_spots_after_qc": 64,
            "n_genes_loaded": 40,
            "n_genes_after_qc": 40,
            "mt_gene_count": 0,
            "median_total_counts": 120.0,
            "median_n_genes_by_counts": 30.0,
            "median_pct_counts_mt": 0.0,
            "max_pct_counts_mt": 0.0,
        },
        qc_table=pd.DataFrame(
            {
                "n_genes_by_counts": np.full(adata.n_obs, 30.0),
                "total_counts": np.arange(adata.n_obs) + 100.0,
                "pct_counts_mt": np.zeros(adata.n_obs),
                "pass_qc": np.ones(adata.n_obs, dtype=bool),
            },
            index=adata.obs_names,
        ),
        warnings=["Moran's I is descriptive; p-values were not calculated."],
        count_source="input.X raw counts",
        nhood_permutations=50,
    )
    # Exercise null encoding for undefined spatial statistics.
    result["nhood_zscore"] = result["nhood_zscore"].astype(float)
    result["nhood_zscore"][0, 0] = np.nan
    result["co_occurrence"] = result["co_occurrence"].astype(float)
    result["co_occurrence"][0, 0, 0] = 1.25
    result["co_occurrence"][0, 1, 0] = np.nan
    result["markers"][0]["pvals_adj"] = 0.0123
    return result


def test_main_generate_report_delivers_scores_heatmap_scope_and_qc(tmp_path):
    """The public renderer must not compute these results and then discard them."""
    result = _result()
    checksum_paths, report_path, result_path = st.generate_report(
        result, tmp_path, source_label="synthetic", demo=True
    )
    assert (tmp_path / "tables" / "co_occurrence.csv").is_file()
    assert (tmp_path / "figures" / "marker_heatmap.png").is_file()
    assert (tmp_path / "figures" / "qc_spot_metrics.png").is_file()
    assert (tmp_path / "tables" / "qc_spot_metrics.csv").is_file()
    assert report_path in checksum_paths
    assert result_path in checksum_paths
    report = report_path.read_text(encoding="utf-8")
    assert "40 genes selected as `all` were tested." in report
    assert "Moran's I is descriptive" in report
    assert "Analysis warnings" in report
    assert "mixed cell populations" in report
    assert "P(*j* | *i*, eligible pair" in report
    assert "0.0123" in report
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    assert payload["co_occurrence"]["axis_order"] == [
        "source_cluster",
        "target_cluster",
        "radius",
    ]
    assert payload["co_occurrence"]["radius_lower_exclusive"] == 0


def test_result_json_exposes_cooccurrence_shape_and_rejects_nonfinite_values(tmp_path):
    result = _result()
    _paths, _report_path, result_path = st.generate_report(
        result, tmp_path, source_label="synthetic", demo=True
    )
    raw = result_path.read_text(encoding="utf-8")
    assert "NaN" not in raw
    payload = json.loads(raw)
    cooc = payload["co_occurrence"]
    assert cooc["clusters"] == result["co_occurrence_clusters"]
    assert cooc["radii"] == pytest.approx(result["co_occurrence_distance"].tolist())
    assert (
        np.asarray(cooc["scores"], dtype=object).shape == result["co_occurrence"].shape
    )
    assert cooc["scores"][0][0][0] == pytest.approx(1.25)
    assert cooc["scores"][0][1][0] is None
    csv_text = (tmp_path / "tables" / "co_occurrence.csv").read_text(encoding="utf-8")
    assert "source_cluster,target_cluster,radius,ratio,status" in csv_text
    assert "undefined_no_eligible_pairs" in csv_text


def test_marker_heatmap_contains_rendered_pixels(tmp_path):
    result = _result()
    st.generate_report(result, tmp_path, source_label="synthetic", demo=True)
    from matplotlib import image as mpimg

    pixels = mpimg.imread(tmp_path / "figures" / "marker_heatmap.png")
    assert pixels.ndim == 3
    assert np.unique(pixels.reshape(-1, pixels.shape[-1]), axis=0).shape[0] > 20


def test_qc_figure_marks_mitochondrial_metric_unavailable_when_no_mt_genes(tmp_path):
    result = _result()
    result["qc_summary"]["median_pct_counts_mt"] = None
    result["qc_summary"]["max_pct_counts_mt"] = None
    result["qc_table"]["pct_counts_mt"] = np.nan
    st.generate_report(result, tmp_path, source_label="synthetic", demo=True)
    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "median mitochondrial fraction: NA%" in report
    from matplotlib import image as mpimg

    assert mpimg.imread(tmp_path / "figures" / "qc_spot_metrics.png").size > 0


def test_report_explains_absent_marker_heatmap_without_claiming_an_artifact(tmp_path):
    result = _result()
    # Reuse an output directory to cover --overwrite and prevent stale biology.
    st.generate_report(result, tmp_path, source_label="synthetic", demo=True)
    assert (tmp_path / "figures" / "marker_heatmap.png").exists()
    result["markers"] = []
    st.generate_report(result, tmp_path, source_label="synthetic", demo=True)
    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "No valid reported markers were available for a heatmap." in report
    assert not (tmp_path / "figures" / "marker_heatmap.png").exists()
    assert "| `figures/marker_heatmap.png`" not in report
