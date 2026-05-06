"""Tests for the rare-disease-rnaseq skill."""

import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

from rare_disease_rnaseq import (
    DEFAULT_PANEL,
    call_outliers,
    cpm_log2,
    generate_demo,
    per_gene_robust_stats,
)


def test_default_panel_loads_and_has_required_columns():
    panel = pd.read_csv(DEFAULT_PANEL)
    assert {"gene", "mechanism", "phenotype"}.issubset(panel.columns)
    assert len(panel) >= 30
    assert "FBN1" in panel["gene"].values
    assert "NF1" in panel["gene"].values


def test_cpm_log2_normalises_library_size():
    counts = pd.DataFrame({"S1": [10, 20, 30], "S2": [100, 200, 300]},
                          index=["G1", "G2", "G3"])
    expr = cpm_log2(counts)
    # After CPM normalisation S1 and S2 should yield identical per-gene values
    np.testing.assert_allclose(expr["S1"].values, expr["S2"].values)


def test_per_gene_robust_stats_zero_mad_for_constant_gene():
    counts = pd.DataFrame(
        {f"C{i}": [100, 200] for i in range(20)},
        index=["constant_gene", "constant_gene_2"],
    )
    expr = cpm_log2(counts)
    stats = per_gene_robust_stats(expr, list(counts.columns))
    assert (stats["mad"] == 0).all()


def test_zero_mad_genes_skipped_in_outlier_call():
    # Construct expr + stats directly so library-size effects do not
    # confound the unit. constant_gene has MAD=0 (must be skipped);
    # noisy_gene has MAD>0 and the case sits 10 SDs above (must be flagged).
    samples = [f"C{i}" for i in range(20)] + ["CASE-1"]
    expr = pd.DataFrame(
        {
            "constant_gene": [10.0] * 21,
            "noisy_gene": [10.0] * 20 + [15.0],
        },
        index=samples,
    ).T
    stats = pd.DataFrame({
        "median": {"constant_gene": 10.0, "noisy_gene": 10.0},
        "mad": {"constant_gene": 0.0, "noisy_gene": 0.3},
    })
    panel = pd.DataFrame([
        {"gene": "constant_gene", "mechanism": "test", "phenotype": "test"},
        {"gene": "noisy_gene", "mechanism": "test", "phenotype": "test"},
    ])
    outliers = call_outliers(expr, ["CASE-1"], stats, panel, z_threshold=3.0)
    assert "constant_gene" not in outliers["gene"].values
    assert "noisy_gene" in outliers["gene"].values
    assert (outliers[outliers["gene"] == "noisy_gene"]["z_score"] > 0).all()


def test_demo_recovers_injected_outliers(tmp_path):
    counts_path, cases_path, controls_path = generate_demo(seed=42, output_dir=tmp_path)
    counts = pd.read_csv(counts_path, index_col=0)
    case_ids = [s.strip() for s in cases_path.read_text().splitlines() if s.strip()]
    control_ids = [s.strip() for s in controls_path.read_text().splitlines() if s.strip()]

    expr = cpm_log2(counts)
    stats = per_gene_robust_stats(expr, control_ids)
    panel = pd.read_csv(DEFAULT_PANEL)
    outliers = call_outliers(expr, case_ids, stats, panel, z_threshold=3.0)

    # FBN1 must be a strong DOWN outlier in CASE-001
    fbn1 = outliers[(outliers["case"] == "CASE-001") & (outliers["gene"] == "FBN1")]
    assert len(fbn1) == 1
    assert fbn1.iloc[0]["direction"] == "down"
    assert fbn1.iloc[0]["z_score"] < -5

    # NF1 must be a strong UP outlier in CASE-002
    nf1 = outliers[(outliers["case"] == "CASE-002") & (outliers["gene"] == "NF1")]
    assert len(nf1) == 1
    assert nf1.iloc[0]["direction"] == "up"
    assert nf1.iloc[0]["z_score"] > 5

    # Both flagged outliers should be in the disease panel
    panel_hits = outliers[outliers["in_panel"]]
    assert "FBN1" in panel_hits["gene"].values
    assert "NF1" in panel_hits["gene"].values


def test_threshold_controls_call_count(tmp_path):
    counts_path, cases_path, controls_path = generate_demo(seed=7, output_dir=tmp_path)
    counts = pd.read_csv(counts_path, index_col=0)
    case_ids = [s.strip() for s in cases_path.read_text().splitlines() if s.strip()]
    control_ids = [s.strip() for s in controls_path.read_text().splitlines() if s.strip()]
    expr = cpm_log2(counts)
    stats = per_gene_robust_stats(expr, control_ids)
    panel = pd.read_csv(DEFAULT_PANEL)
    loose = call_outliers(expr, case_ids, stats, panel, z_threshold=2.0)
    strict = call_outliers(expr, case_ids, stats, panel, z_threshold=5.0)
    assert len(strict) <= len(loose)
