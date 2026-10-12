"""Tests for the claw-ancestry-pca skill."""

from __future__ import annotations

import json
import sys
import textwrap
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent.parent
_SKILL_DIR = Path(__file__).resolve().parent.parent

# Add both project root (for clawbio.common) and skill dir (for ancestry_pca)
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
sys.path.insert(0, str(_SKILL_DIR))

from ancestry_pca import (
    compute_pca,
    generate_report,
    load_population_map,
    run_analysis,
    run_summary,
)

# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

DEMO_VCF = _PROJECT_ROOT / "examples" / "demo_populations.vcf"
DEMO_POP_MAP = _PROJECT_ROOT / "examples" / "demo_population_map.csv"


@pytest.fixture
def small_geno_matrix():
    """Deterministic 10-sample, 20-variant genotype matrix."""
    rng = np.random.RandomState(42)
    mat = rng.choice([0, 1, 2], size=(10, 20))
    return mat


@pytest.fixture
def geno_with_missing():
    """Genotype matrix with some missing values (-1)."""
    rng = np.random.RandomState(42)
    mat = rng.choice([0, 1, 2], size=(10, 20))
    # Sprinkle missing
    mat[0, 0] = -1
    mat[3, 5] = -1
    mat[7, 15] = -1
    return mat


@pytest.fixture
def pop_map_csv(tmp_path):
    """Create a temporary population map CSV."""
    csv_path = tmp_path / "pop_map.csv"
    csv_path.write_text(
        "sample_id,population\n"
        "S1,AFR\nS2,AFR\nS3,EUR\nS4,EUR\nS5,EAS\n"
    )
    return csv_path


@pytest.fixture
def pop_map_tsv(tmp_path):
    """Create a temporary population map TSV."""
    tsv_path = tmp_path / "pop_map.tsv"
    tsv_path.write_text(
        "sample_id\tpopulation\n"
        "S1\tAFR\nS2\tAFR\nS3\tEUR\nS4\tEUR\nS5\tEAS\n"
    )
    return tsv_path


# ---------------------------------------------------------------------------
# Population map parsing
# ---------------------------------------------------------------------------


class TestLoadPopulationMap:
    def test_csv_parsing(self, pop_map_csv):
        samples = ["S1", "S2", "S3", "S4", "S5"]
        result = load_population_map(pop_map_csv, samples)
        assert result["S1"] == "AFR"
        assert result["S3"] == "EUR"
        assert result["S5"] == "EAS"
        assert len(result) == 5

    def test_tsv_parsing(self, pop_map_tsv):
        samples = ["S1", "S2", "S3", "S4", "S5"]
        result = load_population_map(pop_map_tsv, samples)
        assert result["S1"] == "AFR"
        assert result["S5"] == "EAS"

    def test_infer_from_prefix(self):
        samples = ["AFR_001", "EUR_002", "EAS_003", "XYZ_004"]
        result = load_population_map(None, samples)
        assert result["AFR_001"] == "AFR"
        assert result["EUR_002"] == "EUR"
        assert result["EAS_003"] == "EAS"
        assert result["XYZ_004"] == "XYZ"

    def test_missing_file_infers(self):
        samples = ["POP_A", "POP_B"]
        result = load_population_map(Path("/nonexistent.csv"), samples)
        # Falls back to prefix inference
        assert result["POP_A"] == "POP"
        assert result["POP_B"] == "POP"

    def test_flexible_column_names(self, tmp_path):
        csv_path = tmp_path / "flex.csv"
        csv_path.write_text(
            "iid,ancestry\nSAM1,AFR\nSAM2,EUR\n"
        )
        result = load_population_map(csv_path, ["SAM1", "SAM2"])
        assert result["SAM1"] == "AFR"
        assert result["SAM2"] == "EUR"


# ---------------------------------------------------------------------------
# PCA computation
# ---------------------------------------------------------------------------


class TestComputePCA:
    def test_basic_shape(self, small_geno_matrix):
        coords, var = compute_pca(small_geno_matrix, n_components=5)
        assert coords.shape == (10, 5)
        assert len(var) == 5

    def test_deterministic(self, small_geno_matrix):
        c1, v1 = compute_pca(small_geno_matrix, n_components=3)
        c2, v2 = compute_pca(small_geno_matrix, n_components=3)
        np.testing.assert_array_almost_equal(np.abs(c1), np.abs(c2))
        np.testing.assert_array_almost_equal(v1, v2)

    def test_variance_sums_to_one_or_less(self, small_geno_matrix):
        _, var = compute_pca(small_geno_matrix, n_components=9)
        assert var.sum() <= 1.0 + 1e-6

    def test_variance_decreasing(self, small_geno_matrix):
        _, var = compute_pca(small_geno_matrix, n_components=5)
        for i in range(len(var) - 1):
            assert var[i] >= var[i + 1] - 1e-10

    def test_missing_data_handling(self, geno_with_missing):
        coords, var = compute_pca(geno_with_missing, n_components=3)
        assert coords.shape[0] == 10
        assert not np.isnan(coords).any()

    def test_component_clamping(self, small_geno_matrix):
        # Request more components than samples-1
        coords, var = compute_pca(small_geno_matrix, n_components=50)
        assert coords.shape[1] <= small_geno_matrix.shape[0] - 1

    def test_all_missing_column(self):
        mat = np.array([[0, -1], [1, -1], [2, -1]])
        coords, var = compute_pca(mat, n_components=1)
        assert coords.shape == (3, 1)
        assert not np.isnan(coords).any()

    def test_single_component(self, small_geno_matrix):
        coords, var = compute_pca(small_geno_matrix, n_components=1)
        assert coords.shape[1] == 1
        assert len(var) == 1


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


class TestReportGeneration:
    def test_report_has_header(self, tmp_path):
        var = np.array([0.5, 0.2, 0.1])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=50,
            n_variants=500,
            n_components=3,
            explained_var=var,
            pop_counts={"AFR": 20, "EUR": 30},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "Ancestry Decomposition PCA" in report
        assert "claw-ancestry-pca" in report

    def test_report_has_disclaimer(self, tmp_path):
        var = np.array([0.5, 0.2])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=2,
            explained_var=var,
            pop_counts={"AFR": 5, "EUR": 5},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "not a medical device" in report

    def test_report_has_variance_table(self, tmp_path):
        var = np.array([0.4, 0.3, 0.1])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=3,
            explained_var=var,
            pop_counts={"AFR": 5, "EUR": 5},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "PC1" in report
        assert "Variance Explained" in report

    def test_report_has_population_table(self, tmp_path):
        var = np.array([0.5])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=1,
            explained_var=var,
            pop_counts={"AFR": 5, "EAS": 3, "EUR": 2},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "AFR" in report
        assert "EAS" in report
        assert "EUR" in report

    def test_report_figure_section_when_generated(self, tmp_path):
        var = np.array([0.5])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=1,
            explained_var=var,
            pop_counts={"AFR": 5, "EUR": 5},
            output_dir=tmp_path,
            figures_generated=True,
        )
        assert "pca_composite.png" in report

    def test_report_references(self, tmp_path):
        var = np.array([0.5])
        report = generate_report(
            input_path=Path("test.vcf"),
            pop_map_path=None,
            n_samples=10,
            n_variants=100,
            n_components=1,
            explained_var=var,
            pop_counts={"EUR": 10},
            output_dir=tmp_path,
            figures_generated=False,
        )
        assert "Mallick" in report
        assert "ClawBio" in report


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


class TestCSVOutput:
    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_coordinates_csv_written(self, tmp_path):
        run_analysis(DEMO_VCF, DEMO_POP_MAP, tmp_path, n_components=3, no_figures=True)
        csv_path = tmp_path / "tables" / "pc_coordinates.csv"
        assert csv_path.exists()
        df = pd.read_csv(csv_path)
        assert "sample_id" in df.columns
        assert "population" in df.columns
        assert "PC1" in df.columns

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_variance_csv_written(self, tmp_path):
        run_analysis(DEMO_VCF, DEMO_POP_MAP, tmp_path, n_components=3, no_figures=True)
        csv_path = tmp_path / "tables" / "variance_explained.csv"
        assert csv_path.exists()
        df = pd.read_csv(csv_path)
        assert "component" in df.columns
        assert "variance_explained" in df.columns
        assert "cumulative" in df.columns


# ---------------------------------------------------------------------------
# Demo mode (end-to-end)
# ---------------------------------------------------------------------------


class TestDemoMode:
    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_full_report(self, tmp_path):
        result = run_analysis(
            DEMO_VCF, DEMO_POP_MAP, tmp_path,
            n_components=5, no_figures=True,
        )
        assert result["n_samples"] == 50
        assert result["n_variants"] > 0
        assert result["n_components"] <= 5
        assert len(result["explained_var"]) > 0
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_result_json(self, tmp_path):
        run_analysis(DEMO_VCF, DEMO_POP_MAP, tmp_path, n_components=3, no_figures=True)
        result_json = json.loads((tmp_path / "result.json").read_text())
        assert result_json["skill"] == "claw-ancestry-pca"
        assert "summary" in result_json
        assert "n_samples" in result_json["summary"]

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_summary_mode(self):
        text = run_summary(DEMO_VCF, DEMO_POP_MAP, n_components=3)
        assert "ANCESTRY DECOMPOSITION PCA" in text
        assert "PC1" in text
        assert "Samples: 50" in text

    @pytest.mark.skipif(
        not DEMO_VCF.exists(), reason="Demo VCF not available"
    )
    def test_demo_populations_detected(self, tmp_path):
        result = run_analysis(
            DEMO_VCF, DEMO_POP_MAP, tmp_path,
            n_components=3, no_figures=True,
        )
        pops = set(result["pop_counts"].keys())
        assert "AFR" in pops
        assert "EUR" in pops


# ---------------------------------------------------------------------------
# Reproducibility bundle
# ---------------------------------------------------------------------------


@pytest.mark.skipif(not DEMO_VCF.exists(), reason="Demo VCF not available")
def test_demo_writes_verifiable_repro_bundle(tmp_path):
    import hashlib
    import subprocess

    out = tmp_path / "out"
    subprocess.run(
        [sys.executable, str(_SKILL_DIR / "ancestry_pca.py"),
         "--demo", "--output", str(out), "--no-figures"],
        check=True, cwd=_PROJECT_ROOT, capture_output=True,
    )
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).exists(), name

    cmd = (repro / "commands.sh").read_text()
    assert "python skills/claw-ancestry-pca/ancestry_pca.py" in cmd
    assert "--demo" in cmd
    assert str(_PROJECT_ROOT) not in cmd

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "result.json", "reproducibility/environment.yml",
            "tables/pc_coordinates.csv"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert hashlib.sha256((out / label).read_bytes()).hexdigest() == digest, label
