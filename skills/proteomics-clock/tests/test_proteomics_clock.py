import json
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import numpy as np
import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).parent.parent))

SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parents[1]
DEMO_FIXTURE = SKILL_DIR / "data" / "demo_olink_npx.csv.gz"


# ---------------------------------------------------------------------------
# Unit tests
# ---------------------------------------------------------------------------


def test_parse_organ_list_defaults():
    from proteomics_clock import DEFAULT_ORGANS, parse_organ_list

    organs = parse_organ_list(None)
    assert organs == DEFAULT_ORGANS
    assert "Bladder" not in organs
    assert "Heart" in organs
    assert "Brain" in organs


def test_parse_organ_list_custom():
    from proteomics_clock import parse_organ_list

    organs = parse_organ_list("Heart,Brain,Kidney")
    assert organs == ["Heart", "Brain", "Kidney"]


def test_parse_organ_list_empty_raises():
    from proteomics_clock import parse_organ_list

    with pytest.raises(ValueError):
        parse_organ_list("")


def test_mortality_to_years():
    from proteomics_clock import mortality_to_years

    hazards = np.array([0.0, -4.801912])
    years = mortality_to_years(hazards)
    assert len(years) == 2
    assert all(np.isfinite(years))
    # A hazard of 0 should give a reasonable age (positive, < 200)
    assert 0 < years[0] < 200
    # The average hazard (-4.802) maps to ~9.95 per the Gompertz formula
    # (this is the baseline offset, not the population mean age)
    assert 5 < years[1] < 15


def test_predict_gen1_simple():
    from proteomics_clock import predict_organ_ages

    df = pd.DataFrame({
        "sample_id": ["S1", "S2"],
        "ProteinA": [1.0, 2.0],
        "ProteinB": [3.0, 4.0],
    })
    coefficients = {
        "TestOrgan": {"Intercept": 10.0, "ProteinA": 2.0, "ProteinB": 0.5},
    }
    preds, missing = predict_organ_ages(df, coefficients, ["TestOrgan"], "gen1")
    assert "TestOrgan" in preds.columns
    # S1: 10 + 2*1 + 0.5*3 = 13.5
    assert abs(preds.loc[0, "TestOrgan"] - 13.5) < 1e-6
    # S2: 10 + 2*2 + 0.5*4 = 16.0
    assert abs(preds.loc[1, "TestOrgan"] - 16.0) < 1e-6
    assert missing.empty or missing[missing["organ"] == "TestOrgan"].empty


def test_predict_gen2_no_intercept():
    from proteomics_clock import predict_organ_ages

    df = pd.DataFrame({
        "sample_id": ["S1"],
        "ProteinA": [1.0],
        "ProteinB": [3.0],
    })
    # Gen2 coefficients should NOT have Intercept used
    coefficients = {
        "TestOrgan": {"ProteinA": 2.0, "ProteinB": 0.5},
    }
    preds, missing = predict_organ_ages(df, coefficients, ["TestOrgan"], "gen2")
    # S1: 2*1 + 0.5*3 = 3.5
    assert abs(preds.loc[0, "TestOrgan"] - 3.5) < 1e-6


def test_missing_proteins_reported():
    from proteomics_clock import predict_organ_ages

    df = pd.DataFrame({
        "sample_id": ["S1"],
        "ProteinA": [1.0],
        # ProteinB and ProteinC are missing from input
    })
    coefficients = {
        "TestOrgan": {"Intercept": 5.0, "ProteinA": 1.0, "ProteinB": 2.0, "ProteinC": 3.0},
    }
    preds, missing = predict_organ_ages(df, coefficients, ["TestOrgan"], "gen1")
    assert not missing.empty
    missing_for_organ = missing[missing["organ"] == "TestOrgan"]
    missing_proteins = set(missing_for_organ["protein"].tolist())
    assert "ProteinB" in missing_proteins
    assert "ProteinC" in missing_proteins
    assert "ProteinA" not in missing_proteins


def test_load_input_csv(tmp_path):
    from proteomics_clock import load_input

    csv_path = tmp_path / "test.csv"
    df = pd.DataFrame({"sample_id": ["S1"], "ProteinA": [1.5]})
    df.to_csv(csv_path, index=False)
    loaded = load_input(csv_path)
    assert loaded.shape == (1, 2)
    assert "ProteinA" in loaded.columns


# ---------------------------------------------------------------------------
# Integration tests
# ---------------------------------------------------------------------------


def test_run_analysis_with_demo(tmp_path):
    from proteomics_clock import run_analysis

    out_dir = tmp_path / "proteomics_report"
    result = run_analysis(
        output_dir=out_dir,
        organs=["Heart", "Brain", "Kidney"],
        generation="both",
        fold=1,
        convert_mortality_to_years=True,
        age_column="age",
        verbose=False,
        input_path=DEMO_FIXTURE,
    )

    assert result["n_samples"] >= 1
    assert (out_dir / "report.md").exists()
    assert (out_dir / "tables" / "predictions_gen1.csv").exists()
    assert (out_dir / "tables" / "predictions_gen2.csv").exists()
    assert (out_dir / "tables" / "prediction_summary.csv").exists()
    assert (out_dir / "tables" / "missing_proteins.csv").exists()
    assert (out_dir / "tables" / "clock_metadata.json").exists()

    gen1 = pd.read_csv(out_dir / "tables" / "predictions_gen1.csv")
    assert "sample_id" in gen1.columns

    report_text = (out_dir / "report.md").read_text()
    assert "ClawBio is a research and educational tool" in report_text
    assert "Goeminne" in report_text


def test_cli_smoke_with_demo(tmp_path):
    out_dir = tmp_path / "proteomics_cli"
    script = SKILL_DIR / "proteomics_clock.py"

    cmd = [
        sys.executable,
        str(script),
        "--demo",
        "--output",
        str(out_dir),
        "--organs",
        "Heart,Brain,Kidney",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True, timeout=120)

    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    assert data["n_samples"] >= 1
    assert (out_dir / "tables" / "predictions_gen1.csv").exists()
    assert (out_dir / "reproducibility" / "checksums.sha256").exists()
