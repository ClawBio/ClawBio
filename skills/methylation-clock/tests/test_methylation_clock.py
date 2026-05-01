import json
import subprocess
import sys
from pathlib import Path

import pandas as pd
import pytest

pytest.importorskip("pyaging")

sys.path.insert(0, str(Path(__file__).parent.parent))

from methylation_clock import DEFAULT_CLOCKS, parse_clock_list, run_analysis


SKILL_DIR = Path(__file__).resolve().parent.parent
REPO_ROOT = SKILL_DIR.parents[1]
FIXTURE = SKILL_DIR / "data" / "GSE139307_small.csv.gz"


def test_parse_clock_list_defaults():
    clocks = parse_clock_list(None)
    assert clocks == DEFAULT_CLOCKS


def test_run_analysis_with_fixture(tmp_path):
    out_dir = tmp_path / "methylation_report"
    result = run_analysis(
        output_dir=out_dir,
        clocks=["Horvath2013"],
        metadata_cols=["gender", "tissue_type", "dataset"],
        imputer_strategy="mean",
        verbose=False,
        input_path=FIXTURE,
        geo_id=None,
    )

    assert result["n_samples"] >= 1
    assert (out_dir / "report.md").exists()
    assert (out_dir / "tables" / "predictions.csv").exists()
    assert (out_dir / "tables" / "prediction_summary.csv").exists()
    assert (out_dir / "tables" / "missing_features.csv").exists()

    predictions = pd.read_csv(out_dir / "tables" / "predictions.csv")
    assert "sample_id" in predictions.columns
    assert any(col.lower() == "horvath2013" for col in predictions.columns)

    report_text = (out_dir / "report.md").read_text()
    assert "ClawBio is a research and educational tool" in report_text


def test_cli_smoke_with_fixture(tmp_path):
    out_dir = tmp_path / "methylation_cli"
    script = SKILL_DIR / "methylation_clock.py"

    cmd = [
        sys.executable,
        str(script),
        "--input",
        str(FIXTURE),
        "--output",
        str(out_dir),
        "--clocks",
        "Horvath2013",
        "--imputer-strategy",
        "mean",
    ]
    proc = subprocess.run(cmd, capture_output=True, text=True)

    assert proc.returncode == 0, proc.stderr
    data = json.loads(proc.stdout)
    assert data["n_samples"] >= 1
    assert (out_dir / "tables" / "predictions.csv").exists()
    assert (out_dir / "reproducibility" / "checksums.sha256").exists()
