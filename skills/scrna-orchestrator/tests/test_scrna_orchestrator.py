"""Tests for scRNA Orchestrator MVP."""

from __future__ import annotations

import importlib.util
import json
import shlex
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = SKILL_DIR / "scrna_orchestrator.py"
ORCHESTRATOR_PATH = SKILL_DIR.parent / "bio-orchestrator" / "orchestrator.py"


def _run_cmd(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True,
        text=True,
    )


def _require_scanpy() -> None:
    pytest.importorskip("scanpy")
    pytest.importorskip("anndata")


def _load_orchestrator_module():
    spec = importlib.util.spec_from_file_location("bio_orchestrator_module", ORCHESTRATOR_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_demo_end_to_end_outputs(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_output"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    expected = [
        output_dir / "report.md",
        output_dir / "result.json",
        output_dir / "figures" / "qc_violin.png",
        output_dir / "figures" / "umap_leiden.png",
        output_dir / "figures" / "marker_dotplot.png",
        output_dir / "tables" / "cluster_summary.csv",
        output_dir / "tables" / "markers_top.csv",
        output_dir / "tables" / "markers_top.tsv",
        output_dir / "reproducibility" / "commands.sh",
        output_dir / "reproducibility" / "environment.yml",
        output_dir / "reproducibility" / "checksums.sha256",
    ]
    for path in expected:
        assert path.exists(), f"Missing output file: {path}"


def test_demo_summary_in_result_json(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_output"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    payload = json.loads((output_dir / "result.json").read_text())
    summary = payload["summary"]
    assert summary["n_cells_before"] >= summary["n_cells_after"] > 0
    assert summary["n_clusters"] >= 2
    assert summary["n_hvg"] > 0


def test_markers_csv_tsv_schema_match(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_output"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    csv_df = pd.read_csv(output_dir / "tables" / "markers_top.csv")
    tsv_df = pd.read_csv(output_dir / "tables" / "markers_top.tsv", sep="\t")
    assert list(csv_df.columns) == list(tsv_df.columns)
    assert len(csv_df) > 0


def test_report_contains_key_stats(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_output"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    report_text = (output_dir / "report.md").read_text()
    assert "Cells before QC" in report_text
    assert "Cells after QC" in report_text
    assert "Leiden clusters" in report_text
    assert "HVG selected" in report_text
    assert "not a medical device" in report_text


def test_checksums_contains_key_outputs(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_output"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    checksums = (output_dir / "reproducibility" / "checksums.sha256").read_text()
    assert "report.md" in checksums
    assert "tables/markers_top.csv" in checksums
    assert "tables/markers_top.tsv" in checksums
    assert "figures/marker_dotplot.png" in checksums


def test_non_h5ad_input_rejected(tmp_path: Path):
    _require_scanpy()
    input_path = tmp_path / "invalid.csv"
    input_path.write_text("a,b\n1,2\n")

    result = _run_cmd(["--input", str(input_path), "--output", str(tmp_path / "out")])
    assert result.returncode != 0
    assert "Only .h5ad is supported" in result.stderr


def test_tiny_dataset_no_pca_crash(tmp_path: Path):
    _require_scanpy()
    from anndata import AnnData  # type: ignore

    input_path = tmp_path / "tiny.h5ad"
    output_dir = tmp_path / "tiny_output"

    x = np.array(
        [
            [1, 0],
            [2, 1],
            [3, 0],
            [4, 1],
        ],
        dtype=np.int32,
    )
    obs = pd.DataFrame(
        index=pd.Index([f"cell_{i}" for i in range(4)], dtype="object")
    )
    var = pd.DataFrame(index=pd.Index(["GeneA", "GeneB"], dtype="object"))
    AnnData(X=x, obs=obs, var=var).write_h5ad(input_path)

    result = _run_cmd(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_dir),
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--n-top-hvg",
            "2",
            "--n-neighbors",
            "2",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "report.md").exists()


def test_commands_sh_quotes_demo_output_path(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo output (quoted)"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert f"--output {shlex.quote(str(output_dir))}" in commands_sh


def test_orchestrator_no_rds_extension_route():
    module = _load_orchestrator_module()
    assert module.detect_skill_from_file(Path("x.rds")) is None
    assert module.detect_skill_from_file(Path("x.h5ad")) == "scrna-orchestrator"
