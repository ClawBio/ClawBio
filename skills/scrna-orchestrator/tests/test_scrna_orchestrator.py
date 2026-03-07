"""Tests for scRNA Orchestrator MVP."""

from __future__ import annotations

import importlib.util
import json
import os
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
REPO_ROOT = SKILL_DIR.parent.parent
CLAWBIO_PATH = REPO_ROOT / "clawbio.py"


def _run_cmd(
    args: list[str],
    env_overrides: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    run_env = os.environ.copy()
    run_env.setdefault("CLAWBIO_SCRNA_DEMO_SOURCE", "synthetic")
    if env_overrides:
        run_env.update(env_overrides)

    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True,
        text=True,
        env=run_env,
        cwd=str(REPO_ROOT),
    )


def _run_clawbio_scrna_cmd(
    args: list[str],
    env_overrides: dict[str, str] | None = None,
) -> subprocess.CompletedProcess:
    run_env = os.environ.copy()
    run_env.setdefault("CLAWBIO_SCRNA_DEMO_SOURCE", "synthetic")
    if env_overrides:
        run_env.update(env_overrides)

    return subprocess.run(
        [sys.executable, str(CLAWBIO_PATH), "run", "scrna"] + args,
        capture_output=True,
        text=True,
        env=run_env,
        cwd=str(REPO_ROOT),
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


def _load_scrna_module():
    spec = importlib.util.spec_from_file_location("scrna_module", SCRIPT_PATH)
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
    assert payload["data"]["demo_source"] in {"synthetic_forced", "synthetic_fallback"}


def test_demo_prefers_pbmc3k_when_available(monkeypatch: pytest.MonkeyPatch):
    from anndata import AnnData  # type: ignore

    module = _load_scrna_module()
    fake_adata = AnnData(
        X=np.array([[0, 0], [3, 1]], dtype=np.int32),
        obs=pd.DataFrame(index=pd.Index(["cell0", "cell1"], dtype="object")),
        var=pd.DataFrame(index=pd.Index(["GeneA", "GeneB"], dtype="object")),
    )

    class FakeDatasets:
        @staticmethod
        def pbmc3k():
            return fake_adata.copy()

    class FakePP:
        @staticmethod
        def filter_cells(adata, min_counts: int):
            totals = np.ravel(np.asarray(adata.X.sum(axis=1)))
            keep_mask = totals >= min_counts
            adata._inplace_subset_obs(keep_mask)

    fake_scanpy = type("FakeScanpy", (), {"datasets": FakeDatasets, "pp": FakePP})
    monkeypatch.setattr(module, "_import_scanpy", lambda: fake_scanpy)

    adata, source = module.load_demo_adata(random_state=0, demo_source_policy="auto")
    assert source == "pbmc3k_raw"
    assert adata.n_obs == 1


def test_demo_pbmc3k_failure_falls_back_to_synthetic(monkeypatch: pytest.MonkeyPatch, capsys: pytest.CaptureFixture[str]):
    module = _load_scrna_module()

    class FakeDatasets:
        @staticmethod
        def pbmc3k():
            raise RuntimeError("network unavailable")

    class FakePP:
        @staticmethod
        def filter_cells(_adata, min_counts: int):
            return None

    fake_scanpy = type("FakeScanpy", (), {"datasets": FakeDatasets, "pp": FakePP})
    monkeypatch.setattr(module, "_import_scanpy", lambda: fake_scanpy)

    adata, source = module.load_demo_adata(random_state=0, demo_source_policy="auto")
    assert source == "synthetic_fallback"
    assert adata.n_obs > 0
    stderr_text = capsys.readouterr().err
    assert "falling back to synthetic demo data" in stderr_text.lower()


def test_markers_csv_tsv_schema_match(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_output"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    csv_df = pd.read_csv(output_dir / "tables" / "markers_top.csv")
    tsv_df = pd.read_csv(output_dir / "tables" / "markers_top.tsv", sep="\t")
    assert list(csv_df.columns) == list(tsv_df.columns)
    assert len(csv_df) > 0


def test_de_two_group_outputs_and_result_metadata(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "de_demo_output"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "demo_truth",
            "--de-group1",
            "cluster_0",
            "--de-group2",
            "cluster_1",
            "--de-top-genes",
            "7",
        ]
    )
    assert result.returncode == 0, result.stderr

    de_full = output_dir / "tables" / "de_full.csv"
    de_top = output_dir / "tables" / "de_top.csv"
    assert de_full.exists()
    assert de_top.exists()

    de_full_df = pd.read_csv(de_full)
    de_top_df = pd.read_csv(de_top)
    assert len(de_full_df) > 0
    assert 0 < len(de_top_df) <= 7

    payload = json.loads((output_dir / "result.json").read_text())
    de_meta = payload["data"]["de"]
    assert de_meta["enabled"] is True
    assert de_meta["groupby"] == "demo_truth"
    assert de_meta["group1"] == "cluster_0"
    assert de_meta["group2"] == "cluster_1"
    assert de_meta["n_genes_full"] == len(de_full_df)
    assert de_meta["top_table"] == "de_top.csv"
    assert de_meta["volcano_plot"] == ""
    assert "de_full.csv" in payload["data"]["tables"]
    assert "de_top.csv" in payload["data"]["tables"]
    assert "de_volcano.png" not in payload["data"]["figures"]

    report_text = (output_dir / "report.md").read_text()
    assert "Differential Expression (Two-Group)" in report_text
    assert "cluster_0" in report_text
    assert "cluster_1" in report_text
    assert "Volcano plot: not generated" in report_text


def test_de_volcano_generated_when_requested(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "de_volcano_output"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "demo_truth",
            "--de-group1",
            "cluster_0",
            "--de-group2",
            "cluster_1",
            "--de-top-genes",
            "10",
            "--de-volcano",
        ]
    )
    assert result.returncode == 0, result.stderr

    volcano_path = output_dir / "figures" / "de_volcano.png"
    assert volcano_path.exists()

    payload = json.loads((output_dir / "result.json").read_text())
    de_meta = payload["data"]["de"]
    assert de_meta["volcano_plot"] == "de_volcano.png"
    assert "de_volcano.png" in payload["data"]["figures"]

    report_text = (output_dir / "report.md").read_text()
    assert "DE Volcano" in report_text


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


def test_processed_input_rejected_early_with_actionable_message(tmp_path: Path):
    _require_scanpy()
    from anndata import AnnData  # type: ignore

    input_path = tmp_path / "processed_like.h5ad"
    output_dir = tmp_path / "processed_out"

    x = np.array(
        [
            [0.2, -0.1, 1.5],
            [0.4, 0.7, 0.0],
            [1.1, 0.3, 0.2],
        ],
        dtype=np.float32,
    )
    obs = pd.DataFrame(index=pd.Index([f"cell_{i}" for i in range(3)], dtype="object"))
    var = pd.DataFrame(index=pd.Index(["GeneA", "GeneB", "GeneC"], dtype="object"))
    adata = AnnData(X=x, obs=obs, var=var)
    adata.uns["neighbors"] = {"params": {"n_neighbors": 15}}
    adata.write_h5ad(input_path)

    result = _run_cmd(["--input", str(input_path), "--output", str(output_dir)])
    assert result.returncode != 0
    assert "raw-count .h5ad input" in result.stderr
    assert "pbmc3k_processed" in result.stderr


def test_commands_sh_quotes_demo_output_path(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo output (quoted)"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert f"--output {shlex.quote(str(output_dir))}" in commands_sh


def test_commands_sh_contains_de_flags_when_enabled(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_de_commands"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "demo_truth",
            "--de-group1",
            "cluster_0",
            "--de-group2",
            "cluster_1",
            "--de-top-genes",
            "9",
            "--de-volcano",
        ]
    )
    assert result.returncode == 0, result.stderr

    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert "--de-groupby demo_truth" in commands_sh
    assert "--de-group1 cluster_0" in commands_sh
    assert "--de-group2 cluster_1" in commands_sh
    assert "--de-top-genes 9" in commands_sh
    assert "--de-volcano" in commands_sh


def test_de_requires_all_group_flags(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "de_incomplete"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "demo_truth",
        ]
    )
    assert result.returncode != 0
    assert "DE requires --de-groupby, --de-group1, and --de-group2 together" in result.stderr


def test_de_volcano_requires_de_flags(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "de_volcano_requires_flags"
    result = _run_cmd(["--demo", "--output", str(output_dir), "--de-volcano"])
    assert result.returncode != 0
    assert "--de-volcano requires --de-groupby, --de-group1, and --de-group2" in result.stderr


def test_de_groupby_missing_rejected(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "de_missing_groupby"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "not_a_column",
            "--de-group1",
            "cluster_0",
            "--de-group2",
            "cluster_1",
        ]
    )
    assert result.returncode != 0
    assert "DE groupby column not found" in result.stderr


def test_de_groups_missing_rejected(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "de_missing_groups"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "demo_truth",
            "--de-group1",
            "missing_a",
            "--de-group2",
            "missing_b",
        ]
    )
    assert result.returncode != 0
    assert "DE group value(s) not found" in result.stderr


def test_clawbio_run_scrna_accepts_whitelisted_tuning_flags(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "clawbio_scrna_output"

    result = _run_clawbio_scrna_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--n-pcs",
            "20",
            "--n-neighbors",
            "10",
            "--leiden-resolution",
            "0.6",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--top-markers",
            "5",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "result.json").exists()


def test_clawbio_run_scrna_accepts_whitelisted_de_flags(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "clawbio_scrna_de_output"

    result = _run_clawbio_scrna_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "demo_truth",
            "--de-group1",
            "cluster_0",
            "--de-group2",
            "cluster_1",
            "--de-top-genes",
            "8",
            "--de-volcano",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "tables" / "de_full.csv").exists()
    assert (output_dir / "tables" / "de_top.csv").exists()
    assert (output_dir / "figures" / "de_volcano.png").exists()


def test_orchestrator_no_rds_extension_route():
    module = _load_orchestrator_module()
    assert module.detect_skill_from_file(Path("x.rds")) is None
    assert module.detect_skill_from_file(Path("x.h5ad")) == "scrna-orchestrator"
