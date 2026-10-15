"""Tests for scRNA Orchestrator MVP."""

from __future__ import annotations

import gzip
import importlib.util
import json
import os
import shlex
import subprocess
import sys
from pathlib import Path
from types import SimpleNamespace

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


def _build_human_like_input(path: Path) -> None:
    from anndata import AnnData  # type: ignore

    genes = ["CD3D", "TRBC1", "MS4A1", "CD79A", "LYZ", "S100A8", "NKG7", "GNLY"]
    rng = np.random.default_rng(0)
    templates = [
        np.array([18, 16, 1, 1, 1, 1, 4, 3], dtype=np.int32),
        np.array([1, 1, 18, 16, 1, 1, 3, 1], dtype=np.int32),
        np.array([1, 1, 1, 1, 18, 16, 10, 8], dtype=np.int32),
    ]

    rows = []
    for template in templates:
        for _ in range(6):
            rows.append(rng.poisson(lam=template).astype(np.int32) + 1)

    obs = pd.DataFrame(index=pd.Index([f"cell_{i}" for i in range(len(rows))], dtype="object"))
    var = pd.DataFrame(index=pd.Index(genes, dtype="object"))
    AnnData(X=np.vstack(rows), obs=obs, var=var).write_h5ad(path)


def _build_integrated_latent_input(path: Path) -> None:
    from anndata import AnnData  # type: ignore

    genes = ["CD3D", "TRBC1", "MS4A1", "CD79A", "LYZ", "S100A8", "NKG7", "GNLY"]
    rng = np.random.default_rng(2)
    templates = [
        np.array([18, 16, 1, 1, 1, 1, 4, 3], dtype=np.int32),
        np.array([1, 1, 18, 16, 1, 1, 3, 1], dtype=np.int32),
        np.array([1, 1, 1, 1, 18, 16, 10, 8], dtype=np.int32),
    ]

    rows = []
    batches = []
    for batch_idx in range(2):
        for cluster_idx, template in enumerate(templates):
            for _ in range(4):
                rows.append(rng.poisson(lam=template + batch_idx).astype(np.int32) + 1)
                batches.append(f"batch_{batch_idx}")

    counts = np.vstack(rows)
    obs = pd.DataFrame(
        {
            "batch": batches,
            "sample_id": [f"cell_{i}" for i in range(len(rows))],
        },
        index=pd.Index([f"cell_{i}" for i in range(len(rows))], dtype="object"),
    )
    var = pd.DataFrame(index=pd.Index(genes, dtype="object"))
    adata = AnnData(X=counts.astype(np.float32), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.X = np.log1p(counts.astype(np.float32))
    latent = np.column_stack(
        [
            np.repeat([0.0, 3.5, 7.0], repeats=8),
            np.tile(np.repeat([0.0, 2.0], repeats=4), 3),
        ]
    ).astype(np.float32)
    adata.obsm["X_scvi"] = latent
    adata.uns["clawbio_scrna_embedding"] = {
        "source_skill": "scrna-embedding",
        "preferred_rep": "X_scvi",
        "counts_layer": "counts",
    }
    adata.write_h5ad(path)


def _build_within_cluster_latent_input(path: Path) -> None:
    from anndata import AnnData  # type: ignore

    genes = ["CD3D", "TRBC1", "MS4A1", "CD79A", "LYZ", "S100A8", "NKG7", "GNLY"]
    rng = np.random.default_rng(7)
    templates = [
        np.array([18, 16, 1, 1, 1, 1, 4, 3], dtype=np.int32),
        np.array([1, 1, 18, 16, 1, 1, 3, 1], dtype=np.int32),
        np.array([1, 1, 1, 1, 18, 16, 10, 8], dtype=np.int32),
    ]

    rows = []
    conditions = []
    cluster_labels = []
    for cluster_idx, template in enumerate(templates):
        for cell_idx in range(8):
            rows.append(rng.poisson(lam=template).astype(np.int32) + 1)
            cluster_labels.append(f"truth_{cluster_idx}")
            if cluster_idx == 2:
                conditions.append("condition_a")
            else:
                conditions.append("condition_a" if cell_idx < 4 else "condition_b")

    counts = np.vstack(rows)
    obs = pd.DataFrame(
        {
            "condition": conditions,
            "truth_cluster": cluster_labels,
            "sample_id": [f"cell_{i}" for i in range(len(rows))],
        },
        index=pd.Index([f"cell_{i}" for i in range(len(rows))], dtype="object"),
    )
    var = pd.DataFrame(index=pd.Index(genes, dtype="object"))
    adata = AnnData(X=counts.astype(np.float32), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()
    adata.X = np.log1p(counts.astype(np.float32))
    adata.obsm["X_scvi"] = np.column_stack(
        [
            np.repeat([0.0, 3.5, 7.0], repeats=8),
            np.tile(np.linspace(0.0, 1.4, num=8), 3),
        ]
    ).astype(np.float32)
    adata.uns["clawbio_scrna_embedding"] = {
        "source_skill": "scrna-embedding",
        "preferred_rep": "X_scvi",
        "counts_layer": "counts",
    }
    adata.write_h5ad(path)


def _build_10x_mtx_input(matrix_dir: Path, *, compressed: bool = False, prefix: str = "") -> Path:
    from scipy import io, sparse  # type: ignore

    matrix_dir.mkdir(parents=True, exist_ok=True)
    genes = ["CD3D", "TRBC1", "MS4A1", "CD79A", "LYZ", "S100A8", "NKG7", "GNLY"]
    rng = np.random.default_rng(0)
    templates = [
        np.array([18, 16, 1, 1, 1, 1, 4, 3], dtype=np.int32),
        np.array([1, 1, 18, 16, 1, 1, 3, 1], dtype=np.int32),
        np.array([1, 1, 1, 1, 18, 16, 10, 8], dtype=np.int32),
    ]

    rows = []
    for template in templates:
        for _ in range(6):
            rows.append(rng.poisson(lam=template).astype(np.int32) + 1)

    matrix = np.vstack(rows).T
    matrix_path = matrix_dir / f"{prefix}matrix.mtx"
    io.mmwrite(str(matrix_path), sparse.coo_matrix(matrix))

    barcode_path = matrix_dir / f"{prefix}barcodes.tsv"
    barcode_path.write_text(
        "\n".join(f"cell_{i}" for i in range(matrix.shape[1])) + "\n",
        encoding="utf-8",
    )

    features_path = matrix_dir / f"{prefix}features.tsv"
    features_path.write_text(
        "\n".join(f"gene_{i}\t{gene}\tGene Expression" for i, gene in enumerate(genes)) + "\n",
        encoding="utf-8",
    )

    if not compressed:
        return matrix_path

    for plain_path in (matrix_path, barcode_path, features_path):
        gz_path = plain_path.with_suffix(plain_path.suffix + ".gz")
        with plain_path.open("rb") as source_handle, gzip.open(gz_path, "wb") as dest_handle:
            dest_handle.write(source_handle.read())
        plain_path.unlink()

    return matrix_dir / f"{prefix}matrix.mtx.gz"


def _make_args(output_dir: Path, **overrides) -> SimpleNamespace:
    defaults = {
        "input": None,
        "output": str(output_dir),
        "demo": False,
        "min_genes": 200,
        "min_cells": 3,
        "max_mt_pct": 20.0,
        "n_top_hvg": 2000,
        "n_pcs": 50,
        "n_neighbors": 15,
        "use_rep": "auto",
        "leiden_resolution": 1.0,
        "random_state": 0,
        "top_markers": 10,
        "doublet_method": "none",
        "annotate": "none",
        "annotation_model": "Immune_All_Low",
        "contrast_groupby": None,
        "contrast_group1": None,
        "contrast_group2": None,
        "contrast_top_genes": None,
        "contrast_volcano": False,
        "de_groupby": None,
        "de_group1": None,
        "de_group2": None,
        "de_top_genes": None,
        "de_volcano": False,
    }
    defaults.update(overrides)
    return SimpleNamespace(**defaults)


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


def test_dataset_level_contrast_outputs_and_result_metadata(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "dataset_contrast_output"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--contrast-groupby",
            "demo_truth",
            "--contrast-scope",
            "dataset",
            "--contrast-top-genes",
            "7",
        ]
    )
    assert result.returncode == 0, result.stderr

    de_full = output_dir / "tables" / "contrastive_markers_full.csv"
    de_top = output_dir / "tables" / "contrastive_markers_top.csv"
    assert de_full.exists()
    assert de_top.exists()

    de_full_df = pd.read_csv(de_full)
    de_top_df = pd.read_csv(de_top)
    assert len(de_full_df) > 0
    assert set(["scope", "groupby", "group1", "group2", "comparison_id"]).issubset(de_full_df.columns)
    assert de_full_df["comparison_id"].nunique() == 3
    assert de_full_df["scope"].astype(str).eq("dataset").all()
    assert (de_top_df.groupby("comparison_id").size() <= 7).all()

    payload = json.loads((output_dir / "result.json").read_text())
    de_meta = payload["data"]["contrasts"]["dataset"]
    assert de_meta["enabled"] is True
    assert de_meta["groupby"] == "demo_truth"
    assert de_meta["n_contrasts"] == 3
    assert de_meta["n_rows_full"] == len(de_full_df)
    assert de_meta["top_table"] == "contrastive_markers_top.csv"
    assert "contrastive_markers_full.csv" in payload["data"]["tables"]
    assert "contrastive_markers_top.csv" in payload["data"]["tables"]
    assert payload["data"]["contrasts"]["within_cluster"]["enabled"] is False

    report_text = (output_dir / "report.md").read_text()
    assert "Dataset-Level Contrastive Markers" in report_text
    assert "Pairwise comparisons run: **3**" in report_text


def test_within_cluster_contrast_outputs_and_skip_metadata(tmp_path: Path):
    _require_scanpy()
    input_path = tmp_path / "within_cluster_latent.h5ad"
    _build_within_cluster_latent_input(input_path)
    output_dir = tmp_path / "within_cluster_output"
    result = _run_cmd(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_dir),
            "--use-rep",
            "X_scvi",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--n-top-hvg",
            "4",
            "--n-neighbors",
            "4",
            "--contrast-groupby",
            "condition",
            "--contrast-scope",
            "within-cluster",
            "--contrast-clusterby",
            "truth_cluster",
            "--contrast-top-genes",
            "4",
        ]
    )
    assert result.returncode == 0, result.stderr

    within_full = output_dir / "tables" / "within_cluster_contrastive_markers_full.csv"
    within_top = output_dir / "tables" / "within_cluster_contrastive_markers_top.csv"
    assert within_full.exists()
    assert within_top.exists()

    within_full_df = pd.read_csv(within_full)
    within_top_df = pd.read_csv(within_top)
    assert set(["cluster", "scope", "groupby", "group1", "group2", "comparison_id"]).issubset(within_full_df.columns)
    assert within_full_df["scope"].astype(str).eq("within-cluster").all()
    assert within_full_df["cluster"].nunique() >= 2
    assert (within_top_df.groupby(["cluster", "comparison_id"]).size() <= 4).all()

    payload = json.loads((output_dir / "result.json").read_text())
    within_meta = payload["data"]["contrasts"]["within_cluster"]
    assert within_meta["enabled"] is True
    assert within_meta["groupby"] == "condition"
    assert within_meta["clusterby"] == "truth_cluster"
    assert within_meta["n_contrasts"] >= 2
    assert within_meta["skipped_contrasts"] >= 1
    assert within_meta["n_clusters_with_results"] >= 2

    report_text = (output_dir / "report.md").read_text()
    assert "Within-Cluster Contrastive Markers" in report_text
    assert "Skipped cluster/comparison pairs" in report_text


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
    assert "Unsupported input" in result.stderr


def test_10x_mtx_directory_input_runs_end_to_end(tmp_path: Path):
    _require_scanpy()
    matrix_dir = tmp_path / "filtered_feature_bc_matrix"
    _build_10x_mtx_input(matrix_dir, compressed=False)
    output_dir = tmp_path / "mtx_output"

    result = _run_cmd(
        [
            "--input",
            str(matrix_dir),
            "--output",
            str(output_dir),
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--n-top-hvg",
            "8",
            "--n-pcs",
            "5",
            "--n-neighbors",
            "4",
            "--top-markers",
            "4",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "report.md").exists()

    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["data"]["input"]["format"] == "10x_mtx"
    assert set(payload["data"]["input"]["files"]) == {"matrix.mtx", "barcodes.tsv", "features.tsv"}
    assert payload["input_checksum"].startswith("sha256:")

    report_text = (output_dir / "report.md").read_text()
    assert "Input format" in report_text
    assert "matrix.mtx" in report_text
    assert "barcodes.tsv" in report_text
    assert "features.tsv" in report_text

    checksums = (output_dir / "reproducibility" / "checksums.sha256").read_text()
    assert "matrix.mtx" in checksums
    assert "barcodes.tsv" in checksums
    assert "features.tsv" in checksums

    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert f"--input {shlex.quote(str(matrix_dir))}" in commands_sh


def test_load_data_accepts_direct_prefixed_gz_mtx_file(tmp_path: Path):
    _require_scanpy()
    module = _load_scrna_module()
    matrix_dir = tmp_path / "prefixed_matrix"
    matrix_path = _build_10x_mtx_input(matrix_dir, compressed=True, prefix="sample_")

    adata, input_path, is_demo, demo_source, input_source, latent_context = module.load_data(
        str(matrix_path),
        demo=False,
        random_state=0,
        use_rep="auto",
    )

    assert adata.n_obs == 18
    assert adata.n_vars == 8
    assert input_path == matrix_path
    assert is_demo is False
    assert latent_context["resolved_use_rep"] == ""
    assert demo_source is None
    assert input_source["format"] == "10x_mtx"


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
    assert "raw-count .h5ad or 10x single-cell input" in result.stderr
    assert "pbmc3k_processed" in result.stderr


def test_commands_sh_quotes_demo_output_path(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo output (quoted)"
    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode == 0, result.stderr

    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert f"--output {shlex.quote(str(output_dir))}" in commands_sh


def test_commands_sh_contains_contrast_flags_when_enabled(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "demo_de_commands"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--contrast-groupby",
            "demo_truth",
            "--contrast-scope",
            "dataset",
            "--contrast-top-genes",
            "9",
        ]
    )
    assert result.returncode == 0, result.stderr

    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert "--contrast-groupby demo_truth" in commands_sh
    assert "--contrast-scope dataset" in commands_sh
    assert "--contrast-top-genes 9" in commands_sh


def test_demo_doublet_detection_outputs_summary_metadata(tmp_path: Path):
    _require_scanpy()
    pytest.importorskip("scrublet")
    output_dir = tmp_path / "doublet_output"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--doublet-method",
            "scrublet",
        ]
    )
    assert result.returncode == 0, result.stderr

    doublet_table = output_dir / "tables" / "doublet_summary.csv"
    assert doublet_table.exists()
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["n_predicted_doublets"] >= 0
    assert payload["data"]["doublet"]["method"] == "scrublet"
    assert payload["data"]["doublet"]["table"] == "doublet_summary.csv"
    report_text = (output_dir / "report.md").read_text()
    assert "Doublet Detection" in report_text
    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert "--doublet-method scrublet" in commands_sh


def test_doublet_and_dataset_contrasts_can_run_together(tmp_path: Path):
    _require_scanpy()
    pytest.importorskip("scrublet")
    output_dir = tmp_path / "doublet_de_output"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--doublet-method",
            "scrublet",
            "--contrast-groupby",
            "demo_truth",
            "--contrast-scope",
            "dataset",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "tables" / "doublet_summary.csv").exists()
    assert (output_dir / "tables" / "contrastive_markers_full.csv").exists()


def test_contrast_groupby_missing_rejected_when_other_flags_are_used(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "contrast_missing_groupby"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--contrast-scope",
            "both",
        ]
    )
    assert result.returncode != 0
    assert "--contrast-groupby is required" in result.stderr


def test_contrast_top_genes_must_be_positive(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "contrast_top_genes_zero"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--contrast-groupby",
            "demo_truth",
            "--contrast-top-genes",
            "0",
        ]
    )
    assert result.returncode != 0
    assert "--contrast-top-genes must be >= 1" in result.stderr


def test_contrast_groupby_missing_rejected(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "contrast_bad_groupby"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--contrast-groupby",
            "not_a_column",
            "--contrast-scope",
            "dataset",
        ]
    )
    assert result.returncode != 0
    assert "Contrastive marker groupby column not found" in result.stderr


def test_contrast_groupby_requires_multiple_observed_groups(tmp_path: Path):
    _require_scanpy()
    from anndata import AnnData  # type: ignore

    input_path = tmp_path / "single_group.h5ad"
    output_dir = tmp_path / "single_group_output"
    x = np.array(
        [
            [3, 1, 0, 0],
            [4, 2, 0, 0],
            [0, 0, 5, 3],
            [0, 0, 4, 2],
        ],
        dtype=np.int32,
    )
    obs = pd.DataFrame(
        {
            "condition": ["only_group"] * 4,
        },
        index=[f"cell_{i}" for i in range(4)],
    )
    var = pd.DataFrame(index=["GeneA", "GeneB", "GeneC", "GeneD"])
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
            "4",
            "--n-neighbors",
            "2",
            "--contrast-groupby",
            "condition",
            "--contrast-scope",
            "dataset",
        ]
    )
    assert result.returncode != 0
    assert "must contain at least 2 observed groups" in result.stderr


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


def test_clawbio_run_scrna_accepts_whitelisted_feature_flags(tmp_path: Path):
    _require_scanpy()
    pytest.importorskip("scrublet")
    output_dir = tmp_path / "clawbio_scrna_feature_output"

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
            "--doublet-method",
            "scrublet",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "result.json").exists()
    assert (output_dir / "tables" / "doublet_summary.csv").exists()


def test_clawbio_run_scrna_accepts_whitelisted_contrast_flags(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "clawbio_scrna_de_output"

    result = _run_clawbio_scrna_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--contrast-groupby",
            "demo_truth",
            "--contrast-scope",
            "dataset",
            "--contrast-top-genes",
            "8",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "tables" / "contrastive_markers_full.csv").exists()
    assert (output_dir / "tables" / "contrastive_markers_top.csv").exists()


def test_doublet_missing_dependency_message(monkeypatch: pytest.MonkeyPatch):
    from anndata import AnnData  # type: ignore

    module = _load_scrna_module()
    adata = AnnData(
        X=np.array([[1, 0, 3], [2, 1, 0], [1, 1, 2]], dtype=np.int32),
        obs=pd.DataFrame(index=pd.Index(["cell0", "cell1", "cell2"], dtype="object")),
        var=pd.DataFrame(index=pd.Index(["GeneA", "GeneB", "GeneC"], dtype="object")),
    )

    def _raise_missing():
        raise RuntimeError(
            "scrublet is required for --doublet-method scrublet. Install it with: pip install scrublet"
        )

    monkeypatch.setattr(module, "_import_scrublet", _raise_missing)
    with pytest.raises(RuntimeError, match="pip install scrublet"):
        module.run_doublet_detection(adata, method="scrublet", random_state=0)


def test_parse_args_annotation_defaults_and_override(monkeypatch: pytest.MonkeyPatch):
    module = _load_scrna_module()

    monkeypatch.setattr(sys, "argv", ["scrna_orchestrator.py", "--demo"])
    args = module.parse_args()
    assert args.annotate == "none"
    assert args.annotation_model == "Immune_All_Low"

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "scrna_orchestrator.py",
            "--demo",
            "--annotate",
            "celltypist",
            "--annotation-model",
            "CustomModel",
        ],
    )
    args = module.parse_args()
    assert args.annotate == "celltypist"
    assert args.annotation_model == "CustomModel"


def test_celltypist_annotation_pipeline_writes_cluster_summary(
    tmp_path: Path,
    monkeypatch: pytest.MonkeyPatch,
):
    _require_scanpy()
    module = _load_scrna_module()
    input_path = tmp_path / "human_like.h5ad"
    output_dir = tmp_path / "annotated_output"
    _build_human_like_input(input_path)

    model_dir = tmp_path / "models"
    model_dir.mkdir()
    model_path = model_dir / "Immune_All_Low.pkl"
    model_path.write_text("fake model placeholder", encoding="utf-8")

    genes = np.array(["CD3D", "TRBC1", "MS4A1", "CD79A", "LYZ", "S100A8", "NKG7", "GNLY"], dtype=object)

    class FakeModel:
        def __init__(self, features: np.ndarray):
            self.features = features

    class FakeModelLoader:
        @staticmethod
        def load(model: str):
            assert model.endswith("Immune_All_Low.pkl")
            return FakeModel(features=genes)

    def fake_annotate(adata, model, majority_voting: bool = False):
        assert majority_voting is False
        assert np.array_equal(model.features, genes)
        cluster_order = sorted(
            adata.obs["leiden"].astype(str).unique().tolist(),
            key=module._cluster_sort_key,
        )
        label_map = {}
        labels = ["T cell", "B cell", "Myeloid"]
        for index, cluster in enumerate(cluster_order):
            label_map[cluster] = labels[min(index, len(labels) - 1)]

        predicted = []
        probabilities = []
        for cell_name in adata.obs_names:
            cluster = str(adata.obs.loc[cell_name, "leiden"])
            label = label_map[cluster]
            predicted.append(label)
            probabilities.append(
                {
                    "T cell": 0.92 if label == "T cell" else 0.04,
                    "B cell": 0.91 if label == "B cell" else 0.04,
                    "Myeloid": 0.90 if label == "Myeloid" else 0.04,
                }
            )

        return SimpleNamespace(
            predicted_labels=pd.DataFrame({"predicted_labels": predicted}, index=adata.obs_names),
            probability_matrix=pd.DataFrame(probabilities, index=adata.obs_names),
            adata=adata,
        )

    fake_celltypist = SimpleNamespace(
        models=SimpleNamespace(models_path=str(model_dir), Model=FakeModelLoader),
        annotate=fake_annotate,
    )
    monkeypatch.setattr(module, "_import_celltypist", lambda: fake_celltypist)

    args = _make_args(
        output_dir,
        input=str(input_path),
        min_genes=1,
        min_cells=1,
        n_top_hvg=8,
        n_pcs=5,
        n_neighbors=4,
        leiden_resolution=0.7,
        top_markers=4,
        annotate="celltypist",
        annotation_model="Immune_All_Low",
    )
    result = module.run_pipeline(args)
    assert result["n_clusters"] >= 2

    ann = pd.read_csv(output_dir / "tables" / "cluster_annotations.csv")
    assert list(ann.columns) == [
        "cluster",
        "n_cells",
        "predicted_cell_type",
        "support_fraction",
        "mean_confidence",
        "annotation_model",
    ]
    assert set(ann["predicted_cell_type"]).issubset({"T cell", "B cell", "Myeloid"})
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["data"]["annotation"]["backend"] == "celltypist"
    assert payload["data"]["annotation"]["model"] == "Immune_All_Low.pkl"
    assert payload["summary"]["n_clusters_annotated"] == len(ann)
    report_text = (output_dir / "report.md").read_text().lower()
    assert "putative" in report_text
    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert "--annotate celltypist" in commands_sh
    assert "--annotation-model Immune_All_Low" in commands_sh


def test_annotation_missing_local_model_fails_actionably(tmp_path: Path):
    _require_scanpy()
    pytest.importorskip("celltypist")
    output_dir = tmp_path / "annotation_missing_model"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--annotate",
            "celltypist",
            "--annotation-model",
            "DefinitelyMissingModel",
        ]
    )
    assert result.returncode != 0
    assert "Runtime downloads are disabled" in result.stderr
    assert "download_models" in result.stderr


def test_orchestrator_no_rds_extension_route():
    module = _load_orchestrator_module()
    assert module.detect_skill_from_file(Path("x.rds")) is None
    assert module.detect_skill_from_file(Path("x.h5ad")) == "scrna-orchestrator"
    assert module.detect_skill_from_file(Path("matrix.mtx")) == "scrna-orchestrator"
    assert module.detect_skill_from_file(Path("matrix.mtx.gz")) == "scrna-orchestrator"


def test_10x_directory_input_runs(tmp_path: Path):
    _require_scanpy()
    input_dir = _build_10x_mtx_input(tmp_path / "tenx_dir")
    output_dir = tmp_path / "tenx_output"

    result = _run_cmd(
        [
            "--input",
            str(input_dir),
            "--output",
            str(output_dir),
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--n-top-hvg",
            "8",
            "--n-neighbors",
            "4",
        ]
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["data"]["input"]["format"] == "10x_mtx"


def test_integrated_latent_input_auto_mode_runs(tmp_path: Path):
    _require_scanpy()
    input_path = tmp_path / "integrated.h5ad"
    output_dir = tmp_path / "latent_auto_output"
    _build_integrated_latent_input(input_path)

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
            "4",
            "--n-neighbors",
            "4",
            "--contrast-groupby",
            "batch",
            "--contrast-scope",
            "dataset",
        ]
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["graph_basis"] == "X_scvi"
    assert payload["summary"]["use_rep_resolved"] == "X_scvi"
    assert payload["summary"]["counts_layer"] == "counts"
    assert payload["data"]["contrasts"]["dataset"]["enabled"] is True
    assert (output_dir / "tables" / "contrastive_markers_full.csv").exists()


def test_integrated_latent_input_explicit_and_none_modes(tmp_path: Path):
    _require_scanpy()
    input_path = tmp_path / "integrated.h5ad"
    output_dir_explicit = tmp_path / "latent_explicit_output"
    output_dir_none = tmp_path / "latent_none_output"
    _build_integrated_latent_input(input_path)

    explicit_result = _run_cmd(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_dir_explicit),
            "--use-rep",
            "X_scvi",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--n-top-hvg",
            "4",
            "--n-neighbors",
            "4",
        ]
    )
    assert explicit_result.returncode == 0, explicit_result.stderr
    explicit_payload = json.loads((output_dir_explicit / "result.json").read_text())
    assert explicit_payload["summary"]["graph_basis"] == "X_scvi"

    none_result = _run_cmd(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_dir_none),
            "--use-rep",
            "none",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--n-top-hvg",
            "4",
            "--n-neighbors",
            "4",
        ]
    )
    assert none_result.returncode == 0, none_result.stderr
    none_payload = json.loads((output_dir_none / "result.json").read_text())
    assert none_payload["summary"]["graph_basis"] == "pca"
    assert none_payload["summary"]["use_rep_resolved"] == ""


def test_matrix_mtx_gz_input_runs(tmp_path: Path):
    _require_scanpy()
    matrix_path = _build_10x_mtx_input(tmp_path / "tenx_gz", compressed=True)
    output_dir = tmp_path / "tenx_gz_output"

    result = _run_cmd(
        [
            "--input",
            str(matrix_path),
            "--output",
            str(output_dir),
            "--min-genes",
            "1",
            "--min-cells",
            "1",
            "--n-top-hvg",
            "8",
            "--n-neighbors",
            "4",
        ]
    )
    assert result.returncode == 0, result.stderr
    checksums = (output_dir / "reproducibility" / "checksums.sha256").read_text()
    assert "barcodes.tsv.gz" in checksums


def test_invalid_contrast_scope_is_rejected(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "invalid_scope_output"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--contrast-groupby",
            "demo_truth",
            "--contrast-scope",
            "not-a-scope",
        ]
    )
    assert result.returncode != 0
    assert "invalid choice" in result.stderr


def test_clawbio_run_scrna_rejects_removed_de_flags(tmp_path: Path):
    _require_scanpy()
    output_dir = tmp_path / "removed_de_flags_output"
    result = _run_clawbio_scrna_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--de-groupby",
            "demo_truth",
        ]
    )
    assert result.returncode != 0
    assert "unrecognized arguments: --de-groupby" in result.stderr
