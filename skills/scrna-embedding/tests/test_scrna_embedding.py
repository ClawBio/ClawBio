"""Tests for the scRNA Embedding skill."""

from __future__ import annotations

import gzip
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
SCRIPT_PATH = SKILL_DIR / "scrna_embedding.py"
ORCHESTRATOR_PATH = SKILL_DIR.parent / "bio-orchestrator" / "orchestrator.py"
REPO_ROOT = SKILL_DIR.parent.parent
CLAWBIO_PATH = REPO_ROOT / "clawbio.py"


def _run_cmd(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT_PATH)] + args,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )


def _run_clawbio_cmd(args: list[str]) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(CLAWBIO_PATH), "run", "scrna-embedding"] + args,
        capture_output=True,
        text=True,
        cwd=str(REPO_ROOT),
    )


def _require_scvi_stack() -> None:
    pytest.importorskip("scanpy")
    pytest.importorskip("anndata")
    pytest.importorskip("scvi")


def _load_orchestrator_module():
    spec = importlib.util.spec_from_file_location("bio_orchestrator_module", ORCHESTRATOR_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def _build_input_matrix() -> tuple[np.ndarray, list[str], pd.DataFrame]:
    rng = np.random.default_rng(0)
    n_genes = 72
    genes = [f"Gene{i:03d}" for i in range(n_genes)]
    for idx in range(5):
        genes[idx] = f"MT-GENE{idx:02d}"

    templates = []
    for cluster_idx in range(3):
        base = np.full(n_genes, 1.0, dtype=np.float64)
        marker_start = 10 + cluster_idx * 8
        marker_end = marker_start + 6
        base[marker_start:marker_end] = 12.0
        templates.append(base)

    batch_effects = [np.ones(n_genes), np.ones(n_genes)]
    batch_effects[0][40:48] = 2.0
    batch_effects[1][48:56] = 1.8

    rows = []
    truths = []
    batches = []
    for batch_idx, batch_effect in enumerate(batch_effects):
        for cluster_idx, template in enumerate(templates):
            lam = np.clip(template * batch_effect, 0.2, None)
            counts = rng.poisson(lam=lam, size=(12, n_genes))
            counts = np.round(counts * rng.lognormal(0.0, 0.2, size=(12, 1))).astype(np.int32)
            rows.append(counts + 1)
            truths.extend([f"cluster_{cluster_idx}"] * 12)
            batches.extend([f"batch_{batch_idx}"] * 12)

    obs = pd.DataFrame(
        {
            "truth": truths,
            "batch": batches,
        },
        index=pd.Index([f"cell_{idx}" for idx in range(len(truths))], dtype="object"),
    )
    return np.vstack(rows), genes, obs


def _build_h5ad_input(path: Path, *, with_counts_layer: bool = False, processed_x: bool = False) -> None:
    from anndata import AnnData  # type: ignore

    x, genes, obs = _build_input_matrix()
    var = pd.DataFrame(index=pd.Index(genes, dtype="object"))
    matrix = x.astype(np.float32) if processed_x else x.astype(np.int32)
    adata = AnnData(X=matrix, obs=obs, var=var)
    if processed_x:
        adata.X = np.log1p(adata.X)
        adata.uns["neighbors"] = {"params": {"n_neighbors": 10}}
    if with_counts_layer:
        adata.layers["counts"] = x.copy()
    adata.write_h5ad(path)


def _build_10x_input(matrix_dir: Path, *, compressed: bool = False) -> Path:
    from scipy import io, sparse  # type: ignore

    x, genes, _obs = _build_input_matrix()
    matrix_dir.mkdir(parents=True, exist_ok=True)
    matrix_path = matrix_dir / "matrix.mtx"
    io.mmwrite(str(matrix_path), sparse.coo_matrix(x.T))

    (matrix_dir / "barcodes.tsv").write_text(
        "\n".join(f"cell_{idx}" for idx in range(x.shape[0])) + "\n",
        encoding="utf-8",
    )
    (matrix_dir / "features.tsv").write_text(
        "\n".join(f"gene_{idx}\t{gene}\tGene Expression" for idx, gene in enumerate(genes)) + "\n",
        encoding="utf-8",
    )

    if not compressed:
        return matrix_dir

    for plain_path in (matrix_path, matrix_dir / "barcodes.tsv", matrix_dir / "features.tsv"):
        gz_path = plain_path.with_suffix(plain_path.suffix + ".gz")
        with plain_path.open("rb") as source_handle, gzip.open(gz_path, "wb") as dest_handle:
            dest_handle.write(source_handle.read())
        plain_path.unlink()
    return matrix_dir / "matrix.mtx.gz"


def test_demo_end_to_end_outputs(tmp_path: Path):
    _require_scvi_stack()
    output_dir = tmp_path / "demo_output"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr

    expected = [
        output_dir / "report.md",
        output_dir / "result.json",
        output_dir / "integrated.h5ad",
        output_dir / "figures" / "umap_scvi_latent.png",
        output_dir / "tables" / "latent_embeddings.csv",
        output_dir / "reproducibility" / "commands.sh",
        output_dir / "reproducibility" / "environment.yml",
        output_dir / "reproducibility" / "checksums.sha256",
    ]
    for path in expected:
        assert path.exists(), f"Missing output file: {path}"

    from anndata import read_h5ad  # type: ignore

    integrated = read_h5ad(output_dir / "integrated.h5ad")
    assert "X_scvi" in integrated.obsm
    latent_table = pd.read_csv(output_dir / "tables" / "latent_embeddings.csv")
    assert "demo_truth" in latent_table.columns

    report_text = (output_dir / "report.md").read_text(encoding="utf-8")
    assert "Colored by: `demo_truth`" in report_text


def test_demo_batch_outputs_and_result_metadata(tmp_path: Path):
    _require_scvi_stack()
    output_dir = tmp_path / "demo_batch_output"
    result = _run_cmd(
        [
            "--demo",
            "--batch-key",
            "demo_batch",
            "--output",
            str(output_dir),
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr

    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["method"] == "scvi"
    assert payload["summary"]["batch_key"] == "demo_batch"
    assert payload["summary"]["accelerator_used"] == "cpu"
    assert payload["summary"]["latent_plot_color_by"] == "demo_truth"
    assert payload["data"]["integrated_h5ad"] == "integrated.h5ad"
    assert "umap_scvi_batch.png" in payload["data"]["figures"]
    assert "latent_embeddings.csv" in payload["data"]["tables"]
    assert "batch_mixing_metrics.csv" in payload["data"]["tables"]
    assert payload["summary"]["batch_mixing"]["n_batches"] == 2
    assert (output_dir / "figures" / "umap_scvi_batch.png").exists()
    assert (output_dir / "tables" / "latent_embeddings.csv").exists()
    assert (output_dir / "tables" / "batch_mixing_metrics.csv").exists()


def test_h5ad_input_runs_with_batch_key(tmp_path: Path):
    _require_scvi_stack()
    input_path = tmp_path / "input.h5ad"
    output_dir = tmp_path / "h5ad_output"
    _build_h5ad_input(input_path)

    result = _run_cmd(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_dir),
            "--batch-key",
            "batch",
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "integrated.h5ad").exists()
    assert (output_dir / "tables" / "latent_embeddings.csv").exists()
    assert (output_dir / "tables" / "batch_mixing_metrics.csv").exists()

    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["latent_plot_color_by"] == "scvi_latent_group"
    latent_table = pd.read_csv(output_dir / "tables" / "latent_embeddings.csv")
    assert "scvi_latent_group" in latent_table.columns


def test_layer_allows_processed_x_when_raw_counts_layer_is_selected(tmp_path: Path):
    _require_scvi_stack()
    input_path = tmp_path / "layered.h5ad"
    output_dir = tmp_path / "layer_output"
    _build_h5ad_input(input_path, with_counts_layer=True, processed_x=True)

    result = _run_cmd(
        [
            "--input",
            str(input_path),
            "--output",
            str(output_dir),
            "--layer",
            "counts",
            "--batch-key",
            "batch",
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "result.json").exists()


def test_processed_input_rejected_without_raw_layer(tmp_path: Path):
    _require_scvi_stack()
    input_path = tmp_path / "processed_like.h5ad"
    output_dir = tmp_path / "processed_output"
    _build_h5ad_input(input_path, processed_x=True)

    result = _run_cmd(["--input", str(input_path), "--output", str(output_dir)])
    assert result.returncode != 0
    assert "raw-count .h5ad or 10x single-cell input" in result.stderr
    assert "pbmc3k_processed" in result.stderr


def test_10x_directory_input_runs(tmp_path: Path):
    _require_scvi_stack()
    input_dir = _build_10x_input(tmp_path / "tenx_dir")
    output_dir = tmp_path / "tenx_output"

    result = _run_cmd(
        [
            "--input",
            str(input_dir),
            "--output",
            str(output_dir),
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr
    payload = json.loads((output_dir / "result.json").read_text())
    assert payload["summary"]["input_format"] == "10x_mtx"


def test_matrix_mtx_gz_input_runs(tmp_path: Path):
    _require_scvi_stack()
    matrix_path = _build_10x_input(tmp_path / "tenx_gz", compressed=True)
    output_dir = tmp_path / "tenx_gz_output"

    result = _run_cmd(
        [
            "--input",
            str(matrix_path),
            "--output",
            str(output_dir),
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "integrated.h5ad").exists()


def test_non_empty_output_dir_rejected(tmp_path: Path):
    _require_scvi_stack()
    output_dir = tmp_path / "occupied"
    output_dir.mkdir()
    (output_dir / "keep.txt").write_text("present", encoding="utf-8")

    result = _run_cmd(["--demo", "--output", str(output_dir)])
    assert result.returncode != 0
    assert "already exists and is not empty" in result.stderr


def test_commands_sh_quotes_output_path(tmp_path: Path):
    _require_scvi_stack()
    output_dir = tmp_path / "demo output (quoted)"
    result = _run_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr

    commands_sh = (output_dir / "reproducibility" / "commands.sh").read_text()
    assert f"--output {shlex.quote(str(output_dir))}" in commands_sh


def test_clawbio_runner_accepts_whitelisted_embedding_flags(tmp_path: Path):
    _require_scvi_stack()
    output_dir = tmp_path / "clawbio_output"
    result = _run_clawbio_cmd(
        [
            "--demo",
            "--output",
            str(output_dir),
            "--method",
            "scvi",
            "--batch-key",
            "demo_batch",
            "--latent-dim",
            "6",
            "--max-epochs",
            "2",
            "--accelerator",
            "cpu",
            "--min-genes",
            "1",
            "--min-cells",
            "1",
        ]
    )
    assert result.returncode == 0, result.stderr
    assert (output_dir / "result.json").exists()


def test_orchestrator_routes_embedding_keywords_to_new_skill():
    module = _load_orchestrator_module()
    assert module.detect_skill_from_query("Run scvi batch correction on my h5ad") == "scrna-embedding"
    assert module.detect_skill_from_query("Build a latent embedding for this single-cell dataset") == "scrna-embedding"
    assert module.detect_skill_from_query("Cluster my h5ad and find marker genes") == "scrna-orchestrator"
