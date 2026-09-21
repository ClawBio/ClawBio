"""Regression coverage for input semantics, provenance and archive safety."""

from __future__ import annotations

import io
import json
import os
import subprocess
import sys
import tarfile
from pathlib import Path

import numpy as np
import pytest

sc = pytest.importorskip("scanpy")
SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))
import spatial_transcriptomics as st


def _write_h5ad(adata, path):
    import anndata as ad

    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(path)
    return path


def test_normalized_h5ad_requires_explicit_counts_layer(tmp_path):
    adata = st.generate_demo_adata()
    original = adata.X.copy()
    adata.layers["counts"] = original.copy()
    sc.pp.normalize_total(adata, target_sum=1e4)
    sc.pp.log1p(adata)
    path = _write_h5ad(adata, tmp_path / "processed.h5ad")
    with pytest.raises(ValueError, match="counts|transformed"):
        st.load_spatial(path)
    loaded = st.load_spatial(path, counts_layer="counts")
    np.testing.assert_array_equal(loaded.X.toarray(), original.toarray())
    result = st.run_pipeline(loaded)
    np.testing.assert_array_equal(
        result["adata"].layers["counts"].toarray(), original.toarray()
    )
    assert result["count_source"] == "layers[counts]"


@pytest.mark.parametrize("bad", [np.nan, np.inf, -1.0, 0.5])
def test_invalid_counts_are_rejected_before_processing(bad):
    adata = st.generate_demo_adata()
    adata.X = adata.X.astype(float)
    adata.X[0, 0] = bad
    with pytest.raises(ValueError, match="counts"):
        st.run_pipeline(adata)


@pytest.mark.parametrize("coords", [np.zeros((64, 1)), np.full((64, 2), np.nan)])
def test_bad_coordinates_rejected_at_load(tmp_path, coords):
    adata = st.generate_demo_adata()
    adata.obsm["spatial"] = coords
    path = _write_h5ad(adata, tmp_path / "bad-spatial.h5ad")
    with pytest.raises(ValueError, match="spatial|coordinates"):
        st.load_spatial(path)


def test_qc_exports_mito_metrics_and_applies_requested_threshold():
    adata = st.generate_demo_adata()
    names = list(adata.var_names)
    names[-1] = "MT-CO3"
    adata.var_names = names
    matrix = adata.X.tolil()
    matrix[0, -1] = 10000
    adata.X = matrix.tocsr()
    result = st.run_pipeline(
        adata, max_pct_mt=50.0, n_pcs=6, n_neighbors=9, nhood_perms=25
    )
    assert result["adata"].n_obs == 63
    assert result["qc_summary"]["mt_gene_count"] == 1
    assert not result["qc_table"].iloc[0]["pass_qc"]
    assert result["parameters"]["n_pcs"] == 6
    assert result["parameters"]["n_neighbors"] == 9
    assert result["parameters"]["nhood_perms"] == 25


def test_absent_mito_gene_symbols_are_unavailable_not_zero_percent():
    result = st.run_pipeline(st.generate_demo_adata())
    assert result["qc_summary"]["mt_gene_count"] == 0
    assert result["qc_summary"]["median_pct_counts_mt"] is None
    assert result["qc_table"]["pct_counts_mt"].isna().all()


@pytest.mark.parametrize(
    "kind", [tarfile.SYMTYPE, tarfile.LNKTYPE, tarfile.FIFOTYPE, tarfile.CHRTYPE]
)
def test_archive_rejects_nonregular_members(tmp_path, kind):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as archive:
        member = tarfile.TarInfo("spatial/tissue_positions.csv")
        member.type = kind
        member.linkname = "outside"
        archive.addfile(member)
    buf.seek(0)
    with (
        tarfile.open(fileobj=buf) as archive,
        pytest.raises(OSError, match="Refusing|regular|unsafe"),
    ):
        st._safe_extract_visium_tar(archive, tmp_path / "outs")


def test_archive_writer_does_not_use_extractall(tmp_path, monkeypatch):
    buf = io.BytesIO()
    content = b"barcode,in_tissue\nspot,1\n"
    with tarfile.open(fileobj=buf, mode="w") as archive:
        member = tarfile.TarInfo("spatial/tissue_positions.csv")
        member.size = len(content)
        archive.addfile(member, io.BytesIO(content))
    buf.seek(0)
    with tarfile.open(fileobj=buf) as archive:

        def forbidden(*args, **kwargs):
            pytest.fail("downloaded archives must use bounded regular-file extraction")

        monkeypatch.setattr(archive, "extractall", forbidden)
        st._safe_extract_visium_tar(archive, tmp_path / "outs")
    assert (tmp_path / "outs/spatial/tissue_positions.csv").read_bytes() == content


@pytest.mark.parametrize(
    "name",
    [
        "/spatial/tissue_positions.csv",
        "spatial/../../evil",
        "spatial\\..\\evil",
        "unrelated/evil",
    ],
)
def test_archive_rejects_paths_outside_data_allowlist(tmp_path, name):
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as archive:
        archive.addfile(tarfile.TarInfo(name))
    buf.seek(0)
    with tarfile.open(fileobj=buf) as archive, pytest.raises(OSError):
        st._safe_extract_visium_tar(archive, tmp_path / "outs")


def test_archive_rejects_existing_destination_symlink(tmp_path):
    outside = tmp_path / "outside"
    outside.mkdir()
    dest = tmp_path / "outs"
    dest.mkdir()
    (dest / "spatial").symlink_to(outside, target_is_directory=True)
    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w") as archive:
        archive.addfile(tarfile.TarInfo("spatial/tissue_positions.csv"))
    buf.seek(0)
    with tarfile.open(fileobj=buf) as archive, pytest.raises(OSError):
        st._safe_extract_visium_tar(archive, dest)
    assert not list(outside.iterdir())


def test_public_download_hash_failure_does_not_publish_cache(tmp_path, monkeypatch):
    import urllib.request

    def truncated_download(url, path):
        Path(path).write_bytes(b"partial archive")

    monkeypatch.setattr(urllib.request, "urlretrieve", truncated_download)
    dest = tmp_path / "cache/outs"
    with pytest.raises(OSError, match="SHA-256"):
        st.ensure_public_visium_outs(dest)
    assert not dest.exists()
    assert not list((dest.parent / "tarballs").glob("*.tar.gz"))


def test_incomplete_cache_is_not_ready(tmp_path):
    (tmp_path / "filtered_feature_bc_matrix").mkdir()
    (tmp_path / "spatial").mkdir()
    (tmp_path / "filtered_feature_bc_matrix/matrix.mtx.gz").write_bytes(b"matrix")
    (tmp_path / "spatial/tissue_positions_list.csv").write_bytes(b"positions")
    assert not st._visium_outs_ready(tmp_path)


def test_cli_refuses_overwrite_before_computation(tmp_path):
    existing = tmp_path / "report.md"
    existing.write_text("preserve this report", encoding="utf-8")
    proc = subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "spatial_transcriptomics.py"),
            "--demo",
            "--output",
            str(tmp_path),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert "--overwrite" in proc.stderr
    assert existing.read_text() == "preserve this report"


@pytest.mark.parametrize("source_kind", ["demo", "h5ad_counts"])
def test_nondefault_cli_bundle_replays_numeric_results(tmp_path, source_kind):
    output = tmp_path / "original"
    input_path = None
    source_flags = ["--demo"]
    if source_kind == "h5ad_counts":
        adata = st.generate_demo_adata(seed=11)
        adata.layers["raw counts"] = adata.X.copy()
        sc.pp.normalize_total(adata, target_sum=1e4)
        sc.pp.log1p(adata)
        input_path = _write_h5ad(adata, tmp_path / "processed sample.h5ad")
        source_flags = ["--input", str(input_path), "--counts-layer", "raw counts"]
    flags = [
        "--min-genes",
        "6",
        "--min-cells",
        "2",
        "--n-top-hvg",
        "20",
        "--n-pcs",
        "5",
        "--n-neighbors",
        "9",
        "--leiden-resolution",
        "0.7",
        "--random-state",
        "11",
        "--nhood-perms",
        "25",
        "--max-pct-mt",
        "95",
    ]
    proc = subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "spatial_transcriptomics.py"),
            *source_flags,
            *flags,
            "--output",
            str(output),
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    before = json.loads((output / "result.json").read_text())
    commands = (output / "reproducibility/commands.sh").read_text()
    for flag in flags:
        assert flag in commands
    assert "--overwrite" not in commands
    manifest = json.loads((output / "reproducibility/run_manifest.json").read_text())
    assert manifest["parameters"] == before["parameters"]
    assert manifest["environment"]["python"].startswith(
        f"{sys.version_info.major}.{sys.version_info.minor}"
    )
    assert manifest["input"]["kind"] == (
        "synthetic_demo" if source_kind == "demo" else "measured_input"
    )
    if input_path is not None:
        assert "--counts-layer" in commands and "raw counts" in commands
        assert before["count_source"] == "layers[raw counts]"
    for name in ["report.md", "result.json", "reproducibility/run_manifest.json"]:
        assert name in (output / "reproducibility/checksums.sha256").read_text()
    env = (output / "reproducibility/environment.yml").read_text()
    assert "scanpy==" in env and "numpy==" in env
    replay = subprocess.run(
        ["bash", str(output / "reproducibility/commands.sh")],
        cwd=tmp_path,
        env={
            **os.environ,
            "PYTHON": sys.executable,
            "REPLAY_OUTPUT": str(tmp_path / "replay"),
            "INPUT_PATH": str(input_path) if input_path is not None else "",
        },
        capture_output=True,
        text=True,
        check=False,
    )
    assert replay.returncode == 0, replay.stderr
    after = json.loads((tmp_path / "replay/result.json").read_text())
    for key in [
        "parameters",
        "n_spots",
        "n_clusters",
        "top_moran",
        "markers",
        "nhood_zscore",
        "co_occurrence",
        "count_source",
        "input_provenance",
    ]:
        assert after[key] == before[key], key


def test_replay_rejects_modified_measured_input_before_analysis(tmp_path):
    outs = st.write_visium_outs(st.generate_demo_adata(), tmp_path / "outs")
    identity = st._input_identity(outs, seed=7)
    assert len(identity["files"]) == 4
    positions = outs / "spatial/tissue_positions.csv"
    positions.write_text(positions.read_text() + "changed-input\n")
    output = tmp_path / "output"
    proc = subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "spatial_transcriptomics.py"),
            "--input",
            str(outs),
            "--output",
            str(output),
            "--expected-input-sha256",
            identity["sha256"],
        ],
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode != 0
    assert "SHA-256 differs" in proc.stderr
    assert not output.exists()
