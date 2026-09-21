"""Tests for the spatial-transcriptomics skill."""

from __future__ import annotations

import json
import re
import subprocess
import sys
from pathlib import Path

import numpy as np
import pytest

scanpy = pytest.importorskip("scanpy")

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
sys.path.insert(0, str(SKILL_DIR))
sys.path.insert(0, str(PROJECT_ROOT))

import spatial_transcriptomics as st  # noqa: E402

SCRIPT = SKILL_DIR / "spatial_transcriptomics.py"


def run_cli(args: list[str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )


def _parse_output_contract(skill_md: Path) -> list[str]:
    if not skill_md.exists():
        return []
    match = re.search(
        r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```",
        skill_md.read_text(encoding="utf-8"),
        re.S,
    )
    if not match:
        return []
    files: list[str] = []
    parents: dict[int, str] = {}
    for raw in match.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for extra in [k for k in parents if k > depth]:
                del parents[extra]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(rel + "/" + name if rel else name)
    return files


def test_demo_adata_has_two_spatial_domains():
    adata = st.generate_demo_adata()
    assert adata.n_obs == 64
    assert "spatial" in adata.obsm
    assert set(adata.obs["domain"]) == {"epithelium", "stroma"}
    assert adata.obs["domain"].value_counts().min() == 32


def test_moran_i_is_high_for_domain_markers_and_low_for_noise():
    adata = st.generate_demo_adata()
    knn = st.knn_indices(np.asarray(adata.obsm["spatial"]), 6)
    epcam = adata[:, "EPCAM"].X
    noise = adata[:, "GENE20"].X
    if hasattr(epcam, "toarray"):
        epcam = epcam.toarray()
        noise = noise.toarray()
    i_marker = st.moran_i(np.asarray(epcam).reshape(-1), knn)
    i_noise = st.moran_i(np.asarray(noise).reshape(-1), knn)
    assert i_marker > 0.4
    assert i_marker > i_noise + 0.2


def test_nhood_diagonal_is_enriched_for_two_blocks():
    adata = st.generate_demo_adata()
    knn = st.knn_indices(np.asarray(adata.obsm["spatial"]), 6)
    clusters, zscore, _counts = st.nhood_enrichment(
        knn, adata.obs["domain"].to_numpy(), n_perms=40, seed=7
    )
    assert set(clusters) == {"epithelium", "stroma"}
    for i, name in enumerate(clusters):
        assert zscore[i, i] > 2.0, f"{name} should be spatially self-enriched"


def test_h5ad_without_spatial_is_rejected(tmp_path):
    import anndata as ad

    adata = st.generate_demo_adata()
    del adata.obsm["spatial"]
    path = tmp_path / "no_spatial.h5ad"
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(path)
    with pytest.raises(ValueError, match="obsm\\['spatial'\\]"):
        st.load_spatial(path)


def test_visium_outs_roundtrip(tmp_path):
    adata = st.generate_demo_adata()
    outs = st.write_visium_outs(adata, tmp_path / "outs")
    loaded = st.load_spatial(outs)
    assert loaded.n_obs == 64
    assert "spatial" in loaded.obsm
    assert loaded.obsm["spatial"].shape == (64, 2)


def test_too_few_spots_abstain():
    adata = st.generate_demo_adata()
    tiny = adata[:3].copy()
    with pytest.raises(st.InsufficientSpotsError):
        st.run_pipeline(tiny)


def test_pipeline_recovers_two_spatial_programmes():
    adata = st.generate_demo_adata()
    result = st.run_pipeline(adata, random_state=7)
    assert result["adata"].obs["leiden"].nunique() >= 2
    moran = {row["gene"]: row["moran_i"] for row in result["moran"]}
    assert moran["EPCAM"] > moran["GENE20"]
    assert moran["COL1A1"] > moran["GENE20"]
    marker_genes = {row["gene"] for row in result["markers"]}
    assert marker_genes & set(st.MARKER_LEFT + st.MARKER_RIGHT)


def test_demo_cli_writes_contract_and_disclaimer(tmp_path):
    output = tmp_path / "out"
    proc = run_cli(["--demo", "--output", str(output)])
    assert proc.returncode == 0, proc.stderr
    report = (output / "report.md").read_text(encoding="utf-8")
    assert "research and educational tool" in report
    assert "not a medical device" in report.lower()
    payload = json.loads((output / "result.json").read_text(encoding="utf-8"))
    assert payload["demo"] is True
    assert payload["n_spots"] == 64
    assert payload["n_clusters"] >= 2
    assert (output / "figures" / "spatial_leiden.png").exists()
    assert (output / "tables" / "moran_i.csv").exists()
    commands = (output / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--demo" in commands
    assert "CLAWBIO_ROOT" in commands
    env = (output / "reproducibility" / "environment.yml").read_text(encoding="utf-8")
    assert "scanpy" in env
    checksums = (output / "reproducibility" / "checksums.sha256").read_text(encoding="utf-8")
    assert "moran_i.csv" in checksums


def test_bundle_delegates_to_shared_helpers():
    assert st.write_checksums.__module__ == "clawbio.common.reproducibility"
    assert st.write_environment_yml.__module__ == "clawbio.common.reproducibility"
    assert st.write_portable_commands_sh.__module__ == "clawbio.common.reproducibility"


def test_visium_cli_path(tmp_path):
    outs = st.write_visium_outs(st.generate_demo_adata(), tmp_path / "outs")
    output = tmp_path / "run"
    proc = run_cli(["--input", str(outs), "--output", str(output)])
    assert proc.returncode == 0, proc.stderr
    payload = json.loads((output / "result.json").read_text(encoding="utf-8"))
    assert payload["demo"] is False
    assert payload["n_spots"] == 64


def test_safe_extract_rejects_path_traversal(tmp_path):
    import io
    import tarfile

    buf = io.BytesIO()
    with tarfile.open(fileobj=buf, mode="w:gz") as tar:
        info = tarfile.TarInfo(name="../evil.txt")
        info.size = 0
        tar.addfile(info)
    buf.seek(0)
    with tarfile.open(fileobj=buf, mode="r:gz") as tar:
        with pytest.raises(OSError, match="unsafe|Unexpected"):
            st._safe_extract_visium_tar(tar, tmp_path)


def test_public_visium_load_and_run_writes_numeric_spatial_stats(tmp_path):
    """Drive the shipped loader and pipeline on the downloaded 10x lymph-node outs."""
    outs = None
    try:
        outs = st.ensure_public_visium_outs()
    except OSError as exc:
        pytest.skip(f"public Visium download failed: {exc}")
    assert outs is not None
    assert (outs / "filtered_feature_bc_matrix" / "matrix.mtx.gz").is_file()
    assert (outs / "spatial" / "tissue_positions_list.csv").is_file()

    adata = st.load_spatial(outs)
    assert "spatial" in adata.obsm
    assert adata.n_obs >= 1000
    assert adata.obsm["spatial"].shape == (adata.n_obs, 2)

    # Keep the downloaded Visium counts/coords; subset spots only for runtime.
    subset = adata[:400].copy()
    result = st.run_pipeline(subset, n_top_hvg=80, random_state=7)
    st.generate_report(
        result, tmp_path, source_label=str(outs), demo=False
    )
    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert report.splitlines()[0] == "# Spatial Transcriptomics Report"
    assert "V1_Human_Lymph_Node" in report
    assert "**Spots**: 400" in report
    assert "synthetic 8x8" not in report.lower()
    moran_path = tmp_path / "tables" / "moran_i.csv"
    nhood_path = tmp_path / "tables" / "nhood_enrichment.csv"
    assert moran_path.is_file()
    assert nhood_path.is_file()
    moran_lines = moran_path.read_text(encoding="utf-8").strip().splitlines()
    assert len(moran_lines) >= 10
    moran_values = []
    for line in moran_lines[1:]:
        value = line.rsplit(",", 1)[-1]
        moran_values.append(float(value))
    assert any(abs(v) > 0.0 for v in moran_values)
    nhood_text = nhood_path.read_text(encoding="utf-8")
    assert any(ch.isdigit() for ch in nhood_text)


class TestOutputContract:
    def test_documented_outputs_are_produced(self, tmp_path):
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        if not promised:
            pytest.skip("No parseable Output Structure section")
        proc = run_cli(["--demo", "--output", str(tmp_path)])
        assert proc.returncode == 0, proc.stderr
        missing = [path for path in promised if not (tmp_path / path).exists()]
        assert not missing, "SKILL.md promises missing artifacts: " + ", ".join(missing)
