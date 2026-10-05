"""End-to-end tests for the scrna-embedding-audit CLI."""

from __future__ import annotations

import csv
import json
import os
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

import scrna_embedding_audit as sea

SCRIPT = SKILL_DIR / "scrna_embedding_audit.py"
DISCLAIMER_START = "ClawBio is a research and educational tool."
FAST = ["--n-resamples", "10"]  # the minimum the skill accepts; keeps the suite quick
BANNED_PHRASES = (
    "outperform",
    "state of the art",
    "state-of-the-art",
    "superior",
    "beats",
    "wins",
    "significant",
    "best embedding",
)


def run_cli(
    args: list[str], env: dict | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
        env=env,
    )


def _write_h5ad(adata, path: Path) -> Path:
    import anndata as ad
    import pandas as pd

    # With pyarrow installed, pandas 3 backs str data with ArrowStringArray,
    # which anndata 0.12 cannot write to h5ad. Object strings always can;
    # categorical categories need the same cast.
    adata = adata.copy()
    for frame in (adata.obs, adata.var):
        for column in frame.columns:
            if isinstance(frame[column].dtype, pd.CategoricalDtype):
                cat = frame[column].cat
                frame[column] = cat.set_categories(cat.categories.astype(object))
            elif (
                isinstance(frame[column].dtype, pd.StringDtype)
                or frame[column].dtype == "str"
            ):
                frame[column] = frame[column].astype(object)
    adata.obs_names = adata.obs_names.astype(object)
    adata.var_names = adata.var_names.astype(object)
    ad.settings.allow_write_nullable_strings = True
    adata.write_h5ad(path)
    return path


def _parse_output_contract(skill_md: Path) -> list[str]:
    if not skill_md.exists():
        return []
    match = re.search(
        r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```",
        skill_md.read_text(encoding="utf-8"),
        re.DOTALL,
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


def _result(path: Path) -> dict:
    return json.loads((path / "result.json").read_text())


@pytest.fixture(scope="module")
def demo_run(tmp_path_factory):
    """One full demo run (default 50 half-samples) shared by the read-only tests."""
    out = tmp_path_factory.mktemp("demo")
    proc = run_cli(["--demo", "--output", str(out)])
    assert proc.returncode == 0, proc.stderr
    return out, proc


# --- demo data -----------------------------------------------------------------


def test_demo_adata_matches_recipe():
    spec = sea.load_demo_spec()
    adata = sea.generate_demo_adata()
    assert adata.n_obs == spec["n_cells"] == 660 and adata.n_vars == spec["n_genes"]
    assert sorted(set(adata.obs["cell_type"])) == sorted(spec["cell_types"])
    assert set(adata.obs["donor"]) == set(spec["batches"])
    assert set(spec["embedding_keys"]) <= set(adata.obsm)
    single = spec["single_batch_cell_type"]
    nk = adata.obs["cell_type"] == single["name"]
    assert int(nk.sum()) == single["n_cells"] and set(adata.obs.loc[nk, "donor"]) == {
        single["batch"]
    }
    counts = np.asarray(adata.layers["counts"])
    assert np.all(counts >= 0) and np.all(counts == np.round(counts))
    assert np.array_equal(counts, np.asarray(adata.X))
    # pandas 3 stores these as its str dtype; what matters is that no value is missing
    assert (
        not adata.obs["cell_type"].isna().any() and not adata.obs["donor"].isna().any()
    )


def test_demo_embeddings_derive_from_counts_not_labels():
    # the recipe says the embeddings are PCA projections of transformed counts;
    # a label-free check: re-running with another seed changes every embedding
    a = sea.generate_demo_adata(seed=0)
    b = sea.generate_demo_adata(seed=1)
    for key in sea.load_demo_spec()["embedding_keys"]:
        assert a.obsm[key].shape == (660, 10)
        assert not np.allclose(a.obsm[key], b.obsm[key])


def test_demo_generation_is_deterministic():
    a = sea.generate_demo_adata(seed=0)
    b = sea.generate_demo_adata(seed=0)
    assert np.array_equal(
        np.asarray(a.layers["counts"]), np.asarray(b.layers["counts"])
    )
    for key in sea.load_demo_spec()["embedding_keys"]:
        assert np.array_equal(a.obsm[key], b.obsm[key])


# --- demo CLI ------------------------------------------------------------------


class TestOutputContract:
    def test_documented_outputs_are_produced(self, demo_run):
        out, _ = demo_run
        promised = _parse_output_contract(SKILL_DIR / "SKILL.md")
        assert promised, "SKILL.md has no parseable Output Structure section"
        missing = [path for path in promised if not (out / path).exists()]
        assert not missing, "SKILL.md promises missing artifacts: " + ", ".join(missing)

    def test_every_documented_output_is_checksummed(self, demo_run):
        out, _ = demo_run
        checksums = (out / "reproducibility/checksums.sha256").read_text()
        for path in _parse_output_contract(SKILL_DIR / "SKILL.md"):
            if path.startswith("reproducibility/checksums"):
                continue
            assert path in checksums, path


def test_demo_result_json_schema_and_provenance(demo_run):
    out, _ = demo_run
    result = _result(out)
    assert result["skill"] == "scrna-embedding-audit"
    assert result["demo"] is True
    assert result["input_provenance"]["kind"] == "synthetic_demo"
    assert (
        result["n_cells"] == 660
        and result["cells_scored"] == 660
        and result["subsampled"] is False
    )
    assert result["labels"] == {
        "key": "cell_type",
        "n_unique": 4,
        "counts": {"B cell": 200, "Monocyte": 200, "NK cell": 60, "T cell": 200},
    }
    assert result["batch"]["key"] == "donor" and result["batch"]["n_unique"] == 2
    assert result["batch"]["fraction_cells_in_single_batch_labels"] == pytest.approx(
        60 / 660
    )
    assert (
        result["baseline"]["name"] == "pca"
        and result["baseline"]["provenance"] == sea.PROVENANCE_BASELINE
    )
    assert "not scib-metrics' X_pre" in result["baseline"]["recipe"]
    assert result["embedding_keys"] == ["X_corrected", "X_batchy", "X_overcorrected"]
    assert set(result["embeddings"]) == set(result["embedding_keys"]) | {
        sea.BASELINE_KEY,
        sea.CONTROL_KEY,
    }
    for name, entry in result["embeddings"].items():
        assert set(entry["metrics"]) == set(sea.METRIC_NAMES)
        assert set(entry["reasons"]).isdisjoint(
            {k for k, v in entry["metrics"].items() if v is not None}
        )
        expected = {
            sea.BASELINE_KEY: sea.PROVENANCE_BASELINE,
            sea.CONTROL_KEY: sea.PROVENANCE_CONTROL,
        }.get(name, sea.PROVENANCE_USER)
        assert entry["provenance"] == expected
    assert set(result["verdicts"]) == set(
        result["embedding_keys"]
    )  # baseline and control get no verdict
    for entry in result["verdicts"].values():
        assert set(entry["per_metric"]) == set(sea.METRIC_NAMES)
        for c in entry["counts_by_group"].values():
            assert c["n_metrics"] == sum(v for k, v in c.items() if k != "n_metrics")
        for item in entry["per_metric"].values():
            assert item["verdict"] in {
                "higher",
                "lower",
                "not separable",
                "tied",
                "not evaluated",
            }
            assert item["n_resamples"] == result["parameters"]["n_resamples"] == 50
    assert result["references"]["scib_metrics_fixture_version"] == "0.6.1"
    assert set(result["metric_definitions"]) == set(sea.METRIC_NAMES)
    assert set(result["metric_definitions_omitted"]) == {
        "kbet_per_label",
        "pcr_comparison",
        "nmi_ari_cluster_labels_leiden",
    }
    assert (
        result["effective"]["lisi_neighbors"] == 90
        and result["effective"]["perplexity"] == 30
    )
    assert result["effective"]["resample_sizes"] == [330]
    assert (
        result["parameters"]["labels_key"] is None
    )  # demo: not a user flag, so not replayed
    assert (
        "resamples" not in result
    )  # per-resample values live in tables/resamples.csv only
    assert result["timings_seconds"]["total"] > 0


def test_demo_contrast_with_margins(demo_run):
    """Verdicts the recipe is built to produce, asserted with margins, plus their preconditions."""
    out, _ = demo_run
    result = _result(out)
    emb = result["embeddings"]
    base = emb[sea.BASELINE_KEY]["metrics"]
    # preconditions: the baseline carries the batch shift and X_batchy carries more of it
    assert base["ilisi"] - emb["X_batchy"]["metrics"]["ilisi"] > 0.1
    assert emb["X_corrected"]["metrics"]["ilisi"] - base["ilisi"] > 0.1
    verdicts = {
        k: {m: v["verdict"] for m, v in result["verdicts"][k]["per_metric"].items()}
        for k in result["embedding_keys"]
    }
    assert verdicts["X_batchy"]["ilisi"] == "lower"
    assert verdicts["X_corrected"]["ilisi"] == "higher"
    assert result["verdicts"]["X_batchy"]["per_metric"]["ilisi"]["frac_positive"] == 0.0
    assert (
        result["verdicts"]["X_corrected"]["per_metric"]["ilisi"]["frac_positive"] == 1.0
    )
    # removing every marker gene must cost biological conservation
    bio_lower = sum(
        verdicts["X_overcorrected"][m] == "lower" for m in sea.metrics_mod.BIO_METRICS
    )
    assert bio_lower >= 3
    # the random control mixes batches by construction and carries no biology
    ctrl = emb[sea.CONTROL_KEY]["metrics"]
    assert ctrl["ilisi"] > 0.7 and ctrl["clisi"] < 0.5 and ctrl["nmi_kmeans"] < 0.1


def test_demo_tables_have_provenance_and_match_json(demo_run):
    out, _ = demo_run
    result = _result(out)
    with (out / "tables/metrics.csv").open() as fh:
        rows = list(csv.DictReader(fh))
    assert {"embedding", "metric", "group", "value", "reason", "provenance"} <= set(
        rows[0]
    )
    assert len(rows) == len(result["embeddings"]) * len(sea.METRIC_NAMES)
    for row in rows:
        expected = result["embeddings"][row["embedding"]]["metrics"][row["metric"]]
        assert row["provenance"]
        if row["value"] == "":
            assert expected is None
        else:
            assert float(row["value"]) == pytest.approx(expected)
    with (out / "tables/verdicts.csv").open() as fh:
        verdict_rows = list(csv.DictReader(fh))
    assert {
        "embedding",
        "metric",
        "group",
        "verdict",
        "full_value",
        "baseline_value",
        "median_diff",
        "sd_diff",
        "frac_positive",
        "n_resamples",
        "provenance",
    } <= set(verdict_rows[0])
    assert len(verdict_rows) == len(result["embedding_keys"]) * len(sea.METRIC_NAMES)
    for row in verdict_rows:
        assert (
            row["verdict"]
            == result["verdicts"][row["embedding"]]["per_metric"][row["metric"]][
                "verdict"
            ]
        )
    with (out / "tables/resamples.csv").open() as fh:
        resample_rows = list(csv.DictReader(fh))
    assert {"resample", "embedding", "metric", "value", "provenance"} <= set(
        resample_rows[0]
    )
    assert (
        len({r["resample"] for r in resample_rows})
        == result["parameters"]["n_resamples"]
    )
    assert len(resample_rows) == 50 * len(result["embeddings"]) * len(sea.METRIC_NAMES)


def test_demo_report_wording(demo_run):
    out, _ = demo_run
    report = (out / "report.md").read_text()
    assert DISCLAIMER_START in report
    assert "synthetic" in report.lower()
    assert "scib-metrics 0.6.1" in report
    assert re.search(
        r"biological conservation: higher than the PCA baseline on \d of 5 metrics",
        report,
    )
    assert re.search(
        r"batch mixing: higher than the PCA baseline on \d of 4 metrics", report
    )
    assert re.search(
        r"cross-batch label transfer: higher than the PCA baseline on \d of 1 metrics",
        report,
    )
    assert sea.report_mod.SCOPE_SENTENCE in report
    assert sea.report_mod.BASELINE_SENTENCE in report
    assert sea.report_mod.CIRCULARITY_SENTENCE in report
    assert "this skill's own metric" in report
    assert "random control" in report.lower() or "control:random" in report
    lowered = report.lower()
    for banned in BANNED_PHRASES:
        assert banned not in lowered, banned


def test_demo_reproducibility_bundle(demo_run):
    out, _ = demo_run
    commands = (out / "reproducibility/commands.sh").read_text()
    assert "--demo" in commands and "CLAWBIO_ROOT" in commands
    for flag in ("--n-resamples", "--random-state", "--lisi-neighbors"):
        assert flag in commands
    assert "--effective" not in commands and "--labels-key" not in commands
    env_yml = (out / "reproducibility/environment.yml").read_text()
    assert (
        "scanpy==" in env_yml and "scikit-learn==" in env_yml and "numpy==" in env_yml
    )
    checksums = (out / "reproducibility/checksums.sha256").read_text()
    assert (
        "result.json" in checksums
        and "tables/metrics.csv" in checksums
        and "run_manifest.json" in checksums
    )
    manifest = json.loads((out / "reproducibility/run_manifest.json").read_text())
    assert manifest["parameters"] == _result(out)["parameters"]
    assert any(
        f["path"].endswith("scib_metrics_0.6.1_expected.json")
        for f in manifest["source_files"]
    )
    assert sea.write_checksums.__module__ == "clawbio.common.reproducibility"
    assert sea.write_environment_yml.__module__ == "clawbio.common.reproducibility"
    assert sea.write_portable_commands_sh.__module__ == "clawbio.common.reproducibility"


def test_checksums_verify_from_output_dir(demo_run):
    out, _ = demo_run
    proc = subprocess.run(
        ["sha256sum", "-c", "--quiet", "reproducibility/checksums.sha256"],
        cwd=out,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr


def _numbers_equal(a, b, path=""):
    if isinstance(a, dict):
        assert set(a) == set(b), path
        for k in a:
            _numbers_equal(a[k], b[k], f"{path}/{k}")
    elif isinstance(a, list):
        assert len(a) == len(b), path
        for i, (x, y) in enumerate(zip(a, b)):
            _numbers_equal(x, y, f"{path}[{i}]")
    elif isinstance(a, float) and isinstance(b, float):
        assert a == pytest.approx(b, abs=1e-9), path
    else:
        assert a == b, path


def test_replay_from_commands_sh_reproduces_numbers(tmp_path):
    out = tmp_path / "run"
    proc = run_cli(["--demo", "--output", str(out), *FAST])
    assert proc.returncode == 0, proc.stderr
    replay = tmp_path / "replay"
    env = dict(
        os.environ,
        PYTHON=sys.executable,
        REPLAY_OUTPUT=str(replay),
        CLAWBIO_ROOT=str(PROJECT_ROOT),
    )
    proc = subprocess.run(
        ["bash", str(out / "reproducibility/commands.sh")],
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    a, b = _result(out), _result(replay)
    _numbers_equal(a["embeddings"], b["embeddings"])
    _numbers_equal(a["verdicts"], b["verdicts"])


def test_demo_runs_offline(tmp_path):
    """Block every INET socket in the subprocess; the demo must still complete."""
    site = tmp_path / "site"
    site.mkdir()
    (site / "sitecustomize.py").write_text(
        "import socket, sys\n"
        "_real = socket.socket\n"
        "class _Guard(_real):\n"
        "    def __init__(self, family=-1, *a, **k):\n"
        "        if family in (socket.AF_INET, socket.AF_INET6):\n"
        "            raise RuntimeError('network access attempted')\n"
        "        super().__init__(family, *a, **k)\n"
        "socket.socket = _Guard\n"
        "def _blocked(*a, **k):\n"
        "    raise RuntimeError('network access attempted')\n"
        "socket.create_connection = _blocked\n"
        "socket.getaddrinfo = _blocked\n"
        "print('OFFLINE-GUARD-ACTIVE', file=sys.stderr)\n"
    )
    env = {k: v for k, v in os.environ.items() if k != "CLAWBIO_OTLP_ENDPOINT"}
    env["PYTHONPATH"] = str(site) + os.pathsep + env.get("PYTHONPATH", "")
    env["MPLCONFIGDIR"] = str(tmp_path / "mpl")
    env["NUMBA_CACHE_DIR"] = str(tmp_path / "numba")
    proc = run_cli(["--demo", "--output", str(tmp_path / "out"), *FAST], env=env)
    assert "OFFLINE-GUARD-ACTIVE" in proc.stderr
    assert proc.returncode == 0, proc.stderr
    assert "network access attempted" not in proc.stderr


def test_demo_through_clawbio_runner(tmp_path):
    proc = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "clawbio.py"),
            "run",
            "scrna-embedding-audit",
            "--demo",
            "--output",
            str(tmp_path / "out"),
            *FAST,
        ],
        cwd=PROJECT_ROOT,
        capture_output=True,
        text=True,
        check=False,
    )
    assert proc.returncode == 0, proc.stdout + proc.stderr
    assert (tmp_path / "out/result.json").exists()


# --- user input path -----------------------------------------------------------


def _user_h5ad(tmp_path: Path, seed: int = 1) -> Path:
    return _write_h5ad(sea.generate_demo_adata(seed=seed), tmp_path / "cells.h5ad")


def _user_args(path: Path, out: Path, *extra: str) -> list[str]:
    return [
        "--input",
        str(path),
        "--output",
        str(out),
        "--labels-key",
        "cell_type",
        *FAST,
        *extra,
    ]


def test_user_h5ad_with_str_and_categorical_obs_is_scored(tmp_path):
    import anndata as ad
    import pandas as pd

    adata = sea.generate_demo_adata(seed=2)
    adata.obs["cell_type"] = pd.Categorical(adata.obs["cell_type"].to_numpy())
    adata.obs["donor"] = pd.array(adata.obs["donor"].to_numpy(), dtype="string")
    path = _write_h5ad(adata, tmp_path / "typed.h5ad")
    reread = ad.read_h5ad(path)
    assert isinstance(reread.obs["cell_type"].dtype, pd.CategoricalDtype)
    assert (
        reread.obs["donor"].dtype != object
    )  # pandas 3 reads strings back as a str dtype
    proc = run_cli(
        _user_args(
            path,
            tmp_path / "out",
            "--batch-key",
            "donor",
            "--embeddings",
            "X_corrected",
        )
    )
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    assert (
        result["input_provenance"]["kind"] == "measured_input"
        and result["input_provenance"]["sha256"]
    )
    assert result["embedding_keys"] == ["X_corrected"]
    assert result["labels"]["n_unique"] == 4 and result["batch"]["n_unique"] == 2


def test_h5ad_with_uns_and_var_columns_is_scored(tmp_path):
    # files from other tools carry uns dicts and var annotations; neither may break loading
    adata = sea.generate_demo_adata(seed=8)
    adata.uns["prepared_by"] = {"script": "x.py", "n_pcs": 10, "nested": {"a": [1, 2]}}
    adata.var["feature_name"] = np.array(
        [f"GENE{i}" for i in range(adata.n_vars)], dtype=object
    )
    adata.var["highly_variable"] = np.zeros(
        adata.n_vars, dtype=bool
    )  # must not leak into the baseline
    path = _write_h5ad(adata, tmp_path / "annotated.h5ad")
    proc = run_cli(
        _user_args(
            path,
            tmp_path / "out",
            "--batch-key",
            "donor",
            "--embeddings",
            "X_corrected",
        )
    )
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    assert (
        result["baseline"]["n_hvg_used"] == 300
    )  # all genes, not the file's empty HVG flag


def test_default_embedding_selection_skips_layouts_and_warns(tmp_path):
    adata = sea.generate_demo_adata(seed=3)
    adata.obsm["X_umap"] = adata.obsm["X_corrected"][:, :2].copy()
    adata.obsm["X_draw_graph_fa"] = adata.obsm["X_corrected"][:, :3].copy()
    adata.obsm["X_tiny"] = adata.obsm["X_corrected"][:, :2].copy()
    adata.obsm["spatial"] = adata.obsm["X_corrected"][:, :2].copy()
    path = _write_h5ad(adata, tmp_path / "layouts.h5ad")
    proc = run_cli(_user_args(path, tmp_path / "out", "--batch-key", "donor"))
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    assert set(result["embedding_keys"]) == {
        "X_corrected",
        "X_batchy",
        "X_overcorrected",
    }
    skipped = [w for w in result["warnings"] if "skipped by default" in w]
    assert (
        skipped
        and "X_umap" in skipped[0]
        and "X_draw_graph_fa" in skipped[0]
        and "X_tiny" in skipped[0]
    )


def test_counts_layer_is_auto_detected_and_no_batch_is_reported(tmp_path):
    path = _user_h5ad(tmp_path)
    proc = run_cli(_user_args(path, tmp_path / "out", "--embeddings", "X_corrected"))
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    assert result["count_source"] == "layers['counts']"
    assert result["batch"] == {
        "key": None,
        "n_unique": 0,
        "counts": {},
        "fraction_cells_in_single_batch_labels": None,
    }
    entry = result["embeddings"]["X_corrected"]
    for name in (
        "silhouette_batch",
        "bras",
        "ilisi",
        "isolated_labels",
        "knn_transfer_accuracy",
    ):
        assert entry["metrics"][name] is None and "batch" in entry["reasons"][name]
    assert entry["metrics"]["graph_connectivity"] is not None
    report = (tmp_path / "out/report.md").read_text()
    assert sea.report_mod.NO_BATCH_SENTENCE in report
    assert (
        result["verdicts"]["X_corrected"]["per_metric"]["ilisi"]["verdict"]
        == "not evaluated"
    )
    assert result["baseline"]["hvg_batch_aware"] is False


def test_missing_labels_are_dropped_and_counted(tmp_path):
    adata = sea.generate_demo_adata(seed=4)
    labels = adata.obs["cell_type"].to_numpy().astype(object)
    labels[:7] = None
    labels[7:12] = "Unknown"
    adata.obs["cell_type"] = labels
    path = _write_h5ad(adata, tmp_path / "missing.h5ad")
    proc = run_cli(
        _user_args(
            path,
            tmp_path / "out",
            "--batch-key",
            "donor",
            "--exclude-labels",
            "Unknown",
            "--embeddings",
            "X_corrected",
        )
    )
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    assert result["cells_dropped"] == {"missing_label": 7, "excluded_label": 5}
    assert result["cells_scored"] == 660 - 12 and result["labels"]["n_unique"] == 4
    assert "Unknown" not in result["labels"]["counts"]


def test_hvg_path_runs_and_is_batch_aware(tmp_path):
    path = _user_h5ad(tmp_path)
    proc = run_cli(
        _user_args(
            path,
            tmp_path / "out",
            "--batch-key",
            "donor",
            "--n-top-hvg",
            "100",
            "--n-pcs",
            "20",
            "--embeddings",
            "X_corrected",
        )
    )
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    assert (
        result["baseline"]["n_hvg_used"] == 100
        and result["baseline"]["n_pcs_used"] == 20
    )
    assert result["baseline"]["hvg_batch_aware"] is True
    assert result["embeddings"][sea.BASELINE_KEY]["shape"] == [660, 20]
    assert result["embeddings"][sea.CONTROL_KEY]["shape"] == [660, 20]


def test_subsampling_cap_is_applied_and_recorded(tmp_path):
    path = _user_h5ad(tmp_path)
    proc = run_cli(
        _user_args(
            path,
            tmp_path / "out",
            "--batch-key",
            "donor",
            "--max-cells",
            "120",
            "--embeddings",
            "X_corrected",
        )
    )
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    assert (
        result["subsampled"] is True
        and result["cells_scored"] == 120
        and result["n_cells"] == 660
    )
    assert sum(s["n_scored"] for s in result["strata"]) == 120
    assert all(s["n_scored"] >= 1 for s in result["strata"])
    assert (
        result["embeddings"][sea.BASELINE_KEY]["shape"][0] == 120
    )  # baseline fitted on the scored cells
    assert (
        result["effective"]["lisi_neighbors"]
        == result["effective"]["resample_sizes"][0]
        <= 90
    )
    assert "120" in (tmp_path / "out/report.md").read_text()


def test_degenerate_user_embedding_is_reported_not_scored(tmp_path):
    adata = sea.generate_demo_adata(seed=5)
    adata.obsm["X_flat"] = np.full((adata.n_obs, 4), 0.1)
    path = _write_h5ad(adata, tmp_path / "flat.h5ad")
    proc = run_cli(
        _user_args(
            path, tmp_path / "out", "--batch-key", "donor", "--embeddings", "X_flat"
        )
    )
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "out")
    entry = result["embeddings"]["X_flat"]
    assert all(v is None for v in entry["metrics"].values())
    assert all("zero variance" in r for r in entry["reasons"].values())
    assert all(
        v["verdict"] == "not evaluated"
        for v in result["verdicts"]["X_flat"]["per_metric"].values()
    )
    assert all(
        "zero variance" in v["reason"]
        for v in result["verdicts"]["X_flat"]["per_metric"].values()
    )


def test_rotated_and_scaled_copies_of_the_baseline_tie(tmp_path):
    """Independent oracle for the verdict rule: distance-preserving transforms of the baseline."""
    adata = sea.generate_demo_adata(seed=6)
    path = _write_h5ad(adata, tmp_path / "cells.h5ad")
    proc = run_cli(
        _user_args(
            path,
            tmp_path / "first",
            "--batch-key",
            "donor",
            "--embeddings",
            "X_corrected",
            "--n-pcs",
            "10",
        )
    )
    assert proc.returncode == 0, proc.stderr
    # recover the baseline recipe by running the skill's own function, then rotate and scale it
    counts = np.asarray(adata.layers["counts"])
    base, _ = sea.compute_pca_baseline(
        counts,
        adata.var_names,
        batch=adata.obs["donor"].to_numpy().astype(str),
        n_top_hvg=2000,
        n_pcs=10,
        random_state=0,
    )
    rng = np.random.default_rng(0)
    q, _ = np.linalg.qr(rng.normal(size=(10, 10)))
    adata.obsm["X_rotated"] = base @ q
    adata.obsm["X_scaled"] = 3.0 * base
    path2 = _write_h5ad(adata, tmp_path / "copies.h5ad")
    proc = run_cli(
        _user_args(
            path2,
            tmp_path / "second",
            "--batch-key",
            "donor",
            "--embeddings",
            "X_rotated,X_scaled",
            "--n-pcs",
            "10",
        )
    )
    assert proc.returncode == 0, proc.stderr
    result = _result(tmp_path / "second")
    for key in ("X_rotated", "X_scaled"):
        verdicts = {
            m: v["verdict"] for m, v in result["verdicts"][key]["per_metric"].items()
        }
        assert all(v in {"tied", "not separable"} for v in verdicts.values()), verdicts
        # every metric but KMeans (initialisation is scale/rotation sensitive) must tie exactly
        for m, v in verdicts.items():
            if m not in ("nmi_kmeans", "ari_kmeans"):
                assert v == "tied", (key, m, result["verdicts"][key]["per_metric"][m])


# --- rejections ----------------------------------------------------------------


def test_missing_embedding_key_exits_2_and_lists_keys(tmp_path):
    path = _user_h5ad(tmp_path)
    proc = run_cli(_user_args(path, tmp_path / "out", "--embeddings", "X_nope"))
    assert proc.returncode == 2
    assert "X_nope" in proc.stderr and "X_corrected" in proc.stderr
    assert not (tmp_path / "out/result.json").exists()


def test_missing_labels_key_exits_2(tmp_path):
    path = _user_h5ad(tmp_path)
    proc = run_cli(
        [
            "--input",
            str(path),
            "--output",
            str(tmp_path / "out"),
            "--labels-key",
            "nope",
        ]
    )
    assert (
        proc.returncode == 2
        and "nope" in proc.stderr
        and "available columns" in proc.stderr
    )


def test_labels_key_required_without_demo(tmp_path):
    path = _user_h5ad(tmp_path)
    proc = run_cli(["--input", str(path), "--output", str(tmp_path / "out")])
    assert proc.returncode == 2 and "--labels-key is required" in proc.stderr


def test_non_h5ad_input_exits_2(tmp_path):
    bogus = tmp_path / "cells.csv"
    bogus.write_text("a,b\n1,2\n")
    proc = run_cli(
        ["--input", str(bogus), "--output", str(tmp_path / "out"), "--labels-key", "x"]
    )
    assert proc.returncode == 2 and ".h5ad" in proc.stderr


def test_normalised_matrix_without_counts_layer_exits_2(tmp_path):
    adata = sea.generate_demo_adata(seed=3)
    del adata.layers["counts"]
    scanpy.pp.normalize_total(adata, target_sum=1e4)
    scanpy.pp.log1p(adata)
    path = _write_h5ad(adata, tmp_path / "norm.h5ad")
    proc = run_cli(_user_args(path, tmp_path / "out"))
    assert proc.returncode == 2
    assert "normalized/log-transformed" in proc.stderr


def test_negative_values_in_counts_layer_exit_2(tmp_path):
    adata = sea.generate_demo_adata(seed=3)
    adata.layers["counts"] = np.asarray(adata.layers["counts"]) - 1.0
    path = _write_h5ad(adata, tmp_path / "neg.h5ad")
    proc = run_cli(_user_args(path, tmp_path / "out"))
    assert proc.returncode == 2 and "negative" in proc.stderr


def test_too_few_cells_exits_2(tmp_path):
    adata = sea.generate_demo_adata(seed=4)[:40].copy()
    path = _write_h5ad(adata, tmp_path / "tiny.h5ad")
    proc = run_cli(_user_args(path, tmp_path / "out"))
    assert proc.returncode == 2 and f"at least {sea.MIN_CELLS}" in proc.stderr


@pytest.mark.parametrize(
    "flag, value, needle",
    [
        ("--lisi-neighbors", "2", "--lisi-neighbors must be at least 3"),
        ("--n-resamples", "5", "--n-resamples must be at least 10"),
        ("--max-cells", "20", "--max-cells must be at least 50"),
        ("--graph-neighbors", "1", "--graph-neighbors must be at least 2"),
        ("--transfer-k", "0", "--transfer-k must be at least 1"),
    ],
)
def test_parameter_bounds_exit_2(tmp_path, flag, value, needle):
    proc = run_cli(["--demo", "--output", str(tmp_path / "out"), flag, value])
    assert proc.returncode == 2 and needle in proc.stderr


def test_non_empty_output_requires_overwrite(tmp_path):
    out = tmp_path / "out"
    out.mkdir()
    (out / "stale.txt").write_text("x")
    proc = run_cli(["--demo", "--output", str(out), *FAST])
    assert proc.returncode == 2 and "--overwrite" in proc.stderr
    proc = run_cli(["--demo", "--output", str(out), "--overwrite", *FAST])
    assert proc.returncode == 0, proc.stderr


def test_demo_rejects_input_only_flags(tmp_path):
    proc = run_cli(
        ["--demo", "--output", str(tmp_path / "out"), "--counts-layer", "counts"]
    )
    assert proc.returncode == 2 and "require --input" in proc.stderr


def test_expected_sha_mismatch_refuses_before_writing(tmp_path):
    path = _user_h5ad(tmp_path)
    proc = run_cli(
        _user_args(path, tmp_path / "out", "--expected-input-sha256", "0" * 64)
    )
    assert proc.returncode == 2 and "SHA-256 differs" in proc.stderr
    assert not (tmp_path / "out").exists()


# --- registration ----------------------------------------------------------------


def test_cli_registered_with_matching_flags():
    cli_text = (PROJECT_ROOT / "clawbio/cli.py").read_text()
    assert '"scrna-embedding-audit": {' in cli_text
    parser = sea._build_parser()
    option_strings = {
        opt for action in parser._actions for opt in action.option_strings
    }
    for flag in (
        "--embeddings",
        "--labels-key",
        "--batch-key",
        "--exclude-labels",
        "--counts-layer",
        "--n-top-hvg",
        "--n-pcs",
        "--lisi-neighbors",
        "--graph-neighbors",
        "--transfer-k",
        "--max-cells",
        "--n-resamples",
        "--random-state",
        "--expected-input-sha256",
        "--overwrite",
    ):
        assert flag in option_strings, flag


def test_skill_md_under_500_lines_and_has_required_sections():
    text = (SKILL_DIR / "SKILL.md").read_text()
    assert len(text.splitlines()) < 500
    for heading in (
        "## Trigger",
        "## Scope",
        "## Workflow",
        "## Example Output",
        "## Gotchas",
        "## Safety",
        "## Agent Boundary",
        "## Output Structure",
        "## Integration with Bio Orchestrator",
        "## Maintenance",
        "## Citations",
    ):
        assert heading in text, heading
    assert text.count("**You will want to") >= 3
    assert DISCLAIMER_START.split(".")[0] in text
