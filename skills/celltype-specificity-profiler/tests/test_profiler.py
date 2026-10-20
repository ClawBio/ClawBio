"""Tests for the celltype-specificity-profiler skill."""
import json
import os
import subprocess
import sys

import numpy as np
import pytest

SKILL_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, SKILL_DIR)

from profiler import bimodality_coefficient, compute_tau  # noqa: E402


# --------------------------- compute_tau ---------------------------------- #
def test_tau_uniform_is_zero():
    assert compute_tau([3.0, 3.0, 3.0, 3.0]) == pytest.approx(0.0)


def test_tau_single_celltype_is_one():
    # one cell type expresses, the rest are off -> tau == 1
    assert compute_tau([5.0, 0.0, 0.0, 0.0]) == pytest.approx(1.0)


def test_tau_in_unit_interval_for_random_inputs():
    rng = np.random.default_rng(0)
    for _ in range(50):
        means = rng.uniform(0, 10, size=rng.integers(2, 12))
        tau = compute_tau(means)
        assert 0.0 <= tau <= 1.0


def test_tau_all_zero_is_zero():
    assert compute_tau([0.0, 0.0, 0.0]) == 0.0


def test_tau_intermediate_value():
    # x = [2, 1]; tau = (1 - 2/2) + (1 - 1/2) / (n-1=1) = 0.5
    assert compute_tau([2.0, 1.0]) == pytest.approx(0.5)


def test_tau_needs_two_celltypes():
    assert np.isnan(compute_tau([4.0]))


# --------------------- bimodality_coefficient ----------------------------- #
def test_bc_bimodal_high():
    # two well-separated clusters -> strongly bimodal -> BC near/above ~0.8
    vals = [0.0] * 50 + [10.0] * 50
    assert bimodality_coefficient(vals) > 0.8


def test_bc_unimodal_low():
    rng = np.random.default_rng(1)
    vals = rng.normal(0, 1, size=2000)
    # a normal distribution has BC well below the 5/9 bimodality benchmark
    assert bimodality_coefficient(vals) < 0.555


def test_bc_needs_four_values():
    assert np.isnan(bimodality_coefficient([1.0, 2.0, 3.0]))


def test_bc_constant_is_nan():
    assert np.isnan(bimodality_coefficient([2.0, 2.0, 2.0, 2.0, 2.0]))


# --------------------------- end-to-end ----------------------------------- #
# --demo runs on the real, scanpy-bundled pbmc3k dataset (no synthetic fixture).
# Default demo gene is MS4A1, a canonical B-cell marker -> highly cell-type-specific.
def test_end_to_end_demo(tmp_path):
    out = tmp_path / "out"
    result = subprocess.run(
        [
            sys.executable,
            os.path.join(SKILL_DIR, "profiler.py"),
            "--demo",
            "--trial-prior",
            "--out",
            str(out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    # stdout is valid JSON
    stdout_profile = json.loads(result.stdout)
    assert stdout_profile["skill"] == "celltype-specificity-profiler"

    profile_path = out / "profile.json"
    assert profile_path.exists()
    profile = json.loads(profile_path.read_text())

    expected_keys = {
        "skill",
        "gene",
        "atlas",
        "tau",
        "tau_threshold",
        "n_cell_types_used_for_tau",
        "n_cell_types_excluded_small",
        "bimodality_coefficient",
        "interpretation",
        "low_expression",
        "top_cell_types",
        "per_celltype_stats",
        "trial_prior",
    }
    assert expected_keys <= set(profile.keys())

    assert profile["gene"] == "MS4A1"
    assert 0.0 <= profile["tau"] <= 1.0
    assert profile["low_expression"] is False
    assert profile["top_cell_types"][0]["cell_type"] == "B cells"
    assert profile["trial_prior"]["phase_I_to_II_OR"] == 1.27
    # pbmc3k has 8 louvain clusters; Megakaryocytes (~15 cells) is dropped from
    # tau under the Zhang et al. 2026 <20-cell exclusion.
    assert profile["n_cell_types_excluded_small"] == 1
    assert profile["n_cell_types_used_for_tau"] == 7

    # side-car files
    assert (out / "per_celltype.csv").exists()
    repro = out / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "environment.yml").exists()
    assert (repro / "checksums.sha256").exists()


def test_demo_tau_in_expected_range(tmp_path):
    out = tmp_path / "out2"
    subprocess.run(
        [
            sys.executable,
            os.path.join(SKILL_DIR, "profiler.py"),
            "--demo",
            "--out",
            str(out),
        ],
        check=True,
        capture_output=True,
        text=True,
    )
    profile = json.loads((out / "profile.json").read_text())
    # MS4A1 is B-cell restricted in pbmc3k -> very high tau.
    assert profile["tau"] >= 0.85
    assert profile["interpretation"] == "cell-type-specific (tau > 0.69)"


def test_demo_reproducibility_bundle(tmp_path):
    import hashlib

    out = tmp_path / "out3"
    subprocess.run(
        [sys.executable, os.path.join(SKILL_DIR, "profiler.py"), "--demo", "--output", str(out)],
        check=True,
        capture_output=True,
        text=True,
    )
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).exists(), name

    commands = (repro / "commands.sh").read_text()
    assert "python skills/celltype-specificity-profiler/profiler.py" in commands
    assert "--demo" in commands
    assert SKILL_DIR not in commands

    env = (repro / "environment.yml").read_text()
    for dep in ("scanpy", "anndata", "numpy", "scipy", "pandas"):
        assert dep in env

    labels = set()
    for line in (repro / "checksums.sha256").read_text().splitlines():
        digest, label = line.split("  ", 1)
        target = out / label
        assert target.resolve().is_relative_to(out.resolve()), label
        assert hashlib.sha256(target.read_bytes()).hexdigest() == digest, label
        labels.add(label)
    assert {"profile.json", "per_celltype.csv", "reproducibility/environment.yml"} <= labels
