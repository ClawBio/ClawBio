"""Property tests for the metric implementations (Hypothesis).

Hypothesis is not yet in the ClawBio lockfile (ClawBio/ClawBio#537 is open), so
this module is skipped where it is absent and runs locally with
``uv run --with hypothesis pytest skills/scrna-embedding-audit/tests``.
"""

from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest

hypothesis = pytest.importorskip("hypothesis")
from hypothesis import HealthCheck, given, settings
from hypothesis import strategies as st
from hypothesis.extra import numpy as hnp

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import embedding_audit_metrics as eam

SETTINGS = settings(
    max_examples=25, deadline=None, suppress_health_check=[HealthCheck.too_slow]
)


@st.composite
def labelled_embeddings(draw, min_cells: int = 24, max_cells: int = 60):
    n = draw(st.integers(min_cells, max_cells))
    dims = draw(st.integers(2, 6))
    # unique=True keeps every coordinate distinct, so no two cells coincide and
    # k-NN ties cannot make the metrics depend on cell order (duplicate cells are
    # a documented limitation, not a property under test here)
    x = draw(
        hnp.arrays(
            np.float64,
            (n, dims),
            elements=st.floats(
                -50, 50, allow_nan=False, allow_infinity=False, width=32
            ),
            unique=True,
        )
    )
    # guarantee at least two labels and two batches, each with >= 2 cells
    n_labels = draw(st.integers(2, 4))
    labels = np.array([f"L{i % n_labels}" for i in range(n)])
    batches = np.array([f"B{(i // 2) % 2}" for i in range(n)])
    return x, labels, batches


def _all_metrics(x, labels, batches):
    out = eam.score_embedding(
        x,
        labels,
        batches,
        lisi_neighbors=15,
        graph_neighbors=6,
        transfer_k=5,
        random_state=0,
    )
    return out["metrics"], out["reasons"]


@SETTINGS
@given(labelled_embeddings(), st.randoms(use_true_random=False))
def test_metrics_are_invariant_to_cell_order(data, rnd):
    x, labels, batches = data
    perm = np.arange(x.shape[0])
    rnd.shuffle(perm)
    m1, r1 = _all_metrics(x, labels, batches)
    m2, r2 = _all_metrics(x[perm], labels[perm], batches[perm])
    assert set(r1) == set(r2)
    for name in eam.METRIC_NAMES:
        if name in ("nmi_kmeans", "ari_kmeans"):
            continue  # KMeans initialisation depends on row order; covered separately
        if m1[name] is None:
            assert m2[name] is None
        else:
            assert m2[name] == pytest.approx(m1[name], abs=1e-9), name


@SETTINGS
@given(labelled_embeddings())
def test_metrics_are_invariant_to_category_names(data):
    x, labels, batches = data
    m1, _ = _all_metrics(x, labels, batches)
    renamed_labels = np.char.add("cell type ", labels)
    renamed_batches = np.char.add("donor-", batches)
    m2, _ = _all_metrics(x, renamed_labels, renamed_batches)
    for name in eam.METRIC_NAMES:
        if m1[name] is None:
            assert m2[name] is None
        else:
            assert m2[name] == pytest.approx(m1[name], abs=1e-12), name


@SETTINGS
@given(labelled_embeddings())
def test_every_reported_metric_is_within_its_documented_range(data):
    x, labels, batches = data
    metrics, reasons = _all_metrics(x, labels, batches)
    for name, value in metrics.items():
        if value is None:
            assert name in reasons
        elif name == "ari_kmeans":
            assert -1.0 <= value <= 1.0 + 1e-9, (name, value)  # ARI can be negative
        else:
            assert -1e-9 <= value <= 1.0 + 1e-9, (name, value)


@SETTINGS
@given(labelled_embeddings())
def test_lisi_per_cell_lies_between_one_and_number_of_categories(data):
    x, labels, batches = data
    dists, idx = eam.exact_knn(x, min(15, x.shape[0]))
    for values in (labels, batches):
        lisi = eam.lisi_per_cell(dists, idx, values)
        n_cats = len(np.unique(values))
        # -1 marks a cell whose neighbourhood kernel degenerated (scib-metrics convention)
        regular = lisi[lisi > 0]
        assert np.all(regular >= 1.0 - 1e-6)
        assert np.all(regular <= n_cats + 1e-6)
        assert np.all((lisi > 0) | (lisi == -1.0))


@SETTINGS
@given(labelled_embeddings(), st.integers(0, 10_000))
def test_nan_anywhere_makes_the_skill_abstain(data, position):
    x, labels, batches = data
    x = x.copy()
    flat = x.reshape(-1)
    flat[position % flat.size] = np.nan
    metrics, reasons = _all_metrics(x, labels, batches)
    assert all(v is None for v in metrics.values())
    assert all("non-finite" in r for r in reasons.values())


@SETTINGS
@given(
    st.integers(24, 60),
    st.integers(1, 5),
    st.floats(-1e3, 1e3, allow_nan=False, allow_infinity=False),
)
def test_constant_embedding_makes_the_skill_abstain(n, dims, value):
    x = np.full((n, dims), value)
    labels = np.array([f"L{i % 2}" for i in range(n)])
    batches = np.array([f"B{(i // 2) % 2}" for i in range(n)])
    metrics, reasons = _all_metrics(x, labels, batches)
    assert all(v is None for v in metrics.values())
    assert all("zero variance" in r for r in reasons.values())


@SETTINGS
@given(
    st.lists(
        st.floats(-1, 1, allow_nan=False, allow_infinity=False), min_size=1, max_size=30
    )
)
def test_verdict_rule_is_exhaustive_and_consistent(diffs):
    out = eam.verdict(diffs)
    arr = np.asarray(diffs)
    arr = np.where(np.abs(arr) <= eam.DIFF_TOLERANCE, 0.0, arr)
    sd = np.std(arr, ddof=1) if arr.size > 1 else 0.0
    if np.all(arr == 0):
        assert out["verdict"] == "tied"
    elif np.all(arr > 0) and np.median(arr) > eam.SD_FACTOR * sd:
        assert out["verdict"] == "higher"
    elif np.all(arr < 0) and -np.median(arr) > eam.SD_FACTOR * sd:
        assert out["verdict"] == "lower"
    else:
        assert out["verdict"] == "not separable"
    assert out["frac_positive"] == pytest.approx(np.mean(arr > 0))
    assert out["median_diff"] == pytest.approx(np.median(arr))


@SETTINGS
@given(
    st.integers(30, 200),
    st.integers(2, 4),
    st.integers(1, 3),
    st.floats(0.3, 0.9),
    st.integers(0, 1000),
)
def test_stratified_subsample_keeps_every_stratum_and_hits_the_target(
    n, n_labels, n_batches, fraction, seed
):
    labels = np.array([f"L{i % n_labels}" for i in range(n)])
    batches = np.array([f"B{(i // n_labels) % n_batches}" for i in range(n)])
    idx = eam.stratified_subsample(
        labels, batches, fraction=fraction, random_state=seed
    )
    assert len(np.unique(idx)) == len(idx)
    strata_all = set(zip(labels, batches))
    strata_kept = set(zip(labels[idx], batches[idx]))
    assert strata_kept == strata_all
    assert len(idx) == max(len(strata_all), round(fraction * n))
