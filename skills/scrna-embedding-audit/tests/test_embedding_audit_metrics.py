"""Metric definitions for scrna-embedding-audit, checked against scib-metrics 0.6.1.

The fixture under ``fixtures/`` was produced by ``fixtures/generate_scib_fixture.py``
in an isolated environment (scib-metrics 0.6.1, Python 3.12). These tests hold the
numpy/scikit-learn re-implementation to that fixture and to closed-form cases whose
answers do not depend on the implementation.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pytest
from sklearn.metrics import adjusted_rand_score

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import embedding_audit_metrics as eam

FIXTURE_DIR = SKILL_DIR / "fixtures"
EXPECTED = json.loads((FIXTURE_DIR / "scib_metrics_0.6.1_expected.json").read_text())
INPUTS = np.load(FIXTURE_DIR / "scib_fixture_input.npz", allow_pickle=False)
TOL = 1e-5  # observed worst case 2.9e-7 (jax float32 upstream)


def _fixture(key: str):
    return INPUTS[key], INPUTS["labels"], INPUTS["batches"], EXPECTED["embeddings"][key]


@pytest.mark.parametrize("key", ["X_mixed", "X_batchy"])
class TestScibFixtureAgreement:
    def test_silhouette_label(self, key):
        x, labels, _, exp = _fixture(key)
        assert eam.silhouette_label(x, labels) == pytest.approx(
            exp["silhouette_label"], abs=TOL
        )

    def test_silhouette_batch(self, key):
        x, labels, batches, exp = _fixture(key)
        assert eam.silhouette_batch(x, labels, batches) == pytest.approx(
            exp["silhouette_batch"], abs=TOL
        )

    def test_bras(self, key):
        x, labels, batches, exp = _fixture(key)
        assert eam.bras(x, labels, batches) == pytest.approx(exp["bras"], abs=TOL)

    def test_isolated_labels(self, key):
        x, labels, batches, exp = _fixture(key)
        assert eam.isolated_labels(x, labels, batches) == pytest.approx(
            exp["isolated_labels"], abs=TOL
        )

    def test_kmeans_nmi_ari_and_partition(self, key):
        x, labels, _, exp = _fixture(key)
        nmi, ari, pred = eam.kmeans_nmi_ari(x, labels, random_state=0)
        assert nmi == pytest.approx(exp["nmi"], abs=TOL)
        assert ari == pytest.approx(exp["ari"], abs=TOL)
        # same partition as the jax KMeans in scib-metrics, up to label permutation
        assert adjusted_rand_score(pred, exp["kmeans_labels"]) == pytest.approx(1.0)

    def test_lisi_scores(self, key):
        x, labels, batches, exp = _fixture(key)
        k = EXPECTED["n_neighbors"]["lisi"]
        dists, idx = eam.exact_knn(x, k)
        assert eam.clisi(dists, idx, labels) == pytest.approx(exp["clisi"], abs=TOL)
        assert eam.ilisi(dists, idx, batches) == pytest.approx(exp["ilisi"], abs=TOL)

    def test_graph_connectivity(self, key):
        x, labels, _, exp = _fixture(key)
        k = EXPECTED["n_neighbors"]["graph_connectivity"]
        dists, idx = eam.exact_knn(x, k)
        assert eam.graph_connectivity(dists, idx, labels) == pytest.approx(
            exp["graph_connectivity"], abs=TOL
        )


def test_fixture_records_its_provenance():
    assert EXPECTED["packages"]["scib-metrics"] == "0.6.1"
    assert EXPECTED["perplexity"] == EXPECTED["n_neighbors"]["lisi"] // 3
    assert EXPECTED["generator"] == "generate_scib_fixture.py"
    assert (FIXTURE_DIR / "generate_scib_fixture.py").is_file()


def test_exact_knn_includes_self_first():
    rng = np.random.default_rng(0)
    x = rng.normal(size=(40, 3))
    dists, idx = eam.exact_knn(x, 5)
    assert idx.shape == (40, 5) and dists.shape == (40, 5)
    assert np.array_equal(idx[:, 0], np.arange(40))
    assert np.all(dists[:, 0] == 0.0)
    assert np.all(np.diff(dists, axis=1) >= 0)
    # no cell lists itself twice
    assert all(len(set(row)) == 5 for row in idx)


def test_exact_knn_keeps_self_edge_for_duplicate_cells():
    x = np.zeros((6, 2))  # six identical cells
    dists, idx = eam.exact_knn(x, 3)
    assert np.array_equal(idx[:, 0], np.arange(6))
    assert np.all(dists == 0.0)
    assert eam.duplicate_row_fraction(x) == 1.0
    assert eam.duplicate_row_fraction(np.arange(12.0).reshape(6, 2)) == 0.0


def _two_blobs(n_per: int = 40, gap: float = 100.0, seed: int = 1):
    rng = np.random.default_rng(seed)
    x = np.vstack([rng.normal(size=(n_per, 4)), rng.normal(size=(n_per, 4)) + gap])
    labels = np.array(["a"] * n_per + ["b"] * n_per)
    return x, labels


class TestClosedFormCases:
    def test_far_apart_labels_score_perfectly(self):
        x, labels = _two_blobs()
        assert eam.silhouette_label(x, labels) > 0.95
        dists, idx = eam.exact_knn(x, 15)
        assert eam.clisi(dists, idx, labels) == pytest.approx(1.0)
        nmi, ari, _ = eam.kmeans_nmi_ari(x, labels, random_state=0)
        assert nmi == pytest.approx(1.0) and ari == pytest.approx(1.0)

    def test_graph_connectivity_is_one_on_a_chain(self):
        # equally spaced points on a line with k=3 (self + both neighbours):
        # every edge is reciprocated, so each label is one strong component
        x = np.arange(40, dtype=float)[:, None]
        labels = np.where(np.arange(40) < 20, "a", "b")
        dists, idx = eam.exact_knn(x, 3)
        assert eam.graph_connectivity(dists, idx, labels) == pytest.approx(1.0)

    def test_graph_connectivity_counts_unreferenced_cells_as_disconnected(self):
        # scib-metrics uses strongly connected components on the directed k-NN
        # graph, so a cell that nobody lists as a neighbour is its own component
        x, labels = _two_blobs()
        dists, idx = eam.exact_knn(x, 15)
        value = eam.graph_connectivity(dists, idx, labels)
        assert 0.9 <= value <= 1.0

    def test_batch_identical_to_label_gives_zero_ilisi(self):
        x, labels = _two_blobs()
        dists, idx = eam.exact_knn(x, 15)
        assert eam.ilisi(dists, idx, labels) == pytest.approx(0.0, abs=1e-9)

    def test_perfectly_mixed_batches_give_high_ilisi(self):
        rng = np.random.default_rng(3)
        x = rng.normal(size=(300, 5))
        batches = np.where(np.arange(300) % 2 == 0, "b1", "b2")
        dists, idx = eam.exact_knn(x, 90)
        assert eam.ilisi(dists, idx, batches) > 0.8

    def test_lisi_per_cell_is_bounded_by_number_of_categories(self):
        rng = np.random.default_rng(4)
        x = rng.normal(size=(120, 3))
        cats = rng.choice(["p", "q", "r"], size=120)
        dists, idx = eam.exact_knn(x, 30)
        lisi = eam.lisi_per_cell(dists, idx, cats)
        assert lisi.shape == (120,)
        assert np.all(lisi >= 1.0 - 1e-9) and np.all(lisi <= 3.0 + 1e-9)

    def test_silhouette_batch_undefined_when_every_label_has_one_batch(self):
        x, labels = _two_blobs()
        assert eam.silhouette_batch(x, labels, labels) is None
        assert eam.bras(x, labels, labels) is None

    def test_isolated_labels_uses_only_labels_in_fewest_batches(self):
        rng = np.random.default_rng(8)
        # label "a" in two batches, label "b" in one batch -> only "b" is isolated
        x = np.vstack([rng.normal(size=(60, 3)), rng.normal(size=(30, 3)) + 50.0])
        labels = np.array(["a"] * 60 + ["b"] * 30)
        batches = np.array(["d1"] * 30 + ["d2"] * 30 + ["d1"] * 30)
        from sklearn.metrics import silhouette_samples

        per_cell = (silhouette_samples(x, labels) + 1) / 2
        assert eam.isolated_labels(x, labels, batches) == pytest.approx(
            per_cell[labels == "b"].mean()
        )

    def test_bras_is_zero_variance_safe_and_cosine_based(self):
        rng = np.random.default_rng(9)
        x = rng.normal(size=(80, 4))
        labels = np.array(["a"] * 40 + ["b"] * 40)
        batches = np.tile(["d1", "d2"], 40)
        value = eam.bras(x, labels, batches)
        assert 0.0 <= value <= 1.0
        # scaling every cell does not change cosine distances
        assert eam.bras(3.0 * x, labels, batches) == pytest.approx(value)
        with pytest.raises(ValueError):
            eam.bras(
                np.vstack([x, np.zeros((1, 4))]),
                np.append(labels, "a"),
                np.append(batches, "d1"),
            )

    def test_graph_connectivity_detects_split_label(self):
        # label "a" is two far-apart clumps: largest component is half the label
        rng = np.random.default_rng(5)
        x = np.vstack([rng.normal(size=(20, 2)), rng.normal(size=(20, 2)) + 100.0])
        labels = np.array(["a"] * 40)
        dists, idx = eam.exact_knn(x, 5)
        assert eam.graph_connectivity(dists, idx, labels) == pytest.approx(0.5)


class TestKnnTransfer:
    def test_separable_data_transfers_perfectly(self):
        x, labels = _two_blobs()
        batches = np.tile(["d1", "d2"], len(labels) // 2)
        out = eam.knn_transfer_accuracy(x, labels, batches, k=5)
        assert out["accuracy"] == pytest.approx(1.0)
        assert set(out["per_batch_accuracy"]) == {"d1", "d2"}
        assert out["n_cells_with_unseen_label"] == 0

    def test_label_absent_from_training_batches_is_counted_separately(self):
        x, labels = _two_blobs()
        # "b" only in d2 and "a" only in d1: every held-out cell has an unseen label
        batches = np.where(labels == "b", "d2", "d1")
        with pytest.raises(ValueError):
            eam.knn_transfer_accuracy(x, labels, batches, k=5)
        # add a few "a" cells to d2 so the metric is defined; the "b" cells stay unseen
        batches = batches.copy()
        batches[:5] = "d2"
        out = eam.knn_transfer_accuracy(x, labels, batches, k=5)
        assert out["n_cells_with_unseen_label"] == int(np.sum(labels == "b"))
        assert out["accuracy"] == pytest.approx(1.0)  # the seen "a" cells are recovered

    def test_matches_sklearn_leave_one_group_out_oracle(self):
        from sklearn.metrics import balanced_accuracy_score
        from sklearn.model_selection import LeaveOneGroupOut, cross_val_predict
        from sklearn.neighbors import KNeighborsClassifier

        rng = np.random.default_rng(11)
        x = np.vstack(
            [
                rng.normal(size=(60, 3)),
                rng.normal(size=(60, 3)) + 1.5,
                rng.normal(size=(30, 3)) - 1.5,
            ]
        )
        labels = np.array(["a"] * 60 + ["b"] * 60 + ["c"] * 30)
        batches = np.array([f"d{i % 3}" for i in range(150)])
        out = eam.knn_transfer_accuracy(x, labels, batches, k=7)
        pred = cross_val_predict(
            KNeighborsClassifier(n_neighbors=7, weights="distance"),
            x,
            labels,
            groups=batches,
            cv=LeaveOneGroupOut(),
        )
        assert out["n_cells_with_unseen_label"] == 0
        assert out["accuracy"] == pytest.approx(balanced_accuracy_score(labels, pred))

    def test_requires_two_batches(self):
        x, labels = _two_blobs()
        with pytest.raises(ValueError):
            eam.knn_transfer_accuracy(x, labels, np.array(["only"] * len(labels)), k=5)


class TestDegenerateInput:
    def test_nan_is_reported(self):
        x = np.ones((10, 2))
        x[3, 1] = np.nan
        assert eam.degeneracy_reason(x) == "non-finite values in 1 cells"

    def test_constant_embedding_is_reported_on_raw_values(self):
        assert (
            eam.degeneracy_reason(np.full((10, 3), 0.1))
            == "zero variance in every dimension"
        )

    def test_healthy_embedding_has_no_reason(self):
        assert (
            eam.degeneracy_reason(np.random.default_rng(0).normal(size=(10, 3))) is None
        )

    def test_empty_input_raises(self):
        with pytest.raises(ValueError):
            eam.silhouette_label(np.empty((0, 3)), np.array([]))

    def test_single_label_raises_in_primitive(self):
        x = np.random.default_rng(0).normal(size=(20, 3))
        with pytest.raises(ValueError):
            eam.silhouette_label(x, np.array(["same"] * 20))


class TestScoreEmbedding:
    def test_all_metrics_present_with_labels_and_batches(self):
        x, labels, batches, _ = _fixture("X_mixed")
        out = eam.score_embedding(
            x,
            labels,
            batches,
            lisi_neighbors=90,
            graph_neighbors=15,
            transfer_k=15,
            random_state=0,
        )
        assert set(out["metrics"]) == set(eam.METRIC_NAMES)
        assert all(v is not None for v in out["metrics"].values())
        assert out["reasons"] == {}
        assert all(0.0 <= v <= 1.0 for v in out["metrics"].values())
        assert out["details"]["duplicate_row_fraction"] == 0.0
        assert out["details"]["lisi_degenerate_cells"] == 0
        assert set(eam.METRIC_GROUPS.values()) == {"bio", "batch", "transfer"}

    def test_without_batches_batch_metrics_are_null_with_reason(self):
        x, labels, _, _ = _fixture("X_mixed")
        out = eam.score_embedding(
            x,
            labels,
            None,
            lisi_neighbors=90,
            graph_neighbors=15,
            transfer_k=15,
            random_state=0,
        )
        for name in (
            "silhouette_batch",
            "bras",
            "ilisi",
            "isolated_labels",
            "knn_transfer_accuracy",
        ):
            assert out["metrics"][name] is None
            assert "batch" in out["reasons"][name]
        # graph connectivity is grouped with batch mixing in scIB but needs only labels
        for name in (
            "silhouette_label",
            "nmi_kmeans",
            "ari_kmeans",
            "clisi",
            "graph_connectivity",
        ):
            assert out["metrics"][name] is not None
            assert name not in out["reasons"]

    def test_single_label_abstains_on_label_metrics(self):
        x, _, batches, _ = _fixture("X_mixed")
        labels = np.array(["same"] * x.shape[0])
        out = eam.score_embedding(
            x,
            labels,
            batches,
            lisi_neighbors=90,
            graph_neighbors=15,
            transfer_k=15,
            random_state=0,
        )
        assert out["metrics"]["silhouette_label"] is None
        assert "fewer than 2 labels" in out["reasons"]["silhouette_label"]
        assert out["metrics"]["ilisi"] is not None

    def test_degenerate_embedding_abstains_everywhere(self):
        _, labels, batches, _ = _fixture("X_mixed")
        x = np.full((labels.shape[0], 4), 2.5)
        out = eam.score_embedding(
            x,
            labels,
            batches,
            lisi_neighbors=90,
            graph_neighbors=15,
            transfer_k=15,
            random_state=0,
        )
        assert all(v is None for v in out["metrics"].values())
        assert all("zero variance" in r for r in out["reasons"].values())

    def test_k_is_reduced_for_small_inputs_and_recorded(self):
        rng = np.random.default_rng(0)
        x = rng.normal(size=(40, 3))
        labels = np.where(np.arange(40) < 20, "a", "b")
        batches = np.where(np.arange(40) % 2 == 0, "d1", "d2")
        out = eam.score_embedding(
            x,
            labels,
            batches,
            lisi_neighbors=90,
            graph_neighbors=15,
            transfer_k=15,
            random_state=0,
        )
        assert out["effective"]["lisi_neighbors"] == 40  # k counts the cell itself
        assert out["effective"]["graph_neighbors"] == 15


class TestVerdictRule:
    def test_consistent_sign_with_clear_margin_is_higher(self):
        v = eam.verdict([0.10, 0.12, 0.11, 0.09, 0.10])
        assert v["verdict"] == "higher" and v["frac_positive"] == 1.0

    def test_consistent_sign_with_clear_margin_is_lower(self):
        assert eam.verdict([-0.10, -0.12, -0.11, -0.09])["verdict"] == "lower"

    def test_consistent_sign_without_margin_is_not_separable(self):
        # all positive but spread (SD) is large relative to the median
        v = eam.verdict([0.001, 0.2, 0.002, 0.3, 0.001])
        assert v["verdict"] == "not separable"
        assert "SD" in v["reason"]

    def test_mixed_signs_not_separable(self):
        v = eam.verdict([0.1, -0.1, 0.2])
        assert v["verdict"] == "not separable" and "sign" in v["reason"]

    def test_a_zero_difference_breaks_consistency(self):
        assert eam.verdict([0.1, 0.0, 0.2])["verdict"] == "not separable"

    def test_all_zero_differences_are_tied(self):
        assert eam.verdict([0.0, 0.0, 0.0])["verdict"] == "tied"

    def test_differences_within_tolerance_count_as_zero(self):
        assert eam.verdict([1e-6, -3e-7, 5e-6])["verdict"] == "tied"
        assert eam.verdict([0.1, 1e-6, 0.1])["verdict"] == "not separable"

    def test_both_at_ceiling_is_tied_even_with_consistent_differences(self):
        v = eam.verdict(
            [1e-4, 1.1e-4, 0.9e-4], full_embedding=1.0, full_baseline=0.9995
        )
        assert v["verdict"] == "tied" and v["reason"] == "both at ceiling"
        v = eam.verdict([1e-4, 1.1e-4, 0.9e-4], full_embedding=1.0, full_baseline=0.95)
        assert v["verdict"] == "higher"

    def test_any_missing_value_is_not_evaluated(self):
        v = eam.verdict([0.1, None, 0.2])
        assert v["verdict"] == "not evaluated" and v["reason"]

    def test_median_and_sd_reported(self):
        v = eam.verdict([0.1, 0.3, 0.2])
        assert v["median_diff"] == pytest.approx(0.2)
        assert v["sd_diff"] == pytest.approx(np.std([0.1, 0.3, 0.2], ddof=1))


class TestStratifiedSubsample:
    def test_keeps_every_stratum_and_is_seeded(self):
        labels = np.repeat(["a", "b", "c"], 50)
        batches = np.tile(["d1", "d2"], 75)
        idx1 = eam.stratified_subsample(labels, batches, fraction=0.5, random_state=7)
        idx2 = eam.stratified_subsample(labels, batches, fraction=0.5, random_state=7)
        assert np.array_equal(idx1, idx2)
        assert len(idx1) == 75
        strata = {(l, b) for l, b in zip(labels[idx1], batches[idx1])}
        assert len(strata) == 6

    def test_cap_in_cells(self):
        labels = np.repeat(["a", "b"], 500)
        idx = eam.stratified_subsample(labels, None, n_cells=100, random_state=0)
        assert len(idx) == 100
        assert set(labels[idx]) == {"a", "b"}


# --- seeded copies of the Hypothesis properties (run in CI without hypothesis) ---------


def _tie_free_case(seed: int):
    rng = np.random.default_rng(seed)
    n = int(rng.integers(40, 80))
    x = rng.normal(size=(n, int(rng.integers(3, 6))))
    labels = np.array([f"L{i % 3}" for i in range(n)])
    batches = np.array([f"B{(i // 2) % 2}" for i in range(n)])
    return x, labels, batches


def _score(x, labels, batches):
    return eam.score_embedding(
        x,
        labels,
        batches,
        lisi_neighbors=15,
        graph_neighbors=6,
        transfer_k=5,
        random_state=0,
    )


@pytest.mark.parametrize("seed", range(5))
def test_seeded_cell_order_invariance_for_distance_based_metrics(seed):
    x, labels, batches = _tie_free_case(seed)
    perm = np.random.default_rng(100 + seed).permutation(x.shape[0])
    a = _score(x, labels, batches)["metrics"]
    b = _score(x[perm], labels[perm], batches[perm])["metrics"]
    for name in eam.METRIC_NAMES:
        if name in ("nmi_kmeans", "ari_kmeans"):
            continue  # KMeans initialisation reads row order; checked only on separable data
        assert b[name] == pytest.approx(a[name], abs=1e-9), name


@pytest.mark.parametrize("seed", range(5))
def test_seeded_category_name_invariance(seed):
    x, labels, batches = _tie_free_case(seed)
    a = _score(x, labels, batches)["metrics"]
    b = _score(x, np.char.add("type ", labels), np.char.add("donor-", batches))[
        "metrics"
    ]
    for name in eam.METRIC_NAMES:
        assert b[name] == pytest.approx(a[name], abs=1e-12), name


@pytest.mark.parametrize("seed", range(5))
def test_seeded_metric_ranges(seed):
    x, labels, batches = _tie_free_case(seed)
    out = _score(x, labels, batches)
    for name, value in out["metrics"].items():
        assert value is not None, (name, out["reasons"])
        low = -1.0 if name == "ari_kmeans" else -1e-9
        assert low <= value <= 1.0 + 1e-9, (name, value)


@pytest.mark.parametrize("seed", range(3))
def test_seeded_nan_and_constant_abstain(seed):
    x, labels, batches = _tie_free_case(seed)
    bad = x.copy()
    bad[seed, 0] = np.nan
    out = _score(bad, labels, batches)
    assert all(v is None for v in out["metrics"].values())
    assert all("non-finite" in r for r in out["reasons"].values())
    out = _score(np.full_like(x, 0.1 * (seed + 1)), labels, batches)
    assert all(v is None for v in out["metrics"].values())
    assert all("zero variance" in r for r in out["reasons"].values())


def test_kmeans_metrics_are_order_invariant_on_separable_data():
    x, labels = _two_blobs(n_per=60, gap=50.0, seed=2)
    perm = np.random.default_rng(3).permutation(x.shape[0])
    nmi_a, ari_a, _ = eam.kmeans_nmi_ari(x, labels, random_state=0)
    nmi_b, ari_b, _ = eam.kmeans_nmi_ari(x[perm], labels[perm], random_state=0)
    assert (nmi_a, ari_a) == pytest.approx((nmi_b, ari_b))
    assert nmi_a == pytest.approx(1.0) and ari_a == pytest.approx(1.0)


def test_rotation_and_uniform_scaling_leave_distance_based_metrics_unchanged():
    x, labels, batches = _tie_free_case(7)
    q, _ = np.linalg.qr(np.random.default_rng(8).normal(size=(x.shape[1], x.shape[1])))
    a = _score(x, labels, batches)["metrics"]
    rotated = _score(x @ q, labels, batches)["metrics"]
    scaled = _score(2.5 * x, labels, batches)["metrics"]
    for name in eam.METRIC_NAMES:
        if name in ("nmi_kmeans", "ari_kmeans"):
            continue
        assert rotated[name] == pytest.approx(a[name], abs=1e-9), ("rotation", name)
        # LISI is scale-invariant only up to the bandwidth search tolerance (1e-5 in entropy),
        # which is why the verdict dead-band DIFF_TOLERANCE is 1e-5 and not smaller
        assert scaled[name] == pytest.approx(a[name], abs=eam.DIFF_TOLERANCE), (
            "scaling",
            name,
        )
