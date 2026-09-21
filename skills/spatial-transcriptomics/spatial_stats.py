"""Dependency-light spatial statistics used by the Visium skill.

The helpers in this module operate on a coordinate kNN graph rather than the
expression graph used for Leiden clustering.  They deliberately use NumPy only
so importing the command-line module does not require Scanpy or scikit-learn.
"""

from __future__ import annotations

import numpy as np


def _coordinates(coords: np.ndarray) -> np.ndarray:
    """Return validated two-dimensional spatial coordinates."""
    values = np.asarray(coords, dtype=float)
    if values.ndim != 2 or values.shape[1] != 2:
        raise ValueError("coords must be a finite array with shape (n_spots, 2)")
    if values.shape[0] < 2:
        raise ValueError("at least two spatial coordinates are required")
    if not np.isfinite(values).all():
        raise ValueError("coords must contain only finite values")
    return values


def _graph(knn_idx: np.ndarray, n_obs: int | None = None) -> np.ndarray:
    """Validate and return a rectangular directed kNN index matrix."""
    graph = np.asarray(knn_idx)
    if graph.ndim != 2 or graph.shape[0] < 2 or graph.shape[1] < 1:
        raise ValueError(
            "knn_idx must have shape (n_spots, n_neighbors) with n_neighbors >= 1"
        )
    if not np.issubdtype(graph.dtype, np.integer):
        raise ValueError("knn_idx must contain integer indices")
    if n_obs is not None and graph.shape[0] != n_obs:
        raise ValueError("knn_idx and values/labels must have the same number of spots")
    if (graph < 0).any() or (graph >= graph.shape[0]).any():
        raise ValueError("knn_idx contains an out-of-range neighbour index")
    if np.any(graph == np.arange(graph.shape[0])[:, None]):
        raise ValueError("knn_idx must exclude self-neighbours")
    return graph.astype(int, copy=False)


def _labels(labels: np.ndarray, n_obs: int) -> np.ndarray:
    values = np.asarray(labels).astype(str)
    if values.ndim != 1 or len(values) != n_obs:
        raise ValueError("labels must be one-dimensional and match the number of spots")
    return values


def knn_indices(coords: np.ndarray, n_neighbors: int) -> np.ndarray:
    """Build a directed kNN graph, excluding each query point by its identity.

    Excluding ``candidate == query_index`` after the nearest-neighbour query,
    rather than dropping its first result, remains correct when two spots have
    identical coordinates.
    """
    points = _coordinates(coords)
    if not isinstance(n_neighbors, (int, np.integer)) or n_neighbors < 1:
        raise ValueError("n_neighbors must be an integer >= 1")
    from sklearn.neighbors import NearestNeighbors

    n_obs = points.shape[0]
    k = min(int(n_neighbors), n_obs - 1)
    finder = NearestNeighbors(n_neighbors=k + 1, algorithm="auto")
    finder.fit(points)
    candidate_indices = finder.kneighbors(points, return_distance=False)
    indices = np.empty((n_obs, k), dtype=int)
    for row in range(n_obs):
        neighbours = candidate_indices[row]
        neighbours = neighbours[neighbours != row]
        # With tied duplicate coordinates sklearn can omit the query point from
        # the k+1 candidates.  Either way, k non-self candidates are available.
        if len(neighbours) < k:
            raise RuntimeError(
                "nearest-neighbour query did not return enough non-self neighbours"
            )
        indices[row] = neighbours[:k]
    return indices


def moran_i(values: np.ndarray, knn_idx: np.ndarray) -> float:
    """Compute row-standardised kNN Moran's I.

    A constant feature has no variance and therefore no defined Moran's I; this
    function returns ``nan`` rather than reporting it as spatially neutral.
    """
    x = np.asarray(values, dtype=float)
    if x.ndim != 1:
        raise ValueError("values must be one-dimensional")
    if not np.isfinite(x).all():
        raise ValueError("values must contain only finite values")
    graph = _graph(knn_idx, len(x))
    z = x - x.mean()
    denominator = float(np.dot(z, z))
    if denominator == 0.0:
        return float("nan")
    spatial_lag = z[graph].mean(axis=1)
    return float(np.dot(z, spatial_lag) / denominator)


def neighbor_pair_counts(
    knn_idx: np.ndarray, labels: np.ndarray, clusters: list[str]
) -> np.ndarray:
    """Count directed ``(source_cluster, target_cluster)`` graph edges."""
    graph = _graph(knn_idx)
    label_values = _labels(labels, graph.shape[0])
    cluster_values = [str(cluster) for cluster in clusters]
    if len(cluster_values) != len(set(cluster_values)):
        raise ValueError("clusters must not contain duplicates")
    index = {cluster: position for position, cluster in enumerate(cluster_values)}
    unknown = set(label_values) - set(index)
    if unknown:
        raise ValueError(
            f"labels include cluster(s) not listed in clusters: {sorted(unknown)!r}"
        )

    source_codes = np.asarray([index[label] for label in label_values], dtype=int)
    target_codes = source_codes[graph]
    counts = np.zeros((len(cluster_values), len(cluster_values)), dtype=float)
    np.add.at(
        counts,
        (np.repeat(source_codes, graph.shape[1]), target_codes.reshape(-1)),
        1.0,
    )
    return counts


def nhood_enrichment(
    knn_idx: np.ndarray,
    labels: np.ndarray,
    n_perms: int = 50,
    seed: int = 7,
) -> tuple[list[str], np.ndarray, np.ndarray]:
    """Compare observed directed neighbour counts with shuffled-label counts.

    Entries with a zero permutation standard deviation are returned as ``nan``:
    their z-score is undefined, rather than evidence for no enrichment.
    """
    graph = _graph(knn_idx)
    label_values = _labels(labels, graph.shape[0])
    if not isinstance(n_perms, (int, np.integer)) or n_perms < 1:
        raise ValueError("n_perms must be an integer >= 1")
    clusters = sorted(set(label_values.tolist()))
    observed = neighbor_pair_counts(graph, label_values, clusters)
    rng = np.random.default_rng(seed)
    null = np.stack(
        [
            neighbor_pair_counts(graph, rng.permutation(label_values), clusters)
            for _ in range(int(n_perms))
        ]
    )
    mean = null.mean(axis=0)
    std = null.std(axis=0)
    zscores = np.full(observed.shape, np.nan, dtype=float)
    np.divide(observed - mean, std, out=zscores, where=std > 0.0)
    return clusters, zscores, observed


def _positive_pair_distances(points: np.ndarray) -> np.ndarray:
    """Return upper-triangle positive distances without an ``n x n x 2`` array."""
    n_obs = points.shape[0]
    chunks: list[np.ndarray] = []
    for source in range(n_obs - 1):
        distances = np.sqrt(
            np.sum((points[source + 1 :] - points[source]) ** 2, axis=1)
        )
        positive = distances[distances > 0.0]
        if positive.size:
            chunks.append(positive)
    if not chunks:
        return np.empty(0, dtype=float)
    return np.concatenate(chunks)


def _co_occurrence_radii(points: np.ndarray, n_bins: int) -> np.ndarray:
    if not isinstance(n_bins, (int, np.integer)) or n_bins < 1:
        raise ValueError("n_bins must be an integer >= 1")
    distances = _positive_pair_distances(points)
    if distances.size == 0:
        return np.zeros(int(n_bins), dtype=float)
    # These are cumulative distance thresholds, not disjoint bin edges.  They
    # are a useful local default, but are not Squidpy's automatic interval rule.
    return np.quantile(distances, np.linspace(1.0 / n_bins, 1.0, int(n_bins)))


def co_occurrence(
    coords: np.ndarray,
    labels: np.ndarray,
    n_bins: int = 6,
    *,
    radii: np.ndarray | list[float] | None = None,
) -> tuple[np.ndarray, np.ndarray, list[str]]:
    """Compute cumulative cluster co-occurrence probability ratios.

    For each radius ``r``, the score is
    ``P(target | source, 0 < distance <= r) / P(target | eligible pair, r)``.
    This is the cumulative-radius and eligible-pair marginal estimator used by
    Squidpy v1.6.0.  Default radii are six positive-distance quantiles; they are
    intentionally not presented as Squidpy's automatic interval selection.

    The returned tensor is indexed ``[source_cluster, target_cluster, radius]``.
    Duplicate coordinates are excluded because their distance is zero, matching
    Squidpy's ``pw_dist > 0`` eligibility condition.  Scores with no eligible
    source pairs, no eligible target marginal, or an empty radius are ``nan``:
    those ratios are undefined rather than zero co-occurrence.
    """
    points = _coordinates(coords)
    label_values = _labels(labels, points.shape[0])
    clusters = sorted(set(label_values.tolist()))
    codes = np.searchsorted(np.asarray(clusters, dtype=str), label_values)
    n_clusters = len(clusters)

    if radii is None:
        threshold_values = _co_occurrence_radii(points, n_bins)
    else:
        threshold_values = np.asarray(radii, dtype=float)
        if threshold_values.ndim != 1 or threshold_values.size == 0:
            raise ValueError("radii must be a non-empty one-dimensional array")
        if not np.isfinite(threshold_values).all() or (threshold_values <= 0.0).any():
            raise ValueError("radii must contain only finite positive values")
        if np.any(np.diff(threshold_values) < 0.0):
            raise ValueError("radii must be sorted in increasing order")

    scores = np.full(
        (n_clusters, n_clusters, len(threshold_values)), np.nan, dtype=float
    )
    if not np.any(threshold_values > 0.0):
        return scores, threshold_values, clusters

    # A source chunk has shape (chunk, n_spots), avoiding the previous
    # (n_spots, n_spots, 2) coordinate-difference allocation.
    chunk_size = max(1, min(512, points.shape[0]))
    pair_counts = np.zeros_like(scores)
    for start in range(0, points.shape[0], chunk_size):
        stop = min(start + chunk_size, points.shape[0])
        source_points = points[start:stop]
        distances = np.sqrt(
            np.sum((source_points[:, None, :] - points[None, :, :]) ** 2, axis=2)
        )
        source_codes = codes[start:stop]
        for radius_index, radius in enumerate(threshold_values):
            if radius <= 0.0:
                continue
            eligible = (distances > 0.0) & (distances <= radius)
            for source_code in range(n_clusters):
                source_rows = source_codes == source_code
                if not source_rows.any():
                    continue
                target_codes = np.broadcast_to(codes[None, :], eligible.shape)
                targets = target_codes[
                    np.broadcast_to(source_rows[:, None], eligible.shape) & eligible
                ]
                if targets.size:
                    pair_counts[source_code, :, radius_index] += np.bincount(
                        targets, minlength=n_clusters
                    )

    for radius_index in range(len(threshold_values)):
        counts = pair_counts[:, :, radius_index]
        total = counts.sum()
        if total == 0.0:
            continue
        target_marginal = counts.sum(axis=0) / total
        source_totals = counts.sum(axis=1)
        valid = (source_totals[:, None] > 0.0) & (target_marginal[None, :] > 0.0)
        conditional = np.divide(
            counts,
            source_totals[:, None],
            out=np.zeros_like(counts),
            where=source_totals[:, None] > 0.0,
        )
        np.divide(
            conditional,
            target_marginal[None, :],
            out=scores[:, :, radius_index],
            where=valid,
        )
    return scores, threshold_values, clusters
