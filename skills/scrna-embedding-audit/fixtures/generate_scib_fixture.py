"""Generate the scib-metrics reference fixture for scrna-embedding-audit.

Run this in an isolated environment, never in the ClawBio environment, from
the repository root:

    uv run --no-project --python 3.12 --with 'scib-metrics==0.6.1' \
        --with 'jax[cpu]' --with scikit-learn --with numpy \
        python skills/scrna-embedding-audit/fixtures/generate_scib_fixture.py \
        --output skills/scrna-embedding-audit/fixtures

scib-metrics 0.6.1 requires Python >= 3.12 and JAX, which is why it cannot
be a ClawBio dependency. The fixture pins the metric *definitions* of that
release so that the numpy/scikit-learn re-implementation in the skill can be
checked against them.

The kNN graph is computed here with scikit-learn's exact search and handed to
scib-metrics as a ``NeighborsResults`` object. scib-metrics' own
``Benchmarker`` uses pynndescent (approximate), so values from a Benchmarker
run on the same data can differ slightly from this fixture; the fixture tests
the metric definitions, not the neighbour search.
"""

from __future__ import annotations

import argparse
import json
import platform
from datetime import UTC, datetime
from importlib.metadata import version
from pathlib import Path

import numpy as np
import scib_metrics
from scib_metrics.nearest_neighbors import NeighborsResults
from scib_metrics.utils import KMeans as ScibKMeans
from sklearn.neighbors import NearestNeighbors

SEED = 20261005
N_PER_GROUP = 50  # cells per (label, batch) cell
N_LABELS = 3
N_BATCHES = 2
N_DIMS = 10
K_LISI = 90  # scib-metrics Benchmarker: cLISI/iLISI on 90 neighbours
K_GRAPH = 15  # scib-metrics Benchmarker: graph connectivity on 15 neighbours


def make_inputs(rng: np.random.Generator) -> dict[str, np.ndarray]:
    labels = np.repeat(np.arange(N_LABELS), N_PER_GROUP * N_BATCHES)
    batches = np.tile(np.repeat(np.arange(N_BATCHES), N_PER_GROUP), N_LABELS)
    n = labels.shape[0]
    centroids = rng.normal(scale=4.0, size=(N_LABELS, N_DIMS))
    batch_shift = rng.normal(scale=1.0, size=(N_BATCHES, N_DIMS))
    noise_a = rng.normal(size=(n, N_DIMS))
    noise_b = rng.normal(size=(n, N_DIMS))
    noise_c = rng.normal(size=(n, N_DIMS))
    # "mixed": batches overlap inside each label; "batchy": batches separate inside each label;
    # "overlap": labels overlap (centroid spread close to the noise scale) so that NMI, ARI and
    # cLISI land strictly inside (0, 1) and a wrong formula cannot hide behind saturation
    x_mixed = centroids[labels] + 0.2 * batch_shift[batches] + noise_a
    x_batchy = centroids[labels] + 2.5 * batch_shift[batches] + noise_b
    x_overlap = 0.3 * centroids[labels] + 0.6 * batch_shift[batches] + noise_c
    label_names = np.array(["T cell", "B cell", "Monocyte"])[labels]
    batch_names = np.array(["donor_1", "donor_2"])[batches]
    return {
        "X_mixed": x_mixed,
        "X_batchy": x_batchy,
        "X_overlap": x_overlap,
        "labels": label_names,
        "batches": batch_names,
    }


def exact_knn(x: np.ndarray, k: int) -> NeighborsResults:
    nn = NearestNeighbors(n_neighbors=k, metric="euclidean", algorithm="brute").fit(x)
    dists, idx = nn.kneighbors(x)  # includes self at distance 0
    return NeighborsResults(indices=idx, distances=dists)


def score(x: np.ndarray, labels: np.ndarray, batches: np.ndarray) -> dict:
    knn_lisi = exact_knn(x, K_LISI)
    knn_graph = exact_knn(x, K_GRAPH)
    km = ScibKMeans(len(np.unique(labels)))
    km.fit(x)
    nmi_ari = scib_metrics.nmi_ari_cluster_labels_kmeans(x, labels)
    return {
        "silhouette_label": float(scib_metrics.silhouette_label(x, labels)),
        "isolated_labels": float(scib_metrics.isolated_labels(x, labels, batches)),
        "silhouette_batch": float(scib_metrics.silhouette_batch(x, labels, batches)),
        "bras": float(scib_metrics.bras(x, labels, batches)),
        "nmi": float(nmi_ari["nmi"]),
        "ari": float(nmi_ari["ari"]),
        "clisi": float(scib_metrics.clisi_knn(knn_lisi, labels)),
        "ilisi": float(scib_metrics.ilisi_knn(knn_lisi, batches)),
        "graph_connectivity": float(scib_metrics.graph_connectivity(knn_graph, labels)),
        "kmeans_labels": km.labels_.astype(int).tolist(),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--output", required=True, type=Path, help="fixtures directory")
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    rng = np.random.default_rng(SEED)
    inputs = make_inputs(rng)
    np.savez(args.output / "scib_fixture_input.npz", **inputs)

    expected = {
        "generator": "generate_scib_fixture.py",
        "generated_on": datetime.now(UTC).date().isoformat(),
        "python": platform.python_version(),
        "packages": {
            name: version(name)
            for name in (
                "scib-metrics",
                "jax",
                "jaxlib",
                "numpy",
                "scikit-learn",
                "scipy",
            )
        },
        "seed": SEED,
        "n_neighbors": {"lisi": K_LISI, "graph_connectivity": K_GRAPH},
        "perplexity": K_LISI // 3,
        "embeddings": {
            key: score(inputs[key], inputs["labels"], inputs["batches"])
            for key in ("X_mixed", "X_batchy", "X_overlap")
        },
    }
    out = args.output / "scib_metrics_0.6.1_expected.json"
    out.write_text(json.dumps(expected, indent=2) + "\n")
    for key, vals in expected["embeddings"].items():
        print(key, {k: round(v, 6) for k, v in vals.items() if k != "kmeans_labels"})
    print("wrote", out)


if __name__ == "__main__":
    main()
