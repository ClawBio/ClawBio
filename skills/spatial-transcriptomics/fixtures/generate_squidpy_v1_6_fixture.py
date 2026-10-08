"""Regenerate the Squidpy 1.6 cumulative co-occurrence parity fixture.

Run using the exact pinned command recorded in the adjacent JSON file.  This is
kept outside the test path because Squidpy is intentionally not a skill runtime
dependency.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import anndata as ad
import numpy as np
import scanpy
import squidpy as sq
from scipy.sparse import csr_matrix

HERE = Path(__file__).resolve().parent
OUTPUT = HERE / "squidpy_v1_6_co_occurrence.json"
COORDINATES = np.array([[0.0, 0.0], [1.0, 0.0], [4.0, 0.0]])
LABELS = ["A", "B", "B"]
INTERVAL = np.array([0.0, 1.1, 5.0])
KNN_GRAPH = np.array([[1, 2], [0, 2], [1, 3], [1, 2]], dtype=int)
MORAN_VALUES = np.array([1.0, 3.0, 2.0, 7.0])
NHOOD_LABELS = ["A", "B", "B", "B"]


def main() -> None:
    adata = ad.AnnData(np.ones((len(LABELS), 1)))
    adata.obsm["spatial"] = COORDINATES
    adata.obs["cluster"] = LABELS
    adata.obs["cluster"] = adata.obs["cluster"].astype("category")
    sq.gr.co_occurrence(
        adata,
        cluster_key="cluster",
        interval=INTERVAL,
        n_jobs=1,
        show_progress_bar=False,
    )
    graph = np.zeros((4, 4), dtype=float)
    graph[np.arange(4)[:, None], KNN_GRAPH] = 1.0
    graph_adata = ad.AnnData(MORAN_VALUES[:, None])
    graph_adata.var_names = ["gene"]
    graph_adata.obs["cluster"] = NHOOD_LABELS
    graph_adata.obs["cluster"] = graph_adata.obs["cluster"].astype("category")
    graph_adata.obsp["spatial_connectivities"] = csr_matrix(graph)
    moran = sq.gr.spatial_autocorr(
        graph_adata,
        connectivity_key="spatial_connectivities",
        genes=["gene"],
        mode="moran",
        n_perms=None,
        copy=True,
    )
    _zscore, observed = sq.gr.nhood_enrichment(
        graph_adata,
        cluster_key="cluster",
        connectivity_key="spatial",
        n_perms=50,
        seed=7,
        n_jobs=1,
        show_progress_bar=False,
        copy=True,
    )
    payload = {
        "purpose": "Offline numerical parity fixture for the defined entries of cumulative cluster co-occurrence.",
        "generated_with": {
            "python": sys.version.split()[0],
            "squidpy": sq.__version__,
            "anndata": ad.__version__,
            "scanpy": scanpy.__version__,
            "command": "uv run --no-project --python 3.12 --with 'squidpy==1.6.0' --with 'anndata==0.10.9' --with 'scanpy==1.10.4' python skills/spatial-transcriptomics/fixtures/generate_squidpy_v1_6_fixture.py",
        },
        "coordinates": COORDINATES.tolist(),
        "labels": LABELS,
        "clusters": list(adata.obs["cluster"].cat.categories),
        "squidpy_interval": INTERVAL.tolist(),
        "radii": INTERVAL[1:].tolist(),
        "axis_order": ["source_cluster", "target_cluster", "radius"],
        "scores": adata.uns["cluster_co_occurrence"]["occ"].tolist(),
        "moran_and_nhood": {
            "knn_indices": KNN_GRAPH.tolist(),
            "values": MORAN_VALUES.tolist(),
            "labels": NHOOD_LABELS,
            "moran_i": float(moran.loc["gene", "I"]),
            "nhood_clusters": list(graph_adata.obs["cluster"].cat.categories),
            "nhood_observed_counts": observed.tolist(),
            "nhood_permutations": 50,
            "nhood_seed": 7,
            "note": "Only observed counts are parity targets. Squidpy and this skill use different RNG/permutation implementations, so z-scores are not claimed bitwise equivalent.",
        },
    }
    OUTPUT.write_text(json.dumps(payload, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
