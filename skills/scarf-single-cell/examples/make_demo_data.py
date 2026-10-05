#!/usr/bin/env python3
"""Write the synthetic demo dataset for scarf-single-cell as a raw-count H5AD.

2,000 cells x 1,500 genes, Poisson counts with four planted groups (50 marker genes each),
two donors, 10 mitochondrial-like and 20 ribosomal-like genes, and a held-out
`author_cell_type` column that records the planted group. Seeded, so every run writes the
same matrix. No real cells, donors or patients are involved.

Usage:
    python skills/scarf-single-cell/examples/make_demo_data.py OUT.h5ad [--seed 0]

Needs numpy, pandas, scipy and anndata (all installed with scarf[extra]); Scarf itself is
not needed to write the file.
"""

from __future__ import annotations

import argparse
from pathlib import Path

N_CELLS = 2000
N_GENES = 1500
N_GROUPS = 4
MARKERS_PER_GROUP = 50
MARKER_FOLD = 8.0


def make_demo_h5ad(path: Path | str, seed: int = 0) -> Path:
    """Write the demo AnnData (raw integer counts in X) to `path` and return the path."""
    import anndata as ad
    import numpy as np
    import pandas as pd
    import scipy.sparse as sp

    rng = np.random.default_rng(seed)
    genes = (
        [f"GENE{i}" for i in range(N_GENES - 30)]
        + [f"MT-G{i}" for i in range(10)]
        + [f"RPL{i}" for i in range(20)]
    )
    group = rng.integers(0, N_GROUPS, N_CELLS)
    base = rng.gamma(0.3, 1.0, N_GENES)
    mu = np.tile(base, (N_CELLS, 1))
    for g in range(N_GROUPS):
        start = g * MARKERS_PER_GROUP
        mu[group == g, start : start + MARKERS_PER_GROUP] *= MARKER_FOLD
    library = rng.lognormal(0, 0.3, N_CELLS)[:, None]
    counts = sp.csr_matrix(rng.poisson(mu * library).astype(np.int32))
    # Plain object strings: pandas 3 may infer Arrow-backed strings, which older anndata
    # releases cannot write.
    def categorical(values: list[str]) -> pd.Categorical:
        return pd.Categorical(values, categories=pd.Index(sorted(set(values)), dtype=object))

    obs = pd.DataFrame(
        {
            "donor_id": categorical(np.where(np.arange(N_CELLS) % 2, "D1", "D2").tolist()),
            "author_cell_type": categorical([f"type_{g}" for g in group]),
        },
        index=pd.Index([f"cell{i}" for i in range(N_CELLS)], dtype=object),
    )
    var = pd.DataFrame(index=pd.Index(genes, dtype=object))
    path = Path(path)
    path.parent.mkdir(parents=True, exist_ok=True)
    ad.AnnData(X=counts, obs=obs, var=var).write_h5ad(path)
    return path


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("output", type=Path, help="Path of the H5AD file to write")
    parser.add_argument("--seed", type=int, default=0)
    args = parser.parse_args()
    print(make_demo_h5ad(args.output, seed=args.seed))


if __name__ == "__main__":
    main()
