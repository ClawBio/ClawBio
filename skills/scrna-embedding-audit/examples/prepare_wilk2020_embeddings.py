"""Prepare the public validation input for scrna-embedding-audit (Wilk et al. 2020 PBMC).

Input: the CELLxGENE Discover h5ad of Wilk et al., Nat Med 2020 (doi:10.1038/s41591-020-0944-y),
downloaded as described in ``public_validation.md`` (no download happens here). The file keeps
raw UMI counts in ``raw.X`` and author cell-type labels in ``obs["cell_type"]``; donors are in
``obs["donor_id"]``.

Output: an h5ad with raw counts in ``X`` and ``layers["counts"]``, the label and batch columns,
and two embeddings computed here with scanpy only:

- ``X_combat``: 2000 Seurat HVGs (batch-aware), ComBat by donor, scale(max 10), 50-PC PCA.
- ``X_pca_hvg_batchaware``: the same without ComBat (batch-aware HVG selection is the only
  batch information used), so the audit can show what ComBat adds on top of it.

Usage:
    uv run --locked --all-extras python skills/scrna-embedding-audit/examples/prepare_wilk2020_embeddings.py \
        --input wilk2020_covid_pbmc.h5ad --output wilk2020_for_audit.h5ad
"""

from __future__ import annotations

import argparse
from pathlib import Path

import anndata as ad
import numpy as np
import pandas as pd
import scanpy as sc
from scipy import sparse


def _object_strings(adata: ad.AnnData) -> ad.AnnData:
    """Cast string-like obs/var columns and indices to object so anndata 0.12 can write them."""
    for frame in (adata.obs, adata.var):
        for column in frame.columns:
            dtype = frame[column].dtype
            if isinstance(dtype, pd.CategoricalDtype):
                cat = frame[column].cat
                frame[column] = cat.set_categories(cat.categories.astype(object))
            elif isinstance(dtype, pd.StringDtype) or dtype == "str":
                frame[column] = frame[column].astype(object)
    adata.obs_names = adata.obs_names.astype(object)
    adata.var_names = adata.var_names.astype(object)
    return adata


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    parser.add_argument("--n-top-hvg", type=int, default=2000)
    parser.add_argument("--n-pcs", type=int, default=50)
    parser.add_argument("--random-state", type=int, default=0)
    args = parser.parse_args()

    source = ad.read_h5ad(args.input)
    counts = sparse.csr_matrix(source.raw.X).astype(np.float32)
    values = counts.data[:200_000]
    assert np.allclose(values, np.rint(values), rtol=0, atol=1e-6), (
        "raw.X is not integer-valued"
    )
    var = pd.DataFrame(index=pd.Index(source.raw.var_names.astype(str), dtype=object))
    if "feature_name" in source.raw.var.columns:
        var["feature_name"] = source.raw.var["feature_name"].astype(str).to_numpy()
    obs = pd.DataFrame(
        {
            "cell_type": source.obs["cell_type"].astype(str).to_numpy(),
            "donor_id": source.obs["donor_id"].astype(str).to_numpy(),
            "disease": source.obs["disease"].astype(str).to_numpy(),
        },
        index=pd.Index(source.obs_names.astype(str), dtype=object),
    )
    adata = ad.AnnData(X=counts.copy(), obs=obs, var=var)
    adata.layers["counts"] = counts.copy()

    work = ad.AnnData(
        X=counts.copy(), obs=obs.copy(), var=pd.DataFrame(index=var.index)
    )
    sc.pp.filter_genes(work, min_cells=3)
    sc.pp.normalize_total(work, target_sum=1e4)
    sc.pp.log1p(work)
    sc.pp.highly_variable_genes(
        work, n_top_genes=args.n_top_hvg, flavor="seurat", batch_key="donor_id"
    )
    work = work[:, work.var["highly_variable"].to_numpy()].copy()
    work.X = np.asarray(
        work.X.toarray() if sparse.issparse(work.X) else work.X, dtype=np.float64
    )

    plain = work.copy()
    plain.X = np.array(plain.X, dtype=np.float64, copy=True)
    sc.pp.scale(plain, max_value=10)
    sc.pp.pca(
        plain,
        n_comps=args.n_pcs,
        svd_solver="arpack",
        random_state=args.random_state,
        dtype="float64",
    )
    adata.obsm["X_pca_hvg_batchaware"] = np.ascontiguousarray(
        plain.obsm["X_pca"], dtype=np.float64
    )

    combat = work.copy()
    sc.pp.combat(combat, key="donor_id")
    # pandas 3 hands back a read-only array; scale() needs a writable copy
    combat.X = np.array(combat.X, dtype=np.float64, copy=True)
    sc.pp.scale(combat, max_value=10)
    sc.pp.pca(
        combat,
        n_comps=args.n_pcs,
        svd_solver="arpack",
        random_state=args.random_state,
        dtype="float64",
    )
    adata.obsm["X_combat"] = np.ascontiguousarray(
        combat.obsm["X_pca"], dtype=np.float64
    )

    adata.uns["prepared_by"] = {
        "script": "skills/scrna-embedding-audit/examples/prepare_wilk2020_embeddings.py",
        "n_top_hvg": args.n_top_hvg,
        "n_pcs": args.n_pcs,
        "random_state": args.random_state,
        "hvg_genes_used": int(work.n_vars),
        "scanpy": __import__("importlib.metadata").metadata.version("scanpy"),
    }
    ad.settings.allow_write_nullable_strings = True
    # anndata 0.12 + pandas 3 (pyarrow installed): converting string columns to
    # categoricals on write produces Arrow-backed categories that have no h5ad
    # writer, so the string columns are written as plain string arrays instead
    _object_strings(adata).write_h5ad(
        args.output, convert_strings_to_categoricals=False
    )
    print(
        f"wrote {args.output}: {adata.n_obs} cells x {adata.n_vars} genes; obsm {list(adata.obsm.keys())}"
    )


if __name__ == "__main__":
    main()
