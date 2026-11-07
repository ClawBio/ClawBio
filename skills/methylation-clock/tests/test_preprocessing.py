"""Offline preprocessing tests (no clock-weight downloads needed)."""
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

pytest.importorskip("pyaging")

sys.path.insert(0, str(Path(__file__).parent.parent))

from methylation_clock import _aggregate_epicv2_probes, _df_to_adata

METADATA = ["gender", "tissue_type", "dataset"]


def _toy_frame() -> pd.DataFrame:
    return pd.DataFrame(
        {
            "dataset": ["GSE1", "GSE1"],
            "tissue_type": ["blood", "sperm"],
            "gender": ["F", "M"],
            "female": [1, 0],
            "cg00000029": [0.1, 0.2],
            "cg00000108_BC21": [0.4, 0.6],  # EPICv2 replicate probes
            "cg00000108_TC21": [0.6, 0.8],
        },
        index=["s1", "s2"],
    )


def test_aggregation_keeps_metadata_with_underscores():
    # PyAging >=0.5 trims column names at the first "_"; metadata such as
    # "tissue_type" must survive untouched.
    out = _aggregate_epicv2_probes(_toy_frame(), verbose=False)
    assert {"dataset", "tissue_type", "gender", "female"} <= set(out.columns)
    assert "tissue" not in out.columns
    np.testing.assert_allclose(out["cg00000108"].values, [0.5, 0.7])


def test_aggregated_frame_converts_to_anndata():
    out = _aggregate_epicv2_probes(_toy_frame(), verbose=False)
    adata = _df_to_adata(out, METADATA, "mean", verbose=False)
    assert set(METADATA) <= set(adata.obs.columns)
    assert "cg00000108" in adata.var_names
