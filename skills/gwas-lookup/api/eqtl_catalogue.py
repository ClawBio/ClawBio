"""
eqtl_catalogue.py — EBI eQTL Catalogue API.

Endpoint:
  GET https://www.ebi.ac.uk/eqtl/api/v3/associations?rsid={rsid}
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .base_client import BaseClient

BASE_URL = "https://www.ebi.ac.uk/eqtl/api/v3"
RATE_INTERVAL = 0.3


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def get_associations(
    rsid: str,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
) -> dict:
    """Fetch eQTL associations from the EBI eQTL Catalogue."""
    client = _make_client(cache_dir, use_cache)

    try:
        data = client.get("associations", params={"rsid": rsid, "size": 100})
    except Exception as e:
        return {"source": "eqtl_catalogue", "status": "error", "message": str(e)}

    if not isinstance(data, dict):
        return {"source": "eqtl_catalogue", "status": "error", "message": "Unexpected response format"}

    # The API may return _embedded.associations or a direct list
    embedded = data.get("_embedded", {})
    raw_assocs = embedded.get("associations", data.get("associations", []))
    if not isinstance(raw_assocs, list):
        raw_assocs = []

    associations = []
    for a in raw_assocs:
        associations.append({
            "gene_id": a.get("gene_id", ""),
            "gene_name": a.get("gene_name", a.get("molecular_trait_id", "")),
            "tissue": a.get("tissue_label", a.get("qtl_group", "")),
            "study": a.get("study_id", ""),
            "pval": a.get("pvalue", a.get("neg_log10_pvalue")),
            "beta": a.get("beta"),
            "se": a.get("se"),
            "maf": a.get("maf"),
            "dataset": a.get("dataset_id", ""),
        })

    return {
        "source": "eqtl_catalogue",
        "status": "ok",
        "rsid": rsid,
        "total_associations": len(associations),
        "associations": associations,
    }
