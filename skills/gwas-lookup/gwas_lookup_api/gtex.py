"""
gtex.py — GTEx Portal eQTL lookup.

Endpoint:
  GET https://gtexportal.org/api/v2/association/singleTissueEqtl
  Params: variantId=chr{c}_{p}_{r}_{a}_b38&datasetId=gtex_v8
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .base_client import BaseClient

BASE_URL = "https://gtexportal.org/api/v2"
RATE_INTERVAL = 0.5


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def _build_gtex_id(chr: str, pos: int, ref: str, alt: str) -> str:
    """Build GTEx variant ID: chr{c}_{p}_{r}_{a}_b38."""
    return f"chr{chr}_{pos}_{ref}_{alt}_b38"


def get_eqtls(
    chr: str,
    pos: int,
    ref: str,
    alt: str,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
) -> dict:
    """Fetch single-tissue eQTL associations from GTEx."""
    client = _make_client(cache_dir, use_cache)
    variant_id = _build_gtex_id(chr, pos, ref, alt)

    try:
        data = client.get("association/singleTissueEqtl", params={
            "variantId": variant_id,
            "datasetId": "gtex_v8",
        })
    except Exception as e:
        return {"source": "gtex", "status": "error", "message": str(e)}

    if not isinstance(data, dict):
        return {"source": "gtex", "status": "error", "message": "Unexpected response format"}

    raw_eqtls = data.get("data", data.get("singleTissueEqtl", []))
    if not isinstance(raw_eqtls, list):
        raw_eqtls = []

    eqtls = []
    for e in raw_eqtls:
        eqtls.append({
            "gene_symbol": e.get("geneSymbol", e.get("gencodeId", "")),
            "gene_id": e.get("gencodeId", ""),
            "tissue": e.get("tissueSiteDetailId", ""),
            "pval": e.get("pValue", e.get("pval")),
            "nes": e.get("nes"),  # normalized effect size
            "tissue_name": e.get("tissueSiteDetail", ""),
        })

    return {
        "source": "gtex",
        "status": "ok",
        "variant_id": variant_id,
        "total_eqtls": len(eqtls),
        "eqtls": eqtls,
    }
