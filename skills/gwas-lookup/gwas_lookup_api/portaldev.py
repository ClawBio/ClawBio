"""
portaldev.py — LocusZoom / UM PortalDev rsID → position fallback.

Endpoint:
  GET https://portaldev.sph.umich.edu/api/v1/annotation/omnisearch/?q={rsid}&build=GRCh38
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .base_client import BaseClient

BASE_URL = "https://portaldev.sph.umich.edu/api/v1"
RATE_INTERVAL = 0.3


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def resolve_rsid(
    rsid: str,
    build: str = "GRCh38",
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
) -> dict:
    """Resolve an rsID to genomic coordinates via LocusZoom PortalDev."""
    client = _make_client(cache_dir, use_cache)

    try:
        data = client.get("annotation/omnisearch/", params={"q": rsid, "build": build})
    except Exception as e:
        return {"source": "portaldev", "status": "error", "message": str(e)}

    if not isinstance(data, dict):
        return {"source": "portaldev", "status": "error", "message": "Unexpected response format"}

    results = data.get("data", [])
    if not results:
        return {"source": "portaldev", "status": "empty", "message": f"No results for {rsid}"}

    hit = results[0]
    return {
        "source": "portaldev",
        "status": "ok",
        "rsid": rsid,
        "build": build,
        "chr": str(hit.get("chrom", "")),
        "start": hit.get("start"),
        "end": hit.get("end"),
        "term": hit.get("term", ""),
    }
