"""
finngen.py — FinnGen r12 PheWAS lookup.

Endpoint:
  GET https://r12.finngen.fi/api/variant/{chr}:{pos}-{ref}-{alt}
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .base_client import BaseClient

BASE_URL = "https://r12.finngen.fi"
RATE_INTERVAL = 0.5


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def _build_variant_str(chr: str, pos: int, ref: str, alt: str) -> str:
    return f"{chr}:{pos}-{ref}-{alt}"


def get_phewas(
    chr: str,
    pos: int,
    ref: str,
    alt: str,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
) -> dict:
    """Fetch PheWAS associations from FinnGen r12."""
    client = _make_client(cache_dir, use_cache)
    variant_str = _build_variant_str(chr, pos, ref, alt)

    try:
        data = client.get(f"api/variant/{variant_str}")
    except Exception as e:
        return {"source": "finngen", "status": "error", "message": str(e)}

    if not isinstance(data, dict):
        return {"source": "finngen", "status": "error", "message": "Unexpected response format"}

    phenos = data.get("phenos", [])
    associations = []
    for p in phenos:
        associations.append({
            "phenocode": p.get("phenocode", ""),
            "phenostring": p.get("phenostring", ""),
            "category": p.get("category", ""),
            "pval": p.get("pval"),
            "beta": p.get("beta"),
            "sebeta": p.get("sebeta"),
            "maf": p.get("maf"),
            "num_cases": p.get("num_cases"),
            "num_controls": p.get("num_controls"),
        })

    return {
        "source": "finngen",
        "status": "ok",
        "variant": variant_str,
        "total_associations": len(associations),
        "associations": associations,
    }
