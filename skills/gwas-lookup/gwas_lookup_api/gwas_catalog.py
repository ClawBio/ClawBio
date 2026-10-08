"""
gwas_catalog.py — NHGRI-EBI GWAS Catalog REST API v2.

Endpoint:
  GET /v2/associations?rs_id={rsid}&page={n}&size={k}&sort=p_value&direction=asc

The legacy v1 API (/gwas/rest/api/singleNucleotidePolymorphisms/...) was
retired and answers HTTP 410. See
https://www.ebi.ac.uk/gwas/docs/news/rest-api-v2-migration-guide

v2 pages at 20 records by default and reports page.totalElements and
page.totalPages. We request PAGE_SIZE records per page, sorted by p-value
ascending, and keep paging until max_hits records are collected or the
pages run out. Any shortfall against totalElements is reported through
`truncated`, never hidden.
"""

from __future__ import annotations

import re
from pathlib import Path
from typing import Any, Optional

from .base_client import BaseClient

BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api/v2"
RATE_INTERVAL = 0.25  # v2 documents a 15 queries/second limit
PAGE_SIZE = 200
MAX_PAGES = 50  # hard stop against a server that never ends paging

_BETA_RE = re.compile(r"^\s*([-+]?\d*\.?\d+(?:[eE][-+]?\d+)?)\s*(.*?)\s*(increase|decrease)?\s*$")


class MalformedPayload(ValueError):
    """The v2 response did not have the documented shape."""


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def _parse_beta(raw: Any) -> tuple[Optional[float], Optional[str], Optional[str]]:
    """Split v2's free-text beta ("0.03939 unit increase") into number, unit, direction.

    v1 returned betaNum, betaUnit and betaDirection separately; v2 returns one
    string, with "-" when no beta was reported.
    """
    if not isinstance(raw, str) or raw.strip() in ("", "-", "NR"):
        return None, None, None
    m = _BETA_RE.match(raw)
    if not m:
        return None, None, None
    unit = m.group(2) or None
    return float(m.group(1)), unit, m.group(3)


def _page_records(data: Any) -> tuple[list[dict], int, int]:
    """Validate one v2 page; return (records, totalElements, totalPages)."""
    if not isinstance(data, dict) or not isinstance(data.get("page"), dict):
        raise MalformedPayload("malformed GWAS Catalog v2 payload: missing page block")
    page = data["page"]
    total = page.get("totalElements")
    pages = page.get("totalPages", 0)
    if not isinstance(total, int):
        raise MalformedPayload("malformed GWAS Catalog v2 payload: page.totalElements is not an integer")
    embedded = data.get("_embedded", {})
    records = embedded.get("associations", []) if isinstance(embedded, dict) else None
    if not isinstance(records, list) or not all(isinstance(r, dict) for r in records):
        raise MalformedPayload("malformed GWAS Catalog v2 payload: _embedded.associations is not a list of objects")
    return records, total, pages if isinstance(pages, int) else 0


def _to_association(a: dict) -> dict:
    """Map one v2 association onto the schema this skill has always returned."""
    efo = [t for t in a.get("efo_traits") or [] if isinstance(t, dict)]
    effect_alleles = a.get("snp_effect_allele") or []
    beta_num, beta_unit, beta_direction = _parse_beta(a.get("beta"))
    return {
        "pvalue": a.get("p_value"),
        # Historical key name: holds the mantissa, as it did under v1.
        "pvalue_mlog": a.get("pvalue_mantissa"),
        "pvalue_exponent": a.get("pvalue_exponent"),
        "risk_allele": effect_alleles[0] if effect_alleles else "",
        "risk_frequency": a.get("risk_frequency", ""),
        "or_beta": a.get("or_per_copy_num"),
        "beta_num": beta_num,
        "beta_direction": beta_direction,
        "beta_unit": beta_unit,
        "ci": a.get("range", ""),
        "traits": [t.get("efo_trait", "") for t in efo],
        "study_accession": a.get("accession_id", ""),
        # New in v2, additive only.
        "efo_ids": [t.get("efo_id", "") for t in efo],
        "reported_traits": a.get("reported_trait") or [],
        "mapped_genes": a.get("mapped_genes") or [],
        "pubmed_id": a.get("pubmed_id"),
        "association_id": a.get("association_id"),
    }


def get_associations(rsid: str, max_hits: int = 100, cache_dir: Optional[Path] = None, use_cache: bool = True) -> dict:
    """Fetch GWAS associations for a given rsID from the GWAS Catalog v2 API."""
    client = _make_client(cache_dir, use_cache)
    size = max(1, min(max_hits, PAGE_SIZE))
    collected: list[dict] = []
    total_available = 0
    try:
        for page_no in range(MAX_PAGES):
            data = client.get(
                "associations",
                params={
                    "rs_id": rsid,
                    "page": page_no,
                    "size": size,
                    "sort": "p_value",
                    "direction": "asc",
                },
            )
            records, total_available, total_pages = _page_records(data)
            collected.extend(records)
            if len(collected) >= max_hits or not records or page_no + 1 >= total_pages:
                break
    except Exception as e:
        return {"source": "gwas_catalog", "status": "error", "message": str(e)}

    associations = [_to_association(a) for a in collected[:max_hits]]
    return {
        "source": "gwas_catalog",
        "status": "ok",
        "rsid": rsid,
        "total_associations": len(associations),
        "total_available": total_available,
        "truncated": len(associations) < total_available,
        "associations": associations,
    }
