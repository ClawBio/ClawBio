"""
open_targets.py — Open Targets Platform GraphQL API.

Endpoint:
  POST https://api.platform.opentargets.org/api/v4/graphql
  Query: variant(variantId: "chr_pos_ref_alt") → GWAS credible sets
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .base_client import BaseClient

BASE_URL = "https://api.platform.opentargets.org/api/v4"
RATE_INTERVAL = 0.35
# locus(variantIds: ...) narrows each set's rows to the queried variant, so
# posterior, p-value and beta are this variant's, not the set's lead variant's.
# ponytail: first page of 100 only (rs3798220 has 226); page through if the report needs them all.
CREDIBLE_SET_QUERY = """
query CredibleSetQuery($variantId: String!, $variantIds: [String!]) {
  variant(variantId: $variantId) {
    id
    credibleSets(studyTypes: [gwas], page: {index: 0, size: 100}) {
      count
      rows {
        studyId
        study { traitFromSource }
        locus(variantIds: $variantIds) {
          rows {
            posteriorProbability
            pValueMantissa
            pValueExponent
            beta
            is95CredibleSet
            is99CredibleSet
          }
        }
      }
    }
  }
}
"""


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def _build_variant_id(chr: str, pos: int, ref: str, alt: str) -> str:
    """Build Open Targets variant ID: chr_pos_ref_alt."""
    return f"{chr}_{pos}_{ref}_{alt}"


def get_credible_sets(
    chr: str,
    pos: int,
    ref: str,
    alt: str,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
) -> dict:
    """Fetch GWAS credible set membership from Open Targets."""
    client = _make_client(cache_dir, use_cache)
    variant_id = _build_variant_id(chr, pos, ref, alt)

    try:
        data = client.post("graphql", json_body={
            "query": CREDIBLE_SET_QUERY,
            "variables": {"variantId": variant_id, "variantIds": [variant_id]},
        })
    except Exception as e:
        return {"source": "open_targets_credsets", "status": "error", "message": str(e)}

    # GraphQL can answer 200 with errors; a schema change must not read as "no data".
    if data.get("errors"):
        message = "; ".join(err.get("message", "") for err in data["errors"])
        return {"source": "open_targets_credsets", "status": "error", "message": message}

    variant = (data.get("data") or {}).get("variant")
    if not variant:
        return {"source": "open_targets_credsets", "status": "empty", "message": f"No data for {variant_id}"}

    sets = variant.get("credibleSets") or {}
    credible_sets = []
    for cs in sets.get("rows") or []:
        locus = ((cs.get("locus") or {}).get("rows") or [{}])[0]
        mantissa, exponent = locus.get("pValueMantissa"), locus.get("pValueExponent")
        credible_sets.append({
            "study_id": cs.get("studyId", ""),
            "trait": (cs.get("study") or {}).get("traitFromSource", ""),
            "posterior_probability": locus.get("posteriorProbability"),
            "pval": mantissa * 10.0 ** exponent if mantissa is not None and exponent is not None else None,
            "beta": locus.get("beta"),
            "is_95_credible": locus.get("is95CredibleSet", False),
            "is_99_credible": locus.get("is99CredibleSet", False),
        })

    return {
        "source": "open_targets_credsets",
        "status": "ok",
        "variant_id": variant_id,
        "total_credible_sets": sets.get("count", len(credible_sets)),
        "credible_sets": credible_sets,
    }
