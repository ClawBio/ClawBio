"""
open_targets.py — Open Targets Genetics GraphQL API.

Endpoint:
  POST https://api.platform.opentargets.org/api/v4/graphql
  Query: variant(variantId: "chr_pos_ref_alt") → credible sets, V2G scores
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .base_client import BaseClient

BASE_URL = "https://api.platform.opentargets.org/api/v4"
RATE_INTERVAL = 0.35

VARIANT_QUERY = """
query VariantQuery($variantId: String!) {
  variant(variantId: $variantId) {
    id
    rsId
    chromosome
    position
    refAllele
    altAllele
    nearestGene {
      id
      approvedSymbol
    }
    nearestGeneDistance
    mostSevereConsequence
    gnomadNFE
    gnomadAFR
    gnomadEAS
    gnomadAMR
    gnomadFIN
  }
}
"""

CREDIBLE_SET_QUERY = """
query CredibleSetQuery($variantId: String!) {
  variant(variantId: $variantId) {
    id
    credibleSets {
      study {
        studyId
        traitReported
      }
      posteriorProbability
      pval
      beta
      is95CredibleSet
      is99CredibleSet
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


def get_variant(
    chr: str,
    pos: int,
    ref: str,
    alt: str,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
) -> dict:
    """Fetch variant info from Open Targets GraphQL API."""
    client = _make_client(cache_dir, use_cache)
    variant_id = _build_variant_id(chr, pos, ref, alt)

    try:
        data = client.post("graphql", json_body={
            "query": VARIANT_QUERY,
            "variables": {"variantId": variant_id},
        })
    except Exception as e:
        return {"source": "open_targets", "status": "error", "message": str(e)}

    variant = data.get("data", {}).get("variant")
    if not variant:
        return {"source": "open_targets", "status": "empty", "message": f"No data for {variant_id}"}

    nearest = variant.get("nearestGene") or {}
    return {
        "source": "open_targets",
        "status": "ok",
        "variant_id": variant_id,
        "rsid": variant.get("rsId", ""),
        "nearest_gene": nearest.get("approvedSymbol", ""),
        "nearest_gene_distance": variant.get("nearestGeneDistance"),
        "consequence": variant.get("mostSevereConsequence", ""),
        "population_frequencies": {
            "NFE": variant.get("gnomadNFE"),
            "AFR": variant.get("gnomadAFR"),
            "EAS": variant.get("gnomadEAS"),
            "AMR": variant.get("gnomadAMR"),
            "FIN": variant.get("gnomadFIN"),
        },
    }


def get_credible_sets(
    chr: str,
    pos: int,
    ref: str,
    alt: str,
    cache_dir: Optional[Path] = None,
    use_cache: bool = True,
) -> dict:
    """Fetch credible set membership from Open Targets."""
    client = _make_client(cache_dir, use_cache)
    variant_id = _build_variant_id(chr, pos, ref, alt)

    try:
        data = client.post("graphql", json_body={
            "query": CREDIBLE_SET_QUERY,
            "variables": {"variantId": variant_id},
        })
    except Exception as e:
        return {"source": "open_targets_credsets", "status": "error", "message": str(e)}

    variant = data.get("data", {}).get("variant")
    if not variant:
        return {"source": "open_targets_credsets", "status": "empty", "message": f"No data for {variant_id}"}

    raw_sets = variant.get("credibleSets") or []
    credible_sets = []
    for cs in raw_sets:
        study = cs.get("study") or {}
        credible_sets.append({
            "study_id": study.get("studyId", ""),
            "trait": study.get("traitReported", ""),
            "posterior_probability": cs.get("posteriorProbability"),
            "pval": cs.get("pval"),
            "beta": cs.get("beta"),
            "is_95_credible": cs.get("is95CredibleSet", False),
            "is_99_credible": cs.get("is99CredibleSet", False),
        })

    return {
        "source": "open_targets_credsets",
        "status": "ok",
        "variant_id": variant_id,
        "credible_sets": credible_sets,
    }
