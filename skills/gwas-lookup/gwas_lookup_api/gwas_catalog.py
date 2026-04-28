"""
gwas_catalog.py — NHGRI-EBI GWAS Catalog REST API.

Endpoint:
  GET /singleNucleotidePolymorphisms/{rsid}/associations
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional

from .base_client import BaseClient

BASE_URL = "https://www.ebi.ac.uk/gwas/rest/api"
RATE_INTERVAL = 0.25


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def get_associations(rsid: str, max_hits: int = 100, cache_dir: Optional[Path] = None, use_cache: bool = True) -> dict:
    """Fetch GWAS associations for a given rsID from the GWAS Catalog."""
    client = _make_client(cache_dir, use_cache)
    try:
        data = client.get(f"singleNucleotidePolymorphisms/{rsid}/associations")
    except Exception as e:
        return {"source": "gwas_catalog", "status": "error", "message": str(e)}

    embedded = data.get("_embedded", {})
    raw_assocs = embedded.get("associations", [])

    associations = []
    for a in raw_assocs[:max_hits]:
        # Extract trait
        traits = []
        for t in a.get("efoTraits", []):
            traits.append(t.get("trait", ""))

        # Extract risk allele and frequency
        risk_alleles = a.get("riskAlleles", [])
        risk_allele = risk_alleles[0].get("riskAlleleName", "") if risk_alleles else ""
        risk_freq = risk_alleles[0].get("riskFrequency", "") if risk_alleles else ""

        # Extract study info
        study = a.get("study", {})
        # The GWAS Catalog embeds study info via _links
        # We extract what's available in the association payload
        associations.append({
            "pvalue": a.get("pvalue"),
            "pvalue_mlog": a.get("pvalueMantissa"),
            "pvalue_exponent": a.get("pvalueExponent"),
            "risk_allele": risk_allele,
            "risk_frequency": risk_freq,
            "or_beta": a.get("orPerCopyNum"),
            "beta_num": a.get("betaNum"),
            "beta_direction": a.get("betaDirection"),
            "beta_unit": a.get("betaUnit"),
            "ci": a.get("range", ""),
            "traits": traits,
            "study_accession": a.get("studyAccession", ""),
        })

    return {
        "source": "gwas_catalog",
        "status": "ok",
        "rsid": rsid,
        "total_associations": len(raw_assocs),
        "associations": associations,
    }
