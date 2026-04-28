"""
ensembl.py — Ensembl REST API: variant info + VEP annotation.

Endpoints:
  GET /variation/human/{rsid}   → chr, pos, alleles, MAF
  GET /vep/human/id/{rsid}      → consequence, gene, impact
"""

from __future__ import annotations

from pathlib import Path
from typing import Any, Optional

from .base_client import BaseClient

BASE_URL = "https://rest.ensembl.org"
RATE_INTERVAL = 0.15  # Ensembl allows ~15 req/sec


def _make_client(cache_dir: Optional[Path], use_cache: bool) -> BaseClient:
    return BaseClient(
        base_url=BASE_URL,
        rate_interval=RATE_INTERVAL,
        cache_dir=cache_dir,
        use_cache=use_cache,
    )


def get_variant_info(rsid: str, cache_dir: Optional[Path] = None, use_cache: bool = True) -> dict:
    """Fetch variant metadata from Ensembl /variation/human/{rsid}."""
    client = _make_client(cache_dir, use_cache)
    try:
        data = client.get(f"variation/human/{rsid}", params={"content-type": "application/json"})
    except Exception as e:
        return {"source": "ensembl_variation", "status": "error", "message": str(e)}

    mappings = data.get("mappings", [])
    grch38 = None
    grch37 = None
    for m in mappings:
        asm = m.get("assembly_name", "")
        if asm == "GRCh38":
            grch38 = m
        elif asm == "GRCh37":
            grch37 = m

    alleles = data.get("mappings", [{}])[0].get("allele_string", "") if mappings else ""
    maf_data = data.get("MAF", data.get("minor_allele_freq"))
    minor_allele = data.get("minor_allele", "")

    # Population frequencies
    populations = []
    for pop in data.get("populations", []):
        populations.append({
            "population": pop.get("population", ""),
            "allele": pop.get("allele", ""),
            "frequency": pop.get("frequency"),
        })

    result = {
        "source": "ensembl_variation",
        "status": "ok",
        "rsid": rsid,
        "allele_string": alleles,
        "minor_allele": minor_allele,
        "maf": maf_data,
        "var_class": data.get("var_class", ""),
        "most_severe_consequence": data.get("most_severe_consequence", ""),
    }

    if grch38:
        result["chr"] = str(grch38.get("seq_region_name", ""))
        result["pos_grch38"] = grch38.get("start")
        result["ref_allele"] = grch38.get("allele_string", "").split("/")[0] if grch38.get("allele_string") else ""
        result["alt_alleles"] = grch38.get("allele_string", "").split("/")[1:] if grch38.get("allele_string") else []

    if grch37:
        result["pos_grch37"] = grch37.get("start")

    result["populations"] = populations
    return result


def get_vep_annotation(rsid: str, cache_dir: Optional[Path] = None, use_cache: bool = True) -> dict:
    """Fetch VEP annotation from Ensembl /vep/human/id/{rsid}."""
    client = _make_client(cache_dir, use_cache)
    try:
        data = client.get(f"vep/human/id/{rsid}", params={"content-type": "application/json"})
    except Exception as e:
        return {"source": "ensembl_vep", "status": "error", "message": str(e)}

    if not isinstance(data, list) or len(data) == 0:
        return {"source": "ensembl_vep", "status": "empty", "message": "No VEP data"}

    entry = data[0]
    consequences = []
    for tc in entry.get("transcript_consequences", []):
        consequences.append({
            "gene_symbol": tc.get("gene_symbol", ""),
            "gene_id": tc.get("gene_id", ""),
            "consequence_terms": tc.get("consequence_terms", []),
            "impact": tc.get("impact", ""),
            "biotype": tc.get("biotype", ""),
            "sift": tc.get("sift_prediction", ""),
            "polyphen": tc.get("polyphen_prediction", ""),
        })

    return {
        "source": "ensembl_vep",
        "status": "ok",
        "rsid": rsid,
        "most_severe_consequence": entry.get("most_severe_consequence", ""),
        "consequences": consequences,
    }
