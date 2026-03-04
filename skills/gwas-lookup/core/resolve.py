"""
resolve.py — Variant resolution pipeline.

Takes an rsID, queries Ensembl for coordinates and alleles, then constructs
the variant identifiers needed by each downstream API module.
Falls back to PortalDev if Ensembl doesn't return coordinates.
"""

from __future__ import annotations

from pathlib import Path
from typing import Optional


def resolve_variant(rsid: str, cache_dir: Optional[Path] = None, use_cache: bool = True) -> dict:
    """
    Resolve an rsID to full variant metadata.

    Returns a dict with:
      - rsid, chr, pos_grch38, pos_grch37, ref, alt (primary alt allele)
      - variant_ids: pre-built IDs for each API
      - ensembl_info: full Ensembl response
      - vep_info: full VEP response
    """
    from ..api import ensembl, portaldev

    # Step 1: Get variant info from Ensembl
    var_info = ensembl.get_variant_info(rsid, cache_dir=cache_dir, use_cache=use_cache)

    chr_val = var_info.get("chr", "")
    pos_38 = var_info.get("pos_grch38")
    pos_37 = var_info.get("pos_grch37")
    ref = var_info.get("ref_allele", "")
    alt_alleles = var_info.get("alt_alleles", [])
    alt = alt_alleles[0] if alt_alleles else ""

    # Fallback: if Ensembl didn't give us coordinates, try PortalDev
    if not pos_38 and var_info.get("status") != "ok":
        pd_result = portaldev.resolve_rsid(rsid, cache_dir=cache_dir, use_cache=use_cache)
        if pd_result.get("status") == "ok":
            chr_val = pd_result.get("chr", chr_val)
            pos_38 = pd_result.get("start", pos_38)

    # Step 2: Get VEP annotation
    vep_info = ensembl.get_vep_annotation(rsid, cache_dir=cache_dir, use_cache=use_cache)

    # Step 3: Build variant IDs for each API
    variant_ids = {}
    if chr_val and pos_38 and ref and alt:
        variant_ids["open_targets"] = f"{chr_val}_{pos_38}_{ref}_{alt}"
        variant_ids["gtex"] = f"chr{chr_val}_{pos_38}_{ref}_{alt}_b38"
        variant_ids["pheweb"] = f"{chr_val}:{pos_38}-{ref}-{alt}"
    if chr_val and pos_37 and ref and alt:
        variant_ids["bbj"] = f"{chr_val}:{pos_37}-{ref}-{alt}"

    return {
        "rsid": rsid,
        "chr": chr_val,
        "pos_grch38": pos_38,
        "pos_grch37": pos_37,
        "ref": ref,
        "alt": alt,
        "alt_alleles": alt_alleles,
        "allele_string": var_info.get("allele_string", ""),
        "var_class": var_info.get("var_class", ""),
        "most_severe_consequence": var_info.get("most_severe_consequence", ""),
        "minor_allele": var_info.get("minor_allele", ""),
        "maf": var_info.get("maf"),
        "populations": var_info.get("populations", []),
        "variant_ids": variant_ids,
        "ensembl_info": var_info,
        "vep_info": vep_info,
    }
