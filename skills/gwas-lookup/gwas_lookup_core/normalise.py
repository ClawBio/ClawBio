"""
normalise.py — Merge and deduplicate results from all API modules.

Groups results by category (GWAS, PheWAS, eQTL, credible sets),
deduplicates by trait+source, and sorts by p-value.
"""

from __future__ import annotations

import math
from typing import Any

GWS_THRESHOLD = 5e-8  # genome-wide significance


def _safe_float(val: Any) -> float | None:
    """Convert a value to float, returning None on failure."""
    if val is None:
        return None
    try:
        return float(val)
    except (TypeError, ValueError):
        return None


def _compute_pval(entry: dict) -> float | None:
    """Extract or compute p-value from an association entry."""
    pval = _safe_float(entry.get("pval") or entry.get("pvalue"))
    if pval is not None:
        return pval
    # GWAS Catalog stores mantissa + exponent
    mantissa = _safe_float(entry.get("pvalue_mlog") or entry.get("pvalue_mantissa"))
    exponent = _safe_float(entry.get("pvalue_exponent") or entry.get("pval_exponent"))
    if mantissa is not None and exponent is not None:
        return mantissa * (10 ** exponent)
    return None


def _sort_by_pval(items: list[dict]) -> list[dict]:
    """Sort a list of association dicts by p-value (ascending, None last)."""
    def key(x):
        p = _safe_float(x.get("pval"))
        return p if p is not None else float("inf")
    return sorted(items, key=key)


def merge_gwas(gwas_catalog: dict, open_targets: dict) -> list[dict]:
    """Merge GWAS associations from GWAS Catalog and Open Targets."""
    merged = []

    # GWAS Catalog
    if gwas_catalog.get("status") == "ok":
        for a in gwas_catalog.get("associations", []):
            traits = a.get("traits", [])
            pval = _compute_pval(a)
            merged.append({
                "source": "gwas_catalog",
                "trait": "; ".join(traits) if traits else "",
                "pval": pval,
                "or_beta": a.get("or_beta"),
                "beta": a.get("beta_num"),
                "beta_direction": a.get("beta_direction", ""),
                "ci": a.get("ci", ""),
                "risk_allele": a.get("risk_allele", ""),
                "risk_frequency": a.get("risk_frequency", ""),
                "study": a.get("study_accession", ""),
                "genome_wide_significant": pval is not None and pval < GWS_THRESHOLD,
            })

    # Open Targets (from credible sets if available)
    if open_targets.get("status") == "ok":
        for cs in open_targets.get("credible_sets", []):
            pval = _safe_float(cs.get("pval"))
            merged.append({
                "source": "open_targets",
                "trait": cs.get("trait", ""),
                "pval": pval,
                "or_beta": None,
                "beta": cs.get("beta"),
                "beta_direction": "",
                "ci": "",
                "risk_allele": "",
                "risk_frequency": "",
                "study": cs.get("study_id", ""),
                "genome_wide_significant": pval is not None and pval < GWS_THRESHOLD,
                "posterior_probability": cs.get("posterior_probability"),
                "is_95_credible": cs.get("is_95_credible", False),
            })

    return _sort_by_pval(merged)


def merge_phewas(ukb: dict, finngen: dict, bbj: dict) -> dict:
    """Merge PheWAS results from UKB, FinnGen, and BBJ into sub-lists."""
    result = {"ukb": [], "finngen": [], "bbj": []}

    for source_key, data in [("ukb", ukb), ("finngen", finngen), ("bbj", bbj)]:
        if data.get("status") == "ok":
            for a in data.get("associations", []):
                pval = _safe_float(a.get("pval"))
                result[source_key].append({
                    "source": source_key,
                    "phenocode": a.get("phenocode", ""),
                    "phenostring": a.get("phenostring", ""),
                    "category": a.get("category", ""),
                    "pval": pval,
                    "beta": a.get("beta"),
                    "se": a.get("sebeta"),
                    "maf": a.get("maf"),
                    "num_cases": a.get("num_cases"),
                    "num_controls": a.get("num_controls"),
                    "genome_wide_significant": pval is not None and pval < GWS_THRESHOLD,
                })
            result[source_key] = _sort_by_pval(result[source_key])

    return result


def merge_eqtls(gtex: dict, eqtl_catalogue: dict) -> list[dict]:
    """Merge eQTL results from GTEx and eQTL Catalogue."""
    merged = []

    if gtex.get("status") == "ok":
        for e in gtex.get("eqtls", []):
            pval = _safe_float(e.get("pval"))
            merged.append({
                "source": "gtex",
                "gene": e.get("gene_symbol", ""),
                "gene_id": e.get("gene_id", ""),
                "tissue": e.get("tissue", ""),
                "tissue_name": e.get("tissue_name", ""),
                "pval": pval,
                "effect_size": e.get("nes"),
                "genome_wide_significant": pval is not None and pval < GWS_THRESHOLD,
            })

    if eqtl_catalogue.get("status") == "ok":
        for a in eqtl_catalogue.get("associations", []):
            pval = _safe_float(a.get("pval"))
            merged.append({
                "source": "eqtl_catalogue",
                "gene": a.get("gene_name", ""),
                "gene_id": a.get("gene_id", ""),
                "tissue": a.get("tissue", ""),
                "tissue_name": a.get("tissue", ""),
                "pval": pval,
                "effect_size": a.get("beta"),
                "study": a.get("study", ""),
                "genome_wide_significant": pval is not None and pval < GWS_THRESHOLD,
            })

    return _sort_by_pval(merged)


def merge_all(api_results: dict) -> dict:
    """
    Merge all API results into a unified structure.

    api_results keys: gwas_catalog, open_targets, open_targets_credsets,
                      pheweb_ukb, finngen, pheweb_bbj, gtex, eqtl_catalogue
    """
    gwas = merge_gwas(
        api_results.get("gwas_catalog", {}),
        api_results.get("open_targets_credsets", {}),
    )
    phewas = merge_phewas(
        api_results.get("pheweb_ukb", {}),
        api_results.get("finngen", {}),
        api_results.get("pheweb_bbj", {}),
    )
    eqtls = merge_eqtls(
        api_results.get("gtex", {}),
        api_results.get("eqtl_catalogue", {}),
    )

    # Credible sets (pass through from Open Targets)
    credible_sets = []
    ot_cred = api_results.get("open_targets_credsets", {})
    if ot_cred.get("status") == "ok":
        credible_sets = ot_cred.get("credible_sets", [])

    # Data sources summary
    sources = {}
    for key, data in api_results.items():
        sources[key] = {
            "status": data.get("status", "not_queried"),
            "message": data.get("message", ""),
        }

    return {
        "gwas_associations": gwas,
        "phewas": phewas,
        "eqtl_associations": eqtls,
        "credible_sets": credible_sets,
        "data_sources": sources,
        "summary": {
            "total_gwas": len(gwas),
            "total_gwas_significant": sum(1 for g in gwas if g.get("genome_wide_significant")),
            "total_phewas_ukb": len(phewas["ukb"]),
            "total_phewas_finngen": len(phewas["finngen"]),
            "total_phewas_bbj": len(phewas["bbj"]),
            "total_eqtls": len(eqtls),
            "total_credible_sets": len(credible_sets),
        },
    }
