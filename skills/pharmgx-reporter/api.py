"""Importable API for the pharmgx-reporter skill.

Provides a programmatic interface that other skills (e.g. bio-orchestrator)
and the PatientProfile system can call without shelling out to the CLI.

Usage (from project root)::

    import sys; sys.path.insert(0, ".")
    sys.path.insert(0, "skills/pharmgx-reporter")
    from api import run

    result = run({"rs4244285": "AG", "rs3892097": "CT", ...})
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure project root is on sys.path so clawbio.common imports work
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# The skill directory uses a hyphen ("pharmgx-reporter") which is not a valid
# Python package name.  Add the skill dir to sys.path for direct import.
_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

import pharmgx_reporter as _reporter  # noqa: E402  (sibling import)

PGX_SNPS = _reporter.PGX_SNPS
GENE_DEFS = _reporter.GENE_DEFS
call_diplotype = _reporter.call_diplotype
call_phenotype = _reporter.call_phenotype
lookup_drugs = _reporter.lookup_drugs
phenotype_to_key = _reporter.phenotype_to_key


def run(genotypes: dict[str, str], options: dict | None = None) -> dict:
    """Run pharmacogenomics analysis on a genotype dict.

    Args:
        genotypes: Mapping of rsid -> genotype string (e.g. ``{"rs4244285": "AG"}``).
                   Genotypes should be uppercase; they will be normalised internally.
        options: Reserved for future use (e.g. filter to specific genes/drugs).

    Returns:
        dict with keys:
            - pgx_snps: {rsid: {genotype, gene, allele, effect}} for matched PGx variants
            - gene_profiles: {gene: {diplotype, phenotype}} for all profiled genes
            - drug_recommendations: {standard: [...], caution: [...], avoid: [...], indeterminate: [...]}
            - summary: high-level counts (genes, drugs, alerts)
    """
    options = options or {}

    # Normalise genotypes to uppercase, filter empty/missing
    snps = {
        rsid: gt.upper()
        for rsid, gt in genotypes.items()
        if gt and gt not in ("--", "00")
    }

    # Match against PGx panel
    pgx_snps: dict[str, dict] = {}
    for rsid, info in PGX_SNPS.items():
        if rsid in snps:
            pgx_snps[rsid] = {"genotype": snps[rsid], **info}

    # Call diplotypes and phenotypes for every gene
    profiles = _reporter.build_profiles(pgx_snps)

    # Look up drug recommendations
    drug_results = lookup_drugs(profiles)

    n_std = len(drug_results.get("standard", []))
    n_cau = len(drug_results.get("caution", []))
    n_avo = len(drug_results.get("avoid", []))
    n_ind = len(drug_results.get("indeterminate", []))

    return {
        "pgx_snps": pgx_snps,
        "gene_profiles": profiles,
        "marker_coverage": _reporter.marker_coverage(pgx_snps),
        "panel_limitations": _reporter.PANEL_LIMITATIONS,
        "drug_recommendations": drug_results,
        "summary": {
            "pgx_snps_found": len(pgx_snps),
            "pgx_snps_total": len(PGX_SNPS),
            "genes_profiled": len(profiles),
            "drugs_assessed": n_std + n_cau + n_avo + n_ind,
            "drugs_standard": n_std,
            "drugs_caution": n_cau,
            "drugs_avoid": n_avo,
            "drugs_indeterminate": n_ind,
        },
    }
