"""Importable API for the gwas-prs skill.

Allows other skills and the orchestrator to call PRS calculation
programmatically without shelling out to the CLI.

Usage:
    import importlib, sys, pathlib
    _skill_dir = pathlib.Path("<project_root>/skills/gwas-prs")
    if str(_skill_dir) not in sys.path:
        sys.path.insert(0, str(_skill_dir))
    from api import run

    result = run(
        genotypes={"rs7903146": "CT", "rs1801282": "CG", ...},
        options={"curated_panel_id": "CLAWBIO-T2D-8", "build": "GRCh37"},
    )
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure project root is importable
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# The skill directory uses a hyphen ("gwas-prs") which is not a valid
# Python package name, so we load the engine module via importlib.
_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

import gwas_prs as _engine  # noqa: E402  (sibling module in same dir)


def run(genotypes: dict[str, str], options: dict | None = None) -> dict:
    """Run PRS calculation on a genotype dict.

    Args:
        genotypes: Mapping of rsid -> genotype string (e.g. {"rs7903146": "CT"}).
        options: Optional settings dict. Recognised keys:
            - curated_panel_id (str): Specific bundled panel ID (for example
              ``CLAWBIO-T2D-8``). If omitted, all curated panels are used.
            - pgs_id (str): Deprecated benchmark-compatibility alias. Only the
              explicitly pinned alias in ``LEGACY_PGS_PANEL_COMPAT`` is
              accepted; this local API does not download Catalog scores.
            - build (str): Genome build, "GRCh37" (default) or "GRCh38".
            - min_overlap (float): Minimum SNP overlap fraction (default 0.5).

    Returns:
        Dict with keys:
            - results: list of per-score dicts (score_id, pgs_id, trait, raw_score,
              percentile, risk_category, overlap_fraction, ...)
            - scores_calculated: int
            - disclaimer: str
    """
    options = options or {}
    build = options.get("build", "GRCh37")
    curated_panel_id = options.get("curated_panel_id")
    pgs_id = options.get("pgs_id")
    min_overlap = options.get("min_overlap", 0.5)

    # Determine which scoring files to use
    scoring_entries: list[dict] = []

    if curated_panel_id and pgs_id:
        raise ValueError("provide curated_panel_id or pgs_id, not both")

    legacy_compat = False
    output_pgs_id = None
    if curated_panel_id:
        panel_id = curated_panel_id.strip().upper()
        if panel_id not in _engine.CURATED_SCORES:
            raise ValueError(f"unknown curated panel ID: {panel_id}")
        ids_to_score = [panel_id]
    elif pgs_id:
        normalized_pgs_id = pgs_id.strip().upper()
        if not normalized_pgs_id.startswith("PGS"):
            normalized_pgs_id = "PGS" + normalized_pgs_id.lstrip("0")
        panel_id = _engine.LEGACY_PGS_PANEL_COMPAT.get(normalized_pgs_id)
        if panel_id is None:
            raise ValueError(
                "the local gwas-prs API accepts canonical curated_panel_id; "
                f"{normalized_pgs_id} is not a supported compatibility alias"
            )
        ids_to_score = [panel_id]
        legacy_compat = True
        output_pgs_id = normalized_pgs_id
    else:
        # Default: all curated scores
        ids_to_score = list(_engine.CURATED_SCORES.keys())

    for panel_id in ids_to_score:
        fpath = _engine.curated_panel_path(panel_id, build)
        if not fpath.exists():
            continue

        meta = _engine.CURATED_SCORES[panel_id]
        scoring_entries.append({
            "score_id": panel_id,
            "pgs_id": output_pgs_id if legacy_compat else None,
            "trait": meta.get("trait", "Unknown"),
            "filepath": fpath,
            "metadata": {
                "publication": meta.get("publication", ""),
                "variants_count": meta.get("variants_count", 0),
            },
        })

    # Score each entry
    all_results: list[dict] = []
    for sf in scoring_entries:
        scoring_variants = _engine.parse_scoring_file(sf["filepath"])
        if not scoring_variants:
            continue

        prs = _engine.calculate_prs(genotypes, scoring_variants)

        if prs["overlap_fraction"] < min_overlap:
            continue

        pct_info = _engine.estimate_percentile(
            prs["raw_score"], sf["score_id"], scoring_variants
        )

        all_results.append({
            "score_id": sf["score_id"],
            "pgs_id": sf["pgs_id"],
            "curated_demo_panel": True,
            "curated_panel_id": sf["score_id"],
            "legacy_pgs_id": _engine.CURATED_SCORES[
                sf["score_id"]
            ]["legacy_pgs_id"],
            "legacy_pgs_compatibility": legacy_compat,
            "pgs_catalog_id": _engine.CURATED_SCORES[
                sf["score_id"]
            ].get("pgs_catalog_id"),
            "trait": sf["trait"],
            "raw_score": prs["raw_score"],
            "variants_used": prs["variants_used"],
            "variants_total": prs["variants_total"],
            "overlap_fraction": prs["overlap_fraction"],
            "percentile": pct_info["percentile"],
            "risk_category": pct_info["risk_category"],
            "z_score": pct_info["z_score"],
            "method": pct_info["method"],
            "reference_population": pct_info.get("reference_population"),
        })

    return {
        "results": all_results,
        "scores_calculated": len(all_results),
        "disclaimer": _engine.DISCLAIMER,
    }
