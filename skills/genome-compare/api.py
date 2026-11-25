"""Importable API for the genome-compare skill.

Allows other skills and the orchestrator to invoke genome comparison
programmatically without shelling out to the CLI.

Usage::

    from skills.genome_compare.api import run

    result = run(genotypes={"rs1234": "AA", "rs5678": "AG"})
"""

from __future__ import annotations

import sys
from pathlib import Path

# Ensure project root is on sys.path for shared imports
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))


try:
    from .genome_compare import (
        REFERENCE_FILE,
        AIMS_PANEL_FILE,
        compute_ibs,
        estimate_ancestry,
        load_aims_panel,
        _parse_genotype_file,
    )
except ImportError:
    from genome_compare import (  # type: ignore[no-redef]
        REFERENCE_FILE,
        AIMS_PANEL_FILE,
        compute_ibs,
        estimate_ancestry,
        load_aims_panel,
        _parse_genotype_file,
    )


def run(genotypes: dict[str, str], options: dict | None = None) -> dict:
    """Run genome comparison on a genotype dict.

    Args:
        genotypes: Mapping of rsid -> genotype string (e.g. ``{"rs1234": "AA"}``).
        options: Optional settings:
            - ``reference_path``: Path to reference genome file.
            - ``aims_path``: Path to AIMs panel JSON.

    Returns:
        dict with keys: ``ibs_score``, ``n_overlap``, ``n_concordant``,
        ``per_chrom``, ``ancestry``.
    """
    options = options or {}

    reference_path = Path(options.get("reference_path", str(REFERENCE_FILE)))
    aims_path = Path(options.get("aims_path", str(AIMS_PANEL_FILE)))

    # Parse reference genome
    ref_geno, _ = _parse_genotype_file(reference_path)

    # Build a simple positions dict from user genotypes (no position data
    # available when called with a plain dict, so per-chrom breakdown will
    # fall back to "?" chromosome).
    user_pos: dict[str, dict] = {}

    # Compute IBS
    ibs_score, n_overlap, n_concordant, per_chrom = compute_ibs(
        genotypes, ref_geno, user_pos or None
    )

    # Estimate ancestry
    aims_panel, populations = load_aims_panel(aims_path)
    user_ancestry = estimate_ancestry(genotypes, aims_panel, populations)

    return {
        "ibs_score": ibs_score,
        "n_overlap": n_overlap,
        "n_concordant": n_concordant,
        "per_chrom": per_chrom,
        "ancestry": user_ancestry,
    }
