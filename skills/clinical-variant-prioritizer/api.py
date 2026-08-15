"""Importable API for the clinical-variant-prioritizer skill.

Screens a genotype dict against a curated clinical-gene panel and prioritises
carried variants by ClinVar significance, gnomAD frequency, inheritance and
zygosity, following the pathogenicity-screening method of Corpas et al. 2021
(Whole Genome Interpretation for a Family of Five, Front Genet 12:535123).

Usage::

    from api import run
    result = run({"rs28941785": "CT", "rs1800562": "GG"})
"""
from __future__ import annotations

import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from prioritize import load_panel, screen, build_headline  # noqa: E402

VERSION = "0.1.0"
METHOD = (
    "Corpas et al. 2021 (Front Genet 12:535123) pathogenicity screen: "
    "OMIM-morbid + ACMG-SF + Hereditary-Cancer panels -> ClinVar significance "
    "-> gnomAD frequency -> inheritance/zygosity prioritisation"
)
_DEFAULT_PANEL = _SKILL_DIR / "data" / "clinical_panel.json"


def run(genotypes: dict[str, str], options: dict | None = None) -> dict:
    """Run the clinical-variant prioritisation screen on a genotype dict.

    Args:
        genotypes: {rsid_or_variant_id: genotype_str} (e.g. {"rs28941785": "CT"}).
        options: Optional dict. Recognised keys:
            - 'panel_path': custom clinical panel JSON (default: built-in panel).

    Returns:
        dict with keys: skill, version, method, summary, findings (ranked),
        headline, disclaimer.
    """
    options = options or {}
    panel_path = Path(options.get("panel_path") or _DEFAULT_PANEL)
    if not panel_path.exists():
        raise FileNotFoundError(f"Clinical panel not found at {panel_path}")

    panel = load_panel(panel_path)
    findings, summary = screen(genotypes, panel)

    return {
        "skill": "clinical-variant-prioritizer",
        "version": VERSION,
        "method": METHOD,
        "summary": summary,
        "findings": findings,
        "headline": build_headline(summary),
        "disclaimer": (
            "Research and educational use only. Not a clinical diagnosis. "
            "Array-based screening covers only catalogued loci and misses most "
            "rare variants; confirm any finding with an accredited clinical assay."
        ),
    }
