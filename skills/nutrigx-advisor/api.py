"""Importable API for the nutrigx_advisor skill."""
from __future__ import annotations

import json
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# Ensure local skill modules are importable
_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from extract_genotypes import extract_snp_genotypes
from score_variants import compute_nutrient_risk_scores
from generate_report import generate_report


_DEFAULT_PANEL = _SKILL_DIR / "data" / "snp_panel.json"


def run(genotypes: dict[str, str], options: dict | None = None) -> dict:
    """Run nutrigx analysis on a genotype dict.

    Args:
        genotypes: {rsid: genotype_str} dict (e.g. from genotypes_to_simple()).
        options: Optional dict with keys:
            - 'panel_path': str | Path — custom SNP panel JSON (default: built-in panel)
            - 'output_dir': str | Path — write report + figures here (default: no file output)
            - 'no_figures': bool — skip matplotlib figure generation (default: False)
            - 'input_file': str — original input filename for report header (default: "api")

    Returns:
        Result dict compatible with PatientProfile.add_skill_result():
        {
            "skill": "nutrigx",
            "version": "0.2.0",
            "summary": {
                "total_variants": int,
                "panel_snps_tested": int,
                "panel_snps_missing": int,
                "panel_size": int,
                "elevated_domains": list[str],
                "moderate_domains": list[str],
                "domains_assessed": int,
            },
            "risk_scores": { domain: { score, category, ... } },
            "snp_calls": { rsid: { status, genotype, ... } },
            "report_path": str | None,
        }
    """
    options = options or {}

    # Load SNP panel
    panel_path = Path(options.get("panel_path") or _DEFAULT_PANEL)
    if not panel_path.exists():
        raise FileNotFoundError(f"SNP panel not found at {panel_path}")

    with open(panel_path) as f:
        snp_panel = json.load(f)

    # Extract genotypes from panel
    snp_calls = extract_snp_genotypes(genotypes, snp_panel)

    present = sum(1 for v in snp_calls.values() if v["status"] == "found")
    if present == 0:
        raise ValueError(
            f"No panel SNPs found in input genotypes. "
            f"0/{len(snp_panel)} SNPs matched. Cannot generate a report."
        )

    # Compute risk scores
    risk_scores = compute_nutrient_risk_scores(snp_calls, snp_panel)

    # Optionally write report to disk
    report_path = None
    output_dir = options.get("output_dir")
    if output_dir:
        output_dir = Path(output_dir)
        output_dir.mkdir(parents=True, exist_ok=True)
        input_file = options.get("input_file", "api")
        no_figures = options.get("no_figures", False)

        report_path = generate_report(
            snp_calls=snp_calls,
            risk_scores=risk_scores,
            snp_panel=snp_panel,
            output_dir=str(output_dir),
            figures=not no_figures,
            input_file=input_file,
        )

    # Build summary
    elevated_domains = [d for d, v in risk_scores.items() if v["category"] == "Elevated"]
    moderate_domains = [d for d, v in risk_scores.items() if v["category"] == "Moderate"]
    total_tested = sum(v["tested_snps"] for v in risk_scores.values())
    total_missing = sum(v["missing_snps"] for v in risk_scores.values())

    return {
        "skill": "nutrigx",
        "version": "0.2.0",
        "summary": {
            "total_variants": len(genotypes),
            "panel_snps_tested": total_tested,
            "panel_snps_missing": total_missing,
            "panel_size": len(snp_panel),
            "elevated_domains": elevated_domains,
            "moderate_domains": moderate_domains,
            "domains_assessed": len(risk_scores),
        },
        "risk_scores": risk_scores,
        "snp_calls": snp_calls,
        "report_path": report_path,
    }
