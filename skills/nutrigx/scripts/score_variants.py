#!/usr/bin/env python3
"""
score_variants.py — DAG step 3/4: composite per-nutrient-domain risk scores.

Ported verbatim (core compute_nutrient_risk_scores logic unchanged) from
ClawBio 1.x skills/nutrigx/score_variants.py. Reads the snp_calls.json
emitted by extract_genotypes.py plus the SNP panel JSON, and writes a
risk_scores.json (per-domain 0-10 composite score, category, and
contributing SNPs) consumed by generate_report.py.

Usage:
    python3 score_variants.py --snp-calls snp_calls.json --panel snp_panel.json \\
        --out risk_scores.json
"""

import argparse
import json
from collections import defaultdict
from pathlib import Path


def snp_raw_score(risk_count):
    """
    Convert risk allele count to a 0-1 dosage score.
    0 copies -> 0.0 (homozygous reference)
    1 copy   -> 0.5 (heterozygous)
    2 copies -> 1.0 (homozygous risk)
    """
    if risk_count is None:
        return None
    dosage_map = {0: 0.0, 1: 0.5, 2: 1.0}
    if risk_count not in dosage_map:
        raise ValueError(
            f"Unexpected risk_count={risk_count!r}. "
            f"Expected 0, 1, 2, or None."
        )
    return dosage_map[risk_count]


def compute_nutrient_risk_scores(snp_calls: dict, snp_panel: list) -> dict:
    """
    Returns a dict of nutrient_domain -> {
        'score': float (0-10),
        'category': str ('Low' | 'Moderate' | 'Elevated'),
        'contributing_snps': list,
        'tested_snps': int,
        'missing_snps': int,
    }
    """
    # Index panel by rsid
    panel_index = {s["rsid"]: s for s in snp_panel}

    # Group by nutrient domain
    domain_snps = defaultdict(list)
    for snp in snp_panel:
        domain_snps[snp["nutrient_domain"]].append(snp["rsid"])

    results = {}

    for domain, rsids in domain_snps.items():
        weighted_sum = 0.0
        max_possible = 0.0
        contributing = []
        tested = 0
        missing = 0

        for rsid in rsids:
            call = snp_calls.get(rsid)
            panel_entry = panel_index[rsid]
            weight = panel_entry.get("weight", 0.5)

            if call is None or call["status"] not in ("found",):
                missing += 1
                continue

            raw = snp_raw_score(call["risk_count"])
            if raw is None:
                missing += 1
                continue

            tested += 1
            weighted_sum += raw * weight
            max_possible += weight
            contributing.append({
                "rsid": rsid,
                "gene": panel_entry["gene"],
                "genotype": call["normalised"],
                "risk_count": call["risk_count"],
                "raw_score": raw,
                "weight": weight,
                "effect_direction": panel_entry.get("effect_direction", ""),
            })

        # Normalise to 0-10; handle all-missing
        if max_possible > 0:
            score = round((weighted_sum / max_possible) * 10, 2)
        elif tested == 0 and missing > 0:
            score = None  # no data for this domain
        else:
            score = 0.0

        if score is None:
            category = "Unknown"
        elif score < 3.5:
            category = "Low"
        elif score < 6.5:
            category = "Moderate"
        else:
            category = "Elevated"

        total_in_domain = tested + missing
        coverage_str = f"{tested}/{total_in_domain} SNPs tested"

        results[domain] = {
            "score": score,
            "category": category,
            "contributing_snps": contributing,
            "tested_snps": tested,
            "missing_snps": missing,
            "coverage": coverage_str,
        }

    return results


def main():
    ap = argparse.ArgumentParser(description="Compute per-nutrient-domain composite risk scores")
    ap.add_argument("--snp-calls", required=True, help="Path to snp_calls.json from extract_genotypes.py")
    ap.add_argument("--panel", required=True, help="Path to snp_panel.json")
    ap.add_argument("--out", required=True, help="Path to write risk_scores.json")
    args = ap.parse_args()

    with open(args.snp_calls) as f:
        snp_calls = json.load(f)
    with open(args.panel) as f:
        snp_panel = json.load(f)

    print("[score_variants] Computing nutrient risk scores ...")
    risk_scores = compute_nutrient_risk_scores(snp_calls, snp_panel)

    elevated = [d for d, v in risk_scores.items() if v["category"] == "Elevated"]
    moderate = [d for d, v in risk_scores.items() if v["category"] == "Moderate"]
    print(
        f"[score_variants] {len(risk_scores)} domain(s) scored: "
        f"{len(elevated)} elevated, {len(moderate)} moderate"
    )

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(risk_scores, indent=2, sort_keys=True))
    print(f"[score_variants] Wrote risk scores: {out_path}")


if __name__ == "__main__":
    main()
