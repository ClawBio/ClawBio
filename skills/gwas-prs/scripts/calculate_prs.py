#!/usr/bin/env python3
"""Stage 3 of the gwas_prs DAG: compute PRS + percentile for each resolved score.

Port of ClawBio 1.x `gwas-prs` (skills/gwas-prs/gwas_prs.py). Core scoring
math (parse_scoring_file, compute_allele_dosage, calculate_prs,
estimate_percentile) is unchanged from the 1.x source. The I/O boundary
changed: this stage consumes the JSON manifest/genotype emitted by
scripts/resolve_scores.py and scripts/parse_genotype.py rather than doing
file resolution and genotype parsing inline. It also carries the manifest's
curated-demo-panel provenance (curated_demo_panel, curated_panel_id,
legacy_pgs_id) through to prs_results.json so the report can label it.

Standard additive PRS model: PRS = sum(dosage_i * effect_weight_i).
"""

from __future__ import annotations

import argparse
import gzip
import json
import math
import sys
from pathlib import Path

MIN_OVERLAP = 0.5  # skip a score if fewer than this fraction of its variants match

# Percentile thresholds -> risk category. NOTE: this matches the 1.x
# *code* (gwas_prs.py RISK_CATEGORIES), which differs from the 1.x
# SKILL.md prose (20/80/95 there vs 25/75/95 here) — the port is faithful
# to the executed logic, not the marketing description.
RISK_CATEGORIES = [
    (25, "Low"),
    (75, "Average"),
    (95, "Elevated"),
    (100, "High"),
]


def parse_scoring_file(filepath: str | Path) -> list[dict]:
    """Parse a PGS Catalog harmonized scoring file (TSV, may be gzipped)."""
    filepath = Path(filepath)
    if str(filepath).endswith(".gz"):
        with gzip.open(filepath, "rt", errors="replace") as fh:
            raw_lines = fh.readlines()
    else:
        raw_lines = filepath.read_text(errors="replace").splitlines(keepends=True)

    header_line = None
    data_lines = []
    for line in raw_lines:
        stripped = line.strip()
        if stripped.startswith("#"):
            continue
        if header_line is None:
            header_line = stripped
            continue
        if stripped:
            data_lines.append(stripped)

    if header_line is None:
        return []

    headers = header_line.split("\t")
    col_map: dict[str, int] = {h.strip().lower(): i for i, h in enumerate(headers)}

    rsid_col = col_map.get("hm_rsid", col_map.get("rsid"))
    effect_allele_col = col_map.get("effect_allele")
    weight_col = col_map.get("effect_weight")
    other_allele_col = col_map.get("other_allele", col_map.get("hm_inferotherallele"))
    freq_col = col_map.get("allelefrequency_effect")

    if rsid_col is None or effect_allele_col is None or weight_col is None:
        print(f"  WARNING: scoring file missing required columns. Found: {list(col_map.keys())}")
        return []

    variants = []
    for line in data_lines:
        cols = line.split("\t")
        if len(cols) <= max(rsid_col, effect_allele_col, weight_col):
            continue
        rsid = cols[rsid_col].strip()
        if not rsid or not rsid.startswith("rs"):
            continue
        effect_allele = cols[effect_allele_col].strip().upper()
        try:
            weight = float(cols[weight_col].strip())
        except (ValueError, IndexError):
            continue

        variant: dict = {"rsid": rsid, "effect_allele": effect_allele, "effect_weight": weight}
        if other_allele_col is not None and other_allele_col < len(cols):
            oa = cols[other_allele_col].strip().upper()
            if oa:
                variant["other_allele"] = oa
        if freq_col is not None and freq_col < len(cols):
            try:
                variant["allele_freq"] = float(cols[freq_col].strip())
            except (ValueError, IndexError):
                pass
        variants.append(variant)

    return variants


def compute_allele_dosage(genotype: str, effect_allele: str) -> int:
    """Count of effect allele copies (0, 1, or 2; 0 or 1 for hemizygous)."""
    effect_allele = effect_allele.upper()
    genotype = genotype.upper()
    if len(genotype) == 1:
        return 1 if genotype == effect_allele else 0
    return sum(1 for allele in genotype if allele == effect_allele)


def calculate_prs(genotype_dict: dict[str, str], scoring_variants: list[dict]) -> dict:
    """PRS = sum(weight * dosage). Missing genotypes are excluded from the sum."""
    raw_score = 0.0
    variants_used = 0
    variants_missing = 0
    per_variant: list[dict] = []

    for sv in scoring_variants:
        rsid = sv["rsid"]
        effect_allele = sv["effect_allele"]
        weight = sv["effect_weight"]

        if rsid in genotype_dict:
            genotype = genotype_dict[rsid]
            dosage = compute_allele_dosage(genotype, effect_allele)
            contribution = weight * dosage
            raw_score += contribution
            variants_used += 1
            per_variant.append({
                "rsid": rsid, "effect_allele": effect_allele, "genotype": genotype,
                "dosage": dosage, "weight": weight, "contribution": contribution,
                "status": "scored",
            })
        else:
            variants_missing += 1
            per_variant.append({
                "rsid": rsid, "effect_allele": effect_allele, "genotype": None,
                "dosage": None, "weight": weight, "contribution": 0.0, "status": "missing",
            })

    variants_total = len(scoring_variants)
    overlap = variants_used / variants_total if variants_total > 0 else 0.0

    return {
        "raw_score": raw_score,
        "variants_used": variants_used,
        "variants_total": variants_total,
        "variants_missing": variants_missing,
        "overlap_fraction": overlap,
        "per_variant": per_variant,
    }


def _percentile_from_z(z: float) -> float:
    """Standard normal CDF via math.erf (exact closed form, no scipy needed)."""
    return 50.0 * (1.0 + math.erf(z / math.sqrt(2.0)))


def _assign_risk_category(percentile: float) -> str:
    for threshold, label in RISK_CATEGORIES:
        if percentile <= threshold:
            return label
    return "High"


def estimate_percentile(raw_score: float, reference_distribution: dict | None,
                         scoring_variants: list[dict]) -> dict:
    """Tier 1: curated mean/SD reference. Tier 2: allele-frequency HWE estimate.
    Tier 3: unavailable.
    """
    result: dict = {
        "percentile": None, "method": "unavailable", "risk_category": None,
        "z_score": None, "reference_population": None,
    }

    if reference_distribution:
        mean = reference_distribution["mean"]
        sd = reference_distribution["sd"]
        if sd > 0:
            z = (raw_score - mean) / sd
            pct = _percentile_from_z(z)
            result.update(
                percentile=round(pct, 1), method="curated_reference",
                risk_category=_assign_risk_category(pct), z_score=round(z, 3),
                reference_population=reference_distribution.get("population", "EUR"),
            )
            return result

    freqs_available = [
        sv for sv in scoring_variants
        if "allele_freq" in sv and 0 < sv["allele_freq"] < 1
    ]
    if len(freqs_available) >= 3:
        expected = 0.0
        variance = 0.0
        for sv in freqs_available:
            freq = sv["allele_freq"]
            w = sv["effect_weight"]
            expected += 2.0 * freq * w
            variance += 2.0 * freq * (1.0 - freq) * (w ** 2)
        if variance > 0:
            sd = math.sqrt(variance)
            z = (raw_score - expected) / sd
            pct = _percentile_from_z(z)
            result.update(
                percentile=round(pct, 1), method="allele_frequency_estimate",
                risk_category=_assign_risk_category(pct), z_score=round(z, 3),
                reference_population="estimated (allele freq)",
            )
            return result

    return result


def main() -> None:
    parser = argparse.ArgumentParser(description="Compute PRS + percentile for each resolved score")
    parser.add_argument("--genotype-json", required=True, help="Path to genotype.json (from parse_genotype.py)")
    parser.add_argument("--scoring-manifest", required=True, help="Path to scoring_manifest.json (from resolve_scores.py)")
    parser.add_argument("--out-results", required=True, help="Output path for prs_results.json")
    parser.add_argument("--out-variants", required=True, help="Output path for prs_variants.csv")
    parser.add_argument("--min-overlap", type=float, default=MIN_OVERLAP, help="Minimum SNP overlap fraction (default: 0.5)")
    args = parser.parse_args()

    genotype_data = json.loads(Path(args.genotype_json).read_text())
    genotype_dict = genotype_data["genotypes"]
    manifest = json.loads(Path(args.scoring_manifest).read_text())

    all_results: list[dict] = []
    for entry in manifest:
        pgs_id = entry["pgs_id"]
        trait = entry.get("trait", "Unknown")

        if entry["source"] == "error" or not entry.get("filepath"):
            print(f"Skipping {pgs_id}: unresolved ({entry.get('error', 'no scoring file')})")
            continue

        print(f"Scoring {pgs_id} ({trait})...")
        scoring_variants = parse_scoring_file(entry["filepath"])
        if not scoring_variants:
            print("  WARNING: no valid variants in scoring file, skipping")
            continue
        print(f"  Variants in score: {len(scoring_variants)}")

        prs = calculate_prs(genotype_dict, scoring_variants)
        print(f"  Variants matched: {prs['variants_used']}/{prs['variants_total']} "
              f"({prs['overlap_fraction'] * 100:.1f}%)")

        if prs["overlap_fraction"] < args.min_overlap:
            print(f"  Skipping: overlap {prs['overlap_fraction'] * 100:.1f}% "
                  f"below --min-overlap {args.min_overlap * 100:.0f}%")
            continue

        print(f"  Raw PRS: {prs['raw_score']:.6f}")

        pct_info = estimate_percentile(
            prs["raw_score"], entry.get("reference_distribution"), scoring_variants
        )
        if pct_info["percentile"] is not None:
            print(f"  Percentile: {pct_info['percentile']:.1f}% ({pct_info['risk_category']}) "
                  f"[{pct_info['method']}]")
        else:
            print("  Percentile: N/A (no reference distribution available)")

        all_results.append({
            "pgs_id": pgs_id,
            "trait": trait,
            "metadata": entry.get("metadata", {}),
            "curated_demo_panel": bool(entry.get("curated_demo_panel")),
            "curated_panel_id": entry.get("curated_panel_id"),
            "legacy_pgs_id": entry.get("legacy_pgs_id"),
            "prs": prs,
            "percentile_info": pct_info,
        })
        print()

    if not all_results:
        print("No scores produced results. Check genotype file and overlap.")
        sys.exit(1)

    out_results = Path(args.out_results)
    out_results.parent.mkdir(parents=True, exist_ok=True)
    out_results.write_text(json.dumps({
        "input_info": {
            "filepath": genotype_data["filepath"],
            "format": genotype_data["format"],
            "total_snps": genotype_data["total_snps"],
        },
        "results": all_results,
    }, indent=2))
    print(f"Results written to {out_results}")

    csv_lines = ["pgs_id,rsid,effect_allele,genotype,dosage,weight,contribution,status"]
    for r in all_results:
        for v in r["prs"]["per_variant"]:
            gt = v["genotype"] if v["genotype"] else ""
            dosage = str(v["dosage"]) if v["dosage"] is not None else ""
            csv_lines.append(
                f"{r['pgs_id']},{v['rsid']},{v['effect_allele']},{gt},{dosage},"
                f"{v['weight']:.6f},{v['contribution']:.6f},{v['status']}"
            )
    out_variants = Path(args.out_variants)
    out_variants.parent.mkdir(parents=True, exist_ok=True)
    out_variants.write_text("\n".join(csv_lines) + "\n")
    print(f"Variant details written to {out_variants}")


if __name__ == "__main__":
    main()
