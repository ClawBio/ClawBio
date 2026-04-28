#!/usr/bin/env python3
"""
nutrigx_advisor.py — NutriGx Advisor: Personalised Nutrition from Genetic Data
ClawBio Skill v0.2.0

Usage:
    python nutrigx_advisor.py --input genome.csv --output results/
    python nutrigx_advisor.py --input variants.vcf --output results/ --format vcf
"""

import argparse
import json
import sys
from pathlib import Path

# Add project root to path for shared imports
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.parsers import parse_genetic_file, genotypes_to_simple
from clawbio.common.report import write_result_json
from clawbio.common.checksums import sha256_hex

from extract_genotypes import extract_snp_genotypes
from score_variants import compute_nutrient_risk_scores
from generate_report import generate_report
from repro_bundle import create_reproducibility_bundle


def main():
    parser = argparse.ArgumentParser(
        description="NutriGx Advisor — personalised nutrigenomics report from genetic data"
    )
    parser.add_argument(
        "--input", required=True,
        help="Path to genetic data file (23andMe .txt/.csv, AncestryDNA .csv, MyHeritage .csv, or .vcf)"
    )
    parser.add_argument(
        "--output", default="nutrigx_results",
        help="Output directory (created if absent)"
    )
    parser.add_argument(
        "--format", choices=["auto", "23andme", "ancestry", "myheritage", "vcf"], default="auto",
        help="Input file format: 23andme, ancestry, myheritage, vcf (default: auto-detect)"
    )
    parser.add_argument(
        "--panel", default=None,
        help="Path to custom SNP panel JSON (default: data/snp_panel.json)"
    )
    parser.add_argument(
        "--no-figures", action="store_true",
        help="Skip figure generation (useful in headless environments)"
    )
    args = parser.parse_args()

    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)

    # Resolve SNP panel
    panel_path = (
        Path(args.panel) if args.panel
        else Path(__file__).parent / "data" / "snp_panel.json"
    )
    if not panel_path.exists():
        print(f"[ERROR] SNP panel not found at {panel_path}", file=sys.stderr)
        sys.exit(1)

    with open(panel_path) as f:
        snp_panel = json.load(f)

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"[ERROR] Input file not found: {input_path}", file=sys.stderr)
        sys.exit(1)

    # Check for empty file before parsing
    if input_path.stat().st_size == 0:
        print(
            f"[ERROR] Input file is empty: {input_path}. "
            f"Provide a file containing genotype data.",
            file=sys.stderr,
        )
        sys.exit(1)

    print(f"[NutriGx] Parsing input: {input_path}")
    try:
        records = parse_genetic_file(str(input_path), fmt=args.format)
    except ValueError as e:
        print(f"[ERROR] {e}", file=sys.stderr)
        sys.exit(1)
    genotype_table = genotypes_to_simple(records)
    print(f"[NutriGx] Loaded {len(genotype_table):,} variants")

    print("[NutriGx] Extracting SNP genotypes from panel ...")
    snp_calls = extract_snp_genotypes(genotype_table, snp_panel)

    present = sum(1 for v in snp_calls.values() if v["status"] == "found")
    mismatched = sum(1 for v in snp_calls.values() if v["status"] == "allele_mismatch")
    print(f"[NutriGx] Panel coverage: {present}/{len(snp_panel)} SNPs found")
    if mismatched > 0:
        print(f"[NutriGx] Allele mismatches: {mismatched} SNP(s) had unrecognised alleles")

    if present == 0:
        print(
            f"[ERROR] No panel SNPs found in input file. "
            f"0/{len(snp_panel)} SNPs matched. Cannot generate a report. "
            f"Check that the input file contains genotype data for the expected rsIDs.",
            file=sys.stderr,
        )
        sys.exit(1)

    min_coverage = len(snp_panel) * 0.25
    if present < min_coverage:
        print(
            f"[WARNING] *** LOW PANEL COVERAGE: only {present}/{len(snp_panel)} "
            f"SNPs found ({present/len(snp_panel)*100:.0f}%). "
            f"Report will have limited reliability. ***"
        )

    print("[NutriGx] Computing nutrient risk scores ...")
    risk_scores = compute_nutrient_risk_scores(snp_calls, snp_panel)

    print("[NutriGx] Generating report ...")
    report_path = generate_report(
        snp_calls=snp_calls,
        risk_scores=risk_scores,
        snp_panel=snp_panel,
        output_dir=str(output_dir),
        figures=not args.no_figures,
        input_file=str(input_path)
    )

    print("[NutriGx] Creating reproducibility bundle ...")
    create_reproducibility_bundle(
        input_file=str(input_path),
        output_dir=str(output_dir),
        panel_path=str(panel_path),
        args=vars(args)
    )

    # Write structured result.json for programmatic consumption
    print("[NutriGx] Writing result.json ...")
    elevated_domains = [d for d, v in risk_scores.items() if v["category"] == "Elevated"]
    moderate_domains = [d for d, v in risk_scores.items() if v["category"] == "Moderate"]
    total_tested = sum(v["tested_snps"] for v in risk_scores.values())
    total_missing = sum(v["missing_snps"] for v in risk_scores.values())

    write_result_json(
        output_dir=output_dir,
        skill="nutrigx",
        version="0.2.0",
        summary={
            "total_variants_loaded": len(genotype_table),
            "panel_snps_tested": total_tested,
            "panel_snps_missing": total_missing,
            "panel_size": len(snp_panel),
            "elevated_domains": elevated_domains,
            "moderate_domains": moderate_domains,
            "domains_assessed": len(risk_scores),
        },
        data={
            "risk_scores": risk_scores,
            "snp_calls": snp_calls,
        },
        input_checksum=sha256_hex(str(input_path)),
    )

    print(f"\n[NutriGx] Done. Report: {report_path}")
    print(f"[NutriGx] Results in: {output_dir}/")


if __name__ == "__main__":
    main()
