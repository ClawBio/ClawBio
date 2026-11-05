#!/usr/bin/env python3
"""
extract_genotypes.py — DAG step 2/4: SNP panel lookup with forward-strand
normalisation.

Ported verbatim (core extract_snp_genotypes logic unchanged) from ClawBio
1.x skills/nutrigx/extract_genotypes.py. Reads the genotype table JSON
emitted by parse_input.py plus the curated SNP panel JSON, and writes the
per-SNP call table (status/genotype/risk_count) consumed by score_variants.py
and generate_report.py.

The "no panel SNPs found" fatal check that lived in 1.x's nutrigx.py main()
is reproduced here (rather than deferred to score_variants.py), because
panel coverage can only be evaluated once genotypes are extracted — this is
the natural fail-fast point in the DAG, matching where 1.x checked it
(immediately after extraction, before scoring).

Usage:
    python3 extract_genotypes.py --genotype-table genotype_table.json \\
        --panel snp_panel.json --out snp_calls.json
"""

import argparse
import json
import sys
from pathlib import Path

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# SNPs where both alleles are complementary (ambiguous strand)
AMBIGUOUS_PAIRS = {frozenset(["A", "T"]), frozenset(["C", "G"])}


def flip_genotype(genotype: str) -> str:
    """Return the complement strand genotype."""
    return "".join(COMPLEMENT.get(b, b) for b in genotype)


def is_ambiguous(ref: str, alt: str) -> bool:
    return frozenset([ref, alt]) in AMBIGUOUS_PAIRS


def extract_snp_genotypes(genotype_table: dict, snp_panel: list) -> dict:
    """
    For each SNP in the panel, look up the genotype in genotype_table.

    Returns dict keyed by rsid:
    {
      "rsid": "rs1801133",
      "status": "found" | "not_tested",
      "genotype": "CT",           # raw as reported
      "normalised": "CT",         # forward-strand normalised
      "risk_allele": "T",
      "risk_count": 1             # 0, 1, or 2 copies of risk allele
    }
    """
    results = {}

    for snp in snp_panel:
        rsid = snp["rsid"]
        risk_allele = snp["risk_allele"]
        ref_allele = snp.get("ref_allele", "")

        if rsid not in genotype_table:
            results[rsid] = {
                "rsid": rsid,
                "gene": snp["gene"],
                "status": "not_tested",
                "genotype": None,
                "normalised": None,
                "risk_allele": risk_allele,
                "risk_count": None,
                "nutrient_domain": snp["nutrient_domain"],
            }
            continue

        raw_geno = genotype_table[rsid]
        if not raw_geno or len(raw_geno) < 2:
            results[rsid] = {
                "rsid": rsid,
                "gene": snp["gene"],
                "status": "no_call",
                "genotype": raw_geno,
                "normalised": None,
                "risk_allele": risk_allele,
                "risk_count": None,
                "nutrient_domain": snp["nutrient_domain"],
            }
            continue

        # Try direct match first
        norm = raw_geno
        allele_matched = risk_allele in raw_geno
        if not allele_matched:
            # Try strand flip
            flipped = flip_genotype(raw_geno)
            if risk_allele in flipped:
                norm = flipped
                allele_matched = True

        if allele_matched:
            risk_count = norm.count(risk_allele)
            results[rsid] = {
                "rsid": rsid,
                "gene": snp["gene"],
                "status": "found",
                "genotype": raw_geno,
                "normalised": norm,
                "risk_allele": risk_allele,
                "risk_count": risk_count,
                "nutrient_domain": snp["nutrient_domain"],
            }
        else:
            # Genotype does not contain the risk allele on either strand.
            # Check if this is homozygous reference (0 copies of risk allele).
            # A true allele_mismatch is when the genotype contains alleles
            # that are neither ref nor risk (e.g. tri-allelic or data error).
            geno_alleles = set(norm)
            known_alleles = {risk_allele, ref_allele} if ref_allele else {risk_allele}
            flipped_known = {COMPLEMENT.get(a, a) for a in known_alleles}

            if ref_allele and geno_alleles <= {ref_allele}:
                # Homozygous reference: 0 copies of risk allele
                results[rsid] = {
                    "rsid": rsid,
                    "gene": snp["gene"],
                    "status": "found",
                    "genotype": raw_geno,
                    "normalised": norm,
                    "risk_allele": risk_allele,
                    "risk_count": 0,
                    "nutrient_domain": snp["nutrient_domain"],
                }
            elif ref_allele and geno_alleles <= {COMPLEMENT.get(ref_allele, ref_allele)}:
                # Homozygous reference on complement strand
                results[rsid] = {
                    "rsid": rsid,
                    "gene": snp["gene"],
                    "status": "found",
                    "genotype": raw_geno,
                    "normalised": flip_genotype(raw_geno),
                    "risk_allele": risk_allele,
                    "risk_count": 0,
                    "nutrient_domain": snp["nutrient_domain"],
                }
            else:
                # True allele mismatch: genotype contains unknown alleles
                print(
                    f"[WARNING] {rsid} ({snp['gene']}): genotype '{raw_geno}' "
                    f"does not contain risk allele '{risk_allele}' (even after strand flip). "
                    f"Setting allele_mismatch.",
                    file=sys.stderr,
                )
                results[rsid] = {
                    "rsid": rsid,
                    "gene": snp["gene"],
                    "status": "allele_mismatch",
                    "genotype": raw_geno,
                    "normalised": norm,
                    "risk_allele": risk_allele,
                    "risk_count": None,
                    "nutrient_domain": snp["nutrient_domain"],
                    "warning": (
                        f"Genotype '{raw_geno}' does not contain risk allele "
                        f"'{risk_allele}' on either strand"
                    ),
                }

    return results


def main():
    ap = argparse.ArgumentParser(description="Extract curated SNP panel genotypes from a genotype table")
    ap.add_argument("--genotype-table", required=True, help="Path to genotype_table.json from parse_input.py")
    ap.add_argument("--panel", required=True, help="Path to snp_panel.json")
    ap.add_argument("--out", required=True, help="Path to write snp_calls.json")
    args = ap.parse_args()

    with open(args.genotype_table) as f:
        genotype_table = json.load(f)
    with open(args.panel) as f:
        snp_panel = json.load(f)

    print("[extract_genotypes] Extracting SNP genotypes from panel ...")
    snp_calls = extract_snp_genotypes(genotype_table, snp_panel)

    present = sum(1 for v in snp_calls.values() if v["status"] == "found")
    mismatched = sum(1 for v in snp_calls.values() if v["status"] == "allele_mismatch")
    print(f"[extract_genotypes] Panel coverage: {present}/{len(snp_panel)} SNPs found")
    if mismatched > 0:
        print(f"[extract_genotypes] Allele mismatches: {mismatched} SNP(s) had unrecognised alleles")

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

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(snp_calls, indent=2, sort_keys=True))
    print(f"[extract_genotypes] Wrote SNP calls: {out_path}")


if __name__ == "__main__":
    main()
