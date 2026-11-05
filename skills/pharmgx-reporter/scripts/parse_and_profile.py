#!/usr/bin/env python3
"""Step 1 of pharmgx_reporter: parse a 23andMe/AncestryDNA file, extract the
PGx SNP panel, call star-allele diplotypes, and map to CPIC phenotypes.

Ported from the "Parse" + "Extract" + "Call" stages of ClawBio 1.x
`pharmgx_reporter.py`'s `main()` (format detection through gene-profile
construction, including the GRCh37/38 build-detection gate and the UGT1A1
TA-repeat-proxy special case). Splits cleanly here because 1.x's single
script becomes this DAG's first step; `lookup_drugs.py` and `report.py`
consume this step's `profile.json` rather than re-parsing the input file.

Usage:
    python parse_and_profile.py --input demo_patient.txt --out /tmp/out/profile.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pgx_common as pgx


def main() -> int:
    parser = argparse.ArgumentParser(description="Parse genetic data and call PGx gene profiles")
    parser.add_argument("--input", required=True, help="Path to 23andMe/AncestryDNA raw data file")
    parser.add_argument("--out", required=True, help="Output path for profile.json")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        return 1

    print(f"Parsing: {input_path}")
    fmt, total_snps, pgx_snps, ref_genome = pgx.parse_file(input_path)
    print(f"  Format: {fmt}")
    print(f"  Total SNPs: {total_snps}")
    print(f"  PGx SNPs found: {len(pgx_snps)}/{len(pgx.PGX_SNPS)}")
    if ref_genome:
        print(f"  Reference genome: {ref_genome}")

    if fmt == "unknown":
        print("WARNING: Could not detect input file format (expected 23andMe or "
              "AncestryDNA raw data). Results may be unreliable.", file=sys.stderr)

    if ref_genome == "GRCh37":
        print("NOTE: Input coordinates are GRCh37. Diplotype calls are matched "
              "by rsID, so calls are unaffected; positions in this report are "
              "GRCh37.", file=sys.stderr)
    elif ref_genome == "unknown_build":
        print("WARNING: Input coordinates could not be matched to GRCh37 or "
              "GRCh38 (unknown_build). The file may be corrupt or annotated "
              "against an unsupported assembly; results are withheld.",
              file=sys.stderr)

    if len(pgx_snps) == 0:
        print("ERROR: No pharmacogenomic SNPs found in this file.", file=sys.stderr)
        print("Cannot generate a report from zero data. Verify the input file "
              "is a valid 23andMe or AncestryDNA export.", file=sys.stderr)
        return 1

    # Profile genes
    profiles = {}
    for gene in pgx.GENE_DEFS:
        diplotype = pgx.call_diplotype(gene, pgx_snps)
        phenotype = pgx.call_phenotype(gene, diplotype)
        profiles[gene] = {"diplotype": diplotype, "phenotype": phenotype}

    # Genotypes are matched to the PGx panel by rsID, not coordinate, so a
    # GRCh37 file calls identically to the same sample on GRCh38. Only a file
    # whose coordinates match no known build is withheld.
    if ref_genome == "unknown_build":
        for gene in profiles:
            profiles[gene] = {
                "diplotype": "Indeterminate (coordinates match no known build)",
                "phenotype": "Indeterminate (unrecognised reference genome)",
            }

    # UGT1A1: without the TA-repeat marker SNP (rs8175347), *28 cannot be
    # assessed, so the call is unreliable even if the *6 SNP is present.
    if "UGT1A1" in profiles and ref_genome != "unknown_build":
        ugt_sv_rsid = "rs8175347"
        if ugt_sv_rsid not in pgx_snps:
            profiles["UGT1A1"] = {
                "diplotype": "Indeterminate (rs8175347 not tested)",
                "phenotype": "Indeterminate (rs8175347 not tested)",
            }

    not_tested = [g for g, p in profiles.items() if p["diplotype"] == "NOT_TESTED"]
    if not_tested:
        print(f"WARNING: {len(not_tested)} gene(s) not testable from this data: "
              f"{', '.join(not_tested)}", file=sys.stderr)

    print("Gene Profiles:")
    print(f"  {'Gene':<10} {'Diplotype':<20} {'Phenotype'}")
    print(f"  {'-' * 10} {'-' * 20} {'-' * 35}")
    for gene, p in profiles.items():
        print(f"  {gene:<10} {p['diplotype']:<20} {p['phenotype']}")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "input_path": str(input_path.resolve()),
        "input_checksum": pgx.sha256_file(input_path),
        "format": fmt,
        "total_snps": total_snps,
        "pgx_snps_found": len(pgx_snps),
        "pgx_snps_total": len(pgx.PGX_SNPS),
        "ref_genome": ref_genome,
        "pgx_snps": pgx_snps,
        "profiles": profiles,
    }
    out_path.write_text(json.dumps(payload, indent=2))
    print(f"\nProfile written: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
