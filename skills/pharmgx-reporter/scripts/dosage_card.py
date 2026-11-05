#!/usr/bin/env python3
"""Standalone single-drug lookup: resolve a drug name (brand or generic,
fuzzy-matched) against a patient's PGx gene profile and print a dosage card.

Ported from ClawBio 1.x `pharmgx_reporter.py`'s `--drug` mode (the
`resolve_drug_name()` / `lookup_single_drug()` / `format_dosage_card()`
path in `main()`, which exits before the multi-drug report is built).

Why this is a standalone script and not a pipeline.yaml step: it is a
genuinely different *mode* of the underlying tool, not a stage of the
full-report DAG — running it always produces exactly one output (a card
printed to stdout), never a `report.md`/`result.json` pair. `pipeline.yaml`
v1 forbids conditionals, so a DAG cannot branch between "full report" and
"single drug" based on whether `--drug` was passed; encoding that choice in
the yaml would require exactly the `if:`/`when:` machinery CLAUDE.md rules
out. Per CLAUDE.md: "If a step needs branching logic, write a Python script
and make that script the step body" — this script *is* that alternate
body, invoked directly rather than through `clawbio open`/`done`. This is
also the shape `skills/drug_photo/SKILL.md` describes for its own chaining
into this skill: an agent resolves a drug name (via vision) and genotype
lookup separately, outside any single skill's DAG.

Usage:
    python dosage_card.py --input demo_patient.txt --drug Plavix
    python dosage_card.py --input demo_patient.txt --drug codeine --dose 30mg
"""
from __future__ import annotations

import argparse
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pgx_common as pgx


def main() -> int:
    parser = argparse.ArgumentParser(description="Single-drug PGx dosage card lookup")
    parser.add_argument("--input", required=True, help="Path to 23andMe/AncestryDNA raw data file")
    parser.add_argument("--drug", required=True, help="Drug name (brand or generic)")
    parser.add_argument("--dose", default=None, help="Visible dose from packaging (e.g. '50mg')")
    args = parser.parse_args()

    input_path = Path(args.input)
    if not input_path.exists():
        print(f"Error: input file not found: {input_path}", file=sys.stderr)
        return 1

    fmt, total_snps, pgx_snps, ref_genome = pgx.parse_file(input_path)
    if len(pgx_snps) == 0:
        print("ERROR: No pharmacogenomic SNPs found in this file.", file=sys.stderr)
        return 1

    profiles = {}
    for gene in pgx.GENE_DEFS:
        diplotype = pgx.call_diplotype(gene, pgx_snps)
        phenotype = pgx.call_phenotype(gene, diplotype)
        profiles[gene] = {"diplotype": diplotype, "phenotype": phenotype}

    if ref_genome == "unknown_build":
        for gene in profiles:
            profiles[gene] = {
                "diplotype": "Indeterminate (coordinates match no known build)",
                "phenotype": "Indeterminate (unrecognised reference genome)",
            }

    resolved = pgx.resolve_drug_name(args.drug)
    if not resolved:
        print(f"Drug not found: '{args.drug}'. Available drugs: {len(pgx.GUIDELINES)}", file=sys.stderr)
        return 1

    result = pgx.lookup_single_drug(resolved, profiles)
    print(pgx.format_dosage_card(result, visible_dose=args.dose))
    return 0


if __name__ == "__main__":
    sys.exit(main())
