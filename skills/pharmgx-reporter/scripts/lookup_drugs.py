#!/usr/bin/env python3
"""Step 2 of pharmgx_reporter: look up CPIC drug recommendations for every
gene profile produced by `parse_and_profile.py`.

Ported from the "Lookup" stage of ClawBio 1.x `pharmgx_reporter.py`'s
`main()` (the `lookup_drugs()` call and its printed summary/alert block).

Usage:
    python lookup_drugs.py --profile /tmp/out/profile.json --out /tmp/out/drug_recommendations.json
"""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent))
import pgx_common as pgx


def main() -> int:
    parser = argparse.ArgumentParser(description="Look up CPIC drug recommendations from gene profiles")
    parser.add_argument("--profile", required=True, help="Path to profile.json from parse_and_profile.py")
    parser.add_argument("--out", required=True, help="Output path for drug_recommendations.json")
    args = parser.parse_args()

    profile_data = json.loads(Path(args.profile).read_text())
    profiles = profile_data["profiles"]

    drug_results = pgx.lookup_drugs(profiles)
    n_std = len(drug_results["standard"])
    n_cau = len(drug_results["caution"])
    n_avo = len(drug_results["avoid"])
    n_ind = len(drug_results.get("indeterminate", []))
    total_assessed = n_std + n_cau + n_avo + n_ind

    print(f"Drug Recommendations ({total_assessed} drugs):")
    print(f"  Standard:          {n_std}")
    print(f"  Caution:           {n_cau}")
    print(f"  Avoid:             {n_avo}")
    if n_ind > 0:
        print(f"  Insufficient data: {n_ind}")

    if n_avo > 0:
        print("\nALERT - Drugs to AVOID:")
        for d in drug_results["avoid"]:
            print(f"  * {d['drug']} ({d['brand']})")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    payload = {
        "drug_recommendations": drug_results,
        "counts": {
            "standard": n_std,
            "caution": n_cau,
            "avoid": n_avo,
            "indeterminate": n_ind,
            "total": total_assessed,
        },
    }
    out_path.write_text(json.dumps(payload, indent=2))
    print(f"\nDrug recommendations written: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
