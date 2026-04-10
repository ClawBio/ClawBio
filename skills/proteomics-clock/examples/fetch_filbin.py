#!/usr/bin/env python3
"""Download and preprocess Filbin et al. (2021) COVID-19 Olink data.

The Filbin dataset contains longitudinal Olink proteomics from COVID-19
patients at Day 0, 3, and 7. This script converts the raw Excel files
into per-timepoint CSVs ready for the proteomics clock.

Data source: https://doi.org/10.17632/nf853r8xsj.2 (Mendeley Data, free)

Usage:
    # 1. Download manually from Mendeley (3 files):
    #    - Olink_Proteomics.xlsx
    #    - Clinical_Metadata.xlsx
    #    - Supplemental-Table-3-Olink-Models-All.xlsx
    #
    # 2. Run this script:
    python examples/fetch_filbin.py --data-dir /path/to/downloaded/files --output /tmp/filbin_prepared
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

import pandas as pd


def load_and_convert(data_dir: Path, output_dir: Path, verbose: bool = True) -> dict:
    """Load Filbin Excel files, convert OlinkIDs to gene symbols, split by timepoint."""

    olink_path = data_dir / "Olink_Proteomics.xlsx"
    meta_path = data_dir / "Clinical_Metadata.xlsx"
    supp_path = data_dir / "Supplemental-Table-3-Olink-Models-All.xlsx"

    for p in [olink_path, meta_path, supp_path]:
        if not p.exists():
            print(f"ERROR: Missing {p.name}. Download from https://doi.org/10.17632/nf853r8xsj.2")
            sys.exit(1)

    if verbose:
        print("Loading Olink_Proteomics.xlsx ...")
    olink = pd.read_excel(olink_path, sheet_name="Olink Proteomics")

    if verbose:
        print("Loading Clinical_Metadata.xlsx ...")
    clinical = pd.read_excel(meta_path, sheet_name="Subject-level metadata")

    if verbose:
        print("Loading Supplemental-Table-3 for OlinkID -> gene symbol mapping ...")
    supp = pd.read_excel(supp_path, sheet_name="Global F-test")
    conversion = supp[["Assay", "OlinkID"]].drop_duplicates(subset="OlinkID")

    # Convert OlinkID columns to gene symbols
    id_to_gene = dict(zip(conversion["OlinkID"], conversion["Assay"]))
    id_cols = [c for c in olink.columns if c not in ("Public ID", "Public.ID", "day")]
    renamed = {}
    for col in id_cols:
        if col in id_to_gene:
            renamed[col] = id_to_gene[col]
    olink = olink.rename(columns=renamed)

    # Extract sample_id and day
    id_col = "Public ID" if "Public ID" in olink.columns else "Public.ID"
    olink["sample_id"] = olink[id_col].astype(str)

    # Parse day from sample_id (format: "123_D0", "123_D3", "123_D7")
    if "day" not in olink.columns:
        olink["day"] = olink["sample_id"].str.extract(r"_D(\d+)").astype(float)
    olink["subject_id"] = olink["sample_id"].str.replace(r"_D\d+", "", regex=True)

    # Merge clinical metadata
    clinical["subject_id"] = clinical["Public ID"].astype(str) if "Public ID" in clinical.columns else clinical["Public.ID"].astype(str)
    merged = olink.merge(clinical, on="subject_id", how="left", suffixes=("", "_clinical"))

    output_dir.mkdir(parents=True, exist_ok=True)
    results = {"timepoints": {}, "n_proteins": 0, "n_subjects": 0}

    protein_cols = [c for c in olink.columns if c not in (
        "sample_id", "subject_id", "day", id_col,
        "Public ID", "Public.ID",
    ) and not c.startswith("OID")]

    results["n_proteins"] = len(protein_cols)
    results["n_subjects"] = merged["subject_id"].nunique()

    # Split by day and save
    for day_val, group in merged.groupby("day"):
        day_label = f"D{int(day_val)}"
        out_cols = ["sample_id", "subject_id", "day"] + protein_cols
        # Add clinical columns if available
        for cc in ["Acuity", "WHO_Score_Admission", "Age", "Sex", "BMI"]:
            matches = [c for c in group.columns if c.lower().replace(" ", "_") == cc.lower() or c == cc]
            if matches and matches[0] not in out_cols:
                out_cols.append(matches[0])

        available_cols = [c for c in out_cols if c in group.columns]
        subset = group[available_cols].copy()
        out_path = output_dir / f"filbin_olink_{day_label}.csv.gz"
        subset.to_csv(out_path, index=False, compression="gzip")

        results["timepoints"][day_label] = {
            "n_samples": len(subset),
            "file": str(out_path),
        }
        if verbose:
            print(f"  {day_label}: {len(subset)} samples -> {out_path.name}")

    # Save metadata separately
    meta_cols = ["subject_id"]
    for cc in merged.columns:
        if cc not in protein_cols and cc not in ("sample_id", "day", id_col, "subject_id"):
            meta_cols.append(cc)
    meta_out = merged[meta_cols].drop_duplicates(subset="subject_id")
    meta_path_out = output_dir / "filbin_clinical_metadata.csv"
    meta_out.to_csv(meta_path_out, index=False)
    results["metadata_file"] = str(meta_path_out)

    if verbose:
        print(f"\nDone: {results['n_subjects']} subjects, {results['n_proteins']} proteins, "
              f"{len(results['timepoints'])} timepoints")
        print(f"Output: {output_dir}")

    return results


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Preprocess Filbin et al. (2021) COVID Olink data for proteomics clock"
    )
    parser.add_argument(
        "--data-dir", required=True,
        help="Directory containing downloaded Mendeley files (3 xlsx files)"
    )
    parser.add_argument(
        "--output", required=True,
        help="Output directory for per-timepoint CSVs"
    )
    parser.add_argument("--verbose", action="store_true", default=True)
    args = parser.parse_args()

    results = load_and_convert(
        data_dir=Path(args.data_dir),
        output_dir=Path(args.output),
        verbose=args.verbose,
    )

    import json
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
