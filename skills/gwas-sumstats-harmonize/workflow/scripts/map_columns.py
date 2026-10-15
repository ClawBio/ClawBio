#!/usr/bin/env python
"""
scripts/map_columns.py

Stage 1: map a raw summary-statistics file's columns onto canonical names
(SNP CHR BP EA NEA EAF BETA SE P N, plus OR / LOG10P / Z helpers when present).

Column names are detected from a table of known aliases (GWAS-Catalog SSF,
PLINK 1.9/2, REGENIE, METAL, BOLT, SAIGE-style headers); an explicit mapping
always wins. Values are copied verbatim: normalisation and derivation happen in
derive_effects. Two format quirks are resolved here because only the raw
header shows them:
  - PLINK2 REF/ALT/A1: EA = A1, NEA = whichever of REF/ALT A1 is not, per row.
  - No CHR/BP columns (e.g. METAL): parsed from a "chr:pos" SNP id.
  - PLINK TEST column: only TEST == ADD rows are kept (covariate and other
    model rows share the variant's CHR:BP and often have a smaller P); the
    rest are counted as non_additive_test. A TEST column with no ADD rows is
    an error rather than a silently empty output.
This is the only stage besides qc_filter that drops rows.

Example:
    python scripts/map_columns.py \\
        --input input/cohort.glm.logistic \\
        --columns-json results/tmp/cohort.columns.json \\
        --out results/map_columns/cohort.tsv \\
        --summary-json results/map_columns/cohort.summary.json
"""
import argparse
import json
import logging
import os
import sys

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
from harmonize_lib import (  # noqa: E402
    CANONICAL, HELPERS, HarmonizeError, check_derivable, detect_columns, read_raw,
    split_snp_position, write_table,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Raw summary statistics (tsv/csv/whitespace, optional .gz)")
    p.add_argument("--columns-json", help="JSON {canonical: header_name} explicit mapping (optional)")
    p.add_argument("--out", required=True, help="Canonical-column TSV output")
    p.add_argument("--summary-json", required=True, help="Per-run summary JSON")
    return p.parse_args()


def main():
    args = parse_args()
    explicit = {}
    if args.columns_json:
        with open(args.columns_json) as fh:
            explicit = json.load(fh) or {}

    header, rows = read_raw(args.input)
    rows_in = len(rows)
    logger.info("Read %d rows, %d columns from %s", len(rows), len(header), args.input)
    try:
        detected = detect_columns(header, explicit)
        check_derivable(detected)
    except HarmonizeError as exc:
        logger.error("%s", exc)
        sys.exit(2)

    mapping = detected["mapping"]
    idx = {h: i for i, h in enumerate(header)}
    other = detected["plink2_other"]
    parse_pos = not ("CHR" in mapping and "BP" in mapping)
    columns = CANONICAL + [h for h in HELPERS if h in mapping]
    out_rows, short_rows, unparsed_pos = [], 0, 0
    dropped = {}

    test_col = detected["test_column"]
    if test_col is not None:
        tests = [(r[idx[test_col]] if len(r) > idx[test_col] else "").upper() for r in rows]
        if rows and "ADD" not in tests:
            logger.error("TEST column '%s' has no ADD rows (saw: %s); only additive-model results "
                         "are harmonized", test_col, sorted(set(tests))[:10])
            sys.exit(2)
        n_other = sum(t != "ADD" for t in tests)
        rows = [r for r, t in zip(rows, tests) if t == "ADD"]
        if n_other:
            dropped["non_additive_test"] = n_other

    for raw in rows:
        if len(raw) < len(header):
            short_rows += 1
            raw = raw + [""] * (len(header) - len(raw))
        row = {canon: raw[idx[name]] for canon, name in mapping.items()}
        if other:
            ref, alt = raw[idx[other[0]]], raw[idx[other[1]]]
            row["NEA"] = alt if row["EA"].upper() == ref.upper() else ref
        if parse_pos:
            chrom, bp = split_snp_position(row.get("SNP", ""))
            row.setdefault("CHR", chrom or "")
            row.setdefault("BP", bp or "")
            if chrom is None:
                unparsed_pos += 1
        out_rows.append(row)

    notes = list(detected["notes"])
    if parse_pos:
        notes.append(f"CHR/BP parsed from SNP ids ({unparsed_pos} unparseable)")
    if short_rows:
        notes.append(f"{short_rows} rows had fewer fields than the header (padded empty)")
    for n in notes:
        logger.info("note: %s", n)

    write_table(args.out, columns, out_rows)
    summary = {
        "stage": "map_columns",
        "input": args.input,
        "rows_in": rows_in,
        "rows_out": len(out_rows),
        "mapping": mapping,
        "plink2_other": list(other) if other else None,
        "build": detected["build"],
        "dropped": dropped,
        "notes": notes,
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.summary_json)), exist_ok=True)
    with open(args.summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("%d in, %d dropped, %d kept -> %s", rows_in, rows_in - len(out_rows), len(out_rows), args.out)


if __name__ == "__main__":
    main()
