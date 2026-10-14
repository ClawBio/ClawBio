#!/usr/bin/env python
"""
scripts/qc_filter.py

Stage 3: drop invalid or unusable variants, counting every drop by reason, then
sort by chromosome and position (summary statistics are often unsorted).

Reasons, checked in order (a row is counted once, under the first that applies):
  missing_required  CHR, BP, EA, NEA, BETA, SE or P empty
  invalid_p         P not in (0, 1] (checked exactly, so 1e-350 is valid)
  invalid_effect    BETA not numeric, or SE not a positive number
  invalid_allele    alleles not A/C/G/T strings, EA == NEA, or an indel with
                    --keep-indels false
  palindromic       A/T or C/G SNV: with --palindromic all, always; with
                    ambiguous, when EAF is missing or within 0.4-0.6; never
                    with none
  low_maf           min(EAF, 1-EAF) < --min-maf (rows without EAF are kept)
  duplicate         same CHR:BP and allele pair as another row; the smallest P is kept

Example:
    python scripts/qc_filter.py \\
        --input results/derive_effects/cohort.tsv \\
        --out results/qc_filter/cohort.tsv \\
        --summary-json results/qc_filter/cohort.summary.json \\
        --palindromic ambiguous --min-maf 0 --keep-indels true
"""
import argparse
import json
import logging
import os
import re
import sys
from collections import Counter
from decimal import Decimal

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
from harmonize_lib import (  # noqa: E402
    CANONICAL, chr_sort_key, is_palindromic, read_table, to_float, valid_p, write_table,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)

REQUIRED = ["CHR", "BP", "EA", "NEA", "BETA", "SE", "P"]
ALLELE_RE = re.compile(r"^[ACGT]+$")


def str2bool(value):
    v = str(value).strip().lower()
    if v in ("true", "1", "yes"):
        return True
    if v in ("false", "0", "no"):
        return False
    raise argparse.ArgumentTypeError(f"not a boolean: {value!r}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Normalised canonical TSV from derive_effects")
    p.add_argument("--out", required=True, help="Filtered, sorted canonical TSV")
    p.add_argument("--summary-json", required=True, help="Per-run summary JSON")
    p.add_argument("--palindromic", choices=["ambiguous", "all", "none"], required=True,
                   help="Which palindromic A/T, C/G SNVs to drop")
    p.add_argument("--min-maf", type=float, required=True, help="Minimum minor-allele frequency (0 = off)")
    p.add_argument("--keep-indels", type=str2bool, required=True, help="Keep multi-base alleles")
    return p.parse_args()


def drop_reason(row, args):
    if any(not row[c] for c in REQUIRED):
        return "missing_required"
    if not valid_p(row["P"]):
        return "invalid_p"
    se = to_float(row["SE"])
    if to_float(row["BETA"]) is None or se is None or se <= 0:
        return "invalid_effect"
    ea, nea = row["EA"], row["NEA"]
    if not (ALLELE_RE.match(ea) and ALLELE_RE.match(nea)) or ea == nea:
        return "invalid_allele"
    if not args.keep_indels and (len(ea) > 1 or len(nea) > 1):
        return "invalid_allele"
    eaf = to_float(row["EAF"])
    if eaf is not None and not 0 <= eaf <= 1:
        row["EAF"] = ""
        eaf = None
    if is_palindromic(ea, nea):
        if args.palindromic == "all":
            return "palindromic"
        if args.palindromic == "ambiguous" and (eaf is None or 0.4 <= eaf <= 0.6):
            return "palindromic"
    if args.min_maf > 0 and eaf is not None and min(eaf, 1 - eaf) < args.min_maf:
        return "low_maf"
    return None


def main():
    args = parse_args()
    rows = read_table(args.input)
    logger.info("Read %d rows from %s", len(rows), args.input)

    dropped = Counter()
    kept = []
    for row in rows:
        reason = drop_reason(row, args)
        if reason:
            dropped[reason] += 1
        else:
            kept.append(row)

    best = {}
    for i, row in enumerate(kept):
        key = (row["CHR"], row["BP"], tuple(sorted((row["EA"], row["NEA"]))))
        p = Decimal(row["P"])
        if key not in best or p < best[key][0]:
            best[key] = (p, i)
    keep_idx = {i for _, i in best.values()}
    if len(keep_idx) < len(kept):
        dropped["duplicate"] += len(kept) - len(keep_idx)
    kept = [r for i, r in enumerate(kept) if i in keep_idx]
    kept.sort(key=lambda r: chr_sort_key(r["CHR"], r["BP"]))

    for reason, n in sorted(dropped.items()):
        logger.info("dropped %s: %d", reason, n)
    logger.info("%d in, %d dropped, %d kept", len(rows), len(rows) - len(kept), len(kept))

    write_table(args.out, CANONICAL, kept)
    summary = {
        "stage": "qc_filter",
        "input": args.input,
        "rows_in": len(rows),
        "rows_out": len(kept),
        "dropped": dict(sorted(dropped.items())),
        "params": {"palindromic": args.palindromic, "min_maf": args.min_maf, "keep_indels": args.keep_indels},
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.summary_json)), exist_ok=True)
    with open(args.summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
