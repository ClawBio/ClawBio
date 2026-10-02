#!/usr/bin/env python
"""
scripts/align_reference.py

Stage 4: align alleles to a reference panel so that EA = reference ALT.

Reference: tab-separated with header CHR BP REF ALT [AF], AF = ALT frequency,
same genome build as the summary statistics. Per variant (matched on CHR:BP):
  aligned                 EA/NEA = ALT/REF                 unchanged
  swapped                 EA/NEA = REF/ALT                 swap, BETA = -BETA, EAF = 1 - EAF
  strand_flipped          revcomp(EA/NEA) = ALT/REF        reverse-complement alleles
  strand_flipped_swapped  revcomp(EA/NEA) = REF/ALT        reverse-complement + swap
  mismatch                anything else                    dropped
  not_in_reference        no reference row at CHR:BP       dropped with --drop-unmatched true, else kept
Palindromic SNVs (A/T, C/G) cannot be classified from alleles: reverse-strand
A/T looks exactly like a forward swap. They are resolved by frequency (EAF and
the reference ALT frequency on the same side of 0.5 -> EA is ALT) and dropped
as palindromic_unresolved when EAF or AF is missing or within 0.4-0.6.
Missing EAF is filled from the reference AF only after alignment, so the
frequency always belongs to the effect allele.

Without --reference the stage passes rows through unchanged.

Example:
    python scripts/align_reference.py \\
        --input results/qc_filter/cohort.tsv \\
        --reference resources/reference.tsv \\
        --out results/align_reference/cohort.tsv \\
        --summary-json results/align_reference/cohort.summary.json \\
        --drop-unmatched true
"""
import argparse
import json
import logging
import os
import shutil
import sys
from collections import Counter, defaultdict

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
from harmonize_lib import (  # noqa: E402
    CANONICAL, classify_alignment, complement, fmt, normalize_chr, read_table, resolve_palindromic,
    to_float, write_table,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def str2bool(value):
    v = str(value).strip().lower()
    if v in ("true", "1", "yes"):
        return True
    if v in ("false", "0", "no"):
        return False
    raise argparse.ArgumentTypeError(f"not a boolean: {value!r}")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Filtered canonical TSV from qc_filter")
    p.add_argument("--reference", help="Reference TSV (CHR BP REF ALT [AF]); omit to pass through")
    p.add_argument("--out", required=True, help="Aligned canonical TSV")
    p.add_argument("--summary-json", required=True, help="Per-run summary JSON")
    p.add_argument("--drop-unmatched", type=str2bool, required=True,
                   help="Drop variants absent from the reference")
    return p.parse_args()


def load_reference(path):
    ref = defaultdict(list)
    for r in read_table(path):
        ref[(normalize_chr(r["CHR"]), str(int(float(r["BP"]))))].append(
            (r["REF"].upper(), r["ALT"].upper(), r.get("AF", ""))
        )
    return ref


def align(row, candidates, counts):
    """Return the aligned row, or None to drop it."""
    unresolved = False
    for ref, alt, af in candidates:
        case = classify_alignment(row["EA"], row["NEA"], ref, alt)
        if case == "palindromic":
            case = resolve_palindromic(row["EA"], row["NEA"], to_float(row["EAF"]), ref, alt, to_float(af))
            if case is None:
                unresolved = True
                continue
        if case == "mismatch":
            continue
        counts[case] += 1
        if case.startswith("strand_flipped"):
            row["EA"], row["NEA"] = complement(row["EA"]), complement(row["NEA"])
        if case.endswith("swapped"):
            row["EA"], row["NEA"] = row["NEA"], row["EA"]
            row["BETA"] = fmt(-float(row["BETA"]))
            eaf = to_float(row["EAF"])
            if eaf is not None:
                row["EAF"] = fmt(1 - eaf)
        if not row["EAF"] and to_float(af) is not None:
            row["EAF"] = af
            counts["eaf_filled"] += 1
        return row
    counts["palindromic_unresolved" if unresolved else "mismatch"] += 1
    return None


def main():
    args = parse_args()
    rows = read_table(args.input)
    logger.info("Read %d rows from %s", len(rows), args.input)
    counts = Counter()

    if not args.reference:
        os.makedirs(os.path.dirname(os.path.abspath(args.out)), exist_ok=True)
        shutil.copyfile(args.input, args.out)
        kept = rows
        logger.info("No reference given: passing %d rows through", len(rows))
    else:
        reference = load_reference(args.reference)
        logger.info("Loaded %d reference positions from %s", len(reference), args.reference)
        kept = []
        for row in rows:
            candidates = reference.get((row["CHR"], row["BP"]))
            if not candidates:
                counts["not_in_reference"] += 1
                if not args.drop_unmatched:
                    kept.append(row)
                continue
            aligned = align(row, candidates, counts)
            if aligned is not None:
                kept.append(aligned)
        write_table(args.out, CANONICAL, kept)
        for k, v in sorted(counts.items()):
            logger.info("%s: %d", k, v)
        if rows and counts["not_in_reference"] > 0.5 * len(rows):
            logger.warning(
                "%d of %d variants are not in the reference: check that both use the same genome build",
                counts["not_in_reference"], len(rows),
            )

    logger.info("%d in, %d dropped, %d kept", len(rows), len(rows) - len(kept), len(kept))
    summary = {
        "stage": "align_reference",
        "input": args.input,
        "reference": args.reference,
        "rows_in": len(rows),
        "rows_out": len(kept),
        "alignment": dict(sorted(counts.items())),
        "params": {"drop_unmatched": args.drop_unmatched},
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.summary_json)), exist_ok=True)
    with open(args.summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)


if __name__ == "__main__":
    main()
