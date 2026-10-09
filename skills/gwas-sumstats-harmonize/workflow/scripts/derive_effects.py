#!/usr/bin/env python
"""
scripts/derive_effects.py

Stage 2: normalise values and fill derivable fields, writing exactly the
canonical columns (SNP CHR BP EA NEA EAF BETA SE P N).

  - CHR: strip "chr", 23/24/25/26 -> X/Y/XY/MT; alleles upper-cased
  - P:    from LOG10P (exact string, survives float underflow), else from Z,
          else from BETA/SE
  - BETA: ln(OR) when only an odds ratio is given
  - BETA, SE from Z, EAF and N when there is no BETA or OR: the standardised
          (per-SD) scale of Zhu et al. 2016, b = z / sqrt(2p(1-p)(n + z^2)),
          se = 1 / sqrt(2p(1-p)(n + z^2)); rows without EAF or N stay empty
  - SE:   |BETA| / z(P) when missing
  - SNP:  CHR:BP:NEA:EA when missing
NA-like tokens (NA, nan, ., null) become empty. Nothing is dropped here;
invalid rows are removed, with reasons, by qc_filter.

Example:
    python scripts/derive_effects.py \\
        --input results/map_columns/cohort.tsv \\
        --out results/derive_effects/cohort.tsv \\
        --summary-json results/derive_effects/cohort.summary.json
"""
import argparse
import json
import logging
import math
import os
import sys
from collections import Counter
from decimal import InvalidOperation

sys.path.insert(0, os.path.join(os.path.dirname(__file__), "lib"))
from harmonize_lib import (  # noqa: E402
    CANONICAL, beta_se_from_z, fmt, is_missing, normalize_chr, p_from_beta_se, p_from_log10p, p_from_z,
    read_table, se_from_beta_p, to_float, valid_p, write_table,
)

logging.basicConfig(level=logging.INFO, format="%(asctime)s %(levelname)s %(message)s")
logger = logging.getLogger(__name__)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="Canonical-column TSV from map_columns")
    p.add_argument("--out", required=True, help="Normalised canonical TSV")
    p.add_argument("--summary-json", required=True, help="Per-run summary JSON")
    return p.parse_args()


def derive(row, counts):
    out = {c: ("" if is_missing(row.get(c)) else row.get(c).strip()) for c in CANONICAL}
    if out["CHR"]:
        out["CHR"] = normalize_chr(out["CHR"])
    if out["BP"]:
        bp = to_float(out["BP"])
        out["BP"] = str(int(bp)) if bp is not None and bp == int(bp) else out["BP"]
    out["EA"], out["NEA"] = out["EA"].upper(), out["NEA"].upper()

    if not out["P"] and not is_missing(row.get("LOG10P")):
        try:
            out["P"] = p_from_log10p(row["LOG10P"])
            counts["P_from_LOG10P"] += 1
        except (InvalidOperation, ValueError, OverflowError):
            # Non-numeric or infinite LOG10P: leave P empty (qc_filter drops the
            # row as missing_required) and count it so the report says why.
            counts["LOG10P_unparseable"] += 1
    if not out["P"] and to_float(row.get("Z")) is not None:
        out["P"] = fmt(p_from_z(to_float(row["Z"])))
        counts["P_from_Z"] += 1

    beta = to_float(out["BETA"])
    if beta is None:
        or_ = to_float(row.get("OR"))
        if or_ is not None and or_ > 0:
            beta = math.log(or_)
            out["BETA"] = fmt(beta)
            counts["BETA_from_OR"] += 1
    z, eaf, n = to_float(row.get("Z")), to_float(out["EAF"]), to_float(out["N"])
    if beta is None and z is not None and eaf is not None and 0 < eaf < 1 and n is not None and n > 0:
        beta, se_z = beta_se_from_z(z, eaf, n)
        out["BETA"], out["SE"] = fmt(beta), fmt(se_z)
        counts["BETA_SE_from_Z"] += 1

    se = to_float(out["SE"])
    p = to_float(out["P"])
    if se is None and beta is not None and p is not None and 0 < p < 1 and beta != 0:
        se = se_from_beta_p(beta, p)
        if se:
            out["SE"] = fmt(se)
            counts["SE_from_P"] += 1
    if not out["P"] and beta is not None and se:
        out["P"] = fmt(p_from_beta_se(beta, se))
        counts["P_from_BETA_SE"] += 1

    if not out["SNP"] and out["CHR"] and out["BP"]:
        out["SNP"] = f"{out['CHR']}:{out['BP']}:{out['NEA']}:{out['EA']}"
        counts["SNP_filled"] += 1
    if out["P"] and not valid_p(out["P"]):
        counts["P_out_of_range"] += 1
    return out


def main():
    args = parse_args()
    rows = read_table(args.input)
    logger.info("Read %d rows from %s", len(rows), args.input)
    counts = Counter()
    out_rows = [derive(r, counts) for r in rows]
    for k, v in sorted(counts.items()):
        logger.info("derived %s: %d", k, v)

    write_table(args.out, CANONICAL, out_rows)
    summary = {
        "stage": "derive_effects",
        "input": args.input,
        "rows_in": len(rows),
        "rows_out": len(out_rows),
        "derived": dict(sorted(counts.items())),
    }
    os.makedirs(os.path.dirname(os.path.abspath(args.summary_json)), exist_ok=True)
    with open(args.summary_json, "w") as fh:
        json.dump(summary, fh, indent=2)
    logger.info("%d in, 0 dropped, %d kept -> %s", len(rows), len(out_rows), args.out)


if __name__ == "__main__":
    main()
