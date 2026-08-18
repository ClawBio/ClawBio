"""Fetch a build-matched evidence layer for GRCh37 records, and cache it.

Why the GRCh37 endpoint specifically
------------------------------------
The teaching data is GRCh37/b37. gnomAD v4 is GRCh38-native, so joining this
data against v4 requires a liftover, and a failed liftover looks exactly like a
variant with no frequency record. That failure mode is the one the brief singles
out: "Absence of a frequency is not evidence of rarity."

Querying ``grch37.rest.ensembl.org`` removes the liftover from the path
entirely, so an empty frequency result means Ensembl has no record — not that
our coordinates missed.

What this costs in honesty, stated here and repeated in the report: the GRCh37
endpoint serves gnomAD **exomes r2.1.1**, which is exome-only and far smaller
than v4. Absence here is therefore weaker evidence than absence in v4 would be.
We do not use that weakness to imply anything about a variant's frequency.

One batched POST covers the whole file (the endpoint accepts up to 200 variants
per request), so this stays well inside the 15 requests/second limit. The result
is cached to disk and committed, so no downstream step depends on a live call.
"""

from __future__ import annotations

import argparse
import csv
import json
import pathlib
import sys
import urllib.error
import urllib.request

ENDPOINT = "https://grch37.rest.ensembl.org/vep/human/region"
# Kept deliberately minimal, and verified parameter by parameter against the live
# endpoint. The GRCh37 service rejects several parameters the GRCh38 one accepts
# — `mane`, `af_gnomade`, `af_gnomadg`, `af_1kg` — and answers **503 Service
# Unavailable** rather than 400 Bad Request when given them. That is very easy to
# misread as the service being down and to write up as "evidence layer
# unavailable". It was, briefly, in this project. `af=1` alone returns the
# colocated-variant frequencies we use.
PARAMS = "canonical=1&numbers=1&variant_class=1&af=1"
BATCH = 200


def load_variants(tsv: pathlib.Path) -> list[dict]:
    with tsv.open() as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


def to_region_string(row: dict) -> str:
    """VCF-style whitespace-separated record, as the region endpoint expects."""
    ident = row["ID"] if row.get("ID") and row["ID"] != "." else "."
    return f"{row['CHROM']} {row['POS']} {ident} {row['REF']} {row['ALT']} . . ."


def fetch(variants: list[str], *, timeout: int = 180) -> list[dict]:
    out: list[dict] = []
    for i in range(0, len(variants), BATCH):
        chunk = variants[i : i + BATCH]
        req = urllib.request.Request(
            f"{ENDPOINT}?{PARAMS}",
            data=json.dumps({"variants": chunk}).encode(),
            headers={"Content-Type": "application/json", "Accept": "application/json"},
            method="POST",
        )
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            out.extend(json.load(resp))
    return out


def main() -> int:
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--input", required=True, type=pathlib.Path, help="segregation TSV")
    ap.add_argument("--output", required=True, type=pathlib.Path, help="cache JSON path")
    args = ap.parse_args()

    rows = load_variants(args.input)
    variants = [to_region_string(r) for r in rows]
    args.output.parent.mkdir(parents=True, exist_ok=True)

    try:
        data = fetch(variants)
    except (urllib.error.URLError, TimeoutError, OSError) as exc:
        # A failed fetch must not be indistinguishable from an empty result.
        payload = {"status": "unavailable", "reason": str(exc), "requested": len(variants)}
        args.output.write_text(json.dumps(payload, indent=1))
        print(f"evidence layer unavailable: {exc}", file=sys.stderr)
        return 1

    args.output.write_text(
        json.dumps(
            {
                "status": "ok",
                "endpoint": ENDPOINT,
                "params": PARAMS,
                "assembly": "GRCh37",
                "frequency_source": "gnomAD exomes r2.1.1 via Ensembl GRCh37 REST (exome-only)",
                "requested": len(variants),
                "returned": len(data),
                "records": data,
            },
            indent=1,
        )
    )
    print(f"cached {len(data)}/{len(variants)} records -> {args.output}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
