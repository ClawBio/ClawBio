#!/usr/bin/env python3
"""Resolve rsID -> variant coordinates via Ensembl REST.

If {skill_dir}/data/demo_{rsid}.json exists, use its 'variant' block (offline mode).
Otherwise hit https://rest.ensembl.org/variation/human/{rsid}.
"""
from __future__ import annotations

import argparse
import json
import sys
import urllib.request
import urllib.error
from pathlib import Path

ENSEMBL_URL = "https://rest.ensembl.org/variation/human/{rsid}?content-type=application/json"


def from_cache(skill_dir: Path, rsid: str) -> dict | None:
    cache = skill_dir / "data" / f"demo_{rsid}.json"
    if not cache.is_file():
        return None
    bundle = json.loads(cache.read_text())
    variant = bundle.get("variant")
    if not isinstance(variant, dict):
        return None
    return variant


def from_ensembl(rsid: str) -> dict:
    req = urllib.request.Request(ENSEMBL_URL.format(rsid=rsid), headers={"Accept": "application/json"})
    with urllib.request.urlopen(req, timeout=15) as r:
        data = json.load(r)

    grch38 = grch37 = None
    for m in data.get("mappings", []):
        asm = m.get("assembly_name", "")
        if asm == "GRCh38":
            grch38 = m
        elif asm == "GRCh37":
            grch37 = m

    if grch38 is None:
        raise SystemExit(f"Ensembl returned no GRCh38 mapping for {rsid}")

    allele_string = grch38.get("allele_string", "")
    parts = allele_string.split("/") if allele_string else []
    ref = parts[0] if parts else ""
    alt_alleles = parts[1:] if len(parts) > 1 else []

    variant = {
        "rsid": rsid,
        "chr": str(grch38.get("seq_region_name", "")),
        "pos_grch38": grch38.get("start"),
        "pos_grch37": grch37.get("start") if grch37 else None,
        "ref": ref,
        "alt": alt_alleles[0] if alt_alleles else "",
        "alt_alleles": alt_alleles,
        "allele_string": allele_string,
        "var_class": data.get("var_class", ""),
        "most_severe_consequence": data.get("most_severe_consequence", ""),
        "minor_allele": data.get("minor_allele", ""),
        "maf": data.get("MAF", data.get("minor_allele_freq")),
        "populations": [
            {"population": p.get("population", ""), "allele": p.get("allele", ""), "frequency": p.get("frequency")}
            for p in data.get("populations", [])
        ],
    }

    chr_ = variant["chr"]
    pos38 = variant["pos_grch38"]
    pos37 = variant["pos_grch37"]
    alt = variant["alt"]
    variant["variant_ids"] = {
        "open_targets": f"{chr_}_{pos38}_{ref}_{alt}" if pos38 else "",
        "gtex": f"chr{chr_}_{pos38}_{ref}_{alt}_b38" if pos38 else "",
        "pheweb": f"{chr_}:{pos38}-{ref}-{alt}" if pos38 else "",
        "bbj": f"{chr_}:{pos37}-{ref}-{alt}" if pos37 else "",
    }
    return variant


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--rsid", required=True)
    p.add_argument("--skill-dir", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    variant = from_cache(args.skill_dir, args.rsid)
    source = "cache"
    if variant is None:
        try:
            variant = from_ensembl(args.rsid)
            source = "ensembl"
        except urllib.error.URLError as e:
            print(f"[resolve] network error and no cache for {args.rsid}: {e}", file=sys.stderr)
            sys.exit(1)

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(variant, indent=2))
    print(f"[resolve] rsid={args.rsid} source={source} chr={variant.get('chr')} pos38={variant.get('pos_grch38')}")


if __name__ == "__main__":
    main()
