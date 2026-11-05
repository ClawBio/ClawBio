#!/usr/bin/env python3
"""Merge resolve + 8 API source files into a single unified JSON record."""
from __future__ import annotations

import argparse
import json
from pathlib import Path


def main() -> None:
    p = argparse.ArgumentParser()
    p.add_argument("--variant", required=True, type=Path)
    p.add_argument("--gwas-catalog", required=True, type=Path)
    p.add_argument("--open-targets", required=True, type=Path)
    p.add_argument("--open-targets-credsets", required=True, type=Path)
    p.add_argument("--pheweb-ukb", required=True, type=Path)
    p.add_argument("--finngen", required=True, type=Path)
    p.add_argument("--pheweb-bbj", required=True, type=Path)
    p.add_argument("--gtex", required=True, type=Path)
    p.add_argument("--eqtl-catalogue", required=True, type=Path)
    p.add_argument("--out", required=True, type=Path)
    args = p.parse_args()

    variant = json.loads(args.variant.read_text())
    sources = {
        "gwas_catalog": json.loads(args.gwas_catalog.read_text()),
        "open_targets": json.loads(args.open_targets.read_text()),
        "open_targets_credsets": json.loads(args.open_targets_credsets.read_text()),
        "pheweb_ukb": json.loads(args.pheweb_ukb.read_text()),
        "finngen": json.loads(args.finngen.read_text()),
        "pheweb_bbj": json.loads(args.pheweb_bbj.read_text()),
        "gtex": json.loads(args.gtex.read_text()),
        "eqtl_catalogue": json.loads(args.eqtl_catalogue.read_text()),
    }

    statuses = {k: v.get("status", "?") for k, v in sources.items()}
    ok_count = sum(1 for s in statuses.values() if s == "ok")

    merged = {
        "rsid": variant.get("rsid"),
        "variant": variant,
        "sources": sources,
        "summary": {
            "total_sources": len(sources),
            "ok_sources": ok_count,
            "statuses": statuses,
        },
    }

    args.out.parent.mkdir(parents=True, exist_ok=True)
    args.out.write_text(json.dumps(merged, indent=2))
    print(f"[merge] rsid={variant.get('rsid')} ok_sources={ok_count}/{len(sources)}")


if __name__ == "__main__":
    main()
