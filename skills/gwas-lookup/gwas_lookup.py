#!/usr/bin/env python3
"""
GWAS Lookup — Federated variant query across 9 genomic databases.

Inspired by Sasha Gusev's GWAS Lookup (https://sashagusev.github.io/gwas_lookup/).
Queries Ensembl, GWAS Catalog, Open Targets, UKB-TOPMed PheWeb, FinnGen,
Biobank Japan PheWeb, GTEx, EBI eQTL Catalogue, and LocusZoom PortalDev
for a single rsID.

Usage:
    python gwas_lookup.py --rsid rs3798220 --output results/
    python gwas_lookup.py --demo --output /tmp/gwas_demo
    python gwas_lookup.py --rsid rs429358 --skip gtex,bbj --output results/
"""

from __future__ import annotations

import argparse
import json
import sys
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from typing import Optional

# ---------------------------------------------------------------------------
# Shared library imports
# ---------------------------------------------------------------------------

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.report import write_result_json, DISCLAIMER

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

DEMO_RSID = "rs3798220"
SKILL_DIR = Path(__file__).resolve().parent
DEMO_DATA_PATH = SKILL_DIR / "data" / "demo_rs3798220.json"

DEFAULT_CACHE_DIR = Path.home() / ".clawbio" / "gwas_lookup_cache"
MAX_WORKERS = 8

ALL_API_NAMES = [
    "gwas_catalog", "open_targets", "open_targets_credsets",
    "pheweb_ukb", "finngen", "pheweb_bbj",
    "gtex", "eqtl_catalogue",
]

SKIP_ALIASES = {
    "gwas": "gwas_catalog",
    "ot": "open_targets",
    "ukb": "pheweb_ukb",
    "bbj": "pheweb_bbj",
    "eqtl": "eqtl_catalogue",
}


# ---------------------------------------------------------------------------
# API dispatch
# ---------------------------------------------------------------------------


def _fetch_gwas_catalog(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.gwas_catalog import get_associations
    return "gwas_catalog", get_associations(rsid, max_hits=max_hits, cache_dir=cache_dir, use_cache=use_cache)


def _fetch_open_targets(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.open_targets import get_variant
    chr_val = variant.get("chr", "")
    pos = variant.get("pos_grch38")
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")
    if not all([chr_val, pos, ref, alt]):
        return "open_targets", {"source": "open_targets", "status": "skipped", "message": "Missing coordinates"}
    return "open_targets", get_variant(chr_val, pos, ref, alt, cache_dir=cache_dir, use_cache=use_cache)


def _fetch_open_targets_credsets(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.open_targets import get_credible_sets
    chr_val = variant.get("chr", "")
    pos = variant.get("pos_grch38")
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")
    if not all([chr_val, pos, ref, alt]):
        return "open_targets_credsets", {"source": "open_targets_credsets", "status": "skipped", "message": "Missing coordinates"}
    return "open_targets_credsets", get_credible_sets(chr_val, pos, ref, alt, cache_dir=cache_dir, use_cache=use_cache)


def _fetch_pheweb_ukb(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.pheweb_ukb import get_phewas
    chr_val = variant.get("chr", "")
    pos = variant.get("pos_grch38")
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")
    if not all([chr_val, pos, ref, alt]):
        return "pheweb_ukb", {"source": "pheweb_ukb", "status": "skipped", "message": "Missing coordinates"}
    return "pheweb_ukb", get_phewas(chr_val, pos, ref, alt, cache_dir=cache_dir, use_cache=use_cache)


def _fetch_finngen(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.finngen import get_phewas
    chr_val = variant.get("chr", "")
    pos = variant.get("pos_grch38")
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")
    if not all([chr_val, pos, ref, alt]):
        return "finngen", {"source": "finngen", "status": "skipped", "message": "Missing coordinates"}
    return "finngen", get_phewas(chr_val, pos, ref, alt, cache_dir=cache_dir, use_cache=use_cache)


def _fetch_pheweb_bbj(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.pheweb_bbj import get_phewas
    chr_val = variant.get("chr", "")
    pos_37 = variant.get("pos_grch37")
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")
    if not all([chr_val, pos_37, ref, alt]):
        return "pheweb_bbj", {"source": "pheweb_bbj", "status": "skipped", "message": "Missing GRCh37 coordinates"}
    return "pheweb_bbj", get_phewas(chr_val, pos_37, ref, alt, cache_dir=cache_dir, use_cache=use_cache)


def _fetch_gtex(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.gtex import get_eqtls
    chr_val = variant.get("chr", "")
    pos = variant.get("pos_grch38")
    ref = variant.get("ref", "")
    alt = variant.get("alt", "")
    if not all([chr_val, pos, ref, alt]):
        return "gtex", {"source": "gtex", "status": "skipped", "message": "Missing coordinates"}
    return "gtex", get_eqtls(chr_val, pos, ref, alt, cache_dir=cache_dir, use_cache=use_cache)


def _fetch_eqtl_catalogue(rsid, variant, cache_dir, use_cache, max_hits):
    from gwas_lookup_api.eqtl_catalogue import get_associations
    return "eqtl_catalogue", get_associations(rsid, cache_dir=cache_dir, use_cache=use_cache)


API_DISPATCHERS = {
    "gwas_catalog": _fetch_gwas_catalog,
    "open_targets": _fetch_open_targets,
    "open_targets_credsets": _fetch_open_targets_credsets,
    "pheweb_ukb": _fetch_pheweb_ukb,
    "finngen": _fetch_finngen,
    "pheweb_bbj": _fetch_pheweb_bbj,
    "gtex": _fetch_gtex,
    "eqtl_catalogue": _fetch_eqtl_catalogue,
}


# ---------------------------------------------------------------------------
# Main lookup
# ---------------------------------------------------------------------------


def run_lookup(
    rsid: str,
    output_dir: Path,
    cache_dir: Path = DEFAULT_CACHE_DIR,
    use_cache: bool = True,
    skip_apis: list[str] | None = None,
    max_hits: int = 100,
    make_figures: bool = True,
    demo_data: dict | None = None,
) -> dict:
    """
    Run the full GWAS lookup pipeline for one rsID.

    Returns the merged results dict.
    """
    from gwas_lookup_core.resolve import resolve_variant
    from gwas_lookup_core.normalise import merge_all
    from gwas_lookup_core.report import generate_markdown, write_tables, generate_figures, write_reproducibility

    skip_set = set(skip_apis or [])
    # Normalize skip aliases
    skip_set = {SKIP_ALIASES.get(s, s) for s in skip_set}

    output_dir.mkdir(parents=True, exist_ok=True)

    # --- Step 1: Resolve variant ---
    if demo_data:
        variant = demo_data.get("variant", {})
        print(f"  Using pre-fetched demo data for {rsid}")
    else:
        print(f"  Resolving variant {rsid} via Ensembl...")
        variant = resolve_variant(rsid, cache_dir=cache_dir, use_cache=use_cache)

    chr_val = variant.get("chr", "?")
    pos_38 = variant.get("pos_grch38", "?")
    print(f"  Resolved: chr{chr_val}:{pos_38} ({variant.get('allele_string', '?')})")
    print(f"  Consequence: {variant.get('most_severe_consequence', '?')}")
    print()

    # --- Step 2: Parallel API queries ---
    if demo_data:
        api_results = demo_data.get("api_results", {})
        print(f"  Loaded {len(api_results)} pre-fetched API results")
    else:
        api_results = {}
        dispatchers_to_run = {
            name: fn for name, fn in API_DISPATCHERS.items()
            if name not in skip_set
        }

        print(f"  Querying {len(dispatchers_to_run)} APIs in parallel...")
        t0 = time.time()

        with ThreadPoolExecutor(max_workers=MAX_WORKERS) as executor:
            futures = {}
            for name, fn in dispatchers_to_run.items():
                future = executor.submit(fn, rsid, variant, cache_dir, use_cache, max_hits)
                futures[future] = name

            for future in as_completed(futures):
                name = futures[future]
                try:
                    result_name, result_data = future.result()
                    api_results[result_name] = result_data
                    status = result_data.get("status", "unknown")
                    print(f"    {result_name}: {status}")
                except Exception as e:
                    api_results[name] = {"source": name, "status": "error", "message": str(e)}
                    print(f"    {name}: ERROR — {e}")

        elapsed = time.time() - t0
        print(f"  All APIs queried in {elapsed:.1f}s")

    # Mark skipped APIs
    for name in skip_set:
        if name not in api_results:
            api_results[name] = {"source": name, "status": "skipped", "message": "Skipped by user"}

    print()

    # --- Step 3: Merge results ---
    print("  Merging and normalising results...")
    merged = merge_all(api_results)
    summary = merged.get("summary", {})
    print(f"    GWAS: {summary.get('total_gwas', 0)} associations "
          f"({summary.get('total_gwas_significant', 0)} GWS)")
    print(f"    PheWAS: UKB={summary.get('total_phewas_ukb', 0)}, "
          f"FinnGen={summary.get('total_phewas_finngen', 0)}, "
          f"BBJ={summary.get('total_phewas_bbj', 0)}")
    print(f"    eQTLs: {summary.get('total_eqtls', 0)}")
    print(f"    Credible sets: {summary.get('total_credible_sets', 0)}")
    print()

    # --- Step 4: Generate outputs ---
    print("  Writing report...")
    report_md = generate_markdown(variant, merged)
    (output_dir / "report.md").write_text(report_md)

    print("  Writing CSV tables...")
    write_tables(output_dir, merged)

    if make_figures:
        print("  Generating figures...")
        generate_figures(output_dir, merged, variant)

    print("  Writing reproducibility bundle...")
    write_reproducibility(output_dir, variant, list(skip_set))

    # Save raw JSON for debugging
    raw_path = output_dir / "raw_results.json"
    raw_path.write_text(json.dumps({
        "variant": variant,
        "api_results": api_results,
        "merged": merged,
    }, indent=2, default=str))

    # Write standardised result.json envelope
    print("  Writing result.json...")
    write_result_json(
        output_dir=output_dir,
        skill="gwas-lookup",
        version="0.2.0",
        summary={
            "rsid": rsid,
            "chr": chr_val,
            "pos_grch38": pos_38,
            "total_gwas": summary.get("total_gwas", 0),
            "total_gwas_significant": summary.get("total_gwas_significant", 0),
            "total_phewas_ukb": summary.get("total_phewas_ukb", 0),
            "total_phewas_finngen": summary.get("total_phewas_finngen", 0),
            "total_phewas_bbj": summary.get("total_phewas_bbj", 0),
            "total_eqtls": summary.get("total_eqtls", 0),
            "total_credible_sets": summary.get("total_credible_sets", 0),
            "apis_queried": len(api_results),
            "apis_skipped": len(skip_set),
        },
        data={
            "variant": variant,
            "merged": merged,
        },
    )

    print(f"\n  Report: {output_dir / 'report.md'}")
    print(f"  Full output: {output_dir}/")
    print(f"\n  {DISCLAIMER}")

    return merged


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="GWAS Lookup — federated variant query across 9 genomic databases"
    )
    parser.add_argument("--rsid", help="rsID to look up (e.g., rs3798220)")
    parser.add_argument("--demo", action="store_true", help=f"Run with pre-fetched data for {DEMO_RSID}")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--skip", default="", help="Comma-separated API names to skip (e.g., gtex,bbj)")
    parser.add_argument("--no-figures", action="store_true", help="Skip matplotlib figures")
    parser.add_argument("--no-cache", action="store_true", help="Bypass local cache")
    parser.add_argument("--cache-dir", default=str(DEFAULT_CACHE_DIR), help="Cache directory")
    parser.add_argument("--max-hits", type=int, default=100, help="Max GWAS associations (default: 100)")

    args = parser.parse_args()

    if not args.rsid and not args.demo:
        parser.print_help()
        print("\nError: provide --rsid or --demo")
        sys.exit(1)

    rsid = args.rsid or DEMO_RSID
    skip_apis = [s.strip() for s in args.skip.split(",") if s.strip()]
    output_dir = Path(args.output)

    print(f"GWAS Lookup: {rsid}")
    print("=" * 60)
    print()

    # Demo mode: load pre-fetched data
    demo_data = None
    if args.demo:
        if DEMO_DATA_PATH.exists():
            demo_data = json.loads(DEMO_DATA_PATH.read_text())
            print(f"  Demo mode: loading {DEMO_DATA_PATH.name}")
        else:
            print(f"  Demo data not found at {DEMO_DATA_PATH}")
            print(f"  Running live query for {DEMO_RSID} instead")
        print()

    run_lookup(
        rsid=rsid,
        output_dir=output_dir,
        cache_dir=Path(args.cache_dir),
        use_cache=not args.no_cache,
        skip_apis=skip_apis,
        max_hits=args.max_hits,
        make_figures=not args.no_figures,
        demo_data=demo_data,
    )


if __name__ == "__main__":
    main()
