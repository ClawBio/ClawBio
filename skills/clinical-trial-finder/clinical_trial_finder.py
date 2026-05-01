#!/usr/bin/env python3
"""clinical-trial-finder: query ClinicalTrials.gov API v2 for active trials.

Entry point and CLI only.  Business logic lives in:
  constants.py  -- shared constants, FHIR mappings, labels
  api.py        -- ClinicalTrials.gov API client and data normalisation
  writers.py    -- all output generation (report, JSON, FHIR, chart, repro)
  opentargets.py -- OpenTargets gene-to-disease GraphQL client
"""

import argparse
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from api import fetch_trials, parse_input
from constants import ALL_STATUSES, DEFAULT_PAGE_SIZE, DEMO_DATA
from writers import (
    count_recruiting,
    validate_fhir_bundle,
    write_checksums,
    write_commands,
    write_csv,
    write_fhir_bundle,
    write_html,
    write_phase_chart,
    write_report,
    write_summary,
)

# ---------------------------------------------------------------------------
# CLI argument parser
# ---------------------------------------------------------------------------


def _build_parser() -> argparse.ArgumentParser:
    """Build argparse parser with 5 mutually exclusive input modes.

    Modes: --input (file), --query (string), --gene (OpenTargets),
           --rsid (gwas-lookup), --demo.
    """
    p = argparse.ArgumentParser(
        description="Clinical Trial Finder -- ClinicalTrials.gov API v2 + OpenTargets"
    )

    # Input source -- exactly one required
    src = p.add_mutually_exclusive_group()
    src.add_argument("--input", type=Path, help="Query file (one search term per line)")
    src.add_argument("--query", type=str, help="Direct search query string")
    src.add_argument(
        "--gene", type=str, help="Gene symbol (e.g. BRCA1) -- enriched via OpenTargets"
    )
    src.add_argument(
        "--rsid",
        type=str,
        help="rsID (e.g. rs3798220) -- resolves via gwas-lookup to traits + genes",
    )
    src.add_argument(
        "--demo", action="store_true", help="Run with built-in demo data (BRCA1)"
    )

    # Output and filtering
    p.add_argument(
        "--output", type=Path, default=Path("/tmp/clinical_trial_finder_output")
    )
    p.add_argument("--max-results", type=int, default=DEFAULT_PAGE_SIZE)
    p.add_argument(
        "--status",
        type=str,
        default=None,
        choices=ALL_STATUSES,
        help="Filter trials by recruitment status (default: show all)",
    )
    p.add_argument("--fhir", action="store_true", help="Also write fhir_bundle.json")
    p.add_argument(
        "--country",
        type=str,
        default=None,
        help="Filter by country (e.g. 'United States', 'Spain', 'DE')",
    )
    p.add_argument(
        "--euctr",
        action="store_true",
        help="Also search EU Clinical Trials Register (EUCTR)",
    )

    # OpenTargets parameters (only used with --gene)
    p.add_argument(
        "--ot-min-score",
        type=float,
        default=0.6,
        help="OpenTargets association score threshold (default: 0.6)",
    )
    p.add_argument(
        "--ot-max-diseases",
        type=int,
        default=5,
        help="Max diseases to query per gene via OpenTargets (default: 5)",
    )
    return p


# ---------------------------------------------------------------------------
# Main
# ---------------------------------------------------------------------------


def main() -> None:
    """Entry point: resolve input mode, fetch trials, write all outputs."""
    args = _build_parser().parse_args()

    gene_context: dict | None = None
    trials: list[dict] = []
    recruiting: int = 0

    # --- Resolve input source ---

    if args.demo:
        query_info = parse_input(DEMO_DATA)

    elif args.input:
        if not args.input.exists():
            raise SystemExit(f"Input file not found: {args.input}")
        try:
            query_info = parse_input(args.input)
        except ValueError as exc:
            raise SystemExit(str(exc)) from None

    elif args.query:
        query_info = {"query": args.query, "terms": [args.query]}

    elif args.gene:
        # Gene mode: resolve gene via OpenTargets, get associated diseases,
        # then query CT.gov once per disease and deduplicate by NCT ID.
        # This avoids the HTTP 400 that CT.gov returns when you OR too
        # many conditions in a single query.cond string.
        import opentargets  # local module -- only imported when --gene is used

        symbol = args.gene.upper()
        print(f"Resolving {symbol!r} via OpenTargets...")
        try:
            ensembl_id, gene_name = opentargets.resolve_gene(symbol)
        except ValueError as exc:
            raise SystemExit(str(exc)) from None
        diseases = opentargets.get_diseases(
            ensembl_id,
            min_score=args.ot_min_score,
            max_results=args.ot_max_diseases,
        )
        if not diseases:
            raise SystemExit(
                f"No disease associations found for {symbol} "
                f"with score >= {args.ot_min_score}. Try --ot-min-score 0.3"
            )
        disease_names = [d.name for d in diseases]
        print(f"  {symbol} -> {len(diseases)} diseases: {', '.join(disease_names)}")

        # Distribute max_results across diseases, deduplicate by NCT ID
        per_disease = max(1, args.max_results // len(diseases))
        seen: set[str] = set()
        trials = []
        for disease in disease_names:
            print(f"  Querying trials for: {disease!r}")
            for t in fetch_trials(
                disease, max_results=per_disease, country=args.country
            ):
                if t["nct_id"] not in seen:
                    seen.add(t["nct_id"])
                    trials.append(t)

        query_info = {"query": symbol, "terms": disease_names}
        gene_context = {
            "symbol": symbol,
            "name": gene_name,
            "diseases": disease_names,
            "min_score": args.ot_min_score,
        }
        recruiting = count_recruiting(trials)
        print(f"  Found {len(trials)} unique trials ({recruiting} recruiting)")

    elif args.rsid:
        # rsID mode: call gwas-lookup to resolve variant -> traits + genes,
        # then query CT.gov for each trait.  Connects the pipeline:
        # variant -> GWAS trait -> clinical trials
        import gwas_bridge

        rsid = args.rsid.lower()
        print(f"Resolving {rsid!r} via gwas-lookup...")
        try:
            gwas_result = gwas_bridge.resolve_rsid(rsid)
        except ValueError as exc:
            raise SystemExit(str(exc)) from None

        traits = gwas_result["traits"]
        genes = gwas_result["genes"]
        print(f"  {rsid} -> {len(traits)} traits: {', '.join(traits)}")
        if genes:
            print(f"  eQTL genes: {', '.join(genes[:5])}")

        # Query CT.gov for each trait, deduplicate by NCT ID
        per_trait = max(1, args.max_results // max(len(traits), 1))
        seen: set[str] = set()
        trials = []
        for trait in traits:
            print(f"  Querying trials for: {trait!r}")
            for t in fetch_trials(trait, max_results=per_trait, country=args.country):
                if t["nct_id"] not in seen:
                    seen.add(t["nct_id"])
                    trials.append(t)

        query_info = {"query": rsid, "terms": traits}
        gene_context = {
            "symbol": rsid,
            "name": f"GWAS variant ({', '.join(genes[:3])})"
            if genes
            else "GWAS variant",
            "diseases": traits,
            "min_score": 0.0,
        }
        recruiting = count_recruiting(trials)
        print(f"  Found {len(trials)} unique trials ({recruiting} recruiting)")

    else:
        _build_parser().error("Provide --input, --query, --gene, --rsid, or --demo")

    # --- Fetch trials (non-gene/non-rsid modes) ---

    if not args.gene and not args.rsid:
        country_msg = f" in {args.country}" if args.country else ""
        print(f"Querying ClinicalTrials.gov: {query_info['query']!r}{country_msg}")
        trials = fetch_trials(
            query_info["query"], max_results=args.max_results, country=args.country
        )
        recruiting = count_recruiting(trials)
        print(f"Found {len(trials)} trials ({recruiting} recruiting)")

    # --- Optional EUCTR merge (step 14) ---

    if args.euctr:
        import euctr

        query_str = query_info["query"]
        print(f"Querying EU Clinical Trials Register: {query_str!r}")
        eu_trials = euctr.fetch_euctr(query_str)
        if eu_trials:
            seen_ids = {t["nct_id"] for t in trials}
            added = [t for t in eu_trials if t["nct_id"] not in seen_ids]
            trials.extend(added)
            print(f"  EUCTR: {len(eu_trials)} found, {len(added)} new (merged)")
        else:
            print(
                "  EUCTR: no results or API unavailable (continuing with CT.gov only)"
            )

    # --- Optional status filter (applied post-fetch, after merge) ---

    if args.status:
        trials = [t for t in trials if t["status"] == args.status]
        print(f"Filtered by {args.status}: {len(trials)} trials remaining")

    # --- Generate outputs ---

    args.output.mkdir(parents=True, exist_ok=True)

    report = write_report(query_info, trials, args.output, gene_context)
    print(f"Report  -> {report}")

    html = write_html(query_info, trials, args.output, gene_context)
    print(f"HTML    -> {html}")

    summary = write_summary(query_info, trials, args.output)
    print(f"Summary -> {summary}")

    csv_path = write_csv(trials, args.output)
    print(f"CSV     -> {csv_path}")

    if args.fhir:
        bundle = write_fhir_bundle(trials, args.output)
        print(f"FHIR R4 -> {bundle}")

        # Validate the bundle we just wrote
        import json as _json

        bundle_data = _json.loads(bundle.read_text())
        errors = validate_fhir_bundle(bundle_data)
        if errors:
            print(f"  FHIR validation: {len(errors)} issues")
            for e in errors[:5]:
                print(f"    - {e}")
        else:
            print(f"  FHIR validation: PASS ({bundle_data['total']} resources)")

    # Chart title includes gene symbol or query for context
    chart_title = (
        f"Phase Distribution -- {gene_context['symbol']}"
        if gene_context
        else f"Phase Distribution -- {query_info['query'][:40]}"
    )
    chart = write_phase_chart(trials, args.output, title=chart_title)
    if chart:
        print(f"Chart   -> {chart}")

    # Reproducibility outputs -- always last so checksums cover everything
    cmds = write_commands(args, args.output)
    print(f"Repro   -> {cmds}")

    checksums = write_checksums(args.output)
    print(f"SHA-256 -> {checksums}")


if __name__ == "__main__":
    main()
