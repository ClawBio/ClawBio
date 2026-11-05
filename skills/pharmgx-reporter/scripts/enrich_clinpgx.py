#!/usr/bin/env python3
"""Step 3 of pharmgx_reporter: enrich drug recommendations with ClinPGx
evidence levels (evidence grade, CPIC-verified flag, guideline sources).

Ported (trimmed) from ClawBio 1.x `pharmgx_reporter.py`'s
`enrich_with_clinpgx()`. 1.x reached this data by importing a sibling
skill's Python module directly (`sys.path` hack into `skills/clinpgx/`) and
used the `requests` library. Neither carries over cleanly to 2.0:

  - Cross-skill imports don't fit the one-step-one-script DAG model (see
    `skills/drug_photo/SKILL.md`'s "Port status" section for the same
    argument made the other direction). This step is a standalone HTTP
    client instead.
  - `requests` isn't a dependency of any other skill visible from this
    worktree at port time; `urllib.request` (stdlib) is what
    `skills/gwas_lookup/scripts/fetch_api.py` already uses for the same
    "hit a public REST API, degrade gracefully offline" shape, so this
    script follows that precedent instead of adding a new dependency.
  - `extract_phenotype_recs()` (HTML guideline-table parsing to produce
    per-drug prose recommendations via `clawbio.common.rec_shortener`) is
    NOT ported — that package doesn't exist in the 2.0 repo and the
    parsing is a presentational nicety on top of the evidence_level/
    verified/sources fields this script does keep. See pgx_common.py's
    module docstring for the full list of intentionally-omitted 1.x pieces.

This step is declared as a normal (non-skippable) DAG step — this
worktree's `clawbio` CLI does not implement `optional: true` steps /
`clawbio skip` (checked at port time: no `skip` subcommand, `Step` has no
`optional` field). Making network access a hard requirement would strand
any offline run, so instead every per-gene API call is wrapped exactly
like `fetch_api.py`'s pattern: on failure, write a `status: "error"` stub
for that gene and keep going. The step always exits 0 and always emits a
valid (possibly empty/degraded) `clinpgx_enrichment.json`, so `clawbio
done` never blocks on network availability.

Usage:
    python enrich_clinpgx.py --drug-recommendations /tmp/out/drug_recommendations.json --out /tmp/out/clinpgx_enrichment.json
"""
from __future__ import annotations

import argparse
import json
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path
from typing import Any

BASE_URL = "https://api.clinpgx.org/v1"
RATE_LIMIT_INTERVAL = 0.55  # stay under ~2 req/sec
TIMEOUT = 15.0
USER_AGENT = "clawbio2-pharmgx-reporter"

_EVIDENCE_RANK = {"1A": 6, "1B": 5, "2A": 4, "2B": 3, "3": 2, "4": 1}

_last_request_time = 0.0


def _throttle() -> None:
    global _last_request_time
    elapsed = time.time() - _last_request_time
    if elapsed < RATE_LIMIT_INTERVAL:
        time.sleep(RATE_LIMIT_INTERVAL - elapsed)
    _last_request_time = time.time()


def _http_get_json(endpoint: str, params: dict) -> dict:
    _throttle()
    url = f"{BASE_URL}/{endpoint}?{urllib.parse.urlencode(params)}"
    req = urllib.request.Request(url, headers={"Accept": "application/json", "User-Agent": USER_AGENT})
    with urllib.request.urlopen(req, timeout=TIMEOUT) as r:
        return json.load(r)


def _genes_from_drug_recommendations(drug_results: dict) -> list[str]:
    genes: set[str] = set()
    for bucket in drug_results.values():
        for entry in bucket:
            gene_field = entry.get("gene", "")
            for g in gene_field.replace("+", ",").split(","):
                g = g.strip()
                if g:
                    genes.add(g)
    return sorted(genes)


def enrich_gene(gene: str) -> tuple[dict[str, Any], dict[str, Any]]:
    """Return (per_drug_enrichment, gene_status). Raises on network failure —
    caller wraps this per-gene so one bad gene doesn't sink the whole step."""
    enrichment: dict[str, Any] = {}

    gene_resp = _http_get_json("data/gene", {"symbol": gene, "view": "max"})
    gene_records = gene_resp.get("data", [])
    gene_id = gene_records[0].get("id", "") if gene_records else ""

    ann_resp = _http_get_json("data/clinicalAnnotation", {"location.genes.symbol": gene})
    annotations = ann_resp.get("data", [])

    guidelines = []
    if gene_id:
        gl_resp = _http_get_json("data/guidelineAnnotation", {"relatedGenes.accessionId": gene_id})
        guidelines = gl_resp.get("data", [])

    for ann in annotations:
        level_obj = ann.get("levelOfEvidence", {})
        level_term = level_obj.get("term", "") if isinstance(level_obj, dict) else str(level_obj)
        for chem in ann.get("relatedChemicals", []):
            chem_name = chem.get("name", "").lower()
            if not chem_name:
                continue
            existing = enrichment.get(chem_name, {})
            if _EVIDENCE_RANK.get(level_term, 0) > _EVIDENCE_RANK.get(existing.get("evidence_level", ""), 0):
                existing["evidence_level"] = level_term
            existing.setdefault("sources", set())
            enrichment[chem_name] = existing

    for gl in guidelines:
        source = gl.get("source", "")
        has_dosing = gl.get("dosingInformation", False)
        for chem in gl.get("relatedChemicals", []):
            chem_name = chem.get("name", "").lower()
            if not chem_name:
                continue
            existing = enrichment.get(chem_name, {})
            existing.setdefault("sources", set())
            if source:
                existing["sources"].add(source)
            if has_dosing and source == "CPIC":
                existing["verified"] = True
            enrichment[chem_name] = existing

    for entry in enrichment.values():
        entry["sources"] = sorted(entry.get("sources", set()))
        entry.setdefault("verified", False)
        entry.setdefault("evidence_level", "")

    return enrichment, {"gene": gene, "status": "ok", "annotations": len(annotations), "guidelines": len(guidelines)}


def main() -> int:
    parser = argparse.ArgumentParser(
        description="Enrich CPIC drug recommendations with ClinPGx evidence levels (best-effort, network-optional)"
    )
    parser.add_argument("--drug-recommendations", required=True, help="Path to drug_recommendations.json from lookup_drugs.py")
    parser.add_argument("--out", required=True, help="Output path for clinpgx_enrichment.json")
    args = parser.parse_args()

    payload = json.loads(Path(args.drug_recommendations).read_text())
    drug_results = payload["drug_recommendations"]
    genes = _genes_from_drug_recommendations(drug_results)

    enrichment: dict[str, Any] = {}
    gene_statuses: list[dict[str, Any]] = []

    print(f"Querying ClinPGx for evidence data across {len(genes)} gene(s)...")
    for idx, gene in enumerate(genes, 1):
        print(f"  [{idx}/{len(genes)}] {gene}")
        try:
            gene_enrichment, status = enrich_gene(gene)
            for chem_name, entry in gene_enrichment.items():
                existing = enrichment.get(chem_name, {})
                if _EVIDENCE_RANK.get(entry.get("evidence_level", ""), 0) > \
                        _EVIDENCE_RANK.get(existing.get("evidence_level", ""), 0):
                    existing["evidence_level"] = entry["evidence_level"]
                existing["sources"] = sorted(set(existing.get("sources", [])) | set(entry.get("sources", [])))
                existing["verified"] = existing.get("verified", False) or entry.get("verified", False)
                enrichment[chem_name] = existing
            gene_statuses.append(status)
        except (urllib.error.URLError, urllib.error.HTTPError, TimeoutError, json.JSONDecodeError) as e:
            print(f"    WARNING: ClinPGx query failed for {gene}: {type(e).__name__}: {e}", file=sys.stderr)
            gene_statuses.append({"gene": gene, "status": "error", "message": f"{type(e).__name__}: {e}"})

    ok = sum(1 for s in gene_statuses if s["status"] == "ok")
    if ok == 0 and genes:
        print("  ClinPGx enrichment unreachable for all genes (offline or API down); "
              "report will fall back to CPIC classification summaries only.", file=sys.stderr)
    elif enrichment:
        print(f"  Enriched {len(enrichment)} drug(s) with evidence data.")

    out_path = Path(args.out)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps({
        "queried_genes": genes,
        "gene_status": gene_statuses,
        "enrichment": enrichment,
    }, indent=2))
    print(f"\nClinPGx enrichment written: {out_path}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
