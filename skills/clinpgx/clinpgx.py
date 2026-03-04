#!/usr/bin/env python3
"""
ClinPGx: Query the ClinPGx REST API for pharmacogenomic data.

Usage:
    python clinpgx.py --gene CYP2D6 --output report/
    python clinpgx.py --genes "CYP2D6,CYP2C19" --drugs "warfarin" --output report/
    python clinpgx.py --demo --output /tmp/clinpgx_demo
"""

from __future__ import annotations

import argparse
import csv
import hashlib
import json
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional

import requests

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------

BASE_URL = "https://api.clinpgx.org/v1"
RATE_LIMIT_INTERVAL = 0.55  # seconds between requests (stay under 2 req/sec)
CACHE_TTL = 86400  # 24 hours
USER_AGENT = "ClawBio-ClinPGx/0.1.0"
DEMO_GENE = "CYP2D6"

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

# ---------------------------------------------------------------------------
# API Client
# ---------------------------------------------------------------------------


class ClinPGxClient:
    """Rate-limited, caching HTTP client for the ClinPGx REST API."""

    def __init__(self, cache_dir: Path, use_cache: bool = True):
        self.cache_dir = cache_dir
        self.use_cache = use_cache
        self._last_request_time = 0.0
        self.session = requests.Session()
        self.session.headers.update({
            "Accept": "application/json",
            "User-Agent": USER_AGENT,
        })
        if use_cache:
            cache_dir.mkdir(parents=True, exist_ok=True)

    # --- Rate limiting ---

    def _throttle(self):
        elapsed = time.time() - self._last_request_time
        if elapsed < RATE_LIMIT_INTERVAL:
            time.sleep(RATE_LIMIT_INTERVAL - elapsed)
        self._last_request_time = time.time()

    # --- Caching ---

    def _cache_key(self, endpoint: str, params: dict) -> str:
        raw = f"{endpoint}|{json.dumps(params, sort_keys=True)}"
        return hashlib.sha256(raw.encode()).hexdigest()[:16]

    def _get_cached(self, key: str) -> Optional[dict]:
        path = self.cache_dir / f"{key}.json"
        if not path.exists():
            return None
        try:
            data = json.loads(path.read_text())
            if time.time() - data.get("_cached_at", 0) < CACHE_TTL:
                return data.get("response")
        except (json.JSONDecodeError, KeyError):
            pass
        return None

    def _set_cached(self, key: str, response_data):
        path = self.cache_dir / f"{key}.json"
        path.write_text(json.dumps({
            "_cached_at": time.time(),
            "response": response_data,
        }, indent=2))

    # --- Core request ---

    def _request(self, endpoint: str, params: dict | None = None) -> dict:
        params = params or {}
        cache_key = self._cache_key(endpoint, params)

        if self.use_cache:
            cached = self._get_cached(cache_key)
            if cached is not None:
                return cached

        self._throttle()
        url = f"{BASE_URL}/{endpoint}"
        resp = self.session.get(url, params=params, timeout=30)

        if resp.status_code == 429:
            time.sleep(2.0)
            resp = self.session.get(url, params=params, timeout=30)

        resp.raise_for_status()
        data = resp.json()

        if self.use_cache:
            self._set_cached(cache_key, data)

        return data

    # --- Public API methods ---

    def get_gene(self, symbol: str) -> list[dict]:
        """Search gene by HGNC symbol."""
        result = self._request("data/gene", {"symbol": symbol, "view": "max"})
        return result.get("data", [])

    def search_chemical(self, name: str) -> list[dict]:
        """Search chemical/drug by name."""
        result = self._request("data/chemical", {"name": name})
        return result.get("data", [])

    def get_clinical_annotations(
        self, gene_symbol: str | None = None, chemical_name: str | None = None
    ) -> list[dict]:
        """Get clinical annotations filtered by gene and/or drug."""
        params: dict = {}
        if gene_symbol:
            params["location.genes.symbol"] = gene_symbol
        if chemical_name:
            params["relatedChemicals.name"] = chemical_name
        result = self._request("data/clinicalAnnotation", params)
        return result.get("data", [])

    def get_guidelines(
        self, gene_accession_id: str | None = None, source: str | None = None
    ) -> list[dict]:
        """Get guideline annotations."""
        params: dict = {}
        if gene_accession_id:
            params["relatedGenes.accessionId"] = gene_accession_id
        if source:
            params["source"] = source
        result = self._request("data/guidelineAnnotation", params)
        return result.get("data", [])

    def get_drug_labels(
        self, gene_symbol: str | None = None, chemical_name: str | None = None
    ) -> list[dict]:
        """Get FDA/EMA drug labels with PGx info."""
        params: dict = {}
        if gene_symbol:
            params["relatedGenes.symbol"] = gene_symbol
        if chemical_name:
            params["relatedChemicals.name"] = chemical_name
        result = self._request("data/label", params)
        return result.get("data", [])

    def get_variant_annotations(self, gene_symbol: str) -> list[dict]:
        """Get variant annotations for a gene."""
        result = self._request(
            "data/variantAnnotation",
            {"location.genes.symbol": gene_symbol},
        )
        return result.get("data", [])


# ---------------------------------------------------------------------------
# Data extraction helpers
# ---------------------------------------------------------------------------


def extract_gene_summary(gene_data: dict) -> dict:
    """Extract key fields from a gene API response."""
    return {
        "symbol": gene_data.get("symbol", ""),
        "name": gene_data.get("name", ""),
        "id": gene_data.get("id", ""),
        "chr": gene_data.get("chr", {}).get("name", ""),
        "cpic_gene": gene_data.get("cpicGene", False),
        "allele_type": gene_data.get("alleleType", ""),
    }


def extract_annotation_row(ann: dict) -> dict:
    """Extract a row from a clinical annotation."""
    chemicals = ann.get("relatedChemicals", [])
    chem_names = ", ".join(c.get("name", "") for c in chemicals)
    genes = ann.get("relatedGenes", [])
    gene_names = ", ".join(g.get("symbol", "") for g in genes)
    level = ann.get("levelOfEvidence", {})
    level_term = level.get("term", "") if isinstance(level, dict) else str(level)
    return {
        "id": ann.get("accessionId", ann.get("id", "")),
        "gene": gene_names,
        "drug": chem_names,
        "evidence_level": level_term,
        "phenotype_category": ann.get("phenotypeCategory", ""),
    }


def extract_guideline_row(gl: dict) -> dict:
    """Extract a row from a guideline annotation."""
    chemicals = gl.get("relatedChemicals", [])
    chem_names = ", ".join(c.get("name", "") for c in chemicals)
    genes = gl.get("relatedGenes", [])
    gene_names = ", ".join(g.get("symbol", "") for g in genes)
    return {
        "id": gl.get("id", ""),
        "name": gl.get("name", ""),
        "gene": gene_names,
        "drug": chem_names,
        "source": gl.get("source", ""),
        "dosing_info": gl.get("dosingInformation", False),
    }


def extract_label_row(label: dict) -> dict:
    """Extract a row from a drug label."""
    chemicals = label.get("relatedChemicals", [])
    chem_names = ", ".join(c.get("name", "") for c in chemicals)
    genes = label.get("relatedGenes", [])
    gene_names = ", ".join(g.get("symbol", "") for g in genes)
    return {
        "id": label.get("id", ""),
        "name": label.get("name", ""),
        "gene": gene_names,
        "drug": chem_names,
        "source": label.get("source", ""),
        "testing_level": label.get("testingLevel", ""),
    }


# ---------------------------------------------------------------------------
# Query execution
# ---------------------------------------------------------------------------


def query_gene(client: ClinPGxClient, symbol: str) -> dict:
    """Run all queries for a gene and return structured results."""
    print(f"  Querying gene: {symbol}")
    genes = client.get_gene(symbol)
    if not genes:
        print(f"  WARNING: Gene '{symbol}' not found in ClinPGx")
        return {"symbol": symbol, "found": False}

    gene = genes[0]
    gene_summary = extract_gene_summary(gene)
    gene_id = gene.get("id", "")

    print(f"  Fetching clinical annotations...")
    annotations = client.get_clinical_annotations(gene_symbol=symbol)

    print(f"  Fetching guidelines...")
    guidelines = client.get_guidelines(gene_accession_id=gene_id)

    print(f"  Fetching drug labels...")
    labels = client.get_drug_labels(gene_symbol=symbol)

    return {
        "symbol": symbol,
        "found": True,
        "gene": gene_summary,
        "clinical_annotations": [extract_annotation_row(a) for a in annotations],
        "guidelines": [extract_guideline_row(g) for g in guidelines],
        "drug_labels": [extract_label_row(l) for l in labels],
    }


def query_drug(client: ClinPGxClient, name: str) -> dict:
    """Run all queries for a drug and return structured results."""
    print(f"  Querying drug: {name}")
    chemicals = client.search_chemical(name)
    if not chemicals:
        print(f"  WARNING: Drug '{name}' not found in ClinPGx")
        return {"name": name, "found": False}

    chem = chemicals[0]

    print(f"  Fetching clinical annotations...")
    annotations = client.get_clinical_annotations(chemical_name=name)

    print(f"  Fetching drug labels...")
    labels = client.get_drug_labels(chemical_name=name)

    return {
        "name": name,
        "found": True,
        "chemical": {
            "id": chem.get("id", ""),
            "name": chem.get("name", ""),
            "types": chem.get("types", []),
        },
        "clinical_annotations": [extract_annotation_row(a) for a in annotations],
        "drug_labels": [extract_label_row(l) for l in labels],
    }


# ---------------------------------------------------------------------------
# Report generation
# ---------------------------------------------------------------------------


def generate_text_summary(gene_results: list[dict], drug_results: list[dict]) -> str:
    """Generate a concise text summary for stdout."""
    lines = ["CLINPGX REPORT", "=" * 60, ""]

    for gr in gene_results:
        if not gr.get("found"):
            lines.append(f"Gene {gr['symbol']}: NOT FOUND")
            lines.append("")
            continue
        g = gr["gene"]
        lines.append(f"== GENE: {g['symbol']} ==")
        lines.append(f"Name: {g['name']}")
        lines.append(f"ClinPGx ID: {g['id']}")
        lines.append(f"Chromosome: {g['chr']}")
        lines.append(f"CPIC Gene: {'Yes' if g['cpic_gene'] else 'No'}")
        lines.append("")

        anns = gr["clinical_annotations"]
        lines.append(f"Clinical Annotations: {len(anns)} found")
        if anns:
            lines.append(f"  {'Drug':<30s} {'Evidence':<12s} {'Category'}")
            lines.append(f"  {'-'*30} {'-'*12} {'-'*20}")
            seen = set()
            for a in anns[:20]:
                key = (a["drug"], a["evidence_level"])
                if key in seen:
                    continue
                seen.add(key)
                lines.append(
                    f"  {a['drug']:<30s} {a['evidence_level']:<12s} {a['phenotype_category']}"
                )
            if len(anns) > 20:
                lines.append(f"  ... and {len(anns) - 20} more")
        lines.append("")

        gls = gr["guidelines"]
        lines.append(f"CPIC/DPWG Guidelines: {len(gls)} found")
        for gl in gls[:10]:
            lines.append(f"  - {gl['name']} (source: {gl['source']})")
        if len(gls) > 10:
            lines.append(f"  ... and {len(gls) - 10} more")
        lines.append("")

        labels = gr["drug_labels"]
        lines.append(f"Drug Labels: {len(labels)} found")
        for lb in labels[:10]:
            lines.append(f"  - {lb['name']} (source: {lb['source']})")
        if len(labels) > 10:
            lines.append(f"  ... and {len(labels) - 10} more")
        lines.append("")

    for dr in drug_results:
        if not dr.get("found"):
            lines.append(f"Drug {dr['name']}: NOT FOUND")
            lines.append("")
            continue
        c = dr["chemical"]
        lines.append(f"== DRUG: {c['name']} ==")
        lines.append(f"ClinPGx ID: {c['id']}")
        lines.append(f"Types: {', '.join(c['types'])}")
        lines.append("")

        anns = dr["clinical_annotations"]
        lines.append(f"Clinical Annotations: {len(anns)} found")
        if anns:
            lines.append(f"  {'Gene':<15s} {'Evidence':<12s} {'Category'}")
            lines.append(f"  {'-'*15} {'-'*12} {'-'*20}")
            seen = set()
            for a in anns[:20]:
                key = (a["gene"], a["evidence_level"])
                if key in seen:
                    continue
                seen.add(key)
                lines.append(
                    f"  {a['gene']:<15s} {a['evidence_level']:<12s} {a['phenotype_category']}"
                )
        lines.append("")

        labels = dr["drug_labels"]
        lines.append(f"Drug Labels: {len(labels)} found")
        for lb in labels[:10]:
            lines.append(f"  - {lb['name']} (source: {lb['source']})")
        lines.append("")

    lines.append("-" * 60)
    lines.append(f"Source: ClinPGx API ({BASE_URL})")
    lines.append("License: CC BY-SA 4.0")
    lines.append("")
    lines.append(DISCLAIMER)

    return "\n".join(lines)


def generate_markdown_report(
    gene_results: list[dict], drug_results: list[dict], query_desc: str
) -> str:
    """Generate a full markdown report."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines = [
        "# ClinPGx Report",
        "",
        f"**Date**: {now}",
        f"**Query**: {query_desc}",
        f"**API**: ClinPGx REST API v1 ({BASE_URL})",
        "**License**: CC BY-SA 4.0",
        "",
        "---",
        "",
    ]

    for gr in gene_results:
        if not gr.get("found"):
            lines.append(f"## Gene: {gr['symbol']} (not found)")
            lines.append("")
            continue

        g = gr["gene"]
        lines.append(f"## Gene: {g['symbol']}")
        lines.append("")
        lines.append("| Field | Value |")
        lines.append("|-------|-------|")
        lines.append(f"| **Symbol** | {g['symbol']} |")
        lines.append(f"| **Name** | {g['name']} |")
        lines.append(f"| **ClinPGx ID** | {g['id']} |")
        lines.append(f"| **Chromosome** | {g['chr']} |")
        lines.append(f"| **CPIC Gene** | {'Yes' if g['cpic_gene'] else 'No'} |")
        lines.append(f"| **Allele Type** | {g['allele_type']} |")
        lines.append("")

        anns = gr["clinical_annotations"]
        if anns:
            lines.append(f"### Clinical Annotations ({len(anns)})")
            lines.append("")
            lines.append("| ID | Drug | Evidence | Category |")
            lines.append("|----|------|----------|----------|")
            for a in anns:
                lines.append(
                    f"| {a['id']} | {a['drug']} | {a['evidence_level']} | {a['phenotype_category']} |"
                )
            lines.append("")

        gls = gr["guidelines"]
        if gls:
            lines.append(f"### Guidelines ({len(gls)})")
            lines.append("")
            lines.append("| Guideline | Source | Dosing Info |")
            lines.append("|-----------|--------|-------------|")
            for gl in gls:
                dosing = "Yes" if gl["dosing_info"] else "No"
                lines.append(f"| {gl['name']} | {gl['source']} | {dosing} |")
            lines.append("")

        labels = gr["drug_labels"]
        if labels:
            lines.append(f"### Drug Labels ({len(labels)})")
            lines.append("")
            lines.append("| Label | Drug | Source | Testing Level |")
            lines.append("|-------|------|--------|---------------|")
            for lb in labels:
                lines.append(
                    f"| {lb['name']} | {lb['drug']} | {lb['source']} | {lb['testing_level']} |"
                )
            lines.append("")

    for dr in drug_results:
        if not dr.get("found"):
            lines.append(f"## Drug: {dr['name']} (not found)")
            lines.append("")
            continue

        c = dr["chemical"]
        lines.append(f"## Drug: {c['name']}")
        lines.append("")
        lines.append("| Field | Value |")
        lines.append("|-------|-------|")
        lines.append(f"| **Name** | {c['name']} |")
        lines.append(f"| **ClinPGx ID** | {c['id']} |")
        lines.append(f"| **Types** | {', '.join(c['types'])} |")
        lines.append("")

        anns = dr["clinical_annotations"]
        if anns:
            lines.append(f"### Clinical Annotations ({len(anns)})")
            lines.append("")
            lines.append("| ID | Gene | Evidence | Category |")
            lines.append("|----|------|----------|----------|")
            for a in anns:
                lines.append(
                    f"| {a['id']} | {a['gene']} | {a['evidence_level']} | {a['phenotype_category']} |"
                )
            lines.append("")

        labels = dr["drug_labels"]
        if labels:
            lines.append(f"### Drug Labels ({len(labels)})")
            lines.append("")
            lines.append("| Label | Gene | Source | Testing Level |")
            lines.append("|-------|------|--------|---------------|")
            for lb in labels:
                lines.append(
                    f"| {lb['name']} | {lb['gene']} | {lb['source']} | {lb['testing_level']} |"
                )
            lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(f"- **Data source**: ClinPGx REST API v1 ({BASE_URL})")
    lines.append("- **Rate limit**: 2 requests/second observed")
    lines.append("- **Cache**: 24-hour local file cache")
    lines.append("")
    lines.append("## Attribution")
    lines.append("")
    lines.append("Data from ClinPGx (PharmGKB + CPIC + PharmCAT), licensed under CC BY-SA 4.0.")
    lines.append("")
    lines.append("## Disclaimer")
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")

    return "\n".join(lines)


# ---------------------------------------------------------------------------
# CSV output
# ---------------------------------------------------------------------------


def write_csv(filepath: Path, rows: list[dict]):
    """Write a list of dicts as CSV."""
    if not rows:
        return
    filepath.parent.mkdir(parents=True, exist_ok=True)
    with open(filepath, "w", newline="") as f:
        writer = csv.DictWriter(f, fieldnames=rows[0].keys())
        writer.writeheader()
        writer.writerows(rows)


def write_tables(
    output_dir: Path, gene_results: list[dict], drug_results: list[dict]
):
    """Write CSV tables from results."""
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(parents=True, exist_ok=True)

    all_annotations = []
    all_guidelines = []
    all_labels = []

    for gr in gene_results:
        if gr.get("found"):
            all_annotations.extend(gr["clinical_annotations"])
            all_guidelines.extend(gr["guidelines"])
            all_labels.extend(gr["drug_labels"])

    for dr in drug_results:
        if dr.get("found"):
            all_annotations.extend(dr["clinical_annotations"])
            all_labels.extend(dr["drug_labels"])

    if all_annotations:
        write_csv(tables_dir / "clinical_annotations.csv", all_annotations)
        print(f"  Wrote {len(all_annotations)} clinical annotations")

    if all_guidelines:
        write_csv(tables_dir / "guidelines.csv", all_guidelines)
        print(f"  Wrote {len(all_guidelines)} guidelines")

    if all_labels:
        write_csv(tables_dir / "drug_labels.csv", all_labels)
        print(f"  Wrote {len(all_labels)} drug labels")


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="ClinPGx: query the ClinPGx API for pharmacogenomic data"
    )
    parser.add_argument("--gene", help="Single gene symbol (e.g., CYP2D6)")
    parser.add_argument("--drug", help="Single drug name (e.g., warfarin)")
    parser.add_argument(
        "--genes", help="Comma-separated gene symbols (e.g., CYP2D6,CYP2C19)"
    )
    parser.add_argument(
        "--drugs", help="Comma-separated drug names (e.g., warfarin,clopidogrel)"
    )
    parser.add_argument("--output", "-o", help="Output directory for full report")
    parser.add_argument("--demo", action="store_true", help="Run demo with CYP2D6")
    parser.add_argument(
        "--no-cache", action="store_true", help="Skip local cache, always query API"
    )
    parser.add_argument(
        "--cache-dir",
        default=str(Path.home() / ".clawbio" / "clinpgx_cache"),
        help="Cache directory (default: ~/.clawbio/clinpgx_cache/)",
    )
    args = parser.parse_args()

    # Validate
    if not any([args.gene, args.drug, args.genes, args.drugs, args.demo]):
        parser.print_help()
        print("\nError: provide at least one of --gene, --drug, --genes, --drugs, or --demo")
        sys.exit(1)

    # Build query lists
    gene_symbols = []
    drug_names = []

    if args.demo:
        gene_symbols.append(DEMO_GENE)
    if args.gene:
        gene_symbols.append(args.gene.strip().upper())
    if args.genes:
        gene_symbols.extend(g.strip().upper() for g in args.genes.split(",") if g.strip())
    if args.drug:
        drug_names.append(args.drug.strip().lower())
    if args.drugs:
        drug_names.extend(d.strip().lower() for d in args.drugs.split(",") if d.strip())

    # Deduplicate
    gene_symbols = list(dict.fromkeys(gene_symbols))
    drug_names = list(dict.fromkeys(drug_names))

    query_desc = []
    if gene_symbols:
        query_desc.append(f"Genes: {', '.join(gene_symbols)}")
    if drug_names:
        query_desc.append(f"Drugs: {', '.join(drug_names)}")
    query_str = "; ".join(query_desc)

    total_queries = len(gene_symbols) + len(drug_names)
    print(f"ClinPGx: querying {total_queries} item(s) from {BASE_URL}")
    print()

    # Init client
    cache_dir = Path(args.cache_dir)
    client = ClinPGxClient(cache_dir=cache_dir, use_cache=not args.no_cache)

    # Execute queries
    gene_results = []
    for symbol in gene_symbols:
        try:
            result = query_gene(client, symbol)
            gene_results.append(result)
        except requests.HTTPError as e:
            print(f"  ERROR querying gene {symbol}: {e}")
            gene_results.append({"symbol": symbol, "found": False})
        print()

    drug_results = []
    for name in drug_names:
        try:
            result = query_drug(client, name)
            drug_results.append(result)
        except requests.HTTPError as e:
            print(f"  ERROR querying drug {name}: {e}")
            drug_results.append({"name": name, "found": False})
        print()

    # Output
    if args.output:
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)

        # Write report
        report = generate_markdown_report(gene_results, drug_results, query_str)
        report_path = output_dir / "report.md"
        report_path.write_text(report)
        print(f"Report written to {report_path}")

        # Write CSV tables
        write_tables(output_dir, gene_results, drug_results)

        print(f"\nFull output in {output_dir}/")
    else:
        summary = generate_text_summary(gene_results, drug_results)
        print(summary)


if __name__ == "__main__":
    main()
