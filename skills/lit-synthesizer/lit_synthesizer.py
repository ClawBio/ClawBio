#!/usr/bin/env python3
"""
ClawBio — Lit Synthesizer Skill
Searches PubMed and bioRxiv for bioinformatics literature,
synthesizes results, and builds citation graphs.

Usage:
    python lit_synthesizer.py --query "CRISPR off-target effects" --output report/
    python lit_synthesizer.py --demo --output /tmp/demo
"""

import argparse
import csv
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path
from typing import Optional
import urllib.request
import urllib.parse
import xml.etree.ElementTree as ET


# ── Constants ──────────────────────────────────────────────────────────────────
PUBMED_SEARCH_URL = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
PUBMED_FETCH_URL  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
BIORXIV_API_URL   = "https://api.biorxiv.org/details/biorxiv"
TOOL_NAME         = "ClawBio-LitSynthesizer"
TOOL_EMAIL        = "mc@manuelcorpas.com"
MAX_RESULTS       = 10
DISCLAIMER = (
    "\n---\n*ClawBio Lit Synthesizer is a research tool. "
    "Always verify citations independently before use in publications.*\n"
)

# ── Demo data (no network needed) ─────────────────────────────────────────────
# Real PubMed and bioRxiv records, parsed by this module from the saved
# responses in examples/demo_sources/ (tests/test_demo_provenance.py checks
# every field against them). Refresh with
# examples/demo_sources/fetch_demo_sources.py.
DEMO_PAPERS_PATH = Path(__file__).resolve().parent / "examples" / "demo_papers.json"


def load_demo_papers() -> list[dict]:
    """Return the demo papers: real records, no network."""
    return json.loads(DEMO_PAPERS_PATH.read_text(encoding="utf-8"))


# ── PubMed helpers ─────────────────────────────────────────────────────────────

def search_pubmed(query: str, max_results: int = MAX_RESULTS) -> list[str]:
    """Search PubMed and return a list of PMIDs."""
    params = urllib.parse.urlencode({
        "db": "pubmed",
        "term": query,
        "retmax": max_results,
        "retmode": "json",
        "tool": TOOL_NAME,
        "email": TOOL_EMAIL,
    })
    url = f"{PUBMED_SEARCH_URL}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        return data.get("esearchresult", {}).get("idlist", [])
    except Exception as exc:
        print(f"  [warn] PubMed search failed: {exc}", file=sys.stderr)
        return []


def fetch_pubmed_details(pmids: list[str]) -> list[dict]:
    """Fetch title/abstract/authors for a list of PMIDs."""
    if not pmids:
        return []
    params = urllib.parse.urlencode({
        "db": "pubmed",
        "id": ",".join(pmids),
        "retmode": "xml",
        "rettype": "abstract",
        "tool": TOOL_NAME,
        "email": TOOL_EMAIL,
    })
    url = f"{PUBMED_FETCH_URL}?{params}"
    try:
        with urllib.request.urlopen(url, timeout=15) as resp:
            xml_data = resp.read()
        papers = parse_pubmed_xml(xml_data)
        time.sleep(0.34)   # NCBI rate limit: 3 req/s
    except Exception as exc:
        print(f"  [warn] PubMed fetch failed: {exc}", file=sys.stderr)
        papers = []
    return papers


def _all_text(el: Optional[ET.Element]) -> Optional[str]:
    """Element text including inline children such as <sup> and <i>.

    `.text` stops at the first child, which cut abstracts with reference
    superscripts to their first clause.
    """
    if el is None:
        return None
    return " ".join("".join(el.itertext()).split())


def parse_pubmed_xml(xml_data: bytes) -> list[dict]:
    """Parse an efetch PubMed XML response into paper records."""
    papers = []
    root = ET.fromstring(xml_data)
    for article in root.findall(".//PubmedArticle"):
        pmid_el   = article.find(".//PMID")
        title_el  = article.find(".//ArticleTitle")
        abstract_parts = [_all_text(el) for el in article.findall(".//Abstract/AbstractText")]
        journal_el  = article.find(".//Journal/Title")
        year_el     = article.find(".//PubDate/Year")
        doi_el      = article.find(".//ArticleId[@IdType='doi']")
        authors = [
            f"{a.findtext('LastName', '')} {a.findtext('Initials', '')}".strip()
            for a in article.findall(".//Author")
            if a.findtext("LastName")
        ]
        pmid = pmid_el.text if pmid_el is not None else "unknown"
        papers.append({
            "source":   "PubMed",
            "pmid":     pmid,
            "title":    _all_text(title_el) or "No title",
            "authors":  authors,
            "journal":  journal_el.text if journal_el is not None else "Unknown",
            "year":     year_el.text if year_el is not None else "Unknown",
            "abstract": " ".join(p for p in abstract_parts if p) or "No abstract available.",
            "doi":      doi_el.text if doi_el is not None else "",
            "url":      f"https://pubmed.ncbi.nlm.nih.gov/{pmid}",
            "citations": [],
        })
    return papers


# ── bioRxiv helpers ────────────────────────────────────────────────────────────

def biorxiv_item_to_paper(item: dict) -> dict:
    """Map one bioRxiv API `collection` item to a paper record."""
    return {
        "source":   "bioRxiv",
        "pmid":     None,
        "title":    item.get("title", "No title"),
        "authors":  item.get("authors", "").split("; "),
        "journal":  "bioRxiv (preprint)",
        "year":     item.get("date", "")[:4],
        "abstract": item.get("abstract", "No abstract available."),
        "doi":      item.get("doi", ""),
        "url":      f"https://biorxiv.org/content/{item.get('doi', '')}",
        "citations": [],
    }


def search_biorxiv(query: str, max_results: int = 5) -> list[dict]:
    """Search bioRxiv via their public API."""
    # bioRxiv API supports date-range + keyword in server/interval/cursor format.
    # We use the /details endpoint with a keyword search via the publisher's
    # simple query parameter (supported since 2023).
    start_date = "2023-01-01"
    end_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
    url = f"{BIORXIV_API_URL}/{start_date}/{end_date}/0/json"
    papers = []
    try:
        with urllib.request.urlopen(url, timeout=10) as resp:
            data = json.loads(resp.read().decode())
        collection = data.get("collection", [])
        # Filter by keyword in title/abstract
        keywords = query.lower().split()
        for item in collection:
            text = (item.get("title", "") + " " + item.get("abstract", "")).lower()
            if all(kw in text for kw in keywords):
                papers.append(biorxiv_item_to_paper(item))
            if len(papers) >= max_results:
                break
    except Exception as exc:
        print(f"  [warn] bioRxiv search failed: {exc}", file=sys.stderr)
    return papers


# ── Citation graph ─────────────────────────────────────────────────────────────

def build_citation_graph(papers: list[dict]) -> dict:
    """Build a simple citation graph: node per paper, edges from citations list."""
    nodes = []
    edges = []
    pmid_to_title = {p["pmid"]: p["title"] for p in papers if p["pmid"]}
    for paper in papers:
        node_id = paper["pmid"] or paper["doi"] or paper["title"][:40]
        nodes.append({
            "id":     node_id,
            "title":  paper["title"],
            "source": paper["source"],
            "year":   paper["year"],
        })
        for cited_pmid in paper.get("citations", []):
            if cited_pmid in pmid_to_title:
                edges.append({"from": node_id, "to": cited_pmid})
    return {"nodes": nodes, "edges": edges}


# ── Report generation ──────────────────────────────────────────────────────────

def generate_report(query: str, papers: list[dict], graph: dict, output_dir: Path) -> None:
    """Write the full markdown report, JSON results, CSV table, and reproducibility bundle."""
    output_dir.mkdir(parents=True, exist_ok=True)
    (output_dir / "tables").mkdir(exist_ok=True)
    (output_dir / "reproducibility").mkdir(exist_ok=True)

    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    pubmed_count  = sum(1 for p in papers if p["source"] == "PubMed")
    biorxiv_count = sum(1 for p in papers if p["source"] == "bioRxiv")

    # ── report.md ──
    lines = [
        f"# 🦖 ClawBio Lit Synthesizer Report",
        f"",
        f"**Query**: `{query}`  ",
        f"**Date**: {now}  ",
        f"**Sources**: PubMed ({pubmed_count} results) · bioRxiv ({biorxiv_count} results)  ",
        f"**Total papers**: {len(papers)}",
        f"",
        f"---",
        f"",
        f"## Summary",
        f"",
    ]

    if papers:
        # Simple keyword-frequency summary
        all_text = " ".join(p.get("abstract", "") for p in papers).lower()
        bio_terms = [
            "crispr", "genome editing", "off-target", "base editing", "prime editing",
            "cas9", "guide rna", "epigenome", "transcriptome", "variant", "mutation",
            "sequencing", "alignment", "annotation", "pipeline", "workflow",
        ]
        found_terms = [t for t in bio_terms if t in all_text]
        lines.append(
            f"Across {len(papers)} retrieved papers, recurring themes include: "
            + (", ".join(f"**{t}**" for t in found_terms[:6]) if found_terms else "various topics")
            + f". The literature spans {min((p['year'] for p in papers if p['year'].isdigit()), default='N/A')} "
            + f"to {max((p['year'] for p in papers if p['year'].isdigit()), default='N/A')}."
        )
    else:
        lines.append("No papers retrieved for this query.")

    lines += ["", "---", "", "## Papers", ""]

    for i, paper in enumerate(papers, 1):
        lines += [
            f"### {i}. {paper['title']}",
            f"",
            f"| Field | Value |",
            f"|-------|-------|",
            f"| **Source** | {paper['source']} |",
            f"| **Authors** | {', '.join(paper['authors'][:3])}{' et al.' if len(paper['authors']) > 3 else ''} |",
            f"| **Journal** | {paper['journal']} |",
            f"| **Year** | {paper['year']} |",
            f"| **DOI** | {paper['doi'] or 'N/A'} |",
            f"| **URL** | [{paper['url']}]({paper['url']}) |",
            f"",
            f"**Abstract**: {paper['abstract'][:400]}{'...' if len(paper['abstract']) > 400 else ''}",
            f"",
        ]

    lines += [
        "---",
        "",
        "## Citation Graph",
        "",
        f"**Nodes**: {len(graph['nodes'])}  ",
        f"**Edges** (internal citations): {len(graph['edges'])}",
        "",
        "```json",
        json.dumps(graph, indent=2)[:800] + ("\n... (see citation_graph.json for full graph)" if len(json.dumps(graph)) > 800 else ""),
        "```",
        "",
        "---",
        DISCLAIMER,
    ]

    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")

    # ── results.json ──
    results_path = output_dir / "results.json"
    results_path.write_text(json.dumps(papers, indent=2), encoding="utf-8")

    # ── citation_graph.json ──
    graph_path = output_dir / "citation_graph.json"
    graph_path.write_text(json.dumps(graph, indent=2), encoding="utf-8")

    # ── tables/papers.csv ──
    csv_path = output_dir / "tables" / "papers.csv"
    with open(csv_path, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(
            f, fieldnames=["source", "pmid", "title", "authors", "journal", "year", "doi", "url"]
        )
        writer.writeheader()
        for p in papers:
            writer.writerow({
                "source":  p["source"],
                "pmid":    p["pmid"] or "",
                "title":   p["title"],
                "authors": "; ".join(p["authors"]),
                "journal": p["journal"],
                "year":    p["year"],
                "doi":     p["doi"],
                "url":     p["url"],
            })


    # ── reproducibility bundle ──
    commands = f"""#!/usr/bin/env bash
# ClawBio Lit Synthesizer — portable reproducibility bundle
# Generated: {now}
# Query: {query}
#
# How to replay:
#   bash reproducibility/commands.sh
# from anywhere inside the repository clone.

set -euo pipefail

# ── Locate repo root ──────────────────────────────────────────────────────────
SCRIPT_DIR="$(cd "$(dirname "${{BASH_SOURCE[0]}}")" && pwd)"
REPO_ROOT="$SCRIPT_DIR"
while [[ ! -d "$REPO_ROOT/skills" && "$REPO_ROOT" != "/" ]]; do
  REPO_ROOT="$(dirname "$REPO_ROOT")"
done
if [[ ! -d "$REPO_ROOT/skills" ]]; then
  echo "ERROR: Could not locate repo root (no skills/ directory found)" >&2
  exit 1
fi

# ── Replay command ────────────────────────────────────────────────────────────
python "$REPO_ROOT/skills/lit-synthesizer/lit_synthesizer.py" \\
    --query "{query}" \\
    --output "./report"
"""
    (output_dir / "reproducibility" / "commands.sh").write_text(commands)

    env_yml = f"""name: clawbio-lit-synthesizer
channels:
  - defaults
dependencies:
  - python=3.11
  - pip
  - pip:
    - biopython>=1.83
# Generated: {now}
"""
    (output_dir / "reproducibility" / "environment.yml").write_text(env_yml)

    # ── SHA-256 checksums ──
    checksums = {}
    for fpath in output_dir.rglob("*"):
        if fpath.is_file() and "checksums" not in fpath.name:
            h = hashlib.sha256(fpath.read_bytes()).hexdigest()
            checksums[str(fpath.relative_to(output_dir))] = h
    checksum_lines = "\n".join(f"{h}  {f}" for f, h in checksums.items())
    (output_dir / "reproducibility" / "checksums.sha256").write_text(checksum_lines)

    print(f"\n✅ Report written to: {output_dir}/report.md")
    print(f"📄 {len(papers)} papers | 📊 Citation graph: {len(graph['nodes'])} nodes, {len(graph['edges'])} edges")
    print(f"🔒 Checksums: {output_dir}/reproducibility/checksums.sha256")


# ── Main ───────────────────────────────────────────────────────────────────────

def main() -> None:
    parser = argparse.ArgumentParser(
        description="ClawBio Lit Synthesizer — PubMed/bioRxiv search and synthesis"
    )
    parser.add_argument("--query",   type=str, help="Search query (e.g. 'CRISPR off-target effects')")
    parser.add_argument("--output",  type=str, default="lit_report", help="Output directory")
    parser.add_argument("--max",     type=int, default=MAX_RESULTS, help="Max results per source")
    parser.add_argument("--demo",    action="store_true", help="Run with built-in demo data (no network)")
    args = parser.parse_args()

    if not args.demo and not args.query:
        parser.error("Provide --query or use --demo")

    output_dir = Path(args.output)

    if args.demo:
        print("🦖 ClawBio Lit Synthesizer — DEMO MODE")
        print("   Using built-in demo data (no network calls)")
        query  = "CRISPR genome editing"
        papers = load_demo_papers()
    else:
        query = args.query
        print(f"🦖 ClawBio Lit Synthesizer")
        print(f"   Query : {query}")
        print(f"   Output: {output_dir}")
        print()

        print("🔍 Searching PubMed …")
        pmids  = search_pubmed(query, args.max)
        papers = fetch_pubmed_details(pmids)
        print(f"   Found {len(papers)} PubMed results")

        print("🔍 Searching bioRxiv …")
        biorxiv_papers = search_biorxiv(query, max_results=5)
        papers += biorxiv_papers
        print(f"   Found {len(biorxiv_papers)} bioRxiv results")

    print("\n📊 Building citation graph …")
    graph = build_citation_graph(papers)

    print("📝 Generating report …")
    generate_report(query, papers, graph, output_dir)


if __name__ == "__main__":
    main()
