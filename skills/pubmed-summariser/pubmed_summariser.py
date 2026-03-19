#!/usr/bin/env python3
"""
PubMed Summariser — ClawBio skill for literature retrieval.

Queries PubMed for a gene name or disease term and produces:
  - A terminal summary of the top recent papers
  - An HTML report saved to --output/report.html

Usage:
    python pubmed_summariser.py --query BRCA1 --output /tmp/pubmed_demo
    python pubmed_summariser.py --query "type 2 diabetes" --output /tmp/demo
    python pubmed_summariser.py --demo --output /tmp/pubmed_demo
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Project root on sys.path (required to import clawbio.common)
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import pubmed_api
from clawbio.common.html_report import HtmlReportBuilder, write_html_report

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEMO_QUERY = "BRCA1"
MAX_RESULTS_HARD_CAP = 50
SKILL_VERSION = "0.1.0"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clamp_max_results(n: int) -> int:
    """Clamp max_results to 50; print a warning if clamped."""
    if n > MAX_RESULTS_HARD_CAP:
        print(f"[warning] --max-results capped at {MAX_RESULTS_HARD_CAP} (you asked for {n})")
        return MAX_RESULTS_HARD_CAP
    return n


def _format_terminal_summary(query: str, papers: list[dict]) -> str:
    """Build the terminal summary string."""
    lines = []
    header = f"PubMed Research Briefing: {query}"
    lines.append(header)
    lines.append("=" * len(header))

    if not papers:
        lines.append(f"No results found for query: {query}")
        return "\n".join(lines)

    lines.append(f"Found {len(papers)} papers (sorted by date, English only)\n")
    for i, paper in enumerate(papers, 1):
        lines.append(f"{i}. {paper['title']}")
        lines.append(f"   Authors: {paper['authors']}")
        lines.append(f"   Journal: {paper['journal']} | {paper['date']}")
        if paper["abstract"]:
            lines.append(f"   Abstract: {paper['abstract']}")
        lines.append(f"   URL: {paper['url']}")
        lines.append("")

    return "\n".join(lines)


def _build_html_report(query: str, papers: list[dict]) -> str:
    """Build and return the HTML report string using HtmlReportBuilder."""
    builder = HtmlReportBuilder(
        title=f"PubMed Research Briefing: {query}",
        skill="pubmed-summariser",
    )

    builder.add_header_block(
        title=f"PubMed Research Briefing: {query}",
        subtitle=f"{len(papers)} recent English-language papers",
    )
    builder.add_metadata({
        "Query": query,
        "Results": str(len(papers)),
        "Sorted by": "Date (newest first)",
        "Language filter": "English only",
    })
    builder.add_disclaimer()

    builder.add_section("Papers", level=2)
    for i, paper in enumerate(papers, 1):
        card_html = (
            f"<div style='border:1px solid #e0e0e0;border-radius:8px;"
            f"padding:16px;margin:12px 0;background:#fff;'>"
            f"<h3 style='margin:0 0 8px 0;font-size:1em;'>"
            f"<a href='{paper['url']}' target='_blank'>{i}. {paper['title']}</a></h3>"
            f"<p style='margin:4px 0;color:#616161;font-size:0.9em;'>"
            f"<strong>Authors:</strong> {paper['authors']}</p>"
            f"<p style='margin:4px 0;color:#616161;font-size:0.9em;'>"
            f"<strong>Journal:</strong> {paper['journal']} &middot; {paper['date']}</p>"
        )
        if paper["abstract"]:
            card_html += (
                f"<p style='margin:8px 0 0 0;font-size:0.9em;'>"
                f"{paper['abstract']}</p>"
            )
        card_html += "</div>"
        builder.add_raw_html(card_html)

    builder.add_footer_block(skill="pubmed-summariser", version=SKILL_VERSION)
    return builder.render()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="PubMed Summariser — fetch and summarise recent PubMed papers",
    )
    parser.add_argument("--query", help="Gene name or disease term (e.g. BRCA1, type 2 diabetes)")
    parser.add_argument("--output", required=True, help="Directory to save report.html")
    parser.add_argument("--max-results", type=int, default=10, help="Number of results (default 10, max 50)")
    parser.add_argument("--demo", action="store_true", help="Run demo with BRCA1")
    args = parser.parse_args(argv)

    # Resolve query
    if args.demo:
        if args.query:
            print(f"[demo mode] ignoring --query, using {DEMO_QUERY}")
        query = DEMO_QUERY
    elif args.query:
        query = args.query
    else:
        parser.error("--query is required unless --demo is set")

    max_results = _clamp_max_results(args.max_results)

    # Fetch papers
    try:
        papers = pubmed_api.fetch_papers(query, max_results)
    except Exception as exc:
        print(f"Error fetching papers from PubMed: {exc}", file=sys.stderr)
        sys.exit(1)

    # Terminal output
    summary = _format_terminal_summary(query, papers)
    print(summary)

    # HTML report
    if papers:
        html = _build_html_report(query, papers)
        report_path = write_html_report(args.output, "report.html", html)
        print(f"Report saved to: {report_path}")


if __name__ == "__main__":
    main()
