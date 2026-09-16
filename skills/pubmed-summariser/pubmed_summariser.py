#!/usr/bin/env python3
"""
PubMed Summariser — ClawBio skill for literature retrieval.

Queries PubMed for a gene name or disease term and produces:
  - A terminal summary of the top recent papers
  - Complete abstracts and separate summaries in result.json
  - Markdown and HTML reports, plus a reproducibility bundle

Usage:
    python pubmed_summariser.py --query BRCA1 --output /tmp/pubmed_demo
    python pubmed_summariser.py --query "type 2 diabetes" --output /tmp/demo
    python pubmed_summariser.py --demo --output /tmp/pubmed_demo
"""

from __future__ import annotations

import argparse
import html
import importlib.util
import json
import math
import shlex
import sys
from datetime import datetime, timezone
from importlib.metadata import version
from pathlib import Path

# ---------------------------------------------------------------------------
# Project root on sys.path (required to import clawbio.common)
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

def _load_sibling(name):
    spec = importlib.util.spec_from_file_location(name, Path(__file__).resolve().parent / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


pubmed_api = _load_sibling("pubmed_api")
abstract_summary = _load_sibling("abstract_summary")
from clawbio.common.html_report import HtmlReportBuilder, write_html_report
from clawbio.common.report import DISCLAIMER, write_result_json
from clawbio.common.reproducibility import write_commands_sh, write_environment_yml, write_checksums
from clawbio.providers import create_provider, validate_model_params

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEMO_QUERY = "BRCA1"
MAX_RESULTS_HARD_CAP = 50
SKILL_VERSION = "0.2.0"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clamp_max_results(n: int) -> int:
    """Clamp max_results to 50; print a warning if clamped."""
    if n <= 0:
        raise ValueError("--max-results must be positive")
    if n > MAX_RESULTS_HARD_CAP:
        print(f"[warning] --max-results capped at {MAX_RESULTS_HARD_CAP} (you asked for {n})")
        return MAX_RESULTS_HARD_CAP
    return n


def _summary(paper):
    return paper.get("summary") or abstract_summary.summarize_abstract(paper)


def _summary_label(summary: dict) -> str:
    if summary["method"] == "llm":
        return f"AI summary — {summary['provider']} / {summary['model']}"
    if summary["method"] == "unavailable":
        return "Abstract unavailable"
    return "Abstract opening — no AI"


def _format_terminal_summary(query: str, papers: list[dict]) -> str:
    """Build the terminal summary string."""
    lines = []
    header = f"PubMed Research Briefing: {query}"
    lines.append(header)
    lines.append("=" * len(header))

    if not papers:
        lines.append(f"No results found for query: {query}")
        lines.extend(["", DISCLAIMER])
        return "\n".join(lines)

    lines.append(f"Found {len(papers)} papers (sorted by date, English only)\n")
    for i, paper in enumerate(papers, 1):
        lines.append(f"{i}. {paper['title']}")
        lines.append(f"   Authors: {paper['authors']}")
        lines.append(f"   Journal: {paper['journal']} | {paper['date']}")
        summary = _summary(paper)
        lines.append(f"   {_summary_label(summary)}: {summary['text'] or 'No abstract available.'}")
        if summary["fallback_reason"]:
            lines.append("   AI summary could not be generated; showing the abstract opening.")
        lines.append(f"   URL: {paper['url']}")
        lines.append("")

    lines.append(DISCLAIMER)
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
    if not papers:
        builder.add_raw_html("<p>No results found.</p>")
    for i, paper in enumerate(papers, 1):
        card_html = (
            f"<div style='border:1px solid #e0e0e0;border-radius:8px;"
            f"padding:16px;margin:12px 0;background:#fff;'>"
            f"<h3 style='margin:0 0 8px 0;font-size:1em;'>"
            f"<a href=\"{html.escape(paper['url'])}\" target='_blank'>"
            f"{i}. {html.escape(paper['title'])}</a></h3>"
            f"<p style='margin:4px 0;color:#616161;font-size:0.9em;'>"
            f"<strong>Authors:</strong> {html.escape(paper['authors'])}</p>"
            f"<p style='margin:4px 0;color:#616161;font-size:0.9em;'>"
            f"<strong>Journal:</strong> {html.escape(paper['journal'])} &middot; {html.escape(paper['date'])}</p>"
        )
        summary = _summary(paper)
        card_html += (
            f"<p><strong>{html.escape(_summary_label(summary))}:</strong> "
            f"{html.escape(summary['text']) or 'No abstract available.'}</p>"
        )
        if summary["fallback_reason"]:
            card_html += "<p>AI summary could not be generated; showing the abstract opening.</p>"
        if paper["abstract"]:
            card_html += (
                "<details><summary>Full abstract</summary>"
            )
            sections = paper.get("abstract_sections") or [{"text": paper["abstract"], "label": None}]
            for section in sections:
                label = f"<strong>{html.escape(section['label'])}</strong> " if section["label"] else ""
                card_html += f"<p>{label}{html.escape(section['text'])}</p>"
            card_html += "</details>"
        card_html += "</div>"
        builder.add_raw_html(card_html)

    builder.add_footer_block(skill="pubmed-summariser", version=SKILL_VERSION)
    return builder.render()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def _write_outputs(output_dir: Path, query: str, papers: list[dict], config: dict,
                   retrieved_at: str, warnings: list[dict], max_results: int) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    report_html = _build_html_report(query, papers)
    report_path = write_html_report(output_dir, "report.html", report_html)
    lines = [f"# PubMed Research Briefing: {html.escape(query)}", "",
             f"Summary method requested: {config['method']}", "",
             "English-language papers, sorted by date.", ""]
    if not papers:
        lines.extend(["No results found.", ""])
    for paper in papers:
        summary = paper["summary"]
        lines.extend([
            f"## {html.escape(paper['title'])}", "",
            f"Authors: {html.escape(paper['authors'])}", "",
            f"Journal: {html.escape(paper['journal'])} | {paper['date']}", "",
            f"{html.escape(_summary_label(summary))}: {html.escape(summary['text']) or 'No abstract available.'}", "",
            f"[PubMed]({paper['url']})", "",
        ])
        if summary["fallback_reason"]:
            lines.extend(["AI summary could not be generated; showing the abstract opening.", ""])
    lines.extend(["---", "", DISCLAIMER, ""])
    markdown_path = output_dir / "report.md"
    markdown_path.write_text("\n".join(lines), encoding="utf-8")
    result_path = write_result_json(
        output_dir, "pubmed-summariser", SKILL_VERSION,
        summary={"paper_count": len(papers), "fallback_count": len(warnings)},
        data={"query": query, "retrieved_at": retrieved_at,
              "search": {"term": f"{query} AND english[la]", "sort": "date", "max_results": max_results},
              "summary_config": config, "papers": papers, "warnings": warnings},
        datasets={"literature": "PubMed"},
    )
    # Record effective settings rather than relying on model/endpoint env vars
    # or a query file that may no longer be present when replayed.
    command = ["python", Path(__file__).resolve().as_posix(),
               "--query", query, "--output", output_dir.as_posix(),
               "--max-results", str(max_results), "--summary-method", config["method"]]
    if config["method"] == "llm":
        command.extend(["--provider", config["provider"], "--model", config["model"],
                        "--base-url", config["base_url"], "--llm-timeout", str(config["timeout"]),
                        "--model-params", json.dumps(config["model_params"])])
        if config["max_output_tokens"] is not None:
            command.extend(["--summary-max-tokens", str(config["max_output_tokens"])])
    command_path = write_commands_sh(output_dir, shlex.join(command))
    deps = [f"requests=={version('requests')}"]
    if config["method"] == "llm":
        deps.append(f"openai=={version('openai')}")
    environment_path = write_environment_yml(
        output_dir, "clawbio-pubmed", deps,
        python_version=f"{sys.version_info.major}.{sys.version_info.minor}",
    )
    write_checksums([result_path, markdown_path, report_path, command_path, environment_path], output_dir, anchor=output_dir)
    print(f"Reports saved to: {output_dir}")


def main(argv: list[str] | None = None) -> None:
    # The runner captures stdout with the Windows locale encoding. Keep that
    # encoding, but escape unsupported characters instead of losing the run.
    # All saved artifacts retain the original Unicode text.
    for stream in (sys.stdout, sys.stderr):
        if hasattr(stream, "reconfigure"):
            stream.reconfigure(errors="backslashreplace")
    argv = list(sys.argv[1:] if argv is None else argv)
    parser = argparse.ArgumentParser(
        description="PubMed Summariser — fetch and summarise recent PubMed papers",
    )
    query_source = parser.add_mutually_exclusive_group()
    query_source.add_argument("--query", help="Gene name or disease term (e.g. BRCA1, type 2 diabetes)")
    query_source.add_argument("--input", type=Path, help="UTF-8 text file containing one public literature search query")
    parser.add_argument("--output", required=True, help="Directory for JSON, Markdown, HTML and reproducibility files")
    parser.add_argument("--max-results", type=int, default=10, help="Number of results (default 10, max 50)")
    parser.add_argument("--demo", action="store_true", help="Run demo with BRCA1")
    parser.add_argument("--summary-method", choices=["first-sentence", "llm"], default="first-sentence",
                        help="Default: first-sentence. Copies the abstract's opening, up to 300 characters; "
                             "no AI model is used. Choose llm to summarize the complete abstract using a model.")
    parser.add_argument("--provider", choices=["openai", "ollama"], help="Required with --summary-method llm")
    parser.add_argument("--model", help="Model name; otherwise use provider model environment setting")
    parser.add_argument("--base-url", help="Provider endpoint including /v1")
    parser.add_argument("--llm-timeout", type=float, help="Seconds per request (default: 120)")
    parser.add_argument("--summary-max-tokens", type=int, help="Output token budget; default uses the model's limit")
    parser.add_argument("--model-params", help='JSON generation settings, e.g. {"reasoning_effort":"none"}')
    args = parser.parse_args(argv)

    # Resolve query
    if args.demo:
        if args.query or args.input:
            print(f"[demo mode] ignoring query input, using {DEMO_QUERY}")
        query = DEMO_QUERY
    elif args.query:
        query = args.query
    elif args.input:
        try:
            query = args.input.expanduser().read_text(encoding="utf-8").strip()
        except (OSError, UnicodeError) as exc:
            parser.error(f"Could not read query file: {exc}")
    else:
        parser.error("--query or --input is required unless --demo is set")
    query = query.strip()
    if not query:
        parser.error("Query must not be blank")

    provider = None
    try:
        max_results = _clamp_max_results(args.max_results)
        if args.summary_method == "first-sentence" and any(value is not None for value in (
            args.provider, args.model, args.base_url, args.llm_timeout, args.summary_max_tokens, args.model_params,
        )):
            raise ValueError("Provider and model options require --summary-method llm")
        if args.summary_method == "llm" and args.provider is None:
            raise ValueError("--summary-method llm requires --provider")
        params = validate_model_params(json.loads(args.model_params) if args.model_params is not None else {})
        timeout = args.llm_timeout if args.llm_timeout is not None else 120.0
        if not math.isfinite(timeout) or timeout <= 0:
            raise ValueError("--llm-timeout must be finite and positive")
        if args.summary_max_tokens is not None and args.summary_max_tokens <= 0:
            raise ValueError("--summary-max-tokens must be positive")
        if args.summary_method == "llm":
            provider = create_provider(args.provider, args.model, base_url=args.base_url, timeout=timeout)
    except (ValueError, ImportError) as exc:
        parser.error(str(exc))

    output_dir = Path(args.output).expanduser().resolve()
    artifacts = [output_dir / name for name in ("report.html", "report.md", "result.json",
                 "reproducibility/commands.sh", "reproducibility/environment.yml", "reproducibility/checksums.sha256")]
    if any(path.exists() for path in artifacts):
        print(f"[warning] Existing output files will be overwritten in {output_dir}", file=sys.stderr)
    print(f"Summary method: {args.summary_method}")

    # Fetch papers
    try:
        try:
            papers = pubmed_api.fetch_papers(query, max_results)
        except (pubmed_api.requests.RequestException, pubmed_api.ET.ParseError, ValueError) as exc:
            print(f"Error fetching papers from PubMed: {exc}", file=sys.stderr)
            sys.exit(1)
        retrieved_at = datetime.now(timezone.utc).isoformat()
        warnings = []
        # Copy dictionaries so callers' source records are not modified.
        papers = [dict(paper) for paper in papers]
        for index, paper in enumerate(papers, 1):
            if provider is not None and paper["abstract"]:
                print(f"Summarizing paper {index}/{len(papers)} (PMID {paper['pmid']})...")
            paper["summary"] = abstract_summary.summarize_abstract(
                paper, provider, max_output_tokens=args.summary_max_tokens, model_params=params,
            )
            if paper["summary"]["fallback_reason"]:
                warning = {"pmid": paper["pmid"], "reason": paper["summary"]["fallback_reason"]}
                warnings.append(warning)
                print(f"[warning] PMID {paper['pmid']}: {warning['reason']} Showing the abstract opening.", file=sys.stderr)
        config = {"method": args.summary_method, "provider": args.provider,
                  "model": provider.model if provider else None,
                  "base_url": provider.base_url if provider else None,
                  "timeout": timeout if provider else None,
                  "max_output_tokens": args.summary_max_tokens, "model_params": params,
                  "prompt_version": abstract_summary.PROMPT_VERSION if provider else None}
        print(_format_terminal_summary(query, papers))
        _write_outputs(output_dir, query, papers, config, retrieved_at, warnings, max_results)
    finally:
        if provider is not None:
            provider.close()


if __name__ == "__main__":
    main()
