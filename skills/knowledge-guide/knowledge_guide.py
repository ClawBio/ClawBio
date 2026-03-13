#!/usr/bin/env python3
"""
knowledge_guide.py — GTN-backed educational guide for ClawBio
=============================================================
Explains bioinformatics concepts grounded in Galaxy Training Network
tutorials. Supports free-text queries, structured lookups, and
report-embedded "Learn More" sections.

Usage:
    python knowledge_guide.py --query "what is variant calling?" --output /tmp/kg
    python knowledge_guide.py --topic variant-analysis --output /tmp/kg
    python knowledge_guide.py --tool fastqc --output /tmp/kg
    python knowledge_guide.py --concept "polygenic risk" --output /tmp/kg
    python knowledge_guide.py --skill gwas-prs --output /tmp/kg
    python knowledge_guide.py --query "RNA-seq DE" --deep --output /tmp/kg
    python knowledge_guide.py --demo --output /tmp/kg_demo
"""
from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from clawbio.common.report import generate_report_header, generate_report_footer, write_result_json
from clawbio.common.html_report import HtmlReportBuilder

CACHE_PATH = SKILL_DIR / "gtn_cache.json"
RECS_PATH = SKILL_DIR / "skill_recommendations.json"

# Embedded mini-cache for demo mode (no network needed)
DEMO_CACHE = {
    "topics": [
        {
            "name": "variant-analysis",
            "title": "Variant Analysis",
            "summary": "Detecting genetic variants in sequencing data",
            "tutorials": [
                {
                    "title": "Calling very rare variants",
                    "name": "rare-variants",
                    "topic": "variant-analysis",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": [
                        "Process duplex sequencing data",
                        "Identify rare variants without diploid assumptions",
                    ],
                    "tools": ["fastqc", "bwa_mem", "naive_variant_caller"],
                    "url": "https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
                },
                {
                    "title": "Exome sequencing data analysis",
                    "name": "exome-seq",
                    "topic": "variant-analysis",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": [
                        "Perform variant calling on whole-exome data",
                        "Annotate variants with functional impact",
                    ],
                    "tools": ["bwa_mem", "freebayes", "snpeff", "gemini"],
                    "url": "https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/exome-seq/tutorial.html",
                },
            ],
        },
        {
            "name": "transcriptomics",
            "title": "Transcriptomics",
            "summary": "Gene expression analysis approaches",
            "tutorials": [
                {
                    "title": "Reference-based RNA-Seq data analysis",
                    "name": "ref-based",
                    "topic": "transcriptomics",
                    "time_estimation": "4H",
                    "level": "Intermediate",
                    "objectives": [
                        "Perform RNA-seq alignment with HISAT2",
                        "Identify differentially expressed genes with DESeq2",
                    ],
                    "tools": ["hisat2", "featurecounts", "deseq2"],
                    "url": "https://training.galaxyproject.org/training-material/topics/transcriptomics/tutorials/ref-based/tutorial.html",
                },
            ],
        },
        {
            "name": "single-cell",
            "title": "Single Cell",
            "summary": "Single-cell RNA-seq and single-cell analysis",
            "tutorials": [
                {
                    "title": "Clustering 3K PBMCs with Scanpy",
                    "name": "scrna-scanpy-pbmc3k",
                    "topic": "single-cell",
                    "time_estimation": "3H",
                    "level": "Intermediate",
                    "objectives": [
                        "Perform quality control on scRNA-seq data",
                        "Cluster cells and identify marker genes",
                    ],
                    "tools": ["scanpy_filter_cells", "scanpy_normalise", "scanpy_find_cluster"],
                    "url": "https://training.galaxyproject.org/training-material/topics/single-cell/tutorials/scrna-scanpy-pbmc3k/tutorial.html",
                },
            ],
        },
    ],
    "tool_index": {
        "devteam/fastqc/fastqc": [
            {
                "path": "variant-analysis/rare-variants",
                "title": "Calling very rare variants",
                "topic": "Variant Analysis",
                "url": "https://training.galaxyproject.org/training-material/topics/variant-analysis/tutorials/rare-variants/tutorial.html",
            },
        ],
    },
}


def _load_cache(demo: bool = False) -> dict:
    """Load the GTN cache (or demo cache)."""
    if demo:
        return DEMO_CACHE
    if not CACHE_PATH.exists():
        print(f"ERROR: Cache not found at {CACHE_PATH}", file=sys.stderr)
        print("Run: python gtn_cache_builder.py --refresh", file=sys.stderr)
        sys.exit(1)
    return json.loads(CACHE_PATH.read_text())


def _load_recommendations() -> dict:
    """Load skill recommendations if available."""
    if RECS_PATH.exists():
        return json.loads(RECS_PATH.read_text())
    return {}


def _format_tutorial_md(tut: dict, index: int = 1) -> str:
    """Format a single tutorial as markdown."""
    lines = [f"### {index}. {tut.get('title', 'Untitled')}"]

    meta = []
    if tut.get("topic"):
        meta.append(f"**Topic:** {tut['topic']}")
    if tut.get("time_estimation") or tut.get("time"):
        meta.append(f"**Time:** {tut.get('time_estimation') or tut.get('time')}")
    if tut.get("level"):
        meta.append(f"**Level:** {tut['level']}")
    if tut.get("relevance_score"):
        meta.append(f"**Relevance:** {tut['relevance_score']}")
    if meta:
        lines.append(" | ".join(meta))

    objectives = tut.get("objectives", [])
    if isinstance(objectives, list) and objectives:
        lines.append("")
        lines.append("**Learning objectives:**")
        for obj in objectives:
            lines.append(f"- {obj}")

    url = tut.get("url", "")
    if url:
        lines.append("")
        lines.append(f"[Full tutorial]({url})")

    lines.append("")
    return "\n".join(lines)


def _format_deep_content(content: dict) -> str:
    """Format deep tutorial content as markdown."""
    lines = [f"## Deep Content: {content.get('title', '')}"]
    body = content.get("content", "")
    if not body:
        lines.append("*(No content available for this tutorial)*")
        return "\n".join(lines)

    # Extract first ~2000 chars of content as a preview
    text = re.sub(r"<[^>]+>", "", body)  # strip HTML tags
    text = re.sub(r"\n{3,}", "\n\n", text)  # collapse whitespace
    preview = text[:2000]
    if len(text) > 2000:
        preview += "\n\n*(Content truncated — see full tutorial link above)*"

    lines.append("")
    lines.append(preview)
    return "\n".join(lines)


def run_query(args: argparse.Namespace) -> dict:
    """Execute a knowledge-guide query and return results."""
    from query_engine import search_tutorials, lookup_topic, lookup_tool, lookup_concept

    cache = _load_cache(demo=args.demo)
    recs = _load_recommendations()

    mode = "unknown"
    query_text = ""
    results = []

    if args.demo:
        mode = "demo"
        query_text = "variant calling"
        results = search_tutorials("variant calling", cache)
    elif args.query:
        mode = "free-text"
        query_text = args.query
        results = search_tutorials(args.query, cache)
    elif args.topic:
        mode = "topic"
        query_text = args.topic
        topic = lookup_topic(args.topic, cache)
        if topic:
            results = [
                {
                    "title": t["title"],
                    "name": t.get("name", ""),
                    "topic": topic["name"],
                    "time_estimation": t.get("time_estimation", ""),
                    "level": t.get("level", ""),
                    "objectives": t.get("objectives", []),
                    "url": t.get("url", ""),
                }
                for t in topic["tutorials"]
            ]
    elif args.tool:
        mode = "tool"
        query_text = args.tool
        results = lookup_tool(args.tool, cache)
    elif args.concept:
        mode = "concept"
        query_text = args.concept
        results = lookup_concept(args.concept, cache, recs)
    elif args.skill:
        mode = "skill"
        query_text = args.skill
        if args.skill in recs:
            entry = recs[args.skill]
            results = entry.get("tutorials", [])

    # Deep content fetch (live, standalone queries only)
    deep_content = None
    if getattr(args, "deep", False) and results and not args.demo:
        top = results[0]
        topic_id = top.get("topic", "")
        tut_name = top.get("name", top.get("path", "").split("/")[-1] if "/" in top.get("path", "") else "")
        if topic_id and tut_name:
            from gtn_client import fetch_tutorial_content
            deep_content = fetch_tutorial_content(topic_id, tut_name)

    return {
        "mode": mode,
        "query": query_text,
        "matches": len(results),
        "tutorials": results,
        "deep_content": deep_content,
        "deep_content_fetched": deep_content is not None,
    }


def generate_report(result: dict, output_dir: Path) -> None:
    """Generate report.md, report.html, and result.json."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # --- report.md ---
    header = generate_report_header(
        title="Knowledge Guide",
        skill_name="knowledge-guide",
        extra_metadata={
            "Query": result["query"],
            "Mode": result["mode"],
            "Matches": str(result["matches"]),
        },
    )
    body_lines = []
    if not result["tutorials"]:
        body_lines.append("No matching tutorials found for this query.\n")
    else:
        body_lines.append("## Matching Tutorials\n")
        for i, tut in enumerate(result["tutorials"], 1):
            body_lines.append(_format_tutorial_md(tut, i))

    if result.get("deep_content"):
        body_lines.append(_format_deep_content(result["deep_content"]))

    footer = generate_report_footer()
    report_md = header + "\n".join(body_lines) + footer
    (output_dir / "report.md").write_text(report_md)

    # --- report.html ---
    builder = HtmlReportBuilder(title="Knowledge Guide", skill="knowledge-guide")
    builder.add_header_block(
        title="Knowledge Guide",
        subtitle=f"Query: {result['query']} ({result['mode']} mode)",
    )
    if result["tutorials"]:
        rows = []
        for tut in result["tutorials"]:
            url = tut.get("url", "")
            title = tut.get("title", "Tutorial")
            time_est = tut.get("time_estimation", tut.get("time", ""))
            level = tut.get("level", "")
            link = f'<a href="{url}" target="_blank">Open</a>' if url else ""
            rows.append([title, tut.get("topic", ""), time_est, level, link])
        builder.add_section(heading="Matching Tutorials")
        builder.add_table(
            headers=["Tutorial", "Topic", "Time", "Level", "Link"],
            rows=rows,
        )
    else:
        builder.add_section(heading="Results")
        builder.add_paragraph("No matching tutorials found.")
    builder.add_disclaimer()
    builder.add_footer_block(skill="knowledge-guide", version="0.1.0")
    (output_dir / "report.html").write_text(builder.render())

    # --- result.json ---
    data = {k: v for k, v in result.items() if k != "deep_content"}
    write_result_json(
        output_dir=str(output_dir),
        skill="knowledge-guide",
        version="0.1.0",
        summary={
            "query": result["query"],
            "mode": result["mode"],
            "matches": result["matches"],
            "deep_content_fetched": result["deep_content_fetched"],
        },
        data=data,
    )

    # --- Deep content file (if fetched) ---
    if result.get("deep_content"):
        deep_dir = output_dir / "tutorials"
        deep_dir.mkdir(exist_ok=True)
        content = result["deep_content"]
        filename = f"{content.get('name', 'tutorial')}.md"
        (deep_dir / filename).write_text(_format_deep_content(content))

    print(f"Report written to {output_dir}/")


def main():
    parser = argparse.ArgumentParser(
        description="Knowledge Guide — GTN-backed educational guide for ClawBio"
    )
    parser.add_argument("--query", type=str, help="Free-text query")
    parser.add_argument("--topic", type=str, help="GTN topic ID (e.g. variant-analysis)")
    parser.add_argument("--tool", type=str, help="Galaxy tool name (e.g. fastqc)")
    parser.add_argument("--concept", type=str, help="Concept term (e.g. 'polygenic risk')")
    parser.add_argument("--skill", type=str, help="ClawBio skill name for Learn More lookup")
    parser.add_argument("--deep", action="store_true", help="Fetch full tutorial content (live)")
    parser.add_argument("--demo", action="store_true", help="Run demo with embedded data")
    parser.add_argument("--output", type=str, required=True, help="Output directory")

    args = parser.parse_args()

    # Require at least one query mode
    if not any([args.query, args.topic, args.tool, args.concept, args.skill, args.demo]):
        parser.error("Provide one of: --query, --topic, --tool, --concept, --skill, or --demo")

    result = run_query(args)
    generate_report(result, Path(args.output))


if __name__ == "__main__":
    main()
