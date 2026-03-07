#!/usr/bin/env python3
"""
generate_curated_skills.py — Generate curated .md profiles for top Galaxy tools
================================================================================
Reads galaxy_catalog.json and writes structured markdown profiles for the ~200
most important bioinformatics tools, organised by category.

Usage:
    python generate_curated_skills.py                # uses galaxy_catalog.json
    python generate_curated_skills.py --top 50       # only top 50
"""

from __future__ import annotations

import argparse
import json
import re
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
CATALOG_PATH = SKILL_DIR / "galaxy_catalog.json"
OUTPUT_DIR = SKILL_DIR / "galaxy_skills"

# Priority categories and their target counts
CATEGORY_TARGETS = {
    "FASTQ Quality Control": 10,
    "FASTA/FASTQ": 10,
    "Mapping": 15,
    "SAM/BAM": 10,
    "Variant Calling": 15,
    "VCF/BCF": 10,
    "Assembly": 15,
    "Annotation": 15,
    "RNA-seq": 15,
    "Metagenomic Analysis": 20,
    "QIIME2": 10,
    "Proteomics": 10,
    "Single-cell": 15,
    "Epigenetics": 5,
    "ChIP-seq": 5,
    "Nanopore": 10,
    "Statistics": 5,
    "Graph/Display Data": 5,
    "Phenotype Association": 5,
    "Evolution": 5,
    "Machine Learning": 5,
    "Virology": 5,
    "Multiple Alignments": 10,
    "deepTools": 5,
    "Imaging": 5,
    "Metabolomics": 5,
}

# Well-known tools to always include (by partial ID match)
PRIORITY_TOOLS = {
    "fastqc", "trimmomatic", "cutadapt", "fastp",
    "bwa", "bwa_mem2", "bowtie2", "hisat2", "minimap2", "star",
    "samtools", "picard", "bamtools", "bedtools",
    "freebayes", "gatk4", "bcftools", "deepvariant",
    "snpsift", "snpeff", "vep",
    "featurecounts", "htseq_count", "deseq2", "edger", "stringtie", "salmon", "kallisto",
    "kraken2", "metaphlan", "humann", "bracken", "diamond", "megahit",
    "spades", "flye", "unicycler", "quast",
    "prokka", "augustus", "maker", "bakta",
    "iqtree", "raxml", "mafft", "muscle", "clustalw",
    "maxquant", "searchgui", "peptideshaker",
    "scanpy", "cellranger", "seurat",
    "macs2", "deeptools", "diffbind",
    "plink", "regenie",
    "nanoplot", "medaka", "nanopolish",
}


def slugify(name: str) -> str:
    """Convert a tool name to a filename-safe slug."""
    slug = name.lower().strip()
    slug = re.sub(r"[^a-z0-9]+", "-", slug)
    slug = slug.strip("-")
    return slug or "unnamed"


def write_tool_profile(tool: dict, output_dir: Path):
    """Write a markdown profile for a single Galaxy tool."""
    name = tool.get("name", "Unknown")
    tool_id = tool.get("id", "")
    version = tool.get("version", "?")
    desc = tool.get("description", "")
    section = tool.get("section", "")
    edam_topics = tool.get("edam_topics", [])
    edam_ops = tool.get("edam_operations", [])
    inputs = tool.get("inputs", [])
    outputs = tool.get("outputs", [])

    slug = slugify(name)
    path = output_dir / f"{slug}.md"

    lines = [
        f"# {name}",
        "",
        f"**Galaxy Tool ID**: `{tool_id}`",
        f"**Version**: {version}",
        f"**Category**: {section}",
        "",
    ]

    if desc:
        lines.extend([f"> {desc}", ""])

    if edam_topics:
        lines.extend([f"**EDAM Topics**: {', '.join(edam_topics)}", ""])
    if edam_ops:
        lines.extend([f"**EDAM Operations**: {', '.join(edam_ops)}", ""])

    if inputs:
        lines.extend(["## Inputs", ""])
        for inp in inputs[:10]:
            if isinstance(inp, dict):
                label = inp.get("label") or inp.get("name") or "?"
                itype = inp.get("type", "")
                lines.append(f"- **{label}** ({itype})")
        lines.append("")

    if outputs:
        lines.extend(["## Outputs", ""])
        for out in outputs[:10]:
            if isinstance(out, dict):
                oname = out.get("name", "?")
                fmt = out.get("format", "?")
                lines.append(f"- **{oname}** ({fmt})")
        lines.append("")

    lines.extend([
        "## Example Query",
        "",
        f'> "Run {name} on my data"',
        "",
        "## Run via Galaxy Bridge",
        "",
        "```bash",
        f"python galaxy_bridge.py --tool-details {tool_id}",
        f"python galaxy_bridge.py --run {tool_id} --input <file> --output <dir>",
        "```",
        "",
        f"---",
        f"*Auto-generated from usegalaxy.org tool index*",
    ])

    path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return path


def _tool_slug(tool_id: str) -> str:
    """Extract the short tool name from a Galaxy tool ID.

    E.g. 'toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74' → 'fastqc'
    """
    parts = tool_id.strip("/").split("/")
    # ToolShed IDs: .../repos/owner/repo/tool/version — tool is parts[-2]
    if len(parts) >= 2:
        return parts[-2].lower()
    return parts[-1].lower()


def select_tools(catalog: dict, top_n: int) -> list[dict]:
    """Select the most important tools from the catalog."""
    tools = catalog.get("tools", [])
    selected: list[dict] = []
    selected_ids: set[str] = set()
    selected_slugs: set[str] = set()  # prevent duplicates from same tool

    # First pass: priority tools (match on slug, max 1 per priority name)
    priority_matched: set[str] = set()
    for tool in tools:
        tid = tool.get("id") or ""
        slug = _tool_slug(tid)
        name_lower = (tool.get("name") or "").lower()

        for priority in PRIORITY_TOOLS:
            if priority in priority_matched:
                continue
            # Match slug exactly, or name starts with priority
            if slug == priority or name_lower == priority:
                if slug not in selected_slugs:
                    selected.append(tool)
                    selected_ids.add(tid.lower())
                    selected_slugs.add(slug)
                    priority_matched.add(priority)
                break

    # Second pass: fill from category targets
    section_counts: dict[str, int] = {}
    for tool in tools:
        section = tool.get("section", "")
        tid = tool.get("id") or ""
        slug = _tool_slug(tid)

        if tid.lower() in selected_ids or slug in selected_slugs:
            continue

        # Check if this section has quota remaining
        section_lower = section.lower()
        for target_section, target_count in CATEGORY_TARGETS.items():
            if target_section.lower() == section_lower:
                current = section_counts.get(target_section, 0)
                if current < target_count:
                    selected.append(tool)
                    selected_ids.add(tid.lower())
                    selected_slugs.add(slug)
                    section_counts[target_section] = current + 1
                break

        if len(selected) >= top_n:
            break

    return selected[:top_n]


def main():
    parser = argparse.ArgumentParser(description="Generate curated Galaxy tool profiles")
    parser.add_argument("--top", type=int, default=200, help="Number of tools to profile (default: 200)")
    args = parser.parse_args()

    if not CATALOG_PATH.exists():
        print(f"ERROR: {CATALOG_PATH} not found. Run generate_galaxy_catalog.py first.", file=sys.stderr)
        sys.exit(1)

    catalog = json.loads(CATALOG_PATH.read_text(encoding="utf-8"))
    tools = select_tools(catalog, args.top)

    OUTPUT_DIR.mkdir(parents=True, exist_ok=True)

    for tool in tools:
        write_tool_profile(tool, OUTPUT_DIR)

    print(f"Wrote {len(tools)} curated tool profiles to {OUTPUT_DIR}/")

    # Write index
    index_path = OUTPUT_DIR / "INDEX.md"
    lines = [
        "# Galaxy Skills — Curated Tool Profiles",
        "",
        f"**{len(tools)} tools** profiled from usegalaxy.org.",
        "",
        "| Tool | Category | Description |",
        "|------|----------|-------------|",
    ]
    for tool in sorted(tools, key=lambda t: (t.get("section", ""), t.get("name", ""))):
        name = tool.get("name", "?")
        slug = slugify(name)
        section = tool.get("section", "")
        desc = (tool.get("description") or "")[:60]
        lines.append(f"| [{name}]({slug}.md) | {section} | {desc} |")

    index_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(f"Wrote index: {index_path}")


if __name__ == "__main__":
    main()
