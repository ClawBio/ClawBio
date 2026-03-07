#!/usr/bin/env python3
"""
galaxy_bridge.py — Galaxy tool discovery and execution for ClawBio
===================================================================
Search, inspect, and run 1,770+ Galaxy bioinformatics tools from the
command line.  Uses BioBlend (Galaxy's Python SDK) for API calls and
a bundled galaxy_catalog.json for offline discovery.

Usage:
    python galaxy_bridge.py --search "metagenomics"
    python galaxy_bridge.py --list-categories
    python galaxy_bridge.py --tool-details <tool_id>
    python galaxy_bridge.py --run <tool_id> --input <file> --output <dir>
    python galaxy_bridge.py --demo
"""

from __future__ import annotations

import argparse
import hashlib
import json
import os
import sys
import time
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
CATALOG_PATH = SKILL_DIR / "galaxy_catalog.json"
DEMO_DIR = SKILL_DIR / "demo"
CURATED_DIR = SKILL_DIR / "galaxy_skills"

# Galaxy defaults
DEFAULT_GALAXY_URL = "https://usegalaxy.org"


# ---------------------------------------------------------------------------
# Catalog helpers
# ---------------------------------------------------------------------------


def load_catalog() -> dict:
    """Load the bundled galaxy_catalog.json."""
    if not CATALOG_PATH.exists():
        print(f"ERROR: Catalog not found at {CATALOG_PATH}", file=sys.stderr)
        print("Run generate_galaxy_catalog.py first.", file=sys.stderr)
        sys.exit(1)
    return json.loads(CATALOG_PATH.read_text(encoding="utf-8"))


def search_catalog(query: str, catalog: dict, max_results: int = 20) -> list[dict]:
    """Search the local catalog by name, description, section, and EDAM terms."""
    query_lower = query.lower()
    terms = query_lower.split()
    scored: list[tuple[float, dict]] = []

    for tool in catalog.get("tools", []):
        score = 0.0
        name = (tool.get("name") or "").lower()
        desc = (tool.get("description") or "").lower()
        section = (tool.get("section") or "").lower()
        edam_topics = " ".join(tool.get("edam_topics", [])).lower()
        edam_ops = " ".join(tool.get("edam_operations", [])).lower()
        tool_id = (tool.get("id") or "").lower()

        searchable = f"{name} {desc} {section} {edam_topics} {edam_ops} {tool_id}"

        # Exact phrase match (highest weight)
        if query_lower in name:
            score += 10.0
        if query_lower in desc:
            score += 5.0
        if query_lower in section:
            score += 3.0

        # Individual term matches
        for term in terms:
            if term in name:
                score += 3.0
            if term in desc:
                score += 1.5
            if term in searchable:
                score += 0.5

        if score > 0:
            scored.append((score, tool))

    scored.sort(key=lambda x: -x[0])
    return [t for _, t in scored[:max_results]]


def list_categories(catalog: dict) -> dict[str, int]:
    """Count tools per section/category."""
    cats: dict[str, int] = {}
    for tool in catalog.get("tools", []):
        section = tool.get("section") or "Uncategorized"
        cats[section] = cats.get(section, 0) + 1
    return dict(sorted(cats.items(), key=lambda x: -x[1]))


def get_tool_details(tool_id: str, catalog: dict) -> dict | None:
    """Find a specific tool by ID (exact or partial match)."""
    tool_id_lower = tool_id.lower()
    for tool in catalog.get("tools", []):
        tid = (tool.get("id") or "").lower()
        if tid == tool_id_lower or tool_id_lower in tid:
            return tool
    return None


# ---------------------------------------------------------------------------
# Galaxy API (BioBlend)
# ---------------------------------------------------------------------------


def get_galaxy_instance():
    """Create a BioBlend GalaxyInstance from environment variables."""
    url = os.environ.get("GALAXY_URL", DEFAULT_GALAXY_URL)
    key = os.environ.get("GALAXY_API_KEY", "")
    if not key:
        print("WARNING: GALAXY_API_KEY not set. API operations will fail.", file=sys.stderr)
        print("Set it: export GALAXY_API_KEY=your_key_here", file=sys.stderr)
        print(f"Register at: {url}/user/api_key", file=sys.stderr)
        return None

    try:
        from bioblend.galaxy import GalaxyInstance  # type: ignore[import-untyped]
        return GalaxyInstance(url=url, key=key)
    except ImportError:
        print("ERROR: bioblend not installed. Run: pip install bioblend", file=sys.stderr)
        return None


def run_tool_on_galaxy(
    tool_id: str,
    input_path: Path,
    output_dir: Path,
    galaxy_instance=None,
) -> dict:
    """Upload input, run a Galaxy tool, download results."""
    if galaxy_instance is None:
        galaxy_instance = get_galaxy_instance()
        if galaxy_instance is None:
            return {"error": "Cannot connect to Galaxy — check GALAXY_API_KEY"}

    gi = galaxy_instance
    output_dir.mkdir(parents=True, exist_ok=True)

    # Create a new history
    ts = datetime.now(timezone.utc).strftime("%Y%m%d_%H%M%S")
    history_name = f"ClawBio_{tool_id}_{ts}"
    history = gi.histories.create_history(name=history_name)
    history_id = history["id"]

    try:
        # Upload input file
        print(f"  Uploading {input_path.name} to Galaxy...")
        upload = gi.tools.upload_file(str(input_path), history_id)
        dataset_id = upload["outputs"][0]["id"]

        # Wait for upload to complete
        _wait_for_dataset(gi, history_id, dataset_id)

        # Run the tool
        print(f"  Running {tool_id}...")
        inputs = {"input": {"src": "hda", "id": dataset_id}}
        result = gi.tools.run_tool(history_id, tool_id, inputs)

        # Wait for outputs
        output_ids = [o["id"] for o in result.get("outputs", [])]
        for oid in output_ids:
            _wait_for_dataset(gi, history_id, oid)

        # Download outputs
        galaxy_out = output_dir / "galaxy_outputs"
        galaxy_out.mkdir(exist_ok=True)
        downloaded = []
        for oid in output_ids:
            ds = gi.datasets.show_dataset(oid)
            ext = ds.get("extension", "dat")
            fname = f"{ds.get('name', oid)}.{ext}"
            out_path = galaxy_out / fname
            gi.datasets.download_dataset(oid, file_path=str(out_path), use_default_filename=False)
            downloaded.append(str(out_path))
            print(f"  Downloaded: {out_path.name}")

        # Write reproducibility bundle
        _write_reproducibility(output_dir, tool_id, input_path, downloaded, gi.base_url)

        return {
            "status": "success",
            "tool_id": tool_id,
            "history": history_name,
            "outputs": downloaded,
            "output_dir": str(output_dir),
        }

    except Exception as e:
        return {"status": "error", "tool_id": tool_id, "error": str(e)}

    finally:
        # Clean up history on Galaxy
        try:
            gi.histories.delete_history(history_id, purge=True)
        except Exception:
            pass


def _wait_for_dataset(gi, history_id: str, dataset_id: str, timeout: int = 300):
    """Poll until dataset is ready."""
    start = time.time()
    while time.time() - start < timeout:
        ds = gi.datasets.show_dataset(dataset_id)
        state = ds.get("state", "")
        if state == "ok":
            return
        if state in ("error", "discarded", "failed_metadata"):
            raise RuntimeError(f"Dataset {dataset_id} failed: {state}")
        time.sleep(5)
    raise TimeoutError(f"Dataset {dataset_id} timed out after {timeout}s")


def _write_reproducibility(
    output_dir: Path,
    tool_id: str,
    input_path: Path,
    outputs: list[str],
    galaxy_url: str,
):
    """Write commands.sh, environment.yml, and checksums."""
    repro = output_dir / "reproducibility"
    repro.mkdir(exist_ok=True)

    # commands.sh
    (repro / "commands.sh").write_text(
        f"#!/usr/bin/env bash\n"
        f"# Reproduce this analysis on Galaxy\n"
        f"# Galaxy server: {galaxy_url}\n"
        f"# Tool: {tool_id}\n"
        f"# Date: {datetime.now(timezone.utc).isoformat()}\n\n"
        f"python galaxy_bridge.py --run {tool_id} --input {input_path} --output {output_dir}\n",
        encoding="utf-8",
    )

    # environment.yml
    (repro / "environment.yml").write_text(
        f"galaxy_url: {galaxy_url}\n"
        f"tool_id: {tool_id}\n"
        f"date: {datetime.now(timezone.utc).isoformat()}\n"
        f"bioblend_required: true\n",
        encoding="utf-8",
    )

    # checksums
    lines = []
    for fp in [str(input_path)] + outputs:
        p = Path(fp)
        if p.exists():
            sha = hashlib.sha256(p.read_bytes()).hexdigest()
            lines.append(f"{sha}  {p.name}")
    (repro / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ---------------------------------------------------------------------------
# Demo mode
# ---------------------------------------------------------------------------


def run_demo(output_dir: Path | None = None):
    """Run a simulated FastQC demo with pre-cached results."""
    if output_dir is None:
        output_dir = DEMO_DIR

    output_dir.mkdir(parents=True, exist_ok=True)

    print()
    print("Galaxy Bridge — Demo Mode (offline)")
    print("=" * 40)
    print(f"Tool: FastQC v0.74+galaxy1")
    print(f"Input: demo/demo_reads.fq (bundled synthetic FASTQ, 1000 reads)")

    # Write demo FASTQ if not present
    demo_fq = DEMO_DIR / "demo_reads.fq"
    if not demo_fq.exists():
        _generate_demo_fastq(demo_fq)

    # Write demo output
    demo_html = output_dir / "fastqc_demo_output.html"
    if not demo_html.exists():
        _generate_demo_html(demo_html)

    print(f"Output: {demo_html}")
    print()
    print("Result: PASS — Per base sequence quality")
    print("        PASS — Per sequence quality scores")
    print("        WARN — Per base sequence content (normal for Illumina)")
    print("        PASS — Sequence length distribution")
    print()

    # Write reproducibility bundle
    repro = output_dir / "reproducibility"
    repro.mkdir(exist_ok=True)
    (repro / "commands.sh").write_text(
        "#!/usr/bin/env bash\n"
        "# Demo mode — no Galaxy API call made\n"
        "python galaxy_bridge.py --demo\n",
        encoding="utf-8",
    )

    # result.json
    result = {
        "mode": "demo",
        "tool": "fastqc",
        "tool_version": "0.74+galaxy1",
        "input": "demo/demo_reads.fq",
        "output": str(demo_html),
        "qc_modules": {
            "per_base_sequence_quality": "PASS",
            "per_sequence_quality_scores": "PASS",
            "per_base_sequence_content": "WARN",
            "sequence_length_distribution": "PASS",
            "overrepresented_sequences": "PASS",
            "adapter_content": "PASS",
        },
        "timestamp": datetime.now(timezone.utc).isoformat(),
    }
    (output_dir / "result.json").write_text(
        json.dumps(result, indent=2) + "\n", encoding="utf-8"
    )

    print(f"Reproducibility bundle written to {repro}/")
    print()
    return result


def _generate_demo_fastq(path: Path):
    """Generate a small synthetic FASTQ file for demo purposes."""
    import random
    random.seed(42)
    lines = []
    for i in range(1000):
        seq = "".join(random.choices("ACGT", k=150))
        qual = "".join(chr(random.randint(53, 73)) for _ in range(150))  # Phred 20-40
        lines.append(f"@demo_read_{i+1}")
        lines.append(seq)
        lines.append("+")
        lines.append(qual)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def _generate_demo_html(path: Path):
    """Generate a minimal FastQC-like HTML report for demo."""
    html = """<!DOCTYPE html>
<html>
<head><title>FastQC Report — Demo</title></head>
<body>
<h1>FastQC Report — Demo Mode</h1>
<p>Generated by ClawBio Galaxy Bridge (demo mode, no Galaxy API call).</p>
<h2>Summary</h2>
<table border="1" cellpadding="4">
<tr><th>Module</th><th>Status</th></tr>
<tr><td>Per base sequence quality</td><td style="color:green">PASS</td></tr>
<tr><td>Per sequence quality scores</td><td style="color:green">PASS</td></tr>
<tr><td>Per base sequence content</td><td style="color:orange">WARN</td></tr>
<tr><td>Sequence length distribution</td><td style="color:green">PASS</td></tr>
<tr><td>Overrepresented sequences</td><td style="color:green">PASS</td></tr>
<tr><td>Adapter content</td><td style="color:green">PASS</td></tr>
</table>
<h2>Basic Statistics</h2>
<table border="1" cellpadding="4">
<tr><td>Filename</td><td>demo_reads.fq</td></tr>
<tr><td>Total Sequences</td><td>1000</td></tr>
<tr><td>Sequence length</td><td>150</td></tr>
<tr><td>%GC</td><td>50</td></tr>
</table>
<p><em>This is a demo report. Run with --run fastqc for real analysis on usegalaxy.org.</em></p>
</body>
</html>"""
    path.write_text(html, encoding="utf-8")


# ---------------------------------------------------------------------------
# Report generator
# ---------------------------------------------------------------------------


def write_report(output_dir: Path, tool_id: str, results: dict):
    """Write a markdown report for a Galaxy tool run."""
    report_path = output_dir / "report.md"
    ts = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")

    lines = [
        f"# Galaxy Bridge Report: {tool_id}",
        "",
        f"**Date**: {ts}",
        f"**Skill**: galaxy-bridge",
        f"**Galaxy tool**: {tool_id}",
        "",
        "## Results",
        "",
    ]

    if results.get("status") == "success":
        lines.append(f"Tool executed successfully on Galaxy.")
        lines.append("")
        lines.append("### Outputs")
        for o in results.get("outputs", []):
            lines.append(f"- `{Path(o).name}`")
    elif results.get("mode") == "demo":
        lines.append("Demo mode — no Galaxy API call made.")
        lines.append("")
        for module, status in results.get("qc_modules", {}).items():
            lines.append(f"- **{module}**: {status}")
    else:
        lines.append(f"Error: {results.get('error', 'unknown')}")

    lines.extend([
        "",
        "## Reproducibility",
        "",
        "See `reproducibility/commands.sh` to re-run this analysis.",
        "",
        "## Disclaimer",
        "",
        "ClawBio is a research and educational tool. It is not a medical device "
        "and does not provide clinical diagnoses. Consult a healthcare "
        "professional before making any medical decisions.",
    ])

    report_path.write_text("\n".join(lines) + "\n", encoding="utf-8")
    return report_path


# ---------------------------------------------------------------------------
# CLI
# ---------------------------------------------------------------------------


def main():
    parser = argparse.ArgumentParser(
        description="Galaxy Bridge — search, inspect, and run Galaxy tools",
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--search", metavar="QUERY", help="Search Galaxy tools by keyword")
    parser.add_argument("--list-categories", action="store_true", help="List tool categories with counts")
    parser.add_argument("--tool-details", metavar="TOOL_ID", help="Show details for a specific tool")
    parser.add_argument("--run", metavar="TOOL_ID", help="Run a Galaxy tool")
    parser.add_argument("--input", metavar="FILE", help="Input file for --run")
    parser.add_argument("--output", metavar="DIR", help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run FastQC demo (offline, no API key)")
    parser.add_argument("--max-results", type=int, default=20, help="Max search results (default: 20)")
    parser.add_argument("--recommend", metavar="TASK", help="Smart tool recommendation for a task description")
    parser.add_argument("--format", metavar="EXT", help="Input format hint for --recommend (e.g. .fastq, .vcf, .bam)")
    parser.add_argument("--workflow", metavar="TASK", help="Suggest multi-step pipelines for a task")

    args = parser.parse_args()

    # Default: show help
    if not any([args.search, args.list_categories, args.tool_details, args.run, args.demo, args.recommend, args.workflow]):
        parser.print_help()
        sys.exit(0)

    # Demo mode
    if args.demo:
        out = Path(args.output) if args.output else None
        result = run_demo(out)
        if out:
            write_report(out, "fastqc", result)
        return

    # All other modes need the catalog
    catalog = load_catalog()

    # Smart recommendation
    if args.recommend:
        from tool_recommender import recommend_tool, suggest_workflow

        results = recommend_tool(args.recommend, catalog, input_format=args.format)
        if not results:
            print(f"No tools found for: {args.recommend}")
            return

        print(f"\n{'=' * 60}")
        print(f"Galaxy Tool Recommendations")
        print(f"Task: \"{args.recommend}\"")
        if args.format:
            print(f"Input format: {args.format}")
        print(f"{'=' * 60}\n")

        for i, rec in enumerate(results[:5], 1):
            print(f"  {i}. {rec['name']} (v{rec['version']})")
            print(f"     Score: {rec['score']} | {rec['explanation']}")
            print(f"     Category: {rec['section']}")
            if rec["description"]:
                print(f"     {rec['description'][:80]}")
            if rec["edam_labels"]:
                print(f"     EDAM: {', '.join(rec['edam_labels'][:3])}")
            if rec["version_count"] > 1:
                print(f"     Versions: {rec['version_count']} (showing latest)")
            print(f"     ID: {rec['tool_id']}")
            print()

        # Also suggest workflows
        workflows = suggest_workflow(args.recommend, args.format)
        if workflows:
            print(f"{'=' * 60}")
            print("Suggested Pipelines")
            print(f"{'=' * 60}\n")
            for wf in workflows[:2]:
                print(f"  {wf['name']}")
                print(f"  {wf['description']}")
                print()
                for j, step in enumerate(wf["steps"], 1):
                    print(f"    {j}. {step['tool']} — {step['purpose']}")
                print()
        return

    # Workflow suggestion
    if args.workflow:
        from tool_recommender import suggest_workflow

        workflows = suggest_workflow(args.workflow, args.format)
        if not workflows:
            print(f"No workflow templates match: {args.workflow}")
            return
        print(f"\nSuggested Pipelines for: \"{args.workflow}\"\n")
        for wf in workflows:
            print(f"  {wf['name']} (match: {wf['match_score']})")
            print(f"  {wf['description']}")
            print()
            for j, step in enumerate(wf["steps"], 1):
                print(f"    {j}. {step['tool']} — {step['purpose']}")
            print()
        return

    # Search
    if args.search:
        results = search_catalog(args.search, catalog, args.max_results)
        if not results:
            print(f"No tools found for: {args.search}")
            return
        print(f"\nFound {len(results)} Galaxy tools matching '{args.search}':\n")
        for i, tool in enumerate(results, 1):
            name = tool.get("name", "?")
            version = tool.get("version", "?")
            desc = tool.get("description", "")
            section = tool.get("section", "")
            tid = tool.get("id", "")
            print(f"  {i:2d}. {name} (v{version})")
            if desc:
                print(f"      {desc[:80]}")
            if section:
                print(f"      Category: {section}")
            print(f"      ID: {tid}")
            print()
        return

    # List categories
    if args.list_categories:
        cats = list_categories(catalog)
        total = sum(cats.values())
        print(f"\nGalaxy Tool Categories ({len(cats)} categories, {total} tools):\n")
        for cat, count in cats.items():
            print(f"  {cat:40s} {count:4d} tools")
        print()
        return

    # Tool details
    if args.tool_details:
        tool = get_tool_details(args.tool_details, catalog)
        if not tool:
            print(f"Tool not found: {args.tool_details}")
            return
        print(f"\n{'=' * 60}")
        print(f"Tool: {tool.get('name', '?')}")
        print(f"ID:   {tool.get('id', '?')}")
        print(f"Version: {tool.get('version', '?')}")
        print(f"Section: {tool.get('section', '?')}")
        if tool.get("description"):
            print(f"Description: {tool['description']}")
        if tool.get("edam_topics"):
            print(f"EDAM Topics: {', '.join(tool['edam_topics'])}")
        if tool.get("edam_operations"):
            print(f"EDAM Operations: {', '.join(tool['edam_operations'])}")
        if tool.get("inputs"):
            print(f"\nInputs:")
            for inp in tool["inputs"][:10]:
                if isinstance(inp, dict):
                    print(f"  - {inp.get('name', '?')}: {inp.get('label', inp.get('type', '?'))}")
                else:
                    print(f"  - {inp}")
        if tool.get("outputs"):
            print(f"\nOutputs:")
            for out in tool["outputs"][:10]:
                if isinstance(out, dict):
                    print(f"  - {out.get('name', '?')}: {out.get('format', '?')}")
                else:
                    print(f"  - {out}")
        print(f"{'=' * 60}\n")
        return

    # Run tool
    if args.run:
        if not args.input:
            print("ERROR: --input required with --run", file=sys.stderr)
            sys.exit(1)
        input_path = Path(args.input)
        if not input_path.exists():
            print(f"ERROR: Input file not found: {input_path}", file=sys.stderr)
            sys.exit(1)
        output_dir = Path(args.output) if args.output else Path(f"/tmp/galaxy_{args.run}_{int(time.time())}")

        print(f"\nGalaxy Bridge — Running {args.run}")
        print(f"  Input:  {input_path}")
        print(f"  Output: {output_dir}")
        print()

        result = run_tool_on_galaxy(args.run, input_path, output_dir)

        if result.get("status") == "success":
            write_report(output_dir, args.run, result)
            print(f"\nDone. Report: {output_dir / 'report.md'}")
        else:
            print(f"\nError: {result.get('error', 'unknown')}", file=sys.stderr)
            sys.exit(1)


if __name__ == "__main__":
    main()
