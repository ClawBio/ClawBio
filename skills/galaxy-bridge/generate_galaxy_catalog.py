#!/usr/bin/env python3
"""
generate_galaxy_catalog.py — Fetch Galaxy tool index from usegalaxy.org
=======================================================================
Queries the Galaxy API and writes galaxy_catalog.json with metadata for
all available tools.  The catalog is committed to the repo so search
works offline.

Usage:
    python generate_galaxy_catalog.py                     # fetch from usegalaxy.org
    python generate_galaxy_catalog.py --url https://usegalaxy.eu  # alternate server

Requires: requests (or bioblend, which bundles requests)
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
CATALOG_PATH = SKILL_DIR / "galaxy_catalog.json"

DEFAULT_URL = "https://usegalaxy.org"


def fetch_tools(galaxy_url: str) -> list[dict]:
    """Fetch all tools from the Galaxy API."""
    try:
        import requests  # bioblend installs requests
    except ImportError:
        print("ERROR: requests not installed. Run: pip install requests", file=sys.stderr)
        sys.exit(1)

    api_url = f"{galaxy_url.rstrip('/')}/api/tools?in_panel=false"
    print(f"Fetching tools from {api_url} ...")
    resp = requests.get(api_url, timeout=120)
    resp.raise_for_status()
    raw_tools = resp.json()
    print(f"  Received {len(raw_tools)} raw entries")

    tools: list[dict] = []
    seen_ids: set[str] = set()

    for t in raw_tools:
        # Skip non-tool entries (labels, sections, etc.)
        if not isinstance(t, dict):
            continue
        tool_id = t.get("id", "")
        if not tool_id or tool_id in seen_ids:
            continue
        # Skip Galaxy built-in tools (upload, __*)
        if tool_id.startswith("__") or tool_id == "upload1":
            continue

        seen_ids.add(tool_id)

        entry = {
            "id": tool_id,
            "name": t.get("name", ""),
            "description": t.get("description", ""),
            "version": t.get("version", ""),
            "section": t.get("panel_section_name", ""),
            "edam_topics": t.get("edam_topics", []) or [],
            "edam_operations": t.get("edam_operations", []) or [],
        }

        # Compact input/output info (names + types only, skip full schemas)
        inputs = []
        for inp in (t.get("inputs") or []):
            if isinstance(inp, dict):
                inputs.append({
                    "name": inp.get("name", ""),
                    "type": inp.get("type", ""),
                    "label": inp.get("label", ""),
                })
        entry["inputs"] = inputs[:20]  # cap to avoid bloat

        outputs = []
        for out in (t.get("outputs") or []):
            if isinstance(out, dict):
                outputs.append({
                    "name": out.get("name", ""),
                    "format": out.get("format", ""),
                })
        entry["outputs"] = outputs[:10]

        tools.append(entry)

    return tools


def build_catalog(tools: list[dict], galaxy_url: str) -> dict:
    """Assemble the catalog JSON object."""
    # Count by section
    sections: dict[str, int] = {}
    for t in tools:
        sec = t.get("section") or "Uncategorized"
        sections[sec] = sections.get(sec, 0) + 1

    return {
        "version": "1.0.0",
        "generated_by": "generate_galaxy_catalog.py",
        "galaxy_url": galaxy_url,
        "generated_at": datetime.now(timezone.utc).isoformat(),
        "tool_count": len(tools),
        "section_count": len(sections),
        "sections": dict(sorted(sections.items(), key=lambda x: -x[1])),
        "tools": tools,
    }


def main():
    parser = argparse.ArgumentParser(description="Generate Galaxy tool catalog")
    parser.add_argument("--url", default=DEFAULT_URL, help=f"Galaxy server URL (default: {DEFAULT_URL})")
    parser.add_argument("--output", default=str(CATALOG_PATH), help="Output path")
    args = parser.parse_args()

    tools = fetch_tools(args.url)
    catalog = build_catalog(tools, args.url)

    out_path = Path(args.output)
    out_path.write_text(
        json.dumps(catalog, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )
    print(f"Wrote {out_path} — {len(tools)} tools, {catalog['section_count']} sections")


if __name__ == "__main__":
    main()
