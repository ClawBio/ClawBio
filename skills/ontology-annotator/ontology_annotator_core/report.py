"""report.py — report.md generation, including "which ClawBio skills can
consume these ontology IDs", read live from skills/catalog.json.
"""

from __future__ import annotations

import json
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)

ONTOLOGY_NAMES = {
    "uberon": "UBERON (anatomy/tissue)",
    "cl": "CL (Cell Ontology)",
    "mondo": "MONDO (disease)",
    "efo": "EFO (Experimental Factor Ontology / trait)",
}

# Keywords used to search skills/catalog.json for entries relevant to each
# ontology's domain. This drives *discovery* of candidate consuming skills —
# it never stands in for a hardcoded skill list, and every name printed in
# the report is read out of catalog.json at run time, not written here.
ONTOLOGY_HINT_KEYWORDS = {
    "uberon": ["tissue", "uberon", "organ", "anatomy"],
    "cl": ["cell type", "cell-type", "single-cell", "single cell", " cl ", "cell ontology"],
    "mondo": ["disease", "mondo", "diagnosis"],
    "efo": ["trait", "efo", "phenotype", "experimental factor"],
}


def find_consuming_skills(catalog_path: Path, ontologies: list[str]) -> dict[str, list[dict]]:
    """For each ontology prefix this run produced IDs for, find catalog
    entries whose own name/description/tags/trigger_keywords mention that
    ontology's domain. Returns {} per-ontology if the catalog is missing or
    unreadable rather than raising — a missing catalog should never break
    report generation."""
    result: dict[str, list[dict]] = {o: [] for o in ontologies}
    if not catalog_path.exists():
        return result

    try:
        catalog = json.loads(catalog_path.read_text())
    except json.JSONDecodeError:
        return result

    entries = catalog.get("skills", []) if isinstance(catalog, dict) else catalog
    if not isinstance(entries, list):
        return result

    for ontology in ontologies:
        keywords = ONTOLOGY_HINT_KEYWORDS.get(ontology, [ontology])
        seen_names: set[str] = set()
        for entry in entries:
            if not isinstance(entry, dict):
                continue
            name = entry.get("name", "")
            if name == "ontology-annotator" or name in seen_names:
                continue
            haystack = " ".join(
                [
                    str(entry.get("name", "")),
                    str(entry.get("description", "")),
                    " ".join(str(t) for t in entry.get("tags", []) or []),
                    " ".join(str(t) for t in entry.get("trigger_keywords", []) or []),
                ]
            ).lower()
            if any(kw in haystack for kw in keywords):
                result[ontology].append(
                    {"name": name, "description": entry.get("description", "")}
                )
                seen_names.add(name)
    return result


def generate_markdown(
    annotated_columns: dict[str, str],
    row_summaries: dict[str, dict[str, Any]],
    flagged_rows: list[dict[str, Any]],
    n_rows: int,
    threshold: float,
    catalog_path: Path,
) -> str:
    """Build the full report.md text.

    Args:
        annotated_columns: {column_name: ontology_prefix} for this run.
        row_summaries: {column_name: {"n_matched", "n_flagged", "n_no_candidates"}}.
        flagged_rows: list of {row_index, column, value, reason} for every flagged cell.
        n_rows: total input rows processed.
        threshold: confidence threshold used to decide flags.
        catalog_path: path to skills/catalog.json, for the consuming-skills section.
    """
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    lines: list[str] = [
        "# Ontology Annotator Report",
        "",
        f"**Date**: {now}",
        f"**Rows processed**: {n_rows}",
        f"**Confidence threshold**: {threshold}",
        f"**Data source**: EBI OLS4 (`https://www.ebi.ac.uk/ols4/api/search`)",
        "",
        "---",
        "",
        "## Summary",
        "",
        "| Column | Ontology | Matched (≥ threshold) | Flagged | No candidates |",
        "|--------|----------|------------------------|---------|----------------|",
    ]
    for col, ontology in annotated_columns.items():
        s = row_summaries.get(col, {})
        lines.append(
            f"| `{col}` | {ONTOLOGY_NAMES.get(ontology, ontology.upper())} "
            f"| {s.get('n_matched', 0)} | {s.get('n_flagged', 0)} | {s.get('n_no_candidates', 0)} |"
        )
    lines.append("")

    lines.append("## Flagged Rows")
    lines.append("")
    if flagged_rows:
        lines.append(
            "These rows scored below the confidence threshold, or OLS4 returned no "
            "candidates at all. **None of these were auto-picked** — review the "
            "top-3 candidates in `annotated.csv` before trusting them."
        )
        lines.append("")
        lines.append("| Row | Column | Value | Top candidate | Score | Reason |")
        lines.append("|-----|--------|-------|----------------|-------|--------|")
        for f in flagged_rows:
            lines.append(
                f"| {f['row_index']} | `{f['column']}` | {f['value']} | "
                f"{f.get('top_label', '—')} | {f.get('top_score', '—')} | {f['reason']} |"
            )
    else:
        lines.append("None — every annotated value matched with confidence ≥ threshold.")
    lines.append("")

    lines.append("## Ontologies Used")
    lines.append("")
    for ontology in sorted(set(annotated_columns.values())):
        lines.append(f"- **{ONTOLOGY_NAMES.get(ontology, ontology.upper())}**")
    lines.append("")

    lines.append("## ClawBio Skills That Can Consume These IDs")
    lines.append("")
    lines.append(
        "Read live from `skills/catalog.json` by matching each ontology's domain "
        "keywords against every skill's own name/description/tags/trigger_keywords — "
        "not a hardcoded list, so it stays correct as the catalog grows."
    )
    lines.append("")
    consuming = find_consuming_skills(catalog_path, sorted(set(annotated_columns.values())))
    for ontology, skills_list in consuming.items():
        lines.append(f"### {ONTOLOGY_NAMES.get(ontology, ontology.upper())}")
        lines.append("")
        if skills_list:
            for s in skills_list:
                lines.append(f"- **{s['name']}** — {s['description']}")
        else:
            lines.append("- No catalog entries matched this ontology's domain keywords.")
        lines.append("")

    lines.append("---")
    lines.append("")
    lines.append("## Methods")
    lines.append("")
    lines.append(
        "- Each unique (column, ontology) value is looked up once against OLS4 "
        "`/api/search` and cached; identical or differently-cased repeats reuse "
        "the cached response."
    )
    lines.append(
        "- Candidates are re-ranked locally by `difflib.SequenceMatcher` string "
        "similarity between the input value and each candidate's label/synonyms "
        "(deterministic, reproducible by hand)."
    )
    lines.append(
        f"- A row is flagged when its top-scoring candidate is below the "
        f"threshold ({threshold}), or when OLS4 returned zero candidates. "
        "Flagged rows keep their top-3 candidates in the output — nothing is "
        "silently picked."
    )
    lines.append("")
    lines.append(f"*{DISCLAIMER}*")
    lines.append("")

    return "\n".join(lines)
