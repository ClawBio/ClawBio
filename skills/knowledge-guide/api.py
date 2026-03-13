"""Importable API for the knowledge-guide skill."""
from __future__ import annotations
import json
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from clawbio.common.report import write_result_json

RECS_PATH = _SKILL_DIR / "skill_recommendations.json"


def run(genotypes: dict[str, str] | None = None, options: dict | None = None) -> dict:
    """Run a knowledge-guide query.

    Args:
        genotypes: Not used (present for API consistency).
        options: Dict with keys:
            - query: str — free-text query
            - topic: str — GTN topic ID
            - tool: str — Galaxy tool name
            - concept: str — concept term
            - skill: str — ClawBio skill name
            - deep: bool — fetch full tutorial content
            - demo: bool — use embedded demo data
            - output_dir: str — if provided, write report files

    Returns:
        Result dict with summary and data.
    """
    import argparse
    # NOTE: knowledge_guide.py must exist (built in Chunk 2).
    # The argparse.Namespace bridges dict API → CLI internals.
    from knowledge_guide import run_query, generate_report

    options = options or {}

    args = argparse.Namespace(
        query=options.get("query"),
        topic=options.get("topic"),
        tool=options.get("tool"),
        concept=options.get("concept"),
        skill=options.get("skill"),
        deep=options.get("deep", False),
        demo=options.get("demo", False),
    )

    result = run_query(args)

    output_dir = options.get("output_dir")
    if output_dir:
        generate_report(result, Path(output_dir))

    return {
        "summary": {
            "query": result["query"],
            "mode": result["mode"],
            "matches": result["matches"],
        },
        "data": {
            "tutorials": result["tutorials"],
            "deep_content_fetched": result["deep_content_fetched"],
        },
    }


def get_learn_more(
    skill_name: str,
    recommendations_path: Path | None = None,
) -> dict:
    """Get a pre-built 'Learn More' section for a ClawBio skill.

    Reads from skill_recommendations.json (zero network calls).
    Returns structured data + pre-rendered HTML block.

    Args:
        skill_name: ClawBio skill name (e.g. 'gwas-prs').
        recommendations_path: Override path to recommendations JSON.

    Returns:
        Dict with:
            section_title: str
            concepts: list[str]
            tutorials: list[dict]
            html: str (pre-rendered HTML block, empty if no tutorials)
    """
    recs_path = recommendations_path or RECS_PATH

    if recs_path.exists():
        recs = json.loads(recs_path.read_text())
    else:
        recs = {}

    entry = recs.get(skill_name, {})
    concepts = entry.get("concepts", [])
    tutorials = entry.get("tutorials", [])

    # Pre-render HTML block (collapsed by default)
    html = ""
    if tutorials:
        items = []
        for tut in tutorials:
            url = tut.get("url", "#")
            title = tut.get("title", "Tutorial")
            time_est = tut.get("time", tut.get("time_estimation", ""))
            level = tut.get("level", "")
            meta = " | ".join(filter(None, [time_est, level]))
            items.append(
                f'<li><a href="{url}" target="_blank">{title}</a>'
                f'{f" <small>({meta})</small>" if meta else ""}</li>'
            )
        html = (
            '<details class="learn-more">'
            "<summary><strong>Learn More</strong></summary>"
            f'<ul>{"".join(items)}</ul>'
            "</details>"
        )

    return {
        "section_title": "Learn More",
        "concepts": concepts,
        "tutorials": tutorials,
        "html": html,
    }
