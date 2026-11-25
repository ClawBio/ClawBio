"""Robotary — Robot Terry web interface for ClawBio skills.

Run:
    python robotary/server.py
    # → http://localhost:5112
"""
from __future__ import annotations

from pathlib import Path

from fastapi import FastAPI
from fastapi.responses import HTMLResponse
from fastapi.staticfiles import StaticFiles

PROJECT_ROOT = Path(__file__).resolve().parent.parent
SKILLS_DIR = PROJECT_ROOT / "skills"
ROBOTARY_DIR = Path(__file__).resolve().parent

# Manuel Corpas genome — used as default input for genome-based skills
GENOME_PATH = SKILLS_DIR / "genome-compare" / "data" / "manuel_corpas_23andme.txt.gz"

# Maps skill directory names → clawbio.py registry names (core skills only)
SKILL_REGISTRY_MAP: dict[str, str] = {
    "pharmgx-reporter": "pharmgx",
    "nutrigx": "nutrigx",
    "genome-compare": "compare",
    "gwas-prs": "prs",
    "clinpgx": "clinpgx",
    "gwas-lookup": "gwas",
    "profile-report": "profile",
}

# Skills that accept a genome file as --input (vs --demo)
GENOME_SKILLS = {"pharmgx-reporter", "nutrigx", "genome-compare", "gwas-prs"}

# Core skills exposed in Robotary (subset of all ClawBio skills)
CORE_SKILLS = {
    "pharmgx-reporter", "nutrigx", "genome-compare",
    "gwas-prs", "clinpgx", "gwas-lookup", "profile-report",
}


def build_skill_catalog() -> dict[str, str]:
    """Discover skills by scanning skills/*/SKILL.md. Returns {name: description}."""
    catalog = {}
    for skill_dir in sorted(SKILLS_DIR.iterdir()):
        if not skill_dir.is_dir() or skill_dir.name not in CORE_SKILLS:
            continue
        skill_md = skill_dir / "SKILL.md"
        if not skill_md.exists():
            continue
        lines = skill_md.read_text().split("\n")
        desc_lines = []
        in_frontmatter = False
        for line in lines:
            if line.strip() == "---":
                in_frontmatter = not in_frontmatter
                continue
            if not in_frontmatter and line.strip():
                desc_lines.append(line.strip())
                if len(desc_lines) >= 3:
                    break
        catalog[skill_dir.name] = " ".join(desc_lines)[:200]
    return catalog


SKILL_CATALOG = build_skill_catalog()

app = FastAPI(title="Robotary")

# Serve static assets if directory exists
_static_dir = ROBOTARY_DIR / "static"
if _static_dir.is_dir():
    app.mount("/static", StaticFiles(directory=str(_static_dir)), name="static")


@app.get("/")
async def index():
    html_path = ROBOTARY_DIR / "index.html"
    return HTMLResponse(html_path.read_text())


if __name__ == "__main__":
    import uvicorn
    print("=" * 50)
    print("  Robotary — Robot Terry")
    print(f"  Skills loaded: {len(SKILL_CATALOG)}")
    print(f"  Genome: {GENOME_PATH.name}")
    print("  Open http://localhost:5112")
    print("=" * 50)
    uvicorn.run(app, host="0.0.0.0", port=5112)
