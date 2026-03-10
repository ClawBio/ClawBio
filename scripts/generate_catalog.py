#!/usr/bin/env python3
"""
generate_catalog.py — Build skills/catalog.json from SKILL.md + clawbio.py
==========================================================================
Parses YAML frontmatter from each skill's SKILL.md and cross-references the
SKILLS dict in clawbio.py to produce a machine-readable skill index.

Usage:
    python scripts/generate_catalog.py
"""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

CLAWBIO_DIR = Path(__file__).resolve().parents[1]
SKILLS_DIR = CLAWBIO_DIR / "skills"
CATALOG_PATH = SKILLS_DIR / "catalog.json"

sys.path.insert(0, str(CLAWBIO_DIR))

# ---------------------------------------------------------------------------
# YAML frontmatter parser (lightweight, no PyYAML dependency)
# ---------------------------------------------------------------------------


def parse_yaml_frontmatter(text: str) -> dict:
    """Extract YAML frontmatter between --- markers. Returns flat-ish dict."""
    match = re.match(r"^---\s*\n(.*?)\n---", text, re.DOTALL)
    if not match:
        return {}
    raw = match.group(1)
    result: dict = {}
    current_key = None
    for line in raw.split("\n"):
        # Top-level key: value
        m = re.match(r"^(\w[\w-]*):\s*(.*)", line)
        if m:
            key, val = m.group(1), m.group(2).strip()
            if val.startswith("[") and val.endswith("]"):
                result[key] = [v.strip().strip("'\"") for v in val[1:-1].split(",") if v.strip()]
            elif val == "" or val == "|":
                result[key] = ""
                current_key = key
            else:
                result[key] = val.strip("'\"")
                current_key = key
            continue
        # List item under current key
        m2 = re.match(r"^\s+-\s+(.*)", line)
        if m2 and current_key:
            if not isinstance(result.get(current_key), list):
                result[current_key] = []
            result[current_key].append(m2.group(1).strip().strip("'\""))
    return result


# ---------------------------------------------------------------------------
# Gather registered skills from clawbio.py SKILLS dict
# ---------------------------------------------------------------------------


def load_skills_registry() -> set:
    """Parse CLI aliases from clawbio.py SKILLS dict keys without importing.

    clawbio.py requires Python 3.10+ (``str | None`` syntax), so we parse
    the SKILLS dict keys with a regex instead of importing the module.
    """
    source = (CLAWBIO_DIR / "clawbio.py").read_text(encoding="utf-8")
    # Match top-level keys like:  "pharmgx": {
    return set(re.findall(r'^\s{4}"(\w+)":\s*\{', source, re.MULTILINE))


# ---------------------------------------------------------------------------
# Determine skill folder → CLI alias mapping
# ---------------------------------------------------------------------------

# Map skill folder names to their CLI alias in SKILLS dict.
# The SKILLS dict keys are short aliases; script paths reveal the folder name.
FOLDER_TO_ALIAS = {
    "pharmgx-reporter": "pharmgx",
    "equity-scorer": "equity",
    "nutrigx_advisor": "nutrigx",
    "scrna-orchestrator": "scrna",
    "claw-metagenomics": "metagenomics",
    "genome-compare": "compare",
    "drug-photo": "drugphoto",
    "gwas-prs": "prs",
    "clinpgx": "clinpgx",
    "gwas-lookup": "gwas",
    "profile-report": "profile",
    "galaxy-bridge": "galaxy",
}

# Skill folders excluded from the public catalog (local-only / gitignored)
EXCLUDED_FOLDERS = {"pr-audit"}

# Skills that are MVP (have working Python + are in SKILLS dict or are bio-orchestrator)
MVP_FOLDERS = {
    "pharmgx-reporter", "equity-scorer", "nutrigx_advisor", "claw-metagenomics",
    "scrna-orchestrator",
    "genome-compare", "drug-photo", "gwas-prs", "clinpgx", "gwas-lookup",
    "profile-report", "bio-orchestrator", "claw-ancestry-pca", "claw-semantic-sim",
    "ukb-navigator", "galaxy-bridge",
}

# Known trigger keywords for orchestrator routing
TRIGGER_KEYWORDS: dict[str, list[str]] = {
    "pharmgx-reporter": ["pharmacogenomics", "drug interactions", "23andMe medications", "CYP2D6", "CYP2C19", "warfarin", "CPIC"],
    "drug-photo": ["drug photo", "medication photo", "pill photo", "drug image"],
    "clinpgx": ["ClinPGx", "gene-drug", "PharmGKB", "CPIC guideline database", "FDA drug label"],
    "gwas-lookup": ["GWAS", "variant lookup", "rsID", "PheWAS", "eQTL"],
    "gwas-prs": ["polygenic risk", "PRS", "PGS Catalog", "risk score"],
    "profile-report": ["profile report", "unified report", "my profile", "genomic profile"],
    "genome-compare": ["genome comparison", "IBS", "George Church", "Corpasome", "pairwise"],
    "equity-scorer": ["HEIM", "equity", "FST", "heterozygosity", "population representation"],
    "nutrigx_advisor": ["nutrition", "nutrigenomics", "diet genetics", "MTHFR", "caffeine", "lactose"],
    "scrna-orchestrator": ["single-cell", "scrna", "h5ad", "mtx", "10x", "scanpy", "umap", "leiden"],
    "claw-ancestry-pca": ["ancestry", "PCA", "admixture", "SGDP", "population structure"],
    "claw-semantic-sim": ["semantic similarity", "disease neglect", "research gaps", "NTDs", "SII"],
    "claw-metagenomics": ["metagenomics", "Kraken2", "RGI", "CARD", "HUMAnN3", "microbiome"],
    "bio-orchestrator": ["route", "which skill", "orchestrator"],
    "ukb-navigator": ["UK Biobank", "UKB", "biobank schema", "data showcase"],
    "galaxy-bridge": ["galaxy", "usegalaxy", "tool shed", "bioblend", "run on galaxy", "galaxy tool", "galaxy workflow", "NGS pipeline"],
}

# Known chaining partners
CHAINING: dict[str, list[str]] = {
    "pharmgx-reporter": ["drug-photo", "profile-report", "clinpgx"],
    "drug-photo": ["pharmgx-reporter"],
    "clinpgx": ["pharmgx-reporter", "gwas-lookup"],
    "gwas-lookup": ["clinpgx", "gwas-prs", "lit-synthesizer"],
    "gwas-prs": ["profile-report", "gwas-lookup"],
    "profile-report": ["pharmgx-reporter", "nutrigx_advisor", "gwas-prs", "genome-compare"],
    "genome-compare": ["claw-ancestry-pca", "profile-report"],
    "equity-scorer": ["claw-semantic-sim"],
    "nutrigx_advisor": ["profile-report", "pharmgx-reporter"],
    "scrna-orchestrator": [],
    "claw-ancestry-pca": ["genome-compare"],
    "claw-semantic-sim": ["equity-scorer"],
    "claw-metagenomics": [],
    "bio-orchestrator": [],
    "ukb-navigator": [],
    "galaxy-bridge": ["pharmgx-reporter", "claw-metagenomics", "equity-scorer", "vcf-annotator"],
}


# ---------------------------------------------------------------------------
# Build catalog
# ---------------------------------------------------------------------------


def build_catalog() -> list[dict]:
    """Build a list of skill entries for the catalog."""
    registered_aliases = load_skills_registry()
    entries: list[dict] = []

    for skill_dir in sorted(SKILLS_DIR.iterdir()):
        if not skill_dir.is_dir() or skill_dir.name.startswith("."):
            continue
        if skill_dir.name in EXCLUDED_FOLDERS:
            continue

        skill_md = skill_dir / "SKILL.md"
        if not skill_md.exists():
            continue

        folder_name = skill_dir.name
        yaml_data = parse_yaml_frontmatter(skill_md.read_text(encoding="utf-8"))

        # Determine CLI alias
        cli_alias = FOLDER_TO_ALIAS.get(folder_name)

        # Check for Python scripts and tests
        has_script = any(
            f.suffix == ".py" and f.name != "__init__.py" and "test" not in f.name.lower()
            for f in skill_dir.rglob("*.py")
            if "tests" not in str(f.relative_to(skill_dir)).split("/")[0:1]
            and "__pycache__" not in str(f)
        )
        tests_dir = skill_dir / "tests"
        has_tests = tests_dir.exists() and any(tests_dir.glob("test_*.py"))

        # Demo command
        demo_command = None
        if cli_alias and cli_alias in registered_aliases:
            demo_command = f"python clawbio.py run {cli_alias} --demo"
        elif has_script:
            scripts = [f for f in skill_dir.glob("*.py") if f.name != "__init__.py" and f.name != "api.py"]
            if scripts:
                demo_command = f"python {scripts[0].relative_to(CLAWBIO_DIR)} --demo"

        # Status
        status = "mvp" if folder_name in MVP_FOLDERS else "planned"

        # Input types from YAML
        input_types = []
        if isinstance(yaml_data.get("inputs"), list):
            input_types = yaml_data["inputs"] if isinstance(yaml_data["inputs"][0], str) else []

        # Tags
        tags = yaml_data.get("tags", [])
        if isinstance(tags, str):
            tags = [t.strip() for t in tags.split(",")]

        # Dependencies from YAML
        deps = yaml_data.get("dependencies", [])
        if isinstance(deps, str):
            deps = [deps] if deps else []

        entry = {
            "name": folder_name,
            "cli_alias": cli_alias,
            "description": yaml_data.get("description", ""),
            "version": yaml_data.get("version", "0.1.0"),
            "status": status,
            "has_script": has_script,
            "has_tests": has_tests,
            "has_demo": demo_command is not None,
            "demo_command": demo_command,
            "dependencies": deps,
            "tags": tags,
            "trigger_keywords": TRIGGER_KEYWORDS.get(folder_name, []),
            "chaining_partners": CHAINING.get(folder_name, []),
        }
        entries.append(entry)

    return entries


def main() -> None:
    catalog = build_catalog()

    # Inject Galaxy tool count from galaxy_catalog.json if present
    galaxy_tool_count = 0
    galaxy_catalog_path = SKILLS_DIR / "galaxy-bridge" / "galaxy_catalog.json"
    if galaxy_catalog_path.exists():
        try:
            gcat = json.loads(galaxy_catalog_path.read_text(encoding="utf-8"))
            galaxy_tool_count = gcat.get("tool_count", 0)
        except (json.JSONDecodeError, KeyError):
            pass

    catalog_obj = {
        "version": "1.0.0",
        "generated_by": "scripts/generate_catalog.py",
        "skill_count": len(catalog),
        "galaxy_tool_count": galaxy_tool_count,
        "skills": catalog,
    }

    CATALOG_PATH.write_text(
        json.dumps(catalog_obj, indent=2, ensure_ascii=False) + "\n",
        encoding="utf-8",
    )

    mvp = sum(1 for s in catalog if s["status"] == "mvp")
    planned = sum(1 for s in catalog if s["status"] == "planned")
    print(f"Wrote {CATALOG_PATH.relative_to(CLAWBIO_DIR)} — {len(catalog)} skills ({mvp} MVP, {planned} planned)")


if __name__ == "__main__":
    main()
