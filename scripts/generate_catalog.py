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
    """Extract YAML frontmatter between --- markers."""
    match = re.match(r"^---\s*\n(.*?)\n---", text, re.DOTALL)
    if not match:
        return {}
    raw = match.group(1)
    try:
        import yaml

        data = yaml.safe_load(raw) or {}
        return data if isinstance(data, dict) else {}
    except ImportError:
        result: dict = {}

        def _strip(value: str) -> str:
            return value.strip().strip('"').strip("'")

        name_match = re.search(r"^name:\s*(.+)", raw, re.MULTILINE)
        if name_match:
            result["name"] = _strip(name_match.group(1))

        desc_match = re.search(r"^description:\s*(.+)", raw, re.MULTILINE)
        if desc_match:
            result["description"] = _strip(desc_match.group(1))

        for field in ("version", "author", "domain", "license"):
            match = re.search(rf"^{field}:\s*(.+)", raw, re.MULTILINE)
            if match:
                result[field] = _strip(match.group(1))

        tags_match = re.search(r"^tags:\s*\[([^\]]*)\]", raw, re.MULTILINE)
        if tags_match:
            result["tags"] = [_strip(v) for v in tags_match.group(1).split(",") if v.strip()]

        dep_python_match = re.search(r"^dependencies:\s*\n(?:\s+.+\n)*?\s+python:\s*(.+)", raw, re.MULTILINE)
        dep_packages_match = re.search(r"^dependencies:\s*\n(?:\s+.+\n)*?\s+packages:\s*\n((?:\s+-\s+.+\n)*)", raw, re.MULTILINE)
        if dep_python_match or dep_packages_match:
            deps: dict = {}
            if dep_python_match:
                deps["python"] = _strip(dep_python_match.group(1))
            if dep_packages_match:
                deps["packages"] = [
                    _strip(line.strip().lstrip("- "))
                    for line in dep_packages_match.group(1).splitlines()
                    if line.strip()
                ]
            result["dependencies"] = deps

        metadata_match = re.search(r"^metadata:\s*\n((?:\s+.+\n)*)", raw, re.MULTILINE)
        if metadata_match:
            metadata_block = metadata_match.group(1)
            metadata: dict = {}
            for field in ("version", "author", "domain"):
                match = re.search(rf"^\s+{field}:\s*(.+)", metadata_block, re.MULTILINE)
                if match:
                    metadata[field] = _strip(match.group(1))

            tags_block = re.search(r"^\s+tags:\s*\n((?:\s+-\s+.+\n)*)", metadata_block, re.MULTILINE)
            if tags_block:
                metadata["tags"] = [
                    _strip(line.strip().lstrip("- "))
                    for line in tags_block.group(1).splitlines()
                    if line.strip()
                ]

            dep_meta_python = re.search(r"^\s+dependencies:\s*\n(?:\s+.+\n)*?\s+python:\s*(.+)", metadata_block, re.MULTILINE)
            dep_meta_packages = re.search(r"^\s+dependencies:\s*\n(?:\s+.+\n)*?\s+packages:\s*\n((?:\s+-\s+.+\n)*)", metadata_block, re.MULTILINE)
            if dep_meta_python or dep_meta_packages:
                deps: dict = {}
                if dep_meta_python:
                    deps["python"] = _strip(dep_meta_python.group(1))
                if dep_meta_packages:
                    deps["packages"] = [
                        _strip(line.strip().lstrip("- "))
                        for line in dep_meta_packages.group(1).splitlines()
                        if line.strip()
                    ]
                metadata["dependencies"] = deps

            openclaw_block = re.search(r"^\s+openclaw:\s*\n((?:\s+.+\n)*)", metadata_block, re.MULTILINE)
            if openclaw_block:
                oc: dict = {}
                trigger_block = re.search(r"trigger_keywords:\s*\n((?:\s+-\s+.+\n)*)", openclaw_block.group(1))
                if trigger_block:
                    oc["trigger_keywords"] = [
                        _strip(line.strip().lstrip("- "))
                        for line in trigger_block.group(1).splitlines()
                        if line.strip()
                    ]
                metadata["openclaw"] = oc

            result["metadata"] = metadata

        return result


def normalize_skill_metadata(raw: dict) -> dict:
    """Normalize legacy top-level and AgentSkills nested metadata."""
    metadata = raw.get("metadata") if isinstance(raw.get("metadata"), dict) else {}
    openclaw = metadata.get("openclaw") if isinstance(metadata.get("openclaw"), dict) else {}
    return {
        "name": raw.get("name", ""),
        "description": raw.get("description", ""),
        "license": raw.get("license", ""),
        "version": raw.get("version", metadata.get("version", "0.1.0")),
        "author": raw.get("author", metadata.get("author", "")),
        "domain": raw.get("domain", metadata.get("domain", "")),
        "tags": raw.get("tags", metadata.get("tags", [])),
        "inputs": raw.get("inputs", metadata.get("inputs", [])),
        "outputs": raw.get("outputs", metadata.get("outputs", [])),
        "dependencies": raw.get("dependencies", metadata.get("dependencies", [])),
        "demo_data": raw.get("demo_data", metadata.get("demo_data", [])),
        "endpoints": raw.get("endpoints", metadata.get("endpoints", {})),
        "openclaw": openclaw,
    }


def _normalize_dependencies(value: object) -> list[str]:
    if value is None:
        return []
    if isinstance(value, str):
        return [value] if value else []
    if isinstance(value, list):
        return [str(item) for item in value]
    if isinstance(value, dict):
        deps: list[str] = []
        python_req = value.get("python")
        if python_req:
            deps.append(f"python{python_req}")
        packages = value.get("packages")
        if isinstance(packages, list):
            deps.extend(str(pkg) for pkg in packages)
        elif packages:
            deps.append(str(packages))
        return deps
    return [str(value)]


# ---------------------------------------------------------------------------
# Gather registered skills from clawbio.py SKILLS dict
# ---------------------------------------------------------------------------


def load_skills_registry() -> set:
    """Parse CLI aliases from clawbio.py SKILLS dict keys without importing.

    clawbio.py requires Python 3.10+ (``str | None`` syntax), so we parse
    the SKILLS dict keys with a regex instead of importing the module.
    """
    source = (CLAWBIO_DIR / "clawbio.py").read_text(encoding="utf-8")
    # Match top-level keys like:  "pharmgx": { or "scrna-embedding": {
    return set(re.findall(r'^\s{4}"([\w-]+)":\s*\{', source, re.MULTILINE))


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
    "scrna-embedding": "scrna-embedding",
    "claw-metagenomics": "metagenomics",
    "genome-compare": "compare",
    "drug-photo": "drugphoto",
    "gwas-prs": "prs",
    "clinpgx": "clinpgx",
    "gwas-lookup": "gwas",
    "bigquery-public": "bigquery",
    "profile-report": "profile",
    "galaxy-bridge": "galaxy",
    "bioconductor-bridge": "bioc",
    "rnaseq-de": "rnaseq",
    "diff-visualizer": "diffviz",
    "llm-biobank-bench": "llm-bench",
}

# Skill folders excluded from the public catalog (local-only / gitignored)
EXCLUDED_FOLDERS = {"pr-audit", "wes-clinical-report-es"}

# Skills that are MVP (have working Python + are in SKILLS dict or are bio-orchestrator)
MVP_FOLDERS = {
    "pharmgx-reporter", "equity-scorer", "nutrigx_advisor", "claw-metagenomics",
    "scrna-orchestrator", "scrna-embedding",
    "genome-compare", "drug-photo", "gwas-prs", "clinpgx", "gwas-lookup",
    "bigquery-public",
    "profile-report", "bio-orchestrator", "claw-ancestry-pca", "claw-semantic-sim",
    "ukb-navigator", "galaxy-bridge", "rnaseq-de", "diff-visualizer",
    "bioconductor-bridge",
    "llm-biobank-bench",
}

# Known trigger keywords for orchestrator routing
TRIGGER_KEYWORDS: dict[str, list[str]] = {
    "pharmgx-reporter": ["pharmacogenomics", "drug interactions", "23andMe medications", "CYP2D6", "CYP2C19", "warfarin", "CPIC"],
    "drug-photo": ["drug photo", "medication photo", "pill photo", "drug image"],
    "clinpgx": ["ClinPGx", "gene-drug", "PharmGKB", "CPIC guideline database", "FDA drug label"],
    "gwas-lookup": ["GWAS", "variant lookup", "rsID", "PheWAS", "eQTL"],
    "bigquery-public": ["bigquery", "public dataset", "sql", "public data", "cloud query"],
    "gwas-prs": ["polygenic risk", "PRS", "PGS Catalog", "risk score"],
    "profile-report": ["profile report", "unified report", "my profile", "genomic profile"],
    "genome-compare": ["genome comparison", "IBS", "George Church", "Corpasome", "pairwise"],
    "equity-scorer": ["HEIM", "equity", "FST", "heterozygosity", "population representation"],
    "nutrigx_advisor": ["nutrition", "nutrigenomics", "diet genetics", "MTHFR", "caffeine", "lactose"],
    "scrna-orchestrator": ["single-cell", "scrna", "h5ad", "mtx", "10x", "scanpy", "umap", "leiden"],
    "scrna-embedding": ["scvi", "scanvi", "latent", "embedding", "integration", "batch correction", "10x"],
    "rnaseq-de": ["differential expression", "bulk rna", "rna-seq", "count matrix", "deseq2", "pydeseq2"],
    "diff-visualizer": ["visualize de results", "de visualization", "marker heatmap", "marker dotplot", "top genes heatmap"],
    "claw-ancestry-pca": ["ancestry", "PCA", "admixture", "SGDP", "population structure"],
    "claw-semantic-sim": ["semantic similarity", "disease neglect", "research gaps", "NTDs", "SII"],
    "claw-metagenomics": ["metagenomics", "Kraken2", "RGI", "CARD", "HUMAnN3", "microbiome"],
    "bio-orchestrator": ["route", "which skill", "orchestrator"],
    "ukb-navigator": ["UK Biobank", "UKB", "biobank schema", "data showcase"],
    "llm-biobank-bench": ["llm benchmark", "benchmark language models", "biobank knowledge retrieval", "coverage score", "weighted coverage", "model comparison biobank"],
    "galaxy-bridge": ["galaxy", "usegalaxy", "tool shed", "bioblend", "run on galaxy", "galaxy tool", "galaxy workflow", "NGS pipeline"],
    "bioconductor-bridge": ["bioconductor", "bioc", "biocmanager", "summarizedexperiment", "singlecellexperiment", "genomicranges", "variantannotation", "annotationhub", "experimenthub"],
}

# Known chaining partners
CHAINING: dict[str, list[str]] = {
    "pharmgx-reporter": ["drug-photo", "profile-report", "clinpgx"],
    "drug-photo": ["pharmgx-reporter"],
    "clinpgx": ["pharmgx-reporter", "gwas-lookup"],
    "gwas-lookup": ["clinpgx", "gwas-prs", "lit-synthesizer"],
    "bigquery-public": [],
    "gwas-prs": ["profile-report", "gwas-lookup"],
    "profile-report": ["pharmgx-reporter", "nutrigx_advisor", "gwas-prs", "genome-compare"],
    "genome-compare": ["claw-ancestry-pca", "profile-report"],
    "equity-scorer": ["claw-semantic-sim"],
    "nutrigx_advisor": ["profile-report", "pharmgx-reporter"],
    "scrna-orchestrator": [],
    "scrna-embedding": ["scrna-orchestrator"],
    "rnaseq-de": ["diff-visualizer"],
    "diff-visualizer": ["rnaseq-de", "scrna-orchestrator"],
    "claw-ancestry-pca": ["genome-compare"],
    "claw-semantic-sim": ["equity-scorer"],
    "claw-metagenomics": [],
    "bio-orchestrator": [],
    "ukb-navigator": ["llm-biobank-bench"],
    "llm-biobank-bench": ["ukb-navigator", "pubmed-summariser", "lit-synthesizer"],
    "galaxy-bridge": ["pharmgx-reporter", "claw-metagenomics", "equity-scorer", "vcf-annotator"],
    "bioconductor-bridge": ["rnaseq-de", "scrna-orchestrator", "diff-visualizer", "bio-orchestrator"],
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
        skill_meta = normalize_skill_metadata(yaml_data)

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

        tags = [str(tag) for tag in skill_meta.get("tags", [])]
        deps = _normalize_dependencies(skill_meta.get("dependencies"))
        trigger_keywords = skill_meta.get("openclaw", {}).get("trigger_keywords") or TRIGGER_KEYWORDS.get(folder_name, [])

        entry = {
            "name": folder_name,
            "cli_alias": cli_alias,
            "description": skill_meta.get("description", ""),
            "version": str(skill_meta.get("version", "0.1.0")),
            "status": status,
            "has_script": has_script,
            "has_tests": has_tests,
            "has_demo": demo_command is not None,
            "demo_command": demo_command,
            "dependencies": deps,
            "tags": tags,
            "trigger_keywords": trigger_keywords,
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
