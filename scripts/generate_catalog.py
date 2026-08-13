#!/usr/bin/env python3
"""
generate_catalog.py — Build skills/catalog.json from SKILL.md + clawbio/cli.py
==========================================================================
Parses YAML frontmatter from each skill's SKILL.md and cross-references the
SKILLS dict in clawbio/cli.py to produce a machine-readable skill index.

Usage:
    python scripts/generate_catalog.py
"""

from __future__ import annotations

import ast
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
CI_WORKFLOW_PATH = CLAWBIO_DIR / ".github" / "workflows" / "ci.yml"

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

        for field in (
            "version", "author", "domain", "license",
            "model_license", "data_license",
        ):
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
        # Declared separately from `license` on purpose. A permissive
        # wrapper licence grants no rights over third-party weights or
        # reference data, which are frequently non-commercial. Empty means
        # "no such dependency declared", never "permissive".
        "model_license": raw.get("model_license", ""),
        "data_license": raw.get("data_license", ""),
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
# Gather registered skills from clawbio/cli.py SKILLS dict
# ---------------------------------------------------------------------------


CLAWBIO_CLI_PATH = CLAWBIO_DIR / "clawbio" / "cli.py"


def _extract_path_parts(node: ast.AST) -> list[str]:
    """Extract string path fragments from Path-style ``/`` expressions."""
    if isinstance(node, ast.BinOp) and isinstance(node.op, ast.Div):
        return _extract_path_parts(node.left) + _extract_path_parts(node.right)
    if isinstance(node, ast.Constant) and isinstance(node.value, str):
        return [node.value]
    return []


def _iter_static_skill_registry_entries() -> list[tuple[str, str | None]]:
    """Return ``(alias, folder)`` pairs from the static CLI ``SKILLS`` dict.

    The CLI module has optional descriptor augmentation and imports runtime
    helpers, so catalog generation reads the source AST instead of importing it.
    """
    source = CLAWBIO_CLI_PATH.read_text(encoding="utf-8")
    module = ast.parse(source, filename=str(CLAWBIO_CLI_PATH))
    entries: list[tuple[str, str | None]] = []

    for node in module.body:
        if not isinstance(node, ast.Assign):
            continue
        if not any(isinstance(target, ast.Name) and target.id == "SKILLS" for target in node.targets):
            continue
        if not isinstance(node.value, ast.Dict):
            continue

        for key_node, value_node in zip(node.value.keys, node.value.values):
            if not isinstance(key_node, ast.Constant) or not isinstance(key_node.value, str):
                continue
            alias = key_node.value
            folder: str | None = None
            if isinstance(value_node, ast.Dict):
                for entry_key, entry_value in zip(value_node.keys, value_node.values):
                    if isinstance(entry_key, ast.Constant) and entry_key.value == "script":
                        path_parts = _extract_path_parts(entry_value)
                        if path_parts:
                            folder = path_parts[0]
                        break
            entries.append((alias, folder))
        break

    return entries


def load_skills_registry() -> set[str]:
    """Parse CLI aliases from the package CLI ``SKILLS`` dict."""
    return {alias for alias, _folder in _iter_static_skill_registry_entries()}


# ---------------------------------------------------------------------------
# Determine skill folder → CLI alias mapping
# ---------------------------------------------------------------------------

# Aliases that cannot be inferred from CLI script paths.
FOLDER_TO_ALIAS_OVERRIDES = {
    # The photo workflow routes through the PGx reporter after image
    # identification, so its CLI alias points at a different script folder.
    "drug-photo": "drugphoto",
}


def load_folder_to_alias() -> dict[str, str]:
    """Map skill folder names to their primary CLI aliases."""
    folder_to_alias: dict[str, str] = {}
    for alias, folder in _iter_static_skill_registry_entries():
        if folder:
            folder_to_alias.setdefault(folder, alias)
    folder_to_alias.update(FOLDER_TO_ALIAS_OVERRIDES)
    return folder_to_alias


def load_ci_tested_skill_folders() -> set[str]:
    """Return skill folders with explicit pytest coverage in the CI workflow."""
    if not CI_WORKFLOW_PATH.exists():
        return set()
    source = CI_WORKFLOW_PATH.read_text(encoding="utf-8")
    return set(re.findall(r"pytest\s+skills/([^/\s]+)/tests(?:/|\s)", source))


def compute_maturity(
    *,
    has_script: bool,
    has_tests: bool,
    has_demo: bool,
    cli_registered: bool,
    ci_tested: bool,
    benchmark_validated: bool = False,
) -> tuple[str, dict[str, bool]]:
    """Compute an objective maturity tier from observable repository evidence."""
    evidence = {
        "has_skill_md": True,
        "has_script": has_script,
        "has_tests": has_tests,
        "has_demo": has_demo,
        "cli_registered": cli_registered,
        "ci_tested": ci_tested,
        "benchmark_validated": benchmark_validated,
    }

    if benchmark_validated:
        tier = "bench-validated"
    elif ci_tested:
        tier = "ci-validated"
    elif cli_registered:
        tier = "cli-registered"
    elif has_tests:
        tier = "tested"
    elif has_script:
        tier = "scripted"
    else:
        tier = "spec-only"

    return tier, evidence


def select_demo_script(skill_dir: Path, folder_name: str) -> Path | None:
    """Choose a deterministic fallback demo script for non-CLI catalog entries."""
    scripts = sorted(
        f
        for f in skill_dir.glob("*.py")
        if f.name not in {"__init__.py", "api.py"} and not f.name.startswith("test_")
    )
    if not scripts:
        return None

    preferred_names = (
        f"{folder_name.replace('-', '_')}.py",
        "cli.py",
        "__main__.py",
    )
    by_name = {script.name: script for script in scripts}
    for preferred_name in preferred_names:
        if preferred_name in by_name:
            return by_name[preferred_name]

    return scripts[0]

# Skill folders excluded from the public catalog (local-only / gitignored)
EXCLUDED_FOLDERS = {"pr-audit", "wes-clinical-report-es"}

# Skills that are MVP (have working Python + are in SKILLS dict or are bio-orchestrator)
MVP_FOLDERS = {
    "pharmgx-reporter", "equity-scorer", "nutrigx", "claw-metagenomics",
    "nfcore-scrnaseq-wrapper", "nfcore-rnaseq-wrapper", "nfcore-sarek-wrapper", "scrna-orchestrator", "scrna-embedding",
    "genome-compare", "drug-photo", "gwas-prs", "just-prs-mcp", "clinpgx", "gwas-lookup",
    "bigquery-public",
    "profile-report", "bio-orchestrator", "claw-ancestry-pca", "claw-semantic-sim",
    "ukb-navigator", "galaxy-bridge", "rnaseq-de", "diff-visualizer",
    "bioconductor-bridge",
    "sample-qc-triage",
    "crispr-screen-triage",
    "marker-dominance-mapper",
    "llm-biobank-bench",
    "analyze-fasta",
    "phylogenetics-builder",
}

# Known trigger keywords for orchestrator routing
TRIGGER_KEYWORDS: dict[str, list[str]] = {
    "pharmgx-reporter": ["pharmacogenomics", "drug interactions", "23andMe medications", "CYP2D6", "CYP2C19", "warfarin", "CPIC"],
    "drug-photo": ["drug photo", "medication photo", "pill photo", "drug image"],
    "clinpgx": ["ClinPGx", "gene-drug", "PharmGKB", "CPIC guideline database", "FDA drug label"],
    "gwas-lookup": ["GWAS", "variant lookup", "rsID", "PheWAS", "eQTL"],
    "bigquery-public": ["bigquery", "public dataset", "sql", "public data", "cloud query"],
    "gwas-prs": ["polygenic risk", "PRS", "PGS Catalog", "risk score"],
    "just-prs-mcp": ["VCF PRS", "WGS polygenic risk", "evidence-aware PRS", "just-prs"],
    "profile-report": ["profile report", "unified report", "my profile", "genomic profile"],
    "genome-compare": ["genome comparison", "IBS", "George Church", "Corpasome", "pairwise"],
    "equity-scorer": ["HEIM", "equity", "FST", "heterozygosity", "population representation"],
    "nutrigx": ["nutrition", "nutrigenomics", "diet genetics", "MTHFR", "caffeine", "lactose"],
    "nfcore-scrnaseq-wrapper": ["scrnaseq", "nf-core scrnaseq", "single-cell preprocessing", "10x fastq", "generate h5ad from fastq"],
    "nfcore-rnaseq-wrapper": ["bulk RNA-seq preprocessing", "nf-core rnaseq", "run rnaseq from fastq", "FASTQ to count matrix", "STAR Salmon RNA-seq pipeline"],
    "nfcore-sarek-wrapper": ["sarek", "germline variant calling", "somatic variant calling", "tumor-normal pair", "mutect2", "strelka", "haplotypecaller", "ascat", "WES variant calling", "WGS variant calling", "VEP annotation", "nf-core sarek", "tumor-only variant calling", "ControlFREEC", "Manta"],
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
    "sample-qc-triage": ["sample QC triage", "sample identity", "sex mismatch", "fingerprint concordance", "contamination", "batch shift", "low complexity"],
    "crispr-screen-triage": ["CRISPR screen", "guide counts", "rank CRISPR hits", "depleted genes", "knockout screen", "hit triage"],
    "marker-dominance-mapper": ["marker dominance", "map marker spots", "marker-based tissue regions", "tumor core", "immune edge"],
    "analyze-fasta": ["fasta", "analyze fasta", "gc content", "find orfs", "isoelectric point", "gravy index", "protein properties"],
    "phylogenetics-builder": ["phylogeny", "phylogenetic tree", "iqtree", "maximum likelihood tree", "fasta alignment"],
}

# Known chaining partners
CHAINING: dict[str, list[str]] = {
    "pharmgx-reporter": ["drug-photo", "profile-report", "clinpgx"],
    "drug-photo": ["pharmgx-reporter"],
    "clinpgx": ["pharmgx-reporter", "gwas-lookup"],
    "gwas-lookup": ["clinpgx", "gwas-prs", "lit-synthesizer"],
    "bigquery-public": [],
    "gwas-prs": ["profile-report", "gwas-lookup"],
    "just-prs-mcp": ["wgs-prs", "gwas-prs", "profile-report"],
    "profile-report": ["pharmgx-reporter", "nutrigx", "gwas-prs", "genome-compare"],
    "genome-compare": ["claw-ancestry-pca", "profile-report"],
    "equity-scorer": ["claw-semantic-sim"],
    "nutrigx": ["profile-report", "pharmgx-reporter"],
    "nfcore-scrnaseq-wrapper": ["scrna-orchestrator", "scrna-embedding", "bio-orchestrator"],
    "nfcore-rnaseq-wrapper": ["rnaseq-de", "diff-visualizer", "bio-orchestrator"],
    "nfcore-sarek-wrapper": ["variant-annotation", "clinical-variant-reporter", "bio-orchestrator"],
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
    "sample-qc-triage": ["multiqc-reporter", "seq-wrangler"],
    "crispr-screen-triage": ["target-validation-scorer", "omics-target-evidence-mapper"],
    "marker-dominance-mapper": ["scrna-orchestrator", "diff-visualizer"],
    "analyze-fasta": ["struct-predictor", "variant-annotation", "pubmed-summariser"],
    "phylogenetics-builder": ["profile-report"],
}


# ---------------------------------------------------------------------------
# Build catalog
# ---------------------------------------------------------------------------


def build_catalog() -> list[dict]:
    """Build a list of skill entries for the catalog."""
    registered_aliases = load_skills_registry()
    folder_to_alias = load_folder_to_alias()
    ci_tested_folders = load_ci_tested_skill_folders()
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
        cli_alias = folder_to_alias.get(folder_name)

        # Check for Python scripts and tests
        has_script = any(
            f.suffix == ".py" and f.name != "__init__.py" and "test" not in f.name.lower()
            for f in sorted(skill_dir.rglob("*.py"))
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
            script = select_demo_script(skill_dir, folder_name)
            if script is not None:
                demo_command = f"python {script.relative_to(CLAWBIO_DIR)} --demo"
        has_demo = demo_command is not None
        cli_registered = bool(cli_alias and cli_alias in registered_aliases)
        ci_tested = folder_name in ci_tested_folders

        # Status
        status = "mvp" if folder_name in MVP_FOLDERS else "planned"
        maturity_tier, maturity_evidence = compute_maturity(
            has_script=has_script,
            has_tests=has_tests,
            has_demo=has_demo,
            cli_registered=cli_registered,
            ci_tested=ci_tested,
        )

        tags = [str(tag) for tag in skill_meta.get("tags", [])]
        deps = _normalize_dependencies(skill_meta.get("dependencies"))
        trigger_keywords = yaml_data.get("trigger_keywords") or skill_meta.get("openclaw", {}).get("trigger_keywords") or TRIGGER_KEYWORDS.get(folder_name, [])

        entry = {
            "name": folder_name,
            "cli_alias": cli_alias,
            "description": skill_meta.get("description", ""),
            "version": str(skill_meta.get("version", "0.1.0")),
            "status": status,
            "maturity_tier": maturity_tier,
            "maturity_evidence": maturity_evidence,
            "has_script": has_script,
            "has_tests": has_tests,
            "has_demo": has_demo,
            "demo_command": demo_command,
            "dependencies": deps,
            "tags": tags,
            "trigger_keywords": trigger_keywords,
            "chaining_partners": CHAINING.get(folder_name, []),
            # Three separate fields, because they restrict independently and
            # the most restrictive governs. `license` covers the skill's own
            # code; the other two cover third-party weights and reference data
            # the skill pulls at runtime, over which no licence in this repo
            # can grant rights. Empty means "none declared", not "permissive".
            "license": str(skill_meta.get("license", "")),
            "model_license": str(skill_meta.get("model_license", "")),
            "data_license": str(skill_meta.get("data_license", "")),
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
