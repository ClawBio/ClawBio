#!/usr/bin/env python3
"""Bio Orchestrator: routes bioinformatics requests to specialised skills.

Usage:
    python orchestrator.py --input <file_or_query> [--skill <skill_name>] [--output <dir>]
    python orchestrator.py --profile <profile.json> --skills pharmgx,nutrigx --output <dir>

This is the supporting Python code for the Bio Orchestrator skill.
It handles file type detection, skill routing, multi-skill dispatch,
and report assembly.
"""

from __future__ import annotations

from pathlib import Path
import argparse
import hashlib
import json
import sys
from datetime import datetime, timezone
from typing import Dict, List, Optional

# Shared library imports
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file as _shared_sha256
from clawbio.common.report import write_result_json

# ---------------------------------------------------------------------------
# File-type routing
# ---------------------------------------------------------------------------

EXTENSION_MAP: dict[str, str] = {
    ".vcf": "equity-scorer",
    ".vcf.gz": "equity-scorer",
    ".fastq": "seq-wrangler",
    ".fastq.gz": "seq-wrangler",
    ".fq": "seq-wrangler",
    ".fq.gz": "seq-wrangler",
    ".bam": "seq-wrangler",
    ".cram": "seq-wrangler",
    ".pdb": "struct-predictor",
    ".cif": "struct-predictor",
    ".h5ad": "scrna-orchestrator",
    ".mtx": "scrna-orchestrator",
    ".mtx.gz": "scrna-orchestrator",
    ".pkl": "methylation-clock",
    ".pickle": "methylation-clock",
    ".csv": "equity-scorer",
    ".tsv": "equity-scorer",
    ".png": "data-extractor",
    ".jpg": "data-extractor",
    ".jpeg": "data-extractor",
    ".tiff": "data-extractor",
    ".tif": "data-extractor",
}

KEYWORD_MAP: dict[str, str] = {
    "illumina connected analytics": "illumina-bridge",
    "connected analytics": "illumina-bridge",
    "sample sheet": "illumina-bridge",
    "samplesheet": "illumina-bridge",
    "basespace": "illumina-bridge",
    "dragen": "illumina-bridge",
    "illumina": "illumina-bridge",
    "scvi": "scrna-embedding",
    "scanvi": "scrna-embedding",
    "batch correction": "scrna-embedding",
    "batch integration": "scrna-embedding",
    "integration": "scrna-embedding",
    "latent": "scrna-embedding",
    "embedding": "scrna-embedding",
    "x_scvi": "scrna-orchestrator",
    "integrated.h5ad": "scrna-orchestrator",
    "integrated h5ad": "scrna-orchestrator",
    "diversity": "equity-scorer",
    "equity": "equity-scorer",
    "heim": "equity-scorer",
    "heterozygosity": "equity-scorer",
    "fst": "equity-scorer",
    "variant annotation": "vcf-annotator",
    "annotate variant": "vcf-annotator",
    "variant": "vcf-annotator",
    "annotate": "vcf-annotator",
    "vep": "vcf-annotator",
    "structure": "struct-predictor",
    "alphafold": "struct-predictor",
    "fold": "struct-predictor",
    "single-cell": "scrna-orchestrator",
    "scrna": "scrna-orchestrator",
    "cluster": "scrna-orchestrator",
    "literature": "lit-synthesizer",
    "pubmed": "lit-synthesizer",
    "papers": "lit-synthesizer",
    "fastq": "seq-wrangler",
    "alignment": "seq-wrangler",
    "qc": "seq-wrangler",
    "reproducible": "repro-enforcer",
    "nextflow": "repro-enforcer",
    "singularity": "repro-enforcer",
    "conda": "repro-enforcer",
    "labstep": "labstep",
    "clinpgx": "clinpgx",
    "gene-drug pair": "clinpgx",
    "gene drug pair": "clinpgx",
    "cpic guideline": "clinpgx",
    "drug label": "clinpgx",
    "pharmgkb": "clinpgx",
    "clinical annotation": "clinpgx",
    "compare": "genome-compare",
    "corpasome": "genome-compare",
    "ibs": "genome-compare",
    "dna in common": "genome-compare",
    "george church": "genome-compare",
    "genome comparison": "genome-compare",
    "prs": "gwas-prs",
    "polygenic": "gwas-prs",
    "risk score": "gwas-prs",
    "polygenic risk": "gwas-prs",
    "gwas lookup": "gwas-lookup",
    "variant lookup": "gwas-lookup",
    "rs lookup": "gwas-lookup",
    "rsid": "gwas-lookup",
    "look up rs": "gwas-lookup",
    "lookup rs": "gwas-lookup",
    "phewas": "gwas-lookup",
    "gwas": "gwas-lookup",
    "profile report": "profile-report",
    "personal profile": "profile-report",
    "my profile": "profile-report",
    "genomic profile": "profile-report",
    "digitize": "data-extractor",
    "extract data": "data-extractor",
    "plot data": "data-extractor",
    "figure data": "data-extractor",
    "read chart": "data-extractor",
    "bar chart": "data-extractor",
    "scatter plot": "data-extractor",
    "meta-analysis": "data-extractor",
    "bioconductor": "bioconductor-bridge",
    "biocmanager": "bioconductor-bridge",
    "summarizedexperiment": "bioconductor-bridge",
    "singlecellexperiment": "bioconductor-bridge",
    "genomicranges": "bioconductor-bridge",
    "variantannotation": "bioconductor-bridge",
    "annotationhub": "bioconductor-bridge",
    "experimenthub": "bioconductor-bridge",
    "what package should i use": "bioconductor-bridge",
    "which bioconductor package": "bioconductor-bridge",
    "set up bioconductor": "bioconductor-bridge",
    "setup bioconductor": "bioconductor-bridge",
    "visualize de results": "diff-visualizer",
    "visualise de results": "diff-visualizer",
    "flow": "flow-bio",
    "flow.bio": "flow-bio",
    "flow bio": "flow-bio",
    "flow pipeline": "flow-bio",
    "flow sample": "flow-bio",
    "flow execution": "flow-bio",
    "run on flow": "flow-bio",
    "flow upload": "flow-bio",
    "de visualization": "diff-visualizer",
    "differential expression visualization": "diff-visualizer",
    "marker heatmap": "diff-visualizer",
    "marker dotplot": "diff-visualizer",
    "top genes heatmap": "diff-visualizer",
    "differential expression": "rnaseq-de",
    "deseq2": "rnaseq-de",
    "pydeseq2": "rnaseq-de",
    "bulk rna": "rnaseq-de",
    "rna-seq": "rnaseq-de",
    "volcano plot": "rnaseq-de",
    "ma plot": "rnaseq-de",
    "contrast": "rnaseq-de",
    "count matrix": "rnaseq-de",
    "epigenetic age": "methylation-clock",
    "methylation": "methylation-clock",
    "methylation clock": "methylation-clock",
    "dna methylation": "methylation-clock",
    "pyaging": "methylation-clock",
    "horvath": "methylation-clock",
    "altumage": "methylation-clock",
    "grimage": "methylation-clock",
    "dunedinpace": "methylation-clock",
    "geo accession": "methylation-clock",
    "gse": "methylation-clock",
}

SKILLS_DIR = Path(__file__).resolve().parent.parent
SCRNA_LATENT_ARTIFACT_TERMS = (
    "x_scvi",
    "integrated.h5ad",
    "integrated h5ad",
    "after scvi",
    "after scvi embedding",
)
SCRNA_DOWNSTREAM_TERMS = (
    "marker",
    "markers",
    "annotation",
    "annotate",
    "celltypist",
    "contrastive",
    "cluster",
    "clustering",
)
SCRNA_EMBEDDING_TERMS = (
    "scvi",
    "latent",
    "embedding",
    "integration",
    "batch correction",
    "batch integration",
)

ILLUMINA_SAMPLE_SHEET_NAMES = {"samplesheet.csv"}
ILLUMINA_VCF_SUFFIXES = {".vcf", ".vcf.gz"}


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def _looks_like_illumina_bundle(filepath: Path) -> bool:
    """Heuristic detection for DRAGEN-style export directories."""

    if not filepath.exists() or not filepath.is_dir():
        return False
    has_sample_sheet = any(
        candidate.is_file() and candidate.name.lower() in ILLUMINA_SAMPLE_SHEET_NAMES
        for candidate in filepath.rglob("*")
    )
    has_vcf = any(
        candidate.is_file() and "".join(candidate.suffixes).lower() in ILLUMINA_VCF_SUFFIXES
        for candidate in filepath.rglob("*")
    )
    return has_sample_sheet and has_vcf


def detect_skill_from_file(filepath: Path) -> str | None:
    """Determine which skill handles a given file based on extension."""
    if filepath.is_dir():
        if _looks_like_illumina_bundle(filepath):
            return "illumina-bridge"
        return None
    if filepath.name.lower() in ILLUMINA_SAMPLE_SHEET_NAMES:
        return "illumina-bridge"
    suffixes = "".join(filepath.suffixes)  # handles .vcf.gz
    if filepath.suffix.lower() in {".csv", ".tsv"}:
        inferred = detect_skill_from_tabular_header(filepath)
        if inferred:
            return inferred
    if suffixes in EXTENSION_MAP:
        return EXTENSION_MAP[suffixes]
    suffix = filepath.suffix.lower()
    return EXTENSION_MAP.get(suffix)


def detect_skill_from_tabular_header(filepath: Path) -> str | None:
    """Detect skill from tabular headers for CSV/TSV input files."""
    try:
        sep = "\t" if filepath.suffix.lower() == ".tsv" else ","
        with open(filepath, "r", encoding="utf-8") as f:
            first_line = f.readline().strip().lower()
            second_line = f.readline().strip().lower()
    except Exception:
        return None

    if not first_line:
        return None

    headers = [h.strip() for h in first_line.split(sep)]
    header_set = set(headers)

    if {"gene", "log2foldchange"} <= header_set and ({"padj", "pvalue"} & header_set):
        return "diff-visualizer"
    if {"cluster", "names", "scores"} <= header_set:
        return "diff-visualizer"
    if {"names", "scores"} <= header_set:
        return "diff-visualizer"

    equity_markers = {"population", "ancestry", "superpopulation", "ethnicity", "country"}
    if header_set & equity_markers:
        return "equity-scorer"

    rnaseq_metadata_markers = {"condition", "batch", "group", "treatment", "donor", "cell_type"}
    if "sample_id" in header_set and (header_set & rnaseq_metadata_markers):
        return "rnaseq-de"

    gene_like = {"gene", "gene_id", "ensembl_id", "symbol"}
    if headers and headers[0] in gene_like and len(headers) >= 4 and second_line:
        values = [value.strip() for value in second_line.split(sep)]
        numeric_count = 0
        for value in values[1:]:
            try:
                float(value)
                numeric_count += 1
            except ValueError:
                continue
        if numeric_count >= 3:
            return "rnaseq-de"

    methylation_markers = {"gender", "sex", "female", "tissue_type", "dataset"}
    if header_set & methylation_markers:
        cg_like = [h for h in headers if h.startswith("cg")]
        if len(cg_like) >= 10:
            return "methylation-clock"

    return None


def detect_skill_from_query(query: str) -> str | None:
    """Determine which skill matches a natural language query."""
    skill, _ = detect_skill_with_hint_from_query(query)
    return skill


def detect_skill_with_hint_from_query(query: str) -> tuple[str | None, str]:
    """Determine which skill matches a natural language query and explain chain-aware routing."""
    query_lower = query.lower()
    wants_embedding = any(term in query_lower for term in SCRNA_EMBEDDING_TERMS)
    wants_downstream = any(term in query_lower for term in SCRNA_DOWNSTREAM_TERMS)
    has_latent_artifact = any(term in query_lower for term in SCRNA_LATENT_ARTIFACT_TERMS)

    # Chain-aware scRNA routing favors explicit embedding requests unless the
    # user is clearly asking for downstream analysis on an existing latent artifact.
    if has_latent_artifact and wants_downstream:
        return (
            "scrna-orchestrator",
            "Detected a downstream latent-analysis workflow. Use `scrna-orchestrator` "
            "with `--use-rep X_scvi` on `integrated.h5ad` to run clustering, annotation, "
            "and contrastive markers after scVI.",
        )
    if wants_embedding and wants_downstream:
        return (
            "scrna-embedding",
            "Detected a two-step advanced scRNA workflow. First run `scrna-embedding` to "
            "produce `integrated.h5ad`, then run `scrna-orchestrator` with "
            "`--use-rep X_scvi` for downstream clustering, annotation, and contrastive markers.",
        )
    if wants_embedding and has_latent_artifact:
        return (
            "scrna-embedding",
            "Detected an embedding-focused scRNA workflow on an existing latent artifact. "
            "Use `scrna-embedding` to refresh or rebuild the scVI/scANVI latent space "
            "before downstream clustering or annotation.",
        )
    if wants_embedding:
        return (
            "scrna-embedding",
            "Detected an embedding-focused scRNA workflow. Use `scrna-embedding` to "
            "produce `integrated.h5ad` with a scVI/scANVI latent space before "
            "running downstream clustering or annotation.",
        )

    # Prefer longest keyword match to avoid ambiguity (e.g. "variant annotation"
    # should match vcf-annotator, not equity-scorer via "variant" substring)
    best_skill = None
    best_len = 0
    for keyword, skill in KEYWORD_MAP.items():
        if keyword in query_lower and len(keyword) > best_len:
            best_skill = skill
            best_len = len(keyword)
    if best_skill:
        return best_skill, ""
    return None, ""


def detect_routing_hint_for_file(filepath: Path) -> str:
    """Return a routing hint for special-case input files."""
    if filepath.is_dir() and _looks_like_illumina_bundle(filepath):
        return (
            "Detected an Illumina-style export bundle. Use `illumina-bridge` to "
            "normalize SampleSheet, VCF, and QC metrics before downstream analysis."
        )
    if filepath.name == "integrated.h5ad":
        return (
            "Detected `integrated.h5ad`. This is usually the downstream artifact from "
            "`scrna-embedding`; `scrna-orchestrator` can consume it with `--use-rep X_scvi`."
        )
    return ""


def detect_skill_with_flock(query: str) -> tuple[str | None, str]:
    """Use FLock API (open-source LLM) to route ambiguous queries.

    Returns (skill_name, reasoning) or (None, error_message).
    Falls back gracefully if FLock is not configured.
    """
    try:
        from clawbio.providers.flock import FlockRouter
        router = FlockRouter()
        result = router.route_query_safe(query)
        skill = result.get("skill")
        reasoning = result.get("reasoning", "")
        confidence = result.get("confidence", 0.0)
        if skill and confidence >= 0.5:
            return skill, f"FLock LLM routing (confidence={confidence:.1%}): {reasoning}"
        return None, f"FLock LLM low confidence ({confidence:.1%}): {reasoning}"
    except (ImportError, ValueError) as e:
        return None, f"FLock not available: {e}"


def sha256_file(filepath: Path) -> str:
    """Compute SHA-256 checksum of a file (delegates to shared library)."""
    return _shared_sha256(filepath)


def list_available_skills() -> list[str]:
    """List all skill directories that contain a SKILL.md."""
    skills = []
    for d in sorted(SKILLS_DIR.iterdir()):
        if d.is_dir() and (d / "SKILL.md").exists():
            skills.append(d.name)
    return skills


def skill_has_executable(skill_name: str) -> bool:
    """Check if a skill has a runnable Python executable.

    Returns False (stub) if:
      - No .py files exist (SKILL.md only), OR
      - SKILL.md YAML frontmatter declares required bins (``anyBins``) that
        are not found on PATH (e.g. ``boltz`` for struct-predictor).
    """
    import shutil
    skill_dir = SKILLS_DIR / skill_name
    if not skill_dir.is_dir():
        return False
    has_py = any(f.suffix == ".py" and f.name != "__init__.py"
                 for f in skill_dir.iterdir() if f.is_file())
    if not has_py:
        return False

    # Check SKILL.md YAML frontmatter for required external binaries
    skill_md = skill_dir / "SKILL.md"
    if skill_md.exists():
        try:
            text = skill_md.read_text(encoding="utf-8")
            # Quick YAML frontmatter parse (between --- fences)
            if text.startswith("---"):
                end = text.index("---", 3)
                front = text[3:end]
                # Look for anyBins list in openclaw.requires
                import re
                any_bins_match = re.search(r"anyBins:\s*\n((?:\s+-\s+\S+\n?)+)", front)
                if any_bins_match:
                    bins_block = any_bins_match.group(1)
                    bins = [line.strip().lstrip("- ").strip()
                            for line in bins_block.strip().splitlines()
                            if line.strip().startswith("-")]
                    # If any required bin is not on PATH, treat as stub
                    if bins and not any(shutil.which(b) for b in bins):
                        return False
        except Exception:
            pass  # If parsing fails, assume executable

    return True


def generate_report_header(
    title: str,
    skills_used: list[str],
    input_files: list[Path],
) -> str:
    """Generate the standard report header in markdown."""
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    checksums = []
    for f in input_files:
        if f.exists():
            checksums.append(f"- `{f.name}`: `{sha256_file(f)}`")
        else:
            checksums.append(f"- `{f.name}`: (file not found)")

    return f"""# Analysis Report: {title}

**Date**: {now}
**Skills used**: {', '.join(skills_used)}
**Input files**:
{chr(10).join(checksums)}

---
"""


def append_audit_log(output_dir: Path, action: str, details: str = "") -> None:
    """Append an entry to the audit log."""
    log_file = output_dir / "analysis_log.md"
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    entry = f"- **{now}**: {action}"
    if details:
        entry += f" -- {details}"
    entry += "\n"

    with open(log_file, "a") as f:
        if not log_file.exists() or log_file.stat().st_size == 0:
            f.write("# Analysis Audit Log\n\n")
        f.write(entry)


# ---------------------------------------------------------------------------
# Multi-skill routing
# ---------------------------------------------------------------------------

# Maps orchestrator skill names to clawbio.py skill registry names
SKILL_REGISTRY_MAP: dict[str, str] = {
    "pharmgx-reporter": "pharmgx",
    "equity-scorer": "equity",
    "nutrigx_advisor": "nutrigx",
    "scrna-orchestrator": "scrna",
    "scrna-embedding": "scrna-embedding",
    "genome-compare": "compare",
    "gwas-prs": "prs",
    "clinpgx": "clinpgx",
    "gwas-lookup": "gwas",
    "profile-report": "profile",
    "illumina-bridge": "illumina",
    "bioconductor-bridge": "bioc",
    "data-extractor": "data-extract",
    "rnaseq-de": "rnaseq",
    "diff-visualizer": "diffviz",
    "flow-bio": "flow",
}


def detect_multiple_skills(query: str) -> list[str]:
    """Detect all matching skills from a query (not just the first one).

    Returns a list of skill directory names.
    """
    skill, _ = detect_skill_with_hint_from_query(query)
    if skill in {"scrna-embedding", "scrna-orchestrator"}:
        return [skill]

    query_lower = query.lower()
    matched = []
    seen = set()
    for keyword, skill in KEYWORD_MAP.items():
        if keyword in query_lower and skill not in seen:
            matched.append(skill)
            seen.add(skill)
    return matched


def route_to_clawbio(
    skills: list[str],
    input_path: str | None = None,
    profile_path: str | None = None,
    output_dir: str | None = None,
) -> dict:
    """Route to clawbio.py's run_skill for each detected skill.

    Returns a summary dict with per-skill results.
    """
    # Import clawbio.py runner (not the clawbio/ package)
    import importlib.util
    spec = importlib.util.spec_from_file_location("clawbio_runner", _PROJECT_ROOT / "clawbio.py")
    _runner = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(_runner)
    run_skill = _runner.run_skill

    results = {}
    for skill_dir_name in skills:
        # Map orchestrator name to clawbio.py registry name
        registry_name = SKILL_REGISTRY_MAP.get(skill_dir_name, skill_dir_name)

        skill_output = None
        if output_dir:
            skill_output = str(Path(output_dir) / registry_name)

        result = run_skill(
            skill_name=registry_name,
            input_path=input_path,
            output_dir=skill_output,
            profile_path=profile_path,
        )
        results[registry_name] = {
            "success": result["success"],
            "exit_code": result["exit_code"],
            "output_dir": result["output_dir"],
            "files": result["files"],
        }

    return results


# ---------------------------------------------------------------------------
# CLI entry point
# ---------------------------------------------------------------------------

def main() -> None:
    parser = argparse.ArgumentParser(description="Bio Orchestrator: route bioinformatics requests")
    parser.add_argument("--input", "-i", help="Input file path or natural language query")
    parser.add_argument("--skill", "-s", help="Force a specific skill (bypasses auto-detection)")
    parser.add_argument("--skills", help="Comma-separated list of skills to run (multi-skill mode)")
    parser.add_argument("--profile", "-p", help="Path to patient profile JSON (enables profile-aware dispatch)")
    parser.add_argument("--output", "-o", default=".", help="Output directory for reports")
    parser.add_argument("--list-skills", action="store_true", help="List available skills")
    parser.add_argument("--multi", action="store_true", help="Detect and run all matching skills (not just first)")
    parser.add_argument("--provider", choices=["keyword", "flock"], default="keyword",
                        help="Routing strategy: 'keyword' (default, rule-based) or 'flock' (open-source LLM via FLock API)")
    args = parser.parse_args()

    if args.list_skills:
        skills = list_available_skills()
        print("Available skills:")
        for s in skills:
            print(f"  - {s}")
        return

    # Multi-skill mode: explicit skill list
    if args.skills:
        skill_list = [s.strip() for s in args.skills.split(",") if s.strip()]
        print(f"Multi-skill mode: running {skill_list}")
        results = route_to_clawbio(
            skills=skill_list,
            input_path=args.input,
            profile_path=args.profile,
            output_dir=args.output,
        )
        print(json.dumps(results, indent=2))

        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        append_audit_log(output_dir, f"Multi-skill: {skill_list}", f"input={args.input}")

        # Write result.json for orchestration
        write_result_json(
            output_dir=output_dir,
            skill="bio-orchestrator",
            version="0.2.0",
            summary={"skills_run": skill_list, "all_success": all(r["success"] for r in results.values())},
            data=results,
        )
        return

    if not args.input and not args.profile:
        parser.print_help()
        sys.exit(1)

    # Single-skill detection
    if args.input:
        input_path = Path(args.input)
    else:
        input_path = None
    routing_hint = ""

    if args.skill:
        # SEC INT-002: reject path traversal in skill name
        if "/" in args.skill or "\\" in args.skill or ".." in args.skill:
            print(f"Invalid skill name: {args.skill}")
            sys.exit(1)
        skill = args.skill
        method = "user-specified"
    elif input_path and input_path.exists():
        skill = detect_skill_from_file(input_path)
        method = "file-extension"
        routing_hint = detect_routing_hint_for_file(input_path)
    elif args.input:
        # Multi-detect mode: find all matching skills
        if args.multi:
            skills = detect_multiple_skills(args.input)
            if skills:
                print(f"Detected {len(skills)} skills: {skills}")
                results = route_to_clawbio(
                    skills=skills,
                    input_path=args.input if input_path and input_path.exists() else None,
                    profile_path=args.profile,
                    output_dir=args.output,
                )
                print(json.dumps(results, indent=2))
                return
        skill, routing_hint = detect_skill_with_hint_from_query(args.input)
        method = "keyword"
    else:
        skill = None
        method = "none"
        routing_hint = ""

    # Fallback: if keyword matching failed, try FLock LLM routing
    if not skill and args.provider == "flock" and args.input:
        print("Keyword matching failed. Trying FLock LLM routing (open-source model)...")
        try:
            skill, reasoning = detect_skill_with_flock(args.input)
            method = "flock-llm"
            if skill:
                print(f"FLock routed to: {skill} — {reasoning}")
            else:
                print(f"FLock routing: {reasoning}")
        except Exception as exc:
            method = "flock-llm"
            reasoning = f"FLock routing failed: {exc}"
            print(reasoning, file=sys.stderr)
            # Return structured JSON error instead of crashing
            error_result = {
                "input": args.input,
                "detected_skill": None,
                "detection_method": method,
                "error": str(exc),
                "available_skills": list_available_skills(),
            }
            output_dir = Path(args.output)
            output_dir.mkdir(parents=True, exist_ok=True)
            write_result_json(
                output_dir=output_dir,
                skill="bio-orchestrator",
                version="0.2.0",
                summary={"detected_skill": None, "method": method, "error": str(exc)},
                data=error_result,
            )
            print(json.dumps(error_result, indent=2))
            sys.exit(1)

    # FLock fallback removed: keyword mode must not silently send queries
    # to an external API. Use --provider flock explicitly to opt in.

    if not skill:
        print(f"Could not determine skill for input: {args.input}")
        print("Available skills:", ", ".join(list_available_skills()))
        # Emit structured JSON result for harness consumption and exit 0.
        # A no-match is a valid outcome (not a crash), especially when the
        # requested provider (e.g. flock) is unavailable.
        no_match_result = {
            "input": args.input,
            "detected_skill": None,
            "detection_method": method,
            "confidence": 0.0,
            "available_skills": list_available_skills(),
        }
        output_dir = Path(args.output)
        output_dir.mkdir(parents=True, exist_ok=True)
        write_result_json(
            output_dir=output_dir,
            skill="bio-orchestrator",
            version="0.2.0",
            summary={"detected_skill": None, "method": method},
            data=no_match_result,
        )
        print(json.dumps(no_match_result, indent=2))
        sys.exit(1)

    # Check skill exists
    skill_dir = (SKILLS_DIR / skill).resolve()
    # SEC INT-002: ensure resolved path stays within SKILLS_DIR
    if not str(skill_dir).startswith(str(SKILLS_DIR.resolve())):
        print(f"Invalid skill name: {skill}")
        sys.exit(1)
    if not (skill_dir / "SKILL.md").exists():
        print(f"Skill '{skill}' not found")
        sys.exit(1)

    # Warn if skill is a stub (SKILL.md only, no Python executable)
    is_stub = not skill_has_executable(skill)
    if is_stub:
        print(
            f"WARNING: '{skill}' is a SKILL.md-only stub with no Python executable. "
            f"The agent can apply the methodology from SKILL.md but cannot run automated analysis.",
            file=sys.stderr,
        )

    # Output routing decision
    result = {
        "input": args.input,
        "detected_skill": skill,
        "detection_method": method,
        "skill_dir": str(skill_dir),
        "is_stub": is_stub,
        "available_skills": list_available_skills(),
    }
    if routing_hint:
        result["routing_hint"] = routing_hint
    if args.profile:
        result["profile"] = args.profile
    print(json.dumps(result, indent=2))

    # Log the routing
    output_dir = Path(args.output)
    output_dir.mkdir(parents=True, exist_ok=True)
    append_audit_log(output_dir, f"Routed to {skill}", f"input={args.input}, method={method}")

    # Write result.json
    write_result_json(
        output_dir=output_dir,
        skill="bio-orchestrator",
        version="0.2.0",
        summary={"detected_skill": skill, "method": method},
        data=result,
    )


if __name__ == "__main__":
    main()
