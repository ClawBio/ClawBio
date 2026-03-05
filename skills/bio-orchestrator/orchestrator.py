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
    ".csv": "equity-scorer",
    ".tsv": "equity-scorer",
}

KEYWORD_MAP: dict[str, str] = {
    "diversity": "equity-scorer",
    "equity": "equity-scorer",
    "heim": "equity-scorer",
    "heterozygosity": "equity-scorer",
    "fst": "equity-scorer",
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
}

SKILLS_DIR = Path(__file__).resolve().parent.parent


# ---------------------------------------------------------------------------
# Utility functions
# ---------------------------------------------------------------------------

def detect_skill_from_file(filepath: Path) -> str | None:
    """Determine which skill handles a given file based on extension."""
    suffixes = "".join(filepath.suffixes)  # handles .vcf.gz
    if suffixes in EXTENSION_MAP:
        return EXTENSION_MAP[suffixes]
    suffix = filepath.suffix.lower()
    return EXTENSION_MAP.get(suffix)


def detect_skill_from_query(query: str) -> str | None:
    """Determine which skill matches a natural language query."""
    query_lower = query.lower()
    for keyword, skill in KEYWORD_MAP.items():
        if keyword in query_lower:
            return skill
    return None


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
    "genome-compare": "compare",
    "gwas-prs": "prs",
    "clinpgx": "clinpgx",
    "gwas-lookup": "gwas",
    "profile-report": "profile",
}


def detect_multiple_skills(query: str) -> list[str]:
    """Detect all matching skills from a query (not just the first one).

    Returns a list of skill directory names.
    """
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
        skill = detect_skill_from_query(args.input)
        method = "keyword"
    else:
        skill = None
        method = "none"

    if not skill:
        print(f"Could not determine skill for input: {args.input}")
        print("Available skills:", ", ".join(list_available_skills()))
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

    # Output routing decision
    result = {
        "input": args.input,
        "detected_skill": skill,
        "detection_method": method,
        "skill_dir": str(skill_dir),
        "available_skills": list_available_skills(),
    }
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
