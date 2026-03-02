#!/usr/bin/env python3
"""
clawbio.py — ClawBio Bioinformatics Skills Runner
==================================================
Standalone CLI and importable module for running ClawBio skills.

Usage:
    python clawbio.py list
    python clawbio.py run pharmgx --demo
    python clawbio.py run equity --input data.vcf
    python clawbio.py run pharmgx --input patient.txt --output ./results

Importable:
    from clawbio import run_skill, list_skills
    result = run_skill("pharmgx", demo=True)
"""

import argparse
import subprocess
import sys
import time
from datetime import datetime
from pathlib import Path

# --------------------------------------------------------------------------- #
# Paths
# --------------------------------------------------------------------------- #

CLAWBIO_DIR = Path(__file__).resolve().parent
SKILLS_DIR = CLAWBIO_DIR / "skills"
EXAMPLES_DIR = CLAWBIO_DIR / "examples"
DEFAULT_OUTPUT_ROOT = CLAWBIO_DIR / "output"

# Python binary — project standard
PYTHON = "python3.11"

# --------------------------------------------------------------------------- #
# Skills registry
# --------------------------------------------------------------------------- #

SKILLS = {
    "pharmgx": {
        "script": SKILLS_DIR / "pharmgx-reporter" / "pharmgx_reporter.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "pharmgx-reporter" / "demo_patient.txt"),
        ],
        "description": "Pharmacogenomics reporter (12 genes, 31 SNPs, 51 drugs)",
    },
    "equity": {
        "script": SKILLS_DIR / "equity-scorer" / "equity_scorer.py",
        "demo_args": [
            "--input",
            str(EXAMPLES_DIR / "demo_populations.vcf"),
            "--pop-map",
            str(EXAMPLES_DIR / "demo_population_map.csv"),
        ],
        "description": "HEIM equity scorer (FST, heterozygosity, population representation)",
    },
    "nutrigx": {
        "script": SKILLS_DIR / "nutrigx_advisor" / "nutrigx_advisor.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "nutrigx_advisor" / "tests" / "synthetic_patient.csv"),
        ],
        "description": "Nutrigenomics advisor (diet, vitamins, caffeine, lactose)",
    },
    "metagenomics": {
        "script": SKILLS_DIR / "claw-metagenomics" / "metagenomics_profiler.py",
        "demo_args": ["--demo"],
        "description": "Metagenomics profiler (Kraken2, RGI/CARD, HUMAnN3)",
    },
    "compare": {
        "script": SKILLS_DIR / "genome-compare" / "genome_compare.py",
        "demo_args": ["--demo"],
        "description": "Genome comparator (IBS vs George Church + ancestry estimation)",
    },
}

# --------------------------------------------------------------------------- #
# list_skills
# --------------------------------------------------------------------------- #


def list_skills() -> dict:
    """Print available skills and return the registry dict."""
    print("ClawBio Skills")
    print("=" * 55)
    for name, info in SKILLS.items():
        script_exists = info["script"].exists()
        status = "OK" if script_exists else "MISSING"
        print(f"  {name:<15} {info['description']}")
        print(f"  {'':15} script: {info['script'].name} [{status}]")
        print()
    print(f"Run a skill:  python clawbio.py run <skill> --demo")
    print(f"With input:   python clawbio.py run <skill> --input <file>")
    return SKILLS


# --------------------------------------------------------------------------- #
# run_skill
# --------------------------------------------------------------------------- #


def run_skill(
    skill_name: str,
    input_path: str | None = None,
    output_dir: str | None = None,
    demo: bool = False,
    extra_args: list[str] | None = None,
    timeout: int = 300,
) -> dict:
    """
    Run a ClawBio skill as a subprocess.

    Returns a structured dict with success status, output paths, and logs.
    Importable by any agent (RoboTerri, RoboIsaac, Claude Code).
    """
    # Validate skill
    skill_info = SKILLS.get(skill_name)
    if not skill_info:
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": None,
            "files": [],
            "stdout": "",
            "stderr": f"Unknown skill '{skill_name}'. Available: {list(SKILLS.keys())}",
            "duration_seconds": 0,
        }

    script_path = skill_info["script"]
    if not script_path.exists():
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": None,
            "files": [],
            "stdout": "",
            "stderr": f"Script not found: {script_path}",
            "duration_seconds": 0,
        }

    # Build output directory
    if output_dir:
        out_dir = Path(output_dir)
    else:
        ts = datetime.now().strftime("%Y%m%d_%H%M%S")
        out_dir = DEFAULT_OUTPUT_ROOT / f"{skill_name}_{ts}"
    out_dir.mkdir(parents=True, exist_ok=True)

    # Build command
    cmd = [PYTHON, str(script_path)]

    if demo:
        cmd.extend(skill_info["demo_args"])
    elif input_path:
        cmd.extend(["--input", str(input_path)])
    else:
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": str(out_dir),
            "files": [],
            "stdout": "",
            "stderr": "No input provided. Use --demo or --input <file>.",
            "duration_seconds": 0,
        }

    cmd.extend(["--output", str(out_dir)])
    if extra_args:
        cmd.extend(extra_args)

    # Run subprocess
    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=timeout,
            cwd=str(script_path.parent),
        )
        duration = round(time.time() - t0, 2)
    except subprocess.TimeoutExpired:
        duration = round(time.time() - t0, 2)
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": str(out_dir),
            "files": [],
            "stdout": "",
            "stderr": f"Timed out after {timeout} seconds.",
            "duration_seconds": duration,
        }
    except Exception as e:
        duration = round(time.time() - t0, 2)
        return {
            "skill": skill_name,
            "success": False,
            "exit_code": -1,
            "output_dir": str(out_dir),
            "files": [],
            "stdout": "",
            "stderr": str(e),
            "duration_seconds": duration,
        }

    # Collect output files
    output_files = sorted(
        [f.name for f in out_dir.rglob("*") if f.is_file()],
    )

    return {
        "skill": skill_name,
        "success": proc.returncode == 0,
        "exit_code": proc.returncode,
        "output_dir": str(out_dir),
        "files": output_files,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "duration_seconds": duration,
    }


# --------------------------------------------------------------------------- #
# CLI entry point
# --------------------------------------------------------------------------- #


def main():
    parser = argparse.ArgumentParser(
        description="ClawBio — Bioinformatics Skills Runner",
    )
    sub = parser.add_subparsers(dest="command")

    # list
    sub.add_parser("list", help="List available skills")

    # run
    run_parser = sub.add_parser("run", help="Run a skill")
    run_parser.add_argument("skill", help="Skill name (e.g. pharmgx, equity)")
    run_parser.add_argument("--demo", action="store_true", help="Run with demo data")
    run_parser.add_argument("--input", dest="input_path", help="Path to input file")
    run_parser.add_argument("--output", dest="output_dir", help="Output directory")
    run_parser.add_argument(
        "--timeout", type=int, default=300, help="Timeout in seconds (default: 300)"
    )

    args = parser.parse_args()

    if args.command == "list":
        list_skills()
    elif args.command == "run":
        print(f"Running {args.skill}...")
        result = run_skill(
            skill_name=args.skill,
            input_path=args.input_path,
            output_dir=args.output_dir,
            demo=args.demo,
            timeout=args.timeout,
        )
        print()
        if result["success"]:
            print(f"  Status:   OK (exit {result['exit_code']})")
        else:
            print(f"  Status:   FAILED (exit {result['exit_code']})")
        print(f"  Duration: {result['duration_seconds']}s")
        print(f"  Output:   {result['output_dir']}")
        if result["files"]:
            print(f"  Files:    {', '.join(result['files'])}")
        if not result["success"] and result["stderr"]:
            print(f"\n  Error:\n{result['stderr'][-800:]}")
        sys.exit(0 if result["success"] else 1)
    else:
        parser.print_help()


if __name__ == "__main__":
    main()
