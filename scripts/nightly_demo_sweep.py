#!/usr/bin/env python3
"""
nightly_demo_sweep.py - Run every ClawBio skill in demo mode.

Reads skills/catalog.json, runs each skill's demo_command, captures
exit code and stderr. Produces a markdown summary suitable for
GitHub Actions step summary.

Exit code 0 = all demos passed. Non-zero = at least one failed.

Usage:
    python scripts/nightly_demo_sweep.py
    python scripts/nightly_demo_sweep.py --timeout 120
    python scripts/nightly_demo_sweep.py --output /tmp/sweep_results
"""
from __future__ import annotations

import json
import subprocess
import sys
import time
from pathlib import Path

PROJECT_ROOT = Path(__file__).resolve().parents[1]
CATALOG = PROJECT_ROOT / "skills" / "catalog.json"

# Skills that need network access or heavy deps not available in CI
SKIP_IN_CI = {
    "bio-orchestrator",   # meta-router, needs sub-skills configured
    "bigquery-public",    # needs gcloud auth or bq CLI
    "flow-bio",           # needs Flow.bio credentials
    "labstep",            # needs Labstep API token
    "protocols-io",       # needs protocols.io auth for some endpoints
    "cell-detection",     # needs cellpose + torch (heavy)
    "struct-predictor",   # needs boltz (heavy)
    "methylation-clock",  # needs pyaging (heavy)
    "scrna-orchestrator", # needs scanpy + anndata (heavy)
    "galaxy-bridge",      # queries live Galaxy API
    "nfcore-sarek-wrapper",  # needs Nextflow + Java + Docker/Singularity (heavy)
}


def load_demo_skills() -> list[dict]:
    """Return skills from catalog that have demo mode."""
    with open(CATALOG) as f:
        catalog = json.load(f)

    demos = []
    for skill in catalog["skills"]:
        if skill.get("has_demo") and skill.get("demo_command"):
            demos.append(skill)
    return demos


def run_demo(skill: dict, timeout: int, output_dir: Path | None) -> dict:
    """Run a single skill's demo command. Return result dict."""
    name = skill["name"]
    cmd = skill["demo_command"]

    # Use the same Python that's running this script
    python = sys.executable
    cmd = cmd.replace("python clawbio.py", f"{python} clawbio.py")
    cmd = cmd.replace("python skills/", f"{python} skills/")

    # Only append --output if the demo_command already uses --output
    # (some skills like recombinator don't accept it)
    if output_dir and "--output" in skill.get("demo_command", ""):
        skill_out = output_dir / name
        skill_out.mkdir(parents=True, exist_ok=True)
        # Replace existing output path in demo command
        cmd = cmd.replace("--output /tmp/", f"--output {skill_out}/")
    elif output_dir and "--output" not in cmd:
        # Don't force --output on skills that don't support it
        pass

    t0 = time.time()
    try:
        proc = subprocess.run(
            cmd,
            shell=True,
            cwd=str(PROJECT_ROOT),
            capture_output=True,
            text=True,
            timeout=timeout,
        )
        elapsed = time.time() - t0
        return {
            "name": name,
            "command": cmd,
            "exit_code": proc.returncode,
            "elapsed": round(elapsed, 1),
            "stderr_tail": proc.stderr[-500:] if proc.stderr else "",
            "passed": proc.returncode == 0,
        }
    except subprocess.TimeoutExpired:
        elapsed = time.time() - t0
        return {
            "name": name,
            "command": cmd,
            "exit_code": -1,
            "elapsed": round(elapsed, 1),
            "stderr_tail": f"TIMEOUT after {timeout}s",
            "passed": False,
        }
    except Exception as e:
        elapsed = time.time() - t0
        return {
            "name": name,
            "command": cmd,
            "exit_code": -2,
            "elapsed": round(elapsed, 1),
            "stderr_tail": str(e),
            "passed": False,
        }


def collect_output_genes(output_dir: Path) -> dict[str, list[str]]:
    """Scan skill output directories for gene lists in JSON reports.

    Looks for 'gene', 'gene_symbol', or 'genes' fields in any .json file
    in the output directory tree. Returns {skill_name: [gene_list]}.
    """
    skill_genes: dict[str, list[str]] = {}
    if not output_dir or not output_dir.exists():
        return skill_genes

    for skill_dir in sorted(output_dir.iterdir()):
        if not skill_dir.is_dir():
            continue
        genes: set[str] = set()
        for jf in skill_dir.rglob("*.json"):
            try:
                data = json.loads(jf.read_text())
            except (json.JSONDecodeError, OSError):
                continue

            # Extract genes from various output formats
            if isinstance(data, dict):
                # Direct gene list
                if "genes" in data and isinstance(data["genes"], list):
                    genes.update(g for g in data["genes"] if isinstance(g, str))
                # Array of records with gene field
                for key in ("results", "variants", "associations", "credible_sets"):
                    items = data.get(key, [])
                    if isinstance(items, list):
                        for item in items:
                            if isinstance(item, dict):
                                for gf in ("gene", "gene_symbol", "gene_name"):
                                    if gf in item and isinstance(item[gf], str):
                                        genes.add(item[gf])
            elif isinstance(data, list):
                for item in data:
                    if isinstance(item, dict):
                        for gf in ("gene", "gene_symbol", "gene_name"):
                            if gf in item and isinstance(item[gf], str):
                                genes.add(item[gf])

        if genes:
            skill_genes[skill_dir.name] = sorted(genes)
    return skill_genes


def run_benchmark_scoring(output_dir: Path | None) -> dict | None:
    """Score collected genes against AD ground truth if available."""
    if not output_dir:
        return None

    benchmark_dir = PROJECT_ROOT / "tests" / "benchmark"
    gt_path = benchmark_dir / "ad_ground_truth.json"
    if not gt_path.exists():
        return None

    try:
        sys.path.insert(0, str(benchmark_dir))
        from benchmark_scorer import BenchmarkScorer
        scorer = BenchmarkScorer(gt_path)
    except Exception:
        return None

    skill_genes = collect_output_genes(output_dir)
    if not skill_genes:
        return None

    # Merge all discovered genes across skills
    all_genes = set()
    for genes in skill_genes.values():
        all_genes.update(genes)

    result = scorer.score(list(all_genes))
    result["skills_contributing"] = {k: len(v) for k, v in skill_genes.items()}
    return result


def generate_summary(results: list[dict], benchmark: dict | None = None) -> str:
    """Generate markdown summary of sweep results."""
    passed = [r for r in results if r["passed"]]
    failed = [r for r in results if not r["passed"]]
    skipped_names = sorted(SKIP_IN_CI)

    lines = [
        "# ClawBio Nightly Demo Sweep",
        "",
        f"**{len(passed)}/{len(results)} passed** | "
        f"**{len(failed)} failed** | "
        f"**{len(skipped_names)} skipped** (heavy deps / network)",
        "",
    ]

    if failed:
        lines.append("## Failures")
        lines.append("")
        lines.append("| Skill | Exit Code | Time | Error |")
        lines.append("|-------|-----------|------|-------|")
        for r in failed:
            err = r["stderr_tail"].replace("\n", " ")[:100]
            lines.append(
                f"| {r['name']} | {r['exit_code']} | {r['elapsed']}s | {err} |"
            )
        lines.append("")

    lines.append("## Passed")
    lines.append("")
    lines.append("| Skill | Time |")
    lines.append("|-------|------|")
    for r in sorted(passed, key=lambda x: x["name"]):
        lines.append(f"| {r['name']} | {r['elapsed']}s |")

    if skipped_names:
        lines.append("")
        lines.append("## Skipped (CI)")
        lines.append("")
        lines.append(", ".join(f"`{s}`" for s in skipped_names))

    # Benchmark scoring section
    if benchmark:
        status = "PASS" if benchmark.get("passes_minimum") else "FAIL"
        lines.append("")
        lines.append(f"## AD Benchmark [{status}]")
        lines.append("")
        lines.append(f"**Genes evaluated**: {benchmark['pipeline_genes_count']} "
                      f"(from {len(benchmark.get('skills_contributing', {}))} skills)")
        lines.append("")
        lines.append("| Metric | Value | Minimum |")
        lines.append("|--------|-------|---------|")
        mins = benchmark.get("minimums", {})
        lines.append(f"| Gene recovery | {benchmark['gene_recovery_rate']:.4f} | {mins.get('gene_recovery_rate', '-')} |")
        lines.append(f"| Precision | {benchmark['precision']:.4f} | {mins.get('precision', '-')} |")
        lines.append(f"| F1 | {benchmark['f1']:.4f} | {mins.get('f1', '-')} |")
        lines.append(f"| FDR | {benchmark['false_discovery_rate']:.4f} | - |")
        lines.append(f"| Weighted score | {benchmark['weighted_score']:.4f} | - |")

        tb = benchmark.get("tier_breakdown", {})
        lines.append("")
        t1f = len(tb.get("tier1_found", []))
        t1m = len(tb.get("tier1_missed", []))
        t2f = len(tb.get("tier2_found", []))
        t2m = len(tb.get("tier2_missed", []))
        t3f = len(tb.get("tier3_found", []))
        t3m = len(tb.get("tier3_missed", []))
        lines.append(f"Tier 1 (causal): {t1f}/{t1f+t1m} | "
                      f"Tier 2 (GWAS): {t2f}/{t2f+t2m} | "
                      f"Tier 3 (novel): {t3f}/{t3f+t3m}")

        if benchmark.get("false_positive_genes"):
            lines.append(f"\nFalse positives: {', '.join(benchmark['false_positive_genes'])}")

    lines.append("")
    return "\n".join(lines)


def main():
    import argparse

    parser = argparse.ArgumentParser(description="ClawBio nightly demo sweep")
    parser.add_argument(
        "--timeout", type=int, default=180, help="Per-skill timeout in seconds"
    )
    parser.add_argument(
        "--output", type=str, default=None, help="Directory for demo outputs"
    )
    parser.add_argument(
        "--include-heavy",
        action="store_true",
        help="Include skills that need heavy deps (cellpose, torch, etc.)",
    )
    args = parser.parse_args()

    output_dir = Path(args.output) if args.output else None
    if output_dir:
        output_dir.mkdir(parents=True, exist_ok=True)

    skills = load_demo_skills()
    print(f"Found {len(skills)} skills with demo mode")

    skip = set() if args.include_heavy else SKIP_IN_CI
    runnable = [s for s in skills if s["name"] not in skip]
    print(f"Running {len(runnable)} demos (skipping {len(skip)} heavy/network skills)")
    print()

    results = []
    for i, skill in enumerate(runnable, 1):
        name = skill["name"]
        print(f"[{i}/{len(runnable)}] {name} ... ", end="", flush=True)
        result = run_demo(skill, args.timeout, output_dir)
        status = "PASS" if result["passed"] else "FAIL"
        print(f"{status} ({result['elapsed']}s)")
        if not result["passed"] and result["stderr_tail"]:
            for line in result["stderr_tail"].strip().split("\n")[-3:]:
                print(f"         {line}")
        results.append(result)

    # Run benchmark scoring if output dir has gene data
    benchmark = run_benchmark_scoring(output_dir)
    if benchmark:
        print("Benchmark scoring complete.")
    else:
        print("No benchmark data collected (use --output to enable scoring).")

    print()
    summary = generate_summary(results, benchmark=benchmark)
    print(summary)

    # Write summary to file if output dir specified
    if output_dir:
        (output_dir / "sweep_summary.md").write_text(summary)
        (output_dir / "sweep_results.json").write_text(
            json.dumps(results, indent=2)
        )
        if benchmark:
            (output_dir / "benchmark_results.json").write_text(
                json.dumps(benchmark, indent=2)
            )

    # Write to GitHub step summary if available
    summary_file = Path.home() / "GITHUB_STEP_SUMMARY"
    env_summary = __import__("os").environ.get("GITHUB_STEP_SUMMARY")
    if env_summary:
        with open(env_summary, "a") as f:
            f.write(summary)

    failed = [r for r in results if not r["passed"]]
    if failed:
        print(f"\n{len(failed)} skill(s) failed. Exiting with code 1.")
        sys.exit(1)
    else:
        print(f"\nAll {len(results)} demos passed.")
        sys.exit(0)


if __name__ == "__main__":
    main()
