#!/usr/bin/env python3
"""
bioconductor_bridge.py — Bioconductor discovery and setup bridge for ClawBio
=============================================================================
Live recommendation, workflow suggestion, setup inspection, and opt-in
package installation for official Bioconductor workflows.
"""

from __future__ import annotations

import argparse
import csv
import importlib.util
import json
import os
import shlex
import shutil
import subprocess
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

SKILL_DIR = Path(__file__).resolve().parent
PROJECT_ROOT = SKILL_DIR.parent.parent
if str(PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(PROJECT_ROOT))

from bioc_recommender import (
    DOMAIN_DEFINITIONS,
    detect_domain,
    docs_search_catalog,
    fetch_package_documentation,
    get_package_details,
    infer_input_context,
    list_domains,
    load_catalog,
    render_starter_script,
    recommend_packages_with_docs,
    search_catalog_with_docs,
    suggest_workflow,
)


def _load_module(module_name: str, path: Path):
    spec = importlib.util.spec_from_file_location(module_name, path)
    module = importlib.util.module_from_spec(spec)
    assert spec is not None and spec.loader is not None
    spec.loader.exec_module(module)
    return module


_checksums = _load_module("clawbio_common_checksums", PROJECT_ROOT / "clawbio" / "common" / "checksums.py")

sha256_file = _checksums.sha256_file
DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare "
    "professional before making any medical decisions."
)


def generate_report_header(
    title: str,
    skill_name: str,
    input_files: list[Path] | None = None,
    extra_metadata: dict[str, str] | None = None,
) -> str:
    now = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC")
    checksums = []
    if input_files:
        for input_file in input_files:
            path = Path(input_file)
            if path.exists():
                checksums.append(f"- `{path.name}`: `{sha256_file(path)}`")
            else:
                checksums.append(f"- `{path.name}`: (file not found)")

    lines = [
        f"# {title}",
        "",
        f"**Date**: {now}",
        f"**Skill**: {skill_name}",
    ]
    if extra_metadata:
        for key, value in extra_metadata.items():
            lines.append(f"**{key}**: {value}")
    if checksums:
        lines.append("**Input files**:")
        lines.extend(checksums)
    lines.extend(["", "---", ""])
    return "\n".join(lines)


def generate_report_footer() -> str:
    return f"""
---

## Disclaimer

*{DISCLAIMER}*
"""


def write_result_json(
    output_dir: str | Path,
    skill: str,
    version: str,
    summary: dict[str, Any],
    data: dict[str, Any],
    input_checksum: str = "",
) -> Path:
    output_dir = Path(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    envelope = {
        "skill": skill,
        "version": version,
        "completed_at": datetime.now(timezone.utc).isoformat(),
        "input_checksum": f"sha256:{input_checksum}" if input_checksum else "",
        "summary": summary,
        "data": data,
    }
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(envelope, indent=2, default=str), encoding="utf-8")
    return result_path

VERSION = "0.1.0"
DEFAULT_OUTPUT = Path.cwd() / "bioconductor_bridge_output"


def run_process(cmd: list[str], **kwargs) -> subprocess.CompletedProcess[str]:
    """Wrapper for subprocess.run so tests can patch it cleanly."""
    return subprocess.run(cmd, **kwargs)


def warn_if_output_exists(output_dir: Path) -> None:
    """Warn before overwriting an existing non-empty output directory."""
    if output_dir.exists() and any(output_dir.iterdir()):
        print(
            f"WARNING: Output directory already contains files and may be overwritten: {output_dir}",
            file=sys.stderr,
        )


def ensure_output_dir(output_dir: Path) -> None:
    warn_if_output_exists(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)


def _rscript_path() -> str | None:
    return shutil.which("Rscript")


def _run_rscript(expr: str, timeout: int = 120) -> subprocess.CompletedProcess[str] | None:
    rscript = _rscript_path()
    if not rscript:
        return None
    return run_process(
        [rscript, "--vanilla", "-e", expr],
        capture_output=True,
        text=True,
        timeout=timeout,
        cwd=str(SKILL_DIR),
    )


def _r_quote_vector(packages: list[str]) -> str:
    return ", ".join(json.dumps(pkg) for pkg in packages)


def _shell_command_for_expr(expr: str) -> str:
    rscript = _rscript_path() or "Rscript"
    return f"{shlex.quote(rscript)} --vanilla -e {shlex.quote(expr)}"


def _parse_bool_line(text: str) -> bool | None:
    for line in reversed(text.splitlines()):
        stripped = line.strip()
        if stripped in {"TRUE", "FALSE"}:
            return stripped == "TRUE"
    return None


def _parse_key_value_lines(text: str) -> dict[str, str]:
    result: dict[str, str] = {}
    for line in text.splitlines():
        if "=" not in line:
            continue
        key, value = line.split("=", 1)
        result[key.strip()] = value.strip()
    return result


def inspect_setup(target_packages: list[str]) -> dict[str, Any]:
    """Inspect local R/BiocManager/package readiness with live validation."""
    status: dict[str, Any] = {
        "rscript_path": _rscript_path(),
        "r_available": False,
        "r_version": "",
        "r_status": "missing",
        "warnings": [],
        "biocmanager_installed": False,
        "biocmanager_version": "",
        "bioconductor_version": "",
        "installed_packages": {},
        "remote_validation": "not attempted",
    }

    if not status["rscript_path"]:
        status["warnings"].append("Rscript was not found on PATH.")
        return status

    version_proc = _run_rscript('cat(R.version.string)')
    if version_proc and version_proc.returncode == 0:
        version = version_proc.stdout.strip()
        status["r_available"] = True
        status["r_version"] = version
        if "under development" in version.lower() or "unstable" in version.lower():
            status["r_status"] = "devel"
            status["warnings"].append(
                "R appears to be a development/unstable build; Bioconductor release compatibility may differ."
            )
        else:
            status["r_status"] = "release"

    biocmanager_proc = _run_rscript('cat(requireNamespace("BiocManager", quietly = TRUE))')
    if biocmanager_proc and biocmanager_proc.returncode == 0:
        biocmanager_installed = _parse_bool_line(biocmanager_proc.stdout) is True
        status["biocmanager_installed"] = biocmanager_installed

    if status["biocmanager_installed"]:
        version_expr = (
            'cat("biocmanager_version=", as.character(packageVersion("BiocManager")), "\\n", sep="");'
            'cat("bioconductor_version=", as.character(BiocManager::version()), "\\n", sep="")'
        )
        version_proc = _run_rscript(version_expr)
        if version_proc and version_proc.returncode == 0:
            parsed = _parse_key_value_lines(version_proc.stdout)
            status["biocmanager_version"] = parsed.get("biocmanager_version", "")
            status["bioconductor_version"] = parsed.get("bioconductor_version", "")

        valid_expr = (
            'x <- tryCatch(BiocManager::valid(), error=function(e) e);'
            'if (inherits(x, "error")) {'
            'cat("remote_validation=failed\\n");'
            'cat("remote_validation_message=", conditionMessage(x), "\\n", sep="")'
            '} else {'
            'cat("remote_validation=ok\\n")'
            '}'
        )
        valid_proc = _run_rscript(valid_expr, timeout=300)
        if valid_proc and valid_proc.returncode == 0:
            parsed = _parse_key_value_lines(valid_proc.stdout)
            state = parsed.get("remote_validation", "unknown")
            message = parsed.get("remote_validation_message", "")
            if state == "ok":
                status["remote_validation"] = "validated live with BiocManager::valid()"
            elif message:
                status["remote_validation"] = f"live validation failed: {message}"
            else:
                status["remote_validation"] = "live validation failed"
        elif valid_proc:
            status["remote_validation"] = valid_proc.stderr.strip() or "live validation failed"

    if target_packages and status["r_available"]:
        if target_packages:
            pkg_vec = _r_quote_vector(target_packages)
            pkg_expr = (
                f'pkgs <- c({pkg_vec});'
                'for (pkg in pkgs) {'
                'cat(pkg, "=", requireNamespace(pkg, quietly = TRUE), "\\n", sep="")'
                '}'
            )
            package_proc = _run_rscript(pkg_expr)
            if package_proc and package_proc.returncode == 0:
                pkg_status = _parse_key_value_lines(package_proc.stdout)
                status["installed_packages"] = {
                    package: pkg_status.get(package, "FALSE") == "TRUE"
                    for package in target_packages
                }

    return status


def build_install_script(packages: list[str]) -> str:
    """Generate the reproducible R install script."""
    pkg_vec = f"c({_r_quote_vector(packages)})" if packages else "character()"
    return f"""# Generated by bioconductor-bridge
if (!requireNamespace("BiocManager", quietly = TRUE)) {{
  install.packages("BiocManager")
}}

packages <- {pkg_vec}
if (length(packages) > 0) {{
  BiocManager::install(packages, ask = FALSE, update = FALSE)
}} else {{
  message("No target packages were selected for installation.")
}}

BiocManager::valid()

writeLines(capture.output(sessionInfo()), "sessionInfo.txt")
"""


def build_install_expr(packages: list[str]) -> str:
    pkg_vec = _r_quote_vector(packages)
    return (
        'if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager"); '
        f'BiocManager::install(c({pkg_vec}), ask = FALSE, update = FALSE); '
        'writeLines(capture.output(sessionInfo()), stdout())'
    )


def execute_install(packages: list[str]) -> dict[str, Any]:
    """Execute an explicit BiocManager install."""
    expr = build_install_expr(packages)
    cmd = [_rscript_path() or "Rscript", "--vanilla", "-e", expr]

    if not _rscript_path():
        return {
            "requested_packages": packages,
            "executed": False,
            "command": _shell_command_for_expr(expr),
            "returncode": 1,
            "stdout": "",
            "stderr": "Rscript was not found on PATH.",
            "status": "failed",
        }

    proc = run_process(
        cmd,
        capture_output=True,
        text=True,
        timeout=3600,
        cwd=str(SKILL_DIR),
    )
    return {
        "requested_packages": packages,
        "executed": True,
        "command": _shell_command_for_expr(expr),
        "returncode": proc.returncode,
        "stdout": proc.stdout,
        "stderr": proc.stderr,
        "status": "installed" if proc.returncode == 0 else "failed",
    }


def capture_session_info() -> str:
    """Capture sessionInfo() when R is available."""
    proc = _run_rscript('writeLines(capture.output(sessionInfo()), stdout())')
    if proc is None:
        return "Rscript not available; sessionInfo() could not be captured.\n"
    if proc.returncode != 0:
        return (
            "sessionInfo() could not be captured.\n"
            f"stdout:\n{proc.stdout}\n\nstderr:\n{proc.stderr}\n"
        )
    return proc.stdout


def build_environment_yml(
    setup_status: dict[str, Any],
    target_packages: list[str],
    detected_domain: str | None,
) -> str:
    """Write a small environment manifest for reproducibility."""
    python_version = ".".join(str(value) for value in sys.version_info[:3])
    lines = [
        "name: bioconductor-bridge",
        "channels: []",
        "dependencies:",
        f"  - python={python_version}",
    ]
    if setup_status.get("r_version"):
        lines.append(f'  - r-base="{setup_status["r_version"]}"')
    if setup_status.get("bioconductor_version"):
        lines.append(f'  - bioconductor="{setup_status["bioconductor_version"]}"')
    if setup_status.get("biocmanager_version"):
        lines.append(f'  - biocmanager="{setup_status["biocmanager_version"]}"')
    lines.append("variables:")
    lines.append(f'  detected_domain: "{detected_domain or ""}"')
    lines.append(f'  requested_packages: "{",".join(target_packages)}"')
    lines.append(
        f'  remote_validation: "{setup_status.get("remote_validation", "").replace(chr(34), chr(39))}"'
    )
    return "\n".join(lines) + "\n"


def write_checksums(output_dir: Path, files: list[Path]) -> None:
    repro_dir = output_dir / "reproducibility"
    lines = []
    for file_path in files:
        if file_path.exists():
            lines.append(f"{sha256_file(file_path)}  {file_path.relative_to(output_dir)}")
    (repro_dir / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_recommendation_table(output_dir: Path, recommendations: list[dict[str, Any]]) -> Path | None:
    if not recommendations:
        return None
    tables_dir = output_dir / "tables"
    tables_dir.mkdir(exist_ok=True)
    table_path = tables_dir / "recommended_packages.csv"
    with table_path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(
            handle,
            fieldnames=[
                "name",
                "title",
                "domains",
                "workflow_role",
                "score",
                "documentation_score",
                "containers",
                "input_formats",
                "official_url",
                "explanation",
            ],
        )
        writer.writeheader()
        for package in recommendations:
            writer.writerow(
                {
                    "name": package.get("name", ""),
                    "title": package.get("title", ""),
                    "domains": ";".join(package.get("domains", [])),
                    "workflow_role": package.get("workflow_role", ""),
                    "score": package.get("score", ""),
                    "documentation_score": package.get("documentation_score", ""),
                    "containers": ";".join(package.get("containers", [])),
                    "input_formats": ";".join(package.get("input_formats", [])),
                    "official_url": package.get("official_url", ""),
                    "explanation": package.get("explanation", ""),
                }
            )
    return table_path


def collect_official_links(
    recommendations: list[dict[str, Any]],
    workflow: dict[str, Any] | None,
    documentation: dict[str, Any] | None = None,
) -> list[str]:
    links: list[str] = []
    for package in recommendations:
        url = package.get("official_url")
        if url and url not in links:
            links.append(url)
        for link in package.get("documentation_links", []):
            doc_url = link.get("url")
            if doc_url and doc_url not in links:
                links.append(doc_url)
    if workflow:
        for step in workflow.get("steps", []):
            url = step.get("official_url")
            if url and url not in links:
                links.append(url)
    if documentation:
        for link in documentation.get("documentation_links", []):
            doc_url = link.get("url")
            if doc_url and doc_url not in links:
                links.append(doc_url)
    return links


def build_report(
    mode: str,
    query: str,
    detected_domain: str | None,
    input_context: dict[str, str | None],
    recommendations: list[dict[str, Any]],
    workflow: dict[str, Any] | None,
    documentation: dict[str, Any] | None,
    setup_status: dict[str, Any],
    install_status: dict[str, Any] | None,
    domains: list[dict[str, Any]] | None,
    output_dir: Path,
) -> str:
    title_map = {
        "search": "Bioconductor Search Report",
        "recommend": "Bioconductor Recommendation Report",
        "workflow": "Bioconductor Workflow Report",
        "package_details": "Bioconductor Package Detail Report",
        "docs_search": "Bioconductor Documentation Search Report",
        "package_docs": "Bioconductor Package Documentation Report",
        "list_domains": "Bioconductor Supported Domains",
        "setup": "Bioconductor Setup Inspection",
        "install": "Bioconductor Installation Report",
        "demo": "Bioconductor Demo Report",
    }
    header = generate_report_header(
        title_map.get(mode, "Bioconductor Bridge Report"),
        "bioconductor-bridge",
        input_files=[Path(input_context["input_path"])] if input_context.get("input_path") else None,
        extra_metadata={
            "Mode": mode,
            "Detected domain": DOMAIN_DEFINITIONS[detected_domain]["display_name"] if detected_domain else "none",
        },
    )

    lines = [header]
    lines.append("## Query")
    lines.append(f"- `{query}`")
    lines.append("")

    lines.append("## Context")
    lines.append(f"- Input format: `{input_context.get('input_format') or 'none'}`")
    lines.append(f"- Modality: `{input_context.get('modality') or 'none'}`")
    lines.append(f"- Container: `{input_context.get('container') or 'none'}`")
    lines.append("")

    if recommendations:
        lines.append("## Recommended Packages")
        for package in recommendations:
            lines.append(
                f"- `{package['name']}`: {package.get('title', '')} "
                f"({package.get('workflow_role', 'supporting')}) — {package.get('explanation', '')}"
            )
            if package.get("documentation_summary"):
                lines.append(
                    f"  Docs: {' | '.join(package.get('documentation_summary', [])[:3])}"
                )
        lines.append("")

    if documentation:
        lines.append("## Documentation Snapshot")
        lines.append(f"- Package page: `{documentation.get('official_url', '')}`")
        if documentation.get("subtitle"):
            lines.append(f"- Subtitle: {documentation['subtitle']}")
        if documentation.get("documentation_summary"):
            lines.append("- Key documentation lines:")
            for line in documentation["documentation_summary"][:6]:
                lines.append(f"  - {line}")
        if documentation.get("documentation_links"):
            lines.append("- Documentation links:")
            for link in documentation["documentation_links"][:8]:
                lines.append(f"  - {link.get('label', 'link')}: {link.get('url', '')}")
        lines.append("")

    if workflow:
        lines.append("## Suggested Workflow")
        lines.append(
            f"- `{workflow['name']}` using `{workflow['container']}` as the canonical object model."
        )
        for step in workflow.get("steps", []):
            lines.append(
                f"- `{step['package']}` ({step['role']}): {step['description']}"
            )
        lines.append("")

    if domains:
        lines.append("## Supported Domains")
        for domain in domains:
            reps = ", ".join(f"`{pkg}`" for pkg in domain["representative_packages"])
            lines.append(f"- `{domain['display_name']}`: {reps}")
        lines.append("")

    lines.append("## Setup Status")
    lines.append(f"- R available: `{setup_status.get('r_available')}`")
    lines.append(f"- R version: `{setup_status.get('r_version') or 'not detected'}`")
    lines.append(f"- R channel: `{setup_status.get('r_status')}`")
    lines.append(f"- BiocManager installed: `{setup_status.get('biocmanager_installed')}`")
    lines.append(f"- BiocManager version: `{setup_status.get('biocmanager_version') or 'not detected'}`")
    lines.append(f"- Bioconductor version: `{setup_status.get('bioconductor_version') or 'not detected'}`")
    lines.append(f"- Remote validation: {setup_status.get('remote_validation')}")
    if setup_status.get("installed_packages"):
        installed_summary = ", ".join(
            f"{name}={'yes' if installed else 'no'}"
            for name, installed in setup_status["installed_packages"].items()
        )
        lines.append(f"- Target packages installed: {installed_summary}")
    for warning in setup_status.get("warnings", []):
        lines.append(f"- Warning: {warning}")
    lines.append("")

    if install_status:
        lines.append("## Installation Status")
        lines.append(f"- Status: `{install_status.get('status')}`")
        lines.append(f"- Command: `{install_status.get('command')}`")
        if install_status.get("stderr"):
            lines.append("- stderr captured during install execution.")
        lines.append("")

    lines.append("## Reproducibility")
    lines.append("- `reproducibility/commands.sh`")
    lines.append("- `reproducibility/environment.yml`")
    lines.append("- `reproducibility/install_packages.R`")
    lines.append("- `reproducibility/starter_workflow.R`")
    lines.append("- `reproducibility/sessionInfo.txt`")
    lines.append("- `reproducibility/checksums.sha256`")
    lines.append("")
    lines.append(f"Full report directory: `{output_dir}`")
    lines.append(generate_report_footer())
    return "\n".join(lines)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Bioconductor Bridge — live Bioconductor package discovery, workflow recommendation, and setup."
    )
    parser.add_argument("--input", help="Optional input path used for file-aware recommendations")
    parser.add_argument("--output", help="Output directory", default=str(DEFAULT_OUTPUT))
    parser.add_argument("--demo", action="store_true", help="Run the built-in demo recommendation workflow")
    parser.add_argument("--search", metavar="QUERY", help="Search live Bioconductor metadata with documentation-aware reranking")
    parser.add_argument("--recommend", metavar="TASK", help="Recommend Bioconductor packages for a task using metadata plus docs")
    parser.add_argument("--workflow", metavar="TASK", help="Suggest a fixed Bioconductor workflow template")
    parser.add_argument("--package-details", metavar="PACKAGE", help="Show details for a live Bioconductor package")
    parser.add_argument("--docs-search", metavar="QUERY", help="Search package documentation and vignette text for a task")
    parser.add_argument("--package-docs", metavar="PACKAGE", help="Fetch documentation summary and links for a package")
    parser.add_argument("--list-domains", action="store_true", help="List supported domains and representative packages")
    parser.add_argument("--setup", action="store_true", help="Inspect local Bioconductor readiness and write setup artifacts")
    parser.add_argument("--install", metavar="PKG1,PKG2", help="Explicitly install Bioconductor packages with BiocManager")
    parser.add_argument("--format", dest="input_format", metavar="EXT", help="Input format hint such as .vcf or .mtx")
    parser.add_argument("--modality", metavar="NAME", help="Modality hint such as bulk-rnaseq or single-cell")
    parser.add_argument("--container", metavar="NAME", help="Container hint such as SummarizedExperiment or GRanges")
    parser.add_argument("--max-results", type=int, default=5, help="Maximum number of search or recommendation results")
    return parser


def validate_args(args: argparse.Namespace) -> list[str]:
    primary_modes = {
        "demo": args.demo,
        "search": bool(args.search),
        "recommend": bool(args.recommend),
        "workflow": bool(args.workflow),
        "package_details": bool(args.package_details),
        "docs_search": bool(args.docs_search),
        "package_docs": bool(args.package_docs),
        "list_domains": bool(args.list_domains),
        "setup": bool(args.setup),
        "install": bool(args.install),
    }
    active = [name for name, enabled in primary_modes.items() if enabled]
    if len(active) != 1:
        raise SystemExit(
            "ERROR: Choose exactly one primary mode: --demo, --search, --recommend, "
            "--workflow, --package-details, --docs-search, --package-docs, --list-domains, --setup, or --install."
        )

    if args.install is not None:
        packages = [pkg.strip() for pkg in args.install.split(",") if pkg.strip()]
        if not packages:
            raise SystemExit("ERROR: --install requires a comma-separated package list.")

    return active


def mode_query(args: argparse.Namespace) -> str:
    if args.demo:
        return "bulk RNA-seq differential expression from a count matrix"
    if args.search:
        return args.search
    if args.recommend:
        return args.recommend
    if args.workflow:
        return args.workflow
    if args.package_details:
        return args.package_details
    if args.docs_search:
        return args.docs_search
    if args.package_docs:
        return args.package_docs
    if args.list_domains:
        return "supported bioconductor domains"
    if args.setup:
        return "local bioconductor setup inspection"
    return f"install bioconductor packages: {args.install}"


def collect_target_packages(
    recommendations: list[dict[str, Any]],
    workflow: dict[str, Any] | None,
    install_packages: list[str],
) -> list[str]:
    if install_packages:
        return install_packages
    if workflow:
        return [step["package"] for step in workflow.get("steps", [])]
    return [package["name"] for package in recommendations]


def resolve_input_context(args: argparse.Namespace) -> dict[str, str | None]:
    inferred = infer_input_context(args.input)
    return {
        "input_path": args.input,
        "input_format": args.input_format or inferred.get("input_format"),
        "domain": inferred.get("domain"),
        "modality": args.modality or inferred.get("modality"),
        "container": args.container or inferred.get("container"),
    }


def execute_mode(
    args: argparse.Namespace,
    catalog: dict[str, Any],
    input_context: dict[str, str | None],
) -> dict[str, Any]:
    query = mode_query(args)
    input_format = input_context.get("input_format")
    modality = input_context.get("modality")
    container = input_context.get("container")

    mode = (
        "demo" if args.demo else
        "search" if args.search else
        "recommend" if args.recommend else
        "workflow" if args.workflow else
        "package_details" if args.package_details else
        "docs_search" if args.docs_search else
        "package_docs" if args.package_docs else
        "list_domains" if args.list_domains else
        "setup" if args.setup else
        "install"
    )

    detected_domain = detect_domain(query, input_format=input_format, modality=modality, container=container)
    recommendations: list[dict[str, Any]] = []
    workflow: dict[str, Any] | None = None
    documentation: dict[str, Any] | None = None
    domains: list[dict[str, Any]] | None = None
    install_status: dict[str, Any] | None = None
    install_packages: list[str] = []

    if mode == "demo":
        recommendations = recommend_packages_with_docs(
            task=query,
            catalog=catalog,
            input_format=".counts_matrix",
            modality="bulk-rnaseq",
            container="SummarizedExperiment",
            max_results=args.max_results,
        )
        workflow = suggest_workflow(
            task=query,
            catalog=catalog,
            input_format=".counts_matrix",
            modality="bulk-rnaseq",
            container="SummarizedExperiment",
        )
        detected_domain = "bulk_rnaseq"

    elif mode == "search":
        recommendations = search_catalog_with_docs(query, catalog, max_results=args.max_results)
        if not detected_domain and recommendations:
            first_domains = recommendations[0].get("domains", [])
            if first_domains:
                detected_domain = first_domains[0]
        workflow = suggest_workflow(
            task=query,
            catalog=catalog,
            input_format=input_format,
            modality=modality,
            container=container,
        )

    elif mode == "recommend":
        recommendations = recommend_packages_with_docs(
            task=query,
            catalog=catalog,
            input_format=input_format,
            modality=modality,
            container=container,
            max_results=args.max_results,
        )
        if recommendations and not detected_domain:
            detected_domain = recommendations[0].get("detected_domain")
        workflow = suggest_workflow(
            task=query,
            catalog=catalog,
            input_format=input_format,
            modality=modality,
            container=container,
        )

    elif mode == "workflow":
        workflow = suggest_workflow(
            task=query,
            catalog=catalog,
            input_format=input_format,
            modality=modality,
            container=container,
        )
        if workflow is None:
            raise SystemExit(f"ERROR: Could not map workflow request to a supported Bioconductor domain: {query}")
        detected_domain = workflow["domain"]
        for package_name in [step["package"] for step in workflow["steps"]]:
            details = get_package_details(package_name, catalog)
            if details:
                details = dict(details)
                details["score"] = details.get("curated_priority", 0)
                details["explanation"] = f"included in the fixed {workflow['display_name']} workflow"
                recommendations.append(details)

    elif mode == "package_details":
        details = get_package_details(query, catalog)
        if details is None:
            raise SystemExit(f"ERROR: Package not found in live Bioconductor metadata: {query}")
        details = dict(details)
        details["score"] = details.get("curated_priority", 0)
        details["explanation"] = "package details requested explicitly"
        documentation = fetch_package_documentation(details)
        if documentation:
            details["documentation_summary"] = documentation.get("documentation_summary", [])
            details["documentation_links"] = documentation.get("documentation_links", [])
            details["documentation_subtitle"] = documentation.get("subtitle", "")
        recommendations = [details]
        if details.get("domains"):
            detected_domain = details["domains"][0]
            workflow = suggest_workflow(
                task=detected_domain.replace("_", " "),
                catalog=catalog,
                input_format=input_format,
                modality=modality,
                container=container,
            )

    elif mode == "docs_search":
        recommendations = docs_search_catalog(query, catalog, max_results=args.max_results)
        if not detected_domain and recommendations:
            first_domains = recommendations[0].get("domains", [])
            if first_domains:
                detected_domain = first_domains[0]
        workflow = suggest_workflow(
            task=query,
            catalog=catalog,
            input_format=input_format,
            modality=modality,
            container=container,
        )

    elif mode == "package_docs":
        details = get_package_details(query, catalog)
        if details is None:
            raise SystemExit(f"ERROR: Package not found in live Bioconductor metadata: {query}")
        details = dict(details)
        details["score"] = details.get("curated_priority", 0)
        details["explanation"] = "package documentation requested explicitly"
        documentation = fetch_package_documentation(details)
        if documentation:
            details["documentation_summary"] = documentation.get("documentation_summary", [])
            details["documentation_links"] = documentation.get("documentation_links", [])
            details["documentation_subtitle"] = documentation.get("subtitle", "")
        recommendations = [details]
        if details.get("domains"):
            detected_domain = details["domains"][0]
            workflow = suggest_workflow(
                task=detected_domain.replace("_", " "),
                catalog=catalog,
                input_format=input_format,
                modality=modality,
                container=container,
            )

    elif mode == "list_domains":
        domains = list_domains()

    elif mode == "setup":
        workflow = suggest_workflow(
            task=" ".join(
                value for value in [query, input_format or "", modality or "", container or ""] if value
            ),
            catalog=catalog,
            input_format=input_format,
            modality=modality,
            container=container,
        )
        if workflow:
            detected_domain = workflow["domain"]
            for package_name in [step["package"] for step in workflow["steps"]]:
                details = get_package_details(package_name, catalog)
                if details:
                    details = dict(details)
                    details["score"] = details.get("curated_priority", 0)
                    details["explanation"] = "selected from the setup-oriented workflow template"
                    recommendations.append(details)

    elif mode == "install":
        install_packages = [pkg.strip() for pkg in args.install.split(",") if pkg.strip()]
        for package_name in install_packages:
            details = get_package_details(package_name, catalog)
            if details:
                details = dict(details)
                details["score"] = details.get("curated_priority", 0)
                details["explanation"] = "requested explicitly for installation"
                recommendations.append(details)
        if recommendations and not detected_domain:
            domain_counts: dict[str, int] = {}
            for package in recommendations:
                for domain in package.get("domains", []):
                    domain_counts[domain] = domain_counts.get(domain, 0) + 1
            if domain_counts:
                detected_domain = max(domain_counts.items(), key=lambda item: item[1])[0]
        if detected_domain:
            workflow = suggest_workflow(
                task=detected_domain.replace("_", " "),
                catalog=catalog,
                input_format=input_format,
                modality=modality,
                container=container,
            )
        install_status = execute_install(install_packages)

    target_packages = collect_target_packages(recommendations, workflow, install_packages)
    setup_status = inspect_setup(target_packages)

    return {
        "mode": mode,
        "query": query,
        "detected_domain": detected_domain,
        "recommendations": recommendations,
        "workflow": workflow,
        "documentation": documentation,
        "domains": domains,
        "install_status": install_status,
        "setup_status": setup_status,
        "target_packages": target_packages,
    }


def write_reproducibility_bundle(
    output_dir: Path,
    payload: dict[str, Any],
    starter_script: str,
    command_line: str,
    setup_status: dict[str, Any],
) -> list[Path]:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(exist_ok=True)

    install_script_path = repro_dir / "install_packages.R"
    starter_script_path = repro_dir / "starter_workflow.R"
    commands_path = repro_dir / "commands.sh"
    env_path = repro_dir / "environment.yml"
    session_path = repro_dir / "sessionInfo.txt"

    install_script_path.write_text(build_install_script(payload["target_packages"]), encoding="utf-8")
    starter_script_path.write_text(starter_script, encoding="utf-8")

    commands = [
        "#!/usr/bin/env bash",
        "set -euo pipefail",
        "",
        f"# Original invocation",
        command_line,
    ]
    if payload.get("target_packages"):
        install_expr = build_install_expr(payload["target_packages"])
        commands.extend(
            [
                "",
                "# Optional install command",
                _shell_command_for_expr(install_expr),
                "# Live validation",
                "Rscript --vanilla -e 'BiocManager::valid()'",
            ]
        )
    commands_path.write_text("\n".join(commands) + "\n", encoding="utf-8")
    os.chmod(commands_path, 0o755)

    env_path.write_text(
        build_environment_yml(
            setup_status=setup_status,
            target_packages=payload["target_packages"],
            detected_domain=payload["detected_domain"],
        ),
        encoding="utf-8",
    )
    session_path.write_text(capture_session_info(), encoding="utf-8")
    return [commands_path, env_path, install_script_path, starter_script_path, session_path]


def print_summary(payload: dict[str, Any], output_dir: Path) -> None:
    print(f"Bioconductor Bridge — {payload['mode']}")
    if payload.get("detected_domain"):
        print(f"Detected domain: {DOMAIN_DEFINITIONS[payload['detected_domain']]['display_name']}")
    if payload.get("recommendations"):
        print("Top packages:")
        for package in payload["recommendations"][:5]:
            print(f"  - {package['name']}: {package.get('explanation', '')}")
    if payload.get("install_status"):
        print(f"Install status: {payload['install_status']['status']}")
    print(f"Report: {output_dir / 'report.md'}")


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    validate_args(args)

    output_dir = Path(args.output).expanduser().resolve()
    ensure_output_dir(output_dir)

    catalog = load_catalog()
    input_context = resolve_input_context(args)
    payload = execute_mode(args, catalog, input_context)

    starter_script = render_starter_script(payload.get("workflow"))
    table_path = write_recommendation_table(output_dir, payload["recommendations"])
    links = collect_official_links(
        payload["recommendations"],
        payload.get("workflow"),
        payload.get("documentation"),
    )
    report_md = build_report(
        mode=payload["mode"],
        query=payload["query"],
        detected_domain=payload.get("detected_domain"),
        input_context=input_context,
        recommendations=payload["recommendations"],
        workflow=payload.get("workflow"),
        documentation=payload.get("documentation"),
        setup_status=payload["setup_status"],
        install_status=payload.get("install_status"),
        domains=payload.get("domains"),
        output_dir=output_dir,
    )
    report_path = output_dir / "report.md"
    report_path.write_text(report_md, encoding="utf-8")

    summary = {
        "mode": payload["mode"],
        "detected_domain": payload.get("detected_domain", ""),
        "recommendation_count": len(payload["recommendations"]),
        "workflow_name": payload["workflow"]["name"] if payload.get("workflow") else "",
        "install_status": payload["install_status"]["status"] if payload.get("install_status") else "",
    }
    data = {
        "query": payload["query"],
        "detected_domain": payload.get("detected_domain"),
        "input_context": input_context,
        "recommendations": payload["recommendations"],
        "workflow": payload.get("workflow"),
        "documentation": payload.get("documentation"),
        "supported_domains": payload.get("domains"),
        "setup_status": payload["setup_status"],
        "installation_status": payload.get("install_status"),
        "official_links": links,
    }
    result_path = write_result_json(
        output_dir=output_dir,
        skill="bioconductor-bridge",
        version=VERSION,
        summary=summary,
        data=data,
    )

    repro_files = write_reproducibility_bundle(
        output_dir=output_dir,
        payload=payload,
        starter_script=starter_script,
        command_line=shlex.join([sys.executable, str(Path(__file__).resolve()), *sys.argv[1:]]),
        setup_status=payload["setup_status"],
    )
    files_to_checksum = [report_path, result_path, *repro_files]
    if table_path:
        files_to_checksum.append(table_path)
    write_checksums(output_dir, files_to_checksum)

    print_summary(payload, output_dir)

    install_status = payload.get("install_status")
    if install_status and install_status.get("returncode", 0) != 0:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
