#!/usr/bin/env python
"""Query the ucscxenatoolspy TCGA Analysis API — ClawBio skill.

Usage:
    python query_tcga_api.py cancers
    python query_tcga_api.py diff-expr --gene TP53 --cancer LUAD
    python query_tcga_api.py corr --gene TP53 --gene2 EGFR --cancer LUAD
    python query_tcga_api.py survival --gene TP53 --cancer LUAD
    python query_tcga_api.py --demo --output /tmp/xena_demo
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import sys
import time
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)

SCRIPT_REL = Path(__file__).resolve().relative_to(_PROJECT_ROOT)

DEFAULT_BASE_URL = os.getenv("UCSCXENA_API_BASE_URL", "http://biotree.top:38123/ucscxena/")
DEFAULT_API_KEY = os.getenv("UCSCXENA_API_KEY")

# Map hyphenated CLI task names to underscore keys used internally
TASK_KEY_MAP = {
    "cancers": "cancers",
    "diff-expr": "diff_expr",
    "corr": "corr",
    "survival": "survival",
}

CLAWBIO_DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and "
    "does not provide clinical diagnoses. Consult a healthcare professional before "
    "making any medical decisions."
)

PROVENANCE_NOTE = (
    "Results are computed by an author-hosted UCSCXenaToolsPy API service that "
    "queries public TCGA/UCSC Xena-derived datasets. Numerical statistics (fold "
    "change, p-values, survival associations) are calculated server-side and depend "
    "on the hosted service implementation and dataset version."
)


# ──────────────────────────────────────────────
# API helpers
# ──────────────────────────────────────────────

def _request(base_url: str, path: str, params: dict[str, str], api_key: str | None, timeout: int) -> Any:
    url = base_url.rstrip("/") + path
    if params:
        url += "?" + urllib.parse.urlencode(params)
    headers = {"Accept": "application/json"}
    if api_key:
        headers["X-API-Key"] = api_key
    req = urllib.request.Request(url, headers=headers, method="GET")
    try:
        with urllib.request.urlopen(req, timeout=timeout) as resp:
            text = resp.read().decode("utf-8")
    except urllib.error.HTTPError as exc:
        body = exc.read().decode("utf-8", errors="replace")
        raise SystemExit(f"HTTP {exc.code} from {url}: {body}") from exc
    except urllib.error.URLError as exc:
        raise SystemExit(f"Could not reach {url}: {exc.reason}") from exc
    return json.loads(text)


def _p(value: Any, digits: int = 4) -> str:
    if value is None:
        return "NA"
    if isinstance(value, float):
        if value != 0 and abs(value) < 0.001:
            return f"{value:.2e}"
        return f"{value:.{digits}g}"
    return str(value)


# ──────────────────────────────────────────────
# Summarise helpers
# ──────────────────────────────────────────────

def summarize(task: str, data: Any) -> str:
    if task == "cancers":
        cancers = data.get("cancers", [])
        with_normals = sum(1 for item in cancers if item.get("has_normal"))
        examples = ", ".join(item.get("cancer", "?") for item in cancers[:8])
        return (
            f"Available cancer entries: {data.get('count', len(cancers))}; "
            f"with >=3 normals: {with_normals}. Examples: {examples}."
        )

    if task == "diff-expr":
        gene = data.get("gene")
        if data.get("gene_input"):
            gene = f"{data['gene_input']} -> {gene}"
        tumor = data.get("tumor", {})
        normal = data.get("normal", {})
        return (
            f"{gene} in {data.get('cancer')} ({data.get('cancer_full_name')}): "
            f"tumor n={tumor.get('n')}, normal n={normal.get('n')}, "
            f"log2-scale difference={_p(data.get('log2_fold_change'))}, "
            f"Mann-Whitney p={_p(data.get('p_value'))}."
        )

    if task == "corr":
        g1 = data.get("gene1")
        g2 = data.get("gene2")
        if data.get("gene1_input"):
            g1 = f"{data['gene1_input']} -> {g1}"
        if data.get("gene2_input"):
            g2 = f"{data['gene2_input']} -> {g2}"
        return (
            f"{g1} vs {g2} in {data.get('cancer')} primary tumors: "
            f"n={data.get('n')}, Spearman r={_p(data.get('spearman_r'))}, "
            f"p={_p(data.get('p_value'))}."
        )

    if task == "survival":
        lines = [
            f"{data.get('gene')} in {data.get('cancer')} ({data.get('cancer_full_name')}), "
            f"log-rank survival associations:"
        ]
        for endpoint, result in data.get("survival", {}).items():
            if "error" in result:
                lines.append(f"- {endpoint}: {result['error']}")
                continue
            median = result.get("median_cutoff", {})
            optimal = result.get("optimal_cutoff", {})
            lines.append(
                f"- {endpoint}: n={result.get('n_total')}, events={result.get('n_events')}; "
                f"median p={_p(median.get('p_value'))}; "
                f"exploratory optimal p={_p(optimal.get('p_value'))}."
            )
        lines.append(
            "Optimal-cutoff results are exploratory and not adjusted for multiple cutoff testing."
        )
        return "\n".join(lines)

    return json.dumps(data, indent=2, ensure_ascii=False)


# ──────────────────────────────────────────────
# Demo data
# ──────────────────────────────────────────────

DEMO_DIFF_EXPR = {
    "gene": "TP53",
    "gene_input": None,
    "cancer": "LUAD",
    "cancer_full_name": "Lung Adenocarcinoma",
    "expression_scale": "log2(TPM + 0.001)",
    "test": "Mann-Whitney U (two-sided)",
    "tumor": {"n": 515, "mean": 5.12, "median": 5.08, "std": 0.89},
    "normal": {"n": 59, "mean": 4.87, "median": 4.91, "std": 0.76},
    "log2_fold_change": 0.25,
    "p_value": 0.0034,
}

DEMO_CORR = {
    "gene1": "TP53",
    "gene2": "EGFR",
    "gene1_input": None,
    "gene2_input": None,
    "cancer": "LUAD",
    "cancer_full_name": "Lung Adenocarcinoma",
    "expression_scale": "log2(TPM + 0.001)",
    "n": 508,
    "method": "Spearman rank correlation",
    "spearman_r": 0.18,
    "p_value": 4.2e-05,
}

DEMO_SURVIVAL = {
    "gene": "TP53",
    "cancer": "LUAD",
    "cancer_full_name": "Lung Adenocarcinoma",
    "expression_scale": "log2(TPM + 0.001)",
    "survival": {
        "OS": {
            "n_total": 504, "n_events": 189,
            "median_cutoff": {
                "p_value": 0.042, "high_n": 252, "low_n": 252,
                "high_events": 108, "low_events": 81,
                "high_mean_survival_days": 890, "low_mean_survival_days": 1120,
            },
            "optimal_cutoff": {
                "p_value": 0.0081, "cutoff": 5.31, "high_n": 178, "low_n": 326,
                "high_events": 82, "low_events": 107,
                "high_mean_survival_days": 720, "low_mean_survival_days": 1250,
            },
            "optimal_cutoff_note": "Exploratory. p-value is not adjusted for multiple cutoff testing.",
        },
        "DSS": {
            "n_total": 494, "n_events": 142,
            "median_cutoff": {
                "p_value": 0.11, "high_n": 247, "low_n": 247,
                "high_events": 79, "low_events": 63,
            },
            "optimal_cutoff": {
                "p_value": 0.021, "cutoff": 5.31,
            },
            "optimal_cutoff_note": "Exploratory. p-value is not adjusted for multiple cutoff testing.",
        },
        "DFI": {
            "n_total": 312, "n_events": 84,
            "median_cutoff": {"p_value": 0.67},
            "optimal_cutoff": {"p_value": 0.13},
            "optimal_cutoff_note": "Exploratory. p-value is not adjusted for multiple cutoff testing.",
        },
        "PFI": {
            "n_total": 504, "n_events": 218,
            "median_cutoff": {
                "p_value": 0.031, "high_n": 252, "low_n": 252,
                "high_events": 122, "low_events": 96,
            },
            "optimal_cutoff": {
                "p_value": 0.0056, "cutoff": 5.31,
            },
            "optimal_cutoff_note": "Exploratory. p-value is not adjusted for multiple cutoff testing.",
        },
    },
}


def run_demo(output_dir: Path) -> dict[str, Any]:
    """Return the demo results dict for report/JSON writing."""
    return {
        "diff_expr": DEMO_DIFF_EXPR,
        "corr": DEMO_CORR,
        "survival": DEMO_SURVIVAL,
    }


# ──────────────────────────────────────────────
# Report / output writers
# ──────────────────────────────────────────────

def build_report_md(demo_results: dict[str, Any], base_url: str, mode: str = "demo") -> str:
    """Build a Markdown report from demo results (three-task bundle)."""
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    mode_label = "demo (synthetic data — no live API calls)" if mode == "demo" else mode
    lines = [
        "# Xena TCGA Gene Query Report",
        "",
        f"**Date**: {now_utc}",
        f"**API base URL**: {base_url}",
        f"**Mode**: {mode_label}",
        "**Queries**: diff-expr (TP53 in LUAD), corr (TP53 vs EGFR in LUAD), survival (TP53 in LUAD)",
        "",
        f"*{PROVENANCE_NOTE}*",
        "",
        "---",
        "",
        "## 1. Differential Expression — TP53 in LUAD",
        "",
    ]

    de = demo_results["diff_expr"]
    lines.extend([
        f"**Gene**: {de['gene']}",
        f"**Cancer**: {de['cancer']} ({de['cancer_full_name']})",
        f"**Tumor samples**: n = {de['tumor']['n']}",
        f"**Normal samples**: n = {de['normal']['n']}",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Tumor mean (log2) | {_p(de['tumor']['mean'])} |",
        f"| Normal mean (log2) | {_p(de['normal']['mean'])} |",
        f"| log2 Fold Change | {_p(de['log2_fold_change'])} |",
        f"| Mann-Whitney p | {_p(de['p_value'])} |",
        "",
        f"**Interpretation**: TP53 expression is modestly higher in LUAD tumor vs normal "
        f"tissue (p = {_p(de['p_value'])}). The difference (~{_p(de['log2_fold_change'])} "
        f"log2 units) is statistically significant but biologically small.",
        "",
        "---",
        "",
        "## 2. Gene-Gene Correlation — TP53 vs EGFR in LUAD",
        "",
    ])

    cr = demo_results["corr"]
    lines.extend([
        f"**Genes**: {cr['gene1']}, {cr['gene2']}",
        f"**Cancer**: {cr['cancer']} ({cr['cancer_full_name']})",
        f"**Primary tumor samples**: n = {cr['n']}",
        "",
        "| Metric | Value |",
        "|--------|-------|",
        f"| Spearman r | {_p(cr['spearman_r'])} |",
        f"| p-value | {_p(cr['p_value'])} |",
        "",
        f"**Interpretation**: TP53 and EGFR show a weak positive rank correlation in "
        f"LUAD primary tumors (Spearman r = {_p(cr['spearman_r'])}, p = {_p(cr['p_value'])}). "
        f"The correlation is statistically detectable but explains little variance.",
        "",
        "---",
        "",
        "## 3. Survival Association — TP53 in LUAD",
        "",
    ])

    sv = demo_results["survival"]
    lines.extend([
        f"**Gene**: {sv['gene']}",
        f"**Cancer**: {sv['cancer']} ({sv['cancer_full_name']})",
        "",
        "| Endpoint | n (total) | Events | Median-cutoff p | Optimal-cutoff p (exploratory) |",
        "|----------|-----------|--------|-----------------|-------------------------------|",
    ])

    for ep_name in ["OS", "DSS", "DFI", "PFI"]:
        ep = sv["survival"].get(ep_name, {})
        median_p = _p(ep.get("median_cutoff", {}).get("p_value"))
        optimal_p = _p(ep.get("optimal_cutoff", {}).get("p_value"))
        lines.append(
            f"| {ep_name} | {ep.get('n_total', 'NA')} | {ep.get('n_events', 'NA')} | "
            f"{median_p} | {optimal_p} |"
        )

    lines.extend([
        "",
        "Optimal-cutoff results are exploratory and not adjusted for multiple cutoff testing.",
        "",
        "**Interpretation**: Higher TP53 expression is associated with worse overall survival "
        "(OS) and progression-free interval (PFI) at the median split. Disease-specific "
        "survival (DSS) and disease-free interval (DFI) do not reach significance at the "
        "median cutoff. These are statistical associations; they do not prove TP53 is a "
        "causal driver of outcome.",
        "",
        "---",
        "",
        f"*{CLAWBIO_DISCLAIMER}*",
    ])

    return "\n".join(lines) + "\n"


def build_live_report_md(task_key: str, data: dict[str, Any], base_url: str, mode: str = "live") -> str:
    """Build a Markdown report from a single live API query result."""
    now_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
    task_display = task_key.replace("_", "-")
    lines = [
        "# Xena TCGA Gene Query Report",
        "",
        f"**Date**: {now_utc}",
        f"**API base URL**: {base_url}",
        f"**Mode**: {mode}",
        f"**Task**: {task_display}",
        "",
        f"*{PROVENANCE_NOTE}*",
        "",
        "---",
        "",
        "## Results",
        "",
    ]

    if task_key == "cancers":
        cancers = data.get("cancers", [])
        with_normals = sum(1 for item in cancers if item.get("has_normal"))
        lines.append(f"**Available cancer entries**: {data.get('count', len(cancers))}")
        lines.append(f"**With >=3 normal samples**: {with_normals}")
        lines.append("")
        if cancers:
            lines.append("| Cancer | Full Name | Tumor (n) | Normal (n) |")
            lines.append("|--------|-----------|-----------|------------|")
            for item in cancers[:20]:
                lines.append(
                    f"| {item.get('cancer', '?')} | {item.get('cancer_full_name', '?')} | "
                    f"{item.get('tumor_n', '?')} | {item.get('normal_n', '?')} |"
                )

    elif task_key == "diff_expr":
        gene = data.get("gene", "?")
        if data.get("gene_input"):
            gene = f"{data['gene_input']} -> {gene}"
        tumor = data.get("tumor", {})
        normal = data.get("normal", {})
        lines.extend([
            f"**Gene**: {gene}",
            f"**Cancer**: {data.get('cancer')} ({data.get('cancer_full_name')})",
            f"**Tumor samples**: n = {tumor.get('n')}",
            f"**Normal samples**: n = {normal.get('n')}",
            "",
            "| Metric | Value |",
            "|--------|-------|",
            f"| Tumor mean (log2) | {_p(tumor.get('mean'))} |",
            f"| Normal mean (log2) | {_p(normal.get('mean'))} |",
            f"| log2 Fold Change | {_p(data.get('log2_fold_change'))} |",
            f"| Mann-Whitney p | {_p(data.get('p_value'))} |",
        ])

    elif task_key == "corr":
        g1 = data.get("gene1", "?")
        g2 = data.get("gene2", "?")
        if data.get("gene1_input"):
            g1 = f"{data['gene1_input']} -> {g1}"
        if data.get("gene2_input"):
            g2 = f"{data['gene2_input']} -> {g2}"
        lines.extend([
            f"**Genes**: {g1}, {g2}",
            f"**Cancer**: {data.get('cancer')} ({data.get('cancer_full_name')})",
            f"**Primary tumor samples**: n = {data.get('n')}",
            "",
            "| Metric | Value |",
            "|--------|-------|",
            f"| Spearman r | {_p(data.get('spearman_r'))} |",
            f"| p-value | {_p(data.get('p_value'))} |",
        ])

    elif task_key == "survival":
        lines.extend([
            f"**Gene**: {data.get('gene')}",
            f"**Cancer**: {data.get('cancer')} ({data.get('cancer_full_name')})",
            "",
            "| Endpoint | n (total) | Events | Median-cutoff p | Optimal-cutoff p (exploratory) |",
            "|----------|-----------|--------|-----------------|-------------------------------|",
        ])
        for ep_name in ["OS", "DSS", "DFI", "PFI"]:
            ep = data.get("survival", {}).get(ep_name, {})
            if "error" in ep:
                lines.append(f"| {ep_name} | — | — | {ep['error']} | — |")
            else:
                median_p = _p(ep.get("median_cutoff", {}).get("p_value"))
                optimal_p = _p(ep.get("optimal_cutoff", {}).get("p_value"))
                lines.append(
                    f"| {ep_name} | {ep.get('n_total', 'NA')} | {ep.get('n_events', 'NA')} | "
                    f"{median_p} | {optimal_p} |"
                )
        lines.append("")
        lines.append(
            "Optimal-cutoff results are exploratory and not adjusted for multiple cutoff testing."
        )

    lines.extend([
        "",
        "---",
        "",
        f"*{CLAWBIO_DISCLAIMER}*",
    ])
    return "\n".join(lines) + "\n"


def replay_args(args: argparse.Namespace, output_dir: Path) -> list[str | ReproPath]:
    """CLI arguments that reproduce this run. The API key is never recorded;
    replay reads it from UCSCXENA_API_KEY like the original run can."""
    out: list[str | ReproPath] = [
        "--base-url", shlex.quote(args.base_url),
        "--output", ReproPath(output_dir, anchor="output_dir"),
    ]
    if args.demo:
        return out + ["--demo"]
    out.append(args.task)
    for flag in ("gene", "gene2", "cancer"):
        value = getattr(args, flag, None)
        if value is not None:
            out += [f"--{flag}", shlex.quote(value)]
    return out


def write_outputs(
    results: dict[str, Any],
    base_url: str,
    output_dir: Path,
    mode: str = "demo",
    *,
    replay: list[str | ReproPath],
) -> None:
    """Write report.md, result.json, and reproducibility/ to output_dir."""
    output_dir.mkdir(parents=True, exist_ok=True)

    # report.md
    report_path = output_dir / "report.md"
    if mode == "demo":
        report_md = build_report_md(results, base_url, mode=mode)
    else:
        # Live mode: single-task result
        task_key = next(iter(results))
        report_md = build_live_report_md(task_key, results[task_key], base_url, mode=mode)
    report_path.write_text(report_md, encoding="utf-8")

    # result.json
    result_path = output_dir / "result.json"
    result_path.write_text(json.dumps(results, indent=2, ensure_ascii=False), encoding="utf-8")

    # reproducibility/
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    write_portable_commands_sh(
        output_dir,
        ReproCommand(
            script_path=SCRIPT_REL,
            args=replay,
            comment=f"Reproduce this xena-tcga-gene-query run ({mode} mode)",
        ),
        repo_root=_PROJECT_ROOT,
    )
    write_environment_yml(output_dir, "clawbio-xena-tcga-gene-query", pip_deps=[])

    run_json = repro_dir / "run.json"
    run_meta = {
        "skill": "xena-tcga-gene-query",
        "mode": mode,
        "base_url": base_url,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
        "python_version": sys.version,
    }
    run_json.write_text(json.dumps(run_meta, indent=2), encoding="utf-8")

    # Last, so the manifest covers every other output of this run.
    write_checksums(
        [report_path, result_path, repro_dir / "commands.sh", repro_dir / "environment.yml", run_json],
        output_dir,
        anchor=output_dir,
    )


# ──────────────────────────────────────────────
# CLI
# ──────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description=__doc__,
        formatter_class=argparse.RawDescriptionHelpFormatter,
    )
    parser.add_argument("--base-url", default=DEFAULT_BASE_URL)
    parser.add_argument("--api-key", default=DEFAULT_API_KEY)
    parser.add_argument("--timeout", type=int, default=120)
    parser.add_argument("--json", action="store_true", help="Print raw JSON only")
    parser.add_argument("--demo", action="store_true", help="Run with synthetic demo data (no API calls)")
    parser.add_argument("--output", "-o", type=str, default=None, help="Output directory for report.md + result.json + reproducibility/")

    sub = parser.add_subparsers(dest="task")
    sub.add_parser("cancers")

    diff = sub.add_parser("diff-expr")
    diff.add_argument("--gene", required=True)
    diff.add_argument("--cancer", required=True)

    corr = sub.add_parser("corr")
    corr.add_argument("--gene", required=True, help="First gene")
    corr.add_argument("--gene2", required=True, help="Second gene")
    corr.add_argument("--cancer", required=True)

    surv = sub.add_parser("survival")
    surv.add_argument("--gene", required=True)
    surv.add_argument("--cancer", required=True)
    return parser


def main(argv: list[str] | None = None) -> int:
    args = build_parser().parse_args(argv)

    # ── Demo mode ─────────────────────────────
    if args.demo:
        output_dir = Path(args.output) if args.output else Path("/tmp/xena_demo")
        results = run_demo(output_dir)
        write_outputs(results, args.base_url, output_dir, mode="demo", replay=replay_args(args, output_dir))

        # Also print a summary to stdout
        print(f"Demo complete. Output written to {output_dir}")
        print(f"  {output_dir / 'report.md'}")
        print(f"  {output_dir / 'result.json'}")
        print(f"  {output_dir / 'reproducibility'}/")
        print()
        print(summarize("diff-expr", results["diff_expr"]))
        print(summarize("corr", results["corr"]))
        print(summarize("survival", results["survival"]))
        return 0

    # ── Live API mode ─────────────────────────
    if not args.task:
        print("Error: specify a task (cancers, diff-expr, corr, survival) or use --demo", file=sys.stderr)
        return 1

    paths = {
        "cancers": "/api/v1/cancers",
        "diff-expr": "/api/v1/diff-expr",
        "corr": "/api/v1/corr",
        "survival": "/api/v1/survival",
    }
    params: dict[str, str] = {}
    if args.task == "diff-expr":
        params = {"gene": args.gene, "cancer": args.cancer}
    elif args.task == "corr":
        params = {"gene1": args.gene, "gene2": args.gene2, "cancer": args.cancer}
    elif args.task == "survival":
        params = {"gene": args.gene, "cancer": args.cancer}

    data = _request(args.base_url, paths[args.task], params, args.api_key, args.timeout)

    if args.json:
        print(json.dumps(data, indent=2, ensure_ascii=False))
    else:
        print(summarize(args.task, data))
        print("\nRaw JSON:")
        print(json.dumps(data, indent=2, ensure_ascii=False))

    # Write outputs if --output specified
    if args.output:
        output_dir = Path(args.output)
        task_key = TASK_KEY_MAP[args.task]
        results = {task_key: data}
        write_outputs(results, args.base_url, output_dir, mode="live", replay=replay_args(args, output_dir))

    return 0


if __name__ == "__main__":
    raise SystemExit(main())
