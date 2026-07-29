#!/usr/bin/env python3
"""ClawBio bridge to a pinned local just-prs MCP server."""

from __future__ import annotations

import asyncio
import csv
import json
import os
import sys
from importlib.metadata import version as distribution_version
from pathlib import Path
from statistics import mean
from typing import Any, Protocol

import typer

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file  # noqa: E402
from clawbio.common.report import (  # noqa: E402
    generate_report_footer,
    generate_report_header,
    write_result_json,
)
from clawbio.common.textio import write_text_lf  # noqa: E402
from mcp_client import JustPrsMcpClient, McpCallError  # noqa: E402

SKILL_NAME = "just-prs-mcp"
SKILL_VERSION = "0.1.1"
JUST_PRS_MCP_VERSION = "0.3.1"
_SKILL_DIR = Path(__file__).resolve().parent
_DEMO_REPORT = _SKILL_DIR / "tests" / "fixtures" / "demo_trait_report.json"
_DEMO_VCF = _SKILL_DIR / "examples" / "demo_patient.vcf"
_SUPERPOPULATIONS = {"AFR", "AMR", "EAS", "EUR", "SAS", "AUTO"}
_DEFAULT_SUPERPOPULATION = "EUR"
_CHILD_ENV_ALLOWLIST = {
    "PATH",
    "HOME",
    "USER",
    "LOGNAME",
    "TMPDIR",
    "TEMP",
    "TMP",
    "XDG_CACHE_HOME",
    "XDG_DATA_HOME",
    "LANG",
    "LC_ALL",
    "SYSTEMROOT",
    "WINDIR",
    "COMSPEC",
    "PATHEXT",
    "SSL_CERT_FILE",
    "SSL_CERT_DIR",
    "REQUESTS_CA_BUNDLE",
    "CURL_CA_BUNDLE",
    "HTTP_PROXY",
    "HTTPS_PROXY",
    "NO_PROXY",
    "http_proxy",
    "https_proxy",
    "no_proxy",
    "UV_CACHE_DIR",
    "UV_OFFLINE",
    "UV_NATIVE_TLS",
    "UV_HTTP_TIMEOUT",
    "UV_HTTP_RETRIES",
    "PRS_CACHE_DIR",
    "PRS_DUCKDB_MEMORY_LIMIT",
    "PRS_MCP_CACHE_DIR",
    "PRS_MCP_DEFAULT_GENOME_BUILD",
    "PRS_MCP_DUCKDB_MEMORY_LIMIT",
    "PRS_MCP_LOG_LEVEL",
}


class ToolClient(Protocol):
    async def call_tool(self, name: str, arguments: dict[str, Any]) -> Any: ...


class AmbiguousTraitError(ValueError):
    """Trait search returned multiple plausible ontology records."""


def normalize_superpopulation(value: str | None) -> tuple[str, bool]:
    """Resolve ancestry choice; omitted values default to EUR with disclosure."""
    if value is None or not str(value).strip():
        return _DEFAULT_SUPERPOPULATION, True
    resolved = str(value).strip().upper()
    if resolved not in _SUPERPOPULATIONS:
        raise ValueError(
            "--superpopulation must be AFR, AMR, EAS, EUR, SAS, or AUTO."
        )
    return resolved, False


def emit_superpopulation_default_warning(superpopulation: str) -> None:
    """Warn on stderr when EUR was applied because ancestry was omitted."""
    typer.echo(
        "WARNING: --superpopulation was omitted, so percentiles and absolute-risk "
        f"figures use the {_DEFAULT_SUPERPOPULATION} reference panel by default. "
        "Pass --superpopulation explicitly (AFR/AMR/EAS/EUR/SAS) or AUTO to avoid "
        f"silent {superpopulation}-referenced interpretation.",
        err=True,
    )


def server_command() -> list[str]:
    """Return the reproducibly pinned local stdio server command."""
    return [
        "uvx",
        f"just-prs-mcp@{JUST_PRS_MCP_VERSION}",
        "stdio",
        "--mode",
        "essentials",
    ]


def server_environment() -> dict[str, str]:
    """Return the minimum child environment without forwarding unrelated secrets."""
    environment = {
        name: value for name, value in os.environ.items() if name in _CHILD_ENV_ALLOWLIST
    }
    environment["PRS_MCP_MODE"] = "essentials"
    return environment


def resolve_trait(term: str, matches: list[dict[str, Any]]) -> str:
    """Resolve one unambiguous trait search result to its ontology ID."""
    if not matches:
        raise ValueError(f"No PGS Catalog trait matched '{term}'.")

    normalized = term.strip().casefold()
    exact = [
        match
        for match in matches
        if normalized
        in {
            str(match.get("id", "")).casefold(),
            str(match.get("label", "")).casefold(),
        }
    ]
    candidates = exact or matches
    if len(candidates) != 1:
        options = ", ".join(
            f"{match.get('id', 'unknown')} ({match.get('label', 'unlabelled')})"
            for match in candidates[:10]
        )
        raise AmbiguousTraitError(
            f"Trait '{term}' is ambiguous. Re-run with --trait-id using one of: {options}"
        )
    trait_id = str(candidates[0].get("id", "")).strip()
    if not trait_id:
        raise ValueError(f"Trait search for '{term}' returned a record without an ontology ID.")
    return trait_id


def _agreement(scores: list[dict[str, Any]]) -> dict[str, Any]:
    reliable = [
        float(score["percentile"])
        for score in scores
        if score.get("status") == "scored"
        and score.get("percentile_reliable") is True
        and score.get("percentile") is not None
    ]
    if len(reliable) < 2:
        return {
            "reliable_models": len(reliable),
            "verdict": "insufficient_reliable_models",
            "percentile_min": reliable[0] if reliable else None,
            "percentile_max": reliable[0] if reliable else None,
            "percentile_spread": 0.0 if reliable else None,
            "percentile_mean": reliable[0] if reliable else None,
        }
    return {
        "reliable_models": len(reliable),
        "verdict": "descriptive_spread_only",
        "percentile_min": min(reliable),
        "percentile_max": max(reliable),
        "percentile_spread": max(reliable) - min(reliable),
        "percentile_mean": mean(reliable),
    }


def _performance_evidence(performance: Any) -> tuple[str | None, float | None]:
    """Extract display-ready effect-size and AUROC evidence from PRSResult.performance."""
    if not isinstance(performance, dict):
        return None, None

    effect_size: str | None = None
    for metric in performance.get("effect_sizes", []):
        if not isinstance(metric, dict) or metric.get("estimate") is None:
            continue
        name = str(metric.get("name_short") or "effect")
        effect_size = f"{name}={float(metric['estimate']):.2f}"
        if metric.get("ci_lower") is not None and metric.get("ci_upper") is not None:
            effect_size += (
                f" [{float(metric['ci_lower']):.2f}-{float(metric['ci_upper']):.2f}]"
            )
        break

    auroc: float | None = None
    for metric in performance.get("class_acc", []):
        if not isinstance(metric, dict) or metric.get("estimate") is None:
            continue
        if str(metric.get("name_short", "")).casefold() in {"auc", "auroc"}:
            auroc = float(metric["estimate"])
            break
    return effect_size, auroc


def map_trait_report(
    report: dict[str, Any],
    *,
    superpopulation: str,
    absolute_risks: dict[str, dict[str, Any]],
    superpopulation_defaulted: bool = False,
) -> dict[str, Any]:
    """Map the upstream trait report without flattening away uncertainty."""
    scores: list[dict[str, Any]] = []
    for row in report.get("rows", []):
        pgs_id = str(row.get("pgs_id", ""))
        score = {
            "pgs_id": pgs_id,
            "pgs_catalog_url": (
                f"https://www.pgscatalog.org/score/{pgs_id}/" if pgs_id else None
            ),
            "status": row.get("status"),
            "raw_score": row.get("score"),
            "genome_build": row.get("genome_build"),
            "variants_number": row.get("variants_number"),
            "weight_type": row.get("weight_type"),
            "variants_matched": row.get("variants_matched"),
            "variants_total": row.get("variants_total"),
            "match_rate": row.get("match_rate"),
            "weight_mass_coverage": row.get("weight_mass_coverage"),
            "percentile": row.get("percentile"),
            "percentile_method": row.get("percentile_method"),
            "percentile_reliable": row.get("percentile_reliable"),
            "percentile_caveat": row.get("percentile_caveat"),
            "requested_superpopulation": superpopulation,
            "reference_panel_ancestry": row.get("reference_panel_ancestry"),
            "quality_label": row.get("quality_label"),
            "quality_summary": row.get("quality_summary"),
            "effect_size": row.get("effect_size"),
            "auroc_estimate": row.get("auroc_estimate"),
            "absolute_risk": absolute_risks.get(
                pgs_id,
                {
                    "status": "not_attempted" if row.get("status") != "scored" else "unavailable",
                    "reason": (
                        "Score failed, so absolute risk was not attempted."
                        if row.get("status") != "scored"
                        else "No absolute-risk result was returned."
                    ),
                },
            ),
            "error": row.get("error"),
        }
        scores.append(score)

    panel_ancestries = sorted(
        {
            str(score["reference_panel_ancestry"])
            for score in scores
            if score.get("reference_panel_ancestry")
        }
    )
    return {
        "trait_id": report.get("trait_id"),
        "trait": report.get("label"),
        "genome_build": report.get("genome_build"),
        "detected_genome_build": report.get("detected_genome_build"),
        "build_mismatch": bool(report.get("build_mismatch", False)),
        "profile": report.get("profile"),
        "requested_superpopulation": superpopulation,
        "superpopulation_defaulted": superpopulation_defaulted,
        "reference_panel_ancestry": (
            panel_ancestries[0]
            if len(panel_ancestries) == 1
            else (", ".join(panel_ancestries) if panel_ancestries else None)
        ),
        "n_requested": report.get("n_requested", 0),
        "n_scored": report.get("n_scored", 0),
        "n_failed": report.get("n_failed", 0),
        "n_skipped": report.get("n_skipped", 0),
        "n_filtered": report.get("n_filtered", 0),
        "n_omitted": report.get("n_omitted", 0),
        "filter_summary": report.get("filter_summary"),
        "upstream_summary": report.get("summary"),
        "model_agreement": _agreement(scores),
        "scores": scores,
    }


async def _absolute_risks_for_rows(
    client: ToolClient,
    rows: list[dict[str, Any]],
    *,
    superpopulation: str,
) -> dict[str, dict[str, Any]]:
    risks: dict[str, dict[str, Any]] = {}
    for row in rows:
        if row.get("status") != "scored" or row.get("score") is None:
            continue
        pgs_id = str(row["pgs_id"])
        panel = superpopulation
        if panel == "AUTO":
            panel = str(row.get("reference_panel_ancestry") or _DEFAULT_SUPERPOPULATION)
        try:
            percentile = await client.call_tool(
                "percentile",
                {
                    "prs_score": row["score"],
                    "pgs_id": pgs_id,
                    "superpopulation": panel,
                    "weight_mass_coverage": row.get("weight_mass_coverage"),
                },
            )
            row.update(
                {
                    "percentile": percentile.get("percentile"),
                    "percentile_method": percentile.get("method"),
                    "percentile_reliable": percentile.get("reliable"),
                    "percentile_caveat": percentile.get("caveat"),
                    "reference_panel_ancestry": percentile.get(
                        "reference_panel_ancestry"
                    ),
                }
            )
            if not row.get("quality_label"):
                quality = await client.call_tool(
                    "assess_quality",
                    {
                        "match_rate": row.get("match_rate") or 0.0,
                        "auroc": row.get("auroc_estimate"),
                        "percentile": percentile.get("percentile"),
                        "percentile_method": percentile.get("method"),
                        "reliable": percentile.get("reliable") is True,
                        "caveat": percentile.get("caveat") or "",
                    },
                )
                row["quality_label"] = quality.get("quality_label")
                row["quality_summary"] = quality.get("summary")
            if percentile.get("reliable") is not True:
                risks[pgs_id] = {
                    "status": "unavailable",
                    "reason": percentile.get("caveat")
                    or "The percentile was not marked reliable, so absolute risk was not computed.",
                }
                continue
            z_score = percentile.get("z_score")
            if z_score is None:
                raise McpCallError("Percentile result did not include a z-score.")
            risk = await client.call_tool(
                "absolute_risk",
                {"pgs_id": pgs_id, "z_score": z_score},
            )
            risks[pgs_id] = {"status": "available", **risk}
        except McpCallError as exc:
            risks[pgs_id] = {"status": "unavailable", "reason": str(exc)}
    return risks


async def analyze_with_client(
    *,
    client: ToolClient,
    input_path: Path,
    trait: str | None,
    trait_id: str | None,
    pgs_id: str | None,
    superpopulation: str,
    profile: str,
    top_n: int,
    limit: int | None,
    superpopulation_defaulted: bool = False,
) -> dict[str, Any]:
    """Run the live MCP workflow using only the input's resolved local path."""
    resolved_input = input_path.expanduser().resolve()
    if pgs_id:
        computed = await client.call_tool(
            "compute_prs",
            {
                "vcf_path": str(resolved_input),
                "pgs_id": pgs_id,
                "attach_performance": True,
            },
        )
        effect_size, auroc_estimate = _performance_evidence(computed.get("performance"))
        row = {
            "pgs_id": pgs_id,
            "status": "scored",
            "score": computed.get("score"),
            "genome_build": computed.get("genome_build"),
            "variants_number": computed.get("variants_total"),
            "variants_matched": computed.get("variants_matched"),
            "variants_total": computed.get("variants_total"),
            "match_rate": computed.get("match_rate"),
            "weight_mass_coverage": computed.get("weight_mass_coverage"),
            "percentile": computed.get("percentile"),
            "percentile_method": computed.get("percentile_method"),
            "percentile_reliable": None,
            "quality_label": None,
            "quality_summary": None,
            "effect_size": effect_size,
            "auroc_estimate": auroc_estimate,
            "error": None,
        }
        risks = await _absolute_risks_for_rows(
            client, [row], superpopulation=superpopulation
        )
        pseudo_report = {
            "trait_id": None,
            "label": computed.get("trait_reported") or pgs_id,
            "genome_build": computed.get("genome_build"),
            "detected_genome_build": computed.get("detected_genome_build"),
            "build_mismatch": computed.get("build_mismatch", False),
            "profile": "single-score",
            "n_requested": 1,
            "n_scored": 1,
            "n_failed": 0,
            "n_skipped": 0,
            "n_filtered": 0,
            "n_omitted": 0,
            "filter_summary": None,
            "summary": f"Computed {pgs_id}.",
            "rows": [row],
        }
        return map_trait_report(
            pseudo_report,
            superpopulation=superpopulation,
            absolute_risks=risks,
            superpopulation_defaulted=superpopulation_defaulted,
        )

    resolved_trait_id = trait_id
    if resolved_trait_id is None:
        matches = await client.call_tool(
            "search_traits",
            {"term": trait, "limit": 10, "include_pgs_ids": False},
        )
        if not isinstance(matches, list):
            raise McpCallError("search_traits returned an unexpected payload.")
        resolved_trait_id = resolve_trait(str(trait), matches)

    by_trait_arguments: dict[str, Any] = {
        "trait_id": resolved_trait_id,
        "vcf_path": str(resolved_input),
        "interpret": True,
        "superpopulation": superpopulation,
        "profile": profile,
        "top_n": top_n,
    }
    if limit is not None:
        by_trait_arguments["limit"] = limit
    report = await client.call_tool(
        "compute_prs_by_trait",
        by_trait_arguments,
    )
    if not isinstance(report, dict):
        raise McpCallError("compute_prs_by_trait returned an unexpected payload.")
    risks = await _absolute_risks_for_rows(
        client,
        list(report.get("rows", [])),
        superpopulation=superpopulation,
    )
    return map_trait_report(
        report,
        superpopulation=superpopulation,
        absolute_risks=risks,
        superpopulation_defaulted=superpopulation_defaulted,
    )


def _warn_before_overwrite(output_dir: Path) -> None:
    if output_dir.exists() and any(output_dir.iterdir()):
        typer.echo(
            f"WARNING: output directory contains files that may be overwritten: {output_dir}",
            err=True,
        )


def _format_fraction(value: Any) -> str:
    return "n/a" if value is None else f"{float(value):.1%}"


def _format_number(value: Any, digits: int = 2) -> str:
    return "n/a" if value is None else f"{float(value):.{digits}f}"


def _report_markdown(
    data: dict[str, Any],
    *,
    input_path: Path,
    demo: bool,
) -> str:
    agreement = data["model_agreement"]
    requested = str(data.get("requested_superpopulation") or _DEFAULT_SUPERPOPULATION)
    if data.get("superpopulation_defaulted"):
        requested_label = f"{requested} (default — not explicitly chosen)"
    else:
        requested_label = requested
    panel_ancestry = data.get("reference_panel_ancestry") or "unknown"
    header = generate_report_header(
        title="just-prs Polygenic Risk Report",
        skill_name=SKILL_NAME,
        skill_version=SKILL_VERSION,
        input_files=[input_path],
        extra_metadata={
            "Mode": "Synthetic offline demo" if demo else "Live local MCP",
            "Trait": str(data.get("trait") or "unknown"),
            "Trait ontology ID": str(data.get("trait_id") or "single score"),
            "Genome build": str(data.get("genome_build") or "unknown"),
            "Requested ancestry": requested_label,
            "Reference panel ancestry": str(panel_ancestry),
        },
    )
    lines = [header]
    if data.get("superpopulation_defaulted"):
        lines.extend(
            [
                "> **WARNING:** `--superpopulation` was omitted, so percentiles and "
                f"absolute-risk figures are **EUR-referenced** by default. Pass "
                "`--superpopulation` explicitly (AFR/AMR/EAS/EUR/SAS) or `AUTO` if "
                "this ancestry assumption is wrong.",
                "",
            ]
        )
    lines.extend(
        [
            "## Summary",
            "",
            str(data.get("upstream_summary") or "PRS computation complete."),
            "",
            f"- Scored models: **{data.get('n_scored', 0)}**",
            f"- Failed models: **{data.get('n_failed', 0)}**",
            f"- Curated out: **{data.get('n_filtered', 0)}**",
            f"- Build mismatch: **{data.get('build_mismatch', False)}**",
            f"- Requested ancestry: **{requested_label}**",
            f"- Reference panel ancestry: **{panel_ancestry}**",
            "",
            "## Model agreement",
            "",
            f"- Reliable models: **{agreement['reliable_models']}**",
            f"- Verdict: **{agreement['verdict']}**",
            f"- Reliable percentile range: **{_format_number(agreement['percentile_min'])}–"
            f"{_format_number(agreement['percentile_max'])}**",
            f"- Descriptive percentile spread: **{_format_number(agreement['percentile_spread'])}**",
            "",
            "This is a descriptive comparison of reliable models, not a clinical "
            "confidence threshold.",
            "",
            "## Score details",
            "",
            "| PGS ID | Status | Percentile | Reference ancestry | Reliable | C_wt | "
            "Match rate | Quality | Absolute risk |",
            "|---|---|---:|---|---|---:|---:|---|---|",
        ]
    )
    for score in data["scores"]:
        risk = score["absolute_risk"]
        risk_text = (
            _format_fraction(risk.get("absolute_risk"))
            if risk.get("status") == "available"
            else risk.get("status", "unavailable")
        )
        lines.append(
            "| [{pgs}]({url}) | {status} | {pct} | {ancestry} | {reliable} | {coverage} | "
            "{match_rate} | {quality} | {risk} |".format(
                pgs=score["pgs_id"],
                url=score["pgs_catalog_url"],
                status=score.get("status") or "unknown",
                pct=_format_number(score.get("percentile")),
                ancestry=score.get("reference_panel_ancestry") or "n/a",
                reliable=score.get("percentile_reliable"),
                coverage=_format_fraction(score.get("weight_mass_coverage")),
                match_rate=_format_fraction(score.get("match_rate")),
                quality=score.get("quality_label") or "n/a",
                risk=risk_text,
            )
        )
    lines.extend(
        [
            "",
            "## Interpretation boundaries",
            "",
            "- A PRS measures genetic predisposition, not the trait itself and not a diagnosis.",
            "- Percentiles are meaningful only when the reference-panel ancestry is appropriate.",
            "- Absolute-risk percentages may use a single-population prevalence baseline "
            "upstream, so they can be ancestry-inconsistent with an ancestry-stratified "
            "percentile even when both numbers are shown.",
            "- Low C_wt, unreliable percentiles, model disagreement, and build mismatches must not be hidden.",
            "- Lifestyle, environment, age, family history, and genetic factors outside the model also matter.",
            "",
            "## Curation provenance",
            "",
            str(data.get("filter_summary") or "No model-filter summary was returned."),
            "",
            f"Upstream engine: `just-prs-mcp=={JUST_PRS_MCP_VERSION}` via local stdio.",
            "",
            generate_report_footer().strip(),
            "",
        ]
    )
    return "\n".join(lines)


def _write_scores_csv(output_dir: Path, scores: list[dict[str, Any]]) -> Path:
    table_dir = output_dir / "tables"
    table_dir.mkdir(parents=True, exist_ok=True)
    path = table_dir / "scores.csv"
    fieldnames = [
        "pgs_id",
        "status",
        "raw_score",
        "variants_matched",
        "variants_total",
        "match_rate",
        "weight_mass_coverage",
        "percentile",
        "percentile_method",
        "percentile_reliable",
        "reference_panel_ancestry",
        "quality_label",
        "auroc_estimate",
        "absolute_risk_status",
        "absolute_risk",
        "risk_ratio",
        "error",
    ]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        for score in scores:
            risk = score["absolute_risk"]
            writer.writerow(
                {
                    **{key: score.get(key) for key in fieldnames if not key.startswith("absolute_")},
                    "absolute_risk_status": risk.get("status"),
                    "absolute_risk": risk.get("absolute_risk"),
                    "risk_ratio": risk.get("risk_ratio"),
                }
            )
    return path


def _write_reproducibility(
    output_dir: Path,
    *,
    input_path: Path,
    demo: bool,
    trait: str | None = None,
    trait_id: str | None = None,
    pgs_id: str | None = None,
    superpopulation: str = "EUR",
    profile: str = "curated",
    top_n: int = 5,
    limit: int | None = None,
) -> Path:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    args: list[str] = ["--demo"] if demo else ["--input", f'"{input_path}"']
    if not demo:
        if trait:
            args.extend(["--trait", json.dumps(trait)])
        if trait_id:
            args.extend(["--trait-id", trait_id])
        if pgs_id:
            args.extend(["--pgs-id", pgs_id])
        args.extend(
            [
                "--superpopulation",
                superpopulation,
                "--profile",
                profile,
                "--top-n",
                str(top_n),
            ]
        )
        if limit is not None:
            args.extend(["--limit", str(limit)])
    args.extend(["--output", '"$OUTPUT_DIR"'])
    command = " \\\n  ".join(['"${RUNNER[@]}"', *args])
    content = "\n".join(
        [
            "#!/usr/bin/env bash",
            "set -euo pipefail",
            'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
            'OUTPUT_DIR="$(dirname "$SCRIPT_DIR")"',
            f': "${{CLAWBIO_ROOT:={_PROJECT_ROOT}}}"',
            f': "${{CLAWBIO_PACKAGE_SPEC:=clawbio[just-prs]=={distribution_version("clawbio")}}}"',
            'CLAWBIO_SCRIPT="$CLAWBIO_ROOT/skills/just-prs-mcp/just_prs_mcp_bridge.py"',
            'if [[ -f "$CLAWBIO_ROOT/pyproject.toml" ]]; then',
            '  RUNNER=(uv run --project "$CLAWBIO_ROOT" --extra just-prs python "$CLAWBIO_SCRIPT")',
            "else",
            '  RUNNER=(uv run --no-project --with "$CLAWBIO_PACKAGE_SPEC" clawbio run just-prs)',
            "fi",
            command,
            "",
        ]
    )
    path = repro_dir / "commands.sh"
    write_text_lf(path, content)
    path.chmod(path.stat().st_mode | 0o111)
    return path


def write_bundle(
    output_dir: Path,
    data: dict[str, Any],
    *,
    input_path: Path,
    demo: bool,
    trait: str | None = None,
    trait_id: str | None = None,
    pgs_id: str | None = None,
    superpopulation: str = "EUR",
    superpopulation_defaulted: bool = False,
    profile: str = "curated",
    top_n: int = 5,
    limit: int | None = None,
    warn_before_overwrite: bool = True,
) -> dict[str, Any]:
    """Write the stable ClawBio report/result/table/reproducibility contract."""
    if warn_before_overwrite:
        _warn_before_overwrite(output_dir)
    output_dir.mkdir(parents=True, exist_ok=True)
    payload = dict(data)
    payload["requested_superpopulation"] = superpopulation
    payload["superpopulation_defaulted"] = superpopulation_defaulted
    if "reference_panel_ancestry" not in payload:
        panel_ancestries = sorted(
            {
                str(score["reference_panel_ancestry"])
                for score in payload.get("scores", [])
                if score.get("reference_panel_ancestry")
            }
        )
        payload["reference_panel_ancestry"] = (
            panel_ancestries[0]
            if len(panel_ancestries) == 1
            else (", ".join(panel_ancestries) if panel_ancestries else None)
        )
    report_path = output_dir / "report.md"
    write_text_lf(
        report_path,
        _report_markdown(payload, input_path=input_path, demo=demo),
    )
    scores_path = _write_scores_csv(output_dir, payload["scores"])
    input_checksum = sha256_file(input_path)
    ok = int(payload.get("n_scored", 0)) > 0
    result_path = write_result_json(
        output_dir=output_dir,
        skill=SKILL_NAME,
        version=SKILL_VERSION,
        summary={
            "trait": payload.get("trait"),
            "trait_id": payload.get("trait_id"),
            "models_scored": payload.get("n_scored", 0),
            "models_failed": payload.get("n_failed", 0),
            "reliable_models": payload["model_agreement"]["reliable_models"],
            "requested_superpopulation": superpopulation,
            "superpopulation_defaulted": superpopulation_defaulted,
            "reference_panel_ancestry": payload.get("reference_panel_ancestry"),
            "demo": demo,
        },
        data=payload,
        input_checksum=input_checksum,
        datasets={
            "just-prs-mcp": JUST_PRS_MCP_VERSION,
            "PGS Catalog": "live public metadata/scoring files" if not demo else "cached synthetic fixture",
        },
        status="ok" if ok else "failed",
        ok=ok,
    )
    commands_path = _write_reproducibility(
        output_dir,
        input_path=input_path,
        demo=demo,
        trait=trait,
        trait_id=trait_id,
        pgs_id=pgs_id,
        superpopulation=superpopulation,
        profile=profile,
        top_n=top_n,
        limit=limit,
    )
    from clawbio.common.reproducibility import write_checksums

    write_checksums(
        [report_path, result_path, scores_path, commands_path],
        output_dir,
        anchor=output_dir,
    )
    return json.loads(result_path.read_text(encoding="utf-8"))


def run_demo(output_dir: Path) -> dict[str, Any]:
    """Run a deterministic offline demo using an upstream-shaped cached response."""
    report = json.loads(_DEMO_REPORT.read_text(encoding="utf-8"))
    risks = {
        str(row["pgs_id"]): {
            "status": "unavailable",
            "reason": "Offline demo does not call prevalence services.",
        }
        for row in report["rows"]
        if row.get("status") == "scored"
    }
    data = map_trait_report(
        report,
        superpopulation="EUR",
        absolute_risks=risks,
        superpopulation_defaulted=False,
    )
    return write_bundle(
        output_dir,
        data,
        input_path=_DEMO_VCF,
        demo=True,
        superpopulation="EUR",
        superpopulation_defaulted=False,
    )


async def _run_live(
    *,
    input_path: Path,
    output_dir: Path,
    trait: str | None,
    trait_id: str | None,
    pgs_id: str | None,
    superpopulation: str,
    profile: str,
    top_n: int,
    limit: int | None,
    timeout_seconds: int,
    superpopulation_defaulted: bool = False,
) -> dict[str, Any]:
    _warn_before_overwrite(output_dir)
    log_path = output_dir / "reproducibility" / "mcp.stderr.log"
    log_path.parent.mkdir(parents=True, exist_ok=True)
    async with JustPrsMcpClient(
        server_command(),
        server_environment(),
        timeout_seconds=timeout_seconds,
        log_file=log_path,
    ) as client:
        data = await analyze_with_client(
            client=client,
            input_path=input_path,
            trait=trait,
            trait_id=trait_id,
            pgs_id=pgs_id,
            superpopulation=superpopulation,
            profile=profile,
            top_n=top_n,
            limit=limit,
            superpopulation_defaulted=superpopulation_defaulted,
        )
    return write_bundle(
        output_dir,
        data,
        input_path=input_path,
        demo=False,
        trait=trait,
        trait_id=trait_id,
        pgs_id=pgs_id,
        superpopulation=superpopulation,
        superpopulation_defaulted=superpopulation_defaulted,
        profile=profile,
        top_n=top_n,
        limit=limit,
        warn_before_overwrite=False,
    )


def main(
    input_path: Path | None = typer.Option(None, "--input", help="Local VCF/VCF.GZ path."),
    output: Path = typer.Option(Path("output/just-prs"), "--output", help="Report directory."),
    demo: bool = typer.Option(False, "--demo", help="Run the deterministic offline demo."),
    trait: str | None = typer.Option(None, "--trait", help="Trait search term."),
    trait_id: str | None = typer.Option(None, "--trait-id", help="Exact EFO or MONDO ID."),
    pgs_id: str | None = typer.Option(None, "--pgs-id", help="Exact PGS Catalog score ID."),
    superpopulation: str | None = typer.Option(
        None,
        "--superpopulation",
        help=(
            "1000 Genomes ancestry: AFR/AMR/EAS/EUR/SAS, or AUTO to infer. "
            "Defaults to EUR with a visible warning when omitted."
        ),
    ),
    profile: str = typer.Option("curated", "--profile", help="curated or all."),
    top_n: int = typer.Option(5, "--top-n", min=1, help="Maximum returned trait models."),
    limit: int | None = typer.Option(
        None, "--limit", min=1, help="Maximum PGS models computed before curation."
    ),
    timeout_seconds: int = typer.Option(
        3600, "--timeout-seconds", min=30, help="Local MCP call timeout."
    ),
) -> None:
    """Compute evidence-aware PRS reports through a pinned local MCP server."""
    output = output.expanduser().resolve()
    if demo:
        result = run_demo(output)
        typer.echo(json.dumps(result["summary"], indent=2))
        return

    if input_path is None:
        raise typer.BadParameter("--input is required unless --demo is used.")
    input_path = input_path.expanduser().resolve()
    if not input_path.is_file():
        raise typer.BadParameter(f"Input VCF does not exist: {input_path}")
    if not input_path.name.lower().endswith((".vcf", ".vcf.gz", ".vcf.bgz")):
        raise typer.BadParameter("--input must be a VCF, VCF.GZ, or VCF.BGZ file.")
    selectors = [bool(trait), bool(trait_id), bool(pgs_id)]
    if sum(selectors) != 1:
        raise typer.BadParameter("Provide exactly one of --trait, --trait-id, or --pgs-id.")
    try:
        resolved_superpopulation, superpopulation_defaulted = normalize_superpopulation(
            superpopulation
        )
    except ValueError as exc:
        raise typer.BadParameter(str(exc)) from exc
    if superpopulation_defaulted:
        emit_superpopulation_default_warning(resolved_superpopulation)
    if profile not in {"curated", "all"}:
        raise typer.BadParameter("--profile must be curated or all.")

    try:
        result = asyncio.run(
            _run_live(
                input_path=input_path,
                output_dir=output,
                trait=trait,
                trait_id=trait_id,
                pgs_id=pgs_id,
                superpopulation=resolved_superpopulation,
                profile=profile,
                top_n=top_n,
                limit=limit,
                timeout_seconds=timeout_seconds,
                superpopulation_defaulted=superpopulation_defaulted,
            )
        )
    except (AmbiguousTraitError, McpCallError, ValueError) as exc:
        typer.echo(f"ERROR: {exc}", err=True)
        raise typer.Exit(code=1) from exc

    typer.echo(json.dumps(result["summary"], indent=2))
    if not result.get("ok", False):
        raise typer.Exit(code=1)


if __name__ == "__main__":
    typer.run(main)
