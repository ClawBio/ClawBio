from __future__ import annotations

import argparse
import csv
import json
import re
import shlex
import sys
from decimal import Decimal, InvalidOperation
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.report import (  # noqa: E402
    DISCLAIMER,
    generate_report_footer,
    generate_report_header,
)
from clawbio.common.reproducibility import (  # noqa: E402
    write_checksums,
    write_commands_sh,
    write_environment_yml,
)

from bigquery_models import (
    DEFAULT_MAX_BYTES_BILLED,
    DEFAULT_MAX_ROWS,
    DEFAULT_LOCATION,
    DEMO_QUERY_PATH,
    DEMO_RESULT_PATH,
    DISCOVERY_DESCRIBE,
    DISCOVERY_LIST_DATASETS,
    DISCOVERY_LIST_TABLES,
    GOOGLE_AUTH_REQUIREMENT,
    GOOGLE_BIGQUERY_REQUIREMENT,
    QueryExecutionResult,
    QueryParameter,
    QueryValidationError,
    RunPlan,
    DiscoveryRequest,
    SKILL_NAME,
)


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="ClawBio BigQuery Public — SQL-first bridge for public datasets",
    )
    parser.add_argument("--input", help="Path to a SQL file")
    parser.add_argument("--query", help="Inline SQL query string")
    parser.add_argument("--output", required=True, help="Directory to write outputs")
    parser.add_argument("--demo", action="store_true", help="Run offline demo using bundled fixture data")
    parser.add_argument("--dry-run", action="store_true", help="Estimate bytes only; do not execute the query")
    parser.add_argument("--location", default=DEFAULT_LOCATION, help=f"BigQuery location (default: {DEFAULT_LOCATION})")
    parser.add_argument("--max-rows", type=int, default=DEFAULT_MAX_ROWS, help=f"Maximum rows to return (default: {DEFAULT_MAX_ROWS})")
    parser.add_argument(
        "--max-bytes-billed",
        type=int,
        default=DEFAULT_MAX_BYTES_BILLED,
        help=f"Maximum billed bytes safeguard (default: {DEFAULT_MAX_BYTES_BILLED})",
    )
    parser.add_argument(
        "--param",
        action="append",
        default=[],
        help="Scalar query parameter in name=type:value format (repeatable)",
    )
    parser.add_argument("--list-datasets", help="List datasets for a project (project)")
    parser.add_argument("--list-tables", help="List tables for a dataset (project.dataset)")
    parser.add_argument("--describe", help="Describe top-level schema for a table (project.dataset.table)")
    parser.add_argument("--preview", type=int, help="Wrap the SQL query in a preview LIMIT")
    parser.add_argument("--count-only", action="store_true", help="Return only the row count for the SQL query")
    parser.add_argument("--paper", help="Paper reference, DOI, URL, title, or local PDF path")
    parser.add_argument("--note", action="append", default=[], help="Repeatable provenance note")
    return parser


def _mask_sql_literals(text: str) -> str:
    patterns = [
        r"'(?:''|[^'])*'",
        r'"(?:\\"|[^"])*"',
        r"`(?:``|[^`])*`",
    ]
    masked = text
    for pattern in patterns:
        masked = re.sub(pattern, lambda m: " " * len(m.group(0)), masked, flags=re.DOTALL)
    return masked


def _strip_sql_comments(text: str) -> str:
    no_block = re.sub(r"/\*.*?\*/", " ", text, flags=re.DOTALL)
    return re.sub(r"--.*?$", " ", no_block, flags=re.MULTILINE)


def _analysis_sql(query: str) -> str:
    return re.sub(r"\s+", " ", _mask_sql_literals(_strip_sql_comments(query))).strip().upper()


def validate_read_only_sql(query: str) -> str:
    cleaned = query.strip()
    if not cleaned:
        raise QueryValidationError("Query is empty.")

    masked = _mask_sql_literals(_strip_sql_comments(cleaned))
    masked_stripped = masked.strip()
    leading = masked_stripped.upper()
    if not (leading.startswith("SELECT") or leading.startswith("WITH")):
        raise QueryValidationError("Only read-only SELECT/WITH queries are supported.")

    if ";" in masked_stripped[:-1]:
        raise QueryValidationError("Multiple SQL statements are not supported.")
    if masked_stripped.endswith(";"):
        cleaned = cleaned.rstrip()
        cleaned = cleaned[:-1].rstrip()

    analysis = _analysis_sql(cleaned)
    if re.search(r"\bEXECUTE\s+IMMEDIATE\b", analysis, flags=re.IGNORECASE):
        raise QueryValidationError("Unsupported SQL keyword detected: EXECUTE IMMEDIATE")
    if re.search(r"\bBEGIN\b", analysis, flags=re.IGNORECASE):
        raise QueryValidationError("Unsupported SQL keyword detected: BEGIN")
    if re.search(r"\bASSERT\b", analysis, flags=re.IGNORECASE):
        raise QueryValidationError("Unsupported SQL keyword detected: ASSERT")
    if re.search(r"\bEXCEPTION\b", analysis, flags=re.IGNORECASE):
        raise QueryValidationError("Unsupported SQL keyword detected: EXCEPTION")
    if re.search(r"\bWHILE\b", analysis, flags=re.IGNORECASE):
        raise QueryValidationError("Unsupported SQL keyword detected: WHILE")
    # END is only rejected when paired with BEGIN so CASE ... END expressions remain valid.
    if re.search(r"\bEND\b", analysis, flags=re.IGNORECASE) and re.search(r"\bBEGIN\b", analysis, flags=re.IGNORECASE):
        raise QueryValidationError("Unsupported SQL keyword detected: END")

    forbidden = re.compile(
        r"\b(INSERT|UPDATE|DELETE|CREATE|MERGE|EXPORT\s+DATA|DROP|ALTER|TRUNCATE|CALL|DECLARE|SET)\b",
        flags=re.IGNORECASE,
    )
    match = forbidden.search(masked)
    if match:
        raise QueryValidationError(f"Unsupported SQL keyword detected: {match.group(0).strip()}")

    return cleaned


def parse_scalar_param(spec: str) -> QueryParameter:
    if "=" not in spec or ":" not in spec.split("=", 1)[1]:
        raise ValueError(f"Invalid --param value: {spec!r}. Expected name=type:value")

    name, typed_value = spec.split("=", 1)
    type_name, raw_value = typed_value.split(":", 1)
    name = name.strip()
    type_name = type_name.strip().upper()
    raw_value = raw_value.strip()

    if not name:
        raise ValueError(f"Invalid --param value: {spec!r}. Parameter name is empty.")
    if any(char in raw_value for char in ("\n", "\r", "\x00")):
        raise ValueError(f"Invalid --param value: {spec!r}. Parameter values must not contain newlines or NUL bytes.")

    if type_name in {"STRING", "DATE", "DATETIME", "TIMESTAMP"}:
        value: Any = raw_value
    elif type_name in {"INT64", "INTEGER"}:
        value = int(raw_value)
        type_name = "INT64"
    elif type_name in {"FLOAT64", "FLOAT"}:
        value = float(raw_value)
        type_name = "FLOAT64"
    elif type_name == "NUMERIC":
        try:
            value = Decimal(raw_value)
        except InvalidOperation as exc:
            raise ValueError(f"Invalid NUMERIC parameter value: {raw_value!r}") from exc
    elif type_name in {"BOOL", "BOOLEAN"}:
        lowered = raw_value.lower()
        if lowered not in {"true", "false"}:
            raise ValueError(f"Invalid boolean parameter value: {raw_value!r}")
        value = lowered == "true"
        type_name = "BOOL"
    else:
        raise ValueError(
            f"Unsupported parameter type {type_name!r}. "
            "Supported types: STRING, INT64, FLOAT64, BOOL, DATE, DATETIME, TIMESTAMP, NUMERIC."
        )

    return QueryParameter(name=name, type_name=type_name, value=value, original=raw_value)


def parse_scalar_params(specs: list[str]) -> list[QueryParameter]:
    return [parse_scalar_param(spec) for spec in specs]


def parse_project_dataset(spec: str) -> tuple[str, str]:
    parts = spec.strip().split(".")
    if len(parts) != 2 or not all(parts):
        raise ValueError(f"Expected project.dataset, got: {spec!r}")
    return parts[0], parts[1]


def parse_project_dataset_table(spec: str) -> tuple[str, str, str]:
    parts = spec.strip().split(".")
    if len(parts) != 3 or not all(parts):
        raise ValueError(f"Expected project.dataset.table, got: {spec!r}")
    return parts[0], parts[1], parts[2]


def _wrap_preview_query(query: str, preview_rows: int) -> str:
    return f"SELECT * FROM (\n{query}\n) AS clawbio_preview LIMIT {preview_rows}"


def _wrap_count_query(query: str) -> str:
    return f"SELECT COUNT(*) AS row_count FROM (\n{query}\n) AS clawbio_count"


def _has_select_star(query: str) -> bool:
    return bool(re.search(r"\bSELECT\s+(?:ALL\s+|DISTINCT\s+)?\*", _analysis_sql(query), flags=re.IGNORECASE))


def _has_limit_clause(query: str) -> bool:
    return bool(re.search(r"\bLIMIT\b", _analysis_sql(query), flags=re.IGNORECASE))


def collect_query_warnings(
    source_query: str,
    *,
    dry_run: bool,
    preview_rows: int | None,
    count_only: bool,
) -> list[str]:
    warnings: list[str] = []
    if _has_select_star(source_query):
        warnings.append("Query uses SELECT *; consider selecting only the columns you need.")
    if not dry_run and preview_rows is None and not count_only and not _has_limit_clause(source_query):
        warnings.append("Query does not include LIMIT; be sure the result size is intentional.")
    return warnings


def _read_query_from_args(args: argparse.Namespace) -> tuple[str, str]:
    if args.demo:
        return DEMO_QUERY_PATH.read_text(encoding="utf-8"), "demo-query"

    if args.query:
        if args.input:
            print("WARNING: --query provided; ignoring --input SQL file.", file=sys.stderr)
        return args.query, "inline-query"

    if args.input:
        input_path = Path(args.input)
        if not input_path.exists():
            raise FileNotFoundError(f"SQL file not found: {input_path}")
        return input_path.read_text(encoding="utf-8"), str(input_path)

    raise ValueError("Provide --query, --input <sql_file>, --demo, or a discovery option.")


def resolve_run_plan(args: argparse.Namespace) -> RunPlan:
    discovery_options = [
        (DISCOVERY_LIST_DATASETS, args.list_datasets),
        (DISCOVERY_LIST_TABLES, args.list_tables),
        (DISCOVERY_DESCRIBE, args.describe),
    ]
    selected_discovery = [(action, target) for action, target in discovery_options if target]

    if len(selected_discovery) > 1:
        raise ValueError("Choose only one of --list-datasets, --list-tables, or --describe.")
    if args.preview is not None and args.preview <= 0:
        raise ValueError("--preview must be greater than 0.")
    if args.preview is not None and args.count_only:
        raise ValueError("--preview and --count-only cannot be used together.")

    if selected_discovery:
        if args.demo or args.query or args.input:
            raise ValueError("Discovery options are mutually exclusive with --query, --input, and --demo.")
        if args.dry_run:
            raise ValueError("--dry-run is only supported with --query or --input.")
        if args.preview is not None or args.count_only:
            raise ValueError("--preview and --count-only are only supported with --query or --input.")
        action, target = selected_discovery[0]
        if action == DISCOVERY_LIST_TABLES:
            parse_project_dataset(target)
        elif action == DISCOVERY_DESCRIBE:
            parse_project_dataset_table(target)
        return RunPlan(
            mode="discovery",
            query_source=f"discovery:{action}",
            source_query=None,
            effective_query=None,
            discovery=DiscoveryRequest(action=action, target=target),
            warnings=[],
            paper_reference=args.paper,
            notes=args.note or [],
        )

    query_text, query_source = _read_query_from_args(args)
    source_query = validate_read_only_sql(query_text)
    effective_query = source_query

    if args.demo and (args.preview is not None or args.count_only):
        raise ValueError("--preview and --count-only are not supported with --demo.")
    if args.count_only:
        effective_query = validate_read_only_sql(_wrap_count_query(source_query))
    elif args.preview is not None:
        effective_query = validate_read_only_sql(_wrap_preview_query(source_query, args.preview))

    warnings = [] if args.demo else collect_query_warnings(
        source_query,
        dry_run=args.dry_run,
        preview_rows=args.preview,
        count_only=args.count_only,
    )

    return RunPlan(
        mode="demo" if args.demo else "query",
        query_source=query_source,
        source_query=source_query,
        effective_query=effective_query,
        discovery=None,
        warnings=warnings,
        paper_reference=args.paper,
        notes=args.note or [],
    )


def _json_safe_value(value: Any) -> Any:
    if value is None or isinstance(value, (str, int, float, bool)):
        return value
    if isinstance(value, list):
        return [_json_safe_value(item) for item in value]
    if isinstance(value, dict):
        return {str(key): _json_safe_value(val) for key, val in value.items()}
    return str(value)


def _infer_columns(rows: list[dict[str, Any]]) -> list[str]:
    columns: list[str] = []
    seen: set[str] = set()
    for row in rows:
        for key in row.keys():
            if key not in seen:
                seen.add(key)
                columns.append(key)
    return columns


def _write_results_csv(path: Path, rows: list[dict[str, Any]], columns: list[str]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        if not columns:
            return
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        for row in rows:
            writer.writerow(
                {
                    col: json.dumps(_json_safe_value(row.get(col)), ensure_ascii=False)
                    if isinstance(row.get(col), (list, dict))
                    else _json_safe_value(row.get(col))
                    for col in columns
                }
            )


def _build_provenance_payload(plan: RunPlan, result: QueryExecutionResult) -> dict[str, Any]:
    payload = {
        "paper_reference": plan.paper_reference,
        "notes": plan.notes,
        "query_source": plan.query_source,
        "source_query": plan.source_query,
        "effective_query": plan.effective_query,
        "backend": result.backend,
        "project_id": result.project_id,
        "location": result.location,
        "timestamp_utc": datetime.now(timezone.utc).isoformat(),
    }
    if plan.discovery:
        payload["discovery_action"] = plan.discovery.action
        payload["discovery_target"] = plan.discovery.target
    if plan.warnings:
        payload["warnings"] = plan.warnings
    return payload


def _default_query_sql(plan: RunPlan) -> str:
    if plan.effective_query:
        return plan.effective_query.rstrip() + "\n"
    if plan.discovery:
        return (
            "-- No SQL query executed.\n"
            f"-- Discovery action: {plan.discovery.action}\n"
            f"-- Discovery target: {plan.discovery.target}\n"
        )
    return "-- No SQL query executed.\n"


def write_reproducibility_bundle(
    output_dir: Path,
    plan: RunPlan,
    result: QueryExecutionResult,
    cli_args: list[str],
) -> None:
    """Write reproducibility/ last, so checksums cover every file already in output_dir."""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    write_commands_sh(output_dir, shlex.join(["python", f"skills/{SKILL_NAME}/bigquery_public.py", *cli_args]))
    (repro_dir / "query.sql").write_text(_default_query_sql(plan), encoding="utf-8")
    (repro_dir / "job_metadata.json").write_text(
        json.dumps(result.raw_metadata, indent=2, default=str),
        encoding="utf-8",
    )
    (repro_dir / "provenance.json").write_text(
        json.dumps(_build_provenance_payload(plan, result), indent=2, default=str),
        encoding="utf-8",
    )
    write_environment_yml(
        output_dir,
        f"clawbio-{SKILL_NAME}",
        [GOOGLE_AUTH_REQUIREMENT, GOOGLE_BIGQUERY_REQUIREMENT],
    )
    outputs = sorted(p for p in output_dir.rglob("*") if p.is_file() and p.name != "checksums.sha256")
    write_checksums(outputs, output_dir, anchor=output_dir)


def _render_markdown_table(rows: list[dict[str, Any]], columns: list[str], limit: int = 10) -> str:
    if not rows or not columns:
        return "_No rows returned._"

    visible = rows[:limit]
    lines = [
        "| " + " | ".join(columns) + " |",
        "| " + " | ".join(["---"] * len(columns)) + " |",
    ]
    for row in visible:
        values = []
        for col in columns:
            cell = row.get(col, "")
            if isinstance(cell, (dict, list)):
                text = json.dumps(cell, ensure_ascii=False)
            else:
                text = str(cell)
            values.append(text.replace("\n", " ").replace("|", "\\|"))
        lines.append("| " + " | ".join(values) + " |")
    return "\n".join(lines)


def _format_int(value: int | None) -> str:
    return f"{value:,}" if value is not None else "n/a"


def build_report(
    plan: RunPlan,
    result: QueryExecutionResult,
    parameters: list[QueryParameter],
    max_rows: int,
    max_bytes_billed: int | None,
) -> str:
    execution_mode = "Discovery" if plan.mode == "discovery" else ("Dry run" if result.dry_run else "Query")
    metadata = {
        "Execution mode": execution_mode,
        "Mode": plan.mode,
        "Backend": result.backend,
        "Project": result.project_id or "n/a",
        "Location": result.location,
        "Query source": plan.query_source,
        "Rows returned": str(result.row_count),
        "Estimated bytes processed": _format_int(result.estimated_bytes_processed),
        "Actual bytes processed": _format_int(result.total_bytes_processed),
        "Max rows": str(max_rows),
        "Max bytes billed": _format_int(max_bytes_billed),
    }
    if plan.discovery:
        metadata["Discovery action"] = plan.discovery.action
        metadata["Discovery target"] = plan.discovery.target

    lines = [generate_report_header("BigQuery Public Query Report", SKILL_NAME, extra_metadata=metadata)]
    lines.extend(
        [
            "## Summary",
            "",
            f"- Run mode: `{plan.mode}`",
            f"- Execution backend: `{result.backend}`",
            f"- Location: `{result.location}`",
            f"- Rows returned: `{result.row_count}`",
            f"- Estimated bytes processed: `{_format_int(result.estimated_bytes_processed)}`",
            f"- Actual bytes processed: `{_format_int(result.total_bytes_processed)}`",
            "",
        ]
    )

    if plan.warnings:
        lines.append("## Warnings")
        lines.append("")
        for warning in plan.warnings:
            lines.append(f"- {warning}")
        lines.append("")

    if plan.paper_reference or plan.notes:
        lines.append("## Provenance")
        lines.append("")
        if plan.paper_reference:
            lines.append(f"- Paper reference: `{plan.paper_reference}`")
        for note in plan.notes:
            lines.append(f"- Note: {note}")
        lines.append("")

    if plan.discovery:
        lines.extend(
            [
                "## Discovery Request",
                "",
                f"- Action: `{plan.discovery.action}`",
                f"- Target: `{plan.discovery.target}`",
                "",
            ]
        )
    elif plan.effective_query:
        if plan.source_query == plan.effective_query:
            lines.extend(
                [
                    "## Query",
                    "",
                    "```sql",
                    plan.effective_query,
                    "```",
                    "",
                ]
            )
        else:
            lines.extend(
                [
                    "## Source Query",
                    "",
                    "```sql",
                    plan.source_query or "",
                    "```",
                    "",
                    "## Effective Query",
                    "",
                    "```sql",
                    plan.effective_query,
                    "```",
                    "",
                ]
            )

    if parameters:
        lines.append("## Parameters")
        lines.append("")
        for param in parameters:
            lines.append(f"- `{param.name}` ({param.type_name}) = `{param.original}`")
        lines.append("")

    lines.extend(
        [
            "## Results Preview",
            "",
            _render_markdown_table(result.rows, result.columns),
            "",
        ]
    )

    if result.job_id:
        lines.extend(["## Job Metadata", "", f"- Job ID: `{result.job_id}`", ""])

    lines.append(generate_report_footer())
    return "\n".join(lines).strip() + "\n"


def ensure_output_dir_ready(output_dir: Path) -> None:
    if output_dir.exists():
        if any(output_dir.iterdir()):
            raise ValueError(
                f"Output directory already exists and is not empty: {output_dir}. "
                "Choose a new directory to avoid overwriting previous results."
            )
        return
    output_dir.mkdir(parents=True, exist_ok=True)


def load_demo_result(query: str, location: str, max_rows: int, dry_run: bool) -> QueryExecutionResult:
    payload = json.loads(DEMO_RESULT_PATH.read_text(encoding="utf-8"))
    rows = payload["rows"][:max_rows] if not dry_run else []
    columns = payload.get("columns") or _infer_columns(rows)
    bytes_processed = payload.get("total_bytes_processed")
    raw_metadata = {
        "backend": "demo-fixture",
        "project_id": payload.get("project_id"),
        "location": location,
        "job_id": payload.get("job_id"),
        "estimated_bytes_processed": bytes_processed,
        "total_bytes_processed": None if dry_run else bytes_processed,
        "demo_source": str(DEMO_RESULT_PATH),
    }
    return QueryExecutionResult(
        backend="demo-fixture",
        project_id=payload.get("project_id"),
        location=location,
        query=query,
        dry_run=dry_run,
        rows=rows,
        columns=columns,
        estimated_bytes_processed=bytes_processed,
        total_bytes_processed=None if dry_run else bytes_processed,
        row_count=len(rows),
        job_id=payload.get("job_id"),
        raw_metadata=raw_metadata,
    )
