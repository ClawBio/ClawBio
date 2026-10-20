#!/usr/bin/env python3
"""
bigquery_public.py — BigQuery Public bridge for ClawBio
=======================================================
Run read-only SQL queries and lightweight discovery commands against BigQuery
public datasets with a local-first workflow. Demo mode uses a bundled offline
fixture so tests and first-run experience do not require cloud authentication.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.report import DISCLAIMER, write_result_json  # noqa: E402

from bigquery_backends import (  # noqa: E402
    auth_setup_message as _auth_setup_message,
    execute_discovery_with_bq_cli,
    execute_discovery_with_python_client,
    execute_with_bq_cli,
    execute_with_python_client,
    get_gcloud_project as _get_gcloud_project,
)
from bigquery_models import (  # noqa: E402
    BigQuerySetupError,
    DISCOVERY_DESCRIBE,
    DISCOVERY_LIST_DATASETS,
    DISCOVERY_LIST_TABLES,
    DiscoveryRequest,
    QueryValidationError,
    QueryExecutionResult,
    SKILL_NAME,
    SKILL_VERSION,
)
from bigquery_support import (  # noqa: E402
    build_parser,
    build_report,
    ensure_output_dir_ready as _ensure_output_dir_ready,
    load_demo_result,
    parse_scalar_param,
    parse_scalar_params,
    resolve_run_plan as _resolve_run_plan,
    validate_read_only_sql,
    write_reproducibility_bundle,
    _write_results_csv,
)


def execute_query(
    query: str,
    location: str,
    max_rows: int,
    max_bytes_billed: int | None,
    parameters: list,
    dry_run: bool,
) -> QueryExecutionResult:
    failures: list[str] = []
    try:
        return execute_with_python_client(
            query=query,
            location=location,
            max_rows=max_rows,
            max_bytes_billed=max_bytes_billed,
            parameters=parameters,
            dry_run=dry_run,
        )
    except BigQuerySetupError as exc:
        failures.append(f"Python ADC: {exc}")

    try:
        return execute_with_bq_cli(
            query=query,
            location=location,
            max_rows=max_rows,
            max_bytes_billed=max_bytes_billed,
            parameters=parameters,
            dry_run=dry_run,
        )
    except BigQuerySetupError as exc:
        failures.append(f"bq CLI: {exc}")

    raise BigQuerySetupError(_auth_setup_message(failures, _get_gcloud_project()))


def execute_discovery(
    request,
    *,
    max_rows: int,
    location: str,
) -> QueryExecutionResult:
    failures: list[str] = []
    try:
        return execute_discovery_with_python_client(request, max_rows=max_rows, location=location)
    except BigQuerySetupError as exc:
        failures.append(f"Python ADC: {exc}")

    try:
        return execute_discovery_with_bq_cli(request, max_rows=max_rows, location=location)
    except BigQuerySetupError as exc:
        failures.append(f"bq CLI: {exc}")

    raise BigQuerySetupError(_auth_setup_message(failures, _get_gcloud_project()))


def run_plan(
    plan,
    output_dir: Path,
    *,
    parameters: list,
    location: str,
    max_rows: int,
    max_bytes_billed: int | None,
    dry_run: bool,
    cli_args: list[str],
) -> QueryExecutionResult:
    _ensure_output_dir_ready(output_dir)

    if plan.mode == "demo":
        result = load_demo_result(
            query=plan.effective_query or "",
            location=location,
            max_rows=max_rows,
            dry_run=dry_run,
        )
    elif plan.mode == "discovery":
        assert plan.discovery is not None
        result = execute_discovery(plan.discovery, max_rows=max_rows, location=location)
    else:
        assert plan.effective_query is not None
        result = execute_query(
            query=plan.effective_query,
            location=location,
            max_rows=max_rows,
            max_bytes_billed=max_bytes_billed,
            parameters=parameters,
            dry_run=dry_run,
        )

    report_text = build_report(
        plan=plan,
        result=result,
        parameters=parameters,
        max_rows=max_rows,
        max_bytes_billed=max_bytes_billed,
    )
    (output_dir / "report.md").write_text(report_text, encoding="utf-8")
    _write_results_csv(output_dir / "tables" / "results.csv", result.rows, result.columns)

    summary = {
        "mode": plan.mode,
        "dry_run": result.dry_run,
        "backend": result.backend,
        "project_id": result.project_id,
        "location": result.location,
        "row_count": result.row_count,
        "max_rows": max_rows,
        "estimated_bytes_processed": result.estimated_bytes_processed,
        "total_bytes_processed": result.total_bytes_processed,
        "query_source": plan.query_source,
    }
    if plan.discovery:
        summary["discovery_action"] = plan.discovery.action
        summary["discovery_target"] = plan.discovery.target

    data = {
        "query": plan.effective_query or "",
        "source_query": plan.source_query,
        "effective_query": plan.effective_query,
        "columns": result.columns,
        "rows": result.rows,
        "parameters": [
            {"name": param.name, "type": param.type_name, "value": param.original}
            for param in parameters
        ],
        "paper_reference": plan.paper_reference,
        "notes": plan.notes,
        "warnings": plan.warnings,
        "job_metadata": result.raw_metadata,
        "disclaimer": DISCLAIMER,
    }
    if plan.discovery:
        data["discovery"] = {
            "action": plan.discovery.action,
            "target": plan.discovery.target,
        }
    write_result_json(
        output_dir=output_dir,
        skill=SKILL_NAME,
        version=SKILL_VERSION,
        summary=summary,
        data=data,
    )
    write_reproducibility_bundle(output_dir, plan, result, cli_args)
    return result


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    if args.max_rows <= 0:
        parser.error("--max-rows must be greater than 0")
    if args.max_bytes_billed is not None and args.max_bytes_billed <= 0:
        parser.error("--max-bytes-billed must be greater than 0")

    try:
        plan = _resolve_run_plan(args)
        params = parse_scalar_params(args.param) if plan.mode == "query" else []
        output_dir = Path(args.output)
        result = run_plan(
            plan=plan,
            output_dir=output_dir,
            parameters=params,
            location=args.location,
            max_rows=args.max_rows,
            max_bytes_billed=args.max_bytes_billed,
            dry_run=args.dry_run,
            cli_args=sys.argv[1:] if argv is None else argv,
        )
    except (BigQuerySetupError, QueryValidationError, FileNotFoundError, ValueError) as exc:
        print(f"ERROR: {exc}", file=sys.stderr)
        return 1

    for warning in plan.warnings:
        print(f"WARNING: {warning}")
    print(f"Report written to {output_dir / 'report.md'}")
    print(f"Rows returned: {result.row_count}")
    print(f"Backend: {result.backend}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
