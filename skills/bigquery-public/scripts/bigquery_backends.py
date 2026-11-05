from __future__ import annotations

import json
import shutil
import subprocess
from typing import Any

# Port of ClawBio 1.x `bigquery-public/bigquery_backends.py` for ClawBio 2.0.
# Unchanged: this module has no dependency on the 1.x `clawbio.common`
# package, and google.auth / google.cloud.bigquery are imported lazily inside
# _load_python_bigquery_client so importing this file never requires those
# packages to be installed (needed for demo mode, which never calls into
# this module at all — see bigquery_support.run_plan / load_demo_result).

from bigquery_models import (
    BigQuerySetupError,
    DISCOVERY_DESCRIBE,
    DISCOVERY_LIST_DATASETS,
    DISCOVERY_LIST_TABLES,
    DiscoveryRequest,
    QueryExecutionResult,
    QueryParameter,
)
from bigquery_support import _infer_columns, _json_safe_value, parse_project_dataset, parse_project_dataset_table


def get_gcloud_project() -> str | None:
    if not shutil.which("gcloud"):
        return None
    proc = subprocess.run(
        ["gcloud", "config", "get-value", "project"],
        capture_output=True,
        text=True,
        check=False,
    )
    if proc.returncode != 0:
        return None
    value = proc.stdout.strip()
    if not value or value == "(unset)":
        return None
    return value


def auth_setup_message(reasons: list[str], project_id: str | None) -> str:
    lines = [
        "BigQuery authentication is not available for this run.",
        "",
        "Backends tried:",
    ]
    lines.extend(f"- {reason}" for reason in reasons)
    lines.extend(
        [
            "",
            "Suggested setup:",
            "1. gcloud auth login",
            "2. gcloud auth application-default login",
            f"3. gcloud config set project {project_id or 'YOUR_PROJECT_ID'}",
            "4. Re-run the command, or set GOOGLE_APPLICATION_CREDENTIALS for service-account based access.",
        ]
    )
    return "\n".join(lines)


def _extract_named_value(payload: Any, target_key: str) -> Any:
    if isinstance(payload, dict):
        for key, value in payload.items():
            if key.lower() == target_key.lower():
                return value
            nested = _extract_named_value(value, target_key)
            if nested is not None:
                return nested
    elif isinstance(payload, list):
        for item in payload:
            nested = _extract_named_value(item, target_key)
            if nested is not None:
                return nested
    return None


def _try_parse_json(text: str) -> Any:
    if not text.strip():
        return None
    try:
        return json.loads(text)
    except json.JSONDecodeError:
        return None


def _normalize_bq_cli_rows(payload: Any) -> list[dict[str, Any]]:
    if isinstance(payload, list):
        return [_json_safe_value(row) for row in payload if isinstance(row, dict)]

    if isinstance(payload, dict) and isinstance(payload.get("rows"), list):
        schema_fields = payload.get("schema", {}).get("fields", [])
        field_names = [field.get("name", f"col_{idx}") for idx, field in enumerate(schema_fields)]
        rows: list[dict[str, Any]] = []
        for row in payload["rows"]:
            cells = row.get("f", [])
            row_data = {}
            for idx, cell in enumerate(cells):
                if idx >= len(field_names):
                    continue
                row_data[field_names[idx]] = _json_safe_value(cell.get("v"))
            rows.append(row_data)
        return rows

    return []


def _build_python_query_parameters(parameters: list[QueryParameter], bigquery_module: Any) -> list[Any]:
    return [
        bigquery_module.ScalarQueryParameter(param.name, param.type_name, param.value)
        for param in parameters
    ]


def _load_python_bigquery_client(project_id: str | None = None) -> tuple[Any, Any, str]:
    try:
        import google.auth
        from google.auth.exceptions import DefaultCredentialsError
        from google.cloud import bigquery
    except ImportError as exc:
        raise BigQuerySetupError(f"Python BigQuery client unavailable: {exc}") from exc

    try:
        credentials, default_project = google.auth.default(
            scopes=["https://www.googleapis.com/auth/cloud-platform"]
        )
    except DefaultCredentialsError as exc:
        raise BigQuerySetupError(f"ADC unavailable: {exc}") from exc

    active_project = project_id or default_project or get_gcloud_project()
    if not active_project:
        raise BigQuerySetupError("No Google Cloud project configured for the Python client.")

    client = bigquery.Client(project=active_project, credentials=credentials)
    return client, bigquery, active_project


def _execute_with_python_client_once(
    query: str,
    location: str,
    max_rows: int,
    max_bytes_billed: int | None,
    parameters: list[QueryParameter],
    dry_run: bool,
    project_id: str | None,
) -> QueryExecutionResult:
    client, bigquery, active_project = _load_python_bigquery_client(project_id)

    job_config = bigquery.QueryJobConfig(
        dry_run=dry_run,
        use_legacy_sql=False,
        maximum_bytes_billed=max_bytes_billed,
        query_parameters=_build_python_query_parameters(parameters, bigquery),
    )

    try:
        query_job = client.query(query, location=location, job_config=job_config)
    except Exception as exc:
        raise BigQuerySetupError(f"Python BigQuery query failed: {exc}") from exc

    rows: list[dict[str, Any]] = []
    columns: list[str] = []
    if not dry_run:
        try:
            iterator = query_job.result(max_results=max_rows)
        except Exception as exc:
            raise BigQuerySetupError(f"Python BigQuery result fetch failed: {exc}") from exc
        columns = [field.name for field in getattr(iterator, "schema", [])]
        for row in iterator:
            rows.append({key: _json_safe_value(value) for key, value in dict(row).items()})

    total_bytes = getattr(query_job, "total_bytes_processed", None)
    raw_metadata = {
        "backend": "python-adc",
        "project_id": active_project,
        "location": location,
        "job_id": getattr(query_job, "job_id", None),
        "state": getattr(query_job, "state", None),
        "total_bytes_processed": total_bytes,
        "cache_hit": getattr(query_job, "cache_hit", None),
    }
    return QueryExecutionResult(
        backend="python-adc",
        project_id=active_project,
        location=location,
        query=query,
        dry_run=dry_run,
        rows=rows,
        columns=columns or _infer_columns(rows),
        estimated_bytes_processed=total_bytes if dry_run else None,
        total_bytes_processed=None if dry_run else total_bytes,
        row_count=len(rows),
        job_id=getattr(query_job, "job_id", None),
        raw_metadata=raw_metadata,
    )


def execute_with_python_client(
    query: str,
    location: str,
    max_rows: int,
    max_bytes_billed: int | None,
    parameters: list[QueryParameter],
    dry_run: bool,
    project_id: str | None = None,
) -> QueryExecutionResult:
    if dry_run:
        return _execute_with_python_client_once(
            query=query,
            location=location,
            max_rows=max_rows,
            max_bytes_billed=max_bytes_billed,
            parameters=parameters,
            dry_run=True,
            project_id=project_id,
        )

    estimate = _execute_with_python_client_once(
        query=query,
        location=location,
        max_rows=max_rows,
        max_bytes_billed=max_bytes_billed,
        parameters=parameters,
        dry_run=True,
        project_id=project_id,
    )
    actual = _execute_with_python_client_once(
        query=query,
        location=location,
        max_rows=max_rows,
        max_bytes_billed=max_bytes_billed,
        parameters=parameters,
        dry_run=False,
        project_id=estimate.project_id,
    )
    actual.estimated_bytes_processed = estimate.estimated_bytes_processed or estimate.total_bytes_processed
    return actual


def _run_bq_query_command(base_cmd: list[str], query: str) -> tuple[subprocess.CompletedProcess[str], str]:
    # Keep SQL out of argv so long queries do not hit OS argument length limits.
    proc = subprocess.run(base_cmd, input=query, capture_output=True, text=True, check=False)
    return proc, "stdin"


def _execute_with_bq_cli_once(
    query: str,
    location: str,
    max_rows: int,
    max_bytes_billed: int | None,
    parameters: list[QueryParameter],
    dry_run: bool,
    project_id: str | None,
) -> QueryExecutionResult:
    if not shutil.which("bq"):
        raise BigQuerySetupError("bq CLI is not installed.")

    active_project = project_id or get_gcloud_project()
    if not active_project:
        raise BigQuerySetupError("No Google Cloud project configured for the bq CLI.")

    cmd = [
        "bq",
        f"--project_id={active_project}",
        f"--location={location}",
        "query",
        "--use_legacy_sql=false",
        "--format=prettyjson" if dry_run else "--format=json",
        f"--max_rows={max_rows}",
    ]
    if dry_run:
        cmd.append("--dry_run")
    if max_bytes_billed is not None:
        cmd.append(f"--maximum_bytes_billed={max_bytes_billed}")
    for param in parameters:
        cmd.append(f"--parameter={param.to_cli_spec()}")

    proc, query_transport = _run_bq_query_command(cmd, query)
    if proc.returncode != 0:
        detail = proc.stderr.strip() or proc.stdout.strip() or "Unknown bq CLI error."
        raise BigQuerySetupError(f"bq CLI query failed: {detail}")

    parsed = _try_parse_json(proc.stdout)
    rows = [] if dry_run else _normalize_bq_cli_rows(parsed)
    columns = _infer_columns(rows)
    bytes_processed = _extract_named_value(parsed, "totalBytesProcessed")
    try:
        bytes_processed = int(bytes_processed) if bytes_processed is not None else None
    except (TypeError, ValueError):
        bytes_processed = None
    job_id = _extract_named_value(parsed, "jobId")
    raw_metadata = {
        "backend": "bq-cli",
        "project_id": active_project,
        "location": location,
        "job_id": job_id,
        "query_transport": query_transport,
        "total_bytes_processed": bytes_processed,
        "raw_response": parsed if parsed is not None else proc.stdout.strip(),
    }
    return QueryExecutionResult(
        backend="bq-cli",
        project_id=active_project,
        location=location,
        query=query,
        dry_run=dry_run,
        rows=rows,
        columns=columns,
        estimated_bytes_processed=bytes_processed if dry_run else None,
        total_bytes_processed=None if dry_run else bytes_processed,
        row_count=len(rows),
        job_id=str(job_id) if job_id is not None else None,
        raw_metadata=raw_metadata,
    )


def execute_with_bq_cli(
    query: str,
    location: str,
    max_rows: int,
    max_bytes_billed: int | None,
    parameters: list[QueryParameter],
    dry_run: bool,
    project_id: str | None = None,
) -> QueryExecutionResult:
    if dry_run:
        return _execute_with_bq_cli_once(
            query=query,
            location=location,
            max_rows=max_rows,
            max_bytes_billed=max_bytes_billed,
            parameters=parameters,
            dry_run=True,
            project_id=project_id,
        )

    estimate = _execute_with_bq_cli_once(
        query=query,
        location=location,
        max_rows=max_rows,
        max_bytes_billed=max_bytes_billed,
        parameters=parameters,
        dry_run=True,
        project_id=project_id,
    )
    actual = _execute_with_bq_cli_once(
        query=query,
        location=location,
        max_rows=max_rows,
        max_bytes_billed=max_bytes_billed,
        parameters=parameters,
        dry_run=False,
        project_id=estimate.project_id,
    )
    actual.estimated_bytes_processed = estimate.estimated_bytes_processed or estimate.total_bytes_processed
    return actual


def _result_from_discovery_rows(
    *,
    backend: str,
    active_project: str,
    location: str,
    rows: list[dict[str, Any]],
    columns: list[str],
    raw_metadata: dict[str, Any],
) -> QueryExecutionResult:
    return QueryExecutionResult(
        backend=backend,
        project_id=active_project,
        location=location,
        query="",
        dry_run=False,
        rows=rows,
        columns=columns,
        estimated_bytes_processed=None,
        total_bytes_processed=None,
        row_count=len(rows),
        job_id=raw_metadata.get("job_id"),
        raw_metadata=raw_metadata,
    )


def execute_discovery_with_python_client(
    request: DiscoveryRequest,
    *,
    max_rows: int,
    location: str,
) -> QueryExecutionResult:
    client, bigquery, active_project = _load_python_bigquery_client()
    try:
        if request.action == DISCOVERY_LIST_DATASETS:
            target_project = request.target
            rows = []
            for dataset in client.list_datasets(project=target_project, max_results=max_rows):
                rows.append(
                    {
                        "project_id": target_project,
                        "dataset_id": dataset.dataset_id,
                        "location": getattr(dataset, "location", "") or "",
                    }
                )
            result_location = rows[0]["location"] if rows and len({row["location"] for row in rows if row["location"]}) == 1 else location
            raw_metadata = {
                "backend": "python-adc",
                "project_id": active_project,
                "location": result_location,
                "discovery_action": request.action,
                "discovery_target": request.target,
                "target_project": target_project,
            }
            return _result_from_discovery_rows(
                backend="python-adc",
                active_project=active_project,
                location=result_location,
                rows=rows,
                columns=["project_id", "dataset_id", "location"],
                raw_metadata=raw_metadata,
            )

        if request.action == DISCOVERY_LIST_TABLES:
            target_project, dataset_id = parse_project_dataset(request.target)
            dataset_ref = bigquery.DatasetReference(target_project, dataset_id)
            dataset = client.get_dataset(dataset_ref)
            rows = []
            for table in client.list_tables(dataset_ref, max_results=max_rows):
                rows.append(
                    {
                        "project_id": target_project,
                        "dataset_id": dataset_id,
                        "table_id": table.table_id,
                        "table_type": getattr(table, "table_type", "") or "",
                    }
                )
            raw_metadata = {
                "backend": "python-adc",
                "project_id": active_project,
                "location": getattr(dataset, "location", None) or location,
                "discovery_action": request.action,
                "discovery_target": request.target,
                "target_project": target_project,
                "target_dataset": dataset_id,
            }
            return _result_from_discovery_rows(
                backend="python-adc",
                active_project=active_project,
                location=raw_metadata["location"],
                rows=rows,
                columns=["project_id", "dataset_id", "table_id", "table_type"],
                raw_metadata=raw_metadata,
            )

        if request.action == DISCOVERY_DESCRIBE:
            target_project, dataset_id, table_id = parse_project_dataset_table(request.target)
            table_ref = bigquery.TableReference(bigquery.DatasetReference(target_project, dataset_id), table_id)
            table = client.get_table(table_ref)
            rows = [
                {
                    "field_name": field.name,
                    "field_type": field.field_type,
                    "mode": field.mode,
                    "description": field.description or "",
                }
                for field in table.schema
            ]
            raw_metadata = {
                "backend": "python-adc",
                "project_id": active_project,
                "location": getattr(table, "location", None) or location,
                "discovery_action": request.action,
                "discovery_target": request.target,
                "target_project": target_project,
                "target_dataset": dataset_id,
                "target_table": table_id,
                "table_type": getattr(table, "table_type", None),
            }
            return _result_from_discovery_rows(
                backend="python-adc",
                active_project=active_project,
                location=raw_metadata["location"],
                rows=rows,
                columns=["field_name", "field_type", "mode", "description"],
                raw_metadata=raw_metadata,
            )
    except Exception as exc:
        raise BigQuerySetupError(f"Python BigQuery discovery failed: {exc}") from exc

    raise ValueError(f"Unsupported discovery action: {request.action}")


def _run_bq_cli_json(cmd: list[str]) -> Any:
    proc = subprocess.run(cmd, capture_output=True, text=True, check=False)
    if proc.returncode != 0:
        detail = proc.stderr.strip() or proc.stdout.strip() or "Unknown bq CLI error."
        raise BigQuerySetupError(f"bq CLI discovery failed: {detail}")
    parsed = _try_parse_json(proc.stdout)
    return parsed if parsed is not None else proc.stdout.strip()


def execute_discovery_with_bq_cli(
    request: DiscoveryRequest,
    *,
    max_rows: int,
    location: str,
) -> QueryExecutionResult:
    if not shutil.which("bq"):
        raise BigQuerySetupError("bq CLI is not installed.")

    active_project = get_gcloud_project()
    if not active_project:
        raise BigQuerySetupError("No Google Cloud project configured for the bq CLI.")

    if request.action == DISCOVERY_LIST_DATASETS:
        payload = _run_bq_cli_json(
            [
                "bq",
                f"--project_id={active_project}",
                "ls",
                "--datasets",
                "--format=prettyjson",
                request.target,
            ]
        )
        items = payload if isinstance(payload, list) else []
        rows = [
            {
                "project_id": item.get("datasetReference", {}).get("projectId", request.target),
                "dataset_id": item.get("datasetReference", {}).get("datasetId", ""),
                "location": item.get("location", "") or "",
            }
            for item in items[:max_rows]
        ]
        unique_locations = {row["location"] for row in rows if row["location"]}
        result_location = rows[0]["location"] if len(unique_locations) == 1 and rows else location
        raw_metadata = {
            "backend": "bq-cli",
            "project_id": active_project,
            "location": result_location,
            "discovery_action": request.action,
            "discovery_target": request.target,
            "raw_response": payload,
        }
        return _result_from_discovery_rows(
            backend="bq-cli",
            active_project=active_project,
            location=result_location,
            rows=rows,
            columns=["project_id", "dataset_id", "location"],
            raw_metadata=raw_metadata,
        )

    if request.action == DISCOVERY_LIST_TABLES:
        target_project, dataset_id = parse_project_dataset(request.target)
        dataset_meta = _run_bq_cli_json(
            [
                "bq",
                f"--project_id={active_project}",
                "show",
                "--format=prettyjson",
                f"{target_project}:{dataset_id}",
            ]
        )
        payload = _run_bq_cli_json(
            [
                "bq",
                f"--project_id={active_project}",
                "ls",
                "--format=prettyjson",
                f"{target_project}:{dataset_id}",
            ]
        )
        items = payload if isinstance(payload, list) else []
        rows = [
            {
                "project_id": item.get("tableReference", {}).get("projectId", target_project),
                "dataset_id": item.get("tableReference", {}).get("datasetId", dataset_id),
                "table_id": item.get("tableReference", {}).get("tableId", ""),
                "table_type": item.get("type", "") or "",
            }
            for item in items[:max_rows]
        ]
        raw_metadata = {
            "backend": "bq-cli",
            "project_id": active_project,
            "location": dataset_meta.get("location", location) if isinstance(dataset_meta, dict) else location,
            "discovery_action": request.action,
            "discovery_target": request.target,
            "raw_response": payload,
        }
        return _result_from_discovery_rows(
            backend="bq-cli",
            active_project=active_project,
            location=raw_metadata["location"],
            rows=rows,
            columns=["project_id", "dataset_id", "table_id", "table_type"],
            raw_metadata=raw_metadata,
        )

    if request.action == DISCOVERY_DESCRIBE:
        target_project, dataset_id, table_id = parse_project_dataset_table(request.target)
        payload = _run_bq_cli_json(
            [
                "bq",
                f"--project_id={active_project}",
                "show",
                "--format=prettyjson",
                f"{target_project}:{dataset_id}.{table_id}",
            ]
        )
        schema_fields = payload.get("schema", {}).get("fields", []) if isinstance(payload, dict) else []
        rows = [
            {
                "field_name": field.get("name", ""),
                "field_type": field.get("type", ""),
                "mode": field.get("mode", ""),
                "description": field.get("description", "") or "",
            }
            for field in schema_fields
        ]
        raw_metadata = {
            "backend": "bq-cli",
            "project_id": active_project,
            "location": payload.get("location", location) if isinstance(payload, dict) else location,
            "discovery_action": request.action,
            "discovery_target": request.target,
            "raw_response": payload,
        }
        return _result_from_discovery_rows(
            backend="bq-cli",
            active_project=active_project,
            location=raw_metadata["location"],
            rows=rows,
            columns=["field_name", "field_type", "mode", "description"],
            raw_metadata=raw_metadata,
        )

    raise ValueError(f"Unsupported discovery action: {request.action}")
