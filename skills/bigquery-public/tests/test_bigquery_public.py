from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from decimal import Decimal
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
sys.path.insert(0, str(SKILL_DIR))

import bigquery_backends as backends  # noqa: E402
import bigquery_public as skill  # noqa: E402


from clawbio import cli as clawbio_runner  # noqa: E402

_CATALOG_SPEC = importlib.util.spec_from_file_location("generate_catalog_module", PROJECT_ROOT / "scripts" / "generate_catalog.py")
assert _CATALOG_SPEC and _CATALOG_SPEC.loader
generate_catalog = importlib.util.module_from_spec(_CATALOG_SPEC)
_CATALOG_SPEC.loader.exec_module(generate_catalog)


def test_validate_read_only_sql_accepts_select_and_with():
    select_sql = "SELECT 1 AS ok;"
    with_sql = """
    -- comment
    WITH example AS (SELECT 1 AS value)
    SELECT value FROM example
    """
    case_sql = "SELECT CASE WHEN TRUE THEN 1 ELSE 0 END AS ok"
    assert skill.validate_read_only_sql(select_sql) == "SELECT 1 AS ok"
    normalized = skill.validate_read_only_sql(with_sql)
    assert normalized.lstrip().startswith("-- comment")
    assert "WITH example" in normalized
    assert skill.validate_read_only_sql(case_sql) == case_sql


@pytest.mark.parametrize(
    ("query", "expected"),
    [
        ("DELETE FROM dataset.table", "Only read-only SELECT/WITH queries are supported."),
        ("SELECT 1; SELECT 2", "Multiple SQL statements are not supported."),
        ("WITH x AS (SELECT 1) CREATE TABLE t AS SELECT * FROM x", "Unsupported SQL keyword detected"),
        ("SELECT EXECUTE IMMEDIATE 'SELECT 1'", "Unsupported SQL keyword detected: EXECUTE IMMEDIATE"),
        ("WITH x AS (SELECT 1) BEGIN SELECT * FROM x END", "Unsupported SQL keyword detected: BEGIN"),
        ("SELECT 1 AS x ASSERT TRUE", "Unsupported SQL keyword detected: ASSERT"),
        ("SELECT 1 AS x WHILE TRUE DO SELECT 2 END WHILE", "Unsupported SQL keyword detected: WHILE"),
    ],
)
def test_validate_read_only_sql_rejects_unsafe_queries(query: str, expected: str):
    with pytest.raises(skill.QueryValidationError, match=expected):
        skill.validate_read_only_sql(query)


def test_parse_scalar_param_valid_types():
    string_param = skill.parse_scalar_param("gene=STRING:TP53")
    int_param = skill.parse_scalar_param("limit=INT64:5")
    float_param = skill.parse_scalar_param("p=FLOAT64:0.01")
    bool_param = skill.parse_scalar_param("strict=BOOL:true")
    numeric_param = skill.parse_scalar_param("threshold=NUMERIC:12345678901234567890.123456789")

    assert string_param.value == "TP53"
    assert int_param.value == 5
    assert float_param.value == pytest.approx(0.01)
    assert bool_param.value is True
    assert numeric_param.value == Decimal("12345678901234567890.123456789")


def test_parse_scalar_param_invalid_format_and_type():
    with pytest.raises(ValueError, match="Expected name=type:value"):
        skill.parse_scalar_param("broken")
    with pytest.raises(ValueError, match="Unsupported parameter type"):
        skill.parse_scalar_param("gene=ARRAY:TP53")
    with pytest.raises(ValueError, match="must not contain newlines or NUL bytes"):
        skill.parse_scalar_param("gene=STRING:TP53\nBRCA1")


def test_generate_catalog_parser_folds_block_scalar_description():
    payload = """---
name: demo
description: >-
  First line
  second line
tags: [one, two]
---
body
"""
    parsed = generate_catalog.parse_yaml_frontmatter(payload)
    assert parsed["description"] == "First line second line"
    assert parsed["tags"] == ["one", "two"]


def test_ensure_output_dir_ready_rejects_non_empty_dir(tmp_path: Path):
    output_dir = tmp_path / "occupied"
    output_dir.mkdir()
    (output_dir / "stale.txt").write_text("old", encoding="utf-8")
    with pytest.raises(ValueError, match="already exists and is not empty"):
        skill._ensure_output_dir_ready(output_dir)


def test_resolve_run_plan_wraps_preview_and_preserves_source_query(tmp_path: Path):
    args = skill.build_parser().parse_args(
        [
            "--query",
            "SELECT * FROM demo_table",
            "--preview",
            "5",
            "--output",
            str(tmp_path / "preview"),
        ]
    )
    plan = skill._resolve_run_plan(args)
    assert plan.mode == "query"
    assert plan.source_query == "SELECT * FROM demo_table"
    assert plan.effective_query == "SELECT * FROM (\nSELECT * FROM demo_table\n) AS clawbio_preview LIMIT 5"
    assert any("SELECT *" in warning for warning in plan.warnings)


def test_resolve_run_plan_wraps_count_only_query(tmp_path: Path):
    args = skill.build_parser().parse_args(
        [
            "--query",
            "SELECT gene FROM demo_table LIMIT 10",
            "--count-only",
            "--output",
            str(tmp_path / "count"),
        ]
    )
    plan = skill._resolve_run_plan(args)
    assert plan.effective_query == "SELECT COUNT(*) AS row_count FROM (\nSELECT gene FROM demo_table LIMIT 10\n) AS clawbio_count"


@pytest.mark.parametrize(
    "argv",
    [
        ["--query", "SELECT 1", "--list-datasets", "isb-cgc"],
        ["--demo", "--preview", "5"],
        ["--query", "SELECT 1", "--preview", "5", "--count-only"],
        ["--list-datasets", "isb-cgc", "--dry-run"],
    ],
)
def test_resolve_run_plan_rejects_invalid_mode_combinations(tmp_path: Path, argv: list[str]):
    args = skill.build_parser().parse_args([*argv, "--output", str(tmp_path / "invalid")])
    with pytest.raises(ValueError):
        skill._resolve_run_plan(args)


def test_execute_query_prefers_python_backend(monkeypatch: pytest.MonkeyPatch):
    expected = skill.QueryExecutionResult(
        backend="python-adc",
        project_id="demo-project",
        location="US",
        query="SELECT 1",
        dry_run=False,
        rows=[{"example": 1}],
        columns=["example"],
        estimated_bytes_processed=10,
        total_bytes_processed=10,
        row_count=1,
        job_id="python-job",
        raw_metadata={"backend": "python-adc"},
    )

    monkeypatch.setattr(skill, "execute_with_python_client", lambda **_: expected)
    monkeypatch.setattr(skill, "execute_with_bq_cli", lambda **_: (_ for _ in ()).throw(AssertionError("CLI fallback should not run")))

    result = skill.execute_query(
        query="SELECT 1",
        location="US",
        max_rows=10,
        max_bytes_billed=1000,
        parameters=[],
        dry_run=False,
    )
    assert result.backend == "python-adc"


def test_execute_query_falls_back_to_bq_cli(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        skill,
        "execute_with_python_client",
        lambda **_: (_ for _ in ()).throw(skill.BigQuerySetupError("ADC unavailable")),
    )
    monkeypatch.setattr(
        skill,
        "execute_with_bq_cli",
        lambda **_: skill.QueryExecutionResult(
            backend="bq-cli",
            project_id="demo-project",
            location="US",
            query="SELECT 1",
            dry_run=True,
            rows=[],
            columns=[],
            estimated_bytes_processed=123,
            total_bytes_processed=None,
            row_count=0,
            job_id="cli-job",
            raw_metadata={"backend": "bq-cli"},
        ),
    )

    result = skill.execute_query(
        query="SELECT 1",
        location="US",
        max_rows=10,
        max_bytes_billed=1000,
        parameters=[],
        dry_run=True,
    )
    assert result.backend == "bq-cli"
    assert result.estimated_bytes_processed == 123


def test_execute_query_reports_setup_message_when_all_backends_fail(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        skill,
        "execute_with_python_client",
        lambda **_: (_ for _ in ()).throw(skill.BigQuerySetupError("ADC unavailable")),
    )
    monkeypatch.setattr(
        skill,
        "execute_with_bq_cli",
        lambda **_: (_ for _ in ()).throw(skill.BigQuerySetupError("invalid_grant")),
    )
    monkeypatch.setattr(skill, "_get_gcloud_project", lambda: "demo-project")

    with pytest.raises(skill.BigQuerySetupError, match="gcloud auth login"):
        skill.execute_query(
            query="SELECT 1",
            location="US",
            max_rows=10,
            max_bytes_billed=1000,
            parameters=[],
            dry_run=False,
        )


def test_bq_cli_query_prefers_stdin_transport(monkeypatch: pytest.MonkeyPatch):
    calls: list[dict[str, object]] = []

    class Proc:
        def __init__(self, returncode: int, stdout: str, stderr: str = ""):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

    def fake_run(cmd, capture_output, text, check, input=None):
        calls.append({"cmd": cmd, "input": input})
        return Proc(0, '[{"example": 1}]')

    monkeypatch.setattr(backends.shutil, "which", lambda name: "/usr/bin/bq")
    monkeypatch.setattr(backends, "get_gcloud_project", lambda: "demo-project")
    monkeypatch.setattr(backends.subprocess, "run", fake_run)

    result = backends._execute_with_bq_cli_once(
        query="SELECT 1 AS example",
        location="US",
        max_rows=5,
        max_bytes_billed=1000,
        parameters=[],
        dry_run=False,
        project_id=None,
    )
    assert len(calls) == 1
    assert calls[0]["input"] == "SELECT 1 AS example"
    assert result.raw_metadata["query_transport"] == "stdin"


def test_python_backend_preserves_exact_numeric_params():
    class FakeBigQueryModule:
        captured: list[tuple[str, str, object]] = []

        @staticmethod
        def ScalarQueryParameter(name, type_name, value):
            FakeBigQueryModule.captured.append((name, type_name, value))
            return (name, type_name, value)

    params = [skill.parse_scalar_param("threshold=NUMERIC:12345678901234567890.123456789")]
    built = backends._build_python_query_parameters(params, FakeBigQueryModule)
    assert built == [("threshold", "NUMERIC", Decimal("12345678901234567890.123456789"))]
    assert FakeBigQueryModule.captured[0][2] == Decimal("12345678901234567890.123456789")


def test_bq_cli_query_stays_on_stdin_when_bq_returns_an_error(monkeypatch: pytest.MonkeyPatch):
    calls: list[dict[str, object]] = []

    class Proc:
        def __init__(self, returncode: int, stdout: str, stderr: str = ""):
            self.returncode = returncode
            self.stdout = stdout
            self.stderr = stderr

    def fake_run(cmd, capture_output, text, check, input=None):
        calls.append({"cmd": cmd, "input": input})
        return Proc(2, "", "ERROR: query argument is required")

    monkeypatch.setattr(backends.shutil, "which", lambda name: "/usr/bin/bq")
    monkeypatch.setattr(backends, "get_gcloud_project", lambda: "demo-project")
    monkeypatch.setattr(backends.subprocess, "run", fake_run)

    with pytest.raises(backends.BigQuerySetupError, match="query argument is required"):
        backends._execute_with_bq_cli_once(
            query="SELECT 1 AS example",
            location="US",
            max_rows=5,
            max_bytes_billed=1000,
            parameters=[],
            dry_run=False,
            project_id=None,
        )
    assert len(calls) == 1
    assert calls[0]["input"] == "SELECT 1 AS example"


def test_execute_discovery_prefers_python_backend(monkeypatch: pytest.MonkeyPatch):
    expected = skill.QueryExecutionResult(
        backend="python-adc",
        project_id="demo-project",
        location="US",
        query="",
        dry_run=False,
        rows=[{"project_id": "isb-cgc", "dataset_id": "TCGA_bioclin_v0", "location": "US"}],
        columns=["project_id", "dataset_id", "location"],
        estimated_bytes_processed=None,
        total_bytes_processed=None,
        row_count=1,
        job_id=None,
        raw_metadata={"backend": "python-adc"},
    )
    monkeypatch.setattr(skill, "execute_discovery_with_python_client", lambda request, **_: expected)
    monkeypatch.setattr(
        skill,
        "execute_discovery_with_bq_cli",
        lambda request, **_: (_ for _ in ()).throw(AssertionError("CLI fallback should not run")),
    )

    result = skill.execute_discovery(
        skill.DiscoveryRequest(action=skill.DISCOVERY_LIST_DATASETS, target="isb-cgc"),
        max_rows=10,
        location="US",
    )
    assert result.backend == "python-adc"


def test_execute_discovery_falls_back_to_bq_cli(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        skill,
        "execute_discovery_with_python_client",
        lambda request, **_: (_ for _ in ()).throw(skill.BigQuerySetupError("ADC unavailable")),
    )
    monkeypatch.setattr(
        skill,
        "execute_discovery_with_bq_cli",
        lambda request, **_: skill.QueryExecutionResult(
            backend="bq-cli",
            project_id="demo-project",
            location="US",
            query="",
            dry_run=False,
            rows=[{"field_name": "case_barcode", "field_type": "STRING", "mode": "NULLABLE", "description": ""}],
            columns=["field_name", "field_type", "mode", "description"],
            estimated_bytes_processed=None,
            total_bytes_processed=None,
            row_count=1,
            job_id=None,
            raw_metadata={"backend": "bq-cli"},
        ),
    )

    result = skill.execute_discovery(
        skill.DiscoveryRequest(action=skill.DISCOVERY_DESCRIBE, target="isb-cgc.TCGA_bioclin_v0.Clinical"),
        max_rows=10,
        location="US",
    )
    assert result.backend == "bq-cli"


def test_execute_discovery_reports_setup_message_when_all_backends_fail(monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        skill,
        "execute_discovery_with_python_client",
        lambda request, **_: (_ for _ in ()).throw(skill.BigQuerySetupError("ADC unavailable")),
    )
    monkeypatch.setattr(
        skill,
        "execute_discovery_with_bq_cli",
        lambda request, **_: (_ for _ in ()).throw(skill.BigQuerySetupError("invalid_grant")),
    )
    monkeypatch.setattr(skill, "_get_gcloud_project", lambda: "demo-project")

    with pytest.raises(skill.BigQuerySetupError, match="gcloud auth login"):
        skill.execute_discovery(
            skill.DiscoveryRequest(action=skill.DISCOVERY_LIST_DATASETS, target="isb-cgc"),
            max_rows=10,
            location="US",
        )


def test_demo_main_creates_expected_outputs(tmp_path: Path):
    exit_code = skill.main(["--demo", "--output", str(tmp_path)])
    assert exit_code == 0
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "result.json").exists()
    assert (tmp_path / "tables" / "results.csv").exists()
    assert (tmp_path / "reproducibility" / "commands.sh").exists()
    assert (tmp_path / "reproducibility" / "query.sql").exists()
    assert (tmp_path / "reproducibility" / "job_metadata.json").exists()
    assert (tmp_path / "reproducibility" / "provenance.json").exists()

    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["summary"]["mode"] == "demo"
    assert payload["summary"]["backend"] == "demo-fixture"
    assert "query" in payload["data"]
    assert payload["data"]["disclaimer"].startswith("ClawBio is a research and educational tool")

    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "Execution mode" in report
    assert "demo-fixture" in report
    assert "bigquery-public-data.samples.shakespeare" in report


def test_demo_main_respects_max_rows_and_dry_run(tmp_path: Path):
    exit_code = skill.main(
        ["--demo", "--output", str(tmp_path), "--max-rows", "2", "--dry-run", "--location", "EU"]
    )
    assert exit_code == 0
    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["summary"]["location"] == "EU"
    assert payload["summary"]["row_count"] == 0
    assert payload["summary"]["dry_run"] is True


def test_query_main_writes_preview_and_provenance_outputs(
    tmp_path: Path, monkeypatch: pytest.MonkeyPatch
):
    def fake_execute_query(**kwargs):
        return skill.QueryExecutionResult(
            backend="python-adc",
            project_id="demo-project",
            location=kwargs["location"],
            query=kwargs["query"],
            dry_run=kwargs["dry_run"],
            rows=[{"example": 1}],
            columns=["example"],
            estimated_bytes_processed=10,
            total_bytes_processed=10,
            row_count=1,
            job_id="job-123",
            raw_metadata={"backend": "python-adc", "job_id": "job-123"},
        )

    monkeypatch.setattr(skill, "execute_query", fake_execute_query)

    exit_code = skill.main(
        [
            "--query",
            "SELECT * FROM demo_table",
            "--preview",
            "3",
            "--paper",
            "doi:10.1038/example",
            "--note",
            "Use public bulk RNA-seq only",
            "--note",
            "Preview before download",
            "--output",
            str(tmp_path),
        ]
    )
    assert exit_code == 0

    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["data"]["source_query"] == "SELECT * FROM demo_table"
    assert payload["data"]["effective_query"] == "SELECT * FROM (\nSELECT * FROM demo_table\n) AS clawbio_preview LIMIT 3"
    assert payload["data"]["paper_reference"] == "doi:10.1038/example"
    assert payload["data"]["notes"] == ["Use public bulk RNA-seq only", "Preview before download"]
    assert payload["data"]["warnings"]

    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "## Source Query" in report
    assert "## Effective Query" in report
    assert "## Provenance" in report
    assert "## Warnings" in report

    provenance = json.loads((tmp_path / "reproducibility" / "provenance.json").read_text(encoding="utf-8"))
    assert provenance["source_query"] == "SELECT * FROM demo_table"
    assert provenance["effective_query"] == "SELECT * FROM (\nSELECT * FROM demo_table\n) AS clawbio_preview LIMIT 3"
    assert provenance["notes"] == ["Use public bulk RNA-seq only", "Preview before download"]


def test_discovery_main_writes_expected_outputs(tmp_path: Path, monkeypatch: pytest.MonkeyPatch):
    monkeypatch.setattr(
        skill,
        "execute_discovery",
        lambda request, **_: skill.QueryExecutionResult(
            backend="python-adc",
            project_id="demo-project",
            location="US",
            query="",
            dry_run=False,
            rows=[{"project_id": "isb-cgc", "dataset_id": "TCGA_bioclin_v0", "location": "US"}],
            columns=["project_id", "dataset_id", "location"],
            estimated_bytes_processed=None,
            total_bytes_processed=None,
            row_count=1,
            job_id=None,
            raw_metadata={"backend": "python-adc", "discovery_action": request.action},
        ),
    )

    exit_code = skill.main(
        [
            "--list-datasets",
            "isb-cgc",
            "--paper",
            "/tmp/paper.pdf",
            "--note",
            "Inspect datasets before writing SQL",
            "--output",
            str(tmp_path),
        ]
    )
    assert exit_code == 0

    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["summary"]["mode"] == "discovery"
    assert payload["summary"]["discovery_action"] == skill.DISCOVERY_LIST_DATASETS
    assert payload["data"]["discovery"]["target"] == "isb-cgc"

    report = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "## Discovery Request" in report
    assert "isb-cgc" in report


def test_runner_demo_executes_bigquery_skill(tmp_path: Path):
    output_dir = tmp_path / "runner_demo"
    result = clawbio_runner.run_skill(
        skill_name="bigquery",
        demo=True,
        output_dir=str(output_dir),
    )
    assert result["success"] is True
    assert (output_dir / "report.md").exists()
    assert (output_dir / "result.json").exists()


def test_runner_security_filter_allows_bigquery_flags_and_ignores_unknown(tmp_path: Path):
    output_dir = tmp_path / "runner_filtered"
    result = clawbio_runner.run_skill(
        skill_name="bigquery",
        demo=True,
        output_dir=str(output_dir),
        extra_args=["--location", "EU", "--max-rows", "2", "--paper", "demo-paper", "--note", "keep", "--bogus", "nope"],
    )
    assert result["success"] is True
    payload = json.loads((output_dir / "result.json").read_text(encoding="utf-8"))
    assert payload["summary"]["location"] == "EU"
    assert payload["summary"]["max_rows"] == 2
    assert payload["data"]["paper_reference"] == "demo-paper"
    assert payload["data"]["notes"] == ["keep"]


def test_runner_security_filter_preserves_allowed_value_starting_with_dash(
    monkeypatch: pytest.MonkeyPatch,
    tmp_path: Path,
):
    captured: dict[str, object] = {}

    class Proc:
        returncode = 0
        stdout = ""
        stderr = ""

    def fake_run(cmd, capture_output, text, timeout, cwd):
        captured["cmd"] = cmd
        return Proc()

    monkeypatch.setattr(clawbio_runner.subprocess, "run", fake_run)

    result = clawbio_runner.run_skill(
        skill_name="bigquery",
        output_dir=str(tmp_path / "runner_query"),
        extra_args=["--query", "-- comment\nSELECT 1", "--dry-run"],
    )

    assert result["success"] is True
    cmd = captured["cmd"]
    assert "--query" in cmd
    assert "-- comment\nSELECT 1" in cmd


def test_cli_subprocess_demo_round_trip(tmp_path: Path):
    output_dir = tmp_path / "subprocess_demo"
    proc = subprocess.run(
        [
            sys.executable,
            str(PROJECT_ROOT / "clawbio.py"),
            "run",
            "bigquery",
            "--demo",
            "--output",
            str(output_dir),
        ],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
        check=False,
    )
    assert proc.returncode == 0, proc.stderr
    assert (output_dir / "report.md").exists()
