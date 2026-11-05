from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import Any

# Port of ClawBio 1.x `bigquery-public/bigquery_models.py` for ClawBio 2.0.
# Unchanged except SKILL_DIR/DEMO_DIR: 1.x colocated demo/ next to the script
# itself; this port follows the skills/<name>/{scripts,data}/ layout used
# elsewhere in ClawBio 2.0 (see skills/gwas_lookup, skills/claw_ancestry_pca),
# but the fixtures stay at the existing 1.x location, ../demo relative to scripts/.
SCRIPT_DIR = Path(__file__).resolve().parent
SKILL_DIR = SCRIPT_DIR.parent
DEMO_DIR = SKILL_DIR / "demo"
DEMO_QUERY_PATH = DEMO_DIR / "demo_query.sql"
DEMO_RESULT_PATH = DEMO_DIR / "demo_result.json"

SKILL_NAME = "bigquery-public"
SKILL_VERSION = "0.2.1"
DEFAULT_LOCATION = "US"
DEFAULT_MAX_ROWS = 100
DEFAULT_MAX_BYTES_BILLED = 1_000_000_000

GOOGLE_AUTH_REQUIREMENT = "google-auth>=2,<3"
GOOGLE_BIGQUERY_REQUIREMENT = "google-cloud-bigquery>=3,<4"

DISCOVERY_LIST_DATASETS = "list-datasets"
DISCOVERY_LIST_TABLES = "list-tables"
DISCOVERY_DESCRIBE = "describe"


class BigQuerySetupError(RuntimeError):
    """Raised when BigQuery access is unavailable or misconfigured."""


class QueryValidationError(ValueError):
    """Raised when the SQL is unsafe or unsupported."""


@dataclass
class QueryParameter:
    name: str
    type_name: str
    value: Any
    original: str

    def to_cli_spec(self) -> str:
        return f"{self.name}:{self.type_name}:{self.original}"


@dataclass
class QueryExecutionResult:
    backend: str
    project_id: str | None
    location: str
    query: str
    dry_run: bool
    rows: list[dict[str, Any]]
    columns: list[str]
    estimated_bytes_processed: int | None
    total_bytes_processed: int | None
    row_count: int
    job_id: str | None
    raw_metadata: dict[str, Any]


@dataclass
class DiscoveryRequest:
    action: str
    target: str


@dataclass
class RunPlan:
    mode: str
    query_source: str
    source_query: str | None
    effective_query: str | None
    discovery: DiscoveryRequest | None
    warnings: list[str]
    paper_reference: str | None
    notes: list[str]
