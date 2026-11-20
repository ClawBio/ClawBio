"""Smoke tests for turingdb-graph.

These tests do NOT require a running TuringDB server. They verify:

- SKILL.md parses as valid YAML frontmatter + body
- All three demo files exist and are readable
- The CLI parses --help and each subcommand without errors
- The COHORT_QUERIES dict contains the expected analyses
- The HTTP wrapper importable

Live integration tests (actually spinning up TuringDB and running the three demos)
are out of scope for CI; run them locally with:

    python turingdb_graph.py --demo cohort --out /tmp/tg-cohort
    python turingdb_graph.py --demo pathway --out /tmp/tg-pathway
    python turingdb_graph.py --demo antibody --out /tmp/tg-antibody
"""
from __future__ import annotations

import importlib.util
import subprocess
import sys
from pathlib import Path

SKILL_ROOT = Path(__file__).resolve().parent.parent


def test_skill_md_parses() -> None:
    import yaml

    skill_md = (SKILL_ROOT / "SKILL.md").read_text()
    assert skill_md.startswith("---\n"), "SKILL.md must start with YAML frontmatter"
    _, front, body = skill_md.split("---", 2)
    meta = yaml.safe_load(front)

    assert meta["name"] == "turingdb-graph"
    assert meta["license"] == "MIT"
    assert "metadata" in meta
    md = meta["metadata"]
    assert md["domain"] == "graph-analytics"
    assert md["version"] == "0.1.0"
    assert "inputs" in md and "outputs" in md
    assert "openclaw" in md
    assert md["openclaw"]["emoji"] == "\U0001f578"

    for section in ("## Domain Decisions", "## Safety", "## Agent Boundary",
                    "## Core Capabilities", "## CLI Reference", "## Gotchas"):
        assert section in body, f"SKILL.md body missing {section!r}"


def test_demo_datasets_exist() -> None:
    for name in ("cohort.csv", "pathway.gml", "antibody.csv"):
        p = SKILL_ROOT / "demo" / name
        assert p.exists(), f"Missing demo dataset: {p}"
        assert p.stat().st_size > 0, f"Empty demo dataset: {p}"


def test_demo_csv_is_anonymous() -> None:
    """Safety Rule 2: no real names / identifiers in the demo cohort."""
    cohort = (SKILL_ROOT / "demo" / "cohort.csv").read_text()
    # Generic placeholder names only — all names should match 'Patient NNN'
    import csv
    from io import StringIO
    rows = list(csv.DictReader(StringIO(cohort)))
    for row in rows:
        assert row["name"].startswith("Patient "), \
            f"Demo cohort must use 'Patient NNN' placeholders, got: {row['name']!r}"


def test_cli_help() -> None:
    result = subprocess.run(
        [sys.executable, str(SKILL_ROOT / "turingdb_graph.py"), "--help"],
        capture_output=True, text=True, timeout=30,
    )
    assert result.returncode == 0, result.stderr
    out = result.stdout
    for token in ("--build", "--query", "--analyse-cohort", "--demo",
                  "--stop-server", "--host", "--data-dir"):
        assert token in out, f"--help missing flag: {token}"


def test_cohort_queries_present() -> None:
    spec = importlib.util.spec_from_file_location(
        "turingdb_graph", SKILL_ROOT / "turingdb_graph.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["turingdb_graph"] = mod  # required so dataclass annotations resolve
    spec.loader.exec_module(mod)

    expected = {
        "patient_count", "age_stats", "junior_count", "senior_count",
        "top_conditions", "top_medications",
        "top_comorbidities", "top_comedications",
    }
    assert set(mod.COHORT_QUERIES) == expected


def test_http_server_importable() -> None:
    # Importing http_server pulls in fastapi — skip if not installed
    try:
        import fastapi  # noqa: F401
    except ImportError:
        import pytest
        pytest.skip("fastapi not installed")

    # Register turingdb_graph first so http_server's "from turingdb_graph import ..." resolves
    tg_spec = importlib.util.spec_from_file_location(
        "turingdb_graph", SKILL_ROOT / "turingdb_graph.py",
    )
    tg_mod = importlib.util.module_from_spec(tg_spec)
    sys.modules["turingdb_graph"] = tg_mod
    tg_spec.loader.exec_module(tg_mod)

    spec = importlib.util.spec_from_file_location(
        "http_server", SKILL_ROOT / "http_server.py",
    )
    mod = importlib.util.module_from_spec(spec)
    sys.modules["http_server"] = mod
    spec.loader.exec_module(mod)

    routes = {r.path for r in mod.app.routes if hasattr(r, "path")}
    assert "/health" in routes
    assert "/build" in routes
    assert "/query" in routes
    assert "/analyse-cohort" in routes
    assert "/demo" in routes


def test_reference_docs_shipped() -> None:
    ref = SKILL_ROOT / "reference"
    for name in ("startup.md", "querying.md", "writing.md",
                 "algorithms.md", "introspection.md", "biomedical.md"):
        assert (ref / name).exists(), f"Missing reference doc: {name}"


# ---------------------------------------------------------------------------
# Demo output snapshot
# ---------------------------------------------------------------------------
# Fails when a code change alters any demo number. If the change is intended,
# run `uv run pytest skills/turingdb-graph/tests --inline-snapshot=fix`; it
# rewrites the snapshot() values below, and the PR diff shows which numbers moved.

import importlib.util  # noqa: E402

import pytest  # noqa: E402
from inline_snapshot import snapshot  # noqa: E402

from clawbio.common import demo_runner  # noqa: E402


@pytest.mark.demo_snapshot
@pytest.mark.skipif(
    importlib.util.find_spec("turingdb") is None,
    reason="needs the turingdb-graph extra (uv sync --extra turingdb-graph)",
)
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("turingdb-graph", tmp_path) == snapshot({'demo-cohort/analysis/report.md': """\
# Cohort analysis: `demo_cohort`

- **Patients**: 20
- **Ages** (n=20): min 14, max 73, mean 48.6, median 51.0
- **Under 18**: 3
- **Over 65**: 6

## Top 10 conditions

| condition | patients |
|---|---|
| Hypertension | 6 |
| Diabetes Type 2 | 4 |
| Asthma | 3 |
| Arthritis | 3 |
| Cancer | 2 |
| Migraine | 2 |

## Top 10 medications

| medication | patients |
|---|---|
| Metformin | 4 |
| Ibuprofen | 3 |
| Amlodipine | 3 |
| Lisinopril | 3 |
| Salbutamol | 3 |
| Sumatriptan | 2 |
| Tamoxifen | 2 |


---

*ClawBio is a research and educational tool. Not a medical device. This output must not be used for clinical decision-making.*
""", 'demo-cohort/analysis/summary.json': """\
{
  "graph": "demo_cohort",
  "patient_count": 20,
  "age_stats": {
    "n": 20,
    "min": 14,
    "max": 73,
    "mean": 48.6,
    "median": 51.0
  },
  "junior_count": 3,
  "senior_count": 6,
  "top_conditions": [
    {
      "condition": "Hypertension",
      "patients": 6
    },
    {
      "condition": "Diabetes Type 2",
      "patients": 4
    },
    {
      "condition": "Asthma",
      "patients": 3
    },
    {
      "condition": "Arthritis",
      "patients": 3
    },
    {
      "condition": "Cancer",
      "patients": 2
    },
    {
      "condition": "Migraine",
      "patients": 2
    }
  ],
  "top_medications": [
    {
      "medication": "Metformin",
      "patients": 4
    },
    {
      "medication": "Ibuprofen",
      "patients": 3
    },
    {
      "medication": "Amlodipine",
      "patients": 3
    },
    {
      "medication": "Lisinopril",
      "patients": 3
    },
    {
      "medication": "Salbutamol",
      "patients": 3
    },
    {
      "medication": "Sumatriptan",
      "patients": 2
    },
    {
      "medication": "Tamoxifen",
      "patients": 2
    }
  ],
  "top_comorbidities": [],
  "top_comedications": []
}\
""", 'demo-cohort/report.md': """\
# Graph build: `demo_cohort`

- **Source**: `cohort.csv` (csv (cohort demo builder))
- **Loader**: `CREATE (cohort demo builder)`
- **Nodes**: 50
- **Edges**: 120

## Nodes by label

| Label | Count |
|-------|------:|
| `Patient` | 20 |
| `BloodType` | 8 |
| `Medication` | 7 |
| `MedicalCondition` | 6 |
| `Doctor` | 4 |
| `Hospital` | 3 |
| `Gender` | 2 |

## Edge types

`HAS`, `TOOK_MEDICATION`, `IS_TREATED_BY`, `IS_TREATED_IN`, `IS`


---

*ClawBio is a research and educational tool. Not a medical device. This output must not be used for clinical decision-making.*
""", 'demo-cohort/summary.json': """\
{
  "graph": "demo_cohort",
  "source": "cohort.csv",
  "format": "csv (cohort demo builder)",
  "loader": "CREATE (cohort demo builder)",
  "nodes": 50,
  "edges": 120,
  "labels": [
    "MedicalCondition",
    "Medication",
    "Doctor",
    "Hospital",
    "BloodType",
    "Gender",
    "Patient"
  ],
  "edge_types": [
    "HAS",
    "TOOK_MEDICATION",
    "IS_TREATED_BY",
    "IS_TREATED_IN",
    "IS"
  ],
  "properties": [
    {
      "id": 0,
      "propertyType": "displayName",
      "valueType": "String"
    },
    {
      "id": 1,
      "propertyType": "id",
      "valueType": "String"
    },
    {
      "id": 2,
      "propertyType": "age",
      "valueType": "Int64"
    }
  ],
  "nodes_by_label": {
    "MedicalCondition": 6,
    "Medication": 7,
    "Doctor": 4,
    "Hospital": 3,
    "BloodType": 8,
    "Gender": 2,
    "Patient": 20
  }
}\
"""})
