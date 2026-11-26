"""
Tests for pathway-enricher skill.

Covers:
  - Gene-file parsing (comments, blank lines, comma-separated)
  - Demo mode (offline using mocked Enrichr responses)
  - Report and table generation logic
  - result.json structure validation
"""
from __future__ import annotations

import json
import os
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest
from hypothesis import given, settings
from hypothesis import strategies as st

# ---------------------------------------------------------------------------
# Path setup — allow importing the skill module without installing
# ---------------------------------------------------------------------------
SKILL_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = SKILL_DIR.parent.parent
sys.path.insert(0, str(PROJECT_ROOT))

# Import the skill via importlib (no __init__.py, folder name has a hyphen)
import importlib.util as _ilu

_spec = _ilu.spec_from_file_location("pathway_enricher", SKILL_DIR / "pathway_enricher.py")
_mod = _ilu.module_from_spec(_spec)  # type: ignore[arg-type]
_spec.loader.exec_module(_mod)  # type: ignore[union-attr]
pe = _mod


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

MOCK_ENRICHR_RESPONSE = {
    "KEGG_2021_Human": [
        [1, "Alzheimer disease", 0.0001, 3.5, 24.5, ["APP", "APOE", "PSEN1"], 0.002, 0.0001, 0.002],
        [2, "Neurodegeneration", 0.001, 2.8, 15.3, ["APP", "CLU"], 0.01, 0.001, 0.01],
        [3, "Amyloid fiber formation", 0.01, 1.5, 4.2, ["APP"], 0.08, 0.01, 0.08],
    ]
}

FAKE_USER_LIST_RESPONSE = {"userListId": 12345, "shortId": "abc123"}


def _make_mock_requests():
    """Return a mock `requests` module with Enrichr endpoints stubbed."""
    mock_post = MagicMock()
    mock_post.return_value.raise_for_status = MagicMock()
    mock_post.return_value.json.return_value = FAKE_USER_LIST_RESPONSE

    mock_get = MagicMock()
    mock_get.return_value.raise_for_status = MagicMock()
    mock_get.return_value.json.return_value = MOCK_ENRICHR_RESPONSE

    mock_req = MagicMock()
    mock_req.post = mock_post
    mock_req.get = mock_get
    return mock_req


# ---------------------------------------------------------------------------
# Unit tests — parse_gene_file
# ---------------------------------------------------------------------------

def test_parse_gene_file_basic(tmp_path: Path) -> None:
    """Parses plain one-per-line gene file correctly."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("APOE\nBIN1\nCLU\n")
    result = pe.parse_gene_file(gene_file)
    assert result == ["APOE", "BIN1", "CLU"]


def test_parse_gene_file_comments_and_blanks(tmp_path: Path) -> None:
    """Skips comments (#) and blank lines."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("# header\nAPOE\n\nBIN1\n# another comment\nCLU\n")
    result = pe.parse_gene_file(gene_file)
    assert result == ["APOE", "BIN1", "CLU"]


def test_parse_gene_file_comma_separated(tmp_path: Path) -> None:
    """Handles comma-separated genes on one line."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("APOE, BIN1, CLU\n")
    result = pe.parse_gene_file(gene_file)
    assert result == ["APOE", "BIN1", "CLU"]


def test_parse_gene_file_deduplication(tmp_path: Path) -> None:
    """Removes duplicates while preserving order."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("APOE\nBIN1\nAPOE\nCLU\n")
    result = pe.parse_gene_file(gene_file)
    assert result == ["APOE", "BIN1", "CLU"]


def test_parse_gene_file_lowercases_normalised(tmp_path: Path) -> None:
    """Gene symbols are uppercased."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("apoe\nbin1\n")
    result = pe.parse_gene_file(gene_file)
    assert result == ["APOE", "BIN1"]


# ---------------------------------------------------------------------------
# Unit tests — Enrichr API helpers (mocked)
# ---------------------------------------------------------------------------

def test_post_gene_list_uses_multipart_form_data() -> None:
    """_post_gene_list follows Enrichr's multipart upload contract."""
    mock_req = _make_mock_requests()
    with patch.object(pe, "requests", mock_req):
        uid = pe._post_gene_list(["APOE", "BIN1"])
    assert uid == "12345"
    mock_req.post.assert_called_once_with(
        f"{pe.ENRICHR_BASE}/addList",
        files={
            "list": (None, "APOE\nBIN1"),
            "description": (None, "ClawBio gene set"),
        },
        timeout=30,
    )


def test_query_library_parses_rows() -> None:
    """_query_library correctly parses Enrichr result rows."""
    mock_req = _make_mock_requests()
    with patch.object(pe, "requests", mock_req):
        rows = pe._query_library("12345", "KEGG_2021_Human")
    assert len(rows) == 3
    first = rows[0]
    assert first["term"] == "Alzheimer disease"
    assert first["gene_count"] == 3
    assert pytest.approx(first["adj_pvalue"], abs=1e-6) == 0.002


def test_run_enrichr_returns_all_keys(tmp_path: Path) -> None:
    """run_enrichr returns a dict with one key per requested database."""
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = MOCK_ENRICHR_RESPONSE
    with patch.object(pe, "requests", mock_req), patch.object(pe, "time") as mock_time:
        mock_time.sleep = MagicMock()
        result, failures = pe.run_enrichr(["APOE", "BIN1"], ["kegg"])
    assert "kegg" in result
    assert isinstance(result["kegg"], list)
    assert failures == {}


def test_wikipathways_library_id_is_the_name_enrichr_serves() -> None:
    """Enrichr serves 'WikiPathway_2023_Human' (singular); 'WikiPathways_2023_Human'
    does not exist and Enrichr answers it with {} (checked live 2026-11-26)."""
    assert pe.LIBRARIES["wikipathways"][0] == "WikiPathway_2023_Human"


def test_query_library_unknown_library_raises() -> None:
    """Enrichr answers an unknown library name with {}: that is an error, not zero terms."""
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = {}
    with patch.object(pe, "requests", mock_req):
        with pytest.raises(pe.EnrichrLibraryError, match="WikiPathways_2023_Human"):
            pe._query_library("12345", "WikiPathways_2023_Human")


_ROW = st.tuples(
    st.integers(1, 10_000),
    st.text(min_size=1, max_size=40),
    st.floats(0, 1),
    st.floats(-50, 50),
    st.floats(0, 1e5),
    st.lists(st.sampled_from(["APOE", "APP", "BIN1", "CLU", "TREM2"]), min_size=1, max_size=5, unique=True),
    st.floats(0, 1),
    st.floats(0, 1),
    st.floats(0, 1),
).map(list)


@settings(deadline=None)
@given(
    requested=st.sampled_from([lib for lib, _ in pe.LIBRARIES.values()]),
    other=st.dictionaries(st.text(max_size=30), st.lists(_ROW, max_size=3), max_size=4),
)
def test_query_library_never_reads_another_librarys_rows(requested: str, other: dict) -> None:
    """Without the requested library's key the query raises, whatever other libraries'
    rows the response holds: evidence for one library is never read from another."""
    other.pop(requested, None)
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = other
    with patch.object(pe, "requests", mock_req):
        with pytest.raises(pe.EnrichrLibraryError):
            pe._query_library("12345", requested)


@settings(deadline=None)
@given(rows=st.lists(_ROW, max_size=20), noise=st.lists(_ROW, min_size=1, max_size=3))
def test_query_library_parses_exactly_its_own_rows(rows: list, noise: list) -> None:
    """One parsed result per row of the requested library, in order, with its own genes."""
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = {"Reactome_2022": rows, "KEGG_2021_Human": noise}
    with patch.object(pe, "requests", mock_req):
        parsed = pe._query_library("12345", "Reactome_2022")
    assert [r["term"] for r in parsed] == [row[1] for row in rows]
    assert [r["adj_pvalue"] for r in parsed] == [row[6] for row in rows]
    assert all(r["gene_count"] == len(row[5]) for r, row in zip(parsed, rows))


def test_query_library_known_library_with_no_hits_is_empty_not_error() -> None:
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = {"Reactome_2022": []}
    with patch.object(pe, "requests", mock_req):
        assert pe._query_library("12345", "Reactome_2022") == []


def test_run_enrichr_records_failed_library_separately() -> None:
    """A library that fails is reported in failures, not silently as an empty result."""
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = MOCK_ENRICHR_RESPONSE  # KEGG only
    with patch.object(pe, "requests", mock_req), patch.object(pe, "time") as mock_time:
        mock_time.sleep = MagicMock()
        result, failures = pe.run_enrichr(["APOE"], ["kegg", "wikipathways"])
    assert len(result["kegg"]) == 3
    assert "wikipathways" in failures
    assert "WikiPathway_2023_Human" in failures["wikipathways"]


def test_report_and_json_mark_failed_library() -> None:
    failures = {"reactome": "Enrichr returned no entry for library 'Reactome_2022'"}
    all_results = {"kegg": [], "reactome": []}
    report = pe.build_report(["APOE"], all_results, ["kegg", "reactome"], {}, Path("."),
                             "2026-01-01T00:00:00Z", failures=failures)
    assert "Query failed" in report
    assert "Reactome_2022" in report
    result = pe.build_result_json(["APOE"], all_results, ["kegg", "reactome"],
                                  "2026-01-01T00:00:00Z", failures=failures)
    assert result["failed_databases"] == failures
    assert "reactome" not in result["significant_counts"]


def test_main_fails_when_every_library_fails(tmp_path: Path) -> None:
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = {}
    with patch.object(pe, "requests", mock_req), patch.object(pe, "time") as mock_time:
        mock_time.sleep = MagicMock()
        rc = pe.main(["--demo", "--output", str(tmp_path / "out")])
    assert rc == 1


@pytest.mark.network
@pytest.mark.skipif(os.getenv("RUN_LIVE_TESTS") != "1", reason="live tests gated on RUN_LIVE_TESTS=1")
def test_live_every_library_id_exists_in_enrichr() -> None:
    import requests
    stats = requests.get(f"{pe.ENRICHR_BASE}/datasetStatistics", timeout=30).json()["statistics"]
    served = {s["libraryName"] for s in stats}
    missing = [lib for lib, _ in pe.LIBRARIES.values() if lib not in served]
    assert missing == []


# ---------------------------------------------------------------------------
# Unit tests — report & table generation
# ---------------------------------------------------------------------------

def test_table_md_empty() -> None:
    """_table_md returns a graceful message for empty rows."""
    out = pe._table_md([])
    assert "No results" in out


def test_table_md_contains_term() -> None:
    """_table_md includes term names and adj. p-values."""
    rows = [
        {"rank": 1, "term": "Alzheimer disease", "pvalue": 0.001, "adj_pvalue": 0.005,
         "zscore": 3.0, "combined_score": 20.0, "genes": ["APP", "APOE"], "gene_count": 2},
    ]
    out = pe._table_md(rows)
    assert "Alzheimer disease" in out
    assert "5.00e-03" in out or "0.005" in out or "5.0e-03" in out


def test_build_result_json_structure() -> None:
    """build_result_json produces expected keys."""
    genes = ["APOE", "BIN1"]
    rows = [
        {"rank": 1, "term": "Alzheimer disease", "pvalue": 0.001, "adj_pvalue": 0.002,
         "zscore": 3.5, "combined_score": 24.5, "genes": ["APOE"], "gene_count": 1},
    ]
    all_results = {"kegg": rows}
    result = pe.build_result_json(genes, all_results, ["kegg"], "2026-06-01T00:00:00Z")
    assert result["tool"] == "pathway-enricher"
    assert result["input_genes"] == genes
    assert "kegg" in result["top_results"]
    assert result["significant_counts"]["kegg"] == 1


def test_build_report_contains_sections() -> None:
    """build_report includes expected headings and disclaimer."""
    genes = ["APOE", "BIN1"]
    rows = [
        {"rank": 1, "term": "Alzheimer disease", "pvalue": 0.001, "adj_pvalue": 0.002,
         "zscore": 3.5, "combined_score": 24.5, "genes": ["APOE"], "gene_count": 1},
    ]
    all_results = {"kegg": rows}
    report = pe.build_report(genes, all_results, ["kegg"], {}, Path("/tmp"), "2026-06-01T00:00:00Z")
    assert "Pathway Enrichment Report" in report
    assert "KEGG" in report
    assert "Disclaimer" in report or "research and educational" in report


# ---------------------------------------------------------------------------
# Integration-level: demo mode (mocked network)
# ---------------------------------------------------------------------------

def test_demo_mode_creates_outputs(tmp_path: Path) -> None:
    """Demo mode produces report.md, result.json, and table CSVs."""
    mock_req = _make_mock_requests()
    # Return KEGG results for all library queries
    mock_req.get.return_value.json.return_value = MOCK_ENRICHR_RESPONSE

    with patch.object(pe, "requests", mock_req), \
         patch.object(pe, "time") as mock_time, \
         patch.object(pe, "HAS_VIZ", False):  # skip matplotlib in CI
        mock_time.sleep = MagicMock()
        exit_code = pe.main(["--demo", "--output", str(tmp_path), "--databases", "kegg"])

    assert exit_code == 0
    assert (tmp_path / "report.md").exists(), "report.md not created"
    assert (tmp_path / "result.json").exists(), "result.json not created"
    assert (tmp_path / "tables" / "kegg_enrichment.csv").exists(), "kegg CSV not created"

    result = json.loads((tmp_path / "result.json").read_text())
    assert result["tool"] == "pathway-enricher"
    assert len(result["input_genes"]) == 25  # DEMO_GENES


def test_demo_mode_report_not_empty(tmp_path: Path) -> None:
    """Demo mode report.md is non-trivial."""
    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = MOCK_ENRICHR_RESPONSE

    with patch.object(pe, "requests", mock_req), \
         patch.object(pe, "time") as mock_time, \
         patch.object(pe, "HAS_VIZ", False):
        mock_time.sleep = MagicMock()
        pe.main(["--demo", "--output", str(tmp_path), "--databases", "kegg"])

    report = (tmp_path / "report.md").read_text()
    assert len(report) > 200
    assert "KEGG" in report


def test_input_file_mode(tmp_path: Path) -> None:
    """--input mode reads a gene file and produces outputs."""
    gene_file = tmp_path / "genes.txt"
    gene_file.write_text("APOE\nBIN1\nCLU\n")

    mock_req = _make_mock_requests()
    mock_req.get.return_value.json.return_value = MOCK_ENRICHR_RESPONSE

    out_dir = tmp_path / "output"
    with patch.object(pe, "requests", mock_req), \
         patch.object(pe, "time") as mock_time, \
         patch.object(pe, "HAS_VIZ", False):
        mock_time.sleep = MagicMock()
        exit_code = pe.main(["--input", str(gene_file), "--output", str(out_dir), "--databases", "kegg"])

    assert exit_code == 0
    result = json.loads((out_dir / "result.json").read_text())
    assert result["input_genes"] == ["APOE", "BIN1", "CLU"]


def test_missing_input_returns_error(tmp_path: Path) -> None:
    """Passing a non-existent file returns exit code 1."""
    exit_code = pe.main(["--input", str(tmp_path / "nonexistent.txt"),
                         "--output", str(tmp_path / "out")])
    assert exit_code == 1
