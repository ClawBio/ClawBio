"""
Unit tests for pubmed_summariser.py — no network calls required.

Tests the CLI logic, max_results clamping, demo mode behaviour, and
HTML report generation using a mock for pubmed_api.fetch_papers.
"""

import json
import sys
from io import StringIO
from pathlib import Path
from unittest.mock import patch

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

FIXTURES = Path(__file__).parent / "fixtures"

SAMPLE_PAPERS = json.loads((FIXTURES / "sample_papers.json").read_text())


# ── _clamp_max_results ───────────────────────────────────────────────────────


def test_clamp_within_limit():
    """Values <= 50 should be returned unchanged."""
    from pubmed_summariser import _clamp_max_results
    assert _clamp_max_results(10) == 10
    assert _clamp_max_results(50) == 50


def test_clamp_above_limit(capsys):
    """Values > 50 should be clamped to 50 and print a warning."""
    from pubmed_summariser import _clamp_max_results
    result = _clamp_max_results(200)
    assert result == 50
    captured = capsys.readouterr()
    assert "[warning]" in captured.out
    assert "capped at 50" in captured.out


# ── _format_terminal_summary ─────────────────────────────────────────────────


def test_terminal_summary_contains_query():
    """Terminal summary header must include the query string."""
    from pubmed_summariser import _format_terminal_summary
    output = _format_terminal_summary("BRCA1", SAMPLE_PAPERS)
    assert "BRCA1" in output


def test_terminal_summary_contains_paper_title():
    """Each paper title must appear in the terminal output."""
    from pubmed_summariser import _format_terminal_summary
    output = _format_terminal_summary("BRCA1", SAMPLE_PAPERS)
    assert SAMPLE_PAPERS[0]["title"] in output


def test_terminal_summary_empty_results():
    """Zero results should produce a 'No results found' message."""
    from pubmed_summariser import _format_terminal_summary
    output = _format_terminal_summary("UNKNOWN_GENE_XYZ", [])
    assert "No results found" in output


# ── _build_html_report ────────────────────────────────────────────────────────


def test_html_report_contains_query():
    """HTML report must include the query term."""
    from pubmed_summariser import _build_html_report
    html = _build_html_report("BRCA1", SAMPLE_PAPERS)
    assert "BRCA1" in html


def test_html_report_contains_pmid_link():
    """HTML report must contain a link to each paper's PubMed URL."""
    from pubmed_summariser import _build_html_report
    html = _build_html_report("BRCA1", SAMPLE_PAPERS)
    assert SAMPLE_PAPERS[0]["url"] in html


def test_html_report_contains_disclaimer():
    """HTML report must contain the medical disclaimer."""
    from pubmed_summariser import _build_html_report
    html = _build_html_report("BRCA1", SAMPLE_PAPERS)
    assert "research and educational tool" in html.lower() or "disclaimer" in html.lower()


# ── main integration (mocked) ─────────────────────────────────────────────────


def test_main_demo_mode_ignores_query(tmp_path, capsys):
    """When --demo is set alongside --query, demo wins and a note is printed."""
    from pubmed_summariser import main
    with patch("pubmed_summariser.pubmed_api.fetch_papers", return_value=SAMPLE_PAPERS) as mock_fetch:
        main(["--demo", "--query", "TP53", "--output", str(tmp_path)])
        # fetch_papers must have been called with "BRCA1", not "TP53"
        call_query = mock_fetch.call_args[0][0]
        assert call_query == "BRCA1"
    captured = capsys.readouterr()
    assert "[demo mode]" in captured.out


def test_main_writes_report_html(tmp_path):
    """Running with --query must produce report.html in the output directory."""
    from pubmed_summariser import main
    with patch("pubmed_summariser.pubmed_api.fetch_papers", return_value=SAMPLE_PAPERS):
        main(["--query", "BRCA1", "--output", str(tmp_path)])
    assert (tmp_path / "report.html").exists()


def test_main_no_results_exits_cleanly(tmp_path, capsys):
    """Zero API results must not write a report and must print a helpful message."""
    from pubmed_summariser import main
    with patch("pubmed_summariser.pubmed_api.fetch_papers", return_value=[]):
        main(["--query", "UNKNOWNXYZ", "--output", str(tmp_path)])
    assert not (tmp_path / "report.html").exists()
    captured = capsys.readouterr()
    assert "No results found" in captured.out
