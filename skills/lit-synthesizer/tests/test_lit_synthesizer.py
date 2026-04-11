"""
Tests for ClawBio Lit Synthesizer skill.
Run with: pytest skills/lit-synthesizer/tests/test_lit_synthesizer.py -v
"""

import json
import sys
import tempfile
from pathlib import Path

import pytest

# Add skill directory to path
sys.path.insert(0, str(Path(__file__).parent.parent))
from lit_synthesizer import (
    DEMO_PAPERS,
    build_citation_graph,
    generate_report,
)


# ── Citation graph tests ───────────────────────────────────────────────────────

class TestBuildCitationGraph:
    def test_returns_nodes_and_edges_keys(self):
        graph = build_citation_graph(DEMO_PAPERS)
        assert "nodes" in graph
        assert "edges" in graph

    def test_node_count_matches_papers(self):
        graph = build_citation_graph(DEMO_PAPERS)
        assert len(graph["nodes"]) == len(DEMO_PAPERS)

    def test_nodes_have_required_fields(self):
        graph = build_citation_graph(DEMO_PAPERS)
        for node in graph["nodes"]:
            assert "id" in node
            assert "title" in node
            assert "source" in node
            assert "year" in node

    def test_edges_reference_known_nodes(self):
        graph = build_citation_graph(DEMO_PAPERS)
        node_ids = {n["id"] for n in graph["nodes"]}
        for edge in graph["edges"]:
            assert "from" in edge
            assert "to" in edge
            # 'to' should be a PMID from our demo set or an external reference
            assert edge["from"] in node_ids

    def test_empty_papers_returns_empty_graph(self):
        graph = build_citation_graph([])
        assert graph["nodes"] == []
        assert graph["edges"] == []

    def test_paper_without_pmid_uses_doi_as_id(self):
        paper = {
            "source": "bioRxiv",
            "pmid": None,
            "title": "Test Paper",
            "year": "2025",
            "doi": "10.1101/test.001",
            "citations": [],
        }
        graph = build_citation_graph([paper])
        assert graph["nodes"][0]["id"] == "10.1101/test.001"


# ── Report generation tests ────────────────────────────────────────────────────

class TestGenerateReport:
    def test_creates_report_md(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            assert (out / "report.md").exists()

    def test_creates_results_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            assert (out / "results.json").exists()

    def test_creates_citation_graph_json(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            assert (out / "citation_graph.json").exists()

    def test_creates_csv_table(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            assert (out / "tables" / "papers.csv").exists()

    def test_creates_reproducibility_bundle(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            assert (out / "reproducibility" / "commands.sh").exists()
            assert (out / "reproducibility" / "environment.yml").exists()
            assert (out / "reproducibility" / "checksums.sha256").exists()

    def test_results_json_is_valid(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            data = json.loads((out / "results.json").read_text())
            assert isinstance(data, list)
            assert len(data) == len(DEMO_PAPERS)

    def test_report_contains_query(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR genome editing", DEMO_PAPERS, graph, out)
            report_text = (out / "report.md").read_text()
            assert "CRISPR genome editing" in report_text

    def test_report_contains_disclaimer(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            report_text = (out / "report.md").read_text()
            assert "research tool" in report_text.lower()

    def test_report_contains_paper_titles(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            report_text = (out / "report.md").read_text()
            for paper in DEMO_PAPERS:
                assert paper["title"][:30] in report_text

    def test_checksums_cover_all_files(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph(DEMO_PAPERS)
            generate_report("CRISPR", DEMO_PAPERS, graph, out)
            checksums_text = (out / "reproducibility" / "checksums.sha256").read_text()
            # report.md should be checksummed
            assert "report.md" in checksums_text

    def test_empty_papers_does_not_crash(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            out = Path(tmpdir)
            graph = build_citation_graph([])
            generate_report("empty query", [], graph, out)
            assert (out / "report.md").exists()


# ── Demo data integrity tests ──────────────────────────────────────────────────

class TestDemoData:
    def test_demo_papers_have_required_fields(self):
        required = ["source", "title", "authors", "journal", "year", "abstract", "doi", "url", "citations"]
        for paper in DEMO_PAPERS:
            for field in required:
                assert field in paper, f"Demo paper missing field: {field}"

    def test_demo_papers_sources_are_valid(self):
        valid_sources = {"PubMed", "bioRxiv"}
        for paper in DEMO_PAPERS:
            assert paper["source"] in valid_sources

    def test_demo_papers_authors_are_list(self):
        for paper in DEMO_PAPERS:
            assert isinstance(paper["authors"], list)

    def test_demo_papers_citations_are_list(self):
        for paper in DEMO_PAPERS:
            assert isinstance(paper["citations"], list)
