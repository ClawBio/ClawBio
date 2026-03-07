"""
test_galaxy_bridge.py — Tests for Galaxy Bridge skill
======================================================
Covers: catalog loading, local search, demo mode, category listing.
All tests run offline — no Galaxy API key required.
"""

from __future__ import annotations

import json
import shutil
import tempfile
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]

# We need to be able to import the bridge module
import sys
sys.path.insert(0, str(SKILL_DIR))
import galaxy_bridge  # noqa: E402


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

SAMPLE_CATALOG = {
    "version": "1.0.0",
    "generated_by": "test",
    "galaxy_url": "https://usegalaxy.org",
    "tool_count": 5,
    "section_count": 3,
    "sections": {"FASTQ Quality Control": 2, "Metagenomics": 2, "RNA-seq": 1},
    "tools": [
        {
            "id": "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1",
            "name": "FastQC",
            "description": "Read Quality reports",
            "version": "0.74+galaxy1",
            "section": "FASTQ Quality Control",
            "edam_topics": ["Sequencing"],
            "edam_operations": ["Sequence composition calculation"],
            "inputs": [{"name": "input_file", "type": "data", "label": "Raw data"}],
            "outputs": [{"name": "report", "format": "html"}],
        },
        {
            "id": "toolshed.g2.bx.psu.edu/repos/devteam/trimmomatic/trimmomatic/0.39",
            "name": "Trimmomatic",
            "description": "Flexible read trimming tool for Illumina NGS data",
            "version": "0.39",
            "section": "FASTQ Quality Control",
            "edam_topics": [],
            "edam_operations": [],
            "inputs": [],
            "outputs": [],
        },
        {
            "id": "toolshed.g2.bx.psu.edu/repos/iuc/kraken2/kraken2/2.1.3",
            "name": "Kraken2",
            "description": "Assign taxonomic labels to short DNA sequences",
            "version": "2.1.3",
            "section": "Metagenomics",
            "edam_topics": ["Metagenomics"],
            "edam_operations": ["Taxonomic classification"],
            "inputs": [{"name": "input_sequences", "type": "data", "label": "Input sequences"}],
            "outputs": [{"name": "output", "format": "tabular"}],
        },
        {
            "id": "toolshed.g2.bx.psu.edu/repos/iuc/metaphlan/metaphlan/4.0",
            "name": "MetaPhlAn",
            "description": "Metagenomic phylogenetic analysis",
            "version": "4.0",
            "section": "Metagenomics",
            "edam_topics": ["Metagenomics"],
            "edam_operations": [],
            "inputs": [],
            "outputs": [],
        },
        {
            "id": "toolshed.g2.bx.psu.edu/repos/iuc/deseq2/deseq2/1.42.0",
            "name": "DESeq2",
            "description": "Differential gene expression analysis",
            "version": "1.42.0",
            "section": "RNA-seq",
            "edam_topics": ["RNA-Seq"],
            "edam_operations": ["Differential gene expression profiling"],
            "inputs": [],
            "outputs": [],
        },
    ],
}


@pytest.fixture
def catalog_file(tmp_path):
    """Write sample catalog to a temp file and patch galaxy_bridge to use it."""
    cat_path = tmp_path / "galaxy_catalog.json"
    cat_path.write_text(json.dumps(SAMPLE_CATALOG, indent=2), encoding="utf-8")
    original = galaxy_bridge.CATALOG_PATH
    galaxy_bridge.CATALOG_PATH = cat_path
    yield cat_path
    galaxy_bridge.CATALOG_PATH = original


@pytest.fixture
def tmp_output(tmp_path):
    out = tmp_path / "output"
    out.mkdir()
    return out


# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------


class TestCatalogLoading:
    def test_catalog_loads(self, catalog_file):
        """galaxy_catalog.json parses correctly."""
        catalog = galaxy_bridge.load_catalog()
        assert catalog["version"] == "1.0.0"
        assert catalog["tool_count"] == 5
        assert len(catalog["tools"]) == 5

    def test_catalog_missing_exits(self, tmp_path):
        """Missing catalog file causes sys.exit."""
        original = galaxy_bridge.CATALOG_PATH
        galaxy_bridge.CATALOG_PATH = tmp_path / "nonexistent.json"
        with pytest.raises(SystemExit):
            galaxy_bridge.load_catalog()
        galaxy_bridge.CATALOG_PATH = original


class TestSearch:
    def test_search_finds_fastqc(self, catalog_file):
        """Search for 'fastqc' finds FastQC."""
        catalog = galaxy_bridge.load_catalog()
        results = galaxy_bridge.search_catalog("fastqc", catalog)
        assert len(results) >= 1
        assert results[0]["name"] == "FastQC"

    def test_search_metagenomics(self, catalog_file):
        """Search for 'metagenomics' finds Kraken2 and MetaPhlAn."""
        catalog = galaxy_bridge.load_catalog()
        results = galaxy_bridge.search_catalog("metagenomics", catalog)
        names = [r["name"] for r in results]
        assert "Kraken2" in names
        assert "MetaPhlAn" in names

    def test_search_no_results(self, catalog_file):
        """Search for nonexistent term returns empty."""
        catalog = galaxy_bridge.load_catalog()
        results = galaxy_bridge.search_catalog("xyznonexistent123", catalog)
        assert results == []

    def test_search_case_insensitive(self, catalog_file):
        """Search is case-insensitive."""
        catalog = galaxy_bridge.load_catalog()
        r1 = galaxy_bridge.search_catalog("FASTQC", catalog)
        r2 = galaxy_bridge.search_catalog("fastqc", catalog)
        assert len(r1) == len(r2)

    def test_search_by_description(self, catalog_file):
        """Search finds tools by description content."""
        catalog = galaxy_bridge.load_catalog()
        results = galaxy_bridge.search_catalog("differential gene expression", catalog)
        names = [r["name"] for r in results]
        assert "DESeq2" in names


class TestCategories:
    def test_list_categories(self, catalog_file):
        """Categories are counted correctly."""
        catalog = galaxy_bridge.load_catalog()
        cats = galaxy_bridge.list_categories(catalog)
        assert cats["FASTQ Quality Control"] == 2
        assert cats["Metagenomics"] == 2
        assert cats["RNA-seq"] == 1


class TestToolDetails:
    def test_find_by_exact_id(self, catalog_file):
        """Find tool by exact ID."""
        catalog = galaxy_bridge.load_catalog()
        tool = galaxy_bridge.get_tool_details(
            "toolshed.g2.bx.psu.edu/repos/devteam/fastqc/fastqc/0.74+galaxy1",
            catalog,
        )
        assert tool is not None
        assert tool["name"] == "FastQC"

    def test_find_by_partial_id(self, catalog_file):
        """Find tool by partial ID match."""
        catalog = galaxy_bridge.load_catalog()
        tool = galaxy_bridge.get_tool_details("fastqc", catalog)
        assert tool is not None
        assert tool["name"] == "FastQC"

    def test_not_found(self, catalog_file):
        """Returns None for unknown tool."""
        catalog = galaxy_bridge.load_catalog()
        tool = galaxy_bridge.get_tool_details("nonexistent_tool_xyz", catalog)
        assert tool is None


class TestDemoMode:
    def test_demo_runs(self, tmp_output):
        """Demo mode runs without API key and produces output."""
        result = galaxy_bridge.run_demo(tmp_output)
        assert result["mode"] == "demo"
        assert result["tool"] == "fastqc"

    def test_demo_creates_html(self, tmp_output):
        """Demo mode creates HTML output file."""
        galaxy_bridge.run_demo(tmp_output)
        html = tmp_output / "fastqc_demo_output.html"
        assert html.exists()
        content = html.read_text()
        assert "FastQC" in content

    def test_demo_creates_result_json(self, tmp_output):
        """Demo mode creates result.json."""
        galaxy_bridge.run_demo(tmp_output)
        rj = tmp_output / "result.json"
        assert rj.exists()
        data = json.loads(rj.read_text())
        assert data["qc_modules"]["per_base_sequence_quality"] == "PASS"

    def test_demo_creates_reproducibility(self, tmp_output):
        """Demo mode creates reproducibility bundle."""
        galaxy_bridge.run_demo(tmp_output)
        assert (tmp_output / "reproducibility" / "commands.sh").exists()

    def test_demo_creates_fastq(self, tmp_output):
        """Demo mode generates synthetic FASTQ."""
        galaxy_bridge.run_demo(tmp_output)
        fq = galaxy_bridge.DEMO_DIR / "demo_reads.fq"
        assert fq.exists()
        lines = fq.read_text().strip().split("\n")
        # 1000 reads × 4 lines each
        assert len(lines) == 4000


class TestReportWriter:
    def test_write_report_demo(self, tmp_output):
        """Report writer handles demo results."""
        result = {"mode": "demo", "qc_modules": {"quality": "PASS"}}
        rp = galaxy_bridge.write_report(tmp_output, "fastqc", result)
        assert rp.exists()
        content = rp.read_text()
        assert "Galaxy Bridge Report" in content
        assert "Disclaimer" in content

    def test_write_report_success(self, tmp_output):
        """Report writer handles success results."""
        result = {"status": "success", "outputs": ["/tmp/out.html"]}
        rp = galaxy_bridge.write_report(tmp_output, "fastqc", result)
        content = rp.read_text()
        assert "successfully" in content
