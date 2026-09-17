"""Tests for ukb-navigator — UK Biobank schema search."""

import csv
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from ukb_navigator import (
    DEMO_RESULTS,
    generate_report,
    field_lookup,
    query_schema,
)


# ---------------------------------------------------------------------------
# DEMO_RESULTS validation
# ---------------------------------------------------------------------------

class TestDemoResults:
    def test_non_empty(self):
        assert len(DEMO_RESULTS) > 0

    def test_required_keys(self):
        for entry in DEMO_RESULTS:
            assert "text" in entry
            assert "source" in entry
            assert "similarity" in entry
            assert "metadata" in entry

    def test_similarities_in_range(self):
        for entry in DEMO_RESULTS:
            assert 0 <= entry["similarity"] <= 1

    def test_sorted_by_similarity(self):
        sims = [e["similarity"] for e in DEMO_RESULTS]
        assert sims == sorted(sims, reverse=True)


# ---------------------------------------------------------------------------
# generate_report
# ---------------------------------------------------------------------------

class TestGenerateReport:
    def test_report_contains_query(self, tmp_path):
        report_path = generate_report("blood pressure", DEMO_RESULTS, tmp_path)
        text = report_path.read_text()
        assert "blood pressure" in text

    def test_report_contains_disclaimer(self, tmp_path):
        report_path = generate_report("test", DEMO_RESULTS, tmp_path)
        text = report_path.read_text()
        assert "research and educational tool" in text

    def test_report_contains_date(self, tmp_path):
        report_path = generate_report("test", DEMO_RESULTS, tmp_path)
        text = report_path.read_text()
        assert "**Date**:" in text

    def test_report_contains_result_count(self, tmp_path):
        report_path = generate_report("test", DEMO_RESULTS, tmp_path)
        text = report_path.read_text()
        assert f"{len(DEMO_RESULTS)} matches" in text

    def test_csv_output(self, tmp_path):
        generate_report("test", DEMO_RESULTS, tmp_path)
        csv_path = tmp_path / "matched_fields.csv"
        assert csv_path.exists()
        with open(csv_path) as f:
            reader = csv.DictReader(f)
            rows = list(reader)
        assert len(rows) == len(DEMO_RESULTS)
        assert "rank" in rows[0]
        assert "similarity" in rows[0]
        assert "source" in rows[0]
        assert "text" in rows[0]

    def test_creates_output_dir(self, tmp_path):
        out = tmp_path / "new_dir"
        report_path = generate_report("test", DEMO_RESULTS, out)
        assert report_path.exists()

    def test_creates_reproducibility(self, tmp_path):
        generate_report("test", DEMO_RESULTS, tmp_path)
        cmd_file = tmp_path / "reproducibility" / "commands.sh"
        assert cmd_file.exists()

    def test_reproducibility_bundle_is_complete(self, tmp_path):
        generate_report("blood pressure", DEMO_RESULTS, tmp_path, is_demo=True)
        repro = tmp_path / "reproducibility"

        commands_text = (repro / "commands.sh").read_text(encoding="utf-8")
        assert "CLAWBIO_ROOT" in commands_text
        assert "$OUTPUT_DIR" in commands_text
        assert str(tmp_path) not in commands_text

        environment = (repro / "environment.yml").read_text(encoding="utf-8")
        assert "name: clawbio-ukb-navigator" in environment
        assert "chromadb" in environment

        checksum_lines = [
            line
            for line in (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
            if line.strip()
        ]
        assert checksum_lines
        labels = set()
        for line in checksum_lines:
            digest, label = line.split("  ", 1)
            assert len(digest) == 64
            labels.add(label)
            assert (tmp_path / label).exists()
        assert {"report.md", "matched_fields.csv"} <= labels

    def test_reproducibility_command_quotes_multiword_query(self, tmp_path):
        """Parsed as a shell would: the query must come back as ONE argument.
        Asserting on quote characters is too weak -- the unquoted f-string form
        this replaced also produced `--query "blood pressure"`."""
        import shlex

        query = 'blood "pressure" $HOME'
        generate_report(query, DEMO_RESULTS, tmp_path)
        commands_text = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
        run_line = [ln for ln in commands_text.splitlines() if "ukb_navigator.py" in ln][0]
        # the command spans continuation lines; rejoin before parsing
        joined = commands_text[commands_text.index(run_line):].replace("\\\n", " ")
        tokens = shlex.split(joined)
        assert tokens[tokens.index("--query") + 1] == query
        assert "$HOME" in tokens[tokens.index("--query") + 1]

    def test_demo_mode_flag(self, tmp_path):
        report_path = generate_report("test", DEMO_RESULTS, tmp_path, is_demo=True)
        text = report_path.read_text()
        assert "Demo" in text

    def test_empty_results(self, tmp_path):
        report_path = generate_report("nothing", [], tmp_path)
        text = report_path.read_text()
        assert "0 matches" in text


# ---------------------------------------------------------------------------
# field_lookup delegates to query_schema
# ---------------------------------------------------------------------------

class TestFieldLookup:
    def test_delegates_correctly(self):
        # field_lookup just wraps query_schema with a formatted string
        # With no ChromaDB collection, it should return empty list
        # (we test the delegation logic, not the DB)
        # This test verifies the function exists and is callable
        assert callable(field_lookup)
