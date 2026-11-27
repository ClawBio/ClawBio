"""Tests for clawbio.common.report — report generation helpers."""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parent.parent.parent))

from clawbio.common.report import (
    generate_report_header,
    generate_report_footer,
    write_result_json,
)


class TestGenerateReportHeader:
    def test_contains_title_and_skill(self):
        header = generate_report_header(title="Test Report", skill_name="test-skill")
        assert "# Test Report" in header
        assert "test-skill" in header

    def test_contains_date(self):
        header = generate_report_header(title="T", skill_name="s")
        assert "**Date**:" in header

    def test_includes_checksums(self, tmp_path):
        f = tmp_path / "input.txt"
        f.write_text("data")
        header = generate_report_header(title="T", skill_name="s", input_files=[f])
        assert "input.txt" in header
        assert "`" in header  # checksum in backticks

    def test_extra_metadata(self):
        header = generate_report_header(
            title="T", skill_name="s", extra_metadata={"Version": "1.0"}
        )
        assert "**Version**: 1.0" in header


class TestGenerateReportFooter:
    def test_contains_disclaimer(self):
        footer = generate_report_footer()
        assert "Disclaimer" in footer
        assert "research and educational tool" in footer


class TestWriteResultJson:
    def test_creates_valid_envelope(self, tmp_path):
        path = write_result_json(
            output_dir=tmp_path,
            skill="test-skill",
            version="0.1.0",
            summary={"status": "ok"},
            data={"results": [1, 2, 3]},
            input_checksum="abc123",
        )
        assert path.exists()
        envelope = json.loads(path.read_text())
        assert envelope["skill"] == "test-skill"
        assert envelope["version"] == "0.1.0"
        assert "completed_at" in envelope
        assert envelope["input_checksum"] == "sha256:abc123"
        assert envelope["summary"] == {"status": "ok"}
        assert envelope["data"] == {"results": [1, 2, 3]}

    def test_creates_output_dir(self, tmp_path):
        path = write_result_json(
            output_dir=tmp_path / "new_dir",
            skill="s",
            version="1",
            summary={},
            data={},
        )
        assert path.exists()

    def test_status_ok_are_opt_in(self, tmp_path):
        """Backward-compat: without status/ok the envelope has neither key
        (unchanged for the ~40 skills that use this helper)."""
        path = write_result_json(
            output_dir=tmp_path, skill="s", version="1", summary={}, data={}
        )
        envelope = json.loads(path.read_text())
        assert "status" not in envelope
        assert "ok" not in envelope

    def test_status_ok_written_when_provided(self, tmp_path):
        """The nf-core pipeline wrappers opt into the shared status/ok contract."""
        path = write_result_json(
            output_dir=tmp_path,
            skill="s",
            version="1",
            summary={},
            data={},
            status="ok",
            ok=True,
        )
        envelope = json.loads(path.read_text())
        assert envelope["status"] == "ok"
        assert envelope["ok"] is True

    def test_empty_checksum(self, tmp_path):
        path = write_result_json(
            output_dir=tmp_path,
            skill="s",
            version="1",
            summary={},
            data={},
        )
        envelope = json.loads(path.read_text())
        assert envelope["input_checksum"] == ""
