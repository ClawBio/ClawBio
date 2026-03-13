"""Tests for the knowledge-guide CLI."""
from __future__ import annotations
import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "knowledge_guide.py"


class TestDemoMode:
    """Test --demo flag runs without network and produces output."""

    def test_demo_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "report.html").exists()
        assert (tmp_path / "result.json").exists()

    def test_demo_result_json_envelope(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        data = json.loads((tmp_path / "result.json").read_text())
        assert data["skill"] == "knowledge-guide"
        assert "summary" in data
        assert "data" in data


class TestQueryMode:
    """Test --query flag with cache."""

    def test_query_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--query", "variant calling",
             "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()


class TestSkillMode:
    """Test --skill flag for Learn More lookups."""

    def test_skill_lookup(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--skill", "pharmgx-reporter",
             "--output", str(tmp_path)],
            capture_output=True, text=True, timeout=30,
        )
        assert result.returncode == 0, f"STDERR: {result.stderr}"
        assert (tmp_path / "result.json").exists()
