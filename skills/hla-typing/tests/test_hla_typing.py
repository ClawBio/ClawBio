"""Tests for hla-typing. Red/green TDD: these should fail until implementation is complete."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "hla_typing.py"
DEMO_INPUT = SKILL_DIR / "demo_input.txt"


class TestCLI:
    """CLI interface tests."""

    def test_no_args_exits_nonzero(self):
        result = subprocess.run(
            [sys.executable, str(SCRIPT)],
            capture_output=True, text=True
        )
        assert result.returncode != 0

    def test_demo_mode_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (tmp_path / "report.md").exists()
        assert (tmp_path / "result.json").exists()

    def test_input_mode_produces_output(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT),
             "--input", str(DEMO_INPUT),
             "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        assert (tmp_path / "report.md").exists()

    def test_missing_input_exits_nonzero(self, tmp_path):
        result = subprocess.run(
            [sys.executable, str(SCRIPT),
             "--input", str(tmp_path / "nonexistent.txt"),
             "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode != 0


class TestOutputFormat:
    """Output format validation."""

    def test_result_json_is_valid(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert isinstance(result, dict)
        assert "skill" in result
        assert result["skill"] == "hla-typing"

    def test_report_contains_disclaimer(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        report = (tmp_path / "report.md").read_text()
        assert "not a medical device" in report.lower()

    def test_result_has_variants_count(self, tmp_path):
        subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        result = json.loads((tmp_path / "result.json").read_text())
        assert "variants_processed" in result
        assert result["variants_processed"] > 0


class TestReproducibility:
    """Reproducibility bundle written via clawbio.common.reproducibility."""

    def test_demo_writes_verifiable_bundle(self, tmp_path):
        import hashlib

        result = subprocess.run(
            [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
            capture_output=True, text=True
        )
        assert result.returncode == 0, f"stderr: {result.stderr}"
        repro = tmp_path / "reproducibility"
        for name in ("commands.sh", "environment.yml", "checksums.sha256"):
            assert (repro / name).exists(), f"missing {name}"

        commands = (repro / "commands.sh").read_text()
        assert "skills/hla-typing/hla_typing.py" in commands
        assert "--demo" in commands
        assert str(SKILL_DIR) not in commands

        lines = (repro / "checksums.sha256").read_text().splitlines()
        labels = {line.split("  ", 1)[1] for line in lines}
        assert {"report.md", "result.json", "reproducibility/environment.yml"} <= labels
        for line in lines:
            digest, label = line.split("  ", 1)
            target = tmp_path / label
            assert target.is_file(), f"label does not resolve: {label}"
            assert hashlib.sha256(target.read_bytes()).hexdigest() == digest


class TestDemoData:
    """Demo data integrity."""

    def test_demo_input_exists(self):
        assert DEMO_INPUT.exists(), f"Demo data missing: {DEMO_INPUT}"

    def test_demo_input_has_content(self):
        content = DEMO_INPUT.read_text()
        lines = [l for l in content.splitlines() if l.strip() and not l.startswith("#")]
        assert len(lines) > 0, "Demo input has no data lines"
