"""
test_multiqc_reporter.py — Tests for the MultiQC skill
=======================================================
All tests run offline — no MultiQC installation required.
subprocess.run is mocked throughout.
"""

from __future__ import annotations

import json
import sys
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(SKILL_DIR))

import multiqc_reporter  # noqa: E402


def test_multiqc_not_found():
    """check_multiqc() exits with code 1 when multiqc not on PATH."""
    with patch("multiqc_reporter.shutil.which", return_value=None):
        with pytest.raises(SystemExit) as exc:
            multiqc_reporter.check_multiqc()
    assert exc.value.code == 1


def test_missing_input_dir(tmp_path):
    """validate_input_dirs() exits with code 1 when a directory does not exist."""
    missing = str(tmp_path / "does_not_exist")
    with pytest.raises(SystemExit) as exc:
        multiqc_reporter.validate_input_dirs([missing])
    assert exc.value.code == 1


def test_cli_args_forwarded(tmp_path):
    """run_multiqc() invokes multiqc with only input dirs and --outdir (defaults)."""
    with patch("multiqc_reporter.subprocess.run") as mock_run:
        mock_run.return_value = MagicMock(returncode=0)
        multiqc_reporter.run_multiqc(
            input_dirs=[tmp_path],
            output_dir=tmp_path,
        )
    cmd = mock_run.call_args[0][0]
    assert "--export" not in cmd
    assert "--data-format" not in cmd
    assert cmd[:2] == ["multiqc", str(tmp_path)]
    assert cmd[-2:] == ["--outdir", str(tmp_path)]


def test_report_md_written(tmp_path):
    """write_report() creates report.md containing the ClawBio disclaimer."""
    multiqc_reporter.write_report(
        output_dir=tmp_path,
        input_dirs=[tmp_path / "input"],
        samples={},
    )
    report = tmp_path / "report.md"
    assert report.exists()
    assert multiqc_reporter.DISCLAIMER in report.read_text()


def test_per_sample_table(tmp_path):
    """parse_general_stats() reads multiqc_data.json and write_report() renders a table."""
    stats_dir = tmp_path / "multiqc_data"
    stats_dir.mkdir()
    fixture = {
        "report_general_stats_data": {
            "fastqc": {
                "SAMPLE_01": {"percent_duplicates": 5.5, "percent_gc": 49},
                "SAMPLE_02": {"percent_duplicates": 15.0, "percent_gc": 50},
            }
        }
    }
    (stats_dir / "multiqc_data.json").write_text(json.dumps(fixture))

    samples = multiqc_reporter.parse_general_stats(tmp_path)

    assert "SAMPLE_01" in samples
    assert "SAMPLE_02" in samples

    multiqc_reporter.write_report(
        output_dir=tmp_path,
        input_dirs=[tmp_path / "input"],
        samples=samples,
    )
    report_text = (tmp_path / "report.md").read_text()
    assert "SAMPLE_01" in report_text
    assert "SAMPLE_02" in report_text


def test_demo_mode(tmp_path):
    """make_demo_fastqc_files() creates 3 SAMPLE_XX dirs each with fastqc_data.txt."""
    multiqc_reporter.make_demo_fastqc_files(tmp_path)
    sample_dirs = [d for d in tmp_path.iterdir() if d.is_dir()]
    assert len(sample_dirs) == 3
    for d in sample_dirs:
        assert (d / "fastqc_data.txt").exists(), f"fastqc_data.txt missing in {d}"
        assert d.name.startswith("SAMPLE_")


def test_write_reproducibility_bundle(tmp_path):
    """write_reproducibility_bundle writes the three reproducibility files."""
    out = tmp_path / "run"
    out.mkdir()
    (out / "report.md").write_text("# mock\n")
    (out / "multiqc_data").mkdir()
    (out / "multiqc_data" / "multiqc_data.json").write_text("{}")

    inp = tmp_path / "qc_in"
    inp.mkdir()
    multiqc_reporter.write_reproducibility_bundle(out, demo=False, input_dirs=[inp])

    repro = out / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "environment.yml").exists()
    assert (repro / "checksums.sha256").exists()
    assert "report.md" in (repro / "checksums.sha256").read_text()
    assert "multiqc_data/multiqc_data.json" in (repro / "checksums.sha256").read_text()
    cmd = (repro / "commands.sh").read_text()
    assert "multiqc-reporter/multiqc_reporter.py" in cmd
    assert "--input" in cmd
    env_yml = (repro / "environment.yml").read_text()
    assert "clawbio-multiqc-reporter" in env_yml
    assert "python=3.11" in env_yml


def test_write_reproducibility_bundle_demo(tmp_path):
    """Demo runs write --demo in commands.sh and omit --input."""
    out = tmp_path / "run"
    out.mkdir()
    (out / "report.md").write_text("# mock\n")
    multiqc_reporter.write_reproducibility_bundle(out, demo=True, input_dirs=[])
    cmd = (out / "reproducibility" / "commands.sh").read_text()
    assert "--demo" in cmd
    assert "--input" not in cmd


def test_cli_command_for_repro_demo(tmp_path):
    """cli_command_for_repro uses --demo when demo=True."""
    cmd = multiqc_reporter.cli_command_for_repro(
        tmp_path / "out",
        demo=True,
        input_dirs=[],
    )
    assert "--demo" in cmd
    assert "--input" not in cmd
