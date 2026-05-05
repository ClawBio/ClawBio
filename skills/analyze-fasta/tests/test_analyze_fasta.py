"""Tests for analyze-fasta skill (red/green TDD)."""

import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT = SKILL_DIR / "analyze_fasta.py"
EXAMPLE_DIR = SKILL_DIR / "example_data"


def run_cli(args, expect_code=0):
    """Run the CLI and return (returncode, stdout, stderr)."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), *args],
        capture_output=True,
        text=True,
    )
    if expect_code is not None:
        assert result.returncode == expect_code, (
            f"Expected exit {expect_code}, got {result.returncode}\n"
            f"stdout: {result.stdout}\nstderr: {result.stderr}"
        )
    return result


# ──────────────────────────────────────────────
# Demo mode
# ──────────────────────────────────────────────

def test_demo_runs_successfully(tmp_path):
    """--demo --output DIR produces report.md and result.json."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    assert (out / "report.md").exists(), "report.md not generated"
    assert (out / "result.json").exists(), "result.json not generated"


def test_demo_result_json_is_valid(tmp_path):
    """Demo result.json parses and has required top-level keys."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    assert "sequence_type" in data
    assert "sequences" in data
    assert "summary" in data
    assert data["total_sequences"] >= 1


def test_demo_produces_reproducibility_bundle(tmp_path):
    """Output includes reproducibility/ with commands.sh and run.json."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    repro = out / "reproducibility"
    assert repro.is_dir()
    assert (repro / "commands.sh").exists()
    assert (repro / "run.json").exists()


def test_report_includes_clawbio_disclaimer(tmp_path):
    """report.md must reference the ClawBio medical disclaimer."""
    out = tmp_path / "demo_out"
    run_cli(["--demo", "--output", str(out)])
    md = (out / "report.md").read_text()
    assert "ClawBio" in md
    assert "research" in md.lower() and "medical device" in md.lower()


# ──────────────────────────────────────────────
# Nucleotide and protein basic paths
# ──────────────────────────────────────────────

def test_nucleotide_input_basic(tmp_path):
    """Nucleotide FASTA produces gc_content and composition."""
    out = tmp_path / "nuc_out"
    run_cli(["--input", str(EXAMPLE_DIR / "demo_nucleotide.fasta"), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    assert data["sequence_type"] == "nucleotide"
    seq = data["sequences"][0]
    assert "gc_content" in seq
    assert 0 <= seq["gc_content"] <= 100


def test_protein_input_basic(tmp_path):
    """Protein FASTA produces molecular_weight and pI."""
    out = tmp_path / "prot_out"
    run_cli(["--input", str(EXAMPLE_DIR / "demo_protein.fasta"), "--output", str(out)])
    data = json.loads((out / "result.json").read_text())
    assert data["sequence_type"] == "protein"
    seq = data["sequences"][0]
    assert "molecular_weight_da" in seq
    assert "isoelectric_point" in seq
    assert "gravy" in seq


# ──────────────────────────────────────────────
# Loud-failure validations (no silent degradation)
# ──────────────────────────────────────────────

def test_loud_failure_missing_file(tmp_path):
    """Missing input → exit code 1 with explicit error."""
    out = tmp_path / "out"
    result = run_cli(
        ["--input", str(tmp_path / "nope.fasta"), "--output", str(out)],
        expect_code=1,
    )
    assert "no encontrado" in result.stderr.lower() or "not found" in result.stderr.lower()


def test_loud_failure_empty_file(tmp_path):
    """Empty FASTA file → exit 1, not a fake report."""
    empty = tmp_path / "empty.fasta"
    empty.write_text("")
    out = tmp_path / "out"
    result = run_cli(["--input", str(empty), "--output", str(out)], expect_code=1)
    assert "no se encontraron" in result.stderr.lower() or "no sequences" in result.stderr.lower()
    assert not (out / "result.json").exists(), "Should not write fake result.json on failure"


def test_loud_failure_too_short(tmp_path):
    """Sequence shorter than 10 chars → exit 1."""
    tiny = tmp_path / "tiny.fasta"
    tiny.write_text(">tiny\nACGT\n")
    out = tmp_path / "out"
    result = run_cli(["--input", str(tiny), "--output", str(out)], expect_code=1)
    assert "corta" in result.stderr.lower() or "short" in result.stderr.lower()


def test_loud_failure_too_many_ns(tmp_path):
    """Sequence with >50% Ns → exit 1, no silent fake report."""
    ns = tmp_path / "all_n.fasta"
    ns.write_text(">mostly_n\n" + "N" * 100 + "ACGTACGT\n")
    out = tmp_path / "out"
    result = run_cli(["--input", str(ns), "--output", str(out)], expect_code=1)
    assert "n" in result.stderr.lower()  # mentions Ns


# ──────────────────────────────────────────────
# Backward-compat with legacy CLI
# ──────────────────────────────────────────────

def test_legacy_positional_argument_still_works(tmp_path):
    """Old usage `python script.py file.fasta --json` keeps working."""
    result = run_cli([str(EXAMPLE_DIR / "demo_nucleotide.fasta"), "--json"])
    data = json.loads(result.stdout)
    assert data["sequence_type"] == "nucleotide"
