"""Shared fixtures for clawbio tests.

Tests do NOT mock the filesystem or state.json — they run the real commands
against the real demo skill in a per-test tmp directory. 1.x got burned by
mock/prod divergence; 2.0 keeps it real.
"""
from __future__ import annotations

from pathlib import Path

import pytest


REPO_ROOT = Path(__file__).resolve().parents[2]
VARIANT_ANNOTATION = Path(__file__).resolve().parent / "fixtures" / "variant_annotation"


@pytest.fixture
def run_cwd(tmp_path: Path, monkeypatch: pytest.MonkeyPatch) -> Path:
    """Chdir into a fresh tmp dir so each test gets its own .clawbio-run/."""
    monkeypatch.chdir(tmp_path)
    return tmp_path


@pytest.fixture
def variant_inputs(run_cwd: Path) -> dict[str, str]:
    """Create the two input files variant_annotation needs and return paths."""
    vcf = run_cwd / "input.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n")
    ref = run_cwd / "ref.fa"
    ref.write_text(">chr1\nACGT\n")
    return {"vcf_path": str(vcf), "reference": str(ref)}
