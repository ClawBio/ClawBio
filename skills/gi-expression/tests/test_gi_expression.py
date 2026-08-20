"""Tests for gi-expression. Tests marked ``integration`` hit the real GI API."""

from __future__ import annotations

import json
import subprocess
import sys
from pathlib import Path

import pytest

from clawbio.gi.gi_client import read_fasta
from clawbio.gi.gi_runner import (
    EXPRESSION_MAX_BP,
    EXPRESSION_TSS_RADIUS,
    EXPRESSION_WINDOW_BP,
    validate_expression_input,
)

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "gi_expression.py"
DEMO_INPUT = SKILL_DIR / "example_data" / "expression_hbb_k562.fa"


def test_demo_input_exists():
    assert DEMO_INPUT.exists()
    assert DEMO_INPUT.stat().st_size > 1000


def test_demo_fixture_is_exactly_one_window():
    """The bundled fixture must be exactly the 9,198 bp the model scores."""
    _name, seq = read_fasta(DEMO_INPUT)
    assert len(seq) == EXPRESSION_WINDOW_BP


def test_window_constants_match_the_api_contract():
    assert EXPRESSION_WINDOW_BP == 9198
    assert EXPRESSION_TSS_RADIUS == 4599
    assert EXPRESSION_MAX_BP == 500_000


def test_exact_window_without_tss_index_is_accepted():
    assert validate_expression_input("A" * EXPRESSION_WINDOW_BP, None) is None


def test_short_sequence_is_rejected_before_the_api_call():
    problem = validate_expression_input("A" * 600, None)
    assert problem is not None
    assert "600 bp" in problem and "9198 bp" in problem
    assert "TSS" in problem


def test_oversized_sequence_is_rejected():
    problem = validate_expression_input("A" * (EXPRESSION_MAX_BP + 1), 250_000)
    assert problem is not None and "500000 bp" in problem


def test_tss_index_is_required_when_the_sequence_is_not_one_window():
    problem = validate_expression_input("A" * 20_000, None)
    assert problem is not None and "--tss-index" in problem


def test_locus_with_a_valid_tss_index_is_accepted():
    assert validate_expression_input("A" * 20_000, 10_000) is None
    # The inclusive edges of the legal range.
    assert validate_expression_input("A" * 20_000, EXPRESSION_TSS_RADIUS) is None
    assert validate_expression_input("A" * 20_000, 20_000 - EXPRESSION_TSS_RADIUS) is None


def test_tss_index_outside_the_window_is_rejected():
    problem = validate_expression_input("A" * 20_000, 100)
    assert problem is not None
    assert "[4599, 15401]" in problem


def test_cli_rejects_a_short_fasta_without_calling_the_api(tmp_path, monkeypatch):
    """The gate must fire before the client is built, so no key and no network."""
    short = tmp_path / "short.fa"
    short.write_text(">short\n" + "ACGT" * 150 + "\n")
    monkeypatch.delenv("GI_API_KEY", raising=False)
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--input", str(short), "--output", str(tmp_path / "out")],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 1, f"stdout={result.stdout} stderr={result.stderr}"
    assert "9198 bp" in result.stderr
    assert not (tmp_path / "out").exists()


def test_skill_module_imports():
    result = subprocess.run(
        [sys.executable, "-c", "import importlib.util as u, sys; "
         f"s=u.spec_from_file_location('m', r'{SCRIPT}'); m=u.module_from_spec(s); s.loader.exec_module(m); "
         "sys.exit(0 if hasattr(m, 'main') else 1)"],
        capture_output=True, text=True,
    )
    assert result.returncode == 0, f"import failed: {result.stderr}"


@pytest.mark.integration
def test_demo_mode_predicts_high_hbb_k562_expression(tmp_path):
    """HBB in K562 should report HIGH expression (gene-sense FASTA + K562 description).
    A regression here usually means the strand or description wiring broke."""
    result = subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(tmp_path)],
        capture_output=True, text=True, timeout=60,
    )
    assert result.returncode == 0, f"stderr: {result.stderr}"
    body = json.loads((tmp_path / "result.json").read_text())
    log_tpm = body["summary"].get("log_tpm")
    # The floor sits well below the measured gene-sense value and well above
    # the wrong-strand one, because this test exists to catch a strand or
    # plumbing regression, not to pin a prediction. On the bundled fixture
    # gene-sense measured 0.95 and the genomic strand 0.06. Absolute values
    # move when the model checkpoint changes — an earlier `> 2.0` here went
    # stale and failed against production. If this fires, re-measure both
    # strands before touching the number.
    assert log_tpm is not None and log_tpm > 0.3, (
        f"HBB-in-K562 should score well above the wrong-strand level; got {log_tpm}. "
        "Check that the FASTA is RC'd to gene-sense and the default description is reaching the API."
    )
