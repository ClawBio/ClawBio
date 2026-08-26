"""Red/green TDD suite for completeness-aware-caller.

The skill turns (gene search result, completeness estimate, ANI comparisons)
into three-state calls: present / absent / cannot_conclude. Absence may only
be asserted when completeness makes a miss unlikely; ANI-based ranking of two
references may only be asserted when the margin exceeds completeness-induced
drift.
"""

import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
SCRIPT = SKILL_DIR / "completeness_aware_caller.py"

sys.path.insert(0, str(SKILL_DIR))

from completeness_aware_caller import ani_margin_gate, call_gene, parse_completeness


# ---------------------------------------------------------------- gene calls
def test_present_at_any_completeness():
    for c in (0.5, 0.7, 0.99):
        result = call_gene(found=True, completeness=c)
        assert result["status"] == "present"


def test_absent_requires_high_completeness():
    result = call_gene(found=False, completeness=0.98)
    assert result["status"] == "absent"
    assert result["confidence"] == pytest.approx(0.98)


def test_low_completeness_abstains():
    result = call_gene(found=False, completeness=0.70)
    assert result["status"] == "cannot_conclude"
    assert "70" in result["message"]
    assert "cannot" in result["message"].lower()


def test_boundary_exactly_at_threshold_allows_absent():
    result = call_gene(found=False, completeness=0.95)
    assert result["status"] == "absent"


def test_completeness_percent_normalised():
    assert parse_completeness("72.5") == pytest.approx(0.725)
    assert parse_completeness("0.725") == pytest.approx(0.725)


def test_invalid_completeness_raises():
    with pytest.raises(ValueError):
        parse_completeness("150")
    with pytest.raises(ValueError):
        parse_completeness("-0.1")


# ---------------------------------------------------------------- ANI gate
def test_ani_gate_refuses_when_margin_within_noise():
    # Real numbers from the STM815 degradation benchmark: at 50% retention the
    # LB400/J2315 margin (0.70 ANI pts) sits inside 2x the observed drift.
    result = ani_margin_gate(ani_a=81.28, ani_b=80.58, completeness=0.50)
    assert result["decision"] == "cannot_conclude"
    assert result["margin"] == pytest.approx(0.70, abs=0.01)


def test_ani_gate_ranks_when_margin_is_clear():
    result = ani_margin_gate(ani_a=97.0, ani_b=82.0, completeness=0.90)
    assert result["decision"] == "rank"


def test_ani_gate_flags_species_boundary_zone():
    result = ani_margin_gate(ani_a=95.2, ani_b=82.0, completeness=0.70)
    assert result["species_boundary_uncertain"] is True


# ---------------------------------------------------------------- CLI / demo
def test_demo_mode_produces_report_with_refusal(tmp_path):
    out = tmp_path / "demo_out"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    report = (out / "report.md").read_text()
    assert "CANNOT CONCLUDE" in report
    assert "not a medical device" in report
    assert (out / "result.json").exists()
    assert (out / "commands.sh").exists()


def test_demo_result_json_schema(tmp_path):
    out = tmp_path / "demo_json"
    subprocess.run(
        [sys.executable, str(SCRIPT), "--demo", "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    data = json.loads((out / "result.json").read_text())
    assert data["completeness"] == pytest.approx(0.70)
    statuses = {g["gene"]: g["status"] for g in data["gene_calls"]}
    # nodC loses every copy at frag70 -> unprovable absence, must abstain.
    assert statuses["nodC"] == "cannot_conclude"
    assert statuses["nodB"] == "cannot_conclude"
    # nifH survives frag70 only because STM815 carries two copies of it.
    assert statuses["nifH"] == "present"
    assert statuses["nifD"] == "present"
    assert data["ani_gate"]["decision"] in ("rank", "cannot_conclude")


def test_cli_with_real_inputs(tmp_path):
    genes = tmp_path / "genes.json"
    genes.write_text(json.dumps([
        {"gene": "nifH", "found": False},
        {"gene": "nodC", "found": True},
    ]))
    out = tmp_path / "real_out"
    proc = subprocess.run(
        [sys.executable, str(SCRIPT), "--completeness", "0.98",
         "--genes", str(genes), "--output", str(out)],
        capture_output=True, text=True, timeout=120,
    )
    assert proc.returncode == 0, proc.stderr
    data = json.loads((out / "result.json").read_text())
    statuses = {g["gene"]: g["status"] for g in data["gene_calls"]}
    assert statuses["nifH"] == "absent"
    assert statuses["nodC"] == "present"
