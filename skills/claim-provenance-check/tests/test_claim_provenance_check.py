"""
Tests for the Claim Provenance Check skill.

These pin the REFUSAL behaviour, which is the part that matters: the claim this skill makes is
not "it produces a report", it is "it declines correctly" on an invented citation, a missing
citation, or copied text — and that no partial credit is given for a claim that fails any one
of those checks.

Run:
    pytest skills/claim-provenance-check/tests/ -v
or via ClawBio runner:
    python -m pytest skills/claim-provenance-check/tests/ -v
"""

from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest

# Make the skill importable regardless of working directory.
sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from claim_provenance_check import run, run_demo  # noqa: E402
from provenance_core import (  # noqa: E402
    ClaimState,
    CoverageReport,
    EvidenceRow,
    SourceOutcome,
    bind_claim,
    check_no_verbatim,
    extract_tags,
    state_for_claim,
)

ROWS = [
    EvidenceRow("trials:0", "trials", "A phase II study of the agent was terminated for futility.", n=212),
    EvidenceRow("literature:0", "literature", "Elevated wall shear stress correlates with remodelling.", n=48),
    EvidenceRow("ip:0", "ip", "Composition of matter claim granted to a third party in 2019."),
]


@pytest.fixture
def tmp_output(tmp_path: Path) -> Path:
    return tmp_path / "output"


# ---------------------------------------------------------------------------
# Tags
# ---------------------------------------------------------------------------

def test_extract_tags_reads_all_and_lowercases():
    assert extract_tags("a [Trials:0] b [ip:0]") == ["trials:0", "ip:0"]


def test_extract_tags_ignores_non_tags():
    assert extract_tags("see figure [3] and table [x]") == []


# ---------------------------------------------------------------------------
# The binder
# ---------------------------------------------------------------------------

def test_valid_claim_binds_clean():
    assert bind_claim("The trial stopped early [trials:0], and the IP is held elsewhere [ip:0].", ROWS) == []


def test_invented_citation_is_caught():
    issues = bind_claim("Strong signal in cohort [trials:9].", ROWS)
    assert [i.code for i in issues] == ["citation-invalid"]


def test_uncited_claim_is_caught():
    issues = bind_claim("This target is clearly the best option available.", ROWS)
    assert any(i.code == "citation-floor" for i in issues)


def test_verbatim_copying_is_caught():
    # Lifts an 8+ word run straight out of trials:0.
    claim = "Note that a phase II study of the agent was terminated for futility [trials:0]."
    issues = bind_claim(claim, ROWS)
    assert any(i.code == "verbatim-run" for i in issues)


def test_short_claims_do_not_trip_verbatim():
    assert check_no_verbatim("Terminated early.", ROWS) == []


def test_paraphrase_passes():
    claim = "Enrolment halted before completion because the endpoint was not being met [trials:0]."
    assert not any(i.code == "verbatim-run" for i in bind_claim(claim, ROWS))


def test_invented_citation_gets_no_partial_credit_alongside_a_valid_one():
    # One good tag and one bad tag: the whole claim is refused, not half-credited.
    issues = bind_claim("Confirmed in one study [trials:0] and another [trials:9].", ROWS)
    assert any(i.code == "citation-invalid" for i in issues)


# ---------------------------------------------------------------------------
# Claim states
# ---------------------------------------------------------------------------

def test_two_independent_classes_is_supported():
    state, issues = state_for_claim("Halted early [trials:0], and reported elsewhere [literature:0].", ROWS)
    assert state is ClaimState.SUPPORTED and issues == []


def test_single_class_is_contested_not_supported():
    # One source is not corroboration, and pretending otherwise is the quiet failure.
    state, _ = state_for_claim("Halted early [trials:0].", ROWS)
    assert state is ClaimState.CONTESTED


def test_disagreement_is_contested():
    state, _ = state_for_claim(
        "Halted early [trials:0], though the literature reports benefit [literature:0].",
        ROWS,
        disagreeing_tags=["literature:0"],
    )
    assert state is ClaimState.CONTESTED


def test_unbindable_claim_is_no_coverage_not_weak_support():
    state, issues = state_for_claim("A compelling result [trials:9].", ROWS)
    assert state is ClaimState.NO_COVERAGE and issues


# ---------------------------------------------------------------------------
# Coverage report
# ---------------------------------------------------------------------------

def test_coverage_keeps_gaps_visible_with_routes():
    report = CoverageReport(
        "VOL04 synthetic subject",
        [
            SourceOutcome("hemodynamics", True, rows=1),
            SourceOutcome("trials", True, rows=6, cost_usd=0.0012),
            SourceOutcome(
                "device_recalls", False,
                reason="supplier unreachable",
                route="openFDA device recalls, $0.09",
            ),
        ],
    )
    assert report.reached == ["hemodynamics", "trials"]
    assert len(report.gaps) == 1
    rendered = report.render()
    assert "NO COVERAGE" in rendered and "openFDA device recalls" in rendered
    assert report.as_dict()["coverage"] == "2/3"


def test_cost_is_summed_across_reached_classes():
    report = CoverageReport("x", [
        SourceOutcome("a", True, cost_usd=0.005),
        SourceOutcome("b", True, cost_usd=0.0012),
    ])
    assert report.total_cost_usd == 0.0062


# ---------------------------------------------------------------------------
# CLI / end-to-end
# ---------------------------------------------------------------------------

def test_demo_runs(tmp_output: Path) -> None:
    run_demo(tmp_output)


def test_demo_report_generated(tmp_output: Path) -> None:
    run_demo(tmp_output)
    report = tmp_output / "report.md"
    assert report.exists(), "report.md was not created"
    assert report.stat().st_size > 0, "report.md is empty"
    text = report.read_text()
    assert "ClawBio is a research and educational tool" in text
    assert "NO COVERAGE" in text  # the bundled demo case includes a NO COVERAGE claim


def test_demo_result_json_valid(tmp_output: Path) -> None:
    run_demo(tmp_output)
    result_file = tmp_output / "result.json"
    assert result_file.exists(), "result.json was not created"
    data = json.loads(result_file.read_text())
    assert data.get("skill") == "claim-provenance-check"
    assert "coverage" in data and "claims" in data


def test_demo_invented_citation_claim_is_refused(tmp_output: Path) -> None:
    results = run_demo(tmp_output)
    invented = next(c for c in results["claims"] if "population:4" in c["claim"])
    assert invented["state"] == "NO COVERAGE"
    assert any(i["code"] == "citation-invalid" for i in invented["issues"])


def test_demo_uncited_claim_is_refused(tmp_output: Path) -> None:
    results = run_demo(tmp_output)
    uncited = next(c for c in results["claims"] if c["claim"].startswith("This is the strongest"))
    assert uncited["state"] == "NO COVERAGE"
    assert any(i["code"] == "citation-floor" for i in uncited["issues"])


def test_run_rejects_case_missing_required_keys(tmp_path: Path, tmp_output: Path) -> None:
    bad_case = tmp_path / "bad_case.json"
    bad_case.write_text(json.dumps({"subject": "x", "rows": [], "outcomes": []}))
    with pytest.raises(ValueError, match="claims"):
        run(bad_case, tmp_output)
