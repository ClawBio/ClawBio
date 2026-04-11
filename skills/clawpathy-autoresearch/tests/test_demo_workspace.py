"""Tests for the pre-built demo workspace."""
from __future__ import annotations

from pathlib import Path

from skills.clawpathy_autoresearch.workspace import load_workspace, validate_workspace
from skills.clawpathy_autoresearch.scorer import load_scorer


DEMO_DIR = Path(__file__).resolve().parent.parent / "demo_workspace"


def test_demo_workspace_is_valid():
    errors = validate_workspace(DEMO_DIR)
    assert errors == [], f"Demo workspace validation failed: {errors}"


def test_demo_workspace_loads():
    ws = load_workspace(DEMO_DIR)
    assert ws.name == "GWAS Paper Reproduction"
    assert ws.max_iterations == 80
    assert ws.early_stop_n == 5


def test_demo_ground_truth_has_variants():
    ws = load_workspace(DEMO_DIR)
    targets = ws.ground_truth.get("targets", [])
    assert len(targets) >= 3, "Demo should have at least 3 papers"
    for t in targets:
        assert "paper_id" in t
        assert "lead_variants" in t


def test_demo_scorer_runs():
    ws = load_workspace(DEMO_DIR)
    score_fn = load_scorer(ws.scorer_path)

    # Perfect output should score 0
    perfect = {"papers": {}}
    for t in ws.ground_truth["targets"]:
        perfect["papers"][t["paper_id"]] = {
            "variants_found": t["lead_variants"],
            "total_loci_reported": t["total_loci"],
        }
    result = score_fn(perfect, ws.ground_truth)
    assert result == 0.0


def test_demo_scorer_imperfect_output():
    ws = load_workspace(DEMO_DIR)
    score_fn = load_scorer(ws.scorer_path)

    # Empty output should score > 0
    empty = {"papers": {}}
    result = score_fn(empty, ws.ground_truth)
    assert result > 0.0


def test_demo_skill_md_exists():
    ws = load_workspace(DEMO_DIR)
    skill_path = ws.skill_dir / "SKILL.md"
    assert skill_path.exists()
    content = skill_path.read_text()
    assert "## Workflow" in content
