"""Tests for workspace loading and validation."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from skills.clawpathy_autoresearch.workspace import (
    Workspace,
    load_workspace,
    validate_workspace,
)


@pytest.fixture
def valid_workspace(tmp_path: Path) -> Path:
    """Create a minimal valid workspace directory."""
    ws = tmp_path / "my_task"
    ws.mkdir()

    task_json = {
        "name": "Test Task",
        "description": "Reproduce test results",
        "max_iterations": 20,
        "early_stop_n": 5,
    }
    (ws / "task.json").write_text(json.dumps(task_json))

    ground_truth = {
        "targets": [
            {"id": "item_1", "value": 42.0},
            {"id": "item_2", "value": 3.14},
        ]
    }
    (ws / "ground_truth.json").write_text(json.dumps(ground_truth))

    scorer_code = '''
def score(skill_output: dict, ground_truth: dict) -> float:
    """Lower is better. 0 = perfect."""
    targets = {t["id"]: t["value"] for t in ground_truth["targets"]}
    outputs = {o["id"]: o["value"] for o in skill_output.get("results", [])}
    if not targets:
        return 0.0
    errors = []
    for tid, tval in targets.items():
        oval = outputs.get(tid, 0.0)
        errors.append(abs(tval - oval) / max(abs(tval), 1e-9))
    return sum(errors) / len(errors)
'''
    (ws / "scorer.py").write_text(scorer_code)

    skill_dir = ws / "skill"
    skill_dir.mkdir()
    (skill_dir / "SKILL.md").write_text(
        "---\nname: test-skill\nversion: 0.1.0\nauthor: test\n"
        "description: A test skill\n---\n\n# Test Skill\n\n"
        "## Workflow\n\n1. Read the input\n2. Produce output\n"
    )

    sources_dir = ws / "sources"
    sources_dir.mkdir()

    return ws


def test_load_workspace_returns_workspace(valid_workspace: Path):
    ws = load_workspace(valid_workspace)
    assert isinstance(ws, Workspace)
    assert ws.name == "Test Task"
    assert ws.description == "Reproduce test results"
    assert ws.max_iterations == 20
    assert ws.early_stop_n == 5


def test_workspace_has_ground_truth(valid_workspace: Path):
    ws = load_workspace(valid_workspace)
    assert ws.ground_truth is not None
    assert len(ws.ground_truth["targets"]) == 2


def test_workspace_has_paths(valid_workspace: Path):
    ws = load_workspace(valid_workspace)
    assert ws.workspace_dir == valid_workspace
    assert ws.skill_dir == valid_workspace / "skill"
    assert ws.sources_dir == valid_workspace / "sources"
    assert ws.scorer_path == valid_workspace / "scorer.py"


def test_validate_workspace_valid(valid_workspace: Path):
    errors = validate_workspace(valid_workspace)
    assert errors == []


def test_validate_workspace_missing_task_json(valid_workspace: Path):
    (valid_workspace / "task.json").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("task.json" in e for e in errors)


def test_validate_workspace_missing_ground_truth(valid_workspace: Path):
    (valid_workspace / "ground_truth.json").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("ground_truth.json" in e for e in errors)


def test_validate_workspace_missing_scorer(valid_workspace: Path):
    (valid_workspace / "scorer.py").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("scorer.py" in e for e in errors)


def test_validate_workspace_missing_skill_md(valid_workspace: Path):
    (valid_workspace / "skill" / "SKILL.md").unlink()
    errors = validate_workspace(valid_workspace)
    assert any("SKILL.md" in e for e in errors)


def test_load_workspace_invalid_raises(valid_workspace: Path):
    (valid_workspace / "task.json").unlink()
    with pytest.raises(FileNotFoundError):
        load_workspace(valid_workspace)


def test_load_workspace_nonexistent_dir():
    with pytest.raises(FileNotFoundError):
        load_workspace(Path("/nonexistent/workspace"))
