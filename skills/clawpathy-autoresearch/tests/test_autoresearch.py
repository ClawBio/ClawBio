"""Tests for the autoresearch loop and CLI."""
from __future__ import annotations

import json
from pathlib import Path

import pytest

from skills.clawpathy_autoresearch.autoresearch import (
    AutoResearchLoop,
    ExperimentResult,
    load_experiment_log,
    save_experiment_log,
)


@pytest.fixture
def workspace_dir(tmp_path: Path) -> Path:
    """Create a minimal workspace for testing the loop."""
    ws = tmp_path / "test_workspace"
    ws.mkdir()

    task_json = {
        "name": "Test Task",
        "description": "Test",
        "max_iterations": 3,
        "early_stop_n": 2,
    }
    (ws / "task.json").write_text(json.dumps(task_json))

    ground_truth = {
        "targets": [
            {"id": "item_1", "value": 10.0},
        ]
    }
    (ws / "ground_truth.json").write_text(json.dumps(ground_truth))

    scorer_code = '''
def score(skill_output: dict, ground_truth: dict) -> float:
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
        "---\nname: test-skill\n---\n# Test\n\n## Workflow\n\n1. Do thing\n"
    )

    (ws / "sources").mkdir()
    return ws


def test_experiment_result_serialisation():
    result = ExperimentResult(
        experiment=1,
        score=0.35,
        kept=True,
        label="baseline",
        skill_diff={},
        error_breakdown={},
    )
    d = result.to_dict()
    assert d["experiment"] == 1
    assert d["score"] == 0.35
    assert d["kept"] is True


def test_save_and_load_experiment_log(tmp_path: Path):
    results = [
        ExperimentResult(
            experiment=1, score=0.8, kept=True, label="baseline",
            skill_diff={}, error_breakdown={},
        ),
        ExperimentResult(
            experiment=2, score=0.6, kept=True, label="improved",
            skill_diff={"skill": "modified"}, error_breakdown={},
        ),
    ]
    log_path = tmp_path / "experiment_log.json"
    save_experiment_log(results, log_path)
    assert log_path.exists()

    loaded = load_experiment_log(log_path)
    assert len(loaded) == 2
    assert loaded[0]["experiment"] == 1
    assert loaded[1]["score"] == 0.6


def test_loop_init(workspace_dir: Path):
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    assert loop.workspace.name == "Test Task"
    assert loop.max_iterations == 3
    assert loop.early_stop_n == 2
    assert len(loop.history) == 0
    assert loop.best_score == float("inf")


def test_loop_run_iteration_default_agent(workspace_dir: Path):
    """Default run_agent_on_skill returns empty output, scoring > 0."""
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    result = loop.run_iteration(1)
    assert isinstance(result, ExperimentResult)
    assert result.experiment == 1
    assert result.score >= 0.0
    assert result.kept is True  # first iteration is always kept as baseline


def test_loop_early_stop(workspace_dir: Path):
    """Loop stops after early_stop_n consecutive non-improvements."""
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    results = loop.run()
    # With default agent (always returns same output), should early-stop
    # after 1 kept (baseline) + early_stop_n non-improvements = 3 total
    assert len(results) <= loop.max_iterations


def test_loop_produces_plot(workspace_dir: Path):
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    loop.run()
    assert (workspace_dir / "output" / "progress.png").exists()


def test_loop_produces_experiment_log(workspace_dir: Path):
    loop = AutoResearchLoop(
        workspace_dir=workspace_dir,
        output_dir=workspace_dir / "output",
    )
    loop.run()
    assert (workspace_dir / "output" / "experiment_log.json").exists()
