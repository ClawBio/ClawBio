"""Tests for the main autoresearch loop."""
from __future__ import annotations

import json
from pathlib import Path
from unittest.mock import patch

import pytest
import yaml

from skills.clawpathy_autoresearch.autoresearch import (
    AutoResearchLoop,
    ExperimentResult,
    load_experiment_log,
    save_experiment_log,
)


@pytest.fixture
def task_dir(tmp_path: Path) -> Path:
    """Create a complete task directory."""
    gt_dir = tmp_path / "tasks" / "test_task" / "ground_truth"
    gt_dir.mkdir(parents=True)

    task_config = {
        "name": "Test Reproduction",
        "description": "Test task",
        "metric": "reproduction_score",
        "direction": "higher_is_better",
        "scale": [0, 10],
        "skills_dir": str(tmp_path / "skills"),
        "papers": [
            {
                "id": "paper_001",
                "title": "Test GWAS",
                "pmid": "12345678",
                "ground_truth": "ground_truth/paper_001.yaml",
            }
        ],
    }
    task_yaml = tmp_path / "tasks" / "test_task" / "task.yaml"
    task_yaml.write_text(yaml.dump(task_config))

    gt = {
        "lead_variants": [
            {
                "rsid": "rs123456",
                "gene": "TEST1",
                "p_value_order": -20,
                "effect_direction": "risk",
                "or_range": [1.1, 1.3],
            }
        ],
        "qualitative_findings": ["Found test enrichment"],
        "total_loci": 10,
        "ancestry": "European",
    }
    (gt_dir / "paper_001.yaml").write_text(yaml.dump(gt))

    # Create a skills dir with one skill
    skills_dir = tmp_path / "skills" / "test-skill"
    skills_dir.mkdir(parents=True)
    (skills_dir / "SKILL.md").write_text(
        "---\nname: test-skill\n---\n\n# Test Skill\n\n## Workflow\n\n1. Do thing\n"
    )

    return tmp_path


def test_experiment_result_serialisation():
    result = ExperimentResult(
        experiment=1,
        score=5.5,
        kept=True,
        label="baseline",
        skill_changes=["test-skill"],
        skills_created=[],
        per_paper_scores={"paper_001": 5.5},
        automated_score=6.0,
        llm_score=4.5,
    )
    d = result.to_dict()
    assert d["experiment"] == 1
    assert d["score"] == 5.5
    assert d["kept"] is True


def test_save_and_load_experiment_log(tmp_path: Path):
    results = [
        ExperimentResult(
            experiment=1, score=3.0, kept=True, label="baseline",
            skill_changes=[], skills_created=[],
            per_paper_scores={}, automated_score=3.0, llm_score=3.0,
        ),
        ExperimentResult(
            experiment=2, score=4.0, kept=True, label="improved",
            skill_changes=["x"], skills_created=[],
            per_paper_scores={}, automated_score=4.0, llm_score=4.0,
        ),
    ]
    log_path = tmp_path / "experiment_log.json"
    save_experiment_log(results, log_path)
    assert log_path.exists()

    loaded = load_experiment_log(log_path)
    assert len(loaded) == 2
    assert loaded[0]["experiment"] == 1
    assert loaded[1]["score"] == 4.0


def test_autoresearch_loop_init(task_dir: Path):
    loop = AutoResearchLoop(
        task_path=task_dir / "tasks" / "test_task" / "task.yaml",
        output_dir=task_dir / "results",
        max_iterations=3,
    )
    assert loop.task.name == "Test Reproduction"
    assert loop.max_iterations == 3
    assert len(loop.history) == 0
