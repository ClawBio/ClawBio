"""Tests for task loader."""
from __future__ import annotations

import tempfile
from pathlib import Path

import pytest
import yaml

from skills.clawpathy_autoresearch.task import TaskDefinition, load_task


@pytest.fixture
def sample_task_dir(tmp_path: Path) -> Path:
    """Create a minimal task directory with task.yaml and one ground truth."""
    gt_dir = tmp_path / "ground_truth"
    gt_dir.mkdir()

    task_config = {
        "name": "Test Task",
        "description": "Reproduce findings from test papers",
        "metric": "reproduction_score",
        "direction": "higher_is_better",
        "scale": [0, 10],
        "skills_dir": "/tmp/fake_skills",
        "papers": [
            {
                "id": "paper_001",
                "title": "Test GWAS of Test Disease",
                "pmid": "12345678",
                "ground_truth": "ground_truth/paper_001.yaml",
            }
        ],
    }
    (tmp_path / "task.yaml").write_text(yaml.dump(task_config))

    gt = {
        "lead_variants": [
            {
                "rsid": "rs123456",
                "gene": "BRCA1",
                "p_value_order": -20,
                "effect_direction": "risk",
                "or_range": [1.1, 1.3],
            }
        ],
        "qualitative_findings": ["Found immune enrichment"],
        "total_loci": 50,
        "ancestry": "European",
    }
    (gt_dir / "paper_001.yaml").write_text(yaml.dump(gt))

    return tmp_path


def test_load_task_returns_task_definition(sample_task_dir: Path):
    task = load_task(sample_task_dir / "task.yaml")
    assert isinstance(task, TaskDefinition)
    assert task.name == "Test Task"
    assert task.direction == "higher_is_better"
    assert len(task.papers) == 1


def test_load_task_resolves_ground_truth(sample_task_dir: Path):
    task = load_task(sample_task_dir / "task.yaml")
    paper = task.papers[0]
    assert paper.ground_truth is not None
    assert len(paper.ground_truth["lead_variants"]) == 1
    assert paper.ground_truth["lead_variants"][0]["rsid"] == "rs123456"


def test_load_task_missing_file_raises():
    with pytest.raises(FileNotFoundError):
        load_task(Path("/nonexistent/task.yaml"))


def test_task_definition_paper_count(sample_task_dir: Path):
    task = load_task(sample_task_dir / "task.yaml")
    assert task.paper_count == 1
