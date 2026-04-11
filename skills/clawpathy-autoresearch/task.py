"""Task definition loader for clawpathy-autoresearch."""
from __future__ import annotations

from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import yaml


@dataclass
class PaperDefinition:
    """A single paper in the benchmark set."""

    id: str
    title: str
    pmid: str
    ground_truth: dict[str, Any] | None = None


@dataclass
class TaskDefinition:
    """A complete task definition loaded from task.yaml."""

    name: str
    description: str
    metric: str
    direction: str  # "higher_is_better" or "lower_is_better"
    scale: list[float]
    skills_dir: str
    papers: list[PaperDefinition] = field(default_factory=list)

    @property
    def paper_count(self) -> int:
        return len(self.papers)


def load_task(task_yaml: Path) -> TaskDefinition:
    """Load a task definition from a YAML file.

    Resolves ground truth file paths relative to the task.yaml directory.
    """
    task_yaml = Path(task_yaml)
    if not task_yaml.exists():
        raise FileNotFoundError(f"Task file not found: {task_yaml}")

    raw = yaml.safe_load(task_yaml.read_text())
    task_dir = task_yaml.parent

    papers = []
    for p in raw.get("papers", []):
        gt_path = task_dir / p["ground_truth"]
        gt_data = None
        if gt_path.exists():
            gt_data = yaml.safe_load(gt_path.read_text())

        papers.append(
            PaperDefinition(
                id=p["id"],
                title=p["title"],
                pmid=p["pmid"],
                ground_truth=gt_data,
            )
        )

    return TaskDefinition(
        name=raw["name"],
        description=raw["description"],
        metric=raw["metric"],
        direction=raw["direction"],
        scale=raw["scale"],
        skills_dir=raw["skills_dir"],
        papers=papers,
    )
