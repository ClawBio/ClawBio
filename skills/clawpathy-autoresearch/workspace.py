"""Workspace loader for clawpathy-autoresearch.

A workspace is a self-contained directory with task config, ground truth,
scorer, sandbox spec, and the SKILL.md being optimised.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class Workspace:
    name: str
    description: str
    max_iterations: int
    early_stop_n: int
    target_score: float | None
    ground_truth: dict[str, Any]
    workspace_dir: Path
    skill_dir: Path
    scorer_path: Path
    sandbox_path: Path
    data_dir: Path


def validate_workspace(workspace_dir: Path) -> list[str]:
    workspace_dir = Path(workspace_dir)
    errors: list[str] = []
    if not workspace_dir.exists():
        return [f"Workspace directory does not exist: {workspace_dir}"]
    required = [
        ("task.json", workspace_dir / "task.json"),
        ("ground_truth.json", workspace_dir / "ground_truth.json"),
        ("scorer.py", workspace_dir / "scorer.py"),
        ("sandbox.yaml", workspace_dir / "sandbox.yaml"),
        ("skill/SKILL.md", workspace_dir / "skill" / "SKILL.md"),
    ]
    for label, path in required:
        if not path.exists():
            errors.append(f"Missing required file: {label}")
    return errors


def load_workspace(workspace_dir: Path) -> Workspace:
    workspace_dir = Path(workspace_dir)
    if not workspace_dir.exists():
        raise FileNotFoundError(f"Workspace not found: {workspace_dir}")
    errors = validate_workspace(workspace_dir)
    if errors:
        raise FileNotFoundError(f"Invalid workspace: {'; '.join(errors)}")
    task_data = json.loads((workspace_dir / "task.json").read_text())
    ground_truth = json.loads((workspace_dir / "ground_truth.json").read_text())
    return Workspace(
        name=task_data["name"],
        description=task_data.get("description", ""),
        max_iterations=task_data.get("max_iterations", 20),
        early_stop_n=task_data.get("early_stop_n", 5),
        target_score=task_data.get("target_score"),
        ground_truth=ground_truth,
        workspace_dir=workspace_dir,
        skill_dir=workspace_dir / "skill",
        scorer_path=workspace_dir / "scorer.py",
        sandbox_path=workspace_dir / "sandbox.yaml",
        data_dir=workspace_dir / "data",
    )
