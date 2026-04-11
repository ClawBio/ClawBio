"""Workspace loader for clawpathy-autoresearch.

A workspace is a self-contained directory with everything needed to run
the autoresearch optimisation loop: task config, ground truth, scorer,
and the SKILL.md being optimised.
"""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path
from typing import Any


@dataclass
class Workspace:
    """A loaded autoresearch workspace."""

    name: str
    description: str
    max_iterations: int
    early_stop_n: int
    ground_truth: dict[str, Any]
    workspace_dir: Path
    skill_dir: Path
    sources_dir: Path
    scorer_path: Path


def validate_workspace(workspace_dir: Path) -> list[str]:
    """Check a workspace directory for required files. Returns list of errors."""
    workspace_dir = Path(workspace_dir)
    errors = []

    if not workspace_dir.exists():
        return [f"Workspace directory does not exist: {workspace_dir}"]

    required = [
        ("task.json", workspace_dir / "task.json"),
        ("ground_truth.json", workspace_dir / "ground_truth.json"),
        ("scorer.py", workspace_dir / "scorer.py"),
        ("skill/SKILL.md", workspace_dir / "skill" / "SKILL.md"),
    ]
    for label, path in required:
        if not path.exists():
            errors.append(f"Missing required file: {label}")

    return errors


def load_workspace(workspace_dir: Path) -> Workspace:
    """Load a workspace from a directory.

    Raises FileNotFoundError if the directory is missing or invalid.
    """
    workspace_dir = Path(workspace_dir)

    if not workspace_dir.exists():
        raise FileNotFoundError(f"Workspace not found: {workspace_dir}")

    errors = validate_workspace(workspace_dir)
    if errors:
        raise FileNotFoundError(
            f"Invalid workspace: {'; '.join(errors)}"
        )

    task_data = json.loads((workspace_dir / "task.json").read_text())
    ground_truth = json.loads((workspace_dir / "ground_truth.json").read_text())

    return Workspace(
        name=task_data["name"],
        description=task_data.get("description", ""),
        max_iterations=task_data.get("max_iterations", 80),
        early_stop_n=task_data.get("early_stop_n", 5),
        ground_truth=ground_truth,
        workspace_dir=workspace_dir,
        skill_dir=workspace_dir / "skill",
        sources_dir=workspace_dir / "sources",
        scorer_path=workspace_dir / "scorer.py",
    )
