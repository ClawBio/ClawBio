"""Scorer for the hello-world integration task."""
from __future__ import annotations

import json
from pathlib import Path


def score(skill_output: dict, task_dir: Path) -> tuple[float, dict]:
    gt = json.loads((Path(task_dir) / "ground_truth.json").read_text())
    got = skill_output.get("text")
    if got == gt["target"]:
        return 0.0, {"text_match": 0.0}
    return 1.0, {"text_match": 1.0, "got": got, "want": gt["target"]}
