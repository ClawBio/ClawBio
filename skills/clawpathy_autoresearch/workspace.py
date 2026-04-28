"""Workspace layout for an eval-driven skill-tuning task."""
from __future__ import annotations

import json
from dataclasses import dataclass
from pathlib import Path


REQUIRED_FILES = ["task.json", "rubric.md", "skill/SKILL.md"]


@dataclass
class Workspace:
    root: Path
    task: dict

    @property
    def name(self) -> str:
        return self.task.get("name", self.root.name)

    @property
    def description(self) -> str:
        return self.task.get("description", "")

    @property
    def max_iterations(self) -> int:
        return int(self.task.get("max_iterations", 30))

    @property
    def early_stop_n(self) -> int:
        return int(self.task.get("early_stop_n", 6))

    @property
    def target_score(self) -> float | None:
        v = self.task.get("target_score")
        return float(v) if v is not None else None

    @property
    def skill_md(self) -> Path:
        return self.root / "skill" / "SKILL.md"

    @property
    def rubric_md(self) -> Path:
        return self.root / "rubric.md"

    @property
    def reference_dir(self) -> Path:
        return self.root / "reference"

    @property
    def output_dir(self) -> Path:
        return self.root / "output"

    @property
    def runs_dir(self) -> Path:
        return self.root / "executor_runs"

    @property
    def snapshots_dir(self) -> Path:
        return self.root / "snapshots"

    @property
    def history_path(self) -> Path:
        return self.root / "history.jsonl"


def validate_workspace(root: Path) -> list[str]:
    errors = []
    for rel in REQUIRED_FILES:
        if not (root / rel).is_file():
            errors.append(f"missing: {rel}")
    try:
        json.loads((root / "task.json").read_text())
    except Exception as e:
        errors.append(f"task.json not valid JSON: {e}")
    return errors


def load_workspace(root: Path) -> Workspace:
    errors = validate_workspace(root)
    if errors:
        raise ValueError(f"invalid workspace {root}: {errors}")
    task = json.loads((root / "task.json").read_text())
    ws = Workspace(root=root, task=task)
    for d in (ws.output_dir, ws.runs_dir, ws.snapshots_dir):
        d.mkdir(parents=True, exist_ok=True)
    return ws


def recover_state(history_path: Path) -> dict:
    """Rebuild best_score / last_score / last_judgement / recent from history.jsonl.

    Scans every scored entry (serial or parallel format). best_score is the
    minimum score across kept=True entries. last_score/last_judgement come
    from the most recent scored entry. recent is the last 3 scored entries.
    """
    import math
    state: dict = {
        "best_score": math.inf,
        "last_score": None,
        "last_judgement": None,
        "recent": [],
    }
    if not history_path.exists():
        return state
    scored = []
    for line in history_path.read_text().splitlines():
        if not line.strip():
            continue
        try:
            d = json.loads(line)
        except Exception:
            continue
        if not isinstance(d.get("score"), (int, float)):
            continue
        scored.append(d)
        if d.get("kept") and d["score"] < state["best_score"]:
            state["best_score"] = float(d["score"])
    if scored:
        last = scored[-1]
        state["last_score"] = float(last["score"])
        if isinstance(last.get("verdict"), dict):
            state["last_judgement"] = last["verdict"]
        state["recent"] = [
            {"iter": d.get("iter"), "score": d["score"], "kept": d.get("kept", False)}
            for d in scored[-3:]
        ]
    return state
