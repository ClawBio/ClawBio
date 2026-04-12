"""Iteration loop engine for clawpathy-autoresearch."""
from __future__ import annotations

import importlib.util
import json
import math
import shutil
from pathlib import Path
from typing import Any

from .dispatcher import Dispatcher
from .executor import execute_skill
from .proposer import propose_skill
from .sandbox import load_sandbox
from .workspace import load_workspace


def _load_scorer(scorer_path: Path):
    spec = importlib.util.spec_from_file_location("_autoresearch_scorer", scorer_path)
    if spec is None or spec.loader is None:
        raise ImportError(f"Could not load scorer: {scorer_path}")
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    if not hasattr(mod, "score"):
        raise AttributeError(f"Scorer at {scorer_path} has no score() function")
    return mod.score


def _snapshot(skill_md: Path, snap_dir: Path, iter_num: int) -> None:
    snap_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(skill_md, snap_dir / f"iter-{iter_num:03d}.md")


def _revert(skill_md: Path, snap_dir: Path, iter_num: int) -> None:
    snap = snap_dir / f"iter-{iter_num:03d}.md"
    shutil.copy2(snap, skill_md)


def run_loop(workspace_dir: Path, dispatcher: Dispatcher) -> dict[str, Any]:
    ws = load_workspace(Path(workspace_dir))
    sandbox = load_sandbox(ws.sandbox_path)
    scorer = _load_scorer(ws.scorer_path)

    skill_md = ws.skill_dir / "SKILL.md"
    history_path = ws.workspace_dir / "history.jsonl"
    snap_dir = ws.workspace_dir / "snapshots"
    history_path.write_text("")

    task_prompt = ws.description or ws.name
    best_score = math.inf
    last_score: float | None = None
    last_breakdown: dict[str, Any] | None = None
    recent_history: list[dict[str, Any]] = []
    consecutive_non_improvements = 0
    iterations_run = 0

    for i in range(1, ws.max_iterations + 1):
        iterations_run = i
        _snapshot(skill_md, snap_dir, i)

        proposal = propose_skill(
            dispatcher=dispatcher,
            task_prompt=task_prompt,
            current_skill=skill_md.read_text(),
            last_score=last_score,
            last_breakdown=last_breakdown,
            recent_history=recent_history[-3:],
        )
        skill_md.write_text(proposal)

        try:
            output = execute_skill(
                dispatcher=dispatcher,
                skill_content=proposal,
                sandbox=sandbox,
                data_dir_abs=str(ws.data_dir.resolve()),
            )
            total, breakdown = scorer(output, ws.workspace_dir)
        except Exception as exc:
            total = math.inf
            breakdown = {"error": str(exc)}

        kept = total < best_score
        if kept:
            best_score = total
            consecutive_non_improvements = 0
        else:
            _revert(skill_md, snap_dir, i)
            consecutive_non_improvements += 1

        entry = {
            "iter": i,
            "score": total if total != math.inf else None,
            "kept": kept,
            "breakdown": breakdown,
        }
        recent_history.append(entry)
        last_score = total if total != math.inf else None
        last_breakdown = breakdown
        with history_path.open("a") as f:
            f.write(json.dumps(entry) + "\n")

        if ws.target_score is not None and best_score < ws.target_score:
            break
        if consecutive_non_improvements >= max(1, ws.early_stop_n - 1):
            break

    return {
        "best_score": best_score if best_score != math.inf else None,
        "iterations_run": iterations_run,
        "history_path": str(history_path),
    }
