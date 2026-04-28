"""Eval-driven skill-tuning loop. Propose → execute → judge → keep/revert."""
from __future__ import annotations

import json
import math
import shutil
from pathlib import Path
from typing import Any

from .dispatcher import Dispatcher
from .executor import execute_skill
from .judge import judge_run
from .proposer import propose_skill
from .workspace import load_workspace, recover_state


PROMPTS_DIR = Path(__file__).parent / "prompts"


def _read_prompt(name: str) -> str:
    return (PROMPTS_DIR / name).read_text()


def _snapshot(skill_md: Path, snap_dir: Path, i: int) -> None:
    snap_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(skill_md, snap_dir / f"iter-{i:03d}.md")


def _revert(skill_md: Path, snap_dir: Path, i: int) -> None:
    shutil.copy2(snap_dir / f"iter-{i:03d}.md", skill_md)


def run_loop(
    workspace_root: Path,
    dispatcher: Dispatcher,
    proposer_model: str = "sonnet",
    executor_model: str = "sonnet",
    judge_model: str = "opus",
    resume: bool = True,
) -> dict[str, Any]:
    ws = load_workspace(Path(workspace_root))
    proposer_tpl = _read_prompt("proposer.md")
    executor_tpl = _read_prompt("executor.md")
    judge_tpl = _read_prompt("judge.md")
    rubric = ws.rubric_md.read_text()

    if not resume:
        ws.history_path.write_text("")

    st = recover_state(ws.history_path) if resume else {
        "best_score": math.inf, "last_score": None, "last_judgement": None, "recent": []
    }
    best_score = st["best_score"]
    last_score: float | None = st["last_score"]
    last_judgement: dict[str, Any] | None = st["last_judgement"]
    recent: list[dict[str, Any]] = list(st["recent"])
    no_improve = 0

    for i in range(1, ws.max_iterations + 1):
        _snapshot(ws.skill_md, ws.snapshots_dir, i)

        # 1. Propose
        try:
            proposal = propose_skill(
                dispatcher=dispatcher,
                prompt_template=proposer_tpl,
                task_description=ws.description,
                rubric=rubric,
                current_skill=ws.skill_md.read_text(),
                last_score=last_score,
                last_judgement=last_judgement,
                recent_history=recent[-3:],
                model=proposer_model,
            )
            ws.skill_md.write_text(proposal)
        except Exception as exc:
            _log(ws.history_path, {"iter": i, "stage": "propose", "error": str(exc)})
            _revert(ws.skill_md, ws.snapshots_dir, i)
            no_improve += 1
            if no_improve >= ws.early_stop_n:
                break
            continue

        # 2. Execute
        run_log = ws.runs_dir / f"iter-{i:03d}.log"
        try:
            # Clear outputs so the judge sees only this iteration's artifacts.
            if ws.output_dir.exists():
                shutil.rmtree(ws.output_dir)
            ws.output_dir.mkdir(parents=True)
            transcript = execute_skill(
                dispatcher=dispatcher,
                prompt_template=executor_tpl,
                task_description=ws.description,
                skill_content=proposal,
                workspace_root=ws.root,
                model=executor_model,
            )
            run_log.write_text(transcript)
        except Exception as exc:
            run_log.write_text(f"EXECUTOR ERROR: {exc}")
            _log(ws.history_path, {"iter": i, "stage": "execute", "error": str(exc)})
            _revert(ws.skill_md, ws.snapshots_dir, i)
            no_improve += 1
            if no_improve >= ws.early_stop_n:
                break
            continue

        # 3. Judge
        try:
            score, verdict = judge_run(
                dispatcher=dispatcher,
                prompt_template=judge_tpl,
                task_description=ws.description,
                rubric=rubric,
                skill_content=proposal,
                transcript=transcript,
                output_dir=ws.output_dir,
                reference_dir=ws.reference_dir,
                model=judge_model,
            )
        except Exception as exc:
            _log(ws.history_path, {"iter": i, "stage": "judge", "error": str(exc)})
            _revert(ws.skill_md, ws.snapshots_dir, i)
            no_improve += 1
            if no_improve >= ws.early_stop_n:
                break
            continue

        kept = score < best_score
        if kept:
            best_score = score
            no_improve = 0
        else:
            _revert(ws.skill_md, ws.snapshots_dir, i)
            no_improve += 1

        entry = {"iter": i, "score": score, "kept": kept, "verdict": verdict}
        recent.append(entry)
        _log(ws.history_path, entry)
        last_score, last_judgement = score, verdict

        if ws.target_score is not None and best_score <= ws.target_score:
            break
        if no_improve >= ws.early_stop_n:
            break

    return {
        "best_score": None if best_score == math.inf else best_score,
        "iterations_run": len(recent),
        "history_path": str(ws.history_path),
        "final_skill": str(ws.skill_md),
    }


def _log(path: Path, entry: dict[str, Any]) -> None:
    with path.open("a") as f:
        f.write(json.dumps(entry) + "\n")
