"""K-way parallel tuning loop. Per iter: propose K → execute K → judge K → keep best."""
from __future__ import annotations

import json
import math
import shutil
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from typing import Any

from .dispatcher import Dispatcher
from .executor import execute_skill
from .judge import judge_run
from .proposer import propose_skill
from .workspace import load_workspace, recover_state


PROMPTS_DIR = Path(__file__).parent / "prompts"

# Subdirs symlinked into each candidate (shared read-only from executor's POV).
# data/ is included because in this task's current state, data/processed/ is
# pre-built and the executor only reads it; new artifacts go to output/.
_SHARED_SUBDIRS = ["reference", "pipeline", "scripts", "tools", "data"]

DIVERSITY_FOCUSES = [
    "Prioritise reproducibility: environment lockfile, SHA256 provenance, determinism via multi-run averaging when seed unavailable, per-locus reconstructability.",
    "Prioritise numerical fidelity: LD/z-file order verification, harmonization diagnostics, credible-set plausibility gates (no k-saturation, no 100% flip-rate dismissal), lead-SNP sanity checks.",
    "Prioritise simplification: trim the skill to the checks that have caught real bugs; remove redundant gates; make every instruction testable and imperative.",
]

_LOG_LOCK = Lock()


def _read_prompt(name: str) -> str:
    return (PROMPTS_DIR / name).read_text()


def _log(path: Path, entry: dict[str, Any]) -> None:
    with _LOG_LOCK:
        with path.open("a") as f:
            f.write(json.dumps(entry) + "\n")


def _setup_candidate(ws_root: Path, iter_dir: Path, k: int, champion_skill: str) -> Path:
    cand = iter_dir / f"k{k}"
    if cand.exists():
        shutil.rmtree(cand)
    cand.mkdir(parents=True)
    (cand / "skill").mkdir()
    (cand / "output").mkdir()
    (cand / "executor_runs").mkdir()
    (cand / "snapshots").mkdir()
    for f in ("task.json", "rubric.md"):
        shutil.copy2(ws_root / f, cand / f)
    for sub in _SHARED_SUBDIRS:
        src = (ws_root / sub).resolve()
        if src.exists():
            (cand / sub).symlink_to(src)
    (cand / "skill" / "SKILL.md").write_text(champion_skill)
    return cand


def _propose_with_focus(
    *, dispatcher, prompt_tpl, task_description, rubric, champion_skill,
    last_score, last_judgement, recent_history, model, focus,
) -> str:
    augmented = prompt_tpl + f"\n\n# Candidate focus for this proposal\n{focus}\n"
    return propose_skill(
        dispatcher=dispatcher,
        prompt_template=augmented,
        task_description=task_description,
        rubric=rubric,
        current_skill=champion_skill,
        last_score=last_score,
        last_judgement=last_judgement,
        recent_history=recent_history,
        model=model,
    )


def run_parallel_loop(
    workspace_root: Path,
    dispatcher: Dispatcher,
    K: int = 3,
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

    candidates_root = ws.root / "candidates"
    candidates_root.mkdir(exist_ok=True)

    for i in range(1, ws.max_iterations + 1):
        champion_skill = ws.skill_md.read_text()
        iter_dir = candidates_root / f"iter-{i:03d}"
        iter_dir.mkdir(exist_ok=True)

        cand_paths = [_setup_candidate(ws.root, iter_dir, k, champion_skill) for k in range(K)]
        focuses = [DIVERSITY_FOCUSES[k % len(DIVERSITY_FOCUSES)] for k in range(K)]

        # Phase 1: propose K in parallel
        proposals: list[Any] = [None] * K
        with ThreadPoolExecutor(max_workers=K) as pool:
            futs = {pool.submit(
                _propose_with_focus,
                dispatcher=dispatcher, prompt_tpl=proposer_tpl,
                task_description=ws.description, rubric=rubric,
                champion_skill=champion_skill,
                last_score=last_score, last_judgement=last_judgement,
                recent_history=recent[-3:], model=proposer_model,
                focus=focuses[k],
            ): k for k in range(K)}
            for fut in as_completed(futs):
                k = futs[fut]
                try:
                    proposals[k] = fut.result()
                    (cand_paths[k] / "skill" / "SKILL.md").write_text(proposals[k])
                except Exception as exc:
                    proposals[k] = exc
                    _log(ws.history_path, {"iter": i, "candidate": k, "stage": "propose", "error": str(exc)})

        # Phase 2: execute the successful proposals in parallel
        transcripts: list[Any] = [None] * K
        with ThreadPoolExecutor(max_workers=K) as pool:
            futs = {}
            for k in range(K):
                if isinstance(proposals[k], Exception):
                    continue
                futs[pool.submit(
                    execute_skill,
                    dispatcher=dispatcher, prompt_template=executor_tpl,
                    task_description=ws.description,
                    skill_content=proposals[k],
                    workspace_root=cand_paths[k],
                    model=executor_model,
                )] = k
            for fut in as_completed(futs):
                k = futs[fut]
                try:
                    transcripts[k] = fut.result()
                    (cand_paths[k] / "executor_runs" / f"iter-{i:03d}.log").write_text(transcripts[k])
                except Exception as exc:
                    transcripts[k] = exc
                    _log(ws.history_path, {"iter": i, "candidate": k, "stage": "execute", "error": str(exc)})

        # Phase 3: judge in parallel
        verdicts: list[Any] = [None] * K
        scores: list[float | None] = [None] * K
        with ThreadPoolExecutor(max_workers=K) as pool:
            futs = {}
            for k in range(K):
                if isinstance(proposals[k], Exception) or isinstance(transcripts[k], Exception) or transcripts[k] is None:
                    continue
                futs[pool.submit(
                    judge_run,
                    dispatcher=dispatcher, prompt_template=judge_tpl,
                    task_description=ws.description, rubric=rubric,
                    skill_content=proposals[k],
                    transcript=transcripts[k],
                    output_dir=cand_paths[k] / "output",
                    reference_dir=ws.reference_dir,
                    model=judge_model,
                )] = k
            for fut in as_completed(futs):
                k = futs[fut]
                try:
                    score, verdict = fut.result()
                    scores[k] = score
                    verdicts[k] = verdict
                except Exception as exc:
                    verdicts[k] = exc
                    _log(ws.history_path, {"iter": i, "candidate": k, "stage": "judge", "error": str(exc)})

        ranked = sorted([(s, k) for k, s in enumerate(scores) if s is not None])
        if not ranked:
            no_improve += 1
            for cp in cand_paths:
                shutil.rmtree(cp, ignore_errors=True)
            shutil.rmtree(iter_dir, ignore_errors=True)
            if no_improve >= ws.early_stop_n:
                break
            continue

        best_of_iter_score, best_k = ranked[0]
        iter_kept = best_of_iter_score < best_score

        for k in range(K):
            if scores[k] is None:
                continue
            _log(ws.history_path, {
                "iter": i, "candidate": k, "focus": focuses[k],
                "score": scores[k], "kept": (k == best_k) and iter_kept,
                "verdict": verdicts[k],
            })

        if iter_kept:
            shutil.copy2(cand_paths[best_k] / "skill" / "SKILL.md", ws.skill_md)
            (ws.snapshots_dir / f"iter-{i:03d}.md").write_text(
                (cand_paths[best_k] / "skill" / "SKILL.md").read_text()
            )
            best_score = best_of_iter_score
            no_improve = 0
        else:
            no_improve += 1

        last_score = best_of_iter_score
        last_judgement = verdicts[best_k] if isinstance(verdicts[best_k], dict) else None
        recent.append({"iter": i, "score": best_of_iter_score, "kept": iter_kept})

        for cp in cand_paths:
            shutil.rmtree(cp, ignore_errors=True)
        shutil.rmtree(iter_dir, ignore_errors=True)

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
