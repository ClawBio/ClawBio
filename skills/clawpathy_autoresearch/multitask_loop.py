"""Multi-task skill tuner. One central SKILL.md, evaluated across N train traits in parallel each outer iter.

Per outer iter:
  1. Propose a candidate skill (sees central + last round's multi-trait verdicts).
  2. Execute candidate on all train traits in parallel.
  3. Judge each run (opus) in parallel.
  4. Aggregate score = mean of train-trait scores.
  5. If aggregate beats best: promote candidate to central; snapshot.
  6. Every VAL_EVERY iters: execute+judge on val traits (no propose).
  7. Halt on target, max_iters, or val plateau.

Final: evaluate central skill on test traits.
"""
from __future__ import annotations

import json
import shutil
import statistics
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from pathlib import Path
from threading import Lock
from typing import Any

from .dispatcher import ClaudeCLIDispatcher, DispatchRequest
from .executor import execute_skill
from .judge import judge_run


PROMPTS_DIR = Path(__file__).parent / "prompts"
MULTITASK_DIR = Path(__file__).parent / "multitask"
LOG_LOCK = Lock()


def _log(path: Path, entry: dict) -> None:
    with LOG_LOCK:
        with path.open("a") as f:
            f.write(json.dumps(entry) + "\n")


def _load_traits() -> dict:
    return json.loads((MULTITASK_DIR / "traits.json").read_text())["traits"]


def _trait_task_description(name: str, meta: dict) -> str:
    n_note = (
        f"The sumstats file has an N column — pass `--N-col N` to munge_sumstats.py."
        if meta["has_n_column"]
        else f"The sumstats file has NO N column. Total effective sample size from the paper: N={meta['n_effective']}. Pass `--N {meta['n_effective']}` to munge_sumstats.py."
    )
    return (
        f"Reproduce the LDSC SNP-heritability of trait **{name}** ({meta['paper']}). "
        f"Raw sumstats (gzipped TSV): {meta['path']}. Trait type: {meta['trait_type']}. "
        f"Signed-stat column hint: {meta['signed_stat_hint']} (use `--signed-sumstats {meta['signed_stat_hint']},0` in munge). "
        f"{n_note} "
        f"Plausibility range for h²_SNP on this trait: {meta['h2_range'][0]}–{meta['h2_range'][1]} "
        f"(paper's approximate point estimate ≈ {meta['h2_point_approx']}). "
        "Use the LDSC pipeline exactly as the skill prescribes; munge against HapMap3, "
        "run ldsc.py --h2 with 1000G Phase 3 EUR no-MHC weights. Produce output/summary.json, "
        "output/h2.log, output/provenance.json, output/env/. "
        "Do NOT read reference/. Do NOT hardcode trait-specific logic; follow the skill's generic parameterisation."
    )


def _propose_multitask(
    dispatcher,
    central_skill: str,
    last_verdicts: dict[str, dict] | None,
    aggregate_score: float | None,
) -> str:
    """Propose a new central skill. Proposer sees the central skill and the most recent
    multi-trait verdicts so it can write a skill that works on ALL of them."""
    verdicts_block = ""
    if last_verdicts:
        parts = []
        for tname, v in last_verdicts.items():
            vtxt = json.dumps(v, indent=2)[:2000] if isinstance(v, dict) else str(v)[:2000]
            parts.append(f"## {tname}\n{vtxt}")
        verdicts_block = "\n\n# Last round's judge verdicts across training traits\n" + "\n\n".join(parts)
    agg_line = f"\n\nAggregate score last round: {aggregate_score:.3f} (lower = better)" if aggregate_score is not None else ""

    rubric = (MULTITASK_DIR / "rubric.md").read_text()
    prompt = f"""You are a methodology author. Your task is to REWRITE a single SKILL.md that must produce good LDSC h² runs on ANY of several continuous-trait GWAS sumstats — NOT just one trait. The same skill will be executed independently against multiple different raw sumstats files (height, BMI, educational attainment, intelligence, cerebellar volume, anxiety, etc.) and judged separately on each.

# Absolute requirement: the skill must be TRAIT-AGNOSTIC
- The executor is given trait metadata (raw sumstats path, effective N, signed-stat column name, paper citation, plausibility range) in its task description, not in the skill.
- The skill must READ these from the prompt / cwd context, NOT hardcode them.
- Trait-specific hardcoding (paths with 'height_yengo', 'BMI', `if trait == "height"`) → rubric item 1 fails → score floor 0.5.

# Universal rubric (same for every trait)
{rubric}

# Current central SKILL.md
```
{central_skill}
```
{verdicts_block}{agg_line}

Rewrite the SKILL.md so it generalises across the different traits while preserving the reproducibility discipline (env capture, LD preflight, determinism check, demonstrate-don't-assert). Return ONLY the new SKILL.md content. Do not wrap in extra code fences beyond the markdown's own."""
    # Reuse the dispatcher
    out = dispatcher.dispatch(DispatchRequest(
        role="mt_proposer",
        prompt=prompt,
        model="sonnet",
        allowed_tools=[],
    ))
    # Strip outer markdown fences if the model added them
    t = out.strip()
    if t.startswith("```"):
        t = t.split("\n", 1)[1] if "\n" in t else t
        if t.endswith("```"):
            t = t.rsplit("```", 1)[0]
    return t.strip()


def _setup_trait_scratch(scratch_root: Path, trait_name: str) -> Path:
    d = scratch_root / trait_name
    if d.exists():
        shutil.rmtree(d)
    d.mkdir(parents=True)
    (d / "output").mkdir()
    (d / "reference").mkdir()  # empty — judge reads its own reference for ground truth
    return d


def _run_trait(
    *,
    dispatcher,
    executor_tpl: str,
    judge_tpl: str,
    trait_name: str,
    trait_meta: dict,
    skill: str,
    scratch_dir: Path,
    rubric: str,
) -> tuple[float | None, dict | Any, str | Exception]:
    """Execute and judge the skill on one trait. Returns (score, verdict, transcript-or-exc)."""
    task_desc = _trait_task_description(trait_name, trait_meta)
    # Stage ground-truth reference inside the scratch so judge can read it
    gt_src = MULTITASK_DIR / "reference" / f"{trait_name}_ground_truth.json"
    if not gt_src.exists():
        # write ground-truth from meta if we haven't yet
        gt_src.parent.mkdir(parents=True, exist_ok=True)
        gt_src.write_text(json.dumps({
            "trait": trait_name, "paper": trait_meta["paper"],
            "h2_range": trait_meta["h2_range"], "h2_point_approx": trait_meta["h2_point_approx"],
            "notes": "Judge-only. Executor MUST NOT read reference/.",
        }, indent=2))
    (scratch_dir / "reference" / "ground_truth.json").write_text(gt_src.read_text())
    try:
        transcript = execute_skill(
            dispatcher=dispatcher,
            prompt_template=executor_tpl,
            task_description=task_desc,
            skill_content=skill,
            workspace_root=scratch_dir,
            model="sonnet",
        )
    except Exception as exc:
        return None, exc, exc
    try:
        score, verdict = judge_run(
            dispatcher=dispatcher,
            prompt_template=judge_tpl,
            task_description=task_desc,
            rubric=rubric,
            skill_content=skill,
            transcript=transcript,
            output_dir=scratch_dir / "output",
            reference_dir=scratch_dir / "reference",
            model="opus",
        )
        return float(score), verdict, transcript
    except Exception as exc:
        return None, exc, transcript


def _eval_on_traits(
    *,
    dispatcher,
    executor_tpl,
    judge_tpl,
    rubric,
    trait_names: list[str],
    trait_metas: dict,
    skill: str,
    iter_dir: Path,
) -> dict[str, tuple[float | None, Any]]:
    """Run skill on each trait in parallel. Returns {trait: (score_or_None, verdict_or_exc)}."""
    results: dict[str, tuple[float | None, Any]] = {}
    with ThreadPoolExecutor(max_workers=len(trait_names)) as pool:
        futs = {}
        for tn in trait_names:
            scratch = _setup_trait_scratch(iter_dir, tn)
            futs[pool.submit(
                _run_trait,
                dispatcher=dispatcher,
                executor_tpl=executor_tpl,
                judge_tpl=judge_tpl,
                trait_name=tn,
                trait_meta=trait_metas[tn],
                skill=skill,
                scratch_dir=scratch,
                rubric=rubric,
            )] = tn
        for fut in as_completed(futs):
            tn = futs[fut]
            try:
                score, verdict, _tx = fut.result()
                results[tn] = (score, verdict)
            except Exception as exc:
                results[tn] = (None, exc)
    return results


def run_multitask_loop(
    max_outer_iters: int = 10,
    val_every: int = 3,
    val_plateau_patience: int = 2,
    target_aggregate: float = 0.15,
) -> dict:
    dispatcher = ClaudeCLIDispatcher()
    traits = _load_traits()
    train_traits = [n for n, m in traits.items() if m["split"] == "train"]
    val_traits = [n for n, m in traits.items() if m["split"] == "val"]
    test_traits = [n for n, m in traits.items() if m["split"] == "test"]
    print(f"train={train_traits}  val={val_traits}  test={test_traits}", flush=True)

    executor_tpl = (PROMPTS_DIR / "executor.md").read_text()
    judge_tpl = (PROMPTS_DIR / "judge.md").read_text()
    rubric = (MULTITASK_DIR / "rubric.md").read_text()

    central_path = MULTITASK_DIR / "central_skill.md"
    history_path = MULTITASK_DIR / "history.jsonl"
    snapshots_dir = MULTITASK_DIR / "snapshots"
    snapshots_dir.mkdir(exist_ok=True)
    scratch_root = MULTITASK_DIR / "scratch"
    scratch_root.mkdir(exist_ok=True)

    best_aggregate = float("inf")
    last_verdicts: dict[str, dict] | None = None
    last_aggregate: float | None = None
    val_history: list[float] = []
    val_no_improve = 0

    for outer in range(1, max_outer_iters + 1):
        t0 = time.time()
        iter_dir = scratch_root / f"iter-{outer:03d}"
        iter_dir.mkdir(exist_ok=True)

        # 1. Propose
        central_skill = central_path.read_text()
        try:
            candidate_skill = _propose_multitask(
                dispatcher, central_skill, last_verdicts, last_aggregate,
            )
            (iter_dir / "candidate_skill.md").write_text(candidate_skill)
        except Exception as exc:
            _log(history_path, {"outer": outer, "stage": "propose", "error": str(exc)})
            print(f"iter {outer}: proposer failed: {exc}", flush=True)
            continue

        # 2+3. Execute + judge candidate on all train traits (parallel)
        train_results = _eval_on_traits(
            dispatcher=dispatcher,
            executor_tpl=executor_tpl, judge_tpl=judge_tpl, rubric=rubric,
            trait_names=train_traits, trait_metas=traits,
            skill=candidate_skill, iter_dir=iter_dir,
        )

        scores = [s for (s, _) in train_results.values() if isinstance(s, (int, float))]
        if not scores:
            _log(history_path, {"outer": outer, "stage": "eval", "error": "no scores", "results": {
                k: ("err" if isinstance(v[1], Exception) else v[0]) for k, v in train_results.items()
            }})
            print(f"iter {outer}: no train scores produced", flush=True)
            continue

        aggregate = statistics.mean(scores)
        kept = aggregate < best_aggregate
        verdicts_for_log = {
            k: (v[1] if isinstance(v[1], dict) else f"ERR:{type(v[1]).__name__}")
            for k, v in train_results.items()
        }
        _log(history_path, {
            "outer": outer,
            "phase": "train",
            "wall_sec": round(time.time() - t0, 1),
            "train_scores": {k: (v[0] if isinstance(v[0], (int, float)) else None)
                             for k, v in train_results.items()},
            "aggregate": aggregate,
            "kept": kept,
            "best_aggregate_before": None if best_aggregate == float("inf") else best_aggregate,
        })

        if kept:
            central_path.write_text(candidate_skill)
            (snapshots_dir / f"iter-{outer:03d}-agg{aggregate:.3f}.md").write_text(candidate_skill)
            best_aggregate = aggregate
            print(f"iter {outer}: KEPT aggregate={aggregate:.3f} (new central)", flush=True)
        else:
            print(f"iter {outer}: rejected aggregate={aggregate:.3f} (best={best_aggregate:.3f})", flush=True)

        last_verdicts = {k: (v[1] if isinstance(v[1], dict) else None) for k, v in train_results.items()}
        last_aggregate = aggregate

        # 4. Periodic val eval
        if outer % val_every == 0:
            val_results = _eval_on_traits(
                dispatcher=dispatcher,
                executor_tpl=executor_tpl, judge_tpl=judge_tpl, rubric=rubric,
                trait_names=val_traits, trait_metas=traits,
                skill=central_path.read_text(), iter_dir=iter_dir / "_val",
            )
            vscores = [s for (s, _) in val_results.values() if isinstance(s, (int, float))]
            val_agg = statistics.mean(vscores) if vscores else None
            _log(history_path, {
                "outer": outer, "phase": "val",
                "val_scores": {k: (v[0] if isinstance(v[0], (int, float)) else None)
                               for k, v in val_results.items()},
                "val_aggregate": val_agg,
            })
            if val_agg is not None:
                if val_history and val_agg >= min(val_history):
                    val_no_improve += 1
                else:
                    val_no_improve = 0
                val_history.append(val_agg)
                print(f"  val_aggregate={val_agg:.3f}  plateau={val_no_improve}/{val_plateau_patience}", flush=True)
                if val_no_improve >= val_plateau_patience:
                    print(f"  val plateau — halting training", flush=True)
                    break

        if best_aggregate <= target_aggregate:
            print(f"  target aggregate {target_aggregate} reached — halting", flush=True)
            break

    # 5. Final test eval
    print("\n=== FINAL TEST EVAL ===", flush=True)
    test_iter_dir = scratch_root / "final_test"
    test_iter_dir.mkdir(exist_ok=True)
    test_results = _eval_on_traits(
        dispatcher=dispatcher,
        executor_tpl=executor_tpl, judge_tpl=judge_tpl, rubric=rubric,
        trait_names=test_traits, trait_metas=traits,
        skill=central_path.read_text(), iter_dir=test_iter_dir,
    )
    tscores = [s for (s, _) in test_results.values() if isinstance(s, (int, float))]
    test_agg = statistics.mean(tscores) if tscores else None
    _log(history_path, {
        "phase": "test_final",
        "test_scores": {k: (v[0] if isinstance(v[0], (int, float)) else None)
                        for k, v in test_results.items()},
        "test_aggregate": test_agg,
    })
    return {
        "best_train_aggregate": None if best_aggregate == float("inf") else best_aggregate,
        "test_aggregate": test_agg,
        "val_history": val_history,
        "central_skill": str(central_path),
    }


if __name__ == "__main__":
    import argparse
    p = argparse.ArgumentParser()
    p.add_argument("--max-outer-iters", type=int, default=10)
    p.add_argument("--val-every", type=int, default=3)
    p.add_argument("--target-aggregate", type=float, default=0.15)
    args = p.parse_args()
    result = run_multitask_loop(
        max_outer_iters=args.max_outer_iters,
        val_every=args.val_every,
        target_aggregate=args.target_aggregate,
    )
    print(json.dumps(result, indent=2))
