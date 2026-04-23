---
name: clawpathy-autoresearch
description: >-
  Eval-driven skill tuning. Given a task and an LLM-judge rubric, iteratively
  rewrites a SKILL.md until a downstream executor agent performs well against
  the judge. Low-code: all evaluation is LLM-as-judge, not deterministic Python.
version: 1.0.0
author: Jay Moore
license: MIT
tags: [meta, autoresearch, skill-tuning, llm-judge, eval-driven]

inputs:
  - name: paper_query_or_task
    type: string
    description: Paper title/URL/PMID/DOI, or a freeform task description
    required: true

outputs:
  - name: workspace/
    type: directory
    description: Tuned skill/SKILL.md plus history.jsonl, snapshots, executor_runs

metadata:
  openclaw:
    requires:
      bins: [python3, claude]
    always: false
    emoji: "🔁"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    trigger_keywords:
      - auto research
      - autoresearch
      - tune a skill
      - skill tuning
      - improve a skill
      - eval-driven
      - clawpathy
      - replicate paper
      - reproduce paper
---

# clawpathy-autoresearch

Eval-driven skill development. The system iteratively rewrites a `SKILL.md`
so a downstream executor agent performs better at a task class, as judged
by an LLM against a paper/task-specific rubric.

## Core idea

```
  propose (sonnet)  →  execute (sonnet, shell)  →  judge (opus, rubric)
       ↑                                                       │
       └──────── feedback: verdict + recommended edits ────────┘
```

- **Proposer** rewrites SKILL.md based on the last judge verdict.
- **Executor** runs the new SKILL.md end-to-end inside a workspace.
- **Judge** scores methodology (primary) and outputs (secondary) against
  a per-task rubric. Lower is better; 0 = perfect.
- Keep the new SKILL.md only if it strictly beats the best score; else
  revert. Stop on target_score or on `early_stop_n` consecutive regressions.

## You are the orchestrator

You (the agent reading this) don't run the loop yourself. You dispatch
subagents to build the workspace, then hand off to the Python loop.

### Phase 1 — Scout

Dispatch a subagent with `prompts/scout.md` to research the paper/task.
Report key findings to the user in a few lines.

### Phase 2 — Scope (you + user)

Have a conversation. Ask ONE question at a time, multiple-choice where
helpful. Agree on:
- what to reproduce / what success looks like
- which data sources are in-bounds
- what methodology expectations belong in the rubric
- iteration budget and target_score (if any)

Present a summary and get approval.

### Phase 3 — Build

Dispatch a builder subagent with `prompts/builder.md` and the agreed
scope. It writes:
- `task.json`
- `rubric.md` — **the authoritative scoring rubric for the LLM judge**
- `reference/` (optional; judge-only)
- `skill/SKILL.md` — seed

Validate:
```python
from skills.clawpathy_autoresearch import validate_workspace
print(validate_workspace(Path("WORKSPACE")))  # [] means valid
```

### Phase 4 — Loop

```bash
python -m skills.clawpathy_autoresearch WORKSPACE_DIR
# or with custom models:
python -m skills.clawpathy_autoresearch WORKSPACE_DIR \
  --proposer-model sonnet --executor-model sonnet --judge-model opus
```

The loop streams progress to `WORKSPACE/history.jsonl`, snapshots every
iteration's skill to `WORKSPACE/snapshots/iter-NNN.md`, and writes the
executor's full transcript to `WORKSPACE/executor_runs/iter-NNN.log`.

## Workspace layout

```
workspace/
  task.json                  # task metadata + loop knobs
  rubric.md                  # LLM-judge rubric (the heart of the system)
  reference/                 # optional ground truth, judge-only
  skill/SKILL.md             # iterated by the loop
  output/                    # executor outputs (cleared each iter)
  executor_runs/iter-NNN.log # transcripts (judge reads these)
  snapshots/iter-NNN.md      # per-iter SKILL.md snapshots
  history.jsonl              # one row per iter: score, kept, verdict
```

## Key principles

- **LLM judge only.** No deterministic Python scorers. All evaluation goes
  through `judge.md` + opus. This keeps the system low-code and lets the
  rubric carry paper-specific nuance without adding code.
- **Methodology is primary.** The rubric weights "did the agent use sound
  methods?" above "did the numbers match?". Ground-truth match is a signal,
  not the objective — the goal is better SKILL.md files.
- **Never leak ground truth.** `reference/` is judge-only. The executor
  prompt says not to read it, and the judge penalises leakage.
- **No hardcoded answers in SKILL.md.** The proposer prompt and the judge
  both enforce this. The executor must derive results by running methods.
- **Snapshots + strict-better revert.** Score on the first iter becomes the
  floor. Later iters that tie or regress revert to the best.

## Safety

- All processing is local except scout web fetches for public resources.
- ClawBio disclaimer: research/education tool, not a medical device.

## Gotchas

- **Do not skip scoping.** The rubric is paper-specific; a generic rubric
  tunes nothing. Get the user to agree on methodology expectations.
- **Do not write a Python scorer.** Earlier versions of this project did.
  They rewarded API-fetching, not methodology. The judge is the scorer.
- **Do not hand-pick the "best" snapshot yourself.** Trust the loop. If
  the judge is calibrated wrong, fix the rubric, not the history.
