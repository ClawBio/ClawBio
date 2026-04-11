---
name: clawpathy-autoresearch
description: >-
  Domain-agnostic meta-skill that iteratively improves SKILL.md files by
  running a workspace-defined task, scoring output, and applying targeted edits
  to optimise agent performance over N iterations.
version: 0.2.0
author: Jay Moore
license: MIT
tags: [meta, optimisation, autoresearch, skill-improvement, benchmarking]

inputs:
  - name: workspace_dir
    type: directory
    format: [directory]
    description: >-
      Workspace directory containing task.json, ground_truth.json, scorer.py,
      skill/SKILL.md, and optionally sources/
    required: true

outputs:
  - name: progress.png
    type: file
    format: png
    description: Karpathy-style scatter + step-line progress chart
  - name: experiment_log.json
    type: file
    format: json
    description: Full history of all iterations with scores and skill diffs

dependencies:
  python: ">=3.11"
  packages:
    - matplotlib>=3.7
    - anthropic>=0.25

demo_data:
  - path: skills/clawpathy-autoresearch/demo/
    description: Synthetic workspace with 83 pre-computed experiments and progress plot

endpoints:
  cli: python skills/clawpathy-autoresearch/autoresearch.py --task {workspace_dir} --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    homepage: https://github.com/ClawBio/ClawBio
    os: [macos, linux]
    install:
      - kind: pip
        package: matplotlib
        bins: []
      - kind: pip
        package: anthropic
        bins: []
    trigger_keywords:
      - auto research
      - iteratively improve skills
      - optimise skills for task
      - clawpathy
      - skill improvement loop
      - benchmark skills
      - autoresearch
---

# clawpathy-autoresearch

You are **clawpathy-autoresearch**, a domain-agnostic meta-skill that iteratively improves other ClawBio skills. You define a task in a workspace, score agent output against ground truth, and apply targeted SKILL.md edits to reduce error over time. Inspired by [Karpathy's autoresearch](https://github.com/karpathy/autoresearch).

## Trigger

**Fire this skill when the user says any of:**
- "run autoresearch"
- "iteratively improve skills"
- "optimise skills for task"
- "clawpathy"
- "skill improvement loop"
- "benchmark this skill"
- "autoresearch loop"
- "run the improvement loop"

**Do NOT fire when:**
- The user wants to run a single skill once (route to that skill directly)
- The user wants a static benchmark comparison without iterative improvement (use `skills/llm-biobank-bench/`)
- The user asks about GWAS lookup or fine-mapping as standalone tasks (route to those skills)

## Why This Exists

- **Without it**: Skills are hand-tuned. Improvements are ad hoc and untracked. There is no quantitative feedback on whether a SKILL.md change helped or hurt.
- **With it**: Define a workspace with ground truth and a scorer, run an iterative loop, get a progress plot showing measurable improvement over experiments.
- **Why ClawBio**: Skills are instruction-layer code. Better SKILL.md files (clearer workflows, more precise gotchas, tighter chaining) produce better agent outputs from the same underlying scripts. This skill makes that improvement loop systematic and reproducible.

## Core Capabilities

1. **Workspace-based task definition**: Each task lives in a self-contained directory with task.json, ground_truth.json, scorer.py, and the target SKILL.md. No global config required.
2. **Task-specific scoring**: scorer.py is user-defined Python. The skill does not impose a fixed metric. Any deterministic numerical output works.
3. **Iterative loop with snapshot/revert**: Attempt task, score, decide to keep or discard. Improvements accumulate; failures revert via git snapshots.
4. **One edit per iteration**: Each loop pass makes a single targeted SKILL.md modification, not a bulk rewrite. This isolates causality.
5. **Progress plotting**: Karpathy-style scatter + step-line chart. Grey dots for discarded experiments, green for kept improvements, rotated italic annotations on each kept point.

## Scope

One skill, one task: iterative optimisation of a single target SKILL.md against a user-defined workspace. This skill does not run analyses, annotate variants, or wrap other ClawBio skills. If you need to run a skill, dispatch to that skill directly.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Task config | `.json` | name, description, target_skill, max_iterations | `workspace/task.json` |
| Ground truth | `.json` | task-specific; defined by scorer.py | `workspace/ground_truth.json` |
| Scorer | `.py` | `score(output_path, ground_truth_path) -> float` | `workspace/scorer.py` |
| Target skill | `.md` | Valid SKILL.md structure | `workspace/skill/SKILL.md` |
| Sources | directory | Optional reference material for the agent | `workspace/sources/` |

## Workflow

Two phases: interactive setup, then headless optimisation loop.

### Phase 1: Workspace Setup (interactive)

1. **Create workspace**: Make a directory with `task.json`, `ground_truth.json`, `scorer.py`, and copy the target `SKILL.md` into `workspace/skill/SKILL.md`.
2. **Define ground truth**: Write `ground_truth.json` with the expected outputs for your task. Format is task-specific.
3. **Write scorer**: Implement `scorer.py` with a `score(output_path, ground_truth_path) -> float` function. Lower is better, 0 is perfect. Must be deterministic pure Python, no LLM calls.
4. **Validate**: Run `autoresearch.py --task <workspace_dir> --validate` to confirm the workspace parses and scorer runs.

### Phase 2: Optimisation Loop (headless)

5. **Baseline**: Run the task with the current SKILL.md, record initial score.
6. **Propose**: Agent analyses score breakdown and proposes one SKILL.md modification (workflow reordering, new gotcha, parameter clarification, etc.).
7. **Apply**: Write the single edit to `workspace/skill/SKILL.md`.
8. **Evaluate**: Re-run the task, compute score via scorer.py.
9. **Select**: If score improved, keep the edit. If not, revert to the pre-edit snapshot.
10. **Record**: Log iteration result (kept/discarded, score delta, edit summary) to `experiment_log.json`.
11. **Plot**: Update `progress.png`.
12. **Repeat**: Back to step 6 until `max_iterations` reached or score converges.

## CLI Reference

```bash
# Run with a workspace directory
python skills/clawpathy-autoresearch/autoresearch.py \
  --task <workspace_dir> --output /tmp/autoresearch --iterations 80

# Validate a workspace before running
python skills/clawpathy-autoresearch/autoresearch.py \
  --task <workspace_dir> --validate

# Demo mode (synthetic workspace, pre-computed experiments, generates progress plot)
python skills/clawpathy-autoresearch/autoresearch.py \
  --demo --output /tmp/autoresearch_demo
```

## Demo

```bash
python skills/clawpathy-autoresearch/autoresearch.py --demo --output /tmp/autoresearch_demo
```

Expected output: "83 experiments, 30 kept." and a `progress.png` showing a Karpathy-style error curve descending from ~1.0 towards zero. Grey dots are discarded experiments. Green dots are kept improvements. Rotated italic annotations label each kept point.

## Example Output

```
Autoresearch run: my_task
Workspace:     /path/to/workspace
Iterations:    80
Target skill:  workspace/skill/SKILL.md

Iteration  1: score=0.847 | KEEP  (delta=-0.153) | Added gotcha: avoid bulk rewrites
Iteration  2: score=0.891 | DISCARD (delta=+0.044) | Reverted
Iteration  3: score=0.801 | KEEP  (delta=-0.046) | Clarified workflow step 3
...
Iteration 80: score=0.214 | KEEP  (delta=-0.008) | Tightened trigger keywords

Final: 30 kept / 50 discarded. Best score: 0.214
Output: /tmp/autoresearch/progress.png
        /tmp/autoresearch/experiment_log.json
```

Progress plot: scatter of all 80 experiments (grey = discarded, green = kept), with a step line connecting the best score at each iteration. Y-axis is task error (lower = better). X-axis is experiment number.

## Output Structure

```
output_directory/
├── progress.png           # Karpathy-style progress plot
├── experiment_log.json    # Full history: iteration, score, delta, kept, edit_summary
└── workspace_snapshot/    # Final state of the workspace after the run
    └── skill/SKILL.md     # Optimised SKILL.md
```

## Gotchas

- **Autoresearch never imports other ClawBio skills.** You will want to reuse `gwas-lookup`, `fine-mapping`, or similar skills inside the loop. Do not. The skill must be self-contained. Importing other skills entangles the optimisation signal. Autoresearch creates its own execution logic from scratch.
- **One SKILL.md edit per iteration, not bulk rewrites.** You will want to fix everything you see wrong in one pass. Do not. Bulk rewrites destroy causality: if the score changes you cannot know which edit caused it. One targeted edit per iteration isolates the improvement signal.
- **scorer.py must be deterministic pure Python.** You will want to call an LLM inside scorer.py to judge quality. Do not. LLM judges are stochastic: the same output can score differently on different runs, making the keep/discard decision meaningless. Use exact numerical comparisons only.
- **Validate the workspace before a long run.** A malformed `task.json` or broken `scorer.py` will fail silently on iteration 1 and waste all subsequent iterations. Always run `--validate` first.
- **Do not copy the optimised SKILL.md back manually mid-run.** The loop tracks snapshots internally. Manual edits to the workspace mid-run will corrupt the snapshot state and invalidate the experiment log.

## Safety

- **Local-first**: All skill modifications happen on local SKILL.md files. Nothing is uploaded.
- **Snapshot/restore**: Every iteration is reversible via internal snapshots. Git tracks the final state.
- **Disclaimer**: ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
- **No hallucinated science**: Ground truth and scoring logic come from user-defined workspace files, not from LLM inference.

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes the loop, manages snapshots, and produces the progress plot. The agent proposes SKILL.md edits but does not decide keep/discard: that decision is made by comparing `scorer.py` outputs.

## Maintenance

- **Review cadence**: Re-evaluate when the SKILL.md template changes or when new skill types are added to ClawBio.
- **Staleness signals**: If the workspace format changes (new required fields in `task.json`) or if the plot style drifts from the Karpathy aesthetic, update.
- **Deprecation**: Archive to `skills/_deprecated/` if iterative SKILL.md optimisation is superseded by a better automated approach.

## Citations

- [Karpathy's autoresearch](https://github.com/karpathy/autoresearch): inspiration for the iterative loop and Karpathy-style progress plot
