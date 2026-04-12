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

demo_data:
  - path: skills/clawpathy-autoresearch/demo_workspace/
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

Three phases: multi-agent setup, headless optimisation loop, independent audit.

The setup phase is gated by three sequential subagents that prevent the
controller from making assumptions. Grounded in published systems:
task-clarifier (setup gate), prior-art-check (Sakana AI Scientist's
Semantic Scholar novelty check), ground-truth-auditor (Coscientist
reality-gate / OpenAI Deep Research citation pattern). The loop is gated
per iteration by a reflection critic (Google Co-Scientist's Reflection
agent). The final phase is a reproduction auditor on a held-out split
(MLAgentBench pattern).

### Phase 1: Workspace Setup (multi-agent gated)

1. **Clarify task** via `clarifier-prompt.md` subagent. Loop until
   `ready=true`. Relay missing questions to the user.
2. **Verify prior art** via `prior-art-prompt.md` subagent. Confirms the
   paper + target results exist and are cited correctly.
3. **Scaffold**: Run `autoresearch.py --setup <path>` to create
   `task.json`, `ground_truth.json`, `scorer.py`, `skill/SKILL.md`.
4. **Use a dev/held-out split** in `ground_truth.json`: `{"dev": {...},
   "heldout": {...}}`. The loop optimises on `dev`. The reproduction
   auditor verifies on `heldout`.
5. **Audit ground truth** via `ground-truth-auditor-prompt.md` subagent.
   Independently verifies every numerical value against the cited source.
6. **Validate**: `autoresearch.py --task <workspace> --validate`.

### Phase 2: Optimisation Loop (critic-gated)

7. **Baseline**: Run the task with the current SKILL.md, record initial score.
8. **Propose**: LLM proposer rewrites SKILL.md (one targeted edit).
9. **Critic review**: `critic.py` reviews the diff. Rejects edits that
   bake ground-truth values into the skill (memorisation) or are
   cosmetic. Rejected edits are reverted and logged as
   `critic-rejected: <label>` without running the agent.
10. **Evaluate**: Run agent on the new SKILL.md, score via scorer.py
    against the `dev` split.
11. **Select**: If score improved, keep the edit. If not, revert.
12. **Record + plot**: Log to `experiment_log.json`, update `progress.png`.
13. **Repeat**: Back to step 8 until `max_iterations` or early-stop.

### Phase 3: Reproduction Audit (independent)

14. **Run final SKILL.md against `heldout` split.**
15. **Dispatch `reproduction-auditor-prompt.md` subagent.** Independent
    judge never saw `dev`. Decides whether the skill genuinely reproduces
    the paper's claim or overfit the loop's scorer.

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

## Chaining Partners

- **Any ClawBio skill with a Python entrypoint**: autoresearch subclasses `AutoResearchLoop` and implements `run_agent_on_skill` to dispatch the target skill per iteration. It does not import target skills directly — the workspace defines the contract.
- **`final_judge.py`**: runs at the end of the loop, producing a qualitative score and concrete SKILL.md addition proposals via parallel LLM judges.
- **`proposer.py`**: runs every iteration, generating a single targeted SKILL.md edit from the current skill text plus experiment history.

## Maintenance

- **Review cadence**: Re-evaluate when the SKILL.md template changes or when new skill types are added to ClawBio.
- **Staleness signals**: If the workspace format changes (new required fields in `task.json`) or if the plot style drifts from the Karpathy aesthetic, update.
- **Deprecation**: Archive to `skills/_deprecated/` if iterative SKILL.md optimisation is superseded by a better automated approach.

## Citations

- [Karpathy's autoresearch](https://github.com/karpathy/autoresearch): inspiration for the iterative loop and Karpathy-style progress plot
