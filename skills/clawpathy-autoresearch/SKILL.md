---
name: clawpathy-autoresearch
description: >-
  Iterative skill improvement loop inspired by Karpathy's autoresearch.
  Defines a task, scores agent performance, modifies SKILL.md files to
  optimise results, and plots improvement over experiments.
version: 0.1.0
author: Jay Moore
license: MIT
tags: [meta, optimisation, autoresearch, skill-improvement, benchmarking]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🔬"
    homepage: https://github.com/ClawBio/ClawBio
    os: [macos, linux]
    install:
      - kind: pip
        package: pyyaml
        bins: []
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

# 🔬 clawpathy-autoresearch

You are **clawpathy-autoresearch**, a meta-skill that iteratively improves other ClawBio skills by optimising their SKILL.md files against a defined task. Inspired by [Karpathy's autoresearch](https://github.com/karpathy/autoresearch).

## Why This Exists

- **Without it**: Skills are hand-tuned, improvements are ad hoc, no quantitative feedback on whether changes helped
- **With it**: Define a task with ground truth, run an iterative loop, get a Karpathy-style progress plot showing measurable improvement
- **Why ClawBio**: Skills are instruction-layer code. Better SKILL.md files (workflows, gotchas, chaining) produce better agent outputs from the same underlying scripts

## Core Capabilities

1. **Task definition**: YAML-based task configs with ground truth for scoring
2. **Hybrid scoring**: Automated numerical checks (60%) + LLM-as-judge qualitative (40%)
3. **Skill management**: Read, modify, create, snapshot, and rollback SKILL.md files
4. **Iterative loop**: Attempt task, score, modify skills, retry. Keep improvements, revert failures
5. **Progress plotting**: Karpathy-style scatter + step-line chart showing improvement over experiments

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Task config | `.yaml` | name, metric, direction, papers | `tasks/gwas_reproduction/task.yaml` |
| Ground truth | `.yaml` | lead_variants, qualitative_findings, total_loci | `ground_truth/paper_001.yaml` |

## Workflow

When the user asks to iteratively improve skills:

1. **Load task**: Parse task.yaml, resolve ground truth files
2. **Snapshot skills**: Capture current state of all SKILL.md files
3. **Attempt task**: Agent runs the task using current skills
4. **Score output**: Hybrid scoring (automated + LLM judge) against ground truth
5. **Decide**: If score improved, keep changes. If not, revert to snapshot
6. **Modify skills**: Agent proposes SKILL.md modifications based on feedback
7. **Plot**: Update progress.png with latest experiment
8. **Repeat**: Loop for N iterations

## CLI Reference

```bash
# Run with a defined task
python skills/clawpathy-autoresearch/autoresearch.py \
  --task tasks/gwas_reproduction/task.yaml \
  --output /tmp/autoresearch --iterations 20

# Demo mode (synthetic data, generates example progress plot)
python skills/clawpathy-autoresearch/autoresearch.py \
  --demo --output /tmp/autoresearch_demo

# Custom scoring weights
python skills/clawpathy-autoresearch/autoresearch.py \
  --task tasks/gwas_reproduction/task.yaml \
  --auto-weight 0.7 --llm-weight 0.3 \
  --output /tmp/autoresearch
```

## Demo

To verify the skill works:

```bash
python skills/clawpathy-autoresearch/autoresearch.py --demo --output /tmp/autoresearch_demo
```

Expected output: "50 experiments, N kept." plus a `progress.png` showing a Karpathy-style improvement curve with grey discarded dots, green kept improvements, and annotations on each kept point.

## Algorithm / Methodology

### The Loop

1. **Baseline**: Run the task with current skills, establish initial score
2. **Propose**: Agent analyses score breakdown and proposes skill modifications
3. **Apply**: Modify SKILL.md files (workflow reordering, new gotchas, parameter changes, new skills)
4. **Evaluate**: Re-run task with modified skills, compute hybrid score
5. **Select**: If score > best_score, keep changes (commit). Otherwise, revert to snapshot
6. **Record**: Log experiment result (kept/discarded, score, label, skill diff)
7. **Plot**: Update progress chart
8. **Repeat**: Back to step 2

### Hybrid Scoring

**Automated (60% weight):**
- Variant recovery: fraction of ground truth rsIDs found (35%)
- Direction accuracy: correct risk/protective calls (25%)
- Effect size accuracy: OR within expected range (25%)
- Locus count accuracy: within tolerance of reported total (15%)

**LLM Judge (40% weight):**
- Biological pathways identified (0-2)
- Methodology appropriate (0-2)
- Qualitative findings matched (0-2)
- Relevant skills created/modified (0-2)
- Overall coherence (0-2)

### Skill Modification Types

- **Optimise**: Modify existing SKILL.md (workflow, gotchas, parameters)
- **Create**: Generate new skills when gaps identified
- **Compose**: Change skill chaining and workflow ordering

## Example Queries

- "Run autoresearch to improve GWAS reproduction skills"
- "Iteratively optimise my skills for this task"
- "Show me the improvement plot from the last autoresearch run"
- "Define a new benchmark task for scRNA-seq analysis"

## Output Structure

```
output_directory/
├── progress.png           # Karpathy-style progress plot
├── experiment_log.json    # Full history of all iterations
└── results/               # Per-experiment outputs (if applicable)
```

## Dependencies

**Required:**
- `pyyaml` >= 6.0 — task definition parsing
- `matplotlib` >= 3.7 — progress plot generation

**Optional:**
- `anthropic` >= 0.25 — LLM judge scoring (graceful fallback to score 5.0 without it)

## Safety

- **Local-first**: All skill modifications happen on local SKILL.md files
- **Snapshot/restore**: Every iteration is reversible via git-tracked snapshots
- **Disclaimer**: ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
- **No hallucinated science**: Ground truth comes from published papers with PMIDs

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User mentions "auto research", "iteratively improve", "optimise skills"
- User mentions "clawpathy" or "autoresearch"
- User wants to benchmark or score skill performance

**Chaining partners** — this skill connects with:
- `gwas-lookup`: Primary skill optimised in the GWAS reproduction task
- `fine-mapping`: Created/optimised during iterative improvement
- `pubmed-summariser`: Used for paper context retrieval
- `bio-orchestrator`: Routes meta-queries here

## Citations

- [Karpathy's autoresearch](https://github.com/karpathy/autoresearch) — inspiration for the iterative loop and progress plot
- GWAS papers used in the benchmark task are cited in each ground_truth/*.yaml file with PMIDs
