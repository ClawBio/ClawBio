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
2. **Concrete error scoring**: Mean reproduction error from numerical metrics (lower is better, 0 = perfect)
3. **Skill management**: Read, modify, create, snapshot, and rollback SKILL.md files
4. **Iterative loop**: Attempt task, score, modify skills, retry. Keep improvements, revert failures
5. **Progress plotting**: Karpathy-style scatter + step-line chart showing error descending towards zero

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Task config | `.yaml` | name, metric, direction, papers | `tasks/gwas_reproduction/task.yaml` |
| Ground truth | `.yaml` | lead_variants (rsid, gene, neg_log10_p, odds_ratio, effect_allele_freq, effect_direction), total_loci | `ground_truth/paper_001.yaml` |

## Workflow

When the user asks to iteratively improve skills:

1. **Load task**: Parse task.yaml, resolve ground truth files
2. **Snapshot skills**: Capture current state of all SKILL.md files
3. **Attempt task**: Agent runs the task using current skills
4. **Score output**: Compute mean reproduction error against ground truth (lower is better)
5. **Decide**: If error decreased, keep changes. If not, revert to snapshot
6. **Modify skills**: Agent proposes SKILL.md modifications based on feedback
7. **Plot**: Update progress.png with latest experiment
8. **Repeat**: Loop for N iterations

## CLI Reference

```bash
# Run with a defined task
python skills/clawpathy-autoresearch/autoresearch.py \
  --task tasks/gwas_reproduction/task.yaml \
  --output /tmp/autoresearch --iterations 80

# Demo mode (synthetic data, generates example progress plot)
python skills/clawpathy-autoresearch/autoresearch.py \
  --demo --output /tmp/autoresearch_demo
```

## Demo

To verify the skill works:

```bash
python skills/clawpathy-autoresearch/autoresearch.py --demo --output /tmp/autoresearch_demo
```

Expected output: "83 experiments, 30 kept." plus a `progress.png` showing a Karpathy-style error curve descending from ~1.0 towards zero, with grey discarded dots, green kept improvements, and rotated italic annotations on each kept point.

## Algorithm / Methodology

### The Loop

1. **Baseline**: Run the task with current skills, establish initial score
2. **Propose**: Agent analyses score breakdown and proposes skill modifications
3. **Apply**: Modify SKILL.md files (workflow reordering, new gotchas, parameter changes, new skills)
4. **Evaluate**: Re-run task with modified skills, compute mean reproduction error
5. **Select**: If error < best_error, keep changes (commit). Otherwise, revert to snapshot
6. **Record**: Log experiment result (kept/discarded, score, label, skill diff)
7. **Plot**: Update progress chart
8. **Repeat**: Back to step 2

### Error Metrics (Lower is Better, 0 = Perfect)

Score is a weighted mean of six concrete numerical error components:

| Component | Weight | Calculation |
|-----------|--------|-------------|
| **P-value error** | 0.20 | Normalised \|target - reproduced\| / target on neg_log10_p |
| **Odds ratio error** | 0.25 | Absolute \|target - reproduced\| on OR |
| **Effect allele freq error** | 0.10 | Absolute difference in frequency |
| **Locus count error** | 0.15 | Normalised \|target - reproduced\| / target |
| **Variant missing penalty** | 0.20 | 1.0 per missing variant / total variants |
| **Direction error** | 0.10 | 1.0 per wrong risk/protective call / total variants |

Total error = weighted sum. Zero means perfect reproduction. No LLM judge, no vague qualitative rubrics: just numbers.

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
- `anthropic` >= 0.25 — agent integration for real task execution

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
