---
marp: true
theme: default
paginate: true
title: clawpathy-autoresearch — eval-driven skill tuning
---

<!-- _class: lead -->

# clawpathy-autoresearch

### Eval-driven skill tuning for LLM agents

Jay Moore · UK DRI · ClawBio
2026

---

## The bottleneck isn't the model — it's the skill

Agents today are steered by **`SKILL.md`** files: prose instructions that tell
a model *how* to do a domain task (run LDSC, fine-map a locus, call variants…).

- The model is fixed. The **skill** is what makes or breaks a run.
- Skills are hand-written, rarely tested, and silently drift as tools change.
- Every domain expert reinvents the same 200-line playbook.

> If you could *measure* a skill, you could *improve* it.

---

## Idea: treat a SKILL.md like model weights

Borrow the loop that made deep learning work:

1. **Propose** a change to the skill.
2. **Execute** the skill end-to-end on a real task.
3. **Judge** the run against a rubric.
4. **Keep** it only if it beats the champion. Otherwise revert.

Repeat until you plateau.

The skill *is* the learned artifact. The loop is gradient descent, but the
"gradient" is an LLM judge's verdict on a methodology rubric.

---

## The loop

```
  ┌──────────────┐      ┌──────────────┐      ┌──────────────┐
  │  Proposer    │─────▶│  Executor    │─────▶│   Judge      │
  │  (Sonnet)    │      │  (Sonnet +   │      │   (Opus +    │
  │  rewrites    │      │   shell)     │      │   rubric)    │
  │  SKILL.md    │      │  runs skill  │      │   scores run │
  └──────────────┘      └──────────────┘      └──────────────┘
         ▲                                            │
         └───── verdict + recommended edits ──────────┘

  Keep new SKILL.md iff score < best. Else revert.
  Stop on target score or N consecutive regressions.
```

Three modes: **serial**, **K=3 parallel** (diverse candidates per iter),
and **multi-task** (one skill, N traits, held-out val/test).

---

## Why LLM-as-judge

Classical RL needs a programmable reward. Scientific methodology doesn't have one.

- A rubric is prose: *"harmonises alleles correctly," "reports credible sets
  with plausible sizes," "declares provenance of LD panel."*
- Opus-4 reads the run transcript + outputs and scores against the rubric.
- Low-code: no unit tests, no gold data, no task-specific Python.
- The same harness handles GWAS fine-mapping, LDSC heritability, RNA-seq QC…

Trade-off: the judge is noisy. We compensate with K-way parallel candidates
per iter and only keep strict improvements.

---

## What you get out

A **champion SKILL.md** — a frozen, measured, reproducible playbook.

| Task                                   | Judge score (↓ better) |
|----------------------------------------|------------------------|
| Trubetskoy 2022 PGC3 SCZ fine-mapping | 0.235                  |
| Yengo 2022 height LDSC h²              | ~0.12                  |
| Yengo → BMI transfer (no retuning)     | 0.42 → 0.19            |

Champions live in the repo as standalone artifacts — downloadable with
one function, portable to any agent, diffable like code.

---

<!-- _class: lead -->

## Takeaway

You can't optimise what you can't measure.

`clawpathy-autoresearch` measures skills, so skills can compound.

*Repo:* `github.com/ClawBio/ClawBio/tree/main/skills/clawpathy_autoresearch`
