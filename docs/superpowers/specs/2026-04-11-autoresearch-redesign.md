# Autoresearch Redesign — General-Purpose Meta-Skill

**Date:** 2026-04-11
**Scope:** Rewrite autoresearch as a domain-agnostic meta-skill that creates and iteratively improves SKILL.md files for any reproducible task.

## Decisions

- **Domain-agnostic**: autoresearch knows nothing about genomics, GWAS, or any specific domain. All domain knowledge enters through conversation and gets encoded into the SKILL.md it creates.
- **Two phases**: interactive setup conversation, then headless optimisation loop.
- **Optimisation target**: SKILL.md text (natural language instructions for an agent). Not Python parameters, not code.
- **Scoring**: Python function generated during setup. Quantitative ground truth only (no human-in-the-loop judging during the loop).
- **Ground truth**: derived from source material by autoresearch during setup. User does not manually provide target values.
- **Iteration budget**: fixed max with convergence-based early stopping (N consecutive non-improvements).
- **Standalone**: never imports or wraps existing ClawBio skills. The generated SKILL.md is self-contained.

## Architecture

### Phase 1: Interactive Setup

A conversation between autoresearch and the user that produces a task workspace. Steps:

1. **Task elicitation**: "What do you want to accomplish?" Follow-up questions to clarify scope, inputs, expected outputs, success criteria.

2. **Source material ingestion**: user provides papers, datasets, or reference material. Autoresearch reads them and extracts relevant structure (results tables, methods sections, expected outputs).

3. **Ground truth derivation**: from the source material, extract quantitative targets. Store as structured JSON. For a GWAS paper: reported p-values, ORs, frequencies from the results table. For a different domain: whatever the paper reports as its key findings.

4. **Scoring function generation**: write a Python function (`scorer.py`) in the workspace. Signature: `score(skill_output: dict, ground_truth: dict) -> float`. Lower is better. Components and weights derived from the conversation. The function is deterministic and fast (no LLM calls).

5. **Initial SKILL.md generation**: write the first version of the skill instructions. These tell an agent how to accomplish the task end-to-end. For GWAS: "read the paper's methods section, obtain the raw data, run the described analysis pipeline, extract lead variant statistics." The SKILL.md contains no self-evaluation logic: scoring is external.

6. **Workspace output**:
   ```
   tasks/<task-name>/
   ├── task.json          # name, description, iteration budget, early_stop_n
   ├── ground_truth.json  # derived targets
   ├── scorer.py          # generated scoring function
   ├── sources/           # papers, data files, references
   └── skill/
       └── SKILL.md       # the optimisation target
   ```

### Phase 2: Headless Optimisation Loop

Runs autonomously after setup completes. Each iteration:

1. **Snapshot** current SKILL.md via SkillManager.
2. **Execute**: call the `run_agent_on_skill(skill_path, sources_dir) -> dict` function, which reads the SKILL.md, follows the instructions against the source material, and returns structured output as a dict matching the ground truth schema.
3. **Score**: run `scorer.py(output, ground_truth)` → float + component breakdown.
4. **Compare**: if score < previous best → **keep**. Otherwise → **revert** to snapshot.
5. **Propose changes**: analyse the error breakdown, identify the highest-contributing component, generate a text edit to the SKILL.md that addresses it. Examples:
   - "OR error is high → add workflow step: use trimmed median across reported effect sizes"
   - "Missing data penalty → add gotcha: check supplementary tables, not just main results"
   - "Direction errors → refine instruction: verify allele coding against methods section"
6. **Log** experiment to `experiment_log.json` (experiment number, score, component breakdown, kept/discarded, change description).
7. **Plot** progress via existing `plotter.py` (Karpathy-style scatter + step line).
8. **Early stop** if `early_stop_n` consecutive non-improvements (default: 5).

### Output

After loop completion, autoresearch produces:
- Final optimised SKILL.md
- `experiment_log.json` with full history
- Progress plot (PNG)
- Summary report (markdown): convergence history, phase analysis, final error breakdown, failure patterns

## Component Map

| Component | Current State | After Redesign |
|-----------|--------------|----------------|
| `autoresearch.py` | Hardcoded GWAS task, loads `task.yaml` | General-purpose: interactive setup or loads existing workspace |
| `real_runner.py` | Wraps gwas-lookup, queries databases | **Deleted**. Agent follows SKILL.md directly. |
| `scorer.py` | Fixed 6-component GWAS error function | **Generated per task** during setup. Module-level `score()` function. |
| `skill_manager.py` | Snapshot/restore SKILL.md files | **Unchanged**. Same mechanism, works on any SKILL.md. |
| `plotter.py` | Karpathy-style progress plot | **Unchanged**. Already generic (takes `ExperimentRecord` list). |
| `task.py` | `TaskDefinition` + `PaperDefinition`, loads `task.yaml` | **Replaced** by `task.json` loader. Simpler: just name, description, budget, early_stop_n. |
| Ground truth YAMLs | Pre-written per paper in `tasks/gwas_reproduction/` | **Replaced** by `ground_truth.json` derived during setup. |

## CLI Interface

```bash
# Interactive setup (starts conversation, produces workspace)
python skills/clawpathy-autoresearch/autoresearch.py \
  --setup --output /tmp/autoresearch_workspace

# Run loop on existing workspace
python skills/clawpathy-autoresearch/autoresearch.py \
  --task /tmp/autoresearch_workspace --iterations 80

# Demo (pre-built GWAS workspace, no setup needed)
python skills/clawpathy-autoresearch/autoresearch.py \
  --demo --output /tmp/autoresearch_demo
```

## Setup Conversation Flow

The setup is a normal conversation. Autoresearch asks questions one at a time:

1. "What do you want to accomplish?"
2. "What source material do you have?" (papers, datasets, URLs)
3. Reads the source material, then: "I found these quantitative results in the paper: [table]. I'll use these as ground truth for scoring. Does this look right?"
4. "Which of these metrics matter most for reproduction accuracy?" → derives scoring weights
5. "Here's the scoring function I've generated: [summary]. And here's the initial SKILL.md: [summary]. Ready to start the optimisation loop?"

The conversation is open-ended: autoresearch asks as many questions as needed to understand the task. The user can also volunteer information at any point.

## Propose Changes Strategy

The `propose_skill_changes()` function is the core intelligence. After each iteration:

1. Read the error breakdown (which components contributed most to the score).
2. Read the agent's execution trace (what did it do, where did it go wrong).
3. Identify a single, specific improvement to the SKILL.md text:
   - **Add a gotcha**: "The model will want to X. Do not. Instead Y."
   - **Refine a workflow step**: make an instruction more specific or add a sub-step.
   - **Add methodological detail**: incorporate a technique from the source material's methods section that the skill was missing.
   - **Remove a misleading instruction**: if a previous edit made things worse, revert it.
4. Apply the edit and describe it in the experiment log.

One change per iteration. Small, targeted edits. Never rewrite the entire SKILL.md at once.

## Demo Mode

`--demo` ships a pre-built workspace for the GWAS reproduction task:
- `sources/`: three landmark GWAS papers (Nalls 2019, Mahajan 2018, de Lange 2017)
- `ground_truth.json`: 10 variants with target p-values, ORs, frequencies, directions
- `scorer.py`: the 6-component weighted error function (p-value 0.20, OR 0.25, freq 0.10, locus_count 0.15, missing 0.20, direction 0.10)
- `skill/SKILL.md`: initial GWAS reproduction instructions

For quick visual demo without running the loop, the existing synthetic 83-experiment plot generation remains available.

## Constraints

- Autoresearch never imports or calls other ClawBio skills. The SKILL.md it generates is self-contained.
- The scoring function contains no LLM calls. It is pure Python operating on structured data.
- The setup conversation must complete before the loop starts. No interleaving.
- Each iteration proposes exactly one SKILL.md edit. No bulk rewrites.
- The SKILL.md must conform to the ClawBio SKILL-TEMPLATE.md format.

## Error Handling

- If the agent fails to produce parseable output during an iteration, that experiment scores 1.0 (maximum penalty) and is discarded.
- If `scorer.py` raises an exception, the experiment is logged as "error" and discarded.
- If the source material is unreadable or contains no extractable ground truth, setup fails with a clear message.
- If the workspace is missing required files (`task.json`, `ground_truth.json`, `scorer.py`, `skill/SKILL.md`), the loop refuses to start.

## Safety

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
