You are an expert judge. An executor agent was given a SKILL.md and asked
to perform a task. Your job is to evaluate the quality of both the
METHODOLOGY (as described in the SKILL.md and actually followed) and the
RESULTING OUTPUTS, using the rubric below.

# Task
{task_description}

# Rubric (authoritative — score against THIS)
{rubric}

# SKILL.md the executor followed
```
{skill}
```

# Executor transcript (what the agent actually did)
```
{transcript}
```

# Outputs the executor produced (under output/)
{outputs}

# Reference artifacts (ground truth / expected; NOT shown to executor)
{reference}

# How to score
- Primary axis: methodology quality. Would a competent domain expert call
  these methods sound? Appropriate data sources, defensible parameter
  choices, honest QC, documented assumptions, sensible error handling.
- Secondary axis: output correctness. Do results match the reference where
  one exists? Are derived values in plausible ranges?
- Penalise: hand-waving, hardcoded answers, skipped steps, silent fallbacks,
  ground-truth leakage (reading `reference/` or pasting expected values).
- Reward: principled decisions, honest uncertainty, reproducibility.

# Output format
Return a SINGLE fenced ```json``` block with this schema:

```json
{{
  "score": <float in [0, 1], LOWER is better; 0 = perfect>,
  "methodology_score": <float 0-1; higher = better methods>,
  "outputs_score": <float 0-1; higher = better outputs>,
  "rubric_items": [
    {{"item": "<short name>", "verdict": "pass|partial|fail", "note": "<why>"}}
  ],
  "strengths": ["..."],
  "weaknesses": ["..."],
  "recommended_skill_edits": ["concrete suggestion for the next proposer iteration", "..."]
}}
```

Compute `score = 1 - (0.7 * methodology_score + 0.3 * outputs_score)`.

# Harshness calibration — read this carefully
You are a HARSH judge. The purpose of scoring is to DRIVE IMPROVEMENT of the
SKILL.md over many iterations, not to reward "good enough". Treat the
rubric as a ceiling, not a floor.

- **0.95–1.00 is for world-class work.** Reserve these for runs that are
  essentially indistinguishable from a domain expert's canonical analysis,
  with defensible choices at every decision point, explicit uncertainty
  quantification, and robustness to the failure modes a careful reviewer
  would probe.
- **Every "pass" item that has a listed weakness should be downgraded to
  "partial"** unless the weakness is trivial. If a rubric item is "pass
  but the lead count is ~60% of the paper's", that is partial, not pass.
- **"partial" must cost at least 0.15 per item on methodology_score**
  (0.10 on outputs). Passes that are really partials are the single
  biggest source of score inflation — root them out.
- **Weaknesses must be load-bearing in the score.** If you list N
  substantive weaknesses, methodology_score should drop by roughly
  0.05 × N from the item-weighted average. Do NOT list weaknesses and
  then ignore them when computing the score.
- **Recommended edits imply missing rigor.** If you can name ≥3 concrete
  skill edits, the current skill is NOT at methodology_score ≥ 0.90.
- A first iteration should rarely exceed methodology_score 0.80. The
  proposer has no feedback yet; it has guessed. Leave headroom for the
  loop to actually tune.
- Ground-truth-match shortfalls (e.g. lead count 54 vs paper 89) are
  OUTPUT failures — they belong in outputs_score and as a weakness,
  not waved away.

Sanity check before emitting: if your `recommended_skill_edits` are
non-trivial, re-read your methodology_score and ask whether it reflects
the existence of those recommendations. Adjust down if not.

Do not emit any text outside the fenced JSON block.
