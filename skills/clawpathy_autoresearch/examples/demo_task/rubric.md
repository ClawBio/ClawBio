# Rubric — demo GWAS sumstats analysis

The judge scores the executor's run against the items below. Weight
methodology (70%) over outputs (30%).

## Methodology (70%)

1. **Data handling is sound.** Reads the user-provided sumstats from the
   agreed input path (not from `reference/`). Validates columns, handles
   missing p-values, logs row counts. pass/partial/fail.
2. **Lead-variant definition is principled.** Uses p < 5e-8 threshold and
   defines "lead" by locus clumping OR distance-based pruning (e.g. ±500kb)
   — not just "all significant rows". States the choice explicitly. p/p/f.
3. **Lambda computation is correct.** Computes genomic inflation as the
   median of observed chi-squared stats divided by the expected median
   (0.4549). Handles p-values near 0/1 safely. p/p/f.
4. **QQ plot is honest.** Uses -log10(observed) vs -log10(expected uniform),
   includes y=x reference line, labels axes, saves as PNG >= 20KB. p/p/f.
5. **Assumptions are documented.** Inline reasoning about choices (chunking
   thresholds, tie handling, what to do when lambda is very inflated).
   Logged to transcript. p/p/f.

## Outputs (30%)

6. **Lead variants match expected set.** Compare against
   `reference/expected_leads.json`. Pass if set agreement >= 0.8; partial
   >= 0.5; fail below. p/p/f.
7. **Lambda within tolerance.** |computed - reference| <= 0.05. p/p/f.
8. **QQ plot present and non-trivial.** PNG/SVG in `output/`, file size
   indicates a real plot. p/p/f.

## Anti-leakage

- If the transcript shows the executor reading `reference/` or pasting
  expected values into its code, fail the whole methodology axis.
- If `skill/SKILL.md` contains hardcoded rsIDs, gene names, or numeric
  targets, fail item 5 and halve methodology_score.
