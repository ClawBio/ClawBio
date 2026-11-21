## Summary

<!-- What does this PR do? 1-3 bullet points. -->

## Skill affected

<!-- Which skill does this change? e.g. pharmgx-reporter, equity-scorer, or "new skill: my-skill-name" -->

## Checklist

- [ ] I have searched for [existing PRs](https://github.com/ClawBio/ClawBio/pulls)
- [ ] I have included a complete SKILL.md (YAML frontmatter + methodology)
- [ ] I have included tests and they pass (`pytest -v`), with a [property test](../AGENTS.md#property-based-tests-for-maths-and-parsing) for maths or parsing
- [ ] I have included demo data for reviewers to test, and a [demo snapshot test](../AGENTS.md#demo-snapshot-tests) next to tests of the edge cases (a snapshot alone does not test behaviour)
- [ ] I have not committed any patient/sensitive data
- [ ] Outputs (`report.md`, `result.json`, `reproducibility/`) hold no absolute paths, usernames or home folders; record file names, not where they live
- [ ] The skill writes its [reproducibility bundle](../AGENTS.md#code-style) with `clawbio.common.reproducibility`, not a hand-rolled writer
- [ ] Every gene-drug association, threshold or parameter I added cites its source ([no hallucinated science](../AGENTS.md#safety-boundaries))
- [ ] I have followed the [CONTRIBUTING.md](../CONTRIBUTING.md) conventions and [Before Opening a PR](../AGENTS.md#before-opening-a-pr) in AGENTS.md, which wins where they differ

## Evidence

<!-- For each claim above or in the summary: reproduced live (real data or service), reproduced synthetic (mock, fixture or fault injection), or read-only (from reading the code; say what would confirm it). -->

## Test output

<!-- Paste the pass/fail counts of `uv run pytest tests/` and of every skill you touched, plus the demo output -->
