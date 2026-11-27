---
name: analyse
description: Run a ClawBio bioinformatics analysis on user-provided genetic data
---

The user wants to analyse genetic data using ClawBio skills. Follow this workflow:

1. **Identify the input file**: Use the path to the user's data file (23andMe, AncestryDNA, VCF, h5ad, etc.) if they gave one. If they have no input file, do not refuse and do not just ask for one: run the demo in step 3
2. **Route to the right skill**: Match the user's intent against `trigger_keywords`, `tags`, and `description` in `skills/catalog.json`, then confirm with the skill's `SKILL.md` `## Trigger` section
3. **Run the analysis**:
   - With a file: execute the skill's Python script with appropriate arguments
   - With no file: tell the user the run uses the skill's demo data, not theirs, then run its `demo_command` from `skills/catalog.json`, adding `--output <dir>`. If `demo_command` is `null`, the skill is spec-only: walk its SKILL.md methodology instead
4. **Present results**: Show the generated report, open figures, and explain findings
5. **Offer follow-ups**: Suggest related skills (e.g., after pharmgx, suggest genome-compare or prs)

Important:
- All computation runs locally. Genetic data never leaves the machine.
- Always include the ClawBio disclaimer in output.
- If unsure which skill to use, run the bio-orchestrator: `python skills/bio-orchestrator/orchestrator.py --input <file>`
