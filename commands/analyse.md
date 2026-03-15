---
name: analyse
description: Run a ClawBio bioinformatics analysis on user-provided genetic data
---

The user wants to analyse genetic data using ClawBio skills. Follow this workflow:

1. **Identify the input file**: Ask the user for the path to their data file (23andMe, AncestryDNA, VCF, h5ad, etc.)
2. **Route to the right skill**: Read CLAUDE.md's skill routing table to match the user's intent to the correct skill
3. **Run the analysis**: Execute the skill's Python script with appropriate arguments
4. **Present results**: Show the generated report, open figures, and explain findings
5. **Offer follow-ups**: Suggest related skills (e.g., after pharmgx, suggest genome-compare or prs)

Important:
- All computation runs locally. Genetic data never leaves the machine.
- Always include the ClawBio disclaimer in output.
- If unsure which skill to use, run the bio-orchestrator: `python skills/bio-orchestrator/orchestrator.py --input <file>`
