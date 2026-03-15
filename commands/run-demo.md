---
name: run-demo
description: Run a ClawBio skill demo with built-in sample data
---

Run a demo of a ClawBio skill. The user will specify which skill to demo (e.g., "pharmgx", "genome-compare", "scrna").

Steps:
1. Read `skills/catalog.json` to find the skill and its demo command
2. If the skill has `has_demo: true`, run its demo command
3. Show the generated report to the user and explain the results
4. Open any generated figures

If no skill is specified, show the list of skills that have demos available and ask the user to pick one.

Common demo commands:
```bash
python clawbio.py run pharmgx --demo
python clawbio.py run compare --demo
python clawbio.py run scrna --demo
python clawbio.py run clinpgx --demo
python clawbio.py run prs --demo
python clawbio.py run gwas --demo
python clawbio.py run profile --demo
python clawbio.py run ukb --demo
python clawbio.py run galaxy --demo
```
