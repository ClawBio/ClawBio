# demo_task — workspace skeleton

This is a skeleton, not a runnable demo. To make it runnable, a builder
agent (or you) needs to populate `reference/` with:

- `sumstats.tsv` — a synthetic GWAS summary stats table
- `expected_leads.json` — the set of lead rsIDs the judge expects
- `expected_lambda.json` — `{ "lambda": <float> }`

The executor reads `sumstats.tsv` from a path agreed during scoping (NOT
from `reference/`), so you'll also want to copy the sumstats to e.g.
`data/sumstats.tsv` and mention that path in `task.json`'s description.

In real use, Phase 3 (builder) generates all of this for you from the
paper/task. This skeleton is just here to show the file layout.
