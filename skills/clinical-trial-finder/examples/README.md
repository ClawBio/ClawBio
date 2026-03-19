# examples/

Example input and output for clinical-trial-finder.

## Input format

The skill accepts a plain text file with one search term per line.
Lines starting with `#` are comments. See `../demo_input.txt`:

```
# Demo query for clinical-trial-finder
# This is NOT real patient data -- synthetic search terms only
# Searches for trials related to BRCA1 in breast cancer
BRCA1 breast cancer
```

Alternatively, use `--query "EGFR lung cancer"` or `--gene BRCA1` directly.

## Output structure

Running `--demo --fhir` produces:

```
demo_output/
  report.md              # human-readable markdown report
  report.html            # self-contained HTML report (open in browser)
  summary.json           # machine-readable JSON for skill chaining
  fhir_bundle.json       # FHIR R4 Bundle of ResearchStudy resources
  tables/
    trials.csv           # CSV for Excel/R/pandas import
  commands.sh            # exact command to reproduce this run
  checksums.sha256       # SHA-256 of all outputs (verify with sha256sum -c)
```

### report.md

Markdown with a summary table, per-trial sections with NCT links, status labels,
conditions, interventions, and the ClawBio safety disclaimer.

### summary.json

```json
{
  "query": "BRCA1 breast cancer",
  "timestamp": "2026-03-19T...",
  "source": "clinicaltrials.gov/api/v2",
  "total": 20,
  "recruiting": 1,
  "trials": [ { "nct_id": "NCT...", "title": "...", ... } ]
}
```

This is the primary interface for skill chaining -- `profile-report` reads
this file to embed trials in the unified patient report.

### fhir_bundle.json

FHIR R4 Bundle with `ResearchStudy` resources. Status and phase codes follow
the published R4 value sets. Conditions use MeSH coding from CT.gov's
`derivedSection` when available.

### tables/trials.csv

Flat CSV with columns: `nct_id, title, status, phase, study_type, start_date,
completion_date, conditions, interventions`. List fields are pipe-delimited
(`Cancer A | Cancer B`).

## Reproduce

```bash
cd demo_output
sha256sum -c checksums.sha256
# All OK
```

Or regenerate from scratch:

```bash
python skills/clinical-trial-finder/clinical_trial_finder.py \
  --demo --fhir --output examples/demo_output
```
