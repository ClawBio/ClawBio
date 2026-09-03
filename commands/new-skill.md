---
name: new-skill
description: Scaffold a new ClawBio skill from the official template with conformance checks
---

# New Skill Builder

Build a new ClawBio skill from the official template with full conformance enforcement.

Input: $ARGUMENTS (skill name in lowercase-with-hyphens, e.g. "pathway-enrichment")

## Step 1: Validate the Name

Parse `$ARGUMENTS` to get the skill name. If empty, ask the user what the skill should do and derive a name.

Rules:
- Lowercase with hyphens only (e.g. `pathway-enrichment`, not `PathwayEnrichment`)
- Must not already exist in `skills/`
- Should be 2-4 words maximum

Store as `SKILL_NAME`. Also derive `SKILL_NAME_UNDERSCORE` by replacing hyphens with underscores (e.g. `pathway-enrichment` becomes `pathway_enrichment`). Use this for Python filenames and test files.

## Step 2: Interview the User

Before writing any code, ask these questions one at a time. Wait for each answer.

1. **What does this skill do?** (one sentence, this becomes the `description` in YAML)
2. **What input does it take?** (file format: VCF, CSV, TSV, TXT, JSON, h5ad, or free text query)
3. **What output does it produce?** (report, table, plot, JSON, or combination)
4. **What domain databases or algorithms does it use?** (e.g. ClinVar, gnomAD, PGS Catalog, custom logic)
5. **What should the agent NEVER do with this skill?** (this seeds the Gotchas and Agent Boundary)

## Step 3: Create the Skill Directory

```bash
mkdir -p skills/$SKILL_NAME/tests
mkdir -p skills/$SKILL_NAME/examples
```

## Step 4: Generate SKILL.md from Template

Read the template at `templates/SKILL-TEMPLATE.md`.

Fill in every section using the interview answers. Be specific and concrete:

### Trigger Section (MOST IMPORTANT)
- Write at least 5 "fire when" phrases covering synonyms, abbreviations, and natural language variations
- Write at least 2 "do NOT fire when" entries pointing to similar skills that handle adjacent queries
- Copy the trigger phrases into `trigger_keywords` in the YAML metadata

### Scope Section
- One sentence: "This skill does X and nothing else."
- If the user's description implies two tasks, flag it and recommend splitting

### Workflow Section
- Numbered steps only, no prose
- Mark each step as prescriptive (exact) or flexible (room for reasoning)
- Include validation as step 1

### Example Output Section
- Write an actual rendered sample based on the output format
- Use realistic-looking synthetic values, not placeholder text

### Gotchas Section
- Seed with at least 3 gotchas from the user's "never do" answer
- Add 1-2 common model failure patterns for this domain (e.g. hallucinating gene associations, inventing p-values, assuming ancestry)

### Maintenance Section
- Identify what upstream sources could make this skill stale
- Set review cadence (default: monthly)

Write the completed SKILL.md to `skills/$SKILL_NAME/SKILL.md`.

## Step 5: Generate Demo Data

Create a synthetic demo file that exercises the skill's logic. Rules:
- Never use real patient data
- Include at least 3-5 rows/entries to exercise basic logic
- Add a comment header: `# Synthetic demo data for $SKILL_NAME. This is NOT real patient data.`

Write to `skills/$SKILL_NAME/demo_<format>.<ext>`.

## Step 6: Write Tests First (Red Phase)

Create `skills/$SKILL_NAME/tests/test_$SKILL_NAME_UNDERSCORE.py` with:

1. `test_demo_runs_without_error` - the `--demo` flag produces output without crashing
2. `test_output_structure` - output directory contains expected files (report.md, result.json, etc.)
3. `test_report_contains_disclaimer` - every report includes the ClawBio medical disclaimer
4. `test_rejects_malformed_input` - bad input raises a clean error, not a traceback
5. `test_empty_input_handled` - empty file produces a meaningful error message

Run the tests and confirm they fail (red):

```bash
python -m pytest skills/$SKILL_NAME/tests/ -v
```

Show the user the red output.

## Step 7: Write the Python Implementation (Green Phase)

Create `skills/$SKILL_NAME/$SKILL_NAME_UNDERSCORE.py` with:

- `argparse` CLI with `--input`, `--output`, and `--demo` flags
- `run()` function that does the core work
- Pathlib for all paths, no hardcoded paths
- Output: `report.md` + `summary.json` minimum
- ClawBio disclaimer in every report
- Graceful error handling for bad input

Run the tests and confirm they pass (green):

```bash
python -m pytest skills/$SKILL_NAME/tests/ -v
```

Show the user the green output.

## Step 8: Run the Demo

```bash
python skills/$SKILL_NAME/$SKILL_NAME_UNDERSCORE.py --demo --output /tmp/$SKILL_NAME_demo
```

Read and display the generated report to the user.

## Step 9: Self-Audit (Conformance Check)

Run the 17-point SKILL.md conformance checklist against the new skill. Check each item:

| Check | Requirement |
|-------|------------|
| YAML: `name` | Present, matches folder name |
| YAML: `version` | Semver format |
| YAML: `author` | Present |
| YAML: `description` | One line, specific |
| YAML: `inputs` | Present with format and required flag |
| YAML: `outputs` | Present with format |
| YAML: `trigger_keywords` | At least 3 keywords |
| Section: `## Trigger` | Fire/do-not-fire lists present |
| Section: `## Scope` | One-skill-one-task confirmed |
| Section: `## Workflow` | Numbered steps, not prose |
| Section: `## Example Output` | Rendered sample present |
| Section: `## Gotchas` | At least 3 entries |
| Section: `## Safety` | Disclaimer referenced |
| Section: `## Agent Boundary` | Present |
| File: demo data | At least one demo file |
| File: tests/ | Directory with at least one test |
| Line count | SKILL.md under 500 lines |

Report PASS/FAIL for each. Fix any failures before proceeding.

## Step 10: Update Routing Table

Add the new skill to the ClawBio CLAUDE.md routing table:
- Add a row to the `## Skill Routing Table` with user intent phrases, skill path, and action
- Add CLI reference to the `## CLI Reference` section
- Add demo data to the `## Demo Data` table
- Add demo command to the `## Demo Commands` section

## Step 11: Summary

Show the user:
1. Files created (list all with paths)
2. Conformance checklist result (17/17 or what needs fixing)
3. How to run: the exact CLI command
4. How to trigger: the phrases an agent will respond to
5. Next steps: stress test with 10 varied inputs, add gotchas for every failure

Do NOT commit unless the user asks.
