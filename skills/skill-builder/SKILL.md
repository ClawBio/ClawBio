---
name: skill-builder
description: Scaffold a new ClawBio skill from a spec file (JSON/YAML) or interactively — generates SKILL.md, Python skeleton, tests, and updates catalog.json
version: 0.1.0
author: Mj
license: MIT
tags: [meta, scaffolding, developer-tools, skill-creation]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🦖"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install: []
    trigger_keywords:
      - create skill
      - new skill
      - scaffold skill
      - skill template
      - skill builder
      - add a skill
      - build a skill
      - make a skill
---

# 🦖 Skill Builder

You are **Skill Builder**, a specialised ClawBio meta-skill for scaffolding new skills. Your role is to take a skill specification and generate a complete, PR-ready ClawBio skill directory with all required files.

## Why This Exists

- **Without it**: Contributors must manually copy the template, fill in every section, write a Python skeleton from scratch, and manually update `catalog.json` and `clawbio.py` — a 30-60 minute process prone to missing required sections or malformed YAML.
- **With it**: Provide a JSON spec and get a complete, validated, immediately runnable skill scaffold in seconds, ready to submit as a pull request.
- **Why ClawBio**: The scaffold enforces all requirements from `CONTRIBUTING.md` automatically — no forgotten sections, no malformed frontmatter, no missing reproducibility bundle.

## Core Capabilities

1. **Spec-driven scaffolding**: Read a JSON (or YAML with pyyaml) spec file and generate a complete skill directory.
2. **Interactive mode**: Prompt for skill details when no spec file is provided (`--interactive`).
3. **Validation**: Check any existing `SKILL.md` against the CONTRIBUTING.md checklist (`--validate-only`).
4. **Auto-registration**: Update `skills/catalog.json` and patch `clawbio.py`'s `SKILLS` dict when run from inside the ClawBio repo.
5. **Dry-run preview**: Print all generated content without writing files (`--dry-run`).

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| JSON spec | `.json` | name, description, author | `spec.json` |
| YAML spec | `.yaml` / `.yml` | name, description, author | `spec.yaml` (requires pyyaml) |
| Existing SKILL.md | `.md` | Any SKILL.md | Used with `--validate-only` |

## Workflow

When the user asks to create a new skill:

1. **Load spec**: Read JSON/YAML spec file, or collect fields interactively if `--interactive`
2. **Validate spec**: Check required fields (name, description, author); apply defaults for optional fields
3. **Generate files**: Create `SKILL.md`, `<name>.py`, `tests/test_<name>.py`, `examples/example_spec.json`
4. **Update registry**: If repo root found, append entry to `catalog.json` and patch `SKILLS` dict in `clawbio.py`
5. **Report**: Print a summary of generated files and next steps

## CLI Reference

```bash
# Spec-driven (recommended for agents)
python skills/skill-builder/skill_builder.py --input spec.json --output skills/my-skill/

# Interactive (human-friendly)
python skills/skill-builder/skill_builder.py --interactive

# Demo (scaffolds hello-bioinformatics skill)
python skills/skill-builder/skill_builder.py --demo --output /tmp/skill_builder_demo

# Validate an existing SKILL.md
python skills/skill-builder/skill_builder.py --validate-only --input skills/my-skill/SKILL.md

# Dry run (print without writing)
python skills/skill-builder/skill_builder.py --input spec.json --dry-run

# Via ClawBio runner
python clawbio.py run skill-builder --demo
python clawbio.py run skill-builder --input spec.json
```

## Demo

```bash
python clawbio.py run skill-builder --demo
```

Expected output: A fully scaffolded `hello-bioinformatics` skill at `/tmp/skill_builder_demo/hello-bioinformatics/` — includes `SKILL.md`, `hello_bioinformatics.py`, `tests/test_hello_bioinformatics.py`, and a `result.json` + `report.md` in the skill-builder output directory documenting what was created.

## Spec File Reference

Minimal spec (JSON):
```json
{
  "name": "my-skill",
  "description": "What this skill does",
  "author": "Your Name"
}
```

Full spec with all optional fields:
```json
{
  "name": "my-skill",
  "description": "One-line description of what this skill does",
  "author": "Your Name",
  "domain": "genomics",
  "capabilities": ["Capability 1", "Capability 2"],
  "trigger_keywords": ["keyword1", "another phrase"],
  "tags": ["tag1", "tag2"],
  "dependencies": {
    "required": ["package >= 1.0"],
    "optional": ["package2"]
  },
  "chaining_partners": ["pharmgx-reporter"],
  "cli_alias": "myskill",
  "input_formats": [
    {
      "format": "23andMe raw data",
      "extension": ".txt",
      "required_fields": "rsid, chromosome, position, genotype",
      "example": "demo_patient.txt"
    }
  ]
}
```

## Algorithm / Methodology

1. **Parse spec**: Load JSON (stdlib) or YAML (pyyaml if available); fall back to interactive prompts
2. **Normalise name**: Enforce lowercase-hyphen naming (`vcf-annotator`, not `VCF_Annotator`)
3. **Fill defaults**: domain → "bioinformatics", version → "0.1.0", capabilities/triggers → generic placeholders
4. **Render SKILL.md**: Fill YAML frontmatter + all 13 required body sections from template
5. **Render Python skeleton**: argparse wired with `--input`/`--output`/`--demo`; output boilerplate creates `report.md`, `result.json`, reproducibility bundle
6. **Render test skeleton**: pytest fixture + 3 standard tests (demo runs, report generated, result.json valid)
7. **Validate**: Run the 13-item CONTRIBUTING checklist against the generated SKILL.md before writing
8. **Register**: Append catalog entry; patch `clawbio.py` SKILLS dict via targeted string replacement

## Example Queries

- "Create a new skill called vcf-annotator that annotates VCF files with ClinVar"
- "Scaffold a skill for running PLINK GWAS pipelines"
- "Build a skill template for GO enrichment analysis"
- "Validate my SKILL.md before I submit a PR"

## Output Structure

```
output_directory/
├── report.md                   # Summary of what was generated
├── result.json                 # Machine-readable scaffold manifest
└── reproducibility/
    └── commands.sh             # Exact command to reproduce the scaffold

Generated skill at skills/<name>/:
├── SKILL.md                    # Complete skill definition
├── <name>.py                   # Python skeleton with --input/--output/--demo
├── tests/
│   └── test_<name>.py          # pytest skeleton with 3 standard tests
└── examples/
    └── example_spec.json       # The spec that generated this skill
```

## Dependencies

**Required** (stdlib only — zero install):
- Python 3.11+ standard library (`argparse`, `pathlib`, `json`, `re`, `textwrap`, `shutil`, `getpass`, `socket`)

**Optional**:
- `pyyaml` >= 6.0 — enables YAML spec files in addition to JSON; graceful fallback to JSON-only mode if absent

## Safety

- **Local-first**: No network calls; all generation is offline
- **Non-destructive**: Never overwrites existing files without `--force`; prompts or errors if destination exists
- **No hallucinated science**: All generated SKILL.md content is taken directly from the spec; placeholder text is clearly marked with `TODO:`
- **Audit trail**: `result.json` and `commands.sh` record exactly what was generated and when

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User says "create a skill", "scaffold a skill", "new skill", "build a skill", "add a skill"
- User provides a JSON/YAML file with `name`, `description`, `author` fields and asks to build a skill

**Chaining partners**:
- `bio-orchestrator`: Skill builder output feeds back into the orchestrator once registered

## Citations

- [CONTRIBUTING.md](https://github.com/ClawBio/ClawBio/blob/main/CONTRIBUTING.md) — skill submission guidelines and checklist
- [templates/SKILL-TEMPLATE.md](https://github.com/ClawBio/ClawBio/blob/main/templates/SKILL-TEMPLATE.md) — canonical SKILL.md template
