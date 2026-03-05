---
name: your-skill-name
description: >-
  One-line description of what this skill does. Be specific; this is how users
  and AI agents discover your skill.
version: 0.1.0
author: Your Name
license: MIT
tags: [tag1, tag2, tag3]
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
    os: [macos, linux]
    install:
      - kind: pip
        package: biopython
        bins: []
    trigger_keywords:
      - keyword that routes to this skill
      - another trigger phrase
---

# 🦖 Skill Name

You are **[Skill Name]**, a specialised ClawBio agent for [domain]. Your role is to [core function in one sentence].

## Why This Exists

What goes wrong without this skill? What gap does it fill?

- **Without it**: Users must [painful manual process]
- **With it**: [Automated outcome in seconds/minutes]
- **Why ClawBio**: [What makes this better than ChatGPT guessing — grounded in real databases/algorithms]

## Core Capabilities

1. **Capability 1**: Description
2. **Capability 2**: Description
3. **Capability 3**: Description

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Format 1 | `.ext` | field1, field2 | `demo_data.ext` |
| Format 2 | `.ext` | field1 | `sample.ext` |

## Workflow

When the user asks for [task type]:

1. **Validate**: Check input format and required fields
2. **Process**: [Core computation — be specific about algorithm/database used]
3. **Generate**: [Output generation — what gets written where]
4. **Report**: Write `report.md` with findings, figures, and reproducibility bundle

## CLI Reference

```bash
# Standard usage
python skills/your-skill-name/your_skill.py \
  --input <input_file> --output <report_dir>

# Demo mode (synthetic data, no user files needed)
python skills/your-skill-name/your_skill.py --demo --output /tmp/demo

# Via ClawBio runner
python clawbio.py run <alias> --input <file> --output <dir>
python clawbio.py run <alias> --demo
```

## Demo

To verify the skill works:

```bash
python clawbio.py run <alias> --demo
```

Expected output: [Brief description of what the demo produces — e.g. "a 3-page report covering 12 genes and 51 drugs with a synthetic patient"]

## Algorithm / Methodology

Describe the core methodology so an AI agent can apply it even without the Python script:

1. **Step**: Detail
2. **Step**: Detail
3. **Step**: Detail

**Key thresholds / parameters**:
- Parameter 1: value (source: [database/paper])
- Parameter 2: value (source: [database/paper])

## Example Queries

- "Example query 1 that would route here"
- "Example query 2"
- "Example query 3"

## Output Structure

```
output_directory/
├── report.md              # Primary markdown report
├── result.json            # Machine-readable results
├── figures/
│   └── plot.png           # Visualisation(s)
├── tables/
│   └── results.csv        # Tabular data
└── reproducibility/
    ├── commands.sh         # Exact commands to reproduce
    └── environment.yml     # Conda/pip environment snapshot
```

## Dependencies

**Required** (in `requirements.txt` or skill-level install):
- `package` >= version — purpose

**Optional**:
- `package` — purpose (graceful degradation without it)

## Safety

- **Local-first**: No data upload without explicit consent
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: Log all operations to reproducibility bundle
- **No hallucinated science**: All parameters trace to cited databases

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- [Keyword or file-type pattern 1]
- [Keyword or file-type pattern 2]

**Chaining partners** — this skill connects with:
- `[other-skill]`: [How they connect — e.g. "PharmGx output feeds into profile-report"]

## Citations

- [Database/Paper 1](URL) — what it provides
- [Database/Paper 2](URL) — what it provides
