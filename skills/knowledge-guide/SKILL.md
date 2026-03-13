---
name: knowledge-guide
description: >-
  GTN-backed educational guide — explains bioinformatics concepts grounded
  in Galaxy Training Network tutorials. Supports free-text queries, structured
  lookups by topic/tool/concept, and report-embedded "Learn More" sections.
version: 0.1.0
metadata:
  openclaw:
    category: education
    data_source: Galaxy Training Network (training.galaxyproject.org)
---

# Knowledge Guide

Explains bioinformatics concepts by grounding every answer in Galaxy Training
Network (GTN) tutorials. No hallucinated biology — every explanation traces
back to an authoritative training resource.

## Core Capabilities

1. **Free-text search** — "what is variant calling?" → ranked tutorial matches
2. **Topic lookup** — browse all tutorials within a GTN topic
3. **Tool lookup** — find tutorials that use a specific Galaxy tool
4. **Concept lookup** — fuzzy match against skill-seeded concept lists
5. **Skill Learn More** — pre-computed tutorial recommendations per ClawBio skill
6. **Deep content pull** — live-fetch full tutorial content for inline explanation

## Input Modes

| Flag | Example | Description |
|------|---------|-------------|
| `--query` | `"what is variant calling?"` | Free-text keyword search |
| `--topic` | `variant-analysis` | Direct GTN topic lookup |
| `--tool` | `fastqc` | Reverse tool→tutorial index |
| `--concept` | `"polygenic risk"` | Fuzzy concept matching |
| `--skill` | `gwas-prs` | ClawBio skill Learn More |
| `--deep` | (flag) | Fetch full tutorial content live |

## Scoring Weights (Free-Text Mode)

| Signal | Weight | Bonus |
|--------|--------|-------|
| Title match | 4.0/token | +10.0 phrase |
| Objectives match | 3.0/token | +6.0 phrase |
| Topic match | 2.0/token | — |
| Tool match | 1.0/token | — |

## Output Structure

```
output_dir/
├── report.md
├── report.html
├── result.json
└── tutorials/        (only with --deep)
    └── <name>.md
```

## Dependencies

- `requests` (for live GTN API calls and --deep mode)

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device
and does not provide clinical diagnoses. Consult a healthcare professional
before making any medical decisions.*
