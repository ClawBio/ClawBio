---
name: omics-target-evidence-mapper
description: Aggregate public target-level evidence across omics and translational sources for research triage.
version: 0.1.0
tags:
  - omics
  - targets
  - translational-research
  - literature
  - trials
trigger_keywords:
  - target evidence
  - gene disease evidence
  - target triage
  - omics evidence
  - gene disease mapper
---

# Why This Exists

Researchers often need a quick first-pass view of whether a gene or protein target has evidence across multiple public sources. In practice, this usually means checking several websites manually and informally combining results. This skill makes that process reproducible by retrieving and organising public evidence into one structured output.

This skill is for research triage only. It does not infer causality, rank therapeutic value, or make clinical recommendations.

# Core Capabilities

1. Accept a gene or protein target and an optional disease term.
2. Retrieve canonical target information from UniProt.
3. Retrieve disease-target association evidence from Open Targets.
4. Retrieve recent literature hits from PubMed.
5. Optionally retrieve trial records relevant to the target and disease.
6. Produce a machine-readable JSON file and a human-readable Markdown report.

# Input Formats

| Argument | Required | Example | Notes |
|---|---|---|---|
| `--gene` | Yes, unless `--demo` is used | `IL6R` | Gene or target symbol |
| `--disease` | No | `coronary artery disease` | Optional disease context |
| `--output` | Yes | `demo_out` | Output directory |
| `--max-papers` | No | `5` | Number of PubMed hits to include |
| `--max-trials` | No | `5` | Number of trial records to include |
| `--demo` | No | `--demo` | Runs the built-in demo query |

# Workflow

1. Validate CLI inputs.
2. Resolve the query from either user input or demo mode.
3. Query public data sources.
4. Normalise results into a structured evidence object.
5. Write JSON and Markdown outputs.
6. Write file checksums for reproducibility.

# CLI Reference

## Demo mode

```bash
python skills/omics-target-evidence-mapper/omics_target_evidence_mapper.py --demo --output demo_out

