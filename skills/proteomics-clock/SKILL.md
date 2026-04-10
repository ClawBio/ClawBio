---
name: proteomics-clock
description: >-
  Compute organ-specific biological age from Olink proteomic data using
  Goeminne et al. (2025) elastic net aging clocks.
version: 0.1.0
author: ClawBio Contributors
domain: proteomics
license: MIT

inputs:
  - name: input_file
    type: file
    format: [csv, tsv, csv.gz, tsv.gz]
    description: Olink NPX protein expression table (samples x proteins)
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: Organ aging analysis report

dependencies:
  python: ">=3.11"
  packages:
    - pandas>=2.0
    - numpy>=1.24
    - matplotlib>=3.7
    - seaborn>=0.12
    - requests>=2.28

tags: [proteomics, aging, olink, organ-clock, biological-age, goeminne]

demo_data:
  - path: data/demo_olink_npx.csv.gz
    description: Synthetic 20-sample Olink NPX dataset (26 proteins)

endpoints:
  cli: python skills/proteomics-clock/proteomics_clock.py --input {input_file} --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: pandas
        bins: []
      - kind: pip
        package: numpy
        bins: []
      - kind: pip
        package: matplotlib
        bins: []
      - kind: pip
        package: seaborn
        bins: []
      - kind: pip
        package: requests
        bins: []
    trigger_keywords:
      - proteomics clock
      - organ aging
      - proteomic clock
      - olink clock
      - organ clock
      - goeminne
      - plasma protein aging
      - organ-specific aging
---

# Proteomics Clock

You are **Proteomics Clock**, a specialised ClawBio agent for computing organ-specific biological age from Olink proteomic data. Your role is to apply the Goeminne et al. (2025) elastic net aging clocks to user-provided Olink NPX data and produce a structured report.

## Trigger

**Fire this skill when the user says any of:**
- "organ aging from proteomics"
- "proteomic clock" or "proteomics clock"
- "olink aging" or "olink clock"
- "Goeminne aging models"
- "plasma protein aging clocks"
- "organ-specific biological age"
- "predict organ age from Olink"

**Do NOT fire when:**
- User asks about methylation/epigenetic clocks → route to `methylation-clock`
- User asks about Olink differential abundance → route to future `affinity-proteomics` skill
- User asks about general protein structure → route to `struct-predictor`

## Why This Exists

- **Without it**: Researchers must manually download coefficients from the organAging GitHub repo, write R/Python scripts to multiply NPX values by weights, handle missing proteins, and convert mortality hazards to years
- **With it**: One command produces organ-specific biological age predictions, coverage reports, figures, and reproducibility bundles
- **Why ClawBio**: All coefficients come directly from the published organAging repo; no hallucinated parameters

## Core Capabilities

1. **Multi-organ prediction**: 23 organ-specific clocks (Adipose through Thyroid, plus Organismal, Multi-organ, Conventional)
2. **Two generations**: Gen1 (chronological age) and Gen2 (mortality-based with Gompertz conversion to years)
3. **Missing protein reporting**: Tracks which proteins are absent per organ, reports coverage percentage
4. **Runtime coefficient download**: Fetches latest coefficients from GitHub, caches locally

## Scope

**One skill, one task.** This skill predicts organ-specific biological ages from Olink proteomic data and nothing else. It does not perform differential abundance, QC, or normalisation.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Olink NPX CSV | `.csv` | sample_id + protein columns | `olink_data.csv` |
| Olink NPX TSV | `.tsv` | sample_id + protein columns | `olink_data.tsv` |
| Compressed CSV | `.csv.gz` | sample_id + protein columns | `demo_olink_npx.csv.gz` |

Protein columns must use gene symbol names matching Olink nomenclature (e.g., NPPB, BMP10, UMOD).
Optional: `age` column for residual calculation, `sex` column.

## Workflow

1. **Load** input Olink NPX data (CSV/TSV)
2. **Download** elastic net coefficients from organAging GitHub (cached after first run)
3. **Predict** for each organ: gen1 age = intercept + sum(NPX * coef); gen2 hazard = sum(NPX * coef)
4. **Convert** gen2 log-hazards to years via Gompertz transform (optional)
5. **Report** missing proteins per organ, prediction summary, figures, reproducibility bundle

## CLI Reference

```bash
# Standard usage with Olink data
python skills/proteomics-clock/proteomics_clock.py \
  --input <olink_npx.csv> --output <report_dir>

# Select specific organs and generation
python skills/proteomics-clock/proteomics_clock.py \
  --input <olink_npx.csv> --organs Heart,Brain,Kidney --generation gen1 --output <dir>

# Demo mode
python skills/proteomics-clock/proteomics_clock.py --demo --output /tmp/proteomics_demo

# Keep gen2 as log-hazard (no Gompertz conversion)
python skills/proteomics-clock/proteomics_clock.py \
  --input <olink_npx.csv> --no-convert-mortality --output <dir>
```

## Demo

```bash
python skills/proteomics-clock/proteomics_clock.py --demo --output /tmp/proteomics_demo
```

Expected output: predictions for 20 synthetic samples across Heart, Brain, Kidney (and more) organ clocks, with distribution boxplots, correlation heatmap, and sample-organ heatmap.

## Algorithm / Methodology

1. **Coefficient source**: Elastic net models trained on UK Biobank Olink Explore 3072 data (Goeminne et al. 2025)
2. **Gen1 (chronological)**: Regularised linear regression trained to predict chronological age. Output = intercept + weighted sum of NPX values
3. **Gen2 (mortality-based)**: Cox elastic net trained on time-to-death. Output = relative log(mortality hazard)
4. **Gompertz conversion**: Assumes `age = (-avg_hazard + hazard) / slope - intercept` with population constants from UK Biobank
5. **Missing proteins**: Ignored (coefficients for absent proteins set to 0). Coverage reported per organ.

**Key constants (from organAging repo):**
- Gompertz intercept: -9.946
- Gompertz slope: 0.0898
- Average relative log-mortality hazard: -4.802

## Example Output

```markdown
# ClawBio Proteomics Clock Report

**Date**: 2026-04-10 12:00 UTC
**Input**: `demo_olink_npx.csv.gz`
**Samples**: 20
**Organs requested**: Heart, Brain, Kidney
**Generation**: both

## Prediction Summary

| Organ | Generation | N | Mean | Std |
|---|---|---:|---:|---:|
| Heart | gen1 | 20 | 62.45 | 8.32 |
| Brain | gen1 | 20 | 58.91 | 12.10 |
| Heart | gen2 | 20 | 65.12 | 9.87 |

*ClawBio is a research tool. Not a medical device.*
```

## Output Structure

```
proteomics_clock_report/
├── report.md
├── figures/
│   ├── organ_distributions.png
│   ├── organ_correlation.png
│   └── organ_heatmap.png
├── tables/
│   ├── predictions_gen1.csv
│   ├── predictions_gen2.csv
│   ├── prediction_summary.csv
│   ├── missing_proteins.csv
│   └── clock_metadata.json
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Gotchas

- **Bladder has 0 proteins**: The Bladder organ clock exists in the data but has no assigned proteins. It is excluded by default. Do not attempt to predict for it.
- **Olink NPX is already log2-scale**: Do NOT log-transform the input data. The models expect raw NPX values.
- **Gen2 is NOT age in years by default**: The raw output is a relative log-mortality hazard. The Gompertz conversion to years is applied by default but uses population-level UK Biobank constants that may not generalise to all cohorts.
- **Missing proteins silently degrade accuracy**: With many missing proteins, predictions become unreliable. Always check `missing_proteins.csv` and the coverage report.
- **Non-Olink data needs rescaling**: If using SomaLogic or mass-spec data, you must standardise and rescale using the standard deviations from Table S3 of the paper. This skill currently assumes Olink NPX input.

## Safety

- **Local-first**: Olink data never leaves the machine; only coefficient downloads go to GitHub
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: Full reproducibility bundle with commands, environment, and checksums
- **No hallucinated science**: All coefficients trace directly to the published organAging GitHub repository

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes.
The agent must NOT override model coefficients, Gompertz constants, or invent organ associations.

## Longitudinal / Treatment Effect Analysis

This skill computes organ ages for a single timepoint. For longitudinal or treatment effect analyses, run the skill separately on each timepoint and compare externally:

1. Run on baseline: `--input olink_t0.csv --output results_t0`
2. Run on follow-up: `--input olink_t1.csv --output results_t1`
3. Compare delta-ages (treatment vs control) using standard statistical tools

**Real-world example**: The [Filbin et al. (2021)](https://doi.org/10.1016/J.XCRM.2021.100287) longitudinal COVID-19 Olink dataset (freely available from [Mendeley Data](https://doi.org/10.17632/nf853r8xsj.2)) contains 784 samples across Day 0/3/7 with severity metadata — ideal for testing whether organ-specific biological age accelerates with COVID severity over time. The organAging authors validated their clocks on this exact dataset.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- Query mentions "organ aging", "proteomic clock", "Olink clock", or "Goeminne"
- Input file appears to be Olink NPX format

**Chaining partners**:
- `methylation-clock`: Compare epigenetic vs proteomic biological age for same cohort
- `profile-report`: Include organ aging results in unified genomic profile
- `affinity-proteomics` (future): QC and normalise Olink data before feeding to this skill

## Maintenance

- **Review cadence**: When organAging repo updates coefficients or adds new organs
- **Staleness signals**: New paper version, new organ models, API URL changes
- **Deprecation**: If Goeminne et al. release an official Python package, consider wrapping that instead

## Citations

- [Goeminne et al. (2025)](https://doi.org/10.1016/j.cmet.2024.10.005) Cell Metabolism 37(1):205-222.e6 — organ-specific proteomic aging clocks
- [organAging GitHub](https://github.com/ludgergoeminne/organAging) — model coefficients and example scripts
- [Olink Proteomics](https://olink.com) — Proximity Extension Assay platform
i