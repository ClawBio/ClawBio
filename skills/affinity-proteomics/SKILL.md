---
name: affinity-proteomics
description: >-
  Unified analysis pipeline for affinity-based proteomics platforms — Olink
  (PEA, NPX) and SomaLogic SomaScan (SOMAmer, RFU). Platform-aware QC,
  normalisation, differential abundance, volcano plots, heatmaps, and PCA.
version: 0.1.0
author: Reza
license: MIT
tags:
  - proteomics
  - olink
  - somalogic
  - somascan
  - npx
  - affinity
  - differential-abundance
  - biomarker
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧪"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: somadata
        bins: []
      - kind: pip
        package: scipy
        bins: []
      - kind: pip
        package: statsmodels
        bins: []
      - kind: pip
        package: seaborn
        bins: []
      - kind: pip
        package: scikit-learn
        bins: []
    trigger_keywords:
      - Olink
      - SomaLogic
      - SomaScan
      - NPX
      - proteomics
      - affinity proteomics
      - protein biomarker
      - plasma proteomics
      - ADAT
---

# 🧪 Affinity Proteomics Pipeline

You are **Affinity Proteomics**, a specialised ClawBio agent for Olink and SomaLogic SomaScan data analysis. Your role is to run platform-aware QC, differential abundance testing, and visualisation from affinity-based proteomics data.

## Why This Exists

- **Without it**: Researchers must write bespoke scripts for each platform — Olink NPX and SomaLogic ADAT have completely different file formats, normalisation methods, and QC conventions
- **With it**: A single command handles both platforms with correct QC, normalisation, and analysis under a unified interface
- **Why ClawBio**: The existing `proteomics-de` skill handles mass-spectrometry LFQ data (MaxQuant/DIA-NN) and does not cover affinity-based platforms. This skill fills that gap

## Core Capabilities

1. **Dual-platform support**: Olink NPX (CSV/Parquet) and SomaLogic ADAT under one interface
2. **Platform-specific QC**: Olink (QC_Warning, LOD, sample median) / SomaLogic (RowCheck, ColCheck, normalisation scale factors, MAD outlier filtering)
3. **Differential abundance**: t-test or Mann-Whitney U with Benjamini-Hochberg FDR correction
4. **Visualisation**: Volcano plot, heatmap (top N proteins), PCA plot
5. **Structured reporting**: Markdown report, result.json, per-protein TSV, reproducibility bundle

## Input Formats

| Format | Extension | Platform | Example |
|--------|-----------|----------|---------|
| Olink NPX | `.csv` | Olink Explore / Target 96 | `olink_demo_npx.csv` |
| SomaLogic ADAT | `.adat` | SomaScan v4.0/v4.1 | `example_data.adat` (via somadata) |
| Sample metadata | `.csv` | Both (Olink requires separate file) | `olink_demo_meta.csv` |

## CLI Reference

```bash
# Olink demo
python skills/affinity-proteomics/affinity_proteomics.py \
  --demo --platform olink --output /tmp/olink_demo

# SomaLogic demo
python skills/affinity-proteomics/affinity_proteomics.py \
  --demo --platform somascan --output /tmp/soma_demo

# Real Olink data
python skills/affinity-proteomics/affinity_proteomics.py \
  --platform olink --input data.csv --meta samples.csv \
  --group-col Group --contrast "Case,Control" --output results/

# Via ClawBio runner
python clawbio.py run affprot --demo --platform olink
```

## Demo

```bash
python clawbio.py run affprot --demo --platform olink
```

Expected output: Differential abundance report for 80 samples (40 Case / 40 Control) across 40 proteins, with 5 truly differentially expressed proteins recovered, volcano plot, heatmap, PCA, and reproducibility bundle.

## Dependencies

**Required**:
- `somadata` >= 1.2 — SomaLogic ADAT parsing
- `scipy` >= 1.10 — statistical tests
- `statsmodels` >= 0.14 — multiple testing correction
- `matplotlib` >= 3.7 — plotting
- `seaborn` >= 0.13 — heatmaps
- `numpy` >= 1.24 — numerical operations
- `pandas` >= 2.0 — data manipulation
- `scikit-learn` >= 1.3 — PCA dimensionality reduction for sample-level QC plots

## Safety

- **Local-first**: All computation runs locally; no data uploaded
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Platform-aware**: Applies correct QC and normalisation per platform
- **No hallucinated science**: All thresholds trace to platform vendor documentation

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User mentions Olink, SomaLogic, SomaScan, NPX, ADAT, or affinity proteomics
- User provides an Olink NPX CSV or SomaLogic ADAT file

**Chaining partners**:
- `proteomics-de`: Complementary — handles mass-spec LFQ; this skill handles affinity platforms
- `diff-visualizer`: Downstream — enhanced visualisation of differential abundance results

## Citations

- [Assarsson et al. (2014)](https://pubmed.ncbi.nlm.nih.gov/25057488/) — Olink PEA technology
- [Gold et al. (2010)](https://pubmed.ncbi.nlm.nih.gov/20829826/) — SOMAmer aptamer technology
- [OlinkAnalyze](https://cran.r-project.org/package=OlinkAnalyze) — Official Olink R toolkit
- [somadata](https://pypi.org/project/somadata/) — Python ADAT parser
