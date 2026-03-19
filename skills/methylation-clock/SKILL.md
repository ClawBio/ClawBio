---
name: methylation-clock
description: Compute epigenetic age from DNA methylation arrays using PyAging clocks from GEO accessions or local files.
version: 0.1.0
tags: [epigenetics, methylation, aging, clock, pyaging, GEO, illlumina-450k, EPIC]
trigger_keywords: [epigenetic age, methylation clock, pyaging, horvath, grimage, dunedinpace, GEO, GSE]
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
        package: pandas
        bins: []
      - kind: pip
        package: numpy
        bins: []
      - kind: pip
        package: matplotlib
        bins: []
      - kind: pip
        package: pyaging
        bins: []
---

# Methylation Clock

## Domain Decisions

Epigenetic age workflows are difficult to reproduce because preprocessing and clock inputs differ across tools and publications.
This skill standardizes a PyAging-first pipeline from ingestion to report generation, with explicit reproducibility outputs.

### Core Capabilities

1. Accepts exactly one input source: GEO accession (`--geo-id`) or local methylation file (`--input`).
2. Applies notebook-aligned preprocessing (female derivation and EPICv2 aggregation by default).
3. Converts tabular data to AnnData and runs one or more methylation clocks.
4. Exports predictions, missing-feature diagnostics, metadata, figures, and reproducibility artifacts.

### Input Contract

- Exactly one input source:
  - GEO accession with `--geo-id` (example: `GSE139307`)
  - Local file with `--input` (`.pkl`, `.pickle`, `.csv`, `.tsv`, `.csv.gz`, `.tsv.gz`)
- Required output directory via `--output`
- Optional clock list via `--clocks`

### Demo And Usage

Demo fixture provenance and checksum are documented in `skills/methylation-clock/data/PROVENANCE.md`.

Install optional methylation-clock dependency (not part of the global base requirements):

```bash
pip install pyaging>=0.1
```

```bash
# Demo
python skills/methylation-clock/methylation_clock.py \
  --input skills/methylation-clock/data/GSE139307_small.csv.gz \
  --output /tmp/methylation_clock_demo

# GEO input
python skills/methylation-clock/methylation_clock.py \
  --geo-id GSE139307 \
  --output /tmp/methylation_clock_geo

# Local methylation file
python skills/methylation-clock/methylation_clock.py \
  --input my_methylation.pkl \
  --clocks Horvath2013,AltumAge,PCGrimAge,GrimAge2,DunedinPACE \
  --output /tmp/methylation_clock_local
```

### Output Structure

```
methylation_clock_report/
├── report.md
├── figures/
│   ├── clock_distributions.png
│   └── clock_correlation.png
├── tables/
│   ├── predictions.csv
│   ├── prediction_summary.csv
│   ├── missing_features.csv
│   └── clock_metadata.json
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Safety Rules

1. ClawBio is local-first: user methylation data must remain on-device.
2. The skill refuses non-empty output directories to avoid silent overwrite.
3. Reports must include this disclaimer: "ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."

## Agent Boundary

1. Route methylation clock requests to `skills/methylation-clock/methylation_clock.py`.
2. Do not infer clinical diagnosis or treatment from clock estimates.
3. Trigger terms include: epigenetic age, methylation clock, Horvath, GrimAge, DunedinPACE, GEO, GSE.
4. Valid downstream chaining: `rnaseq-de` for transcriptomic-aging contrasts and `equity-scorer` for cohort context.
