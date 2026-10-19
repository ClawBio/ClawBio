---
name: methylation-clock
description: Compute epigenetic age from DNA methylation arrays using PyAging clocks from GEO accessions or local files.
license: MIT
metadata:
  version: 0.1.0
  tags:
  - epigenetics
  - methylation
  - aging
  - clock
  - pyaging
  - GEO
  - illlumina-450k
  - EPIC
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🕰️
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: pandas
    - kind: pip
      package: numpy
    - kind: pip
      package: matplotlib
    - kind: pip
      package: pyaging
    trigger_keywords:
    - epigenetic age
    - methylation clock
    - pyaging
    - Horvath
    - GrimAge
    - DunedinPACE
    - GEO
    - GSE
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
- One row per sample and numeric methylation beta values in `[0, 1]` under `cg`, `ch`, `rs`, or `nv` probe IDs. Missing probe values are allowed and passed to PyAging's selected imputer. Illumina defines beta values on `[0, 1]`; if your export reports methylation percentages on `[0, 100]`, divide the probe values by 100 before running this skill. Do not divide age or other metadata. The skill never guesses the input scale.
- An optional `sample_id` column in CSV/TSV files becomes the sample index and must be present, nonblank, and unique. Otherwise CSV/TSV rows receive numeric IDs. Pickle files retain their DataFrame index, which must also have present, nonblank, unique identifiers.
- Known metadata columns are `gender`, `sex`, `tissue_type`, and `dataset`. List any other nonnumeric columns with `--metadata-cols`; do not list a probe column as metadata. For PCGrimAge or GrimAge2, supply numeric `age` in years and known `gender`/`sex` (`F`/`M`, `Female`/`Male`) or a numeric `female` column (1/0). Do not list `age` or `female` in `--metadata-cols`: PyAging needs them as model features. The skill rejects missing inputs for those two clocks rather than imputing demographic covariates.
- Required output directory via `--output`
- Optional clock list via `--clocks`

### Demo And Usage

The bundled two-sample demo is synthetic. Its fixed-seed generator, public probe-identifier source, and checksums are documented in `skills/methylation-clock/data/PROVENANCE.md`. Its placeholder ages and sex labels have no relationship to the simulated beta values; use its clock predictions only to check that the workflow runs.

Install the optional PyAging dependency in the Python environment used to run ClawBio. The demo below was verified with PyAging 0.5.7; first use of a clock downloads its public model weights. Local input data is not uploaded.

```bash
python -m pip install "pyaging>=0.5.7,<0.6"
```

```bash
# Demo
python skills/methylation-clock/methylation_clock.py \
  --demo --clocks Horvath2013 \
  --output /tmp/methylation_clock_demo

# GEO input
python skills/methylation-clock/methylation_clock.py \
  --geo-id GSE139307 \
  --output /tmp/methylation_clock_geo

# Local methylation file
python skills/methylation-clock/methylation_clock.py \
  --input my_methylation.csv.gz \
  --clocks Horvath2013,DunedinPACE \
  --output /tmp/methylation_clock_local
```

Use pickle input only when you trust its source; deserializing an untrusted pickle can execute code. `--geo-id` downloads a public GEO example into the chosen output directory. The skill refuses a nonempty output directory.

### What The Five Default Clocks Measure

| Clock | PyAging target | Intended tissue | Additional model inputs | Output unit |
|---|---|---|---|---|
| Horvath2013 | Chronological age | Multiple tissues | None | Years |
| AltumAge | Chronological age | Multiple tissues | None | Years |
| PCGrimAge | Mortality-risk proxy reproducing GrimAge | Whole blood | Age and sex | Years |
| GrimAge2 | Mortality risk | Whole blood | Age and sex | Years |
| DunedinPACE | Pace of aging | Whole blood | None | Biological years per chronological year |

The targets, tissues, and units come from the [PyAging clock catalogue](https://pyaging.readthedocs.io/en/latest/clock_glossary.html); its [DNA methylation tutorial](https://pyaging.readthedocs.io/en/latest/tutorials/tutorial_dnam.html) identifies age and sex as features for PCGrimAge and GrimAge2. `tissue_type` is descriptive metadata, not an automatic compatibility check. Compare predictions only within an appropriate tissue and target; a DunedinPACE value is not an age in years. The report reads clock units from PyAging metadata. A difference between an estimate and chronological age is not, by itself, a validated health assessment or an age-acceleration analysis.

### Output Structure

```
methylation_clock_report/
├── report.md
├── result.json
├── figures/
│   ├── clock_distributions.png
│   └── clock_correlation.png  # at least two clocks and ten complete samples
├── tables/
│   ├── predictions.csv
│   ├── prediction_summary.csv
│   ├── missing_features.csv
│   ├── clock_metadata.json
│   └── input_qc.json
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Safety Rules

1. ClawBio is local-first: user methylation data must remain on-device.
2. The skill refuses non-empty output directories to avoid silent overwrite.
3. Reports must include this disclaimer: "ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."
4. Check `tables/input_qc.json` and `tables/missing_features.csv` before interpreting a prediction. Missing-feature counts and successful execution do not establish validity for a tissue or cohort.

## Agent Boundary

1. Route methylation clock requests to `skills/methylation-clock/methylation_clock.py`.
2. Do not infer clinical diagnosis or treatment from clock estimates.
3. Trigger terms include: epigenetic age, methylation clock, Horvath, GrimAge, DunedinPACE, GEO, GSE.
4. Valid downstream chaining: `rnaseq-de` for transcriptomic-aging contrasts and `equity-scorer` for cohort context.
