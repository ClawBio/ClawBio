---
name: organ-aging-studio
description: >-
  Interactive Goeminne proteomic aging clock with organ filters and per-protein
  contribution breakdown (protein NPX × coefficient). Agent- and demo-friendly.
license: MIT
metadata:
  version: 0.1.0
  author: ClawBio hackathon contributor
  domain: proteomics
  tags:
    - aging
    - longevity
    - proteomics
    - biological age
    - organ clock
    - Goeminne
    - interpretability
  inputs:
    - name: input_file
      type: file
      format:
        - csv
        - tsv
        - csv.gz
        - tsv.gz
      description: Olink NPX protein table (samples × proteins)
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Human-readable aging report with per-organ summary
    - name: result
      type: file
      format:
        - json
      description: Machine-readable predictions and protein contributions
  dependencies:
    python: ">=3.11"
    packages:
      - pandas>=2.0
      - numpy>=1.24
      - requests>=2.28
  demo_data:
    - path: ../proteomics-clock/data/demo_olink_npx.csv.gz
      description: Synthetic 20-sample Olink NPX demo (shared with proteomics-clock)
  endpoints:
    cli: python skills/organ-aging-studio/organ_aging_studio.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🕰️"
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
        package: requests
    trigger_keywords:
      - organ aging studio
      - proteomic clock breakdown
      - protein coefficient aging
      - Goeminne clock explain
      - which proteins drive organ age
---

# Organ Aging Studio

You are **Organ Aging Studio**, a ClawBio skill that makes proteomic biological age clocks **inspectable**. Every prediction decomposes into:

```text
predicted_age = intercept + Σ (protein_NPX × coefficient)
```

## Trigger

**Fire this skill when the user says any of:**
- "organ aging studio" or "explain my organ age"
- "which proteins drive biological age"
- "protein breakdown for Goeminne clock"
- "interactive proteomic aging" or "filter proteins by coefficient"

**Do NOT fire when:**
- User only wants batch predictions without breakdown → route to `proteomics-clock`
- User asks about methylation / DNAm clocks → route to `methylation-clock`
- User asks about differential abundance → route to `affinity-proteomics`

## Why This Exists

| Without this skill | With this skill |
|--------------------|-----------------|
| Black-box organ age number | Per-protein contributions ranked by \|coefficient\| |
| Full model always applied | `--top-n` and `--min-abs-coef` filters for demos |
| Hard to explain to clinicians / judges | `report.md` + `protein_contributions.csv` + JSON for agents |

Built on the same pinned [organAging](https://github.com/ludgergoeminne/organAging) coefficients as `proteomics-clock`. **No invented weights.**
Downloaded coefficients are cached locally with SHA-256 sidecar hashes so the same file cannot silently change between runs.

## Core Capabilities

1. **Multi-organ** — any organ supported by Goeminne et al. (2025); default demo set Heart, Brain, Liver, Immune, Organismal
2. **Gen1 / Gen2** — chronological age models or mortality hazard → years (Gompertz)
3. **Protein filters** — `--top-n`, `--min-abs-coef`, single `--sample-id`
4. **Structured outputs** — Markdown report, JSON, contribution table, `reproducibility/` bundle

## Scope

One skill, one task. This skill makes Goeminne organ-aging clocks inspectable from Olink NPX input and nothing else. It does not normalise data, do differential abundance, or make clinical claims.

## Workflow

1. **Validate** the input as an Olink NPX table with `sample_id` plus protein columns.
2. **Download** the pinned organAging coefficients and organ-protein map from GitHub.
3. **Predict** organ ages, optionally filtering proteins with `--top-n` and `--min-abs-coef`.
4. **Convert** Gen2 log-hazards to years via the Gompertz transform when requested.
5. **Write** `report.md`, `result.json`, `protein_contributions.csv`, and the `reproducibility/` bundle (`commands.sh`, `environment.yml`, `checksums.sha256`).

## Input Formats

| Format | Extension | Required columns |
|--------|-----------|------------------|
| Olink NPX CSV | `.csv` | `sample_id` + protein gene symbols |
| Olink NPX TSV | `.tsv` | same |
| Compressed | `.csv.gz` | same |

Optional: `age` (for delta = bio − chrono), `sex`.

## CLI Reference

```bash
# Demo — synthetic Olink data (no download)
python skills/organ-aging-studio/organ_aging_studio.py \
  --demo --output /tmp/studio

# One patient, Heart only, top 5 drivers
python skills/organ-aging-studio/organ_aging_studio.py \
  --input my_olink.csv.gz --output /tmp/studio \
  --organs Heart --sample-id PATIENT_001 --top-n 5

# All demo samples, multiple organs
python skills/organ-aging-studio/organ_aging_studio.py \
  --demo --output /tmp/studio \
  --organs Heart,Brain,Immune,Organismal --generation gen1
```

### Flags

| Flag | Default | Description |
|------|---------|-------------|
| `--demo` | off | Use bundled synthetic Olink table |
| `--organs` | Heart,Brain,Liver,Immune,Organismal | Comma-separated organ list |
| `--generation` | gen1 | `gen1` = years; `gen2` = hazard → years |
| `--sample-id` | all rows | Analyse one sample |
| `--top-n` | all present | Keep top N proteins by \|coef\| |
| `--min-abs-coef` | 0 | Drop small coefficients |

## Demo

```bash
cd ClawBio
uv sync
python skills/organ-aging-studio/organ_aging_studio.py \
  --demo --output /tmp/organ-aging-studio \
  --organs Heart,Brain,Immune,Organismal \
  --sample-id DEMO_000 --top-n 10
```

**Expected outputs** in `/tmp/organ-aging-studio/`:

| File | Contents |
|------|----------|
| `report.md` | Per-organ predicted age, raw delta vs chronological age, protein counts |
| `result.json` | Full nested JSON for agents |
| `tables/protein_contributions.csv` | Long-format NPX × coef × contribution |
| `reproducibility/` | Replay command, environment, checksums |

Example summary row (synthetic demo):

| Organ | Predicted age | Chronological | Raw delta |
|-------|---------------|---------------|-----------|
| Heart | ~67 yr | 66 yr | +1 yr |
| Brain | ~42 yr | 66 yr | −24 yr |

> Demo NPX is **synthetic** — do not use it to validate correlation with age. For real Olink data, see `data/PROVENANCE.md`.
> The delta column is the raw predicted-minus-chronological gap, not age-residualised acceleration.

## Output Structure

```
output_dir/
├── report.md
├── result.json
├── tables/
│   └── protein_contributions.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Gotchas

- **Olink NPX is already log2-scaled**: Do not log-transform the input again.
- **Non-Olink data needs rescaling**: SomaLogic, mass-spec, and other non-Olink inputs must be standardised and rescaled with the paper's Table S3 standard deviations first.
- **Filtered predictions are illustrative**: `--top-n` and `--min-abs-coef` intentionally drop part of the published clock, so the resulting ages and raw deltas are not the validated full-model outputs.
- **Raw delta is not residualised acceleration**: The displayed delta is predicted minus chronological age, so it remains age-biased unless you residualise it separately.
- **Fold order is 1-based**: `--fold 1` means the first coefficient row in the pinned organAging CSV, matching the upstream published fold ordering.

## Real-world data (download separately)

Large cohorts are **not** bundled. See [`data/PROVENANCE.md`](data/PROVENANCE.md) for:

- **Filbin COVID Olink** (real plasma) — Mendeley download + `proteomics-clock/examples/fetch_filbin.py`
- **GEO GSE40279** (blood methylation validation) — for `methylation-clock`, not this skill's input
- **GEO GSE259312** (paired Olink + methylation) — future cross-omics work

## Agent Boundary

- May select organs, filters, and **explain** contributions from `result.json`
- Must **not** invent coefficients or alter the formula
- Must state demo data is synthetic when using `--demo`
- Must refuse clinical diagnosis language

## Safety

- Educational / research use only — **not a medical device**
- Do not run on identifiable patient data without consent
- Do not extrapolate beyond populations represented in clock training (UK Biobank–based models)

## Tests

```bash
pytest skills/organ-aging-studio/tests/ -q
```

## Citation

Goeminne LJE et al. (2025). *Cell Metabolism* 37(1):205-222.e6. DOI: [10.1016/j.cmet.2024.10.005](https://doi.org/10.1016/j.cmet.2024.10.005)
