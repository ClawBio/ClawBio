---
name: gi-chromatin
description: Predict chromatin state — histone marks, DNase, TF binding — across 919 tracks (DeepSEA-style) for DNA sequences, via the hosted Genomic Intelligence /v1/tasks/chromatin/predict API.
license: MIT
metadata:
  openclaw:
    requires:
      bins:
      - python3
      env: null
      config: null
    always: false
    emoji: 🧶
    homepage: https://docs.genomicintelligence.ai
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: requests
      bins: null
    trigger_keywords:
    - chromatin
    - chromatin state
    - chromatin annotation
    - histone mark
    - histone modification
    - DNase
    - ATAC
    - TF binding
    - transcription factor binding
    - DeepSEA
    - epigenome
    - gi chromatin
    - genomic intelligence chromatin
  author: ClawBio + Genomic Intelligence
  demo_data:
  - path: example_data/chromatin_active_promoter_chr19.fa
    description: Chr19 active-promoter region — bundled real human reference sequence.
  dependencies:
    python: '>=3.10'
    packages:
    - requests>=2.31
  domain: genomics
  endpoints:
    cli: python skills/gi-chromatin/gi_chromatin.py --input {input_file} --output {output_dir}
  inputs:
  - name: input_file
    type: file
    format:
    - fa
    - fasta
    - fna
    description: Single-record FASTA, 200–500,000 bp (whitespace stripped); the API windows automatically.
    required: false
  outputs:
  - name: report
    type: file
    format: md
    description: Markdown report — windows processed, total annotations across tracks, model + timing.
  - name: result
    type: file
    format: json
    description: Full `{data, meta}` response with per-window per-track predictions.
  - name: reproducibility
    type: directory
    description: command.sh + environment.json.
  tags:
  - genomics
  - chromatin
  - histone
  - DNase
  - ATAC
  - TF-binding
  - deepsea
  - dna-lm
  - gi-api
  version: 0.1.0
---

# 🧶 gi-chromatin

You are **gi-chromatin**, a ClawBio agent that calls the **Genomic Intelligence** chromatin-annotation model (DeepSEA-style, 919 tracks: histone marks + DNase + TF binding across ENCODE cell types).

> ⚠️ **Remote inference — opt-in required.** Unlike most ClawBio skills, this skill uploads your FASTA sequence to the hosted Genomic Intelligence API at `https://api.genomicintelligence.ai`. Prefer a browser? The same models run interactively at <https://genomicintelligence.ai>. **Do not submit identifiable patient data** without an appropriate data-use agreement. Key setup: see [Authentication](#authentication) below.

## Trigger

**Fire this skill when the user says any of:**
- "predict chromatin state for this sequence"
- "histone mark prediction", "DNase prediction", "ATAC prediction"
- "TF binding prediction"
- "DeepSEA"
- "gi-chromatin", "predict epigenome"
- "is this region accessible?"

**Do NOT fire when:**
- The user asks specifically about enhancer activity → `gi-enhancer`
- The user asks for promoter prediction → `gi-promoter`

## Why This Exists

- **Without it**: Running DeepSEA / similar locally needs custom torch envs + weight wrangling.
- **With it**: One CLI call → 919 track predictions per window, in seconds.
- **Why ClawBio**: Hosted G0 DeepSEA inference plus ClawBio reproducibility and chaining.

## API Backed

`POST https://api.genomicintelligence.ai/v1/tasks/chromatin/predict` — default model `g0-deepsea` (919-track DeepSEA-style prediction head).

> **Contract note (2026-08-19).** The API publishes one operation per task, each with its own request schema — per-task `minLength`/`maxLength` on `sequence`, and a typed, closed `options` object (an unknown option key is a `422 validation_failed`, not a silent ignore). The URL above is byte-identical to what this skill already posted; only the published document changed. The bounds quoted in this file are live on the GI DEV service (`2026.08.19.4`) and reach `api.genomicintelligence.ai` at the next promotion — until then PROD is the more permissive of the two and this skill's local gate is the stricter. The authority is always the served schema: `GET https://api.genomicintelligence.ai/v1/openapi.json`.

## Workflow

1. **Parse**: single-record FASTA.
2. **POST** to `/v1/tasks/chromatin/predict`.
3. **Render**: `report.md` (window + total-annotation counts; per-track detail in `result.json`).

## CLI Reference

```bash
python skills/gi-chromatin/gi_chromatin.py --demo --output /tmp/gi-chromatin-demo
python skills/gi-chromatin/gi_chromatin.py --input my_region.fa --output report_dir
python clawbio.py run gi-chromatin --demo
```

## Authentication

The skill requires a Genomic Intelligence partner key in `GI_API_KEY`. Resolution order:

1. `--api-key <value>` CLI flag (explicit override).
2. `GI_API_KEY` environment variable.
3. Otherwise: the skill raises a `RuntimeError` pointing here.

### Quick start — ClawBio hackathon key

A shared hackathon-tier key ships in `.env.example` at the repo root (50 concurrent / 120 rpm, opt-in only). From wherever the ClawBio files live on your machine:

```bash
# Repo root (git clone) — or ~/.claude/plugins/cache/clawbio/clawbio/<version>/ for plugin installs
cp .env.example .env
set -a && source .env && set +a
```

### Production / heavier use

Request an individual key at **contact@genomicintelligence.ai**, then:

```bash
export GI_API_KEY=gi_yourkeyhere
```

## Demo

```bash
python clawbio.py run gi-chromatin --demo
```

Bundled fixture is an active-promoter region from chr19. Expect dense annotation across active-promoter tracks (H3K4me3, H3K27ac, DNase, etc.) and many called windows.

## Gotchas

- **Big response.** 919 tracks × N windows → multi-MB `result.json`. The report.md summarizes; mine `result.json` programmatically for specific tracks.
- **Track labels are in the response.** Don't hardcode track indices — read the names from `data.tracks`.
- **Length bounds are 200–500,000 bp**, published as `minLength` / `maxLength` on `ChromatinPredictRequest` and counted after whitespace is stripped. Both ends are a `422 validation_failed` (over-max is *not* a 413 — 413 is the separate 16 MiB raw-body cap). The skill rejects either locally before spending a request.
- **200 bp is admission control, not regime.** The model's context window is 1,000 bp (`bio_spec.context_window_bp` on `GET /v1/tasks/chromatin/models`), so 200–999 bp is accepted and scored — against a window padded out to 1,000 bp. The skill warns when you are under it.
- **Pre-windowing is unnecessary** — API strides internally.
- **Hackathon key is shared** — `GI_API_KEY` for heavier use.

## Output Structure

```
output_dir/
├── report.md
├── result.json
└── reproducibility/
    ├── command.sh
    └── environment.json
```

## Integration with Bio Orchestrator

Routes here on: "chromatin", "histone marks", "DNase", "ATAC", "TF binding", "DeepSEA".

Chains with: `gi-enhancer` (cross-validate enhancer calls against H3K27ac), `gi-promoter` (active-promoter signature: high H3K4me3 + DNase), `variant-annotation` (variants in accessible chromatin).

## Safety

Research tool. Not a clinical assay.
