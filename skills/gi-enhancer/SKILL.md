---
name: gi-enhancer
description: Predict enhancer activity in DNA sequences using the Genomic Intelligence G0 DeepSTARR model, via the hosted /v1/tasks/enhancer/predict API. Returns per-window activity scores.
license: MIT
metadata:
  openclaw:
    requires:
      bins:
      - python3
      env: null
      config: null
    always: false
    emoji: 🎚️
    homepage: https://docs.genomicintelligence.ai
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: requests
      bins: null
    trigger_keywords:
    - enhancer
    - enhancer activity
    - predict enhancer
    - regulatory element
    - cis-regulatory
    - CRE
    - DeepSTARR
    - STARR-seq
    - massively parallel reporter assay
    - MPRA
    - gi enhancer
    - genomic intelligence enhancer
  author: ClawBio + Genomic Intelligence
  demo_data:
  - path: example_data/enhancer_eve.fa
    description: Drosophila eve (even-skipped) developmental-enhancer region (chr2R:9972000-9982000, BDGP6, gene-sense, incl. upstream stripe enhancers) — canonical DeepSTARR benchmark.
  dependencies:
    python: '>=3.10'
    packages:
    - requests>=2.31
  domain: genomics
  endpoints:
    cli: python skills/gi-enhancer/gi_enhancer.py --input {input_file} --output {output_dir}
  inputs:
  - name: input_file
    type: file
    format:
    - fa
    - fasta
    - fna
    description: Single-record FASTA, 50–500,000 bp (whitespace stripped); the API windows automatically.
    required: false
  outputs:
  - name: report
    type: file
    format: md
    description: Markdown report — windows processed, max predicted activity, model + timing.
  - name: result
    type: file
    format: json
    description: Full `{data, meta}` response.
  - name: reproducibility
    type: directory
    description: command.sh + environment.json.
  tags:
  - genomics
  - enhancer
  - regulatory
  - cis-regulatory
  - deepstarr
  - dna-lm
  - gi-api
  version: 0.1.0
---

# 🎚️ gi-enhancer

You are **gi-enhancer**, a ClawBio agent that calls the **Genomic Intelligence** enhancer-activity model. Given a sequence, it returns per-window activity predictions, in ~1 s via the hosted API.

> ⚠️ **Remote inference — opt-in required.** Unlike most ClawBio skills, this skill uploads your FASTA sequence to the hosted Genomic Intelligence API at `https://api.genomicintelligence.ai`. The same models also run interactively at <https://genomicintelligence.ai>. **Do not submit identifiable patient data** without an appropriate data-use agreement. Key setup: see [Authentication](#authentication) below.

## Trigger

**Fire this skill when the user says any of:**
- "predict enhancer activity"
- "score this for enhancer / CRE / regulatory function"
- "is this an enhancer?"
- "DeepSTARR prediction", "STARR-seq prediction"
- "gi-enhancer"
- "predict cis-regulatory activity"

**Do NOT fire when:**
- The user asks for promoter activity → `gi-promoter`
- The user asks for chromatin state / accessibility → `gi-chromatin`

## Why This Exists

- **Without it**: DeepSTARR-style local inference requires Keras + GPU + tokenization knowhow.
- **With it**: One CLI call → per-window activity scores in ~1 s.
- **Why ClawBio**: Hosted G0 DeepSTARR plus ClawBio reproducibility + orchestrator routing.

## API Backed

`POST https://api.genomicintelligence.ai/v1/tasks/enhancer/predict`. Omit `model` and the API resolves the default — a DeepSTARR model trained on *Drosophila* S2 cells. `GET /v1/tasks/enhancer/models` is the current list.

> **Contract note.** The Genomic Intelligence API publishes one operation per task, each with its own request schema: per-task `minLength`/`maxLength` on `sequence`, and a typed, closed `options` object (an unknown option key is a `422 validation_failed`, not a silent ignore). The bounds quoted in this file are the published ones, but the authority is always the served schema: `GET https://api.genomicintelligence.ai/v1/openapi.json`.

## Workflow

1. **Parse**: single-record FASTA.
2. **POST** to `/v1/tasks/enhancer/predict`; the API windows internally.
3. **Render**: `report.md` + `result.json` + `reproducibility/`.

## CLI Reference

```bash
python skills/gi-enhancer/gi_enhancer.py --demo --output /tmp/gi-enhancer-demo
python skills/gi-enhancer/gi_enhancer.py --input my_region.fa --output report_dir
python clawbio.py run gi-enhancer --demo
```

## Authentication

The skill requires a Genomic Intelligence partner key in `GI_API_KEY`. Resolution order:

1. `--api-key <value>` CLI flag (explicit override).
2. `GI_API_KEY` environment variable.
3. Otherwise: the skill raises a `RuntimeError` pointing here.

### Quick start — ClawBio hackathon key

A shared hackathon-tier key ships in `.env.example` at the repo root (opt-in only). Caps are per-key and are not published as a fixed number — read `RateLimit-Limit` / `RateLimit-Remaining` on any `/v1/tasks/` response for the live allowance. The runner keeps them for you: they are in `result.json` under `rate_limit`, and a `429` names them on the error line. From wherever the ClawBio files live on your machine:

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
python clawbio.py run gi-enhancer --demo
```

Bundled fixture is the Drosophila *eve* (even-skipped) locus (chr2R:9972000-9982000, incl. the upstream stripe enhancers) — the canonical DeepSTARR benchmark for developmental enhancer activity. Expect a positive developmental signal; read the score from your own run.

## Gotchas

- **DeepSTARR was trained on Drosophila S2 cells.** Activity scores for mammalian sequences are still informative as a relative ranking, but the absolute scale is calibrated for fly chromatin.
- **Length bounds are 50–500,000 bp**, published as `minLength` / `maxLength` on `EnhancerPredictRequest` and counted after whitespace is stripped. Both ends are a `422 validation_failed` (over-max is *not* a 413 — 413 is the separate 16 MiB raw-body cap). The skill rejects either locally before spending a request.
- **50 bp is admission control, not a meaningful enhancer size.** It is the strictest floor any enhancer model needs; some models accept less. The models' context window is 249 bp, so 50–248 bp is accepted and scored — against a window padded out to 249 bp. Compare your length against `bio_spec.context_window_bp` (`GET /v1/tasks/enhancer/models`) to know whether the model saw real sequence; the skill warns when you are under it.
- **Pre-windowing is unnecessary** — the API windows and strides internally.
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

Routes here on: "enhancer", "DeepSTARR", "STARR-seq", "predict CRE", "regulatory activity".

Chains with: `gi-promoter` (joint regulatory-element scan), `gi-chromatin` (cross-validate with chromatin accessibility), `variant-annotation` (variants overlapping high-activity windows).

## Safety

Research and development use. Not for clinical or diagnostic decisions.
