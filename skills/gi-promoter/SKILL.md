---
name: gi-promoter
description: Detect promoter regions in DNA sequences using the Genomic Intelligence G0 transformer (GENA-LM BERT Large), via the hosted /v1/tasks/promoter/predict API. Returns per-window promoter probabilities
  and called regions.
license: MIT
metadata:
  openclaw:
    requires:
      bins:
      - python3
      env: null
      config: null
    always: false
    emoji: 🧬
    homepage: https://docs.genomicintelligence.ai
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: requests
      bins: null
    trigger_keywords:
    - promoter
    - promoter prediction
    - predict promoter
    - find promoter
    - promoter region
    - TSS prediction
    - transcription start site
    - gi promoter
    - genomic intelligence promoter
    - G0 promoter
    - GENA-LM promoter
  author: ClawBio + Genomic Intelligence
  demo_data:
  - path: example_data/promoter_tp53.fa
    description: TP53 locus, gene-sense (chr17:7661779-7687546, GRCh38, 25.8 kbp; TP53 is minus-strand, so this is the reverse complement) — bundled real reference sequence.
  dependencies:
    python: '>=3.10'
    packages:
    - requests>=2.31
  domain: genomics
  endpoints:
    cli: python skills/gi-promoter/gi_promoter.py --input {input_file} --output {output_dir}
  inputs:
  - name: input_file
    type: file
    format:
    - fa
    - fasta
    - fna
    description: Single-record FASTA, 300–500,000 bp (whitespace stripped). The API windows automatically (default model uses 2000 bp context, 1000 bp stride).
    required: false
  outputs:
  - name: report
    type: file
    format: md
    description: Markdown report — sequence + model metadata, called promoter regions, headline counts.
  - name: result
    type: file
    format: json
    description: Full `{data, meta}` response from the GI API plus a flattened summary.
  - name: reproducibility
    type: directory
    description: command.sh + environment.json for exact-rerun reproducibility.
  tags:
  - genomics
  - promoter
  - transcription
  - regulatory
  - dna-lm
  - transformer
  - gi-api
  version: 0.1.0
---

# 🧬 gi-promoter

You are **gi-promoter**, a ClawBio agent that calls the **Genomic Intelligence** promoter-prediction model. Given a DNA sequence of 300–500,000 bp, it returns per-window promoter probabilities and called regions, all in a few hundred milliseconds via the hosted API.

> ⚠️ **Remote inference — opt-in required.** Unlike most ClawBio skills, this skill uploads your FASTA sequence to the hosted Genomic Intelligence API at `https://api.genomicintelligence.ai`. The same models also run interactively at <https://genomicintelligence.ai>. **Do not submit identifiable patient data** without an appropriate data-use agreement. Key setup: see [Authentication](#authentication) below.

## Trigger

**Fire this skill when the user says any of:**
- "predict promoters in this sequence"
- "find promoters in [gene/region]"
- "is this a promoter?"
- "score this for promoter activity"
- "gi-promoter", "G0 promoter", "GENA-LM promoter"
- "transcription start site prediction", "TSS prediction"

**Do NOT fire when:**
- The user asks for splice sites → `gi-splice`
- The user asks for enhancer activity → `gi-enhancer`
- The user asks for chromatin state → `gi-chromatin`
- The user asks for gene/transcript structure → `gi-annotation`

## Why This Exists

- **Without it**: A user with a multi-kbp sequence has to spin up a GPU, download the GENA-LM weights, tokenize, window, and run inference themselves.
- **With it**: One CLI call → annotated report in <1 s for typical sequences. The model is hosted; see [Authentication](#authentication) for key setup.
- **Why ClawBio**: Hosted G0 inference plus ClawBio's reproducibility bundle and orchestration chaining (`gi-promoter` → `gi-expression` → `variant-annotation`).

## API Backed

`POST https://api.genomicintelligence.ai/v1/tasks/promoter/predict`. Omit `model` and the API resolves the default — a GENA-LM BERT Large transformer with a 2000 bp context and a 1000 bp prediction window. Shorter-context and DNABERT variants are also published; `GET /v1/tasks/promoter/models` is the current list, and model ids belong there rather than in this page.

> **Contract note.** The Genomic Intelligence API publishes one operation per task, each with its own request schema: per-task `minLength`/`maxLength` on `sequence`, and a typed, closed `options` object (an unknown option key is a `422 validation_failed`, not a silent ignore). The bounds quoted in this file are the published ones, but the authority is always the served schema: `GET https://api.genomicintelligence.ai/v1/openapi.json`.

## Workflow

1. **Parse**: read single-record FASTA via the shared `clawbio.gi.gi_client.read_fasta` helper (uppercase; refuses multi-record input and any base outside `ACGTN`).
2. **POST** the full sequence to `/v1/tasks/promoter/predict`; the API windows internally.
3. **Render**: write `report.md` (summary + region table), `result.json` (full `{data, meta}` envelope), `reproducibility/`.

## CLI Reference

```bash
# Demo — bundled TP53 region
python skills/gi-promoter/gi_promoter.py --demo --output /tmp/gi-promoter-demo

# Your own FASTA
python skills/gi-promoter/gi_promoter.py --input my_region.fa --output report_dir

# Pick a specific model (ids come from GET /v1/tasks/promoter/models)
python skills/gi-promoter/gi_promoter.py --demo --model <model-id>

# Via ClawBio runner
python clawbio.py run gi-promoter --demo
```

## Demo

```bash
python clawbio.py run gi-promoter --demo
```

Bundled fixture is the TP53 locus (25.8 kbp, GRCh38, gene-sense). Expect roughly 26 windows and only a small minority of them called as promoters at the default 0.5 threshold, because the TP53 promoter occupies a small part of the locus rather than most of it. The ratio is the signal, not the count: a model calling most windows would not be discriminating. Read the counts from your own run.

## Authentication

The skill requires a Genomic Intelligence partner key in `GI_API_KEY`. Resolution order:

1. `--api-key <value>` CLI flag (explicit override).
2. `GI_API_KEY` environment variable.
3. Otherwise: the skill raises a `RuntimeError` pointing here.

### Quick start — ClawBio hackathon key

A shared hackathon-tier key ships in `.env.example` at the repo root (opt-in only). Caps are per-key and are not published as a fixed number — read `RateLimit-Limit` / `RateLimit-Remaining` on any response for the live allowance. From wherever the ClawBio files live on your machine:

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

## Gotchas

- **Length bounds are 300–500,000 bp**, published as `minLength` / `maxLength` on `PromoterPredictRequest` and counted after whitespace is stripped. Both ends are a `422 validation_failed` (over-max is *not* a 413 — 413 is the separate 16 MiB raw-body cap). The skill rejects either locally before spending a request.
- **300 bp is admission control, not regime.** A 400 bp sequence is accepted and scored, but the default model has a 2000 bp context window, so anything shorter is scored against a window padded out to 2000 bp. Compare your length against the model's `bio_spec.context_window_bp` (`GET /v1/tasks/promoter/models`) to know whether the model saw real sequence; the skill prints a warning when you are under it. The 300 bp-context models are in regime at the floor.
- **Do not pre-window the sequence yourself.** Submit the full region; the API windows and strides internally. Pre-windowing inflates rate-limit usage and gives identical results.
- **Strand matters — submit gene-sense.** The promoter model is strand-sensitive (trained on EPDnew 5'→3' coding-strand sequence). For minus-strand genes, reverse-complement to gene-sense before submission. On the bundled TP53 fixture, gene-sense calls several promoter windows above the default 0.5 threshold and the genomic strand calls **none** — the score collapses below threshold across the whole locus. The bundled TP53 fixture is already gene-sense.
- **An empty promoter result is weak evidence of a strand error — and this does not generalise.** Because the promoter score collapses below threshold on the wrong strand (above), an unexpectedly empty result is worth re-checking orientation. Do not carry that heuristic to other tasks: `gi-splice` returns a full set of high-confidence sites on the wrong strand, so there an empty result means no sites, never a strand error.
- **The hackathon key is shared.** If you hit `429`, you are sharing one key's caps with everyone else. Those caps are per-key and can be retuned server-side, so don't hardcode a number — `RateLimit-Limit` on any response is the live per-minute burst allowance and `Retry-After` on a `429` is the wait. Set `GI_API_KEY` to your own key for serious work.
- **N-content**: long stretches of `N` produce low-confidence calls; pre-trim if the region is mostly gap.

## Output Structure

```
output_dir/
├── report.md              # Headline counts, region table, model + timing
├── result.json            # Full {data, meta} envelope from the API
└── reproducibility/
    ├── command.sh         # Exact invocation
    └── environment.json   # API base, model, request_id, timestamp
```

## Integration with Bio Orchestrator

Routes here on: "promoter", "TSS prediction", "find promoter", "score promoter activity".

Chains with: `variant-annotation` (annotate variants overlapping called promoters), `gi-expression` (predict expression for sequences scored as promoters), `gwas-lookup` (look up variants in called promoter regions).

## Safety

Research and development use. Not for clinical or diagnostic decisions. Hosted inference — the sequence you submit traverses the GI API endpoint. Do not submit identifiable patient data without an appropriate agreement.
