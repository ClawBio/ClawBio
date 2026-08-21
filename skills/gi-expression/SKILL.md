---
name: gi-expression
description: Predict tissue / cell-type expression (log TPM + TPM) from a 9,198–500,000 bp TSS-centered DNA sequence (longer than one 9,198 bp window needs --tss-index) using the Genomic Intelligence G0 Expression model, via the hosted /v1/tasks/expression/predict
  API. The model is conditioned on a free-text cell-type / assay description.
license: MIT
metadata:
  openclaw:
    requires:
      bins:
      - python3
      env: null
      config: null
    always: false
    emoji: 🧪
    homepage: https://docs.genomicintelligence.ai
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: requests
      bins: null
    trigger_keywords:
    - expression prediction
    - predict expression
    - sequence to expression
    - TPM prediction
    - cell type expression
    - tissue expression
    - RNA-seq prediction
    - gi expression
    - G0 expression
    - genomic intelligence expression
  author: ClawBio + Genomic Intelligence
  demo_data:
  - path: example_data/expression_hbb_k562.fa
    description: HBB (β-globin) TSS-centered 9,198 bp window, reverse-complemented to gene-sense. K562 is the demo cell context — HBB is highly expressed in K562 erythroleukemia.
  dependencies:
    python: '>=3.10'
    packages:
    - requests>=2.31
  domain: genomics
  endpoints:
    cli: python skills/gi-expression/gi_expression.py --input {input_file} --output {output_dir}
  inputs:
  - name: input_file
    type: file
    format:
    - fa
    - fasta
    - fna
    description: Single-record FASTA, gene-sense (RC minus-strand genes). Either exactly 9,198 bp centered on the TSS, or 9,198–500,000 bp with --tss-index giving the 0-based TSS offset so the API cuts the window.
    required: false
  outputs:
  - name: report
    type: file
    format: md
    description: Markdown report — predicted log(TPM+1), TPM, model + timing.
  - name: result
    type: file
    format: json
    description: Full `{data, meta}` response.
  - name: reproducibility
    type: directory
    description: command.sh + environment.json.
  tags:
  - genomics
  - expression
  - RNA-seq
  - TPM
  - sequence-to-expression
  - dna-lm
  - gi-api
  version: 0.1.0
---

# 🧪 gi-expression

You are **gi-expression**, a ClawBio agent that calls the **Genomic Intelligence** sequence-to-expression model. Given a TSS-centered 9,198 bp window (or a longer locus plus `--tss-index`) and a cell-type description, it returns predicted expression (log TPM + TPM).

> ⚠️ **Remote inference — opt-in required.** Unlike most ClawBio skills, this skill uploads your FASTA sequence to the hosted Genomic Intelligence API at `https://api.genomicintelligence.ai`. The same models also run interactively at <https://genomicintelligence.ai>. **Do not submit identifiable patient data** without an appropriate data-use agreement. Key setup: see [Authentication](#authentication) below.

## Trigger

**Fire this skill when the user says any of:**
- "predict expression for this gene / sequence"
- "what's the expression of this region in [cell type]?"
- "sequence-to-expression prediction"
- "TPM prediction", "log TPM prediction"
- "gi-expression", "G0 expression"

**Do NOT fire when:**
- The user has counts / RNA-seq output and wants differential expression → `rnaseq-de`
- The user wants tissue annotation / GTEx lookup → use external resources

## Why This Exists

- **Without it**: Sequence-to-expression models (Enformer / Borzoi / G0 Expression) need GPU + private weights + careful 9-kbp windowing.
- **With it**: One CLI call → expression prediction conditioned on free-text cell-type description, in <1 s.
- **Why ClawBio**: Private weights, hosted. ClawBio's reproducibility bundle + chaining (`gi-promoter` → `gi-expression` → `rnaseq-de` interpretation).

## API Backed

`POST https://api.genomicintelligence.ai/v1/tasks/expression/predict`. Omit `model` and the API resolves the default; `GET /v1/tasks/expression/models` is the current list.

> **Contract note.** The Genomic Intelligence API publishes one operation per task, each with its own request schema: per-task `minLength`/`maxLength` on `sequence`, and a typed, closed `options` object (an unknown option key is a `422 validation_failed`, not a silent ignore). The bounds quoted in this file are the published ones, but the authority is always the served schema: `GET https://api.genomicintelligence.ai/v1/openapi.json`.

## Workflow

1. **Parse**: single-record FASTA, gene-sense. Either exactly 9,198 bp TSS-centered, or 9,198–500,000 bp with `--tss-index`. Anything else is rejected locally before the request is sent.
2. **Build options**: `{"description": "assay term name is polyA plus RNA-seq. biosample summary is Homo sapiens K562."}` by default; override via `--description "..."`.
3. **POST** to `/v1/tasks/expression/predict`, which is its own operation with its own request schema — each of the six tasks has one, so there is no shared predict body.
4. **Render**: `report.md` (headline log TPM plus the scored window the API actually used) + `result.json` + `reproducibility/`.

## CLI Reference

```bash
# Demo — HBB in K562
python skills/gi-expression/gi_expression.py --demo --output /tmp/gi-expression-demo

# Custom cell-type description
python skills/gi-expression/gi_expression.py \
  --input my_tss_window.fa \
  --description "assay term name is polyA plus RNA-seq. biosample summary is Homo sapiens liver." \
  --output report_dir

# Whole locus — the API cuts the 9,198 bp window around --tss-index
# (0-based offset into the sequence, counted after whitespace is stripped)
python skills/gi-expression/gi_expression.py \
  --input my_locus_50kb.fa --tss-index 24000 \
  --output report_dir

# Via ClawBio runner
python clawbio.py run gi-expression --demo
```

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

## Demo

```bash
python clawbio.py run gi-expression --demo
```

Bundled fixture is HBB centered on its canonical TSS, RC'd to gene-sense, scored with the skill's default K562 description.

> Read the predicted value from your own run rather than from this page. Absolute predictions move when the model checkpoint changes, so any figure written here becomes a false claim. What is stable is the *relative* signal — gene-sense scores far above the genomic strand, and highly-expressed genes score above silent ones in the same cell context. Do not build assertions on an absolute value read from documentation.

## Gotchas

- **9,198 bp is a floor, not a fixed size.** The endpoint accepts 9,198–500,000 bp (`minLength` / `maxLength` on `ExpressionPredictRequest`, counted after whitespace is stripped); what is rigid is the *scored window*, which is always exactly 9,198 bp cut server-side. Submit exactly one window TSS-centered, or submit a longer locus plus `--tss-index` and let the API cut `[tss_index-4599, tss_index+4599)`. Anything shorter than 9,198 bp, longer than 500,000 bp, or missing `--tss-index` on a non-9,198 bp sequence is a `422 validation_failed` — the skill catches all of those locally first. Over-max is a 422, *not* a 413; 413 is the separate 16 MiB raw-body cap. Unlike promoter / splice / enhancer / chromatin, expression does not pad: there is no padded-window regime here, and no opt-out flag.
- **A `tss_index` error reports at `loc: ["body"]`, never `body.tss_index`.** Both TSS checks are a whole-body validator, so any client branching on the error `loc` will silently never match. Match on `error.code` (`validation_failed`) and use `message` for display only — and read `error.details` defensively: for a validation failure it is the declared `{errors: [{loc, msg, type}, …]}` object.
- **A wrong `--tss-index` does not error — it lies.** Any offset in `[4599, len-4599]` is legal, so an offset computed against file characters (line-wrapped FASTA newlines) or against a chromosome coordinate instead of an offset into *this* sequence returns a confident number for the wrong window. Always check the "Scored window" line in `report.md`. The response echoes the applied window in two places, with identical values: `data.input.scored_window` (which also carries `submitted_sequence_length`, and is the pair this skill's report reads) and `meta.task_specific_counts.scored_window`. Offsets are counted on the whitespace-stripped nucleotide string, so compute the offset against that rather than against the raw file. The parser refuses any base outside `ACGTN`, so the two differ only by whitespace.
- **Gene-sense is mandatory.** Minus-strand genes need reverse-complementing. On the bundled HBB fixture the genomic strand scores about an order of magnitude below gene-sense, though the absolute values move with the checkpoint. The wrong strand returns a well-formed low number, not an error.
- **`description` wording changes the answer.** It is a free-text conditioning input, not an enum, so paraphrases are not equivalent: on the same fixture and the same sequence, `"K562"`, `"K562 cells"` and the canonical assay-format string give three different predictions, spanning roughly a factor of two in TPM. Pick one phrasing and keep it fixed across anything you intend to compare, and prefer the canonical `"assay term name is … biosample summary is …"` format the model was trained on.
- **`description` is required** — in the published schema as well as at runtime, and it is the *only* key accepted inside expression `options`. The model is conditioned on it; "assay term name is polyA plus RNA-seq. biosample summary is Homo sapiens [tissue]." is the canonical format.
- **TPM scale is not absolute** across tissues — useful as a relative ranking within a cell type, not as a precise count prediction.
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

Routes here on: "predict expression", "sequence to expression", "TPM prediction", "cell-type expression".

Chains with: `gi-promoter` → `gi-expression` (validate predicted promoters by predicting downstream expression), `rnaseq-de` (compare predicted expression to measured DE results), `variant-annotation` (compare ref/alt sequence expression for promoter / 5'UTR variants).

## Safety

Research and development use. Not for clinical or diagnostic decisions. Predictions are model outputs, not measurements.
