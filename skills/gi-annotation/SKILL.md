---
name: gi-annotation
description: Predict gene and transcript structure (intervals, exons, strand) from a DNA sequence using the Genomic Intelligence DNA Annotation model, via the hosted /v1/tasks/annotation/predict API. Submitted asynchronously
  — the pipeline takes ~20 s for ~20 kbp.
license: MIT
metadata:
  openclaw:
    requires:
      bins:
      - python3
      env: null
      config: null
    always: false
    emoji: 📜
    homepage: https://docs.genomicintelligence.ai
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: requests
      bins: null
    trigger_keywords:
    - gene annotation
    - transcript annotation
    - annotate sequence
    - gene structure prediction
    - predict transcripts
    - de novo gene prediction
    - DNA annotation
    - gene boundaries
    - exon prediction
    - gi annotation
    - genomic intelligence annotation
  author: ClawBio + Genomic Intelligence
  demo_data:
  - path: example_data/annotation_tp53.fa
    description: TP53 locus (chr17:7668402-7687550, GRCh38, 19 kbp) — bundled real reference sequence.
  dependencies:
    python: '>=3.10'
    packages:
    - requests>=2.31
  domain: genomics
  endpoints:
    cli: python skills/gi-annotation/gi_annotation.py --input {input_file} --output {output_dir}
  inputs:
  - name: input_file
    type: file
    format:
    - fa
    - fasta
    - fna
    description: Single-record FASTA (genomic region), 1,000–500,000 bp (whitespace stripped) — tens to hundreds of kbp is the normal case.
    required: false
  outputs:
  - name: report
    type: file
    format: md
    description: Markdown report — predicted transcripts with start / end / strand.
  - name: result
    type: file
    format: json
    description: Full `{data, meta}` response with per-transcript structure.
  - name: reproducibility
    type: directory
    description: command.sh + environment.json.
  tags:
  - genomics
  - annotation
  - gene-prediction
  - transcript-prediction
  - gene-structure
  - dna-lm
  - gi-api
  version: 0.1.0
---

# 📜 gi-annotation

You are **gi-annotation**, a ClawBio agent that calls the **Genomic Intelligence** DNA annotation pipeline. Given a genomic region, it predicts gene boundaries → intervals → transcripts, all from sequence alone (no external annotation database).

> ⚠️ **Remote inference — opt-in required.** Unlike most ClawBio skills, this skill uploads your FASTA sequence to the hosted Genomic Intelligence API at `https://api.genomicintelligence.ai`. The same models also run interactively at <https://genomicintelligence.ai>. **Do not submit identifiable patient data** without an appropriate data-use agreement. Key setup: see [Authentication](#authentication) below.

## Trigger

**Fire this skill when the user says any of:**
- "annotate this DNA sequence"
- "predict genes / transcripts in this region"
- "what genes are encoded here?" (from sequence, not coordinates)
- "de novo gene prediction"
- "gi-annotation"

**Do NOT fire when:**
- The user has a VCF and wants variant consequences → `variant-annotation` (VEP)
- The user wants known gene records by coordinate → external NCBI / Ensembl lookup

## Why This Exists

- **Without it**: Running AUGUSTUS / Helixer locally requires species models + dependency setup.
- **With it**: One CLI call → predicted transcript structures, in ~20 s for ~20 kbp.
- **Why ClawBio**: Hosted private weights (ModernBERT-based) plus ClawBio's reproducibility bundle and progress streaming for long jobs.

## API Backed

`POST https://api.genomicintelligence.ai/v1/tasks/annotation/predict` with `Prefer: respond-async`. The API accepts either delivery mode on every task; this skill always submits async because the pipeline is long-running. The pipeline streams progress through `GET /v1/tasks/jobs/{job_id}` (typically: load → gene-boundaries → gene-intervals → transcripts).

> **Contract note.** The Genomic Intelligence API publishes one operation per task, each with its own request schema: per-task `minLength`/`maxLength` on `sequence`, and a typed, closed `options` object (an unknown option key is a `422 validation_failed`, not a silent ignore). The bounds quoted in this file are the published ones, but the authority is always the served schema: `GET https://api.genomicintelligence.ai/v1/openapi.json`.

## Workflow

1. **Parse**: single-record FASTA.
2. **Submit async**: `POST /v1/tasks/annotation/predict` with `Prefer: respond-async` → 202 + `job_id`.
3. **Poll**: stream progress (`percent`, `message`) until terminal.
4. **Render**: `report.md` (transcripts table) + `result.json` (full response) + `reproducibility/`.

## CLI Reference

```bash
# Demo — bundled TP53 region (~20 s)
python skills/gi-annotation/gi_annotation.py --demo --output /tmp/gi-annotation-demo

# Your own FASTA
python skills/gi-annotation/gi_annotation.py --input my_region.fa --output report_dir

# Via ClawBio runner
python clawbio.py run gi-annotation --demo
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
python clawbio.py run gi-annotation --demo
```

Bundled fixture is the TP53 locus (19 kbp). Expect several transcripts — TP53 has multiple annotated isoforms — and a wall time of roughly 20 s.

## Gotchas

- **Always submitted async.** This skill sends `Prefer: respond-async` and polls, so it never returns a synchronous response — though the API itself serves annotation synchronously when the header is omitted. Delivery mode is a per-request choice, not a property of the task.
- **Length bounds are 1,000–500,000 bp**, published as `minLength` / `maxLength` on `AnnotationPredictRequest` and counted after whitespace is stripped. Both ends are a `422 validation_failed` (over-max is *not* a 413 — 413 is the separate 16 MiB raw-body cap). The skill rejects either locally before spending a request. The 1,000 bp floor is the highest of the six tasks: the gene finder needs a region, not a single exon.
- **Long input is normal.** The model handles tens-to-hundreds of kbp up to the 500 kbp cap; longer regions take proportionally more time. Annotation reports `bio_spec.context_window_bp: null` — there is no sliding-window regime caveat here, unlike promoter / splice / enhancer / chromatin.
- **First-call cold-start.** The annotation pipeline is the heaviest GI model — first request after a cold service takes ~30+ s; subsequent calls are warm.
- **The model is trained on human + a few other vertebrates.** Bacterial / fungal / plant predictions are out of distribution.
- **Hackathon key is shared.** Async jobs count toward concurrent caps too — under heavy hackathon load, you may queue.

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

Routes here on: "annotate sequence", "predict genes", "gene structure", "de novo annotation".

Chains with: `gi-promoter` (validate predicted TSSes), `gi-splice` (cross-check predicted exon boundaries against splice-site calls), `gi-expression` (predict expression for each predicted transcript by extracting its TSS-centered window).

## Safety

Research and development use. Not for clinical or diagnostic decisions. Predicted gene structures are model outputs, not curated reference annotations — for clinical interpretation, anchor to RefSeq / Ensembl.
