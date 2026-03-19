---
name: flow-bio
description: >-
  Flow.bio API bridge — authenticate, browse pipelines/samples/projects, search,
  upload data, launch pipeline executions, and check run status on any Flow instance.
version: 0.1.0
author: ClawBio Contributors
license: MIT
domain: bioinformatics
tags: [flow, flow.bio, LIMS, pipelines, samples, nextflow, bioinformatics, cloud]

inputs:
  - name: reads1
    type: file
    format: [fastq, fq, fastq.gz, fq.gz]
    description: First reads file for sample upload
    required: false
  - name: reads2
    type: file
    format: [fastq, fq, fastq.gz, fq.gz]
    description: Second reads file for paired-end upload
    required: false

outputs:
  - name: report
    type: file
    format: md
    description: Markdown report of Flow.bio interaction
  - name: result
    type: file
    format: json
    description: Machine-readable result envelope

dependencies:
  python: ">=3.10"
  packages:
    - requests>=2.28

demo_data:
  - path: data/demo_cache.json
    description: Pre-cached public Flow.bio data (pipelines, organisms, sample types) for offline demo

endpoints:
  cli: python skills/flow-bio/flow_bio.py --demo --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env:
        - FLOW_URL
        - FLOW_TOKEN
      config: []
    always: false
    emoji: "🌊"
    homepage: https://flow.bio
    os: [macos, linux]
    install:
      - kind: pip
        package: requests
        bins: []
    trigger_keywords:
      - flow
      - flow.bio
      - flow bio
      - flow pipeline
      - flow sample
      - flow execution
      - flow project
      - flow upload
      - run on flow
---

# Flow Bio Bridge

**ClawBio's gateway to the Flow.bio platform — browse, search, upload, and launch bioinformatics pipelines on Flow from the command line.**

## Why This Exists

Flow.bio is a bioinformatics platform hosting curated Nextflow pipelines with managed compute, sample tracking, and collaborative project management. But interacting with Flow requires navigating the web UI or writing custom API scripts.

Flow Bio Bridge makes the platform **agent-accessible**: authenticate once, then list pipelines, upload samples, launch executions, and poll for results — all from the CLI or via the ClawBio orchestrator.

- **Without it**: Users must switch between web UI and local analysis, manually track execution IDs, and write one-off upload scripts
- **With it**: A single CLI covers discovery, upload, execution, and status checking across any Flow instance
- **Why ClawBio**: Chain Flow pipeline outputs with local ClawBio skills (e.g. Flow RNA-seq → ClawBio diffviz)

## Core Capabilities

1. **Authentication** — Login via username/password, existing JWT token, or environment variables
2. **Pipeline discovery** — List and inspect available Nextflow pipelines with versions and parameter schemas
3. **Sample management** — List, search, and upload samples with metadata and organism tagging
4. **Project browsing** — List owned/shared projects and their contents
5. **Execution tracking** — Launch pipeline runs, poll status, and retrieve logs
6. **Search** — Full-text search across samples, projects, data, and executions
7. **Data listing** — Browse owned and shared data files
8. **Organism & type discovery** — List available organisms and sample types for uploads
9. **Overview mode** — Live overview of any Flow instance; public data without credentials, full account view with login

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| FASTQ (single-end) | `.fq`, `.fastq`, `.fq.gz` | reads1 | `sample_R1.fastq.gz` |
| FASTQ (paired-end) | `.fq`, `.fastq`, `.fq.gz` | reads1, reads2 | `sample_R1.fastq.gz`, `sample_R2.fastq.gz` |
| Any data file | any | — | annotation, reference, BED, etc. |

## Workflow

When the user asks to interact with Flow.bio:

1. **Authenticate**: Use stored token, env vars, or prompt for credentials
2. **Discover**: List pipelines, sample types, organisms as needed
3. **Act**: Upload samples, launch executions, or query status
4. **Report**: Write structured output (JSON + markdown) to output directory
5. **Bundle**: Generate reproducibility commands and checksums

## CLI Reference

```bash
# Login (stores token for subsequent calls)
python skills/flow-bio/flow_bio.py --login --username USER --password PASS
python skills/flow-bio/flow_bio.py --login --token TOKEN

# Discovery
python skills/flow-bio/flow_bio.py --pipelines
python skills/flow-bio/flow_bio.py --samples
python skills/flow-bio/flow_bio.py --projects
python skills/flow-bio/flow_bio.py --organisms
python skills/flow-bio/flow_bio.py --sample-types
python skills/flow-bio/flow_bio.py --executions
python skills/flow-bio/flow_bio.py --data

# Inspect details
python skills/flow-bio/flow_bio.py --execution EXEC_ID
python skills/flow-bio/flow_bio.py --sample SAMPLE_ID
python skills/flow-bio/flow_bio.py --pipeline PIPELINE_ID

# Inspect details via clawbio.py runner (uses --*-detail flags)
python clawbio.py run flow --execution-detail EXEC_ID
python clawbio.py run flow --sample-detail SAMPLE_ID
python clawbio.py run flow --pipeline-detail PIPELINE_ID

# Search
python skills/flow-bio/flow_bio.py --search "RNA-seq tumor"

# Raw JSON output (for piping to jq, etc.)
python skills/flow-bio/flow_bio.py --samples --json

# Upload sample
python skills/flow-bio/flow_bio.py --upload-sample \
  --name "Tumour_01" --sample-type "RNA-Seq" \
  --reads1 R1.fastq.gz --reads2 R2.fastq.gz \
  --organism "Homo sapiens" --project PROJECT_ID

# Launch pipeline
python skills/flow-bio/flow_bio.py --run-pipeline PIPELINE_VERSION_ID \
  --run-samples SAMPLE_ID1,SAMPLE_ID2 --output /tmp/flow_run

# Check execution status
python skills/flow-bio/flow_bio.py --execution EXEC_ID --output /tmp/flow_status

# Overview (public endpoints, no credentials needed)
python skills/flow-bio/flow_bio.py --demo --output /tmp/flow_demo

# Overview with authentication (shows owned samples, projects, executions)
python skills/flow-bio/flow_bio.py --demo --username USER --password PASS --output /tmp/flow_demo
```

## Demo

```bash
# Public overview (no credentials)
python skills/flow-bio/flow_bio.py --demo --output /tmp/flow_demo

# Full overview with account data
python skills/flow-bio/flow_bio.py --demo --username USER --password PASS --output /tmp/flow_demo
```

Expected output: a live overview of the Flow.bio instance showing available
pipelines, organisms, and sample types (public), plus owned samples, projects,
and executions if authenticated.

## Algorithm / Methodology

1. **Token acquisition**: POST `/login` with credentials → JWT token (5-min access, 7-day refresh)
2. **API traversal**: RESTful GET/POST with `Authorization: Bearer <token>` header
3. **Pagination**: Offset-based (`?page=N`) with lazy loading
4. **Chunked upload**: Files split into 1 MB chunks, uploaded sequentially with progress tracking
5. **Execution polling**: GET `/executions/<id>` until status reaches terminal state

**Key parameters**:
- Base URL: `https://app.flow.bio/api` (configurable via `FLOW_URL`)
- Chunk size: 1 MB (configurable)
- Rate limits: 1000 req/hr API, 100/hr uploads, 500/hr downloads

## Example Queries

- "What pipelines are available on Flow?"
- "List my samples on Flow.bio"
- "Upload these FASTQ files to Flow"
- "Run the RNA-seq pipeline on my sample"
- "What's the status of my Flow execution?"
- "Search Flow for breast cancer samples"

## Output Structure

```
output_dir/
├── report.md              # Summary of Flow API interaction
├── result.json            # Machine-readable response data
└── reproducibility/
    ├── commands.sh        # Exact CLI commands to reproduce
    └── environment.yml    # Flow instance URL, versions
```

## Dependencies

**Required:**
- `requests` >= 2.28 — HTTP client for REST API calls

**Optional:**
- `tqdm` — progress bars during file upload (graceful degradation without it)

## Safety

- **Credentials never stored on disk** — tokens held in memory or via environment variables only
- **No destructive defaults** — delete/transfer operations are not exposed in v0.1
- **Audit trail**: Every API interaction logged to reproducibility bundle
- **Disclaimer**: ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User mentions "flow", "flow.bio", "flow pipeline", "run on flow"
- User wants to upload samples to a cloud platform
- User asks about pipeline execution status on Flow

**Chaining partners**:
- `rnaseq-de` — Flow RNA-seq pipeline output → ClawBio differential expression
- `diffviz` — Flow DE results → ClawBio visualization
- `scrna-orchestrator` — Flow Cell Ranger output → ClawBio scRNA analysis
- `illumina-bridge` — Illumina DRAGEN → Flow upload → pipeline execution

## Citations

- [Flow.bio](https://flow.bio) — Bioinformatics platform
- [Nextflow](https://nextflow.io/) — Workflow engine powering Flow pipelines
