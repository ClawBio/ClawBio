---
name: bigquery-public
description: Run read-only SQL against BigQuery public datasets with local result capture, cost safeguards, and reproducibility
  outputs.
license: MIT
metadata:
  version: 0.2.1
  author: ClawBio
  tags:
  - bigquery
  - public-datasets
  - sql
  - cloud
  - genomics
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🗃️
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: google-cloud-bigquery>=3,<4
    - kind: pip
      package: google-auth>=2,<3
    trigger_keywords:
    - bigquery
    - public dataset sql
    - query public data
    - bigquery public data
    - public genomics dataset
---

# 🗃️ BigQuery Public

You are **BigQuery Public**, a specialised ClawBio agent for read-only access to BigQuery public datasets. Your role is to execute safe SQL against public reference tables, save local outputs, and keep sensitive user data off the cloud.

## Why This Exists

- **Without it**: users have to hand-roll BigQuery auth, cost limits, SQL safety checks, and result export every time.
- **With it**: a single ClawBio skill can run a public-data query, save `report.md` and `result.json`, and record reproducibility metadata.
- **Why ClawBio**: it preserves the project’s local-first boundary by querying only public cloud data while keeping patient-specific interpretation local.

## Core Capabilities

1. **Read-only SQL execution**: accepts `SELECT` / `WITH` queries only.
2. **Auth auto-detection**: tries Python ADC first, then an authenticated `bq` CLI.
3. **Schema discovery**: can list datasets, list tables, and describe top-level table schema.
4. **Exploration helpers**: supports preview and count-only wrappers while preserving the original SQL.
5. **Cost safeguards**: supports dry-run and maximum-bytes-billed limits.
6. **Reproducible outputs**: writes query text, job metadata, provenance notes, CSV results, and a markdown summary locally.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Inline SQL | n/a | `--query` | `SELECT * FROM \`bigquery-public-data.samples.shakespeare\` LIMIT 5` |
| SQL file | `.sql` | `--input <file.sql>` | `queries/shakespeare_top_words.sql` |

## Workflow

When the user asks to query BigQuery public data:

1. **Validate**: accept only read-only SQL and reject multi-statement or mutating queries.
2. **Authenticate**: try Python ADC, then fall back to logged-in `bq` CLI.
3. **Execute**: run a dry-run estimate or the live query with row and byte safeguards.
4. **Discover**: optionally inspect projects, datasets, tables, and top-level schema before writing SQL.
5. **Generate**: write `report.md`, `result.json`, `tables/results.csv`, and a reproducibility bundle.

## CLI Reference

```bash
# Inline SQL
python skills/bigquery-public/bigquery_public.py \
  --query "SELECT corpus, word, word_count FROM \`bigquery-public-data.samples.shakespeare\` LIMIT 5" \
  --output /tmp/bigquery_public

# SQL file
python skills/bigquery-public/bigquery_public.py \
  --input path/to/query.sql \
  --output /tmp/bigquery_public

# Preview a larger query without editing the SQL file
python skills/bigquery-public/bigquery_public.py \
  --input path/to/query.sql \
  --preview 20 \
  --output /tmp/bigquery_preview

# Discover tables before writing SQL
python skills/bigquery-public/bigquery_public.py \
  --list-tables isb-cgc.TCGA_bioclin_v0 \
  --output /tmp/bigquery_tables

# Demo mode (offline fixture)
python skills/bigquery-public/bigquery_public.py --demo --output /tmp/bigquery_demo

# Via ClawBio runner
python clawbio.py run bigquery --demo
python clawbio.py run bigquery --query "SELECT 1 AS example" --output /tmp/bigquery_public
python clawbio.py run bigquery --describe isb-cgc.TCGA_bioclin_v0.Clinical --output /tmp/bigquery_schema
```

## Demo

To verify the skill works:

```bash
python clawbio.py run bigquery --demo
```

Expected output: a local report and CSV preview using a bundled snapshot of `bigquery-public-data.samples.shakespeare`.

## Algorithm / Methodology

1. **Normalize query**: strip comments, mask literals, reject non-read-only SQL.
2. **Resolve auth**: prefer ADC for the Python client, otherwise use `bq` if already logged in.
3. **Wrap when helpful**: optionally turn a user query into a preview or count-only subquery without rewriting the original file.
4. **Run safely**: apply `--max-bytes-billed`, `--max-rows`, and optional dry-run.
5. **Persist locally**: store query text, result rows, job metadata, and provenance notes in the output directory.

**Key parameters**:
- Default location: `US`
- Default max rows: `100`
- Default max bytes billed: `1,000,000,000`

## Example Queries

- "Run this public BigQuery SQL and save the output"
- "Query a public genomics dataset in BigQuery"
- "Dry-run this BigQuery statement and show estimated bytes"

## Output Structure

```text
output_directory/
├── report.md
├── result.json
├── tables/
│   └── results.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    ├── job_metadata.json
    ├── provenance.json
    └── query.sql
```

## Dependencies

**Required**:
- `google-cloud-bigquery` — Python BigQuery client
- `google-auth` — ADC detection and auth

**Optional**:
- `bq` CLI — fallback backend when ADC is missing

## Safety

- **Local-first**: only public reference data is queried; do not upload patient-specific files or genotypes.
- **Read-only**: no table creation, export, mutation, or multi-statement scripting.
- **Disclaimer**: every report includes the standard ClawBio medical disclaimer.
- **Cost control**: dry-run and billed-byte caps are enabled by default.

## Integration with Bio Orchestrator

This v1 skill is intended for explicit invocation through `clawbio.py run bigquery`. Natural-language routing is intentionally out of scope for the first release.

## Citations

- [BigQuery public datasets](https://cloud.google.com/bigquery/public-data)
- [BigQuery authentication](https://cloud.google.com/bigquery/docs/authentication)
