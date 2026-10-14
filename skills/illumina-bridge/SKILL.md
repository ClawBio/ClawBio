---
name: illumina-bridge
description: Import DRAGEN-exported Illumina result bundles into ClawBio for local tertiary analysis and downstream routing.
license: MIT
metadata:
  version: 0.1.0
  author: ClawBio
  tags:
  - illumina
  - dragen
  - ica
  - tertiary-analysis
  - vcf
  - genomics
  openclaw:
    requires:
      bins:
      - python3
      env:
      - ILLUMINA_ICA_API_KEY
      - ILLUMINA_ICA_BASE_URL
    always: false
    emoji: 🧬
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: requests
    trigger_keywords:
    - illumina
    - dragen
    - ica
    - basespace
    - sample sheet
    - samplesheet
---

# Illumina Bridge

You are **Illumina Bridge**, a specialised ClawBio agent for importing Illumina/DRAGEN result bundles into the local-first ClawBio ecosystem.

## Why This Exists

Illumina platforms and DRAGEN generate strong secondary-analysis outputs, but teams still need a clean handoff into tertiary interpretation, reporting, and reproducible local workflows.

- **Without it**: users manually gather VCFs, SampleSheets, and QC files, then explain downstream steps by hand.
- **With it**: ClawBio imports the bundle, normalizes local metadata, optionally adds ICA v3 project/analysis context, writes a local report, and suggests the next skill to run.
- **Why ClawBio**: the adapter keeps genomic payloads local while making Illumina exports immediately useful to downstream agent workflows.

## Core Capabilities

1. **Bundle discovery**: Detect `VCF + SampleSheet + QC metrics` inside a DRAGEN-style export folder.
2. **Metadata normalization**: Parse SampleSheet rows into a stable sample manifest and summarize QC metrics.
3. **Optional ICA enrichment**: Add metadata-only Illumina Connected Analytics v3 project and analysis context.
4. **ClawBio handoff**: Write `report.md`, `result.json`, `tables/sample_manifest.csv`, and reproducibility artifacts with downstream routing hints.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| DRAGEN bundle directory | directory | `SampleSheet.csv`, one `*.vcf`/`*.vcf.gz`, one QC file | `demo_bundle/` |
| SampleSheet | `.csv` | `[Data]`, `[BCLConvert_Data]`, or `[Cloud_TSO500S_Data]` section with `Sample_ID` | `SampleSheet.csv` |
| QC metrics | `.json`, `.csv`, `.tsv` | run and quality summary metrics | `qc_metrics.json`, `MetricsOutput.tsv` |

## Workflow

1. **Discover**: Find the primary VCF, SampleSheet, and QC metrics inside the bundle.
2. **Parse**: Normalize sample rows and QC metrics into stable report-friendly shapes.
3. **Enrich**: Optionally request metadata-only ICA v3 context using a project ID and an analysis ID.
4. **Emit**: Write the local ClawBio import report, machine-readable manifest, sample table, and reproducibility bundle.

## CLI Reference

```bash
# Standard usage
python skills/illumina-bridge/illumina_bridge.py \
  --input <bundle_dir> --output <report_dir>

# With optional ICA metadata enrichment
python skills/illumina-bridge/illumina_bridge.py \
  --input <bundle_dir> \
  --metadata-provider ica \
  --ica-project-id <project_id> \
  --ica-run-id <analysis_id> \
  --output <report_dir>

# Demo mode
python skills/illumina-bridge/illumina_bridge.py --demo --output /tmp/illumina_demo

# Via ClawBio runner
python clawbio.py run illumina --input <bundle_dir> --output <dir>
python clawbio.py run illumina --demo
```

## Demo

```bash
python clawbio.py run illumina --demo
```

Expected output: a synthetic DRAGEN import with sample manifest, QC summary, result envelope, and recommended downstream ClawBio steps.

Offline ICA demo:

```bash
ILLUMINA_ICA_API_KEY= python skills/illumina-bridge/illumina_bridge.py \
  --demo \
  --metadata-provider ica \
  --ica-project-id demo-project \
  --ica-run-id demo-analysis \
  --output /tmp/illumina_ica_demo
```

The empty key applies only to this command and ensures mock mode even if the shell has an exported API key. `--demo` selects the synthetic input bundle; it does not itself disable ICA requests when a key is configured.

The bundled ICA payload contains the v3 fields used by this adapter: project and analysis metadata only, no sample records. `tables/sample_manifest.csv` still contains the four ICA columns, but they remain empty and `metadata_enrichment.merge.samples_enriched` is `0`.

## Minimal ICA Maintainer Setup

Use environment variables. Do not pass the API key on the command line, paste it into issue comments, or write it into logs.

```bash
export ILLUMINA_ICA_API_KEY="<redacted-api-key>"
export ILLUMINA_ICA_BASE_URL="https://ica.illumina.com/ica/rest"

python clawbio.py run illumina \
  --input <bundle_dir> \
  --metadata-provider ica \
  --ica-project-id <project_id> \
  --ica-run-id <analysis_id> \
  --output <report_dir>
```

`ILLUMINA_ICA_BASE_URL` defaults to `https://ica.illumina.com/ica/rest`. Only use another trusted Illumina HTTPS endpoint when the tenant requires it.

`--ica-run-id` is the ICA analysis ID for this PR. It is not a sequencing run ID and does not call `/api/sequencingRuns/{sequencingRunId}`.

## Algorithm / Methodology

1. **Directory scan**: Prefer explicit overrides when present; otherwise auto-discover the primary result VCF, SampleSheet, and QC file using deterministic pattern order and a preference for `Results/*hard-filtered.vcf`.
2. **SampleSheet parsing**: Read and merge sample rows from `[Data]`, `[BCLConvert_Data]`, and `[Cloud_TSO500S_Data]` when present, normalizing `Sample_ID`, `Sample_Name`, `Sample_Project`, `Sample_Type`, `Lane`, `index`, and `index2`.
3. **QC normalization**: Accept JSON, CSV, or DRAGEN `MetricsOutput.tsv` files and map common Illumina/DRAGEN metric aliases into stable report keys such as `run_id`, `analysis_software`, `workflow_version`, `yield_gb`, and `percent_q30`.
4. **Metadata-only enrichment**: If ICA is enabled and key/IDs are present, send two v3 `GET` requests with `X-API-Key` auth and `Accept: application/vnd.illumina.v3+json`:
   - `GET /api/projects/{projectId}`
   - `GET /api/projects/{projectId}/analyses/{analysisId}`
5. **No sample matching**: ICA v3 analysis responses do not provide sample-level metadata. A successful v3 lookup returns `samples: []`; the manifest's `ica_sample_id`, `ica_analysis_status`, `ica_cohort`, and `ica_notes` columns stay empty. Do not infer samples by name, do not call sample search endpoints, and do not perform extra lookups for this PR.
6. **Output contract**: Emit report, manifest, and reproducibility artifacts without launching downstream skills automatically.

## ICA v3 Result Contract

Successful ICA enrichment reads project and analysis metadata only. The fields below are under `data` in `result.json`; `summary.metadata_status` mirrors `data.metadata_enrichment.status`.

| Field | Meaning |
|---|---|
| `metadata_enrichment.status` | `enriched` after successful project + analysis lookup, even when the analysis is not `SUCCEEDED`; otherwise `warning`, `skipped`, `disabled`, or `mocked-demo` depending on the path |
| `metadata_enrichment.project.active` | Boolean from ICA v3 `Project.active`; if the field is absent, store `null` |
| `metadata_enrichment.project.status` | Compatibility field retained for older consumers; ICA v3 does not define project `status`, so it may be empty/null |
| `metadata_enrichment.run.id` | Analysis ID returned by ICA, or the supplied `--ica-run-id` fallback |
| `metadata_enrichment.run.name` | `AnalysisV3.userReference`; if absent, fall back to `reference` |
| `metadata_enrichment.run.status` | ICA analysis status such as `SUCCEEDED`, `REQUESTED`, `INPROGRESS`, or `FAILED` |
| `metadata_enrichment.run.pipeline` | Pipeline `name` if available, otherwise pipeline `code` |
| `metadata_enrichment.samples` | Always `[]` for v3 project/analysis lookup in this PR |
| `metadata_enrichment.merge.samples_enriched` | `0` for v3 lookup and the offline mock |

Reports must show the analysis status. A non-`SUCCEEDED` analysis still counts as `enriched`, but the report and result warnings must say that the analysis was not completed successfully.

## ICA Error Handling

The local import must finish whenever possible. ICA failures should be warnings, not hard failures.

- Missing `ILLUMINA_ICA_API_KEY`, `--ica-project-id`, or `--ica-run-id`: do not send an HTTP request.
- 401: check the API key.
- 403: check permissions for the tenant/project/analysis.
- 404: check the project ID, analysis ID, and project access.
- 429/5xx: warn and suggest retrying later; keep the local import output. This adapter does not retry automatically.
- Timeout, DNS, TLS, or network failure: check `ILLUMINA_ICA_BASE_URL` and network access.
- Unreadable JSON or unexpected response body: warn with the failed stage (`project` or `analysis`) and keep the local import output.

## Example Queries

- "Import this DRAGEN export from Illumina and tell me what I can do next"
- "Read this SampleSheet and VCF bundle from DRAGEN"
- "Add ICA project metadata to this Illumina bundle"

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── tables/
│   └── sample_manifest.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**:
- `requests` — optional ICA metadata lookup

**Optional**:
- `ILLUMINA_ICA_API_KEY` — enables metadata-only ICA enrichment
- `ILLUMINA_ICA_BASE_URL` — override the ICA API root with a trusted `https://*.illumina.com` endpoint if needed

## Safety

- **Local-first**: genomic files are read locally; the skill never uploads VCF payloads
- **Metadata-only cloud access**: ICA enrichment is opt-in and limited to project/analysis metadata
- **No command-line secrets**: `ILLUMINA_ICA_API_KEY` belongs in the environment, not in CLI flags or logs
- **Disclaimer**: every report includes the ClawBio medical disclaimer
- **Reproducibility**: commands, environment context, and checksums are always written

## Known Limits

- Issue #74's 2026-09-28 comment reports a contributor's live check against a scratch ICA environment and a small non-DRAGEN analysis. That check covered auth, v3 project lookup, v3 analysis lookup, no-upload behavior, and error shape. Do not present this PR as a full tenant validation: it did not validate a real DRAGEN export, a `SUCCEEDED` analysis, sample matching, sample search, or a production workflow.
- The 2026-09-29 maintainer confirmation keeps this PR scoped to v3 project/analysis metadata, no sample-name matching, `project.active` plus compatible `project.status`, and `--ica-run-id` as an analysis ID.
- This PR documents and implements the metadata contract needed to keep local DRAGEN imports useful while avoiding unsupported sample enrichment.

## Testing

Offline checks:

```bash
python skills/illumina-bridge/illumina_bridge.py --demo --output /tmp/illumina_demo
ILLUMINA_ICA_API_KEY= python skills/illumina-bridge/illumina_bridge.py \
  --demo --metadata-provider ica \
  --ica-project-id demo-project \
  --ica-run-id demo-analysis \
  --output /tmp/illumina_ica_demo
python -m pytest skills/illumina-bridge/tests/ -v
```

Maintainer-only live check:

```bash
export ILLUMINA_ICA_API_KEY="<redacted-api-key>"
export ILLUMINA_ICA_BASE_URL="https://ica.illumina.com/ica/rest"

python clawbio.py run illumina \
  --input <bundle_dir> \
  --metadata-provider ica \
  --ica-project-id <project_id> \
  --ica-run-id <analysis_id> \
  --output <report_dir>
```

Use placeholder IDs in docs and logs. Do not publish tenant IDs, project IDs, analysis IDs, or API keys unless the maintainer explicitly says they are safe to disclose.

## Integration with Bio Orchestrator

**Trigger conditions**:
- queries mentioning Illumina, DRAGEN, ICA, BaseSpace, SampleSheet, or sample sheet
- directories that contain a recognizable Illumina bundle (`SampleSheet + VCF`)

**Chaining partners**:
- `equity-scorer`: cohort-level follow-up on imported VCFs
- `clinpgx`: targeted gene-drug follow-up after DRAGEN review
- `gwas-lookup`: per-variant external lookup from imported findings

## Citations

- [DRAGEN secondary analysis](https://www.illumina.com/products/by-type/informatics-products/dragen-secondary-analysis.html)
- [Illumina Connected Analytics](https://www.illumina.com/products/by-type/informatics-products/connected-analytics.html)
- [BCL Convert Sample Sheet](https://support-docs.illumina.com/SW/BCL_Convert/Content/SW/BCLConvert/SampleSheets_swBCL.htm)
- [ICA public OpenAPI v3 spec](https://ica.illumina.com/ica/api/swagger/openapi_public.yaml) (`info.version: "3"`, server path `/ica/rest`)
- [Issue #74 live ICA validation comment, 2026-09-28](https://github.com/ClawBio/ClawBio/issues/74#issuecomment-5876103282)
- [Issue #74 maintainer confirmation, 2026-09-29](https://github.com/ClawBio/ClawBio/issues/74#issuecomment-5885465867)
