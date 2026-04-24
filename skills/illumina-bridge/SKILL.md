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
- **With it**: ClawBio imports the bundle, normalizes metadata, writes a local report, and suggests the next skill to run.
- **Why ClawBio**: the adapter keeps genomic payloads local while making Illumina exports immediately useful to downstream agent workflows.

## Core Capabilities

1. **Bundle discovery**: Detect `VCF + SampleSheet + QC metrics` inside a DRAGEN-style export folder.
2. **Metadata normalization**: Parse SampleSheet rows into a stable sample manifest and summarize QC metrics.
3. **Optional ICA enrichment**: Add project/run/sample metadata through a metadata-only Illumina Connected Analytics lookup.
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
3. **Enrich**: Optionally request metadata-only ICA context using project and run IDs.
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
  --ica-run-id <run_id> \
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

## Algorithm / Methodology

1. **Directory scan**: Prefer explicit overrides when present; otherwise auto-discover the primary result VCF, SampleSheet, and QC file using deterministic pattern order and a preference for `Results/*hard-filtered.vcf`.
2. **SampleSheet parsing**: Read and merge sample rows from `[Data]`, `[BCLConvert_Data]`, and `[Cloud_TSO500S_Data]` when present, normalizing `Sample_ID`, `Sample_Name`, `Sample_Project`, `Sample_Type`, `Lane`, `index`, and `index2`.
3. **QC normalization**: Accept JSON, CSV, or DRAGEN `MetricsOutput.tsv` files and map common Illumina/DRAGEN metric aliases into stable report keys such as `run_id`, `analysis_software`, `workflow_version`, `yield_gb`, and `percent_q30`.
4. **Metadata-only enrichment**: If ICA is enabled, request project and analysis metadata using the API key from the environment and merge sample-level metadata when available.
5. **Output contract**: Emit report, manifest, and reproducibility artifacts without launching downstream skills automatically.

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
- **Metadata-only cloud access**: ICA enrichment is opt-in and limited to project/run metadata
- **Disclaimer**: every report includes the ClawBio medical disclaimer
- **Reproducibility**: commands, environment context, and checksums are always written

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
