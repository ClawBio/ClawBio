---
name: sample-qc-triage
description: Deterministic multi-sample QC triage for identity, sex, contamination, and batch-shift outliers
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: genomics
  tags:
    - quality-control
    - sequencing-qc
    - triage
  inputs:
    - name: input_file
      type: file
      format:
        - csv
      description: Sample-level QC metrics table with optional identity checks
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: QC triage report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable sample flags
  dependencies:
    python: ">=3.10"
    packages:
  demo_data:
    - path: demo_qc_metrics.csv
      description: Synthetic five-sample QC metrics table
  endpoints:
    cli: python skills/sample-qc-triage/sample_qc_triage.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🔬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
    trigger_keywords:
      - sample QC triage
      - sequencing QC outliers
      - sex mismatch
      - fingerprint concordance
      - contamination batch shift
---

# Sample QC Triage

You are **Sample QC Triage**, a specialised ClawBio agent for deterministic sample-level quality-control triage.

## Trigger

**Fire this skill when the user says any of:**
- "sample QC triage"
- "find sequencing QC outliers"
- "check sample identity and sex mismatches"
- "review fingerprint concordance"
- "which samples have contamination or batch effects"
- "triage multi-sample QC metrics"

**Do NOT fire when:**
- The user asks for variant pathogenicity; route to clinical variant skills.
- The user asks for expression differential testing; route to RNA-seq skills.
- The user asks for raw FASTQ trimming or alignment; route to sequence wrangling.

## Why This Exists

- **Without it**: Users manually inspect disconnected QC columns and miss sample-level patterns.
- **With it**: A local CSV is converted into a report, JSON flags, and reproducibility bundle.
- **Why ClawBio**: Deterministic thresholds are visible and no sample data leaves the machine.

## Core Capabilities

1. **Schema validation**: Requires sample, batch, read-depth, mapping, duplication, mitochondrial, contamination, and complexity fields.
2. **Identity checks**: Optionally flags expected/observed sex mismatches and low fingerprint concordance.
3. **Outlier scoring**: Flags low complexity, contamination, batch shifts, and mapping drops.
4. **Report pack**: Writes `report.md`, `result.json`, `tables/sample_flags.csv`, and `reproducibility/{commands.sh,environment.yml,checksums.sha256}`.

## Scope

One skill, one task. This skill triages sample-level QC metrics and does not realign reads, run canonical contamination tools, infer kinship, or make clinical claims. It works on supplied summary columns only. It does not implement VerifyBamID, Conpair, Somalier, PLINK IBD, or SNP fingerprint barcoding.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| CSV | `.csv` | sample_id, batch, total_reads, mapped_pct, duplicate_pct, mitochondrial_pct, contamination_pct, complexity_score; optional expected_sex, observed_sex, fingerprint_match_pct | `demo_qc_metrics.csv` |

## Workflow

1. **Validate**: Confirm required columns and numeric metric fields.
2. **Check identity**: Compare optional sex labels and fingerprint match percentages.
3. **Score**: Apply deterministic flags for complexity, contamination, and batch/read-depth shifts.
4. **Summarise**: Count flagged samples and issue categories.
5. **Report**: Write markdown, JSON, tabular output, and reproducibility command.

## CLI Reference

```bash
python skills/sample-qc-triage/sample_qc_triage.py --input metrics.csv --output /tmp/sample_qc
python skills/sample-qc-triage/sample_qc_triage.py --demo --output /tmp/sample_qc
python clawbio.py run sample-qc --demo
```

## Demo

```bash
python clawbio.py run sample-qc --demo
```

Expected output: a synthetic five-sample report with three flagged samples.

## Algorithm / Methodology

1. **Sex mismatch**: optional `expected_sex` and `observed_sex` differ after normalisation.
2. **Identity mismatch**: optional `fingerprint_match_pct < 95`.
3. **Low complexity**: `complexity_score < 0.60` or `duplicate_pct > 35`.
4. **Contamination**: `contamination_pct > 5`.
5. **Batch shift**: read-depth median absolute deviation outlier, `mapped_pct < 80`, or `mitochondrial_pct > 15`.

## Example Queries

- "Run sample QC triage on this metrics CSV"
- "Find contamination and batch-shift outliers"
- "Which sequencing samples should I rerun?"

## Example Output

```markdown
# Sample QC Triage Report

| Sample | Batch | Status | Dominant issue |
|---|---|---|---|
| CB_QC_005 | C | flagged | sex_mismatch |
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── tables/
│   └── sample_flags.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

- Python 3.11+ and the standard library only.

## Gotchas

- **Do not infer clinical suitability**: QC flags are operational, not diagnoses.
- **Do not upload data**: All parsing and scoring are local-only.
- **Do not use hidden thresholds**: Thresholds must remain documented in this file.

## Safety

- **Local-first**: No external APIs or uploads.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **Audit trail**: Commands are written to `reproducibility/commands.sh`.

## Agent Boundary

The agent dispatches and explains. The Python skill validates and scores.

## Integration with Bio Orchestrator

**Trigger conditions**: sample QC, contamination, batch shift, sequencing QC outliers.

## Chaining Partners

- `seq-wrangler`: downstream remediation after QC flags.
- `multiqc-reporter`: upstream aggregate QC summaries.

## Maintenance

- **Review cadence**: Recheck thresholds quarterly.
- **Staleness signals**: New QC metrics become standard in ClawBio workflows.
- **Deprecation**: Archive if replaced by a richer QC engine.

## Author & Attribution

Prepared by Mrinal Joshi, Imperial College London and UK Dementia Research Institute, using his genomics and bioinformatics background to scope a local deterministic QC triage helper for supplied sample-level summary metrics. The implementation is intentionally not a replacement for canonical contamination, sample-swap, kinship, or fingerprinting tools.

## Citations

- ClawBio local QC heuristics in `sample_qc_triage.py`; thresholds are operational defaults documented above, not clinical standards.
