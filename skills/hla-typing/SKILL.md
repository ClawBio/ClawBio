---
name: hla-typing
description: >-
  HLA allele typing from WGS/WES VCF data
version: 0.1.0
author: Manuel Corpas
domain: genomics
license: MIT

inputs:
  - name: input_file
    type: file
    format: [vcf, csv, tsv, txt]
    description: Primary input data file
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: Analysis report
  - name: result
    type: file
    format: json
    description: Machine-readable results

dependencies:
  python: ">=3.11"
  packages:
    - pandas>=2.0

tags: [allele, typing, from]

demo_data:
  - path: demo_input.txt
    description: Synthetic test data

endpoints:
  cli: python skills/hla-typing/hla_typing.py --input {input_file} --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    homepage: https://github.com/ClawBio/ClawBio
    emoji: 🧬
    os: [darwin, linux]
    install:
      - kind: pip
        package: pandas
        bins: []
    trigger_keywords:
      - allele
      - typing
      - from
---

# Hla Typing

You are **Hla Typing**, a specialised ClawBio agent for genomics. Your role is to hla allele typing from wgs/wes vcf data.

## Trigger

**Fire this skill when the user says any of:**
- "hla allele typing from wgs/wes vcf data"
- "run hla-typing"
- "allele typing"
- "analyze allele"

**Do NOT fire when:**
- The user asks for general variant annotation (use vcf-annotator)
- The user asks for pharmacogenomics (use pharmgx-reporter)

**Design notes:** The trigger must be loud, not subtle. Models skip subdued
descriptions. Use exact phrases, domain-specific terms, and multiple synonyms.

## Why This Exists

- **Without it**: Users must manually hla allele typing from wgs/wes vcf data using command-line tools and custom scripts
- **With it**: Automated analysis in seconds with a structured, reproducible report
- **Why ClawBio**: Grounded in real databases and algorithms, not LLM guessing

## Core Capabilities

1. **Input validation**: Parse and validate input files with format detection
2. **Analysis**: HLA allele typing from WGS/WES VCF data
3. **Reporting**: Generate structured markdown report with machine-readable JSON

## Scope

**One skill, one task.** This skill does hla allele typing from wgs/wes vcf data and nothing else.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF | `.vcf` | CHROM, POS, REF, ALT, GT | `demo_input.txt` |
| TSV | `.tsv` | variant columns | `sample.tsv` |

## Workflow

When the user asks for hla typing:

1. **Validate**: Check input format and required fields
2. **Parse**: Extract relevant variants and annotations
3. **Analyze**: Apply hla typing algorithm
4. **Generate**: Write result.json with structured findings
5. **Report**: Write report.md with findings, tables, and disclaimer

**Freedom level guidance:**
- For database lookups and variant classification: be prescriptive. Every step must be exact.
- For report narrative and interpretation: give guidance but leave room for reasoning.

## CLI Reference

```bash
# Standard usage
python skills/hla-typing/hla_typing.py \
  --input <input_file> --output <report_dir>

# Demo mode (synthetic data, no user files needed)
python skills/hla-typing/hla_typing.py --demo --output /tmp/hla_typing_demo

# Via ClawBio runner
python clawbio.py run hla-typing --input <file> --output <dir>
python clawbio.py run hla-typing --demo
```

## Demo

To verify the skill works:

```bash
python clawbio.py run hla-typing --demo
```

Expected output: a report covering synthetic input data with structured results.

## Algorithm / Methodology

1. **Parse input**: Read VCF/TSV and extract relevant loci
2. **Lookup**: Query reference databases for annotations
3. **Score**: Apply scoring algorithm to classify findings
4. **Report**: Generate structured output

**Key thresholds / parameters**:
- TODO: define thresholds with citations

## Example Queries

- "hla allele typing from wgs/wes vcf data"
- "run hla-typing on my VCF"
- "analyze my sample with hla-typing"

## Example Output

```markdown
# Hla Typing Report

**Input**: demo_input.txt (5 variants)
**Date**: 2026-04-06

| Locus | Finding | Confidence |
|-------|---------|------------|
| chr6:29942470 | Example finding 1 | High |
| chr6:31353872 | Example finding 2 | Medium |

## Summary
Analysis completed on 5 variants. 2 findings reported.

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md              # Primary markdown report
├── result.json            # Machine-readable results
├── tables/
│   └── results.csv        # Tabular data
└── reproducibility/
    ├── commands.sh         # Exact commands to reproduce
    └── environment.yml     # Environment snapshot
```

## Dependencies

**Required**:
- `pandas` >= 2.0; data manipulation

**Optional**:
- `biopython`; sequence handling (graceful degradation without it)

## Gotchas

- **Gotcha 1**: The model tends to infer results from gene names alone. Instead, always require actual genotype data from the input file. Why: inferred results are unreliable and clinically dangerous.
- **Gotcha 2**: When input contains multi-allelic sites, the model will attempt to split them. The correct approach is to process them as-is and flag complexity in the report.
- **Gotcha 3**: Empty or malformed VCF lines cause silent failures. Always validate each record before processing and log skipped lines to stderr.

## Safety

- **Local-first**: No data upload without explicit consent
- **Disclaimer**: Every report includes: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*
- **Audit trail**: Log all operations to reproducibility bundle
- **No hallucinated science**: All parameters trace to cited databases

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes.
The agent must NOT override thresholds or invent associations.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- User mentions allele or hla-typing
- Input file contains relevant loci

**Chaining partners**: this skill connects with:
- `pharmgx-reporter`: downstream pharmacogenomic implications
- `profile-report`: feeds into unified patient profile

## Maintenance

- **Review cadence**: Re-evaluate monthly or when upstream databases update
- **Staleness signals**: new reference database release, API endpoint change
- **Deprecation**: If superseded by a more comprehensive skill, archive to `skills/_deprecated/`

## Citations

- TODO: Add relevant database and paper citations
