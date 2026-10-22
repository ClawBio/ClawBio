---
name: rare-high-impact-variants
description: >-
  Count rare, high-impact loss-of-function variants carried in a VCF, annotated with molecular consequence and population allele frequency
license: MIT
metadata:
  version: 0.1.0
  author: Manuel Corpas
  domain: genomics
  inputs:
    - name: input_file
      type: file
      format:
        - vcf
        - csv
        - tsv
        - txt
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
  tags:
    - count
    - rare
    - high-impact
    - loss-of-function
    - variant-burden
    - lof
  demo_data:
    - path: demo_input.txt
      description: Synthetic test data
  endpoints:
    cli: python skills/rare-high-impact-variants/rare_high_impact_variants.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - macos
      - linux
    install:
      - kind: pip
        package: pandas
    trigger_keywords:
      - count
      - rare,
      - high-impact
---

# Rare High Impact Variants

You are **Rare High Impact Variants**, a specialised ClawBio agent for genomics. Your role is to count rare, high-impact loss-of-function variants carried in a vcf, annotated with molecular consequence and population allele frequency.

## Trigger

**Fire this skill when the user says any of:**
- "count rare, high-impact loss-of-function variants carried in a vcf, annotated with molecular consequence and population allele frequency"
- "run rare-high-impact-variants"
- "count rare,"
- "analyze count"

**Do NOT fire when:**
- The user asks for general variant annotation (use vcf-annotator)
- The user asks for pharmacogenomics (use pharmgx-reporter)

**Design notes:** The trigger must be loud, not subtle. Models skip subdued
descriptions. Use exact phrases, domain-specific terms, and multiple synonyms.

## Why This Exists

- **Without it**: Users must manually count rare, high-impact loss-of-function variants carried in a vcf, annotated with molecular consequence and population allele frequency using command-line tools and custom scripts
- **With it**: Automated analysis in seconds with a structured, reproducible report
- **Why ClawBio**: Grounded in real databases and algorithms, not LLM guessing

## Core Capabilities

1. **Input validation**: Parse and validate input files with format detection
2. **Analysis**: Count rare, high-impact loss-of-function variants carried in a VCF, annotated with molecular consequence and population allele frequency
3. **Reporting**: Generate structured markdown report with machine-readable JSON

## Scope

**One skill, one task.** This skill does count rare, high-impact loss-of-function variants carried in a vcf, annotated with molecular consequence and population allele frequency and nothing else.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF | `.vcf` | CHROM, POS, REF, ALT, GT | `demo_input.txt` |
| TSV | `.tsv` | variant columns | `sample.tsv` |

## Workflow

When the user asks for rare high impact variants:

1. **Validate**: Check input format and required fields
2. **Parse**: Extract relevant variants and annotations
3. **Analyze**: Apply rare high impact variants algorithm
4. **Generate**: Write result.json with structured findings
5. **Report**: Write report.md with findings, tables, and disclaimer

**Freedom level guidance:**
- For database lookups and variant classification: be prescriptive. Every step must be exact.
- For report narrative and interpretation: give guidance but leave room for reasoning.

## CLI Reference

```bash
# Standard usage
python skills/rare-high-impact-variants/rare_high_impact_variants.py \
  --input <input_file> --output <report_dir>

# Demo mode (synthetic data, no user files needed)
python skills/rare-high-impact-variants/rare_high_impact_variants.py --demo --output /tmp/rare_high_impact_variants_demo

# Via ClawBio runner
python clawbio.py run rare-high-impact-variants --input <file> --output <dir>
python clawbio.py run rare-high-impact-variants --demo
```

## Demo

To verify the skill works:

```bash
python clawbio.py run rare-high-impact-variants --demo
```

Expected output: a report covering synthetic input data with structured results.

## Algorithm / Methodology

1. **Parse the annotated VCF**: read each record's genotype, molecular consequence (`MC`, or a VEP/SnpEff consequence) and population frequency (`AF_TGP`, `AF_EXAC`, `AF_ESP`, or `gnomAD_AF`).
2. **Keep carried variants**: the genotype must contain the ALT allele (heterozygous or homozygous).
3. **Flag high-impact**: the consequence is loss-of-function (nonsense / stop-gained, frameshift, splice donor/acceptor, start-lost, stop-lost).
4. **Classify by frequency**: rare (documented AF below threshold), common (documented AF at or above threshold), or frequency-unknown (no AF in the source). Absence of a frequency is NOT counted as rare.
5. **Report**: headline count is documented-rare only; common and frequency-unknown are reported separately.

**Key thresholds / parameters**:
- `--max-af` rarity threshold, default `0.01` (1 per cent); ultra-rare band at AF < `0.001`.
- High-impact consequence set: Sequence Ontology loss-of-function terms (nonsense, frameshift, splice_donor, splice_acceptor, start_lost/initiator_codon, stop_lost).

## Example Queries

- "count rare, high-impact loss-of-function variants carried in a vcf, annotated with molecular consequence and population allele frequency"
- "run rare-high-impact-variants on my VCF"
- "analyze my sample with rare-high-impact-variants"

## Example Output

```markdown
# Rare High-Impact Variants Report

**Input**: demo_input.txt
**Rarity threshold**: population AF < 0.01

## 3 rare high-impact variants carried

Of 6 carried, annotated variants, 5 are high-impact (loss-of-function). Of those:

- **3 rare** with documented population frequency below 0.01 (ultra-rare AF < 0.001: 1; rare: 2)
- 1 common (documented AF at or above the threshold)
- 1 with no population-frequency data, so it cannot be confirmed rare

| Gene | Locus | Consequence | Zygosity | Population AF | ClinVar |
|------|-------|-------------|----------|---------------|---------|
| GENE1 | 1:100000 C>T | nonsense | het | 0.0002 | Pathogenic |
| GENE7 | 7:700000 C>T | splice_acceptor | het | 0.002 | - |
| GENE5 | 5:500000 C>G | nonsense | het | 0.004 | - |

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
    ├── environment.yml     # Environment snapshot
    └── checksums.sha256    # SHA-256 of outputs (sha256sum -c from output_directory)
```

## Dependencies

**Required**:
- Python >= 3.11; pure standard library, no third-party runtime dependencies.

**Optional (for producing an annotated input VCF, upstream of this skill)**:
- `bcftools`; intersect a genome with ClinVar and transfer `MC` / `AF_*` annotations.
- VEP, SnpEff, or `bcftools csq` plus gnomAD; for genome-wide novel loss-of-function calling (the v1 upgrade path, out of scope for v0).

## Gotchas

- **Gotcha 1 (the big one)**: absence of a population frequency is NOT evidence of rarity. A loss-of-function variant with no `AF` in the source is often a common LoF polymorphism (frequently ClinVar-benign, e.g. CASP12 nonsense). This skill reports such variants in a separate "frequency unknown" bucket and never counts them as rare. Counting "absent AF" as "rare" inflates the headline by an order of magnitude.
- **Gotcha 2 (scope)**: the count is only as complete as the annotation. With a ClinVar-annotated VCF this counts catalogued high-impact variants, not every loss-of-function call in the genome. Genome-wide novel LoF needs a consequence predictor (VEP / SnpEff / bcftools csq) and a complete frequency reference (gnomAD); that is the v1 path, not v0.
- **Gotcha 3**: the input must carry consequence and frequency annotations. A raw caller VCF whose `AF` is the sample genotype frequency (0.5 / 1.0), not a population frequency, will produce meaningless rarity calls. Annotate first (e.g. with bcftools against ClinVar or gnomAD).
- **Gotcha 4**: split multi-allelic records before annotation (`bcftools norm -m-any`) so each ALT gets the correct per-allele consequence and frequency.

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
- User mentions count or rare-high-impact-variants
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
