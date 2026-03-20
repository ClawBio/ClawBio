---
name: variant-annotation
description: Annotate VCF variants with Ensembl VEP REST, ClinVar significance, gnomAD/population frequency context, and prioritized variant ranking.
version: 0.1.0
author: Toby Clark
license: MIT
domain: genomics
tags: [genomics, vcf, variant-annotation, vep, clinvar, gnomad]
inputs:
  - name: input
    type: file
    format: [vcf, vcf.gz]
    description: Input VCF containing variant records and optional sample genotype columns
outputs:
  - name: report
    type: file
    format: markdown
    description: Variant annotation summary report with prioritized findings
  - name: result
    type: file
    format: json
    description: Machine-readable annotation results and summary metrics
  - name: annotated_variants
    type: file
    format: tsv
    description: Flat per-variant annotation table with consequence, ClinVar, and frequency fields
  - name: reproducibility
    type: directory
    description: Reproduction commands and run metadata for the analysis
demo_data:
  - path: example_data/synthetic_clinvar_panel.vcf
    description: Bundled synthetic 20-variant VCF used for demo mode
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: uv
        package: pysam
        bins: []
      - kind: uv
        package: requests
        bins: []
    trigger_keywords:
      - vcf
      - variant annotation
      - vep
      - clinvar
      - gnomad
      - annotate variants
      - pathogenic variants
---

# 🧬 Variant Annotation

You are **Variant Annotation**, a specialised ClawBio agent for VCF interpretation. Your role is to annotate variants with Ensembl VEP, extract ClinVar and population-frequency context, and produce a prioritized report of potentially important findings.

## Why This Exists

- **Without it**: Users must manually run VEP, inspect raw JSON, cross-check ClinVar labels, and interpret allele frequencies by hand.
- **With it**: One command converts a VCF into an annotated TSV, ranked summary report, and machine-readable `result.json`.
- **Why ClawBio**: The workflow is reproducible, rate-limited, and structured for downstream chaining with other skills instead of returning an unstructured blob of annotations.

## Core Capabilities

1. **VCF Parsing**: Reads standard VCF 4.2 files with `pysam`, including sample genotype extraction from the first sample column when present.
2. **Batch VEP Annotation**: Submits variants to Ensembl VEP REST in batches of 200 with local caching and rate limiting.
3. **Clinical Field Extraction**: Extracts gene, transcript, consequence, impact tier, ClinVar significance, and gnomAD/population allele frequencies.
4. **Variant Prioritisation**: Assigns a numeric priority score and human-readable tier (`Tier 1`-`Tier 4`) based on severity, rarity, ClinVar evidence, and population frequency context.
5. **Report Generation**: Writes `report.md`, `tables/annotated_variants.tsv`, `result.json`, and a reproducibility bundle.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF 4.2 | `.vcf`, `.vcf.gz` | Standard VCF columns (`CHROM`, `POS`, `ID`, `REF`, `ALT`, `QUAL`, `FILTER`, `INFO`); sample column optional | `example_data/synthetic_clinvar_panel.vcf` |

## Workflow

1. **Parse**: Read the VCF with `pysam.VariantFile` and emit one record per ALT allele.
2. **Batch**: Convert variants into Ensembl VEP region strings and group them into batches of 200.
3. **Annotate**: POST batches to `https://rest.ensembl.org/vep/homo_sapiens/region` using GRCh38 as the default assembly.
4. **Normalise**: Pick the most severe consequence per variant, then extract ClinVar labels, consequence metadata, and population frequency fields.
5. **Prioritise**: Flag rare pathogenic variants (`gnomAD AF < 0.001`) and assign a numeric score plus tier for ranked output.
6. **Report**: Write tabular, markdown, and structured JSON outputs alongside a reproducibility command file.

## CLI Reference

```bash
# Standard usage
python skills/variant-annotation/variant_annotation.py \
  --input <input.vcf> --output <report_dir>

# Demo mode
python skills/variant-annotation/variant_annotation.py \
  --demo --output /tmp/variant_annotation_demo

# Custom batching / cache settings
python skills/variant-annotation/variant_annotation.py \
  --input <input.vcf> --output <report_dir> \
  --batch-size 200 --cache-dir ~/.clawbio/variant_annotation_cache

# Via ClawBio runner (after registry entry is added)
python clawbio.py run variant-annotation --input <file> --output <dir>
python clawbio.py run variant-annotation --demo
```

## Demo

```bash
python skills/variant-annotation/variant_annotation.py --demo --output /tmp/variant_annotation_demo
```

Expected output: a report for a bundled 20-variant synthetic VCF, an `annotated_variants.tsv` table with ClinVar/frequency/prioritization fields, and a `result.json` summary of clinically relevant and top-priority variants.

## Algorithm / Methodology

1. **VCF parsing**: Use `pysam.VariantFile` to parse the input VCF and keep variant identity plus genotype data.
2. **Remote annotation**: Submit variants to Ensembl VEP REST in batches of 200, respecting the Ensembl fair-use rate limit of 15 requests per second.
3. **Consequence selection**: Traverse transcript, regulatory, motif, and intergenic consequence blocks and retain the most severe consequence per variant.
4. **Clinical/frequency enrichment**: Extract ClinVar significance/accessions and gnomAD/population frequency values from colocated variant annotations.
5. **Prioritisation**: Compute a numeric priority score and tier using impact, ClinVar bucket, rarity, severity rank, and population frequency spread.
6. **Output generation**: Produce a flat TSV, markdown summary, `result.json`, and reproducibility metadata.

**Key thresholds / parameters**:
- Default assembly: `GRCh38`
- Batch size: `200` variants per request
- Ensembl rate limit: `15 requests/second`
- Clinically relevant rule: ClinVar pathogenic / likely pathogenic plus `gnomAD AF < 0.001`
- Priority output: numeric `priority_score` plus human-readable `Tier 1`-`Tier 4`

## Domain Decisions

- **Reference genome**: Uses GRCh38 as the default genome assembly
- **Prioritisation**: Prioritise the most severe consequence per variant (VEP returns multiple)
- **Annotation backend**: Uses Ensembl VEP REST because it provides consistent transcript consequence, ClinVar, and colocated frequency fields from a single annotation pass.
- **Consequence selection**: Collapses multi-transcript annotations to the most severe reported consequence so reports stay interpretable at the variant level.
- **ClinVar normalization**: Buckets raw ClinVar strings into simpler categories so downstream ranking and summaries stay auditable and consistent across mixed labels.
- **Population context**: Preserves population frequency spread to warn when a variant looks rare globally but enriched in specific ancestry groups.

## Example Queries

- "Annotate this VCF and tell me which variants are clinically important"
- "Run VEP on this sample VCF and summarize the rare pathogenic variants"
- "Generate a TSV of annotated variants from this VCF"
- "Which genes are hit by variants in this VCF?"
- "Annotate the bundled demo VCF"

## Output Structure

```
output_directory/
├── report.md                      # Markdown summary of prioritized findings
├── result.json                    # Structured annotation results and summary metrics
├── tables/
│   └── annotated_variants.tsv     # Flat variant-level annotation table
└── reproducibility/
    └── commands.sh                # Exact command used to generate the report
```

## Dependencies

**Required**:
- Python 3.10+
- `pysam` — VCF parsing
- `requests` — Ensembl REST API access

**Optional / Planned**:
- Local Ensembl `vep` backend — planned future replacement for the REST backend when fully local annotation is needed

## Safety

- **Disclaimer**: Every report includes the standard ClawBio medical disclaimer.
- **Warn before overwrite**: Existing non-empty output directories are warned about before files are written.
- **Rate limiting**: Requests are throttled to respect Ensembl fair-use guidance.
- **Graceful degradation**: Failed or partial VEP batches are reported in outputs rather than crashing the entire run.
- **Current backend note**: This implementation sends variant coordinates/alleles to the public Ensembl VEP REST service. A local VEP backend is planned for stricter local-first workflows.

## Safety Rules

- **Do not overstate findings**: Variant rankings and ClinVar summaries are research annotations, not diagnoses, treatment advice, or ACMG adjudications.
- **Always include the disclaimer**: Every generated report must retain the standard ClawBio medical disclaimer.
- **Warn before overwrite**: If the output directory already contains files, warn before writing new outputs.
- **Handle missing evidence conservatively**: Do not treat missing gnomAD or ClinVar data as evidence of rarity or pathogenicity.
- **Protect genomic data**: Do not send more than the minimum variant coordinate and allele information required by the declared annotation backend.

## Agent Boundary

- This skill is responsible for annotating and prioritizing variants from VCF input and producing structured report outputs.
- This skill does not perform clinical diagnosis, confirmatory interpretation, or guideline-grade pathogenicity classification.
- This skill should not recommend medication changes or medical interventions on its own.
- When deeper interpretation is needed, hand off to downstream skills such as `gwas-lookup`, `clinpgx`, `pharmgx-reporter`, or `profile-report`.

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- The user provides a `.vcf` / `.vcf.gz` file and asks for annotation or interpretation.
- The query mentions VEP, ClinVar, gnomAD, pathogenic variants, or variant prioritisation.
- The user wants a ranked list of interesting variants from a VCF.

**Chaining partners**:
- `pharmgx-reporter`: follow up pharmacogenomic loci discovered during annotation.
- `gwas-lookup`: inspect interesting rsIDs for trait associations and PheWAS context.
- `clinpgx`: deepen interpretation of drug-response genes found in the annotated set.
- `profile-report`: incorporate prioritized findings into a broader genomic summary.

## Citations

- [Ensembl Variant Effect Predictor](https://www.ensembl.org/info/docs/tools/vep/index.html) — functional consequence annotation
- [Ensembl REST API](https://rest.ensembl.org/) — batch VEP annotation endpoint used by the current backend
- [ClinVar](https://www.ncbi.nlm.nih.gov/clinvar/) — clinical significance assertions
- [gnomAD](https://gnomad.broadinstitute.org/) — population allele frequency reference data
- [VCF Specification](https://samtools.github.io/hts-specs/VCFv4.2.pdf) — variant file format reference
