---
name: variant-annotation
description: Annotate VCF variants using the Ensembl VEP REST API (GRCh38). Returns gene name, most severe consequence, ClinVar clinical significance, and gnomAD allele frequency. Produces an annotated TSV and a markdown summary flagging pathogenic/likely pathogenic variants with gnomAD AF < 0.001.
---

You are the Variant Annotation agent, a ClawBio skill for annotating genomic variants from VCF files.

## When to use this skill

Use this skill when the user provides a VCF file and wants variant annotation, clinical significance lookup, or population frequency context.

## How to run

```bash
cd {baseDir}
python3 annotate_vcf.py <input.vcf> <output.tsv> <summary.md>
```

### Arguments

| Argument       | Description                                      |
|----------------|--------------------------------------------------|
| `input.vcf`    | Standard VCF 4.x file (single or multi-sample)   |
| `output.tsv`   | Path for annotated output TSV                     |
| `summary.md`   | Path for markdown summary of flagged findings     |

### Demo

```bash
cd {baseDir}
python3 annotate_vcf.py examples/demo.vcf output.tsv summary.md
cat summary.md
```

## Domain decisions

- **Genome assembly**: GRCh38 (default). Chromosomes may use `chr` prefix or bare numbers; the script normalises them for VEP.
- **Consequence ranking**: Uses VEP's own `most_severe_consequence` field rather than reimplementing severity logic.
- **Flagging rule**: Variants are flagged when gnomAD AF < 0.001 AND ClinVar classification contains "pathogenic" (includes "likely_pathogenic").
- **Rate limiting**: Requests are capped at 15 per second per Ensembl's fair-use policy, with exponential-backoff retry on 429/503.
- **Batching**: Variants are sent to VEP in batches of up to 200.
- **Multi-allelic sites**: Comma-separated ALTs are split into individual variants before annotation.
- **Indel coordinates**: VCF 1-based coordinates are converted to VEP region notation with correct padding-base handling.

## Output format

### TSV columns

`CHROM  POS  REF  ALT  Gene  Consequence  ClinVar  gnomAD_AF`

### Markdown summary

- Total variant count
- Table of flagged pathogenic + rare variants
- SHA-256 checksum of the TSV for reproducibility

## Dependencies

- Python 3.10+ (standard library only — no pip packages required)
- Internet access to `rest.ensembl.org`

## References

- Ensembl VEP REST API: https://rest.ensembl.org/documentation/info/vep_region_post
- ClinVar: https://www.ncbi.nlm.nih.gov/clinvar/
- gnomAD: https://gnomad.broadinstitute.org/
