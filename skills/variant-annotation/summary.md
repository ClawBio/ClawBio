# Variant Annotation Summary

**Input VCF**: `examples/demo.vcf`

**Total variants annotated**: 20

**Flagged** (ClinVar pathogenic/likely pathogenic + gnomAD AF < 0.001): **2**

## Flagged Variants

| CHROM | POS | REF | ALT | Gene | Consequence | ClinVar | gnomAD_AF |
|-------|-----|-----|-----|------|-------------|---------|----------|
| 13 | 32340300 | G | A | BRCA2 | missense_variant | conflicting_classifications_of_pathogenicity,not_provided,uncertain_significance | 6.841e-07 |
| 17 | 7676154 | G | A | TP53 | missense_variant | benign,conflicting_classifications_of_pathogenicity,likely_benign,pathogenic,uncertain_significance | 6.842e-07 |

## Reproducibility

- **Annotated TSV**: `output.tsv`
- **TSV SHA-256**: `b91adf9e9db99380dee6b6279e1268ca6f475a4b7d6cf6a61bf977ebef372eb1`
- **Genome assembly**: GRCh38
- **VEP endpoint**: `https://rest.ensembl.org/vep/homo_sapiens/region`
- **Batch size**: 200
- **Flagging threshold**: gnomAD AF < 0.001
- **Date**: 2026-03-19 17:46:01 GMT
