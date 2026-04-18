# 🧬 VCF Annotator

**ClawBio skill** — Annotate VCF variants with Ensembl VEP, ClinVar, and gnomAD.

Part of the [ClawBio](https://github.com/ClawBio/ClawBio) bioinformatics AI agent skill library.

## Quick Start

```bash
# Demo (no network, no VCF file needed)
python skills/vcf-annotator/vcf_annotator.py --demo --output /tmp/demo

# Real VCF
python skills/vcf-annotator/vcf_annotator.py \
    --input variants.vcf \
    --output report/
```

## What it produces

```
output/
├── report.md                    # Full annotation report
├── results.json                 # All variants as JSON
├── tables/variants.csv          # Tabular variant data
└── reproducibility/
    ├── commands.sh              # Reproduce this exact run
    ├── environment.yml          # Python environment
    └── checksums.sha256         # SHA-256 of all outputs
```

## Run Tests

```bash
python -m pytest skills/vcf-annotator/tests/ -v
```

## Author

Contributed by [Sooraj](https://github.com/sooraj-codes) — resolves [Issue #6](https://github.com/ClawBio/ClawBio/issues/6).
