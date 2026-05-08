# nfcore-scrnaseq-wrapper

ClawBio wrapper for running the upstream `scrnaseq` Nextflow pipeline from FASTQ inputs with strict preflight, reproducibility artifacts, provenance, and explicit handoff to downstream ClawBio scRNA skills.

## Scope

- Upstream scRNA preprocessing from FASTQ via Nextflow.
- Curated presets for `standard`/simpleaf, `star`, `kallisto`, `cellranger`, `cellrangerarc`, and `cellrangermulti`.
- Validated samplesheet input with absolute normalized FASTQ paths.
- Reproducibility bundle, provenance bundle, `report.md`, and `result.json`.
- Canonical `.h5ad` detection for downstream `scrna` and `scrna-embedding`.

## Out Of Scope

- Clustering, marker detection, normalization, scVI/scANVI, or other downstream analysis.
- Automatic chaining into downstream ClawBio skills.
- Free-form Nextflow passthrough flags.
- Running without preflight.

## Quick Start

```bash
python clawbio.py run scrnaseq-pipeline --demo --output ./outputs/scrnaseq_demo
```

For real data:

```bash
python clawbio.py run scrnaseq-pipeline \
  --input samplesheet.csv \
  --output ./scrnaseq_run \
  --preset star \
  --protocol 10XV3 \
  --genome GRCh38
```

Use `--check` to run preflight without launching Nextflow.
