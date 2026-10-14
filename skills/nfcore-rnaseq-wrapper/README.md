# nfcore-rnaseq-wrapper

ClawBio wrapper for running `nf-core/rnaseq` bulk RNA-seq preprocessing from FASTQ/BAM inputs with strict preflight, reproducibility artifacts, provenance, and explicit handoff to downstream bulk RNA-seq DE skills.

## Scope

- Upstream bulk RNA-seq preprocessing via Nextflow.
- Supported aligners: `star_salmon`, `star_rsem`, `hisat2`, `bowtie2_salmon`.
- Validated samplesheet input with normalized local FASTQ/BAM paths (local-first by default; remote URIs require `--allow-remote-inputs`).
- Reproducibility bundle, provenance bundle, `report.md`, and `result.json`.
- Canonical merged count matrix detection for downstream `rnaseq`.

## Out Of Scope

- Differential expression, contrasts, PCA, volcano plots, or interpretation.
- Automatic execution of `rnaseq-de` without the full flag set (`--run-downstream`, `--metadata`, `--formula`, `--contrast` are all required to trigger the handoff; otherwise only a replay template is written).
- Arbitrary free-form Nextflow passthrough flags (use `--nextflow-config` to supply custom config files).

## Quick Start

```bash
python clawbio.py run rnaseq-pipeline --demo --output ./outputs/rnaseq_demo
```

For real data:

```bash
python clawbio.py run rnaseq-pipeline \
  --input samplesheet.csv \
  --output ./rnaseq_run \
  --aligner star_salmon \
  --genome GRCh38
```

Use `--check` to run preflight without launching Nextflow.

## Validation and resume

Emitted parameters are checked against the pinned upstream scalar contract (types, finite numbers, bounds, enums, and filename patterns). Malformed or ambiguous samplesheets fail before normalized snapshots are written. An atomic launch manifest supports resuming failed launches; resumed `--check` also validates params and samplesheet checksums.

Local Docker runs, including demos, cap both per-process requests and the local executor budget to detected host/Docker CPU and memory capacity. Site/cloud profiles and object-store work directories keep their own resource policy. Explicit `--nextflow-config` files have final precedence. Resource caps cannot make a dataset fit insufficient RAM.

See the [compatibility audit](../../docs/nfcore-wrapper-audit.md) for upstream contracts, fixes, test results, and execution limits.
