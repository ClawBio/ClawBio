# nfcore-rnastructurome-wrapper

SKILL.md-only ClawBio wrapper for running `nf-core/rnastructurome` — chemical-probing RNA structure analysis (SHAPE/DMS, RT-stop/MaP) from FASTQ to per-base reactivity, secondary-structure predictions, and 2D diagrams.

There is no wrapper script here. `SKILL.md` teaches the agent the samplesheet rules, reference-routing logic, and known gotchas needed to construct a correct `nextflow run` invocation directly.

## Scope

- Upstream RNA structure-probing analysis via Nextflow: RNA Framework (rf-count/rf-norm/rf-fold/…), STAR/Bowtie alignment, ViennaRNA, R2DT.
- Genome route (STAR, default) or transcriptome route (Bowtie/Bowtie2, `--transcriptome`).
- Local FASTA/GTF, or automatic Ensembl/NCBI reference download by `organism`.

## Out of Scope

- Cross-condition statistical comparison of reactivity/structure output (the pipeline itself doesn't do this either).
- Summarising or plotting results after the run — no downstream ClawBio skill consumes this output yet.

## Quick Start

```bash
nextflow run nf-core/rnastructurome -r 1.0.0 -profile test,docker --outdir ./rnastructurome_demo
```

For real data:

```bash
nextflow run nf-core/rnastructurome -r 1.0.0 \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results
```

See `SKILL.md` for the samplesheet column rules, reference-routing decision, and full parameter reference.
