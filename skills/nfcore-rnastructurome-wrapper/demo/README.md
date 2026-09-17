# Demo

This skill has no bundled sample data. The demo runs the upstream pipeline's own `test` profile,
which fetches its FASTQs and reference from `nf-core/test-datasets` (`rnastructurome` branch) over HTTPS.

```bash
nextflow run nf-core/rnastructurome -profile test,docker --outdir ./rnastructurome_demo
```

Reference: human mitochondrial chromosome (16,569 bp). Reads: synthetic, from ENST00000389680
(MT-RNR1, MT:648-1601), guaranteed to map. Exercises the full STAR genome-alignment route on
data small enough to finish in minutes.

Other bundled upstream test profiles:
- `test_transcriptome` — exercises the Bowtie/Bowtie2 transcriptome route instead of STAR
- `test_prokaryote` — bacterial reference, NCBI fallback resolution
- `test_full` — full-size dataset (not for routine demoing)

Requires network access (the test data is remote, not bundled) and a container backend
(`docker`, `singularity`, or `conda` in place of `docker` in the profile string above).
