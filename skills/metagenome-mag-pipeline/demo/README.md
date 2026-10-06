# Demo

The demo has no bundled input file. At run time it clones the pinned upstream
checkout (or uses `--pipeline-dir`) and asks upstream's own generator for
synthetic reads:

```
python3 <pipeline>/test/make_test_dataset.py \
  --outdir <output>/pipeline/demo_data --dbdir <output>/pipeline/demo_data/dbs \
  --samples 2 --genomes 4 --genome-length 12000 --pairs 2000 --seed 1
```

It then runs the upstream test configuration (`workflow/config/config_test.yaml`):
adapter trimming, quality filtering, metaSPAdes assembly, QUAST, MultiQC and the
provenance records. Profiling, binning, MAG assessment, contig taxonomy, host
decontamination and sourmash are all off in that configuration, because each of
them needs a reference the test does not ship.

## What it needs

- `snakemake` ≥ 9.0 on PATH (upstream pins 9.11.2)
- `conda` or `mamba` on PATH — every tool environment is built on first run
- network access to conda-forge/bioconda to build those environments
- `git` and network access to clone the workflow, unless `--pipeline-dir` is given

Add `--check` to validate the plan and write a planned report without executing
anything, which is the useful first step on a machine that has no conda yet.

## What it proves

That the wrapper can drive the workflow end to end: validate a samplesheet,
render a config, launch snakemake, and parse the tables it wrote.

It does **not** prove anything about accuracy. The fragments are drawn from
random sequences with locally elevated GC, so they have no taxonomy and no
marker genes. Consensus binning in particular cannot complete on synthetic input,
because MaxBin2 and DAS Tool need single-copy marker genes. Any report from this
demo says so in its caveats.
