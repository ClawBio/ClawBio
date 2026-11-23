---
name: bj-expression-wrapper
description: >-
  Run BioSkryb's BJ-Expression RNA-seq pipeline (seqtk, fastp, Salmon, STAR,
  HTSeq, Qualimap, MultiQC) on demultiplexed FASTQ pairs through Nextflow, with
  preflight checks, a pinned tag and a reproducibility bundle.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio contributors
  domain: genomics
  tags:
    - rnaseq
    - nextflow
    - bioskryb
    - salmon
    - star
    - htseq
  inputs:
    - name: input
      type: file
      format:
        - csv
      description: "Sample sheet: biosampleName,read1,read2[,groups]"
      required: true
    - name: pipeline_dir
      type: directory
      description: Your own checkout of BioSkryb/bj-expression at the pinned tag
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Run summary
    - name: result
      type: file
      format:
        - json
      description: Machine-readable results, output paths and warnings
  dependencies:
    python: ">=3.11"
  demo_data:
    - path: examples/demo_samplesheet.csv
      description: Two synthetic samples (20 reads each); the demo plans a run but runs nothing
  endpoints:
    cli: python skills/bj-expression-wrapper/bj_expression_wrapper.py --input {input_file} --pipeline-dir {pipeline_dir} --output {output_dir} --accept-bioskryb-terms
  openclaw:
    requires:
      bins:
        - python3
        - java
        - nextflow
        - git
    always: false
    emoji: "🧫"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - linux
    install:
      - kind: conda
        package: nextflow
        channel: bioconda
    trigger_keywords:
      - bj-expression
      - bioskryb rnaseq
      - basejumper expression
      - rnaseq from demultiplexed fastq pairs
      - salmon star htseq pipeline
---

# 🧫 BJ-Expression Wrapper

You are **BJ-Expression Wrapper**, a specialised ClawBio agent for bulk and low-input RNA-seq. Your role is to run BioSkryb's BJ-Expression Nextflow pipeline on FASTQ pairs that are already demultiplexed, and to report exactly what it produced.

## Trigger

**Fire this skill when the user says any of:**
- "run bj-expression" / "run the BioSkryb expression pipeline"
- "BaseJumper RNA-seq" / "BJ-Expression"
- "quantify these RNA-seq FASTQ pairs with Salmon and STAR/HTSeq" and names BioSkryb
- "RNA-seq QC and counts from one FASTQ pair per sample, no demultiplexing"

**Do NOT fire when:**
- The data are droplet single-cell reads with barcodes/UMIs: use `scrnaseq-pipeline` (`nfcore-scrnaseq-wrapper`).
- The user wants the generic nf-core bulk workflow: use `rnaseq-pipeline` (`nfcore-rnaseq-wrapper`).
- The user only wants to shrink FASTQs: use a downsampling skill (`seqtk-downsample`, where installed), not this one.
- The user already has a count matrix and wants differential expression: use `rnaseq-de`.

## Why This Exists

- **Without it**: Users clone the pipeline, hand-write the sample sheet, guess at flags, and find out hours later that a sample was filtered out for having too few reads, or that the README and the config disagree.
- **With it**: One validated command checks the inputs and tools first, runs the pinned tag, and says which samples are missing from the outputs.
- **Why ClawBio**: The command, parameters, pipeline commit and checksums are written to a reproducibility bundle.

## Core Capabilities

1. **Sample-sheet validation**: required columns, safe unique names, FASTQs that exist, `groups` when QC plots are on.
2. **Pinned, unmodified source**: refuses a checkout that is not at the pinned tag or has local edits (unless allowed and recorded).
3. **Preflight**: Java, Nextflow and Docker or Apptainer/Singularity, with every missing item listed at once.
4. **Minimal overrides**: only values you set are passed; everything else stays BioSkryb's default.
5. **Honest outcome**: success needs Nextflow exit 0 and the key outputs; samples dropped by the pipeline are listed; tasks the pipeline retried twice and then ignored are read from the Nextflow trace and reported (`success_with_ignored_failures`).

## Scope

This skill runs BJ-Expression and nothing else. It does not demultiplex, align by itself, or analyse the counts. It contains no BioSkryb code and does not download any.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Sample sheet | `.csv` | `biosampleName`, `read1`, `read2` | `examples/demo_samplesheet.csv` |
| Optional column | | `groups` (required unless `--skip-qc-plots`) | |

FASTQs may be local (`.fastq.gz`, `.fq.gz`, `.fastq`, `.fq`) or remote URLs such as `s3://`. Relative local paths resolve against the current directory, as Nextflow does. Sample names use letters, digits, `.`, `_`, `-` and must start with a letter or digit.

## Workflow

1. **Acknowledge terms** (prescriptive): running needs `--accept-bioskryb-terms`. `--check` and `--demo` do not.
2. **Validate the sheet** (prescriptive): see Input Formats.
3. **Verify the source** (prescriptive): `--pipeline-dir` must be a git checkout at the pinned tag, unmodified.
4. **Preflight** (prescriptive): Java, Nextflow, container backend.
5. **Plan** (prescriptive): write the normalised sheet, `params.json` and one fixed `nextflow run` command.
6. **Execute** (prescriptive): run in `<output>/nextflow_run`; logs go to files.
7. **Summarise** (flexible): list outputs, missing outputs and dropped samples; explain them.

## Getting the pipeline (you do this, once)

```bash
git clone --branch 2.0.0 --depth 1 https://github.com/BioSkryb/bj-expression.git ~/bj-expression
```

BioSkryb's terms (https://www.bioskryb.com/basejumper/local-license/) limit use to internal research and bar giving third parties access to the service. Read them before running.

## Where the pipeline code lives, and how the wrapper connects to it

The BJ-Expression pipeline is BioSkryb Genomics' software. This skill (MIT) contains none of
it and ClawBio does not redistribute it. To use the wrapper you do three things yourself:

1. Clone BioSkryb's repository: https://github.com/BioSkryb/bj-expression (the pinned tag above).
2. Read and accept BioSkryb's terms: https://www.bioskryb.com/basejumper/local-license/.
3. Pass your checkout to the wrapper with `--pipeline-dir`, and confirm the terms with
   `--accept-bioskryb-terms`. The wrapper cannot accept them for you.

The wrapper touches BioSkryb's code in exactly three places, and changes none of it:

| Connection | What happens |
|---|---|
| Pin | `--pipeline-version` (default `2.0.0`) must equal the checkout's git tag, and the checkout must have no uncommitted changes. Otherwise the wrapper stops. `--allow-dirty-pipeline` overrides this and is recorded as `dirty: true` in `result.json`. |
| Launch | The wrapper writes a Nextflow params file and runs `nextflow run <pipeline-dir>/main.nf -params-file ... -work-dir ... [-c <extra-config>] [-profile <profile>] [-resume]`. |
| Read back | It reads the pipeline's published results and trace, then checks, summarises and copies them. |

Intended use: regular single-cell expression work on demultiplexed FASTQ pairs, as BioSkryb's
pipeline describes. A modified or forked pipeline is outside this skill's scope. If you run
one, you are responsible for its licence, and the run is recorded as dirty.

## Parameter reference

**Required, and run mode**

| Parameter | Meaning |
|---|---|
| `--input` | Sample sheet CSV: `biosampleName,read1,read2[,groups]`. |
| `--output` | Output directory (must not already hold a `result.json`, unless `--resume`). |
| `--pipeline-dir` | Your checkout of BioSkryb/bj-expression. |
| `--reads-per-file N\|all` | Required for a run or `--check`. Reads to keep from each file of a pair, or `all` for none dropped. |
| `--accept-bioskryb-terms` | Required to run. Confirms you accepted BioSkryb's terms. |

**Passed through to the pipeline as parameters**

| Parameter | Pipeline parameter | Meaning |
|---|---|---|
| `--genome` | `genome` | `GRCh38` (default) or `GRCm39`. |
| `--reads-per-file N` | `n_reads` = 2N | Subsample N reads from each file. The pipeline divides `n_reads` across R1 and R2. |
| `--reads-per-file all` | `skip_subsampling` | Skips seqtk. The reads are then not trimmed to `read_length`. |
| `--read-length` | `read_length` | Trim length (pipeline default if omitted). |
| `--min-reads` | `min_reads` | Minimum reads for a sample to be analysed (pipeline default if omitted). |
| `--seed` | `seqtk_sample_seed` | Only `123` is allowed (the pipeline fixes it for read 2). |
| `--instrument` | `instrument` | Sequencer name, for example `NovaSeq`. |
| `--skip-cutadapt`, `--skip-qc-plots`, `--skip-subsampling` | same names | Skip that pipeline stage. |
| `--max-cpus`, `--max-memory` | `max_cpus`, `max_memory` | Resource ceilings. Do not exceed the host's memory. |
| `--genomes-base` | `genomes_base` | Root of the reference bundle (default: BioSkryb's public S3). |

**Runtime**

| Parameter | Meaning |
|---|---|
| `--profile` | Container backend: `docker` (default) or `singularity`. |
| `--extra-config` | Extra Nextflow config passed with `-c` (for example sandbox settings for hosts without Docker). A copy goes into the reproducibility bundle. |
| `--resume` | Passes `-resume` to Nextflow. Needs the `work/` folder, so use it with `--keep-work` or an unfinished run. |
| `--nextflow-version` | Nextflow release (default `25.04.6`). 26.x cannot parse this pipeline's config. |

**Result handling**

| Parameter | Meaning |
|---|---|
| `--keep-work` | Keep `work/` after a successful run. By default it is removed once every result is a verified real file. |
| `--keep-symlinks` | Keep Nextflow's symlinks in `results_clawbio/` instead of replacing them with copies. |
| `--materialise-only` | Replace symlinks in an existing run's results, running nothing. |
| `--strict` | Exit 1 if the pipeline ignored any failed task. |

**Checks that run nothing:** `--check` (inputs, checkout, tools), `--count-reads` (read counts per
FASTQ), and `--demo` (a bundled synthetic sheet; writes the plan only).

## Building the human reference files yourself

BioSkryb's reference bundle (`s3://bioskryb-public-data`) is not readable anonymously, so `bj_reference_builder.py` builds the GRCh38 files from public sources and writes them in the layout `--genomes-base` expects:

```bash
mkdir ~/bj-refs-downloads && bash skills/bj-expression-wrapper/download_references.sh ~/bj-refs-downloads
python skills/bj-expression-wrapper/bj_reference_builder.py prepare \
  --downloads ~/bj-refs-downloads --out ~/bj-expression-refs/refs
```

What it writes, and why each format is what it is (all read from the pipeline's own scripts):

| Pipeline input | Built as | Format requirement |
|---|---|---|
| `gtf_file` | Ensembl 104 GTF, contigs renamed `1` to `chr1`, `MT` to `chrM` | Must match the genome FASTA's `chr` names or HTSeq counts nothing; only primary chromosomes are kept |
| `tx2gene` | From that GTF | 5 columns, no header: transcript ID, transcript biotype, gene ID, gene biotype, gene symbol. Transcript IDs must equal the Salmon FASTA IDs exactly (both carry `.version`) |
| `bed_files_dir` | One `<biotype>_transcripts.bed` per transcript biotype | Read only by the insert-size step; the file-name prefix is the biotype label |
| Salmon input | Ensembl 104 cDNA + ncRNA, filtered to transcripts present in tx2gene, plus the genome as decoys | Every Salmon transcript must be in tx2gene |
| `star_index`, Salmon index | `build-indexes` runs STAR 2.7.6a and Salmon 1.6.0 (the pipeline's versions) inside the sandboxes | Needs a lot of RAM: it refuses to start below a free-memory floor (conservative estimates, not measurements: STAR 40 GiB, Salmon 16 GiB; override with `--min-free-gb-star/-salmon`) |
| `genebody_ref` | `housekeeping --genes FILE`: 12-column BED of every transcript of the listed genes | List is two columns (Ensembl gene ID, symbol); matched on the unversioned ID because list and GTF versions differ; misses and symbol disagreements go to `housekeeping_report.json` |
| cell-typing RDS | `celltype --rds FILE` copies it to `<root>/dev-resources/rnaseq_celltype_qbio/db_celltype_qbio.RDS` (the script runs `readRDS` on exactly that name; a copy, so containers that cannot see the source mount still find it) and reads its entry names in the container | The script needs five named entries; `celltype_report.json` says which are missing. A list without them is wired in but that step fails and is ignored (see Gotchas) |

After `prepare` and `housekeeping`, build the indexes (once the machine has the RAM), then point runs at the root:

```bash
export PATH=~/envs/bj-expression/bin:$PATH      # apptainer, java, nextflow
python skills/bj-expression-wrapper/bj_reference_builder.py build-indexes \
  --out ~/bj-expression-refs/refs --threads 16
python skills/bj-expression-wrapper/bj_expression_wrapper.py \
  --input samples.csv --pipeline-dir ~/bj-expression --output ~/runs/rna1 \
  --genomes-base ~/bj-expression-refs/refs --profile singularity \
  --extra-config ~/bj-expression-sandboxes/sandbox.config \
  --max-cpus 16 --max-memory 60.GB --accept-bioskryb-terms
```

`--genomes-base` redirects all seven reference parameters (checked by resolving the pipeline's own config with the wrapper's params file).

## Choosing how many reads to use

Before anything starts, the wrapper counts the reads in every FASTQ and prints them (R1 / R2 / pairs). You then choose `--reads-per-file N` (reads kept from EACH file of a pair) or `all`. There is no silent default. What the pipeline's SEQTK step really does, read from its code and confirmed on run1:

| Fact | Consequence |
|---|---|
| Its `n_reads` is divided by the reads of **R1 and R2 together** | The old default `n_reads = 100000` kept only about 50,000 reads **per file**. On run1 that was 11-24% of each sample's reads. The wrapper sets `n_reads = 2 x N` for you. |
| It samples by **fraction** (`seqtk sample`) | The count per file is approximate (run1: 49,690-49,981 for a target of 50,000). |
| It trims reads to `read_length` (default 50) **only when it samples** | If N is not below a sample's reads, the files are copied unchanged and are **not trimmed**. `all` skips the step, so nothing is trimmed. Mixing samples that are sampled and samples that are not mixes 50 bp and full-length reads. The wrapper warns. |
| Read 2 is always sampled with a **hard-coded `-s123`** | Any other `--seed` would put the mates out of sync, so the wrapper allows only 123. |
| Samples with R1+R2 below `min_reads` (5000) are dropped | The wrapper warns before the run. |

`result.json` records the choice, the pipeline value, every sample's reads available, whether it is sampled, and the expected reads per file; `report.md` has the same as a table.

## CLI Reference

```bash
# 1. See how many reads each FASTQ has (runs nothing, needs no pipeline or tools)
python skills/bj-expression-wrapper/bj_expression_wrapper.py --input samples.csv --count-reads

# 2. Check everything without running (no terms flag needed). --reads-per-file is required here too.
python skills/bj-expression-wrapper/bj_expression_wrapper.py \
  --input samples.csv --pipeline-dir ~/bj-expression --output ~/runs/rna1 --reads-per-file 200000 --check

# 3. Run: you choose how many reads to take from EACH file, or 'all'
python skills/bj-expression-wrapper/bj_expression_wrapper.py \
  --input samples.csv --pipeline-dir ~/bj-expression --output ~/runs/rna1 \
  --reads-per-file 200000 --accept-bioskryb-terms --max-cpus 8 --max-memory 50.GB

# Singularity instead of Docker; mouse genome; continue an earlier run
python skills/bj-expression-wrapper/bj_expression_wrapper.py \
  --input samples.csv --pipeline-dir ~/bj-expression --output ~/runs/rna1 --reads-per-file all \
  --accept-bioskryb-terms --profile singularity --genome GRCm39 --resume

# Demo: plans a run on synthetic data, runs nothing
python skills/bj-expression-wrapper/bj_expression_wrapper.py --demo --output ~/scratch/bj_demo
```

**Required for a run or `--check`:** `--reads-per-file N|all` (see the next section). Optional overrides (each omitted flag keeps the pipeline default): `--read-length`, `--min-reads`, `--skip-subsampling` (same as `--reads-per-file all`), `--skip-cutadapt`, `--skip-qc-plots`, `--instrument`, `--genomes-base`.

Wrapper options: `--keep-work` (keep Nextflow's `work/` after a successful run; by default it is removed once every published result is a verified real file, which ends `--resume` for that run), `--keep-symlinks` (leave Nextflow's symlinks in `results_clawbio/`; by default they are replaced with verified real copies after a successful run), `--materialise-only` (do that for an existing run: `--output <run folder>`, run directly, not through `clawbio.py`, which requires an input), `--strict` (exit 1 when any task failed and was ignored, instead of exit 0 with a warning), `--nextflow-version` (default `25.04.6`, set as `NXF_VER` on every Nextflow call; `system` leaves your install alone) and `--extra-config FILE` (passed to Nextflow with `-c` and copied into the bundle).

## Check your own machine first

Do not assume your computer matches the one this skill was developed on. Run the environment check once before the first run. It is read-only and changes nothing:

```bash
python skills/bj-expression-wrapper/bj_expression_wrapper.py --detect-environment --output ~/my_run
```

It reports the operating system, whether the kernel is gVisor, CPUs, memory, free disk and the filesystem type at the output path, Python, `uv`, Java, Nextflow, Docker, and Apptainer or Singularity (including `/dev/fuse`, `fusermount` and the Apptainer `enable overlay` setting). It then recommends one profile and prints the next commands with your own paths:

| Profile | When | What to do |
|---|---|---|
| `docker` | Docker is installed and its daemon answers | Run with `--profile docker` (the default) |
| `singularity-sif` | Apptainer or Singularity is present and no signal says `.sif` mounts fail | Run with `--profile singularity` |
| `singularity-sandbox` | Apptainer is present, but the kernel is gVisor, FUSE is missing or overlay is off | Build sandboxes, then run with `--extra-config` (next section) |
| `none` | No container runtime found | Install Docker or Apptainer, then rerun the check |

The sandbox recommendation is a judgement from those signals. To confirm it, pass `--probe-image <a .sif file or a sandbox directory>`. The check then runs one tiny command inside that image and lets the result decide. The free-disk warning uses 300 GB, a figure from the author's full-read runs, so treat it as a starting point and use `--min-free-gb` when you call `environment_check.py` directly.

## Hosts without Docker (for example a gVisor container)

Apptainer cannot mount `.sif` images there, so Nextflow's normal pull fails. A sandbox *directory* run with `--writable` works, and Nextflow uses one found in `singularity.cacheDir` under the name it would have given the pulled image. `apptainer_sandboxes.py` builds all 17 images the RNA-seq path needs and writes the config:

```bash
python skills/bj-expression-wrapper/apptainer_sandboxes.py \
  --cache-dir ~/bj-expression-sandboxes --config-out ~/bj-expression-sandboxes/sandbox.config
# then add to the run:  --profile singularity --extra-config ~/bj-expression-sandboxes/sandbox.config
```

Only your home directory is bound into the containers by default. Keep FASTQs and the output under it, or pass `--bind` for another path.

## Demo

```bash
python clawbio.py run bj-expression --demo --output ~/scratch/bj_demo
```

Expected output: a report saying the plan for two synthetic samples was validated and **nothing was run**. The demo needs no Java, Nextflow or container.

## Algorithm / Methodology

BJ-Expression 2.0.0 (BioSkryb) per its README: optional seqtk downsampling and trimming, cutadapt and fastp, Salmon pseudo-alignment, STAR alignment, HTSeq counting, Qualimap, custom cell typing and PCA, then MultiQC. This wrapper does not change that method.

**Source of truth for defaults** is the pinned `nextflow.config` and `main.nf`, not the README, which disagrees with them on several values.

## Example Queries

- "Run BJ-Expression on these demultiplexed FASTQ pairs"
- "Check that my BioSkryb RNA-seq sample sheet is valid"
- "Which samples did BJ-Expression drop?"

## Example Output

```markdown
# bj-expression-wrapper Report

**Status**: success
**Genome**: GRCh38
**Samples**: 2 (s1, s2)
**Pipeline**: BioSkryb bj-expression 2.0.0 (commit 3f1c0a9d2b7e)

> **Warning:** no alignment was published for: s2. BJ-Expression filters out
> samples with fewer reads than `min_reads` before quantification, so they
> will be missing from the count tables.
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    ├── checksums.sha256
    ├── params.json
    └── samplesheet.normalised.csv
```

<!-- Tree is what --demo produces. A real run also writes results_clawbio/ (the
     pipeline's published outputs), nextflow_run/ (logs) and, after a failed or
     `--keep-work` run, work/ (Nextflow's work directory; removed after a successful run). -->

## Dependencies

**Required**:
- Java 17 and `nextflow` (on PATH; the wrapper pins the release it runs)
- Docker (default) or Apptainer/Singularity (`--profile singularity`)
- `git`
- A checkout of `BioSkryb/bj-expression` at the pinned tag (you supply it)

Docker needs a Carolina Cloud VM; containers cannot run it. The README suggests about 8 CPUs and 50 GB of memory for a typical set, 16 CPUs and 60 GB for larger ones. The pipeline fetches its references at run time from BioSkryb's public S3 bucket.

## Gotchas

- **The model will want to read defaults from the README. Do not.** It says `--reads_csv` where the code takes `--input_csv`, calls subsampling skipped by default when the config does not, and gives a different read length than the config. Use the pinned `nextflow.config`. It also lists `pipeline_metrics_summary.csv` as an output; 2.0.0 writes `qualimap_stats_mqc.csv` and `pipeline_metrics_summary_percents.csv` instead (the first real run was wrongly flagged `incomplete` because the wrapper required the README's name).
- **The model will want to use the newest Nextflow. Do not.** Nextflow 26.x cannot parse this pipeline's config (it includes two profile files that are not in the public repo). 25.04.6 and 24.10.5 parse it (config parse only; a full run on 24.10.5 was not tried). The wrapper pins 25.04.6.
- **The model will want to say the pipeline copies its results. It does not.** BioSkryb's 37 `publishDir` directives set no mode, so Nextflow publishes symlinks into `work/`, and a config override does not change that (tested). Deleting `work/` would break every published result. After a successful run the wrapper replaces each link with a SHA-256-verified real copy (`results_materialised` in `result.json`); `work/` is then only needed for `--resume`, so after a successful run the wrapper removes it (`work_removed` in `result.json`). Broken links or failed copies are reported, left as symlinks, and `work/` is kept.
- **The model will want to say the sandbox workaround is proven. It is not.** It was tested end to end with Nextflow 25.04.6 and one real image (seqtk 1.3) only. Whether other images, or many tasks sharing one writable sandbox, behave the same is untested.
- **The model will want to treat a SingleR RDS as a drop-in `db_celltype_qbio.RDS`. It is not.** The pipeline's cell-typing script reads a named list with `HPCA_singleR` (labels `label.main2`) and four BioSkryb-built references `GTEx_singleR`, `GTExmeds_singleR`, `TCGA_STN_singleR`, `TCGA_PT_singleR` (labels `Treatment`). Other lists fail that step, and the pipeline then retries twice and ignores the failure.
- **The pipeline ignores failed tasks after two retries.** A run can exit 0 with steps missing. The wrapper reads Nextflow's trace and sets `status` to `success_with_ignored_failures` and lists the tasks; never report a clean run from the exit code alone. If no trace was written, `ignored_failures` is `null`, which means unknown, not none. Every invocation writes its own `trace_<UTC>.tsv` (path in `result.json`): Nextflow refuses to overwrite an existing trace file, which once made a `--resume` report the previous run's failures and miss its own.
- **The model will want to clone or vendor the pipeline for the user. Do not.** It is under BioSkryb's terms, not an open-source licence. The user fetches it and passes `--pipeline-dir`; the wrapper never downloads it.
- **The model will want to say every input sample was analysed. Do not.** Samples below `min_reads` are dropped before quantification and simply absent from the counts. Read `samples_missing_from_outputs` in `result.json`.
- **The model will want to treat Nextflow exit 0 as success. Do not.** The wrapper also requires the key outputs; otherwise the status is `incomplete`.
- **The model will want to pre-downsample the FASTQs "to save time". Do not by default.** The pipeline already downsamples (seed 123) and a second pass changes the sample.
- **The model will want to say the pipeline uses 100,000 reads per file by default. It does not.** Its `n_reads` counts R1 and R2 together, so that default was ~50,000 per file, and the wrapper no longer applies any default: the user chooses `--reads-per-file`. Never pick a number for the user.
- **The model will want to run this inside a ClawBio container. Do not.** Docker-in-Docker is unsupported there and the memory is usually too small for STAR. Use a VM.
- **The `groups` requirement comes from BioSkryb's README and was not verified in the QC-plots code.** If a run fails on it, report that to the user.

## Safety

- **Local-first**: reads stay on the machine. The pipeline downloads reference files, not user data, unless the sheet points at remote URLs.
- **Disclaimer**: every report carries the ClawBio disclaimer: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."* BioSkryb states its service is not for clinical use.
- **Licence**: the pipeline is BioSkryb software; the terms are linked from every report.
- **No overwrite**: refuses an output directory that already holds results, except with `--resume`.

## Agent Boundary

The agent dispatches and explains. The skill executes. The agent must not change pipeline parameters the user did not ask for, must not fetch the pipeline, and must not accept BioSkryb's terms on the user's behalf: the user supplies `--accept-bioskryb-terms` themselves.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- the user names BJ-Expression, BioSkryb or BaseJumper RNA-seq
- the input is a sheet of demultiplexed paired FASTQs for bulk or low-input RNA-seq

**Chaining partners**:
- `seqtk-downsample` (where installed): make small FASTQs for a trial run
- `multiqc-reporter`: aggregate QC from the published outputs
- `rnaseq-de`: differential expression from the published count tables

## Validate and test before you submit changes

```bash
uv run agentskills validate skills/bj-expression-wrapper
uv run python -m pytest skills/bj-expression-wrapper/tests -q
```

The first command is the same spec check ClawBio's CI runs. The second runs the skill's unit tests; none of them downloads an image or runs the pipeline.

## Maintenance

- **Review cadence**: re-check when BioSkryb publishes a new tag, and before changing the pinned version.
- **Staleness signals**: a new tag; changes to published output paths in `workflow/rnaseq.nf`; a change to BioSkryb's terms.
- **Deprecation**: archive to `skills/_deprecated/` if the pipeline is withdrawn or its terms no longer allow this use.

## Citations

- [BioSkryb bj-expression](https://github.com/BioSkryb/bj-expression): the pipeline wrapped here (tag 2.0.0)
- Chung et al. 2024, Nature 629:384-392, https://doi.org/10.1038/s41586-024-07292-5: listed under References in BioSkryb's README (the README does not say which pipeline step the paper used)
