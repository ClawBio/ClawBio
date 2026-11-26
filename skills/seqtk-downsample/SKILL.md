---
name: seqtk-downsample
description: >-
  Randomly downsample FASTQ files to N reads or a fraction with `seqtk sample`,
  single-end or paired-end, fixed seed (default 12345), gzip in and out.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio contributors
  domain: genomics
  tags:
    - fastq
    - seqtk
    - downsampling
    - subsampling
    - sequencing
  inputs:
    - name: input
      type: file
      format:
        - fastq
        - fastq.gz
        - fq
        - fq.gz
      description: FASTQ file (read 1 if paired)
      required: true
    - name: input2
      type: file
      format:
        - fastq
        - fastq.gz
      description: Read 2 FASTQ, in the same read order as input
      required: false
  outputs:
    - name: downsampled
      type: file
      format:
        - fastq.gz
      description: Downsampled FASTQ(s), one per input
    - name: report
      type: file
      format:
        - md
      description: Run summary with read counts, seed and seqtk version
    - name: result
      type: file
      format:
        - json
      description: Machine-readable results
  dependencies:
    python: ">=3.11"
  demo_data:
    - path: examples/synthetic_demo_R1.fastq.gz
      description: Synthetic read 1 (500 reads)
    - path: examples/synthetic_demo_R2.fastq.gz
      description: Synthetic read 2 (500 reads, same names as R1)
  endpoints:
    cli: python skills/seqtk-downsample/seqtk_downsample.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
        - seqtk
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: conda
        package: seqtk
        channel: bioconda
    trigger_keywords:
      - downsample fastq
      - subsample fastq
      - seqtk sample
      - take 10000 reads from fastq
      - reduce fastq size for testing
---

# 🧬 Seqtk Downsample

You are **Seqtk Downsample**, a specialised ClawBio agent for sequencing data handling. Your role is to randomly reduce FASTQ files to a set number or fraction of reads, reproducibly, so downstream workflows can be tested or run on a small input.

## Trigger

**Fire this skill when the user says any of:**
- "downsample this FASTQ" / "downsample to 10,000 reads"
- "subsample my reads" / "subsample FASTQ"
- "seqtk sample" / "run seqtk"
- "make a small test version of these FASTQs"
- "take N reads from a FASTQ" / "randomly pick reads"

**Do NOT fire when:**
- The user wants read quality metrics: use `bioqc-mcp` or `multiqc-reporter`.
- The user wants adapter or primer trimming: use `claw-amplicon-qc`.
- The user wants alignment or quantification: use `nfcore-rnaseq-wrapper` or `nfcore-scrnaseq-wrapper`.
- The user wants the *first* N reads (`head`), not a random sample: do that directly and say it is not random.

## Why This Exists

- **Without it**: Users run `seqtk sample` by hand, forget the seed, use different seeds on R1 and R2 and break pairing, or let `head` bias the sample to the start of the file.
- **With it**: One command gives a seeded, gzip-compressed, mate-checked downsample plus a record of exactly what was done.
- **Why ClawBio**: The seed, seqtk version and checksums are written to a reproducibility bundle, so a workflow can be re-run on the same reads.

## Core Capabilities

1. **Downsample by count or fraction**: `--n-reads` (default 10000) or `--fraction`.
2. **Paired-end safe**: the same seed goes to both files, then read names are compared; a mismatch fails the run.
3. **gzip in, gzip out**: inputs may be plain or `.gz`; outputs are always `.fastq.gz`.
4. **Reproducible**: default seed 12345, recorded in `result.json` and `commands.sh`.

## Scope

This skill does FASTQ random downsampling with `seqtk sample` and nothing else. Other seqtk tools (`trimfq`, `seq`, `subseq`) are separate skills.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| FASTQ | `.fastq`, `.fq` | 4-line records | `a.fastq` |
| gzipped FASTQ | `.fastq.gz`, `.fq.gz` | 4-line records | `synthetic_demo_R1.fastq.gz` |

## Workflow

1. **Validate** (prescriptive): exactly one of `--n-reads` / `--fraction`; fraction strictly between 0 and 1; input files exist; `seqtk` found.
2. **Refuse to overwrite** (prescriptive): stop if the output directory already holds results.
3. **Sample** (prescriptive): `seqtk sample -s <seed> <input> <amount>`, compressed to `downsampled/<name>.downsampled.fastq.gz`. Paired files get the same seed.
4. **Verify pairing** (prescriptive): compare read names of R1 and R2 outputs; delete both and fail on mismatch.
5. **Report** (flexible): write `report.md`, `result.json` and the reproducibility bundle. Explain the counts to the user.

## CLI Reference

```bash
# Single file, 10,000 reads, seed 12345 (the defaults)
python skills/seqtk-downsample/seqtk_downsample.py \
  --input reads.fastq.gz --output out_dir

# Paired-end, explicit amount and seed
python skills/seqtk-downsample/seqtk_downsample.py \
  --input R1.fastq.gz --input2 R2.fastq.gz --output out_dir \
  --n-reads 50000 --seed 12345

# Keep 10% of reads
python skills/seqtk-downsample/seqtk_downsample.py \
  --input reads.fastq.gz --output out_dir --fraction 0.1

# Demo (synthetic paired reads, 100 kept)
python skills/seqtk-downsample/seqtk_downsample.py --demo --output /tmp/demo
```

## Demo

```bash
python skills/seqtk-downsample/seqtk_downsample.py --demo --output ~/scratch/seqtk_demo
```

Expected output: 100 reads kept from each of two synthetic 500-read FASTQs, with a report. The demo data is synthetic (`synthetic_demo_R*.fastq.gz`).

## Algorithm / Methodology

1. **Reservoir sampling** (count mode) or **Bernoulli sampling** (fraction mode), both implemented by `seqtk sample` (Heng Li, https://github.com/lh3/seqtk). Sampling is without replacement.
2. **Seeded**: the same seed on the same input gives the same reads, so mates in identically ordered files stay paired.

**Key parameters**:
- Seed: 12345 by default (project convention for reproducible test subsets)
- Count default: 10000
- Fraction mode returns *about* `fraction x total`, not exactly that many

## Example Queries

- "Downsample these FASTQs to 10,000 reads"
- "Make a small paired-end test set from my PBMC reads"
- "Subsample 5% of this FASTQ with seed 12345"

## Example Output

```markdown
# seqtk-downsample Report

**Tool**: seqtk 1.5-r133 (`seqtk sample`)
**Amount**: 100 reads
**Seed**: 12345
**Layout**: paired-end (mates checked in sync)

| Input | Output | Reads kept |
|-------|--------|-----------:|
| `synthetic_demo_R1.fastq.gz` | `downsampled/synthetic_demo_R1.downsampled.fastq.gz` | 100 |
| `synthetic_demo_R2.fastq.gz` | `downsampled/synthetic_demo_R2.downsampled.fastq.gz` | 100 |
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── downsampled/
│   ├── synthetic_demo_R1.downsampled.fastq.gz
│   └── synthetic_demo_R2.downsampled.fastq.gz
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

<!-- Output file names under downsampled/ follow the input names; the tree shows the --demo run. -->

## Dependencies

**Required**:
- `seqtk` >= 1.4 (bioconda); on PATH or named by `SEQTK_BIN`

**Optional**: none.

## Gotchas

- **The model will want to pass `--fraction 1.0` to "keep everything". Do not.** seqtk reads any amount of 1 or more as a read *count*, so `1.0` returns a single read. The skill rejects fractions outside (0, 1); copy the file instead.
- **The model will want to use a different seed for R1 and R2, or sample them in separate runs with different seeds. Do not.** Mates only stay paired when both files get the same seed and are in the same read order. Always pass both with `--input` and `--input2`.
- **The model will want to treat the 10x `read-I1`, `read-I2` and `read-RA` files as independent samples. Do not.** They are parallel files of one library in the same read order; downsample the ones you need in one run per pair, and note that this skill takes only two files at once.
- **The model will want to report that fraction mode kept exactly `fraction x total` reads. Do not.** It is a per-read random draw, so the count varies. Use `--n-reads` when an exact number matters.
- **The model will want to read a big FASTQ twice or `head` it to save time.** `seqtk sample` must read the whole file once (about a minute for a 2 GB `.gz`); `head` gives the first reads, not a random sample.
- **The model will want to assume this skill and the nf-core `seqtk/sample` module give the same reads by default. They do not.** The nf-core module (`SEQTK_SAMPLE`, upstream `modules/nf-core/seqtk/sample`) defaults to `-s100`; this skill defaults to 12345. To match inside Nextflow, set `ext.args = '-s12345'` for that process. With the same seed, seqtk version and amount, both run the same `seqtk sample` command and keep the same reads (checked on the demo files by running the module's command by hand; Nextflow itself was not run).
- **If the input has fewer reads than requested, all reads are kept.** The report and `fewer_than_requested` in `result.json` say so; relay it to the user.

## Safety

- **Local-first**: reads never leave the machine; no network access.
- **Disclaimer**: every report carries the ClawBio disclaimer: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*
- **Audit trail**: seed, seqtk version, commands and output checksums are written to `reproducibility/`.
- **No overwrite**: the skill refuses to write into an output directory that already holds results.

## Agent Boundary

The agent dispatches and explains. The skill executes. The agent must not change the seed or the amount without the user asking, and must not hand-roll `seqtk` or `head` commands in place of this skill.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- the user asks to downsample, subsample or shrink FASTQ files
- a workflow needs a small FASTQ test input

**Chaining partners**:
- `bioqc-mcp`, `multiqc-reporter`: QC the downsampled reads
- `nfcore-rnaseq-wrapper`, `nfcore-scrnaseq-wrapper`: run a pipeline on the small input first
- `claw-amplicon-qc`: preprocess downsampled amplicon reads

## Maintenance

- **Review cadence**: re-check when seqtk releases a new version.
- **Staleness signals**: a change to how `seqtk sample` parses its amount argument; a new seqtk major version.
- **Deprecation**: archive to `skills/_deprecated/` if seqtk sampling is folded into a broader FASTQ toolkit skill.

## Citations

- [seqtk](https://github.com/lh3/seqtk): the `sample` subcommand used here
