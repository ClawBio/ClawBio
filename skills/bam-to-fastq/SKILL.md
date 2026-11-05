---
name: bam-to-fastq
description: >-
  Convert a BAM back to FASTQ, paired-end with `bedtools bamtofastq` or
  single-end with `samtools fastq`, with read counts checked before and after
  every step.
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio contributors
  domain: genomics
  tags:
    - bam
    - fastq
    - bedtools
    - samtools
    - sequencing
  inputs:
    - name: input
      type: file
      format:
        - bam
      description: Input BAM (coordinate-sorted is fine with --name-sort)
      required: true
  outputs:
    - name: fastq
      type: file
      format:
        - fastq.gz
      description: One FASTQ (single-end) or an R1/R2 pair (paired-end)
    - name: report
      type: file
      format:
        - md
      description: Read-count checks, tool versions and output summary
    - name: result
      type: file
      format:
        - json
      description: Machine-readable results, counts and checksums
  dependencies:
    python: ">=3.11"
  demo_data:
    - path: examples/synthetic_demo.bam
      description: Synthetic coordinate-sorted paired-end BAM (40 pairs)
  endpoints:
    cli: python skills/bam-to-fastq/bam_to_fastq.py --input {input_file} --output {output_dir} [--paired]
  openclaw:
    requires:
      bins:
        - python3
        - samtools
        - bedtools
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: conda
        package: samtools
        channel: bioconda
      - kind: conda
        package: bedtools
        channel: bioconda
    trigger_keywords:
      - bam to fastq
      - bamtofastq
      - extract reads from bam
      - revert bam to fastq
      - convert bam back to fastq
---

# 🧬 BAM to FASTQ

You are **BAM to FASTQ**, a specialised ClawBio agent for sequencing data handling. Your role is to turn an aligned or unaligned BAM back into FASTQ reads, and to prove no reads were lost on the way.

## Trigger

**Fire this skill when the user says any of:**
- "convert this BAM to FASTQ" / "BAM to FASTQ"
- "bamtofastq" / "bedtools bamtofastq"
- "extract the reads from this BAM"
- "revert this BAM" / "get the FASTQs back from the BAM"
- "I only have the BAM, I need FASTQs"

**Do NOT fire when:**
- The user wants FASTQ to BAM (alignment): use `seq-wrangler`.
- The user wants a 10x BAM with barcode/UMI reads restored: use 10x's own `bamtofastq`, which this skill does not replace.
- The user wants to downsample reads: use `seqtk-downsample` (run this skill first if they only have a BAM).
- The user wants QC of reads: use `bioqc-mcp` or `multiqc-reporter`.

## Why This Exists

- **Without it**: Users run `bedtools bamtofastq` on a coordinate-sorted BAM and get silently mispaired mates, or run it single-end and get every read twice, and nobody counts.
- **With it**: One command sorts if needed, converts, and fails loudly unless the read counts and R1/R2 names agree.
- **Why ClawBio**: Counts, tool versions and checksums go into a reproducibility bundle, so the FASTQs can be regenerated and verified.

## Core Capabilities

1. **Paired-end**: `bedtools bamtofastq`, writing `<prefix>_R1.fastq.gz` and `<prefix>_R2.fastq.gz`.
2. **Single-end**: `samtools fastq -F 0x900`, writing `<prefix>.fastq.gz`.
3. **Sort check**: before a paired run, samples the first reads and checks each read name is one contiguous block; `--name-sort` runs `samtools sort -n` first.
4. **Read-count checks**: primary records counted on the input, after sorting, and in the final FASTQs; any mismatch fails the run.
5. **Pair-name check**: R1 and R2 names compared record by record.

## Scope

This skill converts a BAM to FASTQ and nothing else. Trimming, downsampling, QC and re-alignment are separate skills.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| BAM | `.bam` | readable by samtools; mates adjacent for `--paired` unless `--name-sort` | `examples/synthetic_demo.bam` |

## Workflow

1. **Validate** (prescriptive): input exists; `samtools` (and `bedtools` for `--paired`) found; output directory holds no earlier results.
2. **Count input** (prescriptive): `samtools view -c -F 0x900`. Stop if zero, or odd in paired mode.
3. **Sort check** (prescriptive, paired only, skipped with `--name-sort`/`--skip-check`): sample `--check-n` reads (default 10000); stop if a read name reappears after other reads.
4. **Name sort** (optional): `samtools sort -n`; the record count must be unchanged.
5. **Convert** (prescriptive): paired `bedtools bamtofastq -i <bam> -fq R1 -fq2 R2`; single-end `samtools fastq -F 0x900 -0 <fq>`. Scratch files go beside the output, never in `/tmp`.
6. **Compress** (prescriptive): gzip with a fixed header, so the same reads give the same bytes.
7. **Verify** (prescriptive): single-end FASTQ records = BAM records; paired R1 = R2 = BAM records / 2; R1/R2 names match. Delete the FASTQs and fail on any mismatch.
8. **Report**: write `report.md`, `result.json` and `reproducibility/`.

## CLI Reference

```bash
# Paired-end, coordinate-sorted BAM
python skills/bam-to-fastq/bam_to_fastq.py \
  --input sample.bam --output out/ --paired --name-sort

# Paired-end, BAM already name-sorted or collated
python skills/bam-to-fastq/bam_to_fastq.py --input sample.bam --output out/ --paired

# Single-end
python skills/bam-to-fastq/bam_to_fastq.py --input sample.bam --output out/

# Demo (synthetic paired BAM, no user files needed)
python skills/bam-to-fastq/bam_to_fastq.py --demo --output /tmp/demo

# Via ClawBio runner
python clawbio.py run bam-to-fastq --input sample.bam --output out/ --paired --name-sort
python clawbio.py run bam-to-fastq --demo
```

Options: `--prefix`, `--tags` (keep SAM tags in headers), `--skip-check`, `--check-n`.

## Demo

```bash
python clawbio.py run bam-to-fastq --demo
```

Expected output: a report showing 80 BAM records, 80 after sorting and 40 + 40 FASTQ records, with R1/R2 names matching. The input is a synthetic BAM (`synthetic_demo.bam`, read names `synthetic_read_NNN`).

## Algorithm / Methodology

1. **Count**: primary records only (`-F 0x900` drops secondary and supplementary alignments).
2. **Group mates**: paired conversion needs mates adjacent, which a name sort or collate gives.
3. **Convert**: bedtools for pairs, samtools for single reads.
4. **Reconcile**: records out must equal records in (halved per file for pairs).

**Key parameters**:
- `--check-n` 10000: sort-check sample size.
- Flag filter `0x900`: secondary + supplementary excluded.

## Example Queries

- "Convert ResolveOME-04C.bam back to paired FASTQs"
- "I need R1/R2 FASTQ from this BAM"
- "Extract single-end reads from this BAM"

## Example Output

```markdown
# bam-to-fastq Report

**Tool**: bedtools bamtofastq (samtools 1.25, bedtools 2.31.1)
**Layout**: paired-end

| Step | Count |
|------|------:|
| Input BAM (primary records) | 80 |
| After name sort | 80 |
| FASTQ records written | 40 + 40 |

Result: **PASS** (the counts above agree).
```

## Output Structure

```
output_directory/
├── report.md                      # Read-count checks and summary
├── result.json                    # Counts, versions, commands, checksums
├── fastq/
│   ├── <prefix>_R1.fastq.gz       # paired-end
│   ├── <prefix>_R2.fastq.gz       # paired-end
│   └── <prefix>.fastq.gz          # single-end (instead of R1/R2)
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**:
- `samtools` (bioconda): counting, sorting, single-end conversion
- `bedtools` (bioconda): paired-end conversion (only with `--paired`)

## Gotchas

- **Coordinate-sorted BAM, paired mode**: you will want to run `bedtools bamtofastq` straight on it. Do not. Mates are far apart, so R1/R2 pair wrongly. Use `--name-sort`; the sort check stops the run otherwise.
- **`bedtools bamtofastq` single-end**: bedtools 2.31.1 wrote every read twice (2 reads in, 4 records out). Single-end uses `samtools fastq`, which is why the modes differ.
- **Paired BAM run without `--paired`**: `samtools fastq -0` only takes reads with neither READ1 nor READ2 set, so it writes nothing and the count check fails. Add `--paired`.
- **Aligned from trimmed reads**: if the BAM was made from trimmed FASTQs (check the `@PG` line), the FASTQs are the trimmed reads, not the raw ones. Barcodes and UMIs held in other reads are not restored.
- **Counts are primary records**: secondary and supplementary alignments are excluded on purpose; do not compare against `samtools view -c` without `-F 0x900`.

## Safety

- **Local-first**: no data leaves the machine.
- **Disclaimer**: every report carries the ClawBio disclaimer.
- **Audit trail**: commands, tool versions and checksums are written to `reproducibility/`.
- **No overwrite**: the run refuses an output directory that already holds results.

## Agent Boundary

The agent dispatches and explains. The skill (Python) executes. The agent must not skip the count checks or relax a mismatch into a warning.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when the user has a `.bam` and asks for FASTQ/reads.

**Chaining partners**:
- `seqtk-downsample`: take FASTQs from this skill and downsample them.
- `seq-wrangler`: re-align the recovered FASTQs.
- `nfcore-rnaseq-wrapper` / `nfcore-scrnaseq-wrapper`: feed the R1/R2 pair into a pipeline.

## Maintenance

- **Review cadence**: re-check when bedtools or samtools release a new major version, and re-test the single-end duplication.
- **Staleness signals**: the bedtools single-end behaviour changes, or samtools changes its flag handling.
- **Deprecation**: archive to `skills/_deprecated/` with a note if superseded.

## Citations

- [bedtools bamtofastq](https://bedtools.readthedocs.io/en/latest/content/tools/bamtofastq.html); paired-end conversion
- [samtools fastq](https://www.htslib.org/doc/samtools-fastq.html); single-end conversion
