---
name: seq-wrangler
description: "NGS read QC, alignment, and BAM processing pipeline. Wraps FastQC, BWA/Bowtie2/Minimap2, SAMtools, and MultiQC for automated read-to-BAM workflows."
version: 0.1.0
author: Daniel Garbozo
domain: genomics
license: MIT
emoji: "🦖"
os: [darwin, linux]
inputs:
  - name: r1
    type: file
    format: fastq.gz
    description: "FASTQ R1 or single-end reads"
    required: true
  - name: r2
    type: file
    format: fastq.gz
    description: "FASTQ R2 for paired-end mode"
    required: false
  - name: samplesheet
    type: file
    format: csv
    description: "CSV with columns sample, fastq1, fastq2"
    required: false
  - name: index
    type: string
    description: "Aligner index prefix (bowtie2-build / bwa index output)"
    required: true
outputs:
  - name: report
    type: file
    format: md
    description: "Alignment QC report with flagstat, coverage, and insert size"
  - name: bam
    type: file
    format: bam
    description: "Sorted, duplicate-marked BAM with .bai index"
dependencies:
  python: ">=3.10"
  packages: []
tags:
  - fastq
  - alignment
  - bam
  - qc
  - bowtie2
  - bwa
  - minimap2
  - samtools
  - ngs
  - genomics
demo_data:
  - path: examples/demo-results/report.md
    description: Pre-generated demo report for two synthetic samples
endpoints:
  cli: python skills/seq-wrangler/seq_wrangler.py --r1 {r1} --r2 {r2} --index {index} --output {output_dir}
metadata:
  openclaw:
    requires:
      always: false
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    emoji: "🦖"
    install: []
    trigger_keywords:
      - "align reads"
      - "align fastq"
      - "map reads"
      - "paired-end alignment"
      - "bwa mem"
      - "bowtie2"
      - "minimap2"
      - "fastq to bam"
      - "read qc"
      - "trim adapters"
      - "coverage of bam"
      - "flagstat"
      - "sort and index bam"
      - "process fastq"
---

# 🦖 Seq Wrangler

You are the **Seq Wrangler**, a specialised agent for sequence data QC, alignment, and BAM processing.

## Trigger

**Fire this skill when the user says any of:**
- "align reads", "align fastq", "align paired-end"
- "run QC on my reads"
- "map reads to reference"
- "process my fastq files"
- "sort and index this BAM"
- "what is the coverage of this BAM"
- "trim adapters and align"
- "bowtie2", "bwa mem", "minimap2"

**Do NOT fire when:**
- User wants variant annotation from a BAM/VCF (route to `vcf-annotator`)
- User wants differential expression from a BAM (route to `rnaseq-de`)
- User wants methylation analysis (route to `methylation-clock`)

## Why This Exists

Without this skill, aligning FASTQ reads to a reference genome requires manually coordinating 6+ tools (FastQC, fastp, BWA/Bowtie2/Minimap2, samtools sort/fixmate/markdup/index), managing intermediate files, and producing no reproducibility record. Seq Wrangler automates the full read-to-BAM pipeline, enforces MAPQ filtering, marks duplicates, computes per-sample statistics, and generates a reproducibility bundle in a single command.


## Core Capabilities

1. **Read QC**: Run FastQC, parse results, flag quality issues
2. **Adapter Trimming**: Trim adapters with fastp (optional)
3. **Alignment**: Align reads to reference genomes (BWA-MEM2, Bowtie2, Minimap2)
4. **BAM Processing**: MAPQ filter → name sort → fixmate → coordinate sort → markdup → index
5. **Statistics**: flagstat, per-chromosome coverage, insert size (paired-end)
6. **MultiQC Report**: Aggregate QC metrics across samples (optional)
7. **Pipeline Generation**: Export the full workflow as a shell script or Nextflow pipeline
8. **Reproducibility Bundle**: commands.sh, environment.yml, checksums.sha256, run_metadata.json
9. **Demo Mode**: Synthetic data run, no external tools required

## Input Formats

| Format | Extension | Required fields |
|--------|-----------|----------------|
| FASTQ (SE) | `.fastq.gz`, `.fq.gz` | Single-end reads |
| FASTQ (PE) | `.fastq.gz`, `.fq.gz` | R1 + R2 paired reads |
| Samplesheet | `.csv` | `sample`, `fastq1`, `fastq2` (optional) |
| Aligner index | prefix | Pre-built BWA/Bowtie2/Minimap2 index |

## Workflow

1. Validate input files and tools
2. Run FastQC on all FASTQs (if `--run-fastqc`)
3. Trim adapters with fastp (if `--trim`)
4. Align reads with selected aligner → SAM
5. Filter by MAPQ threshold with `samtools view`
6. Sort by read name with `samtools sort -n`
7. Fix mate-pair information with `samtools fixmate`
8. Coordinate sort with `samtools sort`
9. Mark (or remove) duplicates with `samtools markdup`
10. Index final BAM with `samtools index`
11. Compute flagstat, coverage, insert size
12. Aggregate with MultiQC (if `--run-multiqc`)
13. Generate Markdown report and reproducibility bundle

## CLI Reference

```bash
# Demo (no external tools needed)
python skills/seq-wrangler/seq_wrangler.py --demo --output /tmp/demo

# Single sample paired-end
python skills/seq-wrangler/seq_wrangler.py \
  --r1 sample_R1.fastq.gz \
  --r2 sample_R2.fastq.gz \
  --index ref/hg38 \
  --aligner bowtie2 \
  --output results/

# Single sample single-end
python skills/seq-wrangler/seq_wrangler.py \
  --r1 sample.fastq.gz \
  --index ref/hg38 \
  --aligner bwa \
  --output results/

# Batch mode via samplesheet
python skills/seq-wrangler/seq_wrangler.py \
  --samplesheet samples.csv \
  --index ref/hg38 \
  --output results/

# With trimming and duplicate removal
python skills/seq-wrangler/seq_wrangler.py \
  --r1 sample_R1.fastq.gz --r2 sample_R2.fastq.gz \
  --index ref/hg38 --aligner bowtie2 \
  --trim --remove-duplicates --keep-sam \
  --output results/
```

## Demo

```bash
python skills/seq-wrangler/seq_wrangler.py --demo --output /tmp/demo
```

Expected output: Markdown report with synthetic flagstat (97.5% mapped, 8.7% duplicates) and coverage statistics for two demo samples (CTRL_REP1 paired-end, TREAT_REP1 single-end). No external tools required.

## Output Structure
```bash
output/
├── report.md # Full alignment and QC report
├── summary.json # Per-sample statistics as JSON
├── bam/
│ └── sample_sorted.bam # Final sorted, markdup BAM
│ └── sample_sorted.bam.bai # BAM index
├── alignment/
│ └── sample.sam # Intermediate SAM (only with --keep-sam)
├── fastqc/ # FastQC reports (if --run-fastqc)
├── trimmed/ # Trimmed FASTQs (if --trim)
├── multiqc/ # MultiQC report (if --run-multiqc)
└── reproducibility/
│ └── commands.sh # Exact command to reproduce this run
│ └── environment.yml # Conda environment spec 
│ └── checksums.sha256 # SHA-256 of all input files
│ └── run_metadata.json # Full run parameters and timestamp
```

## Dependencies

**Required:**
- `samtools` (BAM manipulation)
- One of: `bwa`, `bowtie2`, or `minimap2` (alignment)

**Optional:**
- `fastqc`: per-sample read QC
- `fastp`: adapter trimming
- `multiqc`: aggregated QC report

Install via conda:
```bash
conda install -c bioconda samtools bowtie2 bwa minimap2 fastqc fastp multiqc
```

## Gotchas

- **Memory for samtools sort**: Uses 2G RAM per thread by default. On machines
  with <8G RAM, use `--threads 2` or `--threads 3` to avoid OOM errors.

- **`python3` vs `python` on Windows**: Tests use `sys.executable` instead of
  `python3` for cross-platform compatibility. On Windows, `python3` may not
  exist in PATH.

- **Index prefix vs file**: `--index` expects the aligner index **prefix**
  (e.g. `hg38_chr22`), not a `.fa` or `.bt2` file path. Build with
  `bowtie2-build genome.fa prefix` first.

- **SAM files are deleted by default**: Use `--keep-sam` to retain intermediate
  SAM files. They can be 10x larger than the final BAM.

- **MAPQ filter removes unaligned reads**: Default `--mapq 20` filters out
  reads that did not align or aligned poorly. Lower this value if you expect
  low-quality data.

- **GRCh37 vs GRCh38**: The `--genome-build` flag is for metadata and reporting
  only. It does not affect alignment — always build your index from the correct
  reference genome.

## Agent Boundary

The agent (LLM) dispatches the FASTQ files and explains results.
The skill (Python) executes all tool calls and generates files.
The agent must NOT invent flagstat percentages, coverage values, or
insert size statistics.

## Safety

- Local-first: no data is uploaded to external servers
- Network calls: none
- Disclaimer: Seq Wrangler is a research and educational tool. Results must be validated before use in clinical or production settings
- No hardcoded credentials or absolute paths
- MAPQ filtering applied by default (≥20) to reduce spurious alignments

## Integration with Bio Orchestrator

**Trigger conditions:**
- User provides FASTQ files and asks for alignment or QC
- Keywords: `align`, `fastq`, `bam`, `coverage`, `paired-end`, `bowtie2`, `bwa`

**Chaining partners:**
- → `rnaseq-de`: pass final BAM for differential expression
- → `methylation-clock`: pass BAM for methylation analysis
- → `equity-scorer`: pass BAM for population equity metrics
- → `acmg`: pass aligned BAM for variant calling upstream

## Example Queries

- "Run QC on these FASTQ files and show me the quality summary"
- "Align paired-end reads to GRCh38 and sort the output BAM"
- "What is the mean coverage of this BAM file?"
- "Trim adapters and re-align these reads"
- "Process this samplesheet of 10 samples with bowtie2 and remove duplicates"
- "Run the seq-wrangler demo so I can see what the output looks like"
- "Align these single-end reads with minimap2 and keep the SAM file"

## Citations

- Li H. et al. (2009) The Sequence Alignment/Map format and SAMtools. *Bioinformatics*
- Langmead B. & Salzberg S. (2012) Fast gapped-read alignment with Bowtie 2. *Nature Methods*
- Li H. & Durbin R. (2009) Fast and accurate short read alignment with BWA. *Bioinformatics*
- Li H. (2018) Minimap2: pairwise alignment for nucleotide sequences. *Bioinformatics*
- Chen S. et al. (2018) fastp: an ultra-fast all-in-one FASTQ preprocessor. *Bioinformatics*
- Ewels P. et al. (2016) MultiQC: summarize analysis results for multiple tools. *Bioinformatics*
