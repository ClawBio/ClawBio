---
name: claw-amplicon-qc
description: 16S rRNA amplicon preprocessing — read stats, N removal, and primer trimming
license: MIT
metadata:
  version: 0.1.0
  author: Zabiulla
  tags:
  - amplicon
  - 16S
  - 18S
  - DADA2
  - cutadapt
  - quality-control
  - environmental
  - microbiome
  inputs:
  - name: raw_folder
    type: directory
    description: Directory containing paired-end FASTQ files (auto-discovers R1/R2 pairs matching *_R1*.fastq.gz and *_R2*.fastq.gz)
  - name: output_folder
    type: directory
    description: Directory where all outputs will be written (created if missing)
  - name: fwd_primer
    type: string
    description: Forward primer sequence in IUPAC nucleotide codes (e.g. CCTACGGGNGGCWGCAG for 341F)
  - name: rev_primer
    type: string
    description: Reverse primer sequence in IUPAC nucleotide codes (e.g. GACTACHVGGGTWTCTAAT for 806R)
  - name: min_length
    type: integer
    description: Minimum read length in bp after primer trimming (reads shorter than this are discarded)
  outputs:
  - name: filtN
    type: directory
    format:
    - fastq.gz
    description: FASTQ files after removal of reads containing N bases (required before Cutadapt)
  - name: cutadapt_trimmed
    type: directory
    format:
    - fastq.gz
    description: Primer-trimmed FASTQ files, ready for downstream DADA2 quality filtering and denoising
  - name: raw_stats
    type: file
    format: tsv
    description: seqkit statistics on raw input FASTQ (baseline snapshot before any transformation)
  - name: filtN_stats
    type: file
    format: tsv
    description: seqkit statistics on N-filtered FASTQ (measures the effect of N removal)
  - name: trimmed_stats
    type: file
    format: tsv
    description: seqkit statistics on primer-trimmed FASTQ (measures the effect of Cutadapt)
  - name: cutadapt_log
    type: file
    format: txt
    description: Full Cutadapt per-sample log capturing primer detection and trimming details
  - name: report
    type: file
    format: md
    description: Human-readable summary of the run — samples processed, retention rates, flagged samples
  - name: qc_summary
    type: file
    format: json
    description: Machine-readable structured summary with per-sample stats, retention rates, and automatic flags for anomalies
  openclaw:
    category: bioinformatics
    emoji: 🧬
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    system_dependencies:
    - R (>=4.4)
    - bioconductor-dada2
    - bioconductor-shortread
    - bioconductor-biostrings
    - r-yaml
    - cutadapt
    - seqkit
    requires:
      bins:
      - Rscript
      - cutadapt
      - seqkit
    always: false
---

# 16S Amplicon Preprocessing (claw-amplicon-qc)

Preprocessing pipeline for 16S/18S rRNA amplicon sequencing data — from raw paired-end FASTQ files through primer-trimmed reads ready for DADA2 quality filtering and denoising.

## What it does

1. Takes a folder of paired-end FASTQ files and auto-discovers R1/R2 sample pairs
2. Runs **seqkit stats** on the raw files to record baseline read counts, lengths, GC content, and quality
3. Removes any read containing an N base using **DADA2 filterAndTrim** with `maxN=0` (Cutadapt cannot detect primers in reads with ambiguous bases)
4. Trims primers using **Cutadapt** with:
   - Forward primer (5') removed from R1
   - Reverse-complement of reverse primer (3') removed from R1
   - Reverse primer (5') removed from R2
   - Reverse-complement of forward primer (3') removed from R2
   - `--nextseq-trim=20` for Illumina two-colour chemistry poly-G trimming
   - `-n 2` to catch reads with primer dimers or read-through
5. Runs **seqkit stats** again on the trimmed files to record post-trimming statistics
6. Calculates retention rates at each stage and flags samples with anomalous drops
7. Writes a human-readable `report.md` and a machine-readable `qc_summary.json`
8. Stops before quality filtering and DADA2 denoising — those decisions require researcher judgment and belong to a separate skill

## Why this exists

If you ask a general AI to "run 16S QC," it will typically:

- Combine N removal and quality filtering into a single step, hiding which reads were lost to which criterion
- Use default DADA2 quality parameters (`maxEE`, `truncQ`) before the researcher has reviewed quality profiles
- Skip the primer sanity check and produce silent failures downstream if primers were missed
- Not distinguish between the forward primer and its reverse-complement, which matters for correct 3' trimming
- Not track per-sample read survival across stages, making it hard to spot problematic samples
- Not produce a machine-readable summary that a downstream tool or AI wrapper can consume

## Scientific decisions encoded

Several methodological choices are baked into this skill. Understanding them helps a researcher (or an AI wrapper) know when the defaults are appropriate and when they are not.

- **N removal happens before primer trimming, not with quality filtering.** Cutadapt cannot detect primer sequences in reads containing ambiguous bases (N). This is not an optimisation — it is a hard prerequisite. Combining N removal with quality filtering (as some tutorials do) makes it impossible to know which reads were lost to which criterion.

- **Quality filtering (`maxEE`, `truncQ`, `truncLen`) is deliberately deferred to a later skill.** These parameters depend on quality profiles the researcher has not yet reviewed. Applying default values here would either be too lenient (wasting downstream compute on unusable reads) or too aggressive (discarding recoverable data). Deferring keeps the human in the loop for a decision that requires their judgment.

- **All four primer positions are trimmed explicitly.** In paired-end amplicon libraries, each read may contain primer sequence at both ends — the primer that primed the read at the 5' end, and the reverse-complement of the opposite primer at the 3' end if the read is long enough to read through the amplicon. Trimming only the 5' primer leaves artificial sequence at the 3' end that would corrupt DADA2 error learning and paired-end merging.

- **Baseline read statistics are captured before any transformation.** This establishes the reference point against which every downstream loss is measured, so retention rates are meaningful and per-sample anomalies are visible from the earliest stage.

- **`min_length` is a required user input, not a hard-coded default.** The appropriate minimum length after primer trimming depends on the amplicon region (V3–V4 vs V4 differ by ~200bp) and the expected paired-end overlap. A generic default would produce silent failures for non-V4 primer pairs.

- **Anomaly flags are informational, not blocking.** When a sample loses more reads than expected, the skill still completes normally and records the flag in `qc_summary.json` and `report.md`. The decision to exclude a flagged sample belongs to the downstream skill or the researcher, not to this preprocessing step.

## Validated On

Environmental 16S samples (aquatic plant roots, fronds, and surrounding water) using 341F/806R primers on Illumina MiSeq 2×300 paired-end sequencing. Sample sizes tested: 20 to 100+ samples per run.

## Pipeline Architecture

```
raw FASTQ files (paired-end)
        │
        ▼
[Stage 1]  seqkit stats  ──►  01_raw_stats.txt
        │
        ▼
[Stage 2]  DADA2 filterAndTrim (maxN=0)
        │                ──►  filtN/ (N-filtered FASTQ)
        │                ──►  02_filtN_stats.txt
        ▼
[Stage 3]  Cutadapt primer trimming
        │                ──►  cutadapt_trimmed/ (trimmed FASTQ)
        │                ──►  03_cutadapt_log.txt
        │                ──►  03_trimmed_stats.txt
        ▼
        Retention analysis + anomaly flagging
        │                ──►  report.md (human-readable)
        │                ──►  qc_summary.json (machine-readable)
        ▼
        Ready for DADA2 quality profiling and denoising
        (separate skill — claw-amplicon-dada2, forthcoming)
```

## Automatic Flags

The skill sets flags in `qc_summary.json` for any sample matching these conditions:

| Flag | Trigger | Action for researcher |
|------|---------|----------------------|
| `extreme_drop` | Sample loses >50% of reads at any single stage | Investigate — likely primer mismatch, low library quality, or contamination |
| `low_read_count` | Sample retains <1000 reads after primer trimming | Consider excluding — insufficient depth for reliable ASV inference |
| `high_overall_loss` | Sample retains <70% of raw reads overall | Investigate — cumulative quality issues |

Flags are informational, not blocking. The skill completes normally and reports flagged samples in `report.md` and `qc_summary.json`. The downstream DADA2 skill (or the researcher) decides how to handle them.

## Usage

```bash
Rscript amplicon_qc.R \
    --raw /path/to/raw_fastq_folder \
    --output /path/to/output_folder \
    --fwd-primer CCTACGGGNGGCWGCAG \
    --rev-primer GACTACHVGGGTWTCTAAT \
    --min-length 200
```

### Arguments

| Argument | Required | Description |
|----------|----------|-------------|
| `--raw` | yes | Folder containing paired-end FASTQ files |
| `--output` | yes | Folder for outputs (created if missing) |
| `--fwd-primer` | yes | Forward primer sequence (IUPAC codes) |
| `--rev-primer` | yes | Reverse primer sequence (IUPAC codes) |
| `--min-length` | yes | Minimum read length in bp after trimming |

### Common primer pairs (for reference)

The skill does not hard-code any primer pairs — the researcher must provide sequences. Common choices include:

| Region | Forward | Sequence | Reverse | Sequence |
|--------|---------|----------|---------|----------|
| V3–V4 | 341F | `CCTACGGGNGGCWGCAG` | 806R | `GACTACHVGGGTWTCTAAT` |
| V4 | 515F | `GTGYCAGCMGCCGCGGTAA` | 806R | `GGACTACNVGGGTWTCTAAT` |
| V1–V2 | 27F | `AGAGTTTGATCMTGGCTCAG` | 338R | `TGCTGCCTCCCGTAGGAGT` |

## Example Output

```
─────────────────────────────────────────────────────────
  claw-amplicon-qc — Stage 1/3: raw read statistics
─────────────────────────────────────────────────────────
  20 sample pairs discovered
  seqkit stats written to: 01_raw_stats.txt

─────────────────────────────────────────────────────────
  claw-amplicon-qc — Stage 2/3: removing reads with N bases
─────────────────────────────────────────────────────────
  Multithreading: enabled (16 cores)
  Reads in:  2,458,932
  Reads out: 2,441,207
  Retention: 99.28%

─────────────────────────────────────────────────────────
  claw-amplicon-qc — Stage 3/3: primer trimming (Cutadapt)
─────────────────────────────────────────────────────────
  Forward primer:  CCTACGGGNGGCWGCAG
  Reverse primer:  GACTACHVGGGTWTCTAAT
  Trimming 20 sample pairs...
  Retention: 87.3%
  Flagged samples: Sample_07 (extreme_drop, 52% loss at Cutadapt)

═════════════════════════════════════════════════════════
  claw-amplicon-qc — Complete
═════════════════════════════════════════════════════════
  Samples processed:  20
  Samples flagged:    1
  Output folder:      /path/to/output_folder

  Next step: quality profiling and DADA2 denoising
  (separate skill — claw-amplicon-dada2, forthcoming)
```

## Testing

A basic end-to-end test lives in `tests/` and can be run with:

```bash
cd tests
bash run_test.sh
```

The test uses a tiny bundled FASTQ pair (~1000 reads per sample) with known primer sequences to verify the pipeline runs to completion and produces expected output files.

## Dependencies

Install via conda:

```bash
conda create -n amplicon-qc -c conda-forge -c bioconda \
    r-base=4.4 \
    bioconductor-dada2 \
    bioconductor-shortread \
    bioconductor-biostrings \
    r-yaml \
    cutadapt \
    seqkit
```

## What Comes Next

This skill deliberately stops before quality filtering and DADA2 denoising. Those steps depend on visual review of quality profiles and researcher judgment about truncation lengths, chimera thresholds, and pooling strategy — decisions that belong in a separate skill with an approval-gate pattern.

The forthcoming `claw-amplicon-dada2` skill will handle:

- Quality profile generation (`plotQualityProfile`)
- Human-in-the-loop recommendation of `truncLen` values
- Quality filtering with user-approved parameters
- Error learning and denoising (`learnErrors`, `dada`)
- Paired-end merging (`mergePairs`)
- ASV table construction and chimera removal

## Citations

If you use this skill in a publication, please cite the underlying tools:

- **DADA2** — Callahan BJ, McMurdie PJ, Rosen MJ, Han AW, Johnson AJ, Holmes SP (2016). *DADA2: High-resolution sample inference from Illumina amplicon data.* Nature Methods 13:581–583. https://doi.org/10.1038/nmeth.3869
- **Cutadapt** — Martin M (2011). *Cutadapt removes adapter sequences from high-throughput sequencing reads.* EMBnet.journal 17(1):10–12. https://doi.org/10.14806/ej.17.1.200
- **seqkit** — Shen W, Le S, Li Y, Hu F (2016). *SeqKit: A cross-platform and ultrafast toolkit for FASTA/Q file manipulation.* PLoS ONE 11(10):e0163962. https://doi.org/10.1371/journal.pone.0163962

## License

MIT — see LICENSE file at the repository root.
