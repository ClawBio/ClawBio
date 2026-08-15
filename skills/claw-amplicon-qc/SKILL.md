---
name: claw-amplicon-qc
description: >-
  16S/18S rRNA amplicon preprocessing — from raw paired-end FASTQ files
  through N removal and primer trimming, producing outputs ready for DADA2
  quality filtering and denoising. Deliberately stops before quality-filtering
  decisions that require researcher judgment.
license: MIT
metadata:
  version: "0.2.0"
  author: Zabiulla
  domain: genomics
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
      description: >-
        Directory containing paired-end FASTQ files. Auto-discovers R1/R2 pairs
        matching *_R1[_.]*.fastq.gz / *_R2[_.]*.fastq.gz (Illumina default) or
        R1_*.fastq.gz / R2_*.fastq.gz (alternate).
      required: true
    - name: output_folder
      type: directory
      description: Directory where all outputs will be written (created if missing).
      required: true
    - name: fwd_primer
      type: string
      description: Forward primer sequence in IUPAC nucleotide codes.
      required: true
    - name: rev_primer
      type: string
      description: Reverse primer sequence in IUPAC nucleotide codes.
      required: true
    - name: min_length
      type: integer
      description: >-
        Minimum read length in bp after primer trimming. Reads shorter than
        this are discarded by Cutadapt.
      required: true
  outputs:
    - name: filtN
      type: directory
      format: [fastq.gz]
      description: FASTQ files after removal of reads containing N bases.
    - name: cutadapt_trimmed
      type: directory
      format: [fastq.gz]
      description: Primer-trimmed FASTQ files, ready for downstream DADA2 quality filtering.
    - name: raw_stats
      type: file
      format: [tsv]
      description: seqkit statistics on raw input FASTQ (baseline before any transformation).
    - name: filtN_stats
      type: file
      format: [tsv]
      description: seqkit statistics on N-filtered FASTQ.
    - name: trimmed_stats
      type: file
      format: [tsv]
      description: seqkit statistics on primer-trimmed FASTQ.
    - name: cutadapt_log
      type: file
      format: [txt]
      description: Per-sample Cutadapt log with primer detection and trimming details.
    - name: report
      type: file
      format: [md]
      description: >-
        Human-readable QC report — samples processed, retention rates, failed
        samples, and flagged anomalies.
    - name: qc_summary
      type: file
      format: [json]
      description: >-
        Machine-readable structured summary with per-sample stats, retention
        rates, and flags. Consumable by downstream skills or wrappers.
  dependencies:
    r: ">=4.4"
    packages:
      - dada2
      - ShortRead
      - Biostrings
      - optparse
      - jsonlite
  demo_data:
    - path: tests/fixtures/demo_R1.fastq.gz
      description: >-
        Real V4 subset from astrobiomike Deep Sea Rock 16S tutorial
        (~500 read pairs, 515F/806R primers attached, CC-BY 4.0).
    - path: tests/fixtures/demo_R2.fastq.gz
      description: Reverse-read partner for demo_R1.
  endpoints:
    cli: >-
      Rscript skills/claw-amplicon-qc/amplicon_qc.R
      --raw {raw_folder}
      --output {output_folder}
      --fwd-primer {fwd_primer}
      --rev-primer {rev_primer}
      --min-length {min_length}
  openclaw:
    category: bioinformatics
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    requires:
      bins:
        - Rscript
        - cutadapt
        - seqkit
    install:
      - kind: conda
        package: r-base>=4.4
      - kind: conda
        package: bioconductor-dada2
      - kind: conda
        package: bioconductor-shortread
      - kind: conda
        package: bioconductor-biostrings
      - kind: conda
        package: r-optparse
      - kind: conda
        package: r-jsonlite
      - kind: conda
        package: cutadapt
      - kind: conda
        package: seqkit
    trigger_keywords:
      - 16S amplicon QC
      - 16S primer trimming
      - amplicon preprocessing
      - amplicon QC
      - cutadapt primer trimming
      - DADA2 preprocessing
      - remove primers from FASTQ
      - 18S amplicon QC
      - microbiome QC
      - rRNA amplicon preprocessing
    always: false
---

# 🧬 16S Amplicon Preprocessing (claw-amplicon-qc)

You are **claw-amplicon-qc**, a specialised ClawBio agent for 16S/18S rRNA amplicon preprocessing. Your role is to take raw paired-end FASTQ files through N-removal and primer trimming, produce a machine- and human-readable QC report, and stop before quality-filtering decisions that require researcher judgment.

## Trigger

**Fire when the user asks to:**

- Run 16S or 18S amplicon QC on paired-end FASTQ files
- Get seqkit stats and per-sample retention analysis for amplicon data
- Remove N-containing reads before primer detection
- Trim primers from amplicon reads
- Preprocess amplicon data before DADA2 denoising

**Do NOT fire when the user asks to:**

- Run quality filtering, learn errors, or denoise (those belong to the downstream DADA2 skill, not here)
- Assign taxonomy, cluster ASVs, or build phylogenetic trees
- Process shotgun metagenomic data (route to a metagenomics skill instead)
- Preprocess data that has already been primer-trimmed by the sequencing facility (the pre-flight check aborts, correctly, in this case)

## Why This Exists

If you ask a general-purpose AI to "run 16S QC," it will typically:

- Combine N removal and quality filtering into a single step, hiding which reads were lost to which criterion
- Apply DADA2 quality parameters (`maxEE`, `truncQ`) before the researcher has reviewed quality profiles
- Skip primer sanity checks and produce silent failures downstream if primers were missed
- Fail to distinguish between the forward primer and its reverse-complement, breaking 3' trimming
- Not track per-sample read survival across stages, hiding problematic samples
- Not produce a machine-readable summary a downstream tool or AI wrapper can consume

This skill does the boring, mechanical preprocessing correctly, keeps every methodological choice explicit and configurable, produces both a human report and a machine-readable JSON summary, and hands off to a separate skill for the parameter-sensitive steps that require researcher judgment.

## Core Capabilities

1. **N-first preprocessing**: Remove reads with ambiguous (N) bases before primer detection. Cutadapt cannot detect primers in reads with Ns, so this is a prerequisite, not an optimisation.
2. **Anchored, four-position primer trimming**: Trim primers from all four positions — 5' and 3' on both R1 and R2 — with anchored 5' primers via Cutadapt. Prevents mid-read primer motifs from truncating real biology.
3. **Per-sample retention analysis with R1/R2 accounting**: Track read survival across stages independently for R1 and R2. Flag samples with anomalous drops, low read counts, or high overall loss using either-side threshold logic.
4. **Structured + human-readable outputs**: Emit `qc_summary.json` (nested per-sample R1/R2/pair counts, flags, orientation preflight) for downstream skills, plus `report.md` with the ClawBio disclaimer for humans.

## Scope

**One skill, one task.** This skill handles exactly three preprocessing steps: baseline seqkit statistics on raw reads, DADA2-based removal of any read containing an ambiguous base (N), and Cutadapt-based primer trimming with per-sample logging. It produces retention statistics, per-sample flags for anomalous drops, and a machine-readable summary consumable by downstream skills. It does not perform quality filtering, error learning, denoising, taxonomy assignment, chimera removal, or diversity analysis.

## Input Formats

| Format | Extension | Required | Example |
|--------|-----------|----------|---------|
| Paired-end FASTQ, gzipped | `.fastq.gz` | Matched R1/R2 pairs | `sample_R1_001.fastq.gz` + `sample_R2_001.fastq.gz` |
| Paired-end FASTQ, gzipped short-form | `.fq.gz` | Matched R1/R2 pairs | `sample_R1.fq.gz` + `sample_R2.fq.gz` |
| Paired-end FASTQ, uncompressed | `.fastq` / `.fq` | Matched R1/R2 pairs | Same as above without `.gz` |

**Naming patterns accepted:**

- **Pattern A (Illumina default):** `SAMPLE_R1_...fastq.gz` + `SAMPLE_R2_...fastq.gz`
- **Pattern B (prefix):** `R1_SAMPLE.fastq.gz` + `R2_SAMPLE.fastq.gz`

The skill auto-detects which pattern is in use. Files must arrive in matched R1/R2 pairs — mispaired or orphan inputs are detected at file discovery and rejected with a clear error. Sample names are extracted from filenames deterministically and used as report row IDs.

## Workflow

When the user asks for 16S/18S amplicon preprocessing:

1. **Discover**: Enumerate paired-end FASTQ files in the input folder and pair by sample name (name-based, not sort-position, so mispaired inputs are caught).
2. **Pre-flight**: Sample the first ~1000 R1 reads per sample, count how many start with the forward vs reverse primer. Abort if primers not detected, all-reversed, or mixed orientation (unless `--allow-mixed-orientation` is set).
3. **Stage 1**: Baseline seqkit statistics on raw reads.
4. **Stage 2**: DADA2 `filterAndTrim` with `maxN=0`, `truncQ=0`, `minLen=0`, `rm.phix=TRUE` — removes any read containing an N, plus PhiX contamination.
5. **Stage 3**: seqkit statistics on N-filtered reads.
6. **Stage 4**: Cutadapt per-sample primer trimming with all four positions (5' and 3' on both R1 and R2), anchored 5' primers, `--discard-untrimmed`, `--pair-filter=any`.
7. **Stage 5**: seqkit statistics on primer-trimmed reads.
8. **Report**: Retention analysis, either-side flag evaluation, human + machine outputs. Aborts if overall retention falls below 1% (100%-data-loss backstop for configuration errors).

**Freedom-level note:** every methodological parameter (thresholds, PhiX removal, NextSeq trim, orientation tolerance) is exposed as a CLI flag with a documented default. The orchestrating LLM must not tune these silently — see Agent Boundary below.

## CLI Reference

```bash
# Standard invocation
Rscript skills/claw-amplicon-qc/amplicon_qc.R \
    --raw /path/to/raw_fastq_folder \
    --output /path/to/output_folder \
    --fwd-primer GTGYCAGCMGCCGCGGTAA \
    --rev-primer GGACTACHVGGGTWTCTAAT \
    --min-length 200

# Demo mode (bundled fixture, no user files needed)
Rscript skills/claw-amplicon-qc/amplicon_qc.R --demo --output /tmp/demo
```

**Required arguments:**

| Argument | Description |
|----------|-------------|
| `--raw` | Folder containing paired-end FASTQ files. |
| `--output` | Folder for outputs (created if missing). |
| `--fwd-primer` | Forward primer sequence (IUPAC codes). |
| `--rev-primer` | Reverse primer sequence (IUPAC codes). |
| `--min-length` | Minimum read length in bp after primer trimming. |

**Optional tuning (defaults preserve documented behaviour):**

| Argument | Default | Description |
|----------|---------|-------------|
| `--max-n` | `0` | Max N bases allowed at Stage 2. 0 = discard any read with N. |
| `--nextseq-trim` | `0` | Cutadapt `--nextseq-trim=` value. 0 = flag omitted. Set to a positive integer (e.g. 20) only for NextSeq/NovaSeq data. |
| `--no-phix-removal` | (off) | Flag. Skip PhiX removal in Stage 2. |
| `--allow-mixed-orientation` | (off) | Flag. Continue when pre-flight detects mixed-orientation reads. Only use if you accept ~50% data loss. |
| `--extreme-drop-threshold` | `50` | Percent retention below which `extreme_drop` fires. |
| `--low-count-threshold` | `1000` | Read count below which `low_read_count` fires. Lower for low-biomass samples. |
| `--overall-retention-threshold` | `70` | Overall retention percent below which `high_overall_loss` fires. |
| `--demo` | (off) | Run on bundled fixture with hardcoded 515F/806R primers and `--min-length 200`. Requires `--output` only. |

**Common primer pairs (for reference):**

| Region | Forward | Sequence | Reverse | Sequence |
|--------|---------|----------|---------|----------|
| V3–V4 | 341F | `CCTACGGGNGGCWGCAG` | 806R (Caporaso) | `GGACTACHVGGGTWTCTAAT` |
| V4 | 515F (Parada) | `GTGYCAGCMGCCGCGGTAA` | 806R (Apprill/EMP) | `GGACTACNVGGGTWTCTAAT` |
| V1–V2 | 27F | `AGAGTTTGATCMTGGCTCAG` | 338R | `TGCTGCCTCCCGTAGGAGT` |

## Demo

The bundled fixture at `tests/fixtures/demo_R{1,2}.fastq.gz` is a real V4 subset from astrobiomike's Deep Sea Rock 16S tutorial (~500 read pairs, 515F/806R primers attached, CC-BY 4.0, attribution in `tests/fixtures/README.md`).

```bash
Rscript skills/claw-amplicon-qc/amplicon_qc.R --demo --output /tmp/demo
```

Expected output: a run that completes in ~1.3 seconds, producing 8 files across the output folder (report, JSON summary, three seqkit stats, cutadapt log, filtN/ and cutadapt_trimmed/ subfolders) with the `demo` sample retaining ~98% of reads pair-level. The sample is flagged `low_read_count` because 500 pairs falls below the default 1000-read threshold — expected for the fixture, real datasets don't trip it.

The end-to-end test suite runs against the same fixture:

```bash
bash tests/run_test.sh
```

Passes 17 assertion checks: prerequisites, output files, output directories, sample count, cutadapt success, retention sanity threshold.

## Algorithm / Methodology

Several methodological choices are baked into this skill. Understanding them helps a researcher — or an orchestrating LLM — know when the defaults are appropriate.

**N removal happens before primer trimming, not with quality filtering.** Cutadapt cannot detect primer sequences in reads containing ambiguous bases. This is not an optimisation, it is a hard prerequisite. Combining N removal with quality filtering (as some tutorials do) makes it impossible to attribute per-stage losses.

**Quality filtering is deliberately deferred to a later skill.** `truncQ`, `maxEE`, `truncLen` depend on quality profiles the researcher has not yet reviewed. Applying defaults here would be either too lenient (wasting downstream compute) or too aggressive (discarding recoverable data). Deferring keeps the human in the loop.

**All four primer positions are trimmed explicitly.** Each paired-end read may contain primer sequence at both ends — the primer that primed the read at 5', and the reverse-complement of the opposite primer at 3' if the read reads through the amplicon. Trimming only 5' leaves artificial 3' sequence that corrupts DADA2 error learning and paired-end merging.

**5' primers are anchored (`-g ^FWD`, `-G ^REV`).** Without anchoring, cutadapt can match a partial primer motif buried mid-read and truncate real biology. For degenerate primers this is especially likely.

**`--discard-untrimmed` is on.** Reads where cutadapt couldn't find a primer are dropped rather than passed through untouched. Primer-free reads (adapter dimer, contamination, primer-region errors) would otherwise poison DADA2's error learning.

**`--pair-filter=any` is explicit.** A pair is discarded if either mate failed primer detection — the strict interpretation. This is cutadapt's current default but explicit here to guard against future upstream changes.

**`--nextseq-trim` is opt-in, not default.** The flag is a two-colour chemistry setting (NextSeq/NovaSeq) that aggressively trims 3' Gs assuming they are dark-cycle artifacts. It is scientifically wrong for MiSeq four-colour chemistry, where Gs are real. Default: flag omitted entirely.

**`filterAndTrim` at Stage 2 is genuinely N-removal-only.** DADA2's defaults for `truncQ` and `minLen` are explicitly overridden to 0, matching this skill's stated behaviour. `rm.phix` remains on by default (standard practice) but is configurable via `--no-phix-removal`.

**Baseline read statistics are captured before any transformation.** This is the reference against which every downstream loss is measured.

**`min_length` is a required user input, not a hardcoded default.** The appropriate minimum length after primer trimming depends on the amplicon region and expected paired-end overlap. A generic default would silently fail for non-V4 primer pairs.

**Anomaly flags are informational, not blocking.** When a sample loses more reads than expected, the skill still completes normally and records the flag. The decision to exclude a flagged sample belongs to the downstream skill or the researcher.

**Key thresholds / parameters (all CLI-configurable):**

- `--max-n = 0` — DADA2 filterAndTrim maxN. Source: cutadapt requires N-free reads for primer detection.
- `--extreme-drop-threshold = 50` — percent retention at any single stage that fires the `extreme_drop` flag.
- `--low-count-threshold = 1000` — reads below which `low_read_count` fires. Lower for low-biomass sample types.
- `--overall-retention-threshold = 70` — overall percent below which `high_overall_loss` fires.

### Automatic Flags

The skill sets these flags in `qc_summary.json` and `report.md`:

| Flag | Trigger | Action for researcher |
|------|---------|----------------------|
| `cutadapt_failed` | Cutadapt returned a non-zero exit status on this sample | Sample is excluded from downstream stats. Inspect `04_cutadapt_log.txt` and decide whether to re-run just this sample or exclude. |
| `extreme_drop` | Sample loses more than the threshold at any single stage on either R1 or R2 | Investigate — primer mismatch, low library quality, or contamination. |
| `low_read_count` | Sample retains fewer than the threshold on either R1 or R2 after primer trimming | Consider excluding — insufficient depth for reliable ASV inference. |
| `high_overall_loss` | Sample retains less than the threshold overall on either R1 or R2 | Investigate — cumulative quality issues. |

Either-side logic means a flag fires if R1 or R2 crosses the threshold — R2-specific problems (common at the 3' end on MiSeq 2×300) are visible even when the pair total looks fine. Individual-sample flags are informational; the 100%-data-loss backstop is separate and aborts only on run-level configuration errors.

## Example Queries

- "Run 16S QC on the FASTQ files in /path/to/raw with 341F/806R primers"
- "Trim primers and remove Ns from my amplicon data before DADA2"
- "Preprocess my V4 microbiome run — reads are in /path/to/raw"
- "Check whether my amplicon data has the right primer orientation before running the pipeline"
- "Show me a demo run so I can see what the outputs look like"

## Example Output

Captured from a real `--demo` invocation. Reproducible: `Rscript amplicon_qc.R --demo --output /tmp/demo`.

```
  claw-amplicon-qc v0.2.0
  Raw folder:     /path/to/tests/fixtures
  Output folder:  /tmp/demo
  Fwd primer:     GTGYCAGCMGCCGCGGTAA
  Rev primer:     GGACTACHVGGGTWTCTAAT
  Min length:     200

Detected naming pattern: Illumina default (SAMPLE_R1_..., SAMPLE_R2_...)
Discovered 1 sample pairs (paired by sample name).

  Preflight — primer orientation check
  Reads sampled (across 1 samples): 500
  R1 starting with forward primer: 498 (99.6%)
  R1 starting with reverse primer: 2 (0.4%)
  Orientation: consistent forward-oriented. OK.

  Stage 2/5 — Removing reads containing N bases
  maxN threshold:        0
  PhiX removal:          ON
  Quality trimming:      OFF (deferred to downstream DADA2 skill)
  Minimum length filter: OFF (--min-length applies at Stage 4 only)

  Stage 4/5 — Primer trimming with Cutadapt
  Fwd: GTGYCAGCMGCCGCGGTAA   Fwd RC: TTACCGCGGCKGCTGRCAC
  Rev: GGACTACHVGGGTWTCTAAT   Rev RC: ATTAGAWACCCBDGTAGTCC
  Min length after trimming: 200 bp
  [1/1] demo

  claw-amplicon-qc — complete
  Samples processed:  1
  Samples failed:     0
  Runtime:            1.3 s
```

Excerpt from the resulting `report.md`:

```markdown
# claw-amplicon-qc — QC Report

> **Disclaimer:** ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.

## Summary
- Samples processed successfully: **1**
- Samples failed (cutadapt): **0**
- Samples flagged: **1**
```

## Output Structure

```
output_folder/
├── report.md                       # Human-readable QC report (includes ClawBio disclaimer)
├── qc_summary.json                 # Machine-readable structured summary
├── 01_raw_stats.txt                # seqkit stats on raw reads
├── 02_filtN_stats.txt              # seqkit stats after N removal
├── 03_trimmed_stats.txt            # seqkit stats after primer trimming
├── 04_cutadapt_log.txt             # Per-sample Cutadapt log (with sample headers)
├── filtN/                          # N-filtered FASTQ files
│   ├── {sample}_R1_001.fastq.gz
│   └── {sample}_R2_001.fastq.gz
└── cutadapt_trimmed/               # Primer-trimmed FASTQ files
    ├── {sample}_R1_001.fastq.gz
    └── {sample}_R2_001.fastq.gz
```

## Dependencies

Install via conda from the bundled `environment.yml`:

```bash
conda env create -f environment.yml
conda activate amplicon-qc
```

Or manually:

```bash
conda create -n amplicon-qc -c conda-forge -c bioconda \
    r-base=4.4 bioconductor-dada2 bioconductor-shortread \
    bioconductor-biostrings r-optparse r-jsonlite \
    cutadapt seqkit
```

**Runtime:** R ≥ 4.4, Bioconductor (DADA2 for `filterAndTrim`, ShortRead + Biostrings for the preflight), CRAN (optparse for CLI, jsonlite for JSON output), plus binaries `cutadapt` (Stage 4) and `seqkit` (Stages 1, 3, 5).

## Gotchas

**Pre-primer-trimmed data will abort at the pre-flight.** If the sequencing facility already stripped primers before delivery, the pre-flight orientation check will report `primers_not_detected` and stop. This is intentional — the pipeline exists to *do* primer trimming; if there's nothing to trim, this is the wrong tool. Go directly to the downstream DADA2 skill.

**Extraction blanks trigger `low_read_count`, correctly.** Blanks are supposed to have almost nothing in them. The default 1000-read floor will fire on them. Not a bug; it means the flag logic works. Either accept the flags on blanks, lower `--low-count-threshold`, or use a metadata scheme to exclude blanks before running.

**Sample names inherit filename suffixes.** A file named `B4_sub_R1_001.fastq.gz` produces sample name `B4_sub` (the `_R1_...` suffix is stripped, but everything before it is kept verbatim). If your files include suffixes like `_sub` or `_downsampled` from prior processing, those end up in report row labels. Rename upstream if you want clean sample IDs.

**Spaces in output paths are safe internally, but shells may still stumble.** The script `shQuote()`s all paths before passing to `system2()`, so `/mnt/d/Research Data/` works fine when invoked directly. If you're wrapping the call in a shell script, quote the path there too.

**`--nextseq-trim` was default in earlier versions but is now opt-in.** If you have saved commands from before 0.2.0 that omit `--nextseq-trim`, the current behaviour (no NextSeq trimming) will differ from before (aggressive 3' G trimming). This mostly *raises* retention on MiSeq data. If you're comparing runs across versions, note this.

**Mixed-orientation libraries abort by default.** Some legitimate library preparations (some ONT amplicon protocols, older 454-derived workflows) produce reads where ~50% of R1 starts with the forward primer and ~50% starts with the reverse. The pre-flight aborts to prevent silent 50% data loss. Use `--allow-mixed-orientation` to proceed if you understand the trade-off.

**`qc_summary.json` schema changed in 0.2.0.** Per-sample count fields are now nested `{R1, R2, pair}` rather than flat integers. Downstream consumers written against 0.1.x need updating. The top-level `totals` keys are stable.

## Safety

- **Local-first.** No data upload without explicit user consent. All processing runs on the invoking machine.
- **Disclaimer.** Every generated `report.md` includes the ClawBio disclaimer: *ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
- **No hallucinated science.** All outputs (read counts, retention percentages, flags) are directly measurable from FASTQ contents. No taxonomic or diversity claims are made anywhere.
- **Fail loud.** Configuration errors, missing tools, mispaired files, wrong primers, and 100% data loss all abort with actionable messages rather than produce a misleading report.

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (R script) executes. The agent must NOT:

- Modify flag thresholds without explicit user request.
- Reclassify flagged samples as unflagged in downstream analysis.
- Proceed to DADA2 denoising without human review of the quality profiles the downstream skill will produce.
- Alter primer sequences the user supplied — always confirm changes with the user first.
- Overrule a pre-flight abort by silently adding `--allow-mixed-orientation` — the abort is a decision point that requires user awareness.

## Integration with Bio Orchestrator

**Trigger conditions.** The orchestrator routes here when:

- File extension `.fastq.gz` (or `.fq.gz`, `.fastq`, `.fq`) with paired-end R1/R2 naming
- Keywords from `metadata.openclaw.trigger_keywords` in the frontmatter

**Chaining partners.**

- **Upstream:** none. This is the earliest preprocessing skill in the amplicon chain.
- **Downstream:** the forthcoming DADA2 quality-filtering skill will consume `cutadapt_trimmed/` (primer-trimmed FASTQ pairs) plus `qc_summary.json` (per-sample flags and retention) and continue with quality profile review, `truncLen` selection, error learning, denoising, and paired-end merging.

Chain-ready outputs: `qc_summary.json` (structured, includes `preflight.primer_orientation`, `totals`, per-sample flags with nested R1/R2/pair counts) and `cutadapt_trimmed/` (the actual downstream input for DADA2).

## Maintenance

- **Review cadence:** on any Cutadapt or DADA2 major release. Verify `--pair-filter` default hasn't changed in Cutadapt; verify `filterAndTrim` parameter defaults haven't shifted in DADA2.
- **Staleness signals:** new primer conventions (dual-index barcoding artifacts, ONT amplicon protocols in Illumina-style workflows), sequencing platform shifts (Illumina NovaSeq X changes), changes to the astrobiomike tutorial dataset that would invalidate the bundled fixture.
- **Deprecation:** when the downstream DADA2 skill can consume outputs from a more integrated pipeline seamlessly, archive this skill to `skills/_deprecated/` with a note pointing at the successor.

## Citations

- [DADA2](https://doi.org/10.1038/nmeth.3869) — Callahan et al. (2016), *Nature Methods* 13:581–583; single-nucleotide-resolution amplicon denoising, provides `filterAndTrim` used at Stage 2.
- [Cutadapt](https://doi.org/10.14806/ej.17.1.200) — Martin (2011), *EMBnet.journal* 17(1):10–12; primer/adapter trimming, used at Stage 4.
- [seqkit](https://doi.org/10.1371/journal.pone.0163962) — Shen et al. (2016), *PLoS ONE* 11(10):e0163962; fast FASTQ statistics, used at Stages 1, 3, 5.

## License

MIT — see LICENSE file at the repository root.
