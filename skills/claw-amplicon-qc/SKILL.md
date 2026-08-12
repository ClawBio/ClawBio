---
name: claw-amplicon-qc
description: >-
  16S/18S rRNA amplicon preprocessing — from raw paired-end FASTQ files
  through N removal and primer trimming, producing outputs ready for DADA2
  quality filtering and denoising. Deliberately stops before quality-filtering
  decisions that require researcher judgment.
license: MIT
metadata:
  version: "0.1.5"
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
  # R-shaped dependency block (this skill is R, not Python)
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
        Tiny synthetic V3-V4 paired-end fixture (~1000 reads, 341F/806R primers
        attached) for end-to-end pipeline verification. Bundled with the skill.
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

<!--
================================================================================
REVIEWER NOTE (Dr. Corpas / audit follow-up)
================================================================================

This SKILL.md is currently mid-migration from the claw-metagenomics sibling-skill
convention (which I originally modeled on) to the templates/SKILL-TEMPLATE.md
convention (which the audit specifies).

Status of this file:
  - FRONTMATTER (above): fully migrated to template format
      - Added: metadata.domain, metadata.dependencies, metadata.demo_data,
        metadata.endpoints.cli, metadata.openclaw.install, and
        metadata.openclaw.trigger_keywords
      - Removed: metadata.openclaw.system_dependencies (superseded by the
        split above)
      - Removed: r-yaml from dependencies (script uses jsonlite, not yaml —
        this was an audit-flagged inconsistency)
      - Version now 0.1.5 (was 0.2.0 in the initial migration). Contributor
        preferred a lower patch-style bump over a minor bump. Happy to
        revert either way.

  - BODY BELOW: NOT YET migrated to template sections. The template's required
    sections (## Trigger, ## Scope, ## Workflow, ## Gotchas, ## Safety,
    ## Agent Boundary, ## Example Queries, ## Integration with Bio Orchestrator,
    ## Maintenance) are not yet present. Existing body content is preserved as-is
    below, with per-section comments flagging what needs to happen.

Two known-inaccurate sections are flagged in-line but preserved rather than
deleted, so you can see exactly what the audit called out:
  - Example Output section: prints "Stage 1/3" and "Multithreading: enabled
    (16 cores)". Script actually prints Stages 1/5 through 5/5 and is
    hardcoded single-threaded. Needs replacing with a captured real run.
  - Testing section: references `bash run_test.sh` and a bundled FASTQ pair.
    Neither currently exists. This is CRITICAL 4 in the audit.

Sections with no direct template equivalent that are preserved below because
they're substantive content I don't want to lose without your call:
  - "Scientific decisions encoded" — could be renamed ## Algorithm / Methodology
  - "Validated On" — no template slot, but useful documentation
  - "Automatic Flags" — no template slot, but useful documentation
  - "What Comes Next" — no template slot, but useful forward-pointer

None of these are wrong; they need a decision on placement.
================================================================================
-->

# 16S Amplicon Preprocessing (claw-amplicon-qc)

Preprocessing pipeline for 16S/18S rRNA amplicon sequencing data — from raw paired-end FASTQ files through primer-trimmed reads ready for DADA2 quality filtering and denoising.

<!--
REVIEWER: The template asks for a persona/role sentence right after the H1,
in the form "You are **[Skill Name]**, a specialised ClawBio agent for [domain].
Your role is to [core function in one sentence]." Not yet added — waiting on
Section B (## Trigger) decisions since they're related.
-->

<!--
================================================================================
SECTION MAPPING TO TEMPLATE (for reviewer reference)
================================================================================
Below is the current SKILL.md body. Section-by-section mapping to what the
template wants:

  "What it does"                → ## Core Capabilities (rename)
  "Why this exists"             → ## Why This Exists (already matches — keep)
  "Scientific decisions encoded"→ ## Algorithm / Methodology (partial, expand)
  "Validated On"                → no template slot (keep or move to Gotchas?)
  "Pipeline Architecture"       → split into ## Workflow + ## Output Structure
  "Automatic Flags"             → no template slot (keep as-is?)
  "Usage" + "Arguments"         → ## CLI Reference (rename, add ClawBio runner)
  "Common primer pairs"         → could stay under CLI Reference or move
  "Example Output"              → ## Example Output (audit: content is wrong)
  "Testing"                     → ## Demo (rewrite when CRITICAL 4 lands)
  "Dependencies"                → ## Dependencies (already matches — keep)
  "What Comes Next"             → no template slot (keep or drop?)
  "Citations"                   → ## Citations (already matches — keep)

MISSING (template requires, not present below):
  ## Trigger
  ## Scope
  ## Workflow (as a named section — content is in Pipeline Architecture)
  ## Gotchas
  ## Safety
  ## Agent Boundary
  ## Example Queries
  ## Output Structure (as a named section)
  ## Integration with Bio Orchestrator
  ## Maintenance
================================================================================
-->

## What it does

<!-- REVIEWER: rename to ## Core Capabilities per template. Content is fine. -->

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

<!--
REVIEWER: Item 4 mentions --nextseq-trim=20 as "for Illumina two-colour chemistry
poly-G trimming". The audit flagged this as scientifically wrong for the stated
platform (MiSeq is four-colour, not two-colour). The description above still
reflects the current script behaviour, but both the code AND this description
should change when we address the method issues. Not touching in this Section-A
pass.
-->

## Why this exists

If you ask a general AI to "run 16S QC," it will typically:

- Combine N removal and quality filtering into a single step, hiding which reads were lost to which criterion
- Use default DADA2 quality parameters (`maxEE`, `truncQ`) before the researcher has reviewed quality profiles
- Skip the primer sanity check and produce silent failures downstream if primers were missed
- Not distinguish between the forward primer and its reverse-complement, which matters for correct 3' trimming
- Not track per-sample read survival across stages, making it hard to spot problematic samples
- Not produce a machine-readable summary that a downstream tool or AI wrapper can consume

<!--
REVIEWER: This section maps cleanly to the template's ## Why This Exists. Content
is well-argued (audit specifically praised the reasoning here). Keep as-is.
-->

## Scientific decisions encoded

<!--
REVIEWER: This has no direct template slot. The template's ## Algorithm /
Methodology section is the closest fit but wants a numbered step-by-step of
the algorithm plus a "Key thresholds / parameters" list with sources.
This existing content is a mix of methodological rationale and design decisions.
Suggest either:
  (a) Rename to ## Algorithm / Methodology and expand with the numbered
      workflow (which currently lives in Pipeline Architecture below)
  (b) Keep as its own subsection under a broader ## Methodology section
Not renaming in this pass — deferring your call.
-->

Several methodological choices are baked into this skill. Understanding them helps a researcher (or an AI wrapper) know when the defaults are appropriate and when they are not.

- **N removal happens before primer trimming, not with quality filtering.** Cutadapt cannot detect primer sequences in reads containing ambiguous bases (N). This is not an optimisation — it is a hard prerequisite. Combining N removal with quality filtering (as some tutorials do) makes it impossible to know which reads were lost to which criterion.

- **Quality filtering (`maxEE`, `truncQ`, `truncLen`) is deliberately deferred to a later skill.** These parameters depend on quality profiles the researcher has not yet reviewed. Applying default values here would either be too lenient (wasting downstream compute on unusable reads) or too aggressive (discarding recoverable data). Deferring keeps the human in the loop for a decision that requires their judgment.

- **All four primer positions are trimmed explicitly.** In paired-end amplicon libraries, each read may contain primer sequence at both ends — the primer that primed the read at the 5' end, and the reverse-complement of the opposite primer at the 3' end if the read is long enough to read through the amplicon. Trimming only the 5' primer leaves artificial sequence at the 3' end that would corrupt DADA2 error learning and paired-end merging.

- **Baseline read statistics are captured before any transformation.** This establishes the reference point against which every downstream loss is measured, so retention rates are meaningful and per-sample anomalies are visible from the earliest stage.

- **`min_length` is a required user input, not a hard-coded default.** The appropriate minimum length after primer trimming depends on the amplicon region (V3–V4 vs V4 differ by ~200bp) and the expected paired-end overlap. A generic default would produce silent failures for non-V4 primer pairs.

- **Anomaly flags are informational, not blocking.** When a sample loses more reads than expected, the skill still completes normally and records the flag in `qc_summary.json` and `report.md`. The decision to exclude a flagged sample belongs to the downstream skill or the researcher, not to this preprocessing step.

<!--
REVIEWER: Note the audit's finding on hardcoded thresholds. The last bullet
above argues persuasively that min_length must not be hardcoded. The audit
correctly extends the same argument to the four thresholds baked into the
flag logic (50% drop, 1000-read floor, 70% overall retention, and the maxN
value implied by filterAndTrim). Those thresholds should also become
configurable in a future iteration. Content preserved as-is here.
-->

## Validated On

<!--
REVIEWER: No template slot. Substantive content — keep, but placement is your
call. Could stay standalone or roll into ## Demo / ## Example Output.
Verification-history could also go under a ## Maintenance section since the
template asks for one.

ALSO — the audit will now be able to note this section can grow. As of the
audit-fix pass, the skill has been verified on:
  - Original environmental 16S set (20 samples, 2×300 MiSeq, 341F/806R)
  - Additional environmental set (25 samples, but pre-trimmed 2×150 — used
    to verify the 100%-data-loss backstop from CRITICAL 3)
  - astrobiomike deep-sea rock dataset (20 samples, 2×300 MiSeq, 515F/806R),
    published tutorial with expected outputs — retention within 5% of
    tutorial's published numbers.
The audit-verified numbers deserve a mention here or in ## Demo.
-->

Environmental 16S samples (aquatic plant roots, fronds, and surrounding water) using 341F/806R primers on Illumina MiSeq 2×300 paired-end sequencing. Sample sizes tested: 20 to 100+ samples per run.

## Pipeline Architecture

<!--
REVIEWER: Template wants this split into ## Workflow (numbered steps) and
## Output Structure (tree). This ASCII diagram serves both roles today.
Suggest keeping the diagram and adding both named sections around it, or
converting to two sections cleanly. Not restructured in this pass.
-->

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
        │                ──►  04_cutadapt_log.txt
        │                ──►  03_trimmed_stats.txt
        ▼
        Retention analysis + anomaly flagging
        │                ──►  report.md (human-readable)
        │                ──►  qc_summary.json (machine-readable)
        ▼
        Ready for DADA2 quality profiling and denoising
        (separate skill — claw-amplicon-dada2, forthcoming)
```

<!--
REVIEWER: DEFECT (re-audit): the diagram previously showed
"03_cutadapt_log.txt" while the script writes "04_cutadapt_log.txt".
Fixed in this pass — the diagram entry now matches. Stage numbering in
the diagram (Stages 1-3 conceptually, showing what happens where) is
kept separate from the terminal Stage 1/5-5/5 which reflects the actual
number of seqkit calls; both are accurate views of the same pipeline.
-->

## Automatic Flags

<!--
REVIEWER: No template slot. Substantive documentation — keep. Could also live
under ## Algorithm / Methodology as "Key thresholds/parameters" per template.

Note the audit's point: these thresholds are hardcoded and uncited. They
should be configurable in a future iteration. Also, as of CRITICAL 3, there
is a NEW flag (`cutadapt_failed`) that fires when a per-sample cutadapt run
returns a non-zero exit status. Not yet added to this table.
-->

The skill sets flags in `qc_summary.json` for any sample matching these conditions:

| Flag | Trigger | Action for researcher |
|------|---------|----------------------|
| `extreme_drop` | Sample loses >50% of reads at any single stage | Investigate — likely primer mismatch, low library quality, or contamination |
| `low_read_count` | Sample retains <1000 reads after primer trimming | Consider excluding — insufficient depth for reliable ASV inference |
| `high_overall_loss` | Sample retains <70% of raw reads overall | Investigate — cumulative quality issues |

<!--
REVIEWER: Missing from table — the `cutadapt_failed` flag added in CRITICAL 3.
Should be added:

| `cutadapt_failed` | Cutadapt returned a non-zero exit status on this sample | Sample is excluded from downstream stats; see 04_cutadapt_log.txt for the error and decide whether to re-run just this sample or exclude entirely |

Not adding in this pass; queuing for the SKILL.md content-review sections.
-->

Flags are informational, not blocking. The skill completes normally and reports flagged samples in `report.md` and `qc_summary.json`. The downstream DADA2 skill (or the researcher) decides how to handle them.

<!--
REVIEWER: Above claim ("skill completes normally") is now partially inaccurate
as of CRITICAL 3. If cutadapt fails on a specific sample, the skill continues
(behaviour (b)) — that part is still true. But if the whole run appears to be
a configuration error (>99% loss across all successful samples), the skill
now aborts with a clear message rather than write a report of zeros. The
claim should be qualified: "individual-sample flags are informational; run-
level aborts fire only on catastrophic configuration errors."
-->

## Usage

<!-- REVIEWER: rename to ## CLI Reference per template. -->

```bash
Rscript amplicon_qc.R \
    --raw /path/to/raw_fastq_folder \
    --output /path/to/output_folder \
    --fwd-primer CCTACGGGNGGCWGCAG \
    --rev-primer GACTACHVGGGTWTCTAAT \
    --min-length 200
```

<!--
REVIEWER: Template also asks for the ClawBio runner invocation, e.g.:
  python clawbio.py run claw-amplicon-qc --raw <dir> --output <dir> ...
And the --demo mode invocation. Both need adding once CRITICAL 4 (--demo)
lands. Placeholder for now.
-->

### Arguments

**Required:**

| Argument | Description |
|----------|-------------|
| `--raw` | Folder containing paired-end FASTQ files |
| `--output` | Folder for outputs (created if missing) |
| `--fwd-primer` | Forward primer sequence (IUPAC codes) |
| `--rev-primer` | Reverse primer sequence (IUPAC codes) |
| `--min-length` | Minimum read length in bp after primer trimming |

**Optional tuning (all with defaults preserving prior behaviour):**

| Argument | Default | Description |
|----------|---------|-------------|
| `--max-n` | `0` | Max N bases allowed at Stage 2 filterAndTrim. Default 0 = discard any read with N. |
| `--nextseq-trim` | `0` | Cutadapt `--nextseq-trim=` value for two-colour chemistry (NextSeq/NovaSeq). Default 0 = flag omitted. Set to a positive integer (e.g. 20) only for NextSeq/NovaSeq data — not for MiSeq. |
| `--no-phix-removal` | (off) | Flag. Skip PhiX removal in Stage 2 (default: PhiX is removed via `rm.phix=TRUE`). |
| `--allow-mixed-orientation` | (off) | Flag. Continue when the pre-flight orientation check detects mixed-orientation reads. Default: abort. Only use if you understand the consequences (~50% data loss). |
| `--extreme-drop-threshold` | `50` | Percent retention below which the `extreme_drop` flag fires. |
| `--low-count-threshold` | `1000` | Read count below which the `low_read_count` flag fires. Lower for low-biomass samples. |
| `--overall-retention-threshold` | `70` | Overall retention percent below which the `high_overall_loss` flag fires. |

**Demo mode:**

| Argument | Description |
|----------|-------------|
| `--demo` | Run on the bundled test fixture (`tests/fixtures/demo_R{1,2}.fastq.gz`) with hardcoded 515F/806R primers and `--min-length 200`. Requires `--output` only; refuses to be combined with content flags. |

### Common primer pairs (for reference)

The skill does not hard-code any primer pairs — the researcher must provide sequences. Common choices include:

| Region | Forward | Sequence | Reverse | Sequence |
|--------|---------|----------|---------|----------|
| V3–V4 | 341F | `CCTACGGGNGGCWGCAG` | 806R (Caporaso) | `GACTACHVGGGTWTCTAAT` |
| V4 | 515F (Parada) | `GTGYCAGCMGCCGCGGTAA` | 806R (Apprill/EMP) | `GGACTACNVGGGTWTCTAAT` |
| V1–V2 | 27F | `AGAGTTTGATCMTGGCTCAG` | 338R | `TGCTGCCTCCCGTAGGAGT` |

## Example Output

<!--
REVIEWER: AUDIT-FLAGGED — content below is stale and does not match the
script's actual behaviour:
  - Says "Stage 1/3" — script actually prints Stages 1/5 through 5/5
  - Says "Multithreading: enabled (16 cores)" — script hardcodes single-
    threaded (multithread <- FALSE at line 213 of amplicon_qc.R)
  - Shows made-up flag messages — real messages differ slightly

Preserved as-is here (rather than deleted) so the audit finding is
inspectable. Should be replaced with a captured real run from
raw_qc_output_v2/ (the CRITICAL-3-verified reference run on 20 samples).
REGENERATED as of the 0.1.5 push — captured from --demo mode against the
bundled fixture. Reproducible by anyone via `Rscript amplicon_qc.R --demo
--output /tmp/demo`.
-->

```
─────────────────────────────────────────────────────────
  claw-amplicon-qc v0.1.5
─────────────────────────────────────────────────────────
  Raw folder:     /path/to/tests/fixtures
  Output folder:  /tmp/demo
  Fwd primer:     GTGYCAGCMGCCGCGGTAA
  Rev primer:     GGACTACHVGGGTWTCTAAT
  Min length:     200
─────────────────────────────────────────────────────────
Detected naming pattern: Illumina default (SAMPLE_R1_..., SAMPLE_R2_...)
Discovered 1 sample pairs (paired by sample name).

─────────────────────────────────────────────────────────
  Preflight — primer orientation check
─────────────────────────────────────────────────────────
  Reads sampled (across 1 samples): 500
  R1 starting with forward primer: 498 (99.6%)
  R1 starting with reverse primer: 2 (0.4%)
  Orientation: consistent forward-oriented. OK.

─────────────────────────────────────────────────────────
  Stage 1/5 — Baseline seqkit stats on raw reads
─────────────────────────────────────────────────────────
  Total raw reads: 1,000
  Stats saved to:  /tmp/demo/01_raw_stats.txt

─────────────────────────────────────────────────────────
  Stage 2/5 — Removing reads containing N bases
─────────────────────────────────────────────────────────
  maxN threshold:        0
  PhiX removal:          ON
  Quality trimming:      OFF (deferred to downstream DADA2 skill)
  Minimum length filter: OFF (--min-length applies at Stage 4 only)
Read in 500 paired-sequences, output 497 (99.4%) filtered paired-sequences.

─────────────────────────────────────────────────────────
  Stage 3/5 — seqkit stats on N-filtered reads
─────────────────────────────────────────────────────────
  Total reads after N removal: 994
  Stats saved to:  /tmp/demo/02_filtN_stats.txt

─────────────────────────────────────────────────────────
  Stage 4/5 — Primer trimming with Cutadapt
─────────────────────────────────────────────────────────
  Fwd: GTGYCAGCMGCCGCGGTAA   Fwd RC: TTACCGCGGCKGCTGRCAC
  Rev: GGACTACHVGGGTWTCTAAT   Rev RC: ATTAGAWACCCBDGTAGTCC
  Min length after trimming: 200 bp
  [1/1] demo

─────────────────────────────────────────────────────────
  Stage 5/5 — seqkit stats on primer-trimmed reads
─────────────────────────────────────────────────────────
  Total reads after primer trimming: 984
  Stats saved to:  /tmp/demo/03_trimmed_stats.txt

─────────────────────────────────────────────────────────
  Retention analysis + anomaly flagging
─────────────────────────────────────────────────────────
  Samples processed successfully:  1
  Samples failed (cutadapt):       0
  Samples flagged:                 1
  Flagged:                         demo

═════════════════════════════════════════════════════════
  claw-amplicon-qc — complete
═════════════════════════════════════════════════════════
  Samples processed:  1
  Samples failed:     0
  Runtime:            1.3 s
  Output folder:      /tmp/demo
  Machine-readable:   /tmp/demo/qc_summary.json
  Human-readable:     /tmp/demo/report.md
═════════════════════════════════════════════════════════
```

(The `demo` sample is flagged `low_read_count` because the bundled fixture
contains only 500 read pairs, below the default 1000-read threshold. This
is expected for the demo — real datasets won't trip that flag.)

## Testing

<!--
REVIEWER: AUDIT-FLAGGED (CRITICAL 4) — content below refers to
`bash run_test.sh` and a bundled FASTQ pair, neither of which currently
exists. The test file that DOES exist (tests/test_amplicon_qc.py) targets
`amplicon_qc.py` which also doesn't exist (only amplicon_qc.R ships).

This section needs to be rewritten AFTER CRITICAL 4 is addressed (bundling
a real fixture at tests/fixtures/ and writing a shell test that drives
Rscript against it). Preserved as-is for now.
-->

A basic end-to-end test lives in `tests/` and can be run with:

```bash
cd tests
bash run_test.sh
```

The test uses a tiny bundled FASTQ pair (~1000 reads per sample) with known primer sequences to verify the pipeline runs to completion and produces expected output files.

## Dependencies

<!--
REVIEWER: The frontmatter now has the canonical dependency information
(metadata.dependencies + metadata.openclaw.install). This body section is
therefore duplicative. The template's ## Dependencies section is meant to
be a human-readable summary rather than a full spec. Content below is fine
but could be simplified. Not touching in this pass.
-->

Install via conda:

```bash
conda create -n amplicon-qc -c conda-forge -c bioconda \
    r-base=4.4 \
    bioconductor-dada2 \
    bioconductor-shortread \
    bioconductor-biostrings \
    r-optparse \
    r-jsonlite \
    cutadapt \
    seqkit
```

Or via the bundled environment file:

```bash
conda env create -f environment.yml
conda activate amplicon-qc
```

<!--
REVIEWER: DEFECT 3 (re-audit): conda command now lists r-optparse and
r-jsonlite explicitly, r-yaml removed. Matches environment.yml exactly.
Frontmatter, environment.yml, and this body block are now all in sync.
-->

## What Comes Next

<!--
REVIEWER: No template slot. Substantive forward-pointer content — keep.
Could roll into ## Integration with Bio Orchestrator (chaining partners) or
## Maintenance. Or keep as-is at the end of the file.
-->

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

<!--
================================================================================
END-OF-FILE REVIEWER SUMMARY
================================================================================

Sections from templates/SKILL-TEMPLATE.md that are NOT YET in this file:

  ## Trigger
      - The routing wiring for the Bio Orchestrator. Frontmatter has
        trigger_keywords (functional wiring); the ## Trigger body section
        wants explicit "Fire when..." and "Do NOT fire when..." lists.
      - Draft candidates for "Fire when":
          * "run 16S QC on my FASTQs"
          * "trim primers from my amplicon data"
          * "preprocess my 16S/18S run for DADA2"
          * "remove Ns and primers from my paired-end FASTQ files"
          * "quality-check my microbiome sequencing"
      - Draft candidates for "Do NOT fire when":
          * User asks for quality filtering / truncLen decisions
            → route to claw-amplicon-dada2 (forthcoming)
          * User asks for taxonomy assignment or ASV clustering
            → route to claw-amplicon-dada2 or claw-amplicon-taxonomy
          * User asks for shotgun metagenomics / assembly
            → route to claw-metagenomics
          * User data is pre-primer-trimmed (skill will detect this at
            run time and abort via 100%-data-loss backstop, but the
            trigger section can pre-empt)

  ## Scope
      - One paragraph: "This skill does X and nothing else." Ready to draft.

  ## Gotchas
      - Only Zabiulla can populate this authoritatively. Candidates from the
        audit-fix debugging session include:
          * Pre-primer-trimmed data will fail the run (raw2 experience)
          * Extraction blanks will trigger low_read_count (astrobiomike)
          * Sample name convention picks up "_sub" suffixes when files are
            named SAMPLE_sub_R1_001.fastq.gz (astrobiomike)
          * Spaces in output paths need double-quoting on the command line
            even though the script now handles them internally
          * --nextseq-trim=20 on MiSeq is scientifically wrong (audit note)
          * `-g` (unanchored 5') means degenerate primers can match mid-read

  ## Safety
      - Boilerplate: "Local-first: No data upload without explicit consent /
        Disclaimer: Every report includes the ClawBio medical disclaimer /
        Audit trail: Log all operations to reproducibility bundle /
        No hallucinated science: All parameters trace to cited databases"
      - Note: this skill's report.md does not currently include the medical
        disclaimer. Should be added when the SKILL.md content-review pass
        happens.

  ## Agent Boundary
      - Draft: "The agent (LLM) explains this skill's outputs and helps the
        researcher interpret retention rates and flags. The agent MUST NOT:
        modify thresholds; re-classify flagged samples as unflagged; proceed
        to DADA2 denoising without human review of quality profiles; alter
        primer sequences without user confirmation."

  ## Example Queries
      - "Run 16S QC on the FASTQ files in /path/to/raw with 341F/806R
        primers"
      - "Trim primers and remove Ns from my amplicon data before DADA2"
      - "Preprocess my V4 microbiome run — reads are in /path/to/raw"

  ## Output Structure
      - The Pipeline Architecture diagram serves this partially. The
        template wants an explicit tree with a testable output contract.
        Should be:

        output_directory/
        ├── report.md                      # Primary markdown report
        ├── qc_summary.json                # Machine-readable summary
        ├── 01_raw_stats.txt               # seqkit stats on raw
        ├── 02_filtN_stats.txt             # seqkit stats after N removal
        ├── 03_trimmed_stats.txt           # seqkit stats after cutadapt
        ├── 04_cutadapt_log.txt            # Cutadapt per-sample log
        ├── filtN/                         # N-filtered FASTQ files
        │   ├── {sample}_R1_001.fastq.gz
        │   └── {sample}_R2_001.fastq.gz
        └── cutadapt_trimmed/              # Primer-trimmed FASTQ files
            ├── {sample}_R1_001.fastq.gz
            └── {sample}_R2_001.fastq.gz

      - Note: no reproducibility/ subdirectory yet. Audit noted this is
        expected by every ClawBio skill ("commands.sh, environment.yml,
        SHA-256 checksums"). Should be added.

  ## Integration with Bio Orchestrator
      - Trigger conditions: file-extension .fastq.gz + keywords listed in
        frontmatter trigger_keywords
      - Chaining partners:
          * upstream: nothing (this is the earliest processing skill)
          * downstream: claw-amplicon-dada2 (forthcoming; consumes
            cutadapt_trimmed/ folder + qc_summary.json)

  ## Maintenance
      - Template asks for review cadence and staleness signals. Draft:
          * Review cadence: on cutadapt or DADA2 major release
          * Staleness signals: new primer conventions (e.g. dual-index
            barcoding artifacts), sequencing platform shifts (e.g. Illumina
            NovaSeq X, ONT for 16S)
          * Deprecation trigger: when claw-amplicon-dada2 is complete and
            can consume the outputs seamlessly

Also missing per audit:
  - Registration in skills/catalog.json (functional gap — skill is
    undiscoverable by the CLI without this)
  - Registration in clawbio/cli.py (same)
  - environment.yml file at repo root or skill root — audit flagged this
    file is missing r-jsonlite. Not visible in the current SKILL.md but
    would live alongside.

================================================================================
-->
