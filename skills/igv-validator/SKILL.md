---
name: igv-validator
description: >-
  Read-level validation of variant calls from any VCF, tumor/normal or tumor-only:
  counts read support in the BAMs (SNVs, indels, SV breakpoints), compares it with
  the caller's reported counts, flags artifact patterns, takes one IGV screenshot
  per variant and writes a report.
license: MIT
metadata:
  version: "0.1.0"
  author: mgdouq-sudo
  domain: genomics
  tags:
    - igv
    - variant-validation
    - somatic
    - tumor-normal
    - bam
    - structural-variants
  inputs:
    - name: vcf
      type: file
      format:
        - vcf
        - vcf.gz
      description: Somatic calls from any caller (Mutect2, Strelka, DRAGEN, Manta tested)
      required: true
    - name: tumor_bam
      type: file
      format:
        - bam
        - cram
      description: Tumor alignments, indexed (.bai/.crai)
      required: true
    - name: normal_bam
      type: file
      format:
        - bam
        - cram
      description: Matched normal alignments, indexed (omit for tumor-only mode)
      required: false
    - name: reference
      type: file
      format:
        - fasta
      description: Indexed reference FASTA the BAMs were aligned to (needed for screenshots and indel clip rescue)
      required: false
    - name: cnv_segments
      type: file
      format:
        - tsv
        - seg
      description: Copy-number mode - GATK called.seg or a table with contig/start/end/log2/cn_call (with --regions)
      required: false
    - name: variants
      type: file
      format:
        - tsv
      description: Optional list of chrom/pos/ref/alt to check (instead of --genes / --pass-only)
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
        - html
      description: Validation report with per-variant counts, flags and screenshots
    - name: result
      type: file
      format:
        - json
      description: Machine-readable per-variant support counts and flags
    - name: support_counts
      type: file
      format:
        - tsv
      description: Tumor/normal read support table
  dependencies:
    python: ">=3.10"
    packages:
      - pysam>=0.22
      - Pillow>=10.1
  demo_data:
    - path: demo/demo_calls.vcf
      description: Synthetic somatic calls (6 variants) on synthetic contigs, with matching synthetic BAMs
  endpoints:
    cli: python skills/igv-validator/igv_validator.py --vcf {vcf} --tumor {tumor_bam} --normal {normal_bam} --reference {reference} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🔬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: pysam
      - kind: pip
        package: Pillow
    trigger_keywords:
      - validate in IGV
      - IGV validation
      - check variant calls in IGV
      - is this mutation real
      - read support for variant
      - tumor normal read counts
      - IGV screenshot of variant
      - validate somatic calls against BAM
      - do the caller calls agree with the reads
      - tumor-only IGV check
---

# 🔬 IGV Validator

You are **IGV Validator**, a specialised ClawBio agent for somatic variant validation. Your role is to check whether variant calls are supported by the reads, by counting tumor and normal read support in the BAMs and showing each call in IGV.

## Trigger

**Fire this skill when the user says any of:**
- "validate these mutations in IGV" / "IGV validation"
- "check my variant calls in IGV" / "look at this variant in IGV"
- "is this mutation real or an artifact?"
- "how many reads support this variant in tumor and normal?"
- "take IGV screenshots of my variants"
- "validate somatic calls against the BAM" / "manual review of calls"
- "check KRAS / <gene> mutations in the reads"
- "do my caller's calls agree with IGV?" / "tumor-only validation"

**Do NOT fire when:**
- The user wants variants *called* from FASTQ/BAM: route to `nfcore-sarek-wrapper`
- The user wants functional annotation (VEP, ClinVar, gnomAD): route to `variant-annotation` or `vcf-annotator`
- The user wants ACMG pathogenicity classification: route to `clinical-variant-reporter` / `cnv-acmg-classifier`
- The user wants sequencing QC (FastQC, coverage summaries): route to `multiqc-reporter` / `seq-wrangler`
- Germline variant review in a single normal sample: the flags assume a tumor sample

## Why This Exists

**Two uses:**
1. **Faster review and reporting.** Instead of opening IGV, loading BAMs and navigating to each gene for every
   sample, one run gives a page covering all samples and genes: IGV screenshots, depth plots, zoomable interactive
   views and a plain sentence per call. It is reproducible, and a PI or collaborator can review it without IGV.
2. **Catching errors in calls and pipelines.** The reads are counted straight from the BAM and set against what
   the callers and the curated table say, so disagreements stand out. Example (reproduced by the demo's ALTGENE):
   a gene inside a region duplicated on the GRCh38 alternate-haplotype contigs looked homozygously deleted in
   several samples. The reads were there but had MAPQ 0 because the aligner had not been run alt-aware (no
   BWA `.alt` file), so every MAPQ-filtering caller saw no coverage. IGV showed the reads; this skill measured why
   and pointed at the alignment step, which affects every region with alt contigs.

**What it isn't:** a variant caller or an automatic truth. Copy number from read depth is a simple measure,
tumor-only data cannot separate inherited from somatic variants, and every verdict is a first pass. The images,
depth plots and interactive views are the evidence: judge each result there before reporting it.

- **Without it**: reviewers open IGV by hand, type each locus, eyeball reads and guess counts from pictures
- **With it**: every call gets exact read counts from the BAM, artifact flags, the caller's own counts next to
  them, and standard images, in minutes
- **Why ClawBio**: counts come from pysam with fixed filters and were cross-checked against samtools and caller
  `AD`; the screenshot shows the evidence, it is never the source of numbers

## Core Capabilities

1. **Variant selection**: by gene (`--genes`, from VEP `CSQ`, SnpEff `ANN`, ANNOVAR `Gene.refGene`/`Gene.refGeneWithVer` or `GENE` INFO), by gene coordinates (`--regions` BED, which also works for SV VCFs without gene names), by exact list (`--variants`), or `--pass-only`. The VCF is streamed, so multi-GB whole-genome files are fine
2. **Read counting**: SNVs (pileup), indels (CIGAR match plus soft-clip rescue at the indel edge), SV breakpoints (split reads plus discordant pairs); each DNA molecule counted once
3. **Artifact flags**: support in normal, fewer than 3 reads, under 5% of reads, single-strand support, alt bases mostly near read ends, low depth, a germline allele already at the site, repeat-like excess depth
4. **Caller comparison**: the caller's own counts (FORMAT `AD`, else `AF` x `DP`) for the tumor's VCF column, next to the BAM counts, with a `caller_disagrees` flag
5. **Two modes**: tumor + matched normal, or tumor-only (omit `--normal`); the normal-based checks are skipped in tumor-only mode
6. **Copy-number mode** (`--cnv`): for each gene in `--regions`, the GATK segment call next to the read depth inside the gene versus its flanks, with a depth plot, a `cnv_disagrees` flag, and IGV screenshots for genes that fit a window
7. **IGV screenshots**: an isolated IGV (own port, own settings folder) on the given reference; one tumor-over-normal (or tumor-only) image per variant, captioned with the BAM and caller counts; closed afterwards
8. **Summary** (`--summarize <folder>`): one page over every run's `result.json` under a folder. With `--curated-calls calls.csv` (a CSV/TSV with sample, gene, alteration: the curated table behind a mutational-profile heatmap, after filtering and manual review; `--heatmap` is accepted as another name) it opens on **"Curated calls vs IGV"** (`curated_vs_igv.tsv`): per sample and gene, what the curated call says, what IGV shows in one plain sentence, *Start here* (*Look first* where the reads and the call differ, *Consistent so far* where they fit; never a verdict), a thumbnail and links to the report, overview and interactive view, with a *details* expander for the comparison, reason and review priority. The comparison follows common curated-table rules (SNV/SV before copy number; DEL/AMP only if under 1 Mb or beyond log2 2), judges copy number from the read depth (per GATK segment, reporting depth steps inside a gene), and explains a WT call over a recurrent germline/artifact variant as consistent (recurrent = the same SNV/indel in 3+ samples, each with at least 2 supporting reads, so low-level copies count). Below it, folded, **"Raw calls vs IGV"** (`igv_agreement.tsv`): each raw caller call before filtering, in one plain sentence, with the read counts and flags under *details*. A row's *report* link opens that check's report showing only that gene (a "show all genes" link returns to the full report). Without `--curated-calls`, the raw calls are the main table. `--overview --regions genes.bed` adds one IGV image per gene and sample; `--interactive` adds zoomable pages; `summary.tsv` keeps every call and field
9. **Report**: `report.md`, self-contained `report.html`, `result.json`, counts TSV and a reproducibility bundle

## Scope

**One skill, one task.** This skill checks existing calls against the reads and nothing else. It does not call, annotate or classify variants.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| VCF | `.vcf`, `.vcf.gz` | CHROM, POS, REF, ALT (BND ALT or SVTYPE/END for SVs) | `demo/demo_calls.vcf` |
| BAM/CRAM | `.bam` + `.bai` | coordinate-sorted, indexed; tumor required, normal optional | `demo/demo_tumor.bam` |
| Reference | `.fa` + `.fai` | same build as the BAMs | `demo/demo_ref.fa` |
| Variant list | `.tsv` | `chrom pos ref alt` (header optional) | `demo/demo_variants.tsv` |
| Regions | `.bed` | `chrom start end [name]`; the name labels the gene | gene bodies from GENCODE or ANNOVAR refGene |

## Workflow

1. **Validate** (prescriptive): VCF, both BAMs and their indexes exist; contig names in the VCF are in the BAM headers; reference (if given) has a `.fai`. Stop with a clear message otherwise.
2. **Select** (prescriptive): apply `--variants` / `--regions` / `--genes` (combined with AND) and `--pass-only`; cap at `--max-variants` (default 50) and say how many were skipped.
3. **Classify** (prescriptive): SNV, insertion, deletion, BND (bracket notation, or symbolic `<TRA>`/`<BND>` with `CHR2`/`END` as SURVIVOR and Delly write) or symbolic `<DEL>/<DUP>/<INV>` of at least 1 kb. Other records (MNVs, small symbolic SVs, CNV-only records) are listed as "not evaluable", never counted. Deletions and insertions written out as sequence (as Manta and SURVIVOR often do) are counted like indels, and from 50 bp they count as SVs in the summary and the curated comparison.
4. **Count** (prescriptive): run the counters with MAPQ >= 20, base quality >= 20, no duplicate/secondary/supplementary/QC-fail reads.
5. **Flag** (prescriptive): apply the rules in *Algorithm*; derive the automatic status.
6. **Screenshot** (prescriptive): unless `--no-igv` or IGV is not found, start the isolated IGV, take one image per variant (two for SVs), check each image rendered, caption it, close IGV.
7. **Report** (flexible): write the outputs. When presenting results, the agent may add its own reading of the screenshots, labelled as a reviewer judgement and kept separate from the automatic status.

## Modes

Ten modes, each with a prompt and a command, are in **[GUIDE.md](GUIDE.md#modes)**: variant checks (tumor + normal or
tumor-only), SVs, copy number, summaries over many samples, the curated-calls comparison, overview images and
interactive views, sharing the report, counts without IGV, the demo, and a whole project in one command
(`--samplesheet`). Agents: read GUIDE.md before running the skill for a user.

## CLI Reference

```bash
# One check: a VCF against a tumor BAM (add --normal for a pair)
python skills/igv-validator/igv_validator.py --vcf calls.vcf --tumor tumor.bam --reference hg38.fa \
  --regions genes.bed --output reports/sampleA/snv

# A whole project in one command (see GUIDE.md)
python skills/igv-validator/igv_validator.py --samplesheet samples.csv --reference hg38.fa \
  --curated-calls curated.csv --annotation gencode.v44.basic.annotation.gtf.gz

# Demo mode (synthetic data, no user files needed)
python skills/igv-validator/igv_validator.py --demo --output /tmp/igv_demo

# Via the ClawBio runner
python clawbio.py run igv-validator --demo
```

| Option | What it does |
|---|---|
| `--vcf` (`--input`) | VCF/VCF.gz of SNV, indel or SV calls (any caller) |
| `--tumor` | tumor BAM/CRAM, indexed (required except for `--summarize`) |
| `--normal` | matched normal BAM/CRAM; leave out for tumor-only |
| `--tumor-only` | ignore any normal (e.g. with `--demo`) |
| `--reference` | FASTA the BAMs were aligned to (screenshots, indel rescue, CRAM) |
| `--genes` | gene symbols. With `--samplesheet`: the genes to check, found in `--annotation` (default: the curated calls' genes). In a single check: selects VCF records by their gene names (VEP, SnpEff, ANNOVAR, GENE) |
| `--regions` | BED of genes (chrom, start, end, name): selects variants/SVs, sets copy-number genes, overview and interactive windows |
| `--variants` | TSV of chrom, pos, [ref, alt] to check exactly (e.g. a filtered table) |
| `--pass-only` | only FILTER=PASS calls |
| `--max-variants` | cap per run (default 50) |
| `--caller-sample` | VCF column whose caller counts to compare (default: the tumor's sample name) |
| `--tumor-name`, `--normal-name` | labels in reports and images (use your sample names) |
| `--cnv`, `--cnv-sample` | copy-number mode: segment file, and the sample in a multi-sample file |
| `--samplesheet` | one command for a whole run: every check for every sample, into a new dated folder with its summary |
| `--skip-unknown-genes` | with `--samplesheet`: run without gene names the annotation does not know (default: stop with suggestions) |
| `--reports-dir` | with `--samplesheet`: where runs are kept (default `igv_reports/`) |
| `--summarize` | build the summary page over every run in a folder (full: with images, interactive pages, zip) |
| `--index` | write `<folder>/index.html` listing every batch folder with a summary; it then keeps itself updated |
| `--project` | project folder of all runs: its `summary.html` is refreshed after this check (automatic for later checks once the folder has a summary) |
| `--curated-calls` (`--heatmap`) | curated table to compare with IGV (sample, gene, alteration) |
| `--overview` | with `--summarize` and `--regions`: one IGV image per gene and sample |
| `--interactive` | with `--summarize` and `--regions`: zoomable pages per sample and gene (needs igv-reports; contain read data) |
| `--annotation` | GTF/GFF3/BED of genes (e.g. GENCODE) drawn as a genes track in every image |
| `--no-bundle` | with `--summarize`: skip `igv_validation_full.zip` (made by default: the whole report, linked from the page) |
| `--no-igv` | counts and reports only |
| `--igv-path`, `--igv-timeout` | IGV to use (found automatically) and how long to wait for it (default 300 s) |
| `--output` | output folder |
| `--demo`, `--demo-cnv` | run on the bundled synthetic data |

## Demo

```bash
python clawbio.py run igv-validator --demo
```

Expected output: a report on 6 synthetic variants on two synthetic contigs: 3 supported (an SNV, a 12 bp deletion, a translocation), 1 single-strand artifact, 1 with too few reads and 1 present in the normal. With IGV installed, 7 captioned screenshots in about 10 seconds; without it, counts and report only.

`examples/` has a ready samplesheet, gene list and curated table for the demo files (see its README): one
command runs the whole-project workflow on synthetic data and shows a planted wrong curated call as *Look first*.

## Algorithm / Methodology

1. **SNV**: pysam pileup at the position; a read supports the call when its base equals ALT. Overlapping mates are counted once (pileup overlap handling).
2. **Indel**: a molecule supports the call when its CIGAR has an indel of the same length starting within ±5 bp. A molecule soft-clipped (>= 5 bp) exactly at the indel edge also supports it when the clipped bases match the ALT haplotype better than the reference (needs `--reference`).
3. **SV breakpoint**: split reads (`SA` tag within 500 bp of the partner breakpoint) plus discordant pairs (mate within 1 kb of the partner), looked up from **both** breakpoints and merged by read name, so each molecule counts once. Mate BND records collapse to one event; either position selects it.
4. **Flags** on the tumor/normal pair: `normal_support` (any supporting read in normal), `low_support` (< 3 reads), `low_vaf` (< 5%), `strand_bias` (>= 4 reads, all one strand; chance 1 in 2^(n-1)), `read_end` (> 50% of alt bases within 10 bp of a read end), `low_depth` (< 10 reads in tumor or normal), `germline_site` (SNV where >= 20% of normal reads carry a third allele), `high_depth` (depth > 2.5x the run median, with >= 3 variants: repeat or mismapping), `caller_disagrees` (caller-reported VAF differs from the BAM by > 10 points, or the caller reports >= 3 alt reads where the BAM has none). In tumor-only mode `normal_support`, `germline_site` and the swap check are skipped.
   - The caller's column is the VCF sample named like the tumor BAM's `SM`, the only sample, the one sample with "tumor" in its name, or `--caller-sample`; otherwise the comparison is skipped and the report says why.
5. **Status**: `supported` (>= 3 reads, no flags), `flagged` (supporting reads with at least one flag), `insufficient` (< 3 reads), `no_coverage` (no reads in either BAM). Records the skill cannot evaluate are listed separately with the reason.
6. **Copy number** (`--cnv`): read starts per bin (MAPQ >= 20, and all MAPQ) in each gene. **Depth log2** = log2(gene depth / sample-wide depth), the sample-wide depth being the median over 300 random 10 kb windows on all autosomes (so a gain or loss of the gene's own chromosome cannot shift the scale); `cnv_disagrees` when it does not reproduce the segment log2 (off by > 0.5, or deep, below -2, on one side only). The **local ratio** compares the gene with flanks 3x its length (10 kb-1 Mb) placed outside any focal (<= 3 Mb) gain/loss and its adjacent same-call segments, to show focal changes; it is not computed ("flanks too low") when the flanks average under 5x, e.g. next to a telomere or an assembly gap. A gene shorter than 1 kb is measured as one bin of its own length; its depth rests on few reads, so when it suggests a change the review notes ask you to judge it on the depth plot or in IGV. `no_segment` when the gene sits between segments; `ambiguous_mapping` (status flagged for a called gain/loss) when most reads in the gene have MAPQ < 20 with alternative hits (XA) or MAPQ 0, and the all-reads depth log2 (as IGV shows) is reported next to the MAPQ >= 20 one; `low_mapq_remaining` when the reads left in a deep loss are mostly MAPQ < 20 without alternative hits.
7. **Run-level checks**: stop if the reference's contigs or lengths differ from the BAM; warn when the normal carries most of the support (samples probably swapped); warn when variants have no reads (region not in the BAMs).

**Key thresholds / parameters**:
- MAPQ >= 20 and base quality >= 20: common somatic-caller defaults (GATK Mutect2 uses MQ 20 / BQ 10-20)
- Minimum 3 supporting reads, 5% VAF: conventional manual-review floors (Barnell et al. 2019 SOP)
- Indel window ±5 bp; SV windows 500 bp (split) / 1 kb (pairs): validated on the osteosarcoma test set (see Maintenance)

## Example Queries

- "Validate the KRAS and BRAF calls from my Mutect2 VCF in IGV"
- "Use igv-validator on samples S1 and S2 for KRAS, BRAF and PTEN; compare with my curated calls in curated.csv"
- "Does the GATK deletion in PTEN for sample S1 hold up in the reads?"
- "Make interactive IGV pages for these genes so my PI can zoom in"
- "Is the ODF1 G>T call real? Check tumor and normal reads"
- "Take IGV screenshots of these 10 variants for my PI"
- "Count read support for my PASS variants, no screenshots, I'm on the cluster"

## Example Output

```markdown
# IGV Validator Report

**Input**: `demo_calls.vcf` (6 variant(s) checked) · tumor: demo_tumor · normal: demo_normal
**Skill**: igv-validator 0.1.0

**3 supported, 2 flagged, 1 insufficient**.

IGV screenshots: 7 image(s) in `figures/igv/`.

| ID | Gene | Variant | Tumor support | Normal support | Caller reported (VCF) | Flags | Status |
|---|---|---|---|---|---|---|---|
| V1 | DEMO1 | demo1:3,000 C>A | 18/64 (28.1%) | 0/60 (0.0%) | 18/64 (28.1%) | none | **supported** |
| V2 | DEMO2 | demo1:5,000 A>T | 6/53 (11.3%) | 0/55 (0.0%) | 6/53 (11.3%) | strand_bias | **flagged** |
| V3 | DEMO3 | demo1:7,000 A>T | 2/46 (4.3%) | 0/60 (0.0%) | 12/56 (21.4%) | low_support; low_vaf; caller_disagrees | **insufficient** |
| V4 | DEMO4 | demo1:9,000 12 bp deletion | 15/60 (25.0%) | 0/63 (0.0%) | 15/60 (25.0%) | none | **supported** |
| V5_1 | DEMO5 | demo1:12,000 <-> demo2:5,001 breakend | 14 (6 split, 8 pairs) / 46 | 0 (0 split, 0 pairs) / 52 | - | none | **supported** |
| V6 | DEMO6 | demo1:15,000 A>G | 31/63 (49.2%) | 25/50 (50.0%) | 31/63 (49.2%) | normal_support | **flagged** |

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md              # Primary markdown report
├── report.html            # Self-contained HTML report (images embedded)
├── result.json            # Machine-readable results
├── tables/
│   └── support_counts.tsv # Tumor/normal read support per variant
├── figures/
│   ├── igv/               # Captioned tumor-over-normal screenshots (optional, needs IGV)
│   └── cnv/               # Copy-number mode: depth plots and IGV views (optional, --cnv only)
└── reproducibility/
    ├── commands.sh         # Exact command to reproduce
    ├── environment.yml     # Conda/pip environment snapshot
    └── checksums.sha256    # Checksums of the outputs
```

## Dependencies

**Required**:
- `pysam` >= 0.22; BAM/VCF/FASTA access and read counting
- `Pillow` >= 10.1; screenshot captions and render checks

**Optional**:
- IGV desktop >= 2.16 (macOS app, or a Linux `igv.sh` / `igv` command such as an HPC module); screenshots. Tested with 2.16.2 and 2.19.7; IGV <= 2.16 has no `currentGenomePath` command, so the skill watches IGV's log for the genome load instead. `--igv-timeout` (default 300 s) covers slow network file systems. Found automatically, or pass `--igv-path` (a path or a command name). Without it the skill writes counts and the report only.
- `igv-reports` >= 1.17 (`pip install igv-reports`); `--interactive` browser pages built on igv.js. Tested with 1.17.0.
- `matplotlib`; depth plots in copy-number mode.

## Step by step: validating a whole project

The full walkthrough is in **[GUIDE.md](GUIDE.md#step-by-step-validating-a-whole-project)**: install and demo,
choosing genes (curated table, `--genes`, a BED), writing `samples.csv` (with the files the curated table was made
from: `snv_list`, `sv_vcf`, `cnv`), checking the reference against the BAMs, the run script, reruns, and runs without a
curated table. In short:

```bash
python skills/igv-validator/igv_validator.py --samplesheet samples.csv --reference hg38.fa \
  --curated-calls curated.csv --annotation gencode.v44.basic.annotation.gtf.gz [--genes KRAS,BRAF]
# -> igv_reports/<date_time>/summary.html (this run) and igv_reports/index.html (all runs)
```

## Reading the summary: what to trust

The page gives no verdicts: *Look first* (the reads and the call differ) and *Consistent so far* (they fit) only
order the review, and every row is checked on its image before it is reported. Details and the key to the images
(IGV colours, depth-plot lines) are in **[GUIDE.md](GUIDE.md#reading-the-summary-what-to-trust)**.

## Gotchas

- **Counting from the picture**: the model will want to count coloured reads in a screenshot. Do not. IGV downsamples and truncates tall panels; numbers must come from `support_counts.tsv`.
- **Status is not a verdict**: the model will want to call `supported` variants "true" and `flagged` ones "false". Do not. The status is a rule-based summary; a germline SNP at the same site, a paralog or a complex allele needs a reviewer. Label any judgement as the agent's reading of the screenshot.
- **Wrong genome**: IGV with a reference other than the BAMs' build shows every read as mismatched. Always pass the same FASTA the BAMs were aligned to; never fall back to a default genome.
- **IGV left running or reused**: never send commands to an IGV the user already has open (it wipes their session). The skill starts its own IGV on a free port and closes it; do not bypass this.
- **Missing ruler / half-drawn images**: IGV can save before it finishes drawing. The skill checks each image and retakes it; do not remove the pause before `snapshot`.
- **Indels undercounted**: reads ending inside a deletion are soft-clipped by aligners. Without `--reference`, clip rescue is off and support is conservative; say so.
- **Restricted data**: with dbGaP or clinical data, do not copy BAMs elsewhere; use `--no-igv` on the cluster where the data lives.
- **IGV start-up failures keep the log**: if IGV does not come up, the report says why and `reproducibility/igv.log` holds IGV's own log; send that, not a guess.
- **HPC screenshots**: IGV needs a display. On a cluster use a remote desktop session (e.g. Open OnDemand "Desktop", not the VS Code app, which has no display), check `igv --version` is at least 2.16, and run the demo there first. On a plain compute node use `--no-igv`.
- **Caller counts differ**: the model will want to treat a mismatch with the caller's `AD` as a bug. Callers realign, cap depth and filter differently (DRAGEN reported 104 reads at a chr4 repeat where 543 MAPQ>=20 reads map). Small SNV differences are normal; large ones usually mean a repeat, which `high_depth` flags.
- **One-sided breakends**: Manta and others write each translocation twice. Never count only the first record's side; the skill merges both ends.
- **Tumor-only cannot see germline**: without a normal, an inherited heterozygous variant looks exactly like a somatic call (about 50% of reads, both strands). Never call a tumor-only `supported` variant somatic; say that a normal (or a population database) is needed to tell.
- **Caller disagreement is a question, not a verdict**: `caller_disagrees` often means the caller capped depth, realigned, or used different filters (a chr4 repeat: DRAGEN 32/104 vs 33/470 MAPQ>=20 reads). Look at `high_depth` and the screenshot before deciding which is right.
- **SV VCFs have no gene names**: SURVIVOR, Manta and Delly write coordinates only, so `--genes` finds nothing there. Use `--regions` with a BED of the genes: an SV matches if a breakpoint is inside a gene or, for deletions/duplications/inversions, if its span covers the gene (as AnnotSV assigns genes). A run where nothing matches writes a report saying so.
- **Copy-number ratios are relative**: the ratio compares the gene with its own flanks, so tumor purity, ploidy (whole-genome doubling makes a one-copy loss look like 0.75) and an altered flank all shift it. Read `cnv_disagrees` as "look again", and report ploidy with the result.
- **Check how the BAMs were aligned**: `samtools view -H tumor.bam | grep '^@PG'` shows the aligner command. With a GRCh38 reference that includes alt contigs (e.g. the Broad `Homo_sapiens_assembly38.fasta`), BWA is alt-aware only when `<index>.alt` sits next to the index; without it, reads in duplicated regions get MAPQ 0 and MAPQ-filtering callers report false losses there. The fix is upstream (alt-aware alignment or a no-alt reference), not in this skill
- **Losses called from MAPQ-filtered depth can be mapping artifacts**: do not trust a deep deletion just because the caller, the MAPQ >= 20 depth and a downstream curated table all agree on it. Check `ambiguous_mapping` first. In regions duplicated on GRCh38 `_alt` contigs, reads get MAPQ 0 with XA hits on the alternate haplotypes, so every MAPQ-filtering tool sees a "homozygous deletion" while IGV shows full coverage (hollow MAPQ-0 reads).
- **Deep deletions in PDX tumors**: reads left inside a homozygous deletion are often mouse or mismapped reads; `low_mapq_remaining` points at that. Do not count them as evidence against the deletion.
- **Whole-VCF runs**: the model will want to run on an entire VCF. With no `--genes`/`--variants`, only the first `--max-variants` (50) records in file order are checked, which are rarely the interesting ones. Ask the user which genes or variants they care about first.
- **No reads is not weak evidence**: `no_coverage` means the BAM has no reads there (a sliced BAM, another sample, another build). Do not report it as "not supported".
- **Swapped samples**: if the report warns that the normal carries most of the support, stop and check the `--tumor`/`--normal` order before interpreting anything.
- **Known gap, clustered mismatches**: alt reads that also carry several other mismatches (misaligned reads; the osteosarcoma ODF1 artifact) are not flagged yet. Look for this in the screenshot and say so.
- **Hallucinated context**: do not add gene-disease claims or VAF interpretations (clonality, purity) that the counts do not show.

## Safety

- **Local-first**: BAMs, VCFs and screenshots stay on the machine; nothing is uploaded
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Audit trail**: command, environment and output checksums in `reproducibility/`
- **No hallucinated science**: thresholds are fixed in code and listed above

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes.
The agent picks the variants with the user, runs the skill, and explains the counts and flags. The agent must NOT change thresholds, invent counts, or present a status as a clinical call.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- The request mentions IGV, read support, manual review, or "is this variant real"
- The inputs are a VCF plus a tumor/normal BAM pair

**Chaining partners**: this skill connects with:
- `nfcore-sarek-wrapper`: its somatic VCFs and recalibrated BAMs are direct inputs here
- `variant-annotation` / `vcf-annotator`: annotate first, then validate the prioritised variants here
- `clinical-variant-reporter`: pass only calls with `supported` status and reviewer agreement

## Maintenance

- **Review cadence**: monthly, and on every IGV desktop release
- **Staleness signals**: IGV batch-command or CLI option changes (`--port`, `--igvDirectory`); pysam API changes; new caller VCF conventions for SVs
- **Validation record**: on a public osteosarcoma WGS dataset (osteosarc.com): 50/50 alt counts identical to `samtools mpileup`; 12/12 tumor/normal counts identical to an earlier validated pipeline using real PURPLE and Manta VCFs; 8 previously unseen DRAGEN calls agreed with DRAGEN `AD` (7 within 1 read); positive/negative controls as expected. Screenshots validated on the demo (3/3 runs) and on real data with the Broad hg38 reference (T0, T1, T2 and two SV callers). Stress-tested on 10 inputs: Mutect2, Strelka, PURPLE, DRAGEN, Manta and ESVEE VCFs; a different sequencing lab; tumor vs another time point's normal; gene selection through VEP `CSQ`; a 34k-record whole-VCF run; wrong reference, missing `chr` prefix and swapped samples; and a run while the user's own IGV was open (it received no commands).
- **Deprecation**: if the skill no longer serves users, archive it to `skills/_deprecated/` with a note explaining why

## Citations

- [IGV desktop](https://igv.org/doc/desktop/); Robinson et al., *Nat Biotechnol* 2011: viewer and batch commands
- [pysam](https://github.com/pysam-developers/pysam): htslib bindings used for counting
- Barnell et al., *Genet Med* 2019, "Standard operating procedure for somatic variant refinement": manual review criteria
- [VCF 4.3 specification](https://samtools.github.io/hts-specs/VCFv4.3.pdf): BND and symbolic SV notation
