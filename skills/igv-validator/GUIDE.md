# IGV Validator: user guide

The detailed companion to `SKILL.md`: every mode with a prompt and a command, the whole-project workflow step by step, and how to read the summary. Agents: read this file when running the skill for a user.

## Modes

Each mode below shows what it does, what to ask Claude (Claude Code or another agent with ClawBio), and the
command to run yourself. Names are placeholders: replace samples, genes and paths with your own.

**1. Check variant calls, tumor + normal.** Counts tumor and normal reads for each SNV, indel or SV in a VCF,
compares them with the caller's own counts, flags artifacts (normal support, strand bias, read ends, low VAF...)
and takes a tumor-over-normal IGV screenshot per call.
- Prompt: *"Use igv-validator to check the KRAS and BRAF calls in my annotated `calls.vcf` against `tumor.bam`
  and `normal.bam` (reference `hg38.fa`)."*
```bash
python skills/igv-validator/igv_validator.py --vcf calls.vcf --tumor tumor.bam --normal normal.bam \
  --reference hg38.fa --genes KRAS,BRAF --output reports/sampleA/snv
```

**2. Check variant calls, tumor-only.** The same without a normal (leave out `--normal`); normal-based checks
are skipped, so the result says whether the reads support the call, not whether it is somatic.
- Prompt: *"Tumor-only: do the Mutect2 calls in `sampleA.vcf` for the genes in `my_genes.bed` agree with
  `sampleA.bam`?"*
```bash
python skills/igv-validator/igv_validator.py --vcf sampleA.vcf --tumor sampleA.bam --reference hg38.fa \
  --regions my_genes.bed --output reports/sampleA/snv
```
`--genes` needs gene names in the VCF (VEP, SnpEff, ANNOVAR); a plain caller VCF is selected with `--regions`
(or `--variants` with a list from your own filtered table).

**3. Check structural variants (SVs).** SV VCFs (Manta, Delly, SvABA, SURVIVOR) have no gene names, so select by
gene coordinates: an SV is kept when a breakpoint falls in a gene, or a DEL/DUP/INV spans it. Split reads and
discordant pairs are counted at both ends.
- Prompt: *"Check the SURVIVOR SVs in `sampleA.sv.vcf` that hit the genes in `my_genes.bed`."*
```bash
python skills/igv-validator/igv_validator.py --vcf sampleA.sv.vcf --tumor sampleA.bam --reference hg38.fa \
  --regions my_genes.bed --output reports/sampleA/sv
```

**4. Check copy-number calls.** For each gene in the BED, sets the caller's segment (GATK `called.seg` or a table
with contig/start/end/log2/cn_call) against the read depth, per segment when several cover the gene, and reports
depth steps inside a gene, ambiguous mapping (MAPQ 0 reads that match alt contigs), a depth plot and an IGV view.
- Prompt: *"Do the GATK copy-number calls for sampleA match the read depth in PTEN and BRAF?"*
```bash
python skills/igv-validator/igv_validator.py --cnv sampleA.called.seg --cnv-sample sampleA --tumor sampleA.bam \
  --reference hg38.fa --regions my_genes.bed --output reports/sampleA/cnv
```

**5. Summarize many samples.** Reads every run under a folder (`reports/<sample>/{snv,sv,cnv}`) and writes one
page, `reports/summary.html`, that links everything: every check of every sample (a *Reports* line at the top,
and each row's own reports, opening at its gene), images and interactive views. Also `igv_agreement.tsv` and
`summary.tsv`. **The page keeps itself up to date:** give the first check `--project reports`, and from then on
every check inside `reports/` refreshes `summary.html` when it finishes (in a second or two: tables and links;
overview images, interactive pages and the zip are rebuilt only by a full `--summarize`, and stay linked).
- Prompt: *"Summarize all igv-validator runs in `reports/`."*
```bash
python skills/igv-validator/igv_validator.py --summarize reports/
```

**6. Compare with your curated calls.** Your final table after filtering and review (for example the one behind a
mutational-profile heatmap): CSV/TSV with `sample, gene, alteration` (WT, SNV, SV, SNV+SV, DEL, AMP, LOH). The page
then opens on "Curated calls vs IGV": one row per sample and gene with a plain sentence of what IGV shows,
*Start here* (*Look first* / *Consistent so far*), a thumbnail and links; raw calls are folded below.
- Prompt: *"Compare my curated calls in `curated.csv` with IGV for all samples in `reports/`."*
```bash
python skills/igv-validator/igv_validator.py --summarize reports/ --curated-calls curated.csv
```

**7. Images and interactive views.** `--overview` draws one IGV image per gene and sample with every call marked;
`--interactive` writes zoomable pages (igv-reports) per sample and gene, where reads can be clicked;
`--annotation` adds a genes track (exons, names) to every image. The last two need `--regions`.
- Prompt: *"Add overview images and interactive views for the genes in `my_genes.bed`, with GENCODE genes drawn."*
```bash
python skills/igv-validator/igv_validator.py --summarize reports/ --curated-calls curated.csv \
  --overview --interactive --regions my_genes.bed --annotation gencode.v44.basic.annotation.gtf.gz
```

**8. Keep and share the report.** Every summary is saved in the folder you summarized (`summary.html`, the
tables, `overview/`, `interactive/`), and the whole report is also zipped there as `igv_validation_full.zip`,
linked from the top of the page as a download, next to the page's own path on the server. Unzip it anywhere and
open `summary.html`: every page, image and interactive view works offline. It contains read data, so keep it
where the BAMs may be. `--no-bundle` skips the zip. A report folder can be moved or copied and summarized
again: relative BAM and reference paths in it are looked up from the folder itself.
- Prompt: *"Summarize `reports/` and give me the download of the whole report for my PI."*
```bash
python skills/igv-validator/igv_validator.py --summarize reports/ --curated-calls curated.csv
# -> reports/summary.html and reports/igv_validation_full.zip
```

**Batches: a separate report each time.** Give every batch its own folder inside one reports folder, for
example `igv_reports/2026-09-30_panel8/` and `igv_reports/2026-10-20_newsamples/`, each with its own
`summary.html`, reports and zip (`--project igv_reports/<batch>`). Nothing is overwritten between batches, and the
same sample can be checked again in a later batch. Run `--index igv_reports/` once: it writes
`igv_reports/index.html`, one row per batch (date, samples, genes, curated agreement, download) linking to each
summary, and it updates itself whenever a batch's summary is built or refreshed. It is only written where you
asked for it.
- Prompt: *"Run the new samples as a separate batch in `igv_reports/2026-10-20_newsamples` and keep the batch index."*
```bash
B=igv_reports/2026-10-20_newsamples
python skills/igv-validator/igv_validator.py --vcf S3.vcf --tumor S3.bam --reference hg38.fa --regions my_genes.bed \
  --output $B/S3/snv --project $B
python skills/igv-validator/igv_validator.py --index igv_reports/     # once; afterwards it keeps itself current
```

**9. Counts only, no IGV.** For compute nodes without a display: counts, flags and reports, no screenshots.
- Prompt: *"Count read support for the PASS calls in `calls.vcf`, no screenshots, I'm on the cluster."*
```bash
python skills/igv-validator/igv_validator.py --vcf calls.vcf --tumor tumor.bam --pass-only --no-igv \
  --output reports/sampleA/snv
```

**10. Demo.** Synthetic data, no files needed: `--demo` (variants) and `--demo-cnv` (copy number, including the
alt-contig case ALTGENE).
- Prompt: *"Run the igv-validator demo."*
```bash
python skills/igv-validator/igv_validator.py --demo --output /tmp/igv_demo
```

**A whole project, in one command** (every sample, every check, a new dated folder, its summary and the index of
runs): `--samplesheet`, see *Step by step* below. Or ask:
*"Use igv-validator on samples S1 and S2. The BAMs, Mutect2, SURVIVOR and GATK outputs are under
`<results folder>`. Tumor-only. Compare with my curated calls in `<table.csv>`."* (add *"only for KRAS, BRAF and
PTEN"* to narrow the genes). Claude then writes `samples.csv`, checks that every file exists and that the reference
matches the BAMs (`@SQ` lines vs the `.fai`), runs `--samplesheet` with `--annotation` and gives the path of the
run's `summary.html`; it asks for anything missing (file locations, reference, tumor-only or paired).

## Step by step: validating a whole project

A worked example with made-up names: samples `sampleA` and `sampleB`, genes KRAS and BRAF, calls from
Mutect2 (SNVs), SURVIVOR (SVs) and GATK (copy number). Replace paths with your own.

**1. Install and try the demo** (no data needed):
```bash
git clone https://github.com/ClawBio/ClawBio.git && cd ClawBio && pip install -e .
python skills/igv-validator/igv_validator.py --demo
```
IGV desktop (2.16 or later) must be installed, or loaded (`module load igv` on an HPC; run from a desktop
session, since IGV needs a display). Add `--no-igv` to skip screenshots.

**2. Choose your genes.** Gene names are enough: the tool finds each gene in the `--annotation` GTF (GENCODE).
- By default, the genes of your curated calls (`--curated-calls`, step 5).
- `--genes KRAS,BRAF` (or a text file with one name per line) checks just these.
- `--regions my_genes.bed` (chrom, start, end, name; start is 0-based) for your own coordinates.

Names must be official symbols (e.g. `CDKN2A`, not `p16`). An unknown or misspelled name stops the run before
anything runs, with suggestions ("did you mean HOMDEL?"); `--skip-unknown-genes` runs without them instead. The run
folder keeps the coordinates it used in `genes.bed`.

**3. List your samples** in a samplesheet (CSV or TSV), one row per sample. `sample` and `tumor` are
required; add whichever calls you have (leave a cell empty to skip that check):

| sample | tumor | normal | snv_vcf | snv_list | sv_vcf | cnv | cnv_sample |
|---|---|---|---|---|---|---|---|
| sampleA | bams/sampleA.bam | | vcf/sampleA.mutect2.vcf | filtered/sampleA.snvs.tsv | vcf/sampleA.survivor.vcf | cnv/sampleA.called.seg | |
| sampleB | bams/sampleB.bam | bams/sampleB_normal.bam | vcf/sampleB.mutect2.vcf | filtered/sampleB.snvs.tsv | | cnv/sampleB.called.seg | |

When file names follow a pattern, a loop writes the samplesheet (check every path exists before running):
```bash
echo "sample,tumor,normal,snv_vcf,snv_list,sv_vcf,cnv,cnv_sample" > samples.csv
for s in S1 S2; do
  echo "$s,/data/bams/$s.bam,,/data/vcf/$s.mutect2.vcf,/data/filtered/$s.snvs.tsv,/data/vcf/$s.survivor.vcf,/data/cnv/$s.called.seg," >> samples.csv
done
tail -n +2 samples.csv | tr ',' '\n' | grep '^/' | while read f; do [ -e "$f" ] || echo "MISSING: $f"; done
```
`normal` empty = tumor-only. `snv_list` = a table of chrom, pos to check instead of every SNV in the genes.
**When you compare with curated calls (step 5), give `snv_list`: the filtered variants your curated table was made
from**, the same file the curated table was built from (it may cover the whole genome; only variants in the
genes are checked; `chr`/`start` columns work). Each listed variant must also be in `snv_vcf`, which supplies the
caller's counts. Without it, every caller SNV in the genes is checked; in tumor-only data most are inherited, so the run warns,
and those variants are listed under *details* without counting against a curated call. `cnv_sample` = the sample's name inside a multi-sample copy-number file.
**The same goes for `sv_vcf` and `cnv`: use the files your curated table was made from.** A pipeline often writes several SV VCFs (raw, merged, repeat-filtered, size-filtered, final); if the table came from an annotated version (e.g. AnnotSV output), give the VCF that annotation was run on. A different file checks other calls, and a curated SV missing from it shows as *Look first* for the wrong reason. To find it, match a few positions from your table against each candidate VCF.

**4. Run everything with one command**:
```bash
python skills/igv-validator/igv_validator.py --samplesheet samples.csv --reference hg38.fa \
  --curated-calls curated.csv --annotation gencode.v44.basic.annotation.gtf.gz [--genes KRAS,BRAF]
```
To rerun later, save the command as a script next to `samples.csv`, once. This block writes it and shows its end
(replace the paths; on an HPC add your `module load` lines, e.g. `module load igv`, above `python`). The script starts
in its own folder, so anyone who copies the folder can run it:
```bash
cat > run_igv_validation.sh <<'EOF'
#!/bin/bash
# IGV validation: each run = one new dated folder in igv_reports/ (never overwrites)
cd "$(dirname "$0")"
python ClawBio/skills/igv-validator/igv_validator.py \
  --samplesheet samples.csv \
  --reference /path/to/hg38.fa \
  --curated-calls curated.csv \
  --annotation gencode.v44.basic.annotation.gtf.gz \
  "$@"
EOF
tail -3 run_igv_validation.sh      # the last line should be "$@"
```
Before the first run, check that the reference is the one the BAMs were aligned to: the BAM's `@SQ` lines must match
the reference `.fai` (names and lengths).
The last line, `"$@"`, passes extra options on, so a run with fewer genes needs no edit:
```bash
bash run_igv_validation.sh --genes KRAS,BRAF,PTEN
```
- `bash run_igv_validation.sh` checks every gene in the curated table.
- `bash run_igv_validation.sh --genes ...` checks only the genes listed (curated rows for other genes are left out
  of that run).
Each run goes into a **new folder named by date and time**, `igv_reports/2026-10-20_14-32/` (never an existing
one, so nothing is overwritten). The tool runs every check for every sample, then builds that run's
`summary.html` (with overview images, interactive views and the whole report zipped) and updates
`igv_reports/index.html`, which lists every run. It prints both paths at the end. A check that fails (say, a
missing file) is reported and skipped; the others still run (`run_log.tsv` in the run folder), and the run folder
keeps a copy of the samplesheet. `--reports-dir` changes where runs are kept; `--no-igv` skips screenshots.

**5. Compare with your curated calls** (optional, `--curated-calls` above): the table your filtering and review
produced (for example the one behind a mutational-profile heatmap), as a CSV/TSV with at least these columns:

| sample | gene | alteration |
|---|---|---|
| sampleA | KRAS | SNV |
| sampleA | BRAF | AMP |
| sampleB | KRAS | WT |

Labels: WT, SNV, SV, SNV+SV, DEL, AMP, LOH. `sample` must match the samplesheet exactly: a near miss (`sample-A`
for `sampleA`) stops the run with a "did you mean"; samples not in this run are listed on the page, not compared.
`gene` must be the official symbol. The page
then opens on "Curated calls vs IGV". The comparison assumes a common convention for curated tables: copy number is
kept only when focal (< 1 Mb) or deep (|log2| > 2), and SNV/SV take priority over copy number for a gene; changes
left out by that rule are explained, not counted as differences. If your table follows other rules, read the
*details* of each row.

**Next run, other samples or genes: edit and run.**
1. **Samples**: edit `samples.csv` (add or remove rows), or make a new file, e.g. `samples_batch2.csv`, and point
   the script's `--samplesheet` line at it.
2. **Genes**: add them to your curated table, or list them when you run: `bash run_igv_validation.sh --genes A,B`.
3. **Run** `the same command again`. The run gets its own dated folder; earlier runs stay as they are, and
   `igv_reports/index.html` lists all of them.

**Without a curated table**, remove the `--curated-calls` line and give the genes with `--genes`, either when you
run (`bash run_igv_validation.sh --genes KRAS,BRAF,PTEN`) or typed in
the script:
```
  --genes KRAS,BRAF,PTEN
```
or from a text file, e.g. `genes.txt` with one name per line:
```
KRAS
BRAF
PTEN
```
and `--genes genes.txt` in the script. The page then compares each raw caller call with IGV instead of curated
calls. With neither a curated table nor `--genes`, the run stops and asks for genes.

**Which genes are checked**, in order:
1. `--regions genes.bed`, if given (your own coordinates).
2. `--genes`, if given. Curated rows for other genes are left out of that run.
3. Otherwise, the genes in the curated table.

**When the curated table and the samples don't line up**:
- New samples with no curated rows are still checked; they appear only under the raw calls.
- Curated rows for samples that aren't in this run are listed as "not compared".

*Running single checks instead* (one sample, one call type) works as in *Modes* above; `--project <folder>` then
keeps that folder's `summary.html` current.

**6. Judge the images yourself**: for each row, open the thumbnail, the report or the interactive view, and use
the "How to read the images" key on the page. Report a difference only after you have seen it yourself.

## Reading the summary: what to trust

The page gives no verdicts. *Start here* only orders your review: *Look first* where the reads and the call
differ, *Consistent so far* where they fit, which still needs your own look before you report it. It is rule-based
and can be wrong, especially for copy number. The IGV screenshots, depth plots and interactive views are the
evidence; the labels only point you to them. Each row's *details* holds the
technical reason and a review priority: *check first* (a difference or a known risk: ambiguous mapping, several
segments over a gene, a depth step, a value near a threshold, no usable flanks, weak support), *quick look* (copy
number left out, not visible or unclear) and *low priority* (a clean SV, a deep deletion with the reads gone, a
clear match). Whatever the priority, look at the image yourself before presenting a result, and never cite the
verdict alone.

The *report* link of a curated row opens the sample's copy-number report (or its only report), showing only that
gene's table row and plots. Tabs at the top (CNV · SNV · SV) switch to the sample's other reports on the same gene;
*show all genes* brings back the rest. In the interactive views, gene-track labels read GENE_TRANSCRIPT, the gene
name and the transcript drawn.

### How to read the images

- **IGV**:
  - coverage track = amount (the reads below are a sample);
  - grey reads = well mapped; hollow = MAPQ 0 (maps equally well elsewhere);
  - coloured ticks = bases that differ (A green, C blue, G orange, T red; a column down many reads is a real variant);
  - purple I = insertion; black line = deletion;
  - coloured tails lining up = split reads at a breakpoint;
  - whole reads coloured = pair problems (red farther apart / deletion, blue closer / insertion, teal wrong
    orientation / inversion or duplication, other colours = mate on another chromosome / translocation);
  - `calls` track = what the callers reported; `genes` track = exons, introns and direction.
- **Depth plots**:
  - blue = MAPQ >= 20 depth (what GATK counts); grey = all reads (what IGV shows);
  - dashed = the sample's normal level; orange = GATK segment; pink = gene;
  - blue at 0 with grey up = ambiguous reads, not a deletion;
  - a sharp jump inside the gene = a copy-number breakpoint.

The same key is on every report page and the summary ("How to read the images").
