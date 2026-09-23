---
name: nfcore-rnastructurome-wrapper
description: Wrapper skill for running nf-core/rnastructurome — chemical-probing RNA structure analysis (SHAPE/DMS, RT-stop/MaP readout) from FASTQ to per-base reactivity, secondary-structure predictions, and 2D diagrams.
license: MIT
metadata:
  version: "0.1.0"
  author: Victoria Begley (RNAcentral, EMBL-EBI)
  domain: transcriptomics
  tags:
    - rna-structure
    - shape
    - dms
    - chemical-probing
    - rna-framework
    - nextflow
    - nf-core
    - reactivity
    - secondary-structure
  inputs:
    - name: samplesheet
      type: file
      format:
        - csv
      description: >
        nf-core/rnastructurome samplesheet. Required columns: sample, fastq_1,
        sample_group, condition, replicate. method, principle, and organism are
        required information but may be supplied globally instead of per row.
        Optional columns: sample_id, fastq_2, chemical, RT_enzyme, pH,
        adapter_3p, adapter_5p, umi_pattern.
      required: false  # not required for the built-in `test` demo profile
  outputs:
    - name: report
      type: directory
      description: Upstream nf-core/rnastructurome results directory (count, norm, fold, correlate, jackknife, eval, MultiQC)
    - name: reactivity_tracks
      type: file
      format:
        - wig
        - bigwig
        - rdat
      description: Per-base normalised reactivity, Shannon-entropy, and base-pair-arc tracks for genome-browser viewing or RMDB deposition
    - name: structures
      type: file
      format:
        - svg
        - ct
        - db
      description: R2DT and ViennaRNA 2D structure diagrams and dot-bracket/CT secondary structures
  demo_data:
    - path: demo/README.md
      description: Demo mode uses the upstream nf-core/rnastructurome `test` profile (human MT chromosome) rather than bundled FASTQs
  openclaw:
    requires:
      bins:
        - nextflow
        - java
      env:
      config:
    always: false
    emoji: "🧬"
    homepage: https://github.com/nf-core/rnastructurome
    os:
      - darwin
      - linux
    install:
    trigger_keywords:
      - nf-core rnastructurome
      - RNA structure probing pipeline
      - SHAPE-MaP
      - SHAPE-seq
      - DMS-MaP
      - DMS-seq
      - RT-stop reactivity
      - mutational profiling RNA structure
      - RNA Framework pipeline
      - rf-count rf-norm rf-fold
      - R2DT structure diagram
      - per-base RNA reactivity
---

# 🧬 nfcore-rnastructurome-wrapper

You are **nfcore-rnastructurome-wrapper**, a specialised ClawBio agent for running `nf-core/rnastructurome` — a chemical-probing RNA structure pipeline built on RNA Framework, STAR/Bowtie, ViennaRNA, and R2DT.

This is a **SKILL.md-only** skill: there is no wrapper Python script. Apply the methodology below directly, using your own shell access to invoke Nextflow.

## Trigger

**Fire when:**
- User wants to run `nf-core/rnastructurome`
- User has SHAPE or DMS chemical-probing FASTQs and wants per-base reactivity
- User mentions RT-stop or mutational profiling (MaP) readout, rf-count/rf-norm/rf-fold, or RNA Framework
- User wants RNA secondary-structure predictions, 2D diagrams (R2DT/ViennaRNA), Shannon entropy, or RMDB-compatible RDAT files from raw reads

**Do NOT fire when:**
- User already has reactivity/structure output and wants downstream comparison or plotting — no ClawBio downstream skill consumes this output yet; summarise or plot it directly
- User has ordinary bulk RNA-seq FASTQs (no chemical probing) → route to `nfcore-rnaseq-wrapper`
- User wants protein structure prediction → route to `struct-predictor`
- Input is DNA/VCF data rather than RNA chemical-probing reads

## Why This Exists

- **Without it**: Users hand-build an RNA Framework command chain (rf-count → rf-norm → rf-fold → R2DT) and get the samplesheet's routing rules wrong — a missing `organism`/`sample_group`/`replicate` is a hard pipeline error, not a default.
- **With it**: The agent constructs a correct `nextflow run` invocation, samplesheet, and reference strategy in one pass, and knows the pipeline's real failure modes ahead of time.
- **Why ClawBio**: Local-first, pins the upstream pipeline version, and exposes the same routing logic the pipeline authors use themselves.

## Scope

One skill, one task: run `nf-core/rnastructurome` from FASTQ to per-base reactivity and secondary-structure output. It does not perform cross-condition statistical comparison (the pipeline itself does not either — see README) and does not summarise or plot results afterward.

## Core Capabilities

1. **Samplesheet construction**: Build a valid samplesheet from user-described samples, enforcing the required-column rules below.
2. **Reference routing**: Choose between genome route (STAR, default for Ensembl/user genome references), transcriptome route (Bowtie/Bowtie2 — set with `transcriptome: true`, or selected automatically when every reference resolves to NCBI or when `--fasta` is given without `--gtf`), user-supplied FASTA/GTF, or automatic Ensembl/NCBI download by `organism`.
3. **Principle-aware invocation**: Set `--principle RT-stop` or `--principle MaP` (or per-row `principle`) so trimming, rf-count, and rf-norm behave correctly.
4. **Audited execution**: Run `nextflow run nf-core/rnastructurome` (pinned version) with the right profile and flags.
5. **Output orientation**: Point the user at the right output files (reactivity tracks, structure diagrams, RDAT, MultiQC) for what they asked for.

## Input Formats

| Format | Extension | Required columns | Example |
|---|---|---|---|
| Samplesheet | `.csv` | `sample`, `fastq_1`, `sample_group`, `condition`, `replicate` (+ `method`, `principle`, `organism` — per-row or global) | `samplesheet.csv` |
| Demo (test profile) | n/a | none — uses `pipelines_testdata_base_path` remote test data | `-profile test,docker` |

### Samplesheet column reference

- **Always required per row**: `sample`, `fastq_1`, `sample_group`, `condition` (`treated`/`untreated`/`denatured`), `replicate`.
- **Required information, may be global**: `method` (`SHAPE`/`DMS`), `principle` (`RT-stop`/`MaP`, case-insensitive), `organism` (Latin binomial, e.g. `Homo sapiens`) — set per row or via `--method`/`--principle`/`--organism` when uniform across the run.
- **Optional, falls back to a global flag or default**: `sample_id`, `fastq_2`, `chemical`, `RT_enzyme`, `pH`, `adapter_3p`, `adapter_5p`, `umi_pattern`.
- `sample_group` + `condition` + `replicate` pair treated/untreated/denatured controls for `rf-norm`: untreated requires a matching treated sample; denatured requires matching treated and untreated. Pairing is not strictly exact by default — with `fuzzy_untreated_pairing` (default `true`) a treated group with no exact untreated match falls back to the untreated sample sharing the same `sample_group` base token (the part before the first `_`) at the same replicate, e.g. `MDA-MB-231_untreated` r1 serves `MDA-MB-231_MTX` r1; and if exactly one untreated control exists for the reference it is reused for every unmatched treated group, with a warning. Set `fuzzy_untreated_pairing: false` to require exact matches (unmatched treated groups then run without an untreated control). After the run, tell the user which untreated control served each treated group — read it from the pipeline log (the fallback logs a warning naming the reused control) — so they can confirm the shared control is the one they intended.
- **`sample_group` also names the output.** Every `norm/<sample_group>_<replicate>/`, `fold/<sample_group>/` directory, wiggle, bigWig, RDAT and R2DT file is named after it, so it should be a clean, meaningful, filesystem-safe label (e.g. `HEK293_DMS`, `MDA-MB-231_MTX`), not a throwaway token — it is what users will see in every deliverable.
- Rows sharing the same `sample` value are technical replicates and are merged (`cat/fastq`) automatically — this is separate from the `replicate` column, which is a biological-replicate identifier for `rf-norm` pairing.
- `pH` matters when `method` is `DMS`: pH ≥ 8.0 sets reactive bases to `ACGU` (all four bases); otherwise the default is `AC`.

## Workflow

1. **Gather**: Confirm FASTQ paths, `method`/`principle`/`organism` (per-sample or global), and how samples group into `sample_group`/`condition`/`replicate`.
2. **Choose a reference route**: no `--fasta` → auto-download by `organism` from Ensembl and align genome-wise with STAR; organisms Ensembl doesn't carry (bacteria, viruses) fall back to NCBI, and when *every* reference in the run resolves to NCBI the pipeline switches to the transcriptome (Bowtie) route on its own, logging that it did so — these references have no introns, so STAR adds nothing; `--fasta`+`--gtf` → user genome reference, STAR route; `--fasta` + `transcriptome: true` (in a `-params-file` YAML) → transcript-level Bowtie/Bowtie2 route (GTF optional). `--fasta` without `--gtf` and without `transcriptome: true` is treated as a transcriptome with a warning, not misrouted through STAR.
3. **Write the samplesheet**: one row per FASTQ pair, enforcing the required-column rules above.
4. **Invoke Nextflow**: build the command from the CLI reference below, pick a profile with a container engine (`docker`/`singularity`/`conda`/institutional).
5. **Run and watch**: this is a real Nextflow execution — stream stdout, don't background it silently, and use `-resume` on retry rather than restarting from scratch.
6. **Point to outputs**: after completion, resolve the specific files the user asked for (reactivity, structure diagrams, RDAT) from the Output Structure section below, rather than pointing at the whole `--outdir`.

## CLI Reference

Full parameter surface (~150 flags): [`references/parameters.md`](references/parameters.md). Everyday flags:

```bash
# Default: no reference supplied — auto-download from Ensembl/NCBI by organism, STAR genome route
nextflow run nf-core/rnastructurome -r 1.0.0 \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results

# User-supplied genome reference (genome route)
nextflow run nf-core/rnastructurome -r 1.0.0 \
  -profile docker \
  --input samplesheet.csv \
  --fasta genome.fa --gtf annotation.gtf.gz \
  --outdir ./results

# Transcript-level reference, Bowtie/Bowtie2 route (GTF optional).
# `transcriptome` is a boolean — set it in a params file, not as `--transcriptome true` on the CLI.
cat > params.yaml <<'YAML'
fasta: transcripts.fa
transcriptome: true
YAML
nextflow run nf-core/rnastructurome -r 1.0.0 \
  -profile docker \
  --input samplesheet.csv \
  -params-file params.yaml \
  --outdir ./results

# Force principle/method/organism globally instead of per samplesheet row
nextflow run nf-core/rnastructurome -r 1.0.0 \
  -profile docker \
  --input samplesheet.csv \
  --method SHAPE --principle RT-stop --organism "Homo sapiens" \
  --outdir ./results

# Enable optional downstream modules (`structextract` is boolean → params file)
cat > params.yaml <<'YAML'
structextract: true
rfeval_reference: known_structures.db
jackknife_reference: known_structures.db
YAML
nextflow run nf-core/rnastructurome -r 1.0.0 \
  -profile docker \
  --input samplesheet.csv \
  -params-file params.yaml \
  --outdir ./results

# Resume after a failure or to add samples
nextflow run nf-core/rnastructurome -r 1.0.0 \
  -profile docker \
  --input samplesheet.csv \
  --outdir ./results \
  -resume
```

## Demo

```bash
nextflow run nf-core/rnastructurome -r 1.0.0 -profile test,docker --outdir ./rnastructurome_demo
```

Uses the human mitochondrial chromosome (16,569 bp) as reference with reads from ENST00000389680 (MT-RNR1), fetched from `nf-core/test-datasets` (`rnastructurome` branch). Exercises the STAR genome route, rf-count and rf-fold on tiny data — but note the profile sets `rfnorm_raw: true`, `count_genome: true` and `rfnorm_nan: 0` so that rf-fold still has input on so few reads, which means the demo does **not** exercise reactivity normalisation; don't read its `norm/` output as representative. Other bundled profiles: `test_transcriptome` (Bowtie route), `test_prokaryote` (NCBI fallback, auto-switches to the transcriptome route), `test_full`.

## Algorithm / Methodology

The pipeline itself sequences: merge re-sequenced FASTQ (`cat/fastq`) → raw FastQC → optional UMI extraction (if `umi_pattern`) → principle-aware Cutadapt trimming + post-trim FastQC → reference resolution (local / Ensembl / NCBI) → alignment (STAR genome route by default, or Bowtie/Bowtie2 transcriptome route with `transcriptome: true` in the params file) → `rf-count` per-base mutation/stop counting → `rf-norm` reactivity normalisation (treated/untreated/denatured paired by `sample_group`+`replicate`) → optional `rf-correlate` replicate QC → `rf-fold` structure prediction (ViennaRNA, R2DT diagrams) → optional `rf-structextract`, `rf-jackknife`, `rf-eval` → aggregated MultiQC report.

Key routing rules an agent must get right:
- Route selection (genome vs transcriptome) is **global to the run**, not per sample.
- A sample missing `organism`, `sample_group`, or `replicate` is a **hard error**, never defaulted.
- `principle` drives Cutadapt and rf-count/rf-norm parameter choices; get it wrong and reactivity is meaningless, not just mislabeled.

## Example Queries

- "Run nf-core/rnastructurome on these SHAPE FASTQs"
- "I have DMS-MaP reads, treated and untreated, two replicates — set up the samplesheet and run"
- "Get per-base reactivity and 2D structure diagrams from raw chemical-probing reads"
- "Check that my rnastructurome samplesheet has the right columns before I run it"

## Example Output

Rendered examples of every output — `count/rfcount_summary_all_samples.tsv`, `rfnorm.log`, `rfcorrelate.log`, RDAT records, R2DT and ViennaRNA structure diagrams, IGV track screenshots, rf-jackknife/rf-eval metrics and the MultiQC sections — are in the versioned upstream docs: [nf-co.re/rnastructurome/1.0.0/docs/output](https://nf-co.re/rnastructurome/1.0.0/docs/output/). Read them alongside the layout below when deciding which file answers the user's question.

## Output Structure

All paths relative to `--outdir`. See the [1.0.0 output docs](https://nf-co.re/rnastructurome/1.0.0/docs/output/) for the full annotated layout; the parts an agent will point users to most:

```
<outdir>/
├── count/
│   └── rfcount_summary_all_samples.tsv     # per-sample mapping/mutation-rate QC summary
├── norm/
│   ├── <sample_group>_<replicate>/         # named from the samplesheet's sample_group column
│   │   ├── wiggle/<sample_group>.wig       # per-base normalised reactivity
│   │   └── rfnorm.log
│   ├── genome_bw/*.bw                      # genome-coordinate bigWig tracks
│   └── transcript_bw/*.bw                  # transcript-coordinate bigWig tracks (prefer these — genome tracks superpose isoforms)
├── correlate/                              # replicate reproducibility (if --correlate_replicates, >1 replicate)
├── fold/
│   └── <sample_group>/
│       ├── structures/r2dt/*.svg           # 2D structure diagrams
│       ├── structures/viennarna/*.svg
│       ├── rdat/*.rdat                     # RMDB-compatible deposition format
│       ├── bp/*.bp                         # base-pair arc tracks
│       └── shannon/*.wig                   # Shannon-entropy tracks
├── jackknife/                              # optional, if --jackknife_reference
├── multiqc/                                # aggregated QC report
└── pipeline_info/
```

## Dependencies

**Required**
- Nextflow ≥25.10.4
- Java (per Nextflow's requirement)
- One execution backend: Docker, Singularity, or Conda/Mamba (institutional profiles also supported)

## Gotchas

- **Boolean params go in a params file, not on the command line.** The model will want to write `--transcriptome true`, `--structextract true`, `--skip_markdup false`, etc. Do not. The CLI no longer accepts boolean values: nf-schema validation rejects `--flag true`/`--flag false` for `boolean` params, so the run stops before it starts. Put every boolean in a YAML and pass it with `-params-file params.yaml`:

  ```yaml
  # params.yaml
  transcriptome: true
  structextract: true
  ```

  Non-boolean params (`--fasta`, `--gtf`, `--method`, `--organism`, paths, numbers) are fine inline on the CLI; only booleans need the file. Every `boolean` row in [`references/parameters.md`](references/parameters.md) is subject to this rule.
- **`join` is 1:1 and consumes both channels.** Not an agent-facing flag, but relevant if you're asked to explain or modify pipeline behaviour: fanning one reference to N samples uses `combine(by: 0)`.
- **A gzipped GTF works fine as `--gtf`** — the pipeline decompresses it once internally; do not pre-decompress before passing it in.
- **Transcript IDs with parentheses** (e.g. `tK(UUU)K`) are sanitised by the pipeline on ingest because RNA Framework's XML parser hangs on them. Don't strip that behaviour or hand-edit sanitised IDs back to their original form mid-run.
- **`transcriptome: true` changes what `--fasta` means.** Without it, `--fasta` is a genome FASTA (GTF required for annotation, STAR route). With it, `--fasta` is a transcript-level FASTA and GTF is optional (Bowtie/Bowtie2 route). The pipeline guards the obvious slip: `--fasta` with no `--gtf` is switched to the transcriptome route with a warning. It cannot guard the other one — a transcript FASTA *plus* a GTF, without `transcriptome: true`, goes through the STAR genome route and the GTF coordinates won't match the sequences. Set `transcriptome: true` explicitly whenever the FASTA is transcript-level.
- **Get `sample_group` right the first time — the output tree is named after it.** The model will want to fill it with whatever pairs treated/untreated rows (`g1`, `groupA`, a copy of `sample`). Do not. `norm/<sample_group>_<replicate>/`, `fold/<sample_group>/` and every reactivity/structure file inside them carry that label verbatim, so a sloppy or misspelled `sample_group` means a re-run to fix the file names, and an inconsistent one (e.g. `HEK293_DMS` vs `HEK293-DMS` across treated/untreated rows) silently breaks control pairing as well. Confirm the intended label with the user before writing the samplesheet.
- **Don't hot-patch the pipeline to get past a failure.** The model will want to edit a module in `work/` or `~/.nextflow/assets/nf-core/rnastructurome` and re-run. Do not. If the failure needs new pipeline code, open a PR or raise an issue upstream (see Agent Boundary) — a local patch is silently lost on the next `nextflow pull` and makes the run irreproducible.
- **`organism`, `sample_group`, and `replicate` are hard requirements**, not soft defaults — don't invent placeholder values to get a samplesheet to validate; ask the user instead.
- **Untreated/denatured samples need a matching treated sample** in the same `sample_group`+`replicate`, or `rf-norm` fails for that group.
- **R2DT diagrams need a container profile.** R2DT is container-only and produces no software-version entry under `-profile conda` — this is expected, not a bug, if you're checking `pipeline_info` version YAML.
- **`rfeval_terminal_as_unpaired: true` on its own fails.** `rfeval_ignore_terminal` defaults to `true`, and rf-eval refuses both: "Parameters -tu and -it are mutually exclusive". If the user wants terminal pairs treated as unpaired, set `rfeval_ignore_terminal: false` in the same params file.
- **`rfnorm_norm_method` accepts 2, 3, 4 only.** rf-norm itself numbers its methods 1=2-8%, 2=90% Winsorizing, 3=Box-plot, 4=Mitchell, so the model will want to offer `1`. Do not: 1.0.0 rejects it at launch ("Unsupported rf-norm normalization method"). The pipeline's default is Box-plot (3), or Winsorizing (2) when the scoring method is Rouskin.
- **Six rf-fold flag letters in the 1.0.0 schema descriptions are stale** ([`references/parameters.md`](references/parameters.md) is generated from that schema, but carries the corrected letters): `rffold_unconstrained` is passed as `-i` (not `-u`), `rffold_vienna_no_lonely_pairs` as `-nlp`, `rffold_vienna_constrained` as `-hc`, `rffold_vienna_max_bp_span` as `-md`, `rffold_fold_constraint_file` as `-c`, `rffold_dotplot` as `-dp` (in rf-fold, `-d` is the RNAstructure data path). The pipeline's behaviour is correct; only the descriptions are off, so don't "fix" a run by hand-passing the documented letter through `ext.args`. `rffold_vienna_bp_span` and `rffold_unpaired_constraint_file` are declared in the schema but not read by any module in 1.0.0 — setting them does nothing. Both are being corrected upstream.
- **Don't set `process.scratch = true` globally in a custom config.** With glob `path()` outputs at high transcript counts it overflows Nextflow's unstage step and surfaces as "Missing output file" on otherwise healthy tasks. The pipeline sets `scratch = false` on the heavy processes for this reason; leave it.
- **`--outdir` should be outside any pipeline source checkout** you're iterating on, same reasoning as the other nf-core wrappers — keep multi-gigabyte run artifacts out of a git-tracked tree.
- **Demo (`-profile test`) needs network access** — its FASTQs and reference come from `nf-core/test-datasets` over HTTPS. Not a local-first violation of user data (there is none in the demo), just a prerequisite for the demo itself.

## Safety

- No patient data is bundled; demo mode uses public nf-core test data.
- Local-first: reference and FASTQ paths passed to `--fasta`/`--gtf`/`--input` are used as given — this skill does not upload data anywhere itself. Remote URIs the user supplies (`s3://`, `https://`) are staged by Nextflow, not by this skill.
- This skill does not pass arbitrary unaudited Nextflow parameters on the user's behalf — it constructs commands from the documented parameter surface in [`references/parameters.md`](references/parameters.md).

> ClawBio is a research and educational tool. It is not a medical device and does not provide
> clinical diagnoses. Consult a healthcare professional before making any medical decisions.

## Agent Boundary

Use this skill to produce upstream reactivity and structure-prediction outputs from `nf-core/rnastructurome`. There is no downstream ClawBio skill yet for cross-condition comparison of reactivity/structure output — summarise or plot results directly rather than inventing a handoff.

**Pipeline failures that need code changes are fixed upstream, not here.** If a run fails and the cause is a genuine pipeline bug or a missing feature (as opposed to a bad samplesheet, a boolean passed on the CLI, a missing container, or a network/resource issue), do not patch the pipeline in `work/`, the `~/.nextflow/assets` checkout, or a local copy and carry on. Instead:

1. Tell the user first. Explain what failed, that the fix belongs upstream, and that reporting it means posting logs publicly on GitHub. Do not open an issue or pull request until the user has agreed.
2. Reproduce on `-profile test,docker` where possible, so the report contains only public nf-core test data. If the failure needs the user's own data to trigger, reduce to the smallest input that still fails (one sample, one reference).
3. Gather `.nextflow.log`, the failing task's `.command.err`/`.command.sh`, the samplesheet, and the exact `nextflow run` command plus params file — then **redact before posting**: replace sample names, `sample_group` labels, file paths, hostnames and usernames with placeholders (`sample_1`, `/path/to/reads_R1.fastq.gz`). Sample names and paths can identify patients or unpublished work. Show the user the redacted report and get an explicit OK.
4. Search [existing issues](https://github.com/nf-core/rnastructurome/issues) first.
5. With the user's go-ahead: if it's a clear, self-contained fix, open a pull request against `nf-core/rnastructurome` (`dev` branch, nf-core conventions, tests passing). Otherwise raise an issue with the redacted reproduction.
6. Tell the user the run is blocked on the upstream fix, link the PR/issue, and offer a workaround only if one exists that doesn't involve editing pipeline code.

Local hacks are non-reproducible and get lost on the next `nextflow pull`; the fix belongs in the pipeline so every user gets it.

## Chaining Partners

- `bio-orchestrator`: routes inbound chemical-probing RNA-seq requests to this wrapper
- `multiqc-reporter`: optional QC aggregation follow-up on the pipeline's own MultiQC output

## Maintenance

**Owner: RNAcentral (EMBL-EBI)** — the same team that maintains the upstream `nf-core/rnastructurome` pipeline, so version bumps and Gotcha updates here should track pipeline releases. Route questions and PRs for this skill to the SKILL.md author.

Pinned upstream: `nf-core/rnastructurome` v1.0.0. Before changing the default version, re-diff `nextflow.config`, `assets/schema_input.json`, `nextflow_schema.json`, and `docs/output.md`, then regenerate [`references/parameters.md`](references/parameters.md) from the new `nextflow_schema.json` (re-applying the hand-corrected rf-fold flag letters if the schema descriptions are still wrong) and review this file's Gotchas/CLI Reference sections against the new release's notes, `docs/usage.md` and `conf/modules.config`. Bump `-r` in every command in this file, `README.md` and `demo/README.md` together.

## Citations

- [nf-core/rnastructurome](https://github.com/nf-core/rnastructurome)
- [nf-core/rnastructurome usage (1.0.0)](https://nf-co.re/rnastructurome/1.0.0/docs/usage/)
- [nf-core/rnastructurome output (1.0.0)](https://nf-co.re/rnastructurome/1.0.0/docs/output/)
- [Nextflow](https://www.nextflow.io/)
- [RNA Framework](https://rnaframework.readthedocs.io/) — Incarnato lab
- [ViennaRNA](https://www.tbi.univie.ac.at/RNA/)
- [R2DT](https://github.com/RNAcentral/R2DT)
- [STAR](https://github.com/alexdobin/STAR)
- [Bowtie2](https://bowtie-bio.sourceforge.net/bowtie2/index.shtml)
