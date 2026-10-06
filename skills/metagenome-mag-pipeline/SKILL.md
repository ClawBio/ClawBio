---
name: metagenome-mag-pipeline
description: >-
  Assemble and bin a shotgun metagenome and assess the recovered MAGs, by
  driving the pinned metagenomics-workflow Snakemake pipeline (metaSPAdes,
  MetaBAT2/MaxBin2/VAMB + DAS Tool, CheckM2, dRep, GTDB-Tk).
license: MIT
metadata:
  version: "0.1.0"
  author: Md Mobashir Rahman
  domain: metagenomics
  tags:
    - metagenomics
    - assembly
    - binning
    - MAG
    - snakemake
  inputs:
    - name: input_path
      type: file
      format:
        - tsv
        - fastq.gz
      description: Samplesheet TSV, or a directory of FASTQ files to discover samples from
      required: false
    - name: pipeline_dir
      type: directory
      format:
        - directory
      description: Existing checkout of mobashirrahman/metagenomics-workflow
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Analysis report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable results
  dependencies:
    python: ">=3.11"
    packages:
      - pyyaml>=6.0
  demo_data:
    - path: demo/README.md
      description: How the synthetic-read demo is generated at run time
  endpoints:
    cli: python skills/metagenome-mag-pipeline/metagenome_mag_pipeline.py --input {input_path} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
        - snakemake
        - conda
      always: false
    emoji: "🦠"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - linux
    install:
      - kind: conda
        package: "snakemake=9.11.2"
        channels:
          - conda-forge
          - bioconda
    trigger_keywords:
      - metagenome assembly
      - MAG
      - metagenome-assembled genomes
      - metagenomic binning
      - metaSPAdes
      - MetaPhlAn
      - CheckM2
      - GTDB-Tk
      - dRep
      - DAS Tool
      - contig taxonomy
      - strain heterogeneity
      - host read removal
      - shotgun metagenomics assembly
---

# 🦠 Metagenome MAG Pipeline

You are **Metagenome MAG Pipeline**, a ClawBio agent that assembles and bins a
shotgun metagenome and assesses the recovered MAGs. Your job is to drive the
upstream Snakemake workflow correctly and report exactly what it produced.

This skill wraps, and does not reimplement, the workflow at
<https://github.com/mobashirrahman/metagenomics-workflow> pinned to commit
`7351a702d29857801963cae7aeff6600d59ebe77` (tag `v0.2.0`, MIT). BBTools QC and host
subtraction, MetaPhlAn profiles, metaSPAdes assembly, MetaBAT2/MaxBin2/VAMB
consensus binning through DAS Tool, CheckM2/dRep/GTDB-Tk genome assessment,
contig taxonomy, functional annotation, strain heterogeneity and a read-loss
diagnostic.

## Trigger

**Fire this skill when the user says any of:**
- "assemble this metagenome" / "metagenome assembly" / "shotgun metagenomics assembly"
- "bin my metagenome" / "metagenomic binning" / "recover MAGs"
- "metagenome-assembled genomes" / "MAG quality" / "CheckM2" / "GTDB-Tk"
- "MetaPhlAn profiles" / "profile my shotgun data"
- "metaSPAdes" / "DAS Tool" / "consensus bins" / "dRep dereplication"
- "contig-level taxonomy" / "strain heterogeneity" / "host read removal from metagenome"

**Do NOT fire when:**
- Read-based profiling or resistome work — Kraken2 taxonomy, RGI/CARD
  resistance genes, HUMAnN3 pathways. That is `claw-metagenomics` (alias
  `metagenomics`), which this skill must not take.
- Amplicon data (16S/18S) → `claw-amplicon-qc`.
- A single isolate genome, or long-read (Nanopore/PacBio) assembly.
- The user only wants an existing QC report aggregated → `multiqc-reporter`.

**Design notes:** every trigger names assembly, binning, MAG recovery or a tool unique to them. "Metagenomics" alone is too broad to claim: it belongs to `claw-metagenomics`.

## Why This Exists

- **Without it**: users hand-edit a dozen upstream YAML keys, then discover at hour six that binning needs assembly, that metaSPAdes refuses single-end libraries, or that GTDB-Tk needs a database nobody downloaded.
- **With it**: one validated samplesheet, one rendered config, a plan to inspect before committing hours of compute, and a report that separates what was measured from what was not.
- **Why ClawBio**: presets, resource caps and database requirements are read from the pinned checkout, not recalled. There are no thresholds in this code — only the ones in the effective config.

## Core Capabilities

1. **Samplesheet validation**: upstream's identifier, lane-count, layout and `sra_run` rules enforced before anything launches, every FASTQ path resolved to an absolute path.
2. **Preset-driven config rendering**: `qc`, `profile`, `assembly` or `full` stage switches on the upstream `config.yaml`, with thread and memory caps and absolute database paths.
3. **Preflight**: pinned source resolution, snakemake ≥ 9 and conda presence, output-directory safety, and an explicit list of what will touch the network.
4. **Execution**: a `--dry-run` plan (`--check`) or a real run with a process-group-aware timeout.
5. **Reporting**: the pipeline's `final/` tables parsed into `report.md` and `result.json`, keeping upstream column names and marking unmeasured cells as unmeasured.

## Scope

**One skill, one task.** This skill runs the upstream shotgun assembly and binning workflow and reports its results. It does not profile reads for resistome or taxonomy (that is `claw-metagenomics`) and invents no analysis of its own. A request needing both routes to both skills.

## Input Formats

| Format | Extension | Required fields | Example |
|--------|-----------|-----------------|---------|
| Samplesheet | `.tsv` | `sample`, `fastq_1`; `fastq_2` for paired | `examples/samplesheet.example.tsv` |
| Reads directory | any | FASTQ pairs named `<id>_R1`/`<id>_R2` or `<id>_1`/`<id>_2`, flat or one folder per sample | `data/reads/S1/…` |
| Host FASTA | `.fa`/`.fasta` | uncompressed | `--host-fasta host.fa` |
| Databases | `.dmnd`, directory | CheckM2 1.1.0 `.dmnd`; extracted GTDB-Tk directory; DIAMOND-indexed protein DBs; unpacked taxdump | `--checkm2-db` |

Multiple lanes of one sample are semicolon-separated in matching order and remain one sample. `layout` is `paired` or `single`, and is required for a row naming a public SRA/ENA/DRA run instead of paths.

## Workflow

1. **Check the output directory**: refuse a non-empty one unless `--force`, and confirm it is writable. Then create `<output>/pipeline/`. *(prescriptive)*
2. **Resolve the pipeline source**: use `--pipeline-dir` when given, else a cached clone at `~/.cache/clawbio/pipelines/metagenomics-workflow-<ref[:12]>`, else clone the upstream repo and check out the pinned commit. Verify `PIPELINE_REQUIRED_FILES` first. *(prescriptive)*
3. **Build the samplesheet**: in `--demo`, run the checkout's own `test/make_test_dataset.py` and read the samplesheet it writes; otherwise load `--input` as a file or discover samples from a directory. Write `<output>/pipeline/samplesheet.tsv` with absolute paths. *(prescriptive)*
4. **Render the config**: load the checkout's `workflow/config/config.yaml` (`config_test.yaml` in demo mode), apply the preset switches, the database paths, `qc.decontaminate.*`, `fetch.sra.enabled` and the `--cores`/`--mem-mb` caps, then write `<output>/pipeline/config.yaml`. *(prescriptive)*
5. **Validate the config** against upstream's own guards — binning needs assembly, MAG needs binning, metaSPAdes needs paired-end reads, isolate mode is incompatible, `taxonomy.sensitivity` must be a DIAMOND setting, every enabled stage must have its databases — and confirm each configured database path exists. *(prescriptive)*
6. **Preflight the host**: snakemake ≥ 9.0 and conda or mamba on PATH. Print what will use the network — conda environments, `sra_run` downloads, the GRCh38 host FASTA, the MetaPhlAn index, the MinPath script — to stderr before launching. *(prescriptive)*
7. **Launch**: the snakemake command upstream's own `run_test.sh` uses, plus `--dry-run` under `--check`. Stream combined output to `<output>/pipeline/snakemake.log`; on timeout terminate the whole process group. *(prescriptive)*
8. **Parse and report**: unless `--check`, read the `final/` tables into findings, then write `result.json`, `report.md` and the `reproducibility/` bundle. *(prescriptive)*

**Freedom level guidance:** all eight steps are prescriptive. Upstream rejects an invalid config at parse time rather than hours into a cohort, and every number in step 8 must come from a table the pipeline wrote.

## CLI Reference

```bash
# Standard usage
python skills/metagenome-mag-pipeline/metagenome_mag_pipeline.py \
  --input <samplesheet.tsv|reads_dir> --output <output_dir>

# Discover samples from a directory of FASTQ files
python skills/metagenome-mag-pipeline/metagenome_mag_pipeline.py \
  --input data/reads --output out

# Inspect the plan without executing anything
python skills/metagenome-mag-pipeline/metagenome_mag_pipeline.py \
  --input data/reads --output out --check

# Binning, MAG assessment, taxonomy and annotation (needs every database)
python skills/metagenome-mag-pipeline/metagenome_mag_pipeline.py \
  --input data/reads --output out --preset full \
  --checkm2-db db/checkm2.dmnd --gtdbtk-db db/gtdbtk \
  --taxonomy-db db/nr.dmnd --taxdump db/taxdump \
  --kegg-db db/kegg.dmnd --cog-db db/cog.dmnd --pfam-db db/Pfam-A.hmm

# Via the ClawBio runner
python clawbio.py run mag-pipeline --input <file> --output <dir>
python clawbio.py run mag-pipeline --demo
python clawbio.py run mag-pipeline --demo --check --output /tmp/plan
```

| Flag | Meaning |
|------|---------|
| `--input PATH` | Samplesheet TSV, or a directory of FASTQ files |
| `--output DIR` | Output directory (required) |
| `--demo` | Synthetic reads, upstream's own test configuration; refuses flags that configure a stage it does not enable |
| `--check` | Validate and dry-run only; nothing is executed |
| `--preset {qc,profile,assembly,full}` | Stage switches (default `assembly`) |
| `--cores INT` | Cores for snakemake; caps every `threads.*` value (default 4) |
| `--mem-mb INT` | Total memory in MiB; caps every `mem_mb.*` value (default 8000) |
| `--coassemble` | Pool cleaned reads by assembly `group` |
| `--binners LIST` | `metabat,maxbin,vamb`; at least two are required |
| `--host-fasta PATH` | Decompressed host FASTA instead of the GRCh38 download |
| `--skip-host-removal` | Disable host read subtraction |
| `--checkm2-db PATH` | CheckM2 1.1.0 DIAMOND `.dmnd` file |
| `--gtdbtk-db PATH` | Extracted GTDB-Tk reference directory |
| `--taxonomy-db PATH` | DIAMOND protein DB whose subjects end in an NCBI taxid |
| `--taxdump PATH` | Unpacked NCBI taxdump directory |
| `--kegg-db PATH` | DIAMOND KEGG-orthology database |
| `--cog-db PATH` | DIAMOND COG database |
| `--pfam-db PATH` | Pfam-A profile database |
| `--pipeline-dir PATH` | Use an existing checkout instead of cloning |
| `--conda-prefix PATH` | Where tool environments are built (default `<cache>/conda-envs`) |
| `--allow-remote-inputs` | Permit `sra_run` rows that download public runs |
| `--timeout-hours FLOAT` | Wall-clock cap (default 12) |
| `--force` | Write into a non-empty output directory |

## Demo

```bash
python clawbio.py run mag-pipeline --demo --check --output /tmp/plan   # plan only
python clawbio.py run mag-pipeline --demo --output /tmp/demo           # real run
```

Expected output: the demo generates synthetic reads with the pinned checkout's
`test/make_test_dataset.py`, then runs QC → assembly → QUAST → MultiQC →
provenance, producing a report with per-sample retained reads, QUAST assembly
metrics and the read-loss diagnostic. It needs snakemake ≥ 9, conda and network
access to build tool environments; with `--check` it needs only the checkout.

The synthetic reads have no biological truth. The demo proves execution, and the
report says so.

A demo run is upstream's `workflow/config/config_test.yaml` with three paths
moved, so `--demo` together with `--input`, `--preset` other than `assembly`,
`--coassemble`, `--binners`, `--host-fasta`, `--skip-host-removal`,
`--allow-remote-inputs` or any database flag is refused with `INVALID_CONFIG`
naming each one. Drop `--demo` to apply those settings to your own reads.
`--cores`, `--mem-mb`, `--timeout-hours`, `--pipeline-dir` and `--conda-prefix`
do apply.

## Algorithm / Methodology

1. **Samplesheet**: upstream's `read_samplesheet()` rules, enforced verbatim —
   identifiers `[A-Za-z0-9_.-]+`, unique samples, equal lane counts on both
   mates, no repeated path, `layout` consistent with the mates present, and
   `sra_run` mutually exclusive with paths and requiring `layout`.
2. **Config**: the pinned `config.yaml`, stages switched per preset, every path
   absolute because upstream's `resolve()` joins relative config paths to the
   pipeline root rather than the working directory.
3. **Launch**: `snakemake -s <pipeline>/workflow/Snakefile --configfile … --use-conda
   --conda-frontend conda --conda-prefix … --cores N --resources mem_mb=M
   --rerun-incomplete --printshellcmds`.
4. **Parse**: every `final/`, `reports/` and `provenance/` table the enabled stages
   produced, keeping upstream column names, converting numeric strings to numbers
   and empty cells to "not measured".

**Key thresholds / parameters** — all read from the effective config:
- CheckM2 gate: `mag.checkm2.min_completeness` and `mag.checkm2.max_contamination`
- dRep gate: `mag.drep.completeness`, `mag.drep.contamination`, `ani_threshold`,
  `cov_threshold`
- Contig taxonomy: `taxonomy.consensus_minimum`, `best_hit_bitscore_fraction`,
  `rank`, `sensitivity`
- Heterogeneity: `heterogeneity.min_depth`, `min_base_quality`,
  `max_dominant_fraction`
- Assembly filter: `assembly.min_contig_length`

## Example Queries

- "assemble my gut metagenome and recover MAGs from it"
- "bin these FASTQs with MetaBAT2 and MaxBin and check the genome quality"
- "which bins pass CheckM2 and what does GTDB-Tk call them"
- "did my reads get assembled and then lost at binning?" (the bottleneck diagnostic)
- "profile these shotgun reads with MetaPhlAn" (use `--preset profile`)

## Example Output

> Illustrative values from a synthetic demo run, shown to give the report its
> shape. The numbers are made up; a real run reads every one of them from the
> pipeline's own tables.

```markdown
# Shotgun metagenomics: MAG pipeline report

**Mode**: demo (synthetic reads) · **Preset**: assembly · **Status**: ok
**Pipeline**: metagenomics-workflow @ 7351a702 (tag `v0.2.0`, local checkout)
**Samples**: 2 · **Stages run**: QC, assembly, QUAST, reporting
**Stages skipped**: profiling, binning, MAG assessment, contig taxonomy,
functional annotation, strain heterogeneity — preset `assembly` does not enable them.

## Read QC

| sample | layout | retained_reads | retained_bases | retained_pairs |
|--------|--------|----------------|----------------|----------------|
| TEST_001 | paired | 3800 | 1140000 | 1900 |
| TEST_002 | paired | 3700 | 1110000 | 1850 |

## Assembly (QUAST) and the read-loss diagnostic

| assembly metric | TEST_001 | TEST_002 |
|-----------------|----------|----------|
| # contigs | 331 | 318 |
| Total length | 1284000 | 1241000 |
| GC (%) | 51.20 | 50.80 |

| sample | reads | assembly_capture_pct | catalogue_capture_pct |
|--------|-------|----------------------|----------------------|
| TEST_001 | 3800 | 77.50 | not measured |
| TEST_002 | 3700 | 81.20 | not measured |

`not measured` means the pipeline could not compute that axis. It is not zero and
must not be plotted as zero.

## Caveats

- The demo reads are synthetic: no taxonomy, no marker genes. The numbers above
  prove execution, not accuracy.
- QUAST counts are for the assembly as filtered by `assembly.min_contig_length`.

*ClawBio is a research and educational tool. It is not a medical device and does
not provide clinical diagnoses. Consult a healthcare professional before making
any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md              # Primary markdown report
├── result.json            # Machine-readable results
├── pipeline/
│   ├── samplesheet.tsv    # Validated samplesheet, absolute paths
│   ├── config.yaml        # Effective config handed to snakemake
│   └── snakemake.log      # Combined snakemake stdout/stderr
└── reproducibility/
    ├── commands.sh        # Exact command to reproduce the run
    ├── environment.yml    # snakemake + pandas + pyyaml environment
    └── checksums.sha256   # sha256 of the report, config and final tables
```

`pipeline/workdir/` is snakemake's `--directory`, not part of the contract. The
pipeline writes the rest, and only when the matching stage ran:

| Path | Stage | Notes |
|------|-------|-------|
| `results/final/qc/read_counts.tsv` | always | `sample, layout, retained_reads, retained_bases, retained_pairs` |
| `results/final/qc/sourmash_similarities.csv` | `qc.sourmash.enabled` | read-sketch similarity matrix |
| `results/final/profile/profile.tsv` | profiler ≠ `none` | MetaPhlAn **percentage relative abundance** |
| `results/final/assembly/quast.tsv` | assembly | QUAST `report.tsv` |
| `results/final/binning/bins.tsv` | MAG | quality joined with taxonomy and source sample |
| `results/final/binning/checkm2/quality_summary.tsv` | MAG | `genome, completeness, contamination, passes_qa` |
| `results/final/binning/drep/GenomeInfo.csv` | MAG | dRep's own table |
| `results/final/binning/gtdb/gtdb_summary.tsv` | MAG | |
| `results/final/diagnostic/bottleneck.tsv` | always | empty cell = **not measured**, never zero |
| `results/final/taxonomy/contig_taxonomy.tsv` | taxonomy + assembly | |
| `results/final/annotation/function_abundance.tsv` | annotation + assembly | `bin, source, function, genes` |
| `results/final/heterogeneity/strain_heterogeneity.tsv` | heterogeneity + binning | blank rate = no position passed the depth filter |
| `results/final/benchmarks/runtime.tsv` | always | |
| `results/reports/multiqc/multiqc_report.html` | `report.multiqc` | |
| `results/provenance/software_versions.tsv` | always | `tool, version, evidence` |
| `results/provenance/pipeline_revision.tsv` | always | `field, value` |
| `results/provenance/reference_db.md5` | always | reference checksums |
| `results/provenance/run_manifest.json` | always | `effective_config`, `source_sha256`, `samplesheet_sha256` |
| `results/reports/assembly/quast.html` | assembly | QUAST's own HTML report |

## Dependencies

**Required**: `python` ≥ 3.11 and `pyyaml` ≥ 6.0 to read and write the config;
`snakemake` ≥ 9.0 on PATH, upstream pinning `snakemake=9.11.2`
(`conda create -n metagenomics -c conda-forge -c bioconda snakemake=9.11.2`);
`conda` or `mamba` on PATH, since every tool environment is built on first run;
and `git`, to resolve the pinned checkout.

**Optional**: an existing checkout via `--pipeline-dir`, which avoids the clone.
Not optional for `--preset full`: CheckM2, GTDB-Tk, contig taxonomy and
functional annotation each need external databases upstream does not bundle.

## Gotchas

- **MetaPhlAn values are not read counts.** The model will want to report
  `final/profile/profile.tsv` values as reads. Do not; upstream profiles with
  relative abundance, so every number is a percentage of the sample.
- **A blank cell is not zero.** The model will want to read the empty
  `assembly_capture_pct` in `bottleneck.tsv` or the empty `heterogeneity_pct`
  as 0%. Do not; blank means the axis could not be measured, which is a different
  claim about the data. Render it as `not measured`.
- **Passing thresholds is not HQ-MAG status.** The model will want to call a bin
  that clears `mag.checkm2.min_completeness` and `max_contamination` a
  "high-quality MAG". Do not; upstream's own docs say threshold passing alone
  does not establish that, and no coverage or marker-gene evidence is here.
- **Demo output is not biology.** The model will want to interpret the demo's
  reads, bins or taxa. Do not; the synthetic fragments carry no taxonomy or
  marker-gene truth. Say the demo proves execution only.
- **`--preset full` without databases fails at config time.** The model will
  want to run `full` to get everything. Do not; `taxonomy`, `annotation` and the
  MAG stage each need external databases, and a half-supplied set is refused
  rather than quietly leaving the stage off.
- **A flag the run cannot honour is an error, never a warning.** `--demo
  --preset full`, or `--preset assembly --checkm2-db X`, exit 1 with
  `INVALID_CONFIG` naming every offending flag at once. The model will want to
  "run it anyway and note the limitation". Do not; the alternative is a user
  believing a binned or annotated analysis ran when no such rule was enabled.
- **Relative config paths resolve against the pipeline root.** The model will
  want to write `results_dir: results` and run from the output directory. Do
  not; upstream's `resolve()` joins relative paths to the checkout, not the
  working directory. Every path this skill writes is absolute.
- **`qc.adapter.fasta` is the one relative path left alone.** It lives inside the
  checkout, so making it absolute would point at a file that does not exist.
- **Two binners are the minimum, not a suggestion.** The model will want to
  run one binner to keep things simple. Do not; DAS Tool consensus refuses
  fewer than two binners at parse time, so `--binners metabat` alone fails.
- **A missing MAG is usually depth, not failure.** The model will want to read
  an absent MAG as a pipeline defect. Do not; at 35x over eight genomes this
  pipeline recovers 4 of 8 bacteria as clean MAGs and misses both 2%-abundance
  yeasts. That is the input's depth, not a bug.
- **No GTDB-Tk database means no MAG tables at all.** The model will want to
  configure binning with `mag.gtdbtk.db` unset and assume only taxonomy is
  lost. Do not; the run fails at classification and `bins.tsv`, `bottleneck.tsv`
  and `runtime.tsv` are never written. Pass `--gtdbtk-db` or turn MAG off.
- **A single-end library cannot be assembled.** The model will want to assemble
  SE data with metaSPAdes. Do not; upstream requires paired-end reads. Use
  `--preset qc` or `--preset profile` for single-end input, and remember a
  co-assembled bin belongs to a `group`, not to one sample.

**How to populate this section:** run the skill 10 times with varied inputs.
Every correction goes here and into a regression test.

## Safety

- **Local-first**: no user or patient data is uploaded; reads, references and
  databases are passed to snakemake as local paths. Network use is limited to
  public reference downloads and is listed on stderr before the run starts: tool
  conda environments on first run, `sra_run` rows (only with
  `--allow-remote-inputs`), the GRCh38 host FASTA (~3 GB) when host subtraction is
  on and no `--host-fasta` is given, the MetaPhlAn index when profiling is on, and
  the MinPath script when its `data_dir` is set.
- **Host subtraction is on by default** because human reads inside a metagenome
  are personal data. Pass `--skip-host-removal` only when the sample is
  demonstrably non-human.
- **Disclaimer**: every report and every `result.json` carries it —
  ClawBio is a research and educational tool. It is not a medical device and does
  not provide clinical diagnoses. Consult a healthcare professional before making
  any medical decisions.
- **No hallucinated science**: every threshold, database and tool version in a
  report is read from the effective config or from a pipeline output.
- **Audit trail**: the exact command, the environment and the checksums of the
  outputs go into `reproducibility/`.

## Agent Boundary

The agent (LLM) dispatches, explains and interprets. The skill executes: it
validates the samplesheet, renders the config, launches snakemake and parses its
tables. The agent must not override thresholds, invent database paths, edit the
rendered `config.yaml` by hand behind the skill's back, or describe a bin as
high-quality on threshold passing alone.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when the user names
assembly, binning, MAG recovery, metaSPAdes, DAS Tool, CheckM2, GTDB-Tk, dRep or
MetaPhlAn; or the input is a shotgun FASTQ pair or an nf-core-style samplesheet
and the user wants contigs or MAGs rather than read-level calls.

**Chaining partners**: this skill connects with:
- `claw-metagenomics`: read-based profiling (Kraken2/RGI/HUMAnN3) before this
  skill, or instead of it when the user only needs read-level calls
- `multiqc-reporter`: consolidate the QC logs when only aggregation is wanted
- `phylogenetics-builder`: build trees from the recovered bin FASTAs
- `bio-orchestrator`: routes here when the request needs assembly or MAG recovery

> `result.json` carries structured stage findings, so it can chain.

## Maintenance

- **Review cadence**: re-check the pinned commit quarterly. The pin is the only
  coupling, so an update is a one-line change plus a re-run of the tests.
- **Staleness signals**: a CheckM2 or GTDB-Tk release the pinned database is
  incompatible with; a binner replacing VAMB or DAS Tool; a snakemake release that
  changes `--conda-frontend`; upstream tagging a release (move to a tag-plus-commit
  pin).
- **Deprecation**: if upstream rewrites its config schema, freeze the pin, publish
  the supported commit range here and archive the skill under
  `skills/_deprecated/` rather than silently reinterpreting the config.

## Citations
- [mobashirrahman/metagenomics-workflow](https://github.com/mobashirrahman/metagenomics-workflow)
  — the workflow this skill wraps, MIT, pinned to `7351a702d29857801963cae7aeff6600d59ebe77` (tag `v0.2.0`)
- QC and fetching: [BBTools](https://jgi.doe.gov/data-and-tools/software-tools/bbtools/),
  [sra-tools](https://github.com/ncbi/sra-tools)
- Profiling and assembly: [MetaPhlAn](https://github.com/biobakery/MetaPhlAn),
  [SPAdes](https://ablab.github.io/spades/), [QUAST](https://quast.sourceforge.net/)
- Binning: [MetaBAT2](https://bitbucket.org/berkeleylab/metabat/),
  [MaxBin2](https://sourceforge.net/projects/maxbin2/),
  [VAMB](https://github.com/RasmussenLab/vamb), [DAS Tool](https://github.com/cmks/DAS_Tool)
- Genomes: [CheckM2](https://github.com/chklovski/CheckM2),
  [dRep](https://drep.readthedocs.io/), [GTDB-Tk](https://ecogenomics.github.io/GTDBTk/)
- Annotation: [Prodigal](https://github.com/hyattpd/Prodigal),
  [DIAMOND](https://github.com/bbuchfink/diamond), [HMMER](http://hmmer.org/),
  [MinPath](https://github.com/mgtools/MinPath)
- Mapping: [Bowtie2](http://bowtie-bio.sourceforge.net/bowtie2/index.shtml),
  [samtools](https://www.htslib.org/), [bwa-mem2](https://github.com/bwa-mem2/bwa-mem2)
