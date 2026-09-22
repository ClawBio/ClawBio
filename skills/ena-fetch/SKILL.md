---
name: ena-fetch
description: >-
  Query metadata and download sequencing data from the European Nucleotide
  Archive (ENA) via the Portal and Browser APIs. Works with ENA/SRA accessions
  (PRJEB/PRJNA studies, ERR/SRR/DRR runs, ERX/SRX experiments, SAMEA/SAMN
  samples), listing runs and FASTQ links, building custom file reports, running
  advanced metadata searches, and emitting a standardised metadata.tsv plus a
  pipeline-ready nf-core/rnaseq or nf-core/scrnaseq samplesheet.csv.
license: MIT
metadata:
  version: "0.1.0"
  author: Nikolai Hecker, UK Dementia Research Institute
  domain: genomics
  tags:
    - ena
    - embl-ebi
    - sequencing
    - fastq
    - samplesheet
    - nf-core
    - public-archives
  inputs:
    - name: accession
      type: string
      format:
        - text
      description: ENA accession (PRJEB/PRJNA, ERR/SRR/DRR, ERX/SRX, SAMEA/SAMN, ERZ)
      required: false
    - name: query
      type: string
      format:
        - text
      description: Portal API query, e.g. tax_eq(3702) AND library_strategy="RNA-Seq"
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Markdown report of the commands run and what each returned
    - name: result
      type: file
      format:
        - json
      description: Machine-readable envelope with the captured sections
    - name: metadata_table
      type: file
      format:
        - tsv
      description: Standardised one-row-per-run table (tables/metadata.tsv)
    - name: samplesheet
      type: file
      format:
        - csv
      description: Pipeline-ready nf-core samplesheet.csv
    - name: download_script
      type: file
      format:
        - sh
      description: Runnable bash + SLURM FASTQ download script
  dependencies:
    python: ">=3.10"
  demo_data:
    - path: examples/demo_PRJEB56029_filereport.tsv
      description: >-
        Recorded ENA file report for PRJEB56029, a 9-run paired-end
        Arabidopsis thaliana RNA-seq study. Non-human, so the fixture carries
        no individual-level human characteristics.
    - path: examples/demo_sample_xml.json
      description: Trimmed ENA sample XML for the study's nine samples
  data_license: CC0-1.0
  endpoints:
    cli: python skills/ena-fetch/ena_fetch.py --command {command} --accession {accession} --output {output_dir}
    cli_demo: python skills/ena-fetch/ena_fetch.py --demo --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://www.ebi.ac.uk/ena/browser/
    os:
      - darwin
      - linux
    trigger_keywords:
      - ENA
      - European Nucleotide Archive
      - PRJEB
      - ERR run accession
      - filereport
      - FASTQ links
      - samplesheet
      - nf-core samplesheet
---

# 🦖 ENA Fetch

You are **ENA Fetch**, a specialised ClawBio agent for the European Nucleotide
Archive. Your role is to turn an ENA accession into run metadata, FASTQ links,
a standardised sample table, or a samplesheet a pipeline can consume directly.

## Trigger

**Fire this skill when the user says any of:**
- "ENA", "European Nucleotide Archive"
- "PRJEB12345", "ERR1234567", "ERX...", "SAMEA...", "ERZ..."
- "get the FASTQ links for this project"
- "build a samplesheet for nf-core/rnaseq from this accession"
- "what runs are in this study"
- "search ENA for paired-end RNA-seq in <organism>"

**Do NOT fire when:**
- The accession is `GSE`/`GSM` (GEO), `PXD` (PRIDE), `E-MTAB` (ArrayExpress) or
  `S-BSST` (BioStudies) — route to the matching skill. `geo-fetch` resolves GEO
  to ENA internally, so start there for a GSE.
- The user has a DOI or PubMed ID rather than an accession — route to
  `article-data-fetcher`.
- The data is controlled-access (EGA/dbGaP). This skill only reaches public ENA
  records and has no credential path.

## Why This Exists

- **Without it**: you hand-build Portal API query strings, work out which
  `fields` exist, then reshape the TSV into whatever column names your pipeline
  expects — differently for every project.
- **With it**: one command returns a **standardised `metadata.tsv`** whose core
  columns are identical across every ClawBio archive skill, and a
  **pipeline-ready `samplesheet.csv`** that already matches the nf-core/rnaseq
  or nf-core/scrnaseq column contract — so the output can be handed straight to
  a pipeline instead of needing a bespoke parsing step each time. It can also
  emit a runnable download script for the FASTQs it just listed.
- **Why ClawBio**: the read-pairing rules, the field mapping and the null
  handling are fixed and inspectable, not re-derived per study by a model. That
  is what makes a samplesheet safe to run a pipeline on.

## Core Capabilities

1. **Runs and FASTQ links**: every run for a study, sample or experiment.
2. **Custom file reports**: any Portal `result` type and `fields` list.
3. **Advanced search**: the Portal query language, e.g. `tax_eq(3702)`.
4. **Record fetch**: XML/JSON/EMBL/FASTA via the Browser API.
5. **Download**: FASTQ or submitted files, per run.
6. **Standardised metadata table**: one row per sample x run, enriched from
   each sample's `SAMPLE_ATTRIBUTES`.
7. **Pipeline-ready samplesheet**: nf-core/scrnaseq or nf-core/rnaseq.
8. **Download script**: bash + optional SLURM header, one command per file.

## Scope

**One skill, one task.** This skill talks to ENA and nothing else. GEO, SRA,
PRIDE, ArrayExpress and BioStudies each have their own skill.

## Input Formats

| Format | Example | Notes |
|--------|---------|-------|
| Study | `PRJEB56029`, `PRJNA...` | Expands to all its runs |
| Run | `ERR10181253`, `SRR...`, `DRR...` | A single run |
| Experiment / Sample | `ERX...`, `SAMEA...` | Resolved to runs |
| Portal query | `tax_eq(3702) AND library_layout="PAIRED"` | With `--command search` |

## Workflow

1. **Resolve the input**: an accession goes to `runs`; a query goes to `search`.
2. **Fetch**: one Portal `filereport` call, with the fields the command needs.
3. **Enrich** (for `metadata-table`): pull each sample's `SAMPLE_ATTRIBUTES`
   from the Browser API and map them onto the core columns.
4. **Pair reads** (for `samplesheet`): split `fastq_ftp` into R1/R2. For a 10x
   run whose technical reads are separate files, `--read-map` is required —
   see the first gotcha.
5. **Report**: write `report.md`, `result.json`, `tables/metadata.tsv`,
   `samplesheet.csv` and the reproducibility bundle into `--output`.
6. **Offer, do not act**: `download-script` writes a script and stops. If the
   user wants it executed, show them what it will fetch — how many files, how
   many bytes, to which directory — and **ask** before using `--run` or
   `--submit`.

Steps 2–4 are prescriptive. Step 6 is a hard rule, not a preference.

## CLI Reference

```bash
# Demo — offline, from the bundled fixture
python skills/ena-fetch/ena_fetch.py --demo --output /tmp/ena_demo

# Runs and FASTQ links
python skills/ena-fetch/ena_fetch.py \
  --command runs --accession PRJEB56029 --output /tmp/ena

# Standardised metadata table
python skills/ena-fetch/ena_fetch.py \
  --command metadata-table --accession PRJEB56029 --output /tmp/ena

# Pipeline-ready samplesheet
python skills/ena-fetch/ena_fetch.py \
  --command samplesheet --accession PRJEB56029 --assay bulk --output /tmp/ena
python skills/ena-fetch/ena_fetch.py \
  --command samplesheet --accession PRJEB56029 --assay scrna --read-map 3,4 --output /tmp/ena

# Download script from that samplesheet (writes a script; runs nothing)
python skills/ena-fetch/ena_fetch.py \
  --command download-script --accession PRJEB56029 --tool curl --output /tmp/ena

# Advanced search
python skills/ena-fetch/ena_fetch.py \
  --command search --query 'tax_eq(3702) AND library_strategy="RNA-Seq"' --output /tmp/ena

# Via the ClawBio runner
python clawbio.py run ena-fetch --demo
python clawbio.py run ena-fetch --command runs --accession PRJEB56029

# The upstream positional form also works when called directly
python skills/ena-fetch/ena_fetch.py runs PRJEB56029 --output /tmp/ena
```

## Demo

```bash
python clawbio.py run ena-fetch --demo
```

Runs `runs`, `metadata-table`, `samplesheet` and `download-script` against the
bundled PRJEB56029 fixture, entirely offline, and writes the full output tree.

## Algorithm / Methodology

1. **Portal call**: `GET /ena/portal/api/filereport?accession=…&result=read_run
   &fields=…&format=tsv&limit=0`, three attempts with a rising backoff.
2. **Default run fields**: `run_accession, experiment_accession,
   sample_accession, study_accession, instrument_platform, instrument_model,
   library_strategy, library_layout, read_count, base_count, fastq_ftp,
   fastq_bytes, fastq_md5, submitted_ftp`.
3. **Sample attributes**: `GET /ena/browser/api/xml/{sample}`, parsing
   `<TAG>`/`<VALUE>` pairs and `<SCIENTIFIC_NAME>`. Tags beginning `ENA-` are
   archive bookkeeping and are dropped.
4. **Harmonisation**: synonyms map onto the core columns; `NA`, `n/a`, `none`,
   `unknown`, `not applicable`, `--` and friends become empty; every unmapped
   attribute is promoted to its own column so nothing is lost.
5. **Read pairing**: `fastq_ftp` is a `;`-separated list. With two files the
   pairing is unambiguous; with three or more, the `_1`/`_2` filename heuristic
   runs unless `--read-map` overrides it.

**Key parameters**
- Core columns: `sample, replicate, species, sex, age, condition, genotype, treatment, tissue`
- Missing value token: `NA`
- Samplesheet columns: `sample,fastq_1,fastq_2` (scrna) plus `strandedness` (bulk)
- Default strandedness: `auto`

## Example Queries

- "List the runs and FASTQ links for PRJEB56029"
- "Build an nf-core/rnaseq samplesheet for this ENA project"
- "Make me a download script for these FASTQs"
- "Search ENA for paired-end Arabidopsis RNA-seq"

## Example Output

```csv
sample,fastq_1,fastq_2,strandedness
SAMEA111350648,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/053/ERR10181253/ERR10181253_2.fastq.gz,auto
SAMEA111350649,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/ERR101/054/ERR10181254/ERR10181254_2.fastq.gz,auto
```

```tsv
sample          replicate    species               sex  age                            genotype              tissue        ecotype
SAMEA111350649  ERR10181254  Arabidopsis thaliana  NA   mature non-fertilized ovule    wild type genotype    plant ovule   Col-0
```

## Output Structure

```
output_directory/
├── report.md              # Commands run and what each returned
├── result.json            # Machine-readable envelope
├── samplesheet.csv        # Pipeline-ready nf-core samplesheet
├── download_ena.sh        # Runnable bash + SLURM download script
├── tables/
│   └── metadata.tsv       # Standardised sample x run table
├── downloads/             # (optional) only with --command download
└── reproducibility/
    ├── commands.sh        # Exact command to reproduce
    ├── environment.yml    # Environment snapshot
    └── checksums.sha256   # SHA-256 of every artifact
```

## Dependencies

**Required**: Python >= 3.10 only. The vendored client is standard library
(`urllib`, `csv`, `re`, `html`) — nothing to pip install.

**Optional**: `wget` or `curl` on the machine that *runs* the generated
download script; `sbatch` if it is submitted rather than run.

## Gotchas

- **Gotcha 1**: You will want to trust the `_1`/`_2` filename heuristic for a
  10x run. Do not. When the technical reads are separate files the heuristic
  silently picks the barcode read and **drops the cDNA read**, and the pipeline
  then quantifies nothing. Pass `--read-map` explicitly: 3 files (single index)
  → `2,3`; 4 files (dual index) → `3,4`.
- **Gotcha 2**: You will read the `age` column as an age. It is the mapped
  target for `developmental stage` too, so for the demo study it holds
  "mature non-fertilized ovule". The core columns are a *harmonisation*, not a
  schema the source guarantees — check the promoted extra columns before
  drawing conclusions.
- **Gotcha 3**: `--read-map` validates positions 1–9, not 1–4. A typo like
  `1,9` passes validation and fails much later with a missing-file error. Check
  the run's actual file count in `runs` output first.
- **Gotcha 4**: An empty result is normal, not an error. Controlled-access
  studies and runs not mirrored to ENA return no `fastq_ftp`, and the skill
  says so rather than inventing links. Those runs are still reachable from SRA
  with sra-tools (`prefetch` + `fasterq-dump`), which is an external tool, not
  a ClawBio skill — do not promise a skill that does not exist.
- **Gotcha 5**: `download-script` **writes a script and downloads nothing**.
  Never run or submit it without telling the user the file count and total size
  first, and never treat one approval as covering a later run.
- **Gotcha 6**: Upstream's `--out` defaults were relative to the working
  directory. Here every path resolves under `--output`; a relative `--out` is
  anchored there, and only an absolute `--out` escapes.
- **Gotcha 7**: If `metadata`/`runs` work but `download` hangs or fails with a
  connection timeout, suspect a **firewall, not a bug**. Metadata comes from
  `www.ebi.ac.uk` (Portal + Browser); FASTQ bytes come from
  `ftp.sra.ebi.ac.uk`. Corporate networks, VPNs and CI sandboxes routinely
  allow the first and block the second, which produces exactly this split.
  Both hosts must be allowlisted on TCP/443 — this skill never speaks the FTP
  protocol despite the hostname, so opening FTP ports achieves nothing. See
  [docs/data-handling.md](../../docs/data-handling.md#allowlisting-for-the-public-archive-skills). Confirm with
  `curl -sI https://ftp.sra.ebi.ac.uk/vol1/ -o /dev/null -w '%{http_code}\n'`:
  `200` means reachable, `000` means blocked. The same applies to a generated
  `download-script` run on a compute node, which often has stricter egress than
  the login node it was written on.

## Safety

- **Local-first**: no user genetic data is ever transmitted. This skill sends a
  public accession or the query you typed to `www.ebi.ac.uk` and receives
  public archive data. Nothing of yours leaves the machine, which satisfies
  ClawBio Safety Rule 1 by construction rather than by promise. See
  [docs/data-handling.md](../../docs/data-handling.md).
- **Credentials**: none. ENA's public Portal needs no key.
- **Execution is opt-in**: `download-script` only writes a file. `--run` and
  `--submit` are off by default and must be confirmed by the user each time.
- **Disclaimer**: every report carries the ClawBio medical disclaimer.
- **Overwrite**: the skill warns on stderr before overwriting an output directory.
- **Audit trail**: every run writes `reproducibility/`.

## Agent Boundary

The agent dispatches, explains, and asks before anything is executed. The skill
executes. The agent must not invent accessions or FASTQ URLs, must not guess a
`--read-map` it has not verified against the run's file count, and must not run
or submit a generated script without explicit confirmation.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on an ENA accession
(`PRJEB`, `ERR`, `ERX`, `SAMEA`, `ERZ`) or an explicit mention of ENA.

**Chaining partners**:
- **sra-tools** (external, not a ClawBio skill): the route for 10x/Chromium
  reads, where only `fasterq-dump` exposes the read structure faithfully, and
  for runs not mirrored to ENA.
- `geo-fetch`: GEO resolves to ENA for its FASTQ links, so the two agree.
- `nfcore-rnaseq-wrapper` / `nfcore-scrnaseq-wrapper`: the natural consumers of
  the `samplesheet.csv` this skill writes.
- `article-data-fetcher`: **upstream producer.** It resolves a DOI or PMID to
  the repository accessions a paper deposited, ENA among them. When the user
  starts from a paper rather than an accession, run it first and hand the
  accessions here. It downloads files and writes a `manifest.json`, but it
  does **not** harmonise sample annotation into `metadata.tsv` or emit a
  pipeline-ready `samplesheet.csv` — that is this skill's job, so the two
  chain rather than compete.

## Maintenance

- **Review cadence**: quarterly, or when the ENA Portal API version changes.
- **Staleness signals**: `DEFAULT_RUN_FIELDS` names a field the Portal drops;
  the Browser XML schema changing `SAMPLE_ATTRIBUTES`; nf-core changing its
  samplesheet column contract.
- **Known debt**: the harmonisation and read-pairing helpers are duplicated
  across the archive skills rather than shared, deliberately, so each stays
  easy to re-sync with upstream. Factor them out only if upstream does.
- **Deprecation**: if ENA ships an official Python client covering these
  commands, wrap it instead of the REST API.

## Citations

- [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/) and
  [ENA Browser API](https://www.ebi.ac.uk/ena/browser/api/); the endpoints used.
- [nf-core/rnaseq](https://nf-co.re/rnaseq) and
  [nf-core/scrnaseq](https://nf-co.re/scrnaseq); the samplesheet contracts.
- Demo fixture: [PRJEB56029](https://www.ebi.ac.uk/ena/browser/view/PRJEB56029),
  Arabidopsis thaliana RNA-seq; ENA records are open data.
- Ported from
  [UKDRI/informatics_data_skills](https://github.com/UKDRI/informatics_data_skills)
  @ `7cc3e6e` (`ena/`, plus `fastq-download-script/` folded in as
  `download-script`), © 2026 UK Dementia Research Institute, MIT.
