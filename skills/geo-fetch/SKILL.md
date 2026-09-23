---
name: geo-fetch
description: >-
  Query metadata and download data from the NCBI Gene Expression Omnibus (GEO).
  Works with GEO accessions (GSE series, GSM samples, GPL platforms, GDS
  datasets) to search GEO DataSets, fetch series and sample metadata, list and
  download series matrix / SOFT / MINiML / supplementary files, fetch the SRA
  Run Selector files, and emit a standardised metadata.tsv plus a pipeline-ready
  nf-core/rnaseq or nf-core/scrnaseq samplesheet.csv.
license: MIT
metadata:
  version: "0.1.0"
  author: Nikolai Hecker, UK Dementia Research Institute
  domain: genomics
  tags:
    - geo
    - ncbi
    - transcriptomics
    - samplesheet
    - nf-core
    - sra
    - public-archives
  inputs:
    - name: accession
      type: string
      format:
        - text
      description: GEO accession (GSE, GSM, GPL, GDS)
      required: false
    - name: query
      type: string
      format:
        - text
      description: Free-text search of GEO DataSets
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
      description: Standardised one-row-per-sample table (tables/metadata.tsv)
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
    - path: examples/demo_GSE30720_http.json.gz
      description: >-
        Recorded URL-to-response map from one real GSE30720 run (42-sample
        Arabidopsis thaliana seedling RNA-seq). Non-human, so the fixture
        carries no individual-level human characteristics.
  endpoints:
    cli: python skills/geo-fetch/geo_fetch.py --command {command} --accession {accession} --output {output_dir}
    cli_demo: python skills/geo-fetch/geo_fetch.py --demo --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "📊"
    homepage: https://www.ncbi.nlm.nih.gov/geo/
    os:
      - darwin
      - linux
    trigger_keywords:
      - GEO
      - GSE
      - GSM
      - gene expression omnibus
      - series matrix
      - GEO supplementary
      - SraRunTable
      - SRR_Acc_List
      - samplesheet
---

# 🦖 GEO Fetch

You are **GEO Fetch**, a specialised ClawBio agent for the NCBI Gene Expression
Omnibus. Your role is to turn a GEO accession into series and sample metadata,
processed or raw files, a standardised sample table, or a samplesheet a
pipeline can consume directly.

## Trigger

**Fire this skill when the user says any of:**
- "GEO", "Gene Expression Omnibus"
- "GSE30720", "GSM762080", "GPL11221", "GDS..."
- "get the series matrix for this accession"
- "what samples are in this GEO series"
- "build a samplesheet for nf-core/rnaseq from this GSE"
- "give me the SraRunTable / SRR_Acc_List for this series"
- "search GEO for <topic>"

**Do NOT fire when:**
- The accession is `E-GEOD-*`. That is ArrayExpress's mirror of a GEO series —
  route to `arrayexpress-fetch` if the MAGE-TAB view is wanted, or translate to
  the `GSE` and use this skill.
- The accession belongs to ENA (`PRJEB`, `ERR`), PRIDE (`PXD`) or BioStudies
  (`S-BSST`) — route to the matching skill.
- The user has a DOI or PubMed ID rather than an accession — route to
  `article-data-fetcher`.

## Why This Exists

- **Without it**: you click through GEO, hand-copy GSM ids, chase the series to
  its SRA project, then reshape the sample characteristics into whatever column
  names your pipeline expects — differently for every series.
- **With it**: one command returns a **standardised `metadata.tsv`** whose core
  columns are identical across every ClawBio archive skill, and a
  **pipeline-ready `samplesheet.csv`** already matching the nf-core/rnaseq or
  nf-core/scrnaseq column contract — so the output can be handed straight to a
  pipeline rather than needing a bespoke parsing step each time. It also emits a
  runnable download script for the FASTQs it resolved.
- **Why ClawBio**: GEO sample characteristics are free text with no schema. The
  mapping rules, the null markers and the read-pairing logic are fixed and
  inspectable, not re-derived per series by a model. That is what makes a
  samplesheet safe to run a pipeline on.

## Core Capabilities

1. **Series metadata**: title, taxon, type, platform, sample count, summary.
2. **Sample listing**: every GSM in a series with its title.
3. **File listing**: series matrix, SOFT, MINiML and supplementary files.
4. **Download**: any of those categories from the GEO FTP mirror.
5. **Search**: free-text GEO DataSets search with organism and type filters.
6. **Standardised metadata table**: one row per sample, from each GSM's SOFT
   record.
7. **Run Selector files**: `SraRunTable.csv` and `SRR_Acc_List.txt`.
8. **Pipeline-ready samplesheet**: resolves GSE → SRA project → ENA FASTQ links.
9. **Download script**: bash + optional SLURM header.

## Scope

**One skill, one task.** This skill talks to GEO (and the SRA/ENA links GEO
publishes) and nothing else.

## Input Formats

| Format | Example | Notes |
|--------|---------|-------|
| Series | `GSE30720` | The normal case |
| Sample | `GSM762080` | A single sample |
| Platform / DataSet | `GPL11221`, `GDS...` | Metadata only |
| Search phrase | `"Arabidopsis seedling transcriptome"` | With `--command search` |

## Workflow

1. **Resolve the input**: an accession goes to `metadata`; a phrase to `search`.
2. **Fetch**: E-utilities for metadata, the GEO FTP listing for files.
3. **Per-sample detail** (for `metadata-table`): fetch each GSM's SOFT record
   and map its characteristics onto the core columns.
4. **Resolve to reads** (for `samplesheet`): GSE → linked SRA project → ENA
   Portal for public FASTQ links. For 10x runs `--read-map` is required — see
   the first gotcha.
5. **Report**: write `report.md`, `result.json`, `tables/metadata.tsv`,
   `samplesheet.csv` and the reproducibility bundle into `--output`.
6. **Offer, do not act**: `download-script` writes a script and stops. If the
   user wants it executed, show them the file count and total size, and **ask**
   before using `--run` or `--submit`.

Steps 2–4 are prescriptive. Step 6 is a hard rule.

## CLI Reference

```bash
# Demo — offline, from the recorded fixture
python skills/geo-fetch/geo_fetch.py --demo --output /tmp/geo_demo

# Series metadata and samples
python skills/geo-fetch/geo_fetch.py \
  --command metadata --accession GSE30720 --output /tmp/geo
python skills/geo-fetch/geo_fetch.py \
  --command samples --accession GSE30720 --output /tmp/geo

# Standardised metadata table
python skills/geo-fetch/geo_fetch.py \
  --command metadata-table --accession GSE30720 --output /tmp/geo

# Pipeline-ready samplesheet
python skills/geo-fetch/geo_fetch.py \
  --command samplesheet --accession GSE30720 --assay bulk --output /tmp/geo

# SRA Run Selector files, then a samplesheet built from them offline
python skills/geo-fetch/geo_fetch.py \
  --command runtable --accession GSE30720 --output /tmp/geo
python skills/geo-fetch/geo_fetch.py \
  --command samplesheet --accession GSE30720 --assay scrna \
  --from-runtable /tmp/geo/runtable/SraRunTable.csv --fastq-dir /data/fastq --output /tmp/geo

# Files and downloads
python skills/geo-fetch/geo_fetch.py --command files --accession GSE30720 --output /tmp/geo
python skills/geo-fetch/geo_fetch.py --command download --accession GSE30720 --matrix --output /tmp/geo

# Opt in to NCBI rate-limit credentials (never sent otherwise)
python skills/geo-fetch/geo_fetch.py \
  --command search --query "Arabidopsis" --use-ncbi-credentials --output /tmp/geo

# Via the ClawBio runner
python clawbio.py run geo-fetch --demo
python clawbio.py run geo-fetch --command metadata --accession GSE30720

# The upstream positional form also works when called directly
python skills/geo-fetch/geo_fetch.py metadata GSE30720 --output /tmp/geo
```

## Demo

```bash
python clawbio.py run geo-fetch --demo
```

Runs `metadata`, `samples`, `files`, `metadata-table`, `samplesheet` and
`download-script` against the recorded GSE30720 fixture, entirely offline,
producing a 42-sample metadata table and samplesheet.

## Algorithm / Methodology

1. **Metadata**: `esearch` on `db=gds` for the accession, then `esummary`.
2. **Samples**: the series summary's GSM list, then each sample's SOFT record
   from `geo/query/acc.cgi`.
3. **Files**: the GEO FTP directory listing, sharded by accession — `GSE30720`
   lives under `series/GSE30nnn/GSE30720/`.
4. **Run Selector**: `esearch db=sra` with history, then `efetch
   rettype=runinfo` for the full table, plus the run accession list.
5. **Samplesheet**: resolve the series to its SRA project, query the ENA Portal
   for that project, and pair `fastq_ftp` into R1/R2.
6. **Harmonisation**: GSM `characteristics` lines are `key: value` free text;
   synonyms map onto the core columns and everything unmapped is promoted to
   its own column so nothing is lost.

**Key parameters**
- Core columns: `sample, replicate, species, sex, age, condition, genotype, treatment, tissue`
- Missing value token: `NA`
- Samplesheet columns: `sample,fastq_1,fastq_2` (scrna) plus `strandedness` (bulk)
- Credentials: `NCBI_EMAIL`, `NCBI_API_KEY` — only with `--use-ncbi-credentials`

## Example Queries

- "Get the metadata for GSE30720"
- "Build an nf-core/rnaseq samplesheet for this GEO series"
- "List the supplementary files for this series"
- "Give me the SRR accession list for this GSE"

## Example Output

```csv
sample,fastq_1,fastq_2,strandedness
SRS243343,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR342/SRR342351/SRR342351_1.fastq.gz,https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR342/SRR342351/SRR342351_2.fastq.gz,auto
```

```text
accession   : GSE30720
title       : Seedling transcriptome sequencing of the Arabidopsis thaliana MAGIC founder accessions
taxon       : Arabidopsis thaliana
gdstype     : Expression profiling by high throughput sequencing
n_samples   : 42
```

## Output Structure

```
output_directory/
├── report.md              # Commands run and what each returned
├── result.json            # Machine-readable envelope
├── samplesheet.csv        # Pipeline-ready nf-core samplesheet
├── download_geo.sh        # Runnable bash + SLURM download script
├── tables/
│   └── metadata.tsv       # Standardised one-row-per-sample table
├── runtable/              # (optional) only with --command runtable
├── downloads/             # (optional) only with --command download
└── reproducibility/
    ├── commands.sh        # Exact command to reproduce
    ├── environment.yml    # Environment snapshot
    └── checksums.sha256   # SHA-256 of every artifact
```

## Dependencies

**Required**: Python >= 3.10 only. The vendored client is standard library.

**Optional**: `wget` or `curl` on the machine that *runs* the generated
download script; `sbatch` if it is submitted.

## Gotchas

- **Gotcha 1**: You will want to trust the `_1`/`_2` filename heuristic for a
  10x run. Do not. When the technical reads are separate files it silently
  picks the barcode read and **drops the cDNA read**. Pass `--read-map`: 3 files
  (single index) → `2,3`; 4 files (dual index) → `3,4`. This skill has no
  `runs` command — that is `ena-fetch`'s — so check the file count in the
  `fastq_ftp` column of `--command runtable`, or ask `ena-fetch` for the run,
  before choosing. Confirm the choice with the user either way.

  The more reliable route for Chromium is `fasterq-dump`, which is the only
  tool that exposes the read structure faithfully. There is no ClawBio skill
  for it yet, so do not promise one: use `--command runtable` here to write
  `SRR_Acc_List.txt`, then run sra-tools directly —
  `prefetch --option-file SRR_Acc_List.txt` followed by
  `fasterq-dump --split-files <SRR>`.
- **Gotcha 2**: Not every GEO series is mirrored to ENA. Recent submissions
  frequently are not, and `samplesheet` then fails with "No public FASTQ found"
  even though the series exists and `metadata` works. That is accurate, not a
  bug. Do not invent links. Use `--command runtable` to get `SRR_Acc_List.txt`
  and fetch from SRA with sra-tools (`prefetch` + `fasterq-dump`).
- **Gotcha 3**: You will want to set `NCBI_EMAIL` and `NCBI_API_KEY` and expect
  higher rate limits. They are **read but never sent** unless
  `--use-ncbi-credentials` is passed. Presence of an environment variable is
  not consent. Without the flag, stay under ~3 requests/second.
- **Gotcha 4**: `metadata-table` fetches **one SOFT record per sample**, so a
  700-sample series is 700 requests. Check `n_samples` from `metadata` first
  and warn the user before running it on a large series.
- **Gotcha 5**: Upstream's `runtable --out` defaulted to `.`, the working
  directory. Here it resolves under `--output`, into `runtable/`. Do not
  reintroduce the cwd default.
- **Gotcha 6**: `download-script` **writes a script and downloads nothing**.
  Never run or submit it without telling the user the file count and total size
  first. Note the FASTQ URLs it writes point at `ftp.sra.ebi.ac.uk`, not at
  NCBI — this skill reads run metadata from ENA. So the script can fail on a
  network where this skill itself worked fine, because the hosts differ. If it
  does, check that host is allowlisted before assuming a bad accession.

## Safety

- **Local-first**: no user genetic data is ever transmitted. This skill sends a
  public accession or the search phrase you typed to `eutils.ncbi.nlm.nih.gov`,
  `ftp.ncbi.nlm.nih.gov` and `www.ebi.ac.uk`, and receives public archive data.
  Nothing of yours leaves the machine, satisfying ClawBio Safety Rule 1 by
  construction. See [docs/data-handling.md](../../docs/data-handling.md).
- **Credentials**: `NCBI_EMAIL` and `NCBI_API_KEY` are sent **only** with
  `--use-ncbi-credentials`, and `result.json` records whether they were.
- **Execution is opt-in**: `download-script` only writes a file.
- **Disclaimer**: every report carries the ClawBio medical disclaimer.
- **Overwrite**: the skill warns on stderr before overwriting an output directory.
- **Audit trail**: every run writes `reproducibility/`.

## Agent Boundary

The agent dispatches, explains, and asks before anything is executed. The skill
executes. The agent must not invent accessions or FASTQ URLs, must not guess a
`--read-map` it has not checked against the run's file count, must not pass
`--use-ncbi-credentials` unless the user asked for it, and must not run or
submit a generated script without explicit confirmation.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on a GEO accession (`GSE`,
`GSM`, `GPL`, `GDS`) or an explicit mention of GEO.

**Chaining partners**:
- **sra-tools** (external, not a ClawBio skill): `--command runtable` writes
  the `SRR_Acc_List.txt` that `prefetch` and `fasterq-dump` consume. This is
  the route for 10x/Chromium reads and for series not mirrored to ENA.
- `ena-fetch`: this skill queries ENA for FASTQ links, so the two agree.
- `nfcore-rnaseq-wrapper` / `nfcore-scrnaseq-wrapper`: the natural consumers of
  the `samplesheet.csv` this skill writes.
- `rnaseq-de`: downstream differential expression once counts exist.
- `article-data-fetcher`: **upstream producer.** It resolves a DOI or PMID to
  the repository accessions a paper deposited, GEO among them. When the user
  starts from a paper rather than an accession, run it first and hand the
  accessions here. It downloads files and writes a `manifest.json`, but it
  does **not** harmonise sample annotation into `metadata.tsv` or emit a
  pipeline-ready `samplesheet.csv` — that is this skill's job, so the two
  chain rather than compete.

## Maintenance

- **Update trigger**: when an archive host changes an endpoint this skill calls.
  Not on a calendar — a fixed cadence either fires when nothing has changed or
  misses a break the week after it lands. The staleness signals below are the
  trigger.
- **How a break surfaces**: from a live call, not from CI. The demo and the tests
  run offline from committed fixtures, so they stay green after an endpoint
  changes. Treat an unexpected HTTP error or an empty result on a real accession
  as the signal, then re-check the fixtures against the live API.
- **Staleness signals**: `esummary` field names changing; the GEO FTP sharding
  scheme changing; NCBI tightening unauthenticated rate limits; nf-core
  changing its samplesheet column contract.
- **Known debt**: the harmonisation and read-pairing helpers are duplicated
  across the archive skills rather than shared, deliberately, so each stays
  easy to re-sync with upstream. Factor them out only if upstream does.
- **Deprecation**: if NCBI ships an official GEO client covering these
  commands, wrap it instead of E-utilities.

## Citations

- [NCBI GEO](https://www.ncbi.nlm.nih.gov/geo/) and
  [E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25501/); the endpoints used.
- [ENA Portal API](https://www.ebi.ac.uk/ena/portal/api/); FASTQ link resolution.
- [nf-core/rnaseq](https://nf-co.re/rnaseq) and
  [nf-core/scrnaseq](https://nf-co.re/scrnaseq); the samplesheet contracts.
- Demo fixture: [GSE30720](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE30720),
  Arabidopsis thaliana seedling RNA-seq; GEO records are public domain.
- Ported from
  [UKDRI/informatics_data_skills](https://github.com/UKDRI/informatics_data_skills)
  @ `7cc3e6e` (`geo/`, plus `fastq-download-script/` folded in as
  `download-script`), © 2026 UK Dementia Research Institute, MIT.
