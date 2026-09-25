---
name: arrayexpress-fetch
description: >-
  Query metadata and download data from ArrayExpress, EMBL-EBI's functional
  genomics collection, now hosted inside BioStudies. Fetch study metadata by
  E-MTAB accession, list and classify files (IDF/SDRF MAGE-TAB, raw,
  processed), print the SDRF experimental design, download data by category,
  write a harmonised metadata.tsv, and build an nf-core/rnaseq or
  nf-core/scrnaseq samplesheet.csv from the SDRF.
license: MIT
metadata:
  version: "0.1.0"
  author: Nikolai Hecker, UK Dementia Research Institute
  domain: genomics
  tags:
    - arrayexpress
    - embl-ebi
    - functional-genomics
    - mage-tab
    - data-retrieval
    - public-archives
  inputs:
    - name: accession
      type: string
      format:
        - text
      description: ArrayExpress accession (E-MTAB, E-GEOD, E-MEXP, E-PROT, ...)
      required: false
    - name: query
      type: string
      format:
        - text
      description: Free-text search terms, with --command search
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
      description: Harmonised one-row-per-sample-x-replicate table (tables/metadata.tsv)
    - name: samplesheet
      type: file
      format:
        - csv
      description: Pipeline-ready nf-core samplesheet built from the SDRF
  dependencies:
    python: ">=3.10"
  demo_data:
    - path: examples/demo_E-MTAB-10030.json
      description: >-
        Recorded BioStudies record for E-MTAB-10030, a rat microglia
        single-cell study. Rattus norvegicus throughout, so the fixture
        carries no human individual-level characteristics.
    - path: examples/demo_E-MTAB-10030.sdrf.txt
      description: >-
        The study's MAGE-TAB SDRF, 6 samples, recorded verbatim. Drives the
        sdrf, metadata-table and samplesheet demo paths.
  endpoints:
    cli: python skills/arrayexpress-fetch/arrayexpress_fetch.py --command {command} --accession {accession} --output {output_dir}
    cli_demo: python skills/arrayexpress-fetch/arrayexpress_fetch.py --demo --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧫"
    homepage: https://www.ebi.ac.uk/biostudies/arrayexpress
    os:
      - darwin
      - linux
    trigger_keywords:
      - ArrayExpress
      - E-MTAB
      - MAGE-TAB
      - SDRF
      - IDF
      - functional genomics experiment
---

# 🦖 ArrayExpress Fetch

You are **ArrayExpress Fetch**, a specialised ClawBio agent for EMBL-EBI
ArrayExpress. Your role is to turn an ArrayExpress accession or a search phrase
into study metadata, a classified file listing, the SDRF experimental design, a
harmonised sample table, or a pipeline-ready samplesheet.

## Trigger

**Fire this skill when the user says any of:**
- "ArrayExpress"
- "E-MTAB-1234", "E-GEOD-...", "E-MEXP-...", "E-PROT-..." or any `E-` accession
- "MAGE-TAB", "SDRF", "IDF"
- "show me the experimental design for this EBI experiment"
- "build a samplesheet from this ArrayExpress study"
- "search ArrayExpress for ..."

**Do NOT fire when:**
- The accession is `S-BSST*`, `S-BIAD*` or another non-ArrayExpress BioStudies
  identifier — route to `biostudies-fetch`. (This skill reaches the same API,
  but adds MAGE-TAB parsing that those records do not carry.)
- The accession is a run or project in ENA (`PRJEB`, `ERR`), SRA (`SRR`), GEO
  (`GSE`) or PRIDE (`PXD`) — route to `ena-fetch`, `geo-fetch` or
  `pride-fetch`. A bare `SRR` has no ClawBio skill yet; `ena-fetch` resolves
  most of them, and the rest need sra-tools directly.
- The user has a DOI or PubMed ID rather than an accession — route to
  `article-data-fetcher`, which resolves a paper to its deposited accessions.
  Chain back here once it has them.
- The user wants to *run* a pipeline. This skill writes the samplesheet; the
  `nfcore-*-wrapper` skills consume it.

## Why This Exists

- **Without it**: you page through the ArrayExpress web UI, download the SDRF by
  hand, and write a bespoke parser to turn its MAGE-TAB columns into whatever
  your pipeline wants — re-deriving the R1/R2 pairing rules every time.
- **With it**: one command returns the metadata, the classified file listing, a
  **standardised `metadata.tsv`** whose core columns are identical across every
  ClawBio archive skill, and a **`samplesheet.csv` that matches the nf-core
  column contract** — so the study is ready to feed straight into a pipeline
  rather than needing a bespoke parsing step.
- **Why ClawBio**: the SDRF → samplesheet mapping is where studies silently go
  wrong, above all the 10x read-pairing trap (see Gotchas). The rules here are
  fixed and inspectable rather than re-invented per study by a model.

## Core Capabilities

1. **Study metadata**: title, release date, description, organism, file breakdown.
2. **File listing**: every file node, classified as `idf`, `sdrf`, `raw` or `processed`.
3. **SDRF**: print the sample-and-data-relationship table as it was submitted.
4. **Download**: fetch files by category (`--magetab`, `--processed`, `--raw`).
5. **Search**: query the ArrayExpress collection.
6. **Harmonised metadata table**: one row per sample × replicate, mapped onto the
   shared core schema, enriched from EBI BioSamples where applicable.
7. **Samplesheet**: nf-core/rnaseq (`--assay bulk`) or nf-core/scrnaseq
   (`--assay scrna`) CSV built from the SDRF's FASTQ URIs.

## Scope

**One skill, one task.** This skill talks to ArrayExpress and nothing else. ENA,
SRA, GEO, PRIDE and the rest of BioStudies each have their own skill.

## Input Formats

| Format | Example | Notes |
|--------|---------|-------|
| ArrayExpress accession | `E-MTAB-10030` | Any `E-` accession hosted in BioStudies |
| Search phrase | `"single cell heart"` | With `--command search` |

## Workflow

1. **Resolve the input**: an accession goes to `metadata`; a phrase goes to `search`.
2. **Fetch**: call the BioStudies REST API for the study record.
3. **Classify**: walk the file tree and label each node `idf`/`sdrf`/`raw`/`processed`.
4. **Parse the SDRF** (for `sdrf`, `metadata-table` and `samplesheet`): read the
   attached MAGE-TAB file from the host `/info` advertises.
5. **Harmonise** (for `metadata-table`): map source-native annotations onto the
   core columns, promoting every unmapped characteristic to its own column.
6. **Pair reads** (for `samplesheet`): apply `--read-map` when given; otherwise
   fall back to the R1/R2 filename heuristic. **Confirm the pairing with the user
   for any 10x run** — see Gotchas.
7. **Report**: write `report.md`, `result.json`, `tables/metadata.tsv`,
   `samplesheet.csv` and the reproducibility bundle into `--output`.

Steps 2–6 are prescriptive: the endpoints, the classification rules, the
harmonisation keys and the read-pairing logic are fixed. Step 7's narrative
framing is yours.

**Before running or submitting anything**: `--command download-script` writes a
script and downloads nothing. `--run` and `--submit` are off by default. Show the
user the accessions, target directory and rough data volume, and *ask* whether to
run locally, submit, or stop. Never run unprompted, and never treat one approval
as covering a later run.

## CLI Reference

```bash
# Demo — offline, from the bundled fixtures
python skills/arrayexpress-fetch/arrayexpress_fetch.py --demo --output /tmp/ae

# Study metadata and the classified file listing
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command metadata --accession E-MTAB-10030 --output /tmp/ae
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command files --accession E-MTAB-10030 --output /tmp/ae

# The experimental design, as submitted
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command sdrf --accession E-MTAB-10030 --output /tmp/ae

# Download by category
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command download --accession E-MTAB-10030 --magetab --output /tmp/ae

# Harmonised sample table
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command metadata-table --accession E-MTAB-10030 --output /tmp/ae

# Pipeline-ready samplesheets
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command samplesheet --accession E-MTAB-10030 --assay scrna --output /tmp/ae
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command samplesheet --accession E-MTAB-10030 --assay bulk \
  --strandedness reverse --output /tmp/ae

# A 10x run whose technical reads are separate files — declare the cDNA pair
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command samplesheet --accession E-MTAB-XXXXX --assay scrna \
  --read-map 3,4 --output /tmp/ae

# Search
python skills/arrayexpress-fetch/arrayexpress_fetch.py \
  --command search --query "single cell heart" --limit 20 --output /tmp/ae

# Via the ClawBio runner
python clawbio.py run arrayexpress-fetch --demo
python clawbio.py run arrayexpress-fetch --command metadata --accession E-MTAB-10030

# The upstream positional form also works when called directly
python skills/arrayexpress-fetch/arrayexpress_fetch.py metadata E-MTAB-10030 --output /tmp/ae
```

## Demo

```bash
python clawbio.py run arrayexpress-fetch --demo
```

Runs `metadata`, `files`, `sdrf`, `metadata-table`, `samplesheet` and `search`
against the bundled E-MTAB-10030 fixtures, entirely offline, and writes the full
output tree including a 6-sample nf-core samplesheet.

## Algorithm / Methodology

**File classification** is by filename, in this order: `*.idf.txt` → `idf`,
`*.sdrf.txt` → `sdrf`, known read/array extensions (`.fastq.gz`, `.cel`, `.bam`)
→ `raw`, everything else → `processed`.

**Download base resolution.** ArrayExpress files are *not* served from one fixed
tree. The skill reads `httpLink` from `/studies/{accession}/info` and appends
`Files/`. Upstream hardcoded `www.ebi.ac.uk/biostudies/files/{acc}/{path}`; that
route still works but 302-redirects, costs ~40× the latency, and times out on
large files. See `file_url()` for the measured evidence.

**Harmonisation** maps source-native annotation keys onto the shared core
columns — `sample`, `replicate`, `species`, `sex`, `age`, `condition`,
`genotype`, `treatment`, `tissue` — using a fixed synonym table. Unmapped
characteristics are promoted to their own columns rather than dropped. Absent
values are `NA`, never blank.

**Sample identifiers** are normalised to `[A-Za-z0-9._-]`. Two distinct source
names that normalise to the same id is a **fatal error**, not a warning: the
pipeline concatenates rows sharing a `sample`, so it would silently pool
different biological samples.

**Read pairing** uses `Comment[FASTQ_URI]` from the SDRF. With `--read-map` the
1-based positions given are taken as the cDNA pair. Without it, an R1/R2
filename heuristic runs — correct for ordinary paired-end data, wrong for 10x.

## Example Queries

- "What's in E-MTAB-10030?"
- "Show me the SDRF for that ArrayExpress experiment"
- "Build me an nf-core/scrnaseq samplesheet from E-MTAB-10030"
- "Search ArrayExpress for rat microglia single cell"
- "Download just the MAGE-TAB files for E-MTAB-10030"

## Example Output

### metadata

```
E-MTAB-10030
Title:        Single-cell RNA-seq of rat microglia in monoculture and in
              coculture with neurons and astrocytes
Released:     2021-03-23
Collection:   ArrayExpress

Files: 26 total
  idf              1
  sdrf             1
  raw             24
```

### files

```
26 file(s) for E-MTAB-10030:

idf                   3.1 kB  E-MTAB-10030.idf.txt
sdrf                  7.2 kB  E-MTAB-10030.sdrf.txt
raw                  1.4 GB  B1_S1_R1.fastq.gz
raw                  1.5 GB  B1_S1_R2.fastq.gz
```

### samplesheet (--assay scrna)

```csv
sample,fastq_1,fastq_2
Sample_1,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S1_R2.fastq.gz
Sample_2,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R1.fastq.gz,https://ftp.ebi.ac.uk/pub/databases/microarray/data/experiment/MTAB/E-MTAB-10030/B1_S2_R2.fastq.gz
```

### tables/metadata.tsv

```
sample     replicate  species            sex    age      condition  tissue
Sample 2   1          Rattus norvegicus  mixed  0 to 2   NA         neocortex
Sample 4   1          Rattus norvegicus  mixed  0 to 2   NA         neocortex
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── samplesheet.csv
├── tables/
│   └── metadata.tsv
├── downloads/                       # (optional) --command download only
├── download_arrayexpress.sh         # (optional) --command download-script only
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

- **Python** ≥ 3.10, standard library only. No pip install required.
- **`curl` or `wget`** — only for a generated `download-script`, not for the
  skill itself.
- **Network access** to `www.ebi.ac.uk` and `ftp.ebi.ac.uk` on TCP/443. See
  Gotcha 5.

## Gotchas

- **Gotcha 1**: You will want to trust the R1/R2 filename heuristic for a 10x
  run. **Do not.** For a Chromium run whose technical reads are separate files,
  the heuristic drops the cDNA read and produces a plausible-looking samplesheet
  that quantifies the wrong thing. Pass `--read-map` explicitly: 3 files → `2,3`,
  4 files (dual index) → `3,4`. Confirm the choice with the user before writing
  the sheet. This is the highest-risk behaviour in the skill.
- **Gotcha 2**: You will want to quote the file `size` as megabytes. It is
  **bytes**, and ArrayExpress raw data is routinely tens of gigabytes. Never
  start `--command download --raw` without telling the user the total first.
- **Gotcha 3**: `--read-map` accepts positions **1–9**, not 1–4. This is upstream
  behaviour, kept deliberately. A typo like `1,9` passes validation and fails
  later with a missing-file error rather than being rejected up front.
- **Gotcha 4**: You will want to hardcode
  `https://www.ebi.ac.uk/biostudies/files/{accession}/{path}`. Do not. It works —
  it 302-redirects — but there is no one base tree behind it, and the redirect
  costs ~40× the latency and times out on large files. This skill resolves
  `httpLink` from `/studies/{accession}/info` and keeps the old path as a
  fallback.
- **Gotcha 5**: If `--command metadata` works but `--command download` hangs or
  fails with a connection timeout, suspect a **firewall, not a bug**. Metadata
  comes from `www.ebi.ac.uk`; file bytes come from `ftp.ebi.ac.uk`. Networks
  routinely allow the first and block the second. Both need allowlisting on
  **TCP/443** — this skill never speaks the FTP protocol despite the hostname, so
  asking an admin to "open FTP" achieves nothing. See
  [docs/data-handling.md](../../docs/data-handling.md#allowlisting-for-the-public-archive-skills).
- **Gotcha 6**: A sample-name collision is **fatal**, by design. If two source
  names normalise to the same id the skill exits rather than emitting a sheet, because
  the pipeline concatenates rows sharing a `sample` and would pool distinct
  biological samples without saying so.
- **Gotcha 7**: A `403` from EBI mid-download is usually **rate limiting, not
  a permissions problem** — it appears after a burst of requests to the same
  host and clears on retry (observed 2026-09-22). Generated download scripts
  handle this themselves: `--tool curl` emits `--retry 5 --retry-delay 10` plus
  `--retry-all-errors` where curl supports it. In-process `--command download`
  does not retry on 403, so just re-run it.
- **Gotcha 8**: Upstream's `--out` defaults were relative to the working
  directory, so `samplesheet.csv` landed wherever you happened to be. Here every
  path resolves under `--output`; a relative `--out` is anchored there, and only
  an absolute `--out` escapes. Do not reintroduce cwd-relative defaults.
- **Gotcha 9**: You will want to reach for `--command download-script` to get the
  FASTQs. There is none, deliberately. **ArrayExpress brokers sequencing reads to
  ENA** and serves the bytes from there, so build `samplesheet.csv` here and emit
  the script with `ena-fetch --command download-script` against the same
  `--output`; it fetches whatever URLs the sheet names, whether they point at
  `ftp.sra.ebi.ac.uk` or `ftp.ebi.ac.uk`. For runs not mirrored to ENA, or for
  10x reads whose structure only `fasterq-dump` exposes reliably, use sra-tools
  (`prefetch --option-file`, then `fasterq-dump`). `--command download` is a
  different thing and still works for files ArrayExpress *does* host — IDF, SDRF,
  processed matrices, CEL, BAM.
- **Gotcha 10**: `tables/metadata.tsv` is one row per sample × replicate, **not
  one per SDRF line**. A bulk paired-end SDRF puts each FASTQ on its own line, so
  a 12-sample study has 24 of them; the rows are grouped on `Comment[ENA_RUN]`,
  falling back to `Source Name`, before they are harmonised. Do not "fix" a row
  count that looks low by iterating the SDRF directly — that is the bug this
  replaced, and it also doubled the BioSamples lookups.
- **Gotcha 11**: `--strandedness` defaults to `auto` and should stay there unless
  the record states the library chemistry. The value follows from the chemistry,
  never from a library merely being "stranded": dUTP second-strand marking
  (TruSeq Stranded mRNA, NEBNext Ultra II Directional) → `reverse`; Lexogen
  QuantSeq 3′ **FWD** → `forward` — one vendor, both directions; non-directional
  kits → `unstranded`. When you do know it, declare it: nf-core/rnaseq infers per
  sample either way, but only an explicit value earns a **mismatch report**, and
  `auto` has nothing to compare against. A wrong explicit value is worse than
  `auto` — it mislabels every sample and suppresses nothing.

## Safety

- **Local-first**: no user genetic data is ever transmitted. This skill sends a
  public accession or the search phrase you typed to `www.ebi.ac.uk` and
  receives public archive data. Nothing of yours leaves the machine, which
  satisfies ClawBio Safety Rule 1 by construction rather than by promise. File
  bytes are fetched from `ftp.ebi.ac.uk` over HTTPS. See
  [docs/data-handling.md](../../docs/data-handling.md).
- **Credentials**: none. ArrayExpress needs no key, and the skill reads no
  credential environment variables.
- **Execution**: `--run` and `--submit` are off by default and must be confirmed
  with the user each time.
- **Disclaimer**: every report carries the ClawBio medical disclaimer —
  *"ClawBio is a research and educational tool. It is not a medical device and
  does not provide clinical diagnoses. Consult a healthcare professional before
  making any medical decisions."*

## Agent Boundary

**The agent dispatches and explains; the skill executes.** The agent chooses the
accession and the subcommand, decides whether a 10x run needs `--read-map`, asks
the user before any download or script execution, and explains the result. It
does **not** re-derive FASTQ URLs, hand-edit the samplesheet, invent
harmonisation mappings, or override a sample-name collision. Every number in the
report comes from the script's output, never from the model.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on an ArrayExpress
accession (`E-MTAB`, `E-GEOD`, `E-MEXP`, `E-PROT`) or an explicit mention of
ArrayExpress, MAGE-TAB, SDRF or IDF.

**Chaining partners**:
- `biostudies-fetch`: the same API for non-ArrayExpress collections; route
  `S-BSST*`/`S-BIAD*` there.
- `ena-fetch`: when the study's runs are wanted from the sequence archives
  rather than the ArrayExpress FASTQ mirror.
- `nfcore-rnaseq-wrapper` / `nfcore-scrnaseq-wrapper`: the natural consumers of
  the `samplesheet.csv` this skill writes.
- `rnaseq-de` / `scrna-orchestrator`: downstream analysis once counts exist.
- `article-data-fetcher`: **upstream producer.** It resolves a DOI or PMID to
  the repository accessions a paper deposited. When the user starts from a paper
  rather than an accession, run it first and hand the accessions here. It
  downloads files and writes a `manifest.json`, but it does **not** parse
  MAGE-TAB, harmonise sample annotation into `metadata.tsv`, or emit a
  pipeline-ready `samplesheet.csv` — that is this skill's job, so the two chain
  rather than compete.

## Maintenance

- **Update trigger**: when an archive host changes an endpoint this skill calls.
  Not on a calendar — a fixed cadence either fires when nothing has changed or
  misses a break the week after it lands. The staleness signals below are the
  trigger.
- **How a break surfaces**: from a live call, not from CI. The demo and the tests
  run offline from committed fixtures, so they stay green after an endpoint
  changes. Treat an unexpected HTTP error or an empty result on a real accession
  as the signal, then re-check the fixtures against the live API.
- **Staleness signals**: `/studies/{acc}/info` stops advertising `httpLink`; the
  MAGE-TAB column vocabulary changes; nf-core changes its samplesheet columns.
- **Known debt**: the metadata-table, sanitiser and read-map helper blocks are
  duplicated across the five archive scripts rather than factored into a shared
  module. This is deliberate — it keeps re-syncing with upstream cheap. Revisit
  only if upstream itself deduplicates.
- **Deprecation criteria**: retire if EMBL-EBI folds ArrayExpress fully into the
  generic BioStudies record shape, at which point `biostudies-fetch` would cover
  it.

## Citations

- **Ported from**
  [`UKDRI/informatics_data_skills`](https://github.com/UKDRI/informatics_data_skills)
  at commit `7cc3e6e` (`arrayexpress/scripts/arrayexpress.py`), MIT licensed.
  Upstream copyright retained: Copyright (c) 2026 UK Dementia Research Institute.
- **ArrayExpress / BioStudies**: Sarkans et al. (2018) *The BioStudies database —
  one stop shop for all data supporting a life sciences study.* Nucleic Acids
  Research 46(D1):D1266–D1270. <https://doi.org/10.1093/nar/gkx965>
- **Demo fixture**: E-MTAB-10030, "Single-cell RNA-seq of rat microglia in
  monoculture and in coculture with neurons and astrocytes", released
  2021-03-23. The record declares **no explicit licence attribute**; it is
  redistributed here as public archive metadata on the terms EMBL-EBI publishes
  it under. *Rattus norvegicus* throughout, so no human individual-level data.
- **MAGE-TAB**: Rayner et al. (2006) *A simple spreadsheet-based, MIAME-supportive
  format for microarray data: MAGE-TAB.* BMC Bioinformatics 7:489.
  <https://doi.org/10.1186/1471-2105-7-489>
