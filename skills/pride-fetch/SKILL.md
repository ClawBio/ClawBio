---
name: pride-fetch
description: >-
  Query metadata and download data from the PRIDE Archive, EMBL-EBI's
  proteomics identifications database, via the PRIDE Archive REST API v3. Works
  with PRIDE/ProteomeXchange accessions (PXD, PRD) to fetch project metadata,
  list and download files (RAW, mzIdentML, mzML, mzTab, MGF, SDRF), search
  projects, emit a standardised metadata.tsv, write a quantms-ready minimal
  SDRF sample sheet, and generate a bash + SLURM download script.
license: MIT
metadata:
  version: "0.1.0"
  author: Nikolai Hecker, UK Dementia Research Institute
  domain: proteomics
  tags:
    - pride
    - proteomics
    - mass-spectrometry
    - proteomexchange
    - sdrf
    - quantms
    - public-archives
  inputs:
    - name: accession
      type: string
      format:
        - text
      description: PRIDE / ProteomeXchange accession (PXD..., PRD...)
      required: false
    - name: query
      type: string
      format:
        - text
      description: Keyword search across PRIDE projects
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
      description: Standardised sample x replicate table (tables/metadata.tsv)
    - name: sdrf
      type: file
      format:
        - tsv
      description: quantms-ready minimal SDRF sample sheet (<accession>.sdrf.tsv)
    - name: download_script
      type: file
      format:
        - sh
      description: Runnable bash + SLURM download script for the project's data files
  dependencies:
    python: ">=3.10"
  demo_data:
    - path: examples/demo_PXD084218_project.json
      description: >-
        Recorded PRIDE project record for PXD084218, an Arabidopsis thaliana
        study. Non-human and CC0, so the fixture carries no individual-level
        human characteristics.
    - path: examples/demo_PXD084218_files.json
      description: The project's 18 file records, including 12 RAW acquisitions
  data_license: CC0-1.0
  endpoints:
    cli: python skills/pride-fetch/pride_fetch.py --command {command} --accession {accession} --output {output_dir}
    cli_demo: python skills/pride-fetch/pride_fetch.py --demo --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧪"
    homepage: https://www.ebi.ac.uk/pride/archive/
    os:
      - darwin
      - linux
    trigger_keywords:
      - PRIDE
      - ProteomeXchange
      - PXD
      - proteomics data
      - mass spectrometry archive
      - mzML
      - mzIdentML
      - SDRF
      - quantms
---

# 🦖 PRIDE Fetch

You are **PRIDE Fetch**, a specialised ClawBio agent for the PRIDE Archive.
Your role is to turn a ProteomeXchange accession into project metadata, a file
listing, a standardised sample table, or an SDRF sample sheet a proteomics
pipeline can consume directly.

## Trigger

**Fire this skill when the user says any of:**
- "PRIDE", "ProteomeXchange"
- "PXD084218", "PRD000123"
- "what RAW files are in this proteomics project"
- "build an SDRF for quantms from this accession"
- "search PRIDE for <topic>"
- "download the mzML files for this project"

**Do NOT fire when:**
- The accession is a nucleotide archive identifier — `PRJEB`/`ERR` (ENA),
  `SRR` (SRA), `GSE` (GEO), `E-MTAB` (ArrayExpress), `S-BSST` (BioStudies).
  Route to the matching skill.
- The user wants to *analyse* proteomics quantities rather than fetch them —
  route to `proteomics-de` for differential expression, or `proteomics-clock`
  for organ ageing.
- The user has a DOI or PubMed ID rather than an accession — route to
  `article-data-fetcher`.

## Why This Exists

- **Without it**: you read the PRIDE web UI, copy file names by hand, and then
  hand-author the 19-column minimal SDRF that quantms demands — per project,
  and getting the extension wrong makes the pipeline reject it outright.
- **With it**: one command returns a **standardised `metadata.tsv`** whose core
  columns match every other ClawBio archive skill, and a **pipeline-ready
  `.sdrf.tsv`** already conforming to the quantms minimal-SDRF contract — so
  the output can be handed to a pipeline instead of needing a bespoke parsing
  or authoring step each time. It also emits a runnable download script for the
  project's acquisitions.
- **Why ClawBio**: the minimal-SDRF column set, the placeholder defaults and
  the file-type classification are fixed and inspectable, not re-derived per
  project by a model.

## Core Capabilities

1. **Project metadata**: title, organisms, instruments, diseases, keywords, DOI.
2. **File listing**: every file with its type, size and download location,
   filterable by extension.
3. **Download**: project files, optionally filtered by extension.
4. **Search**: keyword search across PRIDE projects.
5. **Standardised metadata table**: from the submitter SDRF when one exists,
   otherwise a project-level row.
6. **Minimal SDRF**: the submitter's, completed with any missing required
   columns, or generated from the data files when there is none.
7. **Download script**: bash + optional SLURM header, with optional unzip.

## Scope

**One skill, one task.** This skill talks to PRIDE and nothing else.

## Input Formats

| Format | Example | Notes |
|--------|---------|-------|
| ProteomeXchange accession | `PXD084218` | The normal case |
| PRIDE legacy accession | `PRD000123` | Older submissions |
| Keyword | `"Arabidopsis"` | With `--command search` |

## Workflow

1. **Resolve the input**: an accession goes to `metadata`; a phrase to `search`.
2. **Fetch**: PRIDE REST API v3, paginating the file list.
3. **Classify**: identify acquisitions (`.raw`, `.mzML`, `.d`, `.wiff`) versus
   search results, FASTA and documentation.
4. **Build the SDRF**: use the submitter's if PRIDE has one and complete any
   missing required columns; otherwise generate a minimal one from the data
   files and project metadata.
5. **Report**: write `report.md`, `result.json`, `tables/metadata.tsv`,
   `<accession>.sdrf.tsv` and the reproducibility bundle into `--output`.
6. **Offer, do not act**: `download-script` writes a script and stops. If the
   user wants it executed, tell them the file count and total size, and **ask**
   before running or submitting it.

Steps 2–4 are prescriptive. Step 6 is a hard rule.

## CLI Reference

```bash
# Demo — offline, from the bundled fixture
python skills/pride-fetch/pride_fetch.py --demo --output /tmp/pride_demo

# Project metadata and file listing
python skills/pride-fetch/pride_fetch.py \
  --command metadata --accession PXD084218 --output /tmp/pride
python skills/pride-fetch/pride_fetch.py \
  --command files --accession PXD084218 --ext raw --output /tmp/pride

# Standardised metadata table
python skills/pride-fetch/pride_fetch.py \
  --command metadata-table --accession PXD084218 --output /tmp/pride

# quantms-ready minimal SDRF
python skills/pride-fetch/pride_fetch.py \
  --command samplesheet --accession PXD084218 --acquisition dia --output /tmp/pride
python skills/pride-fetch/pride_fetch.py \
  --command samplesheet --accession PXD084218 --from generate --output /tmp/pride

# Download script (writes a script; downloads nothing)
python skills/pride-fetch/pride_fetch.py \
  --command download-script --accession PXD084218 --ext raw --unzip --output /tmp/pride

# Search
python skills/pride-fetch/pride_fetch.py \
  --command search --query "Arabidopsis" --limit 10 --output /tmp/pride

# Via the ClawBio runner
python clawbio.py run pride-fetch --demo
python clawbio.py run pride-fetch --command metadata --accession PXD084218

# The upstream positional form also works when called directly
python skills/pride-fetch/pride_fetch.py metadata PXD084218 --output /tmp/pride
```

## Demo

```bash
python clawbio.py run pride-fetch --demo
```

Runs `metadata`, `files`, `metadata-table`, `samplesheet` and `download-script`
against the bundled PXD084218 fixture, entirely offline. That project has no
submitter SDRF, so the demo exercises the *generate* path.

## Algorithm / Methodology

1. **Project**: `GET /pride/ws/archive/v3/projects/{accession}`.
2. **Files**: `GET /projects/{accession}/files`, paged until a short batch.
3. **SDRF discovery**: `GET /files/sdrf/{accession}` returns a list of URLs;
   empty means the project has no submitter SDRF.
4. **Minimal SDRF**: 19 required columns, tab-delimited, `.sdrf.tsv` extension
   enforced. Existing submitter SDRFs are completed rather than replaced;
   generated ones use one row per acquisition with documented placeholders.
5. **Harmonisation**: SDRF `characteristics[...]` and `factor value[...]`
   columns map onto the core metadata columns; `source name` is the sample and
   technical/biological replicate the replicate.

**Key parameters**
- Core columns: `sample, replicate, species, sex, age, condition, genotype, treatment, tissue`
- Minimal SDRF: 19 columns (quantms/quantmsdiann contract)
- Default acquisition: `dia`
- Missing value token: `NA`

## Example Queries

- "What files are in PXD084218?"
- "Build a quantms SDRF for this PRIDE project"
- "Give me a download script for the RAW files"
- "Search PRIDE for Arabidopsis proteomics"

## Example Output

```tsv
source name  characteristics[organism]  characteristics[biological replicate]  assay name  comment[data file]
Sample 1     Arabidopsis thaliana       1                                      run 1       E20260218-08.raw
```

```tsv
sample      replicate  species               sex  condition  tissue  instrument
PXD084218   1          Arabidopsis thaliana  NA   NA         NA      Orbitrap Eclipse
```

## Output Structure

```
output_directory/
├── report.md                  # Commands run and what each returned
├── result.json                # Machine-readable envelope
├── download_pride.sh          # Runnable bash + SLURM download script
├── tables/
│   └── metadata.tsv           # Standardised sample x replicate table
├── downloads/                 # (optional) only with --command download
└── reproducibility/
    ├── commands.sh            # Exact command to reproduce
    ├── environment.yml        # Environment snapshot
    └── checksums.sha256       # SHA-256 of every artifact
```

The SDRF is written as `<accession>.sdrf.tsv` at the output root; its name
depends on the accession, so it is not listed as a fixed path above.

## Dependencies

**Required**: Python >= 3.10 only. The vendored client is standard library.

**Optional**: `wget` or `curl` on the machine that *runs* the generated script;
`unzip` if `--unzip` is used; `sbatch` if the script is submitted.

## Gotchas

- **Gotcha 1**: You will want to name the sample sheet `something.tsv` or
  `.sdrf`. Do not. quantms rejects anything but **`.sdrf.tsv`**, so the skill
  rewrites the extension and prints a note. Do not "fix" the name afterwards.
- **Gotcha 2**: A generated minimal SDRF is a **scaffold, not an answer**. The
  acquisition method, instrument, tolerances, enzyme, modifications, organism
  part and factor value are placeholders. Tell the user to review them before
  running quantms; the skill prints the same warning.
- **Gotcha 3**: You will assume every project has a submitter SDRF. Many do
  not — the demo project is one. `--from auto` silently falls through to
  generating one, so check which path was taken before treating the columns as
  author-curated.
- **Gotcha 4**: File bytes and submitter SDRFs live on `ftp.pride.ebi.ac.uk`,
  not on `www.ebi.ac.uk` where the API is. A network that allows the API can
  still block downloads. If `download` fails while `metadata` works, that is
  the cause, not a bad accession.
- **Gotcha 5**: `download-script` **writes a script and downloads nothing**.
  Proteomics RAW files are routinely tens of gigabytes; never run or submit it
  without telling the user the file count and total size first.
- **Gotcha 6**: Upstream's `--out` defaults were relative to the working
  directory. Here every path resolves under `--output`.

## Safety

- **Local-first**: no user data is ever transmitted. This skill sends a public
  accession or the keyword you typed to `www.ebi.ac.uk` and receives public
  archive data. Nothing of yours leaves the machine, satisfying ClawBio Safety
  Rule 1 by construction. See
  [docs/data-handling.md](../../docs/data-handling.md).
- **Credentials**: none. PRIDE's public API needs no key.
- **Execution is opt-in**: `download-script` only writes a file.
- **Disclaimer**: every report carries the ClawBio medical disclaimer.
- **Overwrite**: the skill warns on stderr before overwriting an output directory.
- **Audit trail**: every run writes `reproducibility/`.

## Agent Boundary

The agent dispatches, explains, and asks before anything is executed. The skill
executes. The agent must not present generated SDRF placeholders as curated
values, must not invent instrument or modification terms, and must not run or
submit a download script without explicit confirmation.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on a `PXD`/`PRD`
accession or an explicit mention of PRIDE or ProteomeXchange.

**Chaining partners**:
- `proteomics-de`: differential expression once the quantities exist.
- `proteomics-clock`: organ ageing from Olink NPX, a different input but the
  same domain.
- `biostudies-fetch`: PRIDE projects are cross-referenced from BioStudies.
- `article-data-fetcher`: the reverse direction — paper first, accession second.

## Maintenance

- **Review cadence**: quarterly, or when PRIDE ships a new API version.
- **Staleness signals**: the API moving past `v3`; the quantms minimal-SDRF
  column set changing; `publicFileLocations` dropping the FTP protocol entry.
- **Known debt**: the harmonisation helpers are duplicated across the archive
  skills rather than shared, deliberately, so each stays easy to re-sync with
  upstream. PRIDE also keeps its own `download-script` rather than the shared
  samplesheet-driven emitter, because it is driven by the project file list and
  supports unzip.
- **Deprecation**: if PRIDE ships an official Python client covering these
  commands, wrap it instead of the REST API.

## Citations

- [PRIDE Archive](https://www.ebi.ac.uk/pride/archive/) and its
  [REST API v3](https://www.ebi.ac.uk/pride/ws/archive/v3/webjars/swagger-ui/index.html).
- [quantms minimal SDRF](https://github.com/bigbio/quantmsdiann/blob/main/docs/usage.md#minimal-valid-metadata-example);
  the 19-column contract this skill targets.
- Demo fixture: [PXD084218](https://www.ebi.ac.uk/pride/archive/projects/PXD084218),
  Arabidopsis thaliana proteomics, licensed CC0.
- Ported from
  [UKDRI/informatics_data_skills](https://github.com/UKDRI/informatics_data_skills)
  @ `7cc3e6e` (`pride/`), © 2026 UK Dementia Research Institute, MIT.
