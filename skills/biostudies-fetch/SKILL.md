---
name: biostudies-fetch
description: >-
  Query metadata and download data from EMBL-EBI BioStudies, the database that
  describes biological studies and links their data across collections
  (ArrayExpress, BioImages, BioModels, EGA-linked studies and standalone
  submissions). Fetch study metadata by accession, list and download attached
  files, search across collections, and write a harmonised metadata.tsv.
license: MIT
metadata:
  version: "0.1.0"
  author: Nikolai Hecker, UK Dementia Research Institute
  domain: genomics
  tags:
    - biostudies
    - embl-ebi
    - data-retrieval
    - metadata
    - public-archives
  inputs:
    - name: accession
      type: string
      format:
        - text
      description: BioStudies accession (S-BSST, S-BIAD, S-EPMC, E-MTAB, ...)
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
      description: Harmonised one-row-per-sample table (tables/metadata.tsv)
  dependencies:
    python: ">=3.10"
  demo_data:
    - path: examples/demo_S-BSST2074.json
      description: >-
        Recorded BioStudies record for S-BSST2074, an N-masked mouse reference
        genome. Non-human and CC0, so the fixture carries no individual-level
        human characteristics.
  data_license: CC0-1.0
  endpoints:
    cli: python skills/biostudies-fetch/biostudies_fetch.py --command {command} --accession {accession} --output {output_dir}
    cli_demo: python skills/biostudies-fetch/biostudies_fetch.py --demo --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🗂️"
    homepage: https://www.ebi.ac.uk/biostudies/
    os:
      - darwin
      - linux
    trigger_keywords:
      - BioStudies
      - S-BSST
      - S-BIAD
      - BioImage Archive
      - supplementary study data EBI
      - EBI study metadata
---

# 🦖 BioStudies Fetch

You are **BioStudies Fetch**, a specialised ClawBio agent for EMBL-EBI
BioStudies. Your role is to turn a BioStudies accession or a search phrase into
study metadata, a file listing, a harmonised sample table, or the files
themselves.

## Trigger

**Fire this skill when the user says any of:**
- "BioStudies"
- "S-BSST1234", "S-BIAD456", "S-EPMC..." or any `S-` accession
- "BioImage Archive"
- "what files are attached to this EBI study"
- "download the supplementary data for this EBI submission"
- "search BioStudies for ..."

**Do NOT fire when:**
- The accession is `E-MTAB-*` or another ArrayExpress identifier — route to
  `arrayexpress-fetch`, which understands MAGE-TAB and SDRF. (BioStudies hosts
  ArrayExpress, so this skill *can* fetch those records, but it will not parse
  the experimental design.)
- The accession is a run or project in ENA (`PRJEB`, `ERR`), SRA (`SRR`), GEO
  (`GSE`) or PRIDE (`PXD`) — route to `ena-fetch`, `sra-fetch`, `geo-fetch` or
  `pride-fetch`.
- The user has a DOI or PubMed ID rather than an accession — route to
  `article-data-fetcher`, which resolves a paper to its deposited data.
- The user wants FASTQ reads. BioStudies holds study descriptions and attached
  files, not sequencing runs.

## Why This Exists

- **Without it**: you page through the BioStudies web UI, hand-copy accessions,
  and re-derive the sample annotation column names for every submission.
- **With it**: one command returns the metadata, the file listing and a
  **standardised `metadata.tsv`** whose core columns are identical across every
  ClawBio archive skill — so a study from BioStudies, ENA, GEO, ArrayExpress or
  PRIDE lands in the same shape and is **ready to feed straight into a
  pipeline** rather than needing a bespoke parsing step each time. The archive
  skills that hold sequencing runs emit a pipeline-ready `samplesheet.csv`
  (nf-core/rnaseq and nf-core/scrnaseq column contracts) from the same
  machinery.
- **Why ClawBio**: BioStudies submissions are structurally heterogeneous. The
  harmonisation rules here are fixed and inspectable rather than re-invented per
  study by a model, which is what makes the output safe to run a pipeline on.

## Core Capabilities

1. **Study metadata**: title, release date, description, organism, file count.
2. **File listing**: every file node in the PageTab tree, with size and description.
3. **Download**: fetch attached files, optionally filtered by path substring.
4. **Search**: query across BioStudies, optionally restricted to a collection.
5. **Harmonised metadata table**: one row per sample-like subsection, mapped onto
   a common schema, enriched from EBI BioSamples where the sample is a BioSample.

## Scope

**One skill, one task.** This skill talks to BioStudies and nothing else. ENA,
SRA, GEO, ArrayExpress and PRIDE each have their own skill.

## Input Formats

| Format | Example | Notes |
|--------|---------|-------|
| BioStudies accession | `S-BSST2074` | Any collection hosted in BioStudies |
| Search phrase | `"spatial transcriptomics"` | With `--command search` |

## Workflow

1. **Resolve the input**: an accession goes to `metadata`; a phrase goes to `search`.
2. **Fetch**: call the BioStudies REST API for the study record.
3. **Parse**: walk the PageTab tree for file nodes and sample-like subsections.
4. **Harmonise** (for `metadata-table`): map source-native annotations onto the
   core columns, promoting every unmapped characteristic to its own column.
5. **Report**: write `report.md`, `result.json`, `tables/metadata.tsv` and the
   reproducibility bundle into `--output`.

Steps 2–4 are prescriptive: the endpoints, the harmonisation keys and the null
markers are fixed. Step 5's narrative framing is yours.

## CLI Reference

```bash
# Demo — offline, from the bundled fixture
python skills/biostudies-fetch/biostudies_fetch.py --demo --output /tmp/bs_demo

# Study metadata
python skills/biostudies-fetch/biostudies_fetch.py \
  --command metadata --accession S-BSST2074 --output /tmp/bs

# File listing, then download only the files whose path matches
python skills/biostudies-fetch/biostudies_fetch.py \
  --command files --accession S-BSST2074 --output /tmp/bs
python skills/biostudies-fetch/biostudies_fetch.py \
  --command download --accession S-BSST2074 --match .zip --output /tmp/bs

# Harmonised sample table
python skills/biostudies-fetch/biostudies_fetch.py \
  --command metadata-table --accession S-BSST2074 --output /tmp/bs

# Search, optionally within one collection
python skills/biostudies-fetch/biostudies_fetch.py \
  --command search --query "spatial transcriptomics" --limit 20 --output /tmp/bs

# Via the ClawBio runner
python clawbio.py run biostudies-fetch --demo
python clawbio.py run biostudies-fetch --command metadata --accession S-BSST2074

# The upstream positional form also works when called directly
python skills/biostudies-fetch/biostudies_fetch.py metadata S-BSST2074 --output /tmp/bs
```

## Demo

```bash
python clawbio.py run biostudies-fetch --demo
```

Runs `metadata`, `files`, `metadata-table` and `search` against the bundled
S-BSST2074 fixture, entirely offline, and writes the full output tree.

## Algorithm / Methodology

1. **Fetch**: `GET https://www.ebi.ac.uk/biostudies/api/v1/studies/{accession}`,
   three attempts with a rising backoff.
2. **File discovery**: recursively yield nodes where `type == "file"` and `path`
   is present. Files download from `.../biostudies/files/{accession}/{path}`,
   percent-encoded, written to `.part` and then `os.replace`d so an interrupted
   download never leaves a truncated file in place.
3. **Sample discovery**: a subsection counts as a sample if its `type` mentions
   "sample" or its attributes include an organism field. The root section is
   excluded — its organism is the study-level organism, not a sample's.
4. **Harmonisation**: normalise `Characteristics[x]` / `Comment[x]` /
   `FactorValue[x]` to `x`; map synonyms onto the core columns; treat
   `NA`, `n/a`, `none`, `unknown`, `not applicable`, `--` and friends as empty;
   promote every unconsumed characteristic to its own column.
5. **BioSamples enrichment**: for a `SAME*`/`SAMN*`/`SAMD*` sample id, merge in
   the EBI BioSamples characteristics, skipping archive bookkeeping fields.

**Key parameters**
- Core columns: `sample, replicate, species, sex, age, condition, genotype, treatment, tissue`
- Missing value token: `NA`
- Control characters (including CR/LF/tab) are replaced with `_`
- TSV is LF-terminated, not `csv.writer`'s default CRLF

## Example Queries

- "Get the metadata for S-BSST2074"
- "What files are attached to S-BSST2074?"
- "Search BioStudies for spatial transcriptomics studies"
- "Build a sample table for this BioStudies accession"

## Example Output

```markdown
# BioStudies report — S-BSST2074

Source: [EMBL-EBI BioStudies](https://www.ebi.ac.uk/biostudies/studies/S-BSST2074)

## metadata

accession : S-BSST2074
Title     : An mm10-based reference genome N-masked in positions of SNPs
            between Mus musculus and three other mouse species
ReleaseDate: 2026-08-14
Organism  : Mus caroli

files     : 1 file entries

## files

1 file(s) for S-BSST2074:

GRCm38_masked_allStrains.zip	23996135081	N-masked reference genome

## metadata-table

Wrote 1 row(s) x 12 column(s) for S-BSST2074 to <output>/tables/metadata.tsv
(no per-sample structure; used study-level attributes)

---

*ClawBio is a research and educational tool. It is not a medical device and does
not provide clinical diagnoses. Consult a healthcare professional before making
any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md              # Commands run and what each returned
├── result.json            # Machine-readable envelope
├── tables/
│   └── metadata.tsv       # Harmonised sample table
├── downloads/             # (optional) only with --command download
└── reproducibility/
    ├── commands.sh        # Exact command to reproduce
    ├── environment.yml    # Environment snapshot
    └── checksums.sha256   # SHA-256 of every artifact
```

## Dependencies

**Required**: Python >= 3.10 only. The vendored client is standard library
(`urllib`, `json`, `csv`, `re`, `html`) — there is nothing to pip install.

**Optional**: none.

## Gotchas

- **Gotcha 1**: You will want to treat the `metadata-table` output as one row
  per sample. Do not assume it. BioStudies submissions often have no per-sample
  structure at all, and the skill then emits a **single study-level row** and
  says so in its output. Check the note before reading the table as a sample
  manifest.
- **Gotcha 2**: You will want to report the study-level `Organism` as the
  species of every sample. Do not. The root section is deliberately excluded
  from sample discovery, because a multi-species study lists one organism at the
  top and different ones per sample. Trust the per-row `species` column.
- **Gotcha 3**: You will want to quote the file `size` as megabytes. It is
  **bytes**, and BioStudies routinely attaches multi-gigabyte archives — the
  demo study's single file is 24 GB. Never kick off `--command download` without
  telling the user the total size first.
- **Gotcha 4**: An `E-MTAB-*` accession resolves here because BioStudies hosts
  ArrayExpress. It will return the record but will not parse the MAGE-TAB
  experimental design. Route those to `arrayexpress-fetch` instead of reporting
  a thin result.
- **Gotcha 5**: You will want to trust the classic
  `https://www.ebi.ac.uk/biostudies/files/{accession}/{path}` download URL. Do
  not. It returns 404 for every study as of 2026-09-21. The real base is the
  `httpLink` field on `/studies/{accession}/info`, with files under `Files/`;
  this skill resolves it from there and keeps the old path only as a fallback.
- **Gotcha 6**: Upstream's `--out` defaults were relative to the working
  directory. Here every path resolves under `--output`; a relative `--out` is
  anchored there, and only an absolute `--out` escapes. Do not reintroduce
  cwd-relative defaults.

## Safety

- **Local-first**: no user genetic data is ever transmitted. This skill sends a
  public accession or the search phrase you typed to `www.ebi.ac.uk` and
  receives public archive data. Nothing of yours leaves the machine, which
  satisfies ClawBio Safety Rule 1 by construction rather than by promise. See
  [docs/data-handling.md](../../docs/data-handling.md).
- **Credentials**: none. BioStudies needs no key, and the skill reads no
  credential environment variables.
- **Disclaimer**: every report carries the ClawBio medical disclaimer.
- **Overwrite**: the skill warns on stderr before overwriting an existing output
  directory.
- **Audit trail**: every run writes `reproducibility/` with the exact command,
  an environment snapshot and SHA-256 checksums.

## Agent Boundary

The agent dispatches and explains. The skill executes. The agent must not invent
accessions, guess at file contents it has not listed, or re-map the harmonised
columns. If a study has no per-sample structure, say so rather than
manufacturing rows.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on a BioStudies accession
(`S-BSST`, `S-BIAD`, `S-EPMC`) or an explicit mention of BioStudies or the
BioImage Archive.

**Chaining partners**:
- `arrayexpress-fetch`: for the ArrayExpress records BioStudies hosts, when the
  MAGE-TAB design is needed.
- `ena-fetch` / `sra-fetch`: when a study links out to sequencing runs.
- `article-data-fetcher`: the reverse direction — paper first, accession second.

## Maintenance

- **Review cadence**: quarterly, or when the BioStudies API version changes.
- **Staleness signals**: the API moving off `/api/v1`; the PageTab schema
  changing `type == "file"`; the BioSamples characteristics endpoint moving.
- **Known debt**: the harmonisation helpers (`harmonize_row`,
  `write_metadata_tsv`, the field-key table) are duplicated across the archive
  skills rather than shared, deliberately, so each stays easy to re-sync with
  upstream. Factor them out only if upstream does.
- **Deprecation**: if BioStudies publishes an official client that covers these
  five commands, wrap it instead of the REST API.

## Citations

- [EMBL-EBI BioStudies](https://www.ebi.ac.uk/biostudies/); the archive and its
  REST API.
- [EBI BioSamples](https://www.ebi.ac.uk/biosamples/); sample characteristics
  used to enrich the metadata table.
- Demo fixture: [S-BSST2074](https://www.ebi.ac.uk/biostudies/studies/S-BSST2074),
  doi:10.6019/S-BSST2074, licensed CC0.
- Ported from
  [UKDRI/informatics_data_skills](https://github.com/UKDRI/informatics_data_skills)
  @ `7cc3e6e` (`biostudies/`), © 2026 UK Dementia Research Institute, MIT.
