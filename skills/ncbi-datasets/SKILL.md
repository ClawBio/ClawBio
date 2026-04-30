---
name: ncbi-datasets
description: >-
  Download genomes, genes, virus sequences, and taxonomy data from NCBI using
  the datasets and dataformat CLI tools.
version: 0.1.0
author: nullvoid42
domain: datasets
license: MIT

tags: [
  ncbi,
  genomics,
  bioinformatics,
  genome-download,
  gene,
  virus,
  taxonomy,
  datasets,
  dataformat,
  refseq,
  genbank
]

metadata:
  openclaw:
    requires:
      bins:
        - datasets
        - dataformat
      env: []        
      config: []
    always: false
    homepage: https://www.ncbi.nlm.nih.gov/datasets/
    os: [macos, linux, windows]
    install:
      - kind: conda
        package: ncbi-datasets-cli
        channel: conda-forge
        bins: [datasets, dataformat]
    trigger_keywords:
      - ncbi
      - datasets download
      - datasets summary
      - dataformat
      - download genome
      - download gene
      - reference genome
      - genome accession
      - gene symbol
      - ortholog
      - download virus
      - sars-cov-2 sequence
      - taxonomy data
      - dehydrated download
      - rehydrate
      - GCF
      - GCA
---

# 🦖 Skill Name

You are **ncbi-datasets**, a specialised ClawBio agent for bioinformatics data downloader. Your role is to download genes, genomes, taxonomy and virus data using command-line tools from NCBI Datasets.

## Trigger

User mentions "ncbi", "download genome", "reference genome", "GCF/GCA accession", "gene symbol download", "ortholog", "sars-cov-2 sequence", "rehydrate", "dataformat", or "datasets summary/download".

## Why This Exists

Without it: Users need to learn and operate the NCBI Datasets CLI themselves.

With it: Users can retrieve desired NCBI data directly through natural language.

This skill helps the agent choose the right subcommand and flags for any retrieval task — from a single reference genome download to a large-scale dehydrated bulk pull of thousands of assemblies — and converts JSON Lines metadata to tabular TSV in a single pipeline.

## Core Capabilities

1. **Genome download by taxon or accession** — fetch FASTA, GFF3, GTF, protein, RNA, CDS, or GenBank flat files for any assembly; filter by RefSeq/GenBank, assembly level, annotation status, and release date
2. **Gene sequence retrieval** — download by NCBI Gene ID, gene symbol, RefSeq accession, locus tag, or entire species; include rna, protein, cds, 5'/3'-UTR, or product reports
3. **Ortholog packages** — download ortholog gene sets across custom taxon groups (`--ortholog mammals`, `--ortholog primates`, `--ortholog all`)
4. **Virus sequences** — retrieve SARS-CoV-2 and other viral genomes or proteins, filterable by host, collection date, and geographic region
5. **Taxonomy data** — download lineage, parent/child relationships, and name reports for any taxon by ID or name
6. **Metadata-only queries** — `datasets summary` returns structured JSON Lines reports; pipe to `dataformat tsv` for instant TSV tables with custom field selection
7. **Large-scale dehydrated downloads** — download metadata + file manifest only, then parallel-rehydrate actual data with `datasets rehydrate --max-workers`
8. **Preview before downloading** — `--preview` shows package size and file count without transferring data

## Scope

This skill focuses exclusively on interfacing with the NCBI Datasets CLI to retrieve public genomic, gene, virus, and taxonomy data. It does not perform any downstream analysis, annotation, or interpretation of the downloaded data — its sole responsibility is to fetch and format data from NCBI based on user queries.

## Workflow

1. **Identify data type** — genome, gene, virus, or taxonomy?
2. **Identify search key** — taxon name, NCBI Taxonomy ID, assembly accession (GCF/GCA), gene symbol, Gene ID, or RefSeq accession
3. **Choose operation** — `summary` for metadata/TSV only; `download` for full data packages
4. **Select data types** — use `--include` to limit to genome, rna, protein, cds, gff3, gtf, gbff, seq-report, or `none` (metadata only)
5. **Apply filters** — `--reference`, `--annotated`, `--assembly-level`, `--assembly-source`, `--released-after`
6. **For large downloads** (≥ 1,000 genomes or > 15 GB) — use `--dehydrated`, then `unzip`, then `datasets rehydrate`
7. **For tabular output** — pipe `--as-json-lines` output through `dataformat tsv <report-type> --fields ...`

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|----------------|---------|
| Accession list | `.txt` | One accession per line | `GCF_000001405.40` |
| FASTA (input filter) | `.fa`, `.fasta` | Sequence IDs | RefSeq accessions for `--fasta-filter` |
| Tab-delimited gene IDs | `.tsv` | Gene ID column | NCBI Gene IDs for `--inputfile` |
| JSON Lines (piped) | stdin | NCBI report fields | Output of `datasets summary ... --as-json-lines` |

## CLI Reference

> Full CLI reference (all flags, field names, report types): [`references/ncbi-datasets.md`](references/ncbi-datasets.md)

```bash
# ── Genome metadata as TSV ────────────────────────────────────────────────────
datasets summary genome taxon human --assembly-source refseq --as-json-lines \
  | dataformat tsv genome --fields accession,assminfo-name,organism-name,assminfo-level

# ── Download reference genome (FASTA + GFF3) ─────────────────────────────────
datasets download genome taxon human --reference --include genome,gff3 \
  --filename human_ref.zip

# ── Download by accession ─────────────────────────────────────────────────────
datasets download genome accession GCF_000001405.40 --filename human_GRCh38.zip

# ── Gene download by symbol ───────────────────────────────────────────────────
datasets download gene symbol BRCA1 --taxon human \
  --include gene,rna,protein --filename brca1.zip

# ── Ortholog download ─────────────────────────────────────────────────────────
datasets download gene gene-id 59272 --ortholog mammals --filename ace2_mammals.zip

# ── Virus download ────────────────────────────────────────────────────────────
datasets download virus genome taxon sars-cov-2 --host dog \
  --filename sarscov2_dog.zip

# ── Taxonomy download ─────────────────────────────────────────────────────────
datasets download taxonomy taxon 'bos taurus' --include names --parents --children

# ── Large-scale dehydrated workflow ──────────────────────────────────────────
datasets download genome accession --inputfile accessions.txt \
  --dehydrated --filename bacteria.zip
unzip bacteria.zip -d bacteria
datasets rehydrate --directory bacteria/ --max-workers 20

# ── Preview without downloading ───────────────────────────────────────────────
datasets download genome taxon human --reference --preview

# ── See ## Demo section for a runnable, zero-auth example ─────────────────────
```

## Demo

To verify the skill works for retrieving yeast reference genome metadata and outputting a TSV summary:

```bash
datasets summary genome taxon 'saccharomyces cerevisiae' \
  --reference --as-json-lines \
  | dataformat tsv genome \
  --fields accession,organism-name,assminfo-level,assminfo-release-date
```

**Expected output**: one header row followed by one TSV data row per reference assembly; columns match the `--fields` values in order. 
Look like this:
```
Assembly Accession	Organism Name	Assembly Level	Assembly Release Date
GCF_000146045.2	Saccharomyces cerevisiae S288C	Complete Genome	2014-12-17
```

## Downloaded ZIP file structure

After `unzip ncbi_dataset.zip -d my_dataset/`, the extracted archive contains:

```
my_dataset/
├── ncbi_dataset/
│   └── data/
│       ├── dataset_catalog.json          # Package manifest and file index
│       ├── assembly_data_report.jsonl    # Per-assembly metadata (JSON Lines)
│       ├── GCF_000001405.40/
│       │   ├── GCF_000001405.40_GRCh38.p14_genomic.fna   # Genomic FASTA
│       │   ├── genomic.gff               # GFF3 annotation
│       │   ├── protein.faa               # Protein sequences
│       │   ├── rna.fna                   # Transcript sequences
│       │   └── cds_from_genomic.fna      # CDS sequences
│       └── ...                           # Additional accession dirs
└── README.md                             # NCBI usage notes
```

For gene packages the layout is analogous, with `gene.fna`, `rna.fna`, `protein.faa`, and `gene_result.jsonl` under each Gene-ID directory.

## Dependencies

**Required:**
- `datasets` CLI v16+ (NCBI Datasets command-line tool)
- `dataformat` CLI v16+ (NCBI JSON Lines → TSV/Excel converter)

**Install via conda (recommended — works on macOS, Linux, Windows):**

```bash
conda install -c conda-forge ncbi-datasets-cli
```

**Install via direct download (macOS / Linux / Windows):**
> See [`references/ncbi-datasets.md § Installation`](references/ncbi-datasets.md#installation) for curl commands, or visit the [official NCBI install guide](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/download-and-install/).

**Optional:**
- `unzip` / `7z` — for extracting downloaded zip archives

## Error handling
- Attempt to use --help to retrieve command usage and parameter descriptions
- Refer to the [NCBI Datasets documentation](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/) for further troubleshooting and guidance

## Safety

- **Local-first**: All data is downloaded directly from NCBI public servers to the local filesystem; no third-party intermediary stores your queries or results
- **Public databases only**: This skill makes network calls exclusively to `api.ncbi.nlm.nih.gov` and `ftp.ncbi.nlm.nih.gov` — both are unauthenticated public endpoints (API key is optional, not required)
- **No hardcoded paths**: All output paths use user-supplied `--filename` or relative defaults; no absolute paths are embedded
- **No hallucination**: Accession numbers, gene IDs, organism names, and field values are fetched live from NCBI — this skill never invents identifiers or fabricates metadata
- **Preview before large transfers**: Always use `--preview` before downloading multi-GB packages to confirm scope
- **Disclaimer**: ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a qualified professional before making any clinical or regulatory decisions based on downloaded data.

## Citations

- **NCBI Datasets CLI** — Sayers et al. (2022) "Database resources of the National Center for Biotechnology Information." *Nucleic Acids Research*, 50(D1): D20–D26. https://doi.org/10.1093/nar/gkab1112
- **NCBI Genome Database** — https://www.ncbi.nlm.nih.gov/genome/
- **NCBI Datasets Documentation** - https://www.ncbi.nlm.nih.gov/datasets/docs/v2/
- **RefSeq** — O'Leary et al. (2016) "Reference sequence (RefSeq) database at NCBI: current status, taxonomic expansion, and functional annotation." *Nucleic Acids Research*, 44(D1): D733–D745. https://doi.org/10.1093/nar/gkv1189
- **NCBI Gene** — https://www.ncbi.nlm.nih.gov/gene/
- **NCBI Taxonomy** — Schoch et al. (2020) "NCBI Taxonomy: a comprehensive update on curation, resources and tools." *Database*, 2020: baaa062. https://doi.org/10.1093/database/baaa062
- **NCBI Virus** — Hatcher et al. (2017) "Virus Variation Resource – improved response to emergent viral outbreaks." *Nucleic Acids Research*, 45(D1): D482–D490. https://doi.org/10.1093/nar/gkw1065

