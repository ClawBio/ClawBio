# NCBI Datasets CLI Skill

Use the following knowledge to help users work with the NCBI Datasets CLI tools (`datasets` and `dataformat`).

---

## Overview

NCBI Datasets provides two CLI tools:
- **`datasets`** — download and summarize biological data from NCBI
- **`dataformat`** — convert JSON Lines metadata into TSV or Excel formats

Data is delivered as zip archives; extracted files reside in `ncbi_dataset/data/`.

---

## Installation

### Conda (recommended)
```bash
conda create -n ncbi_datasets
conda activate ncbi_datasets
conda install -c conda-forge ncbi-datasets-cli
```

### curl — macOS (Universal)
```bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/mac/dataformat'
chmod +x datasets dataformat
```

### curl — Linux (AMD64)
```bash
curl -o datasets 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/datasets'
curl -o dataformat 'https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/linux-amd64/dataformat'
chmod +x datasets dataformat
```

### curl — Windows (64-bit)
```bash
curl -o datasets.exe "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/datasets.exe"
curl -o dataformat.exe "https://ftp.ncbi.nlm.nih.gov/pub/datasets/command-line/v2/win64/dataformat.exe"
```

---

## `datasets` — Command Reference

### Top-level subcommands
| Subcommand | Description |
|---|---|
| `summary` | Print metadata report (genome, gene, virus, taxonomy) |
| `download` | Download data as a zip file |
| `rehydrate` | Rehydrate a dehydrated (large) dataset |
| `completion` | Generate shell autocompletion scripts |

### Global flags
| Flag | Description |
|---|---|
| `--api-key string` | NCBI API key (increases rate limits) |
| `--debug` | Emit debugging info |
| `--help` | Print detailed help |
| `--version` | Print version |

---

## Genome Commands

### `datasets summary genome`

```bash
# By taxon
datasets summary genome taxon human
datasets summary genome taxon 'mus musculus'

# By accession (assembly or BioProject)
datasets summary genome accession GCF_000001405.40
datasets summary genome accession PRJEB33226

# Filters
datasets summary genome taxon human --reference
datasets summary genome taxon human --annotated
datasets summary genome taxon human --assembly-level complete
datasets summary genome taxon human --released-after 2020-01-01
datasets summary genome taxon human --search 'T2T Consortium'

# Sequence report
datasets summary genome accession GCF_000006945.2 --report sequence --as-json-lines \
  | dataformat tsv genome-seq --fields accession,genbank-seq-acc,refseq-seq-acc,chr-name

# TSV output with custom fields
datasets summary genome taxon sharks --assembly-source refseq --as-json-lines \
  | dataformat tsv genome --fields accession,assminfo-name,annotinfo-name,annotinfo-release-date,organism-name
```

**`summary genome` flags:**
| Flag | Values / Default | Description |
|---|---|---|
| `--annotated` | — | Limit to annotated genomes |
| `--assembly-level` | `chromosome,complete,contig,scaffold` | Comma-separated assembly levels |
| `--assembly-source` | `RefSeq`, `GenBank`, `all` (default) | Source database |
| `--assembly-version` | `latest` (default), `all` | Version filter |
| `--as-json-lines` | — | JSON Lines output (required for piping to `dataformat`) |
| `--exclude-atypical` | — | Exclude atypical assemblies |
| `--exclude-multi-isolate` | — | Exclude multi-isolate project assemblies |
| `--from-type` | — | Only records with type material |
| `--limit` | `all` (default), number | Number of results |
| `--mag` | `only`, `exclude`, `all` (default) | Metagenome-assembled genomes |
| `--reference` | — | Limit to reference genomes |
| `--released-after` | ISO 8601 date | Released on or after date |
| `--released-before` | ISO 8601 date | Released on or before date |
| `--report` | `genome` (default), `sequence`, `ids_only` | Report type |
| `--search` | string | Text search (repeatable) |

---

### `datasets download genome`

```bash
# By taxon
datasets download genome taxon human --filename human_dataset.zip

# By accession
datasets download genome accession GCF_000001405.40 --filename human_GRCh38.zip

# With specific data files
datasets download genome taxon human --reference --include genome,protein
datasets download genome taxon human --reference --include genome,rna,cds,protein,gff3

# Metadata only
datasets download genome taxon human --reference --include none

# Specific chromosomes
datasets download genome accession GCF_000001405.40 --chromosomes X,Y --include genome,gff3

# Preview without downloading
datasets download genome taxon human --reference --preview
```

**`--include` values for genome:**
| Value | Description |
|---|---|
| `genome` | Genomic sequences (FASTA) — default |
| `rna` | Transcript sequences |
| `protein` | Amino acid sequences |
| `cds` | Nucleotide coding sequences |
| `gff3` | General feature file |
| `gtf` | Gene transfer format |
| `gbff` | GenBank flat file |
| `seq-report` | Sequence report file |
| `all` | All of the above |
| `none` | No sequences (metadata only) |

**`download genome` flags:**
| Flag | Values / Default | Description |
|---|---|---|
| `--annotated` | — | Limit to annotated genomes |
| `--assembly-level` | `chromosome,complete,contig,scaffold` | Assembly level filter |
| `--assembly-source` | `RefSeq`, `GenBank`, `all` (default) | Source database |
| `--chromosomes` | comma-separated list or `all` | Limit to specific chromosomes |
| `--dehydrated` | — | Download metadata only, rehydrate later |
| `--exclude-atypical` | — | Exclude atypical assemblies |
| `--fast-zip-validation` | — | Skip zip checksum validation |
| `--filename` | string | Output filename (default: `ncbi_dataset.zip`) |
| `--include` | see table above | Data file types to include |
| `--no-progressbar` | — | Hide progress bar |
| `--preview` | — | Show package info without downloading |
| `--reference` | — | Limit to reference genomes |
| `--released-after` | ISO 8601 date | Released on or after date |
| `--released-before` | ISO 8601 date | Released on or before date |
| `--search` | string | Text search (repeatable) |

---

### Large Genome Downloads (Dehydrate → Rehydrate)

Use this workflow for ≥ 1,000 genomes or packages > 15 GB:

```bash
# Step 1: Download dehydrated package (metadata + file list only)
datasets download genome accession --inputfile accessions.txt --dehydrated --filename my-genomes.zip

# Step 2: Unzip
unzip my-genomes.zip -d my-genomes

# Step 3: Download actual data files
datasets rehydrate --directory my-genomes/
```

**`rehydrate` flags:**
| Flag | Description |
|---|---|
| `--directory string` | Directory containing unzipped dehydrated package |
| `--gzip` | Rehydrate to gzip format |
| `--list` | Dry run: list files that would be downloaded |
| `--match string` | Only rehydrate files matching this substring |
| `--max-workers int` | Concurrent download workers (1–30, default: 10) |
| `--no-progressbar` | Hide progress bar |

---

## Gene Commands

### `datasets summary gene`

```bash
# By NCBI Gene IDs
datasets summary gene gene-id 1 2 3 9 10

# By gene symbols
datasets summary gene symbol ACRV1 A2M --taxon human
datasets summary gene symbol brca1 --taxon 'mus musculus'

# By RefSeq accession
datasets summary gene accession NM_020107.5 NP_001334352.2

# By taxon (all genes in species)
datasets summary gene taxon human

# Product report
datasets summary gene symbol ACRV1 --report product

# TSV output
datasets summary gene symbol ACRV1 A2M --as-json-lines \
  | dataformat tsv gene --fields symbol,gene-id,synonyms

# Predefined summary template
datasets summary gene symbol ACRV1 A2M --as-json-lines \
  | dataformat tsv gene --template summary

# Gene ontology table
datasets summary gene symbol ACRV1 --as-json-lines \
  | dataformat tsv --template gene-ontology
```

**`summary gene` flags:**
| Flag | Values / Default | Description |
|---|---|---|
| `--as-json-lines` | — | JSON Lines output (required for piping to `dataformat`) |
| `--limit` | `all` (default), number | Number of results |
| `--report` | `gene` (default), `product`, `ids_only` | Report type |

**`summary gene` subcommands:**
- `gene-id` — by NCBI Gene ID
- `symbol` — by gene symbol (requires `--taxon` for disambiguation)
- `accession` — by RefSeq nucleotide or protein accession
- `taxon` — by species (NCBI Taxonomy ID, scientific or common name)
- `locus-tag` — by locus tag

---

### `datasets download gene`

```bash
# By Gene ID
datasets download gene gene-id 672

# By gene symbol
datasets download gene symbol ACRV1 A2M --taxon human
datasets download gene symbol brca1 --taxon 'mus musculus'

# By RefSeq accession
datasets download gene accession NM_020107.5 NP_001334352.2

# All human genes
datasets download gene taxon human

# Include specific data types
datasets download gene gene-id 672 --include gene,protein
datasets download gene gene-id 672 --include gene,rna,cds,protein

# Metadata only
datasets download gene gene-id 672 --include none

# Orthologs
datasets download gene gene-id 59272 --ortholog mammals       # ACE2 in mammals
datasets download gene symbol cftr --ortholog all              # CFTR complete set
datasets download gene accession NM_000492.4 --ortholog primates

# Filter sequences by accession
datasets download gene gene-id 2778 --fasta-filter NC_000020.11,NM_001077490.3
```

**`--include` values for gene:**
| Value | Description |
|---|---|
| `rna` | Transcript sequences — default |
| `protein` | Amino acid sequences — default |
| `gene` | Gene sequences |
| `cds` | Nucleotide coding sequences |
| `5p-utr` | 5'-UTR sequences |
| `3p-utr` | 3'-UTR sequences |
| `product-report` | Transcript and protein locations/metadata |
| `none` | No sequences (metadata only) |

**`download gene` flags:**
| Flag | Description |
|---|---|
| `--fasta-filter strings` | Limit sequences to specified RefSeq accessions (comma-separated) |
| `--fasta-filter-file string` | Read accessions from a file |
| `--filename string` | Output filename (default: `ncbi_dataset.zip`) |
| `--include string` | Data file types to include (default: `rna,protein`) |
| `--no-progressbar` | Hide progress bar |
| `--ortholog string` | Download ortholog set (`all` or a taxon, e.g. `mammals`, `primates`) |
| `--preview` | Show package info without downloading |

---

## Virus Commands

```bash
# Download SARS-CoV-2 genomes from dog hosts
datasets download virus genome taxon sars-cov-2 --host dog

# Download SARS-CoV-2 spike protein from dog hosts
datasets download virus protein S --host dog --filename SARS2-spike-dog.zip

# Get virus metadata
datasets summary virus genome taxon sars-cov-2
```

---

## Taxonomy Commands

```bash
# Download taxonomy for a taxon
datasets download taxonomy taxon 'bos taurus'

# Multiple taxa; include name report
datasets download taxonomy taxon human,'drosophila melanogaster' --include names

# Include parent and child nodes
datasets download taxonomy taxon 10116 --parents --children
```

**Default taxonomy package:** `taxonomy_report.jsonl`, `taxonomy_summary.tsv`, `dataset_catalog.json`

---

## `dataformat` — Command Reference

Converts NCBI JSON Lines metadata into tabular formats.

### Subcommands
| Subcommand | Description |
|---|---|
| `tsv` | Convert to TSV |
| `excel` | Convert to Excel workbook |
| `catalog` | Print the package catalog |
| `completion` | Shell autocompletion |
| `version` | Print version |

### Report types (used with `tsv` or `excel`)
| Report type | Description |
|---|---|
| `genome` | Genome Assembly Data Report |
| `genome-seq` | Genome Assembly Sequence Report |
| `genome-annotations` | Genome Annotation Report |
| `gene` | Gene Report |
| `gene-product` | Gene Product Report |
| `prok-gene` | Prokaryote Gene Report |
| `prok-gene-location` | Prokaryote Gene Location Report |
| `virus-genome` | Virus Data Report |
| `virus-annotation` | Virus Annotation Report |
| `taxonomy` | Taxonomy Report |
| `organelle` | Organelle Report |
| `biocollection` | Biocollection Report |
| `microbigge` | MicroBIGG-E Data Report |

### `dataformat tsv` key flags
| Flag | Description |
|---|---|
| `--fields field1,field2,...` | Specify columns to output |
| `--template name` | Use predefined template (`summary`, `gene-ontology`) |
| `--elide-header` | Omit the header row |
| `--force` | Skip type check prompt |

### Examples
```bash
# Genome TSV with custom fields
datasets summary genome taxon sharks --assembly-source refseq --as-json-lines \
  | dataformat tsv genome --fields accession,assminfo-name,organism-name

# Gene TSV with template
datasets summary gene symbol ACRV1 --as-json-lines \
  | dataformat tsv gene --template summary

# Gene product report
datasets summary gene symbol ACRV1 --report product --as-json-lines \
  | dataformat tsv gene-product --template summary

# Gene ontology
datasets summary gene symbol ACRV1 --as-json-lines \
  | dataformat tsv --template gene-ontology

# Sequence report
datasets summary genome accession GCF_000006945.2 --report sequence --as-json-lines \
  | dataformat tsv genome-seq --fields accession,genbank-seq-acc,refseq-seq-acc,chr-name

# Excel output
datasets summary genome taxon human --reference --as-json-lines \
  | dataformat excel genome --filename human_genomes.xlsx
```

---

## Shell Autocompletion

```bash
# bash
datasets completion bash >> ~/.bashrc

# zsh
datasets completion zsh >> ~/.zshrc

# fish
datasets completion fish > ~/.config/fish/completions/datasets.fish

# PowerShell
datasets completion powershell >> $PROFILE
```

---

## Key Patterns and Tips

| Pattern | Notes |
|---|---|
| Always use `--as-json-lines` when piping to `dataformat` | Required for correct parsing |
| Use `--inputfile accessions.txt` for large ID lists | One accession per line |
| Large downloads (≥1000 genomes or >15 GB): use `--dehydrated` + `rehydrate` | Avoids interrupted large transfers |
| Preview before downloading | `--preview` shows what would be downloaded |
| API rate limits | Use `--api-key <key>` to increase limits |
| Custom output filename | `--filename my_data.zip` |
| Dry run rehydration | `datasets rehydrate --list --directory ...` |

---

## Common Field Names for `dataformat --fields`

**Genome fields:** `accession`, `assminfo-name`, `organism-name`, `organism-tax-id`, `annotinfo-name`, `annotinfo-release-date`, `assminfo-level`, `assminfo-submission-date`, `assminfo-release-date`

**Genome-seq fields:** `accession`, `genbank-seq-acc`, `refseq-seq-acc`, `chr-name`, `seq-length`

**Gene fields:** `symbol`, `gene-id`, `synonyms`, `description`, `tax-name`, `tax-id`, `gene-type`, `transcript-count`

**Gene-product fields:** (use `--template summary` for predefined columns)
