---
name: eqtl-catalogue-region-fetch
description: |
  Fetch a region of cis-eQTL summary statistics from EBI eQTL Catalogue v7+
  via tabix-on-FTP. Use when an agent needs eQTL beta / SE / p-value for
  every variant in a window around a gene's TSS for one specific dataset
  (study × tissue × quantification method). Input: dataset_id, chromosome,
  start, end, optional molecular_trait_id. Output: harmonised TSV slice.
license: MIT
metadata:
  skill-author: Aviv Madar
  version: 0.1.0
  domain: bioinformatics
  tags:
    - eqtl
    - eqtl-catalogue
    - region-fetch
    - tabix
    - summary-statistics
    - cis-eqtl
  inputs:
    - name: dataset_id
      type: string
      description: eQTL Catalogue dataset identifier (e.g. QTD000276 for GTEx minor salivary gland ge-eQTL).
      required: true
    - name: chromosome
      type: string
      description: Chromosome name without `chr` prefix (1, 2, ..., X, Y, MT).
      required: true
    - name: start_bp
      type: integer
      description: Region start, 1-based GRCh38.
      required: true
    - name: end_bp
      type: integer
      description: Region end, 1-based GRCh38 (inclusive).
      required: true
    - name: molecular_trait_id
      type: string
      description: Optional ENSG (versioned or bare) to filter to one gene; required for ge-eQTL datasets where one TSV bundles multiple traits.
      required: false
  outputs:
    - name: variants
      type: list
      description: Per-variant rows with variant_id, chromosome, position, ref, alt, beta, se, p_value, maf, molecular_trait_id, dataset_id.
    - name: release
      type: object
      description: EQTLCatalogueRelease with study_label, tissue_label, condition_label, sample_group, quant_method, dataset_release, fetched_at_utc.
  dependencies:
    - python>=3.10
    - pysam>=0.22
    - pandas>=2.0
    - requests>=2.28
  demo_data:
    - examples/input.json
  endpoints:
    - https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/    # tabix-on-FTP
    - https://www.ebi.ac.uk/eqtl/api/v3/                          # metadata REST
  openclaw:
    requires:
      bins:
        - python3
        - tabix
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install: |
      pip install pysam pandas requests
    trigger_keywords:
      - eqtl region fetch
      - eqtl catalogue tabix
      - eqtl sumstats slice
      - cis-eqtl region pull
      - GTEx eqtl region
---

# eQTL Catalogue Region Fetch

A standalone primitive that pulls a region of cis-eQTL summary statistics from
EBI eQTL Catalogue v7+ via tabix-on-FTP. Returns harmonised per-variant rows
plus dataset-release metadata.

## Overview

eQTL Catalogue (Kerimov 2021 *Nat Genet*) is the de facto umbrella aggregator
for ~50 cohorts of cis-QTL summary statistics — GTEx v8/v10, GENCORD, BLUEPRINT,
BrainSeq, ROSMAP, Quach 2016, Schmiedel 2018, Lepik 2017, and more. Per-dataset
sumstats are bgzip-compressed + tabix-indexed and served from the EBI FTP at
`https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/<QTS>/<QTD>/<QTD>.all.tsv.gz`.

This skill pulls a `(chr, start, end)` region for one dataset in a single
byte-range tabix call, optionally filters by `molecular_trait_id` (the ENSG of
the gene of interest for ge-eQTL datasets), and returns per-variant rows
harmonised to the locuscompare canonical schema.

**Do NOT use the eQTL Catalogue REST API for region-level fetches.** It
silently truncates regional fetches to one side of the TSS (verified May
2026 against the v2 endpoint at `/api/v2/datasets/{id}/associations`).
Tabix-on-FTP is the supported access pattern.

## When to use this skill

- An agent needs eQTL summary statistics for every variant in a 1 Mb window
  around a gene's TSS, for one specific (study × tissue × quantification) dataset.
- An agent is preparing input for `locuscompare-region-render`,
  `fine-mapping`, or `mendelian-randomisation`.
- An agent wants the canonical summary stats for an Open Targets coloc row's
  `left_studyId` after resolving the OT studyId to an eQTL Cat `dataset_id`.

Do NOT use for: single-rsID lookups across many studies (use `gwas-lookup`),
non-eQTL Cat sources (FinnGen direct, Pan-UKBB, UKB-PPP), or reverse-lookup
"which datasets cover this gene" queries (use the eQTL Cat metadata REST
endpoint directly).

## Inputs

See frontmatter `metadata.inputs`. Per-call inputs are passed as keyword
arguments to `EQTLCatalogueClient.fetch_region(...)`. The `dataset_id` is the
canonical eQTL Cat identifier (e.g. `QTD000276`); look up via the eQTL
Catalogue's "Studies" table at <https://www.ebi.ac.uk/eqtl/Studies/> or the
metadata REST endpoint.

`molecular_trait_id` filtering is strongly recommended for ge-eQTL datasets
since one dataset's TSV bundles all genes; without filtering, a 1 Mb window
returns hundreds of thousands of rows.

## Outputs

`RegionResult` dataclass with:
- `variants`: list of `RegionVariant` rows (variant_id, chromosome, position,
  ref, alt, beta, se, p_value, maf, molecular_trait_id, dataset_id).
- `release`: `EQTLCatalogueRelease` metadata block.
- `n_variants`: count.
- `notes`: provenance / soft-failure notes.

A `--output <dir>` invocation also writes a flat `variants.tsv` with the
columns: `variant_id`, `chromosome`, `position_bp`, `allele_a`, `allele_b`,
`beta`, `se`, `p`, `maf`, `molecular_trait_id`, `study_id` — one row per
variant. This shape is consumable by most downstream coloc, fine-mapping,
and regional-plot tooling without further harmonisation.

## Caveats

- **Effect-allele convention**: the `alt` allele is the effect allele; `beta`
  is per-copy of `alt`. Matches the locuscompare canonical schema.
- **Coordinate system**: GRCh38 throughout.
- **Strand handling**: not strand-flipped (forward strand always); tabix
  indexes are coordinate-sorted as expected.
- **License**: data are CC-BY 4.0 (Kerimov 2021 *Nat Genet*); cite the
  original publication for each dataset. The list of per-dataset citations
  is at <https://www.ebi.ac.uk/eqtl/Studies/>.
- **Rate limits**: the EBI FTP rate-limits aggressive parallel access. Limit
  to ≤4 concurrent connections per process.

## Citations

- Kerimov et al. (2021). *A compendium of uniformly processed human gene
  expression and splicing quantitative trait loci.* Nat Genet 53, 1290-1299.
  doi:10.1038/s41588-021-00924-w
- The per-dataset citation list at <https://www.ebi.ac.uk/eqtl/Studies/>.
