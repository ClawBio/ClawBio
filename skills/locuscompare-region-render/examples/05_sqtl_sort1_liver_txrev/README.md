# 06. sQTL × GWAS: SORT1 transcript usage in GTEx liver × cholesterol-VLDL

**v1.3 multi-modality demo, family 2 of 4**. Shows that the regional
LocusCompare orchestrator renders splicing-QTLs (eQTL Catalogue non-ge
quant methods) end-to-end alongside the bulk ge-eQTL renders. No auth
required.

## What this demo proves

> Any OT coloc row referencing an eQTL Catalogue sQTL study (txrev / tx /
> exon / leafcutter quant methods) is renderable today, no auth required.

## Quant-method choice: `txrev`, not `leafcutter`

The folder name on `main` before 2026-05-15 was `06_sqtl_sort1_liver_leafcutter`
and the config pointed at `QTD000270` (GTEx liver leafcutter). Build verified
that the fetcher patch lands rows for QTD000270's region, but those rows
belong to *neighbouring* genes' leafcutter clusters (GSTM1, LAMTOR5-AS1).
SORT1 itself has no leafcutter credible-set hit in GTEx liver, so
filtering by `gene_id=ENSG00000134243` against `QTD000270.cc.tsv.gz`
returns zero rows.

Three other GTEx liver sQTL quant methods do retain a SORT1 credible-set:

| Quant method | eQTL Cat dataset | SORT1 rows in ±500 kb of lead |
|---|---|---|
| `exon` | `QTD000267` | 2,877 (one exon trait) |
| `tx`   | `QTD000268` | 5,754 (two transcripts) |
| `txrev` (this demo) | `QTD000269` | 2,877 (one transcript-usage trait, ENST00000483508) |
| `leafcutter` | `QTD000270` | 0 (no SORT1 cluster in credible-set file) |

`txrev` (proportional transcript usage) is the canonical eQTL Catalogue
splicing-QTL semantic, so it stays the default for the SORT1 demo. The
other three quant methods are available via the same orchestrator path
by swapping `dataset_id` in `config.yaml`.

## Source-file choice: `.cc.tsv.gz`

eQTL Catalogue publishes up to two per-variant FTP files per dataset:

- `<QTD>.all.tsv.gz`: full nominal-pass sumstats, every tested (trait, variant)
  pair. Listed for 306 of the 758 datasets in the catalogue's r7 dataset table
  (every `ge` and `microarray` dataset, plus the one `aptamer` dataset,
  QTD000584). Every `.all` dataset probed also serves a `.cc` (17 of the 306, in
  a 33-dataset probe on 2026-09-13; see eqtl-catalogue-region-fetch gotcha 6);
  the other 452 are listed `.cc` only.
- `<QTD>.cc.tsv.gz`: only the molecular traits with permutation FDR < 1% and
  at least one fine-mapped credible set, each with every tested variant (the
  catalogue's `docs/Columns_parquet.md`, written for its parquet release). For exon / tx / txrevise / leafcutter
  (the only file those methods serve) it further keeps only the top trait per
  independent fine-mapped signal, ~98% smaller (Kerimov 2023 PLoS Genet
  19(9):e1010932, doi:10.1371/journal.pgen.1010932, PMID 37721944).

So a txrevise event that is not the tag of its signal is absent from `.cc`
altogether. Whether the event Open Targets called the coloc on is always the
tag has not been checked; this example's event is present, which is why it
renders.

## OT coloc row resolved

| Field | Value |
|-------|-------|
| OT studyId (exposure) | `gtex_txrev_liver_ensg00000134243` |
| OT studyId (outcome) | `GCST90269602` |
| Pattern | `<study>_<quant>_<sample>_<ensg>` (the canonical OT QTL studyId form) |
| Mapping row | `study_id_mappings.yaml` (caller-supplied, see `load_study_id_mappings`) |

## Upstream studies

- **Exposure**: eQTL Catalogue v7+ dataset `QTD000269` (GTEx liver
  txrev / transcript-usage QTL). Verified live against
  https://www.ebi.ac.uk/eqtl/api/v2/datasets/?size=1000 on 2026-05-15.
- **Outcome**: GWAS Catalog harmonised `GCST90269602`
  (cholesterol in medium VLDL).
- **LD**: 1000G Phase 3 GRCh38, super-pop EUR.

## How to run

```bash
python skills/locuscompare-region-render/cli.py \
    --input skills/locuscompare-region-render/examples/05_sqtl_sort1_liver_txrev/config.yaml \
    --output runs/sqtl_sort1_liver_txrev/
```

No auth required; both exposure (eQTL Catalogue) and outcome
(GWAS Catalog) tabix-on-FTP paths are anonymous. The 1000G region VCF
(~5-50 MB) downloads on first run and caches.

Expected output (committed at `expected_output/`): `manifest.yaml` plus
`1_109274968_G_T_full_locuscompare.png`. Captured 2026-05-15 with
`n_pairs=2648, n_palindromic_excluded=349` (no plink on the dev box,
so the run uses the documented `ld_panel: none` fallback; install plink 1.9
for LD-coloured points).

## Caveats

- `txrev` β is per-ALT-allele change in the *fraction* of expression
  captured by that transcript, not abundance. The β cannot be compared
  directly against a `ge` β at the parent gene; they measure different
  things (isoform switch vs total expression).
- The SORT1 transcript-usage signal at this locus is biologically
  distinct from the ge-eQTL signal in minor salivary gland (example 02).
  A coloc/LocusCompare match between the two is evidence the LDL/CHD
  lead variant alters BOTH expression and splicing, not redundant.
