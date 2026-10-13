# Examples

Seven numbered demos covering the main entry vectors and QTL modalities for
`locuscompare`, one chain pattern, plus four harmonization recipes for
non-bundled data sources.

## When to use which (numbered demos, runnable via `--demo <NAME>`)

| Example | Entry vector | Network required | Use when... |
|---|---|---|---|
| `01_synthetic_demo/` | `--demo` flag | No | Smoke-testing the install; CI; offline use |
| `02_eqtl_catalogue_x_gwas_catalog/` | bundled fetchers | Yes | You have an eQTL Cat dataset_id + GWAS Cat accession |
| `03_open_targets_followup/` | OT coloc-row resolution | Yes | You found an interesting row in an Open Targets coloc query |
| `04_gwas_lookup_followup/` | output of `gwas-lookup` | Yes | You found an interesting variant via `gwas-lookup` |
| `05_sqtl_sort1_liver_txrev/` | bundled fetchers (`.cc.tsv.gz`) | Yes | sQTL exposure (splicing); SORT1 in GTEx liver |
| `06_sceqtl_sort1_onek1k_cd14_mono/` | bundled fetchers | Yes | sceQTL exposure (single-cell eQTL); SORT1 in OneK1K CD14+ monocytes |
| `07_pqtl_sort1_ukbppp_eur/` | bundled fetchers (UKB-PPP dispatch) | Yes | pQTL exposure via `UKB_PPP_*` dispatch; SORT1 plasma sortilin |
| `08_open_targets_credible_set/` | bundled fetchers (eQTL Catalogue aptamer + GWAS Catalog) | Yes | start from an Open Targets credible-set page; IL23R plasma protein × IBD |

## Chain pattern (not in `--list-demos`; requires upstream skill output first)

`chains/finemapping_chain/` documents how to render a regional LocusCompare on top
of a co-localizing pair of credible sets emitted by ClawBio `fine-mapping`. Run
`fine-mapping` on both sides first, then point this config's `sumstats_path:`
at each side's `<output>/inputs/locus.tsv`. Not enumerable via `--demo` because
it has no self-contained data; see `chains/finemapping_chain/README.md`.

## Recipes for non-bundled sources

`recipes/` ships harmonization scripts that convert third-party sumstats into
the canonical TSV schema documented in `../INPUT_SCHEMA.md`. After harmonizing,
point a `sumstats_path` at the produced TSV in your config — locuscompare
treats it identically to a bundled-fetcher result.

| Recipe | Source | License tier | Auth |
|---|---|---|---|
| `finngen_direct/` | FinnGen direct (latest release) | Open-access | None |
| `pan_ukbb_direct/` | Pan-UKBB FTP | Open-access | None |
| `ukb_ppp_pqtl/` | UKB-PPP plasma pQTL | Yellow (registration-gated) | UKB approval |
| `gtex_v10_direct/` | GTEx v10 cis-eQTL bgz | Open-access | None |

Each recipe is a small bash + awk + bgzip + tabix pipeline (~20 lines) plus a
README documenting the source URL, the column-rename map, and any caveats.
