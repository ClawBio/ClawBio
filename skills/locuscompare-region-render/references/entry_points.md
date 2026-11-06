# Where a locuscompare run starts

locuscompare renders one pair: one molecular QTL (the exposure) and one GWAS (the outcome) around one lead variant. Every route into the skill ends in the same config, so the work is turning what the user has in hand into four values:

- `lead.variant_id`: the lead variant, `chrom_pos_ref_alt` on GRCh38
- the exposure: an eQTL Catalogue `dataset_id` plus a gene id, a UKB-PPP protein, or a TSV file
- the outcome: a GWAS Catalog accession, or a TSV file
- the trait labels, for the figure

| The user has | Where the four values come from | Worked example |
|---|---|---|
| An Open Targets credible-set page (a URL) | `ot_credible_set.py options <URL>` lists the colocalising QTLs and checks the GWAS file; `config` writes the run config; see `from_open_targets.md` | `examples/08_open_targets_credible_set/` (one QTL), `examples/09_open_targets_credible_set_many_qtls/` (ten QTLs, GWAS substituted) |
| An Open Targets colocalisation row (two study ids and a lead) | The QTL study's `summarystatsLocation` and `target.id`; the GWAS study id is the accession; see `from_open_targets.md` | `examples/03_open_targets_followup/` |
| An rsID and a `gwas-lookup` result | The GWAS and eQTL tables `gwas-lookup` writes | `examples/04_gwas_lookup_followup/` |
| A protein of interest in UKB-PPP | `source: ukb_ppp` with the protein label and ancestry | `examples/07_pqtl_sort1_ukbppp_eur/` |
| Their own summary statistics, or a source the fetchers do not cover (FinnGen, Pan-UKBB, GTEx v10, ...) | `sumstats_path:` on either side, in the `INPUT_SCHEMA.md` format; `examples/recipes/` converts several sources | `examples/01_synthetic_demo/`, `examples/recipes/` |
| A fine-mapping result | The credible set's lead and the two studies it came from | `examples/chains/finemapping_chain/` |

Two rules hold on every route:

1. **One pair per run.** If the starting point offers several molecular QTLs or several GWAS, list them for the user and ask which pair to render. Do not pick one silently.
2. **The GWAS side must be readable by region.** The GWAS Catalog fetcher reads a window from `harmonised/<accession>.h.tsv.gz` through its tabix index (`.tbi`). Older studies often lack that file: some keep a harmonised file under another name (`<PMID>-<accession>-<EFO>.h.tsv.gz`) with no index, and some have no harmonised folder. For those, choose another study at the same locus or supply the summary statistics as a TSV. (At the IL23R locus below, 6 of the 25 GWAS Catalog studies that colocalise most strongly with the IL23R protein QTL had the indexed file; checked 2026-10-13.)
