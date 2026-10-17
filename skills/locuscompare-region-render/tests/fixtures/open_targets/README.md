# Recorded Open Targets responses

Replayed by `tests/test_ot_credible_set.py` so the helper's tests run offline.

- Source: Open Targets Platform GraphQL API (`https://api.platform.opentargets.org/api/v4/graphql`), retrieved 2026-10-13. Each file keeps the query, its variables and the retrieval time (`_query`, `_variables`, `_recorded`) next to the response.
- `credibleSet__<id>.json`: the credible set (lead variant, GWAS study).
- `coloc_molqtl__<id>.json`: its colocalisations with molecular QTLs.
- `coloc_gwas__<id>.json`: the GWAS colocalisations of a QTL credible set (`227c96a4...` PCSK9 lung expression, `71fea329...` PCSK9 lung exon, `288fda0f...` IL23R plasma protein).
- `ftp_probe_2026-10-13.json`: the file names in the `harmonised/` folder of six GWAS Catalog studies on `ftp.ebi.ac.uk`, retrieved the same day.

Licence: Open Targets Platform data are marked CC0 1.0 (https://platform-docs.opentargets.org/licence, checked 2026-10-13). GWAS Catalog summary statistics are CC0 1.0 per the same page.
