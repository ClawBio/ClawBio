# Starting from an Open Targets page

A user looking at an Open Targets Platform credible-set page, for example

    https://platform.opentargets.org/credible-set/2f6755df0a7737586c9ff36fa2643c61

wants the regional LocusCompare plot for the GWAS signal on that page and one molecular QTL that colocalises with it.

## The page cannot be read; its id can

The Platform is a JavaScript application. Fetching the URL returns an empty shell with none of the ids the user sees. The 32-character id at the end of the URL is the credible set's `studyLocusId`; query the Open Targets GraphQL API with it (POST JSON to `https://api.platform.opentargets.org/api/v4/graphql`) and you get everything the page shows.

## Step 1: the GWAS side and the lead variant

```graphql
{ credibleSet(studyLocusId: "2f6755df0a7737586c9ff36fa2643c61") {
    variant { id rsIds }
    study { id studyType traitFromSource publicationFirstAuthor }
} }
```

`variant.id` is already in locuscompare's form (`1_67240275_G_A`). `study.id` is the GWAS Catalog accession when it starts with `GCST`. A credible set whose `studyType` is not `gwas` is itself a QTL; start from the GWAS side instead.

## Step 2: the molecular QTLs that colocalise with it

```graphql
{ credibleSet(studyLocusId: "2f6755df0a7737586c9ff36fa2643c61") {
    colocalisation(studyTypes: [eqtl, pqtl, sqtl, tuqtl, sceqtl, scpqtl, scsqtl, sctuqtl],
                   page: {index: 0, size: 100}) {
      count
      rows { h4 colocalisationMethod
             otherStudyLocus { study { id studyType projectId
                                       summarystatsLocation
                                       target { id approvedSymbol }
                                       biosample { biosampleName } } } }
} } }
```

Pass `studyTypes` explicitly: without it the list holds GWAS colocalisations only, often many of them, and none of the QTLs. If `count` is more than 1, show the user the QTLs (gene, study, tissue, H4) and ask which one to render.

## Step 3: map the QTL study to the exposure config

| Open Targets field | Config field |
|---|---|
| `summarystatsLocation` ending `.../sumstats/<QTS>/<QTD>` | `exposure.fetch.source: eqtl_catalogue`, `dataset_id: <QTD>` (the last path segment) |
| `target.id` (an Ensembl gene id) | `exposure.fetch.molecular_trait_id` (filters the dataset to that gene's rows) |
| `projectId` `UKB_PPP_EUR` | `exposure.fetch.source: ukb_ppp`; see `examples/07_pqtl_sort1_ukbppp_eur/` |
| any other `summarystatsLocation` | not covered by the bundled fetchers; supply a TSV via `sumstats_path:` |

## Step 4: check the GWAS side is readable

The GWAS Catalog fetcher reads `harmonised/<accession>.h.tsv.gz` through its tabix index. Open Targets' `hasSumstats` does not guarantee that file: older studies keep a harmonised file under another name (`<PMID>-<accession>-<EFO>.h.tsv.gz`) with no index, or have no harmonised folder. List the study's `harmonised/` directory on `ftp.ebi.ac.uk` and confirm both `<accession>.h.tsv.gz` and `<accession>.h.tsv.gz.tbi` are there. If it is not, tell the user, and offer another GWAS that colocalises with the same QTL and has an index (Step 2, run on the QTL's own credible set with `studyTypes: [gwas]`).

## Worked example

`2f6755df0a7737586c9ff36fa2643c61` is inflammatory bowel disease (GCST90225550, Mize 2022) at IL23R, led by rs11209026 (`1_67240275_G_A`, the IL23R R381Q missense variant). Its only molecular-QTL colocalisation is the IL23R plasma protein QTL from Sun 2018 (`sun_2018_aptamer_plasma_il23r_5088_175_3_1`, H4 0.991), whose `summarystatsLocation` ends in `QTD000584` and whose `target.id` is `ENSG00000162594`. The resulting config is `examples/08_open_targets_credible_set/config.yaml`.

The same QTL also colocalises with an IBD credible set from GCST90038683 (`c744d25bd3225993b59f292a2ff225bb`), but that study's harmonised file is `33959723-GCST90038683-EFO_0003767.h.tsv.gz`, with no index, so Step 4 sends you to a study like GCST90225550.
