# Starting from an Open Targets page

A user looking at an Open Targets Platform credible-set page, for example

    https://platform.opentargets.org/credible-set/2f6755df0a7737586c9ff36fa2643c61

wants the regional LocusCompare plot for the GWAS signal on that page and one molecular QTL that colocalises with it. A typical request: "I'm on <that link>; plot LocusCompare for gtex-v10_ge_lung_ensg00000169174", or less precisely "for the PCSK9 lung QTL".

## The page cannot be read; its id can

The Platform is a JavaScript application. Fetching the URL returns an empty shell with none of the ids the user sees. The 32-character id at the end of the URL is the credible set's `studyLocusId`, and the Open Targets GraphQL API returns everything the page shows from it. `ot_credible_set.py` (in this skill's directory, standard library only) runs those queries and makes the two choices below; the raw queries are in the appendix.

## Step 1: list the QTLs and recommend one

```bash
python ot_credible_set.py options https://platform.opentargets.org/credible-set/<id> [--gene SYMBOL]
```

It prints the page's GWAS study and lead variant, then table A: every molecular QTL that colocalises with the signal, one row per gene, tissue, study and QTL kind (rows that measure the same thing are folded and the highest H4 kept; a gene-level and an exon-level eQTL are different kinds and stay separate). Rows are ordered by the gene named with `--gene`, then H4, then sample size; at most 15 are printed and the number not printed is stated. The last column says whether the bundled fetchers can read the QTL (eQTL Catalogue r7 or UKB-PPP), and flags an older release: Open Targets colocalises GTEx v10, eQTL Catalogue r7 holds GTEx v8.

Pass whatever the user said about the QTL to `--qtl`: a row number, an Open Targets study id or credible-set id copied from the page, or a description such as `"PCSK9 lung"` (every word must appear in the row's gene, tissue, QTL kind or study). One match is used as is. Several matches are listed on their own with a recommendation among them: ask the user which one. No match is said, not guessed. With no `--qtl`, show the user table A and the recommendation with its reasons, then ask which QTL to plot. Do not pick silently when there is more than one.

## Step 2: check the GWAS side, and substitute only when it is definitely unreadable

```bash
python ot_credible_set.py options <id> --qtl <row number or study id>
```

With a QTL chosen (or when the page lists only one), the helper checks the page study's `harmonised/` folder on the GWAS Catalog FTP for `<accession>.h.tsv.gz` and its `.tbi`, the files the fetcher reads. Open Targets' `hasSumstats` does not guarantee them: older studies keep a harmonised file under another name (`<PMID>-<accession>-<EFO>.h.tsv.gz`) with no index, or have no harmonised folder. In 29 GWAS Catalog studies listed on 2026-10-13, 18 had only a file under another name and none of those had an index, 8 had `<accession>.h.tsv.gz` and all 8 had its `.tbi`, and 3 had no harmonised folder.

- **Indexed**: the page study is the outcome.
- **Definitely not readable** (no indexed file, no harmonised folder, or not a GWAS Catalog study): table B lists the GWAS studies that colocalise with the chosen QTL, same Open Targets disease first, each with its own file check, and recommends the best indexed same-disease study. If there is none it says so and recommends the best indexed study of a different trait, flagged as such.
- **Check failed** (timeout, network or server error) or skipped: the page study is kept, the output says its file could not be checked, and no substitute is offered unless you pass `--show-substitutes`. A failed check is never treated as "no" and never recommended.

Relay the recommendation and, when the GWAS is substituted, say so in plain words: the plot shows study X, not the page's study.

## Step 3: write the config and render

```bash
python ot_credible_set.py config <id> --qtl <study id> [--gwas <GCST>] --out config.yaml
python cli.py --input config.yaml --output <dir>
```

`config` refuses (exit 2) without `--qtl`, with a description that matches several QTLs, without `--gwas` when the page study is definitely unreadable, and with a `--gwas` that does not colocalise with the chosen QTL or has no indexed file; it exits 3 when Open Targets cannot be read and 4 when the output file cannot be written. The lead is the page's lead variant, or the member of the page's credible set named with `--lead` (recorded as a substitution). `config` also writes every member of the credible set, with its PIP, as `lead.candidates` (PIP high to low; if Open Targets reports more members than it returned, a caveat says so). Whether a lead can carry the plot is only known when the data are read, so the render decides: it uses the configured lead if it is in the exposure data, the outcome data and the 1000 Genomes LD panel, is not strand-ambiguous (A/T or G/C) and lies within a tenth of the window of its centre; otherwise the highest-PIP member that is, preferring presence in both datasets over presence in the LD panel and keeping the configured lead on a tie. The x-range stays on the window fetched around the requested lead, and the caveat names up to three higher-PIP members that were passed over, with the reason. A swap is printed as a caveat and recorded in the manifest (`lead_fallback`, plus one status row per candidate in `lead_candidates`); delete `candidates` from the config to pin the lead. The window is 1 Mb. The exposure is the eQTL Catalogue dataset id from `summarystatsLocation` plus the gene id from `target.id`, or UKB-PPP. Under `provenance:` it records the Open Targets ids, both H4 values and effect directions, and each substitution as `{side, requested, used, reason}`; each substitution is also a plain sentence under `caveats:`, which the render prints.

The eQTL Catalogue dataset id is read from anywhere in `summarystatsLocation`, not from its last path segment: of the 11 molecular-QTL rows on the two worked-example pages (2026-10-13), 9 point to a `gs://` file ending `.../QTDxxxxxx/QTDxxxxxx.cc.parquet` and 2 to an FTP folder ending `/QTSxxxxxx/QTDxxxxxx`.

## Which tissue

The tissue in table A is where Open Targets found the colocalisation, not evidence that the gene acts there. At the PCSK9 locus below, 8 of the 10 molecular-QTL rows are PCSK9, in 5 tissues (lung, tibial nerve, blood, subcutaneous adipose, iPSC), at H4 0.746 to 0.999 (Open Targets, 2026-10-13). Local regulatory effects are often shared across many tissues: "cis-acting genetic variants tend to affect either most tissues or a small number of tissues" (GTEx Consortium. Genetic effects on gene expression across human tissues. Nature. 2017;550(7675):204-213. doi:10.1038/nature24277. PMID: 29022597). Choosing the relevant tissue from expression data is hard: in GTEx v8 "the relatively weak signal indicates that pinpointing the likely relevant tissue GWAS loci is challenging" (GTEx Consortium. The GTEx Consortium atlas of genetic regulatory effects across human tissues. Science. 2020;369(6509):1318-1330. doi:10.1126/science.aaz1776. PMID: 32913098). Say so when you relay a recommendation that rests on the highest H4.

## Worked examples

`2f6755df0a7737586c9ff36fa2643c61` is inflammatory bowel disease (GCST90225550, Mize 2022) at IL23R, led by rs11209026 (`1_67240275_G_A`, the IL23R R381Q missense variant). It lists one molecular QTL, the IL23R plasma protein QTL from Sun 2018 (H4 0.991, eQTL Catalogue QTD000584), and the page study is indexed, so `options` selects it and `config --qtl 1` writes the lead, exposure and outcome of `examples/08_open_targets_credible_set/config.yaml`.

`e723d7b2d797c4c3d9e3889ef3e79997` is high cholesterol (GCST90013932, Mbatchou 2021) near PCSK9, led by rs4500361. It lists 10 molecular QTLs; the helper recommends PCSK9 gene expression in lung (H4 0.9992, GTEx v10 in Open Targets, read as GTEx v8 from eQTL Catalogue r7). The page study has no indexed harmonised file, so `options --qtl 1` recommends GCST90475721 (hypercholesterolemia, Verma 2024), the best indexed same-disease study. `config` keeps the page's lead rs4500361 (`1_55025188_T_C`) and writes the five credible-set members as candidates. At render time that lead is absent from the 1000 Genomes LD panel, the next member by PIP (`1_55021673_C_G`) is a C/G variant, and the render draws its lead at `1_55023869_C_T` (PIP 0.117), the highest-PIP member in the eQTL data, the GWAS and the panel that is not strand-ambiguous (measured 2026-10-13; `1_55026242_C_T`, the substitute study's own lead, is absent from GTEx v8 lung, QTD000271). `examples/09_open_targets_credible_set_many_qtls/` shows the output and the config.

## Appendix: the raw queries

POST JSON (`{"query": ..., "variables": ...}`) to `https://api.platform.opentargets.org/api/v4/graphql`.

```graphql
{ credibleSet(studyLocusId: "2f6755df0a7737586c9ff36fa2643c61") {
    studyType
    variant { id rsIds }
    study { id studyType traitFromSource publicationFirstAuthor diseases { id name } }
} }
```

`variant.id` is already in locuscompare's form (`1_67240275_G_A`). A credible set whose `studyType` is not `gwas` is itself a QTL; start from the GWAS side instead.

```graphql
{ credibleSet(studyLocusId: "2f6755df0a7737586c9ff36fa2643c61") {
    colocalisation(studyTypes: [eqtl, pqtl, sqtl, tuqtl, sceqtl, scpqtl, scsqtl, sctuqtl],
                   page: {index: 0, size: 100}) {
      count
      rows { h4 betaRatioSignAverage colocalisationMethod
             otherStudyLocus { studyLocusId
                               study { id studyType projectId summarystatsLocation
                                       target { id approvedSymbol }
                                       biosample { biosampleId biosampleName } } } }
} } }
```

Pass `studyTypes` explicitly: without it the list holds GWAS colocalisations only, often many of them, and none of the QTLs. The same query on the chosen QTL's `studyLocusId` with `studyTypes: [gwas]` lists the GWAS studies that could stand in for the page's.

| Open Targets field | Config field |
|---|---|
| `summarystatsLocation` containing `QTDxxxxxx` | `exposure.fetch.source: eqtl_catalogue`, `dataset_id: QTDxxxxxx` |
| `target.id` (an Ensembl gene id) | `exposure.fetch.molecular_trait_id` (selects that gene's rows) |
| `projectId` `UKB_PPP_EUR` | `exposure.fetch.source: ukb_ppp`; see `examples/07_pqtl_sort1_ukbppp_eur/` |
| any other `summarystatsLocation` | not covered by the bundled fetchers; supply a TSV via `sumstats_path:` |
