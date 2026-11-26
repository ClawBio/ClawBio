<!-- Output of: python clawbio.py run snp-phewas --demo --output <dir>  (report.md; date and wall time vary) -->

# Variant PheWAS across biobanks: rs12740374

**Date**: 2026-10-30 17:37 UTC
**Skill**: snp-phewas
**Version**: 0.3.0
**Variant (GRCh38)**: 1_109274968_G_T
**Effect allele**: T (other allele G)
**Datasets**: Pan-UK Biobank (European ancestry), FinnGen R12 (Finnish) and Million Veteran Program (European genetic ancestry)

---

## What this answer covers

- **Demo data.** This run reads bundled slices of 3 Pan-UK Biobank outcome files, 3 FinnGen endpoint files and 2 Million Veteran Program study files around rs12740374, offline; the rsID and the GRCh37 position use answers recorded from Open Targets and Ensembl. Every other outcome is listed as not part of this lookup.

## Summary

- **Pan-UK Biobank, European ancestry**: tested in 3 of 3 outcomes read (7,160 in the catalog); 3 below the threshold 0.05 / 3 (over the 3 outcomes looked up).
- **FinnGen R12 (Finnish)**: tested in 3 of 3 outcomes read (2,469 in the catalog); 3 below the threshold 0.05 / 3 (over the 3 outcomes looked up).
- **Million Veteran Program, European genetic ancestry**: tested in 2 of 2 outcomes read (2,070 in the catalog, 849 of them with no readable file); 2 below the threshold 0.05 / 2 (over the 2 outcomes looked up).
- **Replication, Pan-UK Biobank (European ancestry) and FinnGen R12 (Finnish)** (threshold 0.05 / 589 shared outcomes): 2 replicate, 0 same direction significant in one only, 0 opposite directions, 0 null in both.
- **Replication, Pan-UK Biobank (European ancestry) and Million Veteran Program (European genetic ancestry)** (threshold 0.05 / 724 shared outcomes): 2 replicate, 0 same direction significant in one only, 0 opposite directions, 0 null in both.
- **Replication, FinnGen R12 (Finnish) and Million Veteran Program (European genetic ancestry)** (threshold 0.05 / 157 shared outcomes): 1 replicate, 0 same direction significant in one only, 0 opposite directions, 0 null in both.
- **Replication, all three** (threshold 0.05 / 157 outcomes all three share): 1 replicate in all three, 0 in two.

## Notes

- rs12740374 is 1_109274968_G_T (GRCh38), resolved through Open Targets.
- Pan-UK Biobank, European ancestry: tested in 3 of 7,160 outcomes.
- Pan-UK Biobank, European ancestry: 7,157 outcomes have no value: outcome was not part of this lookup.
- Pan-UK Biobank, European ancestry: read from 3 files in 0 s.
- Pan-UK Biobank, European ancestry: only 3 of 7,160 outcomes were requested; the rest are listed as not part of this lookup.
- Pan-UK Biobank, European ancestry uses GRCh37: GRCh38 1:109274968 is GRCh37 1:109817590.
- FinnGen R12 (Finnish): tested in 3 of 2,469 outcomes.
- FinnGen R12 (Finnish): 2,466 outcomes have no value: outcome was not part of this lookup.
- FinnGen R12 (Finnish): read from 0 files plus one web-API request in 0 s.
- FinnGen R12 (Finnish): only 3 of 2,469 outcomes were requested; the rest are listed as not part of this lookup.
- FinnGen R12 (Finnish): the source also returned 1 outcome(s) its published catalog does not list; they are kept separately, not in the counts above.
- Million Veteran Program, European genetic ancestry: tested in 2 of 2,070 outcomes.
- Million Veteran Program, European genetic ancestry: 2,068 outcomes have no value; the most common reason (1,219) is: outcome was not part of this lookup.
- Million Veteran Program, European genetic ancestry: 1,221 of 2,070 outcomes are readable; the other 849 cannot be read from the GWAS Catalog (448 not harmonized by the GWAS Catalog; 401 deposited without summary statistics); they are listed, not counted as nulls, and are outside the significance threshold's denominator.
- Million Veteran Program, European genetic ancestry: read from 2 files in 0 s.
- Million Veteran Program, European genetic ancestry: only 2 of 2,070 outcomes were requested; the rest are listed as not part of this lookup.
- BioBank Japan (NBDC hum0197) is not included: its terms limit use to research and development purposes.

## Strongest associations: Pan-UK Biobank, European ancestry

| Outcome | Outcome id | beta | SE | p | EAF | Status |
|---|---|---|---|---|---|---|
| LDL direct | `biomarkers-30780-both_sexes-irnt` | -0.1345 | 0.002913 | 1.58e-465 | 0.222 | tested |
| E78 Disorders of lipoprotein metabolism and other lipidaemias | `icd10-E78-both_sexes` | -0.1545 | 0.009285 | 3.31e-62 | 0.224 | tested |
| Ischemic Heart Disease | `phecode-411-both_sexes` | -0.1106 | 0.01001 | 2.14e-28 | 0.223 | tested |

## Strongest associations: FinnGen R12 (Finnish)

| Outcome | Outcome id | beta | SE | p | EAF | Status |
|---|---|---|---|---|---|---|
| Disorders of lipoprotein metabolism and other lipidaemias | `E4_LIPOPROT` | -0.1286 | 0.008064 | 3.23e-57 | 0.214 | tested |
| Ischaemic heart disease, wide definition | `I9_IHD` | -0.07619 | 0.007514 | 3.68e-24 | 0.215 | tested |
| Major coronary heart disease event | `I9_CHD` | -0.08011 | 0.00864 | 1.82e-20 | 0.215 | tested |

## Strongest associations: Million Veteran Program, European genetic ancestry

| Outcome | Outcome id | beta | SE | p | EAF | Status |
|---|---|---|---|---|---|---|
| low density lipoprotein cholesterol (LDLC, mean of EHR values, inverse-normal transformed) | `GCST90475416` | -0.1561 | 0.003009 | 6.02e-587 | 0.218 | tested |
| Ischemic heart disease | `GCST90475929` | -0.09349 | 0.005226 | 1.39e-71 | 0.219 | tested |


## Replication across biobanks

1,537 outcomes are harmonized across at least two of the datasets (the same ICD-10 block or phecode, FinnGen's own published match, or a reviewed match of the same lab measurement). An outcome replicates in a pair when both p-values are below 0.05 / n, n = the outcomes both datasets of the pair share, with the same direction.

**Pan-UK Biobank (European ancestry) and FinnGen R12 (Finnish)**: 589 shared outcomes, threshold 0.05 / 589.
- not read in one of the datasets (an outage, outside this lookup, or a study the source holds with no readable file; not a null): 587
- replicates (significant in both, same direction): 2 (Ischemic Heart Disease, E78 Disorders of lipoprotein metabolism and other lipidaemias)
587 harmonized outcomes could not be compared; they are listed apart from the nulls, never counted as nulls.

**Pan-UK Biobank (European ancestry) and Million Veteran Program (European genetic ancestry)**: 724 shared outcomes, threshold 0.05 / 724.
- not read in one of the datasets (an outage, outside this lookup, or a study the source holds with no readable file; not a null): 1,179
- replicates (significant in both, same direction): 2 (LDL direct, Ischemic Heart Disease)
1,179 harmonized outcomes could not be compared; they are listed apart from the nulls, never counted as nulls.

**FinnGen R12 (Finnish) and Million Veteran Program (European genetic ancestry)**: 157 shared outcomes, threshold 0.05 / 157.
- not read in one of the datasets (an outage, outside this lookup, or a study the source holds with no readable file; not a null): 232
- replicates (significant in both, same direction): 1 (Ischemic Heart Disease)
232 harmonized outcomes could not be compared; they are listed apart from the nulls, never counted as nulls.

Different biobanks are treated as independent: they recruit different people (UK Biobank, FinnGen and the Million Veteran Program), which is a matter of study design, not something this lookup measures.

**All three biobanks** (Pan-UKBB EUR, FinnGen and MVP EUR): 157 outcomes with a readable file are shared by all three; each is compared at 0.05 / 157.
- not read in one of the three (an outage, outside this lookup, or a study the source holds with no readable file; not a null): 232
- replicates in all three (significant in each, same direction): 1 (Ischemic Heart Disease)

| Harmonized outcome | Pan-UKBB EUR p | FinnGen p | MVP EUR p | Pan-UKBB EUR × FinnGen | Pan-UKBB EUR × MVP EUR | FinnGen × MVP EUR | All three |
|---|---|---|---|---|---|---|---|
| LDL direct | 1.58e-465 |  | 6.02e-587 |  | replicates (significant in both, same direction) |  |  |
| Ischemic Heart Disease | 2.14e-28 | 3.68e-24 | 1.39e-71 | replicates (significant in both, same direction) | replicates (significant in both, same direction) | replicates (significant in both, same direction) | replicates in all three (significant in each, same direction) |
| E78 Disorders of lipoprotein metabolism and other lipidaemias | 3.31e-62 | 3.23e-57 |  | replicates (significant in both, same direction) |  |  |  |

Full table: `tables/replication.tsv`.

Harmonization basis across all rows:

- no counterpart in another biobank: 7,338
- the outcome is a phecode: 1,191
- same phecode (case definitions may still differ by cohort): 1,155
- the outcome is an ICD-10 code block: 916
- matched by FinnGen's own published phenotype mapping: 546
- only a broader or narrower ICD-10 counterpart exists; not compared: 179
- a one-sex version of an outcome; not compared: 109
- same single ICD-10 code block (case definitions may still differ): 99
- same measurement by a reviewed name match (units and summaries may differ): 63
- the maximum or minimum of a measurement whose mean is the one compared; listed, not compared: 58
- only a near-equivalent measurement exists; not compared: 45

## Other ancestry groups

This lookup read the European-ancestry datasets. Other ancestry groups can be read with --ancestry: African (AFR), admixed American (AMR), East Asian (EAS), Central/South Asian (CSA) or Middle Eastern (MID), or META for the Million Veteran Program's multi-population meta-analysis. For example: `python skills/snp-phewas/snp_phewas.py --variant rs12740374 --ancestry AFR --output DIR`.

## What was read

- Pan-UK Biobank, European ancestry: release `pan_ukbb_2024_06_26` (GRCh37); route tabix; 3 files read; 0 s.
  - Status of every catalog outcome: outcome was not part of this lookup: 7,157, tested: 3.
- FinnGen R12 (Finnish): release `finngen_r12` (GRCh38); route api+tabix; 0 files read; 1 API request(s); 0 s.
  - Status of every catalog outcome: outcome was not part of this lookup: 2,466, tested: 3.
- Million Veteran Program, European genetic ancestry: release `mvp_gwphewas_gwascatalog_2026-10-19` (GRCh38); route tabix; 2 files read; 0 s.
  - Status of every catalog outcome: outcome was not part of this lookup: 1,219, the GWAS Catalog holds this study without a position index (or without summary statistics), so it cannot be read; not a null: 849, tested: 2.
- Tabix indexes: 0 downloaded (0.0 MB), 0 already cached.
- Wall time: 1.5 s.

## Data sources

- Pan-UK Biobank (Karczewski KJ et al., Nat Genet 2025;57:2408-2417, doi:10.1038/s41588-025-02335-7, PMID 40968291); this research has been conducted using the UK Biobank Resource.
- FinnGen R12 (Kurki MI et al., Nature 2023;613:508-518, doi:10.1038/s41586-022-05473-8, PMID 36653562). We want to acknowledge the participants and investigators of the FinnGen study.
- Million Veteran Program (Verma A et al., Science 2024;385:eadj1182, doi:10.1126/science.adj1182, PMID 39024449); summary statistics via the NHGRI-EBI GWAS Catalog.

## Files

- `tables/snp_phewas_rows.tsv`: one row per catalog outcome per dataset, every status.
- `tables/excluded.tsv`: the rows with no value, each with its reason.
- `figures/manhattan.png`: outcome Manhattan, one panel per dataset.
- `tables/replication.tsv`, `figures/replication_matrix.png`: the replication view (every dataset pair that does not share participants, and all three when three datasets are read and none of them do).
- `figures/paired_effects.png`: effect against effect, one row per dataset pair with shared outcomes.

---

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*

