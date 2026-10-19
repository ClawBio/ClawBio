<!-- Output of: python clawbio.py run snp-phewas --demo --output <dir>  (report.md; date and wall time vary) -->

# Variant PheWAS across biobanks: rs12740374

**Date**: 2026-10-19 07:22 UTC
**Skill**: snp-phewas
**Version**: 0.1.0
**Variant (GRCh38)**: 1_109274968_G_T
**Effect allele**: T (other allele G)
**Datasets**: Pan-UK Biobank, European ancestry, FinnGen R12 (Finnish)

---

## What this answer covers

- **Demo data.** This run reads bundled slices of 3 Pan-UK Biobank outcome files and 3 FinnGen endpoint files around rs12740374, offline; the rsID and the GRCh37 position use answers recorded from Open Targets and Ensembl. Every other outcome is listed as not part of this lookup.

## Summary

- **Pan-UK Biobank, European ancestry**: tested in 3 of 3 outcomes read (7,160 in the catalog); 3 below the threshold 0.05 / 3 (over the 3 outcomes looked up).
- **FinnGen R12 (Finnish)**: tested in 3 of 3 outcomes read (2,469 in the catalog); 3 below the threshold 0.05 / 3 (over the 3 outcomes looked up).
- **Replication** (Pan-UK Biobank, European ancestry vs FinnGen R12 (Finnish), threshold 0.05 / 589 shared outcomes): 2 replicate, 0 same direction significant in one only, 0 opposite directions, 0 null in both.

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
- Million Veteran Program (dbGaP phs002453) is not included: its results are openly downloadable but carry no license text. BioBank Japan (NBDC hum0197) is not included: its terms limit use to research and development purposes.

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

## Replication across biobanks (589 harmonized outcomes)

Significance threshold: 0.05 / 589 outcomes both datasets share. Full table: `tables/replication.tsv`.

- replicates (significant in both, same direction): 2
- not read in one of the datasets (an outage, or outside this lookup; not a null): 587

| Harmonized outcome | Pan-UK Biobank, European ancestry p | FinnGen R12 (Finnish) p | Verdict |
|---|---|---|---|
| E78 Disorders of lipoprotein metabolism and other lipidaemias | 3.31e-62 | 3.23e-57 | replicates (significant in both, same direction) |
| Ischemic Heart Disease | 2.14e-28 | 3.68e-24 | replicates (significant in both, same direction) |

Harmonization basis across all rows:

- no counterpart in the other biobank: 7,647
- the outcome is an ICD-10 code block: 916
- matched by FinnGen's own published phenotype mapping: 788
- only a broader or narrower ICD-10 counterpart exists; not compared: 179
- same single ICD-10 code block (case definitions may still differ): 99

## What was read

- Pan-UK Biobank, European ancestry: release `pan_ukbb_2024_06_26` (GRCh37); route tabix; 3 files read; 0 s.
  - Status of every catalog outcome: outcome was not part of this lookup: 7,157, tested: 3.
- FinnGen R12 (Finnish): release `finngen_r12` (GRCh38); route api+tabix; 0 files read; 1 API request(s); 0 s.
  - Status of every catalog outcome: outcome was not part of this lookup: 2,466, tested: 3.
- Tabix indexes: 0 downloaded (0.0 MB), 0 already cached.
- Wall time: 1.1 s.

## Data sources

- Pan-UK Biobank (Karczewski KJ et al., Nat Genet 2025;57:2408-2417, doi:10.1038/s41588-025-02335-7, PMID 40968291); this research has been conducted using the UK Biobank Resource.
- FinnGen R12 (Kurki MI et al., Nature 2023;613:508-518, doi:10.1038/s41586-022-05473-8, PMID 36653562). We want to acknowledge the participants and investigators of the FinnGen study.

## Files

- `tables/snp_phewas_rows.tsv`: one row per catalog outcome per dataset, every status.
- `tables/excluded.tsv`: the rows with no value, each with its reason.
- `figures/manhattan.png`: outcome Manhattan, one panel per dataset.
- `tables/replication.tsv`, `figures/replication_matrix.png`: the replication view.

---

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*

