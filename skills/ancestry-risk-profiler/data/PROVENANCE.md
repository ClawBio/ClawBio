# Data Provenance — ancestry-risk-profiler

## aisnp_panel.csv

**Purpose**: Panel of SNPs with per-allele frequencies across five 1000 Genomes
super-populations (AFR, AMR, EAS, EUR, SAS). Automatic ancestry inference requires
at least 30 matched SNPs with Wright Fst ≥ 0.3. The current bundled panel has five
markers above that floor, so it abstains unless the user supplies `--ancestry`.
Lower-Fst disease/PGx SNPs stay in the file and may contribute to likelihood after
the coverage gate is satisfied, but they cannot pad the high-Fst coverage count.

**Allele frequencies**: Sourced from gnomAD v3.1 population allele frequencies.
- Primary citation: Karczewski et al. (2020) Nature 581:434–443. PMID: 32461654
- Supplementary: 1000 Genomes Project Consortium (2015) Nature 526:68–74. PMID: 26432245

**SNP selection criteria (documented vs enforced)**: Earlier text implied the panel
met an AIM-style Fst threshold even though 67 of 72 rows are below it (T2D, CAD,
PGx SNPs whose frequencies barely differ across super-populations). From v1.4.0 the
runtime computes equal-weighted Wright Fst from the five frequency columns and
counts only Fst ≥ 0.3 markers toward the 30-marker automatic-inference floor. Fst
is a coverage guard from issue #313, not a claim that Nassir et al. or Kosoy et al.
defined or validated this exact software cutoff. The five markers that currently
clear the floor:

Reference context:
- Nassir et al. (2009) BMC Genetics 10:39. [PMID 19630973](https://pubmed.ncbi.nlm.nih.gov/19630973/).
  Evaluated a 93-SNP AIM panel for continental ancestry inference.
- Kosoy et al. (2009) Human Mutation 30(1):69–78. [PMID 18683858](https://pubmed.ncbi.nlm.nih.gov/18683858/).
  Described a 128-marker AIM set and smaller subsets down to 24 markers.

| rsID | Gene | Wright Fst |
|------|------|------------|
| rs1426654 | SLC24A5 | ~0.72 |
| rs2814778 | DARC/ACKR1 | ~0.72 |
| rs16891982 | SLC45A2 | ~0.48 |
| rs3827760 | EDAR | ~0.42 |
| rs4988235 | MCM6/LCT | ~0.32 |

HFE rs1800562 (Fst ~0.01) and ALDH2 rs671 are disease variants, not AIMs; they
remain in the CSV for provenance of what was considered, but they no longer
count toward coverage.

Key markers that actually meet the Fst floor, and their population-specificity:

| rsID | Gene | Specificity |
|------|------|-------------|
| rs1426654 | SLC24A5 | Strong EUR/SAS vs. AFR/EAS differentiator (skin pigmentation locus) |
| rs16891982 | SLC45A2 | EUR enriched |
| rs3827760 | EDAR | EAS enriched (hair follicle morphology) |
| rs4988235 | MCM6/LCT | EUR enriched (lactase persistence) |
| rs2814778 | DARC/ACKR1 | AFR enriched (Duffy antigen) |

**Changelog**:
- v1.1 — fixed trailing whitespace in rsid field of row 16 (rs35205) that caused silent panel miss.
- v1.3.2 — removed 7 near-zero-FST markers (all mean FST < 0.07 across 5 super-populations):
  rs731236, rs2228570 (VDR), rs1801133, rs1801131 (MTHFR C677T/A1298C), rs1800497 (ANKK1/DRD2),
  rs4680 (COMT Val158Met), rs53576 (OXTR). These are candidate-gene pharmacogenomic SNPs, not
  ancestry-informative markers — they provided no population-discriminatory power and were
  incorrectly included. Panel size: 79 → 72 markers. The four high-FST anchors already present
  (rs1426654 SLC24A5, rs16891982 SLC45A2, rs3827760 EDAR, rs2814778 DARC) are retained and
  continue to provide the dominant ancestry signal. The remaining 65 rows were still
  counted toward the 30-marker gate; that is the hole v1.4.0 closes.
- v1.4.0 — **Fst floor is now a runtime coverage property, not an rsID blocklist** (#313).
  Wright Fst is computed from the five super-population frequencies already in the CSV.
  Markers with Fst < 0.3 do not count toward coverage. The conservative software floor
  is back to 30 informative markers; the shipped panel has 5, so automatic inference
  abstains until the AIM panel is expanded. A 30-SNP padding panel of near-zero-Fst
  disease SNPs now abstains. The 7-rsID blocklist from v1.3.2 remains as a second check,
  but it is no longer the thing that enforces the Fst claim.

---

## ancestry_risk_associations.json

**Purpose**: Curated variant-disease-ancestry triples with ancestry-specific odds ratios (ORs)
and matched European reference ORs, enabling ancestry-vs-EUR comparisons via the Ancestry
Elevation Score.

**Coverage**: 5 super-populations × ~11 diseases = ~51 disease-ancestry combinations.
Full panel: ~51 additive association entries + 1 compound recessive group (APOL1 CKD/AFR).

### Changes in v1.3.1 (PR-audit follow-up: dead-data removal + doc reconciliation)

**Removed inert `baseline_prevalences` block from `ancestry_risk_associations.json`**

The top-level `baseline_prevalences` object (5 super-populations, incl. diseases not scored by the
panel such as Prostate Cancer) was no longer read by any code path after the absolute-risk
calculation was dropped in v1.3.0. It was removed entirely so the file cannot be re-wired into a
double-counting absolute-risk calculation, and a CI test (`test_no_baseline_prevalences_block`)
now asserts its absence. No scored output changes.

**rs671 (ALDH2 ESCC) PMID docs reconciled with the data**

The JSON entry carries `pmid: "22960999"` (Wu et al. 2012 Nat Genet), but SKILL.md and this file
still described it as `pmid: null` "pending re-verification". The docs are corrected to match the
data; the entry note continues to flag that the GWAS Catalog accession-to-PMID mapping (GCST001563)
was not independently confirmed.

Version strings aligned to 1.3.1 across SKILL.md, `ancestry_risk_associations.json`, and `catalog.json`.

### Changes in v1.3.0 (PR-audit provenance resolution)

**PMID 17478679 misidentified as "McPherson 2007" — corrected to Helgadottir 2007 and scope narrowed**

PMID 17478679 is Helgadottir A et al. Science 2007;316(5830):1491 — a CAD/9p21 study reporting
rs10757278/rs2383206, not rs1333049. It was cited for hypertension (rs4343/rs699) and
Parkinson's (rs1800497) rows where it is wholly irrelevant.

- rs1333049 CAD EUR → 17634449 (Samani et al. 2007 NEJM, WTCCC — explicitly reports rs1333049)
- rs1333049 CAD SAS → 24246088 (Bhanushali et al. 2013, proxy; OR mismatch noted in entry)
- rs11591147 PCSK9 CAD EUR → 16554528 (Cohen et al. 2006 NEJM, ARIC; HR=0.50 for R46L LOF)
- rs4343 HTN EUR → 20639399 (Niu et al. 2011 meta-analysis; candidate-gene era, no GWAS equivalent)
- rs699 HTN EUR → 1394429 (Jeunemaitre et al. 1992 Cell; seminal M235T linkage study)
- rs1800497 PD EUR → 21663922 (McGuire et al. 2011 PEGASUS; OR=1.5 in non-Hispanic whites)

**PMID 22158537 mislabelled as Kooner et al. (SAS T2D) — is actually Cho et al. (EAS T2D)**

PMID 22158537 is Cho YS et al. Nat Genet 2011;44(1):67, "Meta-analysis of GWAS identifies eight
new loci for T2D in east Asians." The PROVENANCE label was wrong; this is an EAS paper. Correct
Kooner SAS T2D GWAS: PMID 21874001 (Kooner JS et al. Nat Genet 2011;43(10):984).

All 11 entries previously cited under 22158537 re-assigned:
- 7 SAS T2D entries → 21874001; replication of these specific loci in the SAS cohort is pending
  PMC3773920 full-text verification (Kooner reports 6 novel loci; prior-known loci may appear as
  replication covariates in supplementary tables)
- 4 non-T2D SAS entries (HTN rs4343/rs699, CAD rs10455872, AD rs429358) → dropped; no
  SAS-specific replacement paper found for these phenotypes

**PMID 27005778 wholesale mismatch — metabolomics GWAS, not disease GWAS**

PMID 27005778 is Kettunen J et al. Nat Commun 2016;7:11122, "Genome-wide study for circulating
metabolites identifies 62 loci" — NMR metabolomics in ~25,000 mostly EUR/Finnish/Estonian
subjects. Not a disease-association study; not AFR ancestry. All 8 entries citing it were invalid.

- rs11591147 PCSK9 CAD AFR → 28768753 (PCSK9 LOF meta-analysis, 9 Black cohorts, OR=0.51)
- 7 other AFR disease entries → dropped; no replacement papers confirmed (see removed_v1.3.0)
  Notable: rs1800497 PD AFR had a direction conflict (McGuire 2011 shows protective in
  African-Americans, opposite to the stored OR 1.32 risk direction)

**ALDH2 rs671 EAS ESCC — null PMID resolved**

- rs671 ESCC EAS → 22960999 (Wu C et al. 2012 Nat Genet, Chinese Han ESCC GWAS confirming
  ALDH2 12q24 locus). GWAS Catalog accession GCST001563 to-PMID mapping not independently
  verified via Catalog page; flagged in entry note.

---

### Changes in v1.1.0 (domain re-audit — LCT removal + ALDH2 citation fix)

**LCT rs4988235 entries removed (all 4: EAS, SAS, AFR, EUR)**

Root cause: two compounding errors made these entries unshippable:

1. **Direction artifact**: The EUR entry encoded C as a *protective* allele (or=0.45, framing:
   "C is the non-persistence minority — most EUR carry T"). Non-EUR entries encoded C as a
   *risk* allele (ORs 3.2–6.8x, framing: "C causes intolerance"). The AES formula
   `log(OR_ancestry) − log(OR_EUR)` then divides a risk framing by a protective framing,
   manufacturing AES values of 7–15x that have no biological interpretation. The demo patient
   surfaced "Lactose Intolerance | AES 7.11 | N=1" as the top result.

2. **Fabricated PMID**: PMID 14507249 was cited as Enattah et al. (2002) Nat Genet 30:233
   (lactase persistence). PubMed lookup confirms 14507249 = Hashibe et al. (2003)
   (bladder cancer survival study) — an unrelated paper.

3. **Biology**: rs4988235 is near-fixed C/C in EAS populations; there is no meaningful T
   comparison group. An OR computed in this context is not interpretable as a standard GWAS
   risk allele effect.

All four LCT entries are removed rather than repaired. A future version may add lactose
intolerance via a per-population approach that avoids cross-framing comparisons.

**ALDH2 rs671 EAS ESCC — second PMID correction**

The v1.0.0 PMID 20686008 (anaesthesiology case report) was corrected in v1.0.x to 22561518.
The deep domain audit confirmed that 22561518 = Jin et al. (2012) vitiligo GWAS — also wrong.
Two consecutive wrong PMIDs for the same entry break the "no fabricated citations" claim.

Resolution: PMID set to `null`; `gwas_catalog_accession: "GCST001563"` added as the
traceable source. The OR (1.89) and biology (ALDH2 acetaldehyde accumulation → ESCC) are
well-established; only the PMID link to the primary GWAS paper was incorrect.
A verified PMID must be confirmed against PubMed / GWAS Catalog before reinstating.

### Citation table (current — v1.3.0)

| PMID | Confirmed study | Diseases / variants |
|------|----------------|---------------------|
| 1394429 | Jeunemaitre X et al. (1992) Cell 71:169 — M235T angiotensinogen linkage | HTN rs699 EUR (candidate-gene era) |
| 9068472 | Feder et al. (1996) Nat Genet — HFE hereditary haemochromatosis | HFE rs1800562 EUR |
| 16415884 | Grant et al. (2006) Nat Genet — TCF7L2 T2D discovery | T2D rs7903146 EUR |
| 16554528 | Cohen JC et al. (2006) NEJM 354:1264 — PCSK9 LOF and CHD | CAD rs11591147 EUR (HR=0.50 for R46L in ARIC whites) |
| 17463246 | Zeggini et al. (2007) Nat Genet — T2D replication | T2D rs7903146 AFR/EAS/EUR, rs13266634 EUR, rs7756992 EUR, rs5219 EUR |
| 17529973 | Easton et al. (2007) Nature — breast cancer EUR GWAS | Breast cancer rs2981582 EUR, rs3803662 EUR |
| 17603472 | Gudbjartsson et al. (2007) Nat Genet — atrial fibrillation | AF rs10033464 EUR, rs7193343 EUR |
| 17634449 | Samani NJ et al. (2007) NEJM 357:443 — WTCCC CAD GWAS | CAD rs1333049 EUR |
| 20190752 | Dubois et al. (2010) Nat Genet — celiac disease GWAS | Celiac rs2476601 (PTPN22 R620W) EUR |
| 20647424 | Genovese G et al. (2010) Science 329:841 — APOL1 variants and kidney disease in African Americans | CKD rs73885319 (G1), rs60910145 (G2) AFR; compound_recessive_groups APOL1_CKD |
| 20581827 | Voight et al. (2010) Nat Genet — T2D CENTD2/ARAP1 | T2D rs1552224 EUR |
| 20639399 | Niu W et al. (2011) J Renin Angiotensin Aldosterone Syst — ACE meta-analysis | HTN rs4343 EUR (meta-analysis proxy; candidate-gene era) |
| 21347282 | Schunkert et al. (2011) Nat Genet — CAD meta-analysis | CAD rs10455872 EUR |
| 21663922 | McGuire V et al. (2011) J Neurol Sci 307:22 — PEGASUS DRD2/PD | PD rs1800497 EUR |
| 21874001 | Kooner JS et al. (2011) Nat Genet 43:984 — T2D South Asian GWAS | T2D SAS entries (rs7903146, rs13266634, rs2237892, rs7756992, rs1552224, rs8042680, rs5219) — ⚠️ replication of these specific loci pending PMC3773920 full-text verification |
| 22159054 | Naj et al. (2011) Nat Genet — Alzheimer's GWAS | APOE rs429358 AFR/EUR |
| 22960999 | Wu C et al. (2012) Nat Genet — Chinese Han ESCC GWAS | ESCC rs671 EAS — GCST001563 accession-to-PMID mapping unverified via Catalog; flagged in entry |
| 18097733 | Miyake K et al. (2008) J Hum Genet — TCF7L2 T2D Japanese | T2D rs7903146 EAS (2,214 cases/1,873 controls, OR=1.48) |
| 18162508 | Omori S et al. (2008) Diabetes 57:791 — T2D Japanese multi-locus | T2D rs13266634 (SLC30A8), rs7756992 (CDKAL1), rs5219 (KCNJ11) EAS |
| 18711367 | Yasuda K et al. (2008) Nat Genet 40:1092 — KCNQ1 T2D EAS discovery | T2D rs2237892 EAS; original EAS discovery paper |
| 23486537 | Barzan D et al. (2013) Eur J Hum Genet 21:1286 — breast cancer Chinese vs German | Breast cancer rs2981582 (FGFR2), rs3803662 (TOX3) EAS (984 Chinese cases/2,206 controls) |
| 24246088 | Bhanushali AA et al. (2013) Genet Res 95:138 — 9p21 CAD Indian | CAD rs1333049 SAS (proxy; OR 1.564 vs stored 1.34 — mismatch noted) |
| 28768753 | (PCSK9 LOF meta-analysis 2017) Circ Cardiovasc Genet — Black cohorts | CAD rs11591147 AFR (pooled OR=0.51 in 9 Black cohorts) |

**Correction history**:

| Version | Affected entry | Old PMID | Problem | Resolution |
|---------|---------------|----------|---------|------------|
| v1.0.x | APOL1 G1/G2 AFR CKD | 23340282 | Electrochemiluminescence sensor paper | → 20566908 (Genovese 2010) |
| v1.0.x | PTPN22 EUR celiac | 17571276 | Cdk5 neurodegeneration paper | → 20190752 (Dubois 2010) |
| v1.0.x | ALDH2 EAS ESCC | 20686008 | Anaesthesiology case report | → 22561518 (STILL WRONG; see v1.1.0) |
| v1.1.0 | ALDH2 EAS ESCC | 22561518 | Jin 2012 vitiligo GWAS (wrong paper) | → null + GCST001563 |
| v1.1.0 | LCT rs4988235 all ancestries | 14507249 | Hashibe 2003 bladder cancer (wrong paper) + direction artifact | → **Entries removed** |
| v1.3.0 | rs1333049 CAD EUR | 17478679 | Helgadottir 2007 9p21 — wrong variant (reports rs10757278, not rs1333049) | → 17634449 (Samani 2007 NEJM) |
| v1.3.0 | rs1333049 CAD SAS | 17478679 | Same mismatch; also wrong that SAS coverage was assumed | → 24246088 (Bhanushali 2013 proxy) |
| v1.3.0 | rs11591147 CAD EUR | 17478679 | CAD paper but wrong variant reported | → 16554528 (Cohen 2006 NEJM PCSK9) |
| v1.3.0 | rs4343 HTN EUR | 17478679 | CAD paper, wrong phenotype | → 20639399 (meta-analysis proxy) |
| v1.3.0 | rs699 HTN EUR | 17478679 | CAD paper, wrong phenotype | → 1394429 (Jeunemaitre 1992 Cell) |
| v1.3.0 | rs1800497 PD EUR | 17478679 | CAD paper, wrong phenotype entirely | → 21663922 (McGuire 2011 PEGASUS) |
| v1.3.0 | rs11591147 CAD AFR | 27005778 | Kettunen 2016 metabolomics GWAS — wrong phenotype & ancestry | → 28768753 (PCSK9 LOF AFR meta-analysis) |
| v1.3.0 | 7 AFR disease entries | 27005778 | Kettunen 2016 metabolomics GWAS — wrong paper, no replacement | → **Entries removed** (see removed_v1.3.0) |
| v1.3.0 | 11 SAS entries | 22158537 | Cho YS 2011 EAS T2D paper, not Kooner SAS T2D as labelled | → 21874001 (T2D SAS entries); 4 non-T2D SAS entries **removed** |
| v1.3.0 | rs671 ESCC EAS | null | GCST001563 PMID unresolved for 2 versions | → 22960999 (Wu 2012 Nat Genet ESCC GWAS) |
| v1.3.1 | APOL1 G1/G2 AFR CKD (associations + compound_recessive_groups) | 20566908 | Tschiesner 2010 head-and-neck cancer disability survey | → 20647424 (Genovese et al. 2010 Science 329:841, APOL1 variants and kidney disease) |
| v1.3.1 | HFE rs1800562 EUR HH | additive model | Additive model gave heterozygotes spurious 9.8x, homozygotes ~96x | → recessive_compound model; one_allele_or=1.5 (carrier), two_allele_or=9.8 (homozygote) |
| v1.3.2 | rs7903146/rs13266634/rs2237892/rs7756992/rs5219 T2D EAS | 23945395 | Hara 2014 Hum Mol Genet (3-locus Japanese T2D GWAS), not Cho 2012 as labelled | → 18097733 (TCF7L2), 18162508 (SLC30A8/CDKAL1/KCNJ11), 18711367 (KCNQ1) |
| v1.3.2 | rs2981582/rs3803662 Breast Cancer EAS | 23945395 | Same wrong paper; T2D paper cannot source breast cancer ORs | → 23486537 (Barzan 2013, Chinese cohort) |
| v1.3.2 | rs1333049 CAD EAS | 23945395 | Wrong paper; EAS-specific PMID unconfirmed this session | → **Entry removed** |
| v1.3.2 | rs10033464 AF EAS | 23945395 | Wrong paper; rs10033464 not replicated in EAS (OR=1.08 p=0.55 in HK Chinese); EAS 4q25 signal is rs2200733 | → **Entry removed** |
| v1.3.2 | rs429358 AD EAS | 23945395 | Wrong paper; EAS-specific APOE AD PMID unresolved | → **Entry removed** |
| v1.3.2 | AISNP panel (7 markers) | panel inclusion | FST near zero (max 0.06) — candidate-gene SNPs, not AIMs; counted toward 30-marker gate | → **Markers removed**: rs731236, rs2228570 (VDR), rs1801133, rs1801131 (MTHFR), rs1800497 (ANKK1), rs4680 (COMT), rs53576 (OXTR). Panel: 79 → 72 markers |
| v1.4.0 | AISNP coverage gate | 30-marker count of any panel hit | 67 of 72 remaining rows still have Wright Fst < 0.3 (T2D/CAD/PGx) and still padded the gate; the 7-rsID blocklist did not encode the Fst property | → Runtime Fst floor 0.3 for coverage; min 30 informative markers; padding panel of 30 near-zero-Fst SNPs now abstains (#313) |

**Important caveats**:
- ORs are from individual published studies; they are not re-computed here
- EUR reference ORs (`eur_or`) are from matched analyses in the same or companion studies
- Variants are assumed independent within each disease (LD not modelled)
- APOL1 CKD uses a recessive compound model (see `compound_recessive_groups`) — per-allele log-additive OR is biologically wrong for this locus
- This panel is curated for educational demonstration; it is not a validated clinical risk calculator
- ALDH2 rs671: PMID 22960999 (Wu et al. 2012 Nat Genet ESCC GWAS); the GWAS Catalog accession-to-PMID mapping (GCST001563) is flagged for confirmation in the entry note

---

## compound_recessive_groups

**APOL1_CKD**: Two risk alleles (G1+G2, G1/G1, or G2/G2) → ~7x FSGS/CKD risk.
One risk allele → carrier status (OR ~1.3).
Source: Genovese et al. (2010) Science 329:841. PMID 20566908.

---

## demo_patient_south_asian.txt

**Purpose**: Synthetic 23andMe-format file for use with `--demo`. Does not represent any
real individual.

**Construction**: Designed to carry multiple T2D risk alleles (rs7903146 CT,
rs2237892 CT, rs1552224 AC, etc.) so the demo report surfaces a meaningful AES
signal for Type 2 Diabetes when run with `--demo --ancestry SAS`. It includes five
high-Fst AIMs, below the 30-marker floor required for automatic inference.

**Data guarantee**: No real patient data. All genotypes are synthetic.
