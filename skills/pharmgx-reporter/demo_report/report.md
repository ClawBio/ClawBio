# ClawBio PharmGx Report

**Date**: 2026-04-11 19:27 UTC
**Input**: `demo_patient.txt`
**Format detected**: 23andme
**Checksum (SHA-256)**: `ffe44b340edfbb21abf648f00c2ce68715f5c9453a590caef753bc25e316c5cc`
**Total SNPs in file**: 23
**Pharmacogenomic SNPs found**: 23/32
**Genes profiled**: 13
**Drugs assessed**: 59

---

## DATA QUALITY WARNING

**5 gene(s) have unmapped diplotypes**: CYP2C19, DPYD, TPMT, UGT1A1, NUDT15. These diplotypes could not be matched to a known phenotype. Clinical pharmacogenomic testing is recommended.

**Gene-Specific Limitations:**

- CYP2D6: Copy number variation (gene deletion CYP2D6*5, duplication CYP2D6*1xN/*2xN, hybrid alleles CYP2D6*13/*36) cannot be detected from DTC genotyping data. Phenotype assignment may be incomplete. Clinical-grade CYP2D6 testing includes CNV analysis.
- CYP3A5: The CYP3A5*7 allele (rs41303343) is an insertion/deletion variant. DTC genotyping platforms use proxy SNP calls that may not accurately detect this allele.
- UGT1A1: The UGT1A1*28 allele (rs8175347) is a TA-repeat polymorphism that cannot be reliably genotyped by SNP arrays. If rs8175347 is present in the input, the reported genotype is a proxy SNP call and may not reflect true TA repeat count.

---

## Panel Limitations

This report uses **SNP-based genotyping only** from a panel of 32 pharmacogenomic variants. The following cannot be detected:

- **Copy number variants (CNVs)**: Gene deletions (e.g. CYP2D6\*5) and duplications (e.g. CYP2D6\*1xN, \*2xN)
- **Structural variants**: CYP2D6-CYP2D7 hybrid alleles (e.g. \*13, \*36)
- **Repeat polymorphisms**: UGT1A1\*28 (TA7 repeat in rs8175347)
- **HLA typing**: HLA-B\*57:01 (abacavir hypersensitivity)
- **Mitochondrial variants**: MT-RNR1 m.1555A>G (aminoglycoside ototoxicity)
- **G6PD deficiency**: G6PD A- and Mediterranean variants

For CYP2D6, a result of 'Normal Metabolizer' does NOT rule out gene deletions or duplications that alter metabolizer status. Clinical-grade CYP2D6 testing includes CNV analysis.

## Genes Not Assessed by This Panel

The following clinically relevant pharmacogenomic genes are **not included** in this panel. Consider targeted testing if indicated:

| Gene | Clinical Relevance | CPIC Guideline |
|------|-------------------|----------------|
| HLA-B\*57:01 | Abacavir hypersensitivity | CPIC Abacavir (2014, updated 2020) |
| HLA-B\*58:01 | Allopurinol hypersensitivity (SJS/TEN) | CPIC Allopurinol (Hershfield 2013, updated 2015) |
| G6PD | Rasburicase, chloroquine contraindication | CPIC G6PD (2014) |
| MT-RNR1 | Aminoglycoside ototoxicity | FDA label |
| HLA-A\*31:01 | Carbamazepine hypersensitivity | CPIC Carbamazepine (2017) |
| CYP2B6 | Efavirenz metabolism | CPIC Efavirenz (2019) |

---

## Drug Response Summary

| Category | Count |
|----------|-------|
| Standard dosing | 17 |
| Use with caution | 24 |
| Avoid / use alternative | 1 |
| Insufficient data | 17 |

### Actionable Alerts

**AVOID / USE ALTERNATIVE:**

- **Warfarin** (Coumadin) [CYP2C9+VKORC1]

**USE WITH CAUTION:**

- **Codeine** (Tylenol w/ Codeine) [CYP2D6]
- **Tramadol** (Ultram) [CYP2D6]
- **Hydrocodone** (Vicodin) [CYP2D6]
- **Tamoxifen** (Nolvadex) [CYP2D6]
- **Amitriptyline** (Elavil) [CYP2D6]
- **Nortriptyline** (Pamelor) [CYP2D6]
- **Desipramine** (Norpramin) [CYP2D6]
- **Imipramine** (Tofranil) [CYP2D6]
- **Doxepin** (Sinequan) [CYP2D6]
- **Trimipramine** (Surmontil) [CYP2D6]
- **Clomipramine** (Anafranil) [CYP2D6]
- **Risperidone** (Risperdal) [CYP2D6]
- **Haloperidol** (Haldol) [CYP2D6]
- **Phenytoin** (Dilantin) [CYP2C9]
- **Celecoxib** (Celebrex) [CYP2C9]
- **Flurbiprofen** (Ansaid) [CYP2C9]
- **Piroxicam** (Feldene) [CYP2C9]
- **Meloxicam** (Mobic) [CYP2C9]
- **Efavirenz** (Sustiva) [CYP2B6]
- **Methotrexate** (Rheumatrex / Trexall) [MTHFR]
- **Ibuprofen** (Advil / Motrin) [CYP2C9]
- **Naproxen** (Aleve / Naprosyn) [CYP2C9]
- **Propafenone** (Rythmol) [CYP2D6]
- **Methadone** (Dolophine) [CYP2B6]

---

## Gene Profiles

| Gene | Full Name | Diplotype | Phenotype |
|------|-----------|-----------|-----------|
| CYP2C19 | Cytochrome P450 2C19 | Indeterminate (phase ambiguity: *2(rs4244285) + *17(rs12248560)) | Indeterminate (phase ambiguity: *2(rs4244285) + *17(rs12248560)) |
| CYP2D6 | Cytochrome P450 2D6 | *2/*41 | Intermediate Metabolizer |
| CYP2C9 | Cytochrome P450 2C9 | *1/*2 | Intermediate Metabolizer |
| VKORC1 | Vitamin K Epoxide Reductase | TT | High Warfarin Sensitivity |
| SLCO1B1 | Solute Carrier Organic Anion Transporter 1B1 | TT | Normal Function |
| DPYD | Dihydropyrimidine Dehydrogenase | Normal/Normal (2/3 SNPs tested) | Indeterminate (incomplete coverage: 2/3 SNPs tested) |
| TPMT | Thiopurine S-Methyltransferase | *1/*1 (2/3 SNPs tested) | Indeterminate (incomplete coverage: 2/3 SNPs tested) |
| UGT1A1 | UDP-Glucuronosyltransferase 1A1 | Indeterminate (rs8175347 not tested) | Indeterminate (rs8175347 not tested) |
| CYP3A5 | Cytochrome P450 3A5 | *3/*3 | CYP3A5 Non-expressor |
| CYP2B6 | Cytochrome P450 2B6 | *1/*9 | Intermediate Metabolizer |
| NUDT15 | Nudix Hydrolase 15 | *1/*1 (1/2 SNPs tested) | Indeterminate (incomplete coverage: 1/2 SNPs tested) |
| CYP1A2 | Cytochrome P450 1A2 | *1/*1F | Normal Metabolizer |
| MTHFR | Methylenetetrahydrofolate Reductase | 677CC/1298AC | Reduced MTHFR enzyme activity (677CT) |
| HLA-B | HLA-B*57:01 | Not assessed | Indeterminate (not in panel) |
| MT-RNR1 | MT-RNR1 (mitochondrial) | Not assessed | Indeterminate (not in panel) |
| G6PD | Glucose-6-Phosphate Dehydrogenase | Not assessed | Indeterminate (not in panel) |
| HLA-A | HLA-A*31:01 | Not assessed | Indeterminate (not in panel) |

## Interactive ClinPGx Link

**Explore this patient's pharmacogenomic profile** with CPIC, DPWG, and FDA guideline annotations:

[Open in ClinPGx](https://www.clinpgx.org/genotypeResults?q=%7B%22CYP2D6%22%3A%20%5B%22%2A2%22%2C%20%22%2A41%22%2C%20null%5D%2C%20%22CYP2C9%22%3A%20%5B%22%2A1%22%2C%20%22%2A2%22%2C%20null%5D%2C%20%22CYP3A5%22%3A%20%5B%22%2A3%22%2C%20%22%2A3%22%2C%20null%5D%2C%20%22CYP2B6%22%3A%20%5B%22%2A1%22%2C%20%22%2A9%22%2C%20null%5D%2C%20%22CYP1A2%22%3A%20%5B%22%2A1%22%2C%20%22%2A1F%22%2C%20null%5D%2C%20%22MTHFR%22%3A%20%5B%22677CC%22%2C%20%221298AC%22%2C%20null%5D%7D)

This link pre-fills ClinPGx with the patient's genotypes. No data is uploaded; the genotypes are encoded in the URL itself.

### CPIC Gene References

| Gene | Diplotype | Phenotype | CPIC Guideline |
|------|-----------|-----------|----------------|
| CYP2C19 | Indeterminate (phase ambiguity: *2(rs4244285) + *17(rs12248560)) | Indeterminate (phase ambiguity: *2(rs4244285) + *17(rs12248560)) | [CYP2C19 on CPIC](https://cpicpgx.org/gene/cyp2c19/) |
| CYP2D6 | *2/*41 | Intermediate Metabolizer | [CYP2D6 on CPIC](https://cpicpgx.org/gene/cyp2d6/) |
| CYP2C9 | *1/*2 | Intermediate Metabolizer | [CYP2C9 on CPIC](https://cpicpgx.org/gene/cyp2c9/) |
| VKORC1 | TT | High Warfarin Sensitivity | [VKORC1 on CPIC](https://cpicpgx.org/gene/vkorc1/) |
| SLCO1B1 | TT | Normal Function | [SLCO1B1 on CPIC](https://cpicpgx.org/gene/slco1b1/) |
| DPYD | Normal/Normal (2/3 SNPs tested) | Indeterminate (incomplete coverage: 2/3 SNPs tested) | [DPYD on CPIC](https://cpicpgx.org/gene/dpyd/) |
| TPMT | *1/*1 (2/3 SNPs tested) | Indeterminate (incomplete coverage: 2/3 SNPs tested) | [TPMT on CPIC](https://cpicpgx.org/gene/tpmt/) |
| UGT1A1 | Indeterminate (rs8175347 not tested) | Indeterminate (rs8175347 not tested) | [UGT1A1 on CPIC](https://cpicpgx.org/gene/ugt1a1/) |
| CYP3A5 | *3/*3 | CYP3A5 Non-expressor | [CYP3A5 on CPIC](https://cpicpgx.org/gene/cyp3a5/) |
| CYP2B6 | *1/*9 | Intermediate Metabolizer | [CYP2B6 on CPIC](https://cpicpgx.org/gene/cyp2b6/) |
| NUDT15 | *1/*1 (1/2 SNPs tested) | Indeterminate (incomplete coverage: 1/2 SNPs tested) | [NUDT15 on CPIC](https://cpicpgx.org/gene/nudt15/) |
| CYP1A2 | *1/*1F | Normal Metabolizer | [CYP1A2 on CPIC](https://cpicpgx.org/gene/cyp1a2/) |
| MTHFR | 677CC/1298AC | Reduced MTHFR enzyme activity (677CT) | [MTHFR on CPIC](https://cpicpgx.org/gene/mthfr/) |

## Detected Variants

| rsID | Gene | Star Allele | Genotype | Effect |
|------|------|-------------|----------|--------|
| rs762551 | CYP1A2 | *1F | AC | increased_function |
| rs2069514 | CYP1A2 | *1C | GG | decreased_function |
| rs3745274 | CYP2B6 | *9 | GT | decreased_function |
| rs4244285 | CYP2C19 | *2 | AG | no_function |
| rs4986893 | CYP2C19 | *3 | GG | no_function |
| rs12248560 | CYP2C19 | *17 | CT | increased_function |
| rs1799853 | CYP2C9 | *2 | CT | decreased_function |
| rs1057910 | CYP2C9 | *3 | AA | decreased_function |
| rs3892097 | CYP2D6 | *4 | CC | no_function |
| rs16947 | CYP2D6 | *2 | AG | normal_function |
| rs1065852 | CYP2D6 | *10 | CC | decreased_function |
| rs28371725 | CYP2D6 | *41 | CT | decreased_function |
| rs776746 | CYP3A5 | *3 | GG | no_function |
| rs3918290 | DPYD | *2A | CC | no_function |
| rs67376798 | DPYD | D949V | TT | decreased_function |
| rs1801133 | MTHFR | 677T | GG | decreased_function |
| rs1801131 | MTHFR | 1298C | GT | decreased_function |
| rs116855232 | NUDT15 | *3 | CC | no_function |
| rs4149056 | SLCO1B1 | *5 | TT | decreased_function |
| rs1800460 | TPMT | *3B | CC | no_function |
| rs1142345 | TPMT | *3C | AA | no_function |
| rs4148323 | UGT1A1 | *6 | GG | decreased_function |
| rs9923231 | VKORC1 | -1639G>A | TT | decreased_expression |

---

## Complete Drug Recommendations

| Drug | Brand | Class | Gene | Status |
|------|-------|-------|------|--------|
| Warfarin | Coumadin | Anticoagulant | CYP2C9+VKORC1 | AVOID |
| Amitriptyline | Elavil | Tricyclic Antidepressant | CYP2D6 | CAUTION |
| Celecoxib | Celebrex | NSAID | CYP2C9 | CAUTION |
| Clomipramine | Anafranil | Tricyclic Antidepressant | CYP2D6 | CAUTION |
| Codeine | Tylenol w/ Codeine | Opioid Analgesic | CYP2D6 | CAUTION |
| Desipramine | Norpramin | Tricyclic Antidepressant | CYP2D6 | CAUTION |
| Doxepin | Sinequan | Tricyclic Antidepressant | CYP2D6 | CAUTION |
| Efavirenz | Sustiva | Antiretroviral | CYP2B6 | CAUTION |
| Flurbiprofen | Ansaid | NSAID | CYP2C9 | CAUTION |
| Haloperidol | Haldol | Antipsychotic | CYP2D6 | CAUTION |
| Hydrocodone | Vicodin | Opioid Analgesic | CYP2D6 | CAUTION |
| Ibuprofen | Advil / Motrin | NSAID | CYP2C9 | CAUTION |
| Imipramine | Tofranil | Tricyclic Antidepressant | CYP2D6 | CAUTION |
| Meloxicam | Mobic | NSAID | CYP2C9 | CAUTION |
| Methadone | Dolophine | Opioid Analgesic | CYP2B6 | CAUTION |
| Methotrexate | Rheumatrex / Trexall | DMARD / Antineoplastic | MTHFR | CAUTION |
| Naproxen | Aleve / Naprosyn | NSAID | CYP2C9 | CAUTION |
| Nortriptyline | Pamelor | Tricyclic Antidepressant | CYP2D6 | CAUTION |
| Phenytoin | Dilantin | Antiepileptic | CYP2C9 | CAUTION |
| Piroxicam | Feldene | NSAID | CYP2C9 | CAUTION |
| Propafenone | Rythmol | Antiarrhythmic | CYP2D6 | CAUTION |
| Risperidone | Risperdal | Antipsychotic | CYP2D6 | CAUTION |
| Tamoxifen | Nolvadex | SERM (Oncology) | CYP2D6 | CAUTION |
| Tramadol | Ultram | Opioid Analgesic | CYP2D6 | CAUTION |
| Trimipramine | Surmontil | Tricyclic Antidepressant | CYP2D6 | CAUTION |
| Atazanavir | Reyataz | Antiretroviral | UGT1A1 | INDETERMINATE — INSUFFICIENT DATA |
| Azathioprine | Imuran | Immunosuppressant | TPMT | INDETERMINATE — INSUFFICIENT DATA |
| Capecitabine | Xeloda | Antineoplastic | DPYD | INDETERMINATE — INSUFFICIENT DATA |
| Citalopram | Celexa | SSRI Antidepressant | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Clopidogrel | Plavix | Antiplatelet Agent | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Dexlansoprazole | Dexilant | Proton Pump Inhibitor | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Escitalopram | Lexapro | SSRI Antidepressant | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Esomeprazole | Nexium | Proton Pump Inhibitor | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Fluorouracil | 5-FU | Antineoplastic | DPYD | INDETERMINATE — INSUFFICIENT DATA |
| Irinotecan | Camptosar | Antineoplastic | UGT1A1 | INDETERMINATE — INSUFFICIENT DATA |
| Lansoprazole | Prevacid | Proton Pump Inhibitor | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Mercaptopurine | Purinethol | Immunosuppressant | TPMT | INDETERMINATE — INSUFFICIENT DATA |
| Omeprazole | Prilosec | Proton Pump Inhibitor | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Pantoprazole | Protonix | Proton Pump Inhibitor | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Sertraline | Zoloft | SSRI Antidepressant | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Thioguanine | Tabloid | Immunosuppressant | TPMT | INDETERMINATE — INSUFFICIENT DATA |
| Voriconazole | Vfend | Antifungal | CYP2C19 | INDETERMINATE — INSUFFICIENT DATA |
| Aripiprazole | Abilify | Antipsychotic | CYP2D6 | OK |
| Atomoxetine | Strattera | ADHD Medication | CYP2D6 | OK |
| Atorvastatin | Lipitor | Statin | SLCO1B1 | OK |
| Bupropion | Wellbutrin / Zyban | Antidepressant / Smoking Cessation | CYP2B6 | OK |
| Clozapine | Clozaril | Antipsychotic | CYP1A2 | OK |
| Dextromethorphan | Robitussin DM | Antitussive | CYP2D6 | OK |
| Diclofenac | Voltaren | NSAID | CYP2C9 | OK |
| Fluoxetine | Prozac | SSRI Antidepressant | CYP2D6 | OK |
| Metoprolol | Lopressor | Beta-Blocker | CYP2D6 | OK |
| Ondansetron | Zofran | Antiemetic | CYP2D6 | OK |
| Oxycodone | OxyContin | Opioid Analgesic | CYP2D6 | OK |
| Paroxetine | Paxil | SSRI Antidepressant | CYP2D6 | OK |
| Pravastatin | Pravachol | Statin | SLCO1B1 | OK |
| Rosuvastatin | Crestor | Statin | SLCO1B1 | OK |
| Simvastatin | Zocor | Statin | SLCO1B1 | OK |
| Tacrolimus | Prograf | Immunosuppressant | CYP3A5 | OK |
| Venlafaxine | Effexor | SNRI Antidepressant | CYP2D6 | OK |

---

## Disclaimer

This report is for **research and educational purposes only**. It is NOT a diagnostic device and should NOT be used to make medication decisions without consulting a qualified healthcare professional.

Pharmacogenomic recommendations are based on CPIC guidelines (cpicpgx.org). DTC genetic tests have limitations: they may not detect all relevant variants, and results should be confirmed by clinical-grade testing before clinical use.

## Methods

- **Tool**: ClawBio PharmGx Reporter v0.2.0
- **SNP panel**: 31 pharmacogenomic variants across 12 genes
- **Star allele calling**: Simplified DTC-compatible algorithm (single-SNP per allele)
- **Phenotype assignment**: CPIC-based diplotype-to-phenotype mapping
- **Drug guidelines**: 51 drugs from CPIC (cpicpgx.org), simplified for DTC context

## Reproducibility

```bash
python pharmgx_reporter.py --input demo_patient.txt --output report
```

**Input checksum**: `ffe44b340edfbb21abf648f00c2ce68715f5c9453a590caef753bc25e316c5cc`

## References

- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio
- CPIC. Clinical Pharmacogenetics Implementation Consortium. https://cpicpgx.org/
- Caudle, K.E. et al. (2014). Standardizing terms for clinical pharmacogenetic test results. Genet Med, 16(9), 655-663.
- PharmGKB. https://www.pharmgkb.org/
- ClinPGx. Clinical Pharmacogenomics. https://www.clinpgx.org/
