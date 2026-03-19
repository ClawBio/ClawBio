# ClawBio PharmGx Report

**Date**: 2026-03-19 16:40 UTC
**Input**: `genotype_data.txt`
**Format detected**: 23andme
**Checksum (SHA-256)**: `40d540820d4a3dd4db039c800b9133b6fec8e6962a90a53ee66fdc33ca992cc7`
**Total SNPs in file**: 960861
**Pharmacogenomic SNPs found**: 25/30
**Genes profiled**: 12
**Drugs assessed**: 51

---

## DATA QUALITY WARNING

**1 gene(s) could not be assessed** because the relevant SNPs were not found in the input file: NUDT15

Drugs depending on these genes are marked INSUFFICIENT DATA below. Do not assume normal metabolism for untested genes.

**1 gene(s) have unmapped diplotypes**: CYP2C19. These diplotypes could not be matched to a known phenotype. Clinical pharmacogenomic testing is recommended.

---

## Drug Response Summary

| Category | Count |
|----------|-------|
| Standard dosing | 32 |
| Use with caution | 8 |
| Avoid / use alternative | 1 |
| Insufficient data | 10 |

### Actionable Alerts

**AVOID / USE ALTERNATIVE:**

- **Warfarin** (Coumadin) [CYP2C9+VKORC1]

**USE WITH CAUTION:**

- **Phenytoin** (Dilantin) [CYP2C9]
- **Celecoxib** (Celebrex) [CYP2C9]
- **Flurbiprofen** (Ansaid) [CYP2C9]
- **Piroxicam** (Feldene) [CYP2C9]
- **Meloxicam** (Mobic) [CYP2C9]
- **Tacrolimus** (Prograf) [CYP3A5]
- **Efavirenz** (Sustiva) [CYP2B6]
- **Clozapine** (Clozaril) [CYP1A2]

---

## Gene Profiles

| Gene | Full Name | Diplotype | Phenotype |
|------|-----------|-----------|-----------|
| CYP2C19 | Cytochrome P450 2C19 | *17/*4 | Unknown (unmapped diplotype: *17/*4) |
| CYP2D6 | Cytochrome P450 2D6 | *1/*1 (4/5 SNPs tested) | Normal Metabolizer |
| CYP2C9 | Cytochrome P450 2C9 | *1/*2 | Intermediate Metabolizer |
| VKORC1 | Vitamin K Epoxide Reductase | TT | High Warfarin Sensitivity |
| SLCO1B1 | Solute Carrier Organic Anion Transporter 1B1 | TT | Normal Function |
| DPYD | Dihydropyrimidine Dehydrogenase | Normal/Normal | Normal Metabolizer |
| TPMT | Thiopurine S-Methyltransferase | *1/*1 | Normal Metabolizer |
| UGT1A1 | UDP-Glucuronosyltransferase 1A1 | *1/*1 (1/2 SNPs tested) | Normal Metabolizer |
| CYP3A5 | Cytochrome P450 3A5 | *1/*1 | CYP3A5 Expressor |
| CYP2B6 | Cytochrome P450 2B6 | *9/*9 | Poor Metabolizer |
| NUDT15 | Nudix Hydrolase 15 | NOT_TESTED | Indeterminate (not genotyped) |
| CYP1A2 | Cytochrome P450 1A2 | *1F/*1F | Ultrarapid Metabolizer |

## Detected Variants

| rsID | Gene | Star Allele | Genotype | Effect |
|------|------|-------------|----------|--------|
| rs762551 | CYP1A2 | *1F | AA | increased_function |
| rs3745274 | CYP2B6 | *9 | TT | decreased_function |
| rs28399499 | CYP2B6 | *18 | TT | no_function |
| rs4244285 | CYP2C19 | *2 | GG | no_function |
| rs4986893 | CYP2C19 | *3 | GG | no_function |
| rs12248560 | CYP2C19 | *17 | CT | increased_function |
| rs28399504 | CYP2C19 | *4 | AG | no_function |
| rs1799853 | CYP2C9 | *2 | CT | decreased_function |
| rs1057910 | CYP2C9 | *3 | AA | decreased_function |
| rs3892097 | CYP2D6 | *4 | CC | no_function |
| rs5030655 | CYP2D6 | *6 | II | no_function |
| rs1065852 | CYP2D6 | *10 | GG | decreased_function |
| rs28371725 | CYP2D6 | *41 | CC | decreased_function |
| rs776746 | CYP3A5 | *3 | CC | no_function |
| rs10264272 | CYP3A5 | *6 | CC | no_function |
| rs41303343 | CYP3A5 | *7 | II | no_function |
| rs3918290 | DPYD | *2A | CC | no_function |
| rs55886062 | DPYD | *13 | AA | no_function |
| rs67376798 | DPYD | D949V | TT | decreased_function |
| rs4149056 | SLCO1B1 | *5 | TT | decreased_function |
| rs1800460 | TPMT | *3B | CC | no_function |
| rs1142345 | TPMT | *3C | TT | no_function |
| rs1800462 | TPMT | *2 | CC | no_function |
| rs4148323 | UGT1A1 | *6 | GG | decreased_function |
| rs9923231 | VKORC1 | -1639G>A | TT | decreased_expression |

---

## Complete Drug Recommendations

| Drug | Brand | Class | Gene | Status |
|------|-------|-------|------|--------|
| Warfarin | Coumadin | Anticoagulant | CYP2C9+VKORC1 | AVOID |
| Celecoxib | Celebrex | NSAID | CYP2C9 | CAUTION |
| Clozapine | Clozaril | Antipsychotic | CYP1A2 | CAUTION |
| Efavirenz | Sustiva | Antiretroviral | CYP2B6 | CAUTION |
| Flurbiprofen | Ansaid | NSAID | CYP2C9 | CAUTION |
| Meloxicam | Mobic | NSAID | CYP2C9 | CAUTION |
| Phenytoin | Dilantin | Antiepileptic | CYP2C9 | CAUTION |
| Piroxicam | Feldene | NSAID | CYP2C9 | CAUTION |
| Tacrolimus | Prograf | Immunosuppressant | CYP3A5 | CAUTION |
| Citalopram | Celexa | SSRI Antidepressant | CYP2C19 | INSUFFICIENT DATA |
| Clopidogrel | Plavix | Antiplatelet Agent | CYP2C19 | INSUFFICIENT DATA |
| Dexlansoprazole | Dexilant | Proton Pump Inhibitor | CYP2C19 | INSUFFICIENT DATA |
| Escitalopram | Lexapro | SSRI Antidepressant | CYP2C19 | INSUFFICIENT DATA |
| Esomeprazole | Nexium | Proton Pump Inhibitor | CYP2C19 | INSUFFICIENT DATA |
| Lansoprazole | Prevacid | Proton Pump Inhibitor | CYP2C19 | INSUFFICIENT DATA |
| Omeprazole | Prilosec | Proton Pump Inhibitor | CYP2C19 | INSUFFICIENT DATA |
| Pantoprazole | Protonix | Proton Pump Inhibitor | CYP2C19 | INSUFFICIENT DATA |
| Sertraline | Zoloft | SSRI Antidepressant | CYP2C19 | INSUFFICIENT DATA |
| Voriconazole | Vfend | Antifungal | CYP2C19 | INSUFFICIENT DATA |
| Amitriptyline | Elavil | Tricyclic Antidepressant | CYP2D6 | OK |
| Aripiprazole | Abilify | Antipsychotic | CYP2D6 | OK |
| Atazanavir | Reyataz | Antiretroviral | UGT1A1 | OK |
| Atomoxetine | Strattera | ADHD Medication | CYP2D6 | OK |
| Atorvastatin | Lipitor | Statin | SLCO1B1 | OK |
| Azathioprine | Imuran | Immunosuppressant | TPMT | OK |
| Capecitabine | Xeloda | Antineoplastic | DPYD | OK |
| Clomipramine | Anafranil | Tricyclic Antidepressant | CYP2D6 | OK |
| Codeine | Tylenol w/ Codeine | Opioid Analgesic | CYP2D6 | OK |
| Desipramine | Norpramin | Tricyclic Antidepressant | CYP2D6 | OK |
| Doxepin | Sinequan | Tricyclic Antidepressant | CYP2D6 | OK |
| Fluorouracil | 5-FU | Antineoplastic | DPYD | OK |
| Fluoxetine | Prozac | SSRI Antidepressant | CYP2D6 | OK |
| Haloperidol | Haldol | Antipsychotic | CYP2D6 | OK |
| Hydrocodone | Vicodin | Opioid Analgesic | CYP2D6 | OK |
| Imipramine | Tofranil | Tricyclic Antidepressant | CYP2D6 | OK |
| Irinotecan | Camptosar | Antineoplastic | UGT1A1 | OK |
| Mercaptopurine | Purinethol | Immunosuppressant | TPMT | OK |
| Metoprolol | Lopressor | Beta-Blocker | CYP2D6 | OK |
| Nortriptyline | Pamelor | Tricyclic Antidepressant | CYP2D6 | OK |
| Ondansetron | Zofran | Antiemetic | CYP2D6 | OK |
| Oxycodone | OxyContin | Opioid Analgesic | CYP2D6 | OK |
| Paroxetine | Paxil | SSRI Antidepressant | CYP2D6 | OK |
| Pravastatin | Pravachol | Statin | SLCO1B1 | OK |
| Risperidone | Risperdal | Antipsychotic | CYP2D6 | OK |
| Rosuvastatin | Crestor | Statin | SLCO1B1 | OK |
| Simvastatin | Zocor | Statin | SLCO1B1 | OK |
| Tamoxifen | Nolvadex | SERM (Oncology) | CYP2D6 | OK |
| Thioguanine | Tabloid | Immunosuppressant | TPMT | OK |
| Tramadol | Ultram | Opioid Analgesic | CYP2D6 | OK |
| Trimipramine | Surmontil | Tricyclic Antidepressant | CYP2D6 | OK |
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
python pharmgx_reporter.py --input genotype_data.txt --output report
```

**Input checksum**: `40d540820d4a3dd4db039c800b9133b6fec8e6962a90a53ee66fdc33ca992cc7`

## References

- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio
- CPIC. Clinical Pharmacogenetics Implementation Consortium. https://cpicpgx.org/
- Caudle, K.E. et al. (2014). Standardizing terms for clinical pharmacogenetic test results. Genet Med, 16(9), 655-663.
- PharmGKB. https://www.pharmgkb.org/
