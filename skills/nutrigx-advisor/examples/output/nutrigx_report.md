# NutriGx Personalised Nutrition Report

**Generated**: 2026-02-27 12:25 UTC  
**Tool**: ClawBio NutriGx Advisor v0.1.0  
**Input**: `synthetic_patient.csv` (demo — fixed synthetic data)  

> **Disclaimer**: This report is for research and educational purposes only. It does not constitute medical advice. Consult a registered dietitian or clinical geneticist before making significant dietary changes or starting supplements.

## Executive Summary

Analysed **28** of **28** panel SNPs from your genetic data.

**1 nutrient domain(s) at elevated genetic risk:**
- Vitamin D (score: 6.72/10)

**8 domain(s) at moderate genetic risk:**
- Folate / B-Vitamins (score: 5.0/10)
- Omega-3 / LC-PUFA (score: 4.04/10)
- Vitamin A (Beta-carotene) (score: 5.0/10)
- Vitamin B6 (score: 5.0/10)
- Carbohydrate Metabolism (score: 5.0/10)
- Caffeine Metabolism (score: 5.0/10)
- Lactose Tolerance (score: 5.0/10)
- Antioxidant / Detox (score: 5.0/10)

---

## Nutrient Risk Score Overview

| Nutrient Domain | Score (0–10) | Risk Category |
|-----------------|:------------:|:-------------:|
| Vitamin D | 6.72 | 🔴 Elevated |
| Folate / B-Vitamins | 5.0 | 🟡 Moderate |
| Vitamin A (Beta-carotene) | 5.0 | 🟡 Moderate |
| Vitamin B6 | 5.0 | 🟡 Moderate |
| Carbohydrate Metabolism | 5.0 | 🟡 Moderate |
| Caffeine Metabolism | 5.0 | 🟡 Moderate |
| Lactose Tolerance | 5.0 | 🟡 Moderate |
| Antioxidant / Detox | 5.0 | 🟡 Moderate |
| Omega-3 / LC-PUFA | 4.04 | 🟡 Moderate |
| Alcohol Metabolism | 2.35 | 🟢 Low |
| Fat Metabolism | 1.57 | 🟢 Low |
| Vitamin C | 0.0 | 🟢 Low |

![NutriGx Nutrient Risk Profile](nutrigx_radar.png)

![Gene × Nutrient Risk Heatmap](nutrigx_heatmap.png)

---

## Detailed Findings by Nutrient Domain

### Vitamin D

**Risk Score**: 6.72/10 — **Elevated**  
**SNPs tested**: 3 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| VDR | rs2228570 | `TC` | 1/2 | Decreased Vitamin D Response |
| VDR | rs731236 | `CC` | 2/2 | Decreased Vitamin D Response |
| GC | rs4588 | `CA` | 1/2 | Decreased Vitamin D Binding |

**Recommendation**

> Genetic predisposition to low vitamin D binding efficiency. Test and maintain 25(OH)D >100 nmol/L. D3+K2 co-supplementation advisable.

---

### Folate / B-Vitamins

**Risk Score**: 5.0/10 — **Moderate**  
**SNPs tested**: 3 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| MTHFR | rs1801133 | `CT` | 1/2 | Decreased Folate Conversion |
| MTHFR | rs1801131 | `AC` | 1/2 | Decreased Folate Conversion |
| MTR | rs1805087 | `AG` | 1/2 | Increased Homocysteine Risk |

**Recommendation**

> Consider increasing dietary folate. If homocysteine is elevated, 5-MTHF (methylfolate) supplement preferred over folic acid.

---

### Vitamin A (Beta-carotene)

**Risk Score**: 5.0/10 — **Moderate**  
**SNPs tested**: 2 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| BCMO1 | rs7501331 | `CT` | 1/2 | Decreased Beta Carotene Conversion |
| BCMO1 | rs12934922 | `AT` | 1/2 | Decreased Beta Carotene Conversion |

**Recommendation**

> Reduced BCMO1 activity. Include preformed vitamin A (liver, eggs, dairy) alongside carotenoid-rich vegetables.

---

### Vitamin B6

**Risk Score**: 5.0/10 — **Moderate**  
**SNPs tested**: 1 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| ALPL | rs1256335 | `CT` | 1/2 | Decreased B6 Utilisation |

**Recommendation**

> No specific recommendation available.

---

### Carbohydrate Metabolism

**Risk Score**: 5.0/10 — **Moderate**  
**SNPs tested**: 2 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| FTO | rs9939609 | `AT` | 1/2 | Increased Fat Mass Energy Dysregulation |
| TCF7L2 | rs7903146 | `CT` | 1/2 | Impaired Carbohydrate Metabolism |

**Recommendation**

> Glycaemic index awareness recommended. Favour whole-grain, low-GI carbohydrates.

---

### Caffeine Metabolism

**Risk Score**: 5.0/10 — **Moderate**  
**SNPs tested**: 2 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| CYP1A2 | rs762551 | `CA` | 1/2 | Slow Caffeine Metabolism |
| AHR | rs4410790 | `CT` | 1/2 | Reduced Cyp1A2 Induction |

**Recommendation**

> Intermediate metaboliser. Caffeine ≤200 mg/day recommended.

---

### Lactose Tolerance

**Risk Score**: 5.0/10 — **Moderate**  
**SNPs tested**: 1 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| MCM6 | rs4988235 | `GA` | 1/2 | Lactase Non Persistence |

**Recommendation**

> Possible partial lactase non-persistence. Monitor tolerance to high-lactose foods.

---

### Antioxidant / Detox

**Risk Score**: 5.0/10 — **Moderate**  
**SNPs tested**: 4 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| SOD2 | rs4880 | `CT` | 1/2 | Decreased Mitochondrial Antioxidant |
| GPX1 | rs1050450 | `CT` | 1/2 | Decreased Glutathione Peroxidase |
| NQO1 | rs1800566 | `CT` | 1/2 | Decreased Coq10 Regeneration |
| COMT | rs4680 | `GA` | 1/2 | Altered Methylation Catechol Metabolism |

**Recommendation**

> Some oxidative stress pathway variants. Increase dietary antioxidants; selenium and CoQ10 food sources beneficial.

---

### Omega-3 / LC-PUFA

**Risk Score**: 4.04/10 — **Moderate**  
**SNPs tested**: 4 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| FADS1 | rs174546 | `TC` | 1/2 | Decreased Lc Pufa Synthesis |
| FADS2 | rs1535 | `AG` | 1/2 | Altered Omega6 Omega3 Ratio |
| ELOVL2 | rs953413 | `AG` | 1/2 | Decreased Dha Synthesis |
| APOE | rs429358 | `TT` | 0/2 | Increased Ldl On Saturated Fat |

**Recommendation**

> Reduced LC-PUFA synthesis capacity. Increase direct EPA/DHA sources (oily fish, algae oil). Consider omega-6:omega-3 ratio <4:1.

---

### Alcohol Metabolism

**Risk Score**: 2.35/10 — **Low**  
**SNPs tested**: 2 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| ADH1B | rs1229984 | `AG` | 1/2 | Altered Acetaldehyde Accumulation |
| ALDH2 | rs671 | `GG` | 0/2 | Impaired Acetaldehyde Clearance |

**Recommendation**

> Standard alcohol metabolism. General guidelines apply.

---

### Fat Metabolism

**Risk Score**: 1.57/10 — **Low**  
**SNPs tested**: 3 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| PPARG | rs1801282 | `CG` | 1/2 | Reduced Fat Oxidation |
| APOA5 | rs662799 | `TT` | 0/2 | Elevated Triglycerides On Fat |
| APOE | rs7412 | `CC` | 0/2 | Apoe Epsilon2 Lower Ldl |

**Recommendation**

> No specific fat quality concerns from genetics alone.

---

### Vitamin C

**Risk Score**: 0.0/10 — **Low**  
**SNPs tested**: 1 | **Not on chip**: 0

| Gene | rsID | Genotype | Risk Alleles | Effect |
|------|------|----------|:------------:|--------|
| SLC23A1 | rs33972313 | `CC` | 0/2 | Decreased Vitamin C Reabsorption |

**Recommendation**

> Standard dietary vitamin C (citrus, peppers, kiwi) sufficient.

---

## Supplement Interaction Notes

| Supplement | Relevant Genes | Caution |
|------------|---------------|---------|
| Folic acid (synthetic) | MTHFR | Prefer 5-MTHF if MTHFR variant present |
| Vitamin D3 | VDR, GC | K2 co-supplementation improves utilisation |
| Fish oil / EPA+DHA | FADS1, FADS2, ELOVL2 | Higher dose needed if FADS variants present |
| Vitamin A (retinol) | BCMO1 | Do not exceed UL (3000 μg RAE/day) |
| Iron | HFE (not in panel) | Not assessed; request separate iron panel |
| CoQ10 | NQO1 | Ubiquinol form preferred if NQO1 P187S homozygous |

---

## Reproducibility

This report was generated deterministically. See `commands.sh` and `environment.yml`
in the output directory to reproduce this analysis on any machine.

## References

SNP-nutrient associations sourced from GWAS Catalog, ClinVar, and CPIC guidelines.
Full citations available in `skills/nutrigx-advisor/SKILL.md`.
