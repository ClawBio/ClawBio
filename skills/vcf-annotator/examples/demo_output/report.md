# 🦖 ClawBio VCF Annotator Report

**Input**: demo_variants.vcf  
**Date**: 2026-04-18 19:07 UTC  
**Total variants**: 5  
**HIGH impact**: 2 | **MODERATE**: 3 | **LOW**: 0<br>
**ClinVar Pathogenic/Likely Pathogenic**: 3

---

## Summary

Annotated 5 variants against Ensembl VEP, ClinVar, and gnomAD. 2 variant(s) carry HIGH impact consequences. 3 variant(s) are classified as Pathogenic or Likely Pathogenic in ClinVar. Variants are ranked by predicted functional impact below.

---

## Variant Table

| # | Gene | Variant | Consequence | Impact | ClinVar | gnomAD AF |
|---|------|---------|-------------|--------|---------|-----------|
| 1 | BRCA2 | 13:32316461 C>T | stop_gained | **HIGH** | Pathogenic | 0.000004 |
| 2 | CFTR | 7:117548628 CTTT>C | frameshift_variant | **HIGH** | Pathogenic | 0.021000 |
| 3 | BRCA1 | 17:43063931 G>A | missense_variant | **MODERATE** | Pathogenic | 0.000024 |
| 4 | APOE | 19:44908684 T>C | missense_variant | **MODERATE** | Risk factor | 0.147000 |
| 5 | MTHFR | 1:11796321 G>A | missense_variant | **MODERATE** | Benign/Likely benign | 0.312000 |

---

## Detailed Annotations

### 1. BRCA2 — 13:32316461 C>T

| Field | Value |
|-------|-------|
| **rsID** | rs80359550 |
| **HGVS** | NM_000059.4:c.9976A>T |
| **Consequence** | stop_gained |
| **Impact** | HIGH |
| **SIFT** | deleterious |
| **PolyPhen** | probably_damaging |
| **ClinVar** | Pathogenic |
| **ClinVar Condition** | Familial cancer of breast |
| **gnomAD AF (global)** | 0.000004 |
| **gnomAD AF (AFR)** | 0.000000 |
| **gnomAD AF (EUR)** | 0.000006 |

### 2. CFTR — 7:117548628 CTTT>C

| Field | Value |
|-------|-------|
| **rsID** | rs113993960 |
| **HGVS** | NM_000492.4:c.1521_1523delCTT |
| **Consequence** | frameshift_variant |
| **Impact** | HIGH |
| **SIFT** | deleterious |
| **PolyPhen** | N/A |
| **ClinVar** | Pathogenic |
| **ClinVar Condition** | Cystic fibrosis |
| **gnomAD AF (global)** | 0.021000 |
| **gnomAD AF (AFR)** | 0.000300 |
| **gnomAD AF (EUR)** | 0.033000 |

### 3. BRCA1 — 17:43063931 G>A

| Field | Value |
|-------|-------|
| **rsID** | rs55770810 |
| **HGVS** | NM_007294.4:c.5095C>T (p.Arg1699Trp) |
| **Consequence** | missense_variant |
| **Impact** | MODERATE |
| **SIFT** | deleterious |
| **PolyPhen** | probably_damaging |
| **ClinVar** | Pathogenic |
| **ClinVar Condition** | Hereditary breast and ovarian cancer syndrome |
| **gnomAD AF (global)** | 0.000024 |
| **gnomAD AF (AFR)** | N/A |
| **gnomAD AF (EUR)** | N/A |

### 4. APOE — 19:44908684 T>C

| Field | Value |
|-------|-------|
| **rsID** | rs429358 |
| **HGVS** | NM_000041.4:c.388T>C |
| **Consequence** | missense_variant |
| **Impact** | MODERATE |
| **SIFT** | tolerated |
| **PolyPhen** | benign |
| **ClinVar** | Risk factor |
| **ClinVar Condition** | Alzheimer disease |
| **gnomAD AF (global)** | 0.147000 |
| **gnomAD AF (AFR)** | 0.232000 |
| **gnomAD AF (EUR)** | 0.143000 |

### 5. MTHFR — 1:11796321 G>A

| Field | Value |
|-------|-------|
| **rsID** | rs1801133 |
| **HGVS** | NM_005957.5:c.665C>T |
| **Consequence** | missense_variant |
| **Impact** | MODERATE |
| **SIFT** | tolerated |
| **PolyPhen** | benign |
| **ClinVar** | Benign/Likely benign |
| **ClinVar Condition** | Homocystinuria |
| **gnomAD AF (global)** | 0.312000 |
| **gnomAD AF (AFR)** | 0.198000 |
| **gnomAD AF (EUR)** | 0.367000 |

---

---
*ClawBio VCF Annotator is a research tool. Not a clinical diagnostic device. Always consult a qualified geneticist for clinical decisions.*
