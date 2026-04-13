---
name: claw-methylation-cycle
version: 0.1.0
author: Samuel Carmona Aguirre <samuel@unimed-consulting.es>
domain: methylation-genomics
description: >
  Methylation cycle analysis with BH4/neurotransmitter axis interpretation.
  Profiles MTHFR, MTRR, MTR, CBS, BHMT, SHMT1, COMT, and AHCY from raw
  genotype data and estimates Net Methylation Capacity and BH4 production
  capacity (dopamine/serotonin synthesis impact). Designed for clinical
  support in neurodevelopmental disorders (ADHD, depression, anxiety).
license: MIT
inputs:
  - name: genotype_file
    type: file
    format: [txt, csv]
    description: >
      23andMe, AncestryDNA, or ADNTRO raw genotype file.
      Tab-separated with columns: rsid, chromosome, position, genotype.
    required: true
outputs:
  - name: report
    type: file
    format: md
    description: >
      Clinical methylation cycle report with enzymatic activity estimates,
      BH4 axis capacity, compound heterozygosity flags, and
      prioritised supplementation recommendations.
  - name: result
    type: file
    format: json
    description: >
      Structured JSON output for downstream integration with CAPS Digital /
      13MIL or any clinical decision-support system.
    dependencies:
  python: ">=3.11"
  packages:
    - pandas>=2.0
    - numpy>=1.24
tags:
  - methylation
  - MTHFR
  - BH4
  - neurotransmitters
  - ADHD
  - pharmacogenomics
  - nutrigenomics
  - neurodevelopmental
  - Holomedicina
  - CAPS-Digital
demo_data:
  - path: demo_input.txt
    description: >
      Synthetic genotype file with 8 methylation-cycle SNPs
      representing a MTHFR compound heterozygous profile (C677T + A1298C)
      with BHMT and COMT variants. NOT real patient data.
endpoints:
  cli: python skills/claw-methylation-cycle/methylation_cycle.py --input {genotype_file} --output {output_dir}
guideline_authority: PMID
guideline_doi: "10.3390/nu13030768"
validation_tier: community
---
## Domain Decisions

### Genes and Variants Assessed

| Gene   | rsID        | Variant        | Allele Assessed | Effect Direction        |
|--------|-------------|----------------|-----------------|-------------------------|
| MTHFR  | rs1801133   | C677T          | T (risk)        | Decreased MTHFR activity |
| MTHFR  | rs1801131   | A1298C         | C (risk)        | Decreased MTHFR activity |
| MTRR   | rs1801394   | A66G           | G (risk)        | Decreased MTRR activity  |
| MTR    | rs1805087   | A2756G         | G (risk)        | Decreased MTR activity   |
| CBS    | rs234706    | C699T          | T (risk)        | Increased CBS activity   |
| BHMT   | rs3733890   | R239Q          | A (risk)        | Decreased BHMT activity  |
| SHMT1  | rs1979277   | C1420T         | T (risk)        | Decreased SHMT1 activity |
| COMT   | rs4680      | Val158Met      | A/Met (risk)    | Decreased COMT activity  |
| AHCY   | rs819147    | AHCY           | T (risk)        | Decreased AHCY activity  |
### Enzymatic Activity Estimates

Activity is estimated as a percentage of normal function based on homozygous vs.
heterozygous status of risk alleles. These are approximations derived from published
functional studies — they are NOT direct enzyme assays.

| Genotype              | Estimated Activity |
|-----------------------|--------------------|
| 0 risk alleles (WT)   | 100%               |
| 1 risk allele (het)   | 60–80% (gene-specific, see table below) |
| 2 risk alleles (hom)  | 15–40% (gene-specific, see table below) |

Gene-specific activity estimates (homozygous risk):

- MTHFR C677T homozygous: ~30% of normal
- MTHFR A1298C homozygous: ~60% of normal
- MTHFR compound heterozygous (C677T + A1298C): ~15% of normal (combined reduction)
- MTRR A66G homozygous: ~60% of normal
- BHMT R239Q homozygous: ~40% of normal
- COMT Val158Met homozygous (Met/Met): ~25% of normal (slow COMT)
- SHMT1 C1420T homozygous: ~60% of normal

Source: Nazki FH et al. (2014) Gene 533(1):11-20; Ledford AW et al. (2021) Nutrients 13(3):768.
### Net Methylation Capacity (NMC) Score

NMC is a composite index (0–100) derived from weighted enzymatic activities
of the eight assessed genes. Weights reflect clinical significance for
methylation flux:

- MTHFR: weight 0.35 (primary rate-limiting enzyme)
- MTRR: weight 0.15
- BHMT: weight 0.15
- COMT: weight 0.10
- MTR: weight 0.10
- CBS: weight 0.05 (inverse — CBS upregulation diverts homocysteine)
- SHMT1: weight 0.05
- AHCY: weight 0.05

NMC < 40: Severely reduced — clinical intervention indicated
NMC 40–60: Moderately reduced — supplementation strongly recommended
NMC 60–80: Mildly reduced — dietary optimisation recommended
NMC > 80: Normal range

### BH4 Axis Capacity

BH4 (tetrahydrobiopterin) is an essential cofactor for tyrosine hydroxylase
(dopamine synthesis) and tryptophan hydroxylase (serotonin synthesis).
MTHFR activity directly constrains BH4 regeneration via the folate cycle.

BH4 capacity is estimated from MTHFR combined activity and MTRR modifier:
- Base BH4 capacity = MTHFR_activity × MTRR_modifier
- MTRR_modifier: homozygous risk = 0.75; heterozygous = 0.88; WT = 1.0

BH4 Axis thresholds:
- < 40%: Severely reduced — dopamine and serotonin synthesis substantially impaired
- 40–65%: Moderately reduced — neurotransmitter synthesis may be clinically relevant
- > 65%: Within normal range
### Compound Heterozygosity

MTHFR compound heterozygous (C677T + A1298C simultaneously) is the most
clinically significant single-gene methylation finding. When both variants
are present, total MTHFR enzymatic activity is reduced more than either
variant alone. This combination is flagged explicitly in the output.

### Reference Framework

Clinical interpretation is grounded in the Holomedicina® framework
(Samuel Carmona Aguirre, 2014/UNESCO 2016) and the MH-AIAP methodology,
which integrates genomic findings with the subject's biographical history
(FHH — Formulario de Historicidad Holónica) and semantic clinical map
for holonic interpretation. Raw genetic output alone does not constitute
a clinical recommendation.

## Safety Rules

1. Never report a clinical diagnosis. Always include the RUO (Research Use Only) disclaimer.
2. Never recommend specific drug dosages or prescribe medication changes.
3. Always flag MTHFR compound heterozygous status as requiring clinical review.
4. Flag BH4 capacity < 40% with an explicit neurodevelopmental implications warning.
5. Never extrapolate findings to ancestries not represented in the source studies
   (findings validated primarily in European-ancestry populations).
6. Unknown SNPs or SNPs not present in the input file must be reported as
   "Not assessed" — never assume wildtype.
7. All supplementation suggestions are Priority-ranked guidance for a clinician —
   not direct patient instructions.

## Agent Boundary

### In Scope
- Genotype extraction for 9 methylation-cycle genes from raw DTC genotype files
- Enzymatic activity estimation (percentage of normal function)
- Net Methylation Capacity (NMC) composite index calculation
- BH4 axis capacity estimation and neurotransmitter synthesis impact
- MTHFR compound heterozygosity detection and flagging
- Prioritised clinical recommendations for clinician review
- JSON output for integration with downstream clinical decision-support systems

### Out of Scope
- Dosing recommendations (requires clinical context and prescribing authority)
- Diagnosis of methylation disorders or neurodevelopmental conditions
- Drug-drug interaction analysis (see PharmGx Reporter)
- Epigenetic state (methylation is a phenotype; this skill assesses genotype only)
- Whole-genome or whole-exome sequencing data (SNP array input only)
- Direct patient communication (output is for qualified clinician use)
