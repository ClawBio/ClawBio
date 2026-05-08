---
name: methylation-cycle
version: 0.1.1
author: Samuel Carmona Aguirre
license: MIT
trigger_keywords:
  - methylation
  - MTHFR
  - folate cycle
  - BH4
  - homocysteine
  - methylation cycle
  - metilacion
  - neurotransmitter synthesis
input_format: 23andme, adntro, ancestry
output_format: markdown, json
dependencies_python: ">=3.9"
framework: Holomedicina CAPS Digital UNIMED Consulting
tags:
  - genomics
  - methylation
  - MTHFR
  - neurodevelopment
  - CAPS-Digital
demo_data_path: demo_input.txt
guideline_doi: "10.3390/nu13030768"
validation_tier: community
---

# claw-methylation-cycle

Methylation cycle analysis skill for ClawBio. Produces enzymatic activity
profiles, Net Methylation Capacity (NMC), BH4 axis estimates, compound
heterozygosity detection, and clinical recommendations from raw SNP genotype
data. Integrated into the CAPS Digital / Holomedicina clinical framework
(Samuel Carmona Aguirre, 2014/UNESCO 2016).

---

## Trigger

**Fire this skill when:**

- The user asks about methylation, MTHFR variants, folate cycle, or
  homocysteine risk from a genotype file.
- A raw SNP file (23andMe, ADNTRO, or Ancestry format) is provided and the
  clinical question involves methylation, BH4, dopamine or serotonin synthesis
  capacity, or neurodevelopmental contexts such as ADHD, depression, anxiety.
- The upstream workflow (PharmGx Reporter, NutriGx Advisor) has flagged
  MTHFR or MTRR and the clinician needs the full methylation panel.
- Keywords present: methylation, MTHFR, BH4, folate cycle, metilacion,
  homocysteine, 5-MTHF, methylcobalamin, neurotransmitter synthesis.

**Do NOT fire this skill when:**

- The question is purely about folic acid supplementation without a genotype file.
- The user is asking about MTHFR in the context of thrombophilia only.
  Use PharmGx Reporter for warfarin and anticoagulation questions.
- Only N-GENE polygenic risk data is available with no raw SNP file.
  This skill requires genotype-level input. PRS percentiles are not sufficient.
- The SNP file format is VCF, FASTQ, BAM, or PLINK binary.
  These require preprocessing before this skill can run.
- The clinical question is exclusively pharmacogenomic (CYP enzymes).
  Use PharmGx Reporter instead.

---

## Workflow

1. Receive input. Accept a raw genotype file path or a pre-parsed snp_dict.
   Call parse_genotype_file() to extract the rsID to genotype mapping.

2. Panel coverage check. Compare detected rsIDs against the 9-gene methylation
   panel. Mark missing SNPs as not_assessed. Do NOT silently assume normal
   activity. This is Safety Rule 6.

3. Enzymatic activity scoring. For each gene, map the diplotype to an estimated
   activity percentage. Heterozygous risk variants reduce activity by their
   assigned weight. Homozygous variants apply the full reduction.

4. Compound heterozygosity detection. Check MTHFR C677T (rs1801133) and A1298C
   (rs1801131) simultaneously. If both are heterozygous, set
   compound_heterozygosity to True and apply the combined reduction of
   approximately 15% of normal, more severe than either variant alone.

5. Net Methylation Capacity (NMC). Compute the weighted average of all enzyme
   activities. Clamp to 0-100. Expose coverage_pct and snps_missing. Flag NMC
   as partial if key SNPs are missing.

6. BH4 axis capacity. Derive BH4 from MTHFR activity and MTRR modifier. Report
   clinical implications for dopamine and serotonin synthesis.

7. Prioritised recommendations. Generate PRIORITY 1, 2, and 3 recommendations.
   Lead with the highest clinical impact finding.

8. Output. Write report.md and result.json for 13MIL v6.0 integration.

---

## Example Output

```
ClawBio Methylation Cycle Clinical Report
Framework: Holomedicina / CAPS Digital / UNIMED Consulting

Net Methylation Capacity : 53 / 100  REDUCED
BH4 Axis Capacity        : 31 / 100  REDUCED
MTHFR Compound Het.      : YES (C677T + A1298C)
Dopamine Synthesis       : Severely Reduced
Serotonin Synthesis      : Severely Reduced

Gene    Activity  Status               Key Variants
MTHFR    15%     Severely reduced     C677T, A1298C
MTRR     60%     Moderately reduced   A66G
MTR     100%     Normal               -
CBS     100%     Normal               -
BHMT     40%     Moderately reduced   R239Q
SHMT1    80%     Mildly reduced       C1420T
COMT     55%     Moderately reduced   Val158Met
AHCY    100%     Normal               -

PRIORITY 1 - Use 5-MTHF not synthetic folic acid.
PRIORITY 1 - MTHFR compound het: 5-MTHF plus methylcobalamin strongly indicated.
PRIORITY 2 - BH4 at 31%: Riboflavin B2 200-400 mg/day plus Vitamin C 500 mg/day.
PRIORITY 2 - Evaluate ADHD/depression/anxiety re BH4 deficit before pharma.
PRIORITY 2 - MTRR variant: prefer methylcobalamin over cyanocobalamin.
PRIORITY 3 - BHMT R239Q: betaine TMG 500-1000 mg/day plus choline-rich foods.
```

---

## Gotchas

1. Missing SNPs must never be silently normalised. Line 471 assumes normal
   activity for absent SNPs producing artificially high NMC. Always expose
   coverage_pct and snps_missing so downstream consumers know the score is partial.

2. Compound heterozygosity is synergistic not additive. C677T and A1298C affect
   different MTHFR domains. Combined effect is approximately 15% of normal,
   greater than either variant alone. Do not compute as activity(677) multiplied
   by activity(1298).

3. BH4 capacity is an estimate not a measured value. Derived from MTHFR activity
   and literature-based weights. Does not account for DHFR variation or dietary
   cofactors. Always include the RUO disclaimer.

4. COMT Val158Met has a dual role in methylation via SAM consumption and in
   dopamine metabolism. Always note this dual role and do not report in isolation.

5. The pandas import on line 235 is unused. Remove it. Adds approximately 40 MB
   to the dependency footprint with no current function.

6. DTC array coverage varies by platform. ADNTRO covers all 9 panel SNPs for
   most European-ancestry samples. 23andMe v3 and Ancestry v1 may not include
   rs1801394 (MTRR) or rs3733890 (BHMT). Always check snps_missing.

7. This skill does not cover pharmacogenomics. SLCO1B1, CYP enzymes, statin
   risk, and warfarin risk belong to PharmGx Reporter, not this skill.

---

## Safety Rules

1. Never report a clinical diagnosis. Always include the RUO disclaimer.
2. Never recommend specific drug dosages or prescribe medication changes.
3. Always flag MTHFR compound heterozygous status as requiring clinical review.
4. Flag BH4 capacity below 40% with an explicit neurodevelopmental warning.
5. Never extrapolate to ancestries not represented in source studies.
6. Unknown SNPs must be reported as Not assessed. Never assume wildtype.
7. Supplementation suggestions are Priority-ranked guidance for a clinician only.

---

## Agent Boundary

In Scope: genotype extraction for 9 methylation-cycle genes, enzymatic activity
estimation, NMC calculation, BH4 axis capacity, compound heterozygosity detection,
prioritised recommendations, JSON output for clinical decision-support.

Out of Scope: dosing recommendations, diagnosis, drug-drug interactions, epigenetic
state, whole-genome sequencing data, direct patient communication.

---

## References

- Nazki FH et al. (2014). Folate metabolism genes polymorphisms. Gene, 533(1), 11-20.
- Ledford AW et al. (2021). MTHFR and BH4 pathway in neuropsychiatric disorders. Nutrients, 13(3):768.
- Stover PJ (2009). One-carbon metabolism genome interactions. J Nutr, 139(12), 2402-5.
- Carmona Aguirre S. (2014/UNESCO 2016). Holomedicina. UNIMED Consulting.
- ClawBio (2026). https://github.com/ClawBio/ClawBio

---

## Changelog

| Version | Date       | Change |
|---------|------------|--------|
| 0.1.0   | 2026-04-07 | Initial release. Validated on ASES-2307-002. |
| 0.1.1   | 2026-04-14 | Fixed SKILL.md per PR 133. Added Trigger, Workflow, Example Output, Gotchas. Single YAML block. Removed unused pandas. Documented line 471 design decision. |
