---
name: claw-methylation-cycle
license: MIT
metadata:
  openclaw:
    requires:
      bins:
      - python3
      env: null
      config: null
    always: false
    emoji: 🧬
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install: 'pip install -e .

      '
    trigger_keywords:
    - methylation
    - MTHFR
    - folate cycle
    - BH4
    - homocysteine
    - methylation cycle
    - metilación
    - neurotransmitter synthesis
  author: Samuel Carmona Aguirre <samuel@unimed-consulting.es>
  demo_data_path: demo_input.txt
  dependencies_python: '>=3.9'
  guideline_doi: 10.3390/nu13030768
  input_format: 23andme, adntro, ancestry
  output_format: markdown, json
  tags:
  - genomics
  - methylation
  - MTHFR
  - neurodevelopment
  trigger_keywords:
  - methylation
  - MTHFR
  - folate cycle
  - BH4
  - homocysteine
  - methylation cycle
  - metilación
  - neurotransmitter synthesis
  validation_tier: community
  version: 0.1.3
description: Methylation cycle analysis — enzymatic activity profiles, Net Methylation Capacity, BH4 axis estimates, compound heterozygosity detection from SNP genotype data.
---

# claw-methylation-cycle

Methylation cycle analysis skill for ClawBio. Produces enzymatic activity
profiles, Net Methylation Capacity (NMC), BH4 axis estimates, compound
heterozygosity detection, and clinician-review genotype findings from raw
SNP genotype data.

---

## Trigger

**Fire this skill when:**

- The user asks about methylation, MTHFR variants, folate cycle, or
  homocysteine risk from a genotype file.
- A raw SNP file (23andMe / ADNTRO / Ancestry format) is provided and the
  clinical question involves methylation, BH4, dopamine/serotonin synthesis
  capacity, or neurodevelopmental contexts (ADHD, depression, anxiety).
- The upstream workflow (PharmGx Reporter, NutriGx Advisor) has flagged
  MTHFR or MTRR and the clinician needs the full methylation panel.
- Keywords present: `methylation`, `MTHFR`, `BH4`, `folate cycle`,
  `metilación`, `ciclo de metilación`, `homocysteine`, `5-MTHF`,
  `methylcobalamin`, `neurotransmitter synthesis`, `dopamine upstream`.

**Do NOT fire this skill when:**

- The question is purely about folic acid supplementation without a genotype file.
- The user is asking about MTHFR in the context of thrombophilia/clotting
  only — use PharmGx Reporter for warfarin/anticoagulation questions.
- Only N-GENE polygenic risk data is available (no raw SNP file) — this skill
  requires genotype-level input; PRS percentiles are not sufficient.
- The SNP file format is VCF, FASTQ, BAM, or PLINK binary — these require
  preprocessing before this skill can run.
- The clinical question is exclusively pharmacogenomic (CYP enzymes) — use
  PharmGx Reporter instead.

---

## Workflow

1. **Receive input** — Accept either a raw genotype file path or a pre-parsed
   `snp_dict`. If a file is provided, call `parse_genotype_file()` to extract
   the rsID → genotype mapping.

2. **Panel coverage check** — Compare detected rsIDs against the 9-gene
   methylation panel. Log missing SNPs. For any SNP absent from the input,
   mark the corresponding gene as `not_assessed` — do NOT silently assume
   normal activity (Safety Rule 6).

3. **Enzymatic activity scoring** — For each gene, map the diplotype to an
   estimated activity percentage. Heterozygous risk variants reduce activity
   by their assigned weight; homozygous variants apply the full reduction.

4. **Compound heterozygosity detection** — Check MTHFR C677T (rs1801133) and
   A1298C (rs1801131) simultaneously. If both are heterozygous, set
   `compound_heterozygosity = True` and apply the combined activity reduction
   (~15% of normal — more severe than either variant alone).

5. **Net Methylation Capacity (NMC)** — Compute the weighted average of all
   enzyme activities. Clamp to [0, 100]. Expose `coverage_pct` and
   `snps_missing`; flag NMC as partial if key SNPs are absent.

6. **BH4 axis capacity** — Derive BH4 from MTHFR activity and MTRR modifier.
   Report clinical implications for dopamine and serotonin synthesis in
   neurodevelopmental presentations.

7. **Prioritised recommendations** — Generate PRIORITY 1 / 2 / 3
   recommendations based on active findings. Lead with highest clinical
   impact (compound het MTHFR or severely reduced BH4).

8. **Output** — Write `report.md` (human-readable) and `result.json`
   (structured, for downstream integration).

---

## Example Output

```
╔══════════════════════════════════════════════════════════════╗
║  ClawBio · Methylation Cycle Analysis Report                 ║
║  Author: Samuel Carmona Aguirre · RUO — Not a medical device ║
╚══════════════════════════════════════════════════════════════╝

Executive Summary
─────────────────
Net Methylation Capacity : 53 / 100  🔴 Reduced
BH4 Axis Capacity        : 31 / 100  🔴 Reduced
MTHFR Compound Het.      : YES ⚠️   (C677T + A1298C)
Dopamine Synthesis       : Severely Reduced
Serotonin Synthesis      : Severely Reduced

Enzymatic Activity Profile
──────────────────────────
Gene    Activity  Status                Key Variants
MTHFR    15%     🔴 Severely reduced   C677T, A1298C
MTRR     60%     🟡 Moderately reduced A66G
MTR     100%     🟢 Normal             –
CBS     100%     🟢 Normal             –
BHMT     40%     🔴 Moderately reduced R239Q
SHMT1    80%     🟢 Mildly reduced     C1420T
COMT     55%     🟡 Moderately reduced Val158Met
AHCY    100%     🟢 Normal             –

Clinical Recommendations — FOR CLINICIAN REVIEW ONLY
──────────────────────────────────────────────────────
⚠️  The following is genotype-based information for qualified clinician use.
    Do not self-administer. All nutrients listed are reported in the peer-reviewed
    literature for the pathways indicated; dosing and indication require
    individualised clinical assessment.

Genotype findings:
  • 5-MTHF (methylfolate) preferred over synthetic folic acid (MTHFR C677T/A1298C).
    Ref: Lamers Y et al. (2004) Am J Clin Nutr 80(5):1234-41.
  • MTHFR compound het: methylcobalamin co-administration reported in literature.
    Ref: Ledford AW et al. (2021) Nutrients 13(3):768.
  • BH4 capacity at 31%: riboflavin (B2) reported as MTHFR cofactor supporting BH4
    regeneration. Ref: McNulty H et al. (2017) Am J Clin Nutr 106(1):128-36.
  • MTRR A66G: methylcobalamin preferred over cyanocobalamin per functional studies.
    Ref: Olteanu H et al. (2002) Biochemistry 41(45):13378-85.
  • BHMT R239Q: betaine and choline-rich foods reported as alternative methyl donors.
    Ref: Slow S et al. (2004) Clin Chim Acta 340(1-2):57-67.

Some literature reports an association between BH4 deficiency and
ADHD/depression/anxiety phenotypes. This genotype indicates reduced BH4
production capacity. A clinician should contextualise this finding with
the patient's clinical history.
Where clinically relevant, a clinician may evaluate whether neurodevelopmental
symptoms correlate with BH4 capacity for potential non-pharmacological support.
```

---

## Gotchas

1. **Missing SNPs must never be silently normalised.** The current
   implementation (line 471) defaults to assuming normal activity for SNPs
   absent from the input. This produces an artificially high NMC. When key
   SNPs are missing, always expose `coverage_pct` and `snps_missing` so
   downstream consumers know the score is partial.

2. **Compound heterozygosity is synergistic, not additive.** C677T and A1298C
   affect different MTHFR domains. Their combined effect (~15% activity) is
   greater than either variant alone. Do not compute as
   `activity(677) × activity(1298)`.

3. **BH4 capacity is an estimate, not a measured value.** The BH4 score is
   derived from MTHFR activity and literature-based weights. It does not
   account for DHFR variation or dietary cofactor availability. Always
   include the RUO disclaimer.

4. **COMT Val158Met has a dual role.** rs4680 appears in both methylation
   (SAM consumption) and dopamine/catecholamine panels. Always note this —
   do not report it in isolation.

5. **DTC array coverage varies by platform.** ADNTRO covers all 9 panel SNPs
   for most European-ancestry samples. 23andMe v3 and Ancestry v1 may not
   include rs1801394 (MTRR) or rs3733890 (BHMT). Always check `snps_missing`.

6. **This skill does not cover pharmacogenomics.** SLCO1B1, CYP enzymes, and
   statin/warfarin risk belong to PharmGx Reporter, not this skill.

---

## Domain Decisions

### Genes and Variants Assessed

| Gene  | rsID      | Variant   | Allele Assessed | Effect Direction         |
|-------|-----------|-----------|-----------------|--------------------------|
| MTHFR | rs1801133 | C677T     | T (risk)        | Decreased MTHFR activity |
| MTHFR | rs1801131 | A1298C    | C (risk)        | Decreased MTHFR activity |
| MTRR  | rs1801394 | A66G      | G (risk)        | Decreased MTRR activity  |
| MTR   | rs1805087 | A2756G    | G (risk)        | Decreased MTR activity   |
| CBS   | rs234706  | C699T     | T (risk)        | Increased CBS activity   |
| BHMT  | rs3733890 | R239Q     | A (risk)        | Decreased BHMT activity  |
| SHMT1 | rs1979277 | C1420T    | T (risk)        | Decreased SHMT1 activity |
| COMT  | rs4680    | Val158Met | A/Met (risk)    | Decreased COMT activity  |
| AHCY  | rs819147  | AHCY      | T (risk)        | Decreased AHCY activity  |

### Enzymatic Activity Estimates

Activity is estimated as a percentage of normal function based on homozygous vs.
heterozygous status of risk alleles. These are approximations derived from published
functional studies — they are NOT direct enzyme assays.

| Genotype             | Estimated Activity                      |
|----------------------|-----------------------------------------|
| 0 risk alleles (WT)  | 100%                                    |
| 1 risk allele (het)  | 60–80% (gene-specific, see below)       |
| 2 risk alleles (hom) | 15–40% (gene-specific, see below)       |

Gene-specific estimates (homozygous risk):

- MTHFR C677T homozygous: ~30% of normal
- MTHFR A1298C homozygous: ~60% of normal
- MTHFR compound heterozygous (C677T + A1298C): ~15% of normal
- MTRR A66G homozygous: ~60% of normal
- BHMT R239Q homozygous: ~40% of normal
- COMT Val158Met homozygous (Met/Met): ~25% of normal
- SHMT1 C1420T homozygous: ~60% of normal

Source: Nazki FH et al. (2014) Gene 533(1):11-20; Ledford AW et al. (2021) Nutrients 13(3):768.

### Net Methylation Capacity (NMC) Score

NMC is a composite index (0–100) derived from weighted enzymatic activities:

- MTHFR: weight 0.35 (primary rate-limiting enzyme)
- MTRR: weight 0.15
- BHMT: weight 0.15
- COMT: weight 0.10
- MTR: weight 0.10
- CBS: weight 0.05 (inverse — upregulation diverts homocysteine)
- SHMT1: weight 0.05
- AHCY: weight 0.05

NMC < 40: Severely reduced
NMC 40–60: Moderately reduced
NMC 60–80: Mildly reduced
NMC > 80: Within normal range

> **Note (ACMG 2013):** Routine population screening for MTHFR variants is not
> recommended for thrombosis risk assessment. This tool reports genotype facts
> for clinician contextualisation; NMC bands are descriptive outputs, not
> intervention triggers. Clinical decisions require individual patient evaluation.

### BH4 Axis Capacity

BH4 (tetrahydrobiopterin) is an essential cofactor for tyrosine hydroxylase
(dopamine) and tryptophan hydroxylase (serotonin). MTHFR activity directly
constrains BH4 regeneration via the folate cycle.

- Base BH4 capacity = MTHFR_activity × MTRR_modifier
- MTRR_modifier: homozygous risk = 0.75; heterozygous = 0.88; WT = 1.0

BH4 thresholds:
- < 40%: Severely reduced — neurotransmitter synthesis substantially constrained
- 40–65%: Moderately reduced — may be clinically relevant; clinician to contextualise
- > 65%: Within normal range

### Compound Heterozygosity

MTHFR compound heterozygous (C677T + A1298C simultaneously) is the most
clinically significant single-gene methylation finding. Total MTHFR activity
is reduced more than either variant alone. Flagged explicitly in output.

### Author and Attribution

Developed by Samuel Carmona Aguirre (samuel@unimed-consulting.es) as part of a
contribution to the ClawBio open-source bioinformatics library.

**Conflict of Interest (COI):** The author develops clinical genomics workflows
in a private practice context and may use tools derived from this skill as a
component in those workflows. This skill is contributed as a standalone
open-source genotype reporting tool; its output is not specific to any
proprietary clinical platform. All clinical integration decisions rest with
the qualified end-user clinician.

**ACMG caveat on MTHFR testing (Green et al., Genet Med 2013, 15:153–156):**
Routine population screening for MTHFR variants is not currently recommended by
ACMG for assessment of thrombosis risk, neural tube defect risk, or psychiatric
phenotypes. Genotype findings from this tool should be interpreted in the context
of the individual patient's clinical history by a qualified clinician.

---

## Safety Rules

1. Never report a clinical diagnosis. Always include the RUO disclaimer.
2. Never recommend specific drug dosages or prescribe medication changes.
3. Always flag MTHFR compound heterozygous status as requiring clinical review.
4. Flag BH4 capacity < 40% with an explicit neurodevelopmental warning.
5. Never extrapolate to ancestries not represented in source studies
   (validated primarily in European-ancestry populations).
6. Unknown SNPs or SNPs not in the input must be reported as "Not assessed" —
   never assume wildtype.
7. All supplementation suggestions are Priority-ranked guidance for a clinician,
   not direct patient instructions.

---

## Agent Boundary

### In Scope

- Genotype extraction for 9 methylation-cycle genes from raw DTC files
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

---

## Usage

```bash
python skills/claw-methylation-cycle/methylation_cycle.py \
  --input path/to/genotype.txt \
  --output results/
```

---

## Output Structure

```
<output_dir>/
├── report.md
├── result.json
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

---

## References

- Nazki FH et al. (2014). Folate: metabolism, genes, polymorphisms. *Gene*, 533(1), 11-20.
- Ledford AW et al. (2021). MTHFR and BH4 pathway in neuropsychiatric disorders. *Nutrients*, 13(3):768.
- Stover PJ (2009). One-carbon metabolism-genome interactions. *J Nutr*, 139(12), 2402-5.
- Lamers Y et al. (2004). Supplementation with [6S]-5-methyltetrahydrofolate or folic acid equally reduces plasma total homocysteine. *Am J Clin Nutr*, 80(5):1234-41.
- McNulty H et al. (2017). Riboflavin lowers homocysteine in children and adults with common MTHFR polymorphism. *Am J Clin Nutr*, 106(1):128-36.
- Olteanu H et al. (2002). Differences in the efficiency of reductive methylation of cob(II)alamin. *Biochemistry*, 41(45):13378-85.
- Slow S et al. (2004). Plasma betaine and homocysteine. *Clin Chim Acta*, 340(1-2):57-67.
- Esteller M. (2014). Introduccion a la epigenetica. *SEBBM*, 179:4-6.
- Spuch C, Agis-Balboa RC. (2014). Epigenetica en neurociencias. *SEBBM*, 179:18-21.
- Green RC et al. (2013). ACMG recommendations for reporting of incidental findings in clinical exome and genome sequencing. *Genet Med*, 15(7):565-74.
- ClawBio (2026). https://github.com/ClawBio/ClawBio

---

## Changelog

| Version | Date       | Change |
|---------|------------|--------|
| 0.1.0   | 2026-04-07 | Initial release. Validated on ASES-2307-002. |
| 0.1.1   | 2026-04-14 | Fixed SKILL.md per PR #133: single YAML block, added Trigger, Workflow, Example Output, Gotchas. Removed unused pandas. Documented line-471 design decision. |
| 0.1.2   | 2026-05-12 | Framing revisions per Manuel Corpas review: reworded BH4 causal language to association-based; reformatted supplement list as clinician-review block with per-nutrient citations; added 6 DOIs; converted em-dashes to ASCII in test docstrings; removed duplicate top-level test file. |
| 0.1.3   | 2026-05-15 | Per PR #133 third review: (1) removed all specific dosages from recommendations; (2) stripped Holomedicina/CAPS/UNIMED/MH-AIAP branding from YAML, description, trigger, example output, Domain Decisions, report header — kept as author attribution only; (3) reframed NMC and BH4 bands as descriptive outputs, not intervention triggers; (4) added ACMG 2013 caveat on routine MTHFR testing; (5) added author COI disclosure; (6) removed obsolete Gotcha #5 (pandas); (7) softened compound het "strongly indicated" language. |
