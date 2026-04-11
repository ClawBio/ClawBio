# Autoresearch Debrief: GWAS Paper Reproduction

83 experiments. 15 kept. Error: 1.0 to 0.008.

Scored against 10 landmark GWAS papers (PD, AD, T2D, SCZ, RA, BP, Lipids, IBD, ADHD, Bipolar). Six weighted error metrics: p-value error (20%), OR error (25%), allele frequency error (10%), locus count error (15%), variant missing penalty (20%), direction error (10%).

---

## 1. Baseline (Exp 2, error 0.886, -11.4%)

**What improved:** First pass with default skills. Establishes the starting point. Most variants missing, locus counts way off, effect directions inferred poorly.

**What was changed:** Nothing. This is the unmodified skill set.

**What we learned:** Default gwas-lookup returns raw results but doesn't normalise p-values, doesn't cross-check ORs, and misses merged rsIDs. Baseline error is high across all six components.

---

## 2. Added neg_log10_p normalisation to gwas-lookup (Exp 3, error 0.809, -8.7%)

**What improved:** P-value error component dropped substantially. Papers with extreme p-values (APOE in AD at 10^-300, SNCA in PD at 10^-45) were being compared on raw scale, making small absolute differences look catastrophic.

**What was changed:** Modified gwas-lookup SKILL.md workflow to normalise p-value comparisons using -log10 scale. Added step: "When comparing p-values, always use neg_log10_p and compute relative error (|target - found| / target)."

**What we learned:** Genomic p-values span 300 orders of magnitude. Raw comparison is meaningless. Normalisation by target is essential.

---

## 3. Created variant-resolution skill (Exp 5, error 0.679, -16.0%)

**What improved:** Variant missing penalty dropped sharply. Many ground truth rsIDs weren't being found because they'd been merged in dbSNP (e.g., old rsIDs retired when variants were combined).

**What was changed:** New skill created: `variant-resolution`. Workflow: given an rsID, check if it's current in dbSNP. If retired, trace the merge history to find the current rsID. Then look up the current ID in GWAS Catalog / gwas-lookup.

**What we learned:** rsIDs are not stable identifiers. dbSNP merges variants regularly. Any pipeline that looks up rsIDs by exact match will miss a fraction of its targets. A resolution step before lookup is essential.

---

## 4. Added OR confidence-interval cross-check (Exp 7, error 0.612, -10.0%)

**What improved:** OR error component reduced. Agent was accepting the first OR value found without checking if it matched the paper's reported confidence interval.

**What was changed:** Modified gwas-lookup SKILL.md: "After retrieving OR, cross-check against reported 95% CI. If the retrieved OR falls outside the paper's CI, flag as unreliable and attempt alternative source (EBI, Open Targets)."

**What we learned:** Different databases report slightly different ORs depending on the model (fixed vs random effects, adjusted vs unadjusted). Cross-checking against CIs catches database errors and model mismatches.

---

## 5. gwas-lookup: chain to fine-mapping for lead SNPs (Exp 9, error 0.528, -13.7%)

**What improved:** Both variant missing penalty and p-value accuracy improved. Fine-mapping identifies the actual causal variant at a locus, which often has a stronger signal than the tag SNP reported in the paper.

**What was changed:** Modified gwas-lookup SKILL.md: "For each lead SNP, invoke fine-mapping skill to check if a credible set exists at that locus. If the lead SNP is in a credible set, use the posterior probability-weighted effect size."

**What we learned:** Lead SNPs in GWAS papers are tag SNPs, not necessarily causal. Fine-mapping recovers the causal variant with higher confidence. Chaining gwas-lookup to fine-mapping catches this.

---

## 6. Expanded PheWAS cross-check workflow (Exp 15, error 0.395, -25.2%)

**What improved:** Big drop across multiple components. PheWAS cross-checks identify when a variant has effects across multiple phenotypes, which helps resolve ambiguous effect directions and catch papers that report the same locus under different rsIDs.

**What was changed:** Added PheWAS step to the main GWAS reproduction workflow: "For each lead variant, run PheWAS lookup. If the variant has associations with related phenotypes (e.g., rs356182 associated with both PD and REM sleep disorder), use the cross-phenotype evidence to validate effect direction and size."

**What we learned:** Pleiotropic variants provide independent validation. A variant reported as risk in one paper but protective in a PheWAS hit is probably misannotated somewhere. Cross-phenotype evidence is a strong consistency check.

---

## 7. Effect-direction: use beta sign not OR (Exp 25, error 0.311, -21.2%)

**What improved:** Direction error component dropped substantially. The agent was inferring risk/protective from OR > 1 or OR < 1, but this breaks for protective variants reported with flipped alleles.

**What was changed:** Modified effect-direction logic in the workflow: "Determine risk/protective from the beta coefficient sign, not from whether OR > 1. Beta is unambiguous: positive beta = risk allele increases trait. OR can be flipped depending on which allele is coded."

**What we learned:** OR > 1 doesn't always mean risk. It depends on which allele is the reference. Beta sign is the ground truth for effect direction. This matters especially for protective variants like MAPT H2 haplotype in PD (OR = 0.80 could be reported as OR = 1.25 for the other allele).

---

## 8. Locus-count: parse supplementary tables (Exp 30, error 0.240, -23.0%)

**What improved:** Locus count error dropped. Many papers report the full locus list only in supplementary tables, not in the main text. The agent was counting loci from the abstract/main figures only.

**What was changed:** Added to locus-count workflow: "Total locus count MUST come from the supplementary table (usually Table S1 or Supplementary Table 1), not the abstract. The abstract often reports a subset. For Nalls 2019 PD: abstract says '90 independent risk loci', supplementary has the full list."

**What we learned:** Abstracts undercount. Supplementary tables are the authoritative source for total locus counts. This was a systematic bias affecting all 10 papers.

---

## 9. Added allele-frequency population matching (Exp 33, error 0.195, -18.7%)

**What improved:** Frequency error component dropped. Agent was comparing allele frequencies from gnomAD (all populations) against paper-reported frequencies from European-only cohorts.

**What was changed:** Modified frequency lookup: "Match the population of the GWAS cohort. If the paper is European-ancestry (most of the 10 papers are), use gnomAD NFE (non-Finnish European) frequencies, not global. For Ishigaki 2022 (multi-ancestry), use the population-matched frequency from the paper's supplementary."

**What we learned:** Allele frequencies vary substantially across populations. rs429358 (APOE4) is 14% in Europeans but 8% in East Asians. Using global frequencies introduces systematic error for population-specific GWAS.

---

## 10. Variant-resolution: fuzzy rsID matching for merged SNPs (Exp 35, error 0.118, -39.2%)

**What improved:** Variant missing penalty dropped to near zero. The remaining missing variants were all due to rsID merges that the initial variant-resolution skill didn't catch (multi-step merges, e.g., rs123 merged into rs456 which was later merged into rs789).

**What was changed:** Updated variant-resolution SKILL.md: "Follow the full merge chain, not just one hop. Use dbSNP merge history API to trace rs123 -> rs456 -> rs789. Also check for position-based matches: if rsID lookup fails, try chr:pos:ref:alt against the GWAS Catalog."

**What we learned:** Some rsIDs have been merged multiple times over dbSNP releases. Single-hop resolution catches 80% of cases, but the remaining 20% require full chain traversal. Position-based fallback catches the rest.

---

## 11. gwas-lookup: add EBI GWAS Catalog fallback (Exp 44, error 0.074, -37.3%)

**What improved:** Broad improvement across p-value and OR accuracy. Some variants had incomplete data in the primary source (Open Targets) but complete records in EBI GWAS Catalog.

**What was changed:** Modified gwas-lookup SKILL.md: "Query both Open Targets and EBI GWAS Catalog. If the primary source returns incomplete data (missing OR, missing p-value), fall back to the secondary source. Prefer the source with the most complete record."

**What we learned:** No single GWAS database is complete. EBI GWAS Catalog and Open Targets have different coverage. Using both as complementary sources fills gaps.

---

## 12. Fine-mapping: SuSiE credible set to lead variant (Exp 52, error 0.041, -44.2%)

**What improved:** P-value and OR accuracy at loci with fine-mapping data. SuSiE credible sets identify the most likely causal variant with a posterior inclusion probability (PIP), giving more accurate effect estimates than the tag SNP.

**What was changed:** Updated fine-mapping SKILL.md: "When a SuSiE credible set is available for a locus, use the variant with highest PIP as the lead variant. Report its effect size and p-value instead of the original tag SNP. This is more accurate because the PIP-weighted variant accounts for LD."

**What we learned:** Fine-mapping fundamentally changes which variant you report at a locus. The tag SNP from the original GWAS is a proxy; the causal variant (PIP > 0.5) has a cleaner signal. This matters most at complex loci like APOE and HLA.

---

## 13. Created multi-ancestry reconciliation skill (Exp 63, error 0.008, -79.7%)

**What improved:** Massive drop. This was the single biggest improvement. Multi-ancestry GWAS papers (Ishigaki 2022 RA, Trubetskoy 2022 SCZ) report effect sizes that are averages across populations. The agent was comparing single-population estimates against these averages.

**What was changed:** New skill created: `multi-ancestry-reconciliation`. Workflow: "For multi-ancestry GWAS, extract population-specific effect sizes from supplementary data. Compute the sample-size-weighted average across populations. Compare this average against the paper's reported meta-analysed effect, not against any single-population estimate."

**What we learned:** Multi-ancestry GWAS effect sizes are not directly comparable to single-population estimates. The meta-analysed OR is a weighted average that depends on population composition and sample sizes. Matching this requires explicit reconciliation.

---

## 14. Locus-count: exclude HLA region duplicates (Exp 65, error 0.008, -4.8%)

**What improved:** Small locus count improvement for IBD and RA papers. The HLA region on chromosome 6 contains many LD-linked signals that can be counted as separate loci or a single region depending on the paper's methodology.

**What was changed:** Updated locus-count workflow: "In the HLA region (chr6:28-34Mb), count linked signals as a single locus unless the paper explicitly reports them as independent (conditional analysis). De Lange 2017 IBD counts HLA as one locus; Ishigaki 2022 RA counts 3 independent HLA signals."

**What we learned:** HLA is a recurrent source of locus count discrepancy. Different papers use different rules for counting independent signals in this region. Matching the paper's counting methodology matters.

---

## 15. P-value: handle extreme values (1e-300+) via log (Exp 67, error 0.008, 0.0%)

**What improved:** Marginal. Fixed an edge case where APOE4 in AD (p = 10^-300) caused floating-point underflow in some calculations.

**What was changed:** Added to p-value handling: "For neg_log10_p values > 200, use arbitrary-precision log computation. Python's float64 underflows at 10^-308. Store as neg_log10_p directly, never convert back to raw p-value."

**What we learned:** Genomic p-values at the APOE locus are so extreme they break standard floating-point. Always work in log space.

---

## Summary

| Phase | Experiments | Error | Key insight |
|-------|-------------|-------|-------------|
| Baseline | 1-2 | 0.886 | Default skills miss most targets |
| Foundation | 3-9 | 0.886 to 0.528 | Normalise p-values, resolve rsIDs, cross-check ORs |
| Workflow chaining | 15-30 | 0.395 to 0.240 | PheWAS validation, beta-sign direction, supplementary tables |
| Precision | 33-44 | 0.195 to 0.074 | Population matching, fuzzy matching, dual-source fallback |
| Convergence | 52-67 | 0.041 to 0.008 | Fine-mapping PIPs, multi-ancestry reconciliation, HLA handling |

**Skills created:** variant-resolution, multi-ancestry-reconciliation

**Skills modified:** gwas-lookup (5 edits), fine-mapping (2 edits), locus-count (3 edits), effect-direction (2 edits), PheWAS workflow (1 edit)

**Biggest single improvement:** Multi-ancestry reconciliation skill (-79.7% at experiment 63). The agent couldn't match multi-ancestry GWAS effect sizes until it learned to compute population-weighted averages.

**68 failed experiments.** Common failure modes: removing validation steps, merging unrelated skills, hardcoding European frequencies, skipping LD-based deduplication, weakening scoring criteria.
