# Abstention Ledger

Input: `demo_segregation.tsv`  
SHA-256: `1fd4bdf6fefa517c3438ace9b2b9ec267b8b218e6dcb5628417250c10c52e1d3`  
Records: **9**  
Assembly: GRCh37/b37 (contigs without chr prefix)

**5 of 9 records are withheld from the review list.** 4 survived every check we were able to run.

---

## 1. Which sample is which person

Status: **not_checked** — no VCF supplied; role labels taken from TSV column names on trust

## 2. Reproducing the parent-of-origin labels

We re-derive each label from the genotypes rather than reading the supplied column. A record is *paternal* when the father carries and the mother does not, *maternal* when the reverse holds, *ambiguous* when both carry, and *no carrier parent* when neither does.

| Role assignment | paternal | maternal | ambiguous | no carrier parent | disagreements with supplied label |
|---|---|---|---|---|---|
| as labelled | 4 | 5 | 0 | 0 | **0 / 9** |
| sister mother swapped | 3 | 3 | 1 | 2 | **3 / 9** |

The second row is what happens if the sister and mother columns are exchanged. It is not a hypothetical: the published description of this dataset orders the samples in a way that produces exactly that swap. Records with no carrier parent are impossible under the dataset's own stated filter, which is how the mistake announces itself.

## 3. What this data cannot support, whatever the records say

These are properties of the dataset. They bound every claim below.

### `NO_PHENOTYPE`

**Not supported:** Any statement linking a variant to this individual's clinical picture.

**Why:** The input carries no phenotype and no HPO terms. Segregation without a phenotype is arithmetic over genotypes; it cannot be evidence for or against a condition nobody has described. Gene-to-condition links we could look up would be untested against this person.

**Checked against:** Input columns: no phenotype field, no HPO field. Confirmed by the data page.

### `UNPHASED`

**Not supported:** Any statement about two variants being on the same or opposite copies.

**Why:** Parent-of-origin here is Mendelian deduction from genotypes, not molecular phase. Without phase, the two-hit question about a gene is not merely unanswered — it is unaskable from this file.

**Checked against:** All 9 records have son genotype in ['0/1']; origin is inferred from which parent carries, and the data page labels the column UNPHASED.

### `SELECTION_BIAS`

**Not supported:** Any count, burden or per-gene total presented as describing this pedigree.

**Why:** The file contains only records where the son carries and exactly one parent carries. Every variant where both parents carry, and every variant the son does not carry, was removed upstream. Totals computed here describe the filter, not the family.

**Checked against:** 9/9 records have exactly one carrier parent, by construction.

### `HISTORICAL_ANNOTATION`

**Not supported:** Any claim that the supplied effect labels reflect current evidence.

**Why:** The effect annotation is legacy SnpEff against versioned RefSeq transcripts. Transcript sets and consequence rules have both changed since. The supplied labels are a historical record of what a tool said, not current annotation.

**Checked against:** EFF field is classic SnpEff format with pinned NM_ accessions.

### `NON_PROBAND_DISCLOSURE`

**Not supported:** Any per-variant report treated as being only about the designated proband.

**Why:** The son is the teaching proband, but the file states three other people's genotypes at every position. Each parent-of-origin label is a statement about a parent; each shared call is a statement about the sister. Publishing a ranking of the son publishes a partial genotype of three non-probands who are not the subject of the analysis.

**Checked against:** Sample roles resolved from genotypes: SON=SON, FATHER=FATHER, SISTER=SISTER, MOTHER=MOTHER

## 4. Per-record gates

| Reason code | Records |
|---|---|
| `LOW_COMPLEXITY_LOCUS` | 3 |
| `TRANSCRIPT_ARTIFACT` | 2 |
| `SECONDARY_FINDING_NO_CONSENT` | 1 |

### Checks that ran and found nothing

- **ACMG secondary findings:** all 9 records screened against 81 genes — **1 hits**. Source: ACMG_SF_V32_GENES from skills/clinical-variant-reporter/acmg_engine.py (81 genes, ACMG SF v3.2).

## 5. Review list

Ordered by how complete the assembled evidence is, **not** by how severe the consequence appears. Every entry remains subject to section 3.

| # | Variant | Gene(s) | Supplied impact | Current impact | Documented frequency |
|---|---|---|---|---|---|
| 1 | `1:1000100 rsDEMO1 C>T` | DEMOGENE1 | HIGH | not checked | not checked |
| 2 | `1:1000200 rsDEMO2 G>A` | DEMOGENE2 | HIGH | not checked | not checked |
| 3 | `5:5000700 rsDEMO7 A>T` | DEMOGENE7 | HIGH | not checked | not checked |
| 4 | `7:7000800 CG>C` | DEMOGENE8 | HIGH | not checked | not checked |

## 6. The ledger — every record withheld, and why

| Variant | Gene(s) | Codes | Evidence that fired the gate |
|---|---|---|---|
| `2:2000300 rsDEMO3 A>G` | DEMOGENE3 | `TRANSCRIPT_ARTIFACT` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — DEMOGENE3: NM_000003.9=HIGH(START_LOST); NM_000003.1=MODERATE(NON_SYNONYMOUS_CODING) |
| `3:3000400 T>TA` | OR2DEMO1 | `LOW_COMPLEXITY_LOCUS` | `LOW_COMPLEXITY_LOCUS`: OR2DEMO1: olfactory receptor family; largest human paralogous family, high cross-mapping |
| `6:31000500 rsDEMO5 C>A` | DEMOGENE5 | `LOW_COMPLEXITY_LOCUS` | `LOW_COMPLEXITY_LOCUS`: 6:31000500 lies in the extended MHC (6:28477797-33448354, GRCh37); reference-allele bias |
| `4:4000600 rsDEMO6 G>C` | MUCDEMO2 | `TRANSCRIPT_ARTIFACT`, `LOW_COMPLEXITY_LOCUS` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — MUCDEMO2: NM_000006.2=HIGH(FRAME_SHIFT); NM_000006.3=MODIFIER(INTRON)<br>`LOW_COMPLEXITY_LOCUS`: MUCDEMO2: mucin; long VNTR domains, assembly and alignment both unstable |
| `13:13000900 rsDEMO9 C>T` | BRCA2 | `SECONDARY_FINDING_NO_CONSENT` | `SECONDARY_FINDING_NO_CONSENT`: BRCA2 on the ACMG secondary-findings list; returning a finding here would be opportunistic screening of a named person with no phenotype, no clinical relationship and no documented opt-out |

## 7. Limits of the evidence layer we added

- Status: **absent**.  Records carry `NO_EVIDENCE_LAYER` rather than an assumed value.

---

*Any tool can output a ranking. This one ships the check, the value and the source for everything it declined to rank.*

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
