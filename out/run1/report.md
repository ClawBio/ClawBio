# Abstention Ledger

Input: `challenge1-b37-segregation.tsv`  
SHA-256: `495b749ea3271e5b7e17e8d3609b3ad535fddc411ff4c636d74afa701147ac05`  
Records: **68**  
Assembly: GRCh37/b37 (contigs without chr prefix)

**61 of 68 records are withheld from the review list.** 7 survived every check we were able to run.

---

## 1. Which sample is which person

Resolved from genotypes: **1 of 24** possible assignments of sample IDs to family roles reproduces every genotype in the table. The mapping is therefore decided by the data, not taken on trust.

| Role | Sample |
|---|---|
| Son | `ISDBM322015` |
| Father | `ISDBM322016` |
| Sister | `ISDBM322017` |
| Mother | `ISDBM322018` |

## 2. Reproducing the parent-of-origin labels

We re-derive each label from the genotypes rather than reading the supplied column. A record is *paternal* when the father carries and the mother does not, *maternal* when the reverse holds, *ambiguous* when both carry, and *no carrier parent* when neither does.

| Role assignment | paternal | maternal | ambiguous | no carrier parent | disagreements with supplied label |
|---|---|---|---|---|---|
| as labelled | 30 | 38 | 0 | 0 | **0 / 68** |
| sister mother swapped | 11 | 25 | 19 | 13 | **32 / 68** |

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

**Checked against:** All 68 records have son genotype in ['0/1']; origin is inferred from which parent carries, and the data page labels the column UNPHASED.

### `SELECTION_BIAS`

**Not supported:** Any count, burden or per-gene total presented as describing this pedigree.

**Why:** The file contains only records where the son carries and exactly one parent carries. Every variant where both parents carry, and every variant the son does not carry, was removed upstream. Totals computed here describe the filter, not the family.

**Checked against:** 68/68 records have exactly one carrier parent, by construction.

### `HISTORICAL_ANNOTATION`

**Not supported:** Any claim that the supplied effect labels reflect current evidence.

**Why:** The effect annotation is legacy SnpEff against versioned RefSeq transcripts. Transcript sets and consequence rules have both changed since. The supplied labels are a historical record of what a tool said, not current annotation.

**Checked against:** EFF field is classic SnpEff format with pinned NM_ accessions.

### `NON_PROBAND_DISCLOSURE`

**Not supported:** Any per-variant report treated as being only about the designated proband.

**Why:** The son is the teaching proband, but the file states three other people's genotypes at every position. Each parent-of-origin label is a statement about a parent; each shared call is a statement about the sister. Publishing a ranking of the son publishes a partial genotype of three non-probands who are not the subject of the analysis.

**Checked against:** Sample roles resolved from genotypes: SON=ISDBM322015, FATHER=ISDBM322016, SISTER=ISDBM322017, MOTHER=ISDBM322018

## 4. Per-record gates

| Reason code | Records |
|---|---|
| `FREQUENCY_DOCUMENTED_COMMON` | 54 |
| `TRANSCRIPT_ARTIFACT` | 22 |
| `LOW_COMPLEXITY_LOCUS` | 19 |
| `ANNOTATION_SUPERSEDED` | 4 |
| `NO_FREQUENCY_RECORD` | 2 |

## 5. Review list

Ordered by how complete the assembled evidence is, **not** by how severe the consequence appears. Every entry remains subject to section 3.

| # | Variant | Gene(s) | Supplied impact | Current impact | Documented frequency |
|---|---|---|---|---|---|
| 1 | `2:24413293 CA>C` | FAM228A | HIGH | HIGH/MODIFIER | 0.001309 (gnomade_asj) |
| 2 | `4:165878621 G>A` | TRIM61, FAM218A | HIGH | HIGH/MODIFIER | 8.964e-05 (gnomade_afr) |
| 3 | `7:76069617 TG>T` | ZP3 | HIGH | HIGH/MODIFIER | 0.002091 (gnomade_fin) |
| 4 | `8:7718227 C>CT` | SPAG11A | HIGH | HIGH/MODIFIER | 0.0001849 (gnomade_mid) |
| 5 | `14:23371267 CA>C` | RBM23 | HIGH | HIGH/MODIFIER | 0.003401 (gnomadg_mid) |
| 6 | `14:23371269 GCA>G` | RBM23 | HIGH | HIGH/MODIFIER | 0.003401 (gnomadg_mid) |
| 7 | `19:23844937 CA>C` | ZNF675 | HIGH | HIGH/MODIFIER | 4.474e-05 (gnomade_amr) |

## 6. The ledger — every record withheld, and why

| Variant | Gene(s) | Codes | Evidence that fired the gate |
|---|---|---|---|
| `1:11906068 rs5065 A>G` | CLCN6, NPPA-AS1, NPPA | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4183 in afr (>= 0.01) |
| `1:12939546 C>CG` | PRAMEF4 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: PRAMEF4: PRAME family cluster at 1p36, segmental duplication<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1879 in gnomadg_mid (>= 0.01) |
| `1:20501582 rs12139100 G>A` | PLA2G2C | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.3906 in sas (>= 0.01) |
| `1:24710405 rs12090888 A>G` | STPG1 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON`, `ANNOTATION_SUPERSEDED` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — STPG1: NM_178122.4=MODERATE(NON_SYNONYMOUS_CODING); NM_001199012.1=MODERATE(NON_SYNONYMOUS_CODING); NM_001199013.1=MODERATE(NON_SYNONYMOUS_CODING); NM_001199014.1=HIGH(START_LOST)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.2277 in afr (>= 0.01)<br>`ANNOTATION_SUPERSEDED`: supplied impact HIGH is not reproduced by current annotation (now ['MODERATE', 'MODIFIER']) |
| `1:145606274 C>T` | NBPF10, POLR3C, RNF115 | `TRANSCRIPT_ARTIFACT`, `LOW_COMPLEXITY_LOCUS` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — POLR3C: NM_006468.6=MODIFIER(INTRON); NM_006468.6=HIGH(SPLICE_SITE_DONOR)<br>`LOW_COMPLEXITY_LOCUS`: NBPF10: neuroblastoma breakpoint family; DUF1220 tandem repeats |
| `1:152323132 rs12568784 G>T` | FLG2 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: FLG2: filaggrin repeat domains; copy-number variable, repeat-length artefacts<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4554 in eas (>= 0.01) |
| `1:156354347 rs11303415 TC>T` | RHBG | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — RHBG: NR_046115.1=MODIFIER(EXON); NM_001256395.1=MODIFIER(INTRON); NM_001256396.1=MODIFIER(INTRON); NM_020407.4=MODIFIER(INTRON); NM_001256395.1=HIGH(SPLICE_SITE_ACCEPTOR); NM_001256396.1=HIGH(SPLICE_SITE_ACCEPTOR); NM_020407.4=HIGH(SPLICE_SITE_ACCEPTOR); NM_001256395.1=HIGH(SPLICE_SITE_DONOR); NM_001256396.1=HIGH(SPLICE_SITE_DONOR); NM_020407.4=HIGH(SPLICE_SITE_DONOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.6059 in afr (>= 0.01) |
| `1:223285200 rs5744168 G>A` | TLR5 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1033 in sas (>= 0.01) |
| `1:248685789 T>TC` | OR2G6 | `LOW_COMPLEXITY_LOCUS` | `LOW_COMPLEXITY_LOCUS`: OR2G6: olfactory receptor family; largest human paralogous family, high cross-mapping |
| `2:44528267 rs146761937 GT>G` | SLC3A1 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — SLC3A1: NM_000341.3=MODIFIER(INTRON); NM_000341.3=HIGH(SPLICE_SITE_DONOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.2708 in afr (>= 0.01) |
| `2:61361325 rs142269591 TG>T` | KIAA1841 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1382 in gnomadg_ami (>= 0.01) |
| `2:69659126 rs4453725 A>T` | NFU1 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — NFU1: NR_045631.1=MODIFIER(EXON); NR_045632.1=MODIFIER(EXON); NM_001002755.2=MODERATE(NON_SYNONYMOUS_CODING); NM_015700.3=HIGH(START_LOST); NM_001002756.2=MODIFIER(UTR_5_PRIME)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.5539 in gnomade_asj (>= 0.01) |
| `3:56653839 CTG>C` | FAM208A, CCDC66 | `TRANSCRIPT_ARTIFACT` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — CCDC66: NR_024460.1=MODIFIER(EXON); NM_001012506.4=HIGH(FRAME_SHIFT); NM_001141947.1=HIGH(FRAME_SHIFT) |
| `3:195453017 rs144288174 GC>G` | MUC20 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: MUC20: mucin; long VNTR domains, assembly and alignment both unstable<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.191 in gnomade_asj (>= 0.01) |
| `4:88231392 rs35784587 T>TA` | HSD17B13 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON`, `ANNOTATION_SUPERSEDED` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — HSD17B13: NM_001136230.1=MODIFIER(INTRON); NM_178135.3=MODIFIER(INTRON); NM_001136230.1=HIGH(SPLICE_SITE_DONOR); NM_178135.3=HIGH(SPLICE_SITE_DONOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.341 in gnomadg_eas (>= 0.01)<br>`ANNOTATION_SUPERSEDED`: supplied impact HIGH is not reproduced by current annotation (now ['LOW']) |
| `5:134782450 rs12520799 T>A` | TIFAB, C5orf20 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.9319 in afr (>= 0.01) |
| `5:135272373 rs35444976 C>CA` | FBXL21 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.108 in gnomade_mid (>= 0.01) |
| `5:140772898 rs3214276 GC>G` | PCDHGB4, PCDHGA8, PCDHGB1, PCDHGB3, PCDHGA1, PCDHGA4, PCDHGA5, PCDHGB2, PCDHGA2, PCDHGA3, PCDHGA6, PCDHGA7, PCDHGB5 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: PCDHGB4: clustered protocadherin locus, tandem near-identical cassettes<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.184 in gnomade_afr (>= 0.01) |
| `5:180582256 rs140598308 TTGTC>T` | OR2V2 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR2V2: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.2022 in gnomadg_ami (>= 0.01) |
| `6:31106500 rs111966729 T>TC` | CCHCR1, PSORS1C1, PSORS1C2 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: CCHCR1: within the 6p21 MHC block; reference bias<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.3271 in gnomade_asj (>= 0.01) |
| `6:31124849 rs115419545 C>T` | CCHCR1, TCF19 | `TRANSCRIPT_ARTIFACT`, `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — CCHCR1: NM_001105563.1=HIGH(STOP_GAINED); NM_001105564.1.4=HIGH(STOP_GAINED); NM_019052.3=MODIFIER(UTR_5_PRIME)<br>`LOW_COMPLEXITY_LOCUS`: CCHCR1: within the 6p21 MHC block; reference bias<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.6075 in gnomadg_asj (>= 0.01) |
| `6:31125257 rs72856718 C>A` | CCHCR1, TCF19 | `TRANSCRIPT_ARTIFACT`, `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — CCHCR1: NM_019052.3=MODIFIER(INTRON); NM_001105563.1=HIGH(STOP_GAINED); NM_001105564.1.4=HIGH(STOP_GAINED)<br>`LOW_COMPLEXITY_LOCUS`: CCHCR1: within the 6p21 MHC block; reference bias<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1793 in gnomadg_fin (>= 0.01) |
| `6:132203615 G>A` | ENPP1 | `TRANSCRIPT_ARTIFACT` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — ENPP1: NM_006208.2=MODIFIER(INTRON); NM_006208.2=HIGH(SPLICE_SITE_DONOR) |
| `6:154567863 rs34427887 C>T` | OPRM1, IPCEF1 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1153 in gnomadg_asj (>= 0.01) |
| `7:44610376 G>A` | DDX56 | `TRANSCRIPT_ARTIFACT` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — DDX56: NM_001257189.1=MODIFIER(INTRON); NM_019082.3=HIGH(STOP_GAINED) |
| `8:39872935 rs4503083 T>A` | IDO2 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.3783 in gnomadg_ami (>= 0.01) |
| `9:125391241 rs1476860 G>A` | OR1B1 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR1B1: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.5989 in gnomade_eas (>= 0.01) |
| `9:125391770 rs72541201 C>CA` | OR1B1 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR1B1: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.689 in gnomadg_ami (>= 0.01) |
| `10:27702256 rs112067123 G>GC` | PTCHD3 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4275 in gnomadg_ami (>= 0.01) |
| `10:32740798 C>CT` | CCDC7 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.03401 in gnomadg_mid (>= 0.01) |
| `10:33136818 TAA>T` | C10orf68 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1118 in gnomadg_ami (>= 0.01) |
| `10:124214355 rs2736911 C>T` | ARMS2 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1905 in eas (>= 0.01) |
| `11:4389404 rs74427348 AG>A` | OR52B4 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR52B4: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4893 in gnomade_amr (>= 0.01) |
| `11:4790873 rs34672924 CG>C` | OR51F1 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR51F1: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4569 in afr (>= 0.01) |
| `11:5877979 rs12419602 T>A` | OR52E8 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR52E8: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.5122 in gnomadg_ami (>= 0.01) |
| `11:48266736 rs7120775 C>G` | OR4X2 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR4X2: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.2882 in afr (>= 0.01) |
| `11:48286231 rs10838851 T>A` | OR4X1 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR4X1: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.8355 in gnomadg_ami (>= 0.01) |
| `11:61165731 rs11382548 C>CA` | CPSF7, TMEM216 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — TMEM216: NM_001173990.2=MODIFIER(INTRON); NM_001173991.2=MODIFIER(INTRON); NM_016499.5=MODIFIER(INTRON); NM_001173991.2=HIGH(SPLICE_SITE_ACCEPTOR); NM_016499.5=HIGH(SPLICE_SITE_ACCEPTOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.9762 in eas (>= 0.01) |
| `11:61165741 rs10897158 G>C` | CPSF7, TMEM216 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — TMEM216: NM_001173990.2=MODIFIER(INTRON); NM_001173991.2=MODERATE(NON_SYNONYMOUS_CODING); NM_016499.5=MODERATE(NON_SYNONYMOUS_CODING); NM_001173990.2=HIGH(SPLICE_SITE_ACCEPTOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.9714 in gnomade_eas (>= 0.01) |
| `11:118949277 rs36008744 C>G` | VPS11 | `FREQUENCY_DOCUMENTED_COMMON`, `ANNOTATION_SUPERSEDED` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1301 in afr (>= 0.01)<br>`ANNOTATION_SUPERSEDED`: supplied impact HIGH is not reproduced by current annotation (now ['LOW', 'MODIFIER']) |
| `12:51723598 rs60311818 A>AG` | CELA1 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.3914 in gnomadg_ami (>= 0.01) |
| `12:54734289 rs11170877 A>G` | MIR148B, COPZ1 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — COPZ1: NM_001271734.1=MODIFIER(INTRON); NM_001271735.1=MODIFIER(INTRON); NM_016057.2=MODIFIER(INTRON); NR_073424.1=MODIFIER(INTRON); NM_001271736.1=HIGH(START_LOST)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.2554 in gnomade_eas (>= 0.01) |
| `12:55820958 rs57387180 CA>C` | OR6C76 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: OR6C76: olfactory receptor family; largest human paralogous family, high cross-mapping<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.3536 in gnomadg_sas (>= 0.01) |
| `13:31531009 rs12857479 G>A` | TEX26 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — TEX26: NM_152325.1=MODIFIER(INTRON); NM_152325.1=HIGH(SPLICE_SITE_ACCEPTOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4114 in gnomade_mid (>= 0.01) |
| `14:51378590 rs11356035 CT>C` | PYGL | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — PYGL: NM_001163940.1=MODIFIER(INTRON); NM_002863.4=MODIFIER(INTRON); NM_001163940.1=HIGH(SPLICE_SITE_ACCEPTOR); NM_002863.4=HIGH(SPLICE_SITE_ACCEPTOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.7133 in afr (>= 0.01) |
| `14:96730313 G>A` | BDKRB1 | `NO_FREQUENCY_RECORD` | `NO_FREQUENCY_RECORD`: build-matched query returned no population frequency for this position; this is a statement about the reference database, not about the variant |
| `15:31294714 rs3784589 C>A` | TRPM1 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1316 in gnomadg_ami (>= 0.01) |
| `16:1306559 rs3830782 AG>A` | TPSD1 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4191 in afr (>= 0.01) |
| `16:81056441 rs3743503 T>G` | CENPN | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — CENPN: NM_001270474.1=MODIFIER(INTRON); NM_001270473.1=MODIFIER(INTRON); NM_001100624.2=MODIFIER(INTRON); NM_001100625.2=MODIFIER(INTRON); NM_018455.5=HIGH(STOP_LOST)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.2093 in gnomadg_eas (>= 0.01) |
| `17:42979026 rs74349463 T>C` | GFAP, CCDC103, EFTUD2 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON`, `ANNOTATION_SUPERSEDED` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — CCDC103: NM_001258398.1=MODIFIER(INTRON); NM_001258399.1=MODIFIER(INTRON); NM_001258395.1=MODIFIER(INTRON); NM_001258396.1=MODIFIER(INTRON); NM_213607.2=MODIFIER(INTRON); NM_001258397.1=MODIFIER(INTRON); NM_001258398.1=HIGH(SPLICE_SITE_DONOR); NM_001258399.1=HIGH(SPLICE_SITE_DONOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1995 in gnomade_mid (>= 0.01)<br>`ANNOTATION_SUPERSEDED`: supplied impact HIGH is not reproduced by current annotation (now ['LOW', 'MODIFIER']) |
| `17:45468858 rs118004742 T>G` | EFCAB13 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1497 in gnomadg_mid (>= 0.01) |
| `17:72588806 rs545652 C>A` | C17orf77, CD300LD | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.5189 in afr (>= 0.01) |
| `19:6897464 rs330880 C>G` | EMR1 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — EMR1: NM_001256255.1=MODIFIER(INTRON); NM_001256254.1=MODIFIER(INTRON); NM_001256253.1=MODERATE(NON_SYNONYMOUS_CODING); NM_001974.4=MODERATE(NON_SYNONYMOUS_CODING); NM_001256252.1=MODERATE(NON_SYNONYMOUS_CODING); NM_001256255.1=HIGH(SPLICE_SITE_DONOR)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.8328 in afr (>= 0.01) |
| `19:52004791 rs67024588 G>GC` | SIGLEC12 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — SIGLEC12: NM_053003.2=HIGH(FRAME_SHIFT); NM_033329.1=MODIFIER(UPSTREAM)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.7642 in gnomade_amr (>= 0.01) |
| `19:52803669 rs3217319 CTG>C` | ZNF480 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.7836 in gnomadg_fin (>= 0.01) |
| `19:55019261 rs61737751 C>T` | LAIR2 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.05782 in gnomadg_mid (>= 0.01) |
| `20:31756954 rs17124277 G>A` | BPIFA2 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.1068 in gnomade_afr (>= 0.01) |
| `21:11029596 rs138714104 AC>A` | BAGE2, BAGE3, BAGE4, BAGE5 | `TRANSCRIPT_ARTIFACT`, `NO_FREQUENCY_RECORD` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — BAGE2: NM_182482.2=MODIFIER(INTRON); NM_182482.2=HIGH(SPLICE_SITE_DONOR) | BAGE3: NM_182481.1=MODIFIER(INTRON); NM_182481.1=HIGH(SPLICE_SITE_DONOR)<br>`NO_FREQUENCY_RECORD`: build-matched query returned no population frequency for this position; this is a statement about the reference database, not about the variant |
| `21:31744127 rs877346 A>T` | MIR4327, KRTAP13-2 | `LOW_COMPLEXITY_LOCUS`, `FREQUENCY_DOCUMENTED_COMMON` | `LOW_COMPLEXITY_LOCUS`: KRTAP13-2: keratin-associated proteins; clustered, highly similar paralogues<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.6509 in gnomadg_ami (>= 0.01) |
| `22:32875190 rs11107 G>A` | FBXO7 | `TRANSCRIPT_ARTIFACT`, `FREQUENCY_DOCUMENTED_COMMON` | `TRANSCRIPT_ARTIFACT`: supplied max impact HIGH is not reproduced across transcripts of the same gene — FBXO7: NM_012179.3=MODERATE(NON_SYNONYMOUS_CODING); NM_001033024.1=MODERATE(NON_SYNONYMOUS_CODING); NM_001257990.1=HIGH(START_LOST)<br>`FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.7363 in gnomade_eas (>= 0.01) |
| `22:41257834 T>TA` | DNAJB7, XPNPEP3 | `FREQUENCY_DOCUMENTED_COMMON` | `FREQUENCY_DOCUMENTED_COMMON`: documented frequency 0.4475 in gnomade_amr (>= 0.01) |

## 7. Limits of the evidence layer we added

- Source: gnomAD exomes r2.1.1 via Ensembl GRCh37 REST (exome-only)
- Endpoint: `https://grch37.rest.ensembl.org/vep/human/region` (GRCh37), 68 records returned
- The layer is **build-matched**, so an empty frequency result means Ensembl holds no record at that position. It does not mean a liftover failed, because there is no liftover in this path.
- It is nonetheless **exome-only and several releases old**. gnomAD v4 is GRCh38-native and is not used here. `NO_FREQUENCY_RECORD` is therefore a weaker statement than it would be against v4, and we do not treat it as evidence about the variant's frequency in any population.

---

*Any tool can output a ranking. This one ships the check, the value and the source for everything it declined to rank.*
