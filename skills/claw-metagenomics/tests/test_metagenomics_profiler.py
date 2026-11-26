"""Tests for metagenomics_profiler.py."""

from __future__ import annotations

import pytest
from inline_snapshot import snapshot

from clawbio.common import demo_runner

# ── Demo output snapshot ──────────────────────────────────────────────────────
# Fails when a code change alters any demo number. If the change is intended,
# run `uv run pytest skills/claw-metagenomics/tests --inline-snapshot=fix`; it
# rewrites the snapshot() values below, and the PR diff shows which numbers moved.


@pytest.mark.demo_snapshot
def test_demo_output_matches_snapshot(tmp_path):
    assert demo_runner.demo_outputs("claw-metagenomics", tmp_path) == snapshot({'report.md': """\
# Metagenomics Profiling Report

**Date**: <TIMESTAMP> UTC
**Tool**: ClawBio Metagenomics Profiler v0.1.0
**Mode**: Demo (pre-computed Peru sewage data)

---

## Taxonomic Composition

![taxonomy](figures/taxonomy_barplot.png)

### Top Species

  - *Escherichia coli*
  - *Klebsiella pneumoniae*
  - *Pseudomonas aeruginosa*
  - *Acinetobacter baumannii*
  - *Enterococcus faecium*

## Antimicrobial Resistance Profile

![resistome](figures/resistome_heatmap.png)

Total ARG detections: **24**
- WHO-Critical priority: **7**
- WHO-High priority: **7**

### WHO-Critical ARG Summary

![who_critical](figures/who_critical_args.png)

## Functional Pathways (HUMAnN3)

Total pathways detected: **10**

Top pathway: PWY-7219: adenosine ribonucleotides de novo biosynthesis

---

## Methods

- **Taxonomic classification**: Kraken2 (confidence threshold 0.2) with Bracken species-level re-estimation (minimum 10 reads)
- **Resistome profiling**: RGI against the CARD database (Perfect + Strict criteria only)
- **Functional profiling**: HUMAnN3 with MetaCyc pathway stratification
- **WHO Priority mapping**: ARGs classified against the 2024 WHO Bacterial Priority Pathogens List
- **Figures**: 300 dpi, publication-quality

## Input Checksums

  - `demo_peru_sewage.fastq.gz`: demo data (no file)

## Disclaimer

> ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions. Antimicrobial resistance findings must be confirmed by culture-based susceptibility testing.

## Citations

- Wood, D.E., Lu, J. & Langmead, B. (2019). Improved metagenomic analysis with Kraken 2. *Genome Biology*, 20, 257.
- Lu, J. et al. (2017). Bracken: estimating species abundance in metagenomics data. *PeerJ Computer Science*, 3, e104.
- Alcock, B.P. et al. (2023). CARD 2023: expanded curation, support for machine learning, and resistome prediction at the Comprehensive Antibiotic Resistance Database. *Nucleic Acids Research*, 51(D1), D419-D430.
- Beghini, F. et al. (2021). Integrating taxonomic, functional, and strain-level profiling of diverse microbial communities with bioBakery 3. *eLife*, 10, e65088.
- Corpas, M. (2026). ClawBio. https://github.com/ClawBio/ClawBio
""", 'tables/pathway_abundance.tsv': """\
pathway	mean_abundance
PWY-7219: adenosine ribonucleotides de novo biosynthesis	0.0182
PWY-7221: guanosine ribonucleotides de novo biosynthesis	0.0165
PWY-6123: inosine-5'-phosphate biosynthesis I	0.0148
GLYCOLYSIS: glycolysis I (from glucose 6-phosphate)	0.0137
PWY-5100: pyruvate fermentation to acetate and lactate II	0.0121
TCA: TCA cycle I (prokaryotic)	0.0118
PWY-6163: chorismate biosynthesis from 3-dehydroquinate	0.0095
PWY-6122: 5-aminoimidazole ribonucleotide biosynthesis II	0.0089
PANTOSYN-PWY: pantothenate and coenzyme A biosynthesis I	0.0082
PWY-7234: inosine-5'-phosphate biosynthesis III	0.0076
""", 'tables/resistome_profile.tsv': """\
gene	drug_class	criteria	Lima	Cusco	Iquitos
NDM-1	carbapenem	Perfect	18	4	25
OXA-48	carbapenem	Strict	12	2	18
KPC-3	carbapenem	Strict	8	1	14
CTX-M-15	cephalosporin	Perfect	45	18	58
CTX-M-27	cephalosporin	Strict	22	8	31
mecA	methicillin	Strict	5	2	8
VanA	vancomycin	Perfect	3	1	5
VanB	vancomycin	Strict	1	0	3
QnrB	fluoroquinolone	Strict	15	6	22
QnrS	fluoroquinolone	Strict	9	3	14
TEM-1	penicillin	Perfect	52	31	67
SHV-12	cephalosporin	Strict	28	12	35
GyrA_S83L	fluoroquinolone	Strict	11	4	18
ParC_S80I	fluoroquinolone	Strict	7	2	12
aac(6')-Ib	aminoglycoside	Strict	19	8	27
aph(3'')-Ib	aminoglycoside	Strict	14	6	21
sul1	sulfonamide	Perfect	38	22	51
sul2	sulfonamide	Perfect	29	15	38
dfrA12	trimethoprim	Strict	16	9	23
tet(A)	tetracycline	Perfect	33	19	45
erm(B)	macrolide	Strict	12	7	18
mph(A)	macrolide	Strict	8	3	12
catA1	phenicol	Perfect	6	2	9
floR	phenicol	Strict	4	1	7
""", 'tables/taxonomy_species.tsv': """\
species	Lima	Cusco	Iquitos
Escherichia coli	12.3	8.1	15.6
Klebsiella pneumoniae	8.7	5.2	11.3
Pseudomonas aeruginosa	5.1	3.8	7.2
Acinetobacter baumannii	3.9	2.1	5.8
Enterococcus faecium	2.8	1.5	3.9
Bacteroides fragilis	4.5	6.2	3.1
Salmonella enterica	1.9	0.8	3.2
Staphylococcus aureus	1.2	0.5	2.1
Clostridioides difficile	0.8	1.2	0.5
Campylobacter jejuni	0.6	0.9	1.4
Prevotella copri	3.2	5.8	2.1
Faecalibacterium prausnitzii	5.1	7.3	3.8
Bifidobacterium longum	2.1	3.4	1.5
Lactobacillus rhamnosus	1.0	2.1	0.8
Helicobacter pylori	0.4	0.7	1.1
Vibrio cholerae	0.3	0.1	1.8
Shigella flexneri	1.1	0.4	2.3
Proteus mirabilis	0.9	1.3	1.1
Serratia marcescens	0.7	0.6	1.2
Citrobacter freundii	0.5	0.8	0.9
""", 'tables/who_priority_args.tsv': """\
gene	drug_class	who_priority
NDM-1	carbapenem	Critical
OXA-48	carbapenem	Critical
KPC-3	carbapenem	Critical
CTX-M-15	cephalosporin	Critical
CTX-M-27	cephalosporin	Critical
mecA	methicillin	High
VanA	vancomycin	High
VanB	vancomycin	High
QnrB	fluoroquinolone	High
QnrS	fluoroquinolone	High
TEM-1	penicillin	Critical
SHV-12	cephalosporin	Critical
GyrA_S83L	fluoroquinolone	High
ParC_S80I	fluoroquinolone	High
aac(6')-Ib	aminoglycoside	Not classified
aph(3'')-Ib	aminoglycoside	Not classified
sul1	sulfonamide	Not classified
sul2	sulfonamide	Not classified
dfrA12	trimethoprim	Not classified
tet(A)	tetracycline	Not classified
erm(B)	macrolide	Not classified
mph(A)	macrolide	Not classified
catA1	phenicol	Not classified
floR	phenicol	Not classified
"""})
