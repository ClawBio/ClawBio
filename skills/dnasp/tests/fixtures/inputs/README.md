# Test input provenance

These inputs drive the comparisons with DnaSP 6 (Supplementary Tables S1 and S2
of the DnaSP-ClawBio manuscript) and the regression tests built on them.

## DnaSP example data

Example files from the DnaSP distribution, supplied by Julio Rozas (DnaSP 6
authors) for this reimplementation. Credit: Rozas J, Ferrer-Mata A,
Sanchez-DelBarrio JC, et al. (2017) DnaSP 6: DNA Sequence Polymorphism Analysis
of Large Data Sets. Mol Biol Evol 34:3299-3302.

| File | Origin |
|------|--------|
| `rp49_36.nex` | DnaSP example `rp49#36.nex` (*Drosophila subobscura* rp49 region, 36 sequences) |
| `COII_Apes.nex` | DnaSP example `COII_Apes.nex` (hominoid mitochondrial COII) |
| `DmelOSRegion.nex` | DnaSP example `DmelOSRegion.nex` (*D. melanogaster* Os region) |
| `Ex_n1.fas` | DnaSP example `Ex##n1.nbr` converted by DnaSP 3.37 (four sequences) |
| `vcf/Data_Example_DiploidPhased.vcf`, `vcf/Data_Example_DiploidUnphased.vcf`, `vcf/Data_Example_Haploid.vcf`, `vcf/Data_Example_vcf.SG.txt` | DnaSP VCF examples, unchanged |
| `rp49_5regions/Region_1.fas` to `Region_5.fas` | The five regions of DnaSP example `rp49_5regions.fa`, one file per region |

## Derived from the example data

Same sequence content as the source example; only records or columns were selected.

| File | Derived from |
|------|--------------|
| `rp49_ing34_outGUA.fas` | `rp49_36.nex`: 34 ingroup sequences and the outgroup `rp49.gua` |
| `COII_HsaPtr_outPpy.fas` | `COII_Apes.nex`: human and chimpanzee sequences with the outgroup `Ppy1` |
| `COII_Hsa_outPpy.fas` | `COII_Apes.nex`: human sequences with the outgroup `Ppy1` |
| `COII_coding_1_681.fas` | `COII_HsaPtr_outPpy.fas` cropped to positions 1-681, the coding region without the terminal stop codon |
| `rp49_36_pops.txt`, `COII_Apes_HsaPtr_pops.txt` | Population assignments written for these tests from the example sequence names |

## DnaSP help worked examples

| File | Origin |
|------|--------|
| `indel_help.fas` | The worked example in DnaSP's help topic on InDel polymorphism (13 sequences, 18 sites) |
| `ld_help.fas` | The worked example in DnaSP's help topic on linkage disequilibrium |

## Constructed for these tests

| File | Purpose |
|------|---------|
| `stop_family.fas` | Synthetic sequences for the RSCU stop-codon family |
| `hka_synthetic.tsv` | Two synthetic loci built to be neutral (chi-square 0) for the HKA direct mode |

DnaSP output captures used as evidence are described in
`../validation5/README.md`.
