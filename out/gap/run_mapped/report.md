# Rare High-Impact Variants Report

**Input**: ../out/gap/mapped.vcf
**Rarity threshold**: population AF < 0.01

## 12 rare high-impact variants carried

Of 68 carried, annotated variants, 68 are high-impact (loss-of-function). Of those:

- **12 rare** with documented population frequency below 0.01 (ultra-rare AF < 0.001: 7; rare: 5)
- 54 common (documented AF at or above the threshold)
- 2 with no population-frequency data, so they cannot be confirmed rare (absence of a frequency is not evidence of rarity; many are common LoF polymorphisms)

## Variants

| Gene | Locus | Consequence | Zygosity | Population AF | ClinVar |
|------|-------|-------------|----------|---------------|---------|
| ENPP1 | 6:132203615 G>A | splice_donor | het | 9e-07 | - |
| DDX56 | 7:44610376 G>A | nonsense | het | 3e-05 | - |
| ZNF675 | 19:23844937 CA>C | frameshift | het | 4.5e-05 | - |
| FAM218A | 4:165878621 G>A | nonsense | het | 9e-05 | - |
| CCDC66 | 3:56653839 CTG>C | frameshift | het | 0.00018 | - |
| SPAG11A | 8:7718227 C>CT | frameshift | het | 0.00018 | - |
| POLR3C | 1:145606274 C>T | splice_donor | het | 0.00029 | - |
| FAM228A | 2:24413293 CA>C | frameshift | het | 0.0013 | - |
| OR2G6 | 1:248685789 T>TC | frameshift | het | 0.0016 | - |
| ZP3 | 7:76069617 TG>T | frameshift | het | 0.0021 | - |
| RBM23 | 14:23371267 CA>C | frameshift | het | 0.0034 | - |
| RBM23 | 14:23371269 GCA>G | frameshift | het | 0.0034 | - |

## Scope

Counts high-impact (loss-of-function) variants annotated with molecular consequence and population frequency in the input VCF. 'Rare' requires a documented frequency below the threshold; variants with no frequency are reported separately and not called rare. Genome-wide novel LoF calling (VEP / SnpEff / bcftools csq) and a complete frequency reference (gnomAD) are out of scope for v0.1.0.

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
