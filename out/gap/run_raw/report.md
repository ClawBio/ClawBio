# Rare High-Impact Variants Report

**Input**: ../out/gap/raw.vcf
**Rarity threshold**: population AF < 0.01

## 0 rare high-impact variants carried

Of 68 carried, annotated variants, 0 are high-impact (loss-of-function). Of those:

- **0 rare** with documented population frequency below 0.01 (ultra-rare AF < 0.001: 0; rare: 0)
- 0 common (documented AF at or above the threshold)
- 0 with no population-frequency data, so they cannot be confirmed rare (absence of a frequency is not evidence of rarity; many are common LoF polymorphisms)

## Scope

Counts high-impact (loss-of-function) variants annotated with molecular consequence and population frequency in the input VCF. 'Rare' requires a documented frequency below the threshold; variants with no frequency are reported separately and not called rare. Genome-wide novel LoF calling (VEP / SnpEff / bcftools csq) and a complete frequency reference (gnomAD) are out of scope for v0.1.0.

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
