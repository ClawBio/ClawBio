# Genomebook

**Genomebook** is a synthetic-genetics sandbox built into ClawBio. It turns fictional or historical character profiles ("souls") into diploid genomes, scores compatibility, and breeds offspring across generations, complete with Mendelian inheritance, de novo mutations, and clinical evaluation.

## Pipeline

```
SOUL.md  -->  Soul2DNA  -->  .genome.json  -->  GenomeMatch  -->  Recombinator  -->  Gen-N offspring
 (trait       (compiler)     (diploid loci)     (M x F rank)     (meiosis +        (.genome.json
  scores)                                                         mutation)          with clinical
                                                                                     history)
```

## Data

- **20 souls** in `GENOMEBOOK/DATA/SOULS/` (Einstein, Curie, Turing, Hypatia, Da Vinci, ...)
- **20 generation-0 genomes** in `GENOMEBOOK/DATA/GENOMES/`
- **26 traits** across 55 loci (trait_registry.json)
- **Disease registry** with penetrance, fitness costs, and onset probabilities

## Quick Start

```bash
# Compile all souls to genomes (Genomebook sandbox script; not a catalogued skill)
python GENOMEBOOK/PYTHON/01-soul2dna.py

# Score all M x F compatibility pairings
python skills/genome-match/genome_match.py --demo

# Breed Einstein x Anning (3 offspring)
python skills/recombinator/recombinator.py --demo
```

## What It Models

- **Additive, dominant, recessive** inheritance at each locus
- **Heterozygosity advantage** (heterosis) in compatibility scoring
- **Carrier risk** flagging for autosomal recessive conditions
- **De novo mutations** with hotspot categories (cognitive, immune, metabolic)
- **Clinical evaluation** with penetrance, onset probability, and fitness costs
- **Health scores** derived from cumulative condition burden
