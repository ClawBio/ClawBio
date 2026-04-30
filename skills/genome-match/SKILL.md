---
name: genome-match
description: Score genetic compatibility across all male-female pairings in a Genomebook generation
license: MIT
metadata:
  version: 0.1.0
  author: Manuel Corpas
  tags:
  - genomebook
  - compatibility
  - heterozygosity
  - disease-risk
  - mating-selection
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 💞
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    trigger_keywords:
    - genome match
    - compatibility
    - mating pairs
    - genomebook match
    - heterozygosity score
    - breeding pairs
---

# 💞 GenomeMatch

## Purpose

Score genetic compatibility between all male-female pairings in a Genomebook
generation. The engine evaluates heterozygosity advantage, disease carrier risk,
and trait complementarity to rank optimal mating pairs for the next generation.

## How It Works

1. **Load genomes** for a target generation from `GENOMEBOOK/DATA/GENOMES/`.
2. **Compute pairwise compatibility** for every M x F combination:
   - **Heterozygosity score (40%)**: fraction of loci where offspring would be
     heterozygous (genetic diversity advantage).
   - **Trait complementarity (40%)**: reward balanced trait combinations and high
     average trait values across the pair.
   - **Disease risk penalty (20%)**: flag pairs where both parents carry recessive
     disease alleles (25% affected offspring risk per flagged condition).
3. **Rank all pairings** by composite score (0.0 to 1.0).
4. **Select non-overlapping mating pairs** via greedy selection from the top of
   the ranked list (each individual mates at most once per generation).

## Input

- `GENOMEBOOK/DATA/GENOMES/*.genome.json`
- `GENOMEBOOK/DATA/disease_registry.json`

## Output

- Ranked compatibility table (all M x F pairings)
- Selected mating pairs for the next generation

## CLI Usage

```bash
# Score all pairings for generation 0
python skills/genome-match/genome_match.py

# Score a specific generation
python skills/genome-match/genome_match.py --generation 1

# Demo mode
python skills/genome-match/genome_match.py --demo

# Limit output to top N pairings
python skills/genome-match/genome_match.py --top 10
```

## Output Format

```
Rank          Male x Female              Score   Het   Comp   Risk  Flags
   1      einstein-g0 x curie-g0         0.8234  0.650  0.821  0.000  --
   2      darwin-g0   x franklin-g0      0.7891  0.600  0.790  0.000  --
...

SELECTED MATING PAIRS (generation 0 -> 1):
  Albert Einstein x Marie Curie  (compat: 0.8234)
  Charles Darwin x Rosalind Franklin  (compat: 0.7891)
```
