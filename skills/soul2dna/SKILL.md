---
name: soul2dna
description: Compile SOUL.md character profiles into synthetic diploid genomes (.genome.json) via trait-to-allele mapping
license: MIT
metadata:
  version: 0.1.0
  author: Manuel Corpas
  tags:
  - genomebook
  - synthetic-genomics
  - soul-compiler
  - trait-mapping
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🧬
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    trigger_keywords:
    - soul2dna
    - soul compiler
    - soul to genome
    - genomebook compile
    - synthetic genome
    - character genome
---

# 🧬 Soul2DNA Compiler

## Purpose

Compile SOUL.md character profiles into synthetic diploid genomes. Each soul file
describes a historical or fictional figure with trait scores (0.0 to 1.0). The
compiler maps these scores to alleles at defined loci using additive, dominant, or
recessive inheritance models, producing a `.genome.json` file per character.

## How It Works

1. **Parse SOUL.md** files from `GENOMEBOOK/DATA/SOULS/` extracting identity
   metadata (name, sex, ancestry, domain, era) and trait scores.
2. **Load trait registry** (`GENOMEBOOK/DATA/trait_registry.json`) which defines
   loci, alleles, chromosomal positions, dominance models, and effect sizes for
   each trait.
3. **Assign genotypes** at each locus based on trait score thresholds:
   - Additive: <0.33 ref/ref, 0.33-0.66 ref/alt, >0.66 alt/alt
   - Dominant: <0.40 ref/ref, 0.40-0.75 ref/alt, >0.75 alt/alt
   - Recessive: <0.50 ref/ref, 0.50-0.80 ref/alt, >0.80 alt/alt
4. **Write genome** as JSON with full locus detail, trait scores, and metadata.

## Input

- `GENOMEBOOK/DATA/SOULS/*.soul.md` (20 historical figures)
- `GENOMEBOOK/DATA/trait_registry.json`

## Output

- `GENOMEBOOK/DATA/GENOMES/<name>-g0.genome.json` per character

## CLI Usage

```bash
# Compile all souls to genomes
python skills/soul2dna/soul2dna.py

# Demo mode (shows summary without writing files)
python skills/soul2dna/soul2dna.py --demo
```

## Output Format

Each `.genome.json` contains:

```json
{
  "id": "einstein-g0",
  "name": "Albert Einstein",
  "sex": "Male",
  "sex_chromosomes": "XY",
  "ancestry": "...",
  "generation": 0,
  "parents": [null, null],
  "loci": { "<locus_id>": { "chromosome": "...", "alleles": ["A","G"], ... } },
  "trait_scores": { "curiosity": 0.95, ... }
}
```
