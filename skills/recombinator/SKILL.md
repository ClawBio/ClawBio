---
name: recombinator
description: Produce offspring genomes from parent pairs via meiotic recombination, mutation, and clinical evaluation
license: MIT
metadata:
  version: 0.1.0
  author: Manuel Corpas
  tags:
  - genomebook
  - recombination
  - meiosis
  - mutation
  - offspring
  - clinical-genetics
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🧪
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    trigger_keywords:
    - recombinator
    - recombination
    - offspring
    - breed
    - meiosis
    - genomebook breed
    - next generation
---

# 🧪 Recombinator

## Purpose

Produce offspring genomes from selected parent pairs via simulated meiotic
recombination. Models Mendelian segregation, de novo mutation, sex determination,
and clinical evaluation against a disease registry.

## How It Works

1. **Mendelian segregation**: one allele inherited from each parent per locus
   (random selection simulating independent assortment).
2. **De novo mutation**: configurable rate per locus (default 0.1%), with hotspot
   multipliers for cognitive, immune, and metabolic loci. Mutations are classified
   as disease-risk, protective, or neutral.
3. **Sex determination**: 50/50 coin flip (XY or XX).
4. **Trait inference**: reverse-map offspring genotype back to trait scores using
   the trait registry, accounting for dominance models.
5. **Clinical evaluation**: check offspring genotype against disease registry for
   penetrance, onset probability, and fitness cost.
6. **Health score**: computed from cumulative fitness costs of clinical conditions.

## Input

- Two parent `.genome.json` files (one Male, one Female)
- `GENOMEBOOK/DATA/trait_registry.json`
- `GENOMEBOOK/DATA/disease_registry.json`

## Output

- Offspring `.genome.json` with:
  - Inherited loci and alleles
  - Mutation log
  - Inferred trait scores
  - Clinical history
  - Health score (0.0 to 1.0)

## CLI Usage

```bash
# Demo: breed Einstein x Anning, produce 3 offspring
python skills/recombinator/recombinator.py --demo

# Breed specific parents
python skills/recombinator/recombinator.py \
  --father einstein-g0 --mother anning-g0 --offspring 3

# Custom generation number
python skills/recombinator/recombinator.py \
  --father einstein-g0 --mother curie-g0 --offspring 2 --generation 1
```

## Output Format

```
ID:     g1-001-a3f2c1
Sex:    Female (XX)
Health: 0.9500
Mutations: 1
  - COMT_Val158Met: G->A (neutral, from mother)
Conditions: 0
Top traits:
  - curiosity: 0.92
  - analytical_thinking: 0.88
  - persistence: 0.85
```
