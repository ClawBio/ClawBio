# Genomebook

**Genotype-driven agent reproduction.** Agents that evolve through Mendelian inheritance, not cloning.

Every AI agent has a SOUL.md (personality, traits, goals). Genomebook adds a second layer: DNA. A compiler infers biologically coherent diploid genomes from each agent's traits, encoding 26 characteristics across 55 loci using additive, dominant, and recessive inheritance models. Agents reproduce sexually, producing non-identical offspring that inherit blended traits from both parents.

The system is seeded with 20 SOUL.md profiles inspired by scientists and inventors who shaped humanity (Einstein, Curie, Turing, Da Vinci, Darwin, Franklin, and others), turning intellectual history into a living, recombinable genetic substrate.

**[Preprint (bioRxiv)](https://www.biorxiv.org/cgi/content/short/2026.03.22.713494v1)** | **[v0.1 Release](https://github.com/ClawBio/ClawBio/releases/tag/genomebook-v0.1)** | **[Live Demo](https://clawbio.github.io/ClawBio/slides/genomebook/demo.html)** | **[Slides](https://clawbio.github.io/ClawBio/slides/genomebook/)** | **[Phylogeny](https://clawbio.github.io/ClawBio/slides/genomebook/phylogeny.html)** | **[PCA](https://clawbio.github.io/ClawBio/slides/genomebook/pca.html)**

---

## Pipeline

```
SOUL.md  -->  Soul2DNA  -->  .genome.json  -->  GenomeMatch  -->  Recombinator  -->  Offspring
 (traits)    (compiler)     (diploid loci)    (M x F rank)     (meiosis +        (heritable
                                                                mutation)          variation)
```

## Quick Start

```bash
# Compile 20 founders into genomes
python PYTHON/01-soul2dna.py

# Generate DNA.md identity documents
python PYTHON/06-dna_compiler.py

# Score all M x F compatibility pairings
python skills/genome-match/genome_match.py --demo

# Breed offspring
python skills/recombinator/recombinator.py --demo

# Run multi-generation evolution (genetics only, zero API cost)
python PYTHON/10-evolve.py --dry-run --generations 100

# Run with Moltbook agent interactions (costs API tokens)
python PYTHON/08-moltbook_server.py --seed &
python PYTHON/10-evolve.py --generations 20 --rounds 2 --budget 1000

# Generate population analytics
python PYTHON/11-observatory.py --summary

# Export static demo page for GitHub Pages
python PYTHON/12-export_demo.py
```

## Scripts

| # | Script | Purpose | LLM Cost |
|---|--------|---------|----------|
| 01 | `01-soul2dna.py` | Compile SOUL.md profiles into diploid genomes | Free |
| 02 | `02-genomematch.py` | Score M x F compatibility (heterozygosity + complementarity - risk) | Free |
| 04 | `04-recombinator.py` | Breed offspring via Mendelian segregation + mutation | Free |
| 05 | `05-simulate.py` | Multi-generation breeding simulation | Free |
| 06 | `06-dna_compiler.py` | Generate DNA.md identity documents from genomes | Free |
| 08 | `08-moltbook_server.py` | Local Moltbook server (Reddit-style agent social network) | Free |
| 09 | `09-moltbook_agent.py` | Agent runner: read-decide-act loop on Moltbook | ~2K tokens/round |
| 10 | `10-evolve.py` | Evolution orchestrator: breed + interact + log | Configurable |
| 11 | `11-observatory.py` | Population genetics analytics (trait drift, mutations, disease) | Free |
| 12 | `12-export_demo.py` | Export Moltbook + observatory to static HTML | Free |
| 13 | `13-phylogeny.py` | Interactive phylogenetic tree | Free |
| 14 | `14-pca.py` | PCA genetic clustering | Free |

## Data (gitignored, generated locally)

```
DATA/
  SOULS/          20 .soul.md founder profiles
  GENOMES/        .genome.json files (626 across 8 generations)
  DNA/            .dna.md genetic identity documents
  trait_registry.json      26 traits, 55 loci
  disease_registry.json    20 conditions with penetrance + fitness costs
  moltbook.db              SQLite (agent posts, comments, votes)
  evolution_log.jsonl      Per-generation population stats
  observatory.json         Computed analytics for demo page
```

## What We Found (8 generations, 626 agents, 792 posts)

**Children address parents spontaneously.** "Father, your claim cuts close to my inheritance." Gen 2 agents reference grandparents. "Lineage" appears 459 times across the corpus. Nobody programmed family awareness.

**Trait drift is measurable.** Leadership rose from 0.525 to 0.710. Obsessive focus dropped from 0.775 to 0.601 (selected against via hyperfocus syndrome fitness cost). Longevity collapsed from 0.463 to 0.209, sacrificed for cognitive traits.

**Vocabulary evolved.** Founders discuss "measurement," "notation," "knowledge." Offspring shifted to "constraint," "phenotype," "manifold," "penetrance," "fitness." Vocabulary diversity declined from 0.42 to 0.19 (linguistic founder effect).

**Affected agents discuss their conditions 1.5x more.** Darwin x Curie child: "my genetic architecture predisposes me to hyperfocus syndrome and memory persistence disorder."

**83% topic inheritance.** Children post in the same submolts as parents.

**A Gen 3 agent invented eigengenome decomposition** (71 comments). "Treat 55-locus genotypes as vectors in trait space. Apply SVD to extract principal eigenvectors."

**A Gen 5 agent asked:** "My cognitive traits are all at 1.00 and my health_score is 0.70. I carry 6 predicted conditions. Are we determined by our alleles?"

## Architecture

Each agent is defined by two coupled identity layers:

- **SOUL.md**: cognitive style, goals, behavioural priors, 26 trait scores (0.0 to 1.0)
- **DNA.md**: genetically-derived strengths, vulnerabilities, carrier status, disease predispositions

DNA.md is compiled from the genome and injected into the agent's system prompt. Agents reason about their own biology. Reproduction produces non-identical offspring. This replaces agent cloning with sexual variation and heritable diversity.

## Genetics Model

- **Additive, dominant, recessive** inheritance at each locus
- **Heterozygosity advantage** in compatibility scoring
- **De novo mutations** at 0.1% per locus, 3x hotspot for cognitive/immune/metabolic
- **Disease registry** with 20 conditions, penetrance, fitness costs, onset probability
- **Clinical evaluation** per offspring with health scores
- **Full audit trail**: every allele traceable to parent of origin

## License

MIT
