# ClawBio

<p align="center">
  <img src="img/clawbio-logo.jpeg" alt="ClawBio Logo" width="300">
</p>

**The first bioinformatics-native AI agent skill library.** Built on [OpenClaw](https://github.com/openclaw/openclaw) (180k+ GitHub stars), ClawBio gives AI agents domain-specific skills for genomics and computational biology — all running locally on your machine.

## Core Principles

- **Local-first**: No genomic data leaves your laptop. No cloud uploads, no data exfiltration.
- **Reproducible**: Every analysis exports `commands.sh`, `environment.yml`, and SHA-256 checksums — anyone can reproduce it without the agent.
- **Modular**: Each skill is a self-contained directory (`SKILL.md` + Python scripts) that plugs into the orchestrator.

## Why ClawBio?

General-purpose AI agents are powerful but blind to the specific needs of biological research:

- **Privacy**: Genomic data is sensitive. Cloud-first agents risk exposing patient samples, proprietary variants, and unpublished findings.
- **Reproducibility**: Biology demands audit trails. Every analysis step must be logged, versioned, and exportable as a reproducible pipeline.
- **Domain expertise**: A generic agent does not know that a VCF file needs ancestry-aware annotation, or that single-cell data requires doublet removal before clustering.

ClawBio fills this gap with skills that understand biology from the ground up.

## Current Functionality (MVP)

### Bio Orchestrator

A meta-agent that routes natural-language bioinformatics requests to the right specialist skill. Detects file types (`.vcf`, `.fastq`, `.h5ad`) and keywords, chains multi-step analyses, and assembles final reports.

### Equity Scorer

Computes the **HEIM (Health Equity Index for Minorities)** score from VCF or ancestry data:

- Observed & expected **heterozygosity** per population
- Pairwise **FST** (Nei's GST method) between all populations
- **PCA** on the genotype matrix
- **HEIM Equity Score** (0-100) quantifying how well a dataset represents global population diversity
- Outputs: markdown report, 5 publication-quality figures (PCA plot, FST heatmap, heterozygosity comparison, ancestry distribution, HEIM gauge), CSV tables, reproducibility bundle

**Demo result**: 50 samples across 5 populations (AFR, AMR, EAS, EUR, SAS), 500 variants → HEIM Score **76.2/100** (Good). Identified EUR overrepresentation at 44% vs 16% global proportion.

## Skills

| Skill | Status | Description |
|-------|--------|-------------|
| [Bio Orchestrator](skills/bio-orchestrator/) | MVP | Meta-agent that routes bioinformatics requests to specialised sub-agents |
| [Equity Scorer](skills/equity-scorer/) | MVP | HEIM diversity metrics from VCF/ancestry data; heterozygosity, FST, PCA, equity reports |
| [VCF Annotator](skills/vcf-annotator/) | Planned | Variant annotation with VEP, ancestry context, and diversity metrics |
| [Lit Synthesizer](skills/lit-synthesizer/) | Planned | PubMed/bioRxiv search with LLM summarisation and citation graphs |
| [scRNA Orchestrator](skills/scrna-orchestrator/) | Planned | Seurat/Scanpy automation: QC, clustering, DE analysis, visualisation |
| [Struct Predictor](skills/struct-predictor/) | Planned | AlphaFold/Boltz/Chai wrappers for local structure prediction |
| [Seq Wrangler](skills/seq-wrangler/) | Planned | FastQC, alignment, BAM processing, QC reporting |
| [Repro Enforcer](skills/repro-enforcer/) | Planned | Export any analysis as Conda env + Singularity container + Nextflow pipeline |

## Quick Start

### Prerequisites

- [OpenClaw](https://github.com/openclaw/openclaw) installed and configured
- Python 3.11+
- Bioinformatics tools for your skill of choice (see individual SKILL.md files)

### Install a skill

```bash
# Install the Bio Orchestrator (routes to sub-skills automatically)
openclaw install skills/bio-orchestrator

# Install the Equity Scorer for diversity analysis
openclaw install skills/equity-scorer
```

### Use a skill

```bash
# Ask the orchestrator to analyse a VCF file
openclaw "Analyse the diversity metrics in my VCF file at data/samples.vcf"

# Run equity scoring directly
openclaw "Score the population diversity in data/ancestry.csv using HEIM metrics"
```

## Architecture

```
User Request
    |
    v
Bio Orchestrator (routing + file I/O + reporting)
    |
    +---> Equity Scorer (diversity metrics, HEIM index)
    +---> VCF Annotator (variant annotation, VEP)
    +---> Lit Synthesizer (literature search, summarisation)
    +---> scRNA Orchestrator (single-cell pipelines)
    +---> Struct Predictor (protein structure)
    +---> Seq Wrangler (sequence QC, alignment)
    +---> Repro Enforcer (reproducibility export)
    |
    v
Markdown Report + Audit Log + Reproducibility Bundle
```

Each skill is a standalone SKILL.md + supporting scripts. The Bio Orchestrator routes requests to the right skill based on input type and user intent, but every skill also works independently.

See [docs/architecture.md](docs/architecture.md) for the full design.

## Future Skills (Roadmap)

| Skill | What it does | Target |
|-------|-------------|--------|
| **VCF Annotator** | Variant annotation with VEP, ClinVar, gnomAD + ancestry context | Mar 2026 |
| **Lit Synthesizer** | PubMed/bioRxiv search with LLM summarisation and citation graphs | Mar 2026 |
| **scRNA Orchestrator** | Seurat/Scanpy automation: QC, clustering, DE analysis, visualisation | Mar 2026 |
| **Struct Predictor** | Local AlphaFold/Boltz/Chai wrappers for protein structure prediction | Apr 2026 |
| **Seq Wrangler** | FastQC, alignment, BAM processing, QC reporting | Apr 2026 |
| **Repro Enforcer** | Export any analysis as Conda env + Singularity container + Nextflow pipeline | Apr 2026 |

## Community Wanted Skills

We want skills from the bioinformatics community. If you work with genomics, proteomics, metabolomics, imaging, or clinical data, you can contribute a skill.

See [CONTRIBUTING.md](CONTRIBUTING.md) for the submission process and [templates/SKILL-TEMPLATE.md](templates/SKILL-TEMPLATE.md) for the skill skeleton.

- **GWAS Pipeline** — PLINK/REGENIE automation
- **Metagenomics Classifier** — Kraken2/MetaPhlAn wrapper
- **Clinical Variant Reporter** — ACMG classification
- **Pathway Enricher** — GO/KEGG enrichment analysis
- **Phylogenetics Builder** — IQ-TREE/RAxML automation
- **Proteomics Analyser** — MaxQuant/DIA-NN
- **Spatial Transcriptomics** — Visium/MERFISH

## License

MIT — clone it, run it, build a skill, submit a PR.

## Citation

If you use ClawBio in your research, please cite:

```bibtex
@software{clawbio_2026,
  author = {Corpas, Manuel},
  title = {ClawBio: An Open-Source Library of AI Agent Skills for Reproducible Bioinformatics},
  year = {2026},
  url = {https://github.com/manuelcorpas/ClawBio}
}
```

## Links

- **Demo slides**: [manuelcorpas.github.io/ClawBio/slides/](https://manuelcorpas.github.io/ClawBio/slides/)
- [OpenClaw](https://github.com/openclaw/openclaw) — The agent platform
- [ClawHub](https://clawhub.ai) — Skill registry
- [HEIM Index](https://heim-index.org) — Health Equity Index for Minorities
- [Corpus Core](https://github.com/manuelcorpas/corpus-core) — RAG memory system for research
