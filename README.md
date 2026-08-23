<h3 align="center">ClawBio</h3>

<p align="center">
  <strong>The first bioinformatics-native AI agent skill library.</strong><br>
  Built on <a href="https://github.com/openclaw/openclaw">OpenClaw</a> (180k+ GitHub stars). Local-first. Privacy-focused. Reproducible.
</p>

<p align="center">
  <a href="https://github.com/ClawBio/ClawBio/actions/workflows/ci.yml"><img src="https://github.com/ClawBio/ClawBio/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="#quick-start"><img src="https://img.shields.io/badge/python-3.11+-blue?logo=python&logoColor=white" alt="Python 3.11+"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="MIT License"></a>
  <a href="https://clawhub.ai"><img src="https://img.shields.io/badge/ClawHub-97_skills-orange" alt="ClawHub Skills"></a>
  <a href="https://doi.org/10.5281/zenodo.19420648"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.19420648.svg" alt="DOI"></a>
  <a href="https://github.com/ClawBio/ClawBio/issues"><img src="https://img.shields.io/github/issues/ClawBio/ClawBio" alt="Open Issues"></a>
  <a href="https://clawbio.github.io/ClawBio/slides/"><img src="https://img.shields.io/badge/slides-London_Bioinformatics_Meetup-purple" alt="Slides"></a>
</p>

---

### Install in 30 seconds

```bash
pip install clawbio                 # Python 3.11+
clawbio run pharmgx --demo
```

Prefer [conda](https://docs.conda.io/)? `conda install -c bioconda clawbio`.

Or use as a Python library:

```python
from clawbio import run_skill, list_skills
result = run_skill("pharmgx", demo=True)
```

Or install as a [Claude Code](https://claude.ai/claude-code) plugin: `/plugin marketplace add ClawBio/ClawBio`

Using **Cursor, Zed, VS Code or Claude Desktop**? ClawBio ships an [MCP](https://modelcontextprotocol.io)
server. No install needed:

```json
{
  "mcpServers": {
    "clawbio": {
      "command": "uvx",
      "args": ["--from", "clawbio[mcp]", "clawbio", "mcp"]
    }
  }
}
```

It runs locally over stdio and defaults to demo data only, so connecting it cannot by itself
give an agent access to your genome. Full setup: [docs.clawbio.ai/reference/mcp](https://docs.clawbio.ai/reference/mcp/).

**Developing ClawBio or want all skills with full demo data?** Work from a
source checkout instead ([uv](https://docs.astral.sh/uv/) recommended):

```bash
git clone https://github.com/ClawBio/ClawBio.git
cd ClawBio
uv sync                            # installs from pyproject.toml + uv.lock
uv run python clawbio.py run pharmgx --demo
```

---

<p align="center">
  <img src="img/clawbio-demo.gif" alt="ClawBio GWAS Lookup demo — querying 9 genomic databases from the terminal" width="700">
</p>

---

## What ClawBio Does Today

**97 skills (91 with runnable demo data) + 8,182 Galaxy tools + 4,539 tests + benchmark validation. Local-first by default. Reproducible. No guessing.**
> **v0.5.0 released** (4 Apr 2026): Validation and Benchmark Infrastructure. AD ground truth benchmark, mock API server for offline testing, swappable fine-mapping pipeline (SuSiE vs ABF), 74 benchmark tests, red/green TDD mandate. [Release notes](https://github.com/ClawBio/ClawBio/releases/tag/v0.5.0). DOI: [10.5281/zenodo.19420648](https://doi.org/10.5281/zenodo.19420648).

Snap a photo of a medication in Telegram. ClawBio identifies the drug from the packaging, queries your pharmacogenomic profile from [your own genome](docs/demo-genome.md), and returns a personalised dosage card — on your machine, in seconds:

<p align="center">
  <img src="skills/drug-photo/demo-images/00-00-warfarin.jpg" width="360" alt="Warfarin 2mg medication packaging — ClawBio identifies the drug from a photo and returns a personalised pharmacogenomic report">
</p>

> **Warfarin** | CYP2C9 \*1/\*2 Intermediate · VKORC1 High Sensitivity
> **AVOID — DO NOT USE** · Standard dose causes over-anticoagulation in this genotype.

Or take any genetic variant (identified by its rsID — a unique label like [rs9923231](https://www.ncbi.nlm.nih.gov/snp/rs9923231)) and search nine genomic databases at once to find every known disease association, tissue-specific effect, and population frequency. Or estimate your genetic predisposition to conditions like type 2 diabetes by combining thousands of small-effect variants into a single polygenic risk score. Or explore the [UK Biobank](https://www.ukbiobank.ac.uk/) — a half-million-person research dataset — by asking in plain English what fields measure blood pressure, grip strength, or depression, and get back the exact field IDs, descriptions, and linked publications you need.

Many ClawBio analyses write a `reproducibility/` bundle with replay commands, environment metadata, and output checksums. The exact files can vary by skill, and some replays still require the original external inputs to be present. See [docs/reproducibility.md](docs/reproducibility.md).

---

## Reference Genome

ClawBio's demo data is built on a real, fully open human genome: the **Corpasome**. The [23andMe SNP chip](docs/demo-genome.md) (~600K variants) has been available since launch. Now, the project also ships subsets from a **30x Illumina whole-genome sequence** (GRCh37), covering ~4M SNPs, ~600K indels, and structural variants (DEL, DUP, INV, BND, INS, CNVs). All data comes from a single individual (Manuel Corpas), licensed CC0, and published on Zenodo ([doi:10.5281/zenodo.19297389](https://doi.org/10.5281/zenodo.19297389)). This dataset is provided for research and educational purposes only.

See [docs/reference-genome.md](docs/reference-genome.md) for use cases, subsets, and citation details.

---

## Validation & Benchmarking

ClawBio v0.5.0 introduces the platform's first systematic validation infrastructure. Skills are scored against curated ground truth with objective metrics.

**What's in the benchmark suite:**

| Component | Description |
|-----------|-------------|
| **AD Ground Truth** | 34 Alzheimer's disease genes across 3 evidence tiers (Mendelian, GWAS-replicated, novel), 20 negative controls, 10 lead variants |
| **Mock API Server** | Deterministic endpoints for Ensembl, GWAS Catalog, ClinPGx. Offline CI without rate limits |
| **Benchmark Scorer** | Gene recovery rate, FDR, precision, recall, F1, tier-weighted composite score |
| **Swappable Fine-Mapping** | ABF vs SuSiE head-to-head on same data. Method registry pattern for adding FINEMAP, PolyFun |
| **Nightly Sweep** | Runs every skill in demo mode, scores outputs against ground truth |

```bash
# Run the fine-mapping benchmark
python tests/benchmark/finemapping_benchmark.py --output /tmp/fm_bench

# Score a gene list against AD ground truth
python tests/benchmark/benchmark_scorer.py --genes "APP,BIN1,CLU,TREM2,GAPDH"

# Start mock API server for offline testing
python tests/benchmark/mock_api_server.py &
```

**74 benchmark tests** at v0.5.0 baseline, all green. The public leaderboard now tracks **161 / 182 tests passing (88.5%)** across 10 audited skills, up from 80 / 140 (57.1%) at the original audit. See [benchmarks.html](https://clawbio.ai/benchmarks.html) for the live leaderboard and [CHANGELOG.md](CHANGELOG.md) for full details.

---

## The Problem

You read a paper. You want to reproduce Figure 3. So you:

1. Go to GitHub. Clone the repo.
2. Wrong Python version. Fix dependencies.
3. Need the reference data — where is it?
4. Download 2GB from Zenodo. **Link is dead.**
5. Email the first author. **Wait 3 weeks.**
6. Paths are hardcoded to `/home/jsmith/data/`.
7. Two days later: still broken. **You give up.**

Now imagine the same paper published a **skill**:

```bash
python skills/claw-ancestry-pca/ancestry_pca.py --demo --output fig3
# Example demo output regenerated. Verify bundled checksums where provided.
```

**That's ClawBio.** The goal is to make replayable bioinformatics workflows straightforward when a skill ships demo data and helper-backed reproducibility outputs.

---

## What Is ClawBio?

Current agentic bioinformatics systems address either the **reasoning layer** (constraining LLM outputs with knowledge graphs or fine-tuning) or the **connectivity layer** (wrapping tools as MCP servers). Neither addresses the **specification layer**: the encoding of a domain expert's analytical decisions into a machine-readable contract that constrains agent behaviour. Without this layer, the agent must reconstruct expert knowledge from its training distribution, a stochastic, unversioned process.

A **skill** is a self-contained directory comprising a declarative specification (`SKILL.md`), validated Python code, demo data, and reproducibility support. In many skills, that includes a `reproducibility/` bundle with `commands.sh`, `environment.yml`, and SHA-256 checksums, sometimes with extra lock metadata such as `runtime-lock.json`. The specification is a contract, not a prompt: it encodes the domain expert's analytical decisions so the LLM orchestrates but does not improvise.

```
Ad-hoc LLM code generation  = stochastic, unversioned, unverifiable
ClawBio skill                = specification-constrained, versioned, reproducible
```

- **Specification-first**: Domain expertise resides in `SKILL.md`, not in model weights. Specifications are versioned, human-readable, peer-reviewable, and trivially updatable.
- **Agent-agnostic**: Skills execute identically whether invoked by Claude, ChatGPT, or a locally hosted model via Ollama. Reproducibility is decoupled from any specific AI vendor.
- **Local-first by default**: Your genome is analysed on your machine. Skills that send data to external services (hosted inference, public annotation APIs) are individually labelled and run only when you explicitly invoke them.
- **Reproducible**: Many skills export replay metadata such as `commands.sh`, `environment.yml`, and SHA-256 checksums so runs can be rechecked without the original agent session.
- **MIT licensed**: Open-source, free, community-driven.

## Why Not Just Use an LLM?

Ask any LLM to "profile the pharmacogenes in my VCF file." It will produce plausible Python. But the code may use outdated CPIC guidelines, hallucinate star allele functional classifications, or confuse reduced-function and no-function alleles. CYP2D6\*4 is a no-function allele; misclassifying it as reduced-function determines whether a patient receives a dose adjustment or is told to avoid a drug entirely. Approximately 7% of individuals of European ancestry are CYP2D6 poor metabolisers for whom codeine provides zero analgesic effect, and roughly 0.5% carry DPYD variants where standard fluorouracil dosing can be lethal.

An LLM should not be improvising these from training data. ClawBio encodes the correct bioinformatics decisions in versioned, peer-reviewable specifications so the agent executes them correctly every time.

---

## Provenance and Reproducibility

ClawBio's current reproducibility contract centers on a **`reproducibility/` bundle** written inside the output directory for skills that use the shared reproducibility helpers:

```
report/
├── report.md              # Full analysis with figures and tables
├── figures/               # Publication-quality PNGs
├── tables/                # CSV data tables
├── reproducibility/
│   ├── commands.sh        # Replay command for this run
│   ├── environment.yml    # Suggested Conda environment snapshot
│   ├── checksums.sha256   # SHA-256 for selected output files
│   └── runtime-lock.json  # Optional extra lock metadata when supported
```

The exact contents can vary by skill, and some replays also require the original external inputs or tools to be available locally. For a user-facing walkthrough, including a truthful `multiqc-reporter` example based on direct script invocation, see [docs/reproducibility.md](docs/reproducibility.md).

---

## Skills

A curated cross-section of ClawBio's 97 skills. The full machine-readable catalog (with status flags, objective `maturity_tier` evidence, trigger keywords, demo commands, and chaining partners) lives in [`skills/catalog.json`](skills/catalog.json); browse the directory at [`skills/`](skills/) to see every skill folder.

Catalog maturity tiers are computed from repository evidence: `spec-only` (SKILL.md only), `scripted` (has runnable code), `tested` (has skill tests), `cli-registered` (available via `python clawbio.py run`), `ci-validated` (explicitly tested in CI), and `bench-validated` (reserved for skills with blocking scientific benchmark validation).

| Skill | Scale | Description |
|-------|-------|-------------|
| [Bio Orchestrator](skills/bio-orchestrator/) | Infrastructure | Routes free-text or file inputs to the right skill automatically |
| [PharmGx Reporter](skills/pharmgx-reporter/) | Personal | 12 genes, 51 drugs, CPIC guidelines from consumer genetic data |
| [Drug Photo](skills/drug-photo/) | Personal | Snap a medication photo, get a personalised dosage card from your genotype |
| [GWAS Lookup](skills/gwas-lookup/) | Population | Federated variant query across 9 genomic databases (gnomAD, ClinVar, Open Targets, GTEx, LDlink, ...) |
| [GWAS PRS](skills/gwas-prs/) | Population | Polygenic risk scores from the PGS Catalog for 6+ traits |
| [Ancestry PCA](skills/claw-ancestry-pca/) | Population | PCA vs SGDP (345 samples, 164 populations) with confidence ellipses |
| [Fine-Mapping](skills/fine-mapping/) | Population | SuSiE / ABF credible sets with posterior inclusion probabilities from GWAS summary stats |
| [UKB Navigator](skills/ukb-navigator/) | Research | Semantic search across the UK Biobank schema (22,000+ fields) |
| [Galaxy Bridge](skills/galaxy-bridge/) | Research | Search, run, and chain 8,000+ Galaxy bioinformatics tools |
| [Variant Annotation](skills/variant-annotation/) | Clinical | Annotate VCF variants with Ensembl VEP REST, ClinVar significance, gnomAD frequencies |
| [Clinical Variant Reporter](skills/clinical-variant-reporter/) | Clinical | ACMG-guided clinical variant classification from VCF with GiAB validation |
| [nf-core scRNA Wrapper](skills/nfcore-scrnaseq-wrapper/) | Single-cell | Upstream FASTQ → h5ad preprocessing via nf-core/scrnaseq (simpleaf, STARsolo, kallisto, CellRanger) with strict preflight and reproducibility bundle |
| [nf-core RNA-seq Wrapper](skills/nfcore-rnaseq-wrapper/) | Bulk RNA-seq | Upstream FASTQ/BAM → count matrices via nf-core/rnaseq with strict preflight and reproducibility bundle |
| [nf-core Sarek Wrapper](skills/nfcore-sarek-wrapper/) | Variant calling | Upstream FASTQ/BAM/CRAM → VCF via nf-core/sarek (germline, tumor-only, somatic; Strelka/Mutect2/HaplotypeCaller, VEP/snpEff) with strict preflight and reproducibility bundle |
| [scRNA Orchestrator](skills/scrna-orchestrator/) | Single-cell | Scanpy automation: QC, optional doublet detection, clustering, markers, annotation |
| [Equity Scorer](skills/equity-scorer/) | Systemic | HEIM diversity metrics from VCF or ancestry CSV |
| [DnaSP](skills/dnasp/) | Population *(community)* | Python reimplementation of DnaSP 6: 16 population-genetics analyses (Pi, Tajima's D, Fst, Ka/Ks, McDonald-Kreitman) |

For the complete list including pharmacogenomics extensions, single-cell tooling, proteomics, methylation, metagenomics, structure prediction, clinical reporting, and meta-skills, see [`skills/catalog.json`](skills/catalog.json) or run `python clawbio.py list`.

### Hosted inference — Genomic Intelligence

Six skills (`gi-promoter`, `gi-splice`, `gi-enhancer`, `gi-chromatin`, `gi-expression`, `gi-annotation`) wrap [Genomic Intelligence](https://genomicintelligence.ai) — a DNA language-model platform providing transformer-based inference for promoter, splice-site, enhancer, chromatin-accessibility, expression, and annotation prediction. Unlike most ClawBio skills these call a hosted API (`api.genomicintelligence.ai`); a shared hackathon key ships in `.env.example` (opt-in only; caps are per-key — read the `RateLimit-*` response headers for the live allowance).

You can also run the same models interactively at <https://genomicintelligence.ai> — a browser UI with gene search and IGV/JBrowse viewers, no setup required. API reference: <https://docs.genomicintelligence.ai>. Questions, errors, or production-tier keys: **contact@genomicintelligence.ai**.

### Contributing a Skill

Wrap your bioinformatics pipeline as a skill and submit a PR. Several community-contributed skills are now in production (NutriGx Advisor, analyze-fasta, WGS-PRS, DnaSP, ClawPathy Autoresearch, and others). See [skills/catalog.json](skills/catalog.json) for the full list with status flags.

```bash
cp templates/SKILL-TEMPLATE.md skills/<your-skill-name>/SKILL.md
# Edit SKILL.md, add Python implementation, demo data, and tests
```

See [CONTRIBUTING.md](CONTRIBUTING.md) for the full submission process. Join the contributors community on Telegram: [t.me/ClawBioContributors](https://t.me/ClawBioContributors).

---

## Genomebook

**Genomebook** is a synthetic-genetics sandbox built into ClawBio. It turns fictional or historical character profiles ("souls") into diploid genomes, scores compatibility, and breeds offspring across generations, complete with Mendelian inheritance, de novo mutations, and clinical evaluation.

### Pipeline

```
SOUL.md  -->  Soul2DNA  -->  .genome.json  -->  GenomeMatch  -->  Recombinator  -->  Gen-N offspring
 (trait       (compiler)     (diploid loci)     (M x F rank)     (meiosis +        (.genome.json
  scores)                                                         mutation)          with clinical
                                                                                     history)
```

### Data

- **20 souls** in `GENOMEBOOK/DATA/SOULS/` (Einstein, Curie, Turing, Hypatia, Da Vinci, ...)
- **20 generation-0 genomes** in `GENOMEBOOK/DATA/GENOMES/`
- **26 traits** across 55 loci (trait_registry.json)
- **Disease registry** with penetrance, fitness costs, and onset probabilities

### Quick Start

```bash
# Compile all souls to genomes
python skills/soul2dna/soul2dna.py --demo

# Score all M x F compatibility pairings
python skills/genome-match/genome_match.py --demo

# Breed Einstein x Anning (3 offspring)
python skills/recombinator/recombinator.py --demo
```

### What It Models

- **Additive, dominant, recessive** inheritance at each locus
- **Heterozygosity advantage** (heterosis) in compatibility scoring
- **Carrier risk** flagging for autosomal recessive conditions
- **De novo mutations** with hotspot categories (cognitive, immune, metabolic)
- **Clinical evaluation** with penetrance, onset probability, and fitness costs
- **Health scores** derived from cumulative condition burden

---

## Skills in Detail

### PharmGx Reporter — *Personal Scale*

Generates a pharmacogenomic report from consumer genetic data (23andMe, AncestryDNA):

- Parses raw genetic data (auto-detects format, including gzip)
- Extracts **31 pharmacogenomic SNPs** across **12 genes** (CYP2C19, CYP2D6, CYP2C9, VKORC1, SLCO1B1, DPYD, TPMT, UGT1A1, CYP3A5, CYP2B6, NUDT15, CYP1A2)
- Calls star alleles and determines metabolizer phenotypes
- Looks up **CPIC drug recommendations** for **51 medications**
- Zero dependencies. Runs in **< 1 second**.

```bash
python pharmgx_reporter.py --input demo_patient.txt --output report
```

**Demo result**: CYP2D6 \*4/\*4 (Poor Metabolizer) → **10 drugs AVOID** (codeine, tramadol, 7 TCAs, tamoxifen), 20 caution, 21 standard.

> ~7% of people are CYP2D6 Poor Metabolizers — codeine gives them zero pain relief. ~0.5% carry DPYD variants where standard 5-FU dose can be lethal. This skill catches both.

### Drug Photo — *Personal Scale*

Snap a photo of any medication in Telegram. ClawBio identifies the drug from the packaging and returns a personalised dosage card against your own genotype.

- Claude vision extracts drug name and visible dose from the photo
- Cross-references your 23andMe genotype against 31 PGx SNPs
- Four-tier classification: **STANDARD DOSING / USE WITH CAUTION / AVOID / INSUFFICIENT DATA**
- Correct VKORC1 complement-strand handling (23andMe reports minus strand for rs9923231)
- Works for warfarin, clopidogrel, codeine, simvastatin, tamoxifen, sertraline, and 20+ others

```bash
python pharmgx_reporter.py --drug warfarin --dose "5mg" --input my_23andme.txt --output report
```

> No command needed in Telegram — send any medication photo and RoboTerri triggers the skill automatically.

### GWAS Lookup — *Population Scale*

Federated variant query across nine genomic databases in a single command:

| Database | What you get |
|----------|-------------|
| GWAS Catalog | Genome-wide significant associations |
| gnomAD | Allele frequencies across 125,748 exomes |
| ClinVar | Clinical significance and condition links |
| Open Targets | Disease-gene evidence scores |
| Ensembl | Functional annotation, regulatory impact |
| GTEx | eQTL data, tissue-specific expression effects |
| LDlink | Linkage disequilibrium across 26 populations |
| UK Biobank PheWAS | Phenome-wide associations across 4,000+ traits |
| LOVD | Variant pathogenicity database |

```bash
python gwas_lookup.py --rsid rs3798220 --output report
python gwas_lookup.py --demo --output /tmp/gwas_lookup_demo
```

### UKB Navigator — *Research Scale*

Semantic search across the UK Biobank schema. Ask in plain English what UK Biobank measures about any phenotype — get field IDs, descriptions, data types, participant counts, and linked publications back instantly.

```bash
python ukb_navigator.py --query "grip strength"   --output report
python ukb_navigator.py --field 21001              --output report   # BMI
python ukb_navigator.py --demo                     --output /tmp/ukb_demo
```

Built on a ChromaDB embedding of the full UKB Data Showcase (22,000+ fields).

### Ancestry PCA — *Population Scale*

Runs principal component analysis on your cohort against the SGDP reference panel (345 samples, 164 global populations):

- Contig normalisation (chr1 vs 1)
- IBD removal (related individuals filtered)
- Common biallelic SNPs only
- Confidence ellipses per population
- Publication-quality **4-panel figure** generated instantly

```bash
python ancestry_pca.py --demo --output ancestry_report
```

**Demo result**: 736 Peruvian samples across 28 indigenous populations. Amazonian groups (Matzes, Awajun, Candoshi) sit in genetic space that no SGDP population occupies — genuinely underrepresented, not just in GWAS, but in the reference panels themselves.

### Semantic Similarity Index — *Systemic Scale*

Computes a Semantic Isolation Index for diseases using 13.1M PubMed abstracts and PubMedBERT embeddings (768-dim):

- **SII** (Semantic Isolation Index): higher = more isolated in literature
- **KTP** (Knowledge Transfer Potential): higher = more cross-disease spillover
- **RCC** (Research Clustering Coefficient): diversity of research approaches
- **Temporal Drift**: how research focus evolves over time
- Publication-quality **4-panel figure**

**Status: not yet implemented.** `semantic_sim.py` and the pipeline scripts behind it are not
in this repository yet (see `skills/claw-semantic-sim/SKILL.md`); there is no runnable demo.

**Key finding**: Neglected tropical diseases are **+38% more semantically isolated** (P < 0.0001, Cohen's d = 0.84). 14 of the 25 most isolated diseases are Global South priority conditions. Knowledge silos kill innovation — a malaria immunology breakthrough could help leishmaniasis, but the literatures don't talk to each other.

> Corpas et al. (2026). *HEIM: Health Equity Index for Measuring structural bias in biomedical research.* Under review.

---

## Quick Start

ClawBio is first and foremost a local bioinformatics skill library. You can execute it directly from the UNIX command line or through its Python module without any messenger or chat interface at all.

The OpenClaw-borrowed "magic" happens when an external LLM interprets a user's request and translates it into those local library or command-line calls. In other words, the LLM operates at the orchestration layer, while the biological data processing remains local, inspectable, and reproducible.

To experience that orchestration layer conversationally, ClawBio can be used through Telegram or Discord via RoboTerri, as a Claude Code skill, or through a self-hosted OpenClaw gateway with browser-based webchat. Of these options, the self-hosted OpenClaw gateway offers the strongest privacy story together with LLM-assisted interaction, because the LLM stays at the meta/orchestration level rather than operating on the underlying biological data themselves.

### Install and run

```bash
pip install clawbio                      # or: conda install -c bioconda clawbio
clawbio run pharmgx --demo
```

PharmGx demo runs in <2 seconds. Needs Python 3.11+. To develop ClawBio or get
all skills with full demo data, work from a source checkout instead:

```bash
git clone https://github.com/ClawBio/ClawBio.git && cd ClawBio
uv sync                                  # or: pip install -e .
uv run python clawbio.py run pharmgx --demo
```

### Use as a Python library

```python
from clawbio import run_skill, list_skills

# List available skills
skills = list_skills()

# Run pharmacogenomics with demo data
result = run_skill("pharmgx", demo=True)

# Run with your own input
result = run_skill("gwas-lookup", rsid="rs3798220", output="/tmp/report")
```

### Install as a Claude Code plugin

Inside [Claude Code](https://claude.ai/claude-code):

```
/plugin marketplace add ClawBio/ClawBio
/plugin install clawbio
```

All skills are then available as agent-routable commands. Alternatively, clone the repo and open it as your working directory in Claude Code; the `CLAUDE.md` at the repo root teaches Claude how to route requests to skills automatically.

### Slash Commands

The repository also ships reusable slash commands in [`commands/`](commands/) for Claude Code and compatible agents:

| Command | Purpose |
|---------|---------|
| `/analyse` | Analyse a file or input with the appropriate ClawBio skill |
| `/new-skill` | Scaffold a new skill from the official template |
| `/list-skills` | List available skills from `skills/catalog.json` |
| `/run-demo` | Run a skill demo with built-in sample data |

### Try all skills

```bash
python clawbio.py list                           # See available skills
python clawbio.py run pharmgx --demo             # Pharmacogenomics (1s)
python clawbio.py run equity --demo              # Equity scoring (55s)
python clawbio.py run nutrigx --demo             # Nutrigenomics (60s)
python clawbio.py run metagenomics --demo        # Metagenomics (3s)
python clawbio.py run scrna --demo               # scRNA clustering + marker detection (PBMC3k-first demo)
python clawbio.py run scrna --demo --doublet-method scrublet
                                                 # Optional doublet detection before clustering
python clawbio.py run compare --demo             # Manuel Corpas vs George Church (10s)
python clawbio.py run gwas-lookup --demo         # rs3798220 across 9 databases (5s)
python clawbio.py run prs --demo                 # Polygenic risk scores (10s)
python clawbio.py run ukb-navigator --demo       # UK Biobank schema search (5s)
python clawbio.py run profile --demo             # Unified genomic profile (30s)
python clawbio.py run galaxy --demo              # Galaxy Bridge FastQC demo (offline)
python clawbio.py run rnaseq --demo              # RNA-seq DE demo (bulk/pseudo-bulk)
python clawbio.py run methylation --demo        # Epigenetic methylation clocks via PyAging
python clawbio.py run protocols-io --demo       # Protocols.io protocol search
python clawbio.py run variant-annotation --demo # VCF variant annotation (VEP + ClinVar + gnomAD)
python clawbio.py run proteomics --demo         # Proteomics differential expression
python clawbio.py run clinical-trial --demo     # Clinical trial finder
```

### Run with your own data

```bash
python clawbio.py run pharmgx --input my_23andme.txt --output results/
python clawbio.py run rnaseq --input counts.csv,metadata.csv --output results_rnaseq/
python clawbio.py run methylation --geo-id GSE139307 --output results_methylation/
```

### Run tests

```bash
uv run pytest                # or: pip install pytest && python -m pytest
```

### Dependencies

Core dependencies are declared in [`pyproject.toml`](pyproject.toml) and pinned in [`uv.lock`](uv.lock): biopython, pandas, numpy, scikit-learn, matplotlib, openai, pydeseq2, google-cloud-bigquery, google-auth, conda-lock, rocrate. Most skills run with just these.

`uv sync` installs everything in a reproducible virtual environment. To add or update a dependency, run `uv add <package>` (or edit `pyproject.toml` and re-run `uv sync`); commit the resulting `uv.lock` change.

Some skills have additional requirements:

| Skill | Extra dependency | Install |
|-------|-----------------|---------|
| Metagenomics | Kraken2, RGI, HUMAnN3 | Conda (see skill README) |
| Methylation Clock | PyAging | `pip install pyaging` |
| scRNA Embedding | scvi-tools | `pip install scvi-tools` |
| Galaxy Bridge | BioBlend | `pip install bioblend` |

No Docker or Singularity required for core functionality. Skills that need external bioinformatics tools document their setup in their own `SKILL.md`.

---

## Conversational Interfaces

### RoboTerri

<p align="center">
  <img src="img/terri_attwood_avatar_top_left.png" alt="RoboTerri" width="250">
  <br><em>RoboTerri — ClawBio's messaging-facing agent, inspired by <a href="https://en.wikipedia.org/wiki/Teresa_Attwood">Prof. Teresa K. Attwood</a></em>
</p>

**RoboTerri** is the messaging-facing interface to ClawBio. It currently runs over Telegram and Discord, letting an external LLM interpret free-text requests and translate them into local ClawBio skill executions. The public demo bot runs against a real human genome ([Manuel Corpas](https://en.wikipedia.org/wiki/Manuel_Corpas), CC0 public domain). Named after [Prof. Teresa K. Attwood](https://en.wikipedia.org/wiki/Teresa_Attwood) — a pioneer of bioinformatics education, founding Chair of GOBLET, and winner of the 2021 ISCB Outstanding Contributions Award.

<p align="center">
  <a href="https://t.me/RoboTerri_bot">
    <img src="demo/roboterri-preview.gif" alt="RoboTerri Telegram bot demo — querying a real genome" width="300">
  </a>
  <br>
  <a href="https://t.me/RoboTerri_bot"><strong>Try RoboTerri now — no install needed →</strong></a>
</p>

Ask it anything:

- **"Give me my pharmacogenomic summary"** — analyses 12 genes, 51 drugs
- **"What diseases am I at risk for?"** — polygenic risk scores for 6 conditions
- **Send a photo of any medication** — checks CYP2D6/CYP2C19 metaboliser status
- `/demo pharmgx` `/demo prs` `/demo nutrigx` `/demo compare` `/demo profile`

```
You:        [send 23andMe file]
RoboTerri:  Running PharmGx Reporter...
            CYP2D6 *4/*4 — Poor Metabolizer → 10 drugs AVOID
            [report.md attached]
            [3 figures attached]

You:        [send photo of warfarin packet]
RoboTerri:  Warfarin detected. Running Drug Photo skill...
            CYP2C9 *1/*2 · VKORC1 High Sensitivity
            AVOID — DO NOT USE at standard dose.

You:        run gwas-lookup rs3798220
RoboTerri:  Querying 9 databases...
            rs3798220 (LPA) — coronary artery disease, Lp(a) levels.
            eQTL in liver (GTEx). gnomAD MAF 0.07.
```

RoboTerri auto-detects file type (23andMe `.txt`, AncestryDNA `.csv`, VCF, FASTQ) and routes to the right skill via the Bio Orchestrator. Photos of medications trigger the Drug Photo skill automatically — no command needed.

> **[Install your own RoboTerri](docs/tutorial-roboterri-install.md)**: Set up your own Telegram or Discord adapter running ClawBio skills in ~20 minutes.

### OpenClaw Gateway

The **OpenClaw gateway** is the direct, self-hosted browser/webchat interface for ClawBio. Instead of routing through Telegram or Discord, it places a chat client directly in front of a local gateway that orchestrates ClawBio skill calls on your machine.

This is the strongest privacy-preserving conversational setup in the ClawBio stack. It avoids reliance on external messenger platforms, keeps operational control in your hands, and lets the LLM remain at the orchestration layer while the biological analysis itself stays local.

See [docs/custom-domain-webchat.md](docs/custom-domain-webchat.md) for a plain-language walkthrough of the gateway and how to expose a browser-based webchat on your own domain.

---

## Galaxy Integration

ClawBio indexes **8,000+ bioinformatics tools** from [usegalaxy.org](https://usegalaxy.org) via the Galaxy Bridge skill. Search by natural language, inspect tool schemas, and execute remotely — all from the CLI.

```bash
# Search Galaxy tools by keyword
python skills/galaxy-bridge/galaxy_bridge.py --search "metagenomics"

# Browse all 86 tool categories
python skills/galaxy-bridge/galaxy_bridge.py --list-categories

# Run a tool on Galaxy (requires GALAXY_API_KEY)
python skills/galaxy-bridge/galaxy_bridge.py --run fastqc --input reads.fq.gz --output results/

# Demo mode (offline, no API key)
python skills/galaxy-bridge/galaxy_bridge.py --demo
```

**Cross-platform chaining**: Galaxy VEP annotates variants → ClawBio PharmGx generates dosage report. Galaxy Kraken2 classifies reads → ClawBio metagenomics profiler. Neither can do this alone.

Built on [BioBlend](https://bioblend.readthedocs.io/) (Galaxy Python SDK). Developed in collaboration with the Galaxy ML SIG.

---

## Architecture

```
CLI (clawbio.py)     Python (import clawbio)     Claude Code
        │                     │                    │
        └─────────────┬───────┴──────────────┬─────┘
                      │                      │
        OpenClaw Gateway (webchat)    RoboTerri (Telegram/Discord)
                      │                      │
                      └──────────┬───────────┘
                                 │
                          ┌──────▼──────┐
                          │  Bio         │  ← routes by file type + keywords
                          │  Orchestrator│
                          └──────┬──────┘
                                 │
         ┌───────────────────────▼──────────────────────────────────────┐
         │                                                               │
         PharmGx    Equity     NutriGx    Metagenomics   Ancestry
         Reporter   Scorer     Advisor    Profiler        PCA    ...
         │                                                               │
         └───────────────────────┬──────────────────────────────────────┘
                                 │
                          ┌──────▼──────┐
                          │  Markdown    │  ← report + figures + checksums
                          │  Report      │     + reproducibility bundle
                          └─────────────┘
```

Each skill is standalone — the orchestrator routes to the right one, but every skill also works independently. The `clawbio.run_skill()` API is importable by any agent or adapter layer (RoboTerri, Claude Code, or an OpenClaw gateway), and it remains usable directly from local Python without any conversational interface at all.

See [docs/architecture.md](docs/architecture.md) for the full design.

---

## For AI Agents

ClawBio is designed to be discovered and used by AI coding agents, not just humans.

| Resource | Purpose |
|----------|---------|
| [`llms.txt`](llms.txt) | Token-optimized project summary for any LLM ([llmstxt.org](https://llmstxt.org) standard) |
| [`AGENTS.md`](AGENTS.md) | Universal guide for AI coding agents — setup, commands, style, structure, git workflow |
| [`CLAUDE.md`](CLAUDE.md) | Claude-specific routing table, CLI reference, demo commands, safety rules |
| [`skills/catalog.json`](skills/catalog.json) | Machine-readable skill index with trigger keywords, chaining partners, and demo commands |

Agents can also run `python clawbio.py list` to discover available skills programmatically.

---

## Wanted Skills

Open skill requests (PRs welcome):

| Skill | Domain |
|-------|--------|
| **claw-gwas** | PLINK/REGENIE automation (statistical genetics) |
| **claw-acmg** | Clinical variant classification (clinical genomics) |
| **claw-pathway** | GO/KEGG enrichment (functional genomics) |
| **claw-phylogenetics** | IQ-TREE/RAxML automation (evolutionary biology) |
| **claw-spatial** | Visium/MERFISH (spatial transcriptomics) |
| **claw-long-reads** | ONT/PacBio QC and assembly (long-read sequencing) |

See [Contributing a Skill](#contributing-a-skill) above for the submission process.

---

## Presentations and Demos

- **London Bioinformatics Meetup** (26 Feb 2026): project announcement. [Slides](https://clawbio.github.io/ClawBio/slides/).
- **UK AI Agent Hack, Imperial College London** (1 Mar 2026): introduced ClawBio to Peter Steinberger, creator of OpenClaw. [Video](https://www.youtube.com/watch?v=eEEA71qSOmU).
- **DoraHacks Demo Day, Imperial College London** (7 Mar 2026): live demo of pharmacogenomics, intelligent routing, multi-channel agents, and Drug Photo. [Video](https://www.youtube.com/watch?v=vxUDdjXMFwk).

---

## Versioning

ClawBio follows [Semantic Versioning](https://semver.org/). The current release is **v0.5.0**. See [CHANGELOG.md](CHANGELOG.md) for a full history of additions and breaking changes.

---

## FAQ

### What is ClawBio?

ClawBio is the **first bioinformatics-native AI agent skill library**. Built on OpenClaw (180k+ GitHub stars), it provides 96 skills (90 with runnable demo data) for genomics analysis, pharmacogenomics, ancestry profiling, and more. Local-first, privacy-focused, and reproducible.

### What are ClawBio skills?

A **skill** is a self-contained directory comprising:
- Declarative specification (`SKILL.md`)
- Validated Python code
- Demo data
- Reproducibility support (`reproducibility/` bundle with `commands.sh`, `environment.yml`, SHA-256 checksums)

Skills encode domain expert decisions into machine-readable contracts so AI agents orchestrate but do not improvise.

### What can ClawBio do?

| Category | Skills | Example Use Cases |
|----------|--------|-------------------|
| **Personal** | PharmGx, Drug Photo | Pharmacogenomic reports, medication dosage cards |
| **Population** | GWAS Lookup, Ancestry PCA | Variant queries, ancestry analysis |
| **Clinical** | Variant Annotation, Clinical Reporter | ACMG classification, VCF annotation |
| **Research** | UKB Navigator, Galaxy Bridge | UK Biobank search, 8,000+ Galaxy tools |
| **Single-cell** | scRNA Orchestrator | QC, clustering, marker detection |

### How do I install ClawBio?

```bash
git clone https://github.com/ClawBio/ClawBio.git && cd ClawBio
uv sync                            # Python 3.11+ (installs dependencies)
uv run python clawbio.py run pharmgx --demo
```

Or use pip: `pip install -e . && python clawbio.py run pharmgx --demo`

### What is RoboTerri?

**RoboTerri** is ClawBio's messaging-facing agent (Telegram/Discord). It interprets free-text requests and routes them to ClawBio skills:
- Send a 23andMe file → get pharmacogenomic report
- Send a medication photo → get dosage card from your genotype
- Query variants → get disease associations

Try it: [t.me/RoboTerri_bot](https://t.me/RoboTerri_bot)

### What is the Corpasome reference genome?

ClawBio ships demo data from a **real, fully open human genome**: the Corpasome (Manuel Corpas, CC0):
- 23andMe SNP chip (~600K variants)
- 30x Illumina WGS (~4M SNPs, ~600K indels, structural variants)
- Published on Zenodo: [doi:10.5281/zenodo.19297389](https://doi.org/10.5281/zenodo.19297389)

### How does ClawBio ensure reproducibility?

Many skills export a `reproducibility/` bundle:
- `commands.sh` — replay command
- `environment.yml` — Conda environment snapshot
- `checksums.sha256` — SHA-256 checksums for outputs
- `runtime-lock.json` — optional extra lock metadata

### What license does ClawBio use?

MIT License. Open-source, free, community-driven.

---

## Links

- **Try RoboTerri**: [t.me/RoboTerri_bot](https://t.me/RoboTerri_bot) -- query a real genome on Telegram, no install needed
- **Slides**: [clawbio.github.io/ClawBio/slides/](https://clawbio.github.io/ClawBio/slides/)
- **Tutorial**: [Install your own RoboTerri](docs/tutorial-roboterri-install.md)
- **OpenClaw**: [github.com/openclaw/openclaw](https://github.com/openclaw/openclaw) -- the agent platform
- **ClawHub**: [clawhub.ai](https://clawhub.ai) -- skill registry
- **HEIM Index**: [heim-index.org](https://heim-index.org) -- Health Equity Index for Minorities

## License

MIT. See [LICENSE](LICENSE).
