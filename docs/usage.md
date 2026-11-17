# Using ClawBio

The README covers install. This page covers everything else: the Python library, Claude Code, slash commands, every demo, your own data, tests and dependencies.

ClawBio is first and foremost a local bioinformatics skill library. You can execute it directly from the UNIX command line or through its Python module without any messenger or chat interface at all.

The OpenClaw-borrowed "magic" happens when an external LLM interprets a user's request and translates it into those local library or command-line calls. In other words, the LLM operates at the orchestration layer, while the biological data processing remains local, inspectable, and reproducible.

To experience that orchestration layer conversationally, ClawBio can be used through Telegram via RoboTerri, as a Claude Code skill, or through a self-hosted OpenClaw gateway with browser-based webchat. Of these options, the self-hosted OpenClaw gateway offers the strongest privacy story together with LLM-assisted interaction, because the LLM stays at the meta/orchestration level rather than operating on the underlying biological data themselves.

## Use as a Python library

```python
from clawbio import run_skill, list_skills

# List available skills
skills = list_skills()

# Run pharmacogenomics with demo data
result = run_skill("pharmgx", demo=True)

# Run with your own input
result = run_skill("gwas-lookup", rsid="rs3798220", output="/tmp/report")
```

## Install as a Claude Code plugin

Inside [Claude Code](https://claude.ai/claude-code):

```
/plugin marketplace add ClawBio/ClawBio
/plugin install clawbio
```

All skills are then available as agent-routable commands. Alternatively, clone the repo and open it as your working directory in Claude Code; the `AGENTS.md` at the repo root teaches Claude how to route requests to skills automatically.

> **Needs Claude Code 2.1.277+** (`claude --version`, `claude update`): native `AGENTS.md` support landed there, and older versions read only `CLAUDE.md`, which this repo no longer ships. Every other agent reads `AGENTS.md` regardless of version.

## Slash Commands

The repository also ships reusable slash commands in [`commands/`](../commands/) for Claude Code and compatible agents:

| Command | Purpose |
|---------|---------|
| `/analyse` | Analyse a file or input with the appropriate ClawBio skill |
| `/new-skill` | Scaffold a new skill from the official template |
| `/list-skills` | List available skills from `skills/catalog.json` |
| `/run-demo` | Run a skill demo with built-in sample data |

## Try all skills

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

## Run with your own data

```bash
python clawbio.py run pharmgx --input my_23andme.txt --output results/
python clawbio.py run rnaseq --input counts.csv,metadata.csv --output results_rnaseq/
python clawbio.py run methylation --geo-id GSE139307 --output results_methylation/
```

## Run tests

```bash
uv sync --all-extras                            # install every optional extra
uv run pytest                                   # whole repo
uv run pytest skills/pharmgx-reporter/tests/ -v # one skill
```

## Dependencies

Core dependencies are declared in [`pyproject.toml`](../pyproject.toml) and pinned in [`uv.lock`](../uv.lock): biopython, pandas, numpy, scikit-learn, matplotlib, openai, pydeseq2, google-cloud-bigquery, google-auth, conda-lock, rocrate. Most skills run with just these.

`uv sync` installs everything in a reproducible virtual environment. To add or update a dependency, run `uv add <package>` (or edit `pyproject.toml` and re-run `uv sync`); commit the resulting `uv.lock` change.

Some skills have additional requirements, declared in that skill's own `SKILL.md` `install:` block and installed only if you use that skill:

| Skill | Extra dependency | Install |
|-------|-----------------|---------|
| Metagenomics | Kraken2, RGI, HUMAnN3 | Conda (see skill README) |
| Methylation Clock | PyAging | `pip install pyaging` |
| scRNA Embedding | scvi-tools | `pip install scvi-tools` |
| Galaxy Bridge | BioBlend | `pip install bioblend` |
| Affinity Proteomics | somadata, scipy, statsmodels, seaborn, scikit-learn | `pip install somadata scipy statsmodels seaborn scikit-learn` |
| Cell Detection | cellpose, tifffile, czifile, nd2, Pillow, scikit-image | `pip install "cellpose>=4.0" tifffile "czifile>=2019.7.2.2" "nd2>=0.11.1" Pillow scikit-image` |
| Data Extractor | anthropic, opencv-python-headless | `pip install anthropic opencv-python-headless` |
| eQTL Catalogue / GWAS Catalog Region Fetch | pysam | `pip install pysam pandas requests` |
| Variant Annotation | pysam | `uv add pysam requests` |
| Celltype Specificity Profiler | scanpy, anndata | `uv add scanpy anndata numpy scipy pandas` |
| Drug Repurposing Screen | pyarrow (parquet engine) | `pip install numpy pandas scipy pyyaml pyarrow` |
| Proteomics Clock | seaborn | `pip install pandas numpy matplotlib seaborn requests` |
| RoboTerri (`robotary/`) | fastapi | **Undocumented in its own SKILL.md as of this writing.** `uv add fastapi` |

No Docker or Singularity required for core functionality. Skills that need external bioinformatics tools document their setup in their own `SKILL.md`.
