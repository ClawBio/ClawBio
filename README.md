<h3 align="center">ClawBio</h3>

<p align="center">
  <strong>A bioinformatics-native AI agent skill library.</strong><br>
  Built on <a href="https://github.com/openclaw/openclaw">OpenClaw</a> (180k+ GitHub stars). Local-first. Privacy-focused. Reproducible.
</p>

<p align="center">
  <a href="https://github.com/ClawBio/ClawBio/actions/workflows/ci.yml"><img src="https://github.com/ClawBio/ClawBio/actions/workflows/ci.yml/badge.svg" alt="CI"></a>
  <a href="#install"><img src="https://img.shields.io/badge/python-3.12+-blue?logo=python&logoColor=white" alt="Python 3.12+"></a>
  <a href="skills/catalog.json"><img src="https://img.shields.io/badge/dynamic/json?url=https%3A%2F%2Fraw.githubusercontent.com%2FClawBio%2FClawBio%2Fmain%2Fskills%2Fcatalog.json&query=%24.skill_count&label=skills&color=3fb950" alt="Skill count, read live from skills/catalog.json"></a>
  <a href="LICENSE"><img src="https://img.shields.io/badge/license-MIT-green" alt="MIT License"></a>
  <a href="https://luma.com/clawbio"><img src="https://img.shields.io/badge/Events-Follow_on_Luma-7c3aed" alt="Follow ClawBio Events on Luma"></a>
  <a href="https://clawbioworkspace.slack.com/channels/all-clawbio"><img src="https://img.shields.io/badge/slack-%23all--clawbio-4A154B?logo=slack&logoColor=white" alt="Get help on the ClawBio Slack"></a>
  <a href="https://doi.org/10.5281/zenodo.19420648"><img src="https://zenodo.org/badge/DOI/10.5281/zenodo.19420648.svg" alt="DOI"></a>
  <a href="https://github.com/ClawBio/ClawBio/issues"><img src="https://img.shields.io/github/issues/ClawBio/ClawBio" alt="Open Issues"></a>
  <a href="https://clawbio.github.io/ClawBio/slides/"><img src="https://img.shields.io/badge/slides-London_Bioinformatics_Meetup-purple" alt="Slides"></a>
</p>

<p align="center"><strong>📅 <a href="https://luma.com/clawbio">Follow the ClawBio Events Calendar</a></strong> for hackathons, workshops and meetups.</p>

---

## Install

```bash
pip install clawbio                 # Python 3.12+
clawbio run pharmgx --demo          # runs in under 2 seconds
```

- **[Conda](https://docs.conda.io/):** `conda install -c bioconda clawbio`
- **Python library:** `from clawbio import run_skill; run_skill("pharmgx", demo=True)`
- **[Claude Code](https://claude.ai/claude-code) plugin:** `/plugin marketplace add ClawBio/ClawBio` (needs Claude Code 2.1.277+)
- **Cursor, VS Code, Codex, Zed or any editor that reads [Agent Skills](https://agentskills.io):** every skill is a plain Agent Skills folder. Copy or symlink the ones you need from `skills/` into `~/.agents/skills/` or your project's `.agents/skills/`. No server involved.
- **All skills with full demo data, or developing ClawBio:** work from a source checkout ([uv](https://docs.astral.sh/uv/) recommended):

```bash
git clone https://github.com/ClawBio/ClawBio.git && cd ClawBio
uv sync                             # installs from pyproject.toml + uv.lock
uv run python clawbio.py run pharmgx --demo
```

Library calls, slash commands, every demo, running on your own data, tests and per-skill dependencies: [docs/usage.md](docs/usage.md).

<p align="center">
  <img src="img/clawbio-demo.gif" alt="ClawBio GWAS Lookup demo — querying 9 genomic databases from the terminal" width="700">
</p>

---

## What ClawBio does

Snap a photo of a medication in Telegram. ClawBio identifies the drug from the packaging, queries your pharmacogenomic profile from [your own genome](docs/demo-genome.md), and returns a personalised dosage card, on your machine, in seconds:

<p align="center">
  <img src="skills/drug-photo/demo-images/00-00-warfarin.jpg" width="360" alt="Warfarin 2mg medication packaging — ClawBio identifies the drug from a photo and returns a personalised pharmacogenomic report">
</p>

> **Warfarin** | CYP2C9 \*1/\*2 Intermediate · VKORC1 High Sensitivity
> **AVOID — DO NOT USE** · Standard dose causes over-anticoagulation in this genotype.

Or look up a variant (such as [rs9923231](https://www.ncbi.nlm.nih.gov/snp/rs9923231)) across nine genomic databases at once, combine thousands of small-effect variants into a polygenic risk score, or ask the [UK Biobank](https://www.ukbiobank.ac.uk/) schema in plain English which fields measure blood pressure or grip strength. More examples: [docs/skill-highlights.md](docs/skill-highlights.md).

## Why skills, not just an LLM

Ask an LLM to "profile the pharmacogenes in my VCF file" and it writes plausible Python that may use outdated CPIC guidelines or misclassify star alleles. CYP2D6\*4 is a no-function allele; calling it reduced-function decides whether a patient gets a dose adjustment or is told to avoid a drug. About 7% of people of European ancestry are CYP2D6 poor metabolisers, for whom codeine gives no pain relief, and about 0.5% carry DPYD variants where standard fluorouracil dosing can be lethal.

Agentic bioinformatics systems usually address the **reasoning layer** (knowledge graphs, fine-tuning) or the **connectivity layer** (tools wrapped as MCP servers). ClawBio addresses the **specification layer**. A **skill** is a folder holding a declarative specification (`SKILL.md`), validated Python, demo data and tests. The specification is a contract, not a prompt: it encodes a domain expert's analytical decisions, so the LLM orchestrates but does not improvise.

- **Specification-first**: domain expertise lives in versioned, peer-reviewable `SKILL.md` files, not in model weights.
- **Agent-agnostic**: skills run the same whether Claude, ChatGPT or a local model via Ollama invokes them.
- **Local-first by default**: Your genome is analysed on your machine, and the `clawbio` package itself makes no network calls. Some skills do reach public annotation APIs or hosted models, and a few send your variants or sequences to do so. Every one of them is named, with what it sends and to whom, in [docs/data-handling.md](docs/data-handling.md), and a test fails if a networked skill is missing from that page.
- **Reproducible**: many skills write a `reproducibility/` bundle next to their report. Contents vary by skill, and some replays still need the original inputs ([docs/reproducibility.md](docs/reproducibility.md)):

```
report/
├── report.md              # analysis with figures and tables
├── figures/  tables/
└── reproducibility/
    ├── commands.sh        # replay command for this run
    ├── environment.yml    # suggested Conda environment
    ├── checksums.sha256   # SHA-256 of selected outputs
    └── runtime-lock.json  # extra lock metadata, where supported
```

## Reference genome

ClawBio's demo data is built on a real, fully open human genome: the **Corpasome**. The [23andMe SNP chip](docs/demo-genome.md) (~600K variants) has been available since launch. Now, the project also ships subsets from a **30x Illumina whole-genome sequence** (GRCh37), covering ~4M SNPs, ~600K indels, and structural variants (DEL, DUP, INV, BND, INS, CNVs). All data comes from a single individual (Manuel Corpas), licensed CC0, and published on Zenodo ([doi:10.5281/zenodo.19297389](https://doi.org/10.5281/zenodo.19297389)). This dataset is provided for research and educational purposes only.

See [docs/reference-genome.md](docs/reference-genome.md) for use cases, subsets, and citation details.

## Validation

Skills are scored against curated ground truth with objective metrics. The public leaderboard tracks **161 / 182 tests passing (88.5%)** across 10 audited skills, up from 80 / 140 (57.1%) at the original audit: see [benchmarks.html](https://clawbio.ai/benchmarks.html). The suite, its ground truth and how to run it: [docs/benchmarking.md](docs/benchmarking.md).

## Skills

Browse every skill in the [skill directory](https://clawbio.ai/skills.html) or [`skills/catalog.json`](skills/catalog.json), or run `python clawbio.py list`.

- **Hosted inference:** the `gi-*` skills (promoter, splice, enhancer, chromatin, expression, annotation) call the hosted [Genomic Intelligence](https://genomicintelligence.ai) API (`api.genomicintelligence.ai`), unlike most skills. A shared hackathon key ships in `.env.example` (50 concurrent / 120 rpm, opt-in only). API reference: <https://docs.genomicintelligence.ai>; production keys: **contact@genomicintelligence.ai**.
- **Galaxy:** [Galaxy Bridge](skills/galaxy-bridge/) searches and runs 8,000+ [usegalaxy.org](https://usegalaxy.org) tools from the CLI.
- **Genomebook:** a synthetic-genetics sandbox that turns character profiles into diploid genomes and breeds offspring: [docs/genomebook.md](docs/genomebook.md).

## Chat interfaces

**RoboTerri** runs ClawBio over Telegram against a real, CC0-licensed human genome: [try it, no install needed](https://t.me/RoboTerri_bot). For the strongest privacy, a self-hosted **OpenClaw gateway** puts a browser webchat in front of skills running on your own machine. Both, and how to run your own: [docs/conversational-interfaces.md](docs/conversational-interfaces.md). How the pieces fit: [docs/architecture.md](docs/architecture.md).

## Contributing

Wrap your pipeline as a skill and open a PR:

```bash
cp templates/SKILL-TEMPLATE.md skills/<your-skill-name>/SKILL.md
# Edit SKILL.md, add Python implementation, demo data, and tests
```

See [CONTRIBUTING.md](CONTRIBUTING.md). Hackathon work is coordinated on Slack in [`#all-clawbio`](https://clawbioworkspace.slack.com/channels/all-clawbio) ([invite](https://clawbio.ai/join/slack/)).

**Security and data handling.** Report vulnerabilities privately through [SECURITY.md](SECURITY.md), never in a public issue. Before running ClawBio on data you are responsible for, read [docs/data-handling.md](docs/data-handling.md): it lists every skill that can send data off your machine, what it sends, to which host, and which environment variable holds its credential.

## For AI agents

| Resource | Purpose |
|----------|---------|
| [`AGENTS.md`](AGENTS.md) | The single guide for AI coding agents: routing, CLI, demos, setup, development and safety rules |
| [`llms.txt`](llms.txt) | Token-optimized project summary ([llmstxt.org](https://llmstxt.org)) |
| [`skills/catalog.json`](skills/catalog.json) | Machine-readable skill index with trigger keywords and demo commands |

## Versioning

ClawBio follows [Semantic Versioning](https://semver.org/). The current release is **v0.7.1** (5 Sep 2026): tagged releases publish to PyPI on their own again. [Release notes](https://github.com/ClawBio/ClawBio/releases/tag/v0.7.1) · [CHANGELOG.md](CHANGELOG.md) · DOI [10.5281/zenodo.19420648](https://doi.org/10.5281/zenodo.19420648).

## Links

- **OpenClaw**: [github.com/openclaw/openclaw](https://github.com/openclaw/openclaw), the agent platform
- **ClawHub**: [clawhub.ai](https://clawhub.ai), skill registry
- **HEIM Index**: [heim-index.org](https://heim-index.org), Health Equity Index for Minorities
- **Talks**: [docs/talks/](docs/talks/README.md), including the [London Bioinformatics Meetup slides](https://clawbio.github.io/ClawBio/slides/)

## License

MIT. See [LICENSE](LICENSE).
