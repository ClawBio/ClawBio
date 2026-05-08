# Changelog

All notable changes to ClawBio are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

### New Skills
- **nfcore-scrnaseq-wrapper** (`skills/nfcore-scrnaseq-wrapper/`, `scrnaseq-pipeline`): Upstream single-cell RNA-seq preprocessing from FASTQ using nf-core/scrnaseq. Supports six presets (simpleaf/standard, STARsolo/star, kallisto, cellranger, cellrangerarc, cellrangermulti), strict preflight for Java/Nextflow/backend, samplesheet validation, `params.yaml`-driven execution, SHA-256 reproducibility bundle, and automatic handoff to `scrna-orchestrator` (via `--run-downstream`). Includes macOS/Apple Silicon Docker workaround. 282 tests.

## [v0.5.0] — 2026-04-04 — Validation & Benchmark Infrastructure

### Added
- **AD Ground Truth Benchmark Set** (`tests/benchmark/ad_ground_truth.json`): Curated set of 34 positive Alzheimer's disease genes across 3 evidence tiers (4 Mendelian causal, 20 GWAS-replicated from Bellenguez 2022, 10 novel Bellenguez), 20 brain-expressed negative control genes, 10 lead variants with GRCh38 coordinates, and scoring criteria with minimum acceptable thresholds. This is the first disease-specific validation dataset for any agentic bioinformatics platform.
- **Mock API Server** (`tests/benchmark/mock_api_server.py`): Deterministic mock endpoints for Ensembl REST, GWAS Catalog, and ClinPGx APIs. Threaded HTTP server with context manager for test integration. Enables offline CI testing without rate limits or API drift. Inspired by StrongDM's simulated Slack/Jira pattern.
- **Benchmark Scorer** (`tests/benchmark/benchmark_scorer.py`): Scores pipeline outputs against ground truth using gene recovery rate, false discovery rate, precision, recall, F1, and tier-weighted composite score. CLI and Python API. Outputs markdown reports with tier breakdown.
- **Swappable Fine-Mapping Pipeline** (`tests/benchmark/finemapping_benchmark.py`): First autoresearch-style benchmark. Runs ABF and SuSiE fine-mapping on the same synthetic locus with known causal signals, scores each method on recall, precision, PIP concentration, credible set size, and composite score, picks the winner. Method registry pattern: adding FINEMAP or PolyFun requires only a single function. First result: SuSiE wins (composite=0.80) vs ABF (composite=0.65).
- **Nightly Sweep Benchmark Integration** (`scripts/nightly_demo_sweep.py`): Nightly demo sweep now collects gene lists from skill outputs and scores them against the AD ground truth. Reports gene recovery rate, FDR, precision, recall, F1, and tier breakdown in the sweep summary. Benchmark section appears when `--output` is used.
- **Red/Green TDD Mandate** (`CLAUDE.md`): All skill development and modification must use test-driven development. Tests first, watch them fail, implement, watch them pass. Contributing workflow updated to enforce this.
- **74 benchmark tests** across ground truth integrity (8), mock API responses (5), HTTP endpoints (6), benchmark scoring (9), fine-mapping locus generation (7), method runners (6), scoring logic (3), benchmark runner (3), and reference genome (27). All green.

### New Skills (since v0.4.0)
- **struct-predictor** (PR #102, @camlloyd): AlphaFold/Boltz protein structure prediction
- **cell-detection** (PR #101, @camlloyd): CellposeSAM cell segmentation from fluorescence microscopy
- **bigquery-public** (PR #93, @YonghaoZhao722): SQL against BigQuery public genomics datasets
- **clinical-variant-reporter** (PR #89, @RezaJF): ACMG/AMP variant classification
- **fine-mapping** (PR #88, @camlloyd): SuSiE and ABF statistical fine-mapping
- **labstep** (PR #84, @camlloyd): Labstep ELN bridge for experiments, protocols, inventory
- **protocols-io** (PR #83, @camlloyd): protocols.io search, retrieval, authentication

### Community
- **UK AI Agent Hackathon 2026 Winner**: Won the biggest prize at Europe's largest AI hackathon
- **Genomebook 3rd place at AI London hackathon** (20-21 Mar)
- **Bioinformatics Application Note submitted** (2 Apr via ScholarOne)
- **Nature feature interview** (Nicola Jones, 2 Apr): ClawBio as case study for vibe coding in science
- **15 contributors**, 108 forks, 579 GitHub stars
- **PHURI Workshop accepted** (22 Apr, Queen Mary University of London)
- **Google.org AI for Science LOI** in preparation with UKDRI (Nathan Skene PI)

### Workshops & Tutorials
- 5 tutorial tracks with Colab notebooks, slides, and docs pages
- Unified 25-slide deck for live delivery
- 30x WGS workshop with Corpas genome (Zenodo DOI: 10.5281/zenodo.19297389)
- All tutorials tested end-to-end 3 Apr

### Infrastructure
- **Corpas 30x WGS reference genome**: First-class resource with VCF subsets, QC baselines, 28 benchmark tests
- **Nightly demo sweep**: `scripts/nightly_demo_sweep.py` with catalog-driven execution, skip list for heavy deps, GitHub Actions integration
- **Skill catalog**: `scripts/generate_catalog.py` auto-generates `skills/catalog.json` (42 skills indexed)
- **170 new tests** across common library and 4 previously untested skills (PR #85)

### Security
- Token redaction filter for httpx logs
- Structured JSONL audit logging for usage analytics and security events
- Filesystem write restriction to PROJECT_ROOT
- Conversation history sanitisation and global error handler
- Disclaimer enforcement in all Telegram messages

## [v0.4.0] — 2026-03-10 — Galaxy Integration

### Added
- **Galaxy Bridge skill** — search, inspect, and run 8,000+ bioinformatics tools from usegalaxy.org through natural language
- **galaxy_catalog.json** — bundled index of all Galaxy tools for offline discovery (8,182 tools across 86 categories)
- **200 curated tool profiles** — structured markdown profiles for the most important Galaxy tools (FastQC, Kraken2, DESeq2, BWA-MEM2, etc.)
- **BioBlend integration** — remote tool execution on Galaxy via Python SDK with full reproducibility bundles
- **Demo mode** — `python galaxy_bridge.py --demo` runs simulated FastQC analysis offline (no API key needed)
- **Cross-platform chaining** — Galaxy tools chain with ClawBio skills (e.g., Galaxy VEP → PharmGx Reporter)
- **Galaxy tool count in catalog.json** — `galaxy_tool_count` field shows total accessible tools

## [v0.3.1] — 2026-03-05 — Agent-Friendly

### Added
- **llms.txt** — LLM-friendly project summary following the emerging `llms.txt` standard; lists all docs, skills, and entry points in a format optimised for AI agent context windows
- **AGENTS.md** — Universal guide for AI coding agents (Codex, Devin, Cursor, Claude Code, Copilot Workspace); covers setup, commands, code style, project structure, safety boundaries, and contribution workflow
- **Machine-readable skill catalog** — `skills/catalog.json` auto-generated by `scripts/generate_catalog.py`; indexes all 21 skills with name, version, status, dependencies, tags, and trigger keywords
- **Standardised SKILL.md files** — All 21 skill specifications upgraded to consistent YAML frontmatter schema with emoji, OS compatibility, install instructions, and structured methodology sections
- **Upgraded SKILL-TEMPLATE.md** — Best-practice template matching the new standardised format so new contributors start right
- **Agent pointers in README and CONTRIBUTING** — Added references to `llms.txt`, `AGENTS.md`, and `catalog.json` so both human and AI contributors can find agent-specific documentation

## [v0.3.0] — 2026-03-01 — Imperial College AI Agent Hack

### Added
- Video introduction of ClawBio to Peter Steinberger at the UK AI Agent Hack, Imperial College London
- Security audit: 32 fixes for silent degradation across 4 production skills (`SECURITY-AUDIT.md`)
- README overhaul with demo video, provenance section, and architecture diagram

## [v0.2.0] — 2026-02-28 — Tests, CI, and ClawHub

### Added
- Test suites: 57 tests across PharmGx Reporter (24), Equity Scorer (24), NutriGx Advisor (9)
- GitHub Actions CI running on Python 3.10, 3.11, and 3.12 for every push and PR
- ClawHub registry: 3 skills published and installable via `clawhub install pharmgx-reporter`
- Org migration: repo moved to `github.com/ClawBio/ClawBio`
- Community infrastructure: issue templates, PR template, Discussions seeded, 8 open skill issues

[v0.5.0]: https://github.com/ClawBio/ClawBio/compare/v0.3.1...v0.5.0
[v0.4.0]: https://github.com/ClawBio/ClawBio/compare/v0.3.1...v0.4.0
[v0.3.1]: https://github.com/ClawBio/ClawBio/compare/v0.3.0...v0.3.1
[v0.3.0]: https://github.com/ClawBio/ClawBio/compare/v0.2.0...v0.3.0
[v0.2.0]: https://github.com/ClawBio/ClawBio/releases/tag/v0.2.0
