# Changelog

All notable changes to ClawBio are documented in this file.

The format is based on [Keep a Changelog](https://keepachangelog.com/en/1.1.0/),
and this project adheres to [Semantic Versioning](https://semver.org/spec/v2.0.0.html).

## [Unreleased]

## [0.7.1] - 2026-09-05

### Fixed
- **Tagged releases publish themselves again.** Every previous run of
  `publish.yml` failed with `invalid-publisher`, so 0.5.2, 0.6.1 and 0.7.0 were
  uploaded by hand. The workflow was correct: the OIDC claims GitHub sends match
  what `packaging/README.md` documents, and no matching publisher is registered
  on PyPI. Until it is, the job authenticates with an API token held as an
  environment secret on `pypi`, an environment that accepts `v*` tags only.

### Added
- **`SECURITY.md`**: a private disclosure path (GitHub private vulnerability reporting,
  now enabled on the repository, with an email fallback), response commitments, what is
  in and out of scope, and which versions receive fixes.
- **`docs/data-handling.md`**: every skill that can send data off the machine, grouped by
  what it sends (reference downloads only, query terms you typed, your variants, whole
  sequences or files, or chat messages to a hosted LLM), with the host, the trigger, the
  credential variable and the offline behaviour for each. Enforced by
  `tests/test_data_handling_doc.py`, which scans every skill for outbound-call code and
  fails if a networked skill is absent from the page.

### Changed
- README and `docs/architecture.md` no longer claim that "no network calls" are made for
  data processing. Several skills send variants to Ensembl VEP, gnomAD and ClinVar, and the
  six `gi-*` skills upload whole sequences to a hosted model. The privacy statement now says
  exactly that and points at the data-handling page. The README versioning line also read
  `v0.5.0` two releases late; it now says `v0.7.0`.

### Removed
- **`soul2dna` is no longer a ClawBio skill** (#111). The skill folder was a thin
  wrapper around the Genomebook sandbox compiler and presented "compile character
  profiles into synthetic genomes" as a catalogued bioinformatics capability, which
  the external audit rightly called out: there is no scientific basis for mapping
  traits of character to genotypes, and a skill catalogue that hopes to be trusted
  with clinical work should not carry one that implies there is. The compiler itself
  stays where it always lived, `GENOMEBOOK/PYTHON/01-soul2dna.py`, inside the
  sandbox whose README now states its scope: synthetic fixtures for exercising the
  orchestrator, with no biological meaning. `genome-match` and `recombinator`, which
  consume those fixtures, are unchanged. Catalogue count 97 to 96.

### Fixed
- **gwas-prs no longer substitutes the curated T2D panel for PGS000013** (#356,
  remaining half). `--pgs-id PGS000013` now goes to the PGS Catalog like any other
  accession and, when it cannot, fails with a message saying what the accession is
  (Khera 2018, coronary artery disease, 6,630,150 variants) and how to request the
  curated panel instead (`--panel-id CLAWBIO-T2D-8`). The benchmark compatibility
  alias still exists but is opt-in through `CLAWBIO_ALLOW_LEGACY_PGS_ALIAS=1`, which
  only `.github/workflows/bench-leaderboard.yml` sets for the pinned `clawbio_bench`
  revision. A network failure on that path now reports cleanly instead of raising.

## [0.7.0] - 2026-09-03

### Deprecated
- **MCP server** (`clawbio mcp`, the `[mcp]` extra): deprecated as of 0.7.0 and
  scheduled for removal in 0.8.0. Nothing breaks yet: every existing client
  configuration keeps working through the 0.7.x series, and the server prints a
  notice on stderr when it starts (never stdout, which is the stdio protocol channel).
  Why: the skills are plain [Agent Skills](https://agentskills.io) folders, and the
  editors the server was built for in 0.6.0 (Cursor, VS Code, Codex, Zed) now load such
  folders directly from `~/.agents/skills/`, so a remote-call shim in front of them adds a process, a dependency
  pin (`mcp<2`) and a second code path for no capability the skills lack. Migration:
  Claude Code users install the plugin (`/plugin marketplace add ClawBio/ClawBio`,
  `/plugin install clawbio`); other editors copy or symlink the folders they need from
  `skills/` into `~/.agents/skills/`; everything else uses the CLI or the
  Python API. Skills that wrap third-party MCP servers (`bioqc-mcp`, `bgpt-mcp`,
  `just-prs-mcp`) are unaffected; this covers ClawBio's own server only.
- `server.json`, the manifest for the official MCP Registry, is removed. The server was
  never published there and a deprecated server must not be.

### Added
- **Claude Code plugin manifests are current again.** `.claude-plugin/plugin.json` and
  `marketplace.json` were frozen at 0.3.0 and "24 skills" while PyPI shipped 0.6.1 and
  the catalog held 97, so `/plugin install clawbio` reported a version three releases
  old. Both now read 0.7.0 with counts taken from `skills/catalog.json` and the CLI
  registry (97 Agent Skills, 50 with a CLI entry point), and
  `clawbio/tests/test_packaging.py` fails any future release whose manifests or
  `CITATION.cff` disagree with `clawbio.__version__`. `commands/new-skill.md` gained
  the frontmatter `claude plugin validate --strict` requires.
- **Reproducibility bundles as a shared contract.** SHA-256 helpers moved to
  `clawbio.common.checksums` (#351, @camlloyd); `scaffold_skill.py` now generates a
  working bundle for every new skill (#374, @camlloyd); `equity-scorer` (#352,
  @camlloyd), `nutrigx` (#350, @camlloyd) and `gwas-prs` (#379, @krudo-taco) write one.
- **New skills**: `deepspot-m`, virtual spatial transcriptomics from H&E tiles (#331,
  @KalinNonchev); `just-prs-mcp`, evidence-aware VCF/WGS PRS through a pinned local
  just-prs server (#324, @antonkulaga).
- A Token Factory agent example driving ClawBio skills from a Nebius-hosted model, with
  the report-upload behaviour disclosed (#336, #366).
- Artefact licences are stated separately from the wrapper licence (#339).

### Fixed
- `pharmgx-reporter`: three reference positions were GRCh37 in a GRCh38 table, which
  blocked GRCh38 input (#325).
- `clinical-variant-reporter`: fails closed on VEP annotation batch failures (#368,
  @krudo-taco) and reports the live Ensembl release in its Data Sources table (#327,
  @AmirF194).
- `gwas-prs`: corrected panel citations, stopped demo panels shadowing real PGS Catalog
  scores (#357) and gave curated panels honest identifiers (#380, @krudo-taco).
- `nutrigx`: rs4988235 lactose polarity and dominance were inverted (#376,
  @krudo-taco); catalog version synced (#377, @krudo-taco).
- `rnaseq-de`: keeps gene IDs on the nf-core count handoff (#371) and verifies LFC
  shrinkage matches the requested contrast (#378, both @krudo-taco).
- `nfcore-rnaseq-wrapper`: manifest version parse scoped to the manifest block, so
  `custom_config_version` or `vep_version` can no longer be mistaken for it (#334).
- `pathway-enricher`: submits gene lists as multipart form data (#381, @krudo-taco).
- `gi-*` skills: strand guidance and stale expression values corrected (#354, @boldakov).
- `claw-semantic-sim`: dropped a false "works out of the box" demo claim (#358,
  @AmirF194). `vcf-annotator` CLI fixed (#336).

### CI
- The `test` job runs the tests of any skill a PR touches, diffed against the PR's own
  merge-base (#330, #332), and no longer assumes every skill is Python (#361).
- Queued workflow runs are auto-approved for PRs that touch only `skills/` (#362).
- `scientific-audit` can fail: it gates on a committed bench baseline (#343).
- `vcf-annotator` tests added to the main allowlist (#370, @krudo-taco).

## [0.6.1] - 2026-07-29

### Fixed
- The MCP server was unusable when installed from PyPI. `mcp_server.py` computed
  its own `skills/` path relative to the package parent, which is correct only in
  a source checkout; an installed wheel bundles `skills/` *inside* the package.
  Every tool failed with `FileNotFoundError` on `catalog.json`. It now reuses
  `clawbio.cli.SKILLS_DIR`, which already resolves both layouts. v0.6.0 is
  broken for MCP use and should be skipped; the CLI was unaffected.

## [0.6.0] - 2026-07-29

### Added
- **MCP server** (`clawbio mcp`): ClawBio is now usable from any Model Context Protocol
  client, including Cursor, Zed, VS Code agent mode and Claude Desktop, not just Claude
  Code. Runs locally over stdio; no hosted endpoint exists by design, so genomic data
  never traverses a third-party server. Three tools:
  - `clawbio_list_skills(query)` searches the catalog and flags each result `runnable`,
    since agent-readable (`SKILL.md`-only) skills cannot be executed.
  - `clawbio_describe_skill(name)` returns the full `SKILL.md` contract; accepts a skill
    name or CLI alias.
  - `clawbio_run_skill(...)` executes a skill and returns its structured result.

  Install with the new optional extra, or run without installing:

      pip install 'clawbio[mcp]'
      uvx --from 'clawbio[mcp]' clawbio mcp

  **Demo data only by default.** Passing `input_path` or `output_dir` is refused unless
  `CLAWBIO_MCP_ALLOW_LOCAL_FILES=1` is set. Adding an MCP server to a client config is a
  low-friction, easily-forgotten action and must not silently grant an agent read access
  to a patient genome. Documented at
  [docs.clawbio.ai/reference/mcp](https://docs.clawbio.ai/reference/mcp/).
- `server.json` manifest for publication to the official MCP Registry.

### Fixed
- `skills/catalog.json` was stale at 94 skills; `llm-bench` (`skills/llm-biobank-bench/`)
  was CLI-registered but absent from it. Regenerated to the true **95 skills** (29 MVP,
  66 planned), which also fixes the failing `test_checked_in_catalog_is_current`.
- Public skill and test counts were inconsistent across `README.md`, `llms.txt` and
  `index.html` (variously 78, 88 and 94 skills; 2,318 and 4,183 tests). All now derive
  from `skills/catalog.json` and agree: 95 skills, 89 with runnable demo data, 8,182
  Galaxy tools, 4,217 tests. `llms.txt` records the commands that regenerate them.
- README described 88 skills as "production-ready"; that figure is the catalog's
  `has_demo` count, and the `maturity_tier` evidence records 10 `ci-validated` and 38
  `tested`. Reworded to "with runnable demo data", which is what the field measures.

### Changed
- The `[mcp]` extra pins `mcp>=1.9,<2`: `mcp` 2.0.0 removed `mcp.server.fastmcp`.

### New Skills (merged between 0.5.2 and 0.6.0)
- **cnv-acmg-classifier** (`skills/cnv-acmg-classifier/`, `cnv-acmg`): Germline CNV/SV (deletion/duplication) ClinGen/ACMG 2019 (Riggs 2020) point classification. Strictly additive Sections 1-5 (no terminal short-circuit; a complete 2A inherited from an unaffected parent = 0.70 VUS, a de novo 2A gain = 1.45 Pathogenic); partial-overlap 2C/2D/2E sub-calls derived from breakpoint geometry (gene strand + coding boundaries); Section 3 omitted only on a complete 2A/2F; symmetric Benign<=-0.99 boundary; CN1 sex-chromosome guard. VCF or CSV/TSV input, swappable dosage-map/gene-model, full reproducibility bundle. Stdlib-only, local-first. 24 tests.
- **nfcore-scrnaseq-wrapper** (`skills/nfcore-scrnaseq-wrapper/`, `scrnaseq-pipeline`): Upstream single-cell RNA-seq preprocessing from FASTQ using nf-core/scrnaseq. Supports six presets (simpleaf/standard, STARsolo/star, kallisto, cellranger, cellrangerarc, cellrangermulti), strict preflight for Java/Nextflow/backend, samplesheet validation, `params.yaml`-driven execution, SHA-256 reproducibility bundle, and automatic handoff to `scrna-orchestrator` (via `--run-downstream`). Includes macOS/Apple Silicon Docker workaround. Audit round 5 hardening: `--demo` is now fully hermetic (only the four forced essentials reach `params.yaml`; all QC/skip/tuning/save/reporting flags are ignored and warned); FastQC is a required output for **every** aligner — including the Cell Ranger family — matching the 4.1.0 workflow (FASTQC on the shared `ch_fastq` before aligner branching); `feature_type=crispr` now requires `--fb-reference` (shared feature reference with antibody capture); the `-c` params-override lint also catches bracket/whole-map/newline-block forms; a wall-clock cap left active on object-store/institutional runs emits a `--timeout-hours 0` hint; `sample_type`/`feature_type` enum values are validated whenever present under **any** preset (matching `assets/schema_input.json` property-level enums) so an invalid value fails fast in preflight rather than late in Nextflow; and SKILL.md docs corrected (preset-conditional samplesheet columns, `--transcript-fasta`/`--txp2gene` for simpleaf, `--gex-barcode-sample-assignment` is not an OCM selector, STARsolo-velocity example now passes `--protocol`). 387 tests.
- **nfcore-rnaseq-wrapper** (`skills/nfcore-rnaseq-wrapper/`, `rnaseq-pipeline`): Upstream bulk RNA-seq preprocessing from FASTQ/BAM using nf-core/rnaseq v3.26.0. Supports STAR+Salmon, STAR+RSEM, HISAT2, and Bowtie2+Salmon routes; strict preflight for Java/Nextflow/backend, samplesheet strandedness and references; `params.yaml`-driven execution; SHA-256 reproducibility bundle; provenance JSONs; and template handoff to `rnaseq-de`. Hardening round: contaminant screening (`--contaminant-screening`, `--kraken-db`, `--sylph-db`, `--bracken-precision`, BBSplit auto-enable), iGenomes name validation with fast-fail in preflight, GENCODE GTF auto-detect, real `duration_seconds` measurement, auto-handoff to `rnaseq-de` (`--run-downstream --metadata --formula --contrast`), `--prokaryotic` restricted to profile modifier (never standalone backend), `--check` guaranteed to never invoke Nextflow, passthrough flags `--enable-preseq`, `--multiqc-config`, `--multiqc-logo`, `--rsem-extra-args`. 538 tests.
- **nfcore-sarek-wrapper** (`skills/nfcore-sarek-wrapper/`, `sarek-pipeline`): nf-core/sarek v3.8.1 wrapper with step-aware restart sheets, somatic/germline validation, caller and annotation resources, output discovery, and portable reproducibility bundles. Alignment audit hardening includes effective iGenomes resources with the documented `false` sentinel, final `--extra-param` precedence, `--outdir-cache` cache-download preflight, full integrated CLI help/forwarding, and exact portable replay of the captured Nextflow invocation. 321 tests.


## [v0.5.2] - 2026-06-10 - pip / conda packaging

### Packaging
- ClawBio is now installable with `pip install clawbio`, and a bioconda recipe is ready for `conda install -c bioconda clawbio`. The CLI engine moved into the importable `clawbio` package with a `clawbio` console entry point; all skills' logic is bundled into the wheel (full demo data for the headline skills); output and patient profiles route to the working directory when installed; and the version is single-sourced from `clawbio/__init__.py`. First PyPI release was 0.5.1; 0.5.2 additionally bundles the skills below.

### New Skills
- **drug-repurposing-screen**, **pathway-enricher**, **phylogenetics-builder**, **bioqc-mcp**: contributor skills now bundled in the installable package.

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
