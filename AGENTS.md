# AGENTS.md — ClawBio Guide for AI Coding Agents

This guide is for AI coding agents (Codex, Devin, Cursor, Claude Code, Copilot Workspace, etc.) working on the ClawBio codebase. For Claude-specific routing and skill invocation instructions, see [CLAUDE.md](CLAUDE.md).

## Project Overview

ClawBio is a bioinformatics AI agent skill library built on [OpenClaw](https://github.com/OpenClaw/OpenClaw). Each skill is a self-contained module that performs a specific genomic analysis task (pharmacogenomics, ancestry PCA, GWAS lookup, etc.) via CLI or Python API. All processing is local-first — genetic data never leaves the machine.

## Setup

```bash
# Clone
git clone https://github.com/ClawBio/ClawBio.git && cd ClawBio

# Install dependencies
pip install -r requirements.txt

# Verify
python clawbio.py list                      # List all skills
python clawbio.py run pharmgx --demo        # Run flagship skill demo
```

## Commands

| Command | Purpose |
|---------|---------|
| `python clawbio.py list` | List all registered skills and their status |
| `python clawbio.py run <skill> --demo` | Run a skill with built-in demo data |
| `python clawbio.py run <skill> --input <file> --output <dir>` | Run a skill with user data |
| `python -m pytest -v` | Run all tests |
| `python -m pytest skills/<name>/tests/ -v` | Run tests for one skill |
| `python scripts/generate_catalog.py` | Regenerate `skills/catalog.json` |
| `python scripts/lint_skills.py` | Run SKILL.md conformance linter |
| `make test` | Alias for `python -m pytest -v` |
| `make demo` | Run PharmGx demo |

## Testing

Tests live in `skills/<name>/tests/`. The test paths are registered in `pytest.ini`. When adding a new skill with tests, add its test path there.

```bash
python -m pytest -v                                    # All tests
python -m pytest skills/pharmgx-reporter/tests/ -v     # Single skill
python -m pytest -k "test_demo" -v                     # By pattern
```

Tests must pass on Python 3.10, 3.11, and 3.12. CI runs all three via GitHub Actions.

## Code Style

- **Python**: 3.10+ (type hints encouraged, `X | None` syntax over `Optional[X]`)
- **Paths**: Always use `pathlib.Path`. Derive from `Path(__file__).resolve().parent`. Never hardcode absolute paths.
- **Naming**: Skill folders use lowercase-hyphens (`gwas-lookup`). Python files use lowercase_underscores (`gwas_lookup.py`). One exception: `nutrigx_advisor` (legacy).
- **Imports**: Sibling modules loaded via `importlib.util.spec_from_file_location` — no package structure.
- **CLI**: Every skill script accepts `--input`, `--output`, and `--demo`. Use `argparse`.
- **Output**: Skills write to `<output_dir>/report.md` (primary), plus `figures/`, `tables/`, `reproducibility/` subdirectories as needed. Return a `result.json` with structured findings.
- **Dependencies**: Add to `requirements.txt` for core deps. Skill-specific heavy deps go in the skill's SKILL.md YAML `install` section.

## Project Structure

```
ClawBio/
├── clawbio.py              # Main CLI runner + Python API (SKILLS dict here)
├── clawbio/                # Shared utilities package
├── skills/                 # One directory per skill
│   ├── pharmgx-reporter/   # Example MVP skill
│   │   ├── SKILL.md        # Skill specification (YAML frontmatter + methodology)
│   │   ├── pharmgx_reporter.py  # Implementation
│   │   ├── api.py          # Importable API (optional)
│   │   ├── demo_patient.txt     # Demo data
│   │   └── tests/
│   │       └── test_pharmgx.py
│   ├── gwas-lookup/        # Example with subpackages
│   │   ├── SKILL.md
│   │   ├── gwas_lookup.py
│   │   ├── api/            # API clients for 9 databases
│   │   ├── core/           # Normalisation, resolution, reporting
│   │   └── tests/
│   └── catalog.json        # Machine-readable skill index (auto-generated)
├── templates/
│   └── SKILL-TEMPLATE.md   # Template for new skills
├── examples/               # Shared demo data (VCFs, CSVs)
├── profiles/               # PatientProfile JSONs
├── output/                 # Generated reports (timestamped)
├── docs/                   # Tutorials and reference docs
├── bot/                    # RoboTerri Telegram integration
├── requirements.txt        # Core Python dependencies
├── pytest.ini              # Test configuration
├── Makefile                # Convenience targets
├── CLAUDE.md               # Claude-specific agent instructions
├── AGENTS.md               # This file
└── llms.txt                # LLM-friendly project summary
```

## Skill Architecture

Every skill is defined by a `SKILL.md` file with:
1. **YAML frontmatter** (`openclaw` schema) — name, description, version, dependencies, install instructions, emoji, OS compatibility
2. **Markdown body** — Core Capabilities, Workflow, Output Structure, Dependencies, Safety, Integration with Bio Orchestrator

Skills with Python implementations are registered in the `SKILLS` dict in `clawbio.py` (line ~252). Each entry maps a CLI alias to its script path, demo args, description, allowed flags, and capabilities.

The **Bio Orchestrator** (`skills/bio-orchestrator/`) routes user queries to the right skill based on file type and keywords.

### SKILL.md is the primary artifact

Even skills without Python scripts are usable — an AI agent reads the SKILL.md and applies the methodology directly. Python scripts are optional accelerators.

## How to Add a New Skill

1. **Create the directory**: `mkdir skills/<your-skill-name>`
2. **Copy the template**: `cp templates/SKILL-TEMPLATE.md skills/<your-skill-name>/SKILL.md`
3. **Fill in SKILL.md**: Complete all YAML frontmatter fields and every Markdown section
4. **Add Python implementation** (optional): Main script accepting `--input`, `--output`, `--demo`
5. **Add demo data**: Include a small synthetic dataset for `--demo` mode
6. **Add tests**: Create `tests/test_<name>.py` with at least demo-mode coverage
7. **Register in `clawbio.py`**: Add an entry to the `SKILLS` dict with script path, demo_args, description, and allowed_extra_flags
8. **Register in `pytest.ini`**: Add the test path to `testpaths`
9. **Regenerate catalog**: `python scripts/generate_catalog.py`
10. **Verify**: `python clawbio.py list` shows the skill; `python -m pytest` passes

## How to Modify an Existing Skill

1. **Read SKILL.md first** — understand the methodology before changing code
2. **Run existing tests**: `python -m pytest skills/<name>/tests/ -v`
3. **Make changes** — keep the skill's CLI interface stable (`--input`, `--output`, `--demo`)
4. **Update SKILL.md** if you changed capabilities, inputs, or outputs
5. **Run tests again** — all must pass
6. **Regenerate catalog** if SKILL.md YAML changed: `python scripts/generate_catalog.py`

## Git Workflow

- Branch from `main` — name branches `feat/<skill-name>` or `fix/<skill-name>`
- **One skill per PR** — keep changes focused
- Tests must pass on Python 3.10, 3.11, and 3.12
- Include demo output in the PR description so reviewers can verify
- Update `skills/catalog.json` if you changed any SKILL.md YAML frontmatter
- When opening a PR, use `.github/PULL_REQUEST_TEMPLATE.md` as the body — fill in every section, do not leave placeholders

## Safety Boundaries

These are non-negotiable constraints:

1. **Local-first**: Genetic data never leaves the machine. No cloud uploads, no external API calls with patient data.
2. **Disclaimer required**: Every report must include: *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*
3. **No hallucinated science**: Gene-drug associations, thresholds, and parameters must trace back to SKILL.md methodology or cited databases (CPIC, PharmGKB, ClinVar, etc.). Never invent bioinformatics values.
4. **Security filtering**: `clawbio.py` enforces per-skill `allowed_extra_flags` whitelists (INT-001). Do not bypass this.
5. **Warn before overwriting**: Check for existing output before writing to a directory.

## Key Files Reference

| File | Purpose |
|------|---------|
| `clawbio.py` | CLI runner, SKILLS dict, security filtering, profile management |
| `skills/catalog.json` | Machine-readable skill index (auto-generated) |
| `CLAUDE.md` | Claude-specific routing table and demo commands |
| `CONTRIBUTING.md` | Human contributor guide and wanted skills list |
| `requirements.txt` | Core Python dependencies |
| `pytest.ini` | Test path registration |
| `Makefile` | `make test`, `make demo`, `make list` |
