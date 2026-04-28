# 🦖 Contributing to ClawBio

We welcome skills from anyone working in bioinformatics, computational biology, or related fields.

**Join the contributors community: [t.me/ClawBioContributors](https://t.me/ClawBioContributors)**

## How to Contribute a Skill

### 1. Copy the template

```bash
mkdir -p skills/your-skill-name/tests skills/your-skill-name/examples
cp templates/SKILL-TEMPLATE.md skills/your-skill-name/SKILL.md
```

### 2. Define your skill

Use [`templates/SKILL-TEMPLATE.md`](templates/SKILL-TEMPLATE.md) as the structural source of truth.

Edit `SKILL.md` with:
- **Top-level YAML fields**: `name`, `description`, and `license`
- **`metadata` fields**: `version`, `author`, `domain`, `tags`, `inputs`, `outputs`, `dependencies`, `demo_data`, and `endpoints`
- **`metadata.openclaw` fields**: runtime and routing metadata such as `requires`, `always`, `emoji`, `homepage`, `os`, `install`, and `trigger_keywords`
- **Markdown body**: follow the template sections and keep them aligned with the PR audit expectations in [`CLAUDE.md`](CLAUDE.md)

If there is any discrepancy between documents, follow the template for structure and `CLAUDE.md` for audit requirements.

### 3. Add supporting code (optional)

If your skill needs Python/R scripts, add them alongside the SKILL.md:

```
skills/your-skill-name/
├── SKILL.md           # Required
├── your_script.py     # Optional
├── tests/             # Encouraged
│   └── test_script.py
└── examples/          # Encouraged
    ├── input.csv
    └── expected_output.md
```

If your skill requires a conda environment, add an `environment.yml` using `conda-forge` as the sole channel and `nodefaults` to prevent fallback to the Anaconda `defaults` channel (which has commercial licensing restrictions):

```yaml
channels:
  - conda-forge
  - nodefaults
```

### 4. Test locally

```bash
# Confirm the runner is working
python clawbio.py list

# Test a registered skill through the real CLI
python clawbio.py run <registered-alias> --demo

# Test a new skill directly before registration
python skills/your-skill-name/your_skill.py --demo --output /tmp/your-skill-demo
```

Use the direct script path until the skill is registered in `clawbio.py`. Once you add the alias to `clawbio.py`, confirm it appears in `python clawbio.py list` and then validate it through `python clawbio.py run <registered-alias> --demo`.

If the skill includes tests, run:

```bash
python -m pytest skills/your-skill-name/tests/ -v
```

If you changed `SKILL.md` YAML frontmatter, regenerate the catalog:

```bash
python scripts/generate_catalog.py
```

### 5. Submit

**Option A: Pull request to this repo**
```bash
git checkout -b add-your-skill-name
git add skills/your-skill-name/
git commit -m "Add your-skill-name skill"
git push -u origin add-your-skill-name
# Open PR on GitHub
```

**Option B: Submit to ClawHub**
Follow the [ClawHub submission guide](https://clawhub.ai/docs/submit).

## Skill Guidelines

1. **Local-first**: No mandatory cloud data uploads. Network calls only for public databases (PubMed, PDB, UniProt).
2. **Reproducible**: Generate audit logs and reproducibility bundles.
3. **One job well**: Each skill does one thing. Compose via the Bio Orchestrator.
4. **Documented**: Include example queries, expected outputs, and dependency lists.
5. **Safe**: Minimal permissions. Warn before destructive actions. No hardcoded credentials.

## Naming Conventions

- Skill folder: lowercase, hyphens (`vcf-annotator`, not `VCF_Annotator`)
- Python files: lowercase, underscores (`equity_scorer.py`)
- Skill name in YAML: matches folder name exactly

## Code Standards

- Python 3.10+
- Type hints encouraged
- pathlib for all file paths
- No hardcoded absolute paths
- Tests with pytest

## For AI Agents Contributing Skills

AI coding agents (Codex, Devin, Claude Code, Cursor, etc.) should follow the same workflow as human contributors, plus:

1. Read [`AGENTS.md`](AGENTS.md) for setup, commands, code style, and project structure
2. Read the target skill's `SKILL.md` before modifying any code
3. Use `python clawbio.py list` to verify skills still load after changes
4. Run `python -m pytest -v` to confirm all tests pass
5. Regenerate `skills/catalog.json` if you changed any SKILL.md YAML frontmatter: `python scripts/generate_catalog.py`

### SKILL.md Quality Checklist

Treat this as a contributor summary, not as a second source of truth.

- [`templates/SKILL-TEMPLATE.md`](templates/SKILL-TEMPLATE.md) defines the canonical `SKILL.md` structure.
- [`CLAUDE.md`](CLAUDE.md) defines the formal PR audit and conformance checklist.
- Every skill should cover `Trigger`, `Scope`, `Workflow`, `Example Output`, `Gotchas`, `Safety`, and `Agent Boundary`.
- Include synthetic demo data and support `--demo` whenever the skill has executable automation.
- Add tests for demo mode and the main expected path when the skill includes code.
- Keep dependencies, installation steps, inputs, outputs, and safety boundaries accurate and specific.
- Regenerate `skills/catalog.json` whenever you change YAML frontmatter.

If you are unsure whether a `SKILL.md` is PR-ready, defer to the checklist in `CLAUDE.md`.

## 🦖 Skill Ideas We Need

If you are looking for something to build:

- **GWAS Pipeline**: PLINK/REGENIE automation
- **Metagenomics Classifier**: Kraken2/MetaPhlAn wrapper
- **Pathway Enricher**: GO/KEGG enrichment analysis
- **Clinical Variant Reporter**: ACMG classification
- **Phylogenetics Builder**: IQ-TREE/RAxML automation
- **Proteomics Analyser**: MaxQuant/DIA-NN wrapper
- **Spatial Transcriptomics**: Visium/MERFISH analysis

## External Validation

ClawBio skills are independently audited by [clawbio_bench](https://github.com/biostochastics/clawbio_bench), a standalone safety, correctness, and honesty benchmark maintained by Sergey Kornilov (Biostochastics, LLC). The benchmark runs in CI on every PR and produces tamper-evident verdict artifacts.

If you find a scientific error in ClawBio, you can:
- Open an issue on this repo with the finding
- Add a test case to clawbio_bench (the external standard)
- Both (preferred)

We commit to triaging all externally reported findings and referencing the finding ID in fix commits.

## Acknowledgements

- **Sergey Kornilov** ([Biostochastics](https://github.com/biostochastics)) for building the independent [clawbio_bench](https://github.com/biostochastics/clawbio_bench) audit suite and identifying critical correctness and honesty findings across multiple skills.

## Questions?

Open an issue or reach out via the repo discussions.
