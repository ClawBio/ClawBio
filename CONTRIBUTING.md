# 🦖 Contributing to ClawBio

We welcome skills from anyone working in bioinformatics, computational biology, or related fields.

**Join the contributors community: [t.me/ClawBioContributors](https://t.me/ClawBioContributors)**

## How to Contribute a Skill

### 1. Copy the template

```bash
cp -r templates/SKILL-TEMPLATE.md skills/your-skill-name/SKILL.md
```

### 2. Define your skill

Edit `SKILL.md` with:
- **YAML frontmatter**: AgentSkills-compatible metadata. Keep `name`, `description`, and `license` at the top level, then place skill details such as `version`, `author`, `domain`, `tags`, `inputs`, `outputs`, `dependencies`, `demo_data`, and `endpoints` under the top-level `metadata:` block.
- **Markdown body**: Instructions the AI agent follows. Include capabilities, workflow steps, example queries, output format, and safety rules.

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

### 4. Test locally

```bash
# Install your skill
openclaw install skills/your-skill-name

# Test with a sample query
openclaw "Your example query here"
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

- Python 3.11+
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

Every SKILL.md should include these sections (check the template at [`templates/SKILL-TEMPLATE.md`](templates/SKILL-TEMPLATE.md)):

- [ ] **YAML frontmatter** with `openclaw` schema (name, description, version, tags, trigger_keywords)
- [ ] **AgentSkills-compatible layout** with skill metadata nested under `metadata` and OpenClaw routing fields under `metadata.openclaw`
- [ ] **Why This Exists** (what goes wrong without the skill)
- [ ] **Core Capabilities** (numbered list)
- [ ] **Input Formats** (table with format, extension, required fields)
- [ ] **Workflow** (numbered steps)
- [ ] **CLI Reference** (bash examples with `--input`, `--output`, `--demo`)
- [ ] **Demo** section with expected output description
- [ ] **Output Structure** (directory tree)
- [ ] **Dependencies** (required and optional)
- [ ] **Safety** (local-first, disclaimer, no hallucination)
- [ ] **Integration with Bio Orchestrator** (trigger conditions, chaining partners)
- [ ] **Citations** (databases and papers used)

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
