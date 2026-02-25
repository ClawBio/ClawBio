# Contributing to ClawBio

We welcome skills from anyone working in bioinformatics, computational biology, or related fields.

## How to Contribute a Skill

### 1. Copy the template

```bash
cp -r templates/SKILL-TEMPLATE.md skills/your-skill-name/SKILL.md
```

### 2. Define your skill

Edit `SKILL.md` with:
- **YAML frontmatter**: name, description, dependencies (bins, env vars, packages)
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

## Skill Ideas We Need

If you are looking for something to build:

- **GWAS Pipeline**: PLINK/REGENIE automation
- **Metagenomics Classifier**: Kraken2/MetaPhlAn wrapper
- **Pathway Enricher**: GO/KEGG enrichment analysis
- **Clinical Variant Reporter**: ACMG classification
- **Phylogenetics Builder**: IQ-TREE/RAxML automation
- **Proteomics Analyser**: MaxQuant/DIA-NN wrapper
- **Spatial Transcriptomics**: Visium/MERFISH analysis

## Questions?

Open an issue or reach out via the repo discussions.
