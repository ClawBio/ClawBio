# 🦖 ClawBio Skill Development Manual

This manual provides a detailed, step-by-step developer guide on how to build, test, and register a new genomic analysis skill in ClawBio.

---

## Table of Contents
1. [Overview & Architecture](#1-overview--architecture)
2. [Step 1: Planning and Naming Conventions](#2-step-1-planning-and-naming-conventions)
3. [Step 2: Scaffolding the Skill](#3-step-2-scaffolding-the-skill)
4. [Step 3: Documenting the Skill Specification (SKILL.md)](#4-step-3-documenting-the-skill-specification-skillmd)
5. [Step 4: Writing Tests (Red Phase)](#5-step-4-writing-tests-red-phase)
6. [Step 5: Implementing Python CLI Logic (Green Phase)](#6-step-5-implementing-python-cli-logic-green-phase)
7. [Step 6: Registering the Skill in ClawBio](#7-step-6-registering-the-skill-in-clawbio)
8. [Step 7: Regenerating the Catalog and Linting](#8-step-7-regenerating-the-catalog-and-linting)
9. [Step 8: Updating routing in AGENTS.md](#9-step-8-updating-routing-in-agentsmd)
10. [Step 9: Git Workflow and Pull Request Checklist](#10-step-9-git-workflow-and-pull-request-checklist)

---

## 1. Overview & Architecture

Each skill in ClawBio is a self-contained module designed to perform a specific genomic or bioinformatics analysis task.

All processing is **local-first** — genetic data must never leave the machine. There should be no mandatory cloud data uploads, and external network calls should only be made to fetch public databases (e.g., PubMed, PDB, UniProt) without transmitting patient data.

Skills are composed of:
1. **`SKILL.md`**: The primary specification file containing YAML metadata and markdown sections that define the skill's capabilities, inputs, workflow, and gotchas.
2. **`Python Script`**: The execution script accepting `--input`, `--output`, and `--demo` arguments.
3. **`api.py`**: A clean, programmatic entry point for imports and orchestrator chaining.
4. **Demo Data**: Synthetic test data for testing and demo execution.
5. **Unit & Integration Tests**: Executable test cases validating CLI args, output files, and contract adherence.

---

## 2. Step 1: Planning and Naming Conventions

Before coding, plan your skill and ensure it adheres to these naming rules:
* **Folder Name**: Lowercase with hyphens (e.g., `pathway-enrichment`, `hla-typing`).
* **Python File Name**: Lowercase with underscores (e.g., `pathway_enrichment.py`, `hla_typing.py`).
* **Test File Name**: Prefix with `test_` and match the Python file name (e.g., `tests/test_pathway_enrichment.py`).
* **YAML Name**: Must match the folder name exactly.

**Scope Principle**: **One skill, one task.** If your skill tries to perform multiple unrelated jobs, split it into separate, distinct skills.

---

## 3. Step 2: Scaffolding the Skill

ClawBio provides an automated scaffolding script, [scaffold_skill.py](scaffold_skill.py), which generates a benchmark-ready skill directory structure pre-populated with conformant boilerplate files.

Run the scaffolder from the project root:
```bash
python scaffold_skill.py <your-skill-name> "<One-sentence description of the skill>"
```

### Example:
```bash
python scaffold_skill.py pathway-enrichment "Perform GO and KEGG pathway enrichment analysis from gene lists"
```

This creates the following file structure inside `skills/pathway-enrichment/`:
```
skills/pathway-enrichment/
├── SKILL.md                 # Conforming specification file (17-point check ready)
├── pathway_enrichment.py    # Main script skeleton with CLI/demo handling
├── api.py                   # Programmatic API wrapper
├── demo_input.txt           # Synthetic demo data file
├── tests/
│   └── test_pathway_enrichment.py   # Complete pytest test suite
└── bench_test_cases/        # Curated benchmark test cases (happy, empty, malformed)
```

---

## 4. Step 3: Documenting the Skill Specification (SKILL.md)

The `SKILL.md` file is the primary artifact. If a skill does not have a Python script, agents can still interpret `SKILL.md` to run the methodology. It must follow the exact structure defined in [templates/SKILL-TEMPLATE.md](templates/SKILL-TEMPLATE.md).

Every skill's `SKILL.md` must pass validation with `agentskills validate skills/<name>/`.

### YAML Frontmatter Details
Your YAML frontmatter at the top of `SKILL.md` must contain:
```yaml
---
name: pathway-enrichment
description: >-
  Perform GO and KEGG pathway enrichment analysis from gene lists.
version: 0.1.0
author: Your Name
domain: genomics
license: MIT

inputs:
  - name: input_file
    type: file
    format: [csv, txt]
    description: List of gene symbols
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: Analysis report
  - name: result
    type: file
    format: json
    description: Machine-readable results

dependencies:
  python: ">=3.11"
  packages:
    - pandas>=2.0

tags: [pathway, enrichment, go, kegg]

demo_data:
  - path: demo_input.txt
    description: Synthetic test gene list

endpoints:
  cli: python skills/pathway-enrichment/pathway_enrichment.py --input {input_file} --output {output_dir}

metadata:
  openclaw:
    requires:
      bins:
        - python3
    always: false
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]             # MUST use Node.js process.platform values: darwin, linux, win32
    trigger_keywords:
      - pathway enrichment
      - go enrichment
      - kegg analysis
---
```

> [!WARNING]
> **OS Platform Values**: Under `metadata.openclaw.os`, always use process.platform strings: `darwin`, `linux`, or `win32`. Do not use `macos`, `windows`, `ubuntu`, or `mac`.

### Required Markdown Sections
1. **`# 🦖 Skill Title`**: Contains the identity role of the skill.
2. **`## Trigger`**: Loud and explicit "Fire this skill when..." and "Do NOT fire when..." lists.
3. **`## Why This Exists`**: Contrast manual vs. automated workflows, and describe ClawBio's grounding.
4. **`## Core Capabilities`**: Enumeration of capabilities.
5. **`## Scope`**: Reconfirming "One skill, one task".
6. **`## Input Formats`**: Markdown table of supported input formats, extensions, required fields, and examples.
7. **`## Workflow`**: Numbered steps describing verification, processing, and output generation. Specify prescriptive vs. flexible freedom.
8. **`## CLI Reference`**: Standard, demo, and `clawbio.py` usage examples.
9. **`## Demo`**: Instructions on how to run demo mode and what output to expect.
10. **`## Algorithm / Methodology`**: Detailed scientific description. No hallucinated science; ground parameters in literature.
11. **`## Example Queries`**: Example user prompts routing to this skill.
12. **`## Example Output`**: An actual rendered example of what the report file produces.
13. **`## Output Structure`**: The directory tree of generated artifacts.
14. **`## Dependencies`**: Explicit lists of required and optional packages.
15. **`## Gotchas`**: **Minimum 3 gotchas** written in the pattern: *"The model will want to do X. Do not. Here is why."*
16. **`## Safety`**: Emphasize local-first execution and reference the medical disclaimer.
17. **`## Agent Boundary`**: Define what the LLM agent does (explaining/dispatching) vs. what the skill script does (executing).
18. **`## Integration with Bio Orchestrator`**: Trigger conditions and chaining partners.
19. **`## Maintenance`**: Monthly review cadence, deprecation, and staleness signals.
20. **`## Citations`**: Academic or database citations.

---

## 5. Step 4: Writing Tests (Red Phase)

ClawBio mandates **Red/Green Test-Driven Development (TDD)** for all skill changes.

### 1. Register the Test Path
Add the skill's test folder path to `testpaths` inside [pytest.ini](pytest.ini):
```ini
testpaths =
    ...
    skills/pathway-enrichment/tests
```

### 2. Write the Tests First
The scaffolded tests inside `skills/<name>/tests/test_<name_underscore>.py` check:
* No-arguments execution exits with non-zero.
* Demo mode runs successfully and outputs `report.md` and `result.json`.
* Input mode successfully produces reports from user input files.
* Missing input paths cause non-zero exit codes.
* Machine-readable `result.json` exists and is formatted correctly.
* Every generated markdown report contains the required medical disclaimer.
* The output structure matches the `SKILL.md` output contract.

### 3. Verify Test Failures (Red Phase)
Run the test command and confirm that the tests fail since the Python script is just a skeleton:
```bash
python -m pytest skills/pathway-enrichment/tests/ -v
```

---

## 6. Step 5: Implementing Python CLI Logic (Green Phase)

Write the implementation in `skills/<name>/<name_underscore>.py` to satisfy the tests.

### Code Style Guidelines
1. **Python**: Use Python 3.10+ types (e.g. use `X | None` instead of `Optional[X]`).
2. **Paths**: Always use `pathlib.Path` for filesystem operations. Use `Path(__file__).resolve().parent` to resolve paths relatively. Never hardcode absolute paths.
3. **CLI Arguments**: Use `argparse` with `--input`, `--output`, and `--demo` arguments.
4. **Warning on Overwrites**: Always check if output files exist before overwriting them.

### Disclaimer Requirement
Every markdown report generated by a ClawBio skill must end with the following disclaimer text verbatim:

> *"ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions."*

### Standard Output Files
* **Primary Report**: `report.md` detailing human-readable findings.
* **Structured Output**: `result.json` containing machine-readable results (e.g., `{"skill": "pathway-enrichment", "variants_processed": 5, "findings": [...]}`).

### Run the Tests (Green Phase)
Verify that all tests pass after your implementation is complete:
```bash
python -m pytest skills/pathway-enrichment/tests/ -v
```

---

## 7. Step 6: Registering the Skill in ClawBio

To hook your skill up to the CLI runner, you must register it in the `SKILLS` registry inside [clawbio.py](clawbio.py) around line 258:

```python
SKILLS = {
    ...
    "pathway-enrich": {
        "script": SKILLS_DIR / "pathway-enrichment" / "pathway_enrichment.py",
        "demo_args": [
            "--input",
            str(SKILLS_DIR / "pathway-enrichment" / "demo_input.txt"),
        ],
        "description": "Perform GO and KEGG pathway enrichment analysis from gene lists",
        "allowed_extra_flags": set(),  # Flags allowed to pass validation (prevents execution exploits)
        "accepts_genotypes": False,
    },
}
```

> [!IMPORTANT]
> **Allowed Extra Flags**: As a security feature, `clawbio.py` screens inputs and flags. Only arguments in `allowed_extra_flags` will be forwarded to the script when run via `clawbio.py`. Do not bypass this whitelist.

Verify the skill registration:
```bash
python clawbio.py list
```

---

## 8. Step 7: Regenerating the Catalog and Linting

To make the skill indexable by the platform, you must register its attributes and regenerate the machine-readable catalog.

### 1. Update [scripts/generate_catalog.py](scripts/generate_catalog.py)
Open `scripts/generate_catalog.py` and add entries for your skill in:
* **`FOLDER_TO_ALIAS`**: Map the folder to your `clawbio.py` CLI alias.
  ```python
  FOLDER_TO_ALIAS = {
      ...
      "pathway-enrichment": "pathway-enrich",
  }
  ```
* **`MVP_FOLDERS`**: Add your folder to mark the skill as fully functional.
  ```python
  MVP_FOLDERS = {
      ...
      "pathway-enrichment",
  }
  ```
* **`TRIGGER_KEYWORDS`**: Define orchestrator routing keywords.
  ```python
  TRIGGER_KEYWORDS = {
      ...
      "pathway-enrichment": ["pathway enrichment", "go enrichment", "kegg analysis", "gene ontology"],
  }
  ```
* **`CHAINING`**: Add chaining partners.
  ```python
  CHAINING = {
      ...
      "pathway-enrichment": ["profile-report"],
  }
  ```

### 2. Regenerate Catalog
Generate `skills/catalog.json`:
```bash
python scripts/generate_catalog.py
```

### 3. Run the Conformance Validator
Run the validator to verify that `SKILL.md` is discoverable and error-free:
```bash
agentskills validate skills/<your-skill-name>/
```

---

## 9. Step 8: Updating routing in AGENTS.md

Add your skill to the routing table and CLI guides inside [AGENTS.md](AGENTS.md) (which is the source of truth for routing tables and demo commands):
* Add a row to the **Skill Routing Table** under the appropriate categories.
* Add your CLI execution example to the **CLI Reference** section.
* Add your demo input under the **Demo Data** table.
* Add the demo command under the **Demo Commands** list.

---

## 10. Step 9: Git Workflow and Pull Request Checklist

1. Create a clean branch from `main`:
   ```bash
   git checkout -b feat/add-pathway-enrichment
   ```
2. Run the full test suite to make sure no regressions exist:
   ```bash
   python -m pytest -v
   ```
3. Commit and push your changes:
   ```bash
   git add skills/pathway-enrichment/ clawbio.py pytest.ini scripts/generate_catalog.py AGENTS.md
   git commit -m "Add pathway-enrichment skill with TDD, CLI, and registry registration"
   git push -u origin feat/add-pathway-enrichment
   ```
4. Open a Pull Request. Use [.github/PULL_REQUEST_TEMPLATE.md](.github/PULL_REQUEST_TEMPLATE.md) for the description, filling out every section fully. Include demo output snippets to help reviewers.
