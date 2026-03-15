---
name: list-skills
description: List all available ClawBio bioinformatics skills with their status and capabilities
---

List all available ClawBio skills by reading `skills/catalog.json`. For each skill, show:

- Name and CLI alias
- Description (one line)
- Status (MVP or Planned)
- Whether it has a Python script, tests, and demo data

Format as a clean table. If `catalog.json` is not found, run `python scripts/generate_catalog.py` first, then read the generated file.
