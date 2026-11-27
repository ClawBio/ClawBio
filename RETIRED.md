# Retired skills

The record referenced by [docs/skill-lifecycle.md](docs/skill-lifecycle.md): every skill
removed from the catalogue, with the date, the evidence and who decided. Nothing is deleted —
a retired skill's folder moves to `skills/RETIRED/<name>/` with its history intact.

| Skill | Date | Evidence | Decided by |
|---|---|---|---|
| `skill-builder` | 2026-11-27 | Contributor tooling, not a bioinformatics skill, so it was counted and routed as one. Its SKILL.md generator does not read `templates/SKILL-TEMPLATE.md`, its checker tests 13 sections that differ from the AGENTS.md checklist, and it writes `skills/catalog.json` itself. `/new-skill` and the template cover the same job. [#551](https://github.com/ClawBio/ClawBio/issues/551) | Proposed by Cameron Lloyd; lead maintainer approval on the PR |
