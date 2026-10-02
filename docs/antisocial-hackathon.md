# Antisocial — Hackathon Project Note

Team: **Antisocial**. Hackathon repo bootstrap only — no scientific skill code lands
in this commit. This note is the shared reference point for the team as
implementation starts.

## Scientific premise

Explore novel antimicrobial repurposing hypotheses, with a focus on resistant and
severe infection contexts. The angle is comparative: find drug targets that are
conserved across pathogens but diverge from the human host (so hitting them is less
likely to harm the patient), then connect those targets back to chemistry that
already exists — approved drugs, tool compounds, or agrochemicals — that could be
repurposed against them.

## Workflow concept

```
scientific activities
  -> narrow, reusable ClawBio skills   (one job each, testable in isolation)
  -> deterministic composition/workflow (fixed pipeline: skill A's output feeds skill B)
  -> bio-orchestrator                   (adaptive investigation: decides which skills
                                          to run and in what order for a given question)
```

The deterministic layer and the orchestrator are doing different jobs and should not
be blurred:

- **Deterministic skills/workflows**: given the same input, always produce the same
  output. This is where the actual comparative-genomics and cheminformatics logic
  lives — target conservation scoring, ortholog mapping, evidence aggregation. No
  judgment calls, no branching on "what seems interesting."
- **Agent reasoning (bio-orchestrator)**: decides *which* deterministic skills to run,
  in what order, and how to interpret a result well enough to decide the next step.
  It does not itself compute scientific results — it composes the skills that do.

## Intended reusable skill boundaries (not yet built)

Two new skills, scoped narrow on purpose so each is independently testable and
reusable outside this project:

1. **Comparative antimicrobial target discovery** — given a pathogen (or a resistant
   strain), find targets that are conserved across pathogens but divergent in the
   human host, optionally cross-checked against plant/agrochemical targets as
   additional supporting evidence.
2. **Target-to-chemistry evidence mapping** — given a target, return the graded
   chemistry evidence that exists against it: exact-target hits, ortholog hits,
   target-family hits, assay-level evidence, and chemotype-level evidence.

Each should do one job well and compose through the bio-orchestrator rather than
reimplementing pieces of each other.

## Expected to reuse, not rebuild

Resistance profiling, reference retrieval, literature work, and general orchestration
already have a home in the existing ClawBio skill catalog (see the Skill Routing
Table in `AGENTS.md` and `skills/catalog.json`). The two skills above are the actual
gap; everything else should be composed from what already exists.

## Current hackathon outputs

- This bootstrap: a clean `antisocial` branch off `main`, this project note, and the
  team's local dev conventions confirmed against `AGENTS.md` / `CONTRIBUTING.md`.
- No scientific results, code, or tests yet — implementation starts from here.
- Reproducibility and repo-contributable code, once implementation starts, follow the
  same bar as any other ClawBio skill: `templates/SKILL-TEMPLATE.md` structure,
  demo data + `--demo` support, tests under `skills/<name>/tests/`, and
  `python scripts/generate_catalog.py` after any `SKILL.md` frontmatter change.

## Current work / next likely implementation steps

1. Scaffold the comparative-target-discovery skill with
   `python scaffold_skill.py <name> "<description>"`, following
   `CONTRIBUTING.md` → "How to Contribute a Skill".
2. Scaffold the target-to-chemistry evidence-mapping skill the same way.
3. Wire both into `clawbio.py` for CLI registration once each has a working
   `--demo` path, then run `python clawbio.py list` and `make lint` to confirm.
4. Identify the specific existing skills to call for resistance profiling,
   reference retrieval, and literature lookup (via the bio-orchestrator's routing
   table) rather than duplicating that logic here.
