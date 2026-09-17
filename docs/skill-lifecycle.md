# Skill lifecycle

How a skill enters ClawBio, how it is re-measured as models improve, and how it leaves. Version
1.0, 2026-09-17. This page is the editorial process; `GOVERNANCE.md` says who decides.

## Why a lifecycle, and why retirement is normal

A skill that encodes a public table has a shelf life measured in model generations. In September
2026 a general-purpose agent with database access and no skill reached 0.942 phenotype accuracy
on a 110-case pharmacogenomics benchmark, against 0.966 for the authored skill. The knowledge
inside a skill is a wasting asset. What does not waste is the contract around it: what the skill
refuses, what it records, whether it gives the same answer twice, and whether it runs where a
consumer product cannot. A library that only grows is a library nobody trusts, so retirement is
a designed outcome here, not a failure.

## Classes

Every skill is classified by what would be lost if it were deleted and a frontier model with
tool access did the task instead.

- **D1, knowledge-encoding.** The value is a mapping or threshold that already exists publicly
  and authoritatively elsewhere. Depreciates fastest.
- **D2a, local orchestration.** Drives local scientific binaries or cluster pipelines (nextflow,
  PLINK, samtools). Absorbed only when a harness reaches that infrastructure.
- **D2b, hosted orchestration.** Calls hosted services or imports libraries. Absorbed first by
  harnesses with built-in connectors.
- **D3, contract and policy.** Input validation, abstention, provenance, determinism, local
  governance. Survives the model getting smarter.
- **D0, too thin to judge.** A candidate for retirement or completion, not for shipping.

Modifiers: `verifiability` (ground-truth, partial, none), `consequence` (clinical, research,
exploratory), `locality` (portable, local-required). Clinical-consequence skills are hardened by
default and never retired on accuracy parity alone.

## Tiers and admission

Tiers follow the maintainers' decision of 2026-08-07: **research**, **benchmarked**, and
**clinical** (deferred until a regulator's decision exists). A skill is admitted at the tier its
evidence supports, never at the tier its author hopes for.

- **research**: a `SKILL.md` that validates, a script, tests, a demo, a declared code licence.
- **benchmarked**: everything above plus a ClawBench card with a DOI that an outsider can
  recompute in one command, and `benchmark_validated: true` in the catalogue set by that card.
- **clinical**: not available. The word is reserved.

## Re-measurement

On each frontier model release, and at least twice a year, the benchmarked skills are re-run
against a no-skill agent on the same cases, with three replicates and hashed case identifiers.
The endpoints are not accuracy alone:

- accuracy on well-formed inputs
- over-answer rate on ill-formed inputs (missing data, unreported locus, negation, contradiction)
- abstention stability across identical runs
- provenance completeness and re-runnability from inputs alone

Half of every perturbation suite is held back unpublished, so the next model cannot have read it.

## Decisions a re-measurement can produce

| Outcome | Meaning |
|---|---|
| **keep** | The skill still does something the no-skill agent does not, on the endpoints above |
| **harden** | Accuracy parity reached but the contract layer is missing; add validation and abstention before the next round |
| **demote** | Reference material only; catalogued, not shipped in the wheel |
| **retire** | Removed from the catalogue with a dated entry in `RETIRED.md` naming the evidence |
| **promote** | Evidence now supports a higher tier |

Rules written before any measurement, so they cannot be fitted to it: accuracy parity alone
never retires a skill; a skill with no ground truth can only be retired on recorded judgement;
clinical-consequence skills are hardened rather than retired; locality beats capability; any
skill that survives a round acquires a benchmark or is demoted.

## Editors

Named domain reviewers, listed in `MAINTAINERS.md`, review skill PRs in their area against the
classes above and recommend an outcome. Editors recommend; maintainers decide. There is no
steering group and no vote. An editor's recommendation is recorded alongside the decision so
that disagreement is visible later.

## Retirement, step by step

1. A re-measurement or an editor proposes retirement with the evidence.
2. The lead maintainer decides; the decision and the evidence are recorded in `RETIRED.md`.
3. The skill folder moves to `skills/RETIRED/<name>/` with its history intact. Public URLs to
   the old path keep resolving through a redirect note; nothing is deleted.
4. The catalogue regenerates without it. The release notes name it.

## What this page is not

Not a promise that any skill is clinically valid. Not a standard. Not a substitute for the
benchmark card, which is the only evidence that counts for the benchmarked tier.
