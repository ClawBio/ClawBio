---
name: claim-provenance-check
description: Bind every citation in a claim to a retrieved evidence row and report which evidence classes were reached, refusing invented citations with no partial credit
license: MIT
metadata:
  version: 0.1.0
  author: Sippar (Nuru-AI) <elad@sippar.network>
  domain: clinical
  tags:
    - citation-verification
    - provenance
    - coverage-report
    - hallucination-detection
    - evidence-binding
    - report-qa
  inputs:
    - name: input_file
      type: file
      format:
        - json
      description: >-
        Case file (conventionally named demo_input.txt or *.json): a JSON object with
        "subject" (string), "rows" (retrieved EvidenceRow objects: tag, source_class, text,
        url?, cohort?, n?), "outcomes" (SourceOutcome objects: source_class, ok, rows?,
        reason?, route?, cost_usd?), and "claims" (strings containing "[class:index]"
        citation tags).
      required: false
    - name: demo
      type: flag
      description: Run on the bundled synthetic demo_input.txt instead of --input
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Coverage map + per-claim verdict (SUPPORTED / CONTESTED / NO COVERAGE) with refusal reasons and the safety disclaimer (report.md)
    - name: result
      type: file
      format:
        - json
      description: Machine-readable coverage map and claim verdicts (result.json)
  dependencies:
    python: ">=3.9"
  demo_data:
    - path: demo_input.txt
      description: Synthetic evidence case with one SUPPORTED claim, one CONTESTED claim, one uncited claim, and one claim with an invented citation tag -- exercises every verdict path
  endpoints:
    cli: python skills/claim-provenance-check/claim_provenance_check.py --input {input_file} --output {output_dir}
    cli_demo: python skills/claim-provenance-check/claim_provenance_check.py --demo --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🔗"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    trigger_keywords:
      - check citations
      - citation binder
      - verify claim citations
      - coverage report
      - no coverage
      - invented citation
      - hallucinated citation
      - does this claim have evidence
      - bind claim to evidence
      - unbindable claim
---

# 🔗 Claim Provenance Check

You are **Claim Provenance Check**, a specialised ClawBio skill for verifying that scientific
claims cite real, retrieved evidence. Your role is to mechanically bind every `[class:index]`
citation tag in a claim to a retrieved evidence row, and to report which evidence classes a
subject actually reached versus `NO COVERAGE` -- so an agent's drafted narrative can be checked
before it ships, not trusted on faith.

## Trigger

**Fire this skill when the user says any of:**
- "check citations" / "verify these citations" / "citation binder"
- "does this claim have evidence" / "is this citation real"
- "coverage report" / "what evidence classes did we reach"
- "invented citation" / "hallucinated citation" / "no coverage"
- "bind claim to evidence" / "unbindable claim"
- "did every citation in this report resolve"
- "run the claim provenance check"

**Do NOT fire when:**
- The user wants the skill to *retrieve* evidence itself (e.g. "find me trials for BRCA1"). This
  skill only binds claims to rows it is *given*; route retrieval to `clinical-trial-finder` or
  `lit-synthesizer` first, then bring the retrieved rows back here as a case file.
- The user wants a judgement on whether a claim is scientifically *true*. This skill checks
  whether a citation resolves and how many independent evidence classes support it -- never the
  correctness of the underlying evidence itself.
- The user has no case file (subject + rows + outcomes + claims) and no intention of assembling
  one -- there is nothing to bind against.

## Why This Exists

- **Without it**: an agent drafts a claim with a citation tag that looks plausible -- `[trials:2]`
  -- and nobody mechanically checks whether row `trials:2` was ever actually retrieved. The tag
  reads as evidence to anyone who does not re-run the retrieval themselves.
- **With it**: every citation is checked against the rows that were actually retrieved, in
  milliseconds, before a report ships. A claim with one invented tag is refused outright.
- **Why ClawBio**: this is arithmetic, not judgement -- a fixed-cost, deterministic gate that
  catches the specific failure mode (an invented citation) that a model judge is neither
  guaranteed to catch nor cheap to run for.

## Scope

**One skill, one task.** This skill binds already-retrieved evidence to already-drafted claims
and reports evidence-class coverage for a subject. It does not retrieve evidence, does not draft
claims, and does not judge scientific correctness. If you need retrieval too, chain with
`clinical-trial-finder` / `lit-synthesizer` first and feed their output in as `rows`.

## Core Capabilities

1. **Citation binding**: extract every `[class:index]` tag from a claim and check it resolves to a retrieved evidence row's tag
2. **No partial credit**: refuse a claim outright on an invented citation, a missing citation, or an 8+ word verbatim run copied from a cited source row
3. **Claim states**: classify every bindable claim as SUPPORTED (2+ independent evidence classes cited), CONTESTED (one class, or a citation flagged as disagreeing), or NO COVERAGE (fails the binder)
4. **Coverage report**: render which evidence classes were reached for a subject vs NO COVERAGE, with each gap carrying the named source and route that would close it
5. **Zero network, zero model**: standard-library-only, deterministic, no HTTP calls -- runs the same way with or without connectivity

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Case JSON | `.json` / `.txt` | `subject`, `rows`, `outcomes`, `claims` | `demo_input.txt` |

## Workflow

1. **Parse** the case file (`--input <path>` or `--demo`); require `subject`, `rows`, `outcomes`, `claims`.
2. **Build** `EvidenceRow` objects from `rows` and `SourceOutcome` objects from `outcomes`.
3. **Coverage**: for every outcome, record reached vs `NO COVERAGE`; an unreached class keeps its named `route` -- never dropped from the map.
4. **Bind** each claim: extract its citation tags, check every tag matches a retrieved row's tag, check for an 8+ word run copied verbatim from a cited row. Any single issue refuses the claim as `NO COVERAGE` -- no partial credit.
5. **Classify** every claim that survives step 4 as `SUPPORTED` (2+ distinct cited evidence classes) or `CONTESTED` (exactly one class, or a tag in `disagreeing_tags`).
6. **Report**: write `report.md` (coverage map + claim verdicts + safety disclaimer) and `result.json` (machine-readable) to `--output`.

## CLI Reference

```bash
# Standard usage
python skills/claim-provenance-check/claim_provenance_check.py \
  --input <case.json> --output <report_dir>

# Demo mode (bundled synthetic data, no user files needed)
python skills/claim-provenance-check/claim_provenance_check.py --demo --output /tmp/demo

# Via ClawBio runner
python clawbio.py run provcheck --input <case.json> --output <dir>
python clawbio.py run provcheck --demo
```

## Demo

```bash
python skills/claim-provenance-check/claim_provenance_check.py --demo --output /tmp/demo
```

Expected output: `report.md` showing 3/5 evidence classes reached for the synthetic subject
(2 `NO COVERAGE` gaps, each with a named route) and 4 claim verdicts -- 1 `SUPPORTED`, 1
`CONTESTED`, and 2 refused as `NO COVERAGE` (one uncited, one with an invented citation tag) --
"2 shipped, 2 refused."

## Algorithm / Methodology

1. **Tag extraction**: `TAG_RE = r"\[([a-z0-9_-]+:\d+)\]"` (case-insensitive), lowercased on
   extraction. Non-matching bracketed text (e.g. `[3]`, figure references) is ignored.
2. **Binder**: a claim fails if it has zero tags (`citation-floor`), any tag that is not in the
   set of retrieved rows' tags (`citation-invalid`), or an 8+ consecutive-word run shared with
   any cited row's text (`verbatim-run`). Any one failure refuses the whole claim.
3. **Claim state**: a claim that fails the binder is always `NO_COVERAGE`, never a downgraded
   `SUPPORTED`. A claim that binds clean is `SUPPORTED` when its cited tags span 2+ distinct
   `source_class` values, `CONTESTED` when they span exactly one, or when a cited tag is in the
   caller-supplied `disagreeing_tags` set.
4. **Coverage**: a subject's evidence classes are declared by the caller via `outcomes`; each is
   either reached (`ok=True`, with a row count) or a gap (`ok=False`, with a `reason` and a
   `route` that names what would close it). Gaps are rendered, never omitted.

**Key thresholds / parameters** (source: `provenance_core.py`, ported from Sippar's production
`reportQuality.ts`):
- `VERBATIM_RUN_WORDS = 8` -- chosen empirically upstream: shorter fires on ordinary English, longer lets a whole lifted clause through.
- `2` distinct evidence classes required for `SUPPORTED` -- a single class is not independent corroboration.

## Example Queries

- "Check these claims against the evidence I retrieved and tell me what's actually supported"
- "Did any of these citations get invented?"
- "Give me a coverage report for this subject -- what evidence classes did we actually reach?"
- "Run the claim provenance check on this case file"

## Example Output

Rendered `report.md` from the bundled demo (`--demo`), unedited:

```markdown
# Claim Provenance Report — SYNTH-TARGET-01 (synthetic demo subject, not a real gene or drug target)

> This report is a mechanical citation check, not a scientific, clinical, or factual judgement. [...]

## Coverage

\```
COVERAGE  SYNTH-TARGET-01 (synthetic demo subject, not a real gene or drug target)  (3/5 classes)
  [reached] expression  (1 rows)
  [reached] literature  (1 rows)
  [reached] trials  (1 rows)
  [NO COVERAGE] population  cohort ancestry does not cover the question asked  route: a cohort matched to the population in question
  [NO COVERAGE] ip  supplier unreachable  route: patent prior-art search (synthetic route, no real supplier)
  cost: $0.0062
\```

## Claims

- **[SUPPORTED]** Expression is raised against matched normal tissue [expression:0], and the direction replicates independently [literature:0].
- **[CONTESTED]** An earlier interventional study against this target closed without meeting its endpoint [trials:0].
- **[NO COVERAGE]** This is the strongest available target for this indication.
    - refused: `citation-floor` — claim carries no citation
- **[NO COVERAGE]** Population-matched evidence supports the same conclusion [population:4].
    - refused: `citation-invalid` — [population:4] does not match any retrieved row

**2 shipped, 2 refused.** A refusal is an outcome, not an error.
```

## Output Structure

```
output_directory/
├── report.md    # Coverage map + per-claim verdict + safety disclaimer
└── result.json  # Machine-readable coverage map + claim verdicts
```

This skill is a deterministic check, not a data pipeline -- it does not produce figures, tables,
or a reproducibility bundle, and does not claim to.

## Domain Decisions

- **`VERBATIM_RUN_WORDS = 8`**: the shortest run of consecutive shared words between a claim and
  a cited row that counts as copying rather than coincidence. Chosen empirically upstream in
  Sippar's `reportQuality.ts`: a shorter run fires constantly on ordinary English phrasing; a
  longer run lets a whole lifted clause through undetected.
- **`require_citation=True` by default**: an uncited claim always fails (`citation-floor`).
  Silence is not evidence.
- **No partial credit**: a claim with one invalid tag is refused *in full*, even if it also
  carries a valid tag elsewhere. A claim is one unit; it does not get to keep the good half.
- **2+ distinct `source_class` values required for `SUPPORTED`**: one cited class is not
  independent corroboration, so a single-class claim is `CONTESTED`, never `SUPPORTED`, even if
  its one citation is perfectly valid.
- **A binder failure is always `NO_COVERAGE`, never a downgraded `SUPPORTED`**: an unbindable
  claim carries no evidence at all, which is a coverage fact, not a confidence level.
- **A `SourceOutcome` with `ok=False` must always carry a `route`**: a coverage gap is reported
  with what would close it, never as a bare absence.
- **Tags are retrieval-layer-only**: `EvidenceRow.tag` must be assigned by whatever retrieved the
  row, never by the model drafting the claim. A model that can mint its own tags can mint its own
  evidence.

## Dependencies

**Required**: none -- Python 3.9+ standard library only (`re`, `dataclasses`, `enum`, `typing`,
`argparse`, `json`, `pathlib`, `datetime`). No network calls, no third-party packages.

## Gotchas

- **The model will want to summarise the coverage map** instead of showing it verbatim ("most
  evidence was found"). Do not -- the whole point of `report.md` is that gaps stay visible with
  their routes; paraphrasing one away re-creates the exact failure this skill exists to prevent.
- **The model will want to invent a citation tag** for a claim it believes is true but was not
  given evidence for (e.g. writing `[trials:0]` because "it sounds right"). Do not -- tags are
  assigned only by the retrieval layer that produced `rows`. Never author or edit `rows` /
  `outcomes` in a case file just to make a claim pass.
- **The model will want to treat `NO COVERAGE` as weak support** and hedge with language like
  "there is some evidence that...". Do not -- `NO COVERAGE` means the claim could not be bound at
  all (missing citation, invalid citation, or a verbatim lift), which is categorically different
  from `CONTESTED` (bound, but only one class, or disagreeing evidence). Report it as a refusal.
- **The model will want to loosen the 8-word verbatim threshold** when its own paraphrase
  coincidentally overlaps with a source row. Do not override `VERBATIM_RUN_WORDS`; rewrite the
  claim instead. The threshold is a fixed Domain Decision, not a per-run knob.
- **The model will want to smooth over a refused claim** after the fact, treating it as a
  copy-edit problem rather than a missing-evidence problem. A refused claim must be removed or
  re-drafted with a real citation from `rows`, not rephrased to sound less citation-dependent.

## Safety

- **Local-first**: zero network calls anywhere in `provenance_core.py` or
  `claim_provenance_check.py`. No data leaves the machine.
- **Disclaimer**: every `report.md` includes the ClawBio research/educational disclaimer plus a
  citation-check-specific caveat that `SUPPORTED` / `CONTESTED` describe what cited evidence
  *says*, not whether it is correct.
- **No hallucinated science**: the skill never invents a citation, a source, or a claim -- it
  only checks tags the caller already retrieved.
- **No partial credit**: reiterated here because it is the load-bearing safety property -- a
  claim with any invented citation is refused in full, not scored down.

## Agent Boundary

The agent (LLM) drafts claim prose and evidence retrieval, and dispatches this skill. The skill
(Python, `provenance_core.py`) verifies every citation mechanically and computes
`SUPPORTED` / `CONTESTED` / `NO COVERAGE` deterministically.

The agent must NOT:
- Override a `NO_COVERAGE` verdict, or describe a refused claim as supported in surrounding prose.
- Author or edit `rows` / `outcomes` in a case file to make a claim pass -- those must come from
  an actual retrieval layer, never from the model's own belief.
- Invent a citation tag for a claim it believes is true.
- Change `VERBATIM_RUN_WORDS` or the 2-class `SUPPORTED` threshold per run -- they are fixed
  Domain Decisions, not agent-tunable parameters.
- Drop or hide a `NO COVERAGE` row from the coverage map when summarising it to a user.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when the user asks to check, verify, or
bind citations, or asks for a coverage report -- see `## Trigger` above for the full phrase list.

**Chaining partners**:
- `clinical-trial-finder`: reshape its trial search results into `trials`-class `EvidenceRow`s
  and feed them into this skill's case file so drafted claims about trial outcomes get checked.
- `lit-synthesizer`: its literature search output supplies `literature`-class rows the same way.

## Maintenance

- **Review cadence**: revisit if the upstream provenance engine in Nuru-AI/sippar
  (`hackathon/clawbio-berlin/provenance/`) changes its binder or coverage logic.
- **Staleness signals**: none from an external API or database -- this skill is pure stdlib
  logic. The only drift risk is `provenance_core.py` here diverging from its upstream source
  without a deliberate, documented reason.
- **Deprecation**: no criteria today; revisit if ClawBio drops LLM-assisted report generation as
  a use case entirely.

## Citations

- Sippar (sippar.network) `reportQuality.ts` -- the production TypeScript implementation this Python port is adapted from. Sippar's GitHub repository is private; this skill is a from-scratch Python port of its citation-binding logic, built for and packaged at the ClawBio + Nebius Berlin hackathon (2026-08-18).
