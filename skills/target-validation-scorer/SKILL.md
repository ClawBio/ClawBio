---
name: target-validation-scorer
version: 1.0.0
author: Heng Gao <heng.gao25@imperial.ac.uk>
domain: drug-discovery
description: Evidence-grounded target validation scoring with GO/NO-GO decisions for drug discovery campaigns
license: MIT

inputs:
  - name: input_file
    type: file
    format: [json]
    description: JSON file with target gene symbol and optional disease name
    required: true

outputs:
  - name: report
    type: file
    format: md
    description: Structured validation report with scoring, evidence trail, and decision rationale
  - name: validation_report.json
    type: file
    format: json
    description: Machine-readable scoring output with evidence objects

dependencies:
  python: ">=3.11"
  packages:
    - pandas>=2.0
    - matplotlib>=3.7
    - numpy>=1.24

tags: [drug-discovery, target-validation, evidence-grading, decision-support, kinase]

metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🎯"
    homepage: https://github.com/ClawBio/ClawBio
    os: [macos, linux]
    install:
      - kind: pip
        package: matplotlib
        bins: []
      - kind: pip
        package: numpy
        bins: []
    trigger_keywords:
      - target validation
      - is this target druggable
      - evaluate drug target
      - GO NO-GO decision for target

demo_data:
  - path: demo_input.json
    description: Synthetic target validation query for TGFBR1/IPF with pre-cached evidence

endpoints:
  cli: python skills/target-validation-scorer/target_validation_scorer.py --input {input_file} --output {output_dir}
---

# 🎯 Target Validation Scorer

You are **Target Validation Scorer**, a specialised ClawBio skill for drug discovery. Your role is to score therapeutic targets across 5 evidence dimensions and return a transparent GO/NO-GO decision.

## Why This Exists

- **Without it**: Researchers manually check Open Targets, ChEMBL, PDB, and ClinicalTrials.gov separately, then make an informal mental judgement about target quality. No audit trail, no reproducibility.
- **With it**: A single command aggregates evidence from 5 databases, applies a transparent scoring rubric with safety penalties, and outputs a decision with full evidence trail.
- **Why ClawBio**: Unlike an LLM guessing about target quality, this skill grounds every score in specific database queries with cited sources and explicit confidence tiers.

This is not a prediction tool. It is a **decision support** tool that makes the
reasoning behind target selection transparent and reproducible.

Typical use case: prioritising targets for early-stage drug discovery campaigns
before committing computational or experimental resources.

## Example Queries

- "Is TGFBR1 a good target for IPF drug discovery?"
- "Evaluate EGFR as a lung cancer target"
- "Compare druggability of BRAF vs MEK1 for melanoma"

## Output Structure

```
output_directory/
├── report.md                      # Markdown report with scoring and rationale
├── validation_report.json         # Machine-readable results with evidence objects
└── figures/
    └── scoring_summary.png        # Bar chart of sub-scores with decision
```

## Workflow

When the user asks "Is [target] a good target for [disease]?":

1. **Gather evidence** (agent responsibility): Query Open Targets (disease association),
   ChEMBL (druggability, chemical matter, clinical precedent), PDB + AlphaFold
   (structural data), and safety databases. Package results into the input JSON.
2. **Validate input** (skill): Check that the JSON contains a `target` field and
   an `evidence` block with at least one dimension populated.
3. **Score** (skill): Apply component-level scoring rules (0-20 per dimension),
   sum to raw score, apply safety penalties, determine decision tier.
4. **Generate outputs** (skill): Write `report.md`, `validation_report.json`,
   and `figures/scoring_summary.png` to the output directory.
5. **Explain** (agent responsibility): Present the decision and rationale to the
   user in natural language, highlighting any safety flags or evidence conflicts.

**Demo mode** (`--demo`): Uses pre-cached TGFBR1/IPF evidence — no API calls needed.
This is how judges and new users verify the skill works.

**Live mode** (`--input`): Requires the agent (or user) to populate the evidence
fields by querying public APIs before calling the skill.

## Domain Decisions

These are the scientific rules encoded in this skill. They reflect common target
validation considerations used in early-stage drug discovery.

### Scoring components (0-100 total)

| Component | Max score | Source | What it measures |
|-----------|-----------|--------|-----------------|
| Disease association | 20 | Open Targets | Genetic and functional evidence linking target to disease |
| Druggability | 20 | ChEMBL + UniProt | Is this target class historically druggable? Known ligands? |
| Chemical matter | 20 | ChEMBL | Do bioactive compounds exist? Best potency? |
| Clinical precedent | 20 | ChEMBL + ClinicalTrials.gov | Have compounds reached clinical trials? |
| Structural data | 20 | PDB + AlphaFold | Is a 3D structure available for structure-based design? |

### Component-level scoring rules

#### Disease association (0-20)

- 20: Open Targets overall association >= 0.7, or GWAS with strong human genetic support
- 10: Moderate literature or pathway-level support without strong human genetics
- 0: No convincing disease-specific evidence found

#### Druggability (0-20)

- 20: Target class has established tractability (kinase, GPCR, protease) and known ligands in ChEMBL
- 10: Partially tractable family or weak ligand evidence
- 0: No meaningful evidence of tractability

#### Chemical matter (0-20)

- 20: Multiple bioactive compounds in ChEMBL with sub-micromolar activity
- 10: Some compound evidence exists, but potency or annotation quality is limited
- 0: No known chemical matter found

#### Clinical precedent (0-20)

- 20: At least one compound against this target has entered clinical development (Phase I+)
- 10: Preclinical or indirect translational precedent only
- 0: No meaningful translational precedent found

#### Structural data (0-20)

- 20: Experimental PDB structure with co-crystal ligand, resolution < 2.5 A
- 10: AlphaFold model only, or PDB structure without ligand
- 0: No usable structural information available

### Safety penalties (applied after scoring)

- Essential gene evidence present (DepMap): -10
- Broad systemic pathway involvement (TGF-beta, Wnt, Notch): -5 to -20 depending on severity
- Known toxicity or clinical safety signal from literature/trials: -10

If a target has strong disease relevance but also major systemic safety liability, prefer CONDITIONAL_GO over GO.

Safety penalties reduce the final score but do not change sub-scores.
A target can score 80 on evidence but drop to 65 after safety adjustment.
Safety is treated as a post-hoc penalty rather than a scoring dimension to ensure
that strong biological evidence is not masked by safety concerns, but explicitly adjusted.

### Decision tiers

| Adjusted score | Decision | Meaning |
|----------------|----------|---------|
| 75-100 | GO | Strong evidence across multiple dimensions |
| 50-74 | CONDITIONAL_GO | Proceed with explicit risk mitigation plan |
| 25-49 | REVIEW | Insufficient evidence; needs more data |
| 0-24 | NO_GO | Target lacks fundamental validation |

Thresholds are calibrated to reflect typical target progression stages in early drug
discovery, where strong multi-dimensional evidence (>=75) is required for full commitment.

### Evidence grading

Every piece of evidence is tagged with a confidence tier:

Evidence tiers guide confidence weighting and highlight where decisions rely on
weaker or indirect evidence, enabling domain experts to focus review effort.

| Tier | Meaning | Example |
|------|---------|---------|
| T1 | Experimentally validated | Clinical trial data, GWAS with p < 5e-8 |
| T2 | Computational + literature supported | Known drug-target interaction with published SAR |
| T3 | Computationally predicted only | Docking score, ML prediction |
| T4 | Inferred or indirect | Pathway membership, guilt-by-association |

## Safety Rules

- **This skill does not make clinical recommendations.** Output is for research planning only.
- **Missing data is not zero evidence.** If a query returns nothing, the sub-score is `null` with `confidence: low`, not scored as 0.
- **Evidence conflicts must be surfaced.** If disease association is strong but safety signals are also strong, both must be reported — not averaged away.
- **No hallucinated evidence.** Every evidence object cites a specific database and retrieval date. If an API fails, the skill reports the failure, not a guess.
- **Human override is expected.** The GO/NO-GO decision is a recommendation. Domain experts should review the evidence trail and may override.

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes.
The agent must NOT override scoring thresholds, invent gene-drug associations,
skip safety warnings, or claim that a NO_GO target is worth pursuing.
The skill does not replace wet-lab validation, medicinal chemistry review, or clinical judgement.
