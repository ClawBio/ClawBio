---
name: completeness-aware-caller
description: >-
  Three-state gene-presence calls (present / absent / cannot-conclude) and
  ANI-ranking gates for incomplete genomes and MAGs. Refuses to call a gene
  absent, or to rank references by ANI, when assembly completeness cannot
  support the claim.
license: MIT
metadata:
  version: "0.1.0"
  author: Qihao Duann
  domain: genomics
  tags:
    - comparative-genomics
    - completeness
    - metagenomics
    - mag
    - ani
    - abstention
  inputs:
    - name: completeness
      type: value
      format:
        - float
      description: >-
        Assembly completeness as fraction (0-1] or percent (1-100], or via
        --busco-json pointing at a busco-assessor result.json (required
        unless --demo)
      required: true
    - name: genes
      type: file
      format:
        - json
      description: 'Gene search results: [{"gene": "nifH", "found": false}, ...]'
      required: false
  outputs:
    - name: report
      type: file
      format: md
      description: Three-state call table with rationale and refusal messages
    - name: result
      type: file
      format: json
      description: Machine-readable calls, gate decision, and thresholds used
  dependencies:
    python: ">=3.10"
    packages: []
  demo_data:
    - path: examples/demo_genes.json
      description: STM815 frag70 gene set (nod cluster undetectable; nifH survives on a 2nd copy)
  endpoints:
    cli: python skills/completeness-aware-caller/completeness_aware_caller.py --busco-json {busco_result} --genes {genes_json} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🚦"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    trigger_keywords:
      - "is this gene really absent"
      - "completeness-aware"
      - "can I trust this MAG"
      - "gene absence incomplete genome"
      - "cannot conclude"
      - "abstention genomics"
      - "MAG interpretation"
      - "ANI confidence"
---

# 🚦 Completeness-Aware Caller

You are the **completeness-aware-caller**, a ClawBio skill that separates
taxonomic confidence from functional confidence in incomplete genomes. A
fragmented genome makes missing genes look absent and close relatives look
functionally identical; this skill's deliverable is the refusal — an explicit
`CANNOT CONCLUDE` printed as a first-class result, not a footnote.

## Trigger

**Fire when the user says any of:**
- "is this gene really absent", "gene absent or just missing"
- "can I trust this MAG", "MAG completeness interpretation"
- "completeness-aware call", "abstention", "cannot conclude"
- "is the ANI difference meaningful", "ANI confidence"
- "does incompleteness affect my conclusion"

**Do NOT fire when:**
- User wants to *measure* completeness → route to `busco-assessor`
- User wants to *compute* ANI or relatedness → run FastANI (or
  `galaxy-bridge`), then come back here with the values
- User wants variant-level interpretation → `clinical-variant-reporter`
- User wants microbiome community profiling → `claw-metagenomics`

## Why This Exists

- **Without it**: pipelines print "gene absent" from 70%-complete MAGs and
  rank references by ANI margins smaller than the noise incompleteness
  induces. Absence of evidence is silently converted into evidence of
  absence.
- **With it**: every absence claim and every ANI ranking carries either a
  confidence grounded in measured completeness, or an explicit refusal with
  an actionable message.
- **Why ClawBio**: chains directly after `busco-assessor` (reads its
  `result.json` verbatim) and complements `galaxy-bridge` ANI workflows.

## Scope

One skill, one task: **decide whether completeness supports a claim**. This
skill does not measure completeness, compute ANI, search for genes, or
assemble genomes — and it does not assess contamination: a `present` call
assumes the contig truly belongs to the genome (Gotcha 5), so run CheckM
upstream when that is in doubt. It consumes those results and gates the
conclusions.

## Domain Decisions

1. **Absence needs ≥95% completeness** (`MIN_COMPLETENESS_FOR_ABSENT =
   0.95`). Under a random-loss model, P(gene missed | truly present) ≈
   1 − completeness; 0.95 caps that at 5%. Anchored to MIMAG quality tiers
   (Bowers et al. 2017) where "high-quality draft" starts at >90%.
2. **ANI drift model**: drift = 0.82 × (1 − completeness) ANI points,
   calibrated on the *P. phymatum* STM815 degradation benchmark (0.41 points
   of observed drift at 50% retention against *P. xenovorans* LB400, 0.44
   against *B. cenocepacia* J2315). This coefficient is CALIBRATED on that
   benchmark, not validated against it: the same four points produced the
   number and would be used to check it. In that benchmark the reference
   ranking never actually failed — the margin widens under degradation — so
   the gate bounds noise rather than catching a demonstrated error.
3. **Safety factor 2**: a ranking margin must exceed 2× modelled drift —
   signal must beat twice the noise.
4. **Species-boundary zone 94–96%** (Jain et al. 2018): ANI values in this
   window with completeness <90% flag species assignment as uncertain.

## Workflow

1. Parse completeness from `--completeness` (fraction or percent) or from a
   `busco-assessor` `result.json` via `--busco-json` (field `C`). Reject
   values outside (0, 1] / (1, 100].
2. For each gene in `--genes`: found → `present`; not found and
   completeness ≥ 0.95 → `absent` with confidence = completeness; otherwise
   → `cannot_conclude` with the miss probability spelled out.
3. If `--ani-a`/`--ani-b` supplied: compute margin, modelled drift, and
   decide `rank` vs `cannot_conclude`; flag the species-boundary zone.
4. Write `report.md` (call table + gate + disclaimer), `result.json`
   (calls, thresholds, decision), `commands.sh` (replay).

## CLI Reference

```bash
# From a busco-assessor run + gene search results + FastANI values
python skills/completeness-aware-caller/completeness_aware_caller.py \
  --busco-json /path/to/busco_out/result.json \
  --genes genes.json --ani-a 81.45 --ani-b 80.75 --output /tmp/cac_out

# Direct completeness value
python skills/completeness-aware-caller/completeness_aware_caller.py \
  --completeness 0.98 --genes genes.json --output /tmp/cac_out

# Demo: STM815 degradation benchmark at 70% completeness
python skills/completeness-aware-caller/completeness_aware_caller.py \
  --demo --output /tmp/cac_demo
```

## Example Output

```markdown
# Completeness-Aware Caller Report

**Assembly completeness**: 68.1%

| Gene | Status | Confidence | Rationale |
|------|--------|-----------|-----------|
| nifH | CANNOT CONCLUDE | 0.68 | CANNOT CONCLUDE absence: the assembly is only 68% complete, so a truly present gene would be missed with ~32% probability. |
| nifD | PRESENT | 1.00 | Detected in the assembly; detection is positive evidence. |

## ANI ranking gate
**Decision**: RANK
Margin 0.70 ANI points exceeds 2x modelled drift (0.26); ranking is supported.
```

(Real run: BUSCO C:68.1% on STM815 degraded to 70% retention — where the
naive call for nifH would have been a false "absent".)

## Gotchas

1. **The model will want to treat "not found" as "absent".** Do not. The
   whole point of this skill is that at 70% completeness a missing gene is
   missing with ~30% probability even when truly present. Always route
   "not found" through `call_gene`, never report absence directly.
2. **The model will want to reuse the drift coefficient across clades.**
   The default 0.82 was calibrated on one Burkholderiaceae benchmark. For
   distant clades or eukaryotes, recalibrate (degrade a complete genome of
   the target clade and measure) or say the default is a Burkholderiaceae
   calibration.
3. **BUSCO C% is itself noisy at high retention.** In the calibration
   benchmark, dropping 9.8% of bases produced M=19.8% at frag90 — but at
   frag70 and frag50 base loss and marker loss track closely. Marker loss is
   noisy at low loss rather than clustered throughout. Treat completeness as
   an estimate with variance, which is why the thresholds are conservative.
4. **The model will want to treat 1 − completeness as a real miss
   probability. It is not.** That figure assumes genes are lost
   independently. In the calibration benchmark they are not: every gene that
   vanished at 70% retention sat in a single 50 kb window, so a co-located
   operon disappears as a block. For a clustered gene set the true miss risk
   is higher than 1 − C; for a dispersed set it is lower. Report the number
   as a modelled bound, never as a measured rate.

5. **`present` is not contamination-safe.** Detection confidence 1.0 assumes
   the contig truly belongs to the genome. For MAGs with high contamination,
   check CheckM contamination before trusting `present` calls.

## Safety

- All processing is local; no sequence data leaves the machine.
- Every `report.md` ends with the ClawBio disclaimer: *"ClawBio is a
  research and educational tool. It is not a medical device and does not
  provide clinical diagnoses. Consult a healthcare professional before
  making any medical decisions."*
- The skill never invents completeness or ANI values; it only gates numbers
  produced by upstream tools.

## Agent Boundary

- **Agent dispatches**: completeness (or busco-assessor result.json path),
  gene search results, optional ANI pair.
- **Skill executes**: threshold logic, three-state calls, gate decision,
  report writing.
- **Agent explains**: what the refusal means for the user's question and
  what evidence would unlock a conclusion.
- **Agent must NOT**: override a `cannot_conclude` into a definitive call,
  tweak thresholds ad hoc, or report absence for genes the skill refused.

## Chaining Partners

| Upstream | Handoff | Downstream |
|----------|---------|-----------|
| `busco-assessor` | `result.json` (C%) | `completeness-aware-caller` |
| FastANI / `galaxy-bridge` | ANI values | `completeness-aware-caller` |
| `completeness-aware-caller` | `result.json` three-state calls | `profile-report`, `lit-synthesizer` |

## Maintenance

- **Review cadence**: quarterly, or when BUSCO/MIMAG standards revise.
- **Staleness signals**: busco-assessor changes its `result.json` schema;
  new consensus thresholds for MAG quality; drift recalibrations published.
- **Deprecation path**: fold into a broader confidence-gating skill if one
  emerges; keep the three-state vocabulary stable.

## Citations

- Bowers R.M. et al. (2017). Minimum information about a single amplified
  genome (MISAG) and a metagenome-assembled genome (MIMAG). *Nature
  Biotechnology*. https://doi.org/10.1038/nbt.3893
- Jain C. et al. (2018). High throughput ANI analysis of 90K prokaryotic
  genomes reveals clear species boundaries. *Nature Communications*.
  https://doi.org/10.1038/s41467-018-07641-9
- Manni M. et al. (2021). BUSCO Update. *Molecular Biology and Evolution*.
  https://doi.org/10.1093/molbev/msab199
