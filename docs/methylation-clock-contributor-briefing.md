# Methylation clock: contributor briefing

**Meeting:** Monday, 19 October 2026  
**Suggested slot:** 5 minutes explanation, 2 minutes demo, 5 minutes decisions

## One sentence

The skill takes a local table of methylation beta values, prepares it for PyAging, runs selected clocks, and writes predictions, input checks, model metadata, and a reproducibility bundle. It is a research tool; the outputs need clock- and tissue-specific interpretation.

## Suggested five-minute explanation

**1. Why it exists (45 seconds).** Researchers often have array methylation data but need a reproducible path from a table to a clock prediction. This skill gives ClawBio agents a discoverable `SKILL.md` and a CLI implementation. Agents can find it from requests such as “estimate epigenetic age from this methylation file,” then inspect the input and run the script.

**2. What goes in (45 seconds).** One input is required: a local CSV, TSV, compressed version, trusted pickle, or a public GEO accession. Rows are samples; columns include `cg`/`ch`/`rs`/`nv` probes with numeric beta values from 0 to 1. A `sample_id` column is optional. Nonprobe categorical columns belong in `--metadata-cols`.

**3. What it does (75 seconds).** The script derives a numeric `female` covariate when sex labels are present, aggregates EPICv2 duplicate probe IDs, converts the frame to AnnData, and calls PyAging. The draft fix keeps metadata such as `tissue_type` out of the aggregation call. The current improvement branch checks values and identifiers before inference and records basic input QC.

**4. What comes out (60 seconds).** `report.md` points to predictions, input QC, missing-feature diagnostics, clock metadata, and figures. `result.json` has structured findings; `reproducibility/` holds a replay command, environment recipe, and checksums. The report now shows each clock's unit and lists a correlation figure only when one was made.

**5. What the result means (55 seconds).** Horvath2013 is an age estimate in years; DunedinPACE is a pace measure. Running a model successfully does not establish that its estimate is meaningful for every tissue or cohort. The bundled two-sample GEO-derived demo is a software smoke test, not a validation set. The skill does not diagnose health, establish causality, or calculate a validated age-acceleration metric.

## Live demo

From the repository root, in an environment with PyAging installed:

```bash
python skills/methylation-clock/methylation_clock.py \
  --demo --clocks Horvath2013 --imputer-strategy mean \
  --output /tmp/methylation-clock-meeting-demo
```

Use a fresh output path because the skill protects existing results. The first run may download public clock weights; no methylation input is uploaded. Open `report.md`, `tables/input_qc.json`, and `tables/predictions.csv`. The generated `reproducibility/commands.sh` replays into a fresh `replay/` subdirectory. If network is unavailable, show the already generated local demo report and explain when it was run.

Example from the 17 October PyAging 0.5.7 run:

| Demo row | Horvath2013 output |
|---|---:|
| 0 | 26.0016 years |
| 1 | 28.0660 years |

The source fixture lists chronological ages of 72 and 67. Do not interpret this difference as biological age reversal; it illustrates why tissue applicability and independent validation matter.

In the fresh 17 October run, input QC recorded 485,512 probe columns and three missing beta values. The checksum file verified all eight output artifacts, and a replay produced byte-identical predictions. All 18 skill tests passed with PyAging installed; 13 offline checks passed without it.

## Questions to invite

- Which input representations and array platforms should the skill officially support?
- Which clock or clocks should a first-time user see by default?
- What per-clock coverage and model provenance must appear before someone can interpret a result?
- What is an acceptable CI cost for actual PyAging coverage?
- Should the demo remain a public GEO subset or become synthetic?

## Sources

- [PyAging DNA methylation tutorial](https://pyaging.readthedocs.io/en/latest/tutorials/tutorial_dnam.html): example data, preprocessing, and five-clock demonstration.
- [PyAging Horvath2013 metadata](https://pyaging.readthedocs.io/en/stable/clock_notebooks/horvath2013.html): chronological-age target and years unit.
- [PyAging DunedinPACE metadata](https://pyaging.readthedocs.io/en/stable/clock_notebooks/dunedinpace.html): pace-of-aging target.
- [Skill specification](../skills/methylation-clock/SKILL.md) and [fixture provenance](../skills/methylation-clock/data/PROVENANCE.md): ClawBio's input, output, and data-source contract.
