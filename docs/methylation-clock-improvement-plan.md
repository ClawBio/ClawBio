# Methylation clock improvement plan

Prepared for the contributor meeting on **Monday, 19 October 2026**.

## Goal

Make the skill dependable for research use: inputs are checked before inference, reports show what PyAging actually did, and routine CI catches regressions without downloading clock weights on every run. Keep interpretation tied to each clock's documented target, tissue, and unit.

## Status as of 17 October

| Priority | Work | Acceptance evidence | Status |
|---|---|---|---|
| P0 | Keep metadata out of EPICv2 probe aggregation | The unfixed demo fails on `tissue`; the patched Horvath2013 demo completes with PyAging 0.5.7. | In draft fork PR #1; live run verified locally. |
| P0 | Check input beta values, sample IDs, and categorical metadata | Invalid probes and duplicate IDs fail early; a CSV `sample_id` remains in predictions. | Implemented on `codex/methylation-clock-improvements`; offline tests. |
| P0 | Give CI an offline regression check | A mocked aggregator sees only probe columns, including on a machine without PyAging. | Implemented; real PyAging preprocessing checks still run when installed. |
| P1 | Make outputs auditable | Report states clock units and input counts; `result.json`, QC table, and checksums include the report and result. | Implemented; live demo and checksum verification passed. |
| P1 | Keep GEO downloads inside the chosen output directory | A mock download writes under `<output>/input/`; no data file is created in the checkout. | Implemented; live GEO download is not yet verified. |
| P1 | Improve user instructions and meeting demo | `SKILL.md` documents beta input, sample IDs, trust boundary for pickle, and first-run weight download. | Implemented. |

## Decisions for the contributors

1. **CI dependency strategy.** `uv sync --all-extras` runs for three Python versions. Adding PyAging there also pulls Torch and Linux CUDA dependencies into every matrix job. Decide whether to provision a dedicated CPU-only PyAging job, pin a PyAging version for that job, and run one live clock test on a schedule. Keep the no-download regression check in ordinary CI.
2. **Input contract.** Confirm beta values in `[0, 1]` as the supported input. Decide how to handle normalized values, M-values, and named sample columns beyond `sample_id`. Reject unsupported representations with a clear error instead of guessing a conversion.
3. **Metadata and covariates.** Agree which clocks require age or sex and how missing or ambiguous values should be handled. The current change leaves unknown sex as missing rather than treating it as male; this needs clock-specific review.
4. **Default clocks.** Five defaults require five model downloads and mix age estimates with pace measures. Decide whether the default should remain five clocks or become a single quick-start clock, with an explicit comparison option.
5. **Demo data.** The current two-sample fixture is a subset of public GEO accession GSE139307, not synthetic patient data. Confirm that a public fixture is acceptable under the repository's synthetic-demo rule or replace it with clearly labelled synthetic data.

## Next development stages

### Completed before the meeting

- A fresh Horvath2013 demo completed with PyAging 0.5.7. Input QC counted 485,512 probe columns and three missing beta values across two samples. All eight recorded output checksums verified.
- The generated `commands.sh` reproduced the demo in a fresh `replay/` directory; prediction CSV files matched byte for byte.
- The core-only environment ran 13 tests with two PyAging-dependent checks skipped. With PyAging and live clock weights available, all 18 skill tests passed.

### Still to settle at the meeting

- Review the code diff with contributors; keep this branch separate from the narrow metadata fix in draft PR #1 until they agree on how to split or stack the changes.
- Run broader repository contract tests before opening a PR to upstream; this pass verified the touched skill and replay only.

### Following the meeting

- Add a dedicated PyAging CI job using CPU wheels, with documented runtime and cost. Test actual preprocessing on supported Python versions without clock downloads on each PR.
- Add an optional live integration gate or scheduled run that downloads one clock and checks output structure and reproducibility.
- Add per-clock input coverage diagnostics and a documented interpretation table sourced from PyAging metadata and the original clock publications. Do not introduce arbitrary quality cutoffs without validation.
- Review whether model versions and weights should be pinned or recorded more explicitly than the current package version and PyAging metadata.
- Decide whether to offer a separate, properly defined age-acceleration analysis. It would need chronological ages and a chosen reference method; a simple difference is not a substitute.

## Evidence and limits

- The unfixed two-sample demo failed with `Feature columns must be numeric ... ['tissue']` on PyAging 0.5.7. Applying the draft PR patch produced Horvath2013 values of 26.0016 and 28.0660, with report and tables. This verifies execution, not biological validity.
- The fixture contains two sperm samples with `age` values 72 and 67. The gap between those ages and the Horvath estimates must not be presented as rejuvenation or a clinical finding.
- The offline tests use small frames and a mocked predictor where indicated. They catch input and reporting regressions, not changes in downloaded model weights.
- PyAging's [DNA methylation tutorial](https://pyaging.readthedocs.io/en/latest/tutorials/tutorial_dnam.html) provides the GSE139307 workflow and EPICv2 aggregation example. Its [Horvath2013 page](https://pyaging.readthedocs.io/en/stable/clock_notebooks/horvath2013.html) describes the chronological-age target. [DunedinPACE](https://pyaging.readthedocs.io/en/stable/clock_notebooks/dunedinpace.html) has a different target and tissue context.
