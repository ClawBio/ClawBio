# PharmGx reporting and interpretation resolution

- **Date:** 2026-10-14
- **Skill version:** 0.3.0
- **Scope:** Consistency and conservative interpretation of the existing DTC marker panel.

## Problem

The terminal preview and generated reports could disagree about phenotypes and
coverage. Some outputs implied normal function or genotype-guided prescribing
where the available data or reviewed evidence did not support that conclusion.

This resolution does not expand the panel to whole-genome analysis or validate
the report as a clinical prescribing tool.

## Resolved defects

| Area | Previous behavior | Resolution |
|---|---|---|
| Terminal preview | Rapid, indeterminate and other unrecognized phenotypes could appear as Normal; arbitrary inputs received a public-demo attribution. | Preserve phenotype text and actual input identity; remove sample-specific attribution and hard-coded warfarin genotypes, percentages and treatment suggestions. |
| Coverage | A detected variant could hide missing panel markers. | Withhold definitive calls on incomplete coverage and expose exact present, missing, supported and unsupported markers. |
| Genotype handling | Unsupported encodings or opposite-strand calls could silently fall through as reference. | Validate whitelisted nucleotide alphabets; normalize only disjoint, unambiguous complements; never interpret unsupported encodings as reference. |
| Structural/repeat markers | Array encodings could appear to establish CYP2D6*6, CYP3A5*7 or UGT1A1*28. | Explicitly withhold unsupported structural/repeat interpretations. Presence in the input is not equivalent to clinical callability. |
| CLI/API agreement | UGT1A1 withholding was implemented only in the command-line path. | Share gene profiling across CLI and API; publish coverage and limitations in structured results. |
| NUDT15 nomenclature | rs147390019 was labeled *2. | Correct its name to *4. Withhold phenotype translation pending verification of current clinical function rather than inheriting the old mapping. |
| Thiopurines | Lookup used TPMT alone despite unresolved NUDT15. | Review both genes. Uncertainty prevents standard reassurance; an established risk in one gene is retained when the other is unresolved. Do not calculate indication-specific doses. |
| Clozapine | CYP1A2 markers could produce an ultrarapid clinical phenotype and dose alert. | Withhold the clinical phenotype and explain that CYP1A2 genotype alone does not establish a clozapine dose adjustment. |
| Methotrexate | MTHFR annotations implied a validated dose adjustment. | Distinguish descriptive variant annotations from prescribing evidence and remove the genotype-only dose alert. |
| Methadone | CYP2B6 intermediate/poor phenotypes triggered automatic caution. | Align supported phenotypes with CPIC standard prescribing, titration and monitoring; include the rationale and source. |
| Warfarin | High sensitivity or poor clearance was treated as a blanket avoidance rule. | Flag clinical dosing review using validated algorithms and INR monitoring, not a universal contraindication or calculated dose. |
| Report disclosures | CYP2B6 appeared as both profiled and unassessed; method counts were stale; HTML omitted important limitations. | Correct the unassessed-gene list, derive counts from definitions and share coverage/limitation disclosures between report formats. |
| Evidence display | Catalog associations could be presented alongside unresolved patient interpretations as if actionable. | Suppress phenotype-specific recommendation extraction for unresolved results and remove evidence badges from insufficient-data rows. |
| Provenance | JSON hashed the input path rather than file contents; reproduction commands omitted usable paths and options. | Hash file contents; write quoted paths and enrichment options through shared helpers, refuse reproduction-output overwrite, and provide environment and output-checksum files. |

## Output contract

The reporter continues to produce `report.md`, `report.html` and `result.json`,
plus `reproducibility/commands.sh`, `environment.yml` and `checksums.sha256`.

- JSON retains raw panel calls, normalized supported calls, coverage, gene profiles,
  recommendations, reasons and cited sources where reviewed.
- The four existing classification buckets remain compatible.
- STANDARD means no adjustment flagged by this limited lookup, not proven safety
  or efficacy.
- INSUFFICIENT includes unavailable or unsupported patient data and insufficient
  prescribing evidence; it is not a contraindication.
- `--no-enrich` is available through the ClawBio runner for fully offline reporting.
- A local corrected report was regenerated separately from the original. Input
  integrity, cross-format phenotype agreement and output checksums were verified.
  No patient artifacts are included in this document or commit.

## Validation

Validation was performed locally with CPython 3.13.15 on macOS.

- Existing baseline: **71 PharmGx tests passed** before changes.
- Red/green regression cycle: new synthetic cases demonstrated the defects before
  implementation; additional provenance and warfarin regressions also failed before
  their fixes.
- Final focused suite: **116 tests passed**, including 45 added safety test cases
  and a Hypothesis strand-normalization property.
- Broader affected-scope run: **264 passed, 1 failed**. The failure was reproduced
  independently in unchanged shared-report files.
- Final broader run excluding that known baseline failure: **264 passed,
  1 deselected**. These counts include the PharmGx suite; they are not additive.
- Ruff static correctness checks (`E9,F63,F7,F82`), runtime-module compilation and
  `git diff --check` passed.
- Catalog regeneration changed only the PharmGx description and version.
- Pytest emitted temporary-directory cleanup warnings; the full repository suite
  and other Python versions were not run.

### Known unrelated baseline failure

`clawbio/common/tests/test_report.py::test_generate_report_header_requires_skill_version`
expects a `TypeError` when the version is omitted, while the unchanged
`generate_report_header` function supplies a default version argument.

Neither that function nor its test was modified. Resolving their contract mismatch
is a separate task.

### Re-running validation

Test dependencies were supplied through isolated uv invocations; project dependency
and lock files were not changed merely to install test tooling.

```bash
uv run --locked --with pytest --with hypothesis python -m pytest \
  skills/pharmgx-reporter/tests/ -q

uv run --locked --with pytest --with hypothesis python -m pytest \
  skills/pharmgx-reporter/tests/ \
  clawbio/tests/test_package_api.py clawbio/tests/test_parsers.py \
  clawbio/common/tests/test_report.py clawbio/common/tests/test_html_report.py \
  clawbio/common/tests/test_rec_shortener.py \
  tests/test_cli_smoke.py tests/test_core_skill_intents.py \
  tests/test_skill_intents_contract_alerts.py tests/test_generate_catalog.py \
  -q -k 'not test_generate_report_header_requires_skill_version'

uv run --locked --with ruff ruff check \
  clawbio/cli.py skills/pharmgx-reporter/pharmgx_reporter.py \
  skills/pharmgx-reporter/api.py \
  skills/pharmgx-reporter/tests/test_pharmgx.py \
  skills/pharmgx-reporter/tests/test_report_safety.py \
  --select E9,F63,F7,F82
```

## Remaining limitations

- Existing mappings for all drugs have **not** been comprehensively revalidated
  against current guidelines.
- The panel does not establish CYP2D6 copy number or hybrid alleles, UGT1A1 repeat
  length, HLA types, G6PD status or mitochondrial drug-risk variants.
- Complete coverage of the small panel does not exclude rare or untested variants.
- Clinical conditions, medication interactions and SMS deletion coordinates are not
  incorporated. This update does not infer copy-number loss from consumer-array calls.
- The report does not recommend individual doses or replace clinical-grade testing.
- Private genome files, generated reports and the unrelated nested checkout remain
  outside version control for this change.

## Follow-up integration and inventory review

The branch was subsequently merged with current upstream main. This preserved
upstream CPIC table corrections and added DPYD HapB3 support while retaining the
consumer-array safeguards; the skill version is now 0.3.1. The original local
patient report remains a historical 0.3.0 output and was not overwritten.

The inventory generator now uses SKILL.md domains by default and JSON membership
as an optional override. Unclassified skills render under Uncategorised with a
warning; retired overrides warn and are dropped; duplicate assignments still fail.
The `soul2dna` override was removed and `claw-amplicon-qc` assigned to sequence QC.
Catalog counts no longer affect the output. Tests use a synthetic tree and a golden
fixture; a live-repository freshness check is opt-in through
`CLAWBIO_CHECK_SKILLS_SUMMARY=1`.

Follow-up validation on the merged branch:

- Inventory plus PharmGx: 177 passed, one opt-in check skipped.
- Broader core/shared/root/PharmGx scope with the declared `otlp` and `fine-mapping`
  extras: 1,273 passed, 2 skipped, 1 deselected and 2 expected failures.
- `scripts/generate_skills_summary.py --check` passed on the regenerated inventory.
- The earlier shared-report version-argument failure was fixed upstream and is
  no longer excluded from this broader validation.
- These are local results, not a claim that GitHub CI or a maintainer approved the PR.

Use `uv run --locked --extra otlp --extra fine-mapping --with pytest --with hypothesis`
for that broader validation. Missing those extras caused six dependency-related
failures on the first attempt; no dependency or lockfile edits were needed.

## Reviewed evidence

- [CPIC CYP2B6 and methadone guideline, 2024](https://files.cpicpgx.org/data/guideline/publication/methadone/2024/38863207.pdf)
- [CPIC thiopurine guideline, 2025 update published in 2026](https://files.cpicpgx.org/data/guideline/publication/thiopurines/2026/41618934.pdf)
- [Thiopurine supplement, including NUDT15 marker nomenclature](https://files.cpicpgx.org/data/guideline/publication/thiopurines/2026/41618934-supplement.pdf)
- [NCBI: Clozapine therapy and CYP genotype](https://www.ncbi.nlm.nih.gov/books/NBK367795/)
- [CPIC guideline prioritization, including MTHFR/methotrexate](https://www.clinpgx.org/cpic/prioritization/guidelines)
- [CPIC pharmacogenetics-guided warfarin dosing](https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/)

## Disclaimer

ClawBio is a research and educational tool. It is not a medical device and does
not provide clinical diagnoses. Consult a healthcare professional before making
any medical decisions.
