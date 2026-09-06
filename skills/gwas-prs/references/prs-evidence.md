# PRS evidence policy 1.0.0

This gate assesses **research percentile eligibility**, not clinical utility.
It implements deterministic checks against a caller-curated evidence manifest.
It does not verify publications, infer ancestry, measure equity from population
labels, or convert heritability/Mendelian randomisation into predictive accuracy.

## Invocation and migration

Real-input runs now withhold percentile, z-score and risk category unless a
complete manifest is supplied with `--evidence-json manifest.json`.
Raw sums remain available for audit, including partial scores previously skipped
by `--min-overlap`. They must not be interpreted when the assessment is withheld.
The built-in `--demo` remains explicitly `synthetic_demo_only`; providing
`--evidence-json` also gates that demo, and illustrative panels cannot qualify.

The public CLI accepts the flag:
`python clawbio.py run prs --input genotypes.txt --pgs-id PGS... --evidence-json manifest.json --output new-directory`.
Use the registered alias printed by `clawbio.py list` if it differs in your release.
No new dependency, remote API, patient upload or separate skill is introduced.

## Manifest schema

All hashes are lowercase SHA-256 of exact file bytes, including compression.
Schema version is the integer 1. The root has `input` and `scores` objects.
Unknown fields are ignored; missing evidence withholds interpretation.

| Object | Required fields |
|---|---|
| input | sha256, build (GRCh37/GRCh38), strand (forward), population, context, source |
| scores[exact score_id] | sha256, version, source, build, variant_count, validation, reference |
| validation | population, context, independent (boolean true), n (positive integer), source, metric, estimate, ci_lower, ci_upper |
| reference | population, context, source, n (positive integer), mean, sd, distribution (normal), missingness (complete_only) |

Metric is `r2`, `auc` or `correlation`. Confidence limits must contain the
estimate and fall within the metric's domain. All numbers must be finite;
booleans are not numbers. AUC/R2 use [0,1], correlation [-1,1].
No minimum performance threshold is invented: a percentile is a position in a
reference distribution, not a claim of useful disease prediction.

The validation and reference apply to the enclosing exact score version/hash.
Both must match the declared target population **and context**. These are
explicit cohort/context identifiers chosen by the evidence curator, not inferred
identity labels. Record age range, ascertainment and relevant assay/phenotype
conditions in the cited context source. String equality checks the declaration;
it does not scientifically establish applicability.

Input build and forward strand are external declarations bound to the input
hash. The gate checks their consistency with the requested/declared score build;
it does not independently infer build or resolve strand. Only complete,
diploid, biallelic SNV scores with explicit effect and other alleles qualify in
v1. Missing calls, incompatible alleles, unsupported variant types, duplicate
rsIDs, mismatched declared counts, or non-finite weights withhold interpretation.
No strand guessing, lift-over, imputation or missingness extrapolation is done.
A variant count from a trusted score manifest guards against silently dropped
parser rows. The curator must verify that count against the original score.

Optional `heritability` and `causal_evidence` objects are retained separately.
They never satisfy the validation requirement. `causal_claim_supported` is always
false: causal assessment is outside this gate's scope.

## Output

Every real-input result contains `evidence_assessment`, with policy version,
input/score hashes, evidence source/version, performance, separate heritability
and causal evidence, status and machine-readable reasons with remediation.
Statuses are `supported`, `not_established`, and `incompatible`.
Supported means **eligible under supplied evidence**, not independently verified.
A supported normal-reference percentile is 100 * Phi((score - mean)/sd).
Risk category remains null even when the percentile is supported.

Decisions propagate into Markdown, compact JSON and the standard result envelope.
The evidence file hash is in replay provenance; the file itself is not copied.
Replay requires `PRS_EVIDENCE_FILE` pointing to the original evidence JSON.
Evidence and reports can contain sensitive context and should stay local.
The CLI refuses to overwrite existing report artifacts; choose a new output path.

## Scientific motivation and fixtures

Schwaba et al., Nature (2026), DOI:
https://doi.org/10.1038/s41586-026-10992-9
motivates separating population association/heritability, independent prediction,
familial confounding, and population portability. It does not supply this policy's
software thresholds or validate an individual personality report.

Tests include personality and metabolic contexts using **synthetic evidence,
weights and reference distributions**. No published score or validation metric
is fabricated or attributed to that paper. These fixtures establish software
behaviour only; they are not a reproduced experiment, external expert review,
or an independently adjudicated scientific benchmark.

Run `python -m pytest skills/gwas-prs/tests/test_prs_evidence.py -q`.
Before deployment, curate evidence for a real score, independently review it,
and validate performance and refusal rates on relevant external cohorts.

## Runnable offline demonstration

Run from the repository root:

```bash
python skills/gwas-prs/benchmark_evidence.py --demo --output /tmp/prs-evidence-demo
```

This produces result.json and report.md for 14 fixed synthetic cases: two
supported and twelve withheld. Expected decisions are manually specified
software contracts; no external adjudication is claimed. The summary reports
false-supported and false-withheld counts. Use a fresh output directory.

Replay requires `PRS_REPLAY_OUTPUT` naming a fresh output directory, preserving
the original run for comparison. When evidence was supplied, also set
`PRS_EVIDENCE_FILE` to the original manifest.
