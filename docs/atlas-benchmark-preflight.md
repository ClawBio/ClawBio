# Atlas benchmark preflight

Offline check of whether two **declared observation metadata exports** support
a donor-disjoint, cell-disjoint benchmark split. This is a benchmark utility,
not a new atlas connector, analysis skill, model benchmark or clinical tool.

## Why this check exists

The combined Tabula Sapiens 2.0 atlas incorporates earlier samples. A release
called "v2" must not automatically be treated as an unseen test set for a model
trained on v1. The [CZI benchmark data card](https://virtualcellmodels.cziscience.com/dataset/tabula-sapiens)
explicitly describes removing v1 samples previously used in training. This
utility checks the exported identities instead of trusting release names.

The [official data registry](https://registry.opendata.aws/tabula-sapiens/)
distinguishes processed data from controlled raw FASTQ access. This utility
downloads neither. Do not submit donor clinical histories or raw reads to it.

These are data-release claims, not a claim that the 2025 preprint and the final
[2026 Cell article](https://doi.org/10.1016/j.cell.2026.08.010) contain identical
analyses. Keep publication versions separate from dataset versions.

## Run

From the repository root, with Python 3.10 or later and no third-party runtime dependencies:

```bash
python tests/benchmark/atlas_preflight.py \
  --input /path/to/manifest.json --output /path/to/new-report
```

Gate a downstream command on success with `&&`. This utility is **opt-in**; it
does not silently change existing benchmark runners or the ClawBio skill CLI.
Its regression tests are discovered by the existing `tests/benchmark` pytest path.

```bash
python tests/benchmark/atlas_preflight.py --input manifest.json --output new-report && \
  python your_existing_benchmark.py
```

Exit statuses:

| Exit | Result | Meaning |
|---|---|---|
| 0 | PASS | No donor/cell overlap in the declared metadata; warnings may remain |
| 1 | FAIL | Overlapping donors/cells or duplicate cell IDs within a split |
| 2 | ERROR | Missing, malformed, inconsistent or unpinned input, or output cannot be safely written |

An existing output path is always refused, even when empty. Choose a new directory.

## Input contract

Copy the shape of `tests/benchmark/fixtures/atlas_preflight/demo_manifest.json`,
but replace all **synthetic** provenance, identifiers, paths, hashes and counts.
Never rename the toy fixture to suggest it is a real atlas benchmark.

Each of the exactly two splits, `train` and `test`, declares:

- `dataset_id`: source dataset identifier.
- `dataset_version`: explicit release identifier, not `latest`, `stable`, `current`, `main`, `master` or `HEAD`.
- `source_url`: source citation URL, recorded but never fetched.
- `identity_namespace`: the **same stable study-level identity convention** in both splits.
- `observations_path`: local JSONL path, relative to the manifest directory or absolute.
- `observations_sha256`: SHA-256 of the exact JSONL bytes.
- `expected_cells`: positive expected number of observation records.

JSONL has one object per cell, with required string `cell_id` and `donor_id`.
Optional `sex` and `ancestry` fields are donor metadata strings or null.
Other fields are ignored. Empty lines are ignored for cell counting but remain
included in the checksum. Duplicate JSON keys and non-finite JSON numbers are rejected.

Use original globally unique cell identities preserved across releases. A bare
10x barcode is **not** globally unique. An identity may need the original
donor, library/sample and barcode. Do not prepend the release version or split
name: doing so would hide overlapping cells. Donor IDs must also be harmonised
across releases before this audit. Separate identifier namespaces are rejected,
not treated as evidence of independence.

Export all observations actually used by each pipeline from its AnnData `obs`
metadata, preserving the originals. Verify the export count against the matrix
before declaring `expected_cells`. Hash the export with `sha256sum` on Linux or
`shasum -a 256` on macOS. The utility does not verify this export against H5AD;
the caller remains responsible for completeness and faithful identity mapping.

## Outputs and provenance

- `report.md`: status, split counts, finding codes and interpretation limits.
- `result.json`: structured findings, versions, source citations, observed checksums,
  overlap counts and donor-level sex/ancestry coverage.
- `reproducibility/manifest.json`: exact input snapshot, if the input was readable.
- `reproducibility/environment.json`: Python/platform, utility version and script hash.
- `reproducibility/commands.sh`: original shell-quoted invocation and working directory.
- `reproducibility/checksums.sha256`: hashes of every generated bundle file except itself.

Identity lists and expression values are not written to the report. The source
JSONL files remain local and are not copied into the report. Source paths and
demographic distributions may still be sensitive: review a report before sharing.
No shell commands are executed by the utility and no network requests are made.

Checksums bind the files to this run; they are not digital signatures proving
the source declarations true. Retain the original input exports for replay.

## What a PASS does not mean

It does not establish that a model never saw the data, that releases are complete,
that cell annotations are correct, that donors represent global populations, or
that preprocessing and tuning were conducted without test-set leakage. It also
does not perform similarity-based duplicate detection on expression matrices.

Sex and ancestry summaries count **unique donors, not cells**. Missing ancestry
remains unknown. Conflicting labels are flagged and excluded from distributions;
no label is inferred. Warnings do not invalidate a disjoint split, but may limit
the interpretation of downstream results. There is no invented "equity score".

## Verification and bounded scope

```bash
python -m pytest tests/benchmark/test_atlas_preflight.py -q
python tests/benchmark/atlas_preflight.py --demo --output fresh-demo-output
```

The demo deliberately returns **FAIL / exit 1**, with one overlapping donor and
one overlapping cell. It is a successful negative-control check on toy data,
not a failed scientific experiment. Positive-control and malformed-input cases
are included in the regression tests. No real Tabula Sapiens benchmark has been
run by this implementation.

Design decision: extend the existing benchmark utilities, preserving existing
single-cell analysis skills. A new atlas service, separate research programme,
full million-cell download and LLM-based clinical-record interface were excluded.
The immediate deliverable is reproducible evaluation protection, not a claim
of adoption, agent superiority or time saved.

ClawBio is a research and educational tool. It is not a medical device and does
not provide clinical diagnoses. Consult a healthcare professional before making
any medical decisions.
