# Public Visium validation: Human Lymph Node

Validation date: 2026-09-21. This case uses measured public SpaceRanger counts,
not the synthetic demo. No patient data or downloaded matrices are committed.

Dataset: [10x Genomics Human Lymph Node, standard 1.1.0](https://www.10xgenomics.com/datasets/human-lymph-node-1-standard-1-1-0).

| Archive | Bytes | SHA-256 |
|---|---:|---|
| `V1_Human_Lymph_Node_filtered_feature_bc_matrix.tar.gz` | 81075087 | `93f7424de945eb886db17e5184d2112c77855bd54c330b2208a846870595e4e8` |
| `V1_Human_Lymph_Node_spatial.tar.gz` | 8246238 | `812808883366ff9623dc8354847a7211b0d922b2bfc4c9359d6e12e993ea6a73` |

The helper verifies these digests before extracting newly downloaded archives.
Each analysis also records hashes of its actual matrix, barcodes, features and
positions in `reproducibility/run_manifest.json`. A cache hit by itself is not a
claim of biological correctness.

## Commands

From the repository root, provision the public input explicitly:

```bash
uv run --extra spatial python -c "import sys; sys.path.insert(0, 'skills/spatial-transcriptomics'); from spatial_transcriptomics import ensure_public_visium_outs; print(ensure_public_visium_outs())"
```

Use fresh output directories for each run:

```bash
uv run --extra spatial python skills/spatial-transcriptomics/spatial_transcriptomics.py \
  --input ~/.cache/clawbio/visium/V1_Human_Lymph_Node/outs \
  --output /tmp/visium_ln_500 --n-top-hvg 500

uv run --extra spatial python skills/spatial-transcriptomics/spatial_transcriptomics.py \
  --input ~/.cache/clawbio/visium/V1_Human_Lymph_Node/outs \
  --output /tmp/visium_ln_2000

INPUT_PATH=~/.cache/clawbio/visium/V1_Human_Lymph_Node/outs \
REPLAY_OUTPUT=/tmp/visium_ln_replay \
PYTHON="$PWD/.venv/bin/python" \
  bash /tmp/visium_ln_500/reproducibility/commands.sh
```

## Observed full-slide results

Environment: Python 3.14.4, Scanpy 1.12.4, NumPy 2.4.6, scikit-learn 1.9.0.
Embedding defaults: 8 PCs, 8 expression neighbours, Leiden resolution 0.5,
random state 7. Spatial graph: coordinate kNN with k=6. Neighbourhood null:
50 label permutations. No mitochondrial ceiling applied.

| Quantity | 500 HVGs | Default 2000 HVGs |
|---|---:|---:|
| CLI exit status | 0 | 0 |
| Spots loaded / retained | 4035 / 4035 | 4035 / 4035 |
| Genes loaded / after QC | 36601 / 25187 | 36601 / 25187 |
| Genes evaluated by Moran and Wilcoxon | 500 | 2000 |
| Leiden clusters | 8 | 11 |
| Top Moran gene | IGHG2 | IGKC |
| Top Moran statistic | 0.653139 | 0.729871 |
| Co-occurrence tensor shape | 8 × 8 × 6 | 11 × 11 × 6 |
| Exported co-occurrence rows | 384 | 726 |

All exported co-occurrence ratios were finite in both runs. CSV values matched
the full JSON tensor. The 500-HVG run retained the earlier Moran, marker and
neighbourhood values; co-occurrence now follows the corrected cumulative-radius
definition and is actually exported.

QC identified 13 `MT-` gene symbols. Across retained spots, median total counts
were 20,239, median detected genes 5,999, median mitochondrial percentage
0.989313%, and maximum mitochondrial percentage 5.410122%. These are diagnostics,
not evidence that technical spatial confounding has been eliminated.

The generated 500-HVG replay script was executed from another working directory.
Its parameters, effective embedding dimensions, tested-gene list, QC summary,
spot/gene/cluster counts, top Moran values, marker rows, neighbourhood matrix,
co-occurrence tensor and input provenance matched the original run. All 15 files
listed in each original run's checksum manifest verified, including the primary
report, JSON and run manifest. Timestamps/output locations are not numerical
reproducibility targets.

The different 8/11 cluster counts show why the HVG parameter must be preserved in
the replay bundle. Two analyses with different gene selection are not repeated
runs of the same configuration.

## Limits and independent reference checks

This demonstrates full-slide input handling, output delivery and same-environment
replay. It does not validate cell-type annotations, a clinical atlas, a second
tissue or cross-platform bitwise equivalence. Moran values are descriptive and
have no permutation p-values. Marker tests after clustering are exploratory.

The offline fixture in `../fixtures/squidpy_v1_6_co_occurrence.json` was generated
by actual Squidpy 1.6.0 in a separate compatible environment, with explicit
coordinates, labels and radii. It also records Moran's I (-0.13855421686746988)
and observed neighbourhood counts on a four-spot directed graph. Its regeneration
script is provided alongside it. These are small, defined estimator checks;
they do not claim that Squidpy's automatic graph/radius choices or random streams
reproduce the full-slide pipeline. Permutation z-scores are not a bitwise parity
target because the two implementations use different random streams.

The default test suite remains offline. To run the separate 400-spot measured
input integration test, opt in explicitly; download/processing failures then fail:

```bash
CLAWBIO_RUN_PUBLIC_VISIUM=1 uv run --extra spatial --with pytest \
  pytest skills/spatial-transcriptomics/tests/ -m network -q
```
