---
name: spatial-transcriptomics
description: >-
  Analyse 10x Visium spatial transcriptomics: SpaceRanger outs or spatial h5ad
  in, then QC, Leiden clustering, Wilcoxon markers, Moran's I, neighbourhood
  enrichment and co-occurrence in one local report.
license: MIT
metadata:
  version: "0.1.0"
  author: Zhihao Wan
  domain: spatial-transcriptomics
  tags:
    - spatial-transcriptomics
    - visium
    - scanpy
    - moran
    - clustering
  inputs:
    - name: visium_input
      type: file
      format:
        - dir
        - h5ad
      description: SpaceRanger outs/ or h5ad with raw counts and finite two-dimensional obsm['spatial'] coordinates
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Analysis report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable summary
  dependencies:
    python: ">=3.11"
    packages:
      - scanpy>=1.10
      - leidenalg>=0.10
      - numpy>=1.24
      - pandas>=2.0
      - matplotlib>=3.7
      - scikit-learn>=1.3
      - scipy>=1.10
  demo_data:
    - path: examples/demo_spec.json
      description: Recipe for the offline 8x8 two-domain synthetic grid
  endpoints:
    cli: python skills/spatial-transcriptomics/spatial_transcriptomics.py --input {visium_input} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: scanpy
      - kind: pip
        package: leidenalg
    trigger_keywords:
      - visium
      - spatial transcriptomics
      - SpaceRanger
      - Moran's I
      - neighbourhood enrichment
      - spatially variable genes
---

# 🧬 Spatial Transcriptomics (Visium)

You are **spatial-transcriptomics**, a ClawBio agent that analyses measured 10x Visium data. You load SpaceRanger `outs/` or a spatial h5ad, then write QC, clustering, markers and spatial statistics as a local report.

## Trigger

**Fire this skill when the user says any of:**
- "analyse my Visium data"
- "spatial transcriptomics QC and clustering"
- "SpaceRanger outs"
- "spatially variable genes"
- "Moran's I on visium"
- "neighbourhood enrichment"
- "spot co-occurrence"

**Do NOT fire when:**
- The user has an H&E tile and wants predicted expression. That is `deepspot-m`.
- The user has dissociated scRNA-seq (h5ad/mtx with no `obsm['spatial']`). That is `scrna-orchestrator`.
- The user has a marker-by-spot table and wants region labels. That is `marker-dominance-mapper`.

## Why This Exists

- **Without it**: Visium analysis is a Scanpy plus Squidpy notebook stitched by hand, with no `--demo` and no reproducibility bundle.
- **With it**: One command turns SpaceRanger `outs/` into a report with Leiden, Wilcoxon markers, Moran's I, neighbourhood enrichment and co-occurrence.
- **Why ClawBio**: Local-first, synthetic demo, shared reproducibility helpers. Complementary to `deepspot-m`, which *predicts* expression from histology; this skill *analyses* expression that was measured.

## Core Capabilities

1. **Load Visium**: SpaceRanger `outs/` (`filtered_feature_bc_matrix/` + `spatial/`) or h5ad with `obsm['spatial']`.
2. **QC and clustering**: Scanpy filter, normalise, HVG, PCA, UMAP, Leiden.
3. **Markers**: Wilcoxon cluster-vs-rest on log-normalised expression.
4. **Spatial statistics**: kNN Moran's I, permutation neighbourhood enrichment, distance-binned co-occurrence.
5. **Report**: Markdown, JSON, figures, tables, reproducibility bundle.

## Scope

**One skill, one task.** Measured Visium-like spot data in, one analysis report out. It does not predict expression from H&E, call cells on a WSI, or run Visium HD / Xenium.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| SpaceRanger outs | directory | `filtered_feature_bc_matrix/` mtx + `spatial/tissue_positions.csv` | `sample/outs` |
| Spatial AnnData | `.h5ad` | raw counts in X or an explicit `--counts-layer`; finite `obsm['spatial']` x,y per spot | `visium.h5ad` |
| Demo | n/a | none | `--demo` |

HDF5 `filtered_feature_bc_matrix.h5` is not read in v0.1; pass the mtx folder. Tissue images are not required.

Counts must be finite, nonnegative integers. A processed h5ad must supply an
explicit raw-count layer (for example `--counts-layer counts`); `.raw` is not
assumed to contain counts. Normalised/log-transformed X is rejected. Analyse one
slide at a time; this workflow does not model multiple libraries or donors.

## Workflow

1. **Validate (prescriptive)**: Accept `outs/`, raw-count spatial h5ad, or `--demo`. Validate counts and finite two-dimensional coordinates. Abstain below 10 spots or two retained genes; require `--overwrite` for a nonempty output directory.
2. **Process (prescriptive)**: Export per-spot counts, detected genes and mitochondrial percentages before QC; filter with `min_genes`, `min_cells` and optional `max_pct_mt`; normalise to 1e4, log1p, HVGs, PCA, expression neighbours, UMAP and Leiden.
3. **Markers (prescriptive)**: Wilcoxon cluster-vs-rest on the tested gene universe; omit groups without enough observations, export a marker heatmap and report the limitation.
4. **Spatial graph**: k=6 nearest spots on `obsm['spatial']` (not the PCA graph).
5. **Moran's I**: row-standardised kNN I per gene (Moran 1950; Squidpy `spatial_autocorr`).
6. **Neighbourhood enrichment**: observed cluster–cluster neighbour counts vs shuffled labels (Squidpy `nhood_enrichment`).
7. **Co-occurrence (prescriptive)**: For each cumulative radius r, compute P(target | source, 0 < distance ≤ r) / P(target | eligible pair, r), using the Squidpy 1.6.0 estimator and six positive-distance quantile radii. Export all scores, cluster axes and radii.
8. **Generate (prescriptive outputs, flexible narrative)**: Write `report.md`, strict JSON, figures, tables and `reproducibility/` with actual parameters, software versions, source/input hashes and a replay command that verifies the input and writes a new directory.

Steps 1–7 are prescriptive. Report narrative is flexible.

## CLI Reference

```bash
python skills/spatial-transcriptomics/spatial_transcriptomics.py \
  --input sample/outs --output /tmp/visium_out

python skills/spatial-transcriptomics/spatial_transcriptomics.py \
  --input visium.h5ad --output /tmp/visium_out

python skills/spatial-transcriptomics/spatial_transcriptomics.py \
  --demo --output /tmp/spatial_demo

python clawbio.py run spatial --input sample/outs --output /tmp/visium_out
python clawbio.py run spatial --demo

# A processed h5ad with a preserved raw-count layer:
python clawbio.py run spatial --input processed.h5ad --counts-layer counts \
  --n-pcs 30 --n-neighbors 15 --nhood-perms 1000 --output /tmp/visium_review
```

| Flag | Default | Purpose |
|------|---------|---------|
| `--min-genes` | 5 | Drop spots with fewer genes |
| `--min-cells` | 1 | Drop genes in fewer spots |
| `--leiden-resolution` | 0.5 | Leiden resolution |
| `--n-top-hvg` | 2000 | Highly variable genes (capped at the gene count) |
| `--random-state` | 7 | PCA / neighbours / Leiden / permutations |
| `--n-pcs` | 8 | PCA components, capped by spots and selected genes |
| `--n-neighbors` | 8 | Expression graph neighbours; spatial k stays 6 |
| `--nhood-perms` | 50 | Label permutations for exploratory neighbourhood z-scores |
| `--top-markers` | 5 | Reported markers per supported cluster |
| `--max-pct-mt` | none | Optional mitochondrial percentage ceiling (0–100) |
| `--counts-layer` | none | Explicit raw-count layer for h5ad |
| `--overwrite` | false | Explicitly replace report files in a nonempty directory |
| `--expected-input-sha256` | none | Replay integrity check before analysis |

## Demo

```bash
python clawbio.py run spatial --demo
```

Expected output: 64-spot synthetic grid, two spatial domains, Leiden ≥ 2, EPCAM/COL1A1 among high Moran's I genes, figures, tables, reproducibility bundle. No download.

## Algorithm / Methodology

1. **Load**: `scanpy.read_10x_mtx` plus `tissue_positions.csv` (or `tissue_positions_list.csv`); keep `in_tissue==1`.
2. **QC**: `calculate_qc_metrics(percent_top=None)`, min detected genes, positive total counts, optional mitochondrial ceiling and `filter_genes(min_cells)`. Mitochondrial genes match case-insensitive `MT-`; absence of such symbols is reported as unavailable mitochondrial QC.
3. **Normalise**: `normalize_total(1e4)`, `log1p`. Raw counts kept in `layers["counts"]`.
4. **Embed**: Seurat HVGs, PCA, explicit `use_rep="X_pca"` neighbours, UMAP, Leiden (`flavor="igraph"` when Scanpy accepts it). Requested and effective embedding dimensions are recorded.
5. **Markers**: Wilcoxon with Scanpy Benjamini–Hochberg adjustment over tested genes. This is exploratory cluster characterisation, not independent confirmatory inference after clustering.
6. **Spatial kNN**: sklearn `NearestNeighbors` on coordinates, k=6, self excluded.
7. **Moran's I**: I = (zᵀWz)/(zᵀz) with row-standardised W.
8. **Enrichment**: Configurable label permutations on the frozen directed spatial graph; z = (obs − mean_null) / sd_null. Zero null SD is undefined (`null` in JSON, blank in CSV, `NA` in the report).
9. **Co-occurrence**: Cumulative Euclidean radii, excluding zero-distance pairs; conditional target frequency divided by the target marginal over eligible pairs. The six quantile radii differ from Squidpy's automatic radius selection. Tensor axis order is source cluster, target cluster, radius; radius intervals are `(0, r]`. Unsupported ratios are undefined.

Above 80 post-QC genes, **both Moran and Wilcoxon evaluate HVGs only**;
otherwise every retained gene is evaluated. `result.json.analysis_scope` records
the exact tested gene names and count. A gene missing from the Moran table has
not been evaluated, and cannot be called spatially neutral. Moran is a descriptive
statistic with no permutation p-value or multiple-testing correction. Constant
genes have undefined Moran's I, represented as `null`/blank/`NA`.

**Key thresholds**:
- Minimum spots: 10 (below this the kNN graph is not meaningful)
- Spatial k: 6 (hex-like Visium neighbourhood)
- Leiden resolution: 0.5 (demo default; user-overridable)
- Permutations: 50 (exploratory default; configurable with `--nhood-perms`)
- Mitochondrial ceiling: none by default; inspect QC and choose a tissue-appropriate threshold rather than applying a universal cutoff

## Example Queries

- "Run QC and clustering on this Visium outs folder"
- "Which genes are spatially variable in my Visium sample?"
- "Neighbourhood enrichment on my visium h5ad"

## Example Output

```markdown
# Spatial Transcriptomics Report (demo)

**Spots**: 64
**Leiden clusters**: 2

## Spatially variable genes (Moran's I)
| Gene | Moran's I |
|------|-----------|
| DCN | 0.745 |
| VIM | 0.639 |
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── figures/
│   ├── umap_leiden.png
│   ├── spatial_leiden.png
│   ├── marker_heatmap.png  # optional when no valid cluster-vs-rest markers exist
│   └── qc_spot_metrics.png
├── tables/
│   ├── markers_top.csv
│   ├── moran_i.csv
│   ├── nhood_enrichment.csv
│   ├── co_occurrence.csv
│   ├── qc_spot_metrics.csv
│   └── qc_summary.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    ├── checksums.sha256
    └── run_manifest.json
```

## Dependencies

**Required**:
- `scanpy` >= 1.10; QC, HVG, PCA, UMAP, Leiden, Wilcoxon
- `leidenalg` >= 0.10; Leiden
- `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `scipy`; spatial graph, stats, figures

**Not required**:
- `squidpy`. Its SpatialData/Dask/OME-Zarr dependency chain can pull S3 dependencies into every `uv sync --all-extras` job. Runtime estimators are implemented locally. Co-occurrence follows the explicitly cited 1.6.0 cumulative-radius definition, not the annular 1.4 definition; automatic graph/radius construction and random streams are not claimed to be interchangeable.

## Gotchas

- **You will want to route an H&E tile here. Do not.** This skill needs measured spot counts and coordinates. Predicted expression from histology is `deepspot-m`.
- **You will want to cluster on the spatial kNN graph. Do not, unless you mean it.** Leiden uses the PCA neighbour graph. The spatial kNN graph is only for Moran's I and enrichment. Mixing them silently changes what a cluster is.
- **You will want to treat Moran's I as a p-value. Do not.** v0.1 reports the statistic, not a permutation p-value per gene.
- **You will want to quote demo numbers as a tissue result. Do not.** `--demo` is an 8×8 synthetic grid.
- **You will want to pass a Visium HDF5 matrix. Do not in v0.1.** Supply the mtx `filtered_feature_bc_matrix/` directory.
- **You will want to threshold enrichment at |z|>1.96 as a discovery claim. Do not.** 50 permutations make the tails coarse; the diagonal sign is the supported reading.
- **You will want to load already normalised X as counts. Do not.** Select a verified raw-count layer explicitly; a `.raw` attribute is not proof of raw counts.
- **You will want to replay into the original output. Do not by default.** `commands.sh` writes to `replay/` or `$REPLAY_OUTPUT`, verifies `$INPUT_PATH` against the recorded digest, and preserves all analysis parameters.
- **You will want to label a missing Moran gene as non-spatial. Do not.** Check `analysis_scope`; HVG screening omits untested genes.

## Safety

- **Local-first**: Spots are read from disk. No upload. `--demo` does not download a Visium dataset.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **No hallucinated science**: Cluster labels, Moran's I and enrichment come from the matrices above.
- **Audit trail**: `reproducibility/commands.sh`, `environment.yml`, `checksums.sha256` via `clawbio.common.reproducibility`.
- **Archive safety**: The public-data helper verifies pinned SHA-256 digests, stages downloads, rejects traversal/links/special files, and copies only known matrix/position/scalefactor files. Images are not extracted.

## Validation and Replay

```bash
# Default tests and demo are offline. Live-test failures are failures when opted in.
uv run --extra spatial --with pytest pytest skills/spatial-transcriptomics/tests/ -m 'not network'
CLAWBIO_RUN_PUBLIC_VISIUM=1 uv run --extra spatial --with pytest \
  pytest skills/spatial-transcriptomics/tests/ -m network

# On a compatible machine, install the recorded environment and use the same source.
cd /path/to/output
sha256sum -c reproducibility/checksums.sha256
INPUT_PATH=/path/to/original/outs REPLAY_OUTPUT=/tmp/visium_replay \
  bash reproducibility/commands.sh
```

The environment file pins the installed analysis dependency closure and Python
version; the manifest records the platform and source hashes. Cross-platform
bitwise numerical identity is not promised. The public integration test uses 400
measured spots for runtime; full-slide validation is a separate explicit run.

See [the full-slide validation record](examples/public_visium_validation.md).
The offline [Squidpy 1.6 fixture](fixtures/squidpy_v1_6_co_occurrence.json) was
generated in a separate pinned environment using the adjacent regeneration
script. It checks co-occurrence at explicit radii, Moran on an explicit graph,
and observed neighbourhood counts. It does not assert identical permutation
z-scores across different random streams. Undefined ratios use `null` here
instead of Squidpy's zero convention.

## Agent Boundary

The agent dispatches and explains. The Python skill loads data, runs Scanpy and the spatial estimators, and writes files. The agent must not invent Moran's I, relabel clusters, or present demo values as a patient sample.

## Integration with Bio Orchestrator

**Trigger conditions**: Visium, SpaceRanger `outs/`, spatially variable genes, Moran's I, neighbourhood enrichment, spot co-occurrence.

**Chaining partners**:
- `deepspot-m`: complementary. Predicted per-tile expression is not a Visium `outs/` tree; do not pipe it here without building a spatial AnnData first.
- `scrna-orchestrator`: dissociated scRNA-seq without coordinates.
- `marker-dominance-mapper`: downstream if you export a marker-by-spot table.

## Maintenance

- **Review cadence**: Recheck Scanpy Leiden (`flavor="igraph"`) and 10x position CSV headers each quarter.
- **Staleness signals**: SpaceRanger position file rename, Scanpy dropping `rank_genes_groups` Wilcoxon, a request for Visium HD / Xenium.
- **Deprecation**: Archive if a maintained Visium wrapper in this repo supersedes the report contract.

## Citations

- Wolf, Angerer and Theis (2018) Genome Biol 19:15. PMID 29409532. Scanpy.
- Moran (1950) Biometrika 37:17–23. Global Moran's I.
- Palla et al. (2022) Nat Methods 19:171–178. [PMID 35102346](https://pubmed.ncbi.nlm.nih.gov/35102346/). [DOI 10.1038/s41592-021-01358-2](https://doi.org/10.1038/s41592-021-01358-2). Squidpy.
- [Squidpy 1.6.0 co-occurrence source](https://github.com/scverse/squidpy/blob/v1.6.0/src/squidpy/gr/_ppatterns.py): cumulative radii and eligible-pair marginals; explicit thresholds are required for numeric comparisons.
- Traag, Waltman and van Eck (2019) Sci Rep 9:5233. Leiden.
