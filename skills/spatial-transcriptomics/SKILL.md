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
      description: SpaceRanger outs/ directory or an h5ad with obsm['spatial']
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
| Spatial AnnData | `.h5ad` | `obsm['spatial']` with x,y per spot | `visium.h5ad` |
| Demo | n/a | none | `--demo` |

HDF5 `filtered_feature_bc_matrix.h5` is not read in v0.1; pass the mtx folder. Tissue images are not required.

## Workflow

1. **Validate**: Accept `outs/`, spatial h5ad, or `--demo`. Reject h5ad without `obsm['spatial']`. Abstain below 10 in-tissue spots.
2. **Process**: Scanpy QC (`min_genes`, `min_cells`), `normalize_total` to 1e4, `log1p`, HVGs, PCA, neighbours, UMAP, Leiden.
3. **Markers**: `rank_genes_groups(..., method="wilcoxon")`.
4. **Spatial graph**: k=6 nearest spots on `obsm['spatial']` (not the PCA graph).
5. **Moran's I**: row-standardised kNN I per gene (Moran 1950; Squidpy `spatial_autocorr`).
6. **Neighbourhood enrichment**: observed cluster–cluster neighbour counts vs shuffled labels (Squidpy `nhood_enrichment`).
7. **Co-occurrence**: P(cluster j at distance d from i) / P(j) in quantile distance bins (Squidpy `co_occurrence`).
8. **Generate**: `report.md`, `result.json`, figures, tables, `reproducibility/`.

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
```

| Flag | Default | Purpose |
|------|---------|---------|
| `--min-genes` | 5 | Drop spots with fewer genes |
| `--min-cells` | 1 | Drop genes in fewer spots |
| `--leiden-resolution` | 0.5 | Leiden resolution |
| `--random-state` | 7 | PCA / neighbours / Leiden / permutations |

## Demo

```bash
python clawbio.py run spatial --demo
```

Expected output: 64-spot synthetic grid, two spatial domains, Leiden ≥ 2, EPCAM/COL1A1 among high Moran's I genes, figures, tables, reproducibility bundle. No download.

## Algorithm / Methodology

1. **Load**: `scanpy.read_10x_mtx` plus `tissue_positions.csv` (or `tissue_positions_list.csv`); keep `in_tissue==1`.
2. **QC**: `calculate_qc_metrics`, `filter_cells(min_genes)`, `filter_genes(min_cells)`.
3. **Normalise**: `normalize_total(1e4)`, `log1p`. Raw counts kept in `layers["counts"]`.
4. **Embed**: Seurat HVGs, PCA, kNN on PCs, UMAP, Leiden (`flavor="igraph"` when Scanpy accepts it).
5. **Markers**: Wilcoxon, Benjamini–Hochberg adjusted p-values from Scanpy.
6. **Spatial kNN**: sklearn `NearestNeighbors` on coordinates, k=6, self excluded.
7. **Moran's I**: I = (zᵀWz)/(zᵀz) with row-standardised W.
8. **Enrichment**: 50 label permutations on the frozen spatial graph; z = (obs − mean_null) / sd_null.
9. **Co-occurrence**: pairwise Euclidean distances, 6 quantile bins, frequency-normalised.

**Key thresholds**:
- Minimum spots: 10 (below this the kNN graph is not meaningful)
- Spatial k: 6 (hex-like Visium neighbourhood)
- Leiden resolution: 0.5 (demo default; user-overridable)
- Permutations: 50 (speed; raise in a follow-up if a paper needs tighter tails)

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
| EPCAM | 0.82 |
| COL1A1 | 0.80 |
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── figures/
│   ├── umap_leiden.png
│   └── spatial_leiden.png
├── tables/
│   ├── markers_top.csv
│   ├── moran_i.csv
│   └── nhood_enrichment.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**:
- `scanpy` >= 1.10; QC, HVG, PCA, UMAP, Leiden, Wilcoxon
- `leidenalg` >= 0.10; Leiden
- `numpy`, `pandas`, `matplotlib`, `scikit-learn`, `scipy`; spatial graph, stats, figures

**Not required**:
- `squidpy`. Versions 1.4–1.6 depend on spatialdata, dask and s3fs. Those would land in every ClawBio CI job via `uv sync --all-extras`. The estimators above are the ones Squidpy documents; this skill computes them with numpy.

## Gotchas

- **You will want to route an H&E tile here. Do not.** This skill needs measured spot counts and coordinates. Predicted expression from histology is `deepspot-m`.
- **You will want to cluster on the spatial kNN graph. Do not, unless you mean it.** Leiden uses the PCA neighbour graph. The spatial kNN graph is only for Moran's I and enrichment. Mixing them silently changes what a cluster is.
- **You will want to treat Moran's I as a p-value. Do not.** v0.1 reports the statistic, not a permutation p-value per gene.
- **You will want to quote demo numbers as a tissue result. Do not.** `--demo` is an 8×8 synthetic grid.
- **You will want to pass a Visium HDF5 matrix. Do not in v0.1.** Supply the mtx `filtered_feature_bc_matrix/` directory.
- **You will want to threshold enrichment at |z|>1.96 as a discovery claim. Do not.** 50 permutations make the tails coarse; the diagonal sign is the supported reading.

## Safety

- **Local-first**: Spots are read from disk. No upload. `--demo` does not download a Visium dataset.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **No hallucinated science**: Cluster labels, Moran's I and enrichment come from the matrices above.
- **Audit trail**: `reproducibility/commands.sh`, `environment.yml`, `checksums.sha256` via `clawbio.common.reproducibility`.

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
- Palla et al. (2022) Nat Methods 19:171–178. PMID 35165433. Squidpy (estimators this skill reimplements).
- Traag, Waltman and van Eck (2019) Sci Rep 9:5233. Leiden.
