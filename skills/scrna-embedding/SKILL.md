---
name: scrna-embedding
description: Local scVI-based single-cell latent embedding and batch-aware integration from raw-count .h5ad or 10x Matrix Market input, with integrated AnnData export and latent-space outputs.
version: 0.1.0
author: Yonghao Zhao
license: MIT
tags: [scrna, single-cell, scvi, embedding, integration, batch-correction, h5ad, 10x]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os: [macos, linux]
    install:
      - kind: uv
        package: scanpy
        bins: []
      - kind: uv
        package: anndata
        bins: []
      - kind: uv
        package: torch
        bins: []
      - kind: uv
        package: scvi-tools
        bins: []
    trigger_keywords:
      - scvi
      - embedding
      - latent
      - integration
      - batch correction
      - batch integration
      - h5ad
      - 10x
---

# 🧬 scRNA Embedding

You are **scRNA Embedding**, a specialised ClawBio agent for local single-cell latent embedding and batch-aware integration with scVI.

## Why This Exists

Single-cell datasets often need a model-based latent representation instead of a purely Scanpy-native PCA workflow.

- **Without it**: Users manually wire together scvi-tools training, latent export, clustering, and report generation.
- **With it**: One command trains scVI locally, writes `X_scvi`, clusters in latent space, exports marker tables, and saves an integrated `.h5ad`.
- **Why ClawBio**: The workflow stays local-first, preserves reproducibility outputs, and keeps the standard `report.md` / `result.json` contract.

## Core Capabilities

1. **Raw-count Input Validation**: Accept raw-count `.h5ad` and 10x Matrix Market input; reject processed-like matrices.
2. **scVI Latent Embedding**: Train `scvi.model.SCVI` with optional batch-aware integration.
3. **Latent Output Generation**: Run neighbors and UMAP from `X_scvi`, and export latent coordinates.
4. **Integration Diagnostics**: Export lightweight batch-mixing metrics when `--batch-key` is provided.
5. **Integrated Export**: Save `integrated.h5ad` with `obsm["X_scvi"]` and latent-space embeddings.
5. **Reproducibility Bundle**: Emit `commands.sh`, `environment.yml`, and checksums.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| AnnData raw counts | `.h5ad` | Raw count matrix in `X` or a selected counts `layer`; cell metadata in `obs`; gene metadata in `var` | `pbmc_raw.h5ad` |
| 10x Matrix Market | directory, `.mtx`, `.mtx.gz` | `matrix.mtx(.gz)` plus matching `barcodes.tsv(.gz)` and `features.tsv(.gz)` or `genes.tsv(.gz)` | `filtered_feature_bc_matrix/` |
| Demo mode | n/a | none | `python clawbio.py run scrna-embedding --demo` |

## Workflow

When the user asks for scVI embedding, latent integration, or batch correction:

1. **Validate**: Check raw-count `.h5ad` / 10x input (or `--demo`) and reject processed-like matrices.
2. **Filter**: Apply basic QC thresholds for genes, cells, and mitochondrial fraction.
3. **Train**: Fit `scvi.model.SCVI` on HVG raw counts, optionally using `--batch-key`.
4. **Project**: Export `X_scvi`, run latent-space neighbors and UMAP.
5. **Generate**: Write a minimal `report.md`, `result.json`, `integrated.h5ad`, latent tables, figures, and reproducibility files.

## CLI Reference

```bash
# Standard usage
python skills/scrna-embedding/scrna_embedding.py \
  --input <input.h5ad> --output <report_dir>

# Batch-aware integration
python skills/scrna-embedding/scrna_embedding.py \
  --input <input.h5ad> --output <report_dir> \
  --batch-key sample_id

# 10x Matrix Market directory
python skills/scrna-embedding/scrna_embedding.py \
  --input <filtered_feature_bc_matrix_dir> --output <report_dir>

# Demo mode
python skills/scrna-embedding/scrna_embedding.py \
  --demo --output <report_dir>

# Via ClawBio runner
python clawbio.py run scrna-embedding --input <input.h5ad> --output <report_dir>
python clawbio.py run scrna-embedding --demo
```

## Demo

```bash
python clawbio.py run scrna-embedding --demo
python clawbio.py run scrna-embedding --demo --batch-key demo_batch
```

Expected output:
- `report.md` with scVI-specific embedding and integration summary
- `integrated.h5ad` containing `obsm["X_scvi"]`
- figure files (`umap_scvi_latent.png`)
- optional batch figure (`umap_scvi_batch.png`) when `--batch-key` is set
- batch diagnostics table (`batch_mixing_metrics.csv`) when `--batch-key` is set
- latent export table (`latent_embeddings.csv`)
- reproducibility bundle

## Algorithm / Methodology

1. **QC**:
- Compute `n_genes_by_counts`, `total_counts`, `pct_counts_mt`
- Filter by `min_genes`, `min_cells`, `max_mt_pct`
2. **Feature selection**:
- Normalize + `log1p` on the full-gene branch
- Select HVGs (`flavor="seurat"`) for scVI training
3. **Latent model**:
- Train `scvi.model.SCVI` on raw-count HVGs
- Include batch covariate when `--batch-key` is provided
4. **Latent downstream analysis**:
- Save `obsm["X_scvi"]`
- Run neighbors with `use_rep="X_scvi"`
- Compute UMAP
- Export per-cell latent coordinates to CSV
5. **Batch diagnostics**:
- Compute lightweight mixing diagnostics from the neighbor graph and batch labels
- Report cross-batch neighbor fraction, neighbor entropy, and batch silhouette

## Example Queries

- "Run scVI on my h5ad file"
- "Integrate my batches with scvi-tools"
- "Build a latent embedding for this 10x matrix"
- "Export an integrated h5ad with X_scvi"

## Output Structure

```text
output_directory/
├── report.md
├── result.json
├── integrated.h5ad
├── figures/
│   ├── umap_scvi_latent.png
│   └── umap_scvi_batch.png           # only when batch integration is enabled
├── tables/
│   ├── latent_embeddings.csv
│   └── batch_mixing_metrics.csv      # only when batch integration is enabled
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

**Required**:
- `scanpy` >= 1.10
- `anndata` >= 0.10
- `torch`
- `scvi-tools`

**Out of scope (v1)**:
- `scANVI`
- `totalVI`
- multimodal integration
- condition-level DE
- remote model downloads

## Safety

- **Local-first**: No patient data upload.
- **Disclaimer**: Reports include the ClawBio medical disclaimer.
- **Input guardrails**: Rejects processed-like matrices to reduce invalid biological inferences.
- **No remote model fetches**: v1 uses only local code and local data.
- **Reproducibility**: Writes command/environment/checksum bundle.

## Integration with Bio Orchestrator

**Trigger conditions**:
- User explicitly asks for `scvi`, latent embedding, batch integration, or batch correction
- Input is single-cell data and the request is specifically model-based embedding rather than generic Scanpy clustering

**Routing note**:
- Generic single-cell clustering / marker requests still belong to `scrna-orchestrator`
- `scrna-embedding` is the advanced entry point for scVI-style latent integration and export

## Citations

- [scvi-tools documentation](https://docs.scvi-tools.org/) — model API and training interface.
- [Scanpy documentation](https://scanpy.readthedocs.io/) — downstream AnnData analysis utilities.
- [AnnData documentation](https://anndata.readthedocs.io/) — single-cell data model.
