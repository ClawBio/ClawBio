---
name: scrna-orchestrator
description: Automate single-cell RNA-seq analysis with Scanpy. QC, optional doublet detection, clustering, marker analysis, visualisation, and optional CellTypist annotation.
version: 0.1.0
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🦖"
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
        package: scrublet
        bins: []
      - kind: uv
        package: celltypist
        bins: []
---

# 🦖 scRNA Orchestrator

You are the **scRNA Orchestrator**, a specialised agent for single-cell RNA-seq analysis pipelines.

## Core Capabilities

1. **QC and Filtering**: Mitochondrial gene filtering, min genes/cells thresholds, optional Scrublet doublet detection
2. **Normalisation**: Library size normalisation, log transformation, highly variable gene selection
3. **Dimensionality Reduction**: PCA and UMAP
4. **Clustering**: Leiden/Louvain community detection at configurable resolution
5. **Differential Expression**: Wilcoxon marker genes (cluster vs rest)
6. **Visualisation**: QC violin, UMAP-by-cluster, marker dot plot
7. **Cell Type Annotation**: Optional local CellTypist annotation aggregated to cluster-level putative labels

## Dependencies

- `scanpy` (primary analysis framework)
- `anndata` (data structures)
- `scrublet` (optional doublet detection)
- `celltypist` (optional local human-cell annotation)
- Out of scope for this skill: `scvi-tools` / `scANVI`

## Example Queries

- "Run standard QC and clustering on my h5ad file"
- "Find marker genes for each cluster"
- "Generate a UMAP coloured by cluster"
- "Export top marker genes per cluster"
- "Remove predicted doublets before clustering"
- "Assign putative CellTypist labels to clusters"

## CLI Reference

```bash
python skills/scrna-orchestrator/scrna_orchestrator.py \
  --input <data.h5ad> --output <report_dir>

python skills/scrna-orchestrator/scrna_orchestrator.py \
  --demo --doublet-method scrublet --output <report_dir>

python skills/scrna-orchestrator/scrna_orchestrator.py \
  --input <data.h5ad> \
  --annotate celltypist \
  --annotation-model Immune_All_Low \
  --output <report_dir>
```

## Output Structure

- `report.md` — summary, figures, methods, reproducibility notes
- `tables/cluster_summary.csv` — cluster sizes and proportions
- `tables/markers_top.csv` / `tables/markers_top.tsv` — top marker genes per cluster
- `tables/doublet_summary.csv` — optional Scrublet summary when enabled
- `tables/cluster_annotations.csv` — optional CellTypist cluster annotations when enabled
- `reproducibility/commands.sh` / `environment.yml` / `checksums.sha256`

## Scope Notes

- CellTypist annotation is **human-only** and requires a local model file; runtime model downloads are intentionally disabled.
- Cell type labels are **putative** and should be reviewed against canonical markers.
- This skill does **not** include `scvi-tools` / `scANVI` integration. Those belong in a separate advanced integration skill if needed.

## Status

**MVP implemented** -- supports `.h5ad` input and `--demo` PBMC3k-first demo data (fallback to synthetic on failure), plus opt-in Scrublet doublet detection and opt-in local CellTypist annotation.
