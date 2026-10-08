# Expected output of `--demo`

`report.md` from `SCARF_MEM_BUDGET=2G SCARF_WORKERS=2 python skills/scarf-single-cell/scarf_single_cell.py --demo --output <dir>`
with scarf 1.0.0rc19 on Python 3.12. The date, output path and timings differ per run;
cell counts, clusters, markers and the adjusted Rand index are deterministic for this seed.

````markdown
# Scarf single-cell baseline report

**Date**: <run date>
**Skill**: scarf-single-cell 0.3.0
**Scarf**: 1.0.0rc19
**Run**: `clawbio_demo` (completed)
**Resources**: SCARF_MEM_BUDGET=2G, SCARF_WORKERS=2
**Input**: `demo_input.h5ad` (sha256 `469e9f9ce2cc5ae593a3b595b4c4a46e3ab2a8ed95d832d4ec71acb7fb35fc6c`)

---

> Demo mode: a seeded synthetic dataset (2,000 cells, 1,500 genes, four planted
> groups, two donors). Numbers below describe that toy matrix only.

## Summary

- Cells: 1,997 of 2,000 kept by the MAD 5 filter
- Clusters (pipeline silhouette pick): 4
- Handoff H5AD and ClawBio's raw-count guard: accepted

## Matrix check (scripts/inspect_store.py)

- No sign of corrected or pre-filtered counts in these checks.

## QC retention by `donor_id`

| group | n_cells | n_kept | frac_kept |
|---|---|---|---|
| D1 | 1000 | 999 | 0.999 |
| D2 | 1000 | 998 | 0.998 |

## Cluster sizes

| cluster | n_cells |
|---|---|
| 1 | 523 |
| 2 | 518 |
| 3 | 496 |
| 4 | 460 |

## Top markers

- Cluster 1: GENE166, GENE169, GENE163, GENE190, GENE192
- Cluster 2: GENE145, GENE117, GENE108, GENE128, GENE119
- Cluster 3: GENE94, GENE66, GENE77, GENE90, GENE64
- Cluster 4: GENE17, GENE9, GENE18, GENE2, GENE14

## Composition by `donor_id`

| cluster | D1 | D2 |
|---|---|---|
| 1 | 267 | 256 |
| 2 | 259 | 259 |
| 3 | 242 | 254 |
| 4 | 231 | 229 |

## Handoff

`handoff/run_cells.h5ad` holds the run's cells with raw counts in `X`, Scarf `clusters`,
Leiden and Paris labels in `obs`, and `X_umap`. Continue in Scanpy or scVI with:

```bash
python clawbio.py run scrna --input <dir>/handoff/run_cells.h5ad --output <dir>
python clawbio.py run scrna-embedding --input <dir>/handoff/run_cells.h5ad --output <dir>
```

Those skills recompute QC, normalisation and clustering; Scarf labels stay in `obs`.

## Held-out comparison: `author_cell_type`

Read only after the analysis above was finished; nothing upstream used it.

- Adjusted Rand index (clusters vs held-out labels): 1.0
- Crosstab: `tables/holdout_crosstab.csv`

## Caveats

- This is a baseline pass. Audit QC removals, cluster stability and markers before
  naming cell types (SKILL.md, Gotchas).
- Cells are not replicates: condition claims need donor-level aggregation.
- Run details: `run_report.md`; provenance: `lineage.md`, `handoff.json`.

---

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
````
