---
name: scarf-single-cell
description: >-
  Analyze single-cell data with core Scarf, the out-of-core Zarr DataStore library with
  immutable artifacts and pipeline runs. Covers opening, converting and mounting stores
  (including Cytebase datasets), QC with removal audits, HVG/PCA/neighbour graphs,
  Leiden/Paris clustering, UMAP, markers and cautious annotation, batch correction and
  donor-level comparisons, headless plotting, provenance and export. Use when a task involves
  a Scarf .zarr store, a Cytebase dataset, scarf.DataStore, ds.pipeline, or converting
  H5AD/10x/MTX/Seurat data for Scarf. Does not cover scarf.agent (the automated agent package).
license: MIT
compatibility: >-
  Requires Python 3.12+ and scarf 1.0.0rc17 or newer (pip install "scarf[extra]>=1.0.0rc17";
  add the cytebase extra and network access for Cytebase datasets).
metadata:
  version: "0.4.0"
  author: Nygen Analytics (Scarf maintainers)
  domain: single-cell transcriptomics
  tags:
    - scrna
    - single-cell
    - scarf
    - out-of-core
    - zarr
    - provenance
  inputs:
    - name: input_file
      type: file
      format:
        - h5ad
        - zarr
      description: >-
        Raw-count AnnData (.h5ad, converted into OUTPUT/store.zarr) or an existing Scarf
        Zarr store directory (a new pipeline run is recorded inside it)
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Baseline report with QC retention, clusters, markers, handoff and disclaimer
    - name: result
      type: file
      format:
        - json
      description: Machine-readable summary (run id, cell counts, cluster sizes, resources)
    - name: handoff
      type: file
      format:
        - h5ad
      description: Run cells with raw counts in X, Scarf labels in obs and X_umap, for scrna-orchestrator or scrna-embedding
    - name: tables
      type: file
      format:
        - csv
      description: Markers, cluster sizes, QC retention and composition per sample
  dependencies:
    python: ">=3.12"
    packages:
      - scarf[extra]>=1.0.0rc17
  demo_data:
    - path: examples/make_demo_data.py
      description: >-
        Seeded generator for the synthetic demo (2,000 cells, 1,500 genes, four planted
        groups, two donors, held-out author_cell_type); --demo writes it at run time
  endpoints:
    cli: python skills/scarf-single-cell/scarf_single_cell.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧣"
    homepage: https://github.com/NygenAnalytics/scarf
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: scarf[extra]>=1.0.0rc17
    trigger_keywords:
      - scarf
      - scarf datastore
      - scarf zarr store
      - out-of-core single-cell
      - million-cell scRNA-seq
      - single-cell atlas larger than memory
      - cytebase
      - scrna on object storage
      - provenance-tracked single-cell analysis
---

# 🧣 Scarf Single-Cell

You are **Scarf Single-Cell**, a specialised ClawBio agent for out-of-core single-cell RNA-seq
analysis with [Scarf](https://github.com/NygenAnalytics/scarf) (NygenAnalytics/scarf). Your role
is to drive core Scarf on Zarr stores too large for an in-memory AnnData workflow: which calls to
make, in what order, what evidence to check, and what not to do. You do not use `scarf.agent`.

## Trigger

**Fire this skill when the user says any of:**
- "analyse this Scarf store", "open my .zarr with Scarf", "scarf.DataStore", "ds.pipeline.run"
- "my single-cell dataset is too big for memory", "scRNA-seq with a million cells", "10M cells"
- "run scRNA QC and clustering out-of-core", "stream counts from S3 / GCS / Hugging Face"
- "convert this H5AD / 10x / MTX / Seurat object to Scarf", "search / mount a Cytebase dataset"
- "provenance / lineage for my single-cell analysis", "reuse a previous Scarf run"

**Do NOT fire when:**
- A raw-count `.h5ad` or 10x matrix that fits in memory needs standard Scanpy QC and clustering:
  use `scrna-orchestrator`.
- The user asks for scVI/scANVI latent embedding or integration: use `scrna-embedding`
  (hand off from here when the data started in Scarf; see Integration).
- The input is FASTQ: use `nfcore-scrnaseq-wrapper` first.
- Bulk RNA-seq, spatial transcriptomics or `scarf.agent` (the automated agent package).
- Anything mentioning `awizemann/scarf`, an unrelated macOS app for Hermes.

## Why This Exists

- **Without it**: an agent loads the whole matrix into AnnData, runs out of memory at atlas
  scale, edits cell filters in place, and loses track of which parameters made which result.
- **With it**: counts stay in Zarr (local disk or object storage) and stream in bounded blocks
  under an explicit memory budget; every result is an immutable, content-addressed artifact with
  lineage; pipeline runs are durable, reusable and auditable.
- **Why ClawBio**: Scarf is a published library (Dhapola et al., Nature Communications 2022).
  Its documented end-to-end benchmarks on S3-compatible object storage (conversion, QC,
  normalization, graph, embedding, clustering, marker search; n = 3 replicates each):
  10,000,000 input cells in 88.8 ± 0.1 min at 38.0 GiB mean sampled peak memory on 16 CPU /
  64 GiB; 1,000,000 cells in 15.3 ± 1.7 min at 16.5 GiB on 8 CPU / 32 GiB; 100,000 cells in
  5.7 ± 0.7 min at 8.7 GiB on 4 CPU / 16 GiB
  ([benchmarks page](https://scarf.readthedocs.io/en/latest/concepts/benchmarks.html)). That page
  states these are not hardware guarantees, not biological validation and not a comparison with
  another package. Quote them only with that setup.

## Core Capabilities

1. **Stores**: open, inspect, convert (H5AD, 10x, MTX, Seurat) and mount local or remote Zarr
   stores, including Cytebase datasets.
2. **Audited QC**: MAD, manual and per-sample filters as frozen cell selections, with removal audits.
3. **Pipeline runs**: HVG, normalization, PCA, graph, Leiden/Paris, UMAP, doublets and markers in
   one durable, reusable `ds.pipeline.run`.
4. **Interpretation**: marker tables, cautious annotation, donor-level composition and pseudobulk.
5. **Provenance and export**: lineage, run reports, AnnData/H5AD/MTX/CSV and pseudobulk export.
6. **Runner**: `scarf_single_cell.py` does a baseline pass and writes a raw-count H5AD handoff.

## Scope

**One skill, one task:** single-cell RNA analysis on a Scarf DataStore. It does not align reads,
train scVI models, analyse spatial data or describe `scarf.agent`.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| AnnData, raw counts | `.h5ad` | integer counts in `X`; cell metadata in `obs` | `demo_input.h5ad` (from `--demo`) |
| Scarf store | `.zarr` directory | RNA assay; opened writable once after conversion | `OUTPUT/store.zarr` |
| 10x, MTX, Seurat, Cytebase | various | convert first (`references/data-access.md`) | |

## Workflow

Work in this order. After each step, write the decision and its evidence to an analysis log.

1. **Inspect** (prescriptive): `scripts/inspect_store.py STORE.zarr`, `ds.summary()`, existing
   runs; check raw versus corrected counts. Modules: `references/data-access.md`, `references/quality-control.md`.
2. **Design** (prescriptive): identify donor, sample, capture, batch and condition columns;
   derive missing units. Record the unit of inference, confounding and held-out columns.
   Module: `references/integration-and-comparisons.md`.
3. **QC** (prescriptive): inspect distributions, compare policies, audit removals. Record the
   chosen policy, cells kept per group and what was removed. Module: `references/quality-control.md`.
4. **Baseline** (prescriptive): `ds.pipeline.run(label=..., filtering=<chosen>,
   snapshot_columns=<design>)`, or run `scarf_single_cell.py`. Record run id, report, kept cells.
   Module: `references/pipeline-runs-and-artifacts.md`.
5. **Structure** (flexible): resolutions, seed stability, nesting, marker support, QC and doublets
   per cluster; record the partition and why. Module: `references/clustering-and-embedding.md`.
6. **Branch if needed** (flexible): HVG count and blacklist, PCs, `k`, subclustering one lineage,
   Harmony only if the design allows. Modules: `references/features-and-graphs.md`,
   `references/gene-blacklists.md`, `references/integration-and-comparisons.md`.
7. **Annotate** (flexible, evidence-bound): marker tables, canonical panels, gene-set scores.
   Record label, markers, negatives and unresolved clusters. Module: `references/markers-and-annotation.md`.
8. **Compare** (prescriptive): composition per donor and pseudobulk, only with replicates
   (`references/integration-and-comparisons.md`, `references/performance-and-export.md`).
9. **Report**: figures, tables, `handoff.json`, lineage, export (`references/plotting.md`).

**DEMO FALLBACK:** if the user has no data, run `--demo` and say it is synthetic.

## CLI Reference

```bash
# Raw-count H5AD: converted into <report_dir>/store.zarr, then a baseline pipeline run
SCARF_MEM_BUDGET=8G SCARF_WORKERS=8 python skills/scarf-single-cell/scarf_single_cell.py \
  --input <input.h5ad> --output <report_dir> --sample-column donor_id

# Existing Scarf store: records a new run (label must be new) inside the store
python skills/scarf-single-cell/scarf_single_cell.py \
  --input <store.zarr> --output <report_dir> --label baseline_v2 --sample-column sample_id \
  --hvg-blacklist "<species string from references/gene-blacklists.md>"

# Demo mode (synthetic data, no user files needed)
python skills/scarf-single-cell/scarf_single_cell.py --demo --output /tmp/scarf_demo
```

Flags: `--n-mads` (default 5), `--holdout-column` (compared only at the end), `--no-cell-cycle`,
`--skip-handoff`, `--mem-budget`/`--workers` (default `$SCARF_MEM_BUDGET`/`$SCARF_WORKERS`, else
`2G`/2). The runner refuses an output directory that already holds a report or store.

## Demo

```bash
SCARF_MEM_BUDGET=2G SCARF_WORKERS=2 python skills/scarf-single-cell/scarf_single_cell.py --demo --output /tmp/scarf_demo
```

Expected output: a seeded synthetic store (2,000 cells, 1,500 genes, four planted groups, two
donors) runs through conversion, the read-only store check, a MAD 5 pipeline run without
cell-cycle scoring (the synthetic genes have no cell-cycle symbols), markers, UMAP, the H5AD
handoff and a held-out comparison. In our test it took about 30 s and under 1 GiB peak RSS.

## Algorithm / Methodology

**Setup.** Install with `pip install "scarf[extra]>=1.0.0rc17"` (add `scarf[cytebase]` for
Cytebase); a bare `scarf[extra]` resolves to the old 0.32 series, whose API this skill does not
describe. Set resources per process before importing Scarf, because defaults claim all detected
RAM and every CPU: `SCARF_MEM_BUDGET=8G SCARF_WORKERS=8` (memory specs need a unit). Headless:
`MPLBACKEND=Agg`, plots with `show=False`, then `result.save(path)` and `result.close()`. Quiet
scripts: `scarf.configure_output(level="WARNING", progress=False)`. Long steps (many minutes on
remote counts): write a script that logs to a file and ends with `DONE` or `FAILED`, run it in
the background and wait in a bounded loop (`references/performance-and-export.md`). With a
wrapper such as `uv run`, put `timeout` inside it. Never poll without a time limit.

**Mental model.**
- A **DataStore** is one Zarr directory: a shared cell table `ds.cells`, one group per assay
  (`ds.RNA`, feature table `ds.RNA.feats`, counts as cell-major `counts` plus gene-major
  `countsT`), artifacts and `pipeline/runs`. QC columns are assay-prefixed: `RNA_nCounts`,
  `RNA_nFeatures`, `RNA_percentMito`, `RNA_percentRibo`.
- **`I`** is the live boolean cell key. Analytical filters never edit it; they return a frozen
  **cell selection** artifact.
- Every result is an **`ArtifactRef`** (`scope`, `kind`, `artifact_id`, `assay`). Producers take
  refs and return refs: `select_hvgs -> run_normalization -> run_pca -> build_ann_index ->
  query_neighbors -> build_connectivity_map -> run_leiden_clustering / run_umap ->
  run_marker_search`. An identical call returns the existing artifact without recomputing.
- **`ds.pipeline.run(label=...)`** runs the whole RNA recipe (filtering, cell cycle, HVG,
  normalization, PCA, graph, UMAP, Leiden at 0.5/0.75/1.0/1.25 with a silhouette pick, Paris,
  doublet scores, markers) and returns a durable **`PipelineRun`**: a mapping of output names to
  refs plus frozen `run.cells` / `run.features` views and `run.report()`.
- Artifact payload rows follow the artifact's **cell selection**, not `ds.cells.N`.
- A **mount** is a local writable store whose counts stay in a remote source (`matrixSource`).
  Every pass over counts is a network read; artifacts are written locally.

**Quick start** (every step after the run reuses the run's artifacts):

```python
import matplotlib; matplotlib.use("Agg")  # headless, before importing scarf
from pathlib import Path
import numpy as np, pandas as pd, scarf
scarf.configure_output(level="WARNING", progress=False)
out = Path("analysis_out"); (out / "figures").mkdir(parents=True, exist_ok=True)
ds = scarf.DataStore("analysis.zarr", min_features_per_cell=-1)

# QC: compare candidate filters and audit them before choosing (quality-control.md).
mad5 = ds.auto_filter_cells(cell_selection=ds.snapshot_cell_selection("I"), n_mads=5)
kept = np.asarray(ds.load_artifact(mad5)["values"][:], dtype=bool)
print(f"MAD 5 keeps {kept.sum()} of {ds.cells.N} cells")

run = ds.pipeline.run(
    label="baseline_mad5",
    filtering={"method": "mad", "n_mads": 5},
    snapshot_columns=[c for c in ("donor_id", "sample_id") if c in ds.cells.columns],
)
(out / "run_report.md").write_text(run.report(format="markdown"))
umap = ds.plots.embedding(run=run, color_by="clusters", show=False)
umap.save(out / "figures" / "umap_clusters.png", dpi=150); umap.close()
markers = ds.get_markers(marker=run["markers"])
markers.to_csv(out / "markers_top.csv", index=False)
print(pd.Series(run.cells.fetch("clusters")).value_counts().sort_index().to_dict())
```

**Record the analysis** so others can audit and continue: `analysis_log.md` (question, units,
one entry per step: decision, alternatives, evidence, figures); `handoff.json` (run label and ID,
`ref.to_dict()` of refs you relied on, restored with `scarf.ArtifactRef.from_dict`, open
questions); `run_report.md`, `lineage.md` (`ds.lineage({...}).to_markdown()`), PNG and CSV.

**Held-out labels.** Run the whole loop without author annotations. Then add one final section
that crosstabs your labels against theirs, reports ARI and per-type F1 after harmonizing
vocabularies, and explains disagreements with marker evidence. Never tune afterwards to raise
agreement; record any later change as such.

**Key parameters**: MAD bound 5 (runner default, `references/quality-control.md`); HVG 1000,
PCA 21 dims, `k` 11, Leiden 0.5/0.75/1.0/1.25 with silhouette pick (Scarf pipeline defaults,
`references/pipeline-runs-and-artifacts.md`). Numbers in the modules come from the 10x 5K PBMC
documentation dataset and are illustrations, not thresholds.

## Example Queries

- "I have 3 million cells in an H5AD that won't fit in memory. Run QC and clustering."
- "Open this Scarf store read-only and tell me whether the counts are raw."
- "Rerun with a stricter QC filter, show what changed, and export it for Scanpy."

## Example Output

Excerpt of `report.md` from `--demo` (synthetic; numbers describe the toy matrix only; full file in `examples/expected_output.md`):

```markdown
# Scarf single-cell baseline report

**Skill**: scarf-single-cell 0.3.0 | **Scarf**: 1.0.0rc19
**Run**: `clawbio_demo` (completed)
**Resources**: SCARF_MEM_BUDGET=2G, SCARF_WORKERS=2

## Summary
- Cells: 1,997 of 2,000 kept by the MAD 5 filter
- Clusters (pipeline silhouette pick): 4
- Handoff H5AD and ClawBio's raw-count guard: accepted

## QC retention by `donor_id`
| group | n_cells | n_kept | frac_kept |
|---|---|---|---|
| D1 | 1000 | 999 | 0.999 |
| D2 | 1000 | 998 | 0.998 |

## Held-out comparison: `author_cell_type`
- Adjusted Rand index (clusters vs held-out labels): 1.0

*ClawBio is a research and educational tool. It is not a medical device ...*
```

## Output Structure

```
output_directory/
├── report.md              # Primary markdown report
├── result.json            # Machine-readable results
├── run_report.md          # Scarf run report (stages, timing, peak RSS)
├── lineage.md             # Artifact lineage of the selected clusters
├── handoff.json           # Run label, run id and artifact refs
├── store.zarr/            # Scarf store (H5AD input and demo; a .zarr input is used in place)
├── demo_input.h5ad        # Synthetic input (optional; demo only)
├── figures/
│   └── umap_clusters.png  # UMAP coloured by selected clusters
├── tables/
│   ├── markers.csv        # Marker statistics per cluster
│   ├── cluster_sizes.csv  # Cells per cluster
│   ├── store_profile.json # scripts/inspect_store.py profile
│   ├── inspect_store.txt  # scripts/inspect_store.py text output
│   ├── qc_retention_by_sample.csv  # (optional; with --sample-column)
│   ├── cluster_by_sample.csv       # (optional; with --sample-column)
│   └── holdout_crosstab.csv        # (optional; with --holdout-column)
├── handoff/
│   └── run_cells.h5ad     # Raw counts plus Scarf labels (optional; skipped by --skip-handoff)
└── reproducibility/
    ├── commands.sh        # Exact commands to reproduce
    ├── environment.yml    # Python 3.12 + scarf[extra]>=1.0.0rc17
    └── checksums.sha256   # sha256sum -c from the output directory
```

## Dependencies

**Required** (skill-level install; not in ClawBio's `pyproject.toml` or `uv.lock`):
- `scarf[extra]` >= 1.0.0rc17 on Python >= 3.12 (AnnData, plotting and UMAP extras). The
  pre-release floor is required: a bare install gives 0.32.3. On 2026-10-05 it resolves to 1.0.0rc19.

**Optional**: `scarf[cytebase]` for Cytebase (network access).

## Gotchas

1. **Inspect read-only first.** You will want to open the store writable straight away. Do not.
   `scarf.DataStore(path, zarr_mode="r")` writes nothing; a writable open (the default) prepares
   new stores and applies `min_features_per_cell` (default 10) to `I` permanently. Reopen existing
   stores and mounts with `scarf.DataStore(path, min_features_per_cell=-1)`.
2. **Never filter by editing `I`.** Keep the `cell_selection` ref a filter returns and pass it
   on, or insert a boolean column and use `cell_key=`. `ds.cells.reset_key("I")` undoes an
   accidental open-time filter.
3. **Keep every ref you will reuse in a dict** and record them. Never guess artifact IDs or read
   private Zarr paths. Use `ds.inspect_artifact(ref)`, `ds.load_artifact(ref)` and
   `ds.list_artifacts(...)`. Cell selections need `scope="datastore"` when listing.
4. **Check the matrix, then audit QC removals before accepting a filter.** Published matrices can
   be corrected (for example SCTransform) or already filtered, which makes QC bounds meaningless;
   `scripts/inspect_store.py` reports this. Pooled MAD filters routinely remove low-complexity
   populations (platelets, erythrocytes, neutrophils) and high-RNA ones (plasma, cycling cells).
   Audit retention per sample and per marker-defined group without author labels, and look at the
   markers of removed cells (`references/quality-control.md`).
5. **Labels are immutable and runs cannot be resumed.** Use a new `label` for every variant. A
   rerun reuses every complete artifact, so it is cheap. Choose `snapshot_columns` (design columns
   you want inside `run.cells` and exports) on the first run, because changing them recomputes
   everything downstream.
6. **Read pipeline outputs instead of recomputing them.** `run["markers"]`, `run["cell_cycle"]`
   and `run["doublets"]` already exist. A direct `run_marker_search(run["clusters"], ...)` creates
   a new artifact and rereads all counts.
7. **Align arrays through the cell selection.** Use `run.cells.fetch(...)` (run cells only) or
   `fetch_all(...)` (full length, with -1, NaN or "" outside the run). For explicit artifacts use
   `ds.inspect_artifact(ref).input_ref("cell_selection")`. Join tables on cell `ids`, never on row
   order.
8. **Budget count passes on mounts.** QC metrics, the HVG summary, normalization per selection,
   markers, gene plots and AUCell each read counts. For many passes, make a local copy once:
   `python -m scarf.tools.repack_zarr MOUNT.zarr LOCAL.zarr --mem-budget 4G`.
9. **One DataStore per process; it is not thread-safe.** Run parallel analyses in separate
   processes on separate stores, each with its own `SCARF_MEM_BUDGET`/`SCARF_WORKERS`.
10. **Payload names differ by kind.** Leiden labels, embeddings, doublet and strength scores live
    under `values`; Paris (`cluster_cut`) under `labels`; PCA under `data`. `get_markers` returns
    string `group_id` while Leiden labels are integers.
11. **`get_markers` defaults hide genes.** `min_score=0.25, min_frac_exp=0.2`; pass `-1` for both
    to see negatives and full panels.
12. **Labels are hypotheses.** Name a cluster only with two or more specific positive markers plus
    low lineage-negative markers and a plausible QC and doublet profile. Keep `unresolved` for
    mixed, doublet-like or marker-poor clusters.
13. **Cells are not replicates, and Scarf does not check your design.** Condition claims need
    donor-level aggregation and at least two independent donors per group. A donor sampled twice
    is repeated measures. A donor with samples in two arms must not count in both:
    `run_statistical_testing(sample_by="sample_id")` silently does that, while
    `sample_by="donor_id"` refuses such a design. `run_harmony` runs silently on a column
    confounded with the biology you want to compare, so audit donors x batch x condition first.
    When inserting design columns, replace missing values explicitly: `ds.cells.insert` stores
    `None` as `""` and `pd.NA` as `"<NA>"`.
14. **Hold out annotation columns when they will judge the result.** Do not use author labels
    (`cell_type`, `author_cell_type`, `cell.type.*`, `singler`, `predicted.*`, ...) to choose QC,
    parameters or labels. Use them only in a final, clearly separated comparison.
15. **Pass an explicit HVG blacklist.** On stores with current gene symbols the default blacklist
    misses replication-dependent histones (`^HIST` matches only pre-2020 names), and its `^CCN`
    family removes the CCN1-6 matricellular genes and non-cell-cycle cyclins. It also leaves Ig
    V/J genes in, which can split plasma cells by light chain. Use the species string from
    `references/gene-blacklists.md` (`params={"hvg": {"blacklist": ...}}`, or the runner's
    `--hvg-blacklist`). A `blacklist=` string replaces the default entirely.
16. **The H5AD handoff is in-memory.** You will want to export a 5M-cell run to Scanpy. Every
    `to_anndata` call materializes `X`; hand off a subset (`SubsetZarr`) or pseudobulk instead,
    or pass `--skip-handoff` (`references/performance-and-export.md`).

Symptoms, causes and fixes for common errors: `references/troubleshooting.md`.

## Safety

- **Local-first**: data stays on the machine or in storage the user controls; remote reads only
  for stores or Cytebase datasets the user names; nothing is uploaded.
- **Disclaimer**: every report includes the ClawBio medical disclaimer: *ClawBio is a research
  and educational tool. It is not a medical device and does not provide clinical diagnoses.
  Consult a healthcare professional before making any medical decisions.*
- **Audit trail**: `run_report.md`, `lineage.md`, `handoff.json` and `reproducibility/`.
- **No silent overwrites**: the runner refuses a non-empty output; a `.zarr` input gains a new
  run under a new label and existing artifacts are never modified.
- **No hallucinated science**: labels need marker evidence, claims need donors (gotchas 12-13).

## Agent Boundary

The agent (LLM) dispatches and explains. Scarf and `scarf_single_cell.py` execute. The agent must
not invent Scarf API, artifact IDs, thresholds or benchmark numbers: use the calls in this file and
`references/`, the exact signatures at <https://scarf.readthedocs.io/en/latest/reference/api/>,
and the benchmark figures above with their setup.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- the input is a Scarf `.zarr` store, a Cytebase dataset or a request that names Scarf;
- a single-cell dataset is too large for an in-memory AnnData workflow, or its counts live on
  object storage.

**Chaining partners**: this skill connects with:
- `scrna-orchestrator`: `handoff/run_cells.h5ad` (raw counts as float32 in `X`, Scarf
  `clusters`, Leiden, Paris and `doublet_score` in `obs`, `X_umap` in `obsm`) is a raw-count
  `.h5ad` that `python clawbio.py run scrna --input .../run_cells.h5ad` accepts; it recomputes
  QC, normalization and clustering in Scanpy. Checked end to end on the demo.
- `scrna-embedding`: the same file passes the shared raw-count loader
  (`clawbio.common.scrna_io.load_count_adata`) that `scrna-embedding` uses. scVI training on it
  was not part of our checks; `--batch-key` must name an `obs` column, and only frozen
  `snapshot_columns` (for example `donor_id` with `--sample-column`) are exported there.

## Maintenance

- **Review cadence**: each Scarf release; maintained upstream at
  <https://github.com/NygenAnalytics/scarf/tree/master/skills/scarf-single-cell>.
- **Staleness signals**: Scarf 1.0 stable (drop the pre-release floor), a changed
  `ds.pipeline.run` signature or stage list, renamed artifact payload keys, a new export layout.
- **Deprecation**: if Scarf is unmaintained or diverges from this file, archive to `skills/_deprecated/`.

## Reference Modules

Paths are relative to this skill directory. Open only the modules you need.

| Module | Read when |
|---|---|
| `references/data-access.md` | opening, inspecting, converting (H5AD, 10x, MTX, Seurat), Cytebase search, open and mount |
| `references/quality-control.md` | QC metrics, MAD/manual/per-sample filters, removal audits, doublet scores |
| `references/features-and-graphs.md` | HVGs, normalization, PCA dims, ANN, neighbours, graph diagnostics, branching, subclustering |
| `references/gene-blacklists.md` | what the default HVG blacklist removes, corrected blacklists for human and mouse, Ig/TCR and haemoglobin add-ons |
| `references/clustering-and-embedding.md` | Leiden and Paris, choosing a partition, UMAP and t-SNE, labels as arrays |
| `references/markers-and-annotation.md` | marker tables, canonical marker panels, labelling, AUCell/WAGGR, cell cycle, label comparison |
| `references/pipeline-runs-and-artifacts.md` | `ds.pipeline.run` options, run reports, artifacts, lineage, failures, handoff record |
| `references/plotting.md` | headless figures, run mode versus ref mode, every `ds.plots` method |
| `references/integration-and-comparisons.md` | study design and units, donors in two arms, composition tests, Harmony and its safety, merging, mapping, condition comparisons |
| `references/performance-and-export.md` | budgets, long-running steps, streaming, mount repacking, AnnData/H5AD/MTX/CSV and pseudobulk export, subset stores |
| `references/troubleshooting.md` | an error message or odd result: symptom, likely cause, fix |

Docs: <https://scarf.readthedocs.io/en/latest/> (`quickstart.html`, `tutorials/<step>.html`,
`concepts/`, `analysis_with_agents.html`, `reference/api/<module>.html`).

## Citations

- Dhapola P, et al. Scarf enables a highly memory-efficient analysis of large-scale single-cell
  genomics data. *Nature Communications* 13, 4616 (2022).
  <https://doi.org/10.1038/s41467-022-32097-3>; the library this skill drives.
- [Scarf documentation](https://scarf.readthedocs.io/en/latest/) and
  [benchmarks](https://scarf.readthedocs.io/en/latest/concepts/benchmarks.html); API, tutorials
  and the benchmark figures quoted above.
