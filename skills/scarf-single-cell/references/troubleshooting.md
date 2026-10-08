# Troubleshooting

Symptoms seen while driving core Scarf, their usual cause and the fix. Gotcha numbers refer
to `SKILL.md`.

| Symptom | Likely cause | Fix |
|---|---|---|
| `Assay 'RNA' is not prepared` on a read-only open | store just written by a converter or `SubsetZarr` | open once writable, then read-only |
| Fewer active cells than expected after opening | writable open applied `min_features_per_cell` | `ds.cells.reset_key("I")`; reopen with `-1` |
| `ValueError` on `ds.pipeline.run(label=...)` | label already completed | new label; reuse makes it cheap |
| `PermissionError` from a producer | store opened with `zarr_mode="r"` and no matching artifact | reopen writable |
| `KeyError: 'values'` on a Paris ref | Paris stores `labels` | `ds.load_artifact(ref)["labels"]` |
| Array length differs from `ds.cells.N` | payload follows the cell selection | align with `run.cells.fetch_all` or the selection mask |
| Plot raises with `run=` and a gene or live column | run mode accepts one frozen field only | `layout=run["umap"], color_by=[...]` |
| `TypeError` from `distribution(grouping="col")` | grouping needs a ref or `CellField` | `grouping=scarf.plotting.CellField("col")` |
| Very slow steps on a mount | each count pass is a network read | fewer passes; repack locally (gotcha 8) |
| `MemoryError` (`CountLayoutMemoryError` after 1.0.0rc19) from a converter | default count layout does not fit `mem_budget` | larger `mem_budget`; else the `policy=` the message names |
| `ValueError` plotting after reopening the store | a `PipelineRun` is bound to the store object that opened it | reopen the run from the new `ds` |
| `list_artifacts(kind="cell_selection")` is empty | cell selections are datastore-scoped | add `scope="datastore"` |
| `KeyError: 'groups'` in a dot plot table | with `group_by=` the column is named after the grouping column | read `res.tables["aggregate"].columns` first |
| A wait for a long step never ends | unbounded polling, or the process died | bounded wait that checks the process and the log tail |
| QC bounds look odd or nothing is filtered | counts are corrected or already filtered | check the matrix first; prefer flag-only or gentle filters |
| `None of the s_genes match` in the pipeline | feature names are not gene symbols (Ensembl IDs, synthetic names); matching ignores case, so mouse symbols work | `cell_cycle=False` (runner: `--no-cell-cycle`), or pass lists in the store's naming via `params={"cell_cycle": {"s_genes": [...], "g2m_genes": [...]}}` |
