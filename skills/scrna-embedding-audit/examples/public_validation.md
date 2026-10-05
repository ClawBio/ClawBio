# Public-data validation: Wilk et al. 2020 PBMC (CELLxGENE)

Date: 2026-10-05. Skill version 0.1.0 (branch `feat/scrna-embedding-audit`). This record shows the
skill running end to end on measured data with two embeddings produced by documented scanpy
commands. It validates input handling, run time and the shape of the verdicts; it is not a
statement about which integration method is best.

## Dataset

- Wilk et al., "A single-cell atlas of the peripheral immune response in patients with severe
  COVID-19", Nature Medicine 2020, doi:10.1038/s41591-020-0944-y. Raw data GEO GSE150728.
- CELLxGENE Discover collection `a72afd53-ab92-4511-88da-252fb0e26b9a`, dataset
  `456e8b9b-f872-488b-871d-94534090a865`, dataset version `8d7fad46-2f37-4aa0-ac9d-babbc36efb6f`
  (schema 7.1.0). Distributed by CELLxGENE Discover under its open data policy (CC-BY 4.0 / CC0;
  the per-dataset `license` field in the API is null, so the licence is taken from the policy).
- 44,721 Seq-Well PBMCs, 24,505 genes (Ensembl IDs), 13 donors (`donor_id`: 7 COVID-19, 6
  healthy), 19 author cell types (`cell_type`). Raw UMI counts in `raw.X`; `X` is log-normalised.
- Labels come from the authors' annotation of their own clustering, so the circularity caveat
  printed by the skill applies to any embedding built from the same data.

Download (217,094,236 bytes), done by hand, not by the skill:

```bash
mkdir -p data && cd data
curl -L --fail --retry 3 -o wilk2020_covid_pbmc.h5ad \
  https://datasets.cellxgene.cziscience.com/8d7fad46-2f37-4aa0-ac9d-babbc36efb6f.h5ad
echo "50fa832f6ce824a18f1a5d54d2481ada5bee24036630a56bac7f7e98123a7c16  wilk2020_covid_pbmc.h5ad" | sha256sum -c
```

## Embeddings (scanpy only)

`examples/prepare_wilk2020_embeddings.py` builds an h5ad with raw counts in `X` and
`layers["counts"]`, `obs["cell_type"]`, `obs["donor_id"]`, and two embeddings, both on 2000
batch-aware Seurat HVGs with `scale(max_value=10)` and a 50-component arpack PCA (seed 0):

- `X_pca_hvg_batchaware`: no correction (batch information enters only through HVG selection).
- `X_combat`: `sc.pp.combat(key="donor_id")` before scaling and PCA.

```bash
uv run --locked --all-extras python skills/scrna-embedding-audit/examples/prepare_wilk2020_embeddings.py \
  --input data/wilk2020_covid_pbmc.h5ad --output data/wilk2020_for_audit.h5ad
```

Output file sha256: `ed177ecbbeec126f624fa923fe998af0a5f042afa6a781c55d17ee2d9fea8f1d` (44,721 cells x 24,505 genes; `obsm` keys `X_pca_hvg_batchaware`, `X_combat`)

## Audit command

```bash
OMP_NUM_THREADS=8 OPENBLAS_NUM_THREADS=8 \
uv run --locked --all-extras python skills/scrna-embedding-audit/scrna_embedding_audit.py \
  --input data/wilk2020_for_audit.h5ad --output audit/wilk2020 \
  --labels-key cell_type --batch-key donor_id --embeddings X_combat,X_pca_hvg_batchaware
```

Defaults applied: `--max-cells 5000` (one stratified subsample of the 44,721 cells is scored),
`--n-resamples 50`, `--lisi-neighbors 90`, `--graph-neighbors 15`, `--transfer-k 15`,
`--n-top-hvg 2000`, `--n-pcs 50`, `--random-state 0`.

## Environment

- Host: AMD Ryzen Threadripper PRO 7995WX (96 cores), 503 GB RAM, Linux 6.8.0; no GPU used.
- Python 3.12.13 via `uv sync --locked --all-extras`; scanpy 1.12.4, anndata 0.12.6,
  numpy 2.4.6, scipy 1.17.1, scikit-learn 1.9.1, pandas 3.0.6.

## Result

Run: 2026-10-05 23:22 UTC; wall time 1 min 15 s with 8 BLAS threads, max RSS 2.0 GB. 44721 cells in file, 5000 scored (one stratified subsample over 239 label x donor strata), 19 labels, 13 donors, no cells dropped, non-integer count fraction 0%. Baseline: 2000 batch-aware HVGs, 50 PCs, fitted on the 5,000 scored cells. Effective k: LISI 90 (perplexity 30), graph 15; half-samples of 2500 cells.

Summary lines printed by the skill:

- `X_combat`: biological conservation: 0 higher / 5 lower / 0 not separable / 0 tied of 5; batch mixing: 1 higher / 2 lower / 1 not separable / 0 tied of 4; cross-batch label transfer: 0 higher / 1 lower / 0 not separable / 0 tied of 1
- `X_pca_hvg_batchaware`: biological conservation: 2 higher / 0 lower / 3 not separable / 0 tied of 5; batch mixing: 1 higher / 3 lower / 0 not separable / 0 tied of 4; cross-batch label transfer: 1 higher / 0 lower / 0 not separable / 0 tied of 1

Verdicts (`tables/verdicts.csv`):

| embedding | metric | group | verdict | value | baseline | median diff | SD diff | frac > 0 |
|---|---|---|---|---:|---:|---:|---:|---:|
| `X_combat` | silhouette_label | biological conservation | lower | 0.534 | 0.555 | -0.0174 | 0.0011 | 0.00 |
| `X_combat` | isolated_labels | biological conservation | lower | 0.484 | 0.508 | -0.0405 | 0.0073 | 0.00 |
| `X_combat` | nmi_kmeans | biological conservation | lower | 0.589 | 0.655 | -0.0736 | 0.0157 | 0.00 |
| `X_combat` | ari_kmeans | biological conservation | lower | 0.404 | 0.492 | -0.1038 | 0.0383 | 0.00 |
| `X_combat` | clisi | biological conservation | lower | 0.988 | 0.990 | -0.0042 | 0.0010 | 0.00 |
| `X_combat` | bras | batch mixing | lower | 0.810 | 0.829 | -0.0112 | 0.0047 | 0.00 |
| `X_combat` | silhouette_batch | batch mixing | higher | 0.863 | 0.857 | +0.0153 | 0.0045 | 1.00 |
| `X_combat` | ilisi | batch mixing | lower | 0.179 | 0.214 | -0.0287 | 0.0050 | 0.00 |
| `X_combat` | graph_connectivity | batch mixing | not separable | 0.727 | 0.740 | +0.0075 | 0.0207 | 0.62 |
| `X_combat` | knn_transfer_accuracy | cross-batch label transfer | lower | 0.636 | 0.701 | -0.0555 | 0.0100 | 0.00 |
| `X_pca_hvg_batchaware` | silhouette_label | biological conservation | higher | 0.567 | 0.555 | +0.0115 | 0.0008 | 1.00 |
| `X_pca_hvg_batchaware` | isolated_labels | biological conservation | not separable | 0.519 | 0.508 | +0.0101 | 0.0086 | 0.84 |
| `X_pca_hvg_batchaware` | nmi_kmeans | biological conservation | not separable | 0.660 | 0.655 | +0.0083 | 0.0177 | 0.70 |
| `X_pca_hvg_batchaware` | ari_kmeans | biological conservation | not separable | 0.529 | 0.492 | +0.0123 | 0.0435 | 0.58 |
| `X_pca_hvg_batchaware` | clisi | biological conservation | higher | 0.995 | 0.990 | +0.0059 | 0.0007 | 1.00 |
| `X_pca_hvg_batchaware` | bras | batch mixing | lower | 0.783 | 0.829 | -0.0405 | 0.0033 | 0.00 |
| `X_pca_hvg_batchaware` | silhouette_batch | batch mixing | lower | 0.839 | 0.857 | -0.0133 | 0.0031 | 0.00 |
| `X_pca_hvg_batchaware` | ilisi | batch mixing | lower | 0.186 | 0.214 | -0.0220 | 0.0060 | 0.00 |
| `X_pca_hvg_batchaware` | graph_connectivity | batch mixing | higher | 0.838 | 0.740 | +0.0938 | 0.0205 | 1.00 |
| `X_pca_hvg_batchaware` | knn_transfer_accuracy | cross-batch label transfer | higher | 0.727 | 0.701 | +0.0368 | 0.0092 | 1.00 |

Full-data values including the baseline and the random control (`tables/metrics.csv`):

| embedding | silhouette_label | isolated_labels | nmi_kmeans | ari_kmeans | clisi | bras | silhouette_batch | ilisi | graph_connectivity | knn_transfer_accuracy |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| `X_combat` | 0.534 | 0.484 | 0.589 | 0.404 | 0.988 | 0.810 | 0.863 | 0.179 | 0.727 | 0.636 |
| `X_pca_hvg_batchaware` | 0.567 | 0.519 | 0.660 | 0.529 | 0.995 | 0.783 | 0.839 | 0.186 | 0.838 | 0.727 |
| `baseline:pca` | 0.555 | 0.508 | 0.655 | 0.492 | 0.990 | 0.829 | 0.857 | 0.214 | 0.740 | 0.701 |
| `control:random` | 0.492 | 0.493 | 0.012 | -0.000 | 0.753 | 0.967 | 0.920 | 0.409 | 0.070 | 0.045 |

## Reading the result

- `X_combat` scored lower than the in-run PCA baseline on all five biological-conservation metrics (median differences from -0.004 on cLISI to -0.104 on ARI), lower on BRAS and iLISI, higher on the legacy batch silhouette, not separable on graph connectivity, and lower on cross-batch label transfer. On this data ComBat by donor did not improve neighbourhood-level mixing (iLISI 0.179 vs 0.214) while it cost cell-type separation; `disease` is confounded with `donor_id` here (every donor is either COVID-19 or healthy) and cell-type composition differs strongly between donors, two situations in which a location-scale batch model is expected to remove biology. The verdicts describe this 5,000-cell sample and this ComBat run; they are not a statement about ComBat in general.
- `X_pca_hvg_batchaware` is the same recipe as the baseline except that it was fitted on all 44,721 cells while the baseline is fitted on the 5,000 scored cells. The differences are small (|median difference| <= 0.094, most below 0.02) but several reach `higher`/`lower` because half-sampling resolves small, consistent gaps. This is the expected behaviour of the rule and a useful reading of `--max-cells`: when an embedding was fitted on more cells than are scored, differences of this size against the baseline are attributable to the fit population, not to the method. The SKILL.md Gotchas say so.
- The random control behaves as designed: NMI 0.012, ARI 0.000, graph connectivity 0.070, label transfer 0.045, with iLISI 0.409 and BRAS 0.967 (noise mixes donors perfectly). It never enters a verdict.
- No cell was dropped, no label is confined to a single donor, no embedding had duplicate cells or degenerate LISI kernels, so no abstention fired on this input. The `--max-cells` warning is the only warning.
- Label provenance: the authors annotated clusters of their own Seq-Well data, so the circularity sentence printed in every report applies to any embedding built from the same matrix, including both scored here.

Replaying `audit/wilk2020/reproducibility/commands.sh` with `INPUT_PATH` set to the prepared
h5ad reproduces `result.json` up to floating-point noise (verified for the demo in the test suite;
the same bundle format is written here, with the input SHA-256 guard).
