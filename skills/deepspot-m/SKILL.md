---
name: deepspot-m
description: Transcriptome-wide virtual spatial transcriptomics from H&E histology with DeepSpot-M. Scores a 224x224 tile and returns per-gene log1p-CPM values for any HGNC symbols you ask for, with a CSV, a report and a reproducibility bundle.
license: MIT
model_license: cc-by-nc-sa-4.0
metadata:
  version: "0.2.0"
  author: Kalin Nonchev
  domain: spatial-transcriptomics
  tags:
    - spatial-transcriptomics
    - histology
    - gene-expression
    - foundation-model
    - digital-pathology
    - h-and-e
  inputs:
    - name: input_file
      type: file
      format:
        - png
        - jpg
        - jpeg
        - tif
        - tiff
      description: One 224x224 H&E tile cut at native (~20x) resolution
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Per-gene expression report with the upstream limitations attached
    - name: result
      type: file
      format:
        - json
      description: Machine-readable per-gene log1p-CPM values and run parameters
    - name: tables
      type: file
      format:
        - csv
      description: Gene table, one row per gene
    - name: reproducibility
      type: directory
      format:
        - dir
      description: commands.sh, environment.yml and checksums.sha256
  dependencies:
    python: ">=3.11"
    packages:
      - deepspotm>=1.0,<2
      - Pillow>=9.0
  demo_data:
    - path: examples/demo_tile.png
      description: Synthetic 224x224 H&E-like tile
    - path: examples/demo_expression.json
      description: Offline fixture standing in for a gene panel readout
  endpoints:
    cli: python skills/deepspot-m/deepspot_m.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧬"
    homepage: https://github.com/ratschlab/DeepSpotM
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: deepspotm
    trigger_keywords:
      - virtual spatial transcriptomics
      - gene expression from histology
      - spatial transcriptomics from H&E
      - predict gene expression from a tissue image
      - DeepSpot-M
---

# 🧬 DeepSpot-M Virtual Spatial Transcriptomics

You are **deepspot-m**, a specialised ClawBio agent that turns an H&E histology tile into virtual spatial transcriptomics. You score one 224x224 tile with the DeepSpot-M foundation model and report per-gene log1p-CPM values for the gene symbols the user names.

## Trigger

**Fire this skill when the user says any of:**
- "virtual spatial transcriptomics"
- "predict gene expression from histology"
- "spatial transcriptomics from H&E"
- "what genes are expressed in this tissue image"
- "score this tile for BRAF and COL1A1"
- "run DeepSpot-M on this tile"
- "gene expression map from a slide"
- "H&E to transcriptome"

**Do NOT fire when:**
- The user wants cells counted or outlined in an image. That is `cell-detection`.
- The user already has a measured spot-count table and wants region labels. That is `marker-dominance-mapper`.
- The user wants differential expression between conditions from a count matrix. That is `rnaseq-de`.
- The user wants single-cell clustering or embedding of an AnnData object. That is `scrna-orchestrator` or `scrna-embedding`.
- The user asks for TCGA bulk expression lookups. That is `xena-tcga-gene-query`.

## Why This Exists

- **Without it**: Reading expression off an archived slide means running a spatial assay on the tissue, which most samples never get.
- **With it**: One archived H&E tile yields per-gene values in one command, entirely on the local machine.
- **Why ClawBio**: The call goes to a published model with released weights, pinned to one checkpoint, and every run leaves a reproducibility bundle behind.

This is a research tool, not a substitute for measurement. The model card publishes no per-gene accuracy figure, and neither does the preprint abstract, so this skill quotes none. Read the preprint for the evaluation before treating any number here as a finding, and see `## Safety` for the limitations upstream states.

## Core Capabilities

1. **Score a tile**: Map one 224x224 H&E tile to per-gene log1p-CPM values.
2. **Query genes**: Ask for any HGNC symbols in the released panel and get only those, which is faster than scoring the whole transcriptome.
3. **Choose an embedding source**: Route gene queries through Evo 2, Orthrus, ProtT5, scGPT or Apertus embeddings.
4. **Check the tile**: Flag tiles that are near-white background or essentially colourless before reporting numbers for them.
5. **Report**: Write `report.md`, `result.json`, a gene CSV and a reproducibility bundle.

## Scope

One skill, one task. This skill scores a single H&E tile and writes gene values. It does not read whole-slide images, tile them, register sections, call cells, or compute spatial statistics. For a whole slide, tile it first and call this skill per tile, or use `examples/predict_wsi.py` from the upstream repository.

## Input Formats

| Format | Extension | Required Properties | Example |
|--------|-----------|---------------------|---------|
| PNG | `.png` | Exactly 224x224 px, H&E stained | `examples/demo_tile.png` |
| JPEG | `.jpg`, `.jpeg` | Exactly 224x224 px, H&E stained | `tile.jpg` |
| TIFF | `.tif`, `.tiff` | Exactly 224x224 px, H&E stained | `tile.tif` |

Tiles must be exactly 224x224 pixels. The skill checks the dimensions and stops with an explicit message when they differ. Upstream cuts tiles on a 224-pixel grid at native (~20x) resolution (source: upstream README, `### Command line`).

**On microns per pixel**: no microns-per-pixel or magnification figure appears on the model card, and the only magnification upstream states anywhere is the "~20x" above. So the skill never assumes a pixel size. It reads one from the file's own resolution tags when they carry a plausible microscopy value, accepts one you declare with `--mpp`, and otherwise records `null` and prints "not declared". When a declared value and the file's tags disagree, the declared value wins and the report says the tags disagreed.

## Workflow

1. **Validate**: Confirm the tile is exactly 224x224 pixels and load it.
2. **Check the tile**: Measure mean pixel value and mean saturation. Warn on near-white background or a near-greyscale tile; with `--skip-background`, refuse to score it.
3. **Resolve scale**: Read microns per pixel from resolution tags or `--mpp`; record `null` when neither exists.
4. **Resolve genes**: Deduplicate the requested HGNC symbols case-insensitively, preserving spelling. With no `--genes` flag, use the bundled ten gene marker panel.
5. **Load model**: Call `DeepSpotM.from_pretrained("ratschlab/DeepSpotM", source=..., revision=...)` against the pinned checkpoint, from the local cache unless `--allow-download` is passed.
6. **Match to the panel**: Case-fold the requested symbols against `model.gene_names` and carry forward the panel's own spelling.
7. **Predict**: Run `model.predict_genes(image_processor(tile).unsqueeze(0), genes)`.
8. **Report**: Write `report.md`, `result.json`, `tables/gene_expression.csv` and the reproducibility bundle, in requested gene order.

Steps 1, 5, 6 and 7 are prescriptive. Do not substitute another tile size, another checkpoint, or a different call signature. Step 8 narrative is open to the agent.

## CLI Reference

```bash
# Standard usage
python skills/deepspot-m/deepspot_m.py \
  --input tile.png --output /tmp/deepspot_out

# Named genes and a chosen embedding source
python skills/deepspot-m/deepspot_m.py \
  --input tile.png --genes BRAF,CD37,COL1A1 --source evo2 --output /tmp/deepspot_out

# Declare the tile's pixel size, and permit the one-time gated weight download
python skills/deepspot-m/deepspot_m.py \
  --input tile.tif --mpp 0.5 --allow-download --output /tmp/deepspot_out

# Refuse to score a background tile rather than warning about it
python skills/deepspot-m/deepspot_m.py \
  --input tile.png --skip-background --output /tmp/deepspot_out

# Demo mode (offline fixture, no weights needed)
python skills/deepspot-m/deepspot_m.py --demo --output /tmp/deepspot_demo

# Via the ClawBio runner
python clawbio.py run deepspot-m --input tile.png --genes BRAF,CD37
python clawbio.py run deepspot-m --demo
```

| Flag | Default | Purpose |
|------|---------|---------|
| `--genes` | 10 gene marker panel | Comma separated HGNC symbols to score |
| `--source` | `scgpt` | Frozen gene embedding space |
| `--mpp` | unset | Declared microns per pixel; recorded, never assumed |
| `--skip-background` | off | Refuse rather than warn when a tile fails the checks |
| `--white-mean` | `220` | Mean pixel above which a tile counts as background (upstream's default) |
| `--min-saturation` | `0.05` | Mean HSV saturation below which a tile is flagged as not H&E |
| `--allow-download` | off | Permit the one-time gated weight fetch from Hugging Face |

## Demo

```bash
python clawbio.py run deepspot-m --demo
```

Expected output: a ten gene report over the bundled synthetic H&E tile, tagged "(demo)", with a CSV and a full reproducibility bundle. Demo mode reads `examples/demo_expression.json` instead of the model, so it runs with no weights, no GPU and no network.

## Algorithm / Methodology

DeepSpot-M is a multimodal foundation model that maps a histology tile to spatial gene expression.

1. **Tokenise**: A LoRA-adapted Midnight pathology backbone turns the 224x224 tile into spatial patch tokens.
2. **Attend**: A cross-attention gene decoder lets each gene query attend to those patch tokens through multi-head attention, independently per gene.
3. **Route**: A gene router hypernetwork generates gene-specific output projections from frozen biological embeddings drawn from DNA, RNA, protein, single-cell and text foundation models (Evo 2, Orthrus, ProtT5, scGPT, Apertus).
4. **Emit**: Because genes are represented as queryable embeddings rather than fixed output slots, one model spans the protein-coding transcriptome, including genes it never saw during training.

**Key parameters**:
- Tile size: 224x224 px (source: DeepSpot-M model card, and upstream README)
- Magnification: native ~20x (source: upstream README, `### Command line`). No microns-per-pixel figure is published upstream.
- Output unit: log1p-CPM, the scale used by the TCGA virtual spatial transcriptomics atlas (source: atlas dataset card)
- Released panel: roughly 19,000 genes listed in `tokens.csv`, ordered by `model.gene_names` (source: upstream README)
- Embedding sources: `evo2`, `orthrus`, `prott5`, `scgpt`, `apertus`; default `scgpt`
- Pinned checkpoint: Hugging Face revision `86113ee431248c892d25cf55e1f8017cccec2926`

Applied to TCGA, the model produced a virtual spatial transcriptomics atlas of 28,664 slides across 32 cancer types. That atlas was generated with cancer-specific finetuned models. This skill pins the base checkpoint and runs it zero-shot, so it is not the configuration those numbers came from and should not be read as a description of your run.

## Example Queries

- "Run virtual spatial transcriptomics on this H&E tile"
- "What is the predicted EPCAM and PTPRC expression in this tile?"
- "Score tile.png for BRAF, CD37 and COL1A1 using the Evo 2 gene embeddings"

## Example Output

Verbatim `report.md` from `python skills/deepspot-m/deepspot_m.py --demo --output /tmp/deepspot_demo`. These are fixture values, which is why the run is tagged "(demo)" and the table is headed "Fixture Expression". A run against real weights differs only in the tag, the heading and the numbers.

```markdown
# DeepSpot-M Virtual Spatial Transcriptomics Report (demo)

**Date**: 2026-08-09 19:21 UTC
**Tile**: demo_tile.png
**Tile size**: 224x224 px, cut at native (~20x) resolution
**Microns per pixel**: not declared (pass --mpp, or use a tile whose resolution tags carry it)
**Model**: ratschlab/DeepSpotM @ 86113ee43124
**Gene embedding source**: scgpt
**Unit**: log1p-CPM
**Genes scored**: 10

> Demo mode. The values below come from the bundled offline fixture `examples/demo_expression.json`, not from a model run. They exist so the report format, the CSV schema and the reproducibility bundle can be inspected without the model weights.

## Fixture Expression

Genes appear in the order they were requested.

| Gene | Expression (log1p-CPM) |
|------|------------------------|
| EPCAM | 5.82 |
| KRT19 | 5.41 |
| COL1A1 | 4.97 |
| VIM | 4.63 |
| ACTA2 | 3.88 |
| PTPRC | 3.42 |
| CD68 | 2.91 |
| CD3D | 2.14 |
| CD8A | 1.76 |
| MKI67 | 1.35 |

## How to Read These Values

DeepSpot-M predicts relative expression, so a value means something next to the same gene in another tile, not next to a different gene in this one. Ordering the genes in this table by value would largely recover each gene's average abundance in the training data rather than anything specific to this fixture. `tables/gene_expression.csv` carries a `rank` column for convenience; it inherits that caveat.

Upstream states the following limitations, quoted from the "Limitations and biases" section of the model card:

- Trained on a finite set of cancer indications.
- Performance on unseen tissue types, stains, scanners or resolutions may degrade.
- Predicts relative expression rather than absolute counts.
- Under-sequenced genes are predicted less reliably.
- Trained on oncology cohorts, so it is not representative of healthy tissue or non-oncology contexts.
- Not for clinical or diagnostic use.

## Output Files

| File | Description |
|------|-------------|
| `result.json` | Machine-readable per-gene values and run parameters |
| `tables/gene_expression.csv` | Gene table, one row per gene |
| `reproducibility/commands.sh` | Exact command that produced this run |
| `reproducibility/environment.yml` | Conda and pip environment snapshot |
| `reproducibility/checksums.sha256` | SHA-256 digests of the outputs |

---

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md                      # Per-gene report with limitations attached
├── result.json                    # Per-gene values and run parameters
├── tables/
│   └── gene_expression.csv        # gene, expression, unit, rank
└── reproducibility/
    ├── commands.sh                # Exact command to reproduce
    ├── environment.yml            # conda-forge + nodefaults env snapshot
    └── checksums.sha256           # SHA-256 digests of the outputs
```

## Dependencies

**Required** (in `skills/deepspot-m/requirements.txt`, installed per skill rather than repo wide):
- `deepspotm` >= 1.0, < 2; the model, its loader and the image processor
- `Pillow` >= 9.0; tile loading, dimension checks and the tile quality checks

Installing `deepspotm` pulls in `torch`, `lightning`, `timm`, `peft`, `transformers`, `safetensors`, `huggingface_hub`, `pandas` and `numpy`. Both packages import lazily inside the prediction function, so the skill loads and runs its demo without them.

**Licensing and access**, stated plainly because it decides whether you may use this:
- This skill's own wrapper code is MIT. That grants nothing over the weights; the fields below are the ones that restrict you.
- Upstream code is [PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/). Non-commercial use only. You install `deepspotm` yourself and accept that directly; nothing from upstream is vendored here.
- Model weights are CC-BY-NC-SA-4.0. Non-commercial, ShareAlike, with attribution.
- **The NonCommercial term covers the outputs too.** Upstream's `WEIGHTS_LICENSE.md` applies it to "the weights or their outputs", so the numbers this skill writes are themselves non-commercial and require attribution. Real runs stamp that on `report.md` and `result.json`; demo runs do not, because fixture values never touched the weights.
- ShareAlike bites if you fine-tune: derived weights must be redistributed under CC-BY-NC-SA-4.0. This skill only runs inference, so it does not trigger that.
- The restriction comes from DeepSpot-M itself, not its parts. Per upstream's `THIRD_PARTY_LICENSES.md`, the Midnight backbone and all five gene-embedding sources (Evo 2, Orthrus, ProtT5, Apertus, scGPT) are MIT or Apache-2.0.
- Weights are gated on Hugging Face with **manual approval**. Request access on the model page, wait for a human to grant or refuse it, then run `huggingface-cli login`. Approval is not guaranteed and access can be declined.
- **The gate terms are narrower than the licence alone.** Access is granted "only to individuals whose affiliations are exclusively academic or public non-profit research institutions". A concurrent commercial affiliation — employment, consulting, advisory roles, internships or founding roles at a company or startup — makes you ineligible, and research performed at, for, funded by or in collaboration with a commercial entity counts as commercial use. Internal evaluation, benchmarking and proof-of-concept work in a commercial setting are covered by that, so check your own affiliation before requesting access.
- The skill loads from the local Hugging Face cache by default. It resolves `config.json`, `model.safetensors` and `tokens.csv` itself, passing `local_files_only=True` to huggingface_hub, and hands upstream the resulting directory rather than the repo id. The first fetch needs `--allow-download`; nothing reaches the network without it.
- Nothing from upstream is vendored here. ClawBio ships a wrapper; you install `deepspotm` yourself and accept its terms directly.

## Gotchas

- **You will want to compare two genes in the same tile. Do not.** The model predicts relative expression, not absolute counts, so EPCAM scoring above COL1A1 in one tile mostly reflects EPCAM's higher average abundance in the training data. Compare one gene across tiles instead. The report leads with requested order for this reason, and `rank` in the CSV inherits the caveat.
- **You will want to feed a whole slide or an arbitrary crop. Do not.** The model reads exactly 224x224 pixels. A 256x256 crop, a 40x tile or a downsampled thumbnail changes the effective field of view and the prediction with it. Tile on a 224-pixel grid at native 20x resolution first.
- **You will want to upper-case gene symbols. Do not.** HGNC keeps `orf` lower case in roughly 200 symbols, so `C9ORF72` is not in the panel and `C9orf72` is. Pass symbols as HGNC writes them; the skill case-folds to look up and reports the panel's own spelling either way.
- **You will want to run `--demo` and quote the numbers. Do not.** Demo mode reads `examples/demo_expression.json`, an offline fixture that exists to show the report format without the gated weights. The report is tagged "(demo)" and the table is headed "Fixture Expression" for exactly this reason.
- **You will want to score any image you have. Do not.** A blank, non-H&E or non-oncology tile is still scored and still returns numbers. The skill warns on near-white and near-greyscale tiles, but it cannot tell healthy tissue from tumour, and upstream trained on oncology cohorts only.
- **You will want to ask for every gene at once. Do not, unless you need them.** `predict_genes` computes only the queries you pass, so a four gene request is much faster than the full panel.
- **You will want to treat `--source` as cosmetic. It is not.** The five embedding spaces are distinct frozen models, so the same tile scored under `evo2` and under `scgpt` gives different numbers. Record the source alongside the values, which `result.json` does for you.
- **Values are log1p-CPM, not raw counts.** Do not feed them into a tool that expects integer counts, and do not exponentiate them twice.

## Safety

**Upstream limitations, quoted verbatim from the "Limitations and biases" section of the [model card](https://huggingface.co/ratschlab/DeepSpotM):**

> Trained on a finite set of cancer indications. Performance on unseen tissue types, stains, scanners or resolutions may degrade. Predicts relative expression rather than absolute counts. Under-sequenced genes are predicted less reliably. Trained on oncology cohorts, so it is not representative of healthy tissue or non-oncology contexts. Not for clinical or diagnostic use.

Every report reproduces these, so they travel with the numbers rather than staying in this file.

- **Local-first**: Tiles are read from disk and scored on the local machine. Nothing is uploaded. The model loads from the local Hugging Face cache unless `--allow-download` is passed, which permits the one-time gated weight download and nothing else. The gate is enforced by passing `local_files_only` to huggingface_hub on each of the three checkpoint files, not by setting `HF_HUB_OFFLINE`, which the library reads once at import and would already have read by then.
- **Paths**: `report.md` and `result.json` record the tile's file name only, never the directory it came from, because those two files get forwarded. `reproducibility/commands.sh` keeps the full path, since replaying the run is the one thing that needs it. Scrub it before sharing a bundle from a patient directory.
- **Disclaimer**: Every report ends with the ClawBio disclaimer: *ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
- **Research use**: Upstream marks the model research use only, not for clinical or diagnostic use. This skill inherits that.
- **Audit trail**: Every run writes `reproducibility/commands.sh`, `environment.yml` and `checksums.sha256`, and pins the weight revision in `result.json`.
- **No hallucinated science**: Gene values come from the model. The skill never fills in a symbol it could not score, and never prints a pixel size it did not measure or receive.

## Agent Boundary

The agent dispatches, picks genes and explains. The Python skill validates the tile, calls the model and writes the outputs. The agent must not invent expression values, rescale the model output, relax the 224x224 check, report demo fixture numbers as a model run, assert a microns-per-pixel figure the run did not record, or read a cross-gene ordering as tile-specific biology.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on virtual spatial transcriptomics, gene expression from histology, H&E tiles, and named requests for DeepSpot-M.

## Chaining Partners

- `marker-dominance-mapper`: downstream. Per-tile marker values across a tiled slide give the spot table it maps into tissue regions.
- `diff-visualizer`: downstream. The gene CSV feeds heatmaps and dot plots.
- `cell-detection`: complementary. Segment the same tile for cell counts and morphology alongside the expression readout.

## Maintenance

- **Review cadence**: Check the model card and PyPI release each quarter.
- **Staleness signals**: A new `deepspotm` release, a changed `from_pretrained` signature, a new embedding source beyond the current five, an updated `tokens.csv` panel, a new Hugging Face revision, a published accuracy figure worth citing, or a change to the weight licence or gating.
- **Pinned revision**: `MODEL_REVISION` in `deepspot_m.py` pins the Hugging Face checkpoint. Bump it deliberately, re-read the limitations, and re-run the suite; never let it float.
- **Deprecation**: Archive to `skills/_deprecated/` if upstream withdraws the weights or the API diverges beyond a small wrapper fix.

## Citations

- [DeepSpot-M: a multimodal foundation model for transcriptome-wide virtual spatial transcriptomics from histology](https://doi.org/10.64898/2026.06.19.26356060); Nonchev, Dawo, Silina, Koelzer and Rätsch; medRxiv, posted 22 June 2026. Method, architecture and evaluation.
- [ratschlab/DeepSpotM](https://github.com/ratschlab/DeepSpotM); source code, PolyForm Noncommercial 1.0.0.
- [ratschlab/DeepSpotM on Hugging Face](https://huggingface.co/ratschlab/DeepSpotM); gated model weights, CC-BY-NC-SA-4.0.
