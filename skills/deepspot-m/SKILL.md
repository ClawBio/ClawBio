---
name: deepspot-m
description: Transcriptome-wide virtual spatial transcriptomics from H&E histology with DeepSpot-M. Scores a 224x224 tile and returns per-gene log1p-CPM values for any HGNC symbols you ask for, with a ranked CSV, a report and a reproducibility bundle.
license: PolyForm-Noncommercial-1.0.0
metadata:
  version: "0.1.0"
  author: ClawBio Contributors
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
      description: One 224x224 H&E tile cut at roughly 20x (about 0.5 microns per pixel)
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Ranked gene expression report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable per-gene log1p-CPM values and run parameters
    - name: tables
      type: file
      format:
        - csv
      description: Ranked gene table, one row per gene
    - name: reproducibility
      type: directory
      format:
        - dir
      description: commands.sh, environment.yml and checksums.sha256
  dependencies:
    python: ">=3.11"
    packages:
      - deepspotm>=1.0
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

- **Without it**: Reading expression off a slide means sending tissue for a spatial assay, which costs weeks and thousands of pounds per sample.
- **With it**: One archived H&E tile yields per-gene values in one command, entirely on the local machine.
- **Why ClawBio**: The call goes to a published model with released weights, and every run leaves a reproducibility bundle behind.

## Core Capabilities

1. **Score a tile**: Map one 224x224 H&E tile to per-gene log1p-CPM values.
2. **Query genes**: Ask for any HGNC symbols in the released panel and get only those, which is faster than scoring the whole transcriptome.
3. **Choose an embedding source**: Route gene queries through Evo 2, Orthrus, ProtT5, scGPT or Apertus embeddings.
4. **Report**: Write `report.md`, `result.json`, a ranked CSV and a reproducibility bundle.

## Scope

One skill, one task. This skill scores a single H&E tile and writes gene values. It does not read whole-slide images, tile them, drop background, register sections, call cells, or compute spatial statistics. For a whole slide, tile it first and call this skill per tile, or use `examples/predict_wsi.py` from the upstream repository.

## Input Formats

| Format | Extension | Required Properties | Example |
|--------|-----------|---------------------|---------|
| PNG | `.png` | Exactly 224x224 px, H&E stained, roughly 20x | `examples/demo_tile.png` |
| JPEG | `.jpg`, `.jpeg` | Exactly 224x224 px, H&E stained, roughly 20x | `tile.jpg` |
| TIFF | `.tif`, `.tiff` | Exactly 224x224 px, H&E stained, roughly 20x | `tile.tif` |

Tiles must be exactly 224x224 pixels and cut at native 20x resolution, which is about 0.5 microns per pixel. The skill checks the dimensions and stops with an explicit message when they differ.

## Workflow

1. **Validate**: Confirm the tile is exactly 224x224 pixels and load it as RGB.
2. **Resolve genes**: Uppercase, deduplicate and order the requested HGNC symbols. With no `--genes` flag, use the bundled ten gene marker panel.
3. **Load model**: Call `DeepSpotM.from_pretrained("ratschlab/DeepSpotM", source=...)`, which returns the model and its image processor.
4. **Predict**: Run `model.predict_genes(image_processor(tile).unsqueeze(0), genes)` under `torch.no_grad()`.
5. **Rank**: Sort genes by descending log1p-CPM and attach a 1-based rank.
6. **Report**: Write `report.md`, `result.json`, `tables/gene_expression.csv` and the reproducibility bundle.

Steps 1, 3 and 4 are prescriptive. Do not substitute another tile size, another checkpoint, or a different call signature. Step 6 narrative is open to the agent.

## CLI Reference

```bash
# Standard usage
python skills/deepspot-m/deepspot_m.py \
  --input tile.png --output /tmp/deepspot_out

# Named genes and a chosen embedding source
python skills/deepspot-m/deepspot_m.py \
  --input tile.png --genes BRAF,CD37,COL1A1 --source evo2 --output /tmp/deepspot_out

# Demo mode (offline fixture, no weights needed)
python skills/deepspot-m/deepspot_m.py --demo --output /tmp/deepspot_demo

# Via the ClawBio runner
python clawbio.py run deepspot-m --input tile.png --genes BRAF,CD37
python clawbio.py run deepspot-m --demo
```

## Demo

```bash
python clawbio.py run deepspot-m --demo
```

Expected output: a ten gene report over the bundled synthetic H&E tile, tagged "(demo)", with a ranked CSV and a full reproducibility bundle. Demo mode reads `examples/demo_expression.json` instead of the model, so it runs with no weights, no GPU and no network.

## Algorithm / Methodology

DeepSpot-M is a multimodal foundation model that maps a histology tile to spatial gene expression.

1. **Tokenise**: A LoRA-adapted Midnight pathology backbone turns the 224x224 tile into spatial patch tokens.
2. **Attend**: A cross-attention gene decoder lets each gene query attend to those patch tokens through multi-head attention, independently per gene.
3. **Route**: A gene router hypernetwork generates gene-specific output projections from frozen biological embeddings drawn from DNA, RNA, protein, single-cell and text foundation models (Evo 2, Orthrus, ProtT5, scGPT, Apertus).
4. **Emit**: Because genes are represented as queryable embeddings rather than fixed output slots, one model spans the protein-coding transcriptome, including genes it never saw during training.

**Key parameters**:
- Tile size: 224x224 px (source: DeepSpot-M model card)
- Magnification: native 20x, about 0.5 microns per pixel (source: DeepSpot-M model card)
- Output unit: log1p-CPM, the scale used by the TCGA virtual spatial transcriptomics atlas (source: atlas dataset card)
- Released panel: roughly 19,000 genes listed in `tokens.csv`, ordered by `model.gene_names` (source: upstream README)
- Embedding sources: `evo2`, `orthrus`, `prott5`, `scgpt`, `apertus`; default `scgpt`

Applied to TCGA, the model produced a virtual spatial transcriptomics atlas of 28,664 slides across 32 cancer types.

## Example Queries

- "Run virtual spatial transcriptomics on this H&E tile"
- "What is the predicted EPCAM and PTPRC expression in this tile?"
- "Score tile.png for BRAF, CD37 and COL1A1 using the Evo 2 gene embeddings"

## Example Output

```markdown
# DeepSpot-M Virtual Spatial Transcriptomics Report

**Date**: 2026-08-04 09:14 UTC
**Tile**: tile.png
**Tile size**: 224x224 px at roughly 20x (0.5 microns per pixel)
**Model**: ratschlab/DeepSpotM
**Gene embedding source**: scgpt
**Unit**: log1p-CPM
**Genes scored**: 4

## Predicted Expression

| Rank | Gene | Expression (log1p-CPM) |
|------|------|------------------------|
| 1 | EPCAM | 5.82 |
| 2 | COL1A1 | 4.97 |
| 3 | PTPRC | 3.42 |
| 4 | MKI67 | 1.35 |

## Summary

EPCAM holds the highest value in this tile (5.82 log1p-CPM). Values sit on the
log1p-CPM scale, the same scale as the TCGA virtual spatial transcriptomics atlas
of 28,664 slides across 32 cancer types.

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
```

## Output Structure

```
output_directory/
├── report.md                      # Ranked gene report
├── result.json                    # Per-gene values and run parameters
├── tables/
│   └── gene_expression.csv        # rank, gene, expression, unit
└── reproducibility/
    ├── commands.sh                # Exact command to reproduce
    ├── environment.yml            # conda-forge + nodefaults env snapshot
    └── checksums.sha256           # SHA-256 digests of the outputs
```

## Dependencies

**Required** (in `skills/deepspot-m/requirements.txt`, installed per skill rather than repo wide):
- `deepspotm` >= 1.0; the model, its loader and the image processor
- `Pillow` >= 9.0; tile loading and dimension checks

Installing `deepspotm` pulls in `torch`, `lightning`, `timm`, `peft`, `transformers`, `safetensors`, `huggingface_hub`, `pandas` and `numpy`. Both packages import lazily inside the prediction function, so the skill loads and runs its demo without them.

**Licensing and access**, stated plainly because it decides whether you may use this:
- Upstream code is [PolyForm Noncommercial 1.0.0](https://polyformproject.org/licenses/noncommercial/1.0.0/). Non-commercial use only.
- Model weights are CC-BY-NC-SA-4.0. Non-commercial, ShareAlike, with attribution.
- Weights are gated on Hugging Face. Request access on the model page, then run `huggingface-cli login` once. After that the download is automatic.
- Nothing from upstream is vendored here. ClawBio ships a wrapper; you install `deepspotm` yourself and accept its terms directly.

## Gotchas

- **You will want to feed a whole slide or an arbitrary crop. Do not.** The model reads exactly 224x224 pixels at roughly 20x. A 256x256 crop, a 40x tile or a downsampled thumbnail changes the effective field of view and the prediction with it. Tile on a 224-pixel grid at native 20x resolution first.
- **You will want to run `--demo` and quote the numbers. Do not.** Demo mode reads `examples/demo_expression.json`, an offline fixture that exists to show the report format without the gated weights. The report is tagged "(demo)" for exactly this reason. Use `--input` for numbers worth citing.
- **You will want to ask for every gene at once. Do not, unless you need them.** `predict_genes` computes only the queries you pass, so a four gene request is much faster than the full panel. Pass `--genes` whenever you know what you are looking for.
- **You will want to treat `--source` as cosmetic. It is not.** The five embedding spaces are distinct frozen models, so the same tile scored under `evo2` and under `scgpt` gives different numbers. Record the source alongside the values, which `result.json` does for you.
- **The released checkpoint scores the roughly 19,000 gene panel in `tokens.csv`.** Symbols outside that panel need regenerated source gene embeddings, which are not part of the release. Use current HGNC symbols; a retired alias will not resolve.
- **Values are log1p-CPM, not raw counts.** Do not feed them into a tool that expects integer counts, and do not exponentiate them twice.

## Safety

- **Local-first**: Tiles are read from disk and scored on the local machine. Nothing is uploaded. The only network access is the one-time weight download from Hugging Face.
- **Disclaimer**: Every report ends with the ClawBio disclaimer: *ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
- **Research use**: Upstream marks the model research use only, not for clinical or diagnostic use. This skill inherits that.
- **Audit trail**: Every run writes `reproducibility/commands.sh`, `environment.yml` and `checksums.sha256`.
- **No hallucinated science**: Gene values come from the model. The skill never fills in a symbol it could not score.

## Agent Boundary

The agent dispatches, picks genes and explains. The Python skill validates the tile, calls the model and writes the outputs. The agent must not invent expression values, rescale the model output, relax the 224x224 check, or report demo fixture numbers as a model run.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here on virtual spatial transcriptomics, gene expression from histology, H&E tiles, and named requests for DeepSpot-M.

## Chaining Partners

- `marker-dominance-mapper`: downstream. Per-tile marker values across a tiled slide give the spot table it maps into tissue regions.
- `diff-visualizer`: downstream. The ranked CSV feeds heatmaps and dot plots.
- `cell-detection`: complementary. Segment the same tile for cell counts and morphology alongside the expression readout.

## Maintenance

- **Review cadence**: Check the model card and PyPI release each quarter.
- **Staleness signals**: A new `deepspotm` release, a changed `from_pretrained` signature, a new embedding source beyond the current five, an updated `tokens.csv` panel, or a change to the weight licence or gating.
- **Deprecation**: Archive to `skills/_deprecated/` if upstream withdraws the weights or the API diverges beyond a small wrapper fix.

## Citations

- [DeepSpot-M: a multimodal foundation model for transcriptome-wide virtual spatial transcriptomics from histology](https://doi.org/10.64898/2026.06.19.26356060); medRxiv, posted 22 June 2026. Method, architecture and evaluation.
- [ratschlab/DeepSpotM](https://github.com/ratschlab/DeepSpotM); source code, PolyForm Noncommercial 1.0.0.
- [ratschlab/DeepSpotM on Hugging Face](https://huggingface.co/ratschlab/DeepSpotM); gated model weights, CC-BY-NC-SA-4.0.
