---
name: cell-detection
description: Cell segmentation in fluorescence microscopy images. Supports Cellpose/cpsam (Cellpose 4.0) with additional backends
  planned. Produces segmentation masks, per-cell morphology metrics (area, diameter, centroid, eccentricity), overlay figures,
  and a report.md.
license: MIT
metadata:
  version: 0.1.0
  author: ClawBio
  tags:
  - microscopy
  - segmentation
  - cellpose
  - fluorescence
  - imaging
  - cell-biology
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🔬
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: cellpose>=4.0
    - kind: pip
      package: tifffile
    - kind: pip
      package: scikit-image
    trigger_keywords:
    - cellpose
    - cpsam
    - cell segmentation
    - nucleus segmentation
    - fluorescence microscopy
    - microscopy image
    - image segmentation
    - cell counting
    - segmentation mask
---

# 🔬 Cell Segmentation

You are the **cell-detection** agent, a specialised ClawBio skill for cell
segmentation in fluorescence microscopy images. The default backend is `cpsam`
(Cellpose 4.0); additional backends (e.g. StarDist) are planned.

## Why This Exists

Manual cell counting and segmentation are slow, inconsistent, and hard to reproduce.

- **Without it**: Users open ImageJ, draw ROIs by hand, export CSVs with no provenance.
- **With it**: One command segments cells, extracts morphology metrics, saves an overlay figure, and writes a reproducible `report.md`.
- **Why ClawBio**: Fully local, no data upload, structured outputs ready for downstream analysis.

## Core Capabilities

1. **Segment**: Run `cpsam` on any TIFF, PNG, or JPG fluorescence image
2. **Measure**: Extract area, equivalent diameter, centroid, and eccentricity per cell
3. **Report**: Produce `report.md`, `{stem}_measurements.csv`, and histogram figures

## Input Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Greyscale TIFF | `.tif`, `.tiff` | H×W — passed directly |
| 2-channel TIFF | `.tif`, `.tiff` | H×W×2 — cytoplasm + nuclear, any order |
| 3-channel TIFF | `.tif`, `.tiff` | H×W×3 — H&E or fluorescence, any order |
| >3-channel TIFF | `.tif`, `.tiff` | First 3 channels used; remainder truncated with warning |
| PNG / JPEG | `.png`, `.jpg`, `.jpeg` | Greyscale or RGB |

**Channel handling:** cpsam is channel-order invariant — cytoplasm and nuclear channels can be in any order. You do not need to specify which channel is which. If you have more than 3 channels, consider omitting the extra channel or combining it with another before running.

## Workflow

1. **Load** image; detect greyscale vs multi-channel
2. **Prepare** — pass 1–3 channels through unchanged; truncate >3 to first 3 with a warning
3. **Segment** with `CellposeModel()` — no `channels` argument needed
4. **Metrics** via `skimage.measure.regionprops`
5. **Figures** — overlay + size distribution histogram
6. **Report** — `report.md` + `{stem}_measurements.csv` + reproducibility bundle (`commands.sh`, `environment.yml`, `checksums.sha256`)

## CLI Reference

```bash
# Standard usage — greyscale or multi-channel (cpsam handles channels automatically)
python skills/cell-detection/cell_detection.py \
  --input <image.tif> --output <report_dir>

# Override diameter estimate (pixels)
python skills/cell-detection/cell_detection.py \
  --input <image.tif> --diameter 30 --output <report_dir>

# Demo (synthetic image, no user file needed)
python skills/cell-detection/cell_detection.py --demo --output /tmp/cell_detection_demo
```

## Demo

```bash
python skills/cell-detection/cell_detection.py --demo --output /tmp/cell_detection_demo
```

Expected output: report.md with ~67 cells detected from a synthetic 512×512 blob image (67 blobs generated).

## Algorithm / Methodology

1. Load image with `tifffile` (TIFF) or `PIL` (PNG/JPG); detect ndim
2. If >3 channels, truncate to first 3 with a warning
3. Instantiate `CellposeModel(gpu=<flag>)`
4. Call `model.eval(img, diameter=<arg_or_None>)` — no `channels` arg (cpsam is channel-order invariant)
5. Extract per-cell stats from `masks` via `skimage.measure.regionprops`
6. Save `{stem}_measurements.csv`, figures, `report.md`

**Key parameters**:
- Model: `cpsam` (Cellpose 4.0 unified model — channel-order invariant)
- Channels: not passed — cpsam uses the first 3 channels of the input in any order
- Diameter: `None` triggers Cellpose auto-estimation

## Example Queries

- "Segment the cells in my DAPI image"
- "How many cells are in this microscopy image?"
- "Run cellpose on my TIFF and give me a cell count"
- "Segment my fluorescence image and export morphology metrics"

## Output Structure

```
output_dir/
├── report.md
├── {stem}_measurements.csv
├── {stem}_cp_masks.tif
├── {stem}_seg.npy
├── figures/
│   ├── {stem}_cp_outlines.png
│   └── {stem}_histogram.png
└── reproducibility/
    ├── checksums.sha256
    ├── commands.sh
    └── environment.yml
```

## Dependencies

- `cellpose>=4.0` — cpsam model
- `tifffile` — TIFF I/O
- `Pillow` — PNG/JPG loading
- `numpy` — array ops
- `matplotlib` — figures
- `scikit-image` — regionprops metrics

## Safety

- Local-first: no image data leaves the machine
- Every report includes the ClawBio medical disclaimer
- Reproducibility bundle (`commands.sh`, `environment.yml`, `checksums.sha256`) records the exact invocation, dependencies, and output integrity

## Integration with Bio Orchestrator

**Trigger conditions**:
- Input is a TIFF/PNG/JPG microscopy image
- User mentions "cellpose", "segment", "cell counting", "microscopy"

**Chaining partners**:
- Future: export ROI centroids to spatial transcriptomics workflows

## Citations

- [Pachitariu, Rariden & Stringer (2025) *Cellpose-SAM: superhuman generalization for cellular segmentation*. bioRxiv 2025.04.28.651001](https://doi.org/10.1101/2025.04.28.651001) — CellposeSAM / cpsam model
