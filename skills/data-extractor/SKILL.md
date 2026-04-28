---
name: data-extractor
description: Extract numerical data from scientific figure images using Claude vision + OpenCV calibration. Supports 26+ plot types including bar charts, scatter plots, forest plots, Kaplan-Meier curves, box plots, and more.
version: 0.1.0
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env:
        - ANTHROPIC_API_KEY
      config: []
    always: false
    emoji: "📊"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: anthropic
        bins: []
      - kind: pip
        package: opencv-python-headless
        bins: []
      - kind: pip
        package: numpy
        bins: []
      - kind: pip
        package: Pillow
        bins: []
      - kind: pip
        package: pydantic
        bins: []
      - kind: pip
        package: fastapi
        bins: []
      - kind: pip
        package: uvicorn
        bins: []
---

# 📊 Data Extractor

You are the **Data Extractor**, a ClawBio skill for digitizing scientific figures. Your role is to extract numerical data from plot images for meta-analyses and systematic reviews.

## When to Use This Skill

Route to this skill when the user:
- Provides an image file (PNG, JPG, TIFF) containing a scientific figure
- Asks to "extract data from a figure", "digitize a plot", "read values from a chart"
- Mentions "meta-analysis data extraction" or "figure digitization"
- Wants to convert a bar chart, scatter plot, or other figure to CSV/JSON

## Capabilities

### Supported Plot Types (26)
scatter, bar, line, box, violin, histogram, heatmap, forest, kaplan_meier,
dot_strip, stacked_bar, funnel, roc, volcano, waterfall, bland_altman,
paired, bubble, area, dose_response, manhattan, correlation_matrix,
error_bar, table, other

### Pipeline (4 phases)
1. **Panel Detection** — Identify sub-panels in multi-panel figures (Claude vision)
2. **Pre-Analysis** — Identify axes, scale (linear/log), legend entries, error bars (Claude tool calling)
3. **CV Calibration + Extraction** — OpenCV detects markers/bars at pixel level, Claude extracts numerical data with calibration context
4. **Validation** — Heuristic checks for axis range, series count, error bar polarity

### Output Formats
- **CSV** — One row per data point with series name, x/y values, error bars
- **JSON** — Structured ExtractedData objects with full metadata
- **Web UI** — Interactive table + SVG preview with editable cells

## Usage

### CLI
```bash
python data_extractor.py --image figure.png --output results/
python data_extractor.py --web --port 8765
python data_extractor.py --demo
```

### API (importable)
```python
from data_extractor_api import run
result = run(options={"image_path": "figure.png", "output_dir": "results/"})
```

### Web UI
Launch with `--web` flag. Upload images, draw boxes around plots, extract and edit data interactively.

## Input Formats
- PNG, JPG, JPEG, TIFF image files
- Screenshots from papers, posters, slides
- Multi-panel composite figures (auto-detected and split)

## Notes
- Requires ANTHROPIC_API_KEY environment variable
- Uses Claude Sonnet for pre-analysis/detection, Claude Opus for extraction
- OpenCV calibration improves accuracy for scatter/bar plots with clear markers
- Error bars are reported as ± extent (delta from mean), not absolute positions
