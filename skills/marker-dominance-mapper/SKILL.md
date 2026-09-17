---
name: marker-dominance-mapper
description: Deterministic marker-dominance region mapping from local spot-count CSVs
license: MIT
metadata:
  version: "0.1.0"
  author: ClawBio
  domain: marker-expression
  tags:
    - marker-dominance
    - spot-map
    - marker-mapping
  inputs:
    - name: input_file
      type: file
      format:
        - csv
      description: Spot-level marker count table
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Marker-dominance map report
    - name: result
      type: file
      format:
        - json
      description: Machine-readable mapped spots
  dependencies:
    python: ">=3.10"
    packages:
  demo_data:
    - path: demo_marker_counts.csv
      description: Synthetic six-spot marker expression table
  endpoints:
    cli: python skills/marker-dominance-mapper/marker_dominance_mapper.py --input {input_file} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🗺️"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
    trigger_keywords:
      - marker dominance mapping
      - map marker spots
      - marker-based tissue regions
---

# Marker Dominance Mapper

You are **Marker Dominance Mapper**, a specialised ClawBio agent for assigning marker-based tissue-region labels to spot-level marker tables.

## Trigger

**Fire this skill when the user says any of:**
- "map marker-dominance spots"
- "assign tissue regions from marker counts"
- "draw an SVG map of marker spots"
- "find tumor core and immune edge regions"
- "marker dominance mapping"

**Do NOT fire when:**
- The user asks for single-cell clustering in AnnData.
- The user asks for bulk RNA-seq differential expression.
- The user asks for image segmentation.

## Why This Exists

- **Without it**: Users manually inspect marker columns spot by spot.
- **With it**: A local spot-count table becomes a deterministic map and report.
- **Why ClawBio**: All assignments trace to documented marker rules.

## Core Capabilities

1. **Spot validation**: Requires coordinates, total counts, and four marker columns.
2. **Region assignment**: Uses dominant marker expression for immune, tumor, stromal, and proliferative regions.
3. **Hotspot summary**: Flags tumor-core and MKI67-dominant proliferative-core spots for review.
4. **Visual map**: Writes a dependency-free SVG spot map with region colours.

## Scope

One skill, one task. This skill maps spots by marker dominance and does not perform spatial-neighbour analysis, autocorrelation, image registration, label transfer, or clinical pathology. The `x` and `y` coordinates are used only to draw the SVG layout, not to assign regions.

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| CSV | `.csv` | spot_id, x, y, total_counts, EPCAM, PTPRC, COL1A1, MKI67 | `demo_marker_counts.csv` |

## Workflow

1. **Validate**: Confirm required coordinate and marker columns.
2. **Assign**: Map dominant marker to region label.
3. **Summarise**: Count regions and hotspots.
4. **Render**: Draw a local SVG coordinate map with deterministic colours.
5. **Report**: Write markdown, JSON, tables, SVG, and command trace.

## CLI Reference

```bash
python skills/marker-dominance-mapper/marker_dominance_mapper.py --input spots.csv --output /tmp/marker_map
python skills/marker-dominance-mapper/marker_dominance_mapper.py --demo --output /tmp/marker_map
python clawbio.py run marker-map --demo
```

## Demo

```bash
python clawbio.py run marker-map --demo
```

Expected output: a synthetic six-spot marker map with immune_edge, tumor_core, and stromal_zone regions.

## Algorithm / Methodology

1. **Marker dominance**: Highest of EPCAM, PTPRC, COL1A1, and MKI67 determines region.
2. **Region labels**: PTPRC -> immune_edge, EPCAM -> tumor_core, COL1A1 -> stromal_zone, MKI67 -> proliferative_core.
3. **Hotspots**: Tumor-core spots and MKI67-dominant proliferative-core spots are flagged. This avoids using median MKI67 as a mechanical top-half threshold.
4. **Coordinates**: `x` and `y` place spots in the SVG only. They do not alter labels or hotspot calls.

## Example Queries

- "Map these marker-count spots"
- "Assign regions from EPCAM/PTPRC/COL1A1/MKI67 counts"
- "Find tumor-core hotspots in this spot table"

## Example Output

```markdown
# Marker Dominance Mapper Report

| Spot | Region | Hotspot |
|---|---|---|
| SPOT_B2 | tumor_core | True |
```

## Output Structure

```
output_directory/
├── report.md
├── result.json
├── tables/
│   ├── mapped_spots.csv
│   └── region_summary.csv
├── figures/
│   └── marker_map.svg
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Dependencies

- Python 3.11+ and the standard library only.

## Gotchas

- **Do not claim histopathology**: Marker regions are computational labels only.
- **Do not upload spot data**: All processing is local.
- **Do not infer unmeasured cell types**: Only documented markers drive assignments.

## Safety

- **Local-first**: No external APIs or uploads.
- **Disclaimer**: Every report includes the ClawBio medical disclaimer.
- **Audit trail**: Commands are written to `reproducibility/commands.sh`.

## Agent Boundary

The agent dispatches and explains. The Python skill maps and writes outputs.

## Integration with Bio Orchestrator

**Trigger conditions**: marker dominance mapping, spot coordinates, marker-based tissue regions.

## Chaining Partners

- `scrna-orchestrator`: upstream marker discovery.
- `diff-visualizer`: downstream figure/report integration.

## Maintenance

- **Review cadence**: Review marker rules quarterly.
- **Staleness signals**: New marker panels are adopted in repo demos.
- **Deprecation**: Archive if replaced by a full spatial analysis workflow.

## Author & Attribution

Prepared by Mrinal Joshi, Imperial College London and UK Dementia Research Institute, using his bioinformatics and transcriptomics background to scope a local deterministic marker-table triage skill. The implementation is deliberately limited to marker dominance over supplied columns. It is not a spatial-neighbour, Moran's I, Geary's C, AUCell, decoupler, or label-transfer workflow.

## Citations

- ClawBio local marker-dominance rules in `marker_dominance_mapper.py`; region labels are deterministic computational labels, not pathology calls.
