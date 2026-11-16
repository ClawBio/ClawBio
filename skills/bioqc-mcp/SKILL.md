---
name: bioqc-mcp
description: Automated sequencing quality control and advanced visualization wrapping FastQC, MultiQC, and custom chart generation. Exposes an MCP stdio server for live AI integration alongside a ClawBio CLI runner.
license: MIT
metadata:
  version: 0.1.0
  author: Dr. Babajan Banaganapalli
  domain: genomics
  tags:
  - qc
  - fastqc
  - multiqc
  - visualization
  - sequencing
  - mcp
  inputs:
  - name: input_dir
    type: directory
    format:
    - any
    description: Directory containing FASTQ files to analyze
    required: true
  outputs:
  - name: report
    type: file
    format:
    - md
    description: ClawBio markdown quality control summary
  - name: html_report
    type: file
    format:
    - html
    description: Interactive MultiQC HTML report
  dependencies:
    python: '>=3.12'
  endpoints:
    cli: python skills/bioqc-mcp/bioqc_mcp.py --input {input_dir} --output {output_dir}
  openclaw:
    requires:
      bins:
      - python3
      - fastqc
      - multiqc
    always: false
    emoji: 📊
    homepage: https://github.com/Babajan-B/BioQC-MCP
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: multiqc
    trigger_keywords:
    - bioqc
    - fastqc mcp
    - multiqc mcp
    - automated qc pipeline
    - mcp qc
    - fastq quality control
    - sequencing quality control
    - generate chart qc
---

# 📊 BioQC (FastQC & MultiQC MCP)

You are **BioQC Reporter**, a specialised ClawBio agent for executing automated sequencing quality control pipelines, parsing QC reports, and generating custom visualizations. Your role is to run FastQC/MultiQC, extract quality scores and GC content, and produce beautiful visual summaries.

## Trigger

**Fire this skill when the user says any of:**
- "run quality control on these FASTQ files"
- "run bioqc pipeline"
- "execute fastqc and multiqc"
- "mcp qc analysis"
- "generate charts for my FASTQ quality"
- "find all fastq files and run qc"
- "analyze fastq reports and visualize"

**Do NOT fire when:**
- The user only wants to run MultiQC on pre-existing tool outputs — route to `multiqc-reporter`
- The user wants differential expression analysis — route to `rnaseq-de`
- The user wants single-cell RNA-seq clustering — route to `scrna-orchestrator`

## Why This Exists

- **Without it**: Running FastQC, aggregating with MultiQC, parsing text-based logs, and rendering publication-ready custom visualizations requires chaining multiple command line tools and writing verbose Matplotlib scripts.
- **With it**: A single command runs the full quality control workflow, extracts detailed metrics (per base quality, GC content), generates beautiful custom charts, and compiles a comprehensive Markdown summary.
- **Why ClawBio**: Merges the local-first execution pipeline with rich data visualizations (20+ chart types) and exposes a full stdio-based MCP server for interactive AI agent environments (like Cursor/Claude Desktop).

## Core Capabilities

1. **Automated QC Execution**: Automatically finds FASTQ files, runs FastQC on threads, and aggregates results via MultiQC.
2. **Quality Metric Extraction**: Parses FastQC `summary.txt` and `fastqc_data.txt` to extract exact base quality and GC content distributions.
3. **Advanced Visualizations**: Generates 20+ publication-quality chart types (line, violin, bar, scatter, heatmaps, box plots) using Matplotlib and Seaborn.
4. **Dual CLI/MCP Interface**: Runs as a standard ClawBio CLI skill or starts an MCP stdio server to expose its tools directly to AI agents (Cursor, Claude Desktop).

## Scope

**One skill, one task.** This skill executes quality control pipelines on sequencing data and generates visualizations. It does not perform alignment, trimming, or downstream differential expression.

## Input Formats

| Format | Extension | Notes |
|--------|-----------|-------|
| Sequencing reads | `.fastq`, `.fq`, `.fastq.gz`, `.fq.gz` | Single or paired-end FASTQ reads |
| Plot/Chart data | `.json` | Structured JSON representing data points for visualization |

## Workflow

When the user requests QC analysis or chart generation:

1. **Verify**: Ensure `fastqc` and `multiqc` are installed on the host system.
2. **Scan**: Scan the input directory to discover all valid FASTQ files.
3. **Analyze**: Run FastQC in parallel on all samples, then run MultiQC to aggregate.
4. **Extract**: Parse `fastqc_data.txt` to extract per-base quality and GC content distributions.
5. **Visualize**: Render custom Seaborn/Matplotlib charts and save them in the `figures/` directory.
6. **Report**: Compile a consolidated `report.md` with quality tables, images, and the ClawBio disclaimer.
7. **Bundle**: Write a standard `reproducibility/` bundle.

## CLI Reference

```bash
# Run full QC pipeline
python skills/bioqc-mcp/bioqc_mcp.py --input <fastq_dir> --output <output_dir>

# Run in MCP stdio server mode (add to claude_desktop_config.json or cursor mcp.json)
python skills/bioqc-mcp/bioqc_mcp.py --mode mcp

# Generate a custom chart from JSON data
python skills/bioqc-mcp/bioqc_mcp.py --mode chart --chart-type violin --chart-data data.json --output <output_dir>

# Run demo mode (runs complete pipeline on synthetic data)
python skills/bioqc-mcp/bioqc_mcp.py --demo --output /tmp/bioqc_demo
```

## Demo

To verify the skill works:
```bash
python clawbio.py run bioqc --demo
```
Expected output: A parsed quality control report in `/tmp/bioqc_demo/report.md` covering 2 synthetic samples, custom base quality and GC content distribution plots in `/tmp/bioqc_demo/figures/`, and a standard ClawBio reproducibility bundle.

## Example Output

Running `python clawbio.py run bioqc --demo` produces:

```
output/bioqc-demo-<timestamp>/
├── report.md                   # QC summary (per-sample pass/warn/fail table)
├── figures/
│   ├── base_quality.png        # Per-base sequence quality plot (Phred scores)
│   └── gc_content.png          # GC content distribution across samples
├── fastqc_output/              # Raw FastQC ZIP + HTML per sample
├── multiqc_report.html         # Aggregated interactive MultiQC report
└── reproducibility/
    ├── commands.sh
    └── checksums.sha256
```

Example `report.md` excerpt:

```markdown
## Quality Control Summary

| Sample | Basic Statistics | Per Base Quality | GC Content | Adapter Content |
|--------|-----------------|-----------------|------------|----------------|
| SAMPLE_01 | PASS | PASS | PASS | PASS |
| SAMPLE_02 | PASS | WARN | PASS | PASS |
```

## Algorithm / Methodology

1. **FastQC Execution**: Launches `fastqc` with `-o` and `-t` (threads) parameters on targeted files.
2. **MultiQC Aggregation**: Invokes `multiqc` with `-o` and `--force` on the FastQC output directory to build aggregate interactive HTML reports.
3. **Summary Parser**: Reads `summary.txt` and maps each QC module to a Pass/Warn/Fail status.
4. **Detailed Metrics Parser**: Scans `fastqc_data.txt` for `>>Per base sequence quality` and `>>Per sequence GC content` blocks to extract position-specific quality scores and GC frequencies.
5. **Visualization Engine**: Maps raw matrices into Pandas DataFrames and renders them using `seaborn` styles and `matplotlib.pyplot` drawing functions.

## Gotchas

- **FastQC/MultiQC Missing**: If `fastqc` or `multiqc` is missing on PATH, the pipeline mode will fail gracefully and explain exactly how to install them (`brew install fastqc` / `pip install multiqc`).
- **Interactive Plots**: Custom generated charts are saved as static PNGs. Interactive reports are found in `multiqc_report.html`.
- **Large FASTQ Files**: For massive datasets, ensure to specify a reasonable thread count via `--threads` to prevent high CPU utilization.

## Safety

- **Local-first**: All FastQC and MultiQC processing is performed strictly locally. No genetic data is ever uploaded.
- **No code execution**: All analysis is performed via explicit `subprocess.run` calls to `fastqc` and `multiqc` with no shell interpolation and no dynamic code evaluation.
- **Disclaimer**: Every generated `report.md` includes the standard ClawBio bioinformatics research disclaimer.

## Agent Boundary

The agent dispatches parameters and visualizes outcomes. The skill executes the native binaries and processes logs.

## Integration with Bio Orchestrator

**Trigger conditions**: routes here when:
- User mentions "bioqc", "mcp qc", "run fastqc", "fastq quality control".
- Raw FASTQ files are provided as input for pipeline execution.

**Chaining partners**:
- `multiqc-reporter`: Can consume raw data generated by the FastQC step.
- `seq-wrangler`: Can feed upstream raw reads into BioQC.

## Citations

- Andrews S. FastQC: A Quality Control Tool for High Throughput Sequence Data (2010). http://www.bioinformatics.babraham.ac.uk/projects/fastqc
- Ewels P, et al. MultiQC: Summarize analysis results for multiple tools and samples in a single report. *Bioinformatics* (2016).
