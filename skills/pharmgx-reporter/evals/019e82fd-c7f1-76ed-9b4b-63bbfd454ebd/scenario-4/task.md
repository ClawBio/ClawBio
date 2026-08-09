# Warfarin Safety Analysis via Programmatic API

## Problem / Feature Description

The bio-orchestrator team is integrating the PharmGx Reporter into a larger genomic profiling pipeline. Rather than shelling out to the command-line tool each time, they want to call the reporter programmatically from Python code — using an importable interface the skill exposes — so that the results can be inspected in-memory and passed to downstream pipeline stages without spawning subprocesses.

A patient file is provided at `inputs/patient_warfarin.txt`. Your task is to write a Python script called `pgx_analysis.py` that:

1. Uses the PharmGx skill's Python API (not the CLI) to analyse the genotypes from the patient file.
2. Prints the gene profiles and drug recommendations to stdout.
3. Saves a JSON file called `analysis_result.json` with the complete analysis output.

After writing and running the script, also generate a traditional CLI report into `cli_report/` for comparison purposes.

The patient has a genetic profile relevant to anticoagulant therapy. Ensure the analysis captures the drug recommendation for warfarin.

## Output Specification

- `pgx_analysis.py` — the Python script using the programmatic API
- `analysis_result.json` — the JSON output saved by the script
- `cli_report/` — full report directory generated via the CLI tool for comparison

Clean up any temporary files larger than 50 MB.
