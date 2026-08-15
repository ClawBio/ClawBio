# Pharmacogenomic Report from Raw Genetic Data

## Problem/Feature Description

A precision medicine clinic has recently begun offering pharmacogenomic consultation services. A patient has uploaded their raw genetic data file (exported from a consumer DNA testing service) and wants to understand which of their regular medications might require dose adjustments or alternatives based on their genetic profile.

The clinic's bioinformatics team has asked you to build a command-line tool that can process this raw genetic data file and produce a structured pharmacogenomic report. The tool should work with the data format exported directly from consumer DNA services, and produce a self-contained analysis that a pharmacist can review.

The patient's data file is provided at `inputs/patient_data.txt`. Write outputs to a directory named `report_output/`. The tool must run using Python 3 with no additional software to install.

## Output Specification

Write a Python script named `pharmgx_tool.py` that:
- Accepts `--input` and `--output` command-line arguments
- Processes `inputs/patient_data.txt` and writes results to `report_output/`
- Is runnable as: `python pharmgx_tool.py --input inputs/patient_data.txt --output report_output`

Run the script and leave all output files in the workspace. The workspace should contain `pharmgx_tool.py` and everything written to `report_output/` by running the script.
