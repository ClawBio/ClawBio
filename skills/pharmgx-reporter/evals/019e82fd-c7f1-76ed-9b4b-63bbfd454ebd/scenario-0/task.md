# Pharmacogenomic Report from 23andMe Data

## Problem / Feature Description

A clinical informatics team at a regional hospital network is building a patient-facing portal that helps individuals understand how their genes may affect medication safety. The portal ingests raw 23andMe exports uploaded by consenting patients and produces a structured report a pharmacist can review before a prescription is written.

You have been given a sample 23andMe raw data file (`patient_data.txt`) representing a test patient. Your job is to run the PharmGx Reporter tool against this file and produce a complete pharmacogenomic report in the expected output structure. The hospital's portal expects the files to land in a specific directory layout, and the reproducibility log is required by the team's audit process.

The tool is available at `skills/pharmgx-reporter/pharmgx_reporter.py` and can be run with Python 3 from the project root. No external packages should be needed.

## Output Specification

Run the report and place results under a directory named `output/`. The grader will examine every file produced there. Ensure that:

- The output directory contains a human-readable report and a machine-readable data file.
- A reproducibility record is present showing how to re-run the analysis.
- The drug recommendation output is clearly categorised.
- The report includes a medical disclaimer.

Do not modify the input file. Clean up any temporary files larger than 50 MB.
