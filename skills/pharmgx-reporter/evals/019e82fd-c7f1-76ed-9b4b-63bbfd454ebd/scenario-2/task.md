# Generating a Pharmacogenomic Report from AncestryDNA Data

## Problem / Feature Description

A genetic counselling service has been receiving patient files from a variety of direct-to-consumer genetic testing companies, not just 23andMe. A new batch of patients has submitted their raw data downloaded from AncestryDNA. The counselling team wants to verify that the PharmGx Reporter can handle this format and produce a complete, gene-level report covering all the pharmacogenomic genes in the panel.

You have been provided with an AncestryDNA raw data file at `inputs/ancestry_patient.txt`. Use the PharmGx Reporter tool (located at `skills/pharmgx-reporter/pharmgx_reporter.py`) to analyse this file and generate a full pharmacogenomic report into a directory called `pgx_output/`.

The counselling team specifically needs to know the metaboliser phenotype for each gene in the panel, and wants the drug recommendations clearly categorised. They also need a machine-readable data file for downstream processing.

## Output Specification

- A full report under `pgx_output/` including:
  - A human-readable report file
  - A machine-readable data file
  - A reproducibility record
- The report must reflect the AncestryDNA format in its metadata (the detected format should appear in the report header).

Do not paste the file contents into any report section. Clean up any temporary files larger than 50 MB.
