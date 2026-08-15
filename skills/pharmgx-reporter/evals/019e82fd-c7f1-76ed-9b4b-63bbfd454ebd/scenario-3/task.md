# Pharmacogenomic Report with Incomplete Genetic Data

## Problem / Feature Description

A medical informatics researcher is stress-testing the PharmGx Reporter to see how it handles patients whose raw genetic files are incomplete — either because certain SNP positions were not assayed by their chip, or because the patient carries multiple heterozygous variants in certain genes. Knowing how the tool communicates uncertainty in its output is critical before it can be used in a clinical-adjacent setting.

You have a synthetic 23andMe-format file at `inputs/partial_patient.txt` that was deliberately constructed with partial SNP coverage. Run the PharmGx Reporter against this file and save the report to `report_output/`.

After generating the report, write a short summary called `data_quality_notes.txt` covering:
1. Which genes received an uncertain or unknown result, and what the tool reported as the reason
2. Any special warnings about specific gene assay limitations that appeared in the report

## Output Specification

- `report_output/` — full report directory from the PharmGx Reporter
- `data_quality_notes.txt` — your written summary of uncertain gene results and any assay limitation notices found in the generated report

Clean up any temporary files larger than 50 MB.
