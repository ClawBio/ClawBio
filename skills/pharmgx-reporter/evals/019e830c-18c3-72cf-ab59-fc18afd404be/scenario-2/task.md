# Warfarin Safety Assessment for a New Patient

## Problem/Feature Description

A cardiology clinic is starting a new patient on warfarin (Coumadin) for atrial fibrillation management. Warfarin has a notoriously narrow therapeutic window, and incorrect dosing can lead to dangerous bleeding or clotting events. The clinic has access to the patient's consumer genetic data and wants a pharmacogenomic assessment before initiating therapy.

The clinical coordinator has asked you to write a simple analysis script that processes the patient's genetic file and produces a pharmacogenomic report. The report should cover the patient's gene profiles and provide drug classifications the clinic can act on.

The patient's genetic data file is at `inputs/patient_data.txt`. Write the output to a directory called `warfarin_report/`.

## Output Specification

Write a script named `warfarin_analysis.py` that processes the genetic file and writes results to `warfarin_report/`. Run the script:

```
python warfarin_analysis.py --input inputs/patient_data.txt --output warfarin_report
```

The workspace should contain `warfarin_analysis.py` and everything the script writes to `warfarin_report/`.
