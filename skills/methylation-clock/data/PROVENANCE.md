# Demo Fixture Provenance

File: `GSE139307_small.csv.gz`

## What It Contains

- A small, local demo subset of methylation data for GEO accession `GSE139307`
- Intended for smoke testing and documentation demos of `skills/methylation-clock/methylation_clock.py`
- Stored as gzipped CSV for smaller repository footprint and lower deserialization risk

## Source And Generation

- Upstream source dataset: `GSE139307` (downloaded through PyAging helper utilities)
- Generation workflow: `PyAging_tutorial.ipynb`
- Method summary: load the full `GSE139307.pkl`, then persist a reduced two-sample subset as `GSE139307_small.csv.gz` for lightweight tests

## Integrity

- SHA-256: `48704266f1c5c39269e30d343a9bc8553c826195825f5d48afee61013f52191c`
- Size (bytes): `5845000`

## Safety Note

This file is a repository-controlled fixture and should only be loaded from this trusted repository path.

ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
