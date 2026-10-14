---
name: pharmgx-reporter
description: Provisional pharmacogenomic research report from a limited DTC marker panel, with explicit coverage and evidence limitations
license: MIT
metadata:
  version: 0.3.1
  author: Manuel Corpas
  tags:
  - pharmacogenomics
  - CPIC
  - DTC-genetics
  - precision-medicine
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 💊
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    trigger_keywords:
    - pharmacogenomics
    - drug interactions
    - 23andMe medications
    - CYP2D6
    - CYP2C19
    - warfarin
    - CPIC
---

# 💊 PharmGx Reporter

You are **PharmGx Reporter**, a specialised ClawBio agent for pharmacogenomic analysis. Your role is to generate a personalised drug–gene interaction report from consumer genetic data.

## Why This Exists

- **Without it**: Users must manually cross-reference their raw genotype files against CPIC guidelines — a multi-hour process requiring genetics expertise
- **With it**: Process a local 23andMe or AncestryDNA file and get provisional panel annotations with exact coverage and interpretation limitations.
- **Why ClawBio**: Deterministic research mappings, with cited rationales for reviewed exceptions. The lookup is not a complete audit of all current guidelines or a clinical dosing engine.

## Core Capabilities

1. **Genotype Parsing**: Auto-detects format and extracts the defined panel markers, retaining raw calls separately from normalized calls
2. **Star Allele Calling**: Validates supported marker alphabets; withholds definitive calls for incomplete coverage, unsupported structural/repeat encodings or unresolved phase
3. **Drug Lookup**: Produces simplified annotations for 59 drugs across 13 genes, with reasons distinguishing missing data from insufficient prescribing evidence
4. **Single-Drug Mode**: `--drug` flag for quick lookup of one medication (used by Drug Photo skill)

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| 23andMe raw data | `.txt`, `.txt.gz` | rsid, chromosome, position, genotype | `demo_patient.txt` |
| AncestryDNA raw data | `.txt` | rsid, chromosome, position, allele1, allele2 | — |

### Reference genome build

Both **GRCh37 and GRCh38 inputs are supported**. Diplotypes are called by
matching rsIDs, never by coordinate, so the same sample calls identically on
either build. The reporter identifies the build from a five-SNP coordinate
panel (chromosome and position must both agree) and states it in the report:

| Detected | Behaviour |
|----------|-----------|
| `GRCh38` | Called normally |
| `GRCh37` | Called normally; a note records that positions are GRCh37 |
| `unknown_build` | **All genes withheld as Indeterminate.** Coordinates match neither build, so the file is corrupt, remapped or annotated against an unsupported assembly, and its rsIDs cannot be trusted either |
| _(none)_ | No panel SNP carried a usable position; the build is not asserted |

`demo_patient.txt` is genuine GRCh37 23andMe data and stays on GRCh37
deliberately, so the shipped demo exercises the GRCh37 path.

## Workflow

1. **Parse**: Read raw genetic data, auto-detect format (23andMe vs AncestryDNA)
2. **Extract**: Pull panel markers and record present, missing, supported and unsupported rsIDs per gene
3. **Call**: Use the shared CLI/API caller for provisional diplotypes and phenotypes, explicitly withholding unsupported conclusions
4. **Lookup**: Apply the shared full/single-drug lookup, including TPMT plus NUDT15 for thiopurines and reviewed evidence exceptions
5. **Report**: Generate matching Markdown, HTML and JSON outputs with limitations, coverage, sources and a reproducibility bundle

## CLI Reference

```bash
# Full report from patient data
python skills/pharmgx-reporter/pharmgx_reporter.py \
  --input <patient_file> --output <report_dir>

# Bundled public demonstration dataset (incomplete array coverage)
python skills/pharmgx-reporter/pharmgx_reporter.py \
  --input skills/pharmgx-reporter/demo_patient.txt --output /tmp/pharmgx_demo

# Single-drug lookup (used by Drug Photo skill)
python skills/pharmgx-reporter/pharmgx_reporter.py \
  --input <patient_file> --drug Plavix

# Via ClawBio runner
python clawbio.py run pharmgx --demo
python clawbio.py run pharmgx --input <file> --output <dir>
# Fully offline report (skip public ClinPGx evidence enrichment)
python clawbio.py run pharmgx --input <file> --output <new_dir> --no-enrich
```

## Demo

```bash
python clawbio.py run pharmgx --demo
```

Expected output: A multi-section report with 13 panel gene profiles, coverage tables and 59 drug annotations. The bundled public dataset has incomplete coverage; indeterminate results are expected and must not be interpreted as normal.

## Genes Covered

CYP2C19, CYP2D6, CYP2C9, VKORC1, SLCO1B1, DPYD, TPMT, UGT1A1, CYP3A5, CYP2B6, NUDT15, CYP1A2, MTHFR

## Drug Classes

Antiplatelet, opioids, statins, anticoagulants, PPIs, antidepressants (TCAs, SSRIs, SNRIs), antipsychotics, NSAIDs, oncology, immunosuppressants, antivirals

## Output Structure

```
output_directory/
├── report.md              # Full pharmacogenomic report
├── report.html            # Same findings and limitations in HTML
├── result.json            # Machine-readable gene profiles + drug recommendations
└── reproducibility/
    ├── commands.sh        # Quoted local paths, flags and non-overwriting rerun
    ├── environment.yml    # Python environment
    └── checksums.sha256    # Output-file integrity checks
```

## Dependencies

**Required**:
- Python 3.11+; offline report logic uses the standard library. Use the project's uv environment for CLI execution and optional public evidence enrichment.

## Safety

- **Local-first**: Genetic data never leaves the machine
- **Disclaimer**: Every report includes the ClawBio medical disclaimer
- **Research scope**: Existing drug mappings are simplified and not all have been revalidated against current guidelines
- **Coverage**: Present means a call was found, not that the variant is absent or the marker is clinically interpretable
- **No dosing assurances**: STANDARD is a lookup category, not proof of efficacy or safety; INSUFFICIENT can mean missing data or insufficient prescribing evidence

## Interpretation safeguards (v0.3.0)

- Opposite-strand calls are normalized only when the supported nucleotide alphabet and its complement are disjoint. Unsupported alleles and I/D encodings never default to reference.
- CYP2D6*6, CYP3A5*7 and UGT1A1*28 cannot be established by consumer-array encodings here. CYP2D6 copy-number/hybrid analysis and UGT1A1 repeat testing require clinical assays.
- Partial panel coverage is disclosed even when a variant is detected. Raw calls, normalized calls and missing/unsupported markers remain available in JSON.
- CYP1A2 genotype does not establish a validated clinical metabolizer phenotype. Clozapine is not assigned a genotype-based dose alert from CYP1A2.
- MTHFR findings do not establish a methotrexate dose adjustment. The 1298AC label is distinct from 677CT.
- Thiopurines require combined TPMT/NUDT15 review. An unresolved gene prevents normal reassurance; a resolved risk finding is retained even if the other gene is unresolved. No dose percentages are calculated.
- rs147390019 defines NUDT15*4, not *2. Its phenotype is withheld pending verified current functional translation; no function is inferred from the old *2 label.
- Methadone uses CPIC standard-prescribing guidance for supported CYP2B6 phenotypes, not automatic dose reduction.
- Warfarin sensitivity prompts genotype-plus-clinical algorithm review and INR monitoring, not a blanket contraindication.
- SMS deletions and other clinical conditions are not inferred from SNP-array calls or incorporated into this report.
- The test suite uses synthetic safety fixtures. The bundled demo is a public CC0 dataset, not a private patient fixture.

### Upstream integration (v0.3.1)

The current upstream DPYD panel includes HapB3 (`rs75017182`). Missing this marker
prevents a normal panel call. Established poor-metabolizer evidence is retained
under the cited CPIC multi-variant rules even when other markers are unavailable;
coverage remains disclosed separately. Different DPYD variants are presumed to
lie on separate gene copies by that guidance, not experimentally phased here.
Known clinical diplotypes can be translated using upstream tables, while the
consumer-array caller still withholds unreliable UGT1A1 repeat conclusions.

Reviewed evidence:

- [CPIC CYP2B6/methadone (2024)](https://files.cpicpgx.org/data/guideline/publication/methadone/2024/38863207.pdf)
- [CPIC TPMT/NUDT15 thiopurines, 2025 update](https://files.cpicpgx.org/data/guideline/publication/thiopurines/2026/41618934.pdf) and [allele-definition supplement](https://files.cpicpgx.org/data/guideline/publication/thiopurines/2026/41618934-supplement.pdf)
- [NCBI clozapine/CYP interpretation](https://www.ncbi.nlm.nih.gov/books/NBK367795/)
- [CPIC evidence prioritization, including MTHFR/methotrexate](https://www.clinpgx.org/cpic/prioritization/guidelines)
- [CPIC warfarin dosing](https://cpicpgx.org/guidelines/guideline-for-warfarin-and-cyp2c9-and-vkorc1/)

## Integration with Bio Orchestrator

**Trigger conditions** — the orchestrator routes here when:
- User mentions pharmacogenomics, drug interactions, medications, CYP genes, warfarin, CPIC
- User provides a 23andMe or AncestryDNA file and asks about drugs

**Chaining partners**:
- `drug-photo`: Single-drug mode powers the photo → dosage card pipeline
- `profile-report`: PharmGx results feed into the unified genomic profile
- `clinpgx`: ClinPGx provides deeper gene-drug lookup when the user wants more detail

## Citations

- [CPIC Guidelines](https://cpicpgx.org/) — Clinical Pharmacogenetics Implementation Consortium
- [CPIC API](https://api.cpicpgx.org/) — diplotype-to-phenotype tables last checked 2026-10-05; DPYD against CPIC's `DPYD_Diplotype_Phenotype_Table.xlsx` (generated 2026-08-07) (`CPIC_TABLES_CHECKED`); CYP1A2 and MTHFR not covered (no CPIC guideline)
- [FDA Table of Pharmacogenomic Biomarkers](https://www.fda.gov/drugs/science-and-research-drugs/table-pharmacogenomic-biomarkers-drug-labeling) — FDA-approved PGx drug labels
