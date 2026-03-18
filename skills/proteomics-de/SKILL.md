---
name: proteomics-de
description: Differential expression analysis for label-free quantitative (LFQ) intensity data with standard MaxQuant and DIA-NN output. Workflow includes preprocessing, imputation, and statistical testing.
version: 0.1.0
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🥚"
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux, win32]
    install:
      - kind: uv
        package: pandas
        bins: []
      - kind: uv
        package: numpy
        bins: []
      - kind: uv
        package: matplotlib
        bins: []
      - kind: uv
        package: scikit-learn
        bins: []
      - kind: uv
        package: scipy
        bins: []
      - kind: uv
        package: seaborn
        bins: []
---

# 🥚 Proteomics Differential Expression

This skill performs differential expression analysis on label-free quantitative(LFQ) intensity data from MaxQuant and DIA-NN output.

## Core Capabilities

1. **Multi-format Input Support**:
   - MaxQuant `proteinGroups.txt` format (with automatic contaminant/reverse filtering)
   - DIA-NN output format (extracts protein intensities from `.raw` columns)
2. **Preprocessing**:
   - MaxQuant: Filter reverse hits, contaminants, and proteins only identified by site (supports both "Potential contaminant" and "Contaminant" column naming)
   - DIA-NN: Automatic extraction of protein IDs and intensity columns
3. **Matrix extraction**: Extract protein identifiers and intensity columns automatically
4. **Intensity Log2 Transformation**: Log2 transformation of quantitative intensities
5. **Imputation**: Down-shift imputation for missing values, with distribution comparison visualization
6. **Statistical testing**: Two-sample t-test between treatment and control groups
7. **Significance calibration**: s0-based FDR correction for log2FC and p-value validation (s0_based_FDR_correction)
8. **Smart sample matching**: Automatic cleanup of sample names (supports paths and `.raw` suffix in metadata)
9. **Visualization**: PCA, volcano plot with s0 threshold curve, and imputation distribution comparison
10. **Output**: Markdown report with top differentially expressed proteins tables, result tables, and full reproducibility files**

## Input Contract

### Input Formats
1. **MaxQuant proteinGroups.txt** (`.txt` tab-separated file)
2. **DIA-NN output** (`.tsv`/`.txt` file from DIA-NN quantification)

### Common Requirements
- Metadata table (`.csv` or `.tsv`): one row per sample, must include `sample_id` and `group` columns
  - `sample_id` supports both raw sample names and full paths with `.raw` suffix (e.g. `/path/to/sample1.raw`)
- Contrast: `treatment,control` (e.g. `treated,control`)

### MaxQuant proteinGroups.txt column description
| Column header                        | Description                                                                                           |
| ------------------------------------ | ----------------------------------------------------------------------------------------------------- |
| Protein IDs                          | Uniprot identifiers of all sequences in the protein group                                             |
| Majority protein IDs                 | Uniprot identifiers of all sequences containing at least half of the peptides of the leading sequence |
| Protein names                        | Name of the protein                                                                                   |
| Gene names                           | Gene name                                                                                             |
| Fasta headers                        | Header of the fasta entry                                                                             |
| Number of proteins                   | Number of proteins within the protein group                                                           |
| Peptides                             | Total number of identified peptides of the protein group                                              |
| Razor + unique peptides              | Number of unique peptides plus razor peptides in the protein group                                    |
| Unique peptides                      | Number of unique peptides in the protein group                                                        |
| Peptides…01/2/3/(4)                  | Total number of identified peptides in the sample                                                     |
| Razor+unique peptides…01/2/3/(4)     | Number of unique peptides plus razor peptides in the sample                                           |
| Peptides…01/2/3/(4)                  | Number of unique peptides in the sample                                                               |
| Sequence coverage [%]                | Coverage of the protein sequence in %                                                                 |
| Unique + razor sequence coverage [%] | Coverage of the protein sequence in % using only unique and razor peptides                            |
| Unique sequence coverage [%]         | Coverage of the protein sequence in % using only unique peptides                                      |
| Mol. weight [kDa]                    | Molecular weight of the protein                                                                       |
| Sequence length                      | Length of the amino acid sequence                                                                     |
| PEP                                  | Posterior error probability                                                                           |
| Sequence Coverage…01/2/3/(4)         | Coverage of the protein sequence in % in the sample                                                   |
| Intensity                            | Summed raw intensity                                                                                  |
| Intensity …_01/2/3/(4)               | Raw intensity in the sample                                                                           |
| LFQ intensity..._01/2/3/(4)          | Label-free quantification (LFQ) intensity (MaxLFQ algorithm) in the sample                            |
| MS/MS Count..._01/2/3/(4)            | Number of MS/MS events in the sample                                                                  |
| Only identified by site              | "+" indicates proteins only identified with modified peptides                                         |
| Reverse                              | "+" indicates hits to the reverse database                                                            |
| Contaminant                          | "+" indicates proteins belonging to the contaminants fasta file                                       |

### DIA-NN Output Column Description
| Column header        | Description                                                                 |
| -------------------- | --------------------------------------------------------------------------- |
| `Protein.Ids`        | Protein identifiers (UniProt IDs)                                           |
| `*.raw`              | Quantification LFQ intensity columns for each sample (ends with .raw suffix)|
## Output Structure

```
proteomics_de_report/
├── report.md
├── figures/
│   ├── imputation_distribution.png
│   ├── pca.png
│   └── volcano.png
├── tables/
│   ├── imputed_proteinGroups.csv
│   └── de_results.csv
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Usage

### MaxQuant Input
```bash
python proteomics_de.py \
  --input proteinGroups.txt \
  --input-type maxquant \
  --metadata metadata.csv \
  --contrast "treated,control" \
  [--s0 0.1] \
  [--fdr 0.05] \
  [--ttest-df 4] \
  [--imputation-shift 1.8] \
  [--imputation-scale 0.3] \
  --output report_dir
```

### DIA-NN Input
```bash
python proteomics_de.py \
  --input diann_output.tsv \
  --input-type diann \
  --metadata metadata.csv \
  --contrast "treated,control" \
  --output report_dir
```

## Parameters

| Parameter | Description | Default |
|---|---|---|
| `--input` | Path to input protein quantification file (MaxQuant proteinGroups.txt or DIA-NN output, required) | - |
| `--input-type` | Input file type: `maxquant` (default) or `diann` | `maxquant` |
| `--metadata` | Path to sample metadata CSV/TSV file, must contain `sample_id` and `group` columns (required) | - |
| `--contrast` | Comparison groups in format `treatment,control` | `treated,control` |
| `--s0` | s0 parameter for smooth threshold FDR correction | 0.1 |
| `--fdr` | False discovery rate threshold | 0.05 |
| `--ttest-df` | Degree of freedom for t-test (default 4 for 3+3 replicates) | 4 |
| `--imputation-shift` | Down-shift factor for missing value imputation (median - shift*std) | 1.8 |
| `--imputation-scale` | Scale factor for missing value imputation (scale*std) | 0.3 |
| `--output` | Output directory path (required) | - |

## Output Table Schema

### `tables/de_results.csv`
| Column | Description |
|---|---|
| `protein_id` | Protein identifiers (Majority protein IDs for MaxQuant, Protein.Ids for DIA-NN) |
| `mean_intensity` | Mean log2 intensity across all samples |
| `mean_treatment` | Mean log2 intensity in treatment group |
| `mean_control` | Mean log2 intensity in control group |
| `log2FoldChange` | Log2 fold change (treatment / control) |
| `pvalue` | Raw p-value from two-sample t-test |
| `-log10(pvalue)` | -log10 transformed p-value |
| `s0_corrected_-log10(pvalue)` | s0-based FDR corrected -log10(pvalue) threshold |
| `regulation` | Regulation status: `upregulated`, `downregulated`, `non significant` |

### `tables/imputed_proteinGroups.csv`
Imputed log2-transformed LFQ intensity matrix with protein IDs as index and samples as columns.

## Safety

- Local-only processing, no proteomics data leaves the machine
- Warn before overwriting existing output
- Report-level disclaimer required

## Reference
- test_proteinGroups.txt is from: Keilhauer EC, Hein MY, Mann M. Accurate protein complex retrieval by affinity enrichment mass spectrometry (AE-MS) rather than affinity purification mass spectrometry (AP-MS). Mol Cell Proteomics. 2015 Jan;14(1):120-35. doi: 10.1074/mcp.M114.041012. Epub 2014 Nov 2. PMID: 25363814; PMCID: PMC4288248.
- s0 correction algorithm is from: Giai Gianetto Q, Couté Y, Bruley C, Burger T. Uses and misuses of the fudge factor in quantitative discovery proteomics. Proteomics. 2016 Jul;16(14):1955-60. doi: 10.1002/pmic.201600132. PMID: 27272648.
- s0 correction algorithm is cited by: Michaelis AC, Brunner AD, Zwiebel M, Meier F, Strauss MT, Bludau I, Mann M. The social and structural architecture of the yeast protein interactome. Nature. 2023 Dec;624(7990):192-200. doi: 10.1038/s41586-023-06739-5. Epub 2023 Nov 15. PMID: 37968396; PMCID: PMC10700138.
