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
      - kind: pip
        package: pandas
      - kind: pip
        package: numpy
      - kind: pip
        package: matplotlib
      - kind: pip
        package: scikit-learn
      - kind: pip
        package: scipy
      - kind: pip
        package: seaborn
    trigger_keywords:
      - Differential expression analysis of proteomics data from MaxQuant or DIA-NN output
---

# 🥚 Proteomics Differential Expression

This skill performs differential expression analysis on label-free quantitative (LFQ) intensity data from MaxQuant and DIA-NN outputs, including preprocessing, imputation, statistical testing, and visualization.

---

## Domain Decisions

### 1. Multi-format Input Support
- Supports **MaxQuant `proteinGroups.txt`**
  - Automatic filtering of reverse hits, contaminants, and site-only identifications
- Supports **DIA-NN output**
  - Automatically extracts protein IDs and `.raw` intensity columns

---

### 2. Preprocessing Strategy
- MaxQuant:
  - Filters:
    - `Reverse`
    - `Potential contaminant` / `Contaminant`
    - `Only identified by site`
- DIA-NN:
  - Extracts protein identifiers and intensity matrix directly

---

### 3. Intensity Transformation
- LFQ intensities are transformed using **log2 scaling**
- Ensures approximate normality for downstream statistical testing

---

### 4. Missing Value Imputation
- Uses **down-shifted Gaussian imputation**
  - Mean shifted by: `median - shift × std`
  - Default:
    - `shift = 1.8`
    - `scale = 0.3`
- Assumption:
  - Missing values represent **low-abundance proteins**

---

### 5. Statistical Testing
- Two-sample **t-test** between treatment and control groups
- Default degrees of freedom:
  - `df = 4` (for 3 vs 3 replicates)

---

### 6. s0-based FDR Correction
- Uses **s0-based thresholding** to stabilize variance
- Combines:
  - log2 fold change
  - p-value
- Based on:
  - Giai Gianetto et al. (2016)

---

### 7. Significance Thresholding
- Default:
  - `FDR = 0.05`
  - `s0 = 0.1`
- Produces:
  - Adjusted significance boundary (used in volcano plot)

---

### 8. Visualization Outputs
- PCA plot
- Volcano plot (with s0 curve)
- Imputation distribution comparison

---

## Safety Rules

- **Local-first**
  - No data upload without explicit user consent

- **Statistical caution**
  - Statistical results should be interpreted with caution and not overinterpreted
  - Avoid drawing conclusions beyond what the data supports

- **Missing data assumptions**
  - Imputation assumes missing values correspond to low abundance
  - May not hold in all experimental designs

- **Small sample limitations**
  - t-test reliability depends on sufficient replicates

- **Reproducibility**
  - All parameters and commands are logged

- **No hallucinated science**
  - All methods are based on established proteomics workflows

---

## Agent Boundary

### This skill DOES:
- Perform differential expression analysis on LFQ proteomics data
- Handle MaxQuant and DIA-NN outputs
- Generate statistical results and visualizations
- Produce reproducible reports

---

### This skill DOES NOT:
- Process raw mass spectrometry data (e.g. RAW files)
- Perform peptide identification or database search
- Conduct pathway or functional enrichment analysis
- Provide biological interpretation of results

---

## Input Contract

### Supported Input Formats
1. MaxQuant `proteinGroups.txt`
2. DIA-NN output (`.tsv` / `.txt`)

---

### Metadata Requirements
- `.csv` or `.tsv`
- Must include:
  - `sample_id`
  - `group`

Supports:
- raw names
- full paths (e.g. `/path/sample.raw`)

---

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

---

## Usage

### Demo
```bash
python proteomics_de.py \
  --demo \
  --output report_dir
```

### MaxQuant Input
```bash
python proteomics_de.py \
  --input proteinGroups.txt \
  --input-type maxquant \
  --metadata metadata.csv \
  --contrast "treated,control" \
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

### Parameters

| Parameter            | Description           | Default         |
| -------------------- | --------------------- | --------------- |
| `--input`            | Input file path       | -               |
| `--input-type`       | `maxquant` or `diann` | maxquant        |
| `--metadata`         | Metadata file         | -               |
| `--contrast`         | treatment,control     | treated,control |
| `--s0`               | s0 parameter          | 0.1             |
| `--fdr`              | FDR threshold         | 0.05            |
| `--ttest-df`         | Degrees of freedom    | 4               |
| `--imputation-shift` | Imputation shift      | 1.8             |
| `--imputation-scale` | Imputation scale      | 0.3             |
| `--output`           | Output directory      | -               |

## References
- test_proteinGroups.txt is from: Keilhauer EC, Hein MY, Mann M. Accurate protein complex retrieval by affinity enrichment mass spectrometry (AE-MS) rather than affinity purification mass spectrometry (AP-MS). Mol Cell Proteomics. 2015 Jan;14(1):120-35. doi: 10.1074/mcp.M114.041012. Epub 2014 Nov 2. PMID: 25363814; PMCID: PMC4288248. 
- s0 correction algorithm is from: Giai Gianetto Q, Couté Y, Bruley C, Burger T. Uses and misuses of the fudge factor in quantitative discovery proteomics. Proteomics. 2016 Jul;16(14):1955-60. doi: 10.1002/pmic.201600132. PMID: 27272648.
- s0 correction algorithm is cited by: Michaelis AC, Brunner AD, Zwiebel M, Meier F, Strauss MT, Bludau I, Mann M. The social and structural architecture of the yeast protein interactome. Nature. 2023 Dec;624(7990):192-200. doi: 10.1038/s41586-023-06739-5. Epub 2023 Nov 15. PMID: 37968396; PMCID: PMC10700138.

