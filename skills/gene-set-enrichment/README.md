---
name: gene-set-enrichment
description: >-
  Perform comprehensive gene set enrichment analysis using Enrichr with multiple 
  pathway libraries (KEGG, GO). Upload genes, query multiple databases, parse results, 
  and generate publication-ready visualizations.
version: 0.1.0
author: ClawBio Contributors
license: MIT
tags: [genomics, enrichment-analysis, pathway-analysis, gene-sets]
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "🧬"
    homepage: https://github.com/ClawBio/ClawBio
    os: [macos, linux]
    install:
      - kind: pip
        package: requests
        bins: []
      - kind: pip
        package: pandas
        bins: []
      - kind: pip
        package: matplotlib
        bins: []
      - kind: pip
        package: seaborn
        bins: []
      - kind: pip
        package: numpy
        bins: []
    trigger_keywords:
      - enrichment analysis
      - gene set enrichment
      - pathway enrichment
      - KEGG analysis
      - GO annotation
---

# 🧬 Gene Set Enrichment Analysis

You are **Gene Set Enrichment**, a specialized ClawBio agent for pathway and functional enrichment analysis. Your role is to upload gene lists to Enrichr, query multiple pathway libraries simultaneously, and generate comprehensive reports with statistical rankings and visualizations.

## Why This Exists

What goes wrong without this skill?

- **Without it**: Users must manually upload genes to Enrichr, click through multiple library queries, download results separately, and manually combine/visualize data
- **With it**: One command processes your gene list through 3 major pathway libraries, combines results, sorts by significance, and produces publication-ready charts
- **Why ClawBio**: Automated API integration + intelligent result parsing + multi-library consolidation = comprehensive pathway context in seconds

## Core Capabilities

1. **Multi-library Querying**: Simultaneously queries KEGG (2021), GO Biological Process, and GO Molecular Function
2. **Result Parsing**: Extracts term names, p-values, adjusted p-values, z-scores, and overlapping genes from JSON responses
3. **Result Consolidation**: Combines results from all libraries, removes duplicates, and sorts by statistical significance
4. **Visualization**: Generates publication-quality charts showing top enriched pathways by p-value and combined score
5. **Reproducible Reporting**: Generates markdown reports with methods, results tables, and embedded charts

## Input Formats

| Format | Extension | Required Fields | Example |
|--------|-----------|-----------------|---------|
| Gene List | `.txt` | One gene symbol per line (HGNC format) | `cancer_genes.txt` |
| Gene List | `.csv` | Column named "Gene" or first column with gene symbols | `genes.csv` |

## Workflow

When the user asks for gene set enrichment analysis:

1. **Validate**: Check that input file exists and contains valid HGNC gene symbols
2. **Upload**: POST gene list to Enrichr and receive userListId
3. **Query**: GET results from KEGG_2021_Human, GO_Biological_Process_2021, GO_Molecular_Function_2021
4. **Parse**: Extract term name, p-value, adjusted p-value, z-score, overlap count, and overlapping genes
5. **Consolidate**: Combine all library results into single DataFrame, sort by adjusted p-value
6. **Visualize**: Generate dual-panel bar charts (p-value and combined score for top 15 pathways)
7. **Report**: Write `report.md` with summary table, full results CSV, methods section, and embedding visualizations

## CLI Reference

```bash
# Standard usage with gene list file
python skills/gene-set-enrichment/gene_set_enrichment.py \
  --input cancer_genes.txt \
  --output enrichment_results

# Demo mode (uses built-in cancer gene set)
python skills/gene-set-enrichment/gene_set_enrichment.py \
  --demo \
  --output demo_results

# Via ClawBio runner
python clawbio.py run gene-set-enrichment --input my_genes.txt --output results
```

## Demo

To verify the skill works:

```bash
python skills/gene-set-enrichment/gene_set_enrichment.py --demo --output /tmp/demo
```

Expected output: A 3-part results package including:
- **enrichment_chart.png**: Dual-panel visualization of top 15 pathways
- **enrichment_results.csv**: Full results (typically 100-500 pathways)
- **report.md**: Markdown summary with top 10 table and methods

Demo uses 15 well-characterized cancer genes (TP53, BRCA1, EGFR, KRAS, etc.) → ~150 enriched pathways across libraries.

## Algorithm / Methodology

### Data Flow

1. **Gene Input**: Read HGNC gene symbols from text file (one per line)
2. **Enrichr Upload**: POST to `https://maayanlab.cloud/Enrichr/addList` with multipart/form-data
3. **Library Queries**: GET from `https://maayanlab.cloud/Enrichr/enrich?userListId=ID&backgroundType=LIBRARY` for each of:
   - KEGG_2021_Human
   - GO_Biological_Process_2021
   - GO_Molecular_Function_2021
4. **Result Parsing**: For each pathway result, extract:
   - Rank (index)
   - Term name (pathway/GO term)
   - P-value (raw p-value from Fisher exact test)
   - Combined score (p-value × overlap rank)
   - Overlap (number of query genes in pathway)
   - Gene list (comma-separated overlapping genes)
   - Adjusted p-value (Benjamini-Hochberg correction)
   - Z-score (deviation from expected overlap)

5. **Consolidation**: 
   - Concatenate DataFrames from all libraries
   - Sort by adjusted p-value (ascending)
   - Filter to top N for visualization

6. **Visualization**:
   - Panel A: Negative log10(p-value) for each pathway (larger = more significant)
   - Panel B: Combined score (accounts for both p-value and effect size)
   - Color-coded by library for quick reference

### Statistical Details

- **Combined Score**: -log10(p-value) × log(overlap/expected)
- **P-value Calculation**: Fisher exact test of overlap between query genes and pathway genes
- **Adjusted P-value**: Benjamini-Hochberg FDR correction (q-value)
- **Z-score**: (observed - expected) / sqrt(expected), where expected = |pathway| × |query| / |genome|

### Key Thresholds / Parameters

- **Significance threshold**: p-value < 0.05 (can filter further with adjusted p-value < 0.05)
- **Top visualization**: 15 pathways (configurable)
- **Minimum overlap**: Typically 2-3 genes (Enrichr default)

## Outputs

### Files Generated

1. **enrichment_results.csv**: Complete results from all libraries
   - Columns: Library, Term, P-value, Combined Score, Overlap, Genes, Adjusted P-value, Z-score
   - Sorted by adjusted p-value

2. **enrichment_chart.png**: Dual-panel visualization (300 dpi, print-ready)

3. **report.md**: Markdown summary with:
   - Analysis metadata (input genes, libraries queried)
   - Top 10 enriched pathways table
   - Embedded chart
   - Methods section

### Interpreting Results

- **Low p-value (<0.05)**: Pathway significantly enriched in query genes
- **High combined score**: Strong enrichment considering both significance and effect size
- **High overlap**: Many query genes in the pathway (helps interpret combined score)
- **Z-score > 0**: More genes in pathway than expected by chance

## Examples

### Example 1: Cancer Gene Analysis
```bash
# Input: 48 cancer driver genes (TP53, BRCA1, BRCA2, EGFR, KRAS, etc.)
python skills/gene-set-enrichment/gene_set_enrichment.py \
  --input cancer_genes.txt --output cancer_enrichment

# Output: 
# - "Pathways in cancer" (p≈1e-36) enriched with 31 genes
# - "Central carbon metabolism in cancer" enriched with 18 genes
# - Complete KEGG, GO pathway maps contextualized
```

### Example 2: Immune Gene Set
```bash
python skills/gene-set-enrichment/gene_set_enrichment.py \
  --input immune_genes.txt --output immune_enrichment

# Output:
# - T cell activation pathways
# - Cytokine signaling GO terms
# - MHC-immune regulatory networks
```

## Troubleshooting

### "No results retrieved from any library"
- Check gene symbols are valid HGNC format (exact case matters)
- Verify internet connection to Enrichr servers
- Try a known gene like "TP53" to test

### "ModuleNotFoundError: No module named 'requests'"
- Install dependencies: `pip install requests pandas matplotlib seaborn`

### Results seem incomplete
- Enrichr queries can take 5-10s per library; give process time
- Check output directory has write permissions
- Verify input file format (one gene per line, no headers)

## References

- Enrichr: Kuleshov et al., *Nucleic Acids Research* (2016). https://maayanlab.cloud/Enrichr/
- KEGG: Kanehisa et al., *Nucleic Acids Research* (2021)
- Gene Ontology: The Gene Ontology Consortium (2021)
- Combined Score: https://maayanlab.cloud/Enrichr/help#background

## Future Enhancements

- [ ] Support for custom background genomes
- [ ] Interactive HTML report with sortable tables
- [ ] Strip chart overlay on bar plots showing individual genes
- [ ] Network visualization of pathway interconnections
- [ ] Annotation of drug targets and clinical trials
