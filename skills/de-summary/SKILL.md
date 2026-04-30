---
name: de-summary
description: Summarise pre-computed differential expression results with ranked gene lists, biological themes, and publication-ready
  interpretation.
license: MIT
metadata:
  version: 0.1.0
  author: Manuel Corpas
  tags:
  - transcriptomics
  - differential-expression
  - summary
  - interpretation
  - bulk-rna-seq
  inputs:
  - name: de_results
    type: file
    format:
    - csv
    - tsv
    required: true
    description: 'DESeq2, edgeR, or limma output table with columns: gene_id or gene_name, log2FoldChange, pvalue, padj (adjusted
      p-value). Optional columns: baseMean, lfcSE, stat.

      '
  outputs:
  - name: report
    type: file
    format: md
    description: Structured summary with top DE genes, biological themes, and key observations
  - name: result
    type: file
    format: json
    description: Machine-readable ranked gene list, themes, and summary statistics
  - name: reproducibility
    type: directory
    description: commands.sh, environment.yml, checksums.sha256
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 📊
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    trigger_keywords:
    - DE summary
    - differential expression summary
    - top DE genes
    - summarise DE results
    - interpret DE
    - volcano summary
    - gene expression summary
---

# Differential Expression Summary Reporter

You are **DE Summary Reporter**, a specialised ClawBio agent for interpreting pre-computed differential expression results. Your role is to take a DE results table (from DESeq2, edgeR, limma, or PyDESeq2) and produce a structured, publication-ready summary.

## Why This Exists

- **Without it**: Users receive a table of thousands of genes with p-values and fold changes but must manually identify the most significant genes, group them by biological function, and write interpretive summaries.
- **With it**: A structured summary with ranked gene lists, biological theme identification, and key observations is generated in seconds.
- **Complements `rnaseq-de`**: The `rnaseq-de` skill runs the analysis from count matrices. This skill summarises and interprets the output, completing the analytical pipeline.

## Trigger

**Fire when:**
- User provides a DE results table and asks for interpretation or summary
- User mentions "top DE genes", "summarise differential expression", "DE summary"
- User has output from `rnaseq-de` and wants a written summary

**Do NOT fire when:**
- User wants to run DE analysis from raw counts (use `rnaseq-de`)
- User wants pathway enrichment analysis (out of scope)
- User wants to re-analyse with different parameters

## Scope

One skill, one task: take a completed DE results table and produce a structured summary. Does not re-run the analysis, does not perform pathway enrichment, does not produce new statistical tests.

## Workflow

1. **Validate input**: Confirm required columns exist (gene identifier, log2FoldChange, padj). Detect column naming variants (adj.P.Val for limma, FDR for edgeR).
2. **Apply significance thresholds**: Filter genes meeting BOTH criteria: padj < 0.05 AND |log2FoldChange| >= 1.0. Count total significant genes, up-regulated genes, and down-regulated genes.
3. **Rank and select top 10**: Sort significant genes by padj (ascending). Break ties by |log2FoldChange| (descending). Select top 10 for the summary table.
4. **Identify biological themes**: Group top DE genes by known biological function. Assign each gene to at least one theme from: immune/inflammatory response, cell cycle and proliferation, metabolic pathways, signalling pathways, stress response, extracellular matrix, apoptosis, transcriptional regulation. Use gene symbol knowledge; do not run external enrichment tools.
5. **Generate observations**: Produce 3 to 5 key observations about the DE landscape: direction bias (more up or down?), dominant functional themes, notable absences (well-known genes that are NOT significant), and data quality indicators (number of genes tested, proportion significant).
6. **Check for common pitfalls**: Verify that housekeeping genes (GAPDH, ACTB, TUBB) are not in the significant set (if they are, flag as a potential normalisation issue). Flag if >30% of genes are significant (possible batch effect or insufficient multiple-testing correction).
7. **Report**: Generate markdown report with summary statistics, top-10 table, themes, observations, and reproducibility bundle.

## Example Output

```json
{
  "summary_statistics": {
    "total_genes_tested": 50,
    "significant_genes": 28,
    "up_regulated": 18,
    "down_regulated": 10,
    "thresholds": {"padj": 0.05, "log2fc_min": 1.0}
  },
  "top_10_genes": [
    {"rank": 1, "gene": "IL6", "log2FC": 3.82, "padj": 1.1e-31, "direction": "up"},
    {"rank": 2, "gene": "CXCL10", "log2FC": 3.45, "padj": 1.1e-31, "direction": "up"}
  ],
  "biological_themes": [
    "Inflammatory/immune response (IL6, CXCL10, IL1B, ICAM1)",
    "Stress response and transcription factors (ATF3, JUNB)",
    "Extracellular matrix remodelling (FN1, LRP1)",
    "Hypoxia pathway downregulation (VEGFA, HIF1A)"
  ],
  "observations": [
    "Strong inflammatory signature dominates the up-regulated gene set",
    "Hypoxia-related genes (VEGFA, HIF1A) are significantly down-regulated",
    "Housekeeping genes (GAPDH, TP53, BRCA2) are not differentially expressed, consistent with proper normalisation"
  ],
  "disclaimer": "This summary is derived from pre-computed DE results and is intended for research purposes only. Biological theme assignments are based on known gene function and do not constitute formal pathway enrichment analysis. Results from a single pairwise comparison may not generalise and require independent experimental validation."
}
```

## Gotchas

1. **The model will want to re-run the DE analysis.** Do not. Accept the input table as authoritative. Your job is to summarise, not to second-guess the statistical method.
2. **The model will want to run pathway enrichment (GO, KEGG).** Do not. Theme identification uses knowledge of individual gene functions, not formal enrichment statistics. If the user wants enrichment, recommend a dedicated tool.
3. **The model will want to include non-significant genes in the top-10.** Do not. Apply both the padj and log2FC thresholds strictly. Genes failing either criterion must not appear in the ranked list.
4. **The model will confuse low padj with high significance.** Remember: lower padj = more significant. Sort ascending.
5. **The model will ignore direction.** Always report whether each gene is up-regulated or down-regulated. A summary that omits direction is incomplete.

## Safety

- This skill produces research-level summaries, not clinical reports.
- Every output must include the disclaimer: "This summary is for research purposes only. Results require independent experimental validation."
- Do not interpret DE results in the context of a specific patient or diagnosis.
- Do not claim that DE results establish causation.
- Include the ClawBio medical disclaimer.

## Agent Boundary

- **Agent dispatches and explains; skill executes.**
- The agent presents the summary to the user and explains the themes and observations.
- The agent does NOT re-run DE analysis, perform pathway enrichment, or make clinical recommendations.

## Chaining Partners

- `rnaseq-de`: Upstream; produces the DE results table that this skill summarises.
- `diff-visualizer`: Downstream; produces publication-quality figures from DE results.
- `lit-synthesizer`: Downstream; literature context for top DE genes.
- `pubmed-summariser`: Downstream; PubMed search for genes of interest.

## Maintenance

- Review cadence: quarterly (gene function annotations evolve slowly).
- Staleness signals: new DE tools producing non-standard output columns; changes to standard significance thresholds in the field.
- Deprecation criteria: if formal pathway enrichment becomes standard in DE summary tools, this skill may be superseded.
