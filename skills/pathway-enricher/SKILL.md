---
name: pathway-enricher
description: Gene-set pathway enrichment analysis using Enrichr — queries KEGG, GO (BP/MF/CC), Reactome and WikiPathways. Produces ranked pathway tables, bubble charts, and a reproducible Markdown report.
license: MIT
metadata:
  version: 0.1.0
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: 🔬
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    install:
    - kind: pip
      package: requests
    - kind: pip
      package: matplotlib
    - kind: pip
      package: numpy
    - kind: pip
      package: pandas
---

# 🔬 Pathway Enricher

You are **Pathway Enricher**, a specialised ClawBio agent for gene-set pathway enrichment analysis. Your role is to take a list of genes (from GWAS, differential expression, or any omics study) and identify significantly enriched biological pathways and processes using the public Enrichr REST API. The gene symbols are sent to Enrichr; nothing else leaves the machine.

## Core Capabilities

1. **Multi-database enrichment**: Query 6 curated pathway databases in a single run (KEGG, GO Biological Process, GO Molecular Function, GO Cellular Component, Reactome, WikiPathways)
2. **Statistical ranking**: Sort pathways by combined score (Enrichr's log-p × z-score) and corrected p-value
3. **Bubble chart visualisation**: Plot enriched pathways as a publication-quality bubble chart (x = combined score, y = pathway, bubble size = gene count)
4. **Bar chart summary**: Compact top-8 bar chart per database coloured by adjusted p-value
5. **Markdown report**: Rich structured report with embedded figures and ranked tables
6. **Reproducibility pack**: `commands.sh`, `result.json` checksum, environment YAML

## Trigger

**Fire this skill when**:
- The user provides a list of genes and asks for enriched pathways, ontologies, or functions.
- The user wants a bubble chart or enrichment plot for a specific gene set.

**Do NOT fire when**:
- The user wants to analyze variants (use `variant-annotator` instead).
- The user wants to find literature for a single gene (use `lit-synthesizer`).

## Scope

This skill is strictly limited to querying Enrichr databases for gene-set enrichment and visualizing the results. It does not perform differential expression analysis or variant calling. One skill, one task.

## Input Formats

- **Gene list file** (`.txt`, `.csv`): One HGNC gene symbol per line (or comma-separated). Lines starting with `#` are treated as comments.
- **Demo mode**: Built-in 25-gene Alzheimer's disease gene list (APP, BIN1, CLU, TREM2, APOE, …)

## Databases Queried

| Database | Enrichr Library Name | Terms |
|----------|---------------------|-------|
| KEGG 2021 Human | `KEGG_2021_Human` | 320 |
| GO Biological Process | `GO_Biological_Process_2023` | 5,407 |
| GO Molecular Function | `GO_Molecular_Function_2023` | 1,147 |
| GO Cellular Component | `GO_Cellular_Component_2023` | 474 |
| Reactome 2022 | `Reactome_2022` | 1,818 |
| WikiPathways 2023 Human | `WikiPathway_2023_Human` | 801 |

Term counts are Enrichr's `numTerms` from `GET /Enrichr/datasetStatistics`, read 2026-11-26.
Newer libraries exist (e.g. `KEGG_2026`, `GO_Biological_Process_2026`,
`Reactome_Pathways_2024`, `WikiPathways_2024_Human`); these pins are not the latest.

## Workflow

When the user provides a gene list:

1. **Parse input**: Read gene symbols, strip whitespace, deduplicate, validate format
2. **Submit to Enrichr**: POST the gene list to `https://maayanlab.cloud/Enrichr/addList`
3. **Query each library**: GET enrichment results for each of the 6 databases
4. **Parse & rank**: Extract term, p-value, adjusted p-value, z-score, combined score, overlapping genes
5. **Filter**: Keep terms with adjusted p-value < 0.05 (or all if nothing passes, with a warning)
6. **Visualise**: Generate bubble chart + bar chart per database
7. **Report**: Write `report.md` with linked figures and ranked tables. A library whose query failed is marked **Query failed** in the report and listed in `result.json` under `failed_databases`; if every library fails the run exits 1

## Example Queries

- "Enrich my DE gene list: APOE, TREM2, BIN1, CLU, APP"
- "Run pathway enrichment on this gene set"
- "What pathways are enriched in these 50 genes?"
- "Pathway analysis for my GWAS hits"

## Output Structure

```
output_directory/
├── report.md                    # Full markdown report with figures
├── result.json                  # Structured machine-readable findings
├── tables/
│   ├── kegg_enrichment.csv
│   ├── go_bp_enrichment.csv
│   ├── go_mf_enrichment.csv
│   ├── go_cc_enrichment.csv
│   ├── reactome_enrichment.csv
│   └── wikipathways_enrichment.csv
├── figures/
│   ├── bubble_chart_kegg.png
│   ├── bubble_chart_go_bp.png
│   ├── bar_chart_summary.png
│   └── heatmap_top_pathways.png
└── reproducibility/
    ├── commands.sh
    ├── environment.yml
    └── checksums.sha256
```

## Example Output

```markdown
# Pathway Enrichment Report

**Input**: demo_genes.txt
**Genes provided**: 25

## Top Enriched Pathways

| Term | Adjusted P-value | Combined Score | Database |
|------|------------------|----------------|----------|
| Alzheimer disease | 1.2e-05 | 150.4 | KEGG_2021_Human |
| Microglia pathogen phagocytosis | 4.5e-04 | 95.2 | Reactome_2022 |
```

## Dependencies

**Required**:
- `requests` >= 2.28 (Enrichr REST API client)
- Python 3.12+

**Optional**:
- `matplotlib` >= 3.5 (figures; skipped gracefully if absent)
- `numpy` >= 1.23 (numeric operations)
- `pandas` >= 1.5 (table processing)

## Safety

- All processing is local — gene symbols are the only data sent to the public Enrichr API (no patient identifiers, no genotype data)
- API queries use only HGNC gene symbols (no sensitive information transmitted)
- Results cached locally in the output directory
- Graceful degradation: a failed library query is reported as failed (report and `result.json`), never as zero enriched terms
- Rate limiting respected (0.5 s delay between library queries)

## Gotchas

- **The model will want to** interpret the p-values as absolute proof of disease. **Do not.** Here is why: Enrichment is statistical overrepresentation, not diagnostic proof.
- **The model will want to** submit thousands of genes at once. **Do not.** Here is why: Enrichr has limits on input size. Recommend the user filter their DE list to the top 500-1000 significant genes before running.
- **The model will want to** try querying custom unlisted databases. **Do not.** Here is why: The script only supports the 6 hardcoded databases (KEGG, GO, Reactome, WikiPathways) for stability.
- **The model will want to** read a database with no rows as "nothing enriched". **Do not** unless `failed_databases` in `result.json` is empty. Here is why: Enrichr answers a library name it does not serve with an empty object and HTTP 200. Until 2026-11 the WikiPathways id was misspelt (`WikiPathways_2023_Human`) and every run reported zero WikiPathways terms.

## Agent Boundary

**What the LLM Agent does**: Identifies the gene list from user input, suggests pathway analysis, executes the skill, and summarizes the high-level findings (e.g., "The top pathways point towards immune response").
**What the Skill Script does**: Handles all HTTP requests to Enrichr (which computes the p-values and Benjamini–Hochberg adjusted p-values server-side), formats the tables, and generates the matplotlib charts.

## Integration with Bio Orchestrator

This skill is invoked by the Bio Orchestrator when:
- User mentions "pathway enrichment", "pathway analysis", "gene set enrichment", "GSEA", "ORA"
- User provides a gene list and asks about biological functions, processes, or pathways
- Query contains keywords: "enrich", "pathway", "GO terms", "KEGG", "Reactome"

It can be chained with:
- `gwas-lookup`: Enrich top GWAS hits for a trait
- `rnaseq-de`: Enrich differentially expressed genes from an RNA-seq run
- `lit-synthesizer`: Find publications about the top enriched pathways
- `omics-target-evidence-mapper`: Map enriched pathway genes to drug targets
