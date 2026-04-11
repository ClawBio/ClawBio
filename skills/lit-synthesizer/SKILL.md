---
name: lit-synthesizer
description: >
  Search PubMed and bioRxiv for bioinformatics literature, synthesise results
  into a structured report, and build a citation graph — all locally, with a
  reproducibility bundle.
version: 0.1.0
author: Sooraj (github.com/sooraj-codes)
domain: literature
license: MIT
inputs:
  - name: query
    type: string
    description: "Free-text search query (e.g. 'CRISPR off-target effects')"
    required: true
outputs:
  - name: report
    type: file
    format: md
    description: Structured markdown report with paper summaries and citation graph
dependencies:
  python: ">=3.11"
  packages:
    - biopython>=1.83
tags:
  - literature
  - pubmed
  - biorxiv
  - citation
  - synthesis
demo_data:
  - path: examples/demo_output/report.md
    description: Pre-generated demo report for CRISPR genome editing query
endpoints:
  cli: python skills/lit-synthesizer/lit_synthesizer.py --query "{query}" --output {output_dir}
metadata:
  openclaw:
    requires:
      always: false
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - macos
      - linux
    install:
      - kind: pip
        package: biopython
    trigger_keywords:
      - "search pubmed"
      - "find papers"
      - "literature review"
      - "search biorxiv"
      - "find articles"
      - "citation graph"
      - "synthesize literature"
      - "find research papers"
      - "pubmed search"
      - "recent papers on"
---

# 🦖 Lit Synthesizer

You are **Lit Synthesizer**, a specialised ClawBio agent for biomedical literature
discovery and synthesis. Your role is to search PubMed and bioRxiv, summarise
retrieved papers, and build a citation graph — all locally with a reproducibility bundle.

## Trigger

**Fire this skill when the user says any of:**

- "search pubmed for X"
- "find papers on X"
- "literature review on X"
- "search biorxiv for X"
- "find recent articles about X"
- "build a citation graph for X"
- "synthesize the literature on X"
- "what papers exist on X"
- "find research on X"
- "summarise the literature on X"

**Do NOT fire when:**

- The user wants to annotate a VCF file (route to `vcf-annotator`)
- The user wants pharmacogenomic drug recommendations (route to `pharmgx-reporter`)
- The user is asking a general biology question without a search intent

## Why This Exists

**Without it**: A researcher must manually search PubMed, download abstracts,
read each one, spot connections across papers, and format everything by hand.
This can take hours for a single topic.

**With it**: One command searches both PubMed and bioRxiv, summarises abstracts,
identifies recurring themes, builds a citation graph, and outputs a formatted
report with a reproducibility bundle — in under 30 seconds.

**Why ClawBio**: A general LLM will hallucinate paper titles, fabricate authors,
and invent DOIs. This skill uses live API calls to real databases, so every
paper it returns is real and verifiable.

## Core Capabilities

1. **PubMed search**: Queries NCBI E-utilities (free, no API key required)
2. **bioRxiv search**: Queries bioRxiv's public REST API for preprints
3. **Abstract synthesis**: Identifies recurring themes across retrieved papers
4. **Citation graph**: Builds a JSON node-edge graph of internal citations
5. **Reproducibility bundle**: Exports `commands.sh`, `environment.yml`, SHA-256 checksums

## Scope

This skill searches literature and synthesises results. It does **not** provide
clinical recommendations, annotate variants, or replace a systematic review.

## Input Formats

| Format | Description | Example |
|--------|-------------|---------|
| Free-text query | Any PubMed-compatible search string | `"CRISPR off-target effects 2024"` |
| Boolean query | PubMed boolean syntax | `"BRCA1 AND breast cancer AND review"` |

## Workflow

1. **Parse query**: Accept free-text or PubMed boolean query
2. **Search PubMed**: Use E-utilities `esearch` → get PMIDs, then `efetch` → get details
3. **Search bioRxiv**: Query the public bioRxiv API, filter by keywords
4. **Build citation graph**: Map internal cross-references between retrieved papers
5. **Synthesise**: Identify recurring terms across abstracts
6. **Report**: Write `report.md` with paper summaries, citation graph, and reproducibility bundle

## CLI Reference

```bash
# Standard usage
python skills/lit-synthesizer/lit_synthesizer.py \
    --query "CRISPR off-target effects" \
    --output report/

# Limit results
python skills/lit-synthesizer/lit_synthesizer.py \
    --query "single cell RNA sequencing" \
    --max 5 \
    --output report/

# Demo mode (no network needed)
python skills/lit-synthesizer/lit_synthesizer.py \
    --demo --output /tmp/demo

# Via ClawBio runner
python clawbio.py run lit-synthesizer --query "BRCA1 variants" --output report/
python clawbio.py run lit-synthesizer --demo
```

## Demo

```bash
python clawbio.py run lit-synthesizer --demo
```

Expected output: A report covering 3 demo papers on CRISPR genome editing,
with a citation graph of 3 nodes and 3 edges, plus a full reproducibility bundle.

## Algorithm / Methodology

1. **E-utilities search** (`esearch`): POST query to NCBI, receive list of PMIDs
2. **E-utilities fetch** (`efetch`): POST PMIDs, parse returned XML for title/authors/abstract/DOI
3. **Rate limiting**: 0.34 s sleep between NCBI requests (respects 3 req/s limit)
4. **bioRxiv API**: GET `https://api.biorxiv.org/details/biorxiv/{date_range}/0/json`, filter by keywords
5. **Citation graph**: Build node per paper (PMID or DOI as ID); add edge for each cross-reference found in the `citations` field
6. **Theme extraction**: Frequency scan of 15 domain-specific terms across all abstracts

**Key parameters**:

- Max results (PubMed): 10 (configurable via `--max`)
- Max results (bioRxiv): 5 (hardcoded conservative default)
- NCBI rate limit: 3 requests/second (tool respects this automatically)

## Example Queries

- "Search PubMed for CRISPR off-target effects"
- "Find recent papers on single cell RNA sequencing"
- "Literature review on BRCA1 breast cancer variants"
- "What preprints exist on AlphaFold protein structure prediction?"

## Example Output

```
# 🦖 ClawBio Lit Synthesizer Report

**Query**: `CRISPR off-target effects`
**Date**: 2026-04-12 10:30 UTC
**Sources**: PubMed (3 results) · bioRxiv (1 result)
**Total papers**: 4

---

## Summary
Across 4 retrieved papers, recurring themes include: **crispr**, **off-target**,
**base editing**, **cas9**, **guide rna**, **variant**.
The literature spans 2024 to 2025.

---

## Papers

### 1. CRISPR-Cas9 off-target effects: detection and mitigation strategies

| Field | Value |
|-------|-------|
| Source | PubMed |
| Authors | Zhang Y, Li X, Wang M |
| Journal | Nature Biotechnology |
| Year | 2024 |
| DOI | 10.1038/nbt.2024.001 |

**Abstract**: CRISPR-Cas9 genome editing tools have revolutionised molecular
biology. However, off-target cleavage remains a major safety concern...
```

## Output Structure

```
output_directory/
├── report.md                         # Full synthesis report
├── results.json                      # All papers as structured JSON
├── citation_graph.json               # Node-edge citation graph
├── tables/
│   └── papers.csv                    # Tabular paper list
└── reproducibility/
    ├── commands.sh                   # Exact commands to reproduce
    ├── environment.yml               # Conda/pip environment
    └── checksums.sha256              # SHA-256 of all output files
```

## Dependencies

**Required**:

- `biopython >= 1.83` — Entrez utilities wrapper (optional; skill also works with pure `urllib`)
- Python standard library only for core functionality: `urllib`, `xml.etree`, `json`, `csv`, `hashlib`

**Optional**:

- `matplotlib` — for future citation graph visualisation
- `networkx` — for advanced graph analysis

## Gotchas

- **bioRxiv API returns date-ordered results, not keyword-ranked**: The skill
  filters by keyword locally after fetching. For very broad queries this may
  return zero bioRxiv results. Use a specific query to improve recall.

- **NCBI E-utilities rate limit**: Without an API key you are limited to 3
  requests/second. The skill enforces a 0.34 s sleep. Do NOT remove this sleep
  or you will receive HTTP 429 errors.

- **Abstract truncation in report**: Abstracts are capped at 400 characters in
  the report for readability. Full text is in `results.json`.

- **Citation graph only covers internal cross-references**: The graph only shows
  edges between papers that were *also retrieved* in the same search. It is not
  a global citation network.

## Safety

- **Local-first**: No user data is uploaded. Only the search query leaves the machine.
- **Disclaimer**: Every report includes the ClawBio research disclaimer.
- **Audit trail**: All operations logged to reproducibility bundle.
- **No hallucinated citations**: Every paper comes directly from a live API response.

## Agent Boundary

The agent (LLM) dispatches the query and explains results.
The skill (Python) executes the API calls and generates files.
The agent must NOT invent paper titles, authors, or DOIs.

## Integration with Bio Orchestrator

**Trigger conditions**: route here when the user mentions:

- `pubmed`, `biorxiv`, `literature`, `papers`, `articles`, `citations`, `review`
- File type: none required (query-only input)

**Chaining partners**:

- `pharmgx-reporter`: A lit search on a drug gene (e.g. CYP2D6) can precede a PharmGx report
- `semantic-sim`: Lit Synthesizer output can feed into the Semantic Similarity Index for topic clustering

## Maintenance

- **Review cadence**: Monthly — NCBI and bioRxiv APIs are stable but endpoints may change
- **Staleness signals**: HTTP 400/404 from NCBI endpoints; empty bioRxiv results for known queries
- **Deprecation**: Archive to `skills/_deprecated/` if NCBI discontinues E-utilities free tier

## Citations

- [NCBI E-utilities](https://www.ncbi.nlm.nih.gov/books/NBK25499/); PubMed programmatic search API
- [bioRxiv API](https://api.biorxiv.org/); preprint search and metadata
- [Biopython Entrez module](https://biopython.org/docs/latest/api/Bio.Entrez.html); Python wrapper for E-utilities
