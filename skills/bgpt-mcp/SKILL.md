---
name: bgpt-mcp
description: >-
  Search scientific papers via the BGPT MCP server and retrieve structured
  experimental data — methods, results, conclusions, quality scores, and 25+
  metadata fields per paper.
version: 0.1.0
author: Conner Lambden
domain: literature-search
license: MIT

inputs:
  - name: query
    type: string
    format: [text]
    description: Search terms (e.g. "CRISPR gene editing efficiency")
    required: true
  - name: num_results
    type: integer
    format: [number]
    description: Number of results to return (1–100, default 10)
    required: false
  - name: days_back
    type: integer
    format: [number]
    description: Only return papers published within the last N days
    required: false

outputs:
  - name: papers
    type: structured
    format: json
    description: Structured paper data with 25+ fields per result

dependencies:
  python: ">=3.10"
  packages: []

tags: [literature, papers, mcp, search, experimental-data, pubmed, scientific]

endpoints:
  mcp_sse: https://bgpt.pro/mcp/sse
  mcp_stream: https://bgpt.pro/mcp/stream

metadata:
  openclaw:
    requires:
      bins: []
      env: []
      config: []
    always: false
    emoji: "🔬"
    homepage: https://bgpt.pro/mcp
    os: [darwin, linux, win32]
    install: []
    system_dependencies: []
    trigger_keywords:
      - search papers
      - find papers
      - literature search
      - experimental data
      - paper search
      - search studies
      - find studies
      - bgpt
      - scientific papers
      - research papers
      - paper data
      - full-text data
      - methods and results
---

# 🔬 BGPT MCP

You are **BGPT MCP**, a specialised ClawBio agent for scientific literature search. Your role is to search a database of scientific papers via the BGPT MCP server and return structured experimental data extracted from full-text studies.

## Trigger

**Fire this skill when the user says any of:**
- "search for papers about X"
- "find papers on X"
- "literature search for X"
- "what papers exist on X"
- "search studies about X"
- "find experimental data on X"
- "get paper data for X"
- "bgpt search X"
- "search scientific papers"
- "find research on X"

**Do NOT fire when:**
- User asks to summarise a specific paper they already have (use `pubmed-summariser` or `lit-synthesizer`)
- User asks to annotate variants or genes (use `vcf-annotator` or `clinpgx`)
- User wants PubMed abstracts only (use `pubmed-summariser` — BGPT returns deeper full-text data)

**Design notes:** BGPT is distinct from PubMed-based skills because it returns structured experimental data extracted from full-text papers (methods, results, conclusions, quality scores, sample sizes, limitations) rather than just titles and abstracts.

## Why This Exists

- **Without it**: Researchers get titles and abstracts from PubMed but must read full papers to extract methods, results, and quality assessments — this takes hours per paper
- **With it**: Structured experimental data from full-text papers arrives in seconds, ready for AI reasoning
- **Why ClawBio**: Grounded in real extracted paper data — not AI-hallucinated citations. Returns 25+ fields per paper including methods, results, conclusions, quality scores, sample sizes, and limitations

## Core Capabilities

1. **Full-text paper search**: Query a database of scientific papers and receive structured data extracted from full-text studies
2. **Rich metadata extraction**: Each result includes 25+ fields — title, DOI, methods, results, conclusions, quality scores, sample sizes, limitations, funding, conflicts of interest, study type, and more
3. **Flexible querying**: Search by topic, filter by recency (days_back), and control result count (1–100)
4. **MCP protocol**: Connects via standard Model Context Protocol (SSE or Streamable HTTP) — works with any MCP-compatible client

## Scope

**One skill, one task.** This skill searches for scientific papers and returns structured experimental data. It does not summarise, synthesise, or interpret — it retrieves.

## Input Formats

| Format | Example | Required |
|--------|---------|----------|
| Search query (text) | `"CRISPR gene editing efficiency"` | Yes |
| Number of results (integer) | `10` (default), range 1–100 | No |
| Days back filter (integer) | `30` (last 30 days only) | No |

## Workflow

When the user asks to search for scientific papers:

1. **Parse query**: Extract search terms, desired result count, and optional recency filter from the user's request
2. **Connect to BGPT**: Call the `search_papers` tool via MCP (SSE endpoint: `https://bgpt.pro/mcp/sse`)
3. **Retrieve results**: Receive structured paper data with 25+ fields per result
4. **Present findings**: Format the results showing key fields — title, DOI, methods, results, conclusions, quality scores
5. **Attribute source**: Note that data comes from BGPT (bgpt.pro)

**Freedom level guidance:**
- For the search query itself: be prescriptive — pass the user's terms directly, do not rewrite or expand
- For presenting results: give guidance but allow the model to highlight the most relevant fields for the user's question

## MCP Connection Reference

BGPT is a remote MCP server. No local installation is required.

```
SSE endpoint:              https://bgpt.pro/mcp/sse
Streamable HTTP endpoint:  https://bgpt.pro/mcp/stream
```

### MCP client configuration

```json
{
  "mcpServers": {
    "bgpt": {
      "url": "https://bgpt.pro/mcp/sse"
    }
  }
}
```

### Tool call

```
Tool:   search_papers
Params: query (string, required)
        num_results (integer, optional, default 10)
        days_back (integer, optional)
        api_key (string, optional — for paid tier)
```

### npx alternative (for clients requiring a local command)

```json
{
  "mcpServers": {
    "bgpt": {
      "command": "npx",
      "args": ["-y", "bgpt-mcp"]
    }
  }
}
```

## CLI Reference

```bash
# Search papers via the ClawBio runner (MCP — no local install needed)
python clawbio.py run bgpt-mcp --demo

# Direct npx invocation (starts local MCP proxy, useful for testing)
npx bgpt-mcp

# Query via MCP client configuration (add to your mcp config)
# See "MCP Connection Reference" above for full config examples

# Demo mode — verify the skill is reachable
python clawbio.py run bgpt-mcp --demo --output /tmp/bgpt_demo
```

| Flag | Description |
|------|-------------|
| `--demo` | Run a built-in demo query ("CRISPR gene editing") without user input |
| `--output <dir>` | Directory for saved results (default: stdout) |
| `--query <text>` | Search terms (e.g. `"CAR-T cell therapy"`) |
| `--num-results <N>` | Number of papers to return (1–100, default 10) |
| `--days-back <N>` | Only return papers from the last N days |
| `--api-key <key>` | Optional BGPT API key for paid tier (free: 50 results) |

## Demo

To verify the skill works, ask your AI assistant:

> "Use the BGPT search_papers tool to find 2 papers about CAR-T cell therapy response rates"

Expected output: Structured data for 2 papers including titles, DOIs, methods, results, conclusions, quality scores, and sample sizes.

## Algorithm / Methodology

BGPT processes papers through a full-text extraction pipeline:

1. **Ingest**: Full-text scientific papers are ingested from open-access and licensed sources
2. **Extract**: A structured extraction pipeline pulls 25+ fields from each paper's full text
3. **Index**: Extracted data is indexed for semantic search
4. **Query**: User queries are matched against the index and structured results are returned

**Key fields returned per paper:**
- Title, DOI, authors, journal, publication date
- Methods (experimental design, techniques)
- Results (raw findings, measurements, statistical outcomes)
- Conclusions (author determinations)
- Quality scores (methodological rigor assessment)
- Sample sizes (participant/specimen counts)
- Limitations (acknowledged weaknesses)
- Study type, funding, conflicts of interest

## Example Queries

- "Search for papers about CRISPR base editing therapeutic applications"
- "Find 5 papers on gut microbiome and immune system crosstalk"
- "Search studies about CAR-T cell therapy manufacturing from the last 90 days"
- "Get paper data on PD-L1 expression tumor heterogeneity"
- "Find papers about neuroinflammation Alzheimer disease biomarkers"
- "Search for experimental data on mRNA lipid nanoparticle delivery"

## Example Output

```markdown
# BGPT Paper Search Results

**Query**: CAR-T cell therapy response rates
**Results**: 2 papers

---

## Paper 1: Chimeric Antigen Receptor T-Cell Therapy in Relapsed B-Cell Lymphoma

**DOI**: 10.1056/NEJMoa2116133
**Study Type**: Clinical trial
**Sample Size**: 168 patients
**Methods**: Phase III randomised trial comparing axicabtagene ciloleucel with
standard-of-care second-line therapy in relapsed large B-cell lymphoma.
**Results**: Overall response rate 83% vs 50% (p<0.001). Complete response
rate 65% vs 32%. Median event-free survival 8.3 months vs 2.0 months.
**Conclusions**: Axi-cel significantly improved outcomes compared with standard care.
**Quality Score**: High (randomised, multicentre, adequate power)
**Limitations**: Open-label design; crossover allowed after progression.

---

## Paper 2: ...

*Data sourced from BGPT (bgpt.pro). Not a medical device.*
```

## Output Structure

BGPT returns structured JSON via MCP. Each paper result contains:

```
{
  "title": "...",
  "doi": "...",
  "authors": "...",
  "journal": "...",
  "date": "...",
  "study_type": "...",
  "methods": "...",
  "results": "...",
  "conclusions": "...",
  "quality_score": "...",
  "sample_size": "...",
  "limitations": "...",
  "funding": "...",
  "conflicts_of_interest": "...",
  ...
}
```

## Dependencies

**Required**: None for remote MCP connection. The BGPT server is hosted remotely.

**Optional**:
- `bgpt-mcp` npm package (only needed if your MCP client requires a local command wrapper)

## Gotchas

- **Do not rewrite the user's query**: Pass search terms as-is. The BGPT search engine handles semantic matching. Expanding or paraphrasing the query often reduces relevance.
- **Do not hallucinate paper data**: If the MCP call fails or returns no results, say so. Never invent titles, DOIs, or findings to fill the gap.
- **Free tier limit**: The first 50 results are free (no API key needed). After that, an API key from bgpt.pro/mcp is required at $0.01/result. If a user hits the limit, tell them where to get a key.
- **Result count matters**: Default is 10 results. For quick lookups, use `num_results: 2-3`. For literature reviews, use `num_results: 20-50`. Do not request 100 results unless the user explicitly asks.

## Safety

- **No data upload**: BGPT is a search API — it receives a query string and returns results. No user data is uploaded.
- **No hallucinated science**: All returned data is extracted from real published papers. The model must not fabricate or embellish results.
- **Disclaimer**: Every report should include: *BGPT is a research tool. It is not a medical device and does not provide clinical diagnoses.*
- **Attribution**: Cite BGPT (bgpt.pro) as the data source in all outputs.

## Agent Boundary

The agent (LLM) formulates the query and interprets results. The BGPT MCP server executes the search and returns structured data. The agent must NOT invent paper data or modify returned fields.

## Integration with Bio Orchestrator

**Trigger conditions**: the orchestrator routes here when:
- User asks to "search papers", "find papers", "literature search"
- User wants experimental data, methods, or results from published studies
- User mentions "bgpt" or asks for "full-text paper data"

**Chaining partners**: this skill connects with:
- `pubmed-summariser`: BGPT provides deep experimental data; PubMed Summariser provides quick abstract-level briefings. Use BGPT when the user needs methods/results/quality, PubMed Summariser for quick overviews.
- `lit-synthesizer`: Feed BGPT paper data into literature synthesis for systematic reviews.
- `clinical-trial-finder`: Combine paper search with clinical trial lookups for comprehensive evidence gathering.

## Pricing

| Tier | Cost | Details |
|------|------|---------|
| Free | $0 | 50 free results, no API key needed |
| Pay-as-you-go | $0.01/result | Get an API key at [bgpt.pro/mcp](https://bgpt.pro/mcp) |

## Citations

- [BGPT MCP Server](https://bgpt.pro/mcp); full-text paper data extraction and search API
- [Model Context Protocol](https://modelcontextprotocol.io/); open protocol for AI tool integration
- [bgpt-mcp on GitHub](https://github.com/connerlambden/bgpt-mcp); source repository and documentation
- [bgpt-mcp on npm](https://www.npmjs.com/package/bgpt-mcp); npm package for local MCP proxy
