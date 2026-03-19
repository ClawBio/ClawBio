# PubMed Summariser — Design Spec

**Date:** 2026-03-19
**Status:** Approved
**Skill path:** `skills/pubmed-summariser/`

---

## Overview

A ClawBio skill that accepts a gene name or disease term, queries PubMed via the NCBI Entrez API, and delivers a structured research briefing of the top recent papers. Both a terminal summary and an HTML report are produced.

---

## Goals

- Accept a free-text `--query` argument (gene name or disease term)
- Fetch the top 10 most recent English-language PubMed papers matching the query
- Output a readable terminal summary and a styled HTML report using `HtmlReportBuilder`
- Support `--demo` mode using `BRCA1` as the pre-set query (hardcoded — no file dependency)
- Be registered in the CLAUDE.md routing table and `skills/catalog.json`

---

## Non-Goals (MVP)

- MeSH term expansion
- Citation counts (elink API)
- Abstract word frequency visualisation
- Date-range filtering
- Free-text queries beyond gene names and disease terms
- Shell / Makefile alias
- Offline / pre-cached demo mode

---

## File Structure

```
skills/pubmed-summariser/
├── SKILL.md                  ← ClawBio skill metadata + methodology (see below)
├── pubmed_summariser.py      ← CLI entry point + report renderer
└── pubmed_api.py             ← NCBI API calls (esearch + efetch + XML parsing)
```

The demo query `"BRCA1"` is hardcoded in `pubmed_summariser.py` — no `demo_query.txt` file needed.

---

## SKILL.md

### Frontmatter

Follows `templates/SKILL-TEMPLATE.md`. Use `darwin` (not `macos`) to match project convention. No `author` or `license` fields (consistent with `gwas-lookup` convention in this repo).

```yaml
---
name: pubmed-summariser
description: >-
  Search PubMed for a gene name or disease term and generate a structured
  research briefing of the top recent English-language papers.
version: 0.1.0
metadata:
  openclaw:
    requires:
      bins:
        - python3
      env: []
      config: []
    always: false
    emoji: "📄"
    homepage: https://pubmed.ncbi.nlm.nih.gov/
    os: [darwin, linux]
    install:
      - kind: pip
        package: requests
        bins: []
    trigger_keywords:
      - pubmed
      - summarise papers
      - research briefing
      - papers about
      - recent studies
      - literature search pubmed
      - gene papers
      - disease papers
---
```

### Prose Body

The SKILL.md prose body must include these sections (content abbreviated here; expand during implementation):

**Role statement:** "You are **PubMed Summariser**, a specialised ClawBio agent for literature retrieval. Your role is to take a gene name or disease term, query PubMed, and return a structured briefing of the top recent papers."

**Why This Exists:** Without it, researchers manually search PubMed and read each abstract to stay current. With it, they get a formatted briefing in seconds.

**Core Capabilities:** (1) PubMed query via NCBI Entrez API, (2) structured paper extraction (title, authors, journal, date, abstract), (3) terminal summary + HTML report.

**Workflow:** (1) Receive `--query` or `--demo`, (2) call `esearch` for PMIDs, (3) call `efetch` for full records, (4) parse XML, (5) render output.

**Output Structure:** List of papers with title, authors, journal, date, abstract excerpt, PubMed URL.

**Dependencies:** `requests`, Python `xml.etree.ElementTree` (stdlib).

**Safety:** Standard ClawBio medical disclaimer on every report.

---

## Architecture

### pubmed_api.py

Responsible for all network I/O and data extraction. Exposes one public function:

```python
def fetch_papers(query: str, max_results: int = 10) -> list[dict]
```

Internally:
1. `esearch.fcgi` — retrieves PMIDs matching `query AND english[la]`, sorted by date, limited to `max_results`; includes `email=clawbio@example.com` and `tool=clawbio` per NCBI policy
2. `efetch.fcgi` — fetches full XML records for those PMIDs
3. XML parsing via `xml.etree.ElementTree` — extracts: title, authors, journal, publication date, abstract
4. Network timeout: `10` seconds on all requests

**Author formatting:** List up to 3 authors as "Last FM", then append "et al." if more exist. Single-author: just "Last FM".

**Abstract extraction:** Take the first sentence using this heuristic: split on the first occurrence of `. ` (period + space) followed by an uppercase letter, max 300 characters. If no match, use the first 300 characters. Empty string if no abstract.

Returns a list of dicts, one per paper:
```python
{
    "title": str,
    "authors": str,      # e.g. "Smith J, Lee K, et al."
    "journal": str,
    "date": str,         # e.g. "2025-11"
    "abstract": str,     # first sentence, max 300 chars; empty string if none
    "pmid": str,
    "url": str           # https://pubmed.ncbi.nlm.nih.gov/<pmid>/
}
```

### pubmed_summariser.py

Responsible for CLI parsing and rendering. Flow:

1. Parse args (`--query`, `--output`, `--max-results`, `--demo`)
2. If `--demo`: set query to `"BRCA1"` (hardcoded); if `--query` also supplied, print `"[demo mode] ignoring --query, using BRCA1"` and proceed
3. Clamp `--max-results` to 50 silently if exceeded, print `"[warning] --max-results capped at 50"`
4. Call `pubmed_api.fetch_papers(query, max_results)`
5. Print terminal summary
6. Write HTML report to `<output>/report.html` using `HtmlReportBuilder`; `write_html_report` handles `mkdir(parents=True, exist_ok=True)` internally — no manual directory creation needed

---

## CLI Interface

```bash
# Gene or disease query
python skills/pubmed-summariser/pubmed_summariser.py \
  --query BRCA1 --output /tmp/pubmed_demo

python skills/pubmed-summariser/pubmed_summariser.py \
  --query "type 2 diabetes" --output /tmp/pubmed_demo

# Adjust result count (default 10, max 50 — clamped silently)
python skills/pubmed-summariser/pubmed_summariser.py \
  --query BRCA1 --max-results 20 --output /tmp/pubmed_demo

# Demo mode — runs with BRCA1, ignores --query if also supplied
python skills/pubmed-summariser/pubmed_summariser.py \
  --demo --output /tmp/pubmed_demo
```

| Argument | Required | Default | Description |
|---|---|---|---|
| `--query` | yes (unless `--demo`) | — | Gene name or disease term |
| `--output` | yes | — | Directory to save HTML report (created automatically if absent) |
| `--max-results` | no | 10 | Number of papers; silently clamped to 50 with warning |
| `--demo` | no | — | Run with BRCA1; `--query` ignored if both supplied |

---

## Output

### Terminal
```
PubMed Research Briefing: BRCA1
================================
Found 10 papers (sorted by date, English only)

1. BRCA1-mediated DNA repair in triple-negative breast cancer
   Authors: Smith J, Lee K, et al.
   Journal: Nature Genetics | 2025-11
   Abstract: This study investigates the role of BRCA1 in homologous recombination...
   URL: https://pubmed.ncbi.nlm.nih.gov/12345678/
```

### HTML Report
Saved to `<output>/report.html`. Uses `HtmlReportBuilder` from `clawbio.common.html_report`:
- `add_header_block()` — branded green gradient header
- `add_metadata()` — query, result count, date generated
- One `add_section()` + `add_raw_html()` card block per paper (title as link, authors/journal/date, abstract excerpt)
- `add_disclaimer()` — standard ClawBio medical disclaimer
- `add_footer_block()` — branded footer with timestamp

CSS is fully embedded via `HtmlReportBuilder` — no external stylesheet needed.

---

## Error Handling

| Scenario | Behaviour |
|---|---|
| Zero results returned | Print `"No results found for query: <query>"` to terminal; exit code 0; no HTML written |
| Network timeout or HTTP error from NCBI | Print error message to stderr; exit code 1; no HTML written |
| `--query` missing and `--demo` not set | `argparse` exits with usage error |
| `--output` directory does not exist | `write_html_report` creates it automatically via `mkdir(parents=True, exist_ok=True)` |
| `--max-results` > 50 | Silently clamp to 50; print `"[warning] --max-results capped at 50"` |
| Both `--query` and `--demo` supplied | `--demo` wins; print `"[demo mode] ignoring --query, using BRCA1"` |

**Note:** `--demo` makes a live NCBI API call. Demo mode can fail in offline environments — this is acceptable for MVP.

---

## NCBI API Details

- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi` — retrieve PMIDs
- `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi` — fetch XML records
- Always include `tool=clawbio&email=clawbio@example.com` per NCBI E-utilities policy
- No API key required. Optional NCBI key increases rate limit from 3 → 10 req/s
- Network timeout: 10 seconds on all requests

---

## CLAUDE.md Routing

Add **alongside** (not replacing) the existing `lit-synthesizer` entry. `pubmed-summariser` handles PubMed-specific structured briefings; `lit-synthesizer` handles broader literature search methodology.

| User Intent | Skill | Action |
|---|---|---|
| PubMed search, "summarise PubMed papers about X", "recent papers on BRCA1", "research briefing", gene papers, disease papers | `skills/pubmed-summariser/` | Run `pubmed_summariser.py` |

---

## catalog.json Entry

Add this entry to the `skills` array in `skills/catalog.json`:

```json
{
  "name": "pubmed-summariser",
  "cli_alias": null,
  "description": "Search PubMed for a gene name or disease term and generate a structured research briefing of the top recent English-language papers.",
  "version": "0.1.0",
  "status": "mvp",
  "has_script": true,
  "has_tests": false,
  "has_demo": true,
  "demo_command": "python skills/pubmed-summariser/pubmed_summariser.py --demo --output /tmp/pubmed_demo",
  "dependencies": ["requests"],
  "tags": ["pubmed", "literature", "genomics", "research"],
  "trigger_keywords": ["pubmed", "summarise papers", "research briefing", "papers about", "recent studies", "gene papers"],
  "chaining_partners": ["lit-synthesizer", "gwas-lookup", "gwas-prs"]
}
```

---

## Disclaimer

> ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.
