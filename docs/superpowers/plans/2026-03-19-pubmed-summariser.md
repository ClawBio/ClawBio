# PubMed Summariser Implementation Plan

> **For agentic workers:** REQUIRED SUB-SKILL: Use superpowers:subagent-driven-development (recommended) or superpowers:executing-plans to implement this plan task-by-task. Steps use checkbox (`- [ ]`) syntax for tracking.

**Goal:** Build a ClawBio skill that queries PubMed for a gene name or disease term and produces a terminal summary + HTML report of the top recent papers.

**Architecture:** Two-file split — `pubmed_api.py` handles all NCBI Entrez API calls and XML parsing; `pubmed_summariser.py` handles CLI parsing and report rendering via `HtmlReportBuilder`. Clean separation: one file talks to the internet, one file formats output.

**Tech Stack:** Python 3.10+, `requests`, stdlib `xml.etree.ElementTree`, `argparse`, `clawbio.common.html_report.HtmlReportBuilder`

---

## Files to Create or Modify

| Action | Path | Responsibility |
|---|---|---|
| Create | `skills/pubmed-summariser/SKILL.md` | Skill metadata + methodology for AI routing |
| Create | `skills/pubmed-summariser/pubmed_api.py` | NCBI esearch + efetch + XML parsing |
| Create | `skills/pubmed-summariser/pubmed_summariser.py` | CLI entry point + terminal/HTML rendering |
| Create | `skills/pubmed-summariser/tests/__init__.py` | Makes tests a package |
| Create | `skills/pubmed-summariser/tests/fixtures/sample_papers.json` | Pre-fetched fixture for unit tests |
| Create | `skills/pubmed-summariser/tests/test_pubmed_api.py` | Unit tests for API parsing logic |
| Create | `skills/pubmed-summariser/tests/test_pubmed_summariser.py` | Unit tests for CLI + rendering |
| Modify | `pytest.ini` | Register new test path |
| Modify | `CLAUDE.md` | Add routing row for pubmed-summariser |
| Modify | `skills/catalog.json` | Register skill in catalog |

---

## Task 1: Create the skill directory and SKILL.md

**Files:**
- Create: `skills/pubmed-summariser/SKILL.md`

- [ ] **Step 1.1: Create the skill directory**

```bash
mkdir -p skills/pubmed-summariser
```

- [ ] **Step 1.2: Write SKILL.md**

Create `skills/pubmed-summariser/SKILL.md` with this exact content:

```markdown
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

# 📄 PubMed Summariser

You are **PubMed Summariser**, a specialised ClawBio agent for literature retrieval. Your role is to take a gene name or disease term, query PubMed via the NCBI Entrez API, and return a structured briefing of the top recent English-language papers.

## Why This Exists

- **Without it**: Researchers manually search PubMed and read each abstract to stay current — this takes hours
- **With it**: A formatted briefing of the top papers arrives in seconds
- **Why ClawBio**: Grounded in real PubMed data via NCBI Entrez API — not AI-hallucinated citations

## Core Capabilities

1. **PubMed query**: Search by gene name (e.g. `BRCA1`) or disease term (e.g. `type 2 diabetes`)
2. **Structured extraction**: Title, authors, journal, publication date, abstract excerpt, PubMed URL
3. **Dual output**: Terminal summary for quick review + HTML report for sharing

## Input Formats

| Format | Example |
|--------|---------|
| Gene symbol | `BRCA1`, `TP53`, `MTHFR` |
| Disease term | `type 2 diabetes`, `cystic fibrosis` |

## Workflow

When the user asks to summarise PubMed papers about a gene or disease:

1. **Receive query**: `--query <term>` or `--demo` (uses BRCA1)
2. **esearch**: Query `https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi` for PMIDs
3. **efetch**: Fetch full XML records for those PMIDs
4. **Parse XML**: Extract title, authors, journal, date, abstract
5. **Render output**: Print terminal summary and write `report.html`

## Algorithm / Methodology

- Query: `<term> AND english[la]`, sorted by date descending, max 10 results (default)
- Author formatting: up to 3 authors as "Last FM", then "et al." if more exist
- Abstract: first sentence heuristic — split on `. ` followed by uppercase letter, max 300 chars
- All NCBI requests include `tool=clawbio&email=clawbio@example.com` per NCBI E-utilities policy
- Network timeout: 10 seconds

## Output Structure

```
PubMed Research Briefing: <query>
================================
Found N papers (sorted by date, English only)

1. <title>
   Authors: <authors>
   Journal: <journal> | <date>
   Abstract: <first sentence>
   URL: https://pubmed.ncbi.nlm.nih.gov/<pmid>/
```

HTML report saved to `<output>/report.html`.

## Dependencies

- `requests` (HTTP)
- `xml.etree.ElementTree` (stdlib — XML parsing)
- `clawbio.common.html_report.HtmlReportBuilder` (HTML rendering)

## Safety

Every report includes the standard ClawBio medical disclaimer:
> ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.

## Integration with Bio Orchestrator

Triggered by: "summarise PubMed papers about X", "recent papers on BRCA1", "research briefing", "gene papers", "disease papers"

Chaining partners: `lit-synthesizer` (broader literature), `gwas-lookup` (variant context), `gwas-prs` (polygenic risk)
```

- [ ] **Step 1.3: Verify the file exists**

```bash
ls -la skills/pubmed-summariser/SKILL.md
```

Expected: file listed, non-zero size.

- [ ] **Step 1.4: Commit**

```bash
git add skills/pubmed-summariser/SKILL.md
git commit -m "feat(pubmed-summariser): add SKILL.md with frontmatter and methodology"
```

---

## Task 2: Write pubmed_api.py with tests (TDD)

**Files:**
- Create: `skills/pubmed-summariser/tests/__init__.py`
- Create: `skills/pubmed-summariser/tests/fixtures/sample_papers.json`
- Create: `skills/pubmed-summariser/tests/test_pubmed_api.py`
- Create: `skills/pubmed-summariser/pubmed_api.py`

### 2a — Set up test infrastructure

- [ ] **Step 2a.1: Create the tests package**

```bash
mkdir -p skills/pubmed-summariser/tests/fixtures
touch skills/pubmed-summariser/tests/__init__.py
```

- [ ] **Step 2a.2: Create the sample fixture**

This fixture simulates a parsed paper dict — what `_parse_article()` should return. Save as `skills/pubmed-summariser/tests/fixtures/sample_papers.json`:

```json
[
  {
    "title": "BRCA1 and BRCA2: different roles in a common pathway of genome protection",
    "authors": "Prakash R, Zhang Y, et al.",
    "journal": "Nat Rev Mol Cell Biol",
    "date": "2015-08",
    "abstract": "BRCA1 and BRCA2 are involved in the repair of DNA double-strand breaks.",
    "pmid": "26130008",
    "url": "https://pubmed.ncbi.nlm.nih.gov/26130008/"
  },
  {
    "title": "Single-author paper example",
    "authors": "Smith J",
    "journal": "Cell",
    "date": "2024-01",
    "abstract": "This is the first sentence.",
    "pmid": "99999999",
    "url": "https://pubmed.ncbi.nlm.nih.gov/99999999/"
  }
]
```

### 2b — Write the failing tests

- [ ] **Step 2b.1: Write test_pubmed_api.py**

Create `skills/pubmed-summariser/tests/test_pubmed_api.py`:

```python
"""
Unit tests for pubmed_api.py — no network calls required.

All tests operate on pre-built XML strings or fixture data.
"""

import json
import sys
from pathlib import Path
import xml.etree.ElementTree as ET

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

FIXTURES = Path(__file__).parent / "fixtures"


def load_fixture():
    return json.loads((FIXTURES / "sample_papers.json").read_text())


# ── _format_authors ────────────────────────────────────────────────────────────


def test_format_authors_single():
    """Single author: return 'Last FM' with no et al."""
    from pubmed_api import _format_authors
    assert _format_authors(["Smith J"]) == "Smith J"


def test_format_authors_two():
    """Two authors: return both joined by ', '."""
    from pubmed_api import _format_authors
    assert _format_authors(["Smith J", "Lee K"]) == "Smith J, Lee K"


def test_format_authors_three():
    """Three authors: return all three joined."""
    from pubmed_api import _format_authors
    assert _format_authors(["Smith J", "Lee K", "Doe A"]) == "Smith J, Lee K, Doe A"


def test_format_authors_four_or_more():
    """Four+ authors: return first 3 + 'et al.'"""
    from pubmed_api import _format_authors
    result = _format_authors(["Smith J", "Lee K", "Doe A", "Brown B"])
    assert result == "Smith J, Lee K, Doe A, et al."


def test_format_authors_empty():
    """Empty author list: return empty string."""
    from pubmed_api import _format_authors
    assert _format_authors([]) == ""


# ── _first_sentence ─────────────────────────────────────────────────────────


def test_first_sentence_basic():
    """Extracts first sentence from normal abstract."""
    from pubmed_api import _first_sentence
    text = "This is the first sentence. This is the second sentence."
    assert _first_sentence(text) == "This is the first sentence."


def test_first_sentence_no_period():
    """No period found: return up to 300 chars."""
    from pubmed_api import _first_sentence
    text = "A" * 400
    result = _first_sentence(text)
    assert len(result) <= 300


def test_first_sentence_empty():
    """Empty string returns empty string."""
    from pubmed_api import _first_sentence
    assert _first_sentence("") == ""


def test_first_sentence_abbreviation_not_split():
    """Should not split on 'et al. ' (lowercase after period+space)."""
    from pubmed_api import _first_sentence
    text = "Smith et al. demonstrated this. The second sentence follows."
    # 'et al. d...' — 'd' is lowercase, so no split there
    # Should split at '. T' (capital T)
    result = _first_sentence(text)
    assert result == "Smith et al. demonstrated this."


def test_first_sentence_max_300():
    """Result never exceeds 300 chars even if first sentence is long."""
    from pubmed_api import _first_sentence
    text = ("Word " * 100) + ". Next sentence."
    result = _first_sentence(text)
    assert len(result) <= 300


# ── _parse_article ──────────────────────────────────────────────────────────


def _make_article_xml(
    pmid="12345678",
    title="Test Title",
    authors=None,
    journal="Test Journal",
    year="2024",
    month="06",
    abstract="This is the abstract. Second sentence.",
):
    """Build a minimal PubmedArticle XML element for testing."""
    if authors is None:
        authors = [("Smith", "J")]
    author_xml = ""
    for last, fore in authors:
        author_xml += f"""
        <Author>
          <LastName>{last}</LastName>
          <ForeName>{fore}</ForeName>
          <Initials>{fore[0]}</Initials>
        </Author>"""
    return ET.fromstring(f"""
    <PubmedArticle>
      <MedlineCitation>
        <PMID>{pmid}</PMID>
        <Article>
          <ArticleTitle>{title}</ArticleTitle>
          <AuthorList>{author_xml}</AuthorList>
          <Journal>
            <Title>{journal}</Title>
            <JournalIssue>
              <PubDate>
                <Year>{year}</Year>
                <Month>{month}</Month>
              </PubDate>
            </JournalIssue>
          </Journal>
          <Abstract>
            <AbstractText>{abstract}</AbstractText>
          </Abstract>
        </Article>
      </MedlineCitation>
    </PubmedArticle>
    """)


def test_parse_article_returns_all_fields():
    """Parsed article must contain all required dict keys."""
    from pubmed_api import _parse_article
    article = _make_article_xml()
    result = _parse_article(article)
    assert set(result.keys()) == {"title", "authors", "journal", "date", "abstract", "pmid", "url"}


def test_parse_article_pmid_url():
    """URL must be constructed from PMID."""
    from pubmed_api import _parse_article
    article = _make_article_xml(pmid="99887766")
    result = _parse_article(article)
    assert result["pmid"] == "99887766"
    assert result["url"] == "https://pubmed.ncbi.nlm.nih.gov/99887766/"


def test_parse_article_date_format():
    """Date must be formatted as YYYY-MM."""
    from pubmed_api import _parse_article
    article = _make_article_xml(year="2023", month="11")
    result = _parse_article(article)
    assert result["date"] == "2023-11"


def test_parse_article_no_abstract():
    """Missing abstract element returns empty string."""
    from pubmed_api import _parse_article
    xml_str = """
    <PubmedArticle>
      <MedlineCitation>
        <PMID>11111111</PMID>
        <Article>
          <ArticleTitle>No abstract</ArticleTitle>
          <AuthorList>
            <Author><LastName>Doe</LastName><Initials>A</Initials></Author>
          </AuthorList>
          <Journal>
            <Title>Some Journal</Title>
            <JournalIssue><PubDate><Year>2020</Year></PubDate></JournalIssue>
          </Journal>
        </Article>
      </MedlineCitation>
    </PubmedArticle>
    """
    result = _parse_article(ET.fromstring(xml_str))
    assert result["abstract"] == ""
```

- [ ] **Step 2b.2: Run tests to confirm they all FAIL (pubmed_api.py doesn't exist yet)**

```bash
pytest skills/pubmed-summariser/tests/test_pubmed_api.py -v 2>&1 | head -30
```

Expected: `ModuleNotFoundError: No module named 'pubmed_api'` — confirms the tests are wired up correctly.

### 2c — Implement pubmed_api.py

- [ ] **Step 2c.1: Write pubmed_api.py**

Create `skills/pubmed-summariser/pubmed_api.py`:

```python
"""
pubmed_api.py — NCBI Entrez API client for PubMed Summariser.

Exposes one public function:
    fetch_papers(query, max_results=10) -> list[dict]

All network I/O is here. No rendering, no CLI.
"""

from __future__ import annotations

import re
import sys
import xml.etree.ElementTree as ET
from pathlib import Path

import requests

# ---------------------------------------------------------------------------
# Project root on sys.path (required to import clawbio.common)
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

# ---------------------------------------------------------------------------
# NCBI Entrez API constants
# ---------------------------------------------------------------------------
_ESEARCH = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"
_EFETCH  = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi"
_PARAMS  = {"tool": "clawbio", "email": "clawbio@example.com"}
_TIMEOUT = 10  # seconds


# ---------------------------------------------------------------------------
# Internal helpers
# ---------------------------------------------------------------------------

def _format_authors(authors: list[str]) -> str:
    """Format author list: up to 3 names, then 'et al.' if more."""
    if not authors:
        return ""
    if len(authors) <= 3:
        return ", ".join(authors)
    return ", ".join(authors[:3]) + ", et al."


def _first_sentence(text: str) -> str:
    """
    Extract the first sentence from an abstract.

    Heuristic: split on the first '. ' followed by an uppercase letter.
    Falls back to first 300 characters if no match.
    Result is always capped at 300 characters.
    """
    if not text:
        return ""
    match = re.search(r'\.\s+(?=[A-Z])', text)
    if match:
        sentence = text[:match.start() + 1]
    else:
        sentence = text[:300]
    return sentence[:300]


def _parse_article(article: ET.Element) -> dict:
    """
    Parse a <PubmedArticle> XML element into a dict.

    Returns keys: title, authors, journal, date, abstract, pmid, url
    """
    mc = article.find("MedlineCitation")

    # PMID
    pmid_el = mc.find("PMID") if mc is not None else None
    pmid = pmid_el.text.strip() if pmid_el is not None and pmid_el.text else ""

    art = mc.find("Article") if mc is not None else None

    # Title
    title_el = art.find("ArticleTitle") if art is not None else None
    title = (title_el.text or "").strip() if title_el is not None else ""

    # Authors — each as "LastName Initials"
    raw_authors: list[str] = []
    author_list = art.find("AuthorList") if art is not None else None
    if author_list is not None:
        for author in author_list.findall("Author"):
            last = author.findtext("LastName", "").strip()
            initials = author.findtext("Initials", "").strip()
            if last:
                raw_authors.append(f"{last} {initials}".strip())

    # Journal
    journal_el = art.find("Journal") if art is not None else None
    journal = ""
    date = ""
    if journal_el is not None:
        journal = (journal_el.findtext("Title") or "").strip()
        pub_date = journal_el.find(".//PubDate")
        if pub_date is not None:
            year = pub_date.findtext("Year", "").strip()
            month = pub_date.findtext("Month", "").strip()
            date = f"{year}-{month}" if month else year

    # Abstract
    abstract_el = art.find(".//AbstractText") if art is not None else None
    raw_abstract = (abstract_el.text or "").strip() if abstract_el is not None else ""
    abstract = _first_sentence(raw_abstract)

    return {
        "title": title,
        "authors": _format_authors(raw_authors),
        "journal": journal,
        "date": date,
        "abstract": abstract,
        "pmid": pmid,
        "url": f"https://pubmed.ncbi.nlm.nih.gov/{pmid}/",
    }


# ---------------------------------------------------------------------------
# Public API
# ---------------------------------------------------------------------------

def fetch_papers(query: str, max_results: int = 10) -> list[dict]:
    """
    Query PubMed and return structured paper dicts.

    Args:
        query:       Gene name or disease term (e.g. 'BRCA1', 'type 2 diabetes')
        max_results: Number of papers to return (default 10, max 50)

    Returns:
        List of dicts with keys: title, authors, journal, date, abstract, pmid, url
        Empty list if no results.

    Raises:
        requests.RequestException: on network failure or HTTP error
    """
    # Step 1: esearch — get PMIDs
    search_params = {
        **_PARAMS,
        "db": "pubmed",
        "term": f"{query} AND english[la]",
        "retmax": max_results,
        "sort": "date",
        "retmode": "json",
    }
    resp = requests.get(_ESEARCH, params=search_params, timeout=_TIMEOUT)
    resp.raise_for_status()
    data = resp.json()
    pmids = data.get("esearchresult", {}).get("idlist", [])

    if not pmids:
        return []

    # Step 2: efetch — fetch full XML records
    fetch_params = {
        **_PARAMS,
        "db": "pubmed",
        "id": ",".join(pmids),
        "rettype": "xml",
        "retmode": "xml",
    }
    resp = requests.get(_EFETCH, params=fetch_params, timeout=_TIMEOUT)
    resp.raise_for_status()

    # Step 3: parse XML
    root = ET.fromstring(resp.content)
    papers = []
    for article in root.findall("PubmedArticle"):
        parsed = _parse_article(article)
        if parsed["pmid"]:  # skip malformed articles
            papers.append(parsed)

    return papers
```

- [ ] **Step 2c.2: Run the tests — they should all pass now**

```bash
pytest skills/pubmed-summariser/tests/test_pubmed_api.py -v
```

Expected: all tests PASS. If any fail, read the error and fix `pubmed_api.py` before proceeding.

- [ ] **Step 2c.3: Commit**

```bash
git add skills/pubmed-summariser/pubmed_api.py \
        skills/pubmed-summariser/tests/__init__.py \
        skills/pubmed-summariser/tests/fixtures/sample_papers.json \
        skills/pubmed-summariser/tests/test_pubmed_api.py
git commit -m "feat(pubmed-summariser): add pubmed_api.py with full test coverage"
```

---

## Task 3: Write pubmed_summariser.py with tests (TDD)

**Files:**
- Create: `skills/pubmed-summariser/tests/test_pubmed_summariser.py`
- Create: `skills/pubmed-summariser/pubmed_summariser.py`

### 3a — Write the failing tests first

- [ ] **Step 3a.1: Write test_pubmed_summariser.py**

Create `skills/pubmed-summariser/tests/test_pubmed_summariser.py`:

```python
"""
Unit tests for pubmed_summariser.py — no network calls required.

Tests the CLI logic, max_results clamping, demo mode behaviour, and
HTML report generation using a mock for pubmed_api.fetch_papers.
"""

import json
import sys
from io import StringIO
from pathlib import Path
from unittest.mock import patch

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

FIXTURES = Path(__file__).parent / "fixtures"

SAMPLE_PAPERS = json.loads((FIXTURES / "sample_papers.json").read_text())


# ── _clamp_max_results ───────────────────────────────────────────────────────


def test_clamp_within_limit():
    """Values <= 50 should be returned unchanged."""
    from pubmed_summariser import _clamp_max_results
    assert _clamp_max_results(10) == 10
    assert _clamp_max_results(50) == 50


def test_clamp_above_limit(capsys):
    """Values > 50 should be clamped to 50 and print a warning."""
    from pubmed_summariser import _clamp_max_results
    result = _clamp_max_results(200)
    assert result == 50
    captured = capsys.readouterr()
    assert "[warning]" in captured.out
    assert "capped at 50" in captured.out


# ── _format_terminal_summary ─────────────────────────────────────────────────


def test_terminal_summary_contains_query():
    """Terminal summary header must include the query string."""
    from pubmed_summariser import _format_terminal_summary
    output = _format_terminal_summary("BRCA1", SAMPLE_PAPERS)
    assert "BRCA1" in output


def test_terminal_summary_contains_paper_title():
    """Each paper title must appear in the terminal output."""
    from pubmed_summariser import _format_terminal_summary
    output = _format_terminal_summary("BRCA1", SAMPLE_PAPERS)
    assert SAMPLE_PAPERS[0]["title"] in output


def test_terminal_summary_empty_results():
    """Zero results should produce a 'No results found' message."""
    from pubmed_summariser import _format_terminal_summary
    output = _format_terminal_summary("UNKNOWN_GENE_XYZ", [])
    assert "No results found" in output


# ── _build_html_report ────────────────────────────────────────────────────────


def test_html_report_contains_query():
    """HTML report must include the query term."""
    from pubmed_summariser import _build_html_report
    html = _build_html_report("BRCA1", SAMPLE_PAPERS)
    assert "BRCA1" in html


def test_html_report_contains_pmid_link():
    """HTML report must contain a link to each paper's PubMed URL."""
    from pubmed_summariser import _build_html_report
    html = _build_html_report("BRCA1", SAMPLE_PAPERS)
    assert SAMPLE_PAPERS[0]["url"] in html


def test_html_report_contains_disclaimer():
    """HTML report must contain the medical disclaimer."""
    from pubmed_summariser import _build_html_report
    html = _build_html_report("BRCA1", SAMPLE_PAPERS)
    assert "research and educational tool" in html.lower() or "disclaimer" in html.lower()


# ── main integration (mocked) ─────────────────────────────────────────────────


def test_main_demo_mode_ignores_query(tmp_path, capsys):
    """When --demo is set alongside --query, demo wins and a note is printed."""
    from pubmed_summariser import main
    with patch("pubmed_summariser.pubmed_api.fetch_papers", return_value=SAMPLE_PAPERS) as mock_fetch:
        main(["--demo", "--query", "TP53", "--output", str(tmp_path)])
        # fetch_papers must have been called with "BRCA1", not "TP53"
        call_query = mock_fetch.call_args[0][0]
        assert call_query == "BRCA1"
    captured = capsys.readouterr()
    assert "[demo mode]" in captured.out


def test_main_writes_report_html(tmp_path):
    """Running with --query must produce report.html in the output directory."""
    from pubmed_summariser import main
    with patch("pubmed_summariser.pubmed_api.fetch_papers", return_value=SAMPLE_PAPERS):
        main(["--query", "BRCA1", "--output", str(tmp_path)])
    assert (tmp_path / "report.html").exists()


def test_main_no_results_exits_cleanly(tmp_path, capsys):
    """Zero API results must not write a report and must print a helpful message."""
    from pubmed_summariser import main
    with patch("pubmed_summariser.pubmed_api.fetch_papers", return_value=[]):
        main(["--query", "UNKNOWNXYZ", "--output", str(tmp_path)])
    assert not (tmp_path / "report.html").exists()
    captured = capsys.readouterr()
    assert "No results found" in captured.out
```

- [ ] **Step 3a.2: Run tests to confirm they FAIL**

```bash
pytest skills/pubmed-summariser/tests/test_pubmed_summariser.py -v 2>&1 | head -20
```

Expected: `ModuleNotFoundError: No module named 'pubmed_summariser'`

### 3b — Implement pubmed_summariser.py

- [ ] **Step 3b.1: Write pubmed_summariser.py**

Create `skills/pubmed-summariser/pubmed_summariser.py`:

```python
#!/usr/bin/env python3
"""
PubMed Summariser — ClawBio skill for literature retrieval.

Queries PubMed for a gene name or disease term and produces:
  - A terminal summary of the top recent papers
  - An HTML report saved to --output/report.html

Usage:
    python pubmed_summariser.py --query BRCA1 --output /tmp/pubmed_demo
    python pubmed_summariser.py --query "type 2 diabetes" --output /tmp/demo
    python pubmed_summariser.py --demo --output /tmp/pubmed_demo
"""

from __future__ import annotations

import argparse
import sys
from pathlib import Path

# ---------------------------------------------------------------------------
# Project root on sys.path (required to import clawbio.common)
# ---------------------------------------------------------------------------
_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

import pubmed_api
from clawbio.common.html_report import HtmlReportBuilder, write_html_report

# ---------------------------------------------------------------------------
# Constants
# ---------------------------------------------------------------------------
DEMO_QUERY = "BRCA1"
MAX_RESULTS_HARD_CAP = 50
SKILL_VERSION = "0.1.0"


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def _clamp_max_results(n: int) -> int:
    """Clamp max_results to 50; print a warning if clamped."""
    if n > MAX_RESULTS_HARD_CAP:
        print(f"[warning] --max-results capped at {MAX_RESULTS_HARD_CAP} (you asked for {n})")
        return MAX_RESULTS_HARD_CAP
    return n


def _format_terminal_summary(query: str, papers: list[dict]) -> str:
    """Build the terminal summary string."""
    lines = []
    header = f"PubMed Research Briefing: {query}"
    lines.append(header)
    lines.append("=" * len(header))

    if not papers:
        lines.append(f"No results found for query: {query}")
        return "\n".join(lines)

    lines.append(f"Found {len(papers)} papers (sorted by date, English only)\n")
    for i, paper in enumerate(papers, 1):
        lines.append(f"{i}. {paper['title']}")
        lines.append(f"   Authors: {paper['authors']}")
        lines.append(f"   Journal: {paper['journal']} | {paper['date']}")
        if paper["abstract"]:
            lines.append(f"   Abstract: {paper['abstract']}")
        lines.append(f"   URL: {paper['url']}")
        lines.append("")

    return "\n".join(lines)


def _build_html_report(query: str, papers: list[dict]) -> str:
    """Build and return the HTML report string using HtmlReportBuilder."""
    builder = HtmlReportBuilder(
        title=f"PubMed Research Briefing: {query}",
        skill="pubmed-summariser",
    )

    builder.add_header_block(
        title=f"PubMed Research Briefing: {query}",
        subtitle=f"{len(papers)} recent English-language papers",
    )
    builder.add_metadata({
        "Query": query,
        "Results": str(len(papers)),
        "Sorted by": "Date (newest first)",
        "Language filter": "English only",
    })
    builder.add_disclaimer()

    builder.add_section("Papers", level=2)
    for i, paper in enumerate(papers, 1):
        card_html = (
            f"<div style='border:1px solid #e0e0e0;border-radius:8px;"
            f"padding:16px;margin:12px 0;background:#fff;'>"
            f"<h3 style='margin:0 0 8px 0;font-size:1em;'>"
            f"<a href='{paper['url']}' target='_blank'>{i}. {paper['title']}</a></h3>"
            f"<p style='margin:4px 0;color:#616161;font-size:0.9em;'>"
            f"<strong>Authors:</strong> {paper['authors']}</p>"
            f"<p style='margin:4px 0;color:#616161;font-size:0.9em;'>"
            f"<strong>Journal:</strong> {paper['journal']} &middot; {paper['date']}</p>"
        )
        if paper["abstract"]:
            card_html += (
                f"<p style='margin:8px 0 0 0;font-size:0.9em;'>"
                f"{paper['abstract']}</p>"
            )
        card_html += "</div>"
        builder.add_raw_html(card_html)

    builder.add_footer_block(skill="pubmed-summariser", version=SKILL_VERSION)
    return builder.render()


# ---------------------------------------------------------------------------
# Main entry point
# ---------------------------------------------------------------------------

def main(argv: list[str] | None = None) -> None:
    parser = argparse.ArgumentParser(
        description="PubMed Summariser — fetch and summarise recent PubMed papers",
    )
    parser.add_argument("--query", help="Gene name or disease term (e.g. BRCA1, type 2 diabetes)")
    parser.add_argument("--output", required=True, help="Directory to save report.html")
    parser.add_argument("--max-results", type=int, default=10, help="Number of results (default 10, max 50)")
    parser.add_argument("--demo", action="store_true", help="Run demo with BRCA1")
    args = parser.parse_args(argv)

    # Resolve query
    if args.demo:
        if args.query:
            print(f"[demo mode] ignoring --query, using {DEMO_QUERY}")
        query = DEMO_QUERY
    elif args.query:
        query = args.query
    else:
        parser.error("--query is required unless --demo is set")

    max_results = _clamp_max_results(args.max_results)

    # Fetch papers
    try:
        papers = pubmed_api.fetch_papers(query, max_results)
    except Exception as exc:
        print(f"Error fetching papers from PubMed: {exc}", file=sys.stderr)
        sys.exit(1)

    # Terminal output
    summary = _format_terminal_summary(query, papers)
    print(summary)

    # HTML report
    if papers:
        html = _build_html_report(query, papers)
        report_path = write_html_report(args.output, "report.html", html)
        print(f"Report saved to: {report_path}")


if __name__ == "__main__":
    main()
```

- [ ] **Step 3b.2: Run tests — all should pass**

```bash
pytest skills/pubmed-summariser/tests/test_pubmed_summariser.py -v
```

Expected: all tests PASS.

- [ ] **Step 3b.3: Run the full skill test suite to make sure nothing broke**

```bash
pytest skills/pubmed-summariser/tests/ -v
```

Expected: all tests PASS.

- [ ] **Step 3b.4: Commit**

```bash
git add skills/pubmed-summariser/pubmed_summariser.py \
        skills/pubmed-summariser/tests/test_pubmed_summariser.py
git commit -m "feat(pubmed-summariser): add pubmed_summariser.py CLI and renderer with tests"
```

---

## Task 4: Register skill in pytest.ini, CLAUDE.md, and catalog.json

**Files:**
- Modify: `pytest.ini`
- Modify: `CLAUDE.md`
- Modify: `skills/catalog.json`

- [ ] **Step 4.1: Add the test path to pytest.ini**

Open `pytest.ini`. In the `testpaths` block, add this line alongside the others:

```
    skills/pubmed-summariser/tests
```

The full `testpaths` block should look like (new line highlighted):
```ini
testpaths =
    skills/diff-visualizer/tests
    skills/pharmgx-reporter/tests
    skills/equity-scorer/tests
    skills/nutrigx_advisor/tests
    skills/genome-compare/tests
    skills/gwas-prs/tests
    skills/gwas-lookup/tests
    skills/clinpgx/tests
    skills/claw-ancestry-pca/tests
    skills/scrna-orchestrator/tests
    skills/scrna-embedding/tests
    skills/profile-report/tests
    skills/galaxy-bridge/tests
    skills/illumina-bridge/tests
    skills/rnaseq-de/tests
    skills/pubmed-summariser/tests
```

- [ ] **Step 4.2: Verify pytest picks up the new tests**

```bash
pytest skills/pubmed-summariser/tests/ -v --co
```

Expected: lists all test names without errors.

- [ ] **Step 4.3: Update CLAUDE.md (three locations)**

`CLAUDE.md` has three sections that every skill populates. Add to all three:

**A) Routing table** — find the routing table and add this row after the `lit-synthesizer` row:

```markdown
| PubMed search, "summarise PubMed papers about X", "recent papers on gene/disease", research briefing, gene papers, disease papers | `skills/pubmed-summariser/` | Run `pubmed_summariser.py` |
```

**B) CLI Reference block** — find the `## CLI Reference` section and add:

```bash
# PubMed research briefing from gene name or disease term
python skills/pubmed-summariser/pubmed_summariser.py \
  --query <gene_or_disease> --output <report_dir>
python skills/pubmed-summariser/pubmed_summariser.py --demo --output /tmp/pubmed_demo
```

**C) Demo Data table** — find the `## Demo Data` table and add a row:

```markdown
| PubMed Summariser demo (BRCA1, live API) | `--demo` flag | pubmed-summariser |
```

**D) Demo Commands section** — find the `### Demo Commands` code block and add:

```bash
# PubMed Summariser demo
python skills/pubmed-summariser/pubmed_summariser.py --demo --output /tmp/pubmed_demo
```

- [ ] **Step 4.4: Add entry to skills/catalog.json**

Open `skills/catalog.json`. Find the `"skills"` array. Add this entry at the end of the array (before the closing `]`), and increment `"skill_count"` at the top of the file by 1:

```json
{
  "name": "pubmed-summariser",
  "cli_alias": null,
  "description": "Search PubMed for a gene name or disease term and generate a structured research briefing of the top recent English-language papers.",
  "version": "0.1.0",
  "status": "mvp",
  "has_script": true,
  "has_tests": true,
  "has_demo": true,
  "demo_command": "python skills/pubmed-summariser/pubmed_summariser.py --demo --output /tmp/pubmed_demo",
  "dependencies": ["requests"],
  "tags": ["pubmed", "literature", "genomics", "research"],
  "trigger_keywords": ["pubmed", "summarise papers", "research briefing", "papers about", "recent studies", "gene papers"],
  "chaining_partners": ["lit-synthesizer", "gwas-lookup", "gwas-prs"]
}
```

Note: also set `"has_tests": true` (the spec said false — now we have tests). Also set `"skill_count"` at the top of `catalog.json` to **28** (the file currently reads 26, which is already stale — set it to exactly 28, do not just increment).

- [ ] **Step 4.5: Run full test suite to confirm nothing is broken**

```bash
pytest skills/pubmed-summariser/tests/ -v
```

Expected: all tests PASS.

- [ ] **Step 4.6: Commit**

```bash
git add pytest.ini CLAUDE.md skills/catalog.json
git commit -m "feat(pubmed-summariser): register skill in pytest.ini, CLAUDE.md, catalog.json"
```

---

## Task 5: Smoke test end-to-end with live API

This task requires a network connection. It is the final sanity check before declaring the skill done.

- [ ] **Step 5.1: Run the demo**

```bash
python skills/pubmed-summariser/pubmed_summariser.py --demo --output /tmp/pubmed_demo
```

Expected:
- Terminal prints "PubMed Research Briefing: BRCA1"
- Terminal lists 10 papers with titles, authors, journals, dates
- Terminal prints "Report saved to: /tmp/pubmed_demo/report.html"

- [ ] **Step 5.2: Open the HTML report**

```bash
open /tmp/pubmed_demo/report.html
```

Expected: a styled HTML page opens in the browser with a green ClawBio header, paper cards, and the disclaimer.

- [ ] **Step 5.3: Test with a disease query**

```bash
python skills/pubmed-summariser/pubmed_summariser.py \
  --query "type 2 diabetes" --output /tmp/pubmed_disease_demo
```

Expected: 10 papers about type 2 diabetes, report saved.

- [ ] **Step 5.4: Test max-results clamping**

```bash
python skills/pubmed-summariser/pubmed_summariser.py \
  --query BRCA1 --max-results 200 --output /tmp/pubmed_clamp_test
```

Expected: terminal prints `[warning] --max-results capped at 50`, then returns up to 50 results.

- [ ] **Step 5.5: Commit final smoke-test sign-off**

```bash
git add skills/pubmed-summariser/ pytest.ini CLAUDE.md skills/catalog.json
git commit -m "feat(pubmed-summariser): complete MVP — smoke tests pass"
```

---

## Summary

| Task | What it produces |
|---|---|
| Task 1 | `SKILL.md` — AI routing and methodology |
| Task 2 | `pubmed_api.py` + unit tests — NCBI API layer |
| Task 3 | `pubmed_summariser.py` + unit tests — CLI + HTML rendering |
| Task 4 | Registration in `pytest.ini`, `CLAUDE.md`, `catalog.json` |
| Task 5 | Live smoke test against real PubMed API |

After Task 5 completes, the skill is functional and registered. Use `/superpowers:finishing-a-development-branch` to decide on merge/PR strategy.
