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

_MONTH_ABBR = {
    "Jan": 1, "Feb": 2, "Mar": 3, "Apr": 4, "May": 5, "Jun": 6,
    "Jul": 7, "Aug": 8, "Sep": 9, "Oct": 10, "Nov": 11, "Dec": 12,
}


def _date_sort_key(date_str: str) -> tuple:
    """
    Convert a date string to a (year, month, day) tuple for comparison.

    Handles formats: YYYY-MM-DD, YYYY-MM, YYYY-Mon, YYYY.
    Returns (0, 0, 0) for unrecognised strings so they sort last.
    """
    # YYYY-MM-DD
    m = re.match(r'^(\d{4})-(\d{2})-(\d{2})$', date_str)
    if m:
        return (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    # YYYY-MM (numeric)
    m = re.match(r'^(\d{4})-(\d{2})$', date_str)
    if m:
        return (int(m.group(1)), int(m.group(2)), 0)
    # YYYY-Mon (3-letter abbreviation)
    m = re.match(r'^(\d{4})-([A-Za-z]{3})$', date_str)
    if m:
        month = _MONTH_ABBR.get(m.group(2).capitalize(), 0)
        return (int(m.group(1)), month, 0)
    # YYYY only
    m = re.match(r'^(\d{4})$', date_str)
    if m:
        return (int(m.group(1)), 0, 0)
    return (0, 0, 0)


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
    title = "".join(title_el.itertext()).strip() if title_el is not None else ""

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

    # Prefer electronic publication date over journal issue date (avoids
    # future eCollection dates on online-first articles)
    if art is not None:
        epublish = art.find("ArticleDate[@DateType='Electronic']")
        if epublish is not None:
            ep_year = epublish.findtext("Year", "").strip()
            ep_month = epublish.findtext("Month", "").strip()
            ep_day = epublish.findtext("Day", "").strip()
            if ep_year:
                date = f"{ep_year}-{ep_month}-{ep_day}" if ep_month and ep_day else (f"{ep_year}-{ep_month}" if ep_month else ep_year)

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
        max_results: Number of papers to return (default 10)

    Returns:
        List of dicts with keys: title, authors, journal, date, abstract, pmid, url
        Empty list if no results.

    Raises:
        requests.RequestException: on network failure or HTTP error
        xml.etree.ElementTree.ParseError: if the NCBI response XML is malformed
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

    papers.sort(key=lambda p: _date_sort_key(p["date"]), reverse=True)
    return papers
