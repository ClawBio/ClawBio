"""
Unit tests for pubmed_api.py — no network calls required.

All tests operate on pre-built XML strings or fixture data.
"""

import sys
from pathlib import Path
import xml.etree.ElementTree as ET

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))


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


# ── _date_sort_key ──────────────────────────────────────────────────────────


def test_date_sort_key_full_iso():
    from pubmed_api import _date_sort_key
    assert _date_sort_key("2026-03-19") == (2026, 3, 19)


def test_date_sort_key_numeric_month():
    from pubmed_api import _date_sort_key
    assert _date_sort_key("2026-03") == (2026, 3, 0)


def test_date_sort_key_month_abbr():
    from pubmed_api import _date_sort_key
    assert _date_sort_key("2026-Mar") == (2026, 3, 0)
    assert _date_sort_key("2026-Apr") == (2026, 4, 0)


def test_date_sort_key_year_only():
    from pubmed_api import _date_sort_key
    assert _date_sort_key("2025") == (2025, 0, 0)


def test_date_sort_key_unknown():
    from pubmed_api import _date_sort_key
    assert _date_sort_key("") == (0, 0, 0)
    assert _date_sort_key("not-a-date") == (0, 0, 0)


def test_date_sort_key_ordering():
    """Newer dates must produce a larger tuple."""
    from pubmed_api import _date_sort_key
    dates = ["2025-11-24", "2026-Mar", "2026-03-19", "2026-Apr", "2025"]
    sorted_dates = sorted(dates, key=_date_sort_key, reverse=True)
    # Apr > Mar 19 > Mar (month-only) > Nov 2025 > 2025 (year-only)
    assert sorted_dates[0] == "2026-Apr"
    assert sorted_dates[1] == "2026-03-19"
    assert sorted_dates[2] == "2026-Mar"
    assert sorted_dates[3] == "2025-11-24"
    assert sorted_dates[4] == "2025"


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
