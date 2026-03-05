"""rec_shortener — Extract phenotype-specific CPIC recommendations from structured tables.

Parses the HTML tables in ClinPGx guideline textMarkdown to find the exact
recommendation for the patient's phenotype. No NLP or LLM needed — the data
is structured in the CPIC guideline tables.

Usage:
    from clawbio.common.rec_shortener import extract_phenotype_rec

    # rec, strength = extract_phenotype_rec(guideline_html, "Intermediate Metabolizer", "CYP2D6")
    # → ("Use label recommended dosing. If no response, consider non-tramadol opioid", "Moderate")
"""

import re
from html.parser import HTMLParser

__all__ = [
    "extract_phenotype_rec",
    "extract_all_recs_from_guidelines",
    "extract_all_source_recs",
    "shorten_rec",
]


# ── HTML table parser ─────────────────────────────────────────────────────────

class _TableParser(HTMLParser):
    """Parse HTML tables into list-of-lists."""

    def __init__(self):
        super().__init__()
        self.tables: list[list[list[str]]] = []
        self._in_table = False
        self._in_row = False
        self._in_cell = False
        self._current_table: list[list[str]] = []
        self._current_row: list[str] = []
        self._current_cell = ""

    def handle_starttag(self, tag, attrs):
        if tag == "table":
            self._in_table = True
            self._current_table = []
        elif tag == "tr":
            self._in_row = True
            self._current_row = []
        elif tag in ("td", "th"):
            self._in_cell = True
            self._current_cell = ""

    def handle_endtag(self, tag):
        if tag in ("td", "th") and self._in_cell:
            self._current_row.append(self._current_cell.strip())
            self._in_cell = False
        elif tag == "tr" and self._in_row:
            if self._current_row:
                self._current_table.append(self._current_row)
            self._in_row = False
        elif tag == "table" and self._in_table:
            self.tables.append(self._current_table)
            self._in_table = False

    def handle_data(self, data):
        if self._in_cell:
            self._current_cell += data


def _parse_tables(html: str) -> list[list[list[str]]]:
    parser = _TableParser()
    parser.feed(html)
    return parser.tables


# ── Phenotype matching ────────────────────────────────────────────────────────

# Map patient phenotype strings to canonical phenotype keywords for matching
_PHENO_KEYWORDS = {
    "ultrarapid metabolizer": "ultrarapid",
    "rapid metabolizer": "rapid",
    "normal metabolizer": "normal",
    "intermediate metabolizer": "intermediate",
    "poor metabolizer": "poor",
    "extensive metabolizer": "extensive",
    # CYP3A5
    "cyp3a5 expressor": "expressor",
    "cyp3a5 non-expressor": "non-expressor",
    # VKORC1
    "high warfarin sensitivity": "high",
    "moderate warfarin sensitivity": "moderate",
    "low warfarin sensitivity": "low",
    # SLCO1B1
    "normal function": "normal",
    "decreased function": "decreased",
    "poor function": "poor",
    # Generic
    "normal (inferred)": "normal",
    "indeterminate": "indeterminate",
}


def _match_phenotype(row_pheno: str, patient_phenotype: str) -> bool:
    """Check if a table row phenotype matches the patient's phenotype."""
    row_lower = row_pheno.lower().strip()
    patient_lower = patient_phenotype.lower().strip()

    # Exact match
    if patient_lower == row_lower:
        return True

    # Word-boundary substring match (avoids "rapid" matching "ultrarapid")
    if re.search(r'\b' + re.escape(patient_lower) + r'\b', row_lower):
        return True
    if re.search(r'\b' + re.escape(row_lower) + r'\b', patient_lower):
        return True

    # Map patient phenotype to keyword and check (word-boundary match)
    for full, keyword in _PHENO_KEYWORDS.items():
        if full in patient_lower:
            if re.search(r'\b' + re.escape(keyword) + r'\b', row_lower):
                return True

    return False


def _find_rec_column(header: list[str]) -> tuple[int, int]:
    """Find the Recommendation and Classification/Strength column indices."""
    rec_idx = -1
    strength_idx = -1
    for i, h in enumerate(header):
        hl = h.lower()
        if "recommendation" in hl and "classification" not in hl:
            rec_idx = i
        elif "classification" in hl or "strength" in hl:
            strength_idx = i
    return rec_idx, strength_idx


# ── Strength normalization ────────────────────────────────────────────────────

# CPIC classification of recommendations (what they actually mean):
#   Strong   = High confidence, change prescribing
#   Moderate = Moderate confidence, consider change
#   Optional = Weak evidence, prescriber discretion
_STRENGTH_MAP = {
    "strong": "Strong",
    "moderate": "Moderate",
    "optional": "Optional",
}


def _clean_strength(raw: str) -> str:
    """Normalize CPIC strength values (strip footnotes, fix casing, drop 'no recommendation')."""
    if not raw:
        return ""
    stripped = raw.strip()
    lower = stripped.lower()
    if "no recommendation" in lower or "n/a" in lower:
        return ""
    # Try direct match first (avoids stripping letters from clean values)
    if lower in _STRENGTH_MAP:
        return _STRENGTH_MAP[lower]
    # Strip trailing footnote letters (a-g) that CPIC appends, e.g. "Stronge" → "Strong"
    cleaned = re.sub(r'[a-g]+$', '', stripped)
    lower_cleaned = cleaned.lower().strip()
    return _STRENGTH_MAP.get(lower_cleaned, cleaned.title())


# ── Public API ────────────────────────────────────────────────────────────────

def extract_phenotype_rec(
    guideline_html: str,
    patient_phenotype: str,
    gene: str = "",
) -> tuple[str, str]:
    """Extract the recommendation for a specific phenotype from a CPIC guideline table.

    Args:
        guideline_html: The textMarkdown.html content from a ClinPGx guideline.
        patient_phenotype: The patient's phenotype (e.g. "Intermediate Metabolizer").
        gene: Gene symbol (e.g. "CYP2D6") for context.

    Returns:
        (recommendation_text, strength) tuple.
        Returns ("", "") if no match found.
    """
    if not guideline_html or not patient_phenotype:
        return ("", "")

    tables = _parse_tables(guideline_html)

    for table in tables:
        if len(table) < 2:
            continue

        header = table[0]
        rec_idx, strength_idx = _find_rec_column(header)

        # If no explicit "Recommendation" column, try heuristic:
        # first column = phenotype, last-ish columns = rec + strength
        if rec_idx == -1:
            # Check if any header contains "recommendation" substring
            for i, h in enumerate(header):
                if "rec" in h.lower():
                    rec_idx = i
                    break
            if rec_idx == -1:
                continue

        # Search rows for matching phenotype
        for row in table[1:]:
            if not row:
                continue
            row_pheno = row[0]
            if _match_phenotype(row_pheno, patient_phenotype):
                rec = row[rec_idx] if rec_idx < len(row) else ""
                strength = row[strength_idx] if strength_idx != -1 and strength_idx < len(row) else ""
                if rec and rec.lower() != "n/a":
                    return (rec.strip(), _clean_strength(strength))

    return ("", "")


def extract_all_recs_from_guidelines(
    guidelines: list[dict],
    drug_name: str,
    patient_phenotype: str,
    gene: str = "",
) -> tuple[str, str, str]:
    """Extract phenotype-specific recommendation from a list of ClinPGx guidelines.

    Searches DPWG first (most concrete dose/alternative advice), then CPIC.

    Args:
        guidelines: List of guideline dicts from ClinPGxClient.get_guidelines().
        drug_name: Drug name to match.
        patient_phenotype: Patient's phenotype string.
        gene: Gene symbol.

    Returns:
        (recommendation, strength, source) tuple.
        Returns ("", "", "") if no match found.
    """
    drug_lower = drug_name.lower()

    # Priority: DPWG first (most concrete advice), then CPIC, then others
    source_order = ["DPWG", "CPIC", "CPNDS", "RNPGx"]
    guidelines_by_source = {s: [] for s in source_order}
    other_guidelines = []

    for g in guidelines:
        source = g.get("source", "")
        name = g.get("name", "").lower()
        # Match by drug name in guideline name
        if drug_lower not in name and gene.lower() not in name:
            continue
        if source in guidelines_by_source:
            guidelines_by_source[source].append(g)
        else:
            other_guidelines.append(g)

    # Try each source in priority order
    for source in source_order:
        for g in guidelines_by_source[source]:
            html = g.get("textMarkdown", {})
            if isinstance(html, dict):
                html = html.get("html", "")
            if not html:
                continue
            rec, strength = extract_phenotype_rec(html, patient_phenotype, gene)
            if rec:
                return (rec, strength, source)

    # Try others
    for g in other_guidelines:
        html = g.get("textMarkdown", {})
        if isinstance(html, dict):
            html = html.get("html", "")
        if not html:
            continue
        rec, strength = extract_phenotype_rec(html, patient_phenotype, gene)
        if rec:
            return (rec, strength, g.get("source", ""))

    return ("", "", "")


def extract_all_source_recs(
    guidelines: list[dict],
    drug_name: str,
    patient_phenotype: str,
    gene: str = "",
) -> list[dict]:
    """Extract phenotype-specific recommendations from ALL sources.

    Returns a list of {source, rec, strength} dicts — one per source that
    has a matching recommendation. Sources checked: DPWG, CPIC, CPNDS, RNPGx.
    """
    drug_lower = drug_name.lower()

    source_order = ["DPWG", "CPIC", "CPNDS", "RNPGx"]
    guidelines_by_source: dict[str, list[dict]] = {s: [] for s in source_order}
    other_guidelines: list[dict] = []

    for g in guidelines:
        source = g.get("source", "")
        name = g.get("name", "").lower()
        if drug_lower not in name and gene.lower() not in name:
            continue
        if source in guidelines_by_source:
            guidelines_by_source[source].append(g)
        else:
            other_guidelines.append(g)

    results = []
    seen_sources = set()

    for source in source_order:
        for g in guidelines_by_source[source]:
            if source in seen_sources:
                break
            html = g.get("textMarkdown", {})
            if isinstance(html, dict):
                html = html.get("html", "")
            if not html:
                continue
            rec, strength = extract_phenotype_rec(html, patient_phenotype, gene)
            if rec:
                results.append({"source": source, "rec": rec, "strength": strength})
                seen_sources.add(source)
                break

    for g in other_guidelines:
        source = g.get("source", "")
        if source in seen_sources:
            continue
        html = g.get("textMarkdown", {})
        if isinstance(html, dict):
            html = html.get("html", "")
        if not html:
            continue
        rec, strength = extract_phenotype_rec(html, patient_phenotype, gene)
        if rec:
            results.append({"source": source, "rec": rec, "strength": strength})
            seen_sources.add(source)

    return results


def shorten_rec(text: str) -> str:
    """Shorten a recommendation to its first sentence, capped at 120 chars.

    Handles CPIC footnote markers (e.g. "dose.g Utilize") by splitting on
    period followed by whitespace+uppercase, which marks a real sentence boundary.
    """
    if not text:
        return ""
    # Split at period followed by optional footnote chars then whitespace+uppercase
    # This handles "dose.g Utilize" → first sentence = "dose."
    m = re.search(r'\.([a-g,]*)[\s]+(?=[A-Z])', text)
    if m:
        first = text[:m.start() + 1]  # up to and including the period
    else:
        first = text.rstrip(".") + "."
    # Cap length
    if len(first) > 120:
        cut = first[:120].rfind(", ")
        if cut == -1:
            cut = first[:120].rfind("; ")
        if cut > 40:
            first = first[:cut] + "."
        else:
            first = first[:117] + "..."
    return first
