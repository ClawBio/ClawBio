"""EU Clinical Trials Register (EUCTR) client -- secondary European source.

Queries the EUCTR REST API as a complement to ClinicalTrials.gov.
Results are normalised to the same dict schema as api._normalise_trial()
so they can be merged and deduplicated transparently.

The EUCTR API is best-effort: it returns XML, has no versioning guarantees,
and may be unavailable.  All failures degrade gracefully to an empty list.
"""

import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET

from constants import EUCTR_API

# EUCTR uses different status strings than CT.gov. Normalise to CT.gov
# constants so --status filter and FHIR_STATUS mapping work uniformly.
_EUCTR_STATUS_MAP: dict[str, str] = {
    "Ongoing": "RECRUITING",
    "Completed": "COMPLETED",
    "Terminated": "TERMINATED",
    "Prematurely Ended": "TERMINATED",
    "Suspended": "SUSPENDED",
    "Withdrawn": "WITHDRAWN",
    "Not Authorised": "WITHDRAWN",
    "Prohibited by CA": "WITHDRAWN",
    "Restarted": "RECRUITING",
}


def fetch_euctr(query: str, max_results: int = 10) -> list[dict]:
    """Search the EU Clinical Trials Register and return normalised records.

    Returns an empty list on any error (network, parse, API change).
    """
    url = f"{EUCTR_API}{urllib.parse.quote(query)}&page=1&mode=current_page"

    try:
        req = urllib.request.Request(url, headers={"Accept": "application/xml"})
        with urllib.request.urlopen(req, timeout=10) as resp:
            raw = resp.read().decode("utf-8", errors="replace")
    except (urllib.error.URLError, OSError):
        return []

    return _parse_euctr_response(raw, max_results)


def _parse_euctr_response(raw: str, max_results: int) -> list[dict]:
    """Parse EUCTR search results (XML or fallback HTML scraping).

    The EUCTR API is inconsistent -- sometimes returns XML, sometimes HTML.
    We try XML first, then fall back to a minimal HTML parse.
    """
    trials: list[dict] = []

    # Try XML parse first
    try:
        root = ET.fromstring(raw)
        for trial_el in root.findall(".//trial")[:max_results]:
            trial = _euctr_element_to_dict(trial_el)
            if trial.get("nct_id"):  # only include if we got an ID
                trials.append(trial)
        if trials:
            return trials
    except ET.ParseError:
        pass

    # Fallback: regex extraction from raw HTML.  We use regex instead of an
    # HTML parser because EUCTR responses are inconsistent and often malformed;
    # a proper parser would fail more often than a targeted regex.
    import re

    # EudraCT IDs follow format YYYY-NNNNNN-CC (year, 6-digit seq, 2-digit check)
    eudract_pattern = re.compile(r"(\d{4}-\d{6}-\d{2})")
    title_pattern = re.compile(r"<td[^>]*class=\"full\"[^>]*>([^<]+)</td>")

    eudract_ids = eudract_pattern.findall(raw)
    titles = title_pattern.findall(raw)

    for i, eid in enumerate(eudract_ids[:max_results]):
        trials.append(
            {
                "nct_id": f"EUCTR{eid}",
                "title": titles[i].strip() if i < len(titles) else f"EU trial {eid}",
                "status": "UNKNOWN",
                "phase": "",
                "study_type": "",
                "start_date": "",
                "completion_date": "",
                "conditions": [],
                "condition_meshes": [],
                "interventions": [],
                "summary": "",
                "source": "euctr",
            }
        )

    return trials


def _euctr_element_to_dict(el: ET.Element) -> dict:
    """Convert a single EUCTR XML <trial> element to a normalised dict.

    Output schema matches api._normalise_trial() so CT.gov and EUCTR trials
    can be merged in a single list.  Fields not available in EUCTR (e.g.
    condition_meshes, completion_date) default to empty values.
    """

    def _text(tag: str) -> str:
        child = el.find(tag)
        return (child.text or "").strip() if child is not None else ""

    return {
        "nct_id": _text("eudract_number") or _text("id"),
        "title": _text("full_title") or _text("title"),
        "status": _EUCTR_STATUS_MAP.get(_text("trial_status"), "UNKNOWN"),
        "phase": _text("trial_phase"),
        "study_type": _text("trial_type"),
        "start_date": _text("date_on"),
        "completion_date": "",
        "conditions": [c for c in [_text("medical_condition")] if c],
        "condition_meshes": [],
        "interventions": [i for i in [_text("trade_name")] if i],
        "summary": _text("primary_end_point"),
        "source": "euctr",
    }
