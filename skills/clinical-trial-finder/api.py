"""ClinicalTrials.gov API v2 client -- fetch and normalise trial data.

All network communication with ClinicalTrials.gov lives here.
No other module in this skill makes outbound HTTP calls (except opentargets.py).
"""

import json
import re
import time
import urllib.error
import urllib.parse
import urllib.request
from pathlib import Path

from constants import CT_API, CT_FIELDS, DEFAULT_PAGE_SIZE, MONTH_NAMES

# ---------------------------------------------------------------------------
# Network helpers
# ---------------------------------------------------------------------------

# Maximum retries for transient network errors (HTTP 429, 500, 502, 503, 504).
_MAX_RETRIES = 3
# Base delay in seconds for exponential backoff (1s, 2s, 4s).
_BACKOFF_BASE = 1.0


def _request_with_retry(url: str, timeout: int = 15) -> bytes:
    """Fetch a URL with exponential backoff on transient failures.

    Retries on HTTP 429 (rate limit), 5xx (server errors), and network
    timeouts.  Non-retryable errors (4xx except 429) raise immediately.
    """
    last_exc: Exception | None = None
    for attempt in range(_MAX_RETRIES):
        try:
            with urllib.request.urlopen(url, timeout=timeout) as resp:
                return resp.read()
        except urllib.error.HTTPError as exc:
            if exc.code in (429, 500, 502, 503, 504):
                last_exc = exc
                delay = _BACKOFF_BASE * (2**attempt)
                time.sleep(delay)
                continue
            raise  # 4xx (except 429) -- not retryable
        except (urllib.error.URLError, OSError) as exc:
            last_exc = exc
            delay = _BACKOFF_BASE * (2**attempt)
            time.sleep(delay)
            continue
    raise RuntimeError(
        f"ClinicalTrials.gov unreachable after {_MAX_RETRIES} retries: {last_exc}"
    ) from last_exc


# ---------------------------------------------------------------------------
# Date normalisation helpers
# ---------------------------------------------------------------------------

# Precompiled regexes for date normalisation (called once per trial field).
_RE_DATE_ISO = re.compile(r"^\d{4}-\d{2}(-\d{2})?$")
_RE_DATE_MONTH = re.compile(r"^(\w+)\s+(\d{4})$")


def _to_fhir_date(date_str: str) -> str:
    """Normalize a ClinicalTrials.gov date to FHIR-valid YYYY-MM[-DD].

    CT.gov returns three formats: "2024-01-15", "2024-01", "January 2024".
    FHIR R4 dateTime accepts all three when normalised to ISO.
    Returns the original string unchanged if it cannot be parsed.
    """
    if not date_str:
        return ""
    # Already ISO -- pass through
    if _RE_DATE_ISO.match(date_str):
        return date_str
    # "Month YYYY" -- convert to YYYY-MM
    m = _RE_DATE_MONTH.match(date_str)
    if m:
        month_num = MONTH_NAMES.get(m.group(1))
        if month_num:
            return f"{m.group(2)}-{month_num}"
    return date_str


# ---------------------------------------------------------------------------
# Input parsing
# ---------------------------------------------------------------------------


def parse_input(input_path: Path) -> dict:
    """Read a query file with one search term per line.

    Lines starting with # are comments.  All non-empty, non-comment lines
    are joined into a single query string for the CT.gov API.
    """
    lines = input_path.read_text().splitlines()
    terms = [
        line.strip()
        for line in lines
        if line.strip() and not line.lstrip().startswith("#")
    ]
    if not terms:
        raise ValueError(
            f"No search terms in {input_path}. Add at least one non-comment line."
        )
    return {"query": " ".join(terms), "terms": terms}


# ---------------------------------------------------------------------------
# ClinicalTrials.gov API
# ---------------------------------------------------------------------------


def fetch_trials(
    query: str,
    max_results: int = DEFAULT_PAGE_SIZE,
    country: str | None = None,
) -> list[dict]:
    """Query ClinicalTrials.gov API v2 and return normalised trial records.

    Uses query.cond (condition field, MeSH-indexed) for better recall than
    free-text search.  Supports multi-page pagination via the nextPageToken
    cursor returned by the API.  Each page request uses exponential backoff
    via ``_request_with_retry``.

    When *country* is provided, it is passed as query.locn to restrict
    results to trials in that country (ISO 3166-1 name or code).
    """
    # CT.gov API v2 caps pageSize at 1000.  We request in pages of up to
    # 1000 and accumulate until we reach max_results or run out of data.
    page_size = min(max_results, 1000)
    params: dict[str, str | int] = {
        "query.cond": query,
        "pageSize": page_size,
        "fields": CT_FIELDS,
        "format": "json",
    }
    if country:
        params["query.locn"] = country

    all_trials: list[dict] = []
    page_token: str | None = None

    while len(all_trials) < max_results:
        if page_token:
            params["pageToken"] = page_token
        url = f"{CT_API}?{urllib.parse.urlencode(params)}"

        try:
            raw = _request_with_retry(url)
            data = json.loads(raw.decode())
        except json.JSONDecodeError as exc:
            raise RuntimeError(
                f"ClinicalTrials.gov returned malformed JSON: {exc}"
            ) from exc

        studies = data.get("studies", [])
        if not studies:
            break  # No more results

        all_trials.extend(_normalise_trial(s) for s in studies)

        # CT.gov API v2 returns nextPageToken when more pages exist.
        page_token = data.get("nextPageToken")
        if not page_token:
            break  # Last page

    return all_trials[:max_results]


def _normalise_trial(study: dict) -> dict:
    """Flatten a nested CT.gov study object into a simple dict.

    The CT.gov API returns deeply nested JSON (protocolSection ->
    identificationModule -> nctId, etc).  We flatten to a single level
    so downstream functions (report, FHIR, chart) don't need to know
    the API schema.
    """
    proto = study.get("protocolSection", {})
    id_mod = proto.get("identificationModule", {})
    status_mod = proto.get("statusModule", {})
    design_mod = proto.get("designModule", {})
    desc_mod = proto.get("descriptionModule", {})
    cond_mod = proto.get("conditionsModule", {})
    interv_mod = proto.get("armsInterventionsModule", {})
    derived = study.get("derivedSection", {})

    summary = desc_mod.get("briefSummary", "")
    interventions = [
        arm["name"] for arm in interv_mod.get("interventions", []) if arm.get("name")
    ]

    # MeSH codes from CT.gov's derivedSection.conditionBrowseModule --
    # NLM computes these during study indexing, so we get authoritative
    # MeSH IDs without a separate NLM API call.
    condition_meshes = derived.get("conditionBrowseModule", {}).get("meshes", [])

    return {
        "nct_id": id_mod.get("nctId", ""),
        "title": id_mod.get("briefTitle", ""),
        "status": status_mod.get("overallStatus", "UNKNOWN"),
        # CT.gov returns phases as a list, e.g. ["PHASE1", "PHASE2"] -> "PHASE1 / PHASE2"
        "phase": " / ".join(design_mod.get("phases", [])),
        "study_type": design_mod.get("studyType", ""),
        "start_date": _to_fhir_date(
            status_mod.get("startDateStruct", {}).get("date", "")
        ),
        "completion_date": _to_fhir_date(
            status_mod.get("completionDateStruct", {}).get("date", "")
        ),
        "conditions": cond_mod.get("conditions", []),
        "condition_meshes": condition_meshes,
        "interventions": interventions,
        # CT.gov summaries can exceed 2000 chars; truncate for the report while
        # summary.json preserves the full text via the raw API response.
        "summary": summary[:300] + "..." if len(summary) > 300 else summary,
    }
