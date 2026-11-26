"""EU Clinical Trials Register (EUCTR) client -- secondary European source.

Queries the EUCTR HTML search endpoint as a complement to ClinicalTrials.gov.
Results are normalised to the same dict schema as api._normalise_trial()
so they can be merged and deduplicated transparently.

Note: EUCTR now holds legacy EudraCT trials and their results; ongoing
EU/EEA clinical trials are listed in the Clinical Trials Information System
(CTIS) at https://euclinicaltrials.eu/.

Network, HTTP, and unrecognised-response errors raise EUCTRError instead of
silently returning an empty list so callers surface EU source failures.
"""

from __future__ import annotations

import re
import urllib.error
import urllib.parse
import urllib.request
import xml.etree.ElementTree as ET
from html.parser import HTMLParser

from constants import EUCTR_API


class EUCTRError(RuntimeError):
    """Raised when EUCTR cannot be reached or returns an unrecognised payload."""


# EUCTR uses different status strings than CT.gov. Normalise to CT.gov
# constants so --status filter and FHIR_STATUS mapping work uniformly.
# "Ongoing" means active rather than specifically open to recruitment.
_EUCTR_STATUS_MAP: dict[str, str] = {
    "Ongoing": "ACTIVE_NOT_RECRUITING",
    "Completed": "COMPLETED",
    "Terminated": "TERMINATED",
    "Prematurely Ended": "TERMINATED",
    "Suspended": "SUSPENDED",
    "Temporarily Halted": "SUSPENDED",
    "Withdrawn": "WITHDRAWN",
    "Not Authorised": "WITHDRAWN",
    "Prohibited by CA": "WITHDRAWN",
    "Restarted": "ACTIVE_NOT_RECRUITING",
}

# Deterministic precedence when an EUCTR trial reports multiple country statuses.
_STATUS_PRIORITY = (
    "RECRUITING",
    "ENROLLING_BY_INVITATION",
    "ACTIVE_NOT_RECRUITING",
    "NOT_YET_RECRUITING",
    "SUSPENDED",
    "COMPLETED",
    "TERMINATED",
    "WITHDRAWN",
    "APPROVED_FOR_MARKETING",
    "UNKNOWN",
)

_EUDRACT_RE = re.compile(r"\b(\d{4}-\d{6}-\d{2})\b")
_EMPTY_MARKERS = (
    '<div class="nooutcome"',
    "Query did not match any clinical trials",
)


def _collapse_ws(text: str) -> str:
    return " ".join(text.split())


def _resolve_status(raw_statuses: list[str]) -> str:
    """Resolve one canonical status from possibly multiple country-level statuses."""
    mapped = {_EUCTR_STATUS_MAP.get(s.strip(), "UNKNOWN") for s in raw_statuses if s.strip()}
    if not mapped:
        return "UNKNOWN"
    for candidate in _STATUS_PRIORITY:
        if candidate != "UNKNOWN" and candidate in mapped:
            return candidate
    return "UNKNOWN"


class _EUCTRResultParser(HTMLParser):
    """Parse <table class="result"> blocks from EUCTR HTML search pages."""

    def __init__(self) -> None:
        super().__init__(convert_charrefs=True)
        self.results: list[dict] = []
        self._table_depth = 0
        self._current: dict | None = None
        self._in_td = False
        self._td_depth = 0
        self._label_depth = 0
        self._in_info_tip = 0
        self._status_depth = 0
        self._meddra_depth = 0
        self._current_label: list[str] = []
        self._current_value: list[str] = []
        self._current_status: list[str] = []

    def handle_starttag(self, tag: str, attrs: list[tuple[str, str | None]]) -> None:
        attr_map = dict(attrs)
        classes = set((attr_map.get("class") or "").split())

        if tag == "table":
            if self._current is None and "result" in classes:
                self._current = {"fields": {}, "statuses": []}
                self._table_depth = 1
                return
            if self._current is not None:
                self._table_depth += 1
                if "meddra" in classes and self._meddra_depth == 0:
                    self._meddra_depth = self._table_depth
                return

        if self._current is None:
            return

        if self._meddra_depth > 0:
            # Skip nested Disease/MedDRA summary table cells inside the result row.
            return

        if tag == "td":
            self._td_depth += 1
            if self._td_depth == 1:
                self._in_td = True
                self._current_label = []
                self._current_value = []
        elif tag == "span":
            if self._label_depth > 0:
                self._label_depth += 1
                if "info-tip" in classes or self._in_info_tip > 0:
                    self._in_info_tip += 1
            elif "label" in classes and self._in_td:
                self._label_depth = 1
            elif self._status_depth > 0:
                self._status_depth += 1
            elif "status" in classes:
                self._status_depth = 1
                self._current_status = []

    def handle_endtag(self, tag: str) -> None:
        if self._current is None:
            return

        if tag == "table":
            if self._table_depth == self._meddra_depth:
                self._meddra_depth = 0
            self._table_depth -= 1
            if self._table_depth == 0:
                self.results.append(self._current)
                self._current = None
                self._in_td = False
                self._td_depth = 0
                self._label_depth = 0
                self._in_info_tip = 0
                self._status_depth = 0
            return

        if self._meddra_depth > 0:
            return

        if tag == "span":
            if self._label_depth > 0:
                if self._in_info_tip > 0:
                    self._in_info_tip -= 1
                self._label_depth -= 1
            elif self._status_depth > 0:
                self._status_depth -= 1
                if self._status_depth == 0:
                    raw_status = _collapse_ws("".join(self._current_status)).strip("()")
                    if raw_status:
                        self._current["statuses"].append(raw_status)
                    self._current_status = []
        elif tag == "td" and self._in_td:
            self._td_depth -= 1
            if self._td_depth == 0:
                label = _collapse_ws("".join(self._current_label)).rstrip(":* ").strip()
                value = _collapse_ws("".join(self._current_value))
                if label and value and label not in self._current["fields"]:
                    self._current["fields"][label] = value
                self._in_td = False

    def handle_data(self, data: str) -> None:
        if self._current is None or self._meddra_depth > 0:
            return
        if self._status_depth > 0:
            self._current_status.append(data)
        elif self._label_depth > 0:
            if self._in_info_tip == 0:
                self._current_label.append(data)
        elif self._in_td:
            self._current_value.append(data)


def fetch_euctr(query: str, max_results: int = 10) -> list[dict]:
    """Search the EU Clinical Trials Register and return normalised records.

    Raises EUCTRError on HTTP, network, or unrecognised payload failures.
    """
    url = f"{EUCTR_API}{urllib.parse.quote_plus(query)}&page=1"
    req = urllib.request.Request(
        url,
        headers={
            "Accept": "text/html,application/xml;q=0.9,*/*;q=0.8",
            "User-Agent": "ClawBio-clinical-trial-finder/1.0",
        },
    )

    try:
        with urllib.request.urlopen(req, timeout=15) as resp:
            raw = resp.read().decode("utf-8", errors="replace")
    except urllib.error.HTTPError as exc:
        raise EUCTRError(f"EUCTR HTTP {exc.code}: {exc.reason}") from exc
    except (urllib.error.URLError, OSError) as exc:
        raise EUCTRError(f"EUCTR request failed: {exc}") from exc

    return _parse_euctr_response(raw, max_results)


def _html_record_to_dict(raw_record: dict) -> dict | None:
    fields: dict[str, str] = raw_record["fields"]
    raw_eid = fields.get("EudraCT Number", "")
    match = _EUDRACT_RE.search(raw_eid)
    if not match:
        return None
    eid = match.group(1)

    title = fields.get("Full Title") or f"EU trial {eid}"
    condition = fields.get("Medical condition", "")
    sponsor = fields.get("Sponsor Name", "")

    return {
        "nct_id": f"EUCTR{eid}",
        "title": title,
        "status": _resolve_status(raw_record["statuses"]),
        "phase": "",
        "study_type": "",
        "start_date": fields.get("Start Date", ""),
        "completion_date": "",
        "conditions": [condition] if condition else [],
        "condition_meshes": [],
        "interventions": [],
        "summary": f"Sponsor: {sponsor}" if sponsor else "",
        "source": "euctr",
    }


def _dedupe_and_sort(trials: list[dict], max_results: int) -> list[dict]:
    by_id: dict[str, dict] = {}
    for trial in trials:
        tid = trial.get("nct_id", "")
        if tid and tid not in by_id:
            by_id[tid] = trial
    ordered = [by_id[tid] for tid in sorted(by_id)]
    return ordered[:max_results]


def _parse_euctr_response(raw: str, max_results: int) -> list[dict]:
    """Parse EUCTR search results (HTML search tables or legacy <trials> XML)."""
    if not raw or not raw.strip():
        raise EUCTRError("EUCTR returned an empty response")

    stripped = raw.lstrip()
    if stripped.startswith("<?xml") or stripped.startswith("<trials"):
        try:
            root = ET.fromstring(raw)
        except ET.ParseError as exc:
            raise EUCTRError(f"EUCTR returned malformed XML: {exc}") from exc
        if root.tag == "trials":
            trials = []
            for trial_el in root.findall(".//trial"):
                trial = _euctr_element_to_dict(trial_el)
                if trial.get("nct_id"):
                    trials.append(trial)
            return _dedupe_and_sort(trials, max_results)

    parser = _EUCTRResultParser()
    parser.feed(raw)
    parser.close()

    if parser.results:
        parsed = [
            trial
            for raw_rec in parser.results
            if (trial := _html_record_to_dict(raw_rec)) is not None
        ]
        if not parsed:
            raise EUCTRError("EUCTR result tables contained no valid EudraCT numbers")
        return _dedupe_and_sort(parsed, max_results)

    if any(marker in raw for marker in _EMPTY_MARKERS):
        return []

    raise EUCTRError("Unrecognised EUCTR response format (no result tables or empty-results marker)")


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
