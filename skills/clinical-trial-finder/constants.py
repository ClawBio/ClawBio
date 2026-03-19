"""Shared constants, mappings, and labels for clinical-trial-finder.

All static data lives here so that api.py, writers.py, and the CLI can
import what they need without circular dependencies.
"""

from pathlib import Path

# ---------------------------------------------------------------------------
# Paths
# ---------------------------------------------------------------------------

SKILL_DIR = Path(__file__).resolve().parent
DEMO_DATA = SKILL_DIR / "demo_input.txt"

# ---------------------------------------------------------------------------
# ClinicalTrials.gov API
# ---------------------------------------------------------------------------

# ClinicalTrials.gov API v2 -- the authoritative US registry mandated by
# FDAAA 801.  We use query.cond (MeSH-indexed) rather than query.term
# (free-text) for better recall on condition searches.
CT_API = "https://clinicaltrials.gov/api/v2/studies"

# Fields requested from the API.  We ask for the minimum set needed to
# build the report + FHIR bundle.  ConditionMeshId and ConditionMeshTerm
# pull MeSH codes from CT.gov's own derivedSection (no extra NLM call).
CT_FIELDS = ",".join(
    [
        "NCTId",
        "BriefTitle",
        "OverallStatus",
        "Phase",
        "StartDate",
        "CompletionDate",
        "StudyType",
        "BriefSummary",
        "Condition",
        "InterventionName",
        "ConditionMeshId",
        "ConditionMeshTerm",
    ]
)

# 20 trials balances coverage with actionability -- a clinician reviewing
# more than 20 without eligibility pre-screening is unlikely to act on any.
DEFAULT_PAGE_SIZE = 20

# ---------------------------------------------------------------------------
# Safety disclaimer (required in every report by ClawBio policy)
# ---------------------------------------------------------------------------

DISCLAIMER = (
    "*ClawBio is a research and educational tool. It is not a medical device "
    "and does not provide clinical diagnoses. Consult a healthcare professional "
    "before making any medical decisions.*"
)

# ---------------------------------------------------------------------------
# Status constants
# ---------------------------------------------------------------------------

# Every status ClinicalTrials.gov can return, used for --status CLI choices.
ALL_STATUSES = [
    "RECRUITING",
    "ACTIVE_NOT_RECRUITING",
    "NOT_YET_RECRUITING",
    "ENROLLING_BY_INVITATION",
    "COMPLETED",
    "TERMINATED",
    "WITHDRAWN",
    "SUSPENDED",
    "APPROVED_FOR_MARKETING",
    "UNKNOWN",
]

# Text labels for trial status in the markdown report.
# Deliberate choice over emojis: renders correctly in any terminal,
# screen reader, or plain-text viewer.
STATUS_LABEL: dict[str, str] = {
    "RECRUITING": "[RECRUITING]",
    "ACTIVE_NOT_RECRUITING": "[ACTIVE]",
    "COMPLETED": "[COMPLETED]",
    "NOT_YET_RECRUITING": "[PENDING]",
    "ENROLLING_BY_INVITATION": "[BY INVITATION]",
    "TERMINATED": "[TERMINATED]",
    "WITHDRAWN": "[WITHDRAWN]",
    "SUSPENDED": "[SUSPENDED]",
    "APPROVED_FOR_MARKETING": "[APPROVED]",
    "UNKNOWN": "[UNKNOWN]",
}

# ---------------------------------------------------------------------------
# FHIR R4 mappings
#
# Source: hl7.org/fhir/R4/valueset-research-study-status.html
# We use R4 (not R5) because the ONC 21st Century Cures Act mandates R4
# for certified EHR systems -- Epic, Cerner, Oracle Health all use R4.
# ---------------------------------------------------------------------------

FHIR_STATUS: dict[str, str] = {
    "RECRUITING": "active",
    "ACTIVE_NOT_RECRUITING": "active",
    "NOT_YET_RECRUITING": "approved",
    "ENROLLING_BY_INVITATION": "active",
    "COMPLETED": "completed",
    "TERMINATED": "administratively-completed",
    "WITHDRAWN": "withdrawn",
    "SUSPENDED": "temporarily-closed-to-accrual",
    "APPROVED_FOR_MARKETING": "completed",
    # FHIR R4 schema does not include "unknown" -- map to "in-review" as the
    # closest semantically neutral status for trials CT.gov can't confirm.
    "UNKNOWN": "in-review",
}

# Maps CT.gov phase strings to FHIR R4 research-study-phase codes.
FHIR_PHASE: dict[str, str] = {
    "EARLY_PHASE1": "early-phase-1",
    "PHASE1": "phase-1",
    "PHASE2": "phase-2",
    "PHASE3": "phase-3",
    "PHASE4": "phase-4",
    "PHASE1 / PHASE2": "phase-1-phase-2",
    "PHASE2 / PHASE3": "phase-2-phase-3",
    "NA": "n-a",
    "N/A": "n-a",  # CT.gov API v2 sometimes returns "N/A" with slash
}

# Human-readable display for FHIR phase CodeableConcept.
FHIR_PHASE_DISPLAY: dict[str, str] = {
    "EARLY_PHASE1": "Early Phase 1",
    "PHASE1": "Phase 1",
    "PHASE2": "Phase 2",
    "PHASE3": "Phase 3",
    "PHASE4": "Phase 4",
    "PHASE1 / PHASE2": "Phase 1/Phase 2",
    "PHASE2 / PHASE3": "Phase 2/Phase 3",
    "NA": "N/A",
    "N/A": "N/A",
}

# CT.gov returns dates as "January 2024" -- we need YYYY-MM for FHIR.
MONTH_NAMES: dict[str, str] = {
    "January": "01",
    "February": "02",
    "March": "03",
    "April": "04",
    "May": "05",
    "June": "06",
    "July": "07",
    "August": "08",
    "September": "09",
    "October": "10",
    "November": "11",
    "December": "12",
}

# ---------------------------------------------------------------------------
# Visualisation constants
# ---------------------------------------------------------------------------

# Human-readable phase labels for the chart X axis.
PHASE_DISPLAY: dict[str, str] = {
    "EARLY_PHASE1": "Early Phase 1",
    "PHASE1": "Phase 1",
    "PHASE2": "Phase 2",
    "PHASE3": "Phase 3",
    "PHASE4": "Phase 4",
    "PHASE1 / PHASE2": "Phase 1/2",
    "PHASE2 / PHASE3": "Phase 2/3",
    "NA": "N/A",
    "N/A": "N/A",
}

# Colour palette for stacked bars -- ordered from "positive" to "negative"
# statuses so the legend reads naturally.
STATUS_COLOR: dict[str, str] = {
    "RECRUITING": "#2ecc71",
    "NOT_YET_RECRUITING": "#3498db",
    "ENROLLING_BY_INVITATION": "#27ae60",
    "ACTIVE_NOT_RECRUITING": "#f1c40f",
    "SUSPENDED": "#e67e22",
    "TERMINATED": "#e74c3c",
    "WITHDRAWN": "#7f8c8d",
    "COMPLETED": "#95a5a6",
    "APPROVED_FOR_MARKETING": "#8e44ad",
    "UNKNOWN": "#bdc3c7",
}

# Fixed order so phases always appear left-to-right by progression.
PHASE_ORDER = [
    "Early Phase 1",
    "Phase 1",
    "Phase 1/2",
    "Phase 2",
    "Phase 2/3",
    "Phase 3",
    "Phase 4",
    "N/A",
    "Unknown",
]

# ---------------------------------------------------------------------------
# CSV output columns
# ---------------------------------------------------------------------------

CSV_COLUMNS = [
    "nct_id",
    "title",
    "status",
    "phase",
    "study_type",
    "start_date",
    "completion_date",
    "conditions",
    "interventions",
]

# ---------------------------------------------------------------------------
# FHIR R4 validation -- expected value sets
# ---------------------------------------------------------------------------

FHIR_VALID_STATUSES = {
    "active",
    "administratively-completed",
    "approved",
    "closed-to-accrual",
    "closed-to-accrual-and-intervention",
    "completed",
    "disapproved",
    "in-review",
    "temporarily-closed-to-accrual",
    "temporarily-closed-to-accrual-and-intervention",
    "withdrawn",
}

FHIR_VALID_PHASES = {
    "n-a",
    "early-phase-1",
    "phase-1",
    "phase-1-phase-2",
    "phase-2",
    "phase-2-phase-3",
    "phase-3",
    "phase-4",
}

# ---------------------------------------------------------------------------
# EU Clinical Trials Register (EUCTR)
# ---------------------------------------------------------------------------

EUCTR_API = "https://www.clinicaltrialsregister.eu/ctr-search/rest/search?query="
