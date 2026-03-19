"""E2E test — validates FHIR R4 datatype compliance of the generated bundle.

Hits the live ClinicalTrials.gov API (free, no auth), generates a FHIR bundle,
and checks every resource against FHIR R4 ResearchStudy datatype rules.

Run:  python3 -m pytest skills/clinical-trial-finder/tests/test_fhir_e2e.py -v
"""

import importlib.util
import json
import os
import re
import sys
from pathlib import Path

import pytest

# Skip E2E tests in CI to avoid flaky failures from live API calls.
# Run locally with: RUN_E2E=1 python3 -m pytest tests/test_fhir_e2e.py -v
pytestmark = pytest.mark.skipif(
    os.getenv("RUN_E2E") != "1" and os.getenv("CI") == "true",
    reason="E2E tests skipped in CI (set RUN_E2E=1 to force)",
)

SKILL_DIR = Path(__file__).resolve().parent.parent


def _load(name: str):
    spec = importlib.util.spec_from_file_location(name, SKILL_DIR / f"{name}.py")
    mod = importlib.util.module_from_spec(spec)
    sys.modules[name] = mod
    spec.loader.exec_module(mod)
    return mod


const = _load("constants")
api = _load("api")
wr = _load("writers")
ctf = _load("clinical_trial_finder")

# FHIR R4 dateTime regex — accepts YYYY, YYYY-MM, YYYY-MM-DD, or full datetime
# Source: https://hl7.org/fhir/R4/datatypes.html#dateTime
FHIR_DATETIME_RE = re.compile(
    r"^([0-9]{4})"  # YYYY
    r"(-(0[1-9]|1[0-2])"  # -MM
    r"(-(0[1-9]|[12][0-9]|3[01])"  # -DD
    r"(T([01][0-9]|2[0-3]):[0-5][0-9]"  # Thh:mm
    r"(:[0-5][0-9](\.[0-9]+)?)?"  # :ss.sss
    r"(Z|[+-]((0[0-9]|1[0-3]):[0-5][0-9]|14:00))?"  # timezone
    r")?)?)?$"
)

# FHIR R4 instant regex — full datetime with timezone required
FHIR_INSTANT_RE = re.compile(
    r"^[0-9]{4}-(0[1-9]|1[0-2])-(0[1-9]|[12][0-9]|3[01])"
    r"T([01][0-9]|2[0-3]):[0-5][0-9]:[0-5][0-9](\.[0-9]+)?"
    r"(Z|[+-]((0[0-9]|1[0-3]):[0-5][0-9]|14:00))$"
)

# Valid ResearchStudy.status codes (FHIR R4)
VALID_STATUSES = {
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

# Valid ResearchStudy.phase codes (FHIR R4)
VALID_PHASES = {
    "n-a",
    "early-phase-1",
    "phase-1",
    "phase-1-phase-2",
    "phase-2",
    "phase-2-phase-3",
    "phase-3",
    "phase-4",
}

PHASE_SYSTEM = "http://terminology.hl7.org/CodeSystem/research-study-phase"


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


@pytest.fixture(scope="module")
def fhir_bundle(tmp_path_factory):
    """Fetch real trials and generate FHIR bundle (runs once per module)."""
    trials = api.fetch_trials("BRCA1 breast cancer", max_results=10)
    if not trials:
        pytest.skip("ClinicalTrials.gov returned 0 results — network issue?")
    out = tmp_path_factory.mktemp("fhir_e2e")
    path = wr.write_fhir_bundle(trials, out)
    return json.loads(path.read_text())


@pytest.fixture(scope="module")
def resources(fhir_bundle):
    return [e["resource"] for e in fhir_bundle["entry"]]


# ---------------------------------------------------------------------------
# Bundle-level validation
# ---------------------------------------------------------------------------


class TestBundleStructure:
    def test_resource_type(self, fhir_bundle):
        assert fhir_bundle["resourceType"] == "Bundle"

    def test_bundle_type(self, fhir_bundle):
        assert fhir_bundle["type"] == "searchset"

    def test_total_matches_entries(self, fhir_bundle):
        assert fhir_bundle["total"] == len(fhir_bundle["entry"])

    def test_timestamp_is_fhir_instant(self, fhir_bundle):
        ts = fhir_bundle["timestamp"]
        assert FHIR_INSTANT_RE.match(ts), (
            f"Bundle.timestamp '{ts}' is not a valid FHIR instant"
        )

    def test_entries_have_fullurl(self, fhir_bundle):
        for entry in fhir_bundle["entry"]:
            assert "fullUrl" in entry, "Bundle.entry must have fullUrl"
            assert entry["fullUrl"].startswith("http"), (
                f"fullUrl should be an absolute URI: {entry['fullUrl']}"
            )


# ---------------------------------------------------------------------------
# ResearchStudy resource validation
# ---------------------------------------------------------------------------


class TestResearchStudyResource:
    def test_resource_type(self, resources):
        for r in resources:
            assert r["resourceType"] == "ResearchStudy"

    def test_has_id(self, resources):
        for r in resources:
            assert r.get("id"), "ResearchStudy must have an id"

    def test_has_identifier(self, resources):
        for r in resources:
            ids = r.get("identifier", [])
            assert len(ids) >= 1, f"{r['id']}: must have at least one identifier"
            for ident in ids:
                assert "system" in ident, f"{r['id']}: Identifier missing system"
                assert "value" in ident, f"{r['id']}: Identifier missing value"

    def test_has_title(self, resources):
        for r in resources:
            assert isinstance(r.get("title"), str) and r["title"], (
                f"{r['id']}: title must be a non-empty string"
            )


# ---------------------------------------------------------------------------
# FHIR status code validation
# ---------------------------------------------------------------------------


class TestFHIRStatus:
    def test_status_is_valid_code(self, resources):
        for r in resources:
            status = r.get("status")
            assert status in VALID_STATUSES, (
                f"{r['id']}: status '{status}' is not a valid FHIR "
                f"research-study-status code. Valid: {VALID_STATUSES}"
            )

    def test_all_ct_statuses_mapped(self):
        """Every ClinicalTrials.gov status must map to a valid FHIR code."""
        for ct_status, fhir_code in const.FHIR_STATUS.items():
            assert fhir_code in VALID_STATUSES, (
                f"FHIR_STATUS['{ct_status}'] = '{fhir_code}' is not valid"
            )


# ---------------------------------------------------------------------------
# FHIR phase CodeableConcept validation
# ---------------------------------------------------------------------------


class TestFHIRPhase:
    def test_phase_is_codeable_concept(self, resources):
        for r in resources:
            phase = r.get("phase")
            assert isinstance(phase, dict), (
                f"{r['id']}: phase must be a CodeableConcept (dict)"
            )
            codings = phase.get("coding", [])
            assert len(codings) >= 1, (
                f"{r['id']}: phase.coding must have at least one Coding"
            )

    def test_phase_coding_has_system(self, resources):
        for r in resources:
            for coding in r["phase"]["coding"]:
                assert coding.get("system") == PHASE_SYSTEM, (
                    f"{r['id']}: phase coding system must be {PHASE_SYSTEM}"
                )

    def test_phase_code_is_valid(self, resources):
        for r in resources:
            for coding in r["phase"]["coding"]:
                assert coding["code"] in VALID_PHASES, (
                    f"{r['id']}: phase code '{coding['code']}' not in "
                    f"FHIR value set. Valid: {VALID_PHASES}"
                )

    def test_phase_coding_has_display_when_known(self, resources):
        """FHIR best practice: Coding should include display when phase is known.

        Phases with no display (e.g. OBSERVATIONAL studies with no phase)
        correctly omit the field -- FHIR R4 forbids empty strings.
        """
        missing = []
        for r in resources:
            for coding in r["phase"]["coding"]:
                code = coding.get("code", "")
                if code != "n-a" and "display" not in coding:
                    missing.append(r["id"])
        if missing:
            pytest.fail(
                f"phase.coding missing 'display' in {len(missing)} resources: "
                f"{missing[:5]}... — FHIR recommends display for human readability"
            )


# ---------------------------------------------------------------------------
# FHIR Period (dateTime) validation
# ---------------------------------------------------------------------------


class TestFHIRPeriod:
    def test_period_dates_are_valid_fhir_datetime(self, resources):
        """Period.start and Period.end must be valid FHIR dateTime."""
        invalid = []
        for r in resources:
            period = r.get("period", {})
            for field in ("start", "end"):
                val = period.get(field)
                if val and not FHIR_DATETIME_RE.match(val):
                    invalid.append(f"{r['id']}.period.{field}='{val}'")
        if invalid:
            pytest.fail(
                f"Invalid FHIR dateTime values: {invalid}. "
                f"FHIR dateTime must be YYYY, YYYY-MM, YYYY-MM-DD, or full ISO datetime"
            )

    def test_period_end_after_start(self, resources):
        """If both start and end are present, end >= start."""
        for r in resources:
            period = r.get("period", {})
            start = period.get("start", "")
            end = period.get("end", "")
            if start and end:
                assert end >= start, (
                    f"{r['id']}: period.end '{end}' is before period.start '{start}'"
                )


# ---------------------------------------------------------------------------
# FHIR condition CodeableConcept validation
# ---------------------------------------------------------------------------


class TestFHIRCondition:
    def test_condition_is_codeable_concept_list(self, resources):
        for r in resources:
            conds = r.get("condition", [])
            assert isinstance(conds, list), (
                f"{r['id']}: condition must be a list of CodeableConcept"
            )
            for c in conds:
                assert isinstance(c, dict), (
                    f"{r['id']}: each condition must be a CodeableConcept (dict)"
                )
                has_text = "text" in c
                has_coding = "coding" in c
                assert has_text or has_coding, (
                    f"{r['id']}: CodeableConcept must have text or coding"
                )

    def test_condition_has_coding(self, resources):
        """MeSH coding from CT.gov derivedSection — most conditions should be coded."""
        total, coded = 0, 0
        for r in resources:
            for c in r.get("condition", []):
                total += 1
                if "coding" in c:
                    # Verify system is MeSH
                    assert any(
                        "mesh" in coding.get("system", "").lower()
                        for coding in c["coding"]
                    ), f"Unexpected coding system in {r['id']}: {c['coding']}"
                    coded += 1
        if total > 0:
            pct = coded / total * 100
            assert pct >= 50, (
                f"Only {coded}/{total} ({pct:.0f}%) conditions have MeSH coding. "
                f"Expected ≥50% — check derivedSection parsing."
            )


# ---------------------------------------------------------------------------
# Data completeness — available data should not be silently dropped
# ---------------------------------------------------------------------------


class TestDataCompleteness:
    def test_study_type_mapped_to_category(self, resources):
        """study_type (INTERVENTIONAL/OBSERVATIONAL) should be in FHIR category."""
        missing = [r["id"] for r in resources if "category" not in r]
        if missing:
            pytest.fail(
                f"{len(missing)} resources missing 'category' element. "
                f"ClinicalTrials.gov study_type should map to "
                f"ResearchStudy.category (CodeableConcept). IDs: {missing[:5]}"
            )

    def test_interventions_mapped_to_focus(self, resources):
        """Interventions should appear as ResearchStudy.focus."""
        missing = [r["id"] for r in resources if "focus" not in r]
        # Not all trials have interventions, so only flag if we have entries
        # but none have focus
        total_with_focus = len(resources) - len(missing)
        if total_with_focus == 0 and len(resources) > 0:
            pytest.fail(
                "No resources have 'focus' element. "
                "Intervention names from ClinicalTrials.gov should map to "
                "ResearchStudy.focus (CodeableConcept[])"
            )

    def test_identifier_has_use(self, resources):
        """FHIR best practice: Identifier should include 'use' field."""
        missing = []
        for r in resources:
            for ident in r.get("identifier", []):
                if "use" not in ident:
                    missing.append(r["id"])
        if missing:
            pytest.fail(
                f"{len(missing)} identifiers missing 'use' field. "
                f"Recommended: 'official' for NCT IDs"
            )
