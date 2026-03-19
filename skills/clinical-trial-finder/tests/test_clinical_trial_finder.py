"""Tests for clinical-trial-finder -- no network calls."""

import argparse
import importlib.util
import json
import sys
import urllib.error
from pathlib import Path
from unittest.mock import MagicMock, patch

import pytest

# Load modules from skill directory without installing them as packages.
# Order matters: constants first (no deps), then api (needs constants),
# then writers (needs constants), then the CLI entry point, then opentargets.
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
ot = _load("opentargets")


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------

# Three trials covering the key scenarios:
#   A: RECRUITING, PHASE2, has MeSH coding, has interventions (happy path)
#   B: COMPLETED, PHASE3, no MeSH (text fallback), no interventions
#   C: TERMINATED, no phase, OBSERVATIONAL, all fields empty (edge case)
MOCK_TRIALS = [
    {
        "nct_id": "NCT00000001",
        "title": "Synthetic Trial A",
        "status": "RECRUITING",
        "phase": "PHASE2",
        "study_type": "INTERVENTIONAL",
        "start_date": "2024-01",
        "completion_date": "2026-12",
        "conditions": ["Breast Cancer"],
        "condition_meshes": [{"id": "D001943", "term": "Breast Neoplasms"}],
        "interventions": ["Drug X"],
        "summary": "Synthetic trial for testing.",
    },
    {
        "nct_id": "NCT00000002",
        "title": "Synthetic Trial B",
        "status": "COMPLETED",
        "phase": "PHASE3",
        "study_type": "INTERVENTIONAL",
        "start_date": "2020-01",
        "completion_date": "2023-12",
        "conditions": ["Breast Cancer", "BRCA1 Mutation"],
        "condition_meshes": [],  # no MeSH — falls back to text
        "interventions": [],
        "summary": "",
    },
    {
        "nct_id": "NCT00000003",
        "title": "Synthetic Trial C",
        "status": "TERMINATED",
        "phase": "",
        "study_type": "OBSERVATIONAL",
        "start_date": "",
        "completion_date": "",
        "conditions": [],
        "condition_meshes": [],
        "interventions": [],
        "summary": "",
    },
]

MOCK_QUERY = {"query": "BRCA1 breast cancer", "terms": ["BRCA1 breast cancer"]}


# ---------------------------------------------------------------------------
# parse_input
# ---------------------------------------------------------------------------


def test_parse_input_basic(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("# comment\nBRCA1\nbreast cancer\n")
    result = api.parse_input(f)
    assert result["query"] == "BRCA1 breast cancer"
    assert result["terms"] == ["BRCA1", "breast cancer"]


def test_parse_input_empty_raises(tmp_path):
    f = tmp_path / "empty.txt"
    f.write_text("# only comments\n\n")
    with pytest.raises(ValueError, match="No search terms"):
        api.parse_input(f)


def test_parse_input_single_term(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("EGFR\n")
    result = api.parse_input(f)
    assert result["query"] == "EGFR"


# ---------------------------------------------------------------------------
# _normalise_trial
# ---------------------------------------------------------------------------


def test_normalise_trial_truncates_summary():
    long_summary = "x" * 400
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {"overallStatus": "RECRUITING"},
            "designModule": {"phases": ["PHASE1"], "studyType": "INTERVENTIONAL"},
            "descriptionModule": {"briefSummary": long_summary},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = api._normalise_trial(raw)
    assert len(trial["summary"]) <= 303  # 300 chars + "..."
    assert trial["summary"].endswith("...")


def test_normalise_trial_missing_modules():
    """Graceful handling of empty study object."""
    trial = api._normalise_trial({"protocolSection": {}})
    assert trial["nct_id"] == ""
    assert trial["status"] == "UNKNOWN"
    assert trial["conditions"] == []


# ---------------------------------------------------------------------------
# fetch_trials
# ---------------------------------------------------------------------------


def _make_ct_response(trials_raw: list[dict]) -> bytes:
    return json.dumps({"studies": trials_raw}).encode()


def test_fetch_trials_network_error():
    with patch("urllib.request.urlopen", side_effect=urllib.error.URLError("timeout")):
        with pytest.raises(RuntimeError, match="ClinicalTrials.gov unreachable"):
            api.fetch_trials("BRCA1")


def test_fetch_trials_malformed_json():
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = b"not-json{"
    with patch("urllib.request.urlopen", return_value=mock_resp):
        with pytest.raises(RuntimeError, match="malformed JSON"):
            api.fetch_trials("BRCA1")


def test_fetch_trials_empty_results():
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = json.dumps({"studies": []}).encode()
    with patch("urllib.request.urlopen", return_value=mock_resp):
        result = api.fetch_trials("nothing-matches-this-query-xyz")
    assert result == []


# ---------------------------------------------------------------------------
# count_recruiting
# ---------------------------------------------------------------------------


def test_count_recruiting():
    assert wr.count_recruiting(MOCK_TRIALS) == 1


def test_count_recruiting_none():
    assert wr.count_recruiting([MOCK_TRIALS[1], MOCK_TRIALS[2]]) == 0


# ---------------------------------------------------------------------------
# write_report
# ---------------------------------------------------------------------------


def test_write_report_creates_file(tmp_path):
    path = wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    assert path.exists()


def test_write_report_contains_nct_ids(tmp_path):
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    for t in MOCK_TRIALS:
        assert t["nct_id"] in content


def test_write_report_contains_disclaimer(tmp_path):
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "research and educational tool" in content


def test_write_report_gene_context(tmp_path):
    gene_context = {
        "symbol": "BRCA1",
        "name": "BRCA1 DNA repair associated",
        "ensembl": "ENSG00000012048",
        "diseases": ["breast cancer", "ovarian cancer"],
        "min_score": 0.6,
    }
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path, gene_context=gene_context)
    content = (tmp_path / "report.md").read_text()
    assert "BRCA1" in content
    assert "OpenTargets" in content


def test_write_report_terminated_trial_shown(tmp_path):
    """TERMINATED trials must not be hidden."""
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "NCT00000003" in content


def test_write_report_uses_text_labels_not_emojis(tmp_path):
    """Report should use plain-text status labels, not emojis."""
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "[RECRUITING]" in content
    assert "[TERMINATED]" in content


# ---------------------------------------------------------------------------
# write_summary
# ---------------------------------------------------------------------------


def test_write_summary(tmp_path):
    path = wr.write_summary(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    data = json.loads(path.read_text())
    assert data["total"] == 3
    assert data["recruiting"] == 1
    assert len(data["trials"]) == 3


# ---------------------------------------------------------------------------
# write_fhir_bundle
# ---------------------------------------------------------------------------


def test_write_fhir_bundle_structure(tmp_path):
    path = wr.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads(path.read_text())
    assert bundle["resourceType"] == "Bundle"
    assert bundle["type"] == "searchset"
    assert bundle["total"] == 3


def test_fhir_status_mapping(tmp_path):
    wr.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    by_id = {e["resource"]["id"]: e["resource"] for e in bundle["entry"]}
    assert by_id["NCT00000001"]["status"] == "active"
    assert by_id["NCT00000002"]["status"] == "completed"
    assert by_id["NCT00000003"]["status"] == "administratively-completed"


def test_fhir_phase_mapping(tmp_path):
    wr.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    by_id = {e["resource"]["id"]: e["resource"] for e in bundle["entry"]}
    assert by_id["NCT00000001"]["phase"]["coding"][0]["code"] == "phase-2"


def test_fhir_all_statuses_covered():
    known = [
        "RECRUITING",
        "ACTIVE_NOT_RECRUITING",
        "NOT_YET_RECRUITING",
        "COMPLETED",
        "TERMINATED",
        "WITHDRAWN",
        "SUSPENDED",
        "UNKNOWN",
    ]
    for s in known:
        assert s in const.FHIR_STATUS, f"Missing FHIR mapping for: {s}"


def test_fhir_trial_without_conditions(tmp_path):
    """Trial with no conditions should not have 'condition' key."""
    wr.write_fhir_bundle([MOCK_TRIALS[2]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "condition" not in resource


def test_fhir_identifier_use_official(tmp_path):
    wr.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    identifier = bundle["entry"][0]["resource"]["identifier"][0]
    assert identifier["use"] == "official"
    assert identifier["value"] == "NCT00000001"


def test_fhir_phase_display(tmp_path):
    wr.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    coding = bundle["entry"][0]["resource"]["phase"]["coding"][0]
    assert coding["code"] == "phase-2"
    assert coding["display"] == "Phase 2"


def test_fhir_category_from_study_type(tmp_path):
    wr.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "category" in resource
    assert resource["category"][0]["coding"][0]["code"] == "INTERVENTIONAL"


def test_fhir_focus_from_interventions(tmp_path):
    wr.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "focus" in resource
    assert resource["focus"][0]["text"] == "Drug X"


def test_fhir_no_focus_when_no_interventions(tmp_path):
    wr.write_fhir_bundle([MOCK_TRIALS[1]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "focus" not in resource


def test_fhir_condition_uses_mesh_when_available(tmp_path):
    """MeSH-coded conditions preferred over free-text when derivedSection provides them."""
    wr.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    cond = bundle["entry"][0]["resource"]["condition"][0]
    assert "coding" in cond
    assert cond["coding"][0]["system"] == "http://id.nlm.nih.gov/mesh/"
    assert cond["coding"][0]["code"] == "D001943"
    assert cond["text"] == "Breast Neoplasms"


def test_fhir_condition_falls_back_to_text_without_mesh(tmp_path):
    """Trial without MeSH data falls back to text-only CodeableConcept."""
    wr.write_fhir_bundle([MOCK_TRIALS[1]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    cond = bundle["entry"][0]["resource"]["condition"][0]
    assert "coding" not in cond
    assert cond["text"] == "Breast Cancer"


def test_normalise_trial_extracts_condition_meshes():
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {"overallStatus": "RECRUITING"},
            "designModule": {},
            "descriptionModule": {},
            "conditionsModule": {"conditions": ["Breast Cancer"]},
            "armsInterventionsModule": {},
        },
        "derivedSection": {
            "conditionBrowseModule": {
                "meshes": [{"id": "D001943", "term": "Breast Neoplasms"}]
            }
        },
    }
    trial = api._normalise_trial(raw)
    assert trial["condition_meshes"] == [{"id": "D001943", "term": "Breast Neoplasms"}]


def test_normalise_trial_empty_meshes_when_no_derived():
    raw = {
        "protocolSection": {
            "identificationModule": {},
            "statusModule": {},
            "designModule": {},
            "descriptionModule": {},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = api._normalise_trial(raw)
    assert trial["condition_meshes"] == []


# ---------------------------------------------------------------------------
# _to_fhir_date
# ---------------------------------------------------------------------------


def test_to_fhir_date_passthrough_iso():
    assert api._to_fhir_date("2024-01") == "2024-01"
    assert api._to_fhir_date("2024-01-15") == "2024-01-15"
    assert api._to_fhir_date("") == ""


def test_to_fhir_date_normalizes_month_year():
    assert api._to_fhir_date("January 2024") == "2024-01"
    assert api._to_fhir_date("December 2026") == "2026-12"
    assert api._to_fhir_date("March 2020") == "2020-03"


def test_to_fhir_date_unknown_format_passthrough():
    assert api._to_fhir_date("Q1 2024") == "Q1 2024"


# ---------------------------------------------------------------------------
# opentargets module
# ---------------------------------------------------------------------------


def _make_ot_response(data: dict) -> MagicMock:
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = json.dumps({"data": data}).encode()
    return mock_resp


def test_ot_resolve_gene():
    resp = _make_ot_response(
        {
            "search": {
                "hits": [
                    {
                        "id": "ENSG00000012048",
                        "object": {
                            "approvedSymbol": "BRCA1",
                            "approvedName": "BRCA1 DNA repair associated",
                        },
                    }
                ]
            }
        }
    )
    with patch("urllib.request.urlopen", return_value=resp):
        ensembl_id, name = ot.resolve_gene("BRCA1")
    assert ensembl_id == "ENSG00000012048"
    assert "BRCA1" in name


def test_ot_resolve_gene_not_found():
    resp = _make_ot_response({"search": {"hits": []}})
    with patch("urllib.request.urlopen", return_value=resp):
        with pytest.raises(ValueError, match="not found"):
            ot.resolve_gene("FAKEGENE999")


def test_ot_get_diseases_filters_by_score():
    resp = _make_ot_response(
        {
            "target": {
                "associatedDiseases": {
                    "rows": [
                        {
                            "disease": {"id": "MONDO_1", "name": "Breast Cancer"},
                            "score": 0.85,
                        },
                        {
                            "disease": {"id": "MONDO_2", "name": "Low Score Disease"},
                            "score": 0.3,
                        },
                        {
                            "disease": {"id": "MONDO_3", "name": "Ovarian Cancer"},
                            "score": 0.80,
                        },
                    ]
                }
            }
        }
    )
    with patch("urllib.request.urlopen", return_value=resp):
        diseases = ot.get_diseases("ENSG00000012048", min_score=0.6)
    names = [d.name for d in diseases]
    assert "Breast Cancer" in names
    assert "Ovarian Cancer" in names
    assert "Low Score Disease" not in names


def test_ot_network_error():
    with patch("urllib.request.urlopen", side_effect=urllib.error.URLError("refused")):
        with pytest.raises(RuntimeError, match="unreachable"):
            ot.resolve_gene("BRCA1")


def test_ot_graphql_error():
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = json.dumps(
        {"errors": [{"message": "Field not found"}]}
    ).encode()
    with patch("urllib.request.urlopen", return_value=mock_resp):
        with pytest.raises(RuntimeError, match="GraphQL error"):
            ot.resolve_gene("BRCA1")


# ---------------------------------------------------------------------------
# --status filter
# ---------------------------------------------------------------------------


def test_status_filter_recruiting():
    filtered = [t for t in MOCK_TRIALS if t["status"] == "RECRUITING"]
    assert len(filtered) == 1
    assert filtered[0]["nct_id"] == "NCT00000001"


def test_status_filter_none_returns_all():
    filtered = [t for t in MOCK_TRIALS if True]  # no filter
    assert len(filtered) == len(MOCK_TRIALS)


# ---------------------------------------------------------------------------
# write_checksums
# ---------------------------------------------------------------------------


def test_write_checksums(tmp_path):
    (tmp_path / "report.md").write_text("# Test")
    (tmp_path / "summary.json").write_text("{}")
    path = wr.write_checksums(tmp_path)
    assert path.exists()
    content = path.read_text()
    lines = [line for line in content.strip().split("\n") if line]
    assert len(lines) == 2
    for line in lines:
        digest, name = line.split("  ", 1)
        assert len(digest) == 64
        assert name in ("report.md", "summary.json")


def test_write_checksums_skips_missing(tmp_path):
    (tmp_path / "report.md").write_text("# Test")
    path = wr.write_checksums(tmp_path)
    content = path.read_text().strip()
    assert "summary.json" not in content
    assert "report.md" in content


def test_write_checksums_includes_figures(tmp_path):
    (tmp_path / "report.md").write_text("x")
    fig_dir = tmp_path / "figures"
    fig_dir.mkdir()
    (fig_dir / "phase_distribution.png").write_bytes(b"\x89PNG")
    path = wr.write_checksums(tmp_path)
    assert "figures/phase_distribution.png" in path.read_text()


# ---------------------------------------------------------------------------
# write_commands
# ---------------------------------------------------------------------------


def test_write_commands_demo(tmp_path):
    args = argparse.Namespace(
        demo=True,
        input=None,
        query=None,
        gene=None,
        status=None,
        max_results=20,
        fhir=False,
        ot_min_score=0.6,
        ot_max_diseases=5,
    )
    path = wr.write_commands(args, tmp_path)
    assert path.exists()
    content = path.read_text()
    assert "--demo" in content
    assert f"--output '{tmp_path}'" in content


def test_write_commands_query_with_status(tmp_path):
    args = argparse.Namespace(
        demo=False,
        input=None,
        query="lung cancer",
        gene=None,
        status="RECRUITING",
        max_results=20,
        fhir=True,
        ot_min_score=0.6,
        ot_max_diseases=5,
    )
    path = wr.write_commands(args, tmp_path)
    content = path.read_text()
    assert "--query 'lung cancer'" in content
    assert "--status RECRUITING" in content
    assert "--fhir" in content


def test_write_commands_non_default_ot_params(tmp_path):
    args = argparse.Namespace(
        demo=False,
        input=None,
        query=None,
        gene="BRCA1",
        status=None,
        max_results=50,
        fhir=False,
        ot_min_score=0.3,
        ot_max_diseases=10,
    )
    path = wr.write_commands(args, tmp_path)
    content = path.read_text()
    assert "--gene 'BRCA1'" in content
    assert "--max-results 50" in content
    assert "--ot-min-score 0.3" in content
    assert "--ot-max-diseases 10" in content


# ---------------------------------------------------------------------------
# write_phase_chart — figures/ subdir
# ---------------------------------------------------------------------------


def test_chart_in_figures_subdir(tmp_path):
    chart = wr.write_phase_chart(MOCK_TRIALS, tmp_path)
    if chart is None:
        pytest.skip("matplotlib not installed")
    assert chart.parent.name == "figures"
    assert chart.name == "phase_distribution.png"
    assert chart.exists()


# ---------------------------------------------------------------------------
# constants.py -- mapping consistency
# ---------------------------------------------------------------------------


def test_all_statuses_have_labels():
    """Every CT.gov status should have a text label for the report."""
    for s in const.ALL_STATUSES:
        assert s in const.STATUS_LABEL, f"STATUS_LABEL missing: {s}"


def test_all_statuses_have_fhir_mapping():
    """Every CT.gov status should map to a valid FHIR R4 code."""
    for s in const.ALL_STATUSES:
        assert s in const.FHIR_STATUS, f"FHIR_STATUS missing: {s}"


def test_fhir_status_values_are_valid():
    """Every mapped FHIR status code must be in the R4 value set."""
    for ct, fhir in const.FHIR_STATUS.items():
        assert fhir in const.FHIR_VALID_STATUSES, (
            f"FHIR_STATUS['{ct}'] = '{fhir}' not in FHIR_VALID_STATUSES"
        )


def test_fhir_phase_values_are_valid():
    """Every mapped FHIR phase code must be in the R4 value set."""
    for ct, fhir in const.FHIR_PHASE.items():
        assert fhir in const.FHIR_VALID_PHASES, (
            f"FHIR_PHASE['{ct}'] = '{fhir}' not in FHIR_VALID_PHASES"
        )


def test_csv_columns_match_normalised_trial_keys():
    """CSV_COLUMNS should be a subset of the keys _normalise_trial returns."""
    trial = api._normalise_trial({"protocolSection": {}})
    for col in const.CSV_COLUMNS:
        assert col in trial, f"CSV_COLUMNS has '{col}' but _normalise_trial doesn't"


def test_status_color_covers_all_statuses():
    """Chart colours should be defined for every possible status."""
    for s in const.ALL_STATUSES:
        assert s in const.STATUS_COLOR, f"STATUS_COLOR missing: {s}"


# ---------------------------------------------------------------------------
# api.py -- fetch_trials with country parameter
# ---------------------------------------------------------------------------


def test_fetch_trials_country_param():
    """When country is provided, query.locn should appear in the URL."""
    captured_url = {}

    def mock_urlopen(url, **kwargs):
        captured_url["url"] = url
        mock_resp = MagicMock()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_resp.read.return_value = json.dumps({"studies": []}).encode()
        return mock_resp

    with patch("urllib.request.urlopen", side_effect=mock_urlopen):
        api.fetch_trials("BRCA1", country="Spain")
    assert "query.locn=Spain" in captured_url["url"]


def test_fetch_trials_no_country_param():
    """When country is None, query.locn should NOT appear in the URL."""
    captured_url = {}

    def mock_urlopen(url, **kwargs):
        captured_url["url"] = url
        mock_resp = MagicMock()
        mock_resp.__enter__ = lambda s: s
        mock_resp.__exit__ = MagicMock(return_value=False)
        mock_resp.read.return_value = json.dumps({"studies": []}).encode()
        return mock_resp

    with patch("urllib.request.urlopen", side_effect=mock_urlopen):
        api.fetch_trials("BRCA1")
    assert "query.locn" not in captured_url["url"]


# ---------------------------------------------------------------------------
# api.py -- _normalise_trial edge cases
# ---------------------------------------------------------------------------


def test_normalise_trial_summary_exactly_300():
    """Summary of exactly 300 chars should NOT be truncated."""
    summary_300 = "x" * 300
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {"overallStatus": "RECRUITING"},
            "designModule": {},
            "descriptionModule": {"briefSummary": summary_300},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = api._normalise_trial(raw)
    assert trial["summary"] == summary_300
    assert not trial["summary"].endswith("...")


def test_normalise_trial_multiple_phases():
    """Multiple phases should be joined with ' / '."""
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {"overallStatus": "RECRUITING"},
            "designModule": {"phases": ["PHASE1", "PHASE2"]},
            "descriptionModule": {},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = api._normalise_trial(raw)
    assert trial["phase"] == "PHASE1 / PHASE2"


def test_normalise_trial_multiple_interventions():
    """All interventions should be extracted."""
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {},
            "designModule": {},
            "descriptionModule": {},
            "conditionsModule": {},
            "armsInterventionsModule": {
                "interventions": [
                    {"name": "Drug A"},
                    {"name": "Drug B"},
                    {"name": ""},  # empty name should be skipped
                ]
            },
        }
    }
    trial = api._normalise_trial(raw)
    assert trial["interventions"] == ["Drug A", "Drug B"]


# ---------------------------------------------------------------------------
# api.py -- parse_input edge cases
# ---------------------------------------------------------------------------


def test_parse_input_strips_whitespace(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("  BRCA1  \n  breast cancer  \n")
    result = api.parse_input(f)
    assert result["terms"] == ["BRCA1", "breast cancer"]


def test_parse_input_skips_blank_lines(tmp_path):
    f = tmp_path / "q.txt"
    f.write_text("BRCA1\n\n\nbreast cancer\n\n")
    result = api.parse_input(f)
    assert result["terms"] == ["BRCA1", "breast cancer"]


# ---------------------------------------------------------------------------
# writers.py -- write_report edge cases
# ---------------------------------------------------------------------------


def test_write_report_empty_trials(tmp_path):
    """Report with 0 trials should still be valid markdown with disclaimer."""
    path = wr.write_report(MOCK_QUERY, [], tmp_path)
    content = path.read_text()
    assert "Total trials | 0" in content
    assert "Recruiting now | 0" in content
    assert "research and educational tool" in content


def test_write_report_has_reproducibility_section(tmp_path):
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "## Reproducibility" in content
    assert "commands.sh" in content
    assert "checksums.sha256" in content


def test_write_report_nct_links_are_valid_urls(tmp_path):
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    for t in MOCK_TRIALS:
        assert f"https://clinicaltrials.gov/study/{t['nct_id']}" in content


# ---------------------------------------------------------------------------
# writers.py -- write_summary edge cases
# ---------------------------------------------------------------------------


def test_write_summary_has_required_keys(tmp_path):
    path = wr.write_summary(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    data = json.loads(path.read_text())
    for key in ("query", "timestamp", "source", "total", "recruiting", "trials"):
        assert key in data, f"summary.json missing key: {key}"


def test_write_summary_empty_trials(tmp_path):
    path = wr.write_summary(MOCK_QUERY, [], tmp_path)
    data = json.loads(path.read_text())
    assert data["total"] == 0
    assert data["recruiting"] == 0
    assert data["trials"] == []


# ---------------------------------------------------------------------------
# writers.py -- write_fhir_bundle edge cases
# ---------------------------------------------------------------------------


def test_write_fhir_bundle_empty_trials(tmp_path):
    path = wr.write_fhir_bundle([], tmp_path)
    bundle = json.loads(path.read_text())
    assert bundle["total"] == 0
    assert bundle["entry"] == []


def test_fhir_bundle_timestamp_has_timezone(tmp_path):
    """FHIR instant requires timezone in the timestamp."""
    path = wr.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads(path.read_text())
    assert "+00:00" in bundle["timestamp"]


def test_fhir_trial_without_study_type_has_no_category(tmp_path):
    """Trial with empty study_type should not have 'category' key."""
    wr.write_fhir_bundle([MOCK_TRIALS[2]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    # MOCK_TRIALS[2] has study_type "OBSERVATIONAL" so let's test with empty
    trial_no_type = {**MOCK_TRIALS[2], "study_type": ""}
    wr.write_fhir_bundle([trial_no_type], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "category" not in resource


def test_fhir_trial_without_dates_has_no_period(tmp_path):
    """Trial with empty dates should not have 'period' key."""
    wr.write_fhir_bundle([MOCK_TRIALS[2]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "period" not in resource


def test_fhir_fullurl_is_clinicaltrials_link(tmp_path):
    path = wr.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads(path.read_text())
    entry = bundle["entry"][0]
    assert entry["fullUrl"] == "https://clinicaltrials.gov/study/NCT00000001"


def test_fhir_trial_without_summary_has_no_description(tmp_path):
    """Trial with empty summary should not have 'description' key."""
    wr.write_fhir_bundle([MOCK_TRIALS[2]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "description" not in resource


# ---------------------------------------------------------------------------
# writers.py -- _sha256
# ---------------------------------------------------------------------------


def test_sha256_produces_correct_hash(tmp_path):
    """Verify _sha256 matches hashlib directly."""
    import hashlib

    f = tmp_path / "test.txt"
    f.write_text("hello world")
    expected = hashlib.sha256(b"hello world").hexdigest()
    assert wr._sha256(f) == expected


def test_sha256_empty_file(tmp_path):
    import hashlib

    f = tmp_path / "empty.txt"
    f.write_bytes(b"")
    expected = hashlib.sha256(b"").hexdigest()
    assert wr._sha256(f) == expected


# ---------------------------------------------------------------------------
# writers.py -- write_checksums hash correctness
# ---------------------------------------------------------------------------


def test_write_checksums_hashes_are_correct(tmp_path):
    """Verify that the checksums in the file match actual SHA-256 of the outputs."""
    import hashlib

    (tmp_path / "report.md").write_text("# Report")
    (tmp_path / "summary.json").write_text('{"test": true}')
    wr.write_checksums(tmp_path)
    for line in (tmp_path / "checksums.sha256").read_text().strip().split("\n"):
        digest, rel = line.split("  ", 1)
        actual = hashlib.sha256((tmp_path / rel).read_bytes()).hexdigest()
        assert digest == actual, f"Checksum mismatch for {rel}"


# ---------------------------------------------------------------------------
# writers.py -- write_commands edge cases
# ---------------------------------------------------------------------------


def test_write_commands_with_input_path(tmp_path):
    args = argparse.Namespace(
        demo=False,
        input=Path("/data/queries.txt"),
        query=None,
        gene=None,
        status=None,
        max_results=20,
        fhir=False,
        ot_min_score=0.6,
        ot_max_diseases=5,
    )
    path = wr.write_commands(args, tmp_path)
    content = path.read_text()
    assert "--input '/data/queries.txt'" in content
    assert "--demo" not in content


def test_write_commands_default_omits_optional_flags(tmp_path):
    """Default args should only have source + output, no optional flags."""
    args = argparse.Namespace(
        demo=True,
        input=None,
        query=None,
        gene=None,
        status=None,
        max_results=20,
        fhir=False,
        ot_min_score=0.6,
        ot_max_diseases=5,
    )
    path = wr.write_commands(args, tmp_path)
    content = path.read_text()
    assert "--status" not in content
    assert "--max-results" not in content
    assert "--fhir" not in content
    assert "--ot-min-score" not in content
    assert "--ot-max-diseases" not in content


# ---------------------------------------------------------------------------
# _to_fhir_date -- all months
# ---------------------------------------------------------------------------


def test_to_fhir_date_all_months():
    """Verify all 12 months normalise correctly."""
    months = [
        ("January", "01"),
        ("February", "02"),
        ("March", "03"),
        ("April", "04"),
        ("May", "05"),
        ("June", "06"),
        ("July", "07"),
        ("August", "08"),
        ("September", "09"),
        ("October", "10"),
        ("November", "11"),
        ("December", "12"),
    ]
    for name, num in months:
        assert api._to_fhir_date(f"{name} 2024") == f"2024-{num}", f"Failed for {name}"


# ---------------------------------------------------------------------------
# count_recruiting edge cases
# ---------------------------------------------------------------------------


def test_count_recruiting_empty_list():
    assert wr.count_recruiting([]) == 0


def test_count_recruiting_all_recruiting():
    trials = [{"status": "RECRUITING"}, {"status": "RECRUITING"}]
    assert wr.count_recruiting(trials) == 2


# ---------------------------------------------------------------------------
# Integration sanity
# ---------------------------------------------------------------------------


def test_demo_data_exists_and_parseable():
    demo = SKILL_DIR / "demo_input.txt"
    assert demo.exists()
    result = api.parse_input(demo)
    assert result["query"]


def test_modules_import_correctly():
    """Verify the module split didn't break cross-imports."""
    assert hasattr(const, "CT_API")
    assert hasattr(api, "fetch_trials")
    assert hasattr(api, "parse_input")
    assert hasattr(wr, "write_report")
    assert hasattr(wr, "write_fhir_bundle")
    assert hasattr(wr, "count_recruiting")


# ===========================================================================
# euctr.py
# ===========================================================================

eu = _load("euctr")


# ---------------------------------------------------------------------------
# euctr -- XML parsing
# ---------------------------------------------------------------------------

MOCK_EUCTR_XML = """<?xml version="1.0"?>
<trials>
  <trial>
    <eudract_number>2024-000001-99</eudract_number>
    <full_title>EU Trial Alpha</full_title>
    <trial_status>Ongoing</trial_status>
    <trial_phase>Phase II</trial_phase>
    <trial_type>Interventional</trial_type>
    <date_on>2024-03-01</date_on>
    <medical_condition>Breast Cancer</medical_condition>
    <trade_name>Tamoxifen</trade_name>
    <primary_end_point>Overall survival</primary_end_point>
  </trial>
  <trial>
    <eudract_number>2024-000002-99</eudract_number>
    <title>EU Trial Beta</title>
    <trial_status>Completed</trial_status>
  </trial>
</trials>"""


def test_euctr_parse_xml():
    trials = eu._parse_euctr_response(MOCK_EUCTR_XML, max_results=10)
    assert len(trials) == 2
    assert trials[0]["nct_id"] == "2024-000001-99"
    assert trials[0]["title"] == "EU Trial Alpha"
    assert trials[0]["conditions"] == ["Breast Cancer"]
    assert trials[0]["interventions"] == ["Tamoxifen"]


def test_euctr_parse_xml_respects_max_results():
    trials = eu._parse_euctr_response(MOCK_EUCTR_XML, max_results=1)
    assert len(trials) == 1


def test_euctr_parse_html_fallback():
    """When XML parse fails, fallback to regex-based EudraCT ID extraction."""
    html = """
    <table><tr>
      <td>2024-000003-99</td>
      <td class="full">Some EU Trial Title</td>
    </tr></table>
    """
    trials = eu._parse_euctr_response(html, max_results=10)
    assert len(trials) == 1
    assert trials[0]["nct_id"] == "EUCTR2024-000003-99"
    assert trials[0]["title"] == "Some EU Trial Title"


def test_euctr_parse_empty_response():
    trials = eu._parse_euctr_response("", max_results=10)
    assert trials == []


def test_euctr_parse_garbage():
    trials = eu._parse_euctr_response("not xml not html nothing", max_results=10)
    assert trials == []


def test_euctr_fetch_network_error():
    """Network errors should return empty list, not raise."""
    with patch("urllib.request.urlopen", side_effect=urllib.error.URLError("refused")):
        result = eu.fetch_euctr("BRCA1")
    assert result == []


def test_euctr_trial_has_source_field():
    """EUCTR trials should have source='euctr' for dedup/display."""
    trials = eu._parse_euctr_response(MOCK_EUCTR_XML, max_results=1)
    assert trials[0].get("source") == "euctr"


def test_euctr_element_to_dict_missing_fields():
    """Graceful handling when XML elements are missing."""
    import xml.etree.ElementTree as ET

    el = ET.fromstring("<trial><eudract_number>2024-000099-99</eudract_number></trial>")
    trial = eu._euctr_element_to_dict(el)
    assert trial["nct_id"] == "2024-000099-99"
    assert trial["title"] == ""
    assert trial["status"] == "UNKNOWN"
    assert trial["conditions"] == []


# ---------------------------------------------------------------------------
# writers.py -- write_csv
# ---------------------------------------------------------------------------


def test_write_csv_creates_file(tmp_path):
    path = wr.write_csv(MOCK_TRIALS, tmp_path)
    assert path.exists()
    assert path.name == "trials.csv"
    assert path.parent.name == "tables"


def test_write_csv_has_header_row(tmp_path):
    wr.write_csv(MOCK_TRIALS, tmp_path)
    content = (tmp_path / "tables" / "trials.csv").read_text()
    header = content.split("\n")[0]
    for col in const.CSV_COLUMNS:
        assert col in header


def test_write_csv_correct_row_count(tmp_path):
    wr.write_csv(MOCK_TRIALS, tmp_path)
    lines = (tmp_path / "tables" / "trials.csv").read_text().strip().split("\n")
    assert len(lines) == len(MOCK_TRIALS) + 1  # header + data


def test_write_csv_lists_are_pipe_delimited(tmp_path):
    trials = [
        {
            **MOCK_TRIALS[0],
            "conditions": ["Cancer A", "Cancer B"],
            "interventions": ["Drug X", "Drug Y"],
        }
    ]
    wr.write_csv(trials, tmp_path)
    content = (tmp_path / "tables" / "trials.csv").read_text()
    assert "Cancer A | Cancer B" in content
    assert "Drug X | Drug Y" in content


def test_write_csv_empty_trials(tmp_path):
    path = wr.write_csv([], tmp_path)
    lines = path.read_text().strip().split("\n")
    assert len(lines) == 1  # header only


# ---------------------------------------------------------------------------
# writers.py -- write_html
# ---------------------------------------------------------------------------


def test_write_html_creates_file(tmp_path):
    path = wr.write_html(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    assert path.exists()
    assert path.name == "report.html"


def test_write_html_contains_trial_data(tmp_path):
    wr.write_html(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.html").read_text()
    for t in MOCK_TRIALS:
        assert t["nct_id"] in content


def test_write_html_contains_disclaimer(tmp_path):
    wr.write_html(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.html").read_text()
    assert "research and educational tool" in content


def test_write_html_is_valid_html(tmp_path):
    wr.write_html(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.html").read_text()
    assert content.strip().startswith("<!DOCTYPE html>")
    assert "</html>" in content


def test_write_html_contains_nct_links(tmp_path):
    wr.write_html(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.html").read_text()
    assert "https://clinicaltrials.gov/study/NCT00000001" in content


def test_write_html_gene_context(tmp_path):
    gene_context = {
        "symbol": "BRCA1",
        "name": "BRCA1 DNA repair associated",
        "diseases": ["breast cancer"],
        "min_score": 0.6,
    }
    wr.write_html(MOCK_QUERY, MOCK_TRIALS, tmp_path, gene_context=gene_context)
    content = (tmp_path / "report.html").read_text()
    assert "BRCA1" in content


def test_write_html_empty_trials(tmp_path):
    path = wr.write_html(MOCK_QUERY, [], tmp_path)
    content = path.read_text()
    assert "0 (0 recruiting)" in content


# ---------------------------------------------------------------------------
# writers.py -- validate_fhir_bundle
# ---------------------------------------------------------------------------


def test_validate_fhir_bundle_valid(tmp_path):
    """A valid bundle should return no errors."""
    path = wr.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads(path.read_text())
    errors = wr.validate_fhir_bundle(bundle)
    assert errors == [], f"Unexpected errors: {errors}"


def test_validate_fhir_bundle_catches_wrong_type():
    errors = wr.validate_fhir_bundle(
        {
            "resourceType": "NotBundle",
            "type": "searchset",
            "timestamp": "2024-01-01T00:00:00+00:00",
            "total": 0,
            "entry": [],
        }
    )
    assert any("resourceType" in e for e in errors)


def test_validate_fhir_bundle_catches_missing_timestamp():
    errors = wr.validate_fhir_bundle(
        {
            "resourceType": "Bundle",
            "type": "searchset",
            "total": 0,
            "entry": [],
        }
    )
    assert any("timestamp" in e for e in errors)


def test_validate_fhir_bundle_catches_total_mismatch():
    errors = wr.validate_fhir_bundle(
        {
            "resourceType": "Bundle",
            "type": "searchset",
            "timestamp": "2024-01-01T00:00:00+00:00",
            "total": 99,
            "entry": [],
        }
    )
    assert any("total" in e for e in errors)


def test_validate_fhir_bundle_catches_invalid_status():
    errors = wr.validate_fhir_bundle(
        {
            "resourceType": "Bundle",
            "type": "searchset",
            "timestamp": "2024-01-01T00:00:00+00:00",
            "total": 1,
            "entry": [
                {
                    "fullUrl": "https://example.com",
                    "resource": {
                        "resourceType": "ResearchStudy",
                        "id": "NCT1",
                        "title": "Test",
                        "status": "INVALID_STATUS",
                        "phase": {"coding": [{"code": "phase-1"}]},
                    },
                }
            ],
        }
    )
    assert any("invalid status" in e for e in errors)


def test_validate_fhir_bundle_catches_invalid_phase():
    errors = wr.validate_fhir_bundle(
        {
            "resourceType": "Bundle",
            "type": "searchset",
            "timestamp": "2024-01-01T00:00:00+00:00",
            "total": 1,
            "entry": [
                {
                    "fullUrl": "https://example.com",
                    "resource": {
                        "resourceType": "ResearchStudy",
                        "id": "NCT1",
                        "title": "Test",
                        "status": "active",
                        "phase": {"coding": [{"code": "INVALID_PHASE"}]},
                    },
                }
            ],
        }
    )
    assert any("invalid phase" in e for e in errors)


# ---------------------------------------------------------------------------
# writers.py -- write_commands with new flags
# ---------------------------------------------------------------------------


def test_write_commands_with_country(tmp_path):
    args = argparse.Namespace(
        demo=True,
        input=None,
        query=None,
        gene=None,
        status=None,
        country="Spain",
        max_results=20,
        fhir=False,
        euctr=False,
        ot_min_score=0.6,
        ot_max_diseases=5,
    )
    path = wr.write_commands(args, tmp_path)
    assert "--country 'Spain'" in path.read_text()


def test_write_commands_with_euctr(tmp_path):
    args = argparse.Namespace(
        demo=True,
        input=None,
        query=None,
        gene=None,
        status=None,
        country=None,
        max_results=20,
        fhir=False,
        euctr=True,
        ot_min_score=0.6,
        ot_max_diseases=5,
    )
    path = wr.write_commands(args, tmp_path)
    assert "--euctr" in path.read_text()


# ---------------------------------------------------------------------------
# writers.py -- write_checksums includes new targets
# ---------------------------------------------------------------------------


def test_write_checksums_includes_html_and_csv(tmp_path):
    (tmp_path / "report.md").write_text("x")
    (tmp_path / "report.html").write_text("<html></html>")
    tables = tmp_path / "tables"
    tables.mkdir()
    (tables / "trials.csv").write_text("nct_id,title\n")
    path = wr.write_checksums(tmp_path)
    content = path.read_text()
    assert "report.html" in content
    assert "tables/trials.csv" in content


# ---------------------------------------------------------------------------
# constants.py -- structural integrity
# ---------------------------------------------------------------------------


def test_month_names_has_12_entries():
    assert len(const.MONTH_NAMES) == 12


def test_fhir_phase_and_display_have_same_keys():
    """FHIR_PHASE and FHIR_PHASE_DISPLAY must stay in sync."""
    assert set(const.FHIR_PHASE.keys()) == set(const.FHIR_PHASE_DISPLAY.keys())


def test_ct_fields_contains_required_fields():
    for field in (
        "NCTId",
        "BriefTitle",
        "OverallStatus",
        "Phase",
        "Condition",
        "InterventionName",
    ):
        assert field in const.CT_FIELDS


def test_demo_data_path_exists():
    assert const.DEMO_DATA.exists()


def test_euctr_api_is_https():
    assert const.EUCTR_API.startswith("https://")


# ---------------------------------------------------------------------------
# opentargets.py -- get_diseases edge cases
# ---------------------------------------------------------------------------


def test_ot_get_diseases_respects_max_results():
    """get_diseases should truncate to max_results even if more pass the score."""
    resp = _make_ot_response(
        {
            "target": {
                "associatedDiseases": {
                    "rows": [
                        {
                            "disease": {"id": f"D{i}", "name": f"Disease {i}"},
                            "score": 0.9,
                        }
                        for i in range(10)
                    ]
                }
            }
        }
    )
    with patch("urllib.request.urlopen", return_value=resp):
        diseases = ot.get_diseases("ENSG00000012048", min_score=0.5, max_results=3)
    assert len(diseases) == 3


def test_ot_get_diseases_empty_rows():
    resp = _make_ot_response({"target": {"associatedDiseases": {"rows": []}}})
    with patch("urllib.request.urlopen", return_value=resp):
        diseases = ot.get_diseases("ENSG00000012048")
    assert diseases == []


def test_ot_get_diseases_all_below_threshold():
    resp = _make_ot_response(
        {
            "target": {
                "associatedDiseases": {
                    "rows": [
                        {"disease": {"id": "D1", "name": "Weak"}, "score": 0.1},
                        {"disease": {"id": "D2", "name": "Also Weak"}, "score": 0.2},
                    ]
                }
            }
        }
    )
    with patch("urllib.request.urlopen", return_value=resp):
        diseases = ot.get_diseases("ENSG00000012048", min_score=0.6)
    assert diseases == []


def test_ot_disease_namedtuple_fields():
    d = ot.Disease(id="MONDO_1", name="Cancer", score=0.9)
    assert d.id == "MONDO_1"
    assert d.name == "Cancer"
    assert d.score == 0.9


def test_ot_malformed_json():
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = b"not json at all"
    with patch("urllib.request.urlopen", return_value=mock_resp):
        with pytest.raises(RuntimeError, match="malformed JSON"):
            ot.resolve_gene("BRCA1")


# ---------------------------------------------------------------------------
# api.py -- _normalise_trial boundary and edge cases
# ---------------------------------------------------------------------------


def test_normalise_trial_summary_301_chars():
    """Summary of 301 chars should be truncated to 300 + '...'."""
    summary_301 = "x" * 301
    raw = {
        "protocolSection": {
            "identificationModule": {"nctId": "NCT1", "briefTitle": "T"},
            "statusModule": {},
            "designModule": {},
            "descriptionModule": {"briefSummary": summary_301},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = api._normalise_trial(raw)
    assert trial["summary"].endswith("...")
    assert len(trial["summary"]) == 303


def test_normalise_trial_completely_empty():
    """Completely empty dict should not crash."""
    trial = api._normalise_trial({})
    assert trial["nct_id"] == ""
    assert trial["status"] == "UNKNOWN"


def test_normalise_trial_dates_are_normalised():
    """Verify date normalisation is actually applied."""
    raw = {
        "protocolSection": {
            "identificationModule": {},
            "statusModule": {
                "startDateStruct": {"date": "January 2024"},
                "completionDateStruct": {"date": "December 2026"},
            },
            "designModule": {},
            "descriptionModule": {},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = api._normalise_trial(raw)
    assert trial["start_date"] == "2024-01"
    assert trial["completion_date"] == "2026-12"


def test_normalise_trial_no_interventions_key():
    """Missing interventions key should yield empty list."""
    raw = {
        "protocolSection": {
            "identificationModule": {},
            "statusModule": {},
            "designModule": {},
            "descriptionModule": {},
            "conditionsModule": {},
            "armsInterventionsModule": {},
        }
    }
    trial = api._normalise_trial(raw)
    assert trial["interventions"] == []


# ---------------------------------------------------------------------------
# api.py -- fetch_trials with mock response containing real study data
# ---------------------------------------------------------------------------


def test_fetch_trials_returns_normalised_records():
    """Verify fetch_trials calls _normalise_trial on each study."""
    raw_study = {
        "protocolSection": {
            "identificationModule": {
                "nctId": "NCT99999999",
                "briefTitle": "Mock Trial",
            },
            "statusModule": {"overallStatus": "RECRUITING"},
            "designModule": {"phases": ["PHASE3"], "studyType": "INTERVENTIONAL"},
            "descriptionModule": {"briefSummary": "A mock trial."},
            "conditionsModule": {"conditions": ["Test Condition"]},
            "armsInterventionsModule": {"interventions": [{"name": "Mock Drug"}]},
        }
    }
    mock_resp = MagicMock()
    mock_resp.__enter__ = lambda s: s
    mock_resp.__exit__ = MagicMock(return_value=False)
    mock_resp.read.return_value = json.dumps({"studies": [raw_study]}).encode()
    with patch("urllib.request.urlopen", return_value=mock_resp):
        trials = api.fetch_trials("test")
    assert len(trials) == 1
    t = trials[0]
    assert t["nct_id"] == "NCT99999999"
    assert t["status"] == "RECRUITING"
    assert t["phase"] == "PHASE3"
    assert t["conditions"] == ["Test Condition"]
    assert t["interventions"] == ["Mock Drug"]


# ---------------------------------------------------------------------------
# writers.py -- FHIR bundle metadata
# ---------------------------------------------------------------------------


def test_fhir_bundle_has_meta_tag(tmp_path):
    path = wr.write_fhir_bundle(MOCK_TRIALS, tmp_path)
    bundle = json.loads(path.read_text())
    assert "meta" in bundle
    assert bundle["meta"]["tag"][0]["system"] == "https://clawbio.ai"
    assert bundle["meta"]["tag"][0]["code"] == "clinical-trial-finder"


def test_fhir_resource_has_profile(tmp_path):
    wr.write_fhir_bundle([MOCK_TRIALS[0]], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    resource = bundle["entry"][0]["resource"]
    assert "meta" in resource
    assert (
        "http://hl7.org/fhir/StructureDefinition/ResearchStudy"
        in resource["meta"]["profile"]
    )


def test_fhir_period_start_only(tmp_path):
    """Trial with start_date but no completion_date should have period.start only."""
    trial = {**MOCK_TRIALS[0], "completion_date": ""}
    wr.write_fhir_bundle([trial], tmp_path)
    bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
    period = bundle["entry"][0]["resource"]["period"]
    assert "start" in period
    assert "end" not in period


# ---------------------------------------------------------------------------
# writers.py -- write_report content accuracy
# ---------------------------------------------------------------------------


def test_write_report_summary_metrics_match(tmp_path):
    """Summary table metrics should match actual trial data."""
    wr.write_report(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "Total trials | 3" in content
    assert "Recruiting now | 1" in content


def test_write_report_unknown_status_label(tmp_path):
    """Trial with UNKNOWN status should get [UNKNOWN] label."""
    trial_unknown = {**MOCK_TRIALS[0], "status": "UNKNOWN", "nct_id": "NCT_UNK"}
    wr.write_report(MOCK_QUERY, [trial_unknown], tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "[UNKNOWN]" in content


def test_write_report_conditions_capped_at_3(tmp_path):
    """Report should show at most 3 conditions per trial."""
    trial = {
        **MOCK_TRIALS[0],
        "conditions": [
            "AlphaDisease",
            "BetaDisease",
            "GammaDisease",
            "DeltaDisease",
            "EpsilonDisease",
        ],
    }
    wr.write_report(MOCK_QUERY, [trial], tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "AlphaDisease, BetaDisease, GammaDisease" in content
    assert "DeltaDisease" not in content


def test_write_report_interventions_capped_at_3(tmp_path):
    """Report should show at most 3 interventions per trial."""
    trial = {**MOCK_TRIALS[0], "interventions": ["X", "Y", "Z", "W"]}
    wr.write_report(MOCK_QUERY, [trial], tmp_path)
    content = (tmp_path / "report.md").read_text()
    assert "X, Y, Z" in content
    assert "W" not in content


# ---------------------------------------------------------------------------
# writers.py -- write_summary JSON structure
# ---------------------------------------------------------------------------


def test_write_summary_source_is_correct(tmp_path):
    path = wr.write_summary(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    data = json.loads(path.read_text())
    assert data["source"] == "clinicaltrials.gov/api/v2"


def test_write_summary_timestamp_is_iso(tmp_path):
    path = wr.write_summary(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    data = json.loads(path.read_text())
    # ISO format: starts with YYYY-MM-DD
    assert data["timestamp"][:4].isdigit()
    assert data["timestamp"][4] == "-"


def test_write_summary_query_matches_input(tmp_path):
    path = wr.write_summary(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    data = json.loads(path.read_text())
    assert data["query"] == MOCK_QUERY["query"]


def test_write_summary_trials_preserve_nct_ids(tmp_path):
    path = wr.write_summary(MOCK_QUERY, MOCK_TRIALS, tmp_path)
    data = json.loads(path.read_text())
    nct_ids = {t["nct_id"] for t in data["trials"]}
    assert nct_ids == {"NCT00000001", "NCT00000002", "NCT00000003"}


# ---------------------------------------------------------------------------
# writers.py -- write_phase_chart edge cases
# ---------------------------------------------------------------------------


def test_write_phase_chart_empty_trials(tmp_path):
    """Empty trial list should return None (no chart to draw)."""
    result = wr.write_phase_chart([], tmp_path)
    assert result is None


# ---------------------------------------------------------------------------
# writers.py -- write_commands/write_checksums file format
# ---------------------------------------------------------------------------


def test_write_commands_ends_with_newline(tmp_path):
    args = argparse.Namespace(
        demo=True,
        input=None,
        query=None,
        gene=None,
        status=None,
        max_results=20,
        fhir=False,
        ot_min_score=0.6,
        ot_max_diseases=5,
    )
    path = wr.write_commands(args, tmp_path)
    assert path.read_text().endswith("\n")


def test_write_checksums_ends_with_newline(tmp_path):
    (tmp_path / "report.md").write_text("x")
    path = wr.write_checksums(tmp_path)
    assert path.read_text().endswith("\n")


# ---------------------------------------------------------------------------
# FHIR -- all phase mappings produce valid codes
# ---------------------------------------------------------------------------


def test_fhir_all_phase_mappings_individually(tmp_path):
    """Test each CT.gov phase string maps to a valid FHIR code."""
    for ct_phase, expected_code in const.FHIR_PHASE.items():
        trial = {**MOCK_TRIALS[0], "phase": ct_phase}
        wr.write_fhir_bundle([trial], tmp_path)
        bundle = json.loads((tmp_path / "fhir_bundle.json").read_text())
        code = bundle["entry"][0]["resource"]["phase"]["coding"][0]["code"]
        assert code == expected_code, (
            f"Phase '{ct_phase}' mapped to '{code}', expected '{expected_code}'"
        )
        assert code in const.FHIR_VALID_PHASES, (
            f"Code '{code}' not in FHIR_VALID_PHASES"
        )
