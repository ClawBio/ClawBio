"""Request failures must not be rendered as empty evidence searches."""
from __future__ import annotations

import argparse

import pytest

from pathlib import Path
import importlib.util

MODULE_PATH = Path(__file__).parents[1] / "omics_target_evidence_mapper.py"
SPEC = importlib.util.spec_from_file_location("omics_target_evidence_mapper", MODULE_PATH)
mapper = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(mapper)


def test_uniprot_request_failure_is_unavailable_not_no_result(monkeypatch) -> None:
    monkeypatch.setattr(mapper, "safe_request_json", lambda *args, **kwargs: None)
    result = mapper.fetch_uniprot_summary("IL6R")
    assert result["status"] == "unavailable"
    assert result["gene"] == "IL6R"


def test_uniprot_empty_payload_is_no_result(monkeypatch) -> None:
    monkeypatch.setattr(mapper, "safe_request_json", lambda *args, **kwargs: {"results": []})
    result = mapper.fetch_uniprot_summary("IL6R")
    assert result["status"] == "no_result"


def test_pubmed_and_trials_distinguish_failure_from_empty(monkeypatch) -> None:
    monkeypatch.setattr(mapper, "safe_request_json", lambda *args, **kwargs: None)
    literature = mapper.fetch_pubmed_hits("IL6R", "CAD", 5)
    trials = mapper.fetch_trials("IL6R", "CAD", 5)
    assert literature == {"status": "unavailable", "items": []}
    assert trials == {"status": "unavailable", "items": []}

    monkeypatch.setattr(
        mapper,
        "safe_request_json",
        lambda *args, **kwargs: {"esearchresult": {"idlist": []}, "studies": []},
    )
    literature = mapper.fetch_pubmed_hits("IL6R", "CAD", 5)
    trials = mapper.fetch_trials("IL6R", "CAD", 5)
    assert literature["status"] == "no_result"
    assert literature["items"] == []
    assert trials["status"] == "no_result"
    assert trials["items"] == []


def test_report_wording_differs_for_unavailable_and_empty_sources() -> None:
    base = {
        "query": {"gene": "IL6R", "disease": "CAD", "demo_mode": False},
        "target_summary": {"status": "unavailable", "gene": "IL6R"},
        "disease_association": {"status": "unavailable", "gene": "IL6R", "disease": "CAD"},
        "limitations": [],
        "provenance": {"sources": [], "generated_at_utc": "2026-09-15T00:00:00+00:00", "version": "0.1.0"},
    }
    failed = mapper.build_report(
        {
            **base,
            "literature": {"status": "unavailable", "items": []},
            "trials": {"status": "unavailable", "items": []},
        }
    )
    empty = mapper.build_report(
        {
            **base,
            "literature": {"status": "no_result", "items": []},
            "trials": {"status": "no_result", "items": []},
        }
    )
    assert "Literature unavailable; not assessed." in failed
    assert "Trials unavailable; not assessed." in failed
    assert "No literature hits found." not in failed
    assert "No literature hits found." in empty
    assert "No trial hits found." in empty


def test_mixed_success_still_renders_usable_sources(monkeypatch) -> None:
    def fake_json(method, url, **kwargs):
        if "uniprot" in url:
            return {
                "results": [
                    {
                        "primaryAccession": "P40189",
                        "uniProtkbId": "IL6RB_HUMAN",
                        "proteinDescription": {"recommendedName": {"fullName": {"value": "IL-6R"}}},
                        "organism": {"scientificName": "Homo sapiens"},
                    }
                ]
            }
        if "esearch" in url:
            return None
        if "clinicaltrials.gov" in url:
            return {"studies": []}
        if "opentargets" in url:
            return None
        return None

    monkeypatch.setattr(mapper, "safe_request_json", fake_json)
    evidence = mapper.build_evidence(
        argparse.Namespace(demo=False, gene="IL6R", disease="CAD", max_papers=5, max_trials=5)
    )
    assert evidence["target_summary"]["status"] == "ok"
    assert evidence["literature"]["status"] == "unavailable"
    assert evidence["trials"]["status"] == "no_result"
    report = mapper.build_report(evidence)
    assert "P40189" in report
    assert "Literature unavailable; not assessed." in report
    assert "No trial hits found." in report
