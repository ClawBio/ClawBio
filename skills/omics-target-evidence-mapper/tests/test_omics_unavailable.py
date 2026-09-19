"""Request failures must not be rendered as empty evidence searches."""
from __future__ import annotations

import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

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


def _fake_json_mixed_success(method, url, **kwargs):
    if url.startswith("https://rest.uniprot.org"):
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
    if url.startswith("https://eutils.ncbi.nlm.nih.gov/entrez/eutils/esearch.fcgi"):
        return None
    if url.startswith("https://clinicaltrials.gov"):
        return {"studies": []}
    if url.startswith("https://api.platform.opentargets.org"):
        return None
    return None


def test_mixed_success_still_renders_usable_sources(monkeypatch) -> None:
    monkeypatch.setattr(mapper, "safe_request_json", _fake_json_mixed_success)
    evidence = mapper.build_evidence(
        argparse.Namespace(demo=False, gene="IL6R", disease="CAD", max_papers=5, max_trials=5)
    )
    assert evidence["target_summary"]["status"] == "ok"
    assert evidence["literature"]["status"] == "unavailable"
    assert evidence["trials"]["status"] == "no_result"
    assert "PubMed was unavailable; literature was not assessed." in evidence["limitations"]
    assert "ClinicalTrials.gov was unavailable; trials were not assessed." not in evidence["limitations"]
    report = mapper.build_report(evidence)
    assert "P40189" in report
    assert "Literature unavailable; not assessed." in report
    assert "No trial hits found." in report
    assert "PubMed was unavailable; literature was not assessed." in report


def test_main_writes_item_counts_and_source_status(tmp_path, monkeypatch) -> None:
    output_dir = tmp_path / "output"
    monkeypatch.setattr(
        mapper,
        "parse_args",
        lambda: argparse.Namespace(
            demo=False,
            gene="IL6R",
            disease="CAD",
            output=str(output_dir),
            max_papers=5,
            max_trials=5,
        ),
    )
    monkeypatch.setattr(mapper, "safe_request_json", _fake_json_mixed_success)

    def write_ro_crate(output_dir, **_kwargs):
        path = Path(output_dir) / "ro-crate-metadata.json"
        path.write_text("{}\n", encoding="utf-8")
        return path

    monkeypatch.setattr(mapper, "write_ro_crate", write_ro_crate)

    mapper.main()

    metadata = json.loads((output_dir / "metadata.json").read_text(encoding="utf-8"))
    evidence = json.loads((output_dir / "evidence.json").read_text(encoding="utf-8"))
    assert metadata["counts"]["literature"] == 0
    assert metadata["counts"]["trials"] == 0
    assert metadata["counts"] != {
        "literature": len(evidence["literature"]),
        "trials": len(evidence["trials"]),
    }
    assert metadata["status"]["literature"] == "unavailable"
    assert metadata["status"]["trials"] == "no_result"
    assert metadata["counts"]["literature"] == len(evidence["literature"]["items"])
    assert metadata["counts"]["trials"] == len(evidence["trials"]["items"])

    manifest = output_dir / "reproducibility" / "checksums.sha256"
    entries = dict(line.split("  ", 1) for line in manifest.read_text().splitlines())
    report = (output_dir / "report.md").read_text(encoding="utf-8")
    assert "Literature unavailable; not assessed." in report
    assert "PubMed was unavailable; literature was not assessed." in report
    assert hashlib.sha256(report.encode("utf-8")).hexdigest() in entries
