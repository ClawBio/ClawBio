"""Offline contract tests. All cell/donor records in this module are synthetic."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SCRIPT = Path(__file__).with_name("atlas_preflight.py")
SPEC = importlib.util.spec_from_file_location("atlas_preflight", SCRIPT)
audit = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(audit)


def save_json(path, value):
    path.write_text(json.dumps(value), encoding="utf-8")
    return path


def make_manifest(tmp_path, train=None, test=None):
    records = {
        "train": train if train is not None else [
            {"cell_id": "SYN-A:1", "donor_id": "SYN-A", "sex": "female"},
            {"cell_id": "SYN-A:2", "donor_id": "SYN-A", "sex": "female"},
        ],
        "test": test if test is not None else [
            {"cell_id": "SYN-B:1", "donor_id": "SYN-B", "sex": "male"},
        ],
    }
    manifest = {"schema_version": 1, "splits": {}}
    for name, rows in records.items():
        path = tmp_path / (name + ".jsonl")
        path.write_text("".join(json.dumps(r) + "\n" for r in rows), encoding="utf-8")
        manifest["splits"][name] = {
            "dataset_id": "synthetic-atlas",
            "dataset_version": "1.0" if name == "train" else "2.0",
            "source_url": "https://example.org/synthetic-atlas",
            "identity_namespace": "synthetic-study-stable-ids-v1",
            "observations_path": path.name,
            "observations_sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
            "expected_cells": len(rows),
        }
    return save_json(tmp_path / "manifest.json", manifest)


def change(path, edit):
    value = json.loads(path.read_text())
    edit(value)
    save_json(path, value)


def codes(result):
    return {f["code"] for f in result["findings"]}


def test_disjoint_declared_metadata_passes_but_does_not_prove_model_clean(tmp_path):
    path = make_manifest(tmp_path)
    result = audit.audit_manifest(path)
    assert result["status"] == "PASS"
    assert result["scope"] == "declared_observation_metadata_only"
    assert result["model_training_contamination"] == "not_assessed"
    assert result["splits"]["train"]["cells"] == 2
    assert result["splits"]["train"]["donors"] == 1
    assert result["overlap"] == {"cells": 0, "donors": 0}
    assert "ANCESTRY_COVERAGE_INCOMPLETE" in codes(result)
    assert result["provenance"]["manifest_sha256"] == hashlib.sha256(path.read_bytes()).hexdigest()
    assert result == audit.audit_manifest(path)


def test_cell_and_donor_overlap_block_combined_v2(tmp_path):
    path = make_manifest(tmp_path, test=[{"cell_id": "SYN-A:1", "donor_id": "SYN-A"}])
    result = audit.audit_manifest(path)
    assert result["status"] == "FAIL"
    assert {"CELL_OVERLAP", "DONOR_OVERLAP"} <= codes(result)
    assert result["overlap"] == {"cells": 1, "donors": 1}


def test_distinct_cells_from_same_donor_fail_donor_holdout(tmp_path):
    path = make_manifest(tmp_path, test=[{"cell_id": "SYN-A:NEW", "donor_id": "SYN-A"}])
    result = audit.audit_manifest(path)
    assert result["status"] == "FAIL"
    assert "DONOR_OVERLAP" in codes(result)
    assert "CELL_OVERLAP" not in codes(result)


@pytest.mark.parametrize("field", ["cell_id", "donor_id"])
@pytest.mark.parametrize("value", [None, "", "unknown", "NA", "nan", 3, [], " A "])
def test_invalid_or_unknown_identity_fails_closed(tmp_path, field, value):
    row = {"cell_id": "SYN-B:1", "donor_id": "SYN-B", field: value}
    result = audit.audit_manifest(make_manifest(tmp_path, test=[row]))
    assert result["status"] == "ERROR"
    assert "INVALID_INPUT" in codes(result)


def test_absent_identity_rejected(tmp_path):
    result = audit.audit_manifest(make_manifest(tmp_path, test=[{"cell_id": "SYN-B:1"}]))
    assert result["status"] == "ERROR"


def test_namespace_difference_is_unverifiable_not_disjoint(tmp_path):
    path = make_manifest(tmp_path)
    change(path, lambda m: m["splits"]["test"].update(identity_namespace="renamed-ids"))
    result = audit.audit_manifest(path)
    assert result["status"] == "ERROR"
    assert "IDENTITY_NAMESPACE_MISMATCH" in codes(result)
    assert result["overlap"] is None


@pytest.mark.parametrize("version", ["latest", "stable", "main", "HEAD", "current", ""])
def test_floating_or_empty_version_rejected(tmp_path, version):
    path = make_manifest(tmp_path)
    change(path, lambda m: m["splits"]["test"].update(dataset_version=version))
    assert audit.audit_manifest(path)["status"] == "ERROR"


def test_same_version_with_disjoint_ids_is_not_itself_leakage(tmp_path):
    path = make_manifest(tmp_path)
    change(path, lambda m: m["splits"]["test"].update(dataset_version="1.0"))
    assert audit.audit_manifest(path)["status"] == "PASS"


@pytest.mark.parametrize("field", ["expected_cells", "observations_sha256", "source_url", "dataset_id"])
def test_missing_provenance_rejected(tmp_path, field):
    path = make_manifest(tmp_path)
    change(path, lambda m: m["splits"]["test"].pop(field))
    assert audit.audit_manifest(path)["status"] == "ERROR"


def test_tampered_bytes_fail_hash_check(tmp_path):
    path = make_manifest(tmp_path)
    with (tmp_path / "test.jsonl").open("a") as handle:
        handle.write("\n")
    result = audit.audit_manifest(path)
    assert result["status"] == "ERROR"
    assert "CHECKSUM_MISMATCH" in codes(result)


def test_incomplete_export_fails_count_check(tmp_path):
    path = make_manifest(tmp_path)
    change(path, lambda m: m["splits"]["test"].update(expected_cells=100))
    assert "CELL_COUNT_MISMATCH" in codes(audit.audit_manifest(path))


def test_duplicate_cell_within_split_fails(tmp_path):
    row = {"cell_id": "SYN-A:1", "donor_id": "SYN-A"}
    result = audit.audit_manifest(make_manifest(tmp_path, train=[row, row]))
    assert result["status"] == "FAIL"
    assert "DUPLICATE_CELL_ID" in codes(result)


def test_empty_split_fails(tmp_path):
    assert audit.audit_manifest(make_manifest(tmp_path, test=[]))["status"] == "ERROR"


def test_population_summary_counts_donors_not_cells(tmp_path):
    path = make_manifest(tmp_path)
    result = audit.audit_manifest(path)
    coverage = result["splits"]["train"]["donor_metadata"]["sex"]
    assert coverage == {"known_donors": 1, "unknown_donors": 0, "conflicting_donors": 0,
                        "distribution": {"female": 1}}
    assert result["splits"]["train"]["donor_metadata"]["ancestry"]["unknown_donors"] == 1


def test_inconsistent_demographics_warn_and_do_not_invent_consensus(tmp_path):
    rows = [{"cell_id": "SYN-A:1", "donor_id": "SYN-A", "sex": "female"},
            {"cell_id": "SYN-A:2", "donor_id": "SYN-A", "sex": "male"}]
    result = audit.audit_manifest(make_manifest(tmp_path, train=rows))
    assert "SEX_METADATA_CONFLICT" in codes(result)
    assert result["splits"]["train"]["donor_metadata"]["sex"]["distribution"] == {}


@pytest.mark.parametrize("contents", ['{', '[]', '{"schema_version": 2}', '{"schema_version":1,"schema_version":1}'])
def test_malformed_manifest_reports_error(tmp_path, contents):
    path = tmp_path / "bad.json"
    path.write_text(contents)
    assert audit.audit_manifest(path)["status"] == "ERROR"


def test_remote_observations_path_never_fetched(tmp_path):
    path = make_manifest(tmp_path)
    change(path, lambda m: m["splits"]["test"].update(observations_path="https://example.org/cells.jsonl"))
    assert audit.audit_manifest(path)["status"] == "ERROR"


def run_cli(*args):
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], text=True, capture_output=True)


@pytest.mark.parametrize("overlap,exit_code", [(False, 0), (True, 1)])
def test_cli_writes_verifiable_bundle_and_respects_exit_status(tmp_path, overlap, exit_code):
    rows = [{"cell_id": "SYN-A:1", "donor_id": "SYN-A"}] if overlap else None
    path = make_manifest(tmp_path, test=rows)
    out = tmp_path / "report with spaces"
    completed = run_cli("--input", path, "--output", out)
    assert completed.returncode == exit_code, completed.stderr
    result = json.loads((out / "result.json").read_text())
    assert result["status"] == ("FAIL" if overlap else "PASS")
    assert "not a medical device" in (out / "report.md").read_text()
    assert (out / "reproducibility" / "manifest.json").read_bytes() == path.read_bytes()
    assert "'" in (out / "reproducibility" / "commands.sh").read_text()
    for line in (out / "reproducibility" / "checksums.sha256").read_text().splitlines():
        digest, relative = line.split("  ", 1)
        assert hashlib.sha256((out / relative).read_bytes()).hexdigest() == digest
    before = (out / "result.json").read_bytes()
    assert run_cli("--input", path, "--output", out).returncode == 2
    assert (out / "result.json").read_bytes() == before


def test_cli_error_is_reported_and_nonzero(tmp_path):
    path = make_manifest(tmp_path)
    change(path, lambda m: m["splits"]["test"].update(expected_cells=999))
    out = tmp_path / "error_report"
    assert run_cli("--input", path, "--output", out).returncode == 2
    assert json.loads((out / "result.json").read_text())["status"] == "ERROR"


def test_demo_is_synthetic_and_deliberately_detects_leakage(tmp_path):
    out = tmp_path / "demo"
    result = run_cli("--demo", "--output", out)
    assert result.returncode == 1, result.stderr
    payload = json.loads((out / "result.json").read_text())
    assert payload["status"] == "FAIL"
    assert payload["synthetic_demo"] is True
    assert payload["overlap"] == {"cells": 1, "donors": 1}


def test_cli_requires_input_or_demo_and_never_both(tmp_path):
    assert run_cli("--output", tmp_path / "no_input").returncode == 2
    assert run_cli("--input", "x", "--demo", "--output", tmp_path / "both").returncode == 2
