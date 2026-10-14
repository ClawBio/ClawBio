"""Behavior contracts frozen before extracting shared nf-core mechanics."""
from __future__ import annotations
import hashlib
import importlib
import json
from pathlib import Path
from types import SimpleNamespace

import pytest

_PIPELINE = Path(__file__).resolve().parent.parent.name.split("-")[1]


def test_command_contract(tmp_path):
    builder = importlib.import_module("command_builder")
    command, display = builder.build_nextflow_command(
        pipeline_source={"source_kind": "remote_repo", "source_ref": f"nf-core/{_PIPELINE}", "resolved_version": "1.2.3"},
        profile="docker,arm64", params_path=tmp_path / "params with space.yaml",
        resume=True, work_dir="s3://bucket/work", extra_configs=[Path("first.config"), Path("second.config")],
    )
    assert command == ["nextflow", "run", f"nf-core/{_PIPELINE}", "-r", "1.2.3", "-profile", "docker,arm64", "-params-file", str(tmp_path / "params with space.yaml"), "-work-dir", "s3://bucket/work", "-c", "first.config", "-c", "second.config", "-resume"]
    assert display.endswith("-c first.config -c second.config -resume")


def test_error_payload_contract():
    errors = importlib.import_module("errors")
    err = errors.SkillError("preflight", "INVALID_RESUME_STATE", "rejected", "retry", {"old": False})
    assert str(err) == "rejected"
    assert err.to_dict() == {"status": "error", "ok": False, "stage": "preflight", "error_code": "INVALID_RESUME_STATE", "message": "rejected", "fix": "retry", "details": {"old": False}}


def test_bundle_contract(tmp_path):
    params = importlib.import_module("params_builder")
    path = params.write_params_yaml({"second": False, "first": 7}, output_dir=tmp_path)
    expected = ("# 'input' and 'outdir' are relative to the Nextflow launch directory (output_dir).\n"
                "# Use reproducibility/commands.sh for the recommended replay path.\n"
                f"# Manual replay: cd <output_dir> && nextflow run nf-core/{_PIPELINE} -params-file reproducibility/params.yaml\n"
                "second: false\nfirst: 7\n")
    assert path == tmp_path / "reproducibility" / "params.yaml"
    assert path.read_bytes() == expected.encode()


def test_check_mode_artifact_contract(tmp_path):
    wrapper = importlib.import_module(f"nfcore_{_PIPELINE}_wrapper")
    source = {"source_kind": "remote_repo", "source_ref": f"nf-core/{_PIPELINE}", "resolved_version": "1.2.3"}
    if _PIPELINE == "sarek":
        rc = wrapper._write_check_report(output_dir=tmp_path, args=SimpleNamespace(demo=False, profile="docker"), samplesheet_report={"sample_count": 2}, preflight_result=SimpleNamespace(warnings=["notice"], notes=[]), pipeline_source=source)
    else:
        rc = wrapper._write_check_mode_result(tmp_path, preflight_result={"ok": True}, pipeline_source=source)
    assert rc == 0
    assert {p.name for p in tmp_path.iterdir()} == {"check_result.json"}
    payload = json.loads((tmp_path / "check_result.json").read_text())
    assert payload["ok"] is True
    assert payload["pipeline_source"] == source
    assert not (tmp_path / "result.json").exists()


def test_old_manifest_resume_contract(tmp_path):
    wrapper = importlib.import_module(f"nfcore_{_PIPELINE}_wrapper")
    result = tmp_path / "result.json"
    result.write_text('{"status":"success","old":true}')
    wrapper.write_error_result(tmp_path, {"error_code": "INVALID_RESUME_STATE", "message": "new"})
    assert result.read_text() == '{"status":"success","old":true}'


def test_parameter_fingerprint_contract():
    values = {"second": False, "first": 7}
    if _PIPELINE == "sarek":
        api = importlib.import_module("provenance")
        expected = "sha256:" + hashlib.sha256(b'{"first":7,"second":false}').hexdigest()
        assert api.compute_params_checksum(values) == expected
    elif _PIPELINE == "rnaseq":
        api = importlib.import_module("preflight")
        expected = "sha256:" + hashlib.sha256(b'second: false\nfirst: 7\n').hexdigest()
        assert api.params_payload_checksum(values) == expected
    else:
        api = importlib.import_module("nfcore_scrnaseq_wrapper")
        header = ("# 'input' and 'outdir' are relative to the Nextflow launch directory (output_dir).\n"
                  "# Use reproducibility/commands.sh for the recommended replay path.\n"
                  "# Manual replay: cd <output_dir> && nextflow run nf-core/scrnaseq -params-file reproducibility/params.yaml\n")
        expected = hashlib.sha256((header + "second: false\nfirst: 7\n").encode()).hexdigest()
        assert api._params_payload_checksum(values) == expected


@pytest.mark.parametrize("value", ["nan", "inf", "-inf", "1e309"])
def test_non_finite_timeout_is_a_usage_error(value, tmp_path, capsys):
    wrapper = importlib.import_module(f"nfcore_{_PIPELINE}_wrapper")
    with pytest.raises(SystemExit) as exit_info:
        wrapper.build_parser().parse_args(["--output", str(tmp_path), "--timeout-hours", value])
    assert exit_info.value.code == 2
    assert "--timeout-hours" in capsys.readouterr().err
