from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from unittest.mock import patch

from provenance import write_provenance_bundle, write_reproducibility_checksums


def test_write_reproducibility_checksums_normalises_types(tmp_path):
    """All paths in the checksum list must be Path objects before write_checksums is called."""
    samplesheet = tmp_path / "samplesheet.csv"
    params = tmp_path / "params.yaml"
    stdout = tmp_path / "stdout.txt"
    stderr = tmp_path / "stderr.txt"
    for f in (samplesheet, params, stdout, stderr):
        f.write_text("x", encoding="utf-8")

    collected: list = []

    def fake_write_checksums(paths, output_dir, *, anchor):
        collected.extend(paths)
        return output_dir / "checksums.sha256"

    with patch("provenance.write_checksums", fake_write_checksums):
        write_reproducibility_checksums(
            tmp_path,
            normalized_samplesheet=samplesheet,
            params_path=params,
            preflight_result={"references": {}},
            parsed_outputs={"h5ad_candidates": [], "multiqc_report": ""},
            execution_result={
                "stdout_path": str(stdout),
                "stderr_path": str(stderr),
            },
        )

    assert all(isinstance(p, Path) for p in collected), \
        f"All checksum paths must be Path objects, got: {[type(p) for p in collected]}"


def test_write_provenance_bundle(tmp_path):
    output_dir = tmp_path
    params_path = output_dir / "reproducibility" / "params.yaml"
    params_path.parent.mkdir(parents=True)
    params_path.write_text("{}", encoding="utf-8")
    normalized = output_dir / "reproducibility" / "samplesheet.valid.csv"
    normalized.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    stdout = output_dir / "logs" / "stdout.txt"
    stderr = output_dir / "logs" / "stderr.txt"
    stdout.parent.mkdir(parents=True)
    stdout.write_text("", encoding="utf-8")
    stderr.write_text("", encoding="utf-8")
    preferred = output_dir / "upstream" / "results" / "combined_filtered_matrix.h5ad"
    preferred.parent.mkdir(parents=True)
    preferred.write_text("h5ad", encoding="utf-8")
    args = Namespace(demo=False, check=False, preset="star", profile="docker", pipeline_version="4.1.0", resume=False)
    provenance_dir, manifest_path = write_provenance_bundle(
        output_dir,
        args=args,
        pipeline_source={"source_kind": "local_checkout", "source_ref": "scrnaseq-checkout", "resolved_version": "abc123", "branch": "main", "dirty": False},
        preflight_result={"java": {"version": "17.0.8"}, "nextflow": {"version": "25.04.0"}, "references": {}, "profile": {"profile": "docker"}},
        params_path=params_path,
        params_payload={"aligner": "star"},
        normalized_samplesheet=normalized,
        samplesheet_summary={"sample_count": 1, "fastq_paths": [], "unknown_columns": []},
        parsed_outputs={"preferred_h5ad": str(preferred), "multiqc_report": "", "pipeline_info_dir": "", "h5ad_candidates": [str(preferred)], "rds_candidates": [], "handoff_available": True},
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str="nextflow run ...",
    )
    assert (provenance_dir / "runtime.json").exists()
    assert manifest_path.exists()
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    runtime = json.loads((provenance_dir / "runtime.json").read_text(encoding="utf-8"))
    assert manifest["preset"] == "star"
    assert manifest["profile"] == "docker"
    assert manifest["params_checksum"]
    assert "pipeline_source" in manifest
    assert manifest["pipeline_source"]["source_kind"] == "local_checkout"
    assert manifest["pipeline_source"]["resolved_version"] == "abc123"
    assert runtime["work_dir"] == str(output_dir / "upstream" / "work")
    assert manifest["java_version"] == "17.0.8"
    assert manifest["nextflow_version"] == "25.04.0"
    assert manifest["python_version"] == runtime["python_version"]
    assert manifest["environment_yml_mode"] == "install_recipe"


def test_build_upstream_payload_uses_schema_constant():
    from provenance import build_upstream_payload
    from schemas import DEFAULT_REMOTE_PIPELINE
    source = {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/scrnaseq",
        "resolved_version": "4.1.0",
        "branch": "",
        "dirty": False,
    }
    payload = build_upstream_payload(source)
    assert payload["pipeline"] == DEFAULT_REMOTE_PIPELINE


def test_write_reproducibility_environment_uses_minimum_versions(tmp_path):
    """environment.yml must pin schema minimums, not detected runtime versions."""
    from provenance import write_reproducibility_environment
    from schemas import JAVA_MIN_VERSION, NEXTFLOW_MIN_VERSION

    captured: dict = {}

    def fake_write_env(output_dir, *, env_name, pip_deps, conda_deps, python_version):
        captured["conda_deps"] = conda_deps

    with patch("provenance.write_environment_yml", fake_write_env):
        write_reproducibility_environment(
            tmp_path,
            preflight_result={
                "java": {"version": "999.0.0"},
                "nextflow": {"version": "999.9.9"},
            },
        )

    expected_java = f"openjdk>={JAVA_MIN_VERSION}"
    expected_nf = f"nextflow>={'.'.join(map(str, NEXTFLOW_MIN_VERSION))}"
    assert expected_java in captured["conda_deps"]
    assert expected_nf in captured["conda_deps"]
