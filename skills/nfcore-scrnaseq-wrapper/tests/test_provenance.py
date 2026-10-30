from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from unittest.mock import patch

from provenance import write_provenance_bundle


def _run_write_provenance_bundle(tmp_path):
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
    args = Namespace(
        demo=False,
        check=False,
        preset="star",
        profile="docker",
        pipeline_version="4.1.0",
        resume=False,
    )
    return write_provenance_bundle(
        output_dir,
        args=args,
        pipeline_source={
            "source_kind": "local_checkout",
            "source_ref": "scrnaseq-checkout",
            "resolved_version": "abc123",
            "branch": "main",
            "dirty": False,
        },
        preflight_result={
            "java": {"version": "17.0.8"},
            "nextflow": {"version": "25.04.0"},
            "references": {},
            "profile": {"profile": "docker"},
        },
        params_path=params_path,
        params_payload={"aligner": "star"},
        normalized_samplesheet=normalized,
        samplesheet_summary={
            "sample_count": 1,
            "fastq_paths": [],
            "unknown_columns": [],
        },
        parsed_outputs={
            "preferred_h5ad": str(preferred),
            "multiqc_report": "",
            "pipeline_info_dir": "",
            "h5ad_candidates": [str(preferred)],
            "rds_candidates": [],
            "handoff_available": True,
        },
        execution_result={"stdout_path": str(stdout), "stderr_path": str(stderr)},
        command_str="nextflow run ...",
    )


def test_write_provenance_bundle(tmp_path):
    provenance_dir, manifest_path = _run_write_provenance_bundle(tmp_path)
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
    assert runtime["work_dir"] == "upstream/work"
    assert manifest["work_dir"] == "upstream/work"
    assert manifest["java_version"] == "17.0.8"
    assert manifest["nextflow_version"] == "25.04.0"
    assert manifest["python_version"] == runtime["python_version"]
    assert manifest["environment_yml_mode"] == "install_recipe"


def test_bundle_is_consolidated_under_reproducibility(tmp_path):
    _run_write_provenance_bundle(tmp_path)

    repro = tmp_path / "reproducibility"
    # All provenance JSONs live in reproducibility/, not a separate provenance/.
    assert not (tmp_path / "provenance").exists()
    for name in (
        "runtime.json",
        "upstream.json",
        "inputs.json",
        "outputs.json",
        "skill.json",
        "manifest.json",
    ):
        assert (repro / name).is_file(), f"missing {name}"


def _make_args(tmp_path, **overrides):
    values = dict(
        demo=False,
        check=False,
        preset="star",
        profile="docker",
        pipeline_version="4.1.0",
        resume=False,
        allow_dirty_pipeline=False,
        extra_config=[],
    )
    values.update(overrides)
    return Namespace(**values)


def test_runtime_and_inputs_paths_are_relative(tmp_path):
    import provenance

    runtime = provenance.build_runtime_payload(
        tmp_path,
        args=_make_args(tmp_path),
        preflight_result={
            "java": {"version": "17"},
            "nextflow": {"version": "25.04.0"},
        },
        command_str="nextflow run ...",
        timestamp="2026-05-31T00:00:00+00:00",
    )
    assert runtime["work_dir"] == "upstream/work"
    assert runtime["cwd"] == "."

    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    samplesheet.parent.mkdir(parents=True, exist_ok=True)
    samplesheet.write_text("sample,fastq_1\n", encoding="utf-8")
    params_path = tmp_path / "reproducibility" / "params.yaml"
    params_path.write_text("outdir: upstream/results\n", encoding="utf-8")

    inputs = provenance.build_inputs_payload(
        normalized_samplesheet=samplesheet,
        samplesheet_summary={
            "sample_count": 1,
            "fastq_paths": ["/data/s1_R1.fastq.gz"],
        },
        preflight_result={"references": {}},
        params_path=params_path,
        output_dir=tmp_path,
    )
    assert inputs["samplesheet"] == "reproducibility/samplesheet.valid.csv"
    assert inputs["params_path"] == "reproducibility/params.yaml"
    assert inputs["fastq_paths"] == ["/data/s1_R1.fastq.gz"]


def test_runtime_records_custom_remote_work_dir(tmp_path):
    import provenance

    runtime = provenance.build_runtime_payload(
        tmp_path,
        args=_make_args(tmp_path, work_dir="s3://bucket/scrnaseq/work"),
        preflight_result={
            "java": {"version": "17"},
            "nextflow": {"version": "25.04.0"},
        },
        command_str="nextflow run ...",
        timestamp="2026-05-31T00:00:00+00:00",
    )

    assert runtime["work_dir"] == "s3://bucket/scrnaseq/work"


def test_invocation_records_dirty_opt_in_and_extra_configs(tmp_path):
    import provenance

    args = _make_args(
        tmp_path,
        allow_dirty_pipeline=True,
        require_local_pipeline=True,
        allow_conda_cellranger=True,
        work_dir="s3://bucket/scrnaseq/work",
        extra_config=[str(tmp_path / "cluster.config")],
    )
    invocation = provenance.build_invocation_payload(
        args, timestamp="2026-05-31T00:00:00+00:00"
    )

    assert invocation["allow_dirty_pipeline"] is True
    assert invocation["require_local_pipeline"] is True
    assert invocation["allow_conda_cellranger"] is True
    assert invocation["work_dir"] == "s3://bucket/scrnaseq/work"
    assert invocation["extra_config"] == [str(tmp_path / "cluster.config")]


def test_write_reproducibility_environment_pins_installable_versions(tmp_path):
    """environment.yml must be an installable recipe: nextflow pinned exactly
    (it lives on bioconda, so the bioconda channel is required), and openjdk
    pinned to its major series (conda-forge does not publish every JDK patch, so
    an exact patch pin can be unsatisfiable). No '>=' ranges."""
    from provenance import write_reproducibility_environment

    captured: dict = {}

    def fake_write_env(
        output_dir, *, env_name, pip_deps, conda_deps, python_version, channels
    ):
        captured["conda_deps"] = conda_deps
        captured["channels"] = channels

    with patch("provenance.write_environment_yml", fake_write_env):
        write_reproducibility_environment(
            tmp_path,
            preflight_result={
                "java": {"version": "17.0.10"},
                "nextflow": {"version": "999.9.9"},
            },
        )

    assert "openjdk=17" in captured["conda_deps"]  # major series, resolvable
    assert "openjdk=17.0.10" not in captured["conda_deps"]  # not the exact patch
    assert "nextflow=999.9.9" in captured["conda_deps"]  # engine pinned exactly
    assert "bioconda" in captured["channels"]  # nextflow lives there
    assert not any(">=" in dep for dep in captured["conda_deps"])


def test_policy_and_pinned_files_copied_into_bundle(tmp_path):
    _run_write_provenance_bundle(tmp_path)
    repro = tmp_path / "reproducibility"
    assert (repro / "compatibility_policy.json").is_file()
    assert (repro / "pinned_versions.json").is_file()


def test_environment_yml_is_installable_recipe(tmp_path):
    import provenance

    provenance.write_reproducibility_environment(
        tmp_path,
        preflight_result={
            "java": {"version": "17.0.10"},
            "nextflow": {"version": "25.04.2"},
        },
    )
    env = (tmp_path / "reproducibility" / "environment.yml").read_text(encoding="utf-8")
    # nextflow exact + bioconda channel (it is not on conda-forge); openjdk by
    # major series (resolvable on conda-forge); no open-ended ranges.
    assert "nextflow=25.04.2" in env
    assert "bioconda" in env
    assert "openjdk=17" in env
    assert "openjdk=17.0.10" not in env
    assert ">=" not in env


_FASTQC = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
_FASTQC_DIGEST = "sha256:" + "a" * 64


def test_bundle_records_pins_and_checksums_the_containers(tmp_path, monkeypatch):
    import subprocess

    from clawbio.common.nfcore import containers

    seen_env = []

    def run(command, *, cwd=None, env=None, timeout_seconds):
        if command[:2] == ["nextflow", "log"]:
            seen_env.append(env)
            return subprocess.CompletedProcess(command, 0, f"NFCORE_SCRNASEQ:SCRNASEQ:FASTQC\t{_FASTQC}\tCOMPLETED\n", "")
        if command[1:3] == ["image", "inspect"]:
            payload = {"Id": "sha256:" + "f" * 64, "RepoDigests": [f"quay.io/biocontainers/fastqc@{_FASTQC_DIGEST}"],
                       "Os": "linux", "Architecture": "amd64"}
            return subprocess.CompletedProcess(command, 0, json.dumps(payload), "")
        if command[1:2] == ["version"]:
            return subprocess.CompletedProcess(command, 0, "29.1.3 linux/amd64\n", "")
        raise AssertionError(command)

    monkeypatch.setattr(containers, "_run", run)
    (tmp_path / ".nextflow").mkdir()
    (tmp_path / ".nextflow" / "history").write_text("run\n", encoding="utf-8")
    bundle_dir, manifest_path = _run_write_provenance_bundle(tmp_path)
    record = json.loads((bundle_dir / "containers.json").read_text())
    assert record["images"][0]["pinned_reference"] == f"quay.io/biocontainers/fastqc@{_FASTQC_DIGEST}"
    assert seen_env == [{"NXF_SYNTAX_PARSER": "v1"}]
    assert json.loads(manifest_path.read_text())["containers"]["pin_config"] == "reproducibility/containers.config"
    labels = [line.split("  ", 1)[1] for line in (bundle_dir / "checksums.sha256").read_text().splitlines()]
    assert "reproducibility/containers.json" in labels and "reproducibility/containers.config" in labels
