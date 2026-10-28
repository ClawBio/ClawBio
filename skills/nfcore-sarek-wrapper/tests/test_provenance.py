"""Tests for nfcore-sarek-wrapper / provenance.py."""
from __future__ import annotations

import hashlib
import json
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any

import pytest

_SKILL_DIR = Path(__file__).resolve().parent.parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))


def _purge_foreign(*names: str) -> None:
    for name in names:
        module = sys.modules.get(name)
        module_file = Path(getattr(module, "__file__", "") or "")
        if (
            module is not None
            and _SKILL_DIR not in module_file.parents
            and module_file != _SKILL_DIR / f"{name}.py"
        ):
            sys.modules.pop(name, None)


_purge_foreign("provenance", "schemas", "preflight")

from clawbio.common.tests.nfcore_test_imports import load_tested_module  # noqa: E402

provenance = load_tested_module(_SKILL_DIR, "provenance")
ProvenanceBundle = provenance.ProvenanceBundle
compute_params_checksum = provenance.compute_params_checksum
compute_reference_checksums = provenance.compute_reference_checksums
load_manifest = provenance.load_manifest
write_provenance_bundle = provenance.write_provenance_bundle


# ---------------------------------------------------------------------------
# Helpers / fixtures
# ---------------------------------------------------------------------------


_REPRO_DIR = _SKILL_DIR / "reproducibility"


@dataclass
class _FakePipelineSource:
    source_kind: str = "remote_repo"
    resolved_version: str = "3.8.1"
    resolved_uri: str = "https://github.com/nf-core/sarek"
    local_path: str | None = None


@dataclass
class _FakeSamplesheetReport:
    sample_count: int = 2
    patient_count: int = 1
    analysis_mode: str = "germline"
    sample_names: list[str] = field(default_factory=lambda: ["S1", "S2"])
    fastq_paths: list[Path] = field(default_factory=list)


@dataclass
class _FakeOutputsReport:
    step_completed: str = "variant_calling"
    preprocessing: dict = field(default_factory=dict)
    variant_calling: dict = field(default_factory=dict)
    annotation: dict = field(default_factory=dict)
    qc: dict = field(default_factory=dict)
    pipeline_info: Path | None = None
    reference_outputs: dict | None = None
    samples_detected: list[str] = field(default_factory=list)
    pairs_detected: list[str] = field(default_factory=list)
    csv_handoff: dict = field(default_factory=dict)
    handoff_available: bool = True
    warnings: list[str] = field(default_factory=list)
    missing_outputs: list[str] = field(default_factory=list)


def _write_caller_artifacts(tmp_path: Path) -> tuple[Path, Path, Path]:
    samplesheet = tmp_path / "samplesheet.normalized.csv"
    samplesheet.write_text(
        "patient,sample,lane,fastq_1\nP1,S1,L1,/data/S1_R1.fastq.gz\n",
        encoding="utf-8",
    )
    params_yaml = tmp_path / "params.yaml"
    params_yaml.write_text("input: samplesheet.csv\nstep: mapping\n", encoding="utf-8")
    commands_sh = tmp_path / "commands.sh"
    commands_sh.write_text("#!/usr/bin/env bash\nnextflow run nf-core/sarek -r 3.8.1\n", encoding="utf-8")
    return samplesheet, params_yaml, commands_sh


def _base_params(**overrides) -> dict:
    params = {
        "step": "mapping",
        "aligner": "bwa-mem",
        "profile": "docker",
        "tools": ["strelka", "manta"],
        "skip_tools": [],
        "joint_germline": False,
        "joint_mutect2": False,
        "wes": False,
    }
    params.update(overrides)
    return params


def _run_bundle(tmp_path: Path, **overrides) -> ProvenanceBundle:
    samplesheet_csv, params_yaml, commands_sh = _write_caller_artifacts(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    args: dict[str, Any] = dict(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        samplesheet_csv_src=samplesheet_csv,
        params_yaml_src=params_yaml,
        commands_sh_src=commands_sh,
        params=_base_params(),
        samplesheet_report=_FakeSamplesheetReport(),
        pipeline_source=_FakePipelineSource(),
        outputs_report=_FakeOutputsReport(),
        nextflow_version="25.10.2",
        java_version="17.0.10",
        profile="docker",
    )
    args.update(overrides)
    return write_provenance_bundle(**args)


# ---------------------------------------------------------------------------
# compute_params_checksum
# ---------------------------------------------------------------------------



def test_params_checksum_order_independent():
    a = compute_params_checksum({"a": 1, "b": 2, "c": 3})
    b = compute_params_checksum({"c": 3, "b": 2, "a": 1})
    assert a == b


def test_params_checksum_changes_on_value_change():
    a = compute_params_checksum({"step": "mapping"})
    b = compute_params_checksum({"step": "annotate"})
    assert a != b


# ---------------------------------------------------------------------------
# compute_samplesheet_checksum
# ---------------------------------------------------------------------------


# ---------------------------------------------------------------------------
# compute_reference_checksums
# ---------------------------------------------------------------------------


def test_reference_checksums_uri_passthrough():
    params = {"fasta": "s3://bucket/hg38.fa", "dbsnp": "https://example.com/dbsnp.vcf.gz"}
    out = compute_reference_checksums(params)
    assert out["fasta"] == "s3://bucket/hg38.fa"
    assert out["dbsnp"] == "https://example.com/dbsnp.vcf.gz"





def test_reference_checksums_only_recognised_keys(tmp_path: Path):
    fasta = tmp_path / "hg38.fa"
    fasta.write_bytes(b"x")
    params = {"fasta": str(fasta), "unrelated_key": str(fasta)}
    out = compute_reference_checksums(params)
    assert "fasta" in out
    assert "unrelated_key" not in out


# ---------------------------------------------------------------------------
# write_provenance_bundle — happy path
# ---------------------------------------------------------------------------


def test_bundle_writes_all_expected_files(tmp_path: Path):
    bundle = _run_bundle(tmp_path)
    names = {p.name for p in bundle.files_written}
    expected = {
        "samplesheet.valid.csv",
        "params.yaml",
        "commands.sh",
        "manifest.json",
        "checksums.sha256",
        "environment.yml",
        "pipeline_source.json",
        "parameters.json",
        "samplesheet.json",
        "outputs.json",
        "tool_versions.json",
        "compatibility_policy.json",
    }
    assert expected.issubset(names), f"Missing from bundle: {expected - names}"


def test_provenance_bundle_files_use_lf_line_endings(tmp_path: Path):
    """Cross-OS portability invariant: every provenance-bundle text file
    (checksums.sha256, manifest.json, environment.yml, the copied params.yaml /
    commands.sh / samplesheet, and all JSON) is LF-only. CRLF would break the
    bash commands.sh on Windows and make checksums.sha256 OS-dependent."""
    bundle = _run_bundle(tmp_path)
    text_suffixes = {".sh", ".yaml", ".yml", ".json", ".md", ".csv", ".sha256"}
    checked = 0
    for path in bundle.files_written:
        p = Path(path)
        if p.is_file() and (p.suffix in text_suffixes):
            data = p.read_bytes()
            assert b"\r" not in data, f"{p.name} contains a CR byte (not LF-only)"
            checked += 1
    assert checked >= 6  # csv + yaml + sh + sha256 + json + yml


def _run_bundle_with_fake_output(tmp_path: Path) -> ProvenanceBundle:
    """Build a full bundle whose output dir holds a real pipeline output, so the
    output-scoped checksums.sha256 is populated (sarek checksums the results, not
    the bundle; input integrity lives in the manifest's params/samplesheet hashes)."""
    samplesheet_csv, params_yaml, commands_sh = _write_caller_artifacts(tmp_path)
    output_dir = tmp_path / "out"
    results = output_dir / "upstream" / "results" / "variant_calling"
    results.mkdir(parents=True)
    (results / "sample.vcf.gz").write_bytes(b"fake compressed vcf\n")
    return write_provenance_bundle(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        samplesheet_csv_src=samplesheet_csv,
        params_yaml_src=params_yaml,
        commands_sh_src=commands_sh,
        params=_base_params(),
        samplesheet_report=_FakeSamplesheetReport(),
        pipeline_source=_FakePipelineSource(),
        outputs_report=_FakeOutputsReport(),
        nextflow_version="25.10.2",
        java_version="17.0.10",
        profile="docker",
    )


def test_checksums_self_verifiable_in_place(tmp_path: Path):
    """Every label in checksums.sha256 resolves relative to the output dir and
    matches — the manifest is self-verifiable in place on any OS. After the LF
    fixes the digests are byte-stable across operating systems."""
    bundle = _run_bundle_with_fake_output(tmp_path)
    repro = Path(bundle.manifest_path).parent
    output_dir = repro.parent
    raw = (repro / "checksums.sha256").read_bytes()
    assert b"\r" not in raw, "checksums.sha256 must be LF-only"
    lines = [
        ln for ln in raw.decode("utf-8").splitlines()
        if ln.strip() and not ln.startswith("#")
    ]
    assert lines, "checksums.sha256 is empty"
    failures = []
    for ln in lines:
        digest, label = ln.split("  ", 1)
        target = output_dir / label
        if not target.exists():
            failures.append((label, "MISSING"))
        elif hashlib.sha256(target.read_bytes()).hexdigest() != digest:
            failures.append((label, "MISMATCH"))
    assert failures == [], f"checksums.sha256 not verifiable in place: {failures}"


def test_sha256sum_c_passes_from_output_dir(tmp_path: Path):
    """If a sha256 CLI is available, `-c` returns 0 from the output dir."""
    tool = shutil.which("sha256sum")
    if not tool:
        pytest.skip("sha256sum not available")
    bundle = _run_bundle_with_fake_output(tmp_path)
    repro = Path(bundle.manifest_path).parent
    r = subprocess.run(
        [tool, "-c", "reproducibility/checksums.sha256"],
        cwd=str(repro.parent),
        capture_output=True,
        text=True,
    )
    assert r.returncode == 0, f"sha256sum -c failed:\n{r.stdout}\n{r.stderr}"


def test_manifest_has_required_fields(tmp_path: Path):
    bundle = _run_bundle(tmp_path)
    manifest = json.loads(bundle.manifest_path.read_text())
    for key in (
        "schema_version", "skill_name", "skill_version", "generated_at",
        "pipeline_source", "step", "aligner", "tools", "skip_tools",
        "analysis_mode", "joint_germline", "joint_mutect2", "wes",
        "intervals_hash", "profile", "arm", "spark", "gpu", "resume_used",
        "java_version", "nextflow_version", "params_checksum",
        "samplesheet_checksum", "reference_checksums",
    ):
        assert key in manifest, f"manifest is missing key '{key}'"
    assert manifest["params_checksum"].startswith("sha256:")
    assert manifest["samplesheet_checksum"].startswith("sha256:")


def test_manifest_tools_sorted_and_deduped(tmp_path: Path):
    bundle = _run_bundle(
        tmp_path,
        params=_base_params(tools=["manta", "strelka", "manta", "haplotypecaller"]),
    )
    manifest = json.loads(bundle.manifest_path.read_text())
    assert manifest["tools"] == ["haplotypecaller", "manta", "strelka"]



def test_environment_yml_contains_pinned_nextflow(tmp_path: Path):
    bundle = _run_bundle(tmp_path, nextflow_version="25.10.2")
    env_yml = (bundle.bundle_dir / "environment.yml").read_text()
    assert "name: claw-sarek" in env_yml
    assert "nextflow=25.10.2" in env_yml
    assert "conda-forge" in env_yml
    assert "bioconda" in env_yml


# ---------------------------------------------------------------------------
# outputs.json behaviour
# ---------------------------------------------------------------------------



def test_outputs_json_relative_paths(tmp_path: Path):
    samplesheet_csv, params_yaml, commands_sh = _write_caller_artifacts(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    inside_file = output_dir / "preprocessing" / "csv" / "mapped.csv"
    inside_file.parent.mkdir(parents=True)
    inside_file.write_text("sample\nS1\n")
    report = _FakeOutputsReport(
        csv_handoff={"mapped": inside_file},
        samples_detected=["S1"],
    )
    bundle = write_provenance_bundle(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        samplesheet_csv_src=samplesheet_csv,
        params_yaml_src=params_yaml,
        commands_sh_src=commands_sh,
        params=_base_params(),
        samplesheet_report=_FakeSamplesheetReport(),
        pipeline_source=_FakePipelineSource(),
        outputs_report=report,
        nextflow_version="25.10.2",
    )
    outputs = json.loads((bundle.bundle_dir / "outputs.json").read_text())
    assert outputs["csv_handoff"]["mapped"] == "preprocessing/csv/mapped.csv"


# ---------------------------------------------------------------------------
# pipeline_source / samplesheet / parameters JSON
# ---------------------------------------------------------------------------



def test_pipeline_source_json_handles_dict(tmp_path: Path):
    bundle = _run_bundle(
        tmp_path,
        pipeline_source={
            "source_kind": "local_checkout",
            "resolved_version": "abc123",
            "resolved_uri": "",
            "local_path": "/tmp/sarek",
        },
    )
    payload = json.loads((bundle.bundle_dir / "pipeline_source.json").read_text())
    assert payload["source_kind"] == "local_checkout"
    assert payload["local_path"] == "/tmp/sarek"


def test_parameters_json_is_serialisable(tmp_path: Path):
    bundle = _run_bundle(
        tmp_path,
        params=_base_params(fasta=Path("/refs/hg38.fa")),
    )
    payload = json.loads((bundle.bundle_dir / "parameters.json").read_text())
    assert payload["fasta"] == "/refs/hg38.fa"
    assert payload["step"] == "mapping"


# ---------------------------------------------------------------------------
# tool_versions.json
# ---------------------------------------------------------------------------



def test_tool_versions_parsed_from_yml(tmp_path: Path):
    pytest.importorskip("yaml")
    samplesheet_csv, params_yaml, commands_sh = _write_caller_artifacts(tmp_path)
    output_dir = tmp_path / "out"
    (output_dir / "upstream" / "results" / "pipeline_info").mkdir(parents=True)
    (output_dir / "upstream" / "results" / "pipeline_info" / "software_versions.yml").write_text(
        "Workflow:\n  nf-core/sarek: 3.8.1\nbwa:\n  bwa: 0.7.17\n",
        encoding="utf-8",
    )
    bundle = write_provenance_bundle(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        samplesheet_csv_src=samplesheet_csv,
        params_yaml_src=params_yaml,
        commands_sh_src=commands_sh,
        params=_base_params(),
        samplesheet_report=_FakeSamplesheetReport(),
        pipeline_source=_FakePipelineSource(),
        outputs_report=_FakeOutputsReport(),
        nextflow_version="25.10.2",
    )
    payload = json.loads((bundle.bundle_dir / "tool_versions.json").read_text())
    assert payload["Workflow"]["nf-core/sarek"] == "3.8.1"
    assert payload["bwa"]["bwa"] == "0.7.17"


# ---------------------------------------------------------------------------
# checksums.sha256
# ---------------------------------------------------------------------------


def test_checksums_excludes_work_and_nextflow(tmp_path: Path):
    samplesheet_csv, params_yaml, commands_sh = _write_caller_artifacts(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    # Real output
    (output_dir / "preprocessing").mkdir()
    real = output_dir / "preprocessing" / "result.txt"
    real.write_text("hello", encoding="utf-8")
    # Things we MUST exclude
    (output_dir / "work" / "abc").mkdir(parents=True)
    (output_dir / "work" / "abc" / "junk.txt").write_text("junk")
    (output_dir / ".nextflow" / "history").mkdir(parents=True)
    (output_dir / ".nextflow" / "history" / "stuff").write_text("nope")
    (output_dir / "run.log").write_text("log content")
    # Execution logs (stdout/stderr capture) must never enter the checksums.
    (output_dir / "logs").mkdir()
    (output_dir / "logs" / "stdout.txt").write_text("run stdout")
    # A leftover --check summary at the root (check → real run in same dir) is a
    # wrapper artifact, not a pipeline output, so it stays out of the manifest.
    (output_dir / "check_result.json").write_text("{}")

    bundle = write_provenance_bundle(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        samplesheet_csv_src=samplesheet_csv,
        params_yaml_src=params_yaml,
        commands_sh_src=commands_sh,
        params=_base_params(),
        samplesheet_report=_FakeSamplesheetReport(),
        pipeline_source=_FakePipelineSource(),
        outputs_report=_FakeOutputsReport(),
        nextflow_version="25.10.2",
    )
    text = (bundle.bundle_dir / "checksums.sha256").read_text()
    assert "preprocessing/result.txt" in text
    assert "work/" not in text
    assert ".nextflow" not in text
    assert "run.log" not in text
    assert "logs/" not in text
    assert "check_result.json" not in text


def test_checksums_excludes_nested_upstream_work(tmp_path: Path):
    # Real layout: the pipeline writes results under upstream/results/ and its
    # ephemeral Nextflow scratch under upstream/work/ (the default work dir is
    # <output>/upstream/work). The work tree must never enter the manifest even
    # though it is nested one level below the output root.
    samplesheet_csv, params_yaml, commands_sh = _write_caller_artifacts(tmp_path)
    output_dir = tmp_path / "out"
    output_dir.mkdir()
    # Real pipeline output under upstream/results/.
    real = output_dir / "upstream" / "results" / "preprocessing" / "sample.cram"
    real.parent.mkdir(parents=True)
    real.write_text("cram", encoding="utf-8")
    # Nextflow scratch under upstream/work/ — ephemeral, must be excluded.
    junk = output_dir / "upstream" / "work" / "ab" / "cd1234" / "task.bam"
    junk.parent.mkdir(parents=True)
    junk.write_text("ephemeral", encoding="utf-8")

    bundle = write_provenance_bundle(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        samplesheet_csv_src=samplesheet_csv,
        params_yaml_src=params_yaml,
        commands_sh_src=commands_sh,
        params=_base_params(),
        samplesheet_report=_FakeSamplesheetReport(),
        pipeline_source=_FakePipelineSource(),
        outputs_report=_FakeOutputsReport(),
        nextflow_version="25.10.2",
    )
    text = (bundle.bundle_dir / "checksums.sha256").read_text()
    assert "upstream/results/preprocessing/sample.cram" in text
    assert "upstream/work/" not in text
    assert "task.bam" not in text



# ---------------------------------------------------------------------------
# load_manifest
# ---------------------------------------------------------------------------


def test_load_manifest_round_trip(tmp_path: Path):
    bundle = _run_bundle(tmp_path)
    manifest = load_manifest(bundle.bundle_dir)
    assert manifest is not None
    assert manifest["skill_name"] == "nfcore-sarek-wrapper"
    assert manifest["step"] == "mapping"



# ---------------------------------------------------------------------------
# Pinned versions parity test
# ---------------------------------------------------------------------------


def test_pinned_versions_parity_with_schemas():
    """`pinned_versions.json` must mirror schemas and official Sarek 3.8.1 enums."""
    import schemas

    pinned_path = _REPRO_DIR / "pinned_versions.json"
    assert pinned_path.is_file(), "pinned_versions.json is required next to the skill"
    pinned = json.loads(pinned_path.read_text())
    assert pinned["skill"] == schemas.SKILL_NAME
    assert pinned["skill_version"] == schemas.SKILL_VERSION
    assert pinned["pipeline"]["name"] == schemas.DEFAULT_REMOTE_PIPELINE
    assert pinned["pipeline"]["version"] == schemas.DEFAULT_PIPELINE_VERSION

    assert set(pinned["supported_profiles"]) == set(schemas.SUPPORTED_PROFILES)
    assert set(pinned["supported_aligners"]) == set(schemas.SUPPORTED_ALIGNERS)
    assert set(pinned["supported_steps"]) == set(schemas.SUPPORTED_STEPS)
    assert set(pinned["supported_tools"]) == set(schemas.SUPPORTED_TOOLS)
    assert set(pinned["supported_skip_tools"]) == set(schemas.SUPPORTED_SKIP_TOOLS)
    assert set(pinned["supported_sex"]) == set(schemas.SUPPORTED_SEX)
    assert set(pinned["supported_status"]) == set(schemas.SUPPORTED_STATUS)
    assert set(pinned["igenomes_catalogue"]) == set(schemas.SUPPORTED_IGENOMES_NAMES)

    official_profiles = {
        "debug", "conda", "mamba", "docker", "arm64", "emulate_amd64",
        "singularity", "podman", "shifter", "charliecloud", "apptainer",
        "wave", "gpu", "spark", "test", "test_aws", "test_azure",
        "test_full", "test_full_aws", "test_full_azure",
        "test_full_germline", "test_full_germline_aws",
        "test_full_germline_azure", "test_full_germline_ncbench_agilent",
        "mutect",
    }
    official_tools = {
        "ascat", "bbsplit", "bcfann", "cnvkit", "controlfreec",
        "deepvariant", "freebayes", "haplotypecaller", "indexcov",
        "lofreq", "manta", "merge", "mpileup", "msisensor2",
        "msisensorpro", "muse", "mutect2", "ngscheckmate",
        "sentieon_dedup", "sentieon_dnascope", "sentieon_haplotyper",
        "sentieon_tnscope", "snpeff", "snpsift", "strelka", "tiddit",
        "varlociraptor", "vep",
    }
    official_skip_tools = {
        "baserecalibrator", "baserecalibrator_report", "bcftools",
        "dnascope_filter", "documentation", "fastqc", "haplotypecaller_filter",
        "haplotyper_filter", "markduplicates", "markduplicates_report",
        "mosdepth", "multiqc", "samtools", "vcftools", "versions",
    }
    assert set(pinned["supported_profiles"]) == official_profiles
    assert set(pinned["supported_tools"]) == official_tools
    assert set(pinned["supported_skip_tools"]) == official_skip_tools


# ---------------------------------------------------------------------------
# compatibility_policy.json shape
# ---------------------------------------------------------------------------


def test_compatibility_policy_loads_with_expected_fields():
    policy_path = _REPRO_DIR / "compatibility_policy.json"
    assert policy_path.is_file()
    policy = json.loads(policy_path.read_text())
    # These are the fields actually tracked for resume drift, matching
    # _RESUME_TRACKED_FIELDS in preflight.py plus set-valued and pipeline_source
    # fields compared separately in _check_resume_drift.
    expected_drift_fields = {
        "step", "aligner", "tools", "skip_tools", "analysis_mode",
        "joint_germline", "joint_mutect2", "wes", "profile",
        "pipeline_source.source_kind", "pipeline_source.resolved_version",
        "params_checksum", "reference_checksums", "samplesheet_checksum",
    }
    # Exact equality (not issubset): the policy JSON and the code's tracked-field
    # set are the single source of truth for resume drift and must not diverge in
    # either direction — a policy field the code ignored would silently fail to
    # gate resume.
    assert expected_drift_fields == set(policy["resume_drift_fields"])
    # arm/gpu/spark are profile modifiers, not separate drift fields.
    assert {"arm", "gpu", "spark"}.issubset(set(policy["ignored_during_drift"]))
    assert policy["schema_version"] == 1
    assert "legacy_field_defaults" in policy
    assert "ignored_during_drift" in policy


# ---------------------------------------------------------------------------
# REFERENCE_PATH_PARAMS parity with preflight
# ---------------------------------------------------------------------------


def test_reference_path_params_parity_with_preflight():
    import preflight

    assert provenance.REFERENCE_PATH_PARAMS == preflight.REFERENCE_PATH_PARAMS


# ---------------------------------------------------------------------------
# Pre-run (outputs_report=None) full smoke
# ---------------------------------------------------------------------------


def test_pre_run_bundle_still_writes_core_files(tmp_path: Path):
    bundle = _run_bundle(tmp_path, outputs_report=None)
    must_exist = (
        "samplesheet.valid.csv",
        "params.yaml",
        "commands.sh",
        "manifest.json",
        "checksums.sha256",
        "environment.yml",
        "pipeline_source.json",
        "parameters.json",
        "samplesheet.json",
        "tool_versions.json",
        "compatibility_policy.json",
    )
    for name in must_exist:
        assert (bundle.bundle_dir / name).is_file(), f"Missing {name} in pre-run bundle"
    assert not (bundle.bundle_dir / "outputs.json").exists()


def test_intervals_hash_recorded(tmp_path: Path):
    intervals = tmp_path / "intervals.bed"
    intervals.write_text("chr1\t1\t100\n")
    bundle = _run_bundle(tmp_path, params=_base_params(intervals=str(intervals), wes=True))
    manifest = json.loads(bundle.manifest_path.read_text())
    assert manifest["intervals_hash"] is not None
    assert manifest["intervals_hash"].startswith("sha256:")
    assert manifest["wes"] is True


# ---------------------------------------------------------------------------
# containers.json / containers.config
# ---------------------------------------------------------------------------

_FASTQC = "quay.io/biocontainers/fastqc:0.12.1--hdfd78af_0"
_FASTQC_DIGEST = "sha256:" + "a" * 64


def _fake_container_engine(seen_env: list):
    def run(command, *, cwd=None, env=None, timeout_seconds):
        if command[:2] == ["nextflow", "log"]:
            seen_env.append(env)
            return subprocess.CompletedProcess(command, 0, f"NFCORE_SAREK:SAREK:FASTQC\t{_FASTQC}\tCOMPLETED\n", "")
        if command[1:3] == ["image", "inspect"]:
            payload = {"Id": "sha256:" + "f" * 64, "RepoDigests": [f"quay.io/biocontainers/fastqc@{_FASTQC_DIGEST}"],
                       "Os": "linux", "Architecture": "amd64"}
            return subprocess.CompletedProcess(command, 0, json.dumps(payload), "")
        if command[1:2] == ["version"]:
            return subprocess.CompletedProcess(command, 0, "29.1.3 linux/amd64\n", "")
        raise AssertionError(command)
    return run


def _bundle_after_a_run(tmp_path: Path) -> ProvenanceBundle:
    samplesheet_csv, params_yaml, commands_sh = _write_caller_artifacts(tmp_path)
    output_dir = tmp_path / "out"
    (output_dir / ".nextflow").mkdir(parents=True)
    (output_dir / ".nextflow" / "history").write_text("run\n")
    bundle_args: dict[str, Any] = {
        "output_dir": output_dir, "skill_dir": _SKILL_DIR, "samplesheet_csv_src": samplesheet_csv,
        "params_yaml_src": params_yaml, "commands_sh_src": commands_sh, "params": _base_params(),
        "samplesheet_report": _FakeSamplesheetReport(), "pipeline_source": _FakePipelineSource(),
        "outputs_report": _FakeOutputsReport(), "nextflow_version": "26.04.3", "java_version": "17.0.10",
        "profile": "test,docker",
    }
    return write_provenance_bundle(**bundle_args)


def test_bundle_records_and_pins_the_containers_the_run_used(tmp_path: Path, monkeypatch):
    from clawbio.common.nfcore import containers

    seen_env: list = []
    monkeypatch.setattr(containers, "_run", _fake_container_engine(seen_env))
    bundle = _bundle_after_a_run(tmp_path)
    record = json.loads((bundle.bundle_dir / "containers.json").read_text())
    (image,) = record["images"]
    assert image["pinned_reference"] == f"quay.io/biocontainers/fastqc@{_FASTQC_DIGEST}"
    assert image["processes"] == ["NFCORE_SAREK:SAREK:FASTQC"] and record["engine"] == "docker"
    assert "withName: 'NFCORE_SAREK:SAREK:FASTQC'" in (bundle.bundle_dir / "containers.config").read_text()
    # The task log is read with the environment the pinned release runs under.
    assert seen_env == [{"NXF_SYNTAX_PARSER": "v1"}]
    manifest = json.loads(bundle.manifest_path.read_text())
    assert manifest["containers"]["complete"] is True
    assert manifest["containers"]["pin_config"] == "reproducibility/containers.config"
    assert bundle.bundle_dir / "containers.json" in bundle.files_written


def test_bundle_without_a_nextflow_run_says_its_containers_are_unrecorded(tmp_path: Path):
    bundle = _run_bundle(tmp_path)
    record = json.loads((bundle.bundle_dir / "containers.json").read_text())
    assert record["complete"] is False and not (bundle.bundle_dir / "containers.config").exists()
    assert any("containers" in warning for warning in bundle.warnings)
