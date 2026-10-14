"""Shared regression coverage for the three pinned nf-core wrappers."""

from __future__ import annotations

import argparse
import csv
import importlib.util
import io
import sys
from pathlib import Path

import pytest
from hypothesis import HealthCheck, given, settings, strategies as st

ROOT = Path(__file__).resolve().parents[1]


@pytest.fixture(params=["rnaseq", "scrnaseq", "sarek"])
def wrapper(request):
    previous_modules, previous_path = sys.modules.copy(), sys.path[:]
    directory = ROOT / "skills" / f"nfcore-{request.param}-wrapper"
    spec = importlib.util.spec_from_file_location(
        f"audit_{request.param}", directory / f"nfcore_{request.param}_wrapper.py"
    )
    module = importlib.util.module_from_spec(spec)
    try:
        spec.loader.exec_module(module)
        companions = {
            name: sys.modules[name]
            for name in (
                "executor",
                "samplesheet_builder",
                "params_builder",
                "pipeline_source",
                "errors",
            )
        }
        yield request.param, module, companions
    finally:
        for name in set(sys.modules) - previous_modules.keys():
            sys.modules.pop(name, None)
        sys.modules.update(previous_modules)
        sys.path[:] = previous_path


@pytest.mark.parametrize("value", ["nan", "inf", "-inf"])
def test_timeout_rejects_nonfinite(wrapper, value):
    _, module, _ = wrapper
    with pytest.raises(argparse.ArgumentTypeError):
        module._timeout_hours(value)


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.floats(min_value=0, max_value=1000, allow_nan=False, allow_infinity=False))
def test_timeout_round_trips(wrapper, value):
    _, module, _ = wrapper
    assert module._timeout_hours(str(value)) == value


def test_interrupt_terminates_nextflow(wrapper, tmp_path, monkeypatch):
    _, _, companions = wrapper
    executor = companions["executor"]

    class Process:
        def wait(self, timeout=None):
            raise KeyboardInterrupt

    process, terminated = Process(), []
    monkeypatch.setattr(executor.subprocess, "Popen", lambda *a, **kw: process)
    monkeypatch.setattr(executor, "_terminate_process_tree", terminated.append)
    with pytest.raises(KeyboardInterrupt):
        executor.execute_nextflow(
            ["nextflow"], cwd=tmp_path, output_dir=tmp_path, timeout_seconds=None
        )
    assert terminated == [process]


def test_log_tail_has_bounded_reads(wrapper, tmp_path, monkeypatch):
    _, _, companions = wrapper
    path = tmp_path / "large.log"
    path.write_bytes(b"x" * 100000 + b"final error\n")
    original_open, reads = Path.open, []

    class Reader:
        def __init__(self, handle):
            self.handle = handle

        def __enter__(self):
            return self

        def __exit__(self, *args):
            self.handle.close()

        def __getattr__(self, name):
            return getattr(self.handle, name)

        def read(self, size=-1):
            reads.append(size)
            return self.handle.read(size)

    monkeypatch.setattr(
        Path, "open", lambda p, *a, **kw: Reader(original_open(p, *a, **kw))
    )
    assert (
        companions["executor"]._read_log_tail(path, limit=64).endswith("final error\n")
    )
    assert reads and all(0 <= size <= 64 for size in reads)


def _sheet(name):
    if name == "rnaseq":
        return ["sample", "fastq_1", "strandedness"], [
            "S1",
            "s3://b/R1.fastq.gz",
            "auto",
        ]
    if name == "scrnaseq":
        return ["sample", "fastq_1", "fastq_2"], [
            "S1",
            "s3://b/R1.fastq.gz",
            "s3://b/R2.fastq.gz",
        ]
    return ["patient", "sample", "lane", "fastq_1", "fastq_2"], [
        "P1",
        "S1",
        "1",
        "s3://b/R1.fastq.gz",
        "s3://b/R2.fastq.gz",
    ]


@pytest.mark.parametrize(
    "defect",
    [
        "duplicate",
        "empty_header",
        "extra_value",
        "missing_value",
        "unclosed_quote",
        "directory",
        "encoding",
    ],
)
def test_bad_csv_fails_without_output(wrapper, tmp_path, defect):
    name, _, companions = wrapper
    header, row = _sheet(name)
    if defect == "duplicate":
        header.append("sample")
        row.append("S2")
    if defect == "empty_header":
        header.append("")
        row.append("ignored")
    if defect == "extra_value":
        row.append("ignored")
    if defect == "missing_value":
        row.pop()
    buffer = io.StringIO()
    writer = csv.writer(buffer)
    writer.writerow(header)
    writer.writerow(row)
    content = buffer.getvalue()
    if defect == "unclosed_quote":
        content = ",".join(header) + '\n"' + ",".join(row)
    source, dest = tmp_path / "samples.csv", tmp_path / "out.csv"
    if defect == "directory":
        source.mkdir()
    elif defect == "encoding":
        source.write_bytes(b"\xff\xfe")
    else:
        source.write_text(content)
    with pytest.raises(companions["errors"].SkillError) as exc:
        companions["samplesheet_builder"].validate_and_normalize_samplesheet(
            source, dest
        )
    assert exc.value.error_code in {"INVALID_SAMPLESHEET", "MISSING_INPUT"}
    assert not dest.exists()


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.integers(min_value=0, max_value=5))
def test_duplicate_header_any_position(wrapper, tmp_path, position):
    name, _, companions = wrapper
    header, row = _sheet(name)
    index = position % len(header)
    header.insert(index, "sample")
    row.insert(index, "other")
    source = tmp_path / "duplicates.csv"
    with source.open("w", newline="") as handle:
        writer = csv.writer(handle)
        writer.writerow(header)
        writer.writerow(row)
    with pytest.raises(companions["errors"].SkillError):
        companions["samplesheet_builder"].validate_and_normalize_samplesheet(
            source, tmp_path / "out.csv"
        )


@pytest.mark.parametrize("wrapper", ["scrnaseq"], indirect=True)
def test_multiqc_uri_preserved(wrapper, tmp_path):
    _, module, companions = wrapper
    uri = "https://example.org/multiqc.yaml"
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--multiqc-config", uri]
    )
    params = companions["params_builder"].build_effective_params(
        args, normalized_samplesheet=tmp_path / "samples.csv", output_dir=tmp_path
    )
    assert params["multiqc_config"] == uri


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_local_source_version_checked(wrapper, tmp_path, monkeypatch):
    _, module, companions = wrapper
    args = module.build_parser().parse_args(
        ["--demo", "--output", str(tmp_path / "run")]
    )
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: {
            "source_kind": "local_checkout",
            "source_ref": "/local/sarek",
            "resolved_version": "commit",
            "manifest_version": "3.4.0",
            "dirty": False,
        },
    )

    def not_reached(*args):
        pytest.fail("wrong pipeline reached input preparation")

    monkeypatch.setattr(module, "_prepare_samplesheet", not_reached)
    with pytest.raises(companions["errors"].SkillError) as exc:
        module._run_wrapper(args, tmp_path / "run")
    assert exc.value.error_code == "PIPELINE_SOURCE_INVALID"


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
@pytest.mark.parametrize(
    "content",
    [
        "tools: [unterminated",
        "wes: true\nwes: false\n",
        'wes: "false"\n',
        "clip_r1: .nan\n",
        "clip_r1: -1\n",
    ],
)
def test_invalid_params_file(wrapper, tmp_path, content):
    _, _, companions = wrapper
    source = tmp_path / "params.yaml"
    source.write_text(content)
    with pytest.raises(companions["errors"].SkillError) as exc:
        companions["params_builder"]._load_user_params_file(str(source))
    assert exc.value.error_code == "INVALID_FLAG_COMBINATION"


@pytest.mark.parametrize("available", [1, 2, 4, 7])
def test_memory_limit_fits_runtime(wrapper, monkeypatch, available):
    name, module, _ = wrapper
    monkeypatch.setattr(
        module,
        "_detect_host_memory_gb" if name == "rnaseq" else "_host_memory_gb",
        lambda: 64,
    )
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: available)
    assert 1 <= module._resource_limits_memory_gb() <= available


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_macos_vm_memory_and_timeout(wrapper, tmp_path, monkeypatch):
    _, module, _ = wrapper
    monkeypatch.setattr(module.sys, "platform", "darwin")
    monkeypatch.setattr(module, "_host_memory_gb", lambda: 64)
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: 2)
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--timeout-hours", "48"]
    )
    text = module._write_macos_docker_config(tmp_path, args=args).read_text()
    assert "memory: '1.GB'" in text and "time: '48.h'" in text


@pytest.mark.parametrize("wrapper", ["rnaseq", "sarek"], indirect=True)
def test_manifest_ignores_comments_and_quoted_braces(wrapper, tmp_path):
    _, _, companions = wrapper
    source = tmp_path / "pipeline"
    source.mkdir()
    (source / "assets").mkdir()
    (source / "assets/schema_input.json").write_text("{}")
    (source / "main.nf").write_text("// main")
    (source / "nextflow.config").write_text(
        '/*\nmanifest { version = "0.0.0" }\n*/\nmanifest {\n description = "} https://example.org"\n // version = "1.0.0"\n version = "3.8.1"\n}\n'
    )
    result = companions["pipeline_source"].resolve_pipeline_source(
        requested_version="3.8.1", local_pipeline_dir=source
    )
    assert result["manifest_version"] == "3.8.1"


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
@pytest.mark.parametrize(
    "suffix,content",
    [
        (".json", "[{broken]"),
        (".yaml", "samples: [broken"),
        (".yml", "- sample: S1\n  sample: S2\n"),
    ],
)
def test_bad_structured_sheet(wrapper, tmp_path, suffix, content):
    _, _, companions = wrapper
    source = tmp_path / f"samples{suffix}"
    source.write_text(content)
    with pytest.raises(companions["errors"].SkillError) as exc:
        companions["samplesheet_builder"].validate_and_normalize_samplesheet(
            source, tmp_path / "out.csv"
        )
    assert exc.value.error_code == "INVALID_SAMPLESHEET"


def test_failure_retains_resume_manifest(wrapper, tmp_path, monkeypatch):
    name, module, companions = wrapper
    args = module.build_parser().parse_args(
        ["--demo", "--output", str(tmp_path / "run")]
    )
    source = {
        "source_kind": "remote_repo",
        "source_ref": f"nf-core/{name}",
        "resolved_version": {"rnaseq": "3.26.0", "scrnaseq": "4.1.0", "sarek": "3.8.1"}[
            name
        ],
        "dirty": False,
        "branch": "",
    }
    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: source)
    preflight = {
        "java": {"version": "21.0.7"},
        "nextflow": {"version": "25.10.4"},
        "references": {},
        "warnings": [],
    }
    if name == "sarek":
        from types import SimpleNamespace

        preflight = SimpleNamespace(
            warnings=[], notes=[], java_version="21.0.7", nextflow_version="25.10.4"
        )
        monkeypatch.setattr(module, "_write_macos_docker_config", lambda *a, **kw: None)
        monkeypatch.setattr(
            module, "_write_resource_limits_config", lambda *a, **kw: None
        )
    else:
        monkeypatch.setattr(
            module, "_build_extra_nextflow_configs", lambda *a, **kw: []
        )
    monkeypatch.setattr(module, "run_preflight", lambda *a, **kw: preflight)

    def fail(*a, **kw):
        assert (tmp_path / "run/reproducibility/manifest.json").is_file(), (
            "manifest absent when Nextflow starts"
        )
        raise companions["errors"].SkillError(
            stage="execution",
            error_code="EXECUTION_FAILED",
            message="synthetic failure",
            fix="resume",
            details={},
        )

    monkeypatch.setattr(module, "execute_nextflow", fail)
    with pytest.raises(companions["errors"].SkillError):
        module._run_wrapper(args, tmp_path / "run")
    assert (tmp_path / "run/reproducibility/manifest.json").is_file()


@pytest.mark.parametrize("wrapper", ["scrnaseq"], indirect=True)
def test_single_cell_resume_rejects_changed_sheet(wrapper, tmp_path):
    import hashlib, json

    _, module, companions = wrapper
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "manifest.json").write_text(
        json.dumps({"samplesheet_checksum": hashlib.sha256(b"original").hexdigest()})
    )
    changed = tmp_path / "new.csv"
    changed.write_text("changed")
    args = module.build_parser().parse_args(["--output", str(tmp_path), "--resume"])
    with pytest.raises(companions["errors"].SkillError) as exc:
        module._check_resume_samplesheet_checksum(
            args, output_dir=tmp_path, staged_samplesheet=changed
        )
    assert exc.value.error_code == "INVALID_RESUME_STATE"


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_failed_resume_preserves_original_sheet(wrapper, tmp_path, monkeypatch):
    _, module, companions = wrapper
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    original = repro / "samplesheet.demo.csv"
    original.write_text("previous snapshot\n")
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--demo", "--resume"]
    )
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: {
            "source_kind": "remote_repo",
            "source_ref": "nf-core/sarek",
            "resolved_version": "3.8.1",
        },
    )

    def fail(**kw):
        raise companions["errors"].SkillError(
            stage="preflight",
            error_code="INVALID_RESUME_STATE",
            message="drift",
            fix="same flags",
            details={},
        )

    monkeypatch.setattr(module, "run_preflight", fail)
    with pytest.raises(companions["errors"].SkillError):
        module._run_wrapper(args, tmp_path)
    assert original.read_text() == "previous snapshot\n"


@pytest.mark.parametrize("wrapper", ["scrnaseq"], indirect=True)
def test_user_configs_have_last_precedence(wrapper, tmp_path, monkeypatch):
    _, module, _ = wrapper
    user = tmp_path / "site.config"
    user.write_text("process.resourceLimits = [memory: '32.GB']\n")
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--config", str(user)]
    )
    monkeypatch.setattr(module.platform, "system", lambda: "Darwin")
    configs = module._build_extra_nextflow_configs(args, tmp_path)
    assert configs[-1] == user.resolve()


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_duplicate_json_keys_rejected(wrapper, tmp_path):
    _, _, companions = wrapper
    source = tmp_path / "samples.json"
    source.write_text(
        '[{"patient":"P1","sample":"original","sample":"replaced","lane":"1","fastq_1":"s3://b/R1.fastq.gz","fastq_2":"s3://b/R2.fastq.gz"}]'
    )
    with pytest.raises(companions["errors"].SkillError) as exc:
        companions["samplesheet_builder"].validate_and_normalize_samplesheet(
            source, tmp_path / "out.csv"
        )
    assert exc.value.error_code == "INVALID_SAMPLESHEET"


@pytest.mark.parametrize("wrapper", ["rnaseq", "sarek"], indirect=True)
@pytest.mark.parametrize("value", [float("nan"), float("inf")])
def test_emitted_numeric_params_validated(wrapper, tmp_path, value):
    name, module, companions = wrapper
    if name == "scrnaseq":
        pytest.skip("single-cell CLI has no floating scientific params")
    args = module.build_parser().parse_args(["--output", str(tmp_path)])
    if name == "rnaseq":
        args.aligner = "star_salmon"
        args.min_mapped_reads = value
    else:
        args.ascat_ploidy = value
    with pytest.raises(companions["errors"].SkillError):
        companions["params_builder"].build_effective_params(
            args, normalized_samplesheet=tmp_path / "sheet.csv", output_dir=tmp_path
        )


@pytest.mark.parametrize("wrapper", ["rnaseq", "scrnaseq"], indirect=True)
def test_resume_check_rejects_changed_sheet_without_overwrite(
    wrapper, tmp_path, monkeypatch
):
    import hashlib, json

    name, module, companions = wrapper
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    original = repro / "samplesheet.demo.csv"
    original.write_text("previous snapshot\n")
    (repro / "manifest.json").write_text(
        json.dumps(
            {"samplesheet_checksum": hashlib.sha256(b"previous snapshot\n").hexdigest()}
        )
    )
    args = module.build_parser().parse_args(
        ["--demo", "--check", "--resume", "--output", str(tmp_path)]
    )
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: {
            "source_kind": "remote_repo",
            "source_ref": f"nf-core/{name}",
            "resolved_version": "test",
        },
    )
    monkeypatch.setattr(
        module,
        "run_preflight",
        lambda *a, **kw: {
            "java": {"version": "21"},
            "nextflow": {"version": "25.10.4"},
            "warnings": [],
        },
    )
    monkeypatch.setattr(module, "_check_resume_params_checksum", lambda *a, **kw: None)
    with pytest.raises(companions["errors"].SkillError) as exc:
        module._run_wrapper(args, tmp_path)
    assert exc.value.error_code == "INVALID_RESUME_STATE"
    assert original.read_text() == "previous snapshot\n"


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_params_file_custom_fasta_disables_igenomes(wrapper, tmp_path):
    _, module, _ = wrapper
    params = tmp_path / "params.yaml"
    params.write_text("fasta: /refs/custom.fa\n")
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--params-file", str(params)]
    )
    assert module._auto_disable_igenomes_for_custom_fasta(args)
    assert args.igenomes_ignore is True


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_demo_clears_params_file_inputs_and_references(wrapper, tmp_path):
    _, module, companions = wrapper
    source = tmp_path / "params.yaml"
    source.write_text(
        "input: /private/user_samples.csv\nfasta: /private/custom.fa\ngenome: GATK.GRCh38\nigenomes_ignore: true\n"
    )
    args = module.build_parser().parse_args(
        ["--demo", "--output", str(tmp_path), "--params-file", str(source)]
    )
    module._apply_demo_overrides(args)
    params = companions["params_builder"].build_effective_params(
        args, normalized_samplesheet=tmp_path / "samples.csv", output_dir=tmp_path
    )
    preflight = module._build_params_for_preflight(args, composed_profile="test,docker")
    for field in ("input", "fasta", "genome", "igenomes_ignore"):
        assert field not in params and field not in preflight


@pytest.mark.parametrize("wrapper", ["rnaseq"], indirect=True)
def test_bulk_resume_rejects_work_dir_change(wrapper, tmp_path):
    import json

    _, module, companions = wrapper
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--resume", "--work-dir", "s3://new/work"]
    )
    args.aligner = "star_salmon"
    source = {"source_kind": "remote_repo", "resolved_version": "3.26.0"}
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    manifest = {
        "aligner": "star_salmon",
        "pseudo_aligner": None,
        "profile": "docker",
        "prokaryotic": False,
        "pipeline_source": source,
        "work_dir": "s3://old/work",
    }
    (repro / "manifest.json").write_text(json.dumps(manifest))
    preflight = module.run_preflight.__globals__["_check_resume_compatibility"]
    with pytest.raises(companions["errors"].SkillError):
        preflight(args, output_dir=tmp_path, pipeline_source=source)


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_sarek_resume_rejects_work_dir_change(wrapper):
    _, module, companions = wrapper
    preflight = module.run_preflight.__globals__["_check_resume_drift"]
    with pytest.raises(companions["errors"].SkillError):
        preflight(
            params={"work_dir": "s3://new/work"},
            samplesheet={},
            pipeline_source={},
            manifest={"work_dir": "s3://old/work"},
        )


@pytest.mark.parametrize("wrapper", ["rnaseq", "scrnaseq"], indirect=True)
def test_malformed_resume_manifest_is_structured(wrapper, tmp_path):
    _, module, companions = wrapper
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "manifest.json").write_text("[]")
    args = module.build_parser().parse_args(["--output", str(tmp_path), "--resume"])
    preflight = module.run_preflight.__globals__["_check_resume_compatibility"]
    with pytest.raises(companions["errors"].SkillError) as exc:
        preflight(args, output_dir=tmp_path, pipeline_source={})
    assert exc.value.error_code == "INVALID_RESUME_STATE"


@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(st.sampled_from([" ", "\t", "\n", "\r", "\u00a0"]))
def test_local_checkout_whitespace_detection(wrapper, character):
    _, _, companions = wrapper
    assert companions["pipeline_source"]._path_has_whitespace(
        Path("/pipelines") / f"rna{character}seq"
    )


def test_scheduler_capacity_matches_container_runtime(wrapper, tmp_path, monkeypatch):
    name, module, _ = wrapper
    monkeypatch.setattr(module.sys, "platform", "linux")
    monkeypatch.setattr(module, "_runtime_cpu_count", lambda: 2)
    monkeypatch.setattr(module, "_resource_limits_memory_gb", lambda: 1)
    args = module.build_parser().parse_args(["--output", str(tmp_path)])
    config = (
        module._write_resource_limits_config(tmp_path)
        if name == "rnaseq"
        else module._write_resource_limits_config(tmp_path, args=args)
    )
    text = config.read_text()
    assert "executor.cpus = 2" in text
    assert "executor.memory = '1.GB'" in text


def test_site_profile_does_not_use_launcher_resource_capacity(
    wrapper, tmp_path, monkeypatch
):
    name, module, _ = wrapper
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--profile", "docker,my_hpc"]
    )
    monkeypatch.setattr(module.sys, "platform", "linux")
    if hasattr(module, "platform"):
        monkeypatch.setattr(module.platform, "system", lambda: "Linux")
    if name == "sarek":
        assert module._write_resource_limits_config(tmp_path, args=args) is None
    else:
        assert module._build_extra_nextflow_configs(args, tmp_path) == []


@pytest.mark.integration
def test_real_nextflow_honors_generated_scheduler_limits(
    wrapper, tmp_path, monkeypatch
):
    import shutil, subprocess
    from datetime import datetime

    if not shutil.which("nextflow") or not shutil.which("java"):
        pytest.skip("Nextflow and Java are required for the scheduler smoke test")
    name, module, _ = wrapper
    monkeypatch.setattr(module.sys, "platform", "linux")
    monkeypatch.setattr(module, "_runtime_cpu_count", lambda: 2)
    monkeypatch.setattr(module, "_resource_limits_memory_gb", lambda: 1)
    args = module.build_parser().parse_args(["--output", str(tmp_path)])
    config = (
        module._write_resource_limits_config(tmp_path)
        if name == "rnaseq"
        else module._write_resource_limits_config(tmp_path, args=args)
    )
    config.write_text(
        config.read_text()
        + "\ntrace.fields = ['name','status','start','complete','memory','cpus']\n"
    )
    workflow = tmp_path / "main.nf"
    workflow.write_text('''process capped {
    cpus 1
    memory '1.GB'
    input:
    val n
    output:
    val n
    script:
    """
    sleep 1
    """
}
workflow { capped(Channel.of(1, 2, 3)) }
''')
    trace = tmp_path / "trace.tsv"
    result = subprocess.run(
        [
            "nextflow",
            "run",
            str(workflow),
            "-c",
            str(config),
            "-ansi-log",
            "false",
            "-with-trace",
            str(trace),
        ],
        cwd=tmp_path,
        capture_output=True,
        text=True,
        timeout=120,
    )
    assert result.returncode == 0, result.stdout + result.stderr
    with trace.open() as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    assert len(rows) == 3 and all(row["status"] == "COMPLETED" for row in rows)
    intervals = sorted(
        (datetime.fromisoformat(row["start"]), datetime.fromisoformat(row["complete"]))
        for row in rows
    )
    assert all(
        after[0] >= before[1] for before, after in zip(intervals, intervals[1:])
    ), rows


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_sarek_existing_output_file_has_structured_error(wrapper, tmp_path):
    _, module, companions = wrapper
    output = tmp_path / "existing.txt"
    output.write_text("keep this content")
    args = module.build_parser().parse_args(["--demo", "--output", str(output)])
    with pytest.raises(companions["errors"].SkillError) as exc:
        module._run_wrapper(args, output)
    assert exc.value.error_code == "OUTPUT_DIR_NOT_WRITABLE"
    assert output.read_text() == "keep this content"


@pytest.mark.parametrize("wrapper", ["scrnaseq"], indirect=True)
def test_macos_replay_uses_recorded_runtime_limits(wrapper, tmp_path):
    import os
    import subprocess

    _, module, _ = wrapper
    builder = module.write_repro_commands.__globals__["build_nextflow_commands_sh"]
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "params.yaml").write_text("outdir: upstream/results\n")
    (repro / "resource_limits.config").write_text(
        "executor.cpus = 2\nexecutor.memory = '1.GB'\n"
    )
    script = repro / "commands.sh"
    script.write_text(
        builder(
            pipeline_source={
                "source_kind": "remote_repo",
                "source_ref": "nf-core/scrnaseq",
                "resolved_version": "4.1.0",
            },
            profile="docker",
            demo=False,
            macos_docker_config=False,
            resource_limits_config=True,
            nextflow_version=None,
            generated_at="2026-10-14T00:00:00Z",
        )
    )
    binaries = tmp_path / "bin"
    binaries.mkdir()
    (binaries / "uname").write_text('#!/bin/sh\nprintf "Darwin\\n"\n')
    (binaries / "docker").write_text('#!/bin/sh\nprintf "2 2147483648\\n"\n')
    (binaries / "nextflow").write_text('#!/bin/sh\nprintf "%s\\n" "$@" > "$CAPTURE"\n')
    for path in binaries.iterdir():
        path.chmod(0o755)
    capture = tmp_path / "argv.txt"
    result = subprocess.run(
        ["bash", str(script)],
        env={
            **os.environ,
            "PATH": str(binaries) + os.pathsep + os.environ["PATH"],
            "CAPTURE": str(capture),
        },
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
    assert (
        "reproducibility/resource_limits.replay.config"
        in capture.read_text().splitlines()
    )


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_sarek_environment_uses_shared_writer(wrapper, tmp_path, monkeypatch):
    _, module, _ = wrapper
    calls = []

    def write_environment(output_dir, env_name, **kwargs):
        calls.append(kwargs)
        target = output_dir / "reproducibility/environment.yml"
        target.write_text("name: test\n")
        return target

    monkeypatch.setitem(
        module.write_provenance_bundle.__globals__,
        "write_environment_yml",
        write_environment,
    )
    samples = tmp_path / "samples.csv"
    samples.write_text("patient,sample\n")
    params = tmp_path / "params.yaml"
    params.write_text("outdir: upstream/results\n")
    commands = tmp_path / "commands.sh"
    commands.write_text("#!/bin/sh\n")
    module.write_provenance_bundle(
        output_dir=tmp_path,
        skill_dir=ROOT / "skills/nfcore-sarek-wrapper",
        samplesheet_csv_src=samples,
        params_yaml_src=params,
        commands_sh_src=commands,
        params={},
        samplesheet_report={},
        pipeline_source={
            "source_kind": "remote_repo",
            "source_ref": "nf-core/sarek",
            "resolved_version": "3.8.1",
        },
        outputs_report=None,
        nextflow_version="25.10.4",
    )
    assert calls and calls[0]["conda_deps"] == ["nextflow=25.10.4", "openjdk=17"]


@pytest.mark.parametrize("wrapper", ["scrnaseq"], indirect=True)
@settings(suppress_health_check=[HealthCheck.function_scoped_fixture])
@given(
    st.sampled_from(["10xv1", "10x_v2", "10x-v3", "10X V4", "DROP-SEQ", "SMART_SEQ"])
)
def test_known_protocol_aliases_emit_upstream_lookup_keys(wrapper, tmp_path, value):
    _, module, companions = wrapper
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--preset", "star", "--protocol", value]
    )
    params = companions["params_builder"].build_effective_params(
        args, normalized_samplesheet=tmp_path / "samples.csv", output_dir=tmp_path
    )
    assert params["protocol"] in {
        "10XV1",
        "10XV2",
        "10XV3",
        "10XV4",
        "dropseq",
        "smartseq",
    }


@pytest.mark.parametrize("wrapper", ["scrnaseq"], indirect=True)
def test_single_cell_replay_keeps_user_config_last(wrapper):
    _, module, _ = wrapper
    build = module.write_repro_commands.__globals__["build_nextflow_commands_sh"]
    script = build(
        pipeline_source={
            "source_kind": "remote_repo",
            "source_ref": "nf-core/scrnaseq",
            "resolved_version": "4.1.0",
        },
        profile="docker",
        demo=False,
        macos_docker_config=True,
        resource_limits_config=True,
        extra_configs=["reproducibility/nextflow_configs/site_override.config"],
        nextflow_version=None,
        generated_at="2026-10-14T00:00:00Z",
    )
    command = script[script.index("nextflow run ") :]
    assert (
        command.index("$EXTRA_CONFIG")
        < command.index("$RESOURCE_CONFIG")
        < command.index("site_override.config")
    )


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_sarek_replay_keeps_user_config_last(wrapper, tmp_path):
    _, module, _ = wrapper
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "params.yaml").write_text("outdir: upstream/results\n")
    site = tmp_path / "site_override.config"
    site.write_text("process.cpus = 2\n")
    render = module.write_reports.__globals__["_build_commands_sh"]
    command = [
        "nextflow",
        "run",
        "nf-core/sarek",
        "-r",
        "3.8.1",
        "-params-file",
        str(repro / "params.yaml"),
        "-c",
        str(repro / "resource_limits.config"),
        "-c",
        str(site),
    ]
    script, _ = render(
        output_dir=tmp_path,
        generated_at="2026-10-14",
        nextflow_command=command,
        profile="docker",
        pipeline_version="3.8.1",
        pipeline_source_kind="remote_repo",
    )
    run = script[script.index("# ── Replay the captured") :]
    assert run.index("$RESOURCE_CONFIG") < run.index("site_override.config")


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_user_config_with_reserved_basename_is_not_dropped(wrapper):
    _, module, _ = wrapper
    split = module.write_reports.__globals__["_split_macos_docker_config"]
    config = "$SCRIPT_DIR/nextflow_configs/00_resource_limits.config"
    argv, macos, resources = split(["nextflow", "run", "nf-core/sarek", "-c", config])
    assert argv[-2:] == ["-c", config] and not resources and not macos


def test_error_payloads_are_valid_json_with_nonfinite_details(wrapper):
    import json

    _, _, companions = wrapper
    error = companions["errors"].SkillError(
        stage="preflight",
        error_code="INVALID_FLAG_COMBINATION",
        message="invalid numeric input",
        fix="use finite values",
        details={"value": float("nan"), "nested": [float("inf")]},
    )
    json.dumps(error.to_dict(), allow_nan=False)


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_restart_bam_path_matches_upstream_full_path_pattern(wrapper, tmp_path):
    _, _, companions = wrapper
    directory = tmp_path / "reference with spaces"
    directory.mkdir()
    bam = directory / "sample.bam"
    bam.touch()
    bai = directory / "sample.bam.bai"
    bai.touch()
    source = tmp_path / "samples.csv"
    source.write_text(f"patient,sample,bam,bai\nP1,S1,{bam},{bai}\n")
    with pytest.raises(companions["errors"].SkillError) as exc:
        companions["samplesheet_builder"].validate_and_normalize_samplesheet(
            source, tmp_path / "out.csv", step="variant_calling"
        )
    assert exc.value.error_code in {"INVALID_SAMPLESHEET", "INVALID_BAM"}


@pytest.mark.parametrize("wrapper", ["sarek"], indirect=True)
def test_sarek_check_validates_emitted_schema_before_committing_sheet(
    wrapper, tmp_path, monkeypatch
):
    from types import SimpleNamespace

    _, module, companions = wrapper
    output = tmp_path / "out"
    args = module.build_parser().parse_args(
        [
            "--demo",
            "--check",
            "--output",
            str(output),
            "--extra-param",
            "publish_dir_mode=invalid",
        ]
    )
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: {
            "source_kind": "remote_repo",
            "source_ref": "nf-core/sarek",
            "resolved_version": "3.8.1",
        },
    )
    monkeypatch.setattr(
        module, "run_preflight", lambda **kw: SimpleNamespace(warnings=[], notes=[])
    )
    monkeypatch.setattr(module, "_write_check_report", lambda **kw: 0)
    with pytest.raises(companions["errors"].SkillError) as exc:
        module._run_wrapper(args, output)
    assert exc.value.error_code == "INVALID_FLAG_COMBINATION"
    assert "publish_dir_mode" in exc.value.message
    assert not (output / "reproducibility/samplesheet.demo.csv").exists()


def test_generated_cpu_limit_fits_docker(wrapper, tmp_path, monkeypatch):
    name, module, _ = wrapper
    monkeypatch.setattr(module.sys, "platform", "linux")
    monkeypatch.setattr(module.os, "cpu_count", lambda: 64)
    monkeypatch.setattr(module, "_docker_vm_cpu_count", lambda: 2, raising=False)
    args = module.build_parser().parse_args(["--output", str(tmp_path)])
    config = (
        module._write_resource_limits_config(tmp_path)
        if name == "rnaseq"
        else module._write_resource_limits_config(tmp_path, args=args)
    )
    assert "cpus: 2," in config.read_text()


@pytest.mark.parametrize("wrapper", ["scrnaseq"], indirect=True)
def test_macos_real_single_cell_run_gets_runtime_cap(wrapper, tmp_path, monkeypatch):
    _, module, _ = wrapper
    monkeypatch.setattr(module.platform, "system", lambda: "Darwin")
    args = module.build_parser().parse_args(["--output", str(tmp_path)])
    configs = module._build_extra_nextflow_configs(args, tmp_path)
    assert "resource_limits.config" in [path.name for path in configs]
