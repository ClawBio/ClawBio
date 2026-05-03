from __future__ import annotations

import importlib.util
import json
import runpy
import subprocess
import sys
from pathlib import Path
from unittest.mock import patch

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = SKILL_DIR / "nfcore_scrnaseq_wrapper.py"
PROJECT_ROOT = SKILL_DIR.parent.parent
CLAWBIO_PATH = PROJECT_ROOT / "clawbio.py"

sys.path.insert(0, str(SKILL_DIR))

from errors import SkillError


def test_skill_modules_importable_without_running_script(tmp_path):
    """Sibling-module imports must work when loading files directly."""
    skill_root = Path(__file__).resolve().parent.parent
    skill_dir = str(skill_root)
    original = list(sys.path)
    sibling_modules = {
        "errors",
        "schemas",
        "command_builder",
        "executor",
        "outputs_parser",
        "params_builder",
        "pipeline_source",
        "preflight",
        "provenance",
        "reporting",
        "samplesheet_builder",
    }
    original_modules = {name: sys.modules.pop(name) for name in list(sys.modules) if name in sibling_modules}
    sys.path = [p for p in sys.path if p != skill_dir]
    try:
        for module_path in sorted(skill_root.glob("*.py")):
            if module_path.name == "nfcore_scrnaseq_wrapper.py":
                continue
            for name in sibling_modules:
                sys.modules.pop(name, None)
            spec = importlib.util.spec_from_file_location(f"isolated_{module_path.stem}", module_path)
            assert spec is not None
            assert spec.loader is not None
            mod = importlib.util.module_from_spec(spec)
            sys.modules[spec.name] = mod
            try:
                spec.loader.exec_module(mod)
            finally:
                sys.modules.pop(spec.name, None)
    finally:
        for name in sibling_modules:
            sys.modules.pop(name, None)
        sys.modules.update(original_modules)
        sys.path = original


def _load_skill_module():
    spec = importlib.util.spec_from_file_location("nfcore_scrnaseq_wrapper_module", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parser_accepts_expected_flags():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--input", "samplesheet.csv", "--output", "out", "--preset", "star", "--resume"])
    assert args.input == "samplesheet.csv"
    assert args.output == "out"
    assert args.preset == "star"
    assert args.resume is True


def test_parser_star_feature_rejects_invalid_value():
    """--star-feature must reject values not in the nf-core schema enum."""
    module = _load_skill_module()
    parser = module.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["--output", "out", "--star-feature", "Velocyto"])


def test_parser_star_feature_accepts_gene_velocyto_with_space():
    """'Gene Velocyto' (with space) is the correct schema enum value — must be accepted."""
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", "out", "--star-feature", "Gene Velocyto"])
    assert args.star_feature == "Gene Velocyto"


def test_parser_simpleaf_umi_resolution_rejects_invalid_value():
    """--simpleaf-umi-resolution must reject values not in the nf-core schema enum."""
    module = _load_skill_module()
    parser = module.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["--output", "out", "--simpleaf-umi-resolution", "invalid-strategy"])


def test_parser_kb_workflow_rejects_invalid_value():
    """--kb-workflow must reject values not in the nf-core schema enum."""
    module = _load_skill_module()
    parser = module.build_parser()
    with pytest.raises(SystemExit):
        parser.parse_args(["--output", "out", "--kb-workflow", "unknown"])


def test_main_writes_structured_error_when_input_missing(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/scrnaseq",
        "resolved_version": kw.get("requested_version", ""),
        "branch": "",
        "dirty": False,
    })
    monkeypatch.setattr(sys, "argv", ["nfcore_scrnaseq_wrapper.py", "--output", str(tmp_path)])
    rc = module.main()
    assert rc == 1
    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["error_code"] == "MISSING_INPUT"


def test_main_does_not_write_into_rejected_nonempty_output_dir(tmp_path, monkeypatch):
    """OUTPUT_DIR_NOT_EMPTY must fail before writing wrapper artifacts into the directory.

    Ordering dependency: check_output_dir_available() fires in main() *before*
    validate_and_normalize_samplesheet(), so samplesheet.valid.csv is never written.
    The second check inside run_preflight() is never reached in this path.
    """
    module = _load_skill_module()
    out = tmp_path / "out"
    out.mkdir()
    sentinel = out / "result.json"
    sentinel.write_text("do not overwrite", encoding="utf-8")
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")
    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/scrnaseq",
        "resolved_version": kw.get("requested_version", ""),
        "branch": "",
        "dirty": False,
    })
    monkeypatch.setattr(sys, "argv", [
        "nfcore_scrnaseq_wrapper.py",
        "--input", str(samplesheet),
        "--output", str(out),
        "--genome", "GRCh38",
    ])
    rc = module.main()
    assert rc == 1
    assert sentinel.read_text(encoding="utf-8") == "do not overwrite"
    assert not (out / "reproducibility" / "samplesheet.valid.csv").exists()


def test_output_path_existing_file_returns_structured_error(tmp_path, monkeypatch, capsys):
    module = _load_skill_module()
    output_file = tmp_path / "not_a_directory"
    output_file.write_text("already here", encoding="utf-8")
    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: (_ for _ in ()).throw(AssertionError("should not resolve pipeline source")))
    monkeypatch.setattr(sys, "argv", [
        "nfcore_scrnaseq_wrapper.py",
        "--output", str(output_file),
        "--demo",
    ])
    rc = module.main()
    assert rc == 1
    assert output_file.read_text(encoding="utf-8") == "already here"
    payload = json.loads(capsys.readouterr().err)
    assert payload["stage"] == "preflight"
    assert payload["error_code"] == "OUTPUT_DIR_NOT_WRITABLE"


@pytest.mark.integration
def test_clawbio_list_shows_skill():
    proc = subprocess.run(
        [sys.executable, str(CLAWBIO_PATH), "list"],
        capture_output=True,
        text=True,
        cwd=str(PROJECT_ROOT),
    )
    assert proc.returncode == 0
    assert "scrnaseq-pipeline" in proc.stdout


@pytest.mark.integration
def test_clawbio_registry_has_long_timeout_and_file_cap():
    registry = runpy.run_path(str(CLAWBIO_PATH))["SKILLS"]
    info = registry["scrnaseq-pipeline"]
    assert info["default_timeout_seconds"] > 60 * 60 * 12
    assert info["max_output_files_listed"] <= 50
    assert "--run-downstream" in info["allowed_extra_flags"]
    assert "--skip-downstream" in info["allowed_extra_flags_without_values"]


def test_clawbio_runner_rejects_output_path_that_is_file(tmp_path):
    output_file = tmp_path / "not_a_directory"
    output_file.write_text("already here", encoding="utf-8")
    run_skill = runpy.run_path(str(CLAWBIO_PATH))["run_skill"]

    result = run_skill("scrnaseq-pipeline", output_dir=str(output_file), demo=True)

    assert result["success"] is False
    assert result["exit_code"] == -1
    assert "Traceback" not in result["stderr"]
    payload = json.loads(result["stderr"])
    assert payload["stage"] == "preflight"
    assert payload["error_code"] == "OUTPUT_DIR_NOT_WRITABLE"
    assert output_file.read_text(encoding="utf-8") == "already here"


def test_write_macos_docker_config_produces_valid_groovy(tmp_path):
    module = _load_skill_module()
    config_path = module._write_macos_docker_config(tmp_path)
    assert config_path.exists()
    content = config_path.read_text(encoding="utf-8")
    assert "process {" in content
    assert "stageInMode" in content
    assert '"copy"' in content
    assert "docker {" in content
    assert "--platform linux/amd64" in content


def test_write_macos_docker_config_creates_in_reproducibility_dir(tmp_path):
    module = _load_skill_module()
    repro_dir = tmp_path / "reproducibility"
    repro_dir.mkdir()
    config_path = module._write_macos_docker_config(tmp_path)
    assert config_path.parent == repro_dir


def test_resume_params_checksum_mismatch_preserves_previous_run_artifacts(tmp_path, monkeypatch, capsys):
    """Invalid resume must not overwrite artifacts from the previous successful run."""
    module = _load_skill_module()

    fake_source = {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/scrnaseq",
        "resolved_version": "3.14.0",
        "branch": "",
        "dirty": False,
    }
    # Write a manifest with a checksum that will never match the generated params
    repro_dir = tmp_path / "reproducibility"
    repro_dir.mkdir(parents=True)
    manifest = {
        "preset": "star",
        "profile": "docker",
        "pipeline_source": fake_source,
        "params_checksum": "deadbeef" * 8,
    }
    (repro_dir / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")
    previous_samplesheet = repro_dir / "samplesheet.valid.csv"
    previous_params = repro_dir / "params.yaml"
    previous_result = tmp_path / "result.json"
    previous_samplesheet.write_text("previous samplesheet\n", encoding="utf-8")
    previous_params.write_text("previous params\n", encoding="utf-8")
    previous_result.write_text('{"ok": true, "previous": true}', encoding="utf-8")

    # Minimal samplesheet
    r1 = tmp_path / "r1.fastq.gz"
    r2 = tmp_path / "r2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    ss = tmp_path / "samplesheet.csv"
    ss.write_text(f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8")

    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: fake_source)
    monkeypatch.setattr(module, "run_preflight", lambda args, **kw: {
        "ok": True,
        "java": {"version": "21.0.0", "path": ""},
        "nextflow": {"version": "25.4.0", "path": ""},
        "profile": {"profile": "docker"},
        "pipeline_source": fake_source,
        "references": {},
        "samplesheet": {"sample_count": 1, "unknown_columns": []},
    })

    monkeypatch.setattr(sys, "argv", [
        "nfcore_scrnaseq_wrapper.py",
        "--input", str(ss),
        "--output", str(tmp_path),
        "--preset", "star",
        "--profile", "docker",
        "--resume",
        "--fasta", str(r1),
        "--gtf", str(r2),
    ])
    rc = module.main()
    assert rc == 1
    payload = json.loads(capsys.readouterr().err)
    assert payload["error_code"] == "INVALID_RESUME_STATE"
    assert previous_samplesheet.read_text(encoding="utf-8") == "previous samplesheet\n"
    assert previous_params.read_text(encoding="utf-8") == "previous params\n"
    assert json.loads(previous_result.read_text(encoding="utf-8")) == {"ok": True, "previous": True}


def test_main_returns_1_and_writes_result_on_unexpected_exception(tmp_path, monkeypatch):
    """An unexpected exception must produce a structured result.json with a traceback."""
    module = _load_skill_module()

    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: (_ for _ in ()).throw(RuntimeError("boom")))
    monkeypatch.setattr(sys, "argv", [
        "nfcore_scrnaseq_wrapper.py",
        "--output", str(tmp_path),
        "--demo",
    ])
    rc = module.main()
    assert rc == 1
    result_path = tmp_path / "result.json"
    assert result_path.exists()
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    assert payload["error_code"] == "UNEXPECTED_ERROR"
    assert "traceback" in payload["details"]
    assert "RuntimeError" in payload["details"]["traceback"]


def test_build_extra_nextflow_configs_returns_empty_on_linux(tmp_path, monkeypatch):
    """`_build_extra_nextflow_configs` must return [] on Linux regardless of profile."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(profile="docker")
    with patch("platform.system", return_value="Linux"):
        result = module._build_extra_nextflow_configs(args, tmp_path)
    assert result == [], f"Expected [], got {result}"


def test_build_extra_nextflow_configs_writes_config_on_macos_docker(tmp_path):
    """`_build_extra_nextflow_configs` must write the VirtioFS config on macOS with docker."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(profile="docker")
    with patch("platform.system", return_value="Darwin"):
        result = module._build_extra_nextflow_configs(args, tmp_path)
    assert len(result) == 1
    assert result[0].exists()
    assert "stageInMode" in result[0].read_text(encoding="utf-8")


def test_build_extra_nextflow_configs_returns_empty_on_macos_non_docker(tmp_path):
    """`_build_extra_nextflow_configs` must return [] on macOS when profile is not docker."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(profile="singularity")
    with patch("platform.system", return_value="Darwin"):
        result = module._build_extra_nextflow_configs(args, tmp_path)
    assert result == []


def test_prepare_demo_samplesheet_returns_three_tuple(tmp_path):
    """_prepare_demo_samplesheet must return (normalized, staged, summary) — a 3-tuple.

    Previously it returned only 2 values, causing a ValueError when the caller
    unpacked: normalized_samplesheet, staged_samplesheet, samplesheet_summary = ...
    """
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(preset="star", skip_cellbender=False)
    result = module._prepare_demo_samplesheet(args, tmp_path, staging_dir=None)
    assert len(result) == 3, f"Expected 3-tuple, got {len(result)}-tuple"
    normalized, staged, summary = result
    assert isinstance(normalized, __import__("pathlib").Path)
    assert isinstance(staged, __import__("pathlib").Path)
    assert isinstance(summary, dict)
    assert "sample_count" in summary
    assert "normalized_path" in summary


def test_prepare_demo_samplesheet_mutates_preset_before_writing(tmp_path):
    """Preset mutation must happen before file I/O."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(preset="kallisto", skip_cellbender=False)
    _, _, summary = module._prepare_demo_samplesheet(args, tmp_path, staging_dir=None)
    assert args.preset == "star"
    assert args.skip_cellbender is True
    assert summary["sample_count"] == 0


def test_check_resume_params_checksum_does_not_raise_when_manifest_missing(tmp_path):
    """When manifest.json is absent, _check_resume_params_checksum must return silently."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(resume=True)
    params_payload = {"aligner": "star", "outdir": "/tmp/out"}
    module._check_resume_params_checksum(args, output_dir=tmp_path, params_payload=params_payload)


def test_demo_samplesheet_written_as_demo_csv(tmp_path):
    """Demo mode must write samplesheet.demo.csv, not samplesheet.valid.csv."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(preset="star", skip_cellbender=False)
    normalized, staged, _ = module._prepare_demo_samplesheet(args, tmp_path, staging_dir=None)
    assert normalized.name == "samplesheet.demo.csv", \
        f"Expected samplesheet.demo.csv, got {normalized.name}"


def test_nextflow_execution_cwd_is_output_dir(tmp_path):
    """Relative params.input paths must resolve from the wrapper output directory."""
    module = _load_skill_module()

    assert module._nextflow_execution_cwd(tmp_path) == tmp_path


def test_nextflow_replay_command_records_launch_directory(tmp_path):
    module = _load_skill_module()
    output_dir = tmp_path / "output with spaces"

    replay = module._nextflow_replay_command("nextflow run nf-core/scrnaseq", output_dir)

    assert replay.startswith("cd ")
    assert "output with spaces" in replay
    assert replay.endswith("&& nextflow run nf-core/scrnaseq")


# ---------------------------------------------------------------------------
# Task 18: downstream scrna_orchestrator handoff
# ---------------------------------------------------------------------------

def test_parser_accepts_downstream_flags():
    """Downstream handoff must be explicit, with --skip-downstream kept for compatibility."""
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", "out"])
    assert hasattr(args, "run_downstream")
    assert hasattr(args, "skip_downstream")
    assert args.run_downstream is False
    assert args.skip_downstream is False
    run_args = parser.parse_args(["--output", "out", "--run-downstream"])
    assert run_args.run_downstream is True
    skip_args = parser.parse_args(["--output", "out", "--skip-downstream"])
    assert skip_args.skip_downstream is True


def test_run_downstream_handoff_skipped_when_no_h5ad(tmp_path, capsys):
    """When parsed_outputs contains no preferred_h5ad, no subprocess must be launched."""
    module = _load_skill_module()
    import argparse as _argparse
    import subprocess as _sp
    args = _argparse.Namespace(run_downstream=True, skip_downstream=False)
    launched = []
    original_run = _sp.run
    def fake_run(*a, **kw):
        launched.append(a)
        return original_run(["true"])
    with patch("subprocess.run", fake_run):
        module._run_downstream_handoff(args, parsed_outputs={"preferred_h5ad": ""}, output_dir=tmp_path)
    assert not launched, "No subprocess should launch when preferred_h5ad is empty"


def test_run_downstream_handoff_skipped_by_default(tmp_path):
    """Downstream handoff must not launch unless --run-downstream is set."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(run_downstream=False, skip_downstream=False)
    h5ad = tmp_path / "data.h5ad"
    h5ad.write_text("mock", encoding="utf-8")
    launched = []
    with patch("subprocess.run", lambda *a, **kw: launched.append(a)):
        module._run_downstream_handoff(args, parsed_outputs={"preferred_h5ad": str(h5ad)}, output_dir=tmp_path)
    assert not launched, "Downstream handoff should be opt-in"


def test_run_downstream_handoff_skipped_when_flag_set(tmp_path):
    """--skip-downstream must prevent handoff even if --run-downstream is set."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(run_downstream=True, skip_downstream=True)
    h5ad = tmp_path / "data.h5ad"
    h5ad.write_text("mock", encoding="utf-8")
    launched = []
    with patch("subprocess.run", lambda *a, **kw: launched.append(a)):
        module._run_downstream_handoff(args, parsed_outputs={"preferred_h5ad": str(h5ad)}, output_dir=tmp_path)
    assert not launched, "--skip-downstream should prevent subprocess launch"


def test_run_downstream_handoff_launches_scrna_orchestrator(tmp_path):
    """When a preferred_h5ad is found, scrna_orchestrator must be invoked."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(run_downstream=True, skip_downstream=False)
    h5ad = tmp_path / "combined_filtered_matrix.h5ad"
    h5ad.write_text("mock", encoding="utf-8")
    fake_orchestrator = tmp_path / "scrna_orchestrator.py"
    fake_orchestrator.write_text("# stub", encoding="utf-8")
    captured = {}
    def fake_run(cmd, **kw):
        captured["cmd"] = cmd
        class _R:
            returncode = 0
        return _R()
    with patch("subprocess.run", fake_run), \
         patch.object(module, "_resolve_scrna_orchestrator", return_value=fake_orchestrator):
        module._run_downstream_handoff(args, parsed_outputs={"preferred_h5ad": str(h5ad)}, output_dir=tmp_path)
    assert "cmd" in captured, "subprocess.run must be called when preferred_h5ad is set"
    cmd = captured["cmd"]
    assert any("scrna_orchestrator" in str(part) for part in cmd), \
        f"scrna_orchestrator.py must appear in the command, got: {cmd}"


def test_run_downstream_handoff_graceful_on_failure(tmp_path, capsys):
    """A failure in the downstream handoff must NOT raise — it should warn and continue."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(run_downstream=True, skip_downstream=False)
    h5ad = tmp_path / "combined_filtered_matrix.h5ad"
    h5ad.write_text("mock", encoding="utf-8")
    fake_orchestrator = tmp_path / "scrna_orchestrator.py"
    fake_orchestrator.write_text("# stub", encoding="utf-8")
    def fake_run(cmd, **kw):
        class _R:
            returncode = 1
            stderr = "scrna_orchestrator failed"
        return _R()
    with patch("subprocess.run", fake_run), \
         patch.object(module, "_resolve_scrna_orchestrator", return_value=fake_orchestrator):
        module._run_downstream_handoff(args, parsed_outputs={"preferred_h5ad": str(h5ad)}, output_dir=tmp_path)
    # Must not raise; a warning should appear on stderr
    err = capsys.readouterr().err
    assert "WARNING" in err or "warn" in err.lower() or "downstream" in err.lower()


def test_run_downstream_handoff_graceful_on_timeout(tmp_path, capsys):
    """A downstream timeout must warn and preserve the successful upstream result."""
    module = _load_skill_module()
    import argparse as _argparse
    import subprocess as _sp
    args = _argparse.Namespace(run_downstream=True, skip_downstream=False)
    h5ad = tmp_path / "combined_filtered_matrix.h5ad"
    h5ad.write_text("mock", encoding="utf-8")
    fake_orchestrator = tmp_path / "scrna_orchestrator.py"
    fake_orchestrator.write_text("# stub", encoding="utf-8")

    def fake_run(cmd, **kw):
        raise _sp.TimeoutExpired(cmd=cmd, timeout=kw.get("timeout"))

    with patch("subprocess.run", fake_run), \
         patch.object(module, "_resolve_scrna_orchestrator", return_value=fake_orchestrator):
        module._run_downstream_handoff(args, parsed_outputs={"preferred_h5ad": str(h5ad)}, output_dir=tmp_path)

    err = capsys.readouterr().err
    assert "timed out" in err


def test_run_downstream_handoff_warns_when_orchestrator_missing(tmp_path, capsys):
    """When scrna_orchestrator is not found, a WARNING must be printed to stderr."""
    module = _load_skill_module()
    import argparse as _argparse
    args = _argparse.Namespace(run_downstream=True, skip_downstream=False)
    h5ad = tmp_path / "combined_filtered_matrix.h5ad"
    h5ad.write_text("mock", encoding="utf-8")
    with patch.object(module, "_resolve_scrna_orchestrator", return_value=None):
        module._run_downstream_handoff(args, parsed_outputs={"preferred_h5ad": str(h5ad)}, output_dir=tmp_path)
    err = capsys.readouterr().err
    assert "WARNING" in err
    assert "CLAWBIO_SCRNA_ORCHESTRATOR" in err


def test_resolve_scrna_orchestrator_uses_env_var(tmp_path, monkeypatch):
    """CLAWBIO_SCRNA_ORCHESTRATOR env var must override the default path."""
    module = _load_skill_module()
    custom_path = tmp_path / "custom_scrna_orchestrator.py"
    custom_path.write_text("# stub", encoding="utf-8")
    monkeypatch.setenv("CLAWBIO_SCRNA_ORCHESTRATOR", str(custom_path))
    result = module._resolve_scrna_orchestrator()
    assert result == custom_path


def test_resolve_scrna_orchestrator_returns_none_when_env_path_missing(tmp_path, monkeypatch):
    """CLAWBIO_SCRNA_ORCHESTRATOR pointing to a non-existent file must return None."""
    module = _load_skill_module()
    monkeypatch.setenv("CLAWBIO_SCRNA_ORCHESTRATOR", str(tmp_path / "nonexistent.py"))
    result = module._resolve_scrna_orchestrator()
    assert result is None
