from __future__ import annotations

import importlib.util
import json
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
    original_modules = {
        name: sys.modules.pop(name)
        for name in list(sys.modules)
        if name in sibling_modules
    }
    sys.path = [p for p in sys.path if p != skill_dir]
    try:
        for module_path in sorted(skill_root.glob("*.py")):
            if module_path.name == "nfcore_scrnaseq_wrapper.py":
                continue
            for name in sibling_modules:
                sys.modules.pop(name, None)
            spec = importlib.util.spec_from_file_location(
                f"isolated_{module_path.stem}", module_path
            )
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
    spec = importlib.util.spec_from_file_location(
        "nfcore_scrnaseq_wrapper_module", SCRIPT_PATH
    )
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


def test_parser_accepts_expected_flags():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(
        [
            "--input",
            "samplesheet.csv",
            "--output",
            "out",
            "--preset",
            "star",
            "--resume",
        ]
    )
    assert args.input == "samplesheet.csv"
    assert args.output == "out"
    assert args.preset == "star"
    assert args.resume is True


def test_aligner_alias_maps_to_wrapper_preset():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", "out", "--aligner", "simpleaf"])

    module._apply_aligner_alias(args)

    assert args.preset == "standard"


def test_missing_preset_defaults_after_alias_normalization():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", "out"])

    module._apply_aligner_alias(args)

    assert args.preset == "standard"
    assert args.preset_explicit is False


def test_aligner_alias_overrides_implicit_default_preset():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", "out", "--aligner", "star"])

    module._apply_aligner_alias(args)

    assert args.preset == "star"
    assert args.preset_explicit is False


def test_aligner_alias_rejects_conflicting_explicit_preset():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(
        ["--output", "out", "--preset", "star", "--aligner", "simpleaf"]
    )

    with pytest.raises(SkillError) as exc:
        module._apply_aligner_alias(args)

    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"


def test_aligner_alias_rejects_conflict_with_explicit_default_preset():
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(
        ["--output", "out", "--preset", "standard", "--aligner", "star"]
    )

    with pytest.raises(SkillError) as exc:
        module._apply_aligner_alias(args)

    assert exc.value.details["preset"] == "standard"
    assert exc.value.details["preset_from_aligner"] == "star"


def test_demo_does_not_warn_about_preset_when_none_was_requested(
    tmp_path, capsys
):
    # --demo without --preset should not claim the user "requested" the
    # default preset; the warning is only meaningful for an explicit choice.
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(["--output", str(tmp_path), "--demo"])

    module._prepare_demo_samplesheet(args, tmp_path, staging_dir=None)

    stderr = capsys.readouterr().err
    assert "forces preset=star" not in stderr
    assert args.preset == "star"


def test_demo_warns_when_explicit_non_star_preset_is_overridden(
    tmp_path, capsys
):
    # An explicit --preset that --demo overrides still warns.
    module = _load_skill_module()
    parser = module.build_parser()
    args = parser.parse_args(
        ["--output", str(tmp_path), "--demo", "--preset", "standard"]
    )

    module._prepare_demo_samplesheet(args, tmp_path, staging_dir=None)

    stderr = capsys.readouterr().err
    assert "forces preset=star (requested: 'standard')" in stderr
    assert args.preset == "star"


def test_main_writes_structured_error_for_aligner_preset_conflict(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "nfcore_scrnaseq_wrapper.py",
            "--output",
            str(tmp_path),
            "--preset",
            "standard",
            "--aligner",
            "star",
        ],
    )

    rc = module.main()

    assert rc == 1
    payload = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert payload["error_code"] == "INVALID_PRESET_CONFIGURATION"


def test_main_accepts_explicit_argv(tmp_path, monkeypatch):
    """main() must accept an explicit argv list (parity with nfcore-sarek/rnaseq)
    so the entrypoint is unit-testable without mutating sys.argv."""
    module = _load_skill_module()

    def _boom(*_args, **_kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(module, "_run_wrapper", _boom)
    # Reaching _run_wrapper proves argv was parsed by main(argv=...).
    assert module.main(["--output", str(tmp_path), "--demo"]) == 130


def test_main_keyboard_interrupt_returns_130(tmp_path, monkeypatch):
    """Ctrl+C during a long-running pipeline must exit 130 (SIGINT convention),
    not dump a traceback — parity with nfcore-sarek/rnaseq."""
    module = _load_skill_module()

    def _boom(*_args, **_kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(module, "_run_wrapper", _boom)
    rc = module.main(["--output", str(tmp_path), "--demo"])
    assert rc == 130


def test_parser_save_align_intermeds_is_two_way_switch():
    """--save-align-intermeds must be a two-way switch so the BAM-saving default
    can be explicitly turned OFF (--no-save-align-intermeds → False), turned ON
    (--save-align-intermeds → True), or left unset (None = defer to nf-core)."""
    module = _load_skill_module()
    parser = module.build_parser()
    assert parser.parse_args(["--output", "out"]).save_align_intermeds is None
    assert (
        parser.parse_args(
            ["--output", "out", "--save-align-intermeds"]
        ).save_align_intermeds
        is True
    )
    assert (
        parser.parse_args(
            ["--output", "out", "--no-save-align-intermeds"]
        ).save_align_intermeds
        is False
    )


def test_parser_accepts_nfcore_profiles_and_composite_profiles():
    """nf-core profiles can be combined as comma-separated Nextflow profile strings."""
    module = _load_skill_module()
    parser = module.build_parser()
    args_wave = parser.parse_args(["--output", "out", "--profile", "wave"])
    assert args_wave.profile == "wave"
    args_gpu = parser.parse_args(["--output", "out", "--profile", "gpu"])
    assert args_gpu.profile == "gpu"
    args_composite = parser.parse_args(
        ["--output", "out", "--profile", "docker,emulate_amd64"]
    )
    assert args_composite.profile == "docker,emulate_amd64"


def test_main_writes_structured_error_when_input_missing(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: {
            "source_kind": "remote_repo",
            "source_ref": "nf-core/scrnaseq",
            "resolved_version": kw.get("requested_version", ""),
            "branch": "",
            "dirty": False,
        },
    )
    monkeypatch.setattr(
        sys, "argv", ["nfcore_scrnaseq_wrapper.py", "--output", str(tmp_path)]
    )
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
    samplesheet.write_text(
        f"sample,fastq_1,fastq_2\nsampleA,{r1},{r2}\n", encoding="utf-8"
    )
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: {
            "source_kind": "remote_repo",
            "source_ref": "nf-core/scrnaseq",
            "resolved_version": kw.get("requested_version", ""),
            "branch": "",
            "dirty": False,
        },
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "nfcore_scrnaseq_wrapper.py",
            "--input",
            str(samplesheet),
            "--output",
            str(out),
            "--genome",
            "GRCh38",
        ],
    )
    rc = module.main()
    assert rc == 1
    assert sentinel.read_text(encoding="utf-8") == "do not overwrite"
    assert not (out / "reproducibility" / "samplesheet.valid.csv").exists()


def test_output_path_existing_file_returns_structured_error(
    tmp_path, monkeypatch, capsys
):
    module = _load_skill_module()
    output_file = tmp_path / "not_a_directory"
    output_file.write_text("already here", encoding="utf-8")
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: (_ for _ in ()).throw(
            AssertionError("should not resolve pipeline source")
        ),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "nfcore_scrnaseq_wrapper.py",
            "--output",
            str(output_file),
            "--demo",
        ],
    )
    rc = module.main()
    assert rc == 1
    assert output_file.read_text(encoding="utf-8") == "already here"
    payload = json.loads(capsys.readouterr().err)
    assert payload["stage"] == "preflight"
    assert payload["error_code"] == "OUTPUT_DIR_NOT_WRITABLE"



def test_resume_accepts_a_bundle_fingerprinted_before_params_stopped_wrapping(tmp_path):
    """params.yaml used to be written with PyYAML's 80-column wrapping, and the
    resume fingerprint is the hash of that text. A bundle made then must still
    resume, while new fingerprints match the unwrapped params.yaml on disk."""
    import argparse
    import hashlib

    import yaml

    module = _load_skill_module()
    params = {"fasta": "/data/" + "x" * 80 + " release/genome.fa", "outdir": "upstream/results"}
    header = module.serialize_params_yaml({})[: -len("{}\n")]
    legacy = hashlib.sha256((header + yaml.dump(params, allow_unicode=True, sort_keys=False)).encode()).hexdigest()
    repro = tmp_path / "reproducibility"
    repro.mkdir()
    (repro / "manifest.json").write_text(json.dumps({"params_checksum": legacy}), encoding="utf-8")
    module._check_resume_params_checksum(argparse.Namespace(resume=True), output_dir=tmp_path, params_payload=params)

    written = module.write_params_yaml(params, output_dir=tmp_path)
    assert module._params_payload_checksum(params) == hashlib.sha256(written.read_bytes()).hexdigest()

def test_resume_params_checksum_mismatch_preserves_previous_run_artifacts(
    tmp_path, monkeypatch, capsys
):
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
    monkeypatch.setattr(
        module,
        "run_preflight",
        lambda args, **kw: {
            "ok": True,
            "java": {"version": "21.0.0", "path": ""},
            "nextflow": {"version": "25.4.0", "path": ""},
            "profile": {"profile": "docker"},
            "pipeline_source": fake_source,
            "references": {},
            "samplesheet": {"sample_count": 1, "unknown_columns": []},
        },
    )

    monkeypatch.setattr(
        sys,
        "argv",
        [
            "nfcore_scrnaseq_wrapper.py",
            "--input",
            str(ss),
            "--output",
            str(tmp_path),
            "--preset",
            "star",
            "--profile",
            "docker",
            "--resume",
            "--fasta",
            str(r1),
            "--gtf",
            str(r2),
        ],
    )
    rc = module.main()
    assert rc == 1
    payload = json.loads(capsys.readouterr().err)
    assert payload["error_code"] == "INVALID_RESUME_STATE"
    assert previous_samplesheet.read_text(encoding="utf-8") == "previous samplesheet\n"
    assert previous_params.read_text(encoding="utf-8") == "previous params\n"
    assert json.loads(previous_result.read_text(encoding="utf-8")) == {
        "ok": True,
        "previous": True,
    }


def test_main_returns_1_and_writes_result_on_unexpected_exception(
    tmp_path, monkeypatch
):
    """An unexpected exception must produce a structured result.json with a traceback."""
    module = _load_skill_module()

    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: (_ for _ in ()).throw(RuntimeError("boom")),
    )
    monkeypatch.setattr(
        sys,
        "argv",
        [
            "nfcore_scrnaseq_wrapper.py",
            "--output",
            str(tmp_path),
            "--demo",
        ],
    )
    rc = module.main()
    assert rc == 1
    result_path = tmp_path / "result.json"
    assert result_path.exists()
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    assert payload["error_code"] == "UNEXPECTED_ERROR"
    assert "traceback" in payload["details"]
    assert "RuntimeError" in payload["details"]["traceback"]


def test_parse_outputs_with_effective_aligner_attaches_validation(tmp_path):
    module = _load_skill_module()
    upstream = tmp_path / "upstream" / "results"
    (upstream / "pipeline_info").mkdir(parents=True)
    # Documented core files nf-core always writes (audit finding #6).
    (upstream / "pipeline_info" / "software_versions.yml").write_text(
        "x", encoding="utf-8"
    )
    (upstream / "pipeline_info" / "params_2026-01-01.json").write_text(
        "{}", encoding="utf-8"
    )
    (upstream / "multiqc").mkdir()
    (upstream / "multiqc" / "multiqc_report.html").write_text(
        "report", encoding="utf-8"
    )
    (upstream / "multiqc" / "multiqc_data").mkdir()
    fastqc = upstream / "fastqc"
    fastqc.mkdir()
    (fastqc / "sample_fastqc.html").write_text("html", encoding="utf-8")
    (fastqc / "sample_fastqc.zip").write_text("zip", encoding="utf-8")
    h5ad = upstream / "simpleaf" / "mtx_conversions" / "combined_matrix.h5ad"
    h5ad.parent.mkdir(parents=True)
    h5ad.write_text("h5ad", encoding="utf-8")
    args = module.argparse.Namespace(preset="standard", skip_multiqc=False)

    result = module._parse_outputs_with_effective_aligner(tmp_path, args)

    assert result["aligner_effective"] == "simpleaf"
    assert result["preferred_h5ad"] == str(h5ad)
    assert result["output_validation"] == {
        "missing_required": [],
        "missing_optional": [],
    }


def test_raise_if_expected_outputs_missing_raises_structured_error(tmp_path):
    module = _load_skill_module()
    parsed_outputs = {
        "aligner_effective": "simpleaf",
        "output_validation": {
            "missing_required": ["simpleaf"],
            "missing_optional": ["simpleaf/mtx_conversions/*.h5ad"],
        },
    }

    with pytest.raises(SkillError) as exc:
        module._raise_if_expected_outputs_missing(parsed_outputs, output_dir=tmp_path)

    assert exc.value.stage == "parsing"
    assert exc.value.error_code == "EXPECTED_OUTPUTS_NOT_FOUND"
    assert exc.value.details["missing_required"] == ["simpleaf"]


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


def test_build_extra_configs_emits_resource_limits_on_linux_real_run(tmp_path):
    """A non-macOS docker real run must get a host-scaled resourceLimits cap so the
    local executor does not abort when STARsolo's default request exceeds this host."""
    module = _load_skill_module()
    import argparse as _argparse

    args = _argparse.Namespace(profile="docker", demo=False)
    with patch("platform.system", return_value="Linux"):
        result = module._build_extra_nextflow_configs(args, tmp_path)
    names = [p.name for p in result]
    assert "resource_limits.config" in names
    # Shipped inside the bundle so commands.sh can re-apply it on replay.
    shipped = tmp_path / "reproducibility" / "resource_limits.config"
    assert shipped.is_file()
    text = shipped.read_text(encoding="utf-8")
    assert "resourceLimits" in text
    assert "stageInMode" not in text and "--platform" not in text


def test_build_extra_configs_skips_resource_limits_on_linux_demo(tmp_path):
    module = _load_skill_module()
    import argparse as _argparse

    args = _argparse.Namespace(profile="docker", demo=True)
    with patch("platform.system", return_value="Linux"):
        result = module._build_extra_nextflow_configs(args, tmp_path)
    assert "resource_limits.config" not in [p.name for p in result]


def test_resource_limits_memory_scales_with_host(monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "_host_memory_gb", lambda: 64)
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: None)
    assert module._resource_limits_memory_gb() > 15


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

    with (
        patch("subprocess.run", fake_run),
        patch.object(
            module, "_resolve_scrna_orchestrator", return_value=fake_orchestrator
        ),
    ):
        module._run_downstream_handoff(
            args, parsed_outputs={"preferred_h5ad": str(h5ad)}, output_dir=tmp_path
        )
    assert "cmd" in captured, "subprocess.run must be called when preferred_h5ad is set"
    cmd = captured["cmd"]
    assert any("scrna_orchestrator" in str(part) for part in cmd), (
        f"scrna_orchestrator.py must appear in the command, got: {cmd}"
    )


def test_run_downstream_handoff_warns_when_h5ad_selection_is_ambiguous(
    tmp_path, capsys
):
    module = _load_skill_module()
    import argparse as _argparse

    args = _argparse.Namespace(run_downstream=True, skip_downstream=False)
    parsed_outputs = {
        "preferred_h5ad": "",
        "h5ad_candidates": [
            str(tmp_path / "sampleA_filtered_matrix.h5ad"),
            str(tmp_path / "sampleB_filtered_matrix.h5ad"),
        ],
    }

    module._run_downstream_handoff(
        args, parsed_outputs=parsed_outputs, output_dir=tmp_path
    )

    captured = capsys.readouterr()
    assert "ambiguous" in captured.err.lower()
    assert "h5ad" in captured.err.lower()
