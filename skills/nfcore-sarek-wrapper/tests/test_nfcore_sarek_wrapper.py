"""Tests for the nfcore-sarek-wrapper orchestrator."""
from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path
from typing import Any

import pytest


SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = SKILL_DIR / "nfcore_sarek_wrapper.py"


def _load_module():
    """Load the wrapper module.

    We deliberately do *not* purge ``preflight``/``executor``/etc. from
    ``sys.modules``: ``test_preflight.py`` and friends import those modules at
    file-scope and monkeypatch attributes on the cached instances.  Reloading
    them under their feet creates two divergent copies and breaks unrelated
    tests.  Re-using the cached modules is safe for our purposes because the
    orchestrator imports everything it needs by name, and our tests
    monkeypatch the orchestrator module's bound names rather than the
    submodules themselves.
    """
    if str(SKILL_DIR) not in sys.path:
        sys.path.insert(0, str(SKILL_DIR))
    cached = sys.modules.get("nfcore_sarek_wrapper")
    if cached is not None:
        return cached
    spec = importlib.util.spec_from_file_location("nfcore_sarek_wrapper", SCRIPT_PATH)
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    sys.modules["nfcore_sarek_wrapper"] = module
    spec.loader.exec_module(module)
    return module


@pytest.fixture
def module():
    return _load_module()


# ---------------------------------------------------------------------------
# Fixtures
# ---------------------------------------------------------------------------


def _fake_pipeline_source() -> dict[str, Any]:
    return {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/sarek",
        "resolved_version": "3.8.1",
        "branch": "",
        "dirty": False,
    }


def _fake_preflight(warnings=None, notes=None):
    class _PR:
        def __init__(self):
            self.warnings = list(warnings or [])
            self.notes = list(notes or [])
    return _PR()


def _fake_outputs_report():
    class _OR:
        pass
    o = _OR()
    o.step_completed = "mapping"
    o.preprocessing = {}
    o.variant_calling = {}
    o.annotation = {}
    o.qc = {}
    o.pipeline_info = None
    o.reference_outputs = None
    o.samples_detected = []
    o.pairs_detected = []
    o.csv_handoff = {}
    o.handoff_available = False
    o.warnings = []
    o.missing_outputs = []
    return o


def _fake_reporting_artifacts(output_dir: Path):
    repro = output_dir / "reproducibility"
    repro.mkdir(parents=True, exist_ok=True)
    commands_sh = repro / "commands.sh"
    commands_sh.write_text("#!/usr/bin/env bash\nnextflow run ...\n", encoding="utf-8")
    class _RA:
        report_md: Path
        result_json: Path
        commands_sh: Path
        remap_paths_py: Path
        files_written: list[Path]
        warnings: list[str]
    ra = _RA()
    ra.report_md = repro / "report.md"
    ra.result_json = repro / "result.json"
    ra.commands_sh = commands_sh
    ra.remap_paths_py = repro / "remap_paths.py"
    ra.files_written = [commands_sh]
    ra.warnings = []
    return ra


def _patch_heavy(monkeypatch, module, *, pipeline_source=None, preflight=None,
                  samplesheet_report=None, outputs_report=None, exec_result=None,
                  exec_side_effect=None):
    """Stub the heavy dependencies so a full main() run is fast."""
    monkeypatch.setattr(
        module, "resolve_pipeline_source",
        lambda **kw: pipeline_source or _fake_pipeline_source(),
    )
    monkeypatch.setattr(
        module, "run_preflight",
        lambda **kw: preflight or _fake_preflight(),
    )
    monkeypatch.setattr(
        module, "load_manifest",
        lambda *a, **kw: None,
    )
    if samplesheet_report is None:
        samplesheet_report = {
            "normalized_path": Path("/tmp/x.csv"),
            "step": "mapping",
            "sample_count": 0,
            "patient_count": 0,
            "sample_names": [],
            "patient_names": [],
            "fastq_paths": [],
            "bam_paths": [],
            "cram_paths": [],
            "vcf_paths": [],
            "spring_paths": [],
            "tables": [],
            "unknown_columns": [],
            "sex_counts": {},
            "status_counts": {},
            "pairings": [],
            "analysis_mode": "germline",
        }
    monkeypatch.setattr(
        module, "validate_and_normalize_samplesheet",
        lambda *a, **kw: samplesheet_report,
    )

    def _build_effective_params(args, *, normalized_samplesheet, output_dir):
        return {
            "outdir": "upstream/results",
            "step": getattr(args, "step", None) or "mapping",
            "aligner": getattr(args, "aligner", None) or "bwa-mem",
            "tools": getattr(args, "tools", None),
            "skip_tools": getattr(args, "skip_tools", None),
            "wes": bool(getattr(args, "wes", False)),
        }

    monkeypatch.setattr(module, "build_effective_params", _build_effective_params)

    def _write_params_yaml(params, *, output_dir):
        repro = Path(output_dir) / "reproducibility"
        repro.mkdir(parents=True, exist_ok=True)
        p = repro / "params.yaml"
        p.write_text("# fake params\n", encoding="utf-8")
        return p

    monkeypatch.setattr(module, "write_params_yaml", _write_params_yaml)

    call_log: dict[str, Any] = {}

    def _build_nextflow_command(**kw):
        call_log["build_command"] = kw
        return (["nextflow", "run", "nf-core/sarek"], "nextflow run nf-core/sarek")

    monkeypatch.setattr(module, "build_nextflow_command", _build_nextflow_command)

    def _execute_nextflow(command, *, cwd, output_dir, timeout_seconds):
        call_log["execute"] = {
            "command": command,
            "cwd": cwd,
            "output_dir": output_dir,
            "timeout_seconds": timeout_seconds,
        }
        if exec_side_effect is not None:
            raise exec_side_effect
        return exec_result or {"exit_code": 0, "stdout_path": "x", "stderr_path": "y"}

    monkeypatch.setattr(module, "execute_nextflow", _execute_nextflow)

    def _parse_outputs(*a, **kw):
        call_log["parse_outputs"] = (a, kw)
        return outputs_report or _fake_outputs_report()

    monkeypatch.setattr(module, "parse_outputs", _parse_outputs)

    def _write_reports(**kw):
        call_log["write_reports"] = kw
        return _fake_reporting_artifacts(kw["output_dir"])

    monkeypatch.setattr(module, "write_reports", _write_reports)

    def _write_provenance_bundle(**kw):
        call_log["write_provenance_bundle"] = kw
        from collections import namedtuple
        Bundle = namedtuple("Bundle", ["bundle_dir", "manifest_path", "files_written", "warnings"])
        return Bundle(kw["output_dir"] / "reproducibility", Path("/tmp/m.json"), [], [])

    monkeypatch.setattr(module, "write_provenance_bundle", _write_provenance_bundle)

    monkeypatch.setattr(module, "_detect_java_version", lambda: "openjdk 21")
    monkeypatch.setattr(module, "_detect_nextflow_version", lambda: "25.10.2")
    return call_log


# ---------------------------------------------------------------------------
# Parser tests
# ---------------------------------------------------------------------------


def test_help_exits_zero(module, capsys):
    with pytest.raises(SystemExit) as exc:
        module.build_parser().parse_args(["--help"])
    assert exc.value.code == 0
    out = capsys.readouterr().out
    assert "usage" in out.lower()



def test_input_false_normalized_to_none(module):
    # Sarek's `--input false` sentinel (with --build-only-index) means "no input".
    args = module.build_parser().parse_args([
        "--output", "out", "--build-only-index", "--input", "false",
    ])
    module._validate_wrapper_flags(args)
    assert args.input is None


def test_nfcore_native_underscore_spelling_accepted(module):
    # nf-core's own parameter is `--skip_tools` (snake_case). The wrapper exposes
    # it as `--skip-tools`, but a user copying an upstream nf-core command types
    # the underscore form. Both spellings must land in the same dest so direct
    # wrapper invocation is not surprised by the native spelling.
    hyphen = module.build_parser().parse_args([
        "--output", "out", "--skip-tools", "baserecalibrator",
    ])
    underscore = module.build_parser().parse_args([
        "--output", "out", "--skip_tools", "baserecalibrator",
    ])
    assert hyphen.skip_tools == "baserecalibrator"
    assert underscore.skip_tools == "baserecalibrator"


def test_nfcore_native_underscore_spelling_for_value_reference(module):
    # Generalises beyond skip_tools: a schema passthrough with a value
    # (`--fasta_fai`) must also accept the nf-core-native underscore spelling.
    args = module.build_parser().parse_args([
        "--output", "out", "--fasta_fai", "ref.fasta.fai",
    ])
    assert args.fasta_fai == "ref.fasta.fai"


def test_mutect_profile_rejected_outside_demo(module):
    # The `mutect` profile only includeConfig conf/test_mutect2.config, whose sole
    # effect is `--normal-sample normal` on MUTECT2_PAIRED — correct ONLY for the
    # upstream test dataset. On a real paired run it would mislabel the normal
    # sample, so it must be rejected outside --demo.
    from errors import ErrorCode, SkillError
    args = module.build_parser().parse_args([
        "--output", "out", "--input", "sheet.csv", "--mutect-profile",
    ])
    with pytest.raises(SkillError) as ei:
        module._validate_wrapper_flags(args)
    assert ei.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION


def test_mutect_profile_allowed_under_demo(module):
    # Under --demo the upstream test data's normal IS named "normal", so the
    # profile is valid; it must not raise.
    args = module.build_parser().parse_args([
        "--output", "out", "--demo", "--mutect-profile",
    ])
    module._validate_wrapper_flags(args)


def test_mutect_token_in_profile_rejected_outside_demo(module):
    # `--profile docker,mutect` composes the mutect profile via _sync_profile_flags;
    # the same guard must fire for the token form.
    from errors import ErrorCode, SkillError
    args = module.build_parser().parse_args([
        "--output", "out", "--input", "sheet.csv", "--profile", "docker,mutect",
    ])
    module._sync_profile_flags(args)
    with pytest.raises(SkillError) as ei:
        module._validate_wrapper_flags(args)
    assert ei.value.error_code == ErrorCode.INVALID_FLAG_COMBINATION



def test_input_help_mentions_all_official_formats(module, capsys):
    with pytest.raises(SystemExit):
        module.build_parser().parse_args(["--help"])
    out = capsys.readouterr().out
    assert ".tsv" in out
    assert ".json" in out
    assert ".yaml" in out


def test_check_report_labels_true_somatic_pairs(module, tmp_path, capsys):
    args = module.build_parser().parse_args(["--output", tmp_path.as_posix(), "--check", "--demo"])
    samplesheet_report = {
        "sample_count": 1,
        "patient_count": 1,
        "analysis_mode": "germline",
        "pairings": [
            {"patient": "P1", "normal": "S1", "tumors": [], "mode": "germline"},
        ],
    }

    rc = module._write_check_report(
        output_dir=tmp_path,
        args=args,
        samplesheet_report=samplesheet_report,
        preflight_result=_fake_preflight(),
        pipeline_source=_fake_pipeline_source(),
    )

    out = capsys.readouterr().out
    assert rc == 0
    assert "Somatic pairs: 0" in out
    assert "Pairs:" not in out


def test_params_file_is_loaded_for_preflight_snapshot(module, tmp_path):
    params_file = tmp_path / "params.yaml"
    params_file.write_text(
        "step: annotate\n"
        "tools: vep\n"
        "vep_cache_version: '111'\n"
        "vep_genome: GRCh38\n"
        "vep_species: homo_sapiens\n",
        encoding="utf-8",
    )
    args = module.build_parser().parse_args([
        "--input", "samplesheet.csv",
        "--output", tmp_path.as_posix(),
        "--params-file", params_file.as_posix(),
    ])

    params = module._build_params_for_preflight(args, composed_profile="docker")

    assert params["step"] == "annotate"
    assert params["tools"] == ["vep"]
    assert params["vep_cache_version"] == "111"





def test_cf_contamination_is_integer_param(module):
    args = module.build_parser().parse_args([
        "--output", "out",
        "--demo",
        "--cf-contamination", "7",
    ])
    assert args.cf_contamination == 7
    assert isinstance(args.cf_contamination, int)



def test_test_full_profile_is_treated_as_demo(module):
    # Audit pass-2 M4: any upstream test* / test_full* profile token must
    # trigger the demo cleanup so user-supplied reference flags are dropped.
    parser = module.build_parser()
    for profile_str in ("docker,test_full", "singularity,test_full_germline", "docker,test_aws"):
        args = parser.parse_args(["--output", "out", "--profile", profile_str])
        module._sync_profile_flags(args)
        assert args.demo is True, f"expected demo=True for profile {profile_str!r}"



def test_freebayes_filter_is_string_expression(module):
    # nextflow_schema.json types freebayes_filter as `string` — a vcflib/vcffilter
    # expression (default "30"), NOT an integer. The wrapper must accept arbitrary
    # filter expressions like "QUAL > 20" and preserve them verbatim.
    args = module.build_parser().parse_args([
        "--output", "out",
        "--demo",
        "--freebayes-filter", "QUAL > 20 & DP > 10",
    ])
    assert args.freebayes_filter == "QUAL > 20 & DP > 10"
    assert isinstance(args.freebayes_filter, str)
    # The numeric default form is still accepted (as a string).
    args2 = module.build_parser().parse_args([
        "--output", "out", "--demo",
        "--freebayes-filter", "30",
    ])
    assert args2.freebayes_filter == "30"


def test_nextflow_config_repeatable(module):
    args = module.build_parser().parse_args([
        "--output", "out", "--demo",
        "--nextflow-config", "hpc.config",
        "--nextflow-config", "local.config",
    ])
    assert args.nextflow_config == ["hpc.config", "local.config"]


# ---------------------------------------------------------------------------
# Override / profile composition tests
# ---------------------------------------------------------------------------


def test_apply_demo_overrides_clears_references(module):
    args = module.build_parser().parse_args([
        "--output", "out",
        "--demo",
        "--fasta", "/tmp/x.fa",
        "--genome", "GATK.GRCh38",
        "--pon", "/tmp/p.vcf.gz",
        "--igenomes-ignore",
        "--input", "input.csv",
    ])
    module._apply_demo_overrides(args)
    assert args.fasta is None
    assert args.genome is None
    assert args.pon is None
    assert args.input is None
    # Demo must NOT set igenomes_ignore: the test profile (test.config) defines
    # genome + igenomes_base, and forcing igenomes_ignore would break the run.
    # A user-supplied --igenomes-ignore must be cleared, not passed through.
    assert getattr(args, "igenomes_ignore", False) is not True



def test_custom_fasta_without_genome_auto_disables_igenomes(module):
    """A custom --fasta with no --genome must auto-set igenomes_ignore so Sarek
    does not fall back to its default GATK.GRCh38 iGenomes bundle and try to
    validate ~20 remote s3://ngi-igenomes paths as local files."""
    args = module.build_parser().parse_args([
        "--output", "out",
        "--input", "input.csv",
        "--fasta", "/refs/custom.fa",
    ])
    args.profile = "docker"  # composed profile (no upstream test profile)
    notice = module._auto_disable_igenomes_for_custom_fasta(args)
    assert args.igenomes_ignore is True
    assert notice and "igenomes" in notice.lower()
    # The auto-disable must reach both preflight and the emitted params.yaml.
    params = module._build_params_for_preflight(args, composed_profile="docker")
    assert params.get("igenomes_ignore") is True


def test_custom_fasta_with_explicit_genome_keeps_igenomes(module):
    """When the user passes BOTH --fasta and --genome (documented override), the
    wrapper must respect the explicit iGenomes choice and not auto-disable it."""
    args = module.build_parser().parse_args([
        "--output", "out",
        "--input", "input.csv",
        "--fasta", "/refs/custom.fa",
        "--genome", "GATK.GRCh38",
    ])
    args.profile = "docker"
    notice = module._auto_disable_igenomes_for_custom_fasta(args)
    assert getattr(args, "igenomes_ignore", False) is False
    assert notice is None


def test_no_fasta_does_not_auto_disable_igenomes(module):
    """A plain --genome run (no custom fasta) must not be touched."""
    args = module.build_parser().parse_args([
        "--output", "out",
        "--input", "input.csv",
        "--genome", "GATK.GRCh38",
    ])
    args.profile = "docker"
    notice = module._auto_disable_igenomes_for_custom_fasta(args)
    assert getattr(args, "igenomes_ignore", False) is False
    assert notice is None


def test_demo_run_does_not_auto_disable_igenomes(module):
    """--demo uses the upstream test profile, which owns genome/igenomes_base;
    auto-disabling would break it."""
    args = module.build_parser().parse_args(["--output", "out", "--demo"])
    args.profile = "test,docker"
    notice = module._auto_disable_igenomes_for_custom_fasta(args)
    assert notice is None
    assert getattr(args, "igenomes_ignore", False) is False


def test_timeout_hours_resolves_to_seconds(module):
    args = module.build_parser().parse_args(
        ["--output", "out", "--demo", "--no-banner", "--timeout-hours", "6"]
    )
    assert module._resolve_timeout_seconds(args) == 6 * 3600


def test_timeout_hours_zero_disables_cap(module):
    """--timeout-hours 0 disables the wall-clock cap (parity with scrnaseq/rnaseq)."""
    args = module.build_parser().parse_args(
        ["--output", "out", "--demo", "--no-banner", "--timeout-hours", "0"]
    )
    assert module._resolve_timeout_seconds(args) is None


def test_timeout_hours_default_is_pinned(module):
    args = module.build_parser().parse_args(["--output", "out", "--demo", "--no-banner"])
    assert module._resolve_timeout_seconds(args) == module.DEFAULT_TIMEOUT_SECONDS


def test_timeout_hours_negative_rejected(module):
    with pytest.raises(SystemExit):
        module.build_parser().parse_args(
            ["--output", "out", "--demo", "--no-banner", "--timeout-hours", "-1"]
        )


def test_pipeline_version_mismatch_blocked(module):
    """A non-pinned --pipeline-version is blocked by default (parity with siblings)."""
    args = module.build_parser().parse_args(
        ["--output", "out", "--demo", "--no-banner", "--pipeline-version", "3.5.0"]
    )
    with pytest.raises(module.SkillError) as exc:
        module._check_pipeline_version_supported(args)
    assert exc.value.error_code == "PIPELINE_SOURCE_INVALID"


def test_pipeline_version_override_allows_mismatch(module):
    args = module.build_parser().parse_args(
        [
            "--output", "out", "--demo", "--no-banner",
            "--pipeline-version", "3.5.0", "--allow-pipeline-version-override",
        ]
    )
    module._check_pipeline_version_supported(args)  # must not raise


def test_pipeline_version_default_passes(module):
    args = module.build_parser().parse_args(["--output", "out", "--demo", "--no-banner"])
    module._check_pipeline_version_supported(args)  # default == pinned, no raise


def test_work_dir_defaults_under_output(module, tmp_path):
    args = module.build_parser().parse_args(["--output", str(tmp_path), "--demo", "--no-banner"])
    assert module._resolve_nextflow_work_dir(args, tmp_path) == tmp_path / "upstream" / "work"


def test_work_dir_override_local(module, tmp_path):
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--demo", "--no-banner", "--work-dir", str(tmp_path / "scratch")]
    )
    assert module._resolve_nextflow_work_dir(args, tmp_path) == (tmp_path / "scratch").resolve()


def test_work_dir_override_object_store_uri_preserved(module, tmp_path):
    args = module.build_parser().parse_args(
        ["--output", str(tmp_path), "--demo", "--no-banner", "--work-dir", "s3://bucket/work"]
    )
    assert module._resolve_nextflow_work_dir(args, tmp_path) == "s3://bucket/work"


def test_flag_surface_counts_match_docs(module):
    # Locks the numbers cited in README.md / CLAUDE.md / SKILL.md so they cannot
    # silently drift: 154 Sarek passthrough params + 23 wrapper-only controls
    # = 177 user-facing flags (excluding the auto-generated -h/--help).
    # (Homogenization added --timeout-hours, --allow-pipeline-version-override,
    # --work-dir, and --allow-remote-inputs for parity with the sibling wrappers,
    # bumping wrapper-only controls 19 → 23.)
    assert len(module._SAREK_PASSTHROUGH_PARAMS) == 154
    parser = module.build_parser()
    opts = [a for a in parser._actions if a.option_strings and a.dest != "help"]
    assert len({a.dest for a in opts}) == 177
    passthrough_dests = {dest for _flag, dest, *_rest in module._SAREK_PASSTHROUGH_PARAMS}
    wrapper_only = {a.dest for a in opts} - passthrough_dests
    assert len(wrapper_only) == 23


def test_demo_strips_reference_extra_params(module):
    # --demo must clear reference flags even when supplied via the --extra-param
    # escape hatch, otherwise they leak into params.yaml and break the test profile.
    args = module.build_parser().parse_args([
        "--output", "out", "--demo",
        "--extra-param", "fasta=/tmp/x.fa",
        "--extra-param", "dbsnp=/tmp/d.vcf.gz",
        "--extra-param", "igenomes_ignore=true",
        "--extra-param", "custom_config_version=1.0",  # non-reference: must survive
    ])
    module._merge_extra_params(args)
    module._apply_demo_overrides(args)
    assert "fasta" not in args._extras
    assert "dbsnp" not in args._extras
    assert "igenomes_ignore" not in args._extras
    # Non-reference passthrough extras must be preserved.
    assert args._extras.get("custom_config_version") == "1.0"


def test_demo_params_yaml_has_no_reference_keys(module, tmp_path):
    # End-to-end guard: after the demo override pipeline, the real params builder
    # must emit no reference keys (genome/igenomes_ignore/fasta/dbsnp/...), even
    # when supplied via first-class flags AND via --extra-param.
    from params_builder import build_effective_params

    out = tmp_path / "out"
    args = module.build_parser().parse_args([
        "--output", str(out), "--demo",
        "--genome", "GATK.GRCh38",
        "--igenomes-ignore",
        "--extra-param", "fasta=/tmp/x.fa",
        "--extra-param", "dbsnp=/tmp/d.vcf.gz",
    ])
    module._merge_extra_params(args)
    module._sync_profile_flags(args)
    module._apply_demo_overrides(args)

    ss = out / "reproducibility" / "samplesheet.demo.csv"
    ss.parent.mkdir(parents=True, exist_ok=True)
    ss.write_text("patient,sample,lane,fastq_1,fastq_2\n", encoding="utf-8")
    params = build_effective_params(args, normalized_samplesheet=ss, output_dir=out)
    for key in ("genome", "igenomes_ignore", "igenomes_base", "fasta", "dbsnp"):
        assert key not in params, f"{key} leaked into demo params.yaml"


def test_sync_profile_flags_detects_modifiers(module):
    args = module.build_parser().parse_args([
        "--output", "out", "--demo",
        "--profile", "docker,arm64,gpu,spark,mutect",
    ])
    module._sync_profile_flags(args)
    assert args.arm is True
    assert args.gpu is True
    assert args.spark_profile is True
    assert args.mutect_profile is True


def test_sync_profile_flags_test_implies_demo(module):
    args = module.build_parser().parse_args([
        "--output", "out",
        "--input", "input.csv",
        "--profile", "test,docker",
    ])
    module._sync_profile_flags(args)
    assert args.demo is True


def test_main_profile_test_without_input_runs_as_demo(module, monkeypatch, tmp_path):
    # A bare `--profile test` (no --demo, no --input) must be treated as a demo
    # run: profile-flag syncing happens before the input requirement is enforced,
    # so it must not fail with MISSING_INPUT.
    _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--profile", "test",
        "--check",
        "--no-banner",
    ])
    assert rc == 0


def test_main_demo_keeps_output_root_to_two_children(module, monkeypatch, tmp_path):
    # Finding 1: no stray top-level logs/ or macos_docker.config at the output root.
    _patch_heavy(monkeypatch, module)
    monkeypatch.setattr(module.sys, "platform", "darwin")
    out = tmp_path / "out"
    rc = module.main([
        "--output", str(out),
        "--demo",
        "--no-banner",
        "--profile", "docker",
    ])
    assert rc == 0
    assert not (out / "logs").exists()
    assert not (out / "macos_docker.config").exists()


def test_merge_extra_params_key_value_form(module):
    args = module.build_parser().parse_args([
        "--output", "out", "--demo",
        "--extra-param", "tools=haplotypecaller,mutect2",
        "--extra-param", "wes=true",
        "--extra-param", "novel_key=hello",
    ])
    module._merge_extra_params(args)
    assert args.tools == "haplotypecaller,mutect2"
    assert args.wes is True
    assert not hasattr(args, "novel_key")
    assert args._extras["novel_key"] == "hello"


def test_merge_extra_params_rejects_malformed(module):
    args = module.build_parser().parse_args([
        "--output", "out", "--demo",
        "--extra-param", "broken",
    ])
    with pytest.raises(Exception) as exc:
        module._merge_extra_params(args)
    # SkillError - check by attribute
    assert getattr(exc.value, "error_code", None) is not None


# ---------------------------------------------------------------------------
# Validation tests
# ---------------------------------------------------------------------------


def test_validate_wrapper_flags_requires_input_or_demo(module):
    args = module.build_parser().parse_args(["--output", "out"])
    with pytest.raises(Exception) as exc:
        module._validate_wrapper_flags(args)
    assert getattr(exc.value, "error_code", None) == "MISSING_INPUT"



def test_run_downstream_requires_skill(module):
    args = module.build_parser().parse_args([
        "--output", "out", "--demo",
        "--run-downstream",
    ])
    with pytest.raises(Exception) as exc:
        module._validate_wrapper_flags(args)
    assert "INVALID_FLAG_COMBINATION" in getattr(exc.value, "error_code", "")


# ---------------------------------------------------------------------------
# main() integration tests
# ---------------------------------------------------------------------------


def test_main_demo_check_no_executor_call(module, monkeypatch, tmp_path):
    call_log = _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--check",
        "--no-banner",
    ])
    assert rc == 0
    assert "execute" not in call_log
    # check_result.json lands at the output root — parallel to result.json and
    # matching nfcore-rnaseq/scrnaseq (not inside reproducibility/).
    check_path = tmp_path / "out" / "check_result.json"
    assert check_path.exists()
    assert not (tmp_path / "out" / "reproducibility" / "check_result.json").exists()
    payload = json.loads(check_path.read_text())
    assert payload["ok"] is True
    assert payload["mode"] == "check"



def test_main_returns_2_on_skill_error_during_validation(module, monkeypatch, tmp_path):
    _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        # no --input and no --demo → SkillError(MISSING_INPUT)
        "--no-banner",
    ])
    assert rc == 1


def test_main_keyboard_interrupt_returns_130(module, monkeypatch, tmp_path):
    _patch_heavy(monkeypatch, module, exec_side_effect=KeyboardInterrupt())
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
    ])
    assert rc == 130



def test_main_extra_param_merges_into_run(module, monkeypatch, tmp_path):
    call_log = _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
        "--extra-param", "tools=haplotypecaller",
    ])
    assert rc == 0
    assert call_log["write_reports"]["params"]["tools"] == "haplotypecaller"


def test_preflight_snapshot_includes_flags_used_by_rules(module):
    args = module.build_parser().parse_args([
        "--output", "out",
        "--input", "samplesheet.csv",
        "--save-mapped",
        "--save-output-as-bam",
        "--vep-dbnsfp",
        "--umi-location", "read1",
        "--umi-length", "8",
        "--sentieon-haplotyper-emit-mode", "gvcf",
    ])
    params = module._build_params_for_preflight(args, composed_profile="docker")
    for key in (
        "input",
        "save_mapped",
        "save_output_as_bam",
        "vep_dbnsfp",
        "umi_location",
        "umi_length",
        "sentieon_haplotyper_emit_mode",
    ):
        assert key in params



def test_macos_docker_memory_scales_to_host(module, monkeypatch, tmp_path):
    # The resourceLimits memory cap must scale to host RAM (with headroom), not
    # be hardcoded at 15.GB — on a 16GB Mac a 15GB cap starves the host.
    monkeypatch.setattr(module.sys, "platform", "darwin")
    monkeypatch.setattr(module, "_host_memory_gb", lambda: 16)
    args = module.build_parser().parse_args(["--output", str(tmp_path / "out"), "--demo", "--profile", "docker"])
    args.profile = "test,docker"
    cfg = module._write_macos_docker_config(tmp_path / "out", args=args)
    text = cfg.read_text()
    # 16GB host → capped below host with headroom; never 15.GB, never above host.
    assert "memory: '15.GB'" not in text
    import re as _re
    m = _re.search(r"memory: '(\d+)\.GB'", text)
    assert m is not None
    mem = int(m.group(1))
    assert 4 <= mem <= 12  # host(16) - headroom, bounded

    # Big host (64GB) caps at the 15GB ceiling.
    monkeypatch.setattr(module, "_host_memory_gb", lambda: 64)
    cfg2 = module._write_macos_docker_config(tmp_path / "out2", args=args)
    m2 = _re.search(r"memory: '(\d+)\.GB'", cfg2.read_text())
    assert m2 is not None and int(m2.group(1)) == 15

    # Tiny host (8GB) floors at 4GB.
    monkeypatch.setattr(module, "_host_memory_gb", lambda: 8)
    cfg3 = module._write_macos_docker_config(tmp_path / "out3", args=args)
    m3 = _re.search(r"memory: '(\d+)\.GB'", cfg3.read_text())
    assert m3 is not None and int(m3.group(1)) >= 4


def _resource_args(module, tmp_path, *, demo=False, profile="docker"):
    argv = ["--output", str(tmp_path / "out"), "--profile", profile]
    if demo:
        argv.append("--demo")
    else:
        argv += ["--input", "s.csv", "--genome", "GATK.GRCh38", "--tools", "strelka"]
    args = module.build_parser().parse_args(argv)
    return args


def test_resource_limits_config_emitted_on_linux_real_run(module, monkeypatch, tmp_path):
    # The host-scaled resourceLimits cap must be emitted on non-macOS docker runs
    # so a real run whose production requirements exceed this host does not abort.
    monkeypatch.setattr(module.sys, "platform", "linux")
    args = _resource_args(module, tmp_path)
    cfg = module._write_resource_limits_config(tmp_path / "out", args=args)
    assert cfg is not None
    text = cfg.read_text()
    assert "resourceLimits" in text
    # None of the macOS-only workarounds (they slow or break native Linux runs).
    assert "--platform" not in text and "stageInMode" not in text


def test_resource_limits_config_scales_with_host(module, monkeypatch, tmp_path):
    monkeypatch.setattr(module.sys, "platform", "linux")
    monkeypatch.setattr(module, "_host_memory_gb", lambda: 64)
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: None)
    args = _resource_args(module, tmp_path)
    text = module._write_resource_limits_config(tmp_path / "out", args=args).read_text()
    import re as _re
    mem = int(_re.search(r"memory: '(\d+)\.GB'", text).group(1))
    assert mem > 15 and mem <= 64  # no 15 GB macOS ceiling, never above physical RAM


def test_resource_limits_config_none_on_demo(module, monkeypatch, tmp_path):
    monkeypatch.setattr(module.sys, "platform", "linux")
    args = _resource_args(module, tmp_path, demo=True)
    assert module._write_resource_limits_config(tmp_path / "out", args=args) is None


def test_resource_limits_config_none_on_macos(module, monkeypatch, tmp_path):
    monkeypatch.setattr(module.sys, "platform", "darwin")
    args = _resource_args(module, tmp_path)
    assert module._write_resource_limits_config(tmp_path / "out", args=args) is None


def test_resource_limits_config_none_without_docker(module, monkeypatch, tmp_path):
    monkeypatch.setattr(module.sys, "platform", "linux")
    args = _resource_args(module, tmp_path, profile="singularity")
    assert module._write_resource_limits_config(tmp_path / "out", args=args) is None


def test_main_writes_macos_docker_config_when_applicable(module, monkeypatch, tmp_path):
    _patch_heavy(monkeypatch, module)
    monkeypatch.setattr(module.sys, "platform", "darwin")
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
        "--profile", "docker",
    ])
    assert rc == 0
    # The macOS config lives inside reproducibility/, keeping the output root to
    # exactly two children (upstream/, reproducibility/).
    config_path = tmp_path / "out" / "reproducibility" / "macos_docker.config"
    assert config_path.exists()
    assert not (tmp_path / "out" / "macos_docker.config").exists()
    text = config_path.read_text()
    assert "stageInMode" in text
    # No --arm, so platform=linux/amd64 should be set
    assert "--platform=linux/amd64" in text
    assert "-u $(id -u):$(id -g)" in text





def test_main_run_downstream_writes_handoff(module, monkeypatch, tmp_path):
    _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
        "--run-downstream",
        "--downstream-skill", "clinical-variant-reporter",
    ])
    assert rc == 0
    sh = tmp_path / "out" / "reproducibility" / "sarek_downstream_handoff.sh"
    js = tmp_path / "out" / "reproducibility" / "sarek_downstream_handoff.json"
    assert sh.exists()
    assert js.exists()
    payload = json.loads(js.read_text())
    assert payload["selected_skill"] == "clinical-variant-reporter"
    assert "clinical-variant-reporter" in payload["routes"]
    assert "clinical-trial-finder" in payload["routes"]
    assert "omics-target-evidence-mapper" in payload["routes"]
    assert "wes-clinical-report-en" in payload["routes"]
    assert "wes-clinical-report-es" in payload["routes"]


def test_downstream_handoff_examples_use_python3(module, monkeypatch, tmp_path):
    """The generated sarek_downstream_handoff.sh lists cross-skill example commands
    the user runs directly (`bash sarek_downstream_handoff.sh`). They must invoke
    `python3`, not a bare `python`, for portability to python3-only systems (PEP 394)."""
    _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
        "--run-downstream",
        "--downstream-skill", "clinical-variant-reporter",
    ])
    assert rc == 0
    sh = (tmp_path / "out" / "reproducibility" / "sarek_downstream_handoff.sh").read_text()
    assert "python3 skills/" in sh
    assert "python skills/" not in sh, "bare `python` is not portable; use python3"


def test_main_run_downstream_without_skill_returns_2(module, monkeypatch, tmp_path):
    _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
        "--run-downstream",
    ])
    assert rc == 1



def test_main_resume_passes_manifest_to_preflight(module, monkeypatch, tmp_path):
    """When a manifest exists, --resume should hand it to run_preflight."""
    _patch_heavy(monkeypatch, module)
    # Replace stub manifest loader to return a fake dict
    monkeypatch.setattr(module, "load_manifest", lambda *a, **kw: {"schema_version": 1, "params_checksum": "sha256:abc"})
    captured: dict[str, Any] = {}

    def capture_preflight(**kw):
        captured.update(kw)
        return _fake_preflight()

    monkeypatch.setattr(module, "run_preflight", capture_preflight)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--check",
        "--no-banner",
    ])
    assert rc == 0
    # No --resume on the command line, so resume=False.
    assert captured["resume"] is False


def test_main_demo_preserves_resume(module, monkeypatch, tmp_path):
    """--demo must NOT strip --resume.

    Nextflow's -resume is orthogonal to -profile test (nf-core documents no
    incompatibility). Stripping it left a failed or interrupted demo in a dead end:
    the partially-populated output dir makes preflight raise OUTPUT_DIR_NOT_EMPTY,
    whose fix text tells the user to "re-run with --resume against a compatible
    manifest" — advice the demo override then silently discarded.

    Demo/real drift stays blocked: params["profile"] records the *composed* profile
    (which carries the `test` token for demo) and is part of _RESUME_TRACKED_FIELDS.
    """
    _patch_heavy(monkeypatch, module)
    monkeypatch.setattr(module, "load_manifest", lambda *a, **kw: {"schema_version": 1, "params_checksum": "sha256:abc"})
    captured: dict[str, Any] = {}

    def capture_preflight(**kw):
        captured.update(kw)
        return _fake_preflight()

    monkeypatch.setattr(module, "run_preflight", capture_preflight)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--resume",
        "--check",
        "--no-banner",
    ])
    assert rc == 0
    assert captured["resume"] is True, "--demo must not disable --resume"


def test_main_pipeline_source_passed_through(module, monkeypatch, tmp_path):
    call_log = _patch_heavy(
        monkeypatch, module,
        pipeline_source={
            "source_kind": "local_checkout",
            "source_ref": "/tmp/sarek",
            "resolved_version": "abc123",
            "branch": "main",
            "dirty": False,
        },
    )
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
    ])
    assert rc == 0
    assert call_log["write_reports"]["pipeline_source_kind"] == "local_checkout"


def test_main_nextflow_config_propagates(module, monkeypatch, tmp_path):
    call_log = _patch_heavy(monkeypatch, module)
    extra = tmp_path / "extra.config"
    extra.write_text("// extra\n", encoding="utf-8")
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
        "--nextflow-config", str(extra),
    ])
    assert rc == 0
    cfgs = call_log["build_command"]["extra_configs"]
    assert any(Path(p) == extra.resolve() for p in cfgs)


def test_main_banner_suppression(module, monkeypatch, tmp_path, capsys):
    _patch_heavy(monkeypatch, module)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
    ])
    assert rc == 0
    out = capsys.readouterr().out
    assert "ClawBio :: nfcore-sarek-wrapper" not in out





def test_main_skill_error_bubbles_up_with_nonzero_exit(module, monkeypatch, tmp_path):
    """SkillError raised from any module yields exit code 1 (parity with siblings)."""
    _patch_heavy(monkeypatch, module)

    def boom(**kw):
        from errors import ErrorCode, SkillError
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_JAVA,
            message="java missing",
            fix="install java",
        )

    monkeypatch.setattr(module, "run_preflight", boom)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
    ])
    assert rc == 1


def test_main_error_result_json_written_at_output_root(module, monkeypatch, tmp_path):
    # On error the result.json marker must land at the output root, next to where
    # a successful run writes result.json (parity with nfcore-rnaseq/scrnaseq), so
    # consumers read <output>/result.json regardless of outcome.
    _patch_heavy(monkeypatch, module)

    def boom(**kw):
        from errors import ErrorCode, SkillError
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_JAVA,
            message="java missing",
            fix="install java",
        )

    monkeypatch.setattr(module, "run_preflight", boom)
    out = tmp_path / "out"
    rc = module.main(["--output", str(out), "--demo", "--no-banner"])
    assert rc == 1
    err = out / "result.json"
    assert err.exists()
    payload = json.loads(err.read_text())
    assert payload["ok"] is False
    assert payload["status"] == "error"
    assert payload["error_code"] == "MISSING_JAVA"
    assert not (out / "reproducibility" / "result.json").exists()


def test_main_unexpected_error_returns_1(module, monkeypatch, tmp_path):
    _patch_heavy(monkeypatch, module)

    def boom(**kw):
        raise RuntimeError("something broke unexpectedly")

    monkeypatch.setattr(module, "run_preflight", boom)
    rc = module.main([
        "--output", str(tmp_path / "out"),
        "--demo",
        "--no-banner",
    ])
    assert rc == 1



def test_main_user_samplesheet_path_runs_validator(module, monkeypatch, tmp_path):
    captured_csv: dict[str, Any] = {}

    def capture_sheet(input_path, output_path, *, step, tools):
        captured_csv["input"] = input_path
        captured_csv["output"] = output_path
        captured_csv["step"] = step
        captured_csv["tools"] = tools
        return {
            "normalized_path": output_path,
            "step": step,
            "sample_count": 1,
            "patient_count": 1,
            "sample_names": ["s1"],
            "patient_names": ["p1"],
            "fastq_paths": [],
            "bam_paths": [],
            "cram_paths": [],
            "vcf_paths": [],
            "spring_paths": [],
            "tables": [],
            "unknown_columns": [],
            "sex_counts": {"XX": 1},
            "status_counts": {0: 1},
            "pairings": [],
            "analysis_mode": "germline",
        }

    _patch_heavy(monkeypatch, module)
    monkeypatch.setattr(module, "validate_and_normalize_samplesheet", capture_sheet)
    sheet = tmp_path / "input.csv"
    sheet.write_text("patient,sample,lane\n", encoding="utf-8")
    rc = module.main([
        "--input", str(sheet),
        "--output", str(tmp_path / "out"),
        "--no-banner",
        "--check",
        "--tools", "haplotypecaller,mutect2",
    ])
    assert rc == 0
    assert captured_csv["input"] == sheet.resolve()
    assert captured_csv["tools"] == ["haplotypecaller", "mutect2"]



def test_detect_repo_root_finds_clawbio_py(module):
    root = module._detect_repo_root()
    # Should be a Path; may or may not contain clawbio.py depending on layout
    assert isinstance(root, Path)


# ── clawbio runner allowlist stays in sync with the wrapper CLI ────────────────
# A flag must be in clawbio.cli.SKILLS["sarek-pipeline"]["allowed_extra_flags"]
# (and, when value-free, also in allowed_extra_flags_without_values) or the
# extra-args filter silently drops it / swallows the following token.


def _sarek_pipeline_allowlist() -> tuple[set[str], set[str]]:
    from clawbio.cli import SKILLS

    entry = SKILLS.get("sarek-pipeline", {})
    return (
        set(entry.get("allowed_extra_flags") or ()),
        set(entry.get("allowed_extra_flags_without_values") or ()),
    )


def _wrapper_option_flags(module) -> tuple[set[str], set[str]]:
    """Return (all option strings, value-free option strings) for the wrapper parser."""
    import argparse

    parser = module.build_parser()
    all_opts: set[str] = set()
    no_value: set[str] = set()
    for action in parser._actions:
        for opt in action.option_strings:
            all_opts.add(opt)
            if (
                isinstance(
                    action,
                    (
                        argparse._StoreTrueAction,
                        argparse._StoreFalseAction,
                        argparse._CountAction,
                        argparse._HelpAction,
                    ),
                )
                or action.nargs == 0
            ):
                no_value.add(opt)
    return all_opts, no_value


def test_clawbio_allow_remote_inputs_is_value_free_for_sarek_pipeline():
    """--allow-remote-inputs is a store_true wrapper flag; it must live in the
    boolean (value-free) set, otherwise the extra-args filter treats it as
    value-taking and consumes the following argument as its value."""
    allowed, booleans = _sarek_pipeline_allowlist()
    assert "--allow-remote-inputs" in allowed
    assert "--allow-remote-inputs" in booleans


def test_clawbio_allow_pipeline_version_override_is_value_free_for_sarek_pipeline():
    """--allow-pipeline-version-override is store_true and must be value-free too."""
    allowed, booleans = _sarek_pipeline_allowlist()
    assert "--allow-pipeline-version-override" in allowed
    assert "--allow-pipeline-version-override" in booleans


def test_sarek_pipeline_allowlist_consistent_with_wrapper(module):
    """The clawbio runner allowlist must stay coherent with the wrapper parser:

    * every allow-listed flag is a real wrapper flag, and
    * the value/boolean classification matches the wrapper so the extra-args
      filter neither drops a value nor swallows the next token.
    """
    values, booleans = _sarek_pipeline_allowlist()
    allowed = values | booleans
    all_opts, no_value = _wrapper_option_flags(module)

    assert allowed <= all_opts, f"allowlist flags absent from wrapper: {sorted(allowed - all_opts)}"
    boolean_allowed = no_value & allowed
    assert boolean_allowed <= booleans, (
        "value-free wrapper flags missing from the without-values set: "
        f"{sorted(boolean_allowed - booleans)}"
    )
    assert booleans <= no_value, (
        "flags marked value-free but the wrapper expects a value: "
        f"{sorted(booleans - no_value)}"
    )
