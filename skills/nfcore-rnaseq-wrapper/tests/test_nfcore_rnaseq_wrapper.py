from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from argparse import Namespace
from pathlib import Path
from unittest.mock import patch

import pytest


SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = SKILL_DIR / "nfcore_rnaseq_wrapper.py"
PROJECT_ROOT = SKILL_DIR.parent.parent


def _rnaseq_pipeline_allowlist() -> tuple[set[str], set[str]]:
    """Return (value flags, boolean flags) allow-listed for `clawbio run rnaseq-pipeline`.

    The clawbio skill registry lives in ``clawbio.cli.SKILLS``. A flag must be in the
    allowlist or it is dropped by the extra-args filter before reaching the wrapper.
    """
    from clawbio.cli import SKILLS

    entry = SKILLS.get("rnaseq-pipeline", {})
    return (
        set(entry.get("allowed_extra_flags") or ()),
        set(entry.get("allowed_extra_flags_without_values") or ()),
    )


def _load_skill_module():
    spec = importlib.util.spec_from_file_location("nfcore_rnaseq_wrapper_module", SCRIPT_PATH)
    assert spec is not None
    assert spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in (
        "command_builder",
        "errors",
        "executor",
        "outputs_parser",
        "params_builder",
        "pipeline_source",
        "preflight",
        "provenance",
        "reporting",
        "samplesheet_builder",
        "schemas",
    ):
        loaded = sys.modules.get(name)
        loaded_file = Path(getattr(loaded, "__file__", "") or "")
        if loaded is not None and (loaded_file == SKILL_DIR / f"{name}.py" or SKILL_DIR in loaded_file.parents):
            sys.modules.pop(name, None)
    if str(SKILL_DIR) in sys.path:
        sys.path.remove(str(SKILL_DIR))
    return module


def _fake_pipeline_source(version: str = "3.26.0") -> dict[str, object]:
    return {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/rnaseq",
        "resolved_version": version,
        "branch": "",
        "dirty": False,
    }


def _fake_preflight() -> dict[str, object]:
    return {
        "ok": True,
        "java": {"version": "21.0.0", "path": "/usr/bin/java"},
        "nextflow": {"version": "25.04.3", "path": "/usr/local/bin/nextflow"},
        "profile": {"profile": "docker", "backend_ready": True},
        "pipeline_source": _fake_pipeline_source(),
        "references": {},
        "aligner_effective": "star_salmon",
        "handoff_available": True,
        "warnings": [],
        "samplesheet": {"sample_count": 1, "sample_names": ["sampleA"], "unknown_columns": []},
    }


def _write_samplesheet(tmp_path: Path) -> Path:
    r1 = tmp_path / "sampleA_R1.fastq.gz"
    r2 = tmp_path / "sampleA_R2.fastq.gz"
    r1.write_text("x", encoding="utf-8")
    r2.write_text("x", encoding="utf-8")
    samplesheet = tmp_path / "samplesheet.csv"
    samplesheet.write_text(
        f"sample,fastq_1,fastq_2,strandedness\nsampleA,{r1},{r2},auto\n",
        encoding="utf-8",
    )
    return samplesheet


def test_entrypoint_file_exists():
    assert SCRIPT_PATH.exists()


def test_parser_accepts_core_flags():
    module = _load_skill_module()
    args = module.build_parser().parse_args(
        [
            "--input",
            "samplesheet.csv",
            "--output",
            "out",
            "--demo",
            "--check",
            "--profile",
            "docker",
            "--pipeline-version",
            "3.26.0",
            "--pipeline-local",
            "/tmp/rnaseq",
            "--resume",
        ]
    )
    assert args.input == "samplesheet.csv"
    assert args.output == "out"
    assert args.demo is True
    assert args.check is True
    assert args.profile == "docker"
    assert args.pipeline_version == "3.26.0"
    assert args.pipeline_local == "/tmp/rnaseq"
    assert args.resume is True


@pytest.mark.parametrize(
    ("flag", "dest"),
    [
        ("--fasta", "fasta"),
        ("--genome", "genome"),
    ],
)
def test_parser_exposes_reference_flags(flag, dest):
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", flag, "/ref/file"])
    assert getattr(args, dest) == "/ref/file"


@pytest.mark.parametrize(
    ("flag", "value"),
    [
        ("--aligner", "bad_aligner"),
    ],
)
def test_parser_rejects_invalid_choices(flag, value):
    module = _load_skill_module()
    with pytest.raises(SystemExit):
        module.build_parser().parse_args(["--output", "out", flag, value])


def test_parser_defaults_match_plan():
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out"])
    assert args.profile == "docker"
    assert args.pipeline_version == "3.26.0"
    # --aligner defaults to None (sentinel); star_salmon is applied after all
    # override functions run via _apply_aligner_default in _run_wrapper.
    assert args.aligner is None
    assert args.trimmer == "trimgalore"
    assert args.pseudo_aligner is None
    assert args.pseudo_aligner_kmer_size is None
    assert args.ribo_removal_tool is None
    assert args.stranded_threshold is None
    assert args.unstranded_threshold is None
    assert args.deseq2_vst is None


def test_parser_accepts_nextflow_config_flag():
    """--nextflow-config must be accepted by the argument parser and stored as a list."""
    module = _load_skill_module()
    args = module.build_parser().parse_args([
        "--output", "out",
        "--nextflow-config", "/path/to/hpc.config",
    ])
    assert hasattr(args, "nextflow_config"), "Parser must expose --nextflow-config as args.nextflow_config"
    assert "/path/to/hpc.config" in (args.nextflow_config or [])


def test_prokaryotic_coerces_default_aligner_to_bowtie2_salmon():
    """--prokaryotic with no explicit aligner must coerce to bowtie2_salmon."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--prokaryotic"])
    assert args.aligner is None  # sentinel default — no explicit aligner given
    module._apply_prokaryotic_overrides(args)
    assert args.aligner == "bowtie2_salmon"


def test_prokaryotic_does_not_coerce_explicit_aligner():
    """--prokaryotic with an explicit non-default aligner leaves the aligner unchanged."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--prokaryotic", "--aligner", "hisat2"])
    module._apply_prokaryotic_overrides(args)
    assert args.aligner == "hisat2"


def test_sync_profile_flags_sets_prokaryotic_from_profile_string():
    """_sync_profile_flags must set args.prokaryotic=True when 'prokaryotic' is in --profile."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--profile", "prokaryotic"])
    assert args.prokaryotic is False  # before sync
    module._sync_profile_flags(args)
    assert args.prokaryotic is True


@pytest.mark.parametrize("profile", [
    "docker,test",
])
def test_sync_profile_flags_sets_noinput_for_self_contained_profiles(profile):
    """_sync_profile_flags must set args._noinput=True for self-contained nf-core profiles.

    Official nf-core/rnaseq profiles (test, test_full, test_prokaryotic, etc.) ship
    with params.input in their profile config.  Running them without --input is valid
    upstream; the wrapper must not block them with MISSING_INPUT.
    _sync_profile_flags detects these profiles and sets a _noinput sentinel so
    _prepare_samplesheet and params_builder skip the input requirement.
    """
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--profile", profile])
    module._sync_profile_flags(args)
    assert getattr(args, "_noinput", False) is True, (
        f"_sync_profile_flags did not set _noinput=True for --profile {profile!r}"
    )


def test_sync_profile_flags_test_prokaryotic_sets_both_prokaryotic_and_noinput():
    """test_prokaryotic must trigger both _noinput (no input required) and prokaryotic (aligner coercion)."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--profile", "test_prokaryotic"])
    module._sync_profile_flags(args)
    assert getattr(args, "_noinput", False) is True, "test_prokaryotic must set _noinput"
    assert getattr(args, "prokaryotic", False) is True, "test_prokaryotic must set prokaryotic"


def test_sync_profile_flags_sets_rapid_quant_from_profile_string():
    """_sync_profile_flags must set args.rapid_quant=True when 'rapid_quant' is in --profile."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--profile", "rapid_quant"])
    assert args.rapid_quant is False  # before sync
    module._sync_profile_flags(args)
    assert args.rapid_quant is True


def test_sync_profile_flags_sets_arm_from_arm64_token():
    """--profile arm64 must set args.arm=True so params.yaml gets arm: true."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--profile", "arm64"])
    assert args.arm is False  # before sync
    module._sync_profile_flags(args)
    assert args.arm is True


def test_apply_demo_overrides_clears_all_reference_flags(capsys):
    """--demo must clear every reference/index arg before they reach params.yaml.

    The upstream `test` profile bundles sample FASTQs paired with its own
    reference data. Leaving even one reference flag (e.g. --fasta) would let
    params-file override the profile's matched ref and silently desynchronise
    samples from references — producing garbage counts with no error. This
    test pins the contract: every flag in _DEMO_CLEARED_REFERENCE_FIELDS that
    the user set must be None after _apply_demo_overrides, and the user must
    receive a single structured warning naming the cleared flags.
    """
    module = _load_skill_module()
    args = module.build_parser().parse_args([
        "--output", "out",
        "--demo",
        "--fasta", "/tmp/custom.fa",
        "--gtf", "/tmp/custom.gtf",
        "--genome", "GRCh38",
        "--star-index", "/tmp/star",
        "--salmon-index", "/tmp/salmon",
        "--transcript-fasta", "/tmp/tx.fa",
    ])
    module._apply_demo_overrides(args)
    for field in module._DEMO_CLEARED_REFERENCE_FIELDS:
        assert getattr(args, field, None) is None, (
            f"--demo did not clear {field!r}; partial overrides desync samples from refs"
        )
    captured = capsys.readouterr()
    assert "--demo ignores reference flags" in captured.err
    # The warning should name the specific flags the user actually set.
    for cleared in ("--fasta", "--gtf", "--genome", "--star-index", "--salmon-index", "--transcript-fasta"):
        assert cleared in captured.err, f"warning should list {cleared}"


def test_apply_demo_overrides_silent_when_no_refs_set(capsys):
    """--demo without any reference flags must not print the ref-clearing warning."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--demo"])
    module._apply_demo_overrides(args)
    captured = capsys.readouterr()
    assert "--demo ignores reference flags" not in captured.err


def test_apply_demo_overrides_preserves_resume(capsys):
    """--demo must NOT strip --resume.

    Nextflow's -resume is orthogonal to -profile test: nf-core documents no
    incompatibility, and a demo run keeps its work dir (<output>/upstream/work) and
    Nextflow session cache (<output>/.nextflow) in the output dir, so resuming it in
    place is well-defined. Stripping --resume here was what made a demo bundle's
    reproducibility/commands.sh unreplayable: the replay hit OUTPUT_DIR_NOT_EMPTY and
    the error's suggested --resume fix was silently discarded.
    """
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--demo", "--resume"])
    module._apply_demo_overrides(args)
    assert args.resume is True, "--demo must not disable --resume"
    captured = capsys.readouterr()
    assert "disables --resume" not in captured.err


def test_debug_profile_does_not_set_noinput():
    """--profile debug must NOT set _noinput=True.

    The debug profile only enables debug logging (dumpHashes, cleanup=false).
    It does not include params.input and still requires --input from the user.
    """
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--profile", "debug"])
    module._sync_profile_flags(args)
    assert getattr(args, "_noinput", False) is False, (
        "--profile debug must not set _noinput — it carries no params.input"
    )


def test_prepare_samplesheet_skips_input_requirement_when_noinput(tmp_path):
    """_prepare_samplesheet must not raise MISSING_INPUT when args._noinput is True.

    This is the core of the P2 fix: when a self-contained test profile is used
    without --input, _prepare_samplesheet should behave like demo mode (no user
    samplesheet) instead of raising SkillError(MISSING_INPUT).
    """
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", str(tmp_path), "--profile", "test_full"])
    args._noinput = True
    args.demo = False
    staging_dir = tmp_path / "staging"
    staging_dir.mkdir()
    # Must not raise SkillError
    normalized, staged, summary = module._prepare_samplesheet(args, tmp_path, staging_dir=staging_dir)
    assert normalized is not None
    assert summary["sample_count"] == 0  # no user rows, same as demo mode


def test_noinput_samplesheet_uses_noinput_filename(tmp_path):
    """_prepare_samplesheet with _noinput=True must NOT produce samplesheet.demo.csv.

    The stub written for self-contained test profiles should be named
    samplesheet.noinput.csv so that provenance audits can distinguish it
    from a real --demo run.
    """
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", str(tmp_path), "--profile", "test_full"])
    args._noinput = True
    args.demo = False
    staging_dir = tmp_path / "staging"
    staging_dir.mkdir()

    normalized, staged, summary = module._prepare_samplesheet(args, tmp_path, staging_dir=staging_dir)
    assert "demo" not in normalized.name, (
        f"samplesheet for _noinput must not be named samplesheet.demo.csv; got {normalized.name}"
    )
    assert "noinput" in normalized.name, (
        f"samplesheet for _noinput should contain 'noinput' in the name; got {normalized.name}"
    )


def test_demo_samplesheet_still_uses_demo_filename(tmp_path):
    """--demo must still produce samplesheet.demo.csv (no regression)."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", str(tmp_path), "--demo"])
    staging_dir = tmp_path / "staging"
    staging_dir.mkdir()

    normalized, staged, summary = module._prepare_samplesheet(args, tmp_path, staging_dir=staging_dir)
    assert "demo" in normalized.name, (
        f"--demo samplesheet should be named samplesheet.demo.csv; got {normalized.name}"
    )


def test_noinput_does_not_force_igenomes_ignore_in_params_builder(tmp_path):
    """params_builder must NOT write igenomes_ignore=True solely because _noinput is set.

    test_full* profiles use params.genome='GRCh37' (iGenomes); setting
    igenomes_ignore=True would prevent nf-core from resolving that genome.
    The profile itself controls igenomes_ignore — the wrapper must stay silent.
    """
    import importlib.util as _ilu
    pb_spec = _ilu.spec_from_file_location("params_builder_mod", SKILL_DIR / "params_builder.py")
    pb = _ilu.module_from_spec(pb_spec)
    pb_spec.loader.exec_module(pb)

    module = _load_skill_module()
    fake_samplesheet = tmp_path / "reproducibility" / "samplesheet.noinput.csv"
    fake_samplesheet.parent.mkdir(parents=True, exist_ok=True)
    fake_samplesheet.write_text("", encoding="utf-8")

    args = module.build_parser().parse_args(["--output", str(tmp_path), "--profile", "test_full"])
    args._noinput = True
    args.aligner = "star_salmon"  # _apply_aligner_default has not run yet

    params = pb._build_base_params(args, normalized_samplesheet=fake_samplesheet, output_dir=tmp_path)
    assert "igenomes_ignore" not in params, (
        f"_build_base_params must not set igenomes_ignore for _noinput (non-demo) runs; got {params}"
    )
    assert "input" not in params, (
        "_build_base_params must not set params.input for _noinput runs"
    )


# ── Audit F1: self-contained test profiles own their refs/tuning (like --demo) ─


def test_noinput_profile_clears_reference_overrides(capsys):
    """A self-contained test profile (--profile test_full) ships its own params.input
    AND bundled reference data. User --genome/--fasta/index flags must be cleared and
    warned about — exactly like --demo — so they cannot silently override the profile's
    matched reference and desync samples from refs (audit F1)."""
    module = _load_skill_module()
    args = module.build_parser().parse_args([
        "--output", "out",
        "--profile", "test_full",
        "--genome", "GRCh38",
        "--fasta", "/tmp/custom.fa",
        "--gtf", "/tmp/custom.gtf",
        "--salmon-index", "/tmp/salmon",
    ])
    module._record_aligner_explicit(args)
    module._apply_demo_overrides(args)   # not demo -> no-op
    module._sync_profile_flags(args)     # sets _noinput=True for test_full
    module._apply_noinput_overrides(args)
    for field in ("genome", "fasta", "gtf", "salmon_index"):
        assert getattr(args, field) is None, (
            f"a self-contained test profile must clear {field!r}; partial overrides desync refs"
        )
    err = capsys.readouterr().err
    assert "test profile" in err.lower()
    assert "--genome" in err and "--fasta" in err


def test_noinput_params_builder_omits_reference_and_tuning_overrides(tmp_path):
    """params.yaml for a self-contained test profile must carry neither user references
    nor user tuning flags; the hermetic profile owns every pipeline parameter (audit F1)."""
    import importlib.util as _ilu
    pb_spec = _ilu.spec_from_file_location("params_builder_mod_f1", SKILL_DIR / "params_builder.py")
    pb = _ilu.module_from_spec(pb_spec)
    pb_spec.loader.exec_module(pb)

    module = _load_skill_module()
    args = module.build_parser().parse_args([
        "--output", str(tmp_path),
        "--profile", "test_full",
        "--genome", "GRCh38",
        "--with-umi", "--umitools-bc-pattern", "NNNN",
    ])
    module._record_aligner_explicit(args)
    module._apply_demo_overrides(args)
    module._sync_profile_flags(args)
    module._apply_noinput_overrides(args)
    module._apply_aligner_default(args)

    fake = tmp_path / "reproducibility" / "samplesheet.noinput.csv"
    fake.parent.mkdir(parents=True, exist_ok=True)
    fake.write_text("", encoding="utf-8")

    params = pb.build_effective_params(args, normalized_samplesheet=fake, output_dir=tmp_path)
    assert "genome" not in params, "test profile owns references; --genome must not reach params.yaml"
    assert "with_umi" not in params, "test profile owns tuning; --with-umi must not reach params.yaml"


def test_demo_still_clears_refs_after_noinput_helper_added(capsys):
    """Regression: --demo (which does not set _noinput) must keep clearing references."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--demo", "--genome", "GRCh38"])
    module._record_aligner_explicit(args)
    module._apply_demo_overrides(args)
    module._sync_profile_flags(args)
    module._apply_noinput_overrides(args)
    assert args.genome is None
    assert "--demo ignores reference flags" in capsys.readouterr().err


def test_check_mode_writes_check_result_and_does_not_execute(tmp_path, monkeypatch):
    module = _load_skill_module()
    executed = []
    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: _fake_pipeline_source(kw["requested_version"]))
    monkeypatch.setattr(module, "run_preflight", lambda *a, **kw: _fake_preflight())
    monkeypatch.setattr(module, "execute_nextflow", lambda *a, **kw: executed.append(True))
    rc = module.main(["--output", str(tmp_path), "--demo", "--check"])
    assert rc == 0
    assert not executed
    payload = json.loads((tmp_path / "check_result.json").read_text(encoding="utf-8"))
    assert payload["ok"] is True
    assert payload["skill"] == "nfcore-rnaseq-wrapper"


def test_check_demo_passes_with_mocked_environment(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: _fake_pipeline_source())
    monkeypatch.setattr(module, "run_preflight", lambda *a, **kw: _fake_preflight())
    rc = module.main(["--output", str(tmp_path), "--demo", "--check", "--profile", "docker"])
    assert rc == 0
    assert (tmp_path / "reproducibility" / "samplesheet.demo.csv").exists()


def test_run_execution_mode_builds_params_before_command(tmp_path, monkeypatch):
    module = _load_skill_module()
    order = []
    captured_success = {}
    params_path = tmp_path / "reproducibility" / "params.yaml"
    params_path.parent.mkdir()

    monkeypatch.setattr(module, "build_effective_params", lambda *a, **kw: order.append("params") or {"outdir": "x", "aligner": "star_salmon"})
    monkeypatch.setattr(module, "check_resume_params_checksum", lambda *a, **kw: order.append("params_checksum"))
    monkeypatch.setattr(module, "check_resume_samplesheet_checksum", lambda *a, **kw: order.append("samplesheet_checksum"))
    monkeypatch.setattr(module, "write_params_yaml", lambda *a, **kw: order.append("write_params") or params_path)
    monkeypatch.setattr(module, "execute_nextflow", lambda *a, **kw: order.append("execute") or {"returncode": 0})
    monkeypatch.setattr(module, "parse_outputs", lambda *a, **kw: {"preferred_counts_tsv": "counts.tsv", "pipeline_info_dir": "pipeline_info", "handoff_available": True})
    monkeypatch.setattr(
        module,
        "_write_success_outputs",
        lambda *a, **kw: order.append("success_outputs") or captured_success.update(parsed_outputs=kw["parsed_outputs"]),
    )

    args = Namespace(
        aligner="star_salmon",
        demo=False,
        profile="docker",
        prokaryotic=False,
        resume=True,
        run_downstream=True,
        skip_downstream=False,
    )
    staged = tmp_path / "stage.csv"
    staged.write_text("ok", encoding="utf-8")
    normalized = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    rc = module._run_execution_mode(
        args,
        output_dir=tmp_path,
        pipeline_source=_fake_pipeline_source(),
        preflight_result=_fake_preflight(),
        normalized_samplesheet=normalized,
        staged_samplesheet=staged,
        samplesheet_summary={"sample_count": 1},
    )
    assert rc == 0
    assert order[:5] == ["params", "params_checksum", "samplesheet_checksum", "write_params", "execute"]
    assert (tmp_path / "reproducibility" / "rnaseq_de_handoff.sh").exists()
    assert captured_success["parsed_outputs"]["rnaseq_de_handoff"].endswith("rnaseq_de_handoff.sh")
    assert "rnaseq_de_output_dir" not in captured_success["parsed_outputs"]
    assert captured_success["parsed_outputs"]["rnaseq_de_status"] == "template_only"


# ── Task 9: Mocked Nextflow end-to-end integration test ───────────────────


def _write_fake_star_salmon_artifacts(output_dir: Path) -> Path:
    counts = output_dir / "upstream" / "results" / "star_salmon" / "salmon.merged.gene_counts.tsv"
    counts.parent.mkdir(parents=True, exist_ok=True)
    counts.write_text("gene_id\tsampleA\nGENE1\t100\n", encoding="utf-8")
    multiqc = output_dir / "upstream" / "results" / "multiqc" / "star_salmon" / "multiqc_report.html"
    multiqc.parent.mkdir(parents=True, exist_ok=True)
    multiqc.write_text("<html>multiqc</html>", encoding="utf-8")
    pipeline_info = output_dir / "upstream" / "results" / "pipeline_info"
    pipeline_info.mkdir(parents=True, exist_ok=True)
    (pipeline_info / "execution_report.html").write_text("<html>exec</html>", encoding="utf-8")
    return counts


def test_mocked_nextflow_end_to_end_writes_all_expected_outputs(tmp_path, monkeypatch):
    module = _load_skill_module()
    counts_path = _write_fake_star_salmon_artifacts(tmp_path)

    def fake_execute_nextflow(command, cwd, output_dir, timeout_seconds):
        logs = Path(output_dir) / "logs"
        logs.mkdir(parents=True, exist_ok=True)
        stdout_log = logs / "nextflow.stdout.log"
        stderr_log = logs / "nextflow.stderr.log"
        stdout_log.write_text("Completed successfully", encoding="utf-8")
        stderr_log.write_text("", encoding="utf-8")
        return {"returncode": 0, "stdout_log": str(stdout_log), "stderr_log": str(stderr_log)}

    monkeypatch.setattr(module, "execute_nextflow", fake_execute_nextflow)
    monkeypatch.setattr(module, "time", Namespace(monotonic=iter([0.0, 5.123]).__next__), raising=False)

    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    samplesheet.parent.mkdir(parents=True, exist_ok=True)
    samplesheet.write_text("sample,fastq_1,strandedness\nsampleA,/data/R1.fastq.gz,auto\n", encoding="utf-8")

    args = Namespace(
        aligner="star_salmon",
        demo=True,
        profile="docker",
        prokaryotic=False,
        resume=False,
        run_downstream=False,
        skip_downstream=False,
        pseudo_aligner=None,
        pipeline_version="3.26.0",
        pipeline_local=None,
    )

    rc = module._run_execution_mode(
        args,
        output_dir=tmp_path,
        pipeline_source=_fake_pipeline_source(),
        preflight_result=_fake_preflight(),
        normalized_samplesheet=samplesheet,
        staged_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 1, "sample_names": ["sampleA"], "unknown_columns": []},
    )

    assert rc == 0
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "result.json").exists()
    result = json.loads((tmp_path / "result.json").read_text(encoding="utf-8"))
    assert result["summary"]["preferred_counts_tsv"] == str(counts_path)
    assert (tmp_path / "provenance" / "runtime.json").exists()
    runtime = json.loads((tmp_path / "provenance" / "runtime.json").read_text(encoding="utf-8"))
    assert runtime["duration_seconds"] == 5.123
    assert (tmp_path / "reproducibility" / "commands.sh").exists()


def test_run_downstream_handoff_template_only_without_required_flags(tmp_path):
    module = _load_skill_module()
    args = Namespace(run_downstream=True, skip_downstream=False, metadata=None, formula=None, contrast=None, downstream_output=None)
    with patch("subprocess.run") as run:
        result = module._run_downstream_handoff(
            args,
            parsed_outputs={"preferred_counts_tsv": "counts.tsv", "handoff_available": True},
            output_dir=tmp_path,
        )
    run.assert_not_called()
    assert result == {
        "template_path": str(tmp_path / "reproducibility" / "rnaseq_de_handoff.sh"),
        "downstream_output_dir": None,
        "downstream_status": "template_only",
        "downstream_returncode": None,
    }
    assert (tmp_path / "reproducibility" / "rnaseq_de_handoff.sh").exists()


def test_downstream_handoff_uses_portable_python_interpreter(tmp_path):
    """The generated `rnaseq_de_handoff.sh` is executed on possibly-fresh machines
    (`bash rnaseq_de_handoff.sh`) where only `python3` exists (PEP 394). It must invoke
    the interpreter portably as `${PYTHON:-python3}`, never a bare `python`, mirroring
    the `commands.sh` replay patch — a bare `python` fails with `python: command not found`."""
    module = _load_skill_module()
    args = Namespace(run_downstream=True, skip_downstream=False, metadata=None, formula=None, contrast=None, downstream_output=None)
    with patch("subprocess.run"):
        module._run_downstream_handoff(
            args,
            parsed_outputs={"preferred_counts_tsv": "counts.tsv", "handoff_available": True},
            output_dir=tmp_path,
        )
    script = (tmp_path / "reproducibility" / "rnaseq_de_handoff.sh").read_text(encoding="utf-8")
    assert "\npython " not in script, "bare `python` is not portable; use ${PYTHON:-python3}"
    assert '"${PYTHON:-python3}"' in script


def test_run_downstream_handoff_returns_none_without_run_downstream(tmp_path):
    """Without --run-downstream, NO handoff template is written and None is returned,
    even when counts are available (this is the demo/default path). The report's Next
    Steps still shows the suggested command; the rnaseq_de_handoff.sh file requires
    --run-downstream. Locks the behaviour SKILL.md documents."""
    module = _load_skill_module()
    args = Namespace(run_downstream=False, skip_downstream=False, metadata=None, formula=None, contrast=None, downstream_output=None)
    result = module._run_downstream_handoff(
        args, parsed_outputs={"preferred_counts_tsv": "counts.tsv", "handoff_available": True}, output_dir=tmp_path
    )
    assert result is None
    assert not (tmp_path / "reproducibility" / "rnaseq_de_handoff.sh").exists()


def test_run_downstream_handoff_returns_none_when_skip_downstream(tmp_path):
    """--skip-downstream suppresses the handoff template even with --run-downstream set."""
    module = _load_skill_module()
    args = Namespace(run_downstream=True, skip_downstream=True, metadata=None, formula=None, contrast=None, downstream_output=None)
    result = module._run_downstream_handoff(
        args, parsed_outputs={"preferred_counts_tsv": "counts.tsv", "handoff_available": True}, output_dir=tmp_path
    )
    assert result is None
    assert not (tmp_path / "reproducibility" / "rnaseq_de_handoff.sh").exists()


def test_run_downstream_handoff_launches_rnaseq_de_when_required_flags_present(tmp_path, monkeypatch):
    module = _load_skill_module()
    launched = {}
    args = Namespace(
        run_downstream=True,
        skip_downstream=False,
        metadata=str(tmp_path / "meta with spaces.csv"),
        formula="~ batch + condition",
        contrast="condition,treated,control",
        downstream_output=None,
    )

    def fake_run(cmd, capture_output, text, timeout, cwd):
        launched["cmd"] = cmd
        launched["capture_output"] = capture_output
        launched["text"] = text
        launched["timeout"] = timeout
        launched["cwd"] = cwd
        return subprocess.CompletedProcess(cmd, 0, stdout="ok", stderr="")

    monkeypatch.setattr(module.subprocess, "run", fake_run)
    result = module._run_downstream_handoff(
        args,
        parsed_outputs={"preferred_counts_tsv": str(tmp_path / "counts with spaces.tsv"), "handoff_available": True},
        output_dir=tmp_path,
    )
    assert result == {
        "template_path": str(tmp_path / "reproducibility" / "rnaseq_de_handoff.sh"),
        "downstream_output_dir": str(tmp_path / "rnaseq_de"),
        "downstream_status": "completed",
        "downstream_returncode": 0,
    }
    assert launched["cmd"][:4] == [sys.executable, str(PROJECT_ROOT / "clawbio.py"), "run", "rnaseq"]
    assert launched["cmd"] == [
        sys.executable,
        str(PROJECT_ROOT / "clawbio.py"),
        "run",
        "rnaseq",
        "--counts",
        str(tmp_path / "counts with spaces.tsv"),
        "--metadata",
        str(tmp_path / "meta with spaces.csv"),
        "--formula",
        "~ batch + condition",
        "--contrast",
        "condition,treated,control",
        "--output",
        str(tmp_path / "rnaseq_de"),
    ]
    assert launched["capture_output"] is True
    assert launched["text"] is True
    assert launched["timeout"] == 60 * 60 * 2
    assert launched["cwd"] == str(PROJECT_ROOT)


def test_downstream_handoff_uses_absolute_clawbio_path(tmp_path, monkeypatch):
    """_run_downstream_handoff must invoke clawbio.py via its absolute path."""
    module = _load_skill_module()
    import pathlib

    captured_commands = []

    def fake_run(cmd, **kwargs):
        captured_commands.append(list(cmd))
        return subprocess.CompletedProcess(cmd, returncode=0, stdout="", stderr="")

    monkeypatch.setattr(subprocess, "run", fake_run)

    # Write a stub clawbio.py so the existence check passes
    fake_clawbio = module._PROJECT_ROOT / "clawbio.py"
    already_exists = fake_clawbio.exists()
    if not already_exists:
        fake_clawbio.write_text("# stub\n", encoding="utf-8")

    try:
        counts = str(tmp_path / "counts.tsv")
        pathlib.Path(counts).write_text("gene\ts1\n", encoding="utf-8")
        meta = str(tmp_path / "meta.csv")
        pathlib.Path(meta).write_text("sample,condition\n", encoding="utf-8")

        args = Namespace(
            run_downstream=True,
            skip_downstream=False,
            metadata=meta,
            formula="~ condition",
            contrast="condition,treated,control",
            downstream_output=None,
        )
        parsed_outputs = {"preferred_counts_tsv": counts, "handoff_available": True}

        module._run_downstream_handoff(args, parsed_outputs=parsed_outputs, output_dir=tmp_path)
    finally:
        if not already_exists:
            fake_clawbio.unlink(missing_ok=True)

    assert captured_commands, "subprocess.run was never called"
    invoked = pathlib.Path(captured_commands[0][1])
    assert invoked.is_absolute(), (
        f"clawbio.py was invoked via a relative path: {captured_commands[0][1]!r}"
    )


def test_downstream_handoff_warns_when_clawbio_missing(tmp_path, monkeypatch, capsys):
    """_run_downstream_handoff must warn and return error status if clawbio.py is absent."""
    module = _load_skill_module()
    import pathlib

    fake_root = tmp_path / "fake_project"
    fake_root.mkdir()
    monkeypatch.setattr(module, "_PROJECT_ROOT", fake_root)
    # clawbio.py does NOT exist in fake_root

    counts = str(tmp_path / "counts.tsv")
    pathlib.Path(counts).write_text("gene\ts1\n", encoding="utf-8")
    meta = str(tmp_path / "meta.csv")
    pathlib.Path(meta).write_text("sample,condition\n", encoding="utf-8")

    args = Namespace(
        run_downstream=True,
        skip_downstream=False,
        metadata=meta,
        formula="~ condition",
        contrast="condition,treated,control",
        downstream_output=None,
    )
    parsed_outputs = {"preferred_counts_tsv": counts, "handoff_available": True}

    result = module._run_downstream_handoff(args, parsed_outputs=parsed_outputs, output_dir=tmp_path)

    captured = capsys.readouterr()
    assert "WARNING" in captured.err, f"Expected WARNING in stderr, got: {captured.err!r}"
    assert result["downstream_status"] == "error"


def test_resume_uses_staging_to_preserve_previous_artifacts(tmp_path, monkeypatch):
    module = _load_skill_module()
    previous = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    previous.parent.mkdir(parents=True)
    previous.write_text("previous\n", encoding="utf-8")
    samplesheet = _write_samplesheet(tmp_path)

    def fake_preflight(*a, **kw):
        raise module.SkillError("preflight", "INVALID_RESUME_STATE", "bad resume", "fix", {})

    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: _fake_pipeline_source())
    monkeypatch.setattr(module, "run_preflight", fake_preflight)
    rc = module.main(["--input", str(samplesheet), "--output", str(tmp_path), "--genome", "GRCh38", "--resume"])
    assert rc == 1
    assert previous.read_text(encoding="utf-8") == "previous\n"


def test_main_writes_json_skill_error_to_stderr(tmp_path, capsys, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(
        module,
        "resolve_pipeline_source",
        lambda **kw: (_ for _ in ()).throw(module.SkillError("preflight", "PIPELINE_SOURCE_INVALID", "bad", "fix", {})),
    )
    rc = module.main(["--output", str(tmp_path), "--demo"])
    assert rc == 1
    payload = json.loads(capsys.readouterr().err)
    assert payload["ok"] is False
    assert payload["error_code"] == "PIPELINE_SOURCE_INVALID"


def test_main_writes_json_unexpected_error_to_stderr(tmp_path, capsys, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "resolve_pipeline_source", lambda **kw: (_ for _ in ()).throw(RuntimeError("boom")))
    rc = module.main(["--output", str(tmp_path), "--demo"])
    assert rc == 1
    payload = json.loads(capsys.readouterr().err)
    assert payload["error_code"] == "UNEXPECTED_ERROR"
    assert payload["details"]["exception_type"] == "RuntimeError"


def test_macos_docker_config_has_audited_rnaseq_content(tmp_path):
    module = _load_skill_module()
    config_path = module._write_macos_docker_config(tmp_path)
    assert config_path == tmp_path / ".nextflow_macos_docker.config"
    content = config_path.read_text(encoding="utf-8")
    assert "stageInMode = 'copy'" in content
    assert "--platform linux/amd64" in content
    assert "containerOptions" in content
    assert "resourceLimits" in content
    assert "STARsolo" not in content
    assert "ext.args" not in content
    assert not (tmp_path / "reproducibility" / "macos_docker.config").exists()


# ── Audit follow-up F-4: the macOS Docker resourceLimits time ceiling must track
# --timeout-hours so a large-cohort run raised above the 12h default is not capped
# back to 12h by the generated config.


def test_macos_docker_config_time_defaults_to_twelve_hours(tmp_path):
    module = _load_skill_module()
    content = module._write_macos_docker_config(tmp_path).read_text(encoding="utf-8")
    assert "time: '12.h'" in content


def test_macos_docker_config_time_tracks_timeout_hours(tmp_path):
    module = _load_skill_module()
    content = module._write_macos_docker_config(tmp_path, timeout_hours=48).read_text(encoding="utf-8")
    assert "time: '48.h'" in content
    assert "time: '12.h'" not in content


def test_macos_docker_config_time_floored_at_one_hour(tmp_path):
    module = _load_skill_module()
    content = module._write_macos_docker_config(tmp_path, timeout_hours=0.5).read_text(encoding="utf-8")
    assert "time: '1.h'" in content


# ── Audit follow-up F-8: the per-process memory ceiling must never exceed the
# actual Docker VM memory (Docker Desktop's VM is usually smaller than the host,
# so a host-RAM-derived budget could be OOM-killed).


def test_macos_docker_memory_capped_to_docker_vm(monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "_detect_host_memory_gb", lambda: 64)  # host budget would hit the 15 GB cap
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: 6)     # but the VM only has 6 GB
    assert module._macos_docker_memory_gb() <= 6


def test_macos_docker_memory_falls_back_when_vm_unknown(monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "_detect_host_memory_gb", lambda: 64)
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: None)
    # Unknown VM size → preserve the prior host-derived behaviour (15 GB ceiling).
    assert module._macos_docker_memory_gb() == 15


def test_docker_vm_memory_none_when_docker_absent(monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module.shutil, "which", lambda name: None)
    assert module._docker_vm_memory_gb() is None


# ── Linux resourceLimits gap: the host-scaled resourceLimits config must be
# emitted on non-macOS docker runs too, not only on macOS. Without it a real
# (non-demo) run on a Linux host smaller than the pipeline's production process
# requirements aborts with "Process requirement exceeds available memory".


def test_resource_limits_config_has_portable_content(tmp_path):
    module = _load_skill_module()
    config_path = module._write_resource_limits_config(tmp_path)
    assert config_path == tmp_path / ".nextflow_resource_limits.config"
    content = config_path.read_text(encoding="utf-8")
    assert "resourceLimits" in content
    assert "memory:" in content and "cpus:" in content and "time:" in content
    # Must NOT carry the macOS-only workarounds (they slow or break Linux runs).
    assert "--platform linux/amd64" not in content
    assert "stageInMode" not in content
    assert "containerOptions" not in content


def test_resource_limits_memory_scales_with_host_no_macos_ceiling(monkeypatch):
    module = _load_skill_module()
    # A large Linux host must NOT be clamped to the 15 GB macOS Docker ceiling.
    monkeypatch.setattr(module, "_detect_host_memory_gb", lambda: 64)
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: None)
    mem = module._resource_limits_memory_gb()
    assert mem > module._MACOS_DOCKER_MEMORY_CEILING_GB
    assert mem <= 64  # never above physical RAM


def test_resource_limits_memory_capped_to_docker_vm(monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "_detect_host_memory_gb", lambda: 64)
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: 32)
    assert module._resource_limits_memory_gb() <= 32


def test_resource_limits_memory_fallback_when_undetectable(monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module, "_detect_host_memory_gb", lambda: None)
    monkeypatch.setattr(module, "_docker_vm_memory_gb", lambda: None)
    assert module._resource_limits_memory_gb() == module._MACOS_DOCKER_MEMORY_CEILING_GB


def _extra_config_args(module, **overrides):
    defaults = dict(profile="docker", arm=False, nextflow_config=None, demo=False, timeout_hours=12.0)
    defaults.update(overrides)
    return Namespace(**defaults)


def test_extra_configs_emits_resource_limits_on_linux_real_run(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module.platform, "system", lambda: "Linux")
    args = _extra_config_args(module)
    configs = module._build_extra_nextflow_configs(args, tmp_path)
    names = [p.name for p in configs]
    assert ".nextflow_resource_limits.config" in names
    assert ".nextflow_macos_docker.config" not in names


def test_extra_configs_skips_resource_limits_on_linux_demo(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module.platform, "system", lambda: "Linux")
    args = _extra_config_args(module, demo=True)
    configs = module._build_extra_nextflow_configs(args, tmp_path)
    # -profile test carries its own resourceLimits; do not override them.
    assert ".nextflow_resource_limits.config" not in [p.name for p in configs]


def test_extra_configs_macos_unchanged(tmp_path, monkeypatch):
    module = _load_skill_module()
    monkeypatch.setattr(module.platform, "system", lambda: "Darwin")
    args = _extra_config_args(module)
    names = [p.name for p in module._build_extra_nextflow_configs(args, tmp_path)]
    assert ".nextflow_macos_docker.config" in names
    assert ".nextflow_resource_limits.config" not in names


def test_rapid_quant_override_sets_profile_implied_args():
    """--rapid-quant must sync skip_alignment, pseudo_aligner, skip_quantification_merge."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--rapid-quant"])
    assert args.pseudo_aligner is None
    assert args.skip_alignment is False
    assert args.skip_quantification_merge is False
    module._apply_rapid_quant_overrides(args)
    assert args.pseudo_aligner == "salmon"
    assert args.skip_alignment is True
    assert args.skip_quantification_merge is True


def test_rapid_quant_override_respects_explicit_pseudo_aligner():
    """--rapid-quant must not override an explicit --pseudo-aligner."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out", "--rapid-quant", "--pseudo-aligner", "kallisto"])
    module._apply_rapid_quant_overrides(args)
    assert args.pseudo_aligner == "kallisto"


def test_rapid_quant_override_no_op_when_not_set():
    """_apply_rapid_quant_overrides must not touch args when --rapid-quant is not passed."""
    module = _load_skill_module()
    args = module.build_parser().parse_args(["--output", "out"])
    module._apply_rapid_quant_overrides(args)
    assert args.pseudo_aligner is None
    assert args.skip_alignment is False
    assert args.skip_quantification_merge is False


# ── Step 2: _write_success_outputs provenance isolation ──────────────────────


def test_write_success_outputs_calls_write_report_even_if_provenance_raises(tmp_path, monkeypatch):
    """If write_provenance_bundle raises, write_report and write_result must still be called."""
    module = _load_skill_module()
    calls = []
    monkeypatch.setattr(module, "write_repro_commands", lambda *a, **kw: None)
    monkeypatch.setattr(module, "write_provenance_bundle", lambda *a, **kw: (_ for _ in ()).throw(RuntimeError("disk full")))
    monkeypatch.setattr(module, "write_report", lambda *a, **kw: calls.append("write_report"))
    monkeypatch.setattr(module, "write_result", lambda *a, **kw: calls.append("write_result"))

    params = tmp_path / "reproducibility" / "params.yaml"
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    params.parent.mkdir(parents=True)
    params.write_text("{}", encoding="utf-8")
    samplesheet.write_text("x", encoding="utf-8")

    module._write_success_outputs(
        tmp_path,
        args=Namespace(aligner="star_salmon", profile="docker", demo=False,
                       resume=False, pseudo_aligner=None, prokaryotic=False),
        pipeline_source=_fake_pipeline_source(),
        preflight_result=_fake_preflight(),
        params_path=params,
        params_payload={"aligner": "star_salmon"},
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 1, "fastq_paths": []},
        parsed_outputs={"preferred_counts_tsv": "/counts.tsv"},
        execution_result={},
        command_str="nextflow run nf-core/rnaseq",
    )
    assert "write_report" in calls, "write_report must be called even when provenance fails"
    assert "write_result" in calls, "write_result must be called even when provenance fails"


def test_write_success_outputs_passes_provenance_warnings_to_write_report(tmp_path, monkeypatch):
    """A provenance failure must produce a post_run_warnings entry passed to write_report."""
    module = _load_skill_module()
    captured_warnings = {}
    monkeypatch.setattr(module, "write_repro_commands", lambda *a, **kw: None)
    monkeypatch.setattr(module, "write_provenance_bundle", lambda *a, **kw: (_ for _ in ()).throw(OSError("no space left")))
    monkeypatch.setattr(module, "write_report", lambda *a, **kw: captured_warnings.update(post_run_warnings=kw.get("post_run_warnings", [])))
    monkeypatch.setattr(module, "write_result", lambda *a, **kw: None)

    params = tmp_path / "reproducibility" / "params.yaml"
    samplesheet = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    params.parent.mkdir(parents=True)
    params.write_text("{}", encoding="utf-8")
    samplesheet.write_text("x", encoding="utf-8")

    module._write_success_outputs(
        tmp_path,
        args=Namespace(aligner="star_salmon", profile="docker", demo=False,
                       resume=False, pseudo_aligner=None, prokaryotic=False),
        pipeline_source=_fake_pipeline_source(),
        preflight_result=_fake_preflight(),
        params_path=params,
        params_payload={"aligner": "star_salmon"},
        normalized_samplesheet=samplesheet,
        samplesheet_summary={"sample_count": 1, "fastq_paths": []},
        parsed_outputs={"preferred_counts_tsv": "/counts.tsv"},
        execution_result={},
        command_str="nextflow run nf-core/rnaseq",
    )
    warnings = captured_warnings.get("post_run_warnings", [])
    assert len(warnings) >= 1, "Expected at least one post_run_warning from the provenance failure"
    assert any("provenance" in w.lower() or "OSError" in w or "no space" in w.lower() for w in warnings)
def test_record_aligner_explicit_true_when_user_supplied():
    module = _load_skill_module()
    ns = Namespace(aligner="star_rsem")
    module._record_aligner_explicit(ns)
    assert ns._aligner_explicit is True


def test_record_aligner_explicit_false_when_aligner_defaulted_none():
    module = _load_skill_module()
    ns = Namespace(aligner=None)
    module._record_aligner_explicit(ns)
    assert ns._aligner_explicit is False


def test_contaminant_screening_choices_derive_from_schemas_constant():
    """--contaminant-screening choices must come from the centralised schemas
    constant, not a hand-maintained inline literal (single source of truth)."""
    module = _load_skill_module()
    parser = module.build_parser()
    action = next(a for a in parser._actions if "--contaminant-screening" in a.option_strings)
    spec = importlib.util.spec_from_file_location("schemas_cs_check", SKILL_DIR / "schemas.py")
    schemas_mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(schemas_mod)
    assert set(action.choices) == schemas_mod.SUPPORTED_CONTAMINANT_SCREENING


# ── clawbio CLI integration: rnaseq-pipeline flags must reach the skill ────────
# `clawbio run rnaseq-pipeline --<flag>` forwards a flag only when it is allow-listed
# in clawbio.cli.SKILLS["rnaseq-pipeline"]; otherwise the extra-args filter drops it.
# These regression tests guard the flags the wrapper relies on.


def test_clawbio_nextflow_config_in_rnaseq_pipeline_allowlist():
    """--nextflow-config must be in the rnaseq-pipeline allowlist so it is not
    silently dropped before reaching the wrapper."""
    values, _ = _rnaseq_pipeline_allowlist()
    assert "--nextflow-config" in values, f"--nextflow-config missing from allowlist: {sorted(values)!r}"


def test_clawbio_nextflow_config_forwarded_to_rnaseq_pipeline():
    """--nextflow-config is a value-taking flag (action='append'); it must be in the
    value allowlist (not the boolean set) for every occurrence to reach the wrapper."""
    values, booleans = _rnaseq_pipeline_allowlist()
    assert "--nextflow-config" in values
    assert "--nextflow-config" not in booleans


def test_clawbio_star_index_forwarded_to_rnaseq_pipeline():
    """--star-index must be forwardable (regression: declared in the parser but not forwarded)."""
    values, _ = _rnaseq_pipeline_allowlist()
    assert "--star-index" in values, f"--star-index missing from allowlist: {sorted(values)!r}"


def test_clawbio_kallisto_index_forwarded_to_rnaseq_pipeline():
    """--kallisto-index must be forwardable (regression: absent from the forwarding loop)."""
    values, _ = _rnaseq_pipeline_allowlist()
    assert "--kallisto-index" in values, f"--kallisto-index missing from allowlist: {sorted(values)!r}"


def test_clawbio_timeout_hours_in_rnaseq_pipeline_allowlist():
    """--timeout-hours is a real wrapper flag (and is exposed for scrnaseq-pipeline and
    sarek-pipeline); it must also be forwardable for rnaseq-pipeline rather than being
    silently dropped by the extra-args filter. It takes a value, so it belongs to the
    value allowlist, not the boolean set."""
    values, booleans = _rnaseq_pipeline_allowlist()
    assert "--timeout-hours" in values, f"--timeout-hours missing from allowlist: {sorted(values)!r}"
    assert "--timeout-hours" not in booleans


def test_clawbio_allow_pipeline_version_override_in_rnaseq_pipeline_allowlist():
    """--allow-pipeline-version-override is a store_true wrapper flag (exposed for
    scrnaseq-pipeline and sarek-pipeline); it must be forwardable for rnaseq-pipeline
    and must live in the boolean (value-free) set so the filter does not consume the
    following token as its value."""
    values, booleans = _rnaseq_pipeline_allowlist()
    assert "--allow-pipeline-version-override" in (values | booleans), (
        f"--allow-pipeline-version-override missing from allowlist: {sorted(values | booleans)!r}"
    )
    assert "--allow-pipeline-version-override" in booleans


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


def test_rnaseq_pipeline_allowlist_consistent_with_wrapper():
    """The clawbio runner allowlist must stay coherent with the wrapper parser:

    * every allow-listed flag is a real wrapper flag (else it would error as
      ``unrecognized arguments`` when forwarded), and
    * the value/boolean classification matches the wrapper so the extra-args
      filter neither drops a value nor swallows the next token.
    """
    values, booleans = _rnaseq_pipeline_allowlist()
    allowed = values | booleans
    all_opts, no_value = _wrapper_option_flags(_load_skill_module())

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


def test_main_keyboard_interrupt_returns_130(tmp_path, monkeypatch):
    """Ctrl+C during a long-running pipeline must exit 130 (SIGINT convention),
    not dump a traceback — parity with nfcore-sarek/scrnaseq."""
    module = _load_skill_module()

    def _boom(*_args, **_kwargs):
        raise KeyboardInterrupt

    monkeypatch.setattr(module, "_run_wrapper", _boom)
    rc = module.main(["--output", str(tmp_path), "--demo"])
    assert rc == 130
