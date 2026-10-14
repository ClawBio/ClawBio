"""Audit hardening for nfcore-rnaseq-wrapper — parity with the scrnaseq wrapper.

Three confirmed observations from the cross-skill audit are fixed here, each
pinned by a regression test:

* OBS-2  --pipeline-version is gated to the pinned contract (3.26.0): a different
         remote version is rejected unless --allow-pipeline-version-override
         (the wrapper's parameter/validation logic is 3.26.0-specific).
* OBS-3  --demo is FULLY hermetic: no QC/skip/tuning/contaminant/ribo/umi/
         reporting flag leaks into params.yaml (the upstream test profile owns
         every pipeline parameter; a -params-file value overrides profile config).
* OBS-1  a completed run that produced NO merged count matrix when one was
         expected fails with EXPECTED_OUTPUTS_NOT_FOUND instead of a silent
         success — while hisat2-no-quant and --skip_quantification_merge stay OK.

Module loading mirrors the canonical pattern used by the other wrapper test files
(load the skill under a unique module name and purge its bare submodules from
sys.modules) so this file does not perturb cross-file module-identity isolation.
"""

from __future__ import annotations

import importlib.util
from argparse import Namespace
import json
import os
from pathlib import Path
import sys

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
SCRIPT_PATH = SKILL_DIR / "nfcore_rnaseq_wrapper.py"

_SKILL_BARE_MODULES = (
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
)


def _load_skill_module():
    """Load the wrapper under a unique name, then purge its bare submodules and
    drop SKILL_DIR from sys.path so the next test file imports its own skill's
    modules cleanly (same isolation dance as test_nfcore_rnaseq_wrapper.py)."""
    sys.path.insert(0, str(SKILL_DIR))
    spec = importlib.util.spec_from_file_location(
        "nfcore_rnaseq_wrapper_hardening_module", SCRIPT_PATH
    )
    assert spec is not None and spec.loader is not None
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    for name in _SKILL_BARE_MODULES:
        loaded = sys.modules.get(name)
        loaded_file = Path(getattr(loaded, "__file__", "") or "")
        if loaded is not None and (
            loaded_file == SKILL_DIR / f"{name}.py" or SKILL_DIR in loaded_file.parents
        ):
            sys.modules.pop(name, None)
    if str(SKILL_DIR) in sys.path:
        sys.path.remove(str(SKILL_DIR))
    return module


wrapper = _load_skill_module()
# Bind the symbols from the loaded module so SkillError identity matches what the
# wrapper actually raises (avoids cross-module class-identity mismatches).
SkillError = wrapper.SkillError
build_effective_params = wrapper.build_effective_params
DEFAULT_PIPELINE_VERSION = wrapper.DEFAULT_PIPELINE_VERSION


def _samplesheet(tmp_path: Path) -> Path:
    ss = tmp_path / "reproducibility" / "samplesheet.valid.csv"
    ss.parent.mkdir(parents=True, exist_ok=True)
    ss.write_text("sample,fastq_1,strandedness\n", encoding="utf-8")
    return ss


# ── OBS-2: pipeline-version gate ──────────────────────────────────────────────


def test_pinned_version_passes():
    wrapper._check_pipeline_version_supported(
        Namespace(
            pipeline_version=DEFAULT_PIPELINE_VERSION,
            allow_pipeline_version_override=False,
        )
    )  # must not raise


def test_non_pinned_version_rejected_by_default():
    with pytest.raises(SkillError) as exc:
        wrapper._check_pipeline_version_supported(
            Namespace(pipeline_version="3.20.0", allow_pipeline_version_override=False)
        )
    assert exc.value.error_code == "PIPELINE_SOURCE_INVALID"
    assert exc.value.details["contract_version"] == DEFAULT_PIPELINE_VERSION
    assert exc.value.details["requested"] == "3.20.0"


def test_override_allows_non_pinned(capsys):
    wrapper._check_pipeline_version_supported(
        Namespace(pipeline_version="dev", allow_pipeline_version_override=True)
    )  # must not raise
    assert DEFAULT_PIPELINE_VERSION in capsys.readouterr().err


def test_allow_pipeline_version_override_flag_exists():
    parsed = wrapper.build_parser().parse_args(
        ["--output", "/tmp/x", "--allow-pipeline-version-override"]
    )
    assert parsed.allow_pipeline_version_override is True


# ── OBS-3: fully-hermetic demo ────────────────────────────────────────────────


def test_demo_params_omit_tuning_skip_contaminant_flags(tmp_path):
    ss = _samplesheet(tmp_path)
    args = wrapper.build_parser().parse_args(
        [
            "--output", str(tmp_path), "--demo", "--aligner", "star_salmon",
            "--pseudo-aligner", "kallisto",
            "--trimmer", "fastp",
            "--contaminant-screening", "kraken2",
            "--with-umi",
            "--gencode",
            "--skip-umi-extract",
            "--remove-ribo-rna",
            "--save-bbsplit-reads",
        ]
    )
    params = build_effective_params(args, normalized_samplesheet=ss, output_dir=tmp_path)
    # Only the wrapper-forced essentials remain.
    assert params["aligner"] == "star_salmon"
    assert params.get("igenomes_ignore") is True
    assert "input" not in params
    # Nothing else may leak (a -params-file value overrides the test profile).
    for leaked in (
        "pseudo_aligner", "trimmer", "contaminant_screening", "with_umi",
        "gencode", "skip_umi_extract", "remove_ribo_rna", "save_bbsplit_reads",
    ):
        assert leaked not in params, f"--demo must not write {leaked!r}"


def test_non_demo_still_writes_tuning_flags(tmp_path):
    # Guard: the hermetic-demo change must not affect normal runs.
    ss = _samplesheet(tmp_path)
    fa = tmp_path / "g.fa"
    fa.write_text(">c\nA\n", encoding="utf-8")
    gtf = tmp_path / "g.gtf"
    gtf.write_text("x\n", encoding="utf-8")
    args = wrapper.build_parser().parse_args(
        [
            "--output", str(tmp_path), "--input", str(ss), "--aligner", "star_salmon",
            "--fasta", str(fa), "--gtf", str(gtf),
            "--pseudo-aligner", "kallisto", "--trimmer", "fastp",
        ]
    )
    params = build_effective_params(args, normalized_samplesheet=ss, output_dir=tmp_path)
    assert params.get("pseudo_aligner") == "kallisto"
    assert params.get("trimmer") == "fastp"


# ── OBS-1: required-output gate (no silent empty success) ─────────────────────


def _parsed(**over):
    base = {
        "aligner_effective": "star_salmon",
        "pseudo_aligner_effective": "",
        "skip_quantification_merge": False,
        "pipeline_info_dir": "/out/upstream/results/pipeline_info",
        "preferred_counts_tsv": "",
        "raw_counts_tsv": "",
    }
    base.update(over)
    return base


def _gate_args(**over):
    base = dict(skip_alignment=False, skip_pseudo_alignment=False)
    base.update(over)
    return Namespace(**base)


def test_missing_count_matrix_for_star_salmon_raises(tmp_path):
    with pytest.raises(SkillError) as exc:
        wrapper._raise_if_expected_outputs_missing(
            _parsed(), args=_gate_args(), output_dir=tmp_path
        )
    assert exc.value.error_code == "EXPECTED_OUTPUTS_NOT_FOUND"


def test_present_count_matrix_passes(tmp_path):
    wrapper._raise_if_expected_outputs_missing(
        _parsed(preferred_counts_tsv="/out/star_salmon/salmon.merged.gene_counts.tsv"),
        args=_gate_args(),
        output_dir=tmp_path,
    )  # must not raise


def test_hisat2_without_pseudo_aligner_does_not_raise(tmp_path):
    # hisat2 has no quantifier; a merged count matrix is legitimately absent.
    wrapper._raise_if_expected_outputs_missing(
        _parsed(aligner_effective="hisat2"), args=_gate_args(), output_dir=tmp_path
    )  # must not raise


def test_skip_quantification_merge_does_not_raise(tmp_path):
    # --skip-quantification-merge produces per-sample quant, no merged matrix.
    wrapper._raise_if_expected_outputs_missing(
        _parsed(skip_quantification_merge=True),
        args=_gate_args(),
        output_dir=tmp_path,
    )  # must not raise


def test_missing_pipeline_info_raises(tmp_path):
    with pytest.raises(SkillError) as exc:
        wrapper._raise_if_expected_outputs_missing(
            _parsed(pipeline_info_dir="", preferred_counts_tsv="/x/counts.tsv"),
            args=_gate_args(),
            output_dir=tmp_path,
        )
    assert exc.value.error_code == "EXPECTED_OUTPUTS_NOT_FOUND"


def test_pseudo_only_run_still_requires_counts(tmp_path):
    # --skip-alignment with a pseudo-aligner: counts still expected (from salmon).
    with pytest.raises(SkillError) as exc:
        wrapper._raise_if_expected_outputs_missing(
            _parsed(pseudo_aligner_effective="salmon"),
            args=_gate_args(skip_alignment=True),
            output_dir=tmp_path,
        )
    assert exc.value.error_code == "EXPECTED_OUTPUTS_NOT_FOUND"


# ── OBS-1 integration: real parse_outputs tree → gate (closes unit↔real gap) ──


def _make_results_tree(output_dir: Path, *, with_counts: bool, aligner: str = "star_salmon") -> None:
    results = output_dir / "upstream" / "results"
    (results / "pipeline_info").mkdir(parents=True)
    if with_counts:
        aligner_dir = results / aligner
        aligner_dir.mkdir(parents=True)
        (aligner_dir / "salmon.merged.gene_counts.tsv").write_text(
            "gene_id\tSAMPLE_1\n", encoding="utf-8"
        )


def test_gate_integration_real_tree_with_counts_passes(tmp_path):
    out = tmp_path / "run"
    _make_results_tree(out, with_counts=True)
    args = wrapper.build_parser().parse_args(
        ["--output", str(out), "--input", "x.csv", "--aligner", "star_salmon"]
    )
    parsed = wrapper._parse_outputs_with_effective_aligner(out, args)
    wrapper._raise_if_expected_outputs_missing(parsed, args=args, output_dir=out)  # no raise


def test_gate_integration_real_tree_missing_counts_raises(tmp_path):
    out = tmp_path / "run"
    _make_results_tree(out, with_counts=False)
    args = wrapper.build_parser().parse_args(
        ["--output", str(out), "--input", "x.csv", "--aligner", "star_salmon"]
    )
    parsed = wrapper._parse_outputs_with_effective_aligner(out, args)
    with pytest.raises(SkillError) as exc:
        wrapper._raise_if_expected_outputs_missing(parsed, args=args, output_dir=out)
    assert exc.value.error_code == "EXPECTED_OUTPUTS_NOT_FOUND"


def test_gate_integration_skip_quant_merge_real_tree_passes(tmp_path):
    out = tmp_path / "run"
    _make_results_tree(out, with_counts=False)  # no merged matrix on disk
    args = wrapper.build_parser().parse_args(
        ["--output", str(out), "--input", "x.csv", "--aligner", "star_salmon",
         "--skip-quantification-merge"]
    )
    parsed = wrapper._parse_outputs_with_effective_aligner(out, args)
    wrapper._raise_if_expected_outputs_missing(parsed, args=args, output_dir=out)  # no raise


def test_gate_integration_hisat2_no_quant_real_tree_passes(tmp_path):
    out = tmp_path / "run"
    _make_results_tree(out, with_counts=False, aligner="hisat2")
    args = wrapper.build_parser().parse_args(
        ["--output", str(out), "--input", "x.csv", "--aligner", "hisat2"]
    )
    parsed = wrapper._parse_outputs_with_effective_aligner(out, args)
    wrapper._raise_if_expected_outputs_missing(parsed, args=args, output_dir=out)  # no raise


# ── Audit follow-up F-01: pinned-version contract on a local checkout ─────────


def _local_source(manifest_version):
    return {
        "source_kind": "local_checkout",
        "source_ref": "/repo/rnaseq",
        "resolved_version": "abc123",
        "manifest_version": manifest_version,
        "branch": "main",
        "dirty": False,
    }


def test_local_checkout_matching_version_passes():
    wrapper._check_pipeline_source_version(
        Namespace(allow_pipeline_version_override=False),
        _local_source(DEFAULT_PIPELINE_VERSION),
    )  # must not raise


def test_local_checkout_mismatched_version_rejected():
    with pytest.raises(SkillError) as exc:
        wrapper._check_pipeline_source_version(
            Namespace(allow_pipeline_version_override=False),
            _local_source("3.20.0"),
        )
    assert exc.value.error_code == "PIPELINE_SOURCE_INVALID"
    assert exc.value.details["manifest_version"] == "3.20.0"
    assert exc.value.details["contract_version"] == DEFAULT_PIPELINE_VERSION


def test_local_checkout_mismatch_allowed_with_override(capsys):
    wrapper._check_pipeline_source_version(
        Namespace(allow_pipeline_version_override=True),
        _local_source("3.20.0"),
    )  # must not raise
    assert "3.20.0" in capsys.readouterr().err


def test_local_checkout_unknown_version_warns_not_blocks(capsys):
    wrapper._check_pipeline_source_version(
        Namespace(allow_pipeline_version_override=False),
        _local_source(""),
    )  # must not raise
    assert DEFAULT_PIPELINE_VERSION in capsys.readouterr().err


def test_remote_source_is_not_subject_to_local_version_gate():
    wrapper._check_pipeline_source_version(
        Namespace(allow_pipeline_version_override=False),
        {"source_kind": "remote_repo", "manifest_version": "anything"},
    )  # must not raise


# ── Audit follow-up F-02: --timeout-hours wiring ──────────────────────────────


def test_timeout_hours_flag_default():
    parsed = wrapper.build_parser().parse_args(["--output", "/tmp/x"])
    assert parsed.timeout_hours == wrapper.DEFAULT_TIMEOUT_HOURS


def test_resolve_timeout_seconds_uses_arg():
    assert wrapper._resolve_timeout_seconds(Namespace(timeout_hours=48)) == 48 * 3600


def test_resolve_timeout_seconds_zero_disables_cap():
    """--timeout-hours 0 disables the wall-clock cap (returns None) — parity with
    nfcore-scrnaseq/sarek."""
    assert wrapper._resolve_timeout_seconds(Namespace(timeout_hours=0)) is None


def test_resolve_work_dir_default(tmp_path):
    assert wrapper._resolve_nextflow_work_dir(Namespace(work_dir=None), tmp_path) == (
        tmp_path / "upstream" / "work"
    )


def test_resolve_work_dir_object_store_uri_preserved(tmp_path):
    """Remote object-store work dirs are preserved verbatim (cloud executors) —
    parity with nfcore-scrnaseq/sarek."""
    assert (
        wrapper._resolve_nextflow_work_dir(Namespace(work_dir="s3://bucket/work"), tmp_path)
        == "s3://bucket/work"
    )


def test_resolve_timeout_seconds_falls_back_to_default():
    assert (
        wrapper._resolve_timeout_seconds(Namespace(timeout_hours=None))
        == wrapper.DEFAULT_TIMEOUT_HOURS * 3600
    )


# ── Audit follow-up F-07: macOS Docker memory ceiling derives from host RAM ────


def test_macos_docker_memory_never_exceeds_ceiling(monkeypatch):
    # Isolate the host-RAM ceiling logic: force the Docker-VM probe to "unknown" so the
    # test is deterministic regardless of whether a real Docker daemon is running on the
    # host (the VM cap is exercised separately). Without this, a live Docker daemon makes
    # the (correct) 90%-of-VM cap leak in and the test becomes environment-dependent.
    monkeypatch.setattr(wrapper, "_docker_vm_memory_gb", lambda: None)
    monkeypatch.setattr(wrapper, "_detect_host_memory_gb", lambda: 256)
    assert wrapper._macos_docker_memory_gb() == wrapper._MACOS_DOCKER_MEMORY_CEILING_GB


def test_macos_docker_memory_lowers_on_small_host(monkeypatch):
    monkeypatch.setattr(wrapper, "_docker_vm_memory_gb", lambda: None)
    monkeypatch.setattr(wrapper, "_detect_host_memory_gb", lambda: 8)
    mem = wrapper._macos_docker_memory_gb()
    assert mem <= 8
    assert mem >= wrapper._MACOS_DOCKER_MEMORY_FLOOR_GB


def test_macos_docker_memory_falls_back_when_undetectable(monkeypatch):
    monkeypatch.setattr(wrapper, "_docker_vm_memory_gb", lambda: None)
    monkeypatch.setattr(wrapper, "_detect_host_memory_gb", lambda: None)
    assert wrapper._macos_docker_memory_gb() == wrapper._MACOS_DOCKER_MEMORY_CEILING_GB


def test_macos_docker_memory_capped_to_90pct_of_live_vm(monkeypatch):
    """When the Docker VM's memory IS known, the per-process ceiling is capped to 90% of
    it (audit F-8) — even on a huge host — so a container is never OOM-killed by requesting
    more than the VM has. Mirrors the real-environment behaviour observed with colima."""
    monkeypatch.setattr(wrapper, "_detect_host_memory_gb", lambda: 256)
    monkeypatch.setattr(wrapper, "_docker_vm_memory_gb", lambda: 12)
    assert wrapper._macos_docker_memory_gb() == int(12 * 0.9)


# ── Audit follow-up F-02: demo warns about silently-ignored tuning flags ──────


def test_demo_warns_when_tuning_flags_set(capsys):
    args = wrapper.build_parser().parse_args(
        ["--output", "/tmp/x", "--demo", "--with-umi", "--contaminant-screening", "kraken2"]
    )
    wrapper._apply_demo_overrides(args)
    err = capsys.readouterr().err
    assert "--with-umi" in err
    assert "--contaminant-screening" in err


def test_demo_without_tuning_flags_emits_no_tuning_warning(capsys):
    args = wrapper.build_parser().parse_args(["--output", "/tmp/x", "--demo"])
    wrapper._apply_demo_overrides(args)
    err = capsys.readouterr().err
    assert "ignores tuning flags" not in err


def test_non_demo_run_never_warns_about_tuning(capsys):
    args = wrapper.build_parser().parse_args(
        ["--output", "/tmp/x", "--input", "s.csv", "--with-umi"]
    )
    wrapper._apply_demo_overrides(args)
    assert capsys.readouterr().err == ""


# ── Audit follow-up: CLI enum choices use the schemas single source (no drift) ─


def _parser_choices(dest):
    parser = wrapper.build_parser()
    for action in parser._actions:
        if action.dest == dest:
            return action.choices
    raise AssertionError(f"no argparse action with dest {dest!r}")


def test_salmon_libtype_choices_use_schema_constant():
    assert set(_parser_choices("salmon_quant_libtype")) == wrapper.SUPPORTED_SALMON_LIB_TYPES


def test_umi_extract_method_choices_use_schema_constant():
    assert set(_parser_choices("umitools_extract_method")) == wrapper.SUPPORTED_UMI_EXTRACT_METHODS


def test_umi_grouping_method_choices_use_schema_constant():
    assert set(_parser_choices("umitools_grouping_method")) == wrapper.SUPPORTED_UMI_GROUPING_METHODS


def test_publish_dir_mode_choices_use_schema_constant():
    assert set(_parser_choices("publish_dir_mode")) == wrapper.SUPPORTED_PUBLISH_DIR_MODES


def test_bracken_precision_choices_preserve_schema_order():
    # SUPPORTED_BRACKEN_PRECISION is an ordered tuple (domain→species); --help must
    # preserve that order, so the parser choices must equal it exactly (not a set).
    assert tuple(_parser_choices("bracken_precision")) == wrapper.SUPPORTED_BRACKEN_PRECISION


def test_skipped_samples_hint_added_when_log_signals_it(tmp_path):
    """When nf-core reports it 'completed with skipped sample(s)' and no count
    matrix was produced, the fix must point at --min-trimmed-reads (nf-core drops
    samples below its default 10000-trimmed-read threshold)."""
    logs = tmp_path / "logs"
    logs.mkdir()
    (logs / "stdout.txt").write_text(
        "-[nf-core/rnaseq] Pipeline completed successfully with skipped sample(s)-\n",
        encoding="utf-8",
    )
    with pytest.raises(SkillError) as exc:
        wrapper._raise_if_expected_outputs_missing(
            _parsed(), args=_gate_args(), output_dir=tmp_path
        )
    assert exc.value.error_code == "EXPECTED_OUTPUTS_NOT_FOUND"
    assert "min-trimmed-reads" in exc.value.fix


def test_no_skipped_samples_hint_without_signal(tmp_path):
    """No skipped-samples signal in the logs -> no min-trimmed-reads hint."""
    with pytest.raises(SkillError) as exc:
        wrapper._raise_if_expected_outputs_missing(
            _parsed(), args=_gate_args(), output_dir=tmp_path
        )
    assert "min-trimmed-reads" not in exc.value.fix


# ── Hermetic test profiles: ignored flags are named, what ran is reported ─────

_REAL_RUN = ["--output", "/tmp/x", "--input", "s.csv"]
_SAMPLE_VALUES = (
    "3", "0.5", "31", "x", "reference.fasta", "annotation.gtf",
    "annotation.gff", "genes.bed", "user@example.org", "3.GB",
)


def _real_params(parser, extra):
    args = parser.parse_args([*_REAL_RUN, *extra])
    return build_effective_params(
        args,
        normalized_samplesheet=Path("/tmp/x/reproducibility/samplesheet.valid.csv"),
        output_dir=Path("/tmp/x"),
    )


def _sample_argv(parser, action):
    """A value accepted by both the parser and the pinned scalar contract."""
    flag = next(option for option in action.option_strings if option.startswith("--"))
    if action.nargs == 0:
        return [flag]
    if action.choices:
        values = [choice for choice in action.choices if choice != action.default][:1]
    else:
        values = list(_SAMPLE_VALUES)
    for value in values:
        try:
            parser.parse_args([*_REAL_RUN, flag, str(value)])
        except SystemExit:
            continue
        try:
            _real_params(parser, [flag, str(value)])
        except SkillError:
            continue
        return [flag, str(value)]
    assert not values, f"No schema-valid sample for {flag} among {values}"
    return None


def test_demo_names_every_pipeline_flag_it_ignores(capsys):
    """--demo runs nf-core's test profile, which owns every pipeline parameter. Each
    flag that would change a real run's params.yaml has no effect there and must be
    named in a warning (--skip-qc used to vanish without one); flags that compose a
    Nextflow profile (--arm, --prokaryotic) still apply and must not be named."""
    parser = wrapper.build_parser()
    baseline = _real_params(parser, [])
    named: list[str] = []
    for action in parser._actions:
        if not action.option_strings or action.dest in {"help", "output", "input", "demo"}:
            continue
        argv = _sample_argv(parser, action)
        if argv is None:
            continue
        changes_params = _real_params(parser, argv) != baseline
        wrapper._apply_demo_overrides(parser.parse_args(["--output", "/tmp/x", "--demo", *argv]))
        err = capsys.readouterr().err
        if changes_params and action.dest not in {"arm", "prokaryotic"}:
            assert argv[0] in err, f"{argv[0]} is ignored under --demo without a warning"
            named.append(argv[0])
        else:
            tuning = [line for line in err.splitlines() if "ignores tuning flags" in line]
            assert not any(argv[0] + ")" in line or argv[0] + "," in line for line in tuning), argv[0]
    assert "--skip-qc" in named and "--pseudo-aligner" in named


def _hermetic_results(output_dir: Path, *launches: dict) -> None:
    """A demo results tree: STAR+Salmon and Salmon pseudo-alignment counts, plus one
    nf-core params record per launch, oldest first."""
    results = output_dir / "upstream" / "results"
    info = results / "pipeline_info"
    info.mkdir(parents=True)
    for second, recorded in enumerate(launches):
        record = info / f"params_2026-10-14_15-57-{second:02d}.json"
        record.write_text(json.dumps(recorded), encoding="utf-8")
        os.utime(record, (1_700_000_000 + second, 1_700_000_000 + second))
    for quant in ("star_salmon", "salmon"):
        (results / quant).mkdir()
        (results / quant / "salmon.merged.gene_counts.tsv").write_text(
            "gene_id\tgene_name\tS1\nG1\tG1\t5\n", encoding="utf-8"
        )


def test_demo_reports_the_pseudo_aligner_nf_core_recorded(tmp_path):
    """nf-core/rnaseq 3.26.0's test profile also runs Salmon pseudo-alignment. The demo
    report said 'Pseudo-aligner: none' and left those counts out; it must state what
    ran, even when the user asked for a pseudo-aligner the profile ignored."""
    _hermetic_results(tmp_path, {"aligner": "star_salmon", "pseudo_aligner": "salmon"})
    args = wrapper.build_parser().parse_args(
        ["--output", str(tmp_path), "--demo", "--pseudo-aligner", "kallisto"]
    )
    wrapper._apply_demo_overrides(args)
    parsed = wrapper._parse_outputs_with_effective_aligner(tmp_path, args)
    results = tmp_path / "upstream" / "results"
    assert parsed["pseudo_aligner_effective"] == "salmon"
    assert parsed["salmon_pseudo_counts_tsv"] == str(results / "salmon" / "salmon.merged.gene_counts.tsv")
    assert parsed["raw_counts_tsv"] == str(results / "star_salmon" / "salmon.merged.gene_counts.tsv")
    report = wrapper.write_report(
        tmp_path, args=args, pipeline_source={"source_kind": "remote_repo", "resolved_version": "3.26.0"},
        preflight_result={"warnings": []}, parsed_outputs=parsed, command_str="nextflow run",
    ).read_text(encoding="utf-8")
    assert "- Pseudo-aligner: `salmon` (set by the test profile)" in report


def test_demo_reads_the_latest_launch_record(tmp_path):
    _hermetic_results(tmp_path, {"aligner": "star_salmon", "pseudo_aligner": None},
                      {"aligner": "star_salmon", "pseudo_aligner": "salmon"})
    args = wrapper.build_parser().parse_args(["--output", str(tmp_path), "--demo"])
    wrapper._apply_demo_overrides(args)
    assert wrapper._parse_outputs_with_effective_aligner(tmp_path, args)["pseudo_aligner_effective"] == "salmon"


def test_real_run_reports_the_requested_pseudo_aligner(tmp_path):
    """Outside a test profile the wrapper writes pseudo_aligner itself, so the request
    is what ran and the report keeps saying so."""
    _hermetic_results(tmp_path, {"aligner": "star_salmon", "pseudo_aligner": "salmon"})
    args = wrapper.build_parser().parse_args(
        ["--output", str(tmp_path), "--input", "s.csv", "--pseudo-aligner", "salmon"]
    )
    parsed = wrapper._parse_outputs_with_effective_aligner(tmp_path, args)
    assert parsed["pseudo_aligner_effective"] == "salmon"
    assert not parsed.get("pseudo_aligner_from_pipeline")
