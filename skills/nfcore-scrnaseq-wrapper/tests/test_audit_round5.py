"""Audit round 5 — regression tests for the second external-audit pass.

Each test pins one corrected finding so it cannot silently regress:

* H-01  ``--demo`` is FULLY hermetic: not only references/protocol but every
        pipeline-tuning/skip/save/reporting flag is kept out of params.yaml (they
        would override the upstream ``test`` profile, since ``-params-file`` beats
        profile config), and the user is warned about ALL ignored flags. Output
        validation in demo also treats FastQC/MultiQC as not-skipped.
* H-02  FastQC is a hard gate for EVERY aligner unless ``--skip-fastqc``: the 4.1.0
        workflow runs FASTQC on the shared ``ch_fastq`` consumed by all aligner
        branches (incl. cellrangermulti) and output.md states "FastQC is applied to
        all aligners' input reads".
* H-03  ``feature_type=crispr`` requires ``--fb-reference``: CRISPR Guide Capture
        shares the same feature reference as Antibody Capture in 4.1.0 (no separate
        CRISPR reference param exists).
* H-04  The ``-c`` params-override lint also catches bracket-notation, whole-map
        assignment, and a ``params`` block whose ``{`` is on the next line.
* H-05  A wall-clock cap left active on an object-store work-dir or an institutional
        (HPC/site) profile emits a hint to pass ``--timeout-hours 0``.
"""

from __future__ import annotations

from argparse import Namespace
from pathlib import Path
import re
import sys

import pytest

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR.parent.parent))

import nfcore_scrnaseq_wrapper as wrapper
import preflight
import samplesheet_builder
from errors import SkillError
from nfcore_4_1_0_contract import (
    FASTQC_GATED_ALIGNERS,
)
from outputs_parser import validate_expected_outputs
from params_builder import build_effective_params
from schemas import PRESET_ALIGNERS


# ── H-01: fully-hermetic demo ─────────────────────────────────────────────────


def _full_args(**overrides) -> Namespace:
    defaults = dict(
        demo=False,
        preset="star",
        protocol=None,
        email=None,
        email_on_fail=None,
        multiqc_title=None,
        multiqc_config=None,
        multiqc_logo=None,
        multiqc_methods_description=None,
        publish_dir_mode=None,
        trace_report_suffix=None,
        monochrome_logs=False,
        skip_cellbender=False,
        skip_fastqc=False,
        skip_emptydrops=False,
        skip_multiqc=False,
        skip_cellranger_renaming=False,
        skip_cellrangermulti_vdjref=False,
        fasta=None,
        gtf=None,
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        star_feature=None,
        star_ignore_sjdbgtf=False,
        seq_center=None,
        simpleaf_umi_resolution=None,
        kb_workflow=None,
        kb_t1c=None,
        kb_t2c=None,
        motifs=None,
        cellrangerarc_config=None,
        cellrangerarc_reference=None,
        cellranger_vdj_index=None,
        gex_frna_probe_set=None,
        gex_target_panel=None,
        gex_cmo_set=None,
        fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
        genome=None,
        igenomes_base=None,
        igenomes_ignore=False,
        save_reference=False,
        save_align_intermeds=None,
    )
    defaults.update(overrides)
    return Namespace(**defaults)


def test_demo_params_are_hermetic_no_tuning_skip_or_reporting_leakage(tmp_path):
    ss = tmp_path / "reproducibility" / "samplesheet.demo.csv"
    ss.parent.mkdir(parents=True)
    ss.write_text("sample,fastq_1,fastq_2\n", encoding="utf-8")
    # demo flow forces star + skip_cellbender before params are built.
    args = _full_args(
        demo=True,
        preset="star",
        skip_cellbender=True,
        star_feature="GeneFull",
        star_ignore_sjdbgtf=True,
        seq_center="MYCENTER",
        skip_fastqc=True,
        skip_multiqc=True,
        save_reference=True,
        save_align_intermeds=False,
        email="a@b.co",
        multiqc_title="X",
        publish_dir_mode="copy",
        monochrome_logs=True,
    )
    params = build_effective_params(
        args, normalized_samplesheet=ss, output_dir=tmp_path
    )
    # Wrapper-forced essentials remain.
    assert params["aligner"] == "star"
    assert params.get("igenomes_ignore") is True
    assert params.get("skip_cellbender") is True
    assert "input" not in params
    # Nothing else may leak (it would override the test profile).
    for leaked in (
        "star_feature",
        "star_ignore_sjdbgtf",
        "seq_center",
        "skip_fastqc",
        "skip_multiqc",
        "save_reference",
        "save_align_intermeds",
        "email",
        "multiqc_title",
        "publish_dir_mode",
        "monochrome_logs",
    ):
        assert leaked not in params, f"--demo must not write {leaked!r}"


def test_demo_ignored_flags_helper_lists_non_reference_flags():
    args = _full_args(
        demo=True,
        preset="star",
        genome="GRCh38",
        protocol="10XV3",
        star_feature="GeneFull",
        seq_center="C",
        skip_fastqc=True,
        save_reference=True,
    )
    ignored = wrapper._demo_ignored_flags(args)
    for flag in (
        "--genome",
        "--protocol",
        "--star-feature",
        "--seq-center",
        "--skip-fastqc",
        "--save-reference",
    ):
        assert flag in ignored, flag
    # skip_cellbender is forced by demo, not ignored — never listed.
    assert "--skip-cellbender" not in ignored
    # Nothing supplied → empty list.
    assert wrapper._demo_ignored_flags(_full_args(demo=True, preset="star")) == []


def test_demo_names_every_pipeline_flag_it_ignores():
    """Each flag that would change a real run's params.yaml is owned by the test
    profile under --demo and must be named, so a new flag cannot be dropped silently.
    The exceptions are the settings --demo enforces itself (CellBender off, iGenomes
    ignored), which the user gets regardless."""
    parser = wrapper.build_parser()
    real = ["--output", "/tmp/x", "--input", "s.csv"]

    def real_params(extra):
        return build_effective_params(
            parser.parse_args([*real, *extra]),
            normalized_samplesheet=Path("/tmp/x/reproducibility/samplesheet.valid.csv"),
            output_dir=Path("/tmp/x"),
        )

    baseline = real_params([])
    enforced_by_demo = {"--preset", "--skip-cellbender", "--skip-emptydrops", "--igenomes-ignore"}
    named = []
    for action in parser._actions:
        if not action.option_strings or action.dest in {"help", "output", "input", "demo"}:
            continue
        flag = next(option for option in action.option_strings if option.startswith("--"))
        candidates = [[flag]] if action.nargs == 0 else [
            [flag, str(value)]
            for value in ([c for c in action.choices if c != action.default][:1] if action.choices else [
                "3", "0.5", "x", "reference.fasta", "user@example.org", "3.GB",
            ])
        ]
        for argv in candidates:
            try:
                parser.parse_args([*real, *argv])
            except SystemExit:
                continue
            try:
                params = real_params(argv)
            except SkillError:
                continue
            if params != baseline and flag not in enforced_by_demo:
                ignored = wrapper._demo_ignored_flags(parser.parse_args(["--output", "/tmp/x", "--demo", *argv]))
                assert flag in ignored, f"{flag} is ignored under --demo without a warning"
                named.append(flag)
            break
        else:
            assert not candidates, f"No schema-valid sample for {flag} among {candidates}"
    assert {"--protocol", "--skip-fastqc"} <= set(named)


def test_demo_validation_requires_fastqc_and_multiqc_even_if_user_passed_skip(tmp_path):
    # Build a results tree WITHOUT multiqc/fastqc; the hermetic demo must still
    # require them because the skip flags were not written to params.yaml.
    results = tmp_path / "upstream" / "results"
    (results / "pipeline_info").mkdir(parents=True)
    (results / "pipeline_info" / "software_versions.yml").write_text(
        "x", encoding="utf-8"
    )
    (results / "pipeline_info" / "params_x.json").write_text("{}", encoding="utf-8")
    mtx = results / "star" / "mtx_conversions"
    mtx.mkdir(parents=True)
    (mtx / "combined_filtered_matrix.h5ad").write_text("x", encoding="utf-8")
    args = _full_args(demo=True, preset="star", skip_multiqc=True, skip_fastqc=True)

    parsed = wrapper._parse_outputs_with_effective_aligner(tmp_path, args)
    missing = parsed["output_validation"]["missing_required"]

    assert "multiqc/multiqc_report.html" in missing
    assert "fastqc" in missing


# ── H-02: FastQC is a hard gate for every aligner ─────────────────────────────


def _outputs_without_fastqc(aligner: str) -> dict:
    base = f"/out/upstream/results/{aligner}"
    return {
        "pipeline_info_dir": "/out/upstream/results/pipeline_info",
        "multiqc_report": "/out/upstream/results/multiqc/multiqc_report.html",
        "h5ad_candidates": [f"{base}/mtx_conversions/combined_filtered_matrix.h5ad"],
        "upstream_dir": "/out/upstream/results",
        "top_level_entries": ["pipeline_info", aligner, "multiqc"],
        "official_outputs": {"fastqc": {"present": False}},
    }


@pytest.mark.parametrize("aligner", ["cellranger", "cellrangerarc", "cellrangermulti"])
def test_missing_fastqc_is_required_failure_for_cellranger_family(aligner):
    validation = validate_expected_outputs(
        _outputs_without_fastqc(aligner), aligner=aligner, skip_multiqc=False
    )
    assert "fastqc" in validation["missing_required"]
    assert "fastqc" not in validation["missing_optional"]


def test_fastqc_gated_aligners_are_all_six_aligners():
    assert FASTQC_GATED_ALIGNERS == set(PRESET_ALIGNERS.values())


def test_skip_fastqc_still_silences_the_cellranger_family_gate():
    validation = validate_expected_outputs(
        _outputs_without_fastqc("cellrangermulti"),
        aligner="cellrangermulti",
        skip_multiqc=False,
        skip_fastqc=True,
    )
    assert "fastqc" not in validation["missing_required"]


# ── H-03: crispr requires a feature-barcode reference ─────────────────────────


def test_cellrangermulti_crispr_requires_feature_barcode_reference():
    args = Namespace(preset="cellrangermulti", fb_reference=None)
    with pytest.raises(SkillError) as exc:
        preflight._check_cellrangermulti_feature_references(args, {"crispr"})
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "fb_reference"


def test_cellrangermulti_crispr_with_fb_reference_is_accepted():
    args = Namespace(preset="cellrangermulti", fb_reference="/refs/features.csv")
    preflight._check_cellrangermulti_feature_references(args, {"crispr"})  # no raise


def test_cellrangermulti_antibody_still_requires_fb_reference():
    args = Namespace(preset="cellrangermulti", fb_reference=None)
    with pytest.raises(SkillError) as exc:
        preflight._check_cellrangermulti_feature_references(args, {"ab"})
    assert exc.value.details["missing_field"] == "fb_reference"


def test_cellrangermulti_gex_only_needs_no_fb_reference():
    args = Namespace(preset="cellrangermulti", fb_reference=None)
    preflight._check_cellrangermulti_feature_references(args, {"gex"})  # no raise


# ── H-04: stricter -c params-override lint ─────────────────────────────────────


def _cfg(tmp_path: Path, body: str) -> Namespace:
    p = tmp_path / "site.config"
    p.write_text(body, encoding="utf-8")
    return Namespace(extra_config=[str(p)], trust_config_params=False)


def test_config_bracket_notation_params_override_is_blocked(tmp_path):
    args = _cfg(tmp_path, "params['aligner'] = 'star'\n")
    with pytest.raises(SkillError):
        preflight._check_config_param_overrides(args)


def test_config_double_quoted_bracket_params_override_is_blocked(tmp_path):
    args = _cfg(tmp_path, 'params["outdir"] = "/elsewhere"\n')
    with pytest.raises(SkillError):
        preflight._check_config_param_overrides(args)


def test_config_whole_params_map_assignment_is_blocked(tmp_path):
    args = _cfg(tmp_path, "params = [aligner: 'star']\n")
    with pytest.raises(SkillError):
        preflight._check_config_param_overrides(args)


def test_config_params_block_with_newline_brace_is_blocked(tmp_path):
    args = _cfg(tmp_path, "params\n{\n  outdir = '/elsewhere'\n}\n")
    with pytest.raises(SkillError):
        preflight._check_config_param_overrides(args)


def test_config_reading_bracket_params_on_rhs_is_not_blocked(tmp_path):
    args = _cfg(tmp_path, "process.cpus = params['max_cpus']\n")
    preflight._check_config_param_overrides(args)  # must not raise


# ── H-05: cap-active-on-HPC/cloud warning ─────────────────────────────────────


def test_warns_when_cap_active_on_object_store_work_dir(capsys):
    args = Namespace(timeout_hours=12.0, work_dir="s3://bucket/work", profile="docker")
    preflight._warn_if_capped_remote_run(args)
    assert "timeout-hours 0" in capsys.readouterr().err


def test_warns_when_cap_active_on_institutional_profile(capsys):
    args = Namespace(timeout_hours=12.0, work_dir=None, profile="crick")
    preflight._warn_if_capped_remote_run(args)
    assert "timeout-hours 0" in capsys.readouterr().err


def test_no_timeout_warning_when_cap_disabled(capsys):
    args = Namespace(timeout_hours=0, work_dir="s3://bucket/work", profile="docker")
    preflight._warn_if_capped_remote_run(args)
    assert capsys.readouterr().err == ""


def test_no_timeout_warning_for_local_docker_run(capsys):
    args = Namespace(timeout_hours=12.0, work_dir=None, profile="docker")
    preflight._warn_if_capped_remote_run(args)
    assert capsys.readouterr().err == ""


# ── H-07: doc regression — protocol-required CLI examples must pass --protocol ──


def _skill_md_wrapper_invocations() -> list[str]:
    """Each SKILL.md shell example that invokes the wrapper/runner, line-joined."""
    text = (_SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    joined = re.sub(r"\\\n\s*", " ", text)  # fold shell line-continuations
    return [
        line.strip()
        for line in joined.splitlines()
        if "nfcore_scrnaseq_wrapper.py" in line or "run scrnaseq-pipeline" in line
    ]


def test_skill_md_cli_examples_pass_protocol_for_protocol_required_presets():
    offenders = []
    for command in _skill_md_wrapper_invocations():
        if "--demo" in command:
            continue  # demo is hermetic — protocol is owned by the test profile
        match = re.search(r"--preset\s+(\S+)", command)
        if not match:
            continue
        if (
            match.group(1) in {"standard", "star", "kallisto"}
            and "--protocol" not in command
        ):
            offenders.append(command)
    assert offenders == [], (
        "SKILL.md CLI examples for standard/star/kallisto must include --protocol "
        f"(nf-core/scrnaseq 4.1.0 refuses its default `auto` for them): {offenders}"
    )


# ── H-08 superseded: sample_type / feature_type values are nf-schema's ─────────
# assets/schema_input.json declares both enums, so nf-schema names an invalid value
# when Nextflow starts; the builder only lower-cases the value. A cellrangermulti
# row still needs a feature_type (align_cellrangermulti drops rows without one).


@pytest.mark.parametrize("preset", ["star", "cellrangerarc", "cellrangermulti", None])
def test_sample_and_feature_type_values_are_left_to_nf_schema(preset):
    assert samplesheet_builder._sample_type({"sample_type": "BEAM"}) == "beam"
    assert samplesheet_builder._feature_type({"feature_type": "Beam"}, 2, preset=preset) == "beam"


def test_empty_types_pass_except_feature_type_under_cellrangermulti():
    assert samplesheet_builder._sample_type({"sample_type": ""}) == ""
    assert samplesheet_builder._feature_type({"feature_type": ""}, 2, preset="star") == ""
    with pytest.raises(SkillError) as exc:
        samplesheet_builder._feature_type({"feature_type": ""}, 2, preset="cellrangermulti")
    assert exc.value.details["column"] == "feature_type"
