"""Audit round 3 — regression tests for the production-readiness findings.

Each test pins one finding from the round-3 audit so the fix cannot silently
regress:

* F-01  Cell Ranger family + conda/mamba profile fails by default (Cell Ranger is
        not distributed via bioconda; only trusted site configs may opt in).
* F-02  FastQC is a hard gate for EVERY aligner unless --skip-fastqc (audit H-02
        corrected the earlier Cell-Ranger-family leniency): the 4.1.0 workflow runs
        FASTQC on the shared ch_fastq consumed by all aligner branches and output.md
        states FastQC is applied to all aligners' input reads.
* F-03  argparse ``choices`` and the version/profile/preset policy files stay in
        lockstep with the 4.1.0 contract (no manual-sync drift).
* F-05  the exact 4.1.0 parameter-name surface is frozen, so a future edit that
        drops/adds a contract key fails loudly.
* F-08  a LOCAL ``--igenomes-base`` mirror is existence-checked in preflight when
        ``--genome`` is used, matching the fail-fast guarantee for every other
        reference path.
"""

from __future__ import annotations

import json
from argparse import Namespace
from pathlib import Path
import sys

import pytest

_SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(_SKILL_DIR))

import pipeline_source
import preflight
from command_builder import build_nextflow_command
from errors import SkillError
from nfcore_4_1_0_contract import (
    CELLRANGER_FAMILY_PRESETS,
    FASTQC_GATED_ALIGNERS,
    OFFICIAL_PARAMS,
    INTENTIONALLY_UNSUPPORTED_PARAMS,
    WRAPPER_SUPPORTED_UPSTREAM_PARAMS,
)
from nfcore_scrnaseq_wrapper import build_parser
from nfcore_scrnaseq_wrapper import _build_parameter_support_summary
from outputs_parser import validate_expected_outputs
from samplesheet_builder import validate_and_normalize_samplesheet
from schemas import (
    PRESET_ALIGNERS,
    SUPPORTED_PRESETS,
    SUPPORTED_PROFILES,
)


# ── F-01: Cell Ranger family + conda/mamba runtime policy ─────────────────────


def _cellranger_conda_args() -> Namespace:
    return Namespace(preset="cellranger", profile="conda", demo=False)


def test_conda_profile_with_cellranger_preset_warns(capsys):
    preflight._warn_if_conda_cellranger(_cellranger_conda_args())
    err = capsys.readouterr().err
    assert "Cell Ranger" in err
    assert "conda" in err.lower()


def test_conda_profile_with_cellranger_preset_is_rejected_by_default():
    with pytest.raises(SkillError) as exc:
        preflight._check_cellranger_runtime_policy(_cellranger_conda_args())
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["preset"] == "cellranger"
    assert exc.value.details["profile"] == "conda"


def test_conda_profile_with_cellranger_preset_can_be_explicitly_allowed(capsys):
    args = _cellranger_conda_args()
    args.allow_conda_cellranger = True
    preflight._check_cellranger_runtime_policy(args)
    assert "Cell Ranger" in capsys.readouterr().err


def test_mamba_profile_with_cellrangermulti_preset_warns(capsys):
    preflight._warn_if_conda_cellranger(
        Namespace(preset="cellrangermulti", profile="mamba,arm64", demo=False)
    )
    assert "Cell Ranger" in capsys.readouterr().err


def test_docker_profile_with_cellranger_preset_does_not_warn(capsys):
    preflight._warn_if_conda_cellranger(
        Namespace(preset="cellranger", profile="docker", demo=False)
    )
    assert capsys.readouterr().err == ""


def test_conda_profile_with_star_preset_does_not_warn(capsys):
    preflight._warn_if_conda_cellranger(
        Namespace(preset="star", profile="conda", demo=False)
    )
    assert capsys.readouterr().err == ""


def test_cellranger_family_constant_is_the_three_cellranger_presets():
    assert CELLRANGER_FAMILY_PRESETS == {
        "cellranger",
        "cellrangerarc",
        "cellrangermulti",
    }


# ── F-02: FastQC is a hard gate only for simpleaf/star/kallisto ───────────────


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


@pytest.mark.parametrize("aligner", sorted(FASTQC_GATED_ALIGNERS))
def test_missing_fastqc_is_required_failure_for_gated_aligners(aligner):
    validation = validate_expected_outputs(
        _outputs_without_fastqc(aligner), aligner=aligner, skip_multiqc=False
    )
    assert "fastqc" in validation["missing_required"]


@pytest.mark.parametrize(
    "aligner", ["cellranger", "cellrangerarc", "cellrangermulti"]
)
def test_missing_fastqc_is_required_failure_for_cellranger_family(aligner):
    validation = validate_expected_outputs(
        _outputs_without_fastqc(aligner), aligner=aligner, skip_multiqc=False
    )
    # nf-core/scrnaseq 4.1.0 runs FASTQC on the shared ch_fastq for EVERY aligner
    # (incl. the Cell Ranger family) and publishes results/fastqc/, so a missing
    # fastqc/ tree is a genuine failure, not optional (audit H-02 corrects F-02).
    assert "fastqc" in validation["missing_required"]
    assert "fastqc" not in validation["missing_optional"]


def test_cellranger_family_still_fails_on_missing_h5ad():
    outputs = _outputs_without_fastqc("cellrangermulti")
    outputs["h5ad_candidates"] = []
    validation = validate_expected_outputs(
        outputs, aligner="cellrangermulti", skip_multiqc=False
    )
    assert any("mtx_conversions" in entry for entry in validation["missing_required"])


def test_fastqc_gated_aligners_are_every_aligner():
    # FastQC runs (and publishes results/fastqc/) for every aligner in 4.1.0, so the
    # hard gate covers all of them — including the Cell Ranger family (audit H-02).
    assert FASTQC_GATED_ALIGNERS == set(PRESET_ALIGNERS.values())


# ── F-03: argparse choices and policy files mirror the contract ───────────────


def _parser_choices(flag: str) -> list[str]:
    parser = build_parser()
    for action in parser._actions:
        if flag in action.option_strings:
            return list(action.choices)
    raise AssertionError(f"flag {flag!r} not found in parser")


@pytest.mark.parametrize(
    "flag,param",
    [
        ("--star-feature", "star_feature"),
        ("--simpleaf-umi-resolution", "simpleaf_umi_resolution"),
        ("--kb-workflow", "kb_workflow"),
        ("--publish-dir-mode", "publish_dir_mode"),
    ],
)
def test_argparse_choices_match_official_enum(flag, param):
    assert set(_parser_choices(flag)) == set(OFFICIAL_PARAMS[param]["enum"])


def test_aligner_enum_matches_preset_aligners():
    assert set(OFFICIAL_PARAMS["aligner"]["enum"]) == set(PRESET_ALIGNERS.values())


def _pinned_versions() -> dict:
    return json.loads(
        (_SKILL_DIR / "reproducibility" / "pinned_versions.json").read_text(
            encoding="utf-8"
        )
    )


def test_pinned_versions_profiles_match_supported_profiles():
    assert set(_pinned_versions()["supported_profiles"]) == SUPPORTED_PROFILES


def test_pinned_versions_presets_match_supported_presets():
    assert set(_pinned_versions()["supported_presets"]) == SUPPORTED_PRESETS


# ── F-05: freeze the exact 4.1.0 parameter-name surface ───────────────────────

# Verbatim parameter names exposed by nf-core/scrnaseq 4.1.0
# (nextflow_schema.json / parameters page). Editing OFFICIAL_PARAMS without
# updating this frozen set — or vice versa — fails this test, so contract drift
# from the real pipeline cannot pass silently (audit F-05).
_NFCORE_4_1_0_PARAM_NAMES = frozenset(
    {
        "input",
        "outdir",
        "email",
        "multiqc_title",
        "barcode_whitelist",
        "aligner",
        "protocol",
        "skip_multiqc",
        "skip_fastqc",
        "skip_cellbender",
        "skip_emptydrops",
        "genome",
        "fasta",
        "igenomes_ignore",
        "transcript_fasta",
        "gtf",
        "save_reference",
        "save_align_intermeds",
        "igenomes_base",
        "txp2gene",
        "simpleaf_index",
        "simpleaf_umi_resolution",
        "star_index",
        "star_ignore_sjdbgtf",
        "seq_center",
        "star_feature",
        "kallisto_index",
        "kb_t1c",
        "kb_t2c",
        "kb_workflow",
        "cellranger_index",
        "skip_cellranger_renaming",
        "motifs",
        "cellrangerarc_config",
        "cellrangerarc_reference",
        "cellranger_vdj_index",
        "skip_cellrangermulti_vdjref",
        "gex_frna_probe_set",
        "gex_target_panel",
        "gex_cmo_set",
        "fb_reference",
        "vdj_inner_enrichment_primers",
        "gex_barcode_sample_assignment",
        "cellranger_multi_barcodes",
        "custom_config_version",
        "custom_config_base",
        "config_profile_name",
        "config_profile_description",
        "config_profile_contact",
        "config_profile_url",
        "version",
        "publish_dir_mode",
        "email_on_fail",
        "plaintext_email",
        "max_multiqc_email_size",
        "monochrome_logs",
        "hook_url",
        "multiqc_config",
        "multiqc_logo",
        "multiqc_methods_description",
        "validate_params",
        "pipelines_testdata_base_path",
        "trace_report_suffix",
        "help",
        "help_full",
        "show_hidden",
    }
)


def test_official_params_keys_match_frozen_4_1_0_surface():
    assert set(OFFICIAL_PARAMS) == _NFCORE_4_1_0_PARAM_NAMES


def test_check_mode_parameter_support_summary_is_complete():
    summary = _build_parameter_support_summary()
    assert set(summary["supported_upstream"]) == WRAPPER_SUPPORTED_UPSTREAM_PARAMS
    assert set(summary["intentionally_unsupported"]) == INTENTIONALLY_UNSUPPORTED_PARAMS
    assert set(summary["official_params"]) == set(OFFICIAL_PARAMS)
    assert summary["nfcore_version"] == "4.1.0"


# ── F-08: local --igenomes-base is existence-checked in preflight ─────────────


def test_local_igenomes_base_missing_dir_is_rejected(tmp_path):
    missing = tmp_path / "no_such_mirror"
    args = Namespace(genome="GRCh38", igenomes_base=str(missing), demo=False)
    with pytest.raises(SkillError) as exc:
        preflight._check_igenomes_base(args)
    assert exc.value.error_code == "MISSING_REFERENCE"
    assert exc.value.details["field"] == "igenomes_base"


def test_local_igenomes_base_existing_dir_is_accepted(tmp_path):
    mirror = tmp_path / "igenomes_mirror"
    mirror.mkdir()
    args = Namespace(genome="GRCh38", igenomes_base=str(mirror), demo=False)
    preflight._check_igenomes_base(args)  # must not raise


def test_remote_igenomes_base_is_not_existence_checked():
    args = Namespace(
        genome="GRCh38", igenomes_base="s3://ngi-igenomes/igenomes/", demo=False
    )
    preflight._check_igenomes_base(args)  # must not raise


# ── F-09/F-10: cloud work-dir and require-local safeguards ────────────────────


def test_nextflow_command_preserves_remote_work_dir_uri(tmp_path):
    params_path = tmp_path / "params.yaml"
    params_path.write_text("aligner: star\n", encoding="utf-8")
    command, command_str = build_nextflow_command(
        pipeline_source={
            "source_kind": "remote_repo",
            "source_ref": "nf-core/scrnaseq",
            "resolved_version": "4.1.0",
        },
        profile="docker,awsbatch",
        params_path=params_path,
        resume=False,
        work_dir="s3://bucket/scrnaseq/work",
    )
    assert command[command.index("-work-dir") + 1] == "s3://bucket/scrnaseq/work"
    assert "s3://bucket/scrnaseq/work" in command_str


def test_require_local_pipeline_rejects_whitespace_fallback(tmp_path):
    local = tmp_path / "my scrnaseq checkout"
    local.mkdir()
    with pytest.raises(SkillError) as exc:
        pipeline_source.resolve_pipeline_source(
            requested_version="4.1.0",
            local_pipeline_dir=local,
            require_local=True,
        )
    assert exc.value.error_code == "PIPELINE_SOURCE_INVALID"
    assert exc.value.details["require_local"] is True


# ── F-11: RNA velocity pairings fail fast ─────────────────────────────────────


def test_star_velocity_requires_star_ignore_sjdbgtf():
    args = Namespace(preset="star", star_feature="Gene Velocyto", star_ignore_sjdbgtf=False)
    with pytest.raises(SkillError) as exc:
        preflight._check_rna_velocity_pairings(args)
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_field"] == "star_ignore_sjdbgtf"


def test_kallisto_velocity_partial_capture_files_is_rejected():
    # Providing exactly one capture file is always an inconsistent/partial spec,
    # regardless of whether an index is prebuilt — reject it.
    args = Namespace(preset="kallisto", kb_workflow="nac", kb_t1c="/refs/t1c.txt", kb_t2c=None)
    with pytest.raises(SkillError) as exc:
        preflight._check_rna_velocity_pairings(args)
    assert exc.value.error_code == "INVALID_PRESET_CONFIGURATION"
    assert exc.value.details["missing_fields"] == ["kb_t2c"]


def test_kallisto_velocity_from_fasta_gtf_does_not_require_capture_files():
    # The documented `--kb-workflow nac --fasta --gtf` example builds the index
    # from scratch; `kb ref` generates the capture files, so neither is required.
    # The 4.1.0 parameter docs never mark kb_t1c/kb_t2c as required (audit F-11
    # regression: enforcing them unconditionally rejected this documented run).
    args = Namespace(
        preset="kallisto", kb_workflow="nac", kb_t1c=None, kb_t2c=None,
        kallisto_index=None,
    )
    preflight._check_rna_velocity_pairings(args)  # must not raise


def test_kallisto_velocity_with_prebuilt_index_requires_capture_files():
    # A prebuilt --kallisto-index is NOT regenerated, so the capture lists must be
    # supplied alongside it for lamanno/nac.
    args = Namespace(
        preset="kallisto", kb_workflow="nac", kb_t1c=None, kb_t2c=None,
        kallisto_index="/refs/kallisto_idx",
    )
    with pytest.raises(SkillError) as exc:
        preflight._check_rna_velocity_pairings(args)
    assert exc.value.details["missing_fields"] == ["kb_t1c", "kb_t2c"]


def test_standard_kallisto_workflow_does_not_require_capture_files():
    args = Namespace(preset="kallisto", kb_workflow="standard", kb_t1c=None, kb_t2c=None)
    preflight._check_rna_velocity_pairings(args)


# ── F-12: docs must not retain stale audit-era limitations ────────────────────


def test_skill_docs_do_not_claim_work_dir_is_always_pinned():
    text = (_SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert "-work-dir <output>/upstream/work" not in text
    assert "not supported even via `-c`" not in text
    assert "does **not** validate this pairing" not in text


def test_igenomes_base_ignored_without_genome(tmp_path):
    # No --genome → iGenomes is not consumed, so a local base path is irrelevant
    # and must not be checked (avoids false positives).
    args = Namespace(
        genome=None, igenomes_base=str(tmp_path / "missing"), demo=False
    )
    preflight._check_igenomes_base(args)  # must not raise


# ── I-2: expected_cells advisory covers the whole Cell Ranger family ───────────


def _touch(path: Path) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_bytes(b"")


def test_expected_cells_advisory_covers_cellrangerarc(tmp_path, capsys):
    # cellranger-arc count also auto-estimates cells (no --expect-cells), so the
    # "since Cell Ranger v7 don't supply expected_cells" advisory applies to it
    # too — the whole CELLRANGER_FAMILY_PRESETS, not just cellranger/multi (I-2).
    r1 = tmp_path / "s_S1_L001_R1_001.fastq.gz"
    r2 = tmp_path / "s_S1_L001_R2_001.fastq.gz"
    _touch(r1)
    _touch(r2)
    src = tmp_path / "sheet.csv"
    src.write_text(
        "sample,fastq_1,fastq_2,sample_type,fastq_barcode,expected_cells\n"
        f"s,{r1},{r2},gex,,1000\n",
        encoding="utf-8",
    )
    validate_and_normalize_samplesheet(
        src, tmp_path / "out.csv", preset="cellrangerarc"
    )
    out = capsys.readouterr().out
    assert "expected_cells" in out
    assert "auto-estimate" in out.lower()


def test_expected_cells_advisory_not_emitted_for_simpleaf(tmp_path, capsys):
    r1 = tmp_path / "a_R1.fastq.gz"
    r2 = tmp_path / "a_R2.fastq.gz"
    _touch(r1)
    _touch(r2)
    src = tmp_path / "sheet.csv"
    src.write_text(
        f"sample,fastq_1,fastq_2,expected_cells\na,{r1},{r2},1000\n",
        encoding="utf-8",
    )
    validate_and_normalize_samplesheet(src, tmp_path / "out.csv", preset="standard")
    out = capsys.readouterr().out
    assert "auto-estimate" not in out.lower()
