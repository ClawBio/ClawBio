from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from reporting import write_check_result, write_report, write_repro_commands, write_result


def _make_args(tmp_path: Path, **kwargs) -> argparse.Namespace:
    defaults = dict(
        input=str(tmp_path / "samplesheet.csv"),
        output=str(tmp_path / "out"),
        preset="star",
        profile="docker",
        pipeline_version="3.14.0",
        protocol=None,
        skip_cellbender=False,
        skip_fastqc=False,
        skip_emptydrops=False,
        skip_multiqc=False,
        resume=False,
        demo=False,
        check=False,
        fasta=None,
        gtf=None,
        transcript_fasta=None,
        txp2gene=None,
        simpleaf_index=None,
        kallisto_index=None,
        star_index=None,
        cellranger_index=None,
        barcode_whitelist=None,
        expected_cells=None,
        star_feature=None,
        star_ignore_sjdbgtf=False,
        seq_center=None,
        simpleaf_umi_resolution=None,
        kb_workflow=None,
        kb_t1c=None,
        kb_t2c=None,
        skip_cellranger_renaming=False,
        motifs=None,
        cellrangerarc_config=None,
        cellrangerarc_reference=None,
        cellranger_vdj_index=None,
        skip_cellrangermulti_vdjref=False,
        gex_frna_probe_set=None,
        gex_target_panel=None,
        gex_cmo_set=None,
        fb_reference=None,
        vdj_inner_enrichment_primers=None,
        gex_barcode_sample_assignment=None,
        cellranger_multi_barcodes=None,
        genome=None,
        save_reference=False,
        save_align_intermeds=False,
        run_downstream=False,
        email=None,
        multiqc_title=None,
    )
    defaults.update(kwargs)
    return argparse.Namespace(**defaults)


def _pipeline_source() -> dict:
    return {
        "source_kind": "remote_repo",
        "source_ref": "nf-core/scrnaseq",
        "resolved_version": "3.14.0",
        "branch": "",
        "dirty": False,
    }


def _preflight_result() -> dict:
    return {
        "ok": True,
        "java": {"version": "21.0.0", "path": "/usr/bin/java"},
        "nextflow": {"version": "25.4.0", "path": "/usr/bin/nextflow"},
        "profile": {"profile": "docker", "backend_path": "/usr/bin/docker", "backend_ready": True},
        "pipeline_source": _pipeline_source(),
        "references": {},
        "samplesheet": {"sample_count": 2, "unknown_columns": []},
    }


def _parsed_outputs(preferred_h5ad: str = "") -> dict:
    return {
        "preferred_h5ad": preferred_h5ad,
        "multiqc_report": "",
        "pipeline_info_dir": "",
        "aligner_effective": "star",
        "cellbender_used": False,
        "handoff_available": bool(preferred_h5ad),
        "samples_detected": ["sampleA", "sampleB"],
        "h5ad_candidates": [],
        "rds_candidates": [],
    }


# ── write_report ──────────────────────────────────────────────────────────────

def test_write_report_creates_report_md(tmp_path):
    args = _make_args(tmp_path)
    report_path = write_report(
        tmp_path,
        args=args,
        pipeline_source=_pipeline_source(),
        preflight_result=_preflight_result(),
        parsed_outputs=_parsed_outputs(),
        command_str="nextflow run nf-core/scrnaseq",
    )
    assert report_path == tmp_path / "report.md"
    assert report_path.exists()


def test_write_report_contains_preset_and_profile(tmp_path):
    args = _make_args(tmp_path, preset="cellranger", profile="singularity")
    write_report(
        tmp_path,
        args=args,
        pipeline_source=_pipeline_source(),
        preflight_result=_preflight_result(),
        parsed_outputs=_parsed_outputs(),
        command_str="nextflow run nf-core/scrnaseq",
    )
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "cellranger" in content
    assert "singularity" in content


def test_write_report_handoff_with_preferred_h5ad(tmp_path):
    preferred = "/output/outs/cellbender.h5ad"
    write_report(
        tmp_path,
        args=_make_args(tmp_path),
        pipeline_source=_pipeline_source(),
        preflight_result=_preflight_result(),
        parsed_outputs=_parsed_outputs(preferred_h5ad=preferred),
        command_str="nextflow run nf-core/scrnaseq",
    )
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert preferred in content
    assert "scrna" in content


def test_write_report_handoff_without_preferred_h5ad(tmp_path):
    write_report(
        tmp_path,
        args=_make_args(tmp_path),
        pipeline_source=_pipeline_source(),
        preflight_result=_preflight_result(),
        parsed_outputs=_parsed_outputs(preferred_h5ad=""),
        command_str="nextflow run nf-core/scrnaseq",
    )
    content = (tmp_path / "report.md").read_text(encoding="utf-8")
    assert "not available" in content or "No canonical" in content


# ── write_repro_commands ──────────────────────────────────────────────────────

def test_write_repro_commands_creates_commands_sh(tmp_path):
    args = _make_args(tmp_path)
    write_repro_commands(tmp_path, args=args)
    assert (tmp_path / "reproducibility" / "commands.sh").exists()


def test_write_repro_commands_contains_output_path(tmp_path):
    args = _make_args(tmp_path, output=str(tmp_path))
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    # write_repro_commands uses the first (already-resolved) output_dir, not args.output
    assert tmp_path.as_posix() in content


def test_write_repro_commands_paths_have_no_backslashes(tmp_path):
    args = _make_args(tmp_path)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    # Exclude comment lines from the backslash check
    non_comment_lines = [l for l in content.splitlines() if not l.strip().startswith("#")]
    for line in non_comment_lines:
        assert "\\" not in line.replace(" \\", ""), f"Backslash in arg value: {line!r}"


def test_write_repro_commands_demo_omits_input(tmp_path):
    args = _make_args(tmp_path, demo=True)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--demo" in content
    assert "--input" not in content


def test_write_repro_commands_includes_expected_cells_when_set(tmp_path):
    """--expected-cells must appear in the repro script so replays are identical."""
    args = _make_args(tmp_path, expected_cells=5000)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--expected-cells" in content
    assert "5000" in content


def test_write_repro_commands_omits_expected_cells_when_none(tmp_path):
    args = _make_args(tmp_path, expected_cells=None)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--expected-cells" not in content


def test_write_repro_commands_uses_canonical_skip_flags_when_set(tmp_path):
    """Deprecated --skip-emptydrops replays as the canonical --skip-cellbender flag."""
    args = _make_args(tmp_path, skip_fastqc=True, skip_emptydrops=True, skip_multiqc=True)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--skip-fastqc" in content
    assert "--skip-cellbender" in content
    assert "--skip-emptydrops" not in content
    assert "--skip-multiqc" in content


def test_write_repro_commands_omits_skip_flags_when_not_set(tmp_path):
    args = _make_args(tmp_path, skip_fastqc=False, skip_emptydrops=False, skip_multiqc=False)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--skip-fastqc" not in content
    assert "--skip-emptydrops" not in content
    assert "--skip-multiqc" not in content


def test_write_repro_commands_includes_aligner_tuning_when_set(tmp_path):
    """Aligner-specific tuning flags must appear in the repro script for exact replay."""
    args = _make_args(tmp_path, star_feature="GeneFull", simpleaf_umi_resolution="cr-like-em", kb_workflow="lamanno")
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--star-feature" in content
    assert "GeneFull" in content
    assert "--simpleaf-umi-resolution" in content
    assert "cr-like-em" in content
    assert "--kb-workflow" in content
    assert "lamanno" in content


def test_write_repro_commands_includes_genome_when_set(tmp_path):
    args = _make_args(tmp_path, genome="GRCh38")
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--genome" in content
    assert "GRCh38" in content


def test_write_repro_commands_includes_save_flags_when_set(tmp_path):
    args = _make_args(tmp_path, save_reference=True, save_align_intermeds=True)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--save-reference" in content
    assert "--save-align-intermeds" in content


def test_write_repro_commands_includes_downstream_opt_in_when_set(tmp_path):
    args = _make_args(tmp_path, run_downstream=True)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--run-downstream" in content


def test_write_repro_commands_includes_starsolo_extras(tmp_path):
    args = _make_args(tmp_path, star_ignore_sjdbgtf=True, seq_center="CoreLab")
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--star-ignore-sjdbgtf" in content
    assert "--seq-center" in content
    assert "CoreLab" in content


def test_write_repro_commands_includes_cellrangerarc_reference(tmp_path):
    args = _make_args(tmp_path, cellrangerarc_reference="GRCh38-2024-A")
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--cellrangerarc-reference" in content
    assert "GRCh38-2024-A" in content


def test_write_repro_commands_includes_skip_cellranger_renaming(tmp_path):
    args = _make_args(tmp_path, skip_cellranger_renaming=True)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--skip-cellranger-renaming" in content


def test_write_repro_commands_includes_skip_cellrangermulti_vdjref(tmp_path):
    args = _make_args(tmp_path, skip_cellrangermulti_vdjref=True)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--skip-cellrangermulti-vdjref" in content


def test_write_repro_commands_includes_email_and_multiqc_title(tmp_path):
    args = _make_args(tmp_path, email="lab@example.com", multiqc_title="My QC")
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--email" in content
    assert "lab@example.com" in content
    assert "--multiqc-title" in content
    assert "My QC" in content


# ── write_result ──────────────────────────────────────────────────────────────

def test_write_result_creates_result_json(tmp_path):
    args = _make_args(tmp_path)
    preferred = tmp_path / "preferred.h5ad"
    result_path = write_result(
        tmp_path,
        args=args,
        pipeline_source=_pipeline_source(),
        parsed_outputs=_parsed_outputs(preferred_h5ad=str(preferred)),
        command_str="nextflow run nf-core/scrnaseq",
    )
    assert result_path.exists()
    payload = json.loads(result_path.read_text(encoding="utf-8"))
    assert "summary" in payload
    assert payload["summary"]["preset"] == "star"
    assert payload["summary"]["profile"] == "docker"
    assert payload["summary"]["output_artifacts"]["preferred_h5ad"] == str(preferred)
    assert payload["data"]["output_artifacts"] == payload["summary"]["output_artifacts"]


# ── write_check_result ────────────────────────────────────────────────────────

def test_write_check_result_creates_check_result_json(tmp_path):
    payload = {"ok": True, "skill": "nfcore-scrnaseq-pipeline", "preflight": {}}
    path = write_check_result(tmp_path, payload)
    assert path == tmp_path / "check_result.json"
    assert path.exists()
    loaded = json.loads(path.read_text(encoding="utf-8"))
    assert loaded["ok"] is True


def test_check_result_json_not_counted_as_materialized_entry(tmp_path):
    """check_result.json must not block a subsequent full run (_check_output_dir ignores it)."""
    import sys
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent))
    from preflight import _check_output_dir
    (tmp_path / "check_result.json").write_text('{"ok": true}', encoding="utf-8")
    _check_output_dir(tmp_path, resume=False)


def test_check_output_dir_rejects_existing_file(tmp_path):
    from preflight import _check_output_dir

    output_file = tmp_path / "not_a_directory"
    output_file.write_text("x", encoding="utf-8")
    with pytest.raises(Exception) as exc:
        _check_output_dir(output_file, resume=False)
    assert getattr(exc.value, "error_code", "") == "OUTPUT_DIR_NOT_WRITABLE"


def test_write_repro_commands_creates_remap_script(tmp_path):
    args = _make_args(tmp_path)
    write_repro_commands(tmp_path, args=args)
    assert (tmp_path / "reproducibility" / "remap_paths.py").exists()


def test_remap_script_is_executable_python(tmp_path):
    args = _make_args(tmp_path)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "remap_paths.py").read_text(encoding="utf-8")
    assert "#!/usr/bin/env python3" in content
    assert "remap_csv" in content
    assert "verify_paths" in content


def test_write_repro_commands_portability_notice_in_commands_sh(tmp_path):
    args = _make_args(tmp_path, demo=False)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "remap_paths.py" in content
    assert "portab" in content.lower() or "absolute" in content.lower() or "FASTQ" in content


def test_write_repro_commands_demo_skips_portability_notice(tmp_path):
    args = _make_args(tmp_path, demo=True)
    write_repro_commands(tmp_path, args=args)
    content = (tmp_path / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "remap_paths.py" not in content


def test_handoff_lines_use_absolute_clawbio_path():
    import sys as _sys
    _sys.path.insert(0, str(__import__("pathlib").Path(__file__).resolve().parent.parent))
    from reporting import _build_handoff_lines
    from schemas import SKILL_DIR
    expected_script = (SKILL_DIR.parent.parent / "clawbio.py").as_posix()
    lines = _build_handoff_lines("/data/combined_filtered_matrix.h5ad")
    handoff_line = next((l for l in lines if "clawbio.py" in l), None)
    assert handoff_line is not None, "No handoff line found"
    assert expected_script in handoff_line, \
        f"Expected absolute path {expected_script!r} in line: {handoff_line!r}"
