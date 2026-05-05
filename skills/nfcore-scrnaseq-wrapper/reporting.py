from __future__ import annotations

import json
import shutil
import sys
from datetime import datetime, timezone
from pathlib import Path

from clawbio.common.portable_commands import write_portable_commands_sh
from clawbio.common.report import generate_report_footer, generate_report_header, write_result_json

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from schemas import SKILL_ALIAS, SKILL_DIR, SKILL_NAME, SKILL_VERSION

_CLAWBIO_SCRIPT = (SKILL_DIR.parent.parent / "clawbio.py").as_posix()

_REPRO_PATH_FLAGS = (
    "fasta",
    "gtf",
    "transcript_fasta",
    "txp2gene",
    "simpleaf_index",
    "kallisto_index",
    "star_index",
    "cellranger_index",
    "barcode_whitelist",
    "kb_t1c",
    "kb_t2c",
    "motifs",
    "cellrangerarc_config",
    "cellranger_vdj_index",
    "gex_frna_probe_set",
    "gex_target_panel",
    "gex_cmo_set",
    "fb_reference",
    "vdj_inner_enrichment_primers",
    "gex_barcode_sample_assignment",
    "cellranger_multi_barcodes",
)


def write_report(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, object],
    parsed_outputs: dict[str, object],
    command_str: str,
) -> Path:
    lines = build_report_lines(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
    )
    report_path = output_dir / "report.md"
    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


def build_report_lines(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, object],
    parsed_outputs: dict[str, object],
    command_str: str,
) -> list[str]:
    header = generate_report_header(
        "nf-core/scrnaseq Wrapper Report",
        SKILL_NAME,
        extra_metadata={
            "Preset": args.preset,
            "Profile": args.profile,
            "Pipeline source": str(pipeline_source["source_kind"]),
            "Pipeline ref": str(pipeline_source["resolved_version"]),
        },
    )
    preferred_h5ad = str(parsed_outputs.get("preferred_h5ad", ""))
    return [
        header,
        "## Summary",
        "",
        f"- Preset: `{args.preset}`",
        f"- Effective aligner: `{parsed_outputs.get('aligner_effective', '')}`",
        f"- Pipeline source: `{pipeline_source['source_kind']}`",
        f"- Pipeline ref: `{pipeline_source['resolved_version']}`",
        f"- Output root: `{output_dir}`",
        "",
        "## Preflight",
        "",
        f"- Java: `{preflight_result['java']['version']}`",
        f"- Nextflow: `{preflight_result['nextflow']['version']}`",
        f"- Backend profile: `{args.profile}`",
        f"- Samples: `{len(parsed_outputs.get('samples_detected', []))}`",
        "",
        "## Outputs",
        "",
        f"- Preferred h5ad: `{preferred_h5ad or 'not available'}`",
        f"- MultiQC report: `{parsed_outputs.get('multiqc_report', '') or 'not found'}`",
        f"- Pipeline info: `{parsed_outputs.get('pipeline_info_dir', '') or 'not found'}`",
        f"- CellBender detected: `{parsed_outputs.get('cellbender_used', False)}`",
        "",
        "## Reproducibility",
        "",
        f"- Command: `{command_str}`",
        f"- Repro bundle: `{output_dir / 'reproducibility'}`",
        "",
        *_build_handoff_lines(preferred_h5ad),
        generate_report_footer(),
    ]


def _build_handoff_lines(preferred_h5ad: str) -> list[str]:
    if preferred_h5ad:
        return [
            "## Next Steps",
            "",
            f"- `python {_CLAWBIO_SCRIPT} run scrna --input {preferred_h5ad} --output <dir>`",
            f"- `python {_CLAWBIO_SCRIPT} run scrna-embedding --input {preferred_h5ad} --output <dir>`",
            "",
        ]
    return [
        "## Next Steps",
        "",
        "- No canonical `.h5ad` was selected automatically. Inspect `result.json` and the upstream outputs before chaining downstream skills.",
        "",
    ]


_PORTABILITY_NOTICE = """\

# ── Portability notice (FASTQ paths) ─────────────────────────────────────────
# FASTQ paths in samplesheet.valid.csv are absolute paths (required by
# Nextflow). Before replaying on a different machine, update them with:
#
#   python reproducibility/remap_paths.py --old /original/prefix --new /new/prefix
#
# Preview changes without modifying files:
#   python reproducibility/remap_paths.py --old /original/prefix --new /new/prefix --dry-run
#
# Verify all current paths exist on this machine:
#   python reproducibility/remap_paths.py --verify
"""

_REMAP_SCRIPT_SRC = SKILL_DIR / "remap_paths.py"


def write_repro_commands(
    output_dir: Path,
    *,
    args,
) -> None:
    repro_dir = output_dir / "reproducibility"
    command_args = build_repro_command_args(output_dir, args=args)
    write_portable_commands_sh(
        repro_dir,
        skill_name=SKILL_NAME,
        script_name="nfcore_scrnaseq_wrapper.py",
        args=command_args,
        generated_at=datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M UTC"),
    )
    if not getattr(args, "demo", False):
        commands_sh = repro_dir / "commands.sh"
        with commands_sh.open("a", encoding="utf-8") as fh:
            fh.write(_PORTABILITY_NOTICE)
    _write_remap_script(repro_dir)


def _write_remap_script(repro_dir: Path) -> None:
    repro_dir.mkdir(parents=True, exist_ok=True)
    shutil.copy2(_REMAP_SCRIPT_SRC, repro_dir / "remap_paths.py")


def build_repro_command_args(output_dir: Path, *, args) -> dict[str, str | None]:
    command_args = _base_repro_command_args(output_dir, args=args)
    _add_optional_value_flags(command_args, args)
    _add_optional_boolean_flags(command_args, args)
    _add_repro_path_flags(command_args, args)
    return command_args


def _base_repro_command_args(output_dir: Path, *, args) -> dict[str, str | None]:
    command_args: dict[str, str | None] = {
        "--output": output_dir.as_posix(),
        "--preset": args.preset,
        "--profile": args.profile,
        "--pipeline-version": args.pipeline_version,
    }
    if args.demo:
        command_args["--demo"] = None
    else:
        # Resolve to absolute + forward slashes so the repro script works
        # when pasted into bash on any platform (including Windows Git Bash / WSL2).
        command_args["--input"] = Path(args.input).expanduser().resolve().as_posix()
    return command_args


def _add_optional_value_flags(command_args: dict[str, str | None], args) -> None:
    if args.protocol:
        command_args["--protocol"] = args.protocol
    if getattr(args, "email", None):
        command_args["--email"] = args.email
    if getattr(args, "multiqc_title", None):
        command_args["--multiqc-title"] = args.multiqc_title
    if getattr(args, "expected_cells", None) is not None and not args.demo:
        command_args["--expected-cells"] = str(args.expected_cells)
    if getattr(args, "skip_cellbender", False) or getattr(args, "skip_emptydrops", False):
        command_args["--skip-cellbender"] = None
    if getattr(args, "skip_fastqc", False):
        command_args["--skip-fastqc"] = None
    if getattr(args, "skip_multiqc", False):
        command_args["--skip-multiqc"] = None
    if getattr(args, "skip_cellranger_renaming", False):
        command_args["--skip-cellranger-renaming"] = None
    if getattr(args, "skip_cellrangermulti_vdjref", False):
        command_args["--skip-cellrangermulti-vdjref"] = None
    if getattr(args, "star_ignore_sjdbgtf", False):
        command_args["--star-ignore-sjdbgtf"] = None
    if getattr(args, "seq_center", None):
        command_args["--seq-center"] = args.seq_center
    if getattr(args, "star_feature", None):
        command_args["--star-feature"] = args.star_feature
    if getattr(args, "simpleaf_umi_resolution", None):
        command_args["--simpleaf-umi-resolution"] = args.simpleaf_umi_resolution
    if getattr(args, "kb_workflow", None):
        command_args["--kb-workflow"] = args.kb_workflow
    if getattr(args, "genome", None):
        command_args["--genome"] = args.genome
    if getattr(args, "cellrangerarc_reference", None):
        command_args["--cellrangerarc-reference"] = args.cellrangerarc_reference


def _add_optional_boolean_flags(command_args: dict[str, str | None], args) -> None:
    if getattr(args, "save_reference", False):
        command_args["--save-reference"] = None
    if getattr(args, "save_align_intermeds", False):
        command_args["--save-align-intermeds"] = None
    if args.resume:
        command_args["--resume"] = None
    if getattr(args, "run_downstream", False):
        command_args["--run-downstream"] = None


def _add_repro_path_flags(command_args: dict[str, str | None], args) -> None:
    for flag_name in _REPRO_PATH_FLAGS:
        value = getattr(args, flag_name, None)
        if value:
            command_args[f"--{flag_name.replace('_', '-')}"] = Path(value).expanduser().resolve().as_posix()


def write_result(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, object],
    parsed_outputs: dict[str, object],
    command_str: str,
) -> Path:
    output_artifacts = build_output_artifacts(parsed_outputs)
    summary = build_result_summary(
        args=args,
        pipeline_source=pipeline_source,
        parsed_outputs=parsed_outputs,
        output_artifacts=output_artifacts,
    )
    data = build_result_data(parsed_outputs=parsed_outputs, output_artifacts=output_artifacts, command_str=command_str)
    return write_result_json(
        output_dir,
        skill=SKILL_ALIAS,
        version=SKILL_VERSION,
        summary=summary,
        data=data,
    )


def build_output_artifacts(parsed_outputs: dict[str, object]) -> dict[str, object]:
    return {
        "preferred_h5ad": parsed_outputs.get("preferred_h5ad", ""),
        "multiqc_report": parsed_outputs.get("multiqc_report", ""),
        "pipeline_info_dir": parsed_outputs.get("pipeline_info_dir", ""),
        "h5ad_candidates": parsed_outputs.get("h5ad_candidates", []),
        "rds_candidates": parsed_outputs.get("rds_candidates", []),
    }


def build_result_summary(
    *,
    args,
    pipeline_source: dict[str, object],
    parsed_outputs: dict[str, object],
    output_artifacts: dict[str, object],
) -> dict[str, object]:
    return {
        "preset": args.preset,
        "aligner_effective": parsed_outputs.get("aligner_effective", ""),
        "pipeline_source_kind": pipeline_source["source_kind"],
        "pipeline_version_or_commit": pipeline_source["resolved_version"],
        "profile": args.profile,
        "resume_used": bool(args.resume),
        "multiqc_report": parsed_outputs.get("multiqc_report", ""),
        "pipeline_info_dir": parsed_outputs.get("pipeline_info_dir", ""),
        "preferred_h5ad": parsed_outputs.get("preferred_h5ad", ""),
        "handoff_available": parsed_outputs.get("handoff_available", False),
        "samples_detected": len(parsed_outputs.get("samples_detected", [])),
        "cellbender_used": parsed_outputs.get("cellbender_used", False),
        "output_artifacts": output_artifacts,
    }


def build_result_data(
    *,
    parsed_outputs: dict[str, object],
    output_artifacts: dict[str, object],
    command_str: str,
) -> dict[str, object]:
    return {
        "canonical_skill_name": SKILL_NAME,
        "cli_alias": SKILL_ALIAS,
        "command": command_str,
        "output_artifacts": output_artifacts,
        "outputs": parsed_outputs,
    }


def write_check_result(output_dir: Path, payload: dict[str, object]) -> Path:
    path = output_dir / "check_result.json"
    path.write_text(json.dumps(payload, indent=2), encoding="utf-8")
    return path
