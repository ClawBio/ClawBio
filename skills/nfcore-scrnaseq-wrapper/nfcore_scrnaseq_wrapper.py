#!/usr/bin/env python3
from __future__ import annotations

import argparse
import hashlib
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import traceback
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from command_builder import build_nextflow_command
from errors import ErrorCode, SkillError
from executor import execute_nextflow
from outputs_parser import parse_outputs
from params_builder import build_effective_params, serialize_params_yaml, write_params_yaml
from pipeline_source import resolve_pipeline_source
from preflight import check_output_dir_available, run_preflight
from provenance import write_provenance_bundle
from reporting import write_check_result, write_report, write_repro_commands, write_result
from samplesheet_builder import validate_and_normalize_samplesheet
from schemas import (
    DEFAULT_PIPELINE_VERSION,
    DEFAULT_PRESET,
    DEFAULT_PROFILE,
    DEFAULT_TIMEOUT_SECONDS,
    SKILL_NAME,
    PRESET_ALIGNERS,
    SUPPORTED_PRESETS,
    SUPPORTED_PROFILES,
)

_DOWNSTREAM_HANDOFF_TIMEOUT_SECONDS = 60 * 60


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(description="Run the scrnaseq pipeline through a ClawBio wrapper.")
    parser.add_argument("--input", help="Path to a valid samplesheet.csv")
    parser.add_argument("--output", required=True, help="Output directory for the wrapper results")
    parser.add_argument("--demo", action="store_true", help="Run the pipeline demo with the upstream test profile")
    parser.add_argument("--check", action="store_true", help="Run preflight only and exit")
    parser.add_argument("--profile", default=DEFAULT_PROFILE, choices=sorted(SUPPORTED_PROFILES), help="Execution backend profile")
    parser.add_argument("--pipeline-version", default=DEFAULT_PIPELINE_VERSION, help="Remote pipeline tag/commit fallback")
    parser.add_argument("--preset", default=DEFAULT_PRESET, choices=sorted(SUPPORTED_PRESETS), help="Curated pipeline preset (standard = simpleaf/alevin-fry, star = STARsolo, kallisto = kb-python, cellranger/cellrangerarc/cellrangermulti = CellRanger variants)")
    parser.add_argument("--protocol", default=None, help="Protocol forwarded to the aligner (e.g. 10XV3, auto, dropseq, or any custom string). Use 'auto' only for cellranger/cellrangerarc/cellrangermulti. For standard/star/kallisto an explicit non-auto value is required.")
    parser.add_argument("--email", default=None, help="Email address for pipeline completion notification")
    parser.add_argument("--multiqc-title", default=None, help="Custom title for the MultiQC report")
    parser.add_argument("--expected-cells", type=int, default=None, help="Optional override for expected_cells")
    parser.add_argument("--resume", action="store_true", help="Attempt a compatible Nextflow resume")
    # Skip flags
    parser.add_argument("--skip-cellbender", action="store_true", help="Disable the cellbender subworkflow")
    parser.add_argument("--skip-fastqc", action="store_true", help="Skip the FastQC quality control step")
    parser.add_argument(
        "--skip-emptydrops",
        action="store_true",
        help="Deprecated alias for --skip-cellbender; kept for old commands.",
    )
    parser.add_argument("--skip-multiqc", action="store_true", help="Skip MultiQC report generation")
    parser.add_argument("--skip-cellranger-renaming", action="store_true", help="Skip automatic sample renaming in CellRanger modules")
    parser.add_argument("--skip-cellrangermulti-vdjref", action="store_true", help="Skip mkvdjref step in cellrangermulti (when no VDJ data)")
    # Reference genome
    parser.add_argument("--genome", default=None, help="iGenomes reference shortcut (e.g. GRCh38, mm10). Mutually exclusive with --fasta/--gtf.")
    parser.add_argument("--save-reference", action="store_true", help="Save the built reference index for future reuse")
    parser.add_argument("--save-align-intermeds", action="store_true", help="Save alignment intermediate files (BAMs)")
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--gtf", default=None)
    parser.add_argument("--transcript-fasta", default=None)
    parser.add_argument("--txp2gene", default=None)
    parser.add_argument("--simpleaf-index", default=None)
    parser.add_argument("--kallisto-index", default=None)
    parser.add_argument("--star-index", default=None)
    parser.add_argument("--cellranger-index", default=None)
    parser.add_argument("--barcode-whitelist", default=None)
    # STARsolo extras
    parser.add_argument("--star-feature", default=None, choices=["Gene", "GeneFull", "Gene Velocyto"], help="STARsolo feature type. 'Gene Velocyto' generates RNA velocity matrices.")
    parser.add_argument("--star-ignore-sjdbgtf", action="store_true", help="Do not use GTF for SJDB construction (use with --star-feature 'Gene Velocyto')")
    parser.add_argument("--seq-center", default=None, help="Sequencing center name for BAM read group tag")
    # Simpleaf extras
    parser.add_argument("--simpleaf-umi-resolution", default=None, choices=["cr-like", "cr-like-em", "parsimony", "parsimony-em", "parsimony-gene", "parsimony-gene-em"], help="UMI resolution strategy for alevin-fry")
    # Kallisto/BUS extras
    parser.add_argument("--kb-workflow", default=None, choices=["standard", "lamanno", "nac"], help="Kallisto workflow type")
    parser.add_argument("--kb-t1c", default=None, help="cDNA transcripts-to-capture file for RNA velocity (lamanno/nac workflows)")
    parser.add_argument("--kb-t2c", default=None, help="Intron transcripts-to-capture file for RNA velocity (lamanno/nac workflows)")
    # CellRanger ARC
    parser.add_argument("--motifs", default=None, help="Motif file (e.g. JASPAR) for CellRanger ARC index construction")
    parser.add_argument("--cellrangerarc-config", default=None, help="Config file for CellRanger ARC index construction")
    parser.add_argument("--cellrangerarc-reference", default=None, help="Reference genome name used inside the CellRanger ARC config file")
    # CellRanger Multi
    parser.add_argument("--cellranger-vdj-index", default=None, help="Pre-built CellRanger VDJ reference index")
    parser.add_argument("--gex-frna-probe-set", default=None, help="Probe set CSV for fixed RNA profiling (FFPE samples)")
    parser.add_argument("--gex-target-panel", default=None, help="Target panel CSV for targeted gene expression")
    parser.add_argument("--gex-cmo-set", default=None, help="Cell Multiplexing Oligo (CMO) reference CSV for multiplexed samples")
    parser.add_argument("--fb-reference", default=None, help="Feature barcoding reference CSV (e.g. antibody capture)")
    parser.add_argument("--vdj-inner-enrichment-primers", default=None, help="Text file with V(D)J cDNA enrichment primer sequences")
    parser.add_argument("--gex-barcode-sample-assignment", default=None, help="Barcode-to-sample assignment CSV to override CellRanger defaults")
    parser.add_argument("--cellranger-multi-barcodes", default=None, help="Additional samplesheet CSV with multiplexed sample information for cellrangermulti")
    parser.add_argument("--run-downstream", action="store_true", help="Opt in to running scrna_orchestrator after a canonical h5ad is detected")
    parser.add_argument("--skip-downstream", action="store_true", help="Compatibility flag; downstream handoff is skipped unless --run-downstream is set")
    return parser


def _write_demo_samplesheet(samplesheet_path: Path) -> None:
    samplesheet_path.parent.mkdir(parents=True, exist_ok=True)
    samplesheet_path.write_text("sample,fastq_1,fastq_2,expected_cells,seq_center\n", encoding="utf-8")


def _write_macos_docker_config(output_dir: Path) -> Path:
    """Write a Nextflow config with macOS + Apple Silicon Docker workarounds.

    Three issues on macOS with Docker (Colima / Docker Desktop):

    1. EDEADLK (errno 35): QEMU-emulated amd64 containers cannot exec files from
       VirtioFS bind-mounted paths. stageInMode="copy" forces Nextflow to copy
       inputs to the host work dir before launching each container, giving the
       guest kernel a clean page-cache state.

    2. ARM64 host / amd64 images: nf-core containers ship only linux/amd64 images.
       docker.runOptions = "--platform linux/amd64" tells Docker to pull and run
       the amd64 image under Rosetta / QEMU emulation on Apple Silicon.

    3. STAR FIFO (named-pipe) limitation: STAR creates named pipes under its work
       dir for streaming reads. VirtioFS does not support FIFOs, so STAR fails with
       "could not create FIFO file". --outTmpDir /tmp/star_tmp routes _STARtmp to
       the container's own /tmp (Linux tmpfs) where FIFOs are supported. The full
       ext.args from conf/modules.config is reproduced here so no flags are lost.
       Do NOT use a self-referencing closure (task.ext.args ?: "") — that causes a
       StackOverflowError.

    4. Time limit: macOS Docker (VirtioFS) adds significant I/O overhead. The nf-core
       test profile caps all tasks at 1 h, which is too short for STAR genome generation
       and alignment under emulation. resourceLimits.time is raised to 4 h here so the
       test profile limit is overridden (later configs take precedence in Nextflow).

    IMPORTANT — work directory must be under HOME:
    Colima (QEMU or VZ backend) only mounts the macOS HOME directory into the VM.
    /tmp and /private/tmp are NOT mounted; they map to the VM's own filesystem.
    Use --output under your home directory to avoid ".command.run: No such file"
    errors. Preflight emits a WARNING when --output is under /tmp.
    """
    config_path = output_dir / "reproducibility" / "macos_docker.config"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    config_path.write_text(
        '// macOS + Docker workaround for VirtioFS EDEADLK (errno 35), ARM64 hosts,\n'
        '// STAR FIFOs, and the nf-core test profile 1 h time cap.\n'
        'process {\n'
        '    stageInMode = "copy"\n'
        '    // VirtioFS I/O overhead makes STAR genome generation and alignment\n'
        '    // exceed the test profile\'s 1 h cap. Raise the ceiling to 4 h.\n'
        '    resourceLimits = [\n'
        '        cpus: 4,\n'
        '        memory: \'15.GB\',\n'
        '        time: \'4.h\'\n'
        '    ]\n'
        '    // STAR creates FIFOs in _STARtmp; VirtioFS does not support FIFOs.\n'
        '    // --outTmpDir /tmp/star_tmp routes _STARtmp to the container\'s /tmp\n'
        '    // (Linux tmpfs), which does support FIFOs.\n'
        '    // All other flags mirror nf-core/scrnaseq conf/modules.config verbatim\n'
        '    // so this override does not silently drop --readFilesCommand zcat.\n'
        '    // Note: do NOT use a closure that references task.ext.args — that\n'
        '    // causes a StackOverflowError. Instead we reproduce the flags here.\n'
        '    withName: \'.*STAR_ALIGN.*\' {\n'
        '        ext.args = { "--readFilesCommand zcat --runDirPerm All_RWX --outWigType bedGraph --twopassMode Basic --outSAMtype BAM SortedByCoordinate --limitBAMsortRAM ${task.memory.toBytes()} --outTmpDir /tmp/star_tmp" }\n'
        '    }\n'
        '}\n'
        'docker {\n'
        '    runOptions = "--platform linux/amd64"\n'
        '}\n',
        encoding="utf-8",
    )
    return config_path


def _write_error_result_if_safe(output_dir: Path, payload: dict[str, object]) -> None:
    """Persist structured errors unless doing so would overwrite a rejected output dir."""
    error_code = payload.get("error_code")
    if error_code in {"OUTPUT_DIR_NOT_EMPTY", "OUTPUT_DIR_NOT_WRITABLE"}:
        return
    if error_code == "INVALID_RESUME_STATE" and (output_dir / "result.json").exists():
        return
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        (output_dir / "result.json").write_text(json.dumps(payload, indent=2), encoding="utf-8")
    except OSError:
        # stderr still receives the structured payload; avoid masking the root error.
        return


def main() -> int:
    parser = build_parser()
    args = parser.parse_args()
    output_dir = Path(args.output).expanduser().resolve()

    try:
        return _run_wrapper(args, output_dir)
    except SkillError as exc:
        return _handle_skill_error(output_dir, exc)
    except Exception as exc:
        return _handle_unexpected_error(output_dir, exc)


def _run_wrapper(args: argparse.Namespace, output_dir: Path) -> int:
    if args.resume:
        with tempfile.TemporaryDirectory(prefix="clawbio-scrnaseq-resume-") as staging_dir:
            return _run_wrapper_with_staging(args, output_dir, staging_dir=Path(staging_dir))
    return _run_wrapper_with_staging(args, output_dir, staging_dir=None)


def _run_wrapper_with_staging(args: argparse.Namespace, output_dir: Path, *, staging_dir: Path | None) -> int:
    check_output_dir_available(output_dir, resume=args.resume)
    pipeline_source = resolve_pipeline_source(requested_version=args.pipeline_version)
    normalized_samplesheet, staged_samplesheet, samplesheet_summary = _prepare_samplesheet(
        args,
        output_dir,
        staging_dir=staging_dir,
    )
    preflight_result = run_preflight(args, pipeline_source=pipeline_source, samplesheet_summary=samplesheet_summary)
    if args.check:
        return _write_check_mode_result(output_dir, preflight_result=preflight_result, pipeline_source=pipeline_source)
    return _run_execution_mode(
        args,
        output_dir=output_dir,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        normalized_samplesheet=normalized_samplesheet,
        staged_samplesheet=staged_samplesheet,
        samplesheet_summary=samplesheet_summary,
    )


def _prepare_samplesheet(
    args: argparse.Namespace,
    output_dir: Path,
    *,
    staging_dir: Path | None,
) -> tuple[Path, Path, dict[str, object]]:
    if args.demo:
        return _prepare_demo_samplesheet(args, output_dir, staging_dir=staging_dir)
    return _prepare_user_samplesheet(args, output_dir, staging_dir=staging_dir)


def _prepare_demo_samplesheet(
    args: argparse.Namespace,
    output_dir: Path,
    *,
    staging_dir: Path | None,
) -> tuple[Path, Path, dict[str, object]]:
    # Apply demo overrides first so all subsequent callers see star.
    if args.preset != "star":
        print(
            f"WARNING: --demo forces preset=star (requested: {args.preset!r}). "
            "The nf-core test profile ships STAR-compatible data.",
            file=sys.stderr,
        )
    args.preset = "star"
    args.skip_cellbender = True

    normalized_samplesheet = _final_samplesheet_path(output_dir, demo=True)
    staged_samplesheet = _staged_samplesheet_path(output_dir, staging_dir=staging_dir, demo=True)
    _write_demo_samplesheet(staged_samplesheet)
    return normalized_samplesheet, staged_samplesheet, {
        "normalized_path": normalized_samplesheet,
        "sample_count": 0,
        "sample_names": [],
        "fastq_paths": [],
        "unknown_columns": [],
    }


def _prepare_user_samplesheet(
    args: argparse.Namespace,
    output_dir: Path,
    *,
    staging_dir: Path | None,
) -> tuple[Path, Path, dict[str, object]]:
    if not args.input:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_INPUT,
            message="An input samplesheet is required unless --demo is used.",
            fix="Provide --input <samplesheet.csv> or run with --demo.",
            details={},
        )
    normalized_samplesheet = _final_samplesheet_path(output_dir)
    staged_samplesheet = _staged_samplesheet_path(output_dir, staging_dir=staging_dir)
    samplesheet_summary = validate_and_normalize_samplesheet(
        Path(args.input).expanduser().resolve(),
        staged_samplesheet,
        expected_cells_override=args.expected_cells,
        preset=args.preset,
    )
    samplesheet_summary["normalized_path"] = normalized_samplesheet
    _warn_about_preserved_unknown_columns(samplesheet_summary)
    return normalized_samplesheet, staged_samplesheet, samplesheet_summary


def _final_samplesheet_path(output_dir: Path, *, demo: bool = False) -> Path:
    filename = "samplesheet.demo.csv" if demo else "samplesheet.valid.csv"
    return output_dir / "reproducibility" / filename


def _staged_samplesheet_path(output_dir: Path, *, staging_dir: Path | None, demo: bool = False) -> Path:
    if staging_dir is not None:
        filename = "samplesheet.demo.csv" if demo else "samplesheet.valid.csv"
        return staging_dir / filename
    return _final_samplesheet_path(output_dir, demo=demo)


def _warn_about_preserved_unknown_columns(samplesheet_summary: dict[str, object]) -> None:
    unknown_columns = samplesheet_summary.get("unknown_columns", [])
    if unknown_columns:
        print(
            f"WARNING: samplesheet contains unrecognised columns that will be preserved: {unknown_columns}",
            file=sys.stderr,
        )


def _write_check_mode_result(
    output_dir: Path,
    *,
    preflight_result: dict[str, object],
    pipeline_source: dict[str, object],
) -> int:
    payload = {
        "ok": True,
        "skill": SKILL_NAME,
        "preflight": preflight_result,
        "pipeline_source": pipeline_source,
    }
    write_check_result(output_dir, payload)
    print(json.dumps(payload, indent=2))
    return 0


def _run_execution_mode(
    args: argparse.Namespace,
    *,
    output_dir: Path,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, object],
    normalized_samplesheet: Path,
    staged_samplesheet: Path,
    samplesheet_summary: dict[str, object],
) -> int:
    params_payload = build_effective_params(
        args,
        normalized_samplesheet=normalized_samplesheet,
        output_dir=output_dir,
    )
    _check_resume_params_checksum(args, output_dir=output_dir, params_payload=params_payload)
    _commit_validated_samplesheet(staged_samplesheet, normalized_samplesheet)
    params_path = write_params_yaml(params_payload, output_dir=output_dir)
    command, command_str = _build_nextflow_invocation(args, output_dir, pipeline_source, params_path)
    nextflow_cwd = _nextflow_execution_cwd(output_dir)
    execution_result = execute_nextflow(
        command,
        cwd=nextflow_cwd,
        output_dir=output_dir,
        timeout_seconds=DEFAULT_TIMEOUT_SECONDS,
    )
    parsed_outputs = _parse_outputs_with_effective_aligner(output_dir, args)
    _write_success_outputs(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        params_path=params_path,
        params_payload=params_payload,
        normalized_samplesheet=normalized_samplesheet,
        samplesheet_summary=samplesheet_summary,
        parsed_outputs=parsed_outputs,
        execution_result=execution_result,
        command_str=_nextflow_replay_command(command_str, nextflow_cwd),
    )
    _run_downstream_handoff(args, parsed_outputs=parsed_outputs, output_dir=output_dir)
    print(f"Wrapper completed successfully. Output: {output_dir}")
    return 0


def _commit_validated_samplesheet(staged_samplesheet: Path, normalized_samplesheet: Path) -> None:
    if staged_samplesheet == normalized_samplesheet:
        return
    normalized_samplesheet.parent.mkdir(parents=True, exist_ok=True)
    shutil.copy2(staged_samplesheet, normalized_samplesheet)


def _check_resume_params_checksum(
    args: argparse.Namespace,
    *,
    output_dir: Path,
    params_payload: dict[str, object],
) -> None:
    if not args.resume:
        return
    manifest_path = output_dir / "reproducibility" / "manifest.json"
    if not manifest_path.exists():
        return
    params_checksum = _params_payload_checksum(params_payload)
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("params_checksum") == params_checksum:
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_RESUME_STATE,
        message="Resume state does not match the effective params.yaml for this run.",
        fix="Use the same arguments as the original run or remove --resume.",
        details={
            "previous_params_checksum": manifest.get("params_checksum"),
            "requested_params_checksum": params_checksum,
        },
    )


def _params_payload_checksum(params_payload: dict[str, object]) -> str:
    return hashlib.sha256(serialize_params_yaml(params_payload).encode("utf-8")).hexdigest()


def _build_nextflow_invocation(
    args: argparse.Namespace,
    output_dir: Path,
    pipeline_source: dict[str, object],
    params_path: Path,
) -> tuple[list[str], str]:
    return build_nextflow_command(
        pipeline_source=pipeline_source,
        profile=f"test,{args.profile}" if args.demo else args.profile,
        params_path=params_path,
        resume=args.resume,
        work_dir=output_dir / "upstream" / "work",
        extra_configs=_build_extra_nextflow_configs(args, output_dir),
    )


def _build_extra_nextflow_configs(args: argparse.Namespace, output_dir: Path) -> list[Path]:
    # The VirtioFS EDEADLK fix is macOS-specific:
    #   Linux native: no VirtioFS; Docker bind-mounts are direct overlayfs.
    #   Windows/WSL2: Nextflow runs inside WSL2 (reports "Linux"); work dir is on ext4.
    #   Windows native: not officially supported by Nextflow; users should use WSL2.
    if platform.system() == "Darwin" and args.profile == "docker":
        return [_write_macos_docker_config(output_dir)]
    return []


def _nextflow_execution_cwd(output_dir: Path) -> Path:
    # params.input is written relative to output_dir so it stays schema-safe even
    # when the absolute workspace path contains spaces.
    return output_dir


def _nextflow_replay_command(command_str: str, cwd: Path) -> str:
    return f"cd {shlex.quote(cwd.as_posix())} && {command_str}"


def _parse_outputs_with_effective_aligner(output_dir: Path, args: argparse.Namespace) -> dict[str, object]:
    return {
        **parse_outputs(output_dir),
        "aligner_effective": PRESET_ALIGNERS[args.preset],
    }


def _write_success_outputs(
    output_dir: Path,
    *,
    args: argparse.Namespace,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, object],
    params_path: Path,
    params_payload: dict[str, object],
    normalized_samplesheet: Path,
    samplesheet_summary: dict[str, object],
    parsed_outputs: dict[str, object],
    execution_result: dict[str, object],
    command_str: str,
) -> None:
    write_repro_commands(output_dir, args=args)
    write_provenance_bundle(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        params_path=params_path,
        params_payload=params_payload,
        normalized_samplesheet=normalized_samplesheet,
        samplesheet_summary=samplesheet_summary,
        parsed_outputs=parsed_outputs,
        execution_result=execution_result,
        command_str=command_str,
    )
    write_report(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
    )
    write_result(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
    )


_SCRNA_ORCHESTRATOR_DEFAULT = _SKILL_DIR.parent / "scrna-orchestrator" / "scrna_orchestrator.py"


def _resolve_scrna_orchestrator() -> Path | None:
    env_path = os.environ.get("CLAWBIO_SCRNA_ORCHESTRATOR")
    if env_path:
        p = Path(env_path)
        return p if p.exists() else None
    return _SCRNA_ORCHESTRATOR_DEFAULT if _SCRNA_ORCHESTRATOR_DEFAULT.exists() else None


def _run_downstream_handoff(
    args: argparse.Namespace,
    *,
    parsed_outputs: dict[str, object],
    output_dir: Path,
) -> None:
    if not getattr(args, "run_downstream", False) or getattr(args, "skip_downstream", False):
        return
    preferred_h5ad = str(parsed_outputs.get("preferred_h5ad", "")).strip()
    if not preferred_h5ad:
        return
    orchestrator = _resolve_scrna_orchestrator()
    if orchestrator is None:
        print(
            "WARNING: scrna_orchestrator not found. "
            "Set CLAWBIO_SCRNA_ORCHESTRATOR=/path/to/scrna_orchestrator.py to enable automatic handoff.",
            file=sys.stderr,
        )
        return
    downstream_output = output_dir / "scrna_analysis"
    cmd = [
        sys.executable,
        str(orchestrator),
        "--input", preferred_h5ad,
        "--output", str(downstream_output),
    ]
    print("Handing off to scrna_orchestrator...")
    try:
        result = subprocess.run(
            cmd,
            capture_output=True,
            text=True,
            timeout=_DOWNSTREAM_HANDOFF_TIMEOUT_SECONDS,
            cwd=str(_PROJECT_ROOT),
        )
    except subprocess.TimeoutExpired:
        print(
            "WARNING: downstream scrna_orchestrator timed out. "
            "The nf-core pipeline completed successfully; run downstream analysis manually if needed.",
            file=sys.stderr,
        )
        return
    if result.returncode != 0:
        print(
            f"WARNING: downstream scrna_orchestrator exited with code {result.returncode}. "
            "The nf-core pipeline completed successfully; inspect scrna_analysis/ manually.",
            file=sys.stderr,
        )


def _handle_skill_error(output_dir: Path, exc: SkillError) -> int:
    payload = exc.to_dict()
    _write_error_result_if_safe(output_dir, payload)
    print(json.dumps(payload, indent=2), file=sys.stderr)
    return 1


def _handle_unexpected_error(output_dir: Path, exc: Exception) -> int:
    payload = {
        "ok": False,
        "stage": "internal",
        "error_code": "UNEXPECTED_ERROR",
        "message": str(exc),
        "fix": "Report this as a bug. Include the full traceback and your command arguments.",
        "details": {
            "exception_type": type(exc).__name__,
            "traceback": traceback.format_exc(),
        },
    }
    _write_error_result_if_safe(output_dir, payload)
    print(json.dumps(payload, indent=2), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
