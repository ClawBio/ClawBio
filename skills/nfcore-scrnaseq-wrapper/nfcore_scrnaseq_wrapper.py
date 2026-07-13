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
from typing import Any

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules(
    "command_builder",
    "errors",
    "executor",
    "nfcore_4_1_0_contract",
    "outputs_parser",
    "params_builder",
    "pipeline_source",
    "preflight",
    "provenance",
    "reporting",
    "repro_commands",
    "samplesheet_builder",
    "schemas",
)

from clawbio.common.textio import write_text_lf
from command_builder import build_nextflow_command
from errors import ErrorCode, SkillError
from executor import execute_nextflow
from nfcore_4_1_0_contract import (
    INTENTIONALLY_UNSUPPORTED_PARAMS,
    INTENTIONALLY_UNSUPPORTED_REASONS,
    NFCORE_SCRNASEQ_VERSION,
    OFFICIAL_PARAMS,
    WRAPPER_SUPPORTED_UPSTREAM_PARAMS,
)
from outputs_parser import parse_outputs, validate_expected_outputs
from params_builder import (
    build_effective_params,
    serialize_params_yaml,
    write_params_yaml,
)
from pipeline_source import resolve_pipeline_source
from preflight import check_output_dir_available, run_preflight
from provenance import write_provenance_bundle
from reporting import (
    write_check_result,
    write_report,
    write_repro_commands,
    write_result,
)
from repro_commands import write_macos_docker_config
from samplesheet_builder import validate_and_normalize_samplesheet
from schemas import (
    ALL_REFERENCE_PATH_FIELDS,
    DEFAULT_PIPELINE_VERSION,
    DEFAULT_PRESET,
    DEFAULT_PROFILE,
    DEFAULT_TIMEOUT_SECONDS,
    SKILL_NAME,
    PRESET_ALIGNERS,
    SUPPORTED_PRESETS,
    SYMBOLIC_REFERENCE_FIELDS,
    profile_includes,
)

_DOWNSTREAM_HANDOFF_TIMEOUT_SECONDS = 60 * 60
_DEFAULT_TIMEOUT_HOURS = DEFAULT_TIMEOUT_SECONDS / 3600

# Pipeline-parameter flags the upstream `test` profile owns; supplying any with
# --demo is a no-op (the wrapper never writes them in demo mode, audit H-01) and the
# user is warned. This is the full set the hermetic demo ignores: references +
# protocol + every QC/skip/tuning/save/reporting knob. ``skip_cellbender`` is NOT
# here — demo FORCES it on (the test data is too small for CellBender), so it is a
# wrapper-set essential, not an ignored flag.
_DEMO_IGNORED_FIELDS = (
    *SYMBOLIC_REFERENCE_FIELDS,
    "igenomes_base",
    *ALL_REFERENCE_PATH_FIELDS,
    "protocol",
    # aligner tuning
    "star_feature",
    "star_ignore_sjdbgtf",
    "seq_center",
    "simpleaf_umi_resolution",
    "kb_workflow",
    # skip flags (skip_cellbender is forced by demo, never ignored)
    "skip_fastqc",
    "skip_multiqc",
    "skip_cellranger_renaming",
    "skip_cellrangermulti_vdjref",
    # save flags
    "save_reference",
    "save_align_intermeds",
    # input metadata + generic reporting
    "email",
    "email_on_fail",
    "multiqc_title",
    "multiqc_config",
    "multiqc_logo",
    "multiqc_methods_description",
    "publish_dir_mode",
    "trace_report_suffix",
    "monochrome_logs",
)


def _timeout_hours(value: str) -> float:
    try:
        hours = float(value)
    except (TypeError, ValueError):
        raise argparse.ArgumentTypeError(
            f"--timeout-hours must be a number, got {value!r}"
        )
    if hours < 0:
        raise argparse.ArgumentTypeError(
            "--timeout-hours must be >= 0 (0 disables the timeout)"
        )
    return hours


def _resolve_timeout_seconds(args: argparse.Namespace) -> int | None:
    """Translate --timeout-hours into seconds; 0 disables the cap (returns None)."""
    hours = getattr(args, "timeout_hours", _DEFAULT_TIMEOUT_HOURS)
    if hours is None or hours == 0:
        return None
    return int(round(hours * 3600))


class _PresetAction(argparse.Action):
    def __call__(self, parser, namespace, values, option_string=None):
        setattr(namespace, self.dest, values)
        setattr(namespace, "preset_explicit", True)


def _demo_flag_supplied(args: argparse.Namespace, name: str) -> bool:
    """True when a demo-ignored flag was actually supplied on the CLI."""
    value = getattr(args, name, None)
    if name == "save_align_intermeds":
        # BooleanOptionalAction: None means unset; both --save-… and --no-save-…
        # (True/False) count as supplied and are equally ignored under --demo.
        return value is not None
    return bool(value)


def _demo_ignored_flags(args: argparse.Namespace) -> list[str]:
    """List the pipeline-parameter CLI flags supplied under --demo (all ignored)."""
    return [
        f"--{name.replace('_', '-')}"
        for name in _DEMO_IGNORED_FIELDS
        if _demo_flag_supplied(args, name)
    ]


def build_parser() -> argparse.ArgumentParser:
    parser = argparse.ArgumentParser(
        description="Run the scrnaseq pipeline through a ClawBio wrapper."
    )
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    parser.add_argument(
        "--no-banner", action="store_true", help="Suppress the startup banner"
    )
    parser.add_argument("--input", help="Path to a valid samplesheet.csv")
    parser.add_argument(
        "--output", required=True, help="Output directory for the wrapper results"
    )
    parser.add_argument(
        "--demo",
        action="store_true",
        help="Run the pipeline demo with the upstream test profile",
    )
    parser.add_argument(
        "--check", action="store_true", help="Run preflight only and exit"
    )
    parser.add_argument(
        "--profile",
        default=DEFAULT_PROFILE,
        help="Execution backend profile or comma-separated nf-core profile list",
    )
    parser.add_argument(
        "--pipeline-version",
        default=DEFAULT_PIPELINE_VERSION,
        help="Remote pipeline tag/commit fallback",
    )
    parser.add_argument(
        "--allow-dirty-pipeline",
        action="store_true",
        help="Allow an intentionally modified local nf-core/scrnaseq checkout; disabled by default for reproducibility",
    )
    parser.add_argument(
        "--require-local-pipeline",
        action="store_true",
        help="Require a verifiable local sibling nf-core/scrnaseq checkout instead of falling back to the remote pipeline",
    )
    parser.add_argument(
        "--allow-pipeline-version-override",
        action="store_true",
        help=(
            f"Allow --pipeline-version other than the pinned {NFCORE_SCRNASEQ_VERSION} contract. "
            "The wrapper's parameter/protocol/output validations remain 4.1.0-specific, so this is at your own risk."
        ),
    )
    parser.add_argument(
        "--trust-config-params",
        action="store_true",
        help=(
            "Allow -c/--config files that set params.* (otherwise blocked). nf-core advises against "
            "setting params in configs; detected overrides are recorded in provenance."
        ),
    )
    parser.add_argument(
        "-c",
        "--config",
        "--nextflow-config",
        dest="extra_config",
        action="append",
        default=[],
        help="Additional Nextflow config file (alias: --nextflow-config); may be supplied multiple times for HPC/cloud profiles",
    )
    parser.add_argument(
        "--preset",
        default=DEFAULT_PRESET,
        action=_PresetAction,
        choices=sorted(SUPPORTED_PRESETS),
        help="Curated pipeline preset (standard = simpleaf/alevin-fry, star = STARsolo, kallisto = kb-python, cellranger/cellrangerarc/cellrangermulti = CellRanger variants)",
    )
    parser.set_defaults(preset_explicit=False)
    parser.add_argument(
        "--aligner",
        default=None,
        choices=sorted(PRESET_ALIGNERS.values()),
        help="nf-core/scrnaseq aligner name. Alias for --preset; simpleaf maps to preset=standard.",
    )
    parser.add_argument(
        "--protocol",
        default=None,
        help="Protocol forwarded to the aligner. Required for standard/star/kallisto; auto is CellRanger-only; other values are mapped when known or passed through by nf-core.",
    )
    parser.add_argument(
        "--email",
        default=None,
        help="Email address for pipeline completion notification",
    )
    parser.add_argument(
        "--email-on-fail",
        default=None,
        help="Email address for pipeline failure notification only",
    )
    parser.add_argument(
        "--multiqc-title", default=None, help="Custom title for the MultiQC report"
    )
    parser.add_argument("--multiqc-config", default=None, help="MultiQC config file")
    parser.add_argument("--multiqc-logo", default=None, help="MultiQC logo image")
    parser.add_argument(
        "--multiqc-methods-description",
        default=None,
        help="MultiQC methods description YAML",
    )
    parser.add_argument(
        "--publish-dir-mode",
        default=None,
        choices=["symlink", "rellink", "link", "copy", "copyNoFollow", "move"],
        help="Nextflow publishDir mode",
    )
    parser.add_argument(
        "--trace-report-suffix",
        default=None,
        help="Suffix appended to Nextflow trace report filenames",
    )
    parser.add_argument(
        "--monochrome-logs",
        action="store_true",
        help="Disable ANSI colors in Nextflow logs",
    )
    parser.add_argument(
        "--expected-cells",
        type=int,
        default=None,
        help="Global override for expected_cells. Single-sample only: for multi-sample sheets set expected_cells per-row instead (a global value is rejected to avoid applying one count to heterogeneous samples).",
    )
    parser.add_argument(
        "--timeout-hours",
        type=_timeout_hours,
        default=_DEFAULT_TIMEOUT_HOURS,
        help="Wall-clock cap for the Nextflow run in hours (default 12). Use 0 to disable the cap for long HPC/cloud runs whose walltime is enforced by the scheduler.",
    )
    parser.add_argument(
        "--work-dir",
        default=None,
        help="Nextflow work directory. Defaults to <output>/upstream/work; may be an object-store URI for cloud executors.",
    )
    parser.add_argument(
        "--allow-remote-inputs",
        action="store_true",
        help=(
            "Opt in to remote samplesheet inputs and reference paths (s3://, gs://, "
            "https://, ftp://, …). Default is local-first: remote URIs are rejected so "
            "genetic data and references stay on the local machine. When enabled, a "
            "runtime warning lists the paths fetched over the network."
        ),
    )
    parser.add_argument(
        "--allow-conda-cellranger",
        action="store_true",
        help="Allow Cell Ranger presets with conda/mamba only when a trusted site config supplies Cell Ranger.",
    )
    parser.add_argument(
        "--resume", action="store_true", help="Attempt a compatible Nextflow resume"
    )
    # Skip flags
    parser.add_argument(
        "--skip-cellbender",
        action="store_true",
        help="Disable the cellbender subworkflow",
    )
    parser.add_argument(
        "--skip-fastqc",
        action="store_true",
        help="Skip the FastQC quality control step",
    )
    parser.add_argument(
        "--skip-emptydrops",
        action="store_true",
        help="Deprecated alias for --skip-cellbender; kept for old commands.",
    )
    parser.add_argument(
        "--skip-multiqc", action="store_true", help="Skip MultiQC report generation"
    )
    parser.add_argument(
        "--skip-cellranger-renaming",
        action="store_true",
        help="Skip automatic sample renaming in CellRanger modules",
    )
    parser.add_argument(
        "--skip-cellrangermulti-vdjref",
        action="store_true",
        help="Skip mkvdjref step in cellrangermulti (when no VDJ data)",
    )
    # Reference genome
    parser.add_argument(
        "--genome",
        default=None,
        help="iGenomes reference shortcut (e.g. GRCh38, mm10). Mutually exclusive with --fasta/--gtf.",
    )
    parser.add_argument(
        "--igenomes-base",
        default=None,
        help="Base path/URL for iGenomes references (e.g. a local mirror); overrides the default S3 location.",
    )
    parser.add_argument(
        "--igenomes-ignore",
        action="store_true",
        help="Do not load the iGenomes reference config.",
    )
    parser.add_argument(
        "--save-reference",
        action="store_true",
        help="Save the built reference index for future reuse",
    )
    parser.add_argument(
        "--save-align-intermeds",
        action=argparse.BooleanOptionalAction,
        default=None,
        help="Save alignment intermediate files (BAMs). Use --no-save-align-intermeds to force off; unset defers to nf-core.",
    )
    parser.add_argument(
        "--fasta",
        default=None,
        help="Genome FASTA. Path must contain NO whitespace and end in .fa/.fna/.fasta(.gz) (nf-core/scrnaseq 4.1.0 schema). Mutually exclusive with --genome.",
    )
    parser.add_argument(
        "--gtf",
        default=None,
        help="Gene annotation GTF (paired with --fasta to build the index). Mutually exclusive with --genome.",
    )
    parser.add_argument(
        "--transcript-fasta",
        default=None,
        help="Transcriptome FASTA for simpleaf (used with --txp2gene).",
    )
    parser.add_argument(
        "--txp2gene",
        default=None,
        help="Transcript-to-gene mapping for simpleaf (used with --transcript-fasta).",
    )
    parser.add_argument(
        "--simpleaf-index",
        default=None,
        help="Prebuilt simpleaf/alevin-fry index directory.",
    )
    parser.add_argument(
        "--kallisto-index", default=None, help="Prebuilt kallisto index."
    )
    parser.add_argument(
        "--star-index", default=None, help="Prebuilt STAR genome index directory."
    )
    parser.add_argument(
        "--cellranger-index",
        default=None,
        help="Prebuilt CellRanger / CellRanger ARC reference directory.",
    )
    parser.add_argument(
        "--barcode-whitelist",
        default=None,
        help="Custom barcode whitelist (per-aligner format). Auxiliary file; does not replace the genome reference.",
    )
    # STARsolo extras
    parser.add_argument(
        "--star-feature",
        default=None,
        choices=["Gene", "GeneFull", "Gene Velocyto"],
        help="STARsolo feature type. 'Gene Velocyto' generates RNA velocity matrices.",
    )
    parser.add_argument(
        "--star-ignore-sjdbgtf",
        action="store_true",
        help="Do not use GTF for SJDB construction (use with --star-feature 'Gene Velocyto')",
    )
    parser.add_argument(
        "--seq-center",
        default=None,
        help="Sequencing center name for BAM read group tag",
    )
    # Simpleaf extras
    parser.add_argument(
        "--simpleaf-umi-resolution",
        default=None,
        choices=[
            "cr-like",
            "cr-like-em",
            "parsimony",
            "parsimony-em",
            "parsimony-gene",
            "parsimony-gene-em",
        ],
        help="UMI resolution strategy for alevin-fry",
    )
    # Kallisto/BUS extras
    parser.add_argument(
        "--kb-workflow",
        default=None,
        choices=["standard", "lamanno", "nac"],
        help="Kallisto workflow type",
    )
    parser.add_argument(
        "--kb-t1c",
        default=None,
        help="cDNA transcripts-to-capture file for RNA velocity (lamanno/nac workflows)",
    )
    parser.add_argument(
        "--kb-t2c",
        default=None,
        help="Intron transcripts-to-capture file for RNA velocity (lamanno/nac workflows)",
    )
    # CellRanger ARC
    parser.add_argument(
        "--motifs",
        default=None,
        help="Motif file (e.g. JASPAR) for CellRanger ARC index construction",
    )
    parser.add_argument(
        "--cellrangerarc-config",
        default=None,
        help="Config file for CellRanger ARC index construction",
    )
    parser.add_argument(
        "--cellrangerarc-reference",
        default=None,
        help="Reference genome name used inside the CellRanger ARC config file",
    )
    # CellRanger Multi
    parser.add_argument(
        "--cellranger-vdj-index",
        default=None,
        help="Pre-built CellRanger VDJ reference index",
    )
    parser.add_argument(
        "--gex-frna-probe-set",
        default=None,
        help="Probe set CSV for fixed RNA profiling (FFPE samples)",
    )
    parser.add_argument(
        "--gex-target-panel",
        default=None,
        help="Target panel CSV for targeted gene expression",
    )
    parser.add_argument(
        "--gex-cmo-set",
        default=None,
        help="Cell Multiplexing Oligo (CMO) reference CSV for multiplexed samples",
    )
    parser.add_argument(
        "--fb-reference",
        default=None,
        help="Feature barcoding reference CSV (e.g. antibody capture)",
    )
    parser.add_argument(
        "--vdj-inner-enrichment-primers",
        default=None,
        help="Text file with V(D)J cDNA enrichment primer sequences",
    )
    parser.add_argument(
        "--gex-barcode-sample-assignment",
        default=None,
        help="Barcode-to-sample assignment CSV to override CellRanger defaults",
    )
    parser.add_argument(
        "--cellranger-multi-barcodes",
        default=None,
        help="Additional samplesheet CSV with multiplexed sample information for cellrangermulti",
    )
    parser.add_argument(
        "--run-downstream",
        action="store_true",
        help="Opt in to running scrna_orchestrator after a canonical h5ad is detected",
    )
    parser.add_argument(
        "--skip-downstream",
        action="store_true",
        help="Force-skip the downstream handoff even if --run-downstream is given (handoff is already off by default)",
    )
    return parser


def _write_demo_samplesheet(samplesheet_path: Path) -> None:
    samplesheet_path.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(
        samplesheet_path, "sample,fastq_1,fastq_2,expected_cells,seq_center\n"
    )


def _write_error_result_if_safe(output_dir: Path, payload: dict[str, object]) -> None:
    """Persist structured errors unless doing so would overwrite a rejected output dir."""
    error_code = payload.get("error_code")
    if error_code in {"OUTPUT_DIR_NOT_EMPTY", "OUTPUT_DIR_NOT_WRITABLE"}:
        return
    if error_code == "INVALID_RESUME_STATE" and (output_dir / "result.json").exists():
        return
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_text_lf(output_dir / "result.json", json.dumps(payload, indent=2))
    except OSError:
        # stderr still receives the structured payload; avoid masking the root error.
        return


_BANNER = (
    "==============================================\n"
    f"  ClawBio :: {SKILL_NAME}\n"
    f"  nf-core/scrnaseq {NFCORE_SCRNASEQ_VERSION} orchestrator\n"
    "=============================================="
)


def _print(msg: str) -> None:
    """Log a human-readable status line to stdout (progress/traceability)."""
    print(msg, flush=True)


def _indent(text: str, n: int) -> str:
    pad = " " * n
    return "\n".join(pad + line for line in text.splitlines())


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not getattr(args, "no_banner", False):
        _print(_BANNER)
    output_dir = Path(args.output).expanduser().resolve()

    try:
        _apply_aligner_alias(args)
        _check_pipeline_version_supported(args)
        return _run_wrapper(args, output_dir)
    except SkillError as exc:
        return _handle_skill_error(output_dir, exc, verbose=getattr(args, "verbose", False))
    except KeyboardInterrupt:
        _print("[abort] Interrupted by user.")
        return 130
    except Exception as exc:
        return _handle_unexpected_error(
            output_dir, exc, verbose=getattr(args, "verbose", False)
        )


def _check_pipeline_version_supported(args: argparse.Namespace) -> None:
    """Keep execution on the pinned 4.1.0 contract unless explicitly overridden.

    The wrapper's parameter set, protocol matrix and output validations are
    hardcoded for nf-core/scrnaseq 4.1.0. Running a different remote tag/commit
    would apply 4.1.0 rules to a pipeline that may differ, so any non-4.1.0
    `--pipeline-version` is blocked by default; `--allow-pipeline-version-override`
    is a recorded, warned opt-in for advanced use (audit finding #1).
    """
    requested = str(getattr(args, "pipeline_version", "") or "").strip()
    if requested == NFCORE_SCRNASEQ_VERSION:
        return
    if getattr(args, "allow_pipeline_version_override", False):
        print(
            f"WARNING: --pipeline-version {requested!r} differs from the wrapper's pinned "
            f"nf-core/scrnaseq {NFCORE_SCRNASEQ_VERSION} contract. Parameter, protocol and "
            f"output validations remain {NFCORE_SCRNASEQ_VERSION}-specific and may not match.",
            file=sys.stderr,
        )
        return
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
        message=(
            f"--pipeline-version must be {NFCORE_SCRNASEQ_VERSION} "
            "(the version this wrapper's validations are pinned to)."
        ),
        fix=(
            f"Use --pipeline-version {NFCORE_SCRNASEQ_VERSION}, or pass "
            "--allow-pipeline-version-override to run a different version at your own risk "
            f"(validations stay {NFCORE_SCRNASEQ_VERSION})."
        ),
        details={"requested": requested, "contract_version": NFCORE_SCRNASEQ_VERSION},
    )


def _apply_aligner_alias(
    args: argparse.Namespace, parser: argparse.ArgumentParser | None = None
) -> None:
    """Normalize nf-core's --aligner spelling into the wrapper's --preset spelling."""
    aligner = getattr(args, "aligner", None)
    if not aligner:
        return
    aligner_to_preset = {
        aligner_name: preset_name
        for preset_name, aligner_name in PRESET_ALIGNERS.items()
    }
    requested_preset = aligner_to_preset[aligner]
    supplied_preset = getattr(args, "preset", DEFAULT_PRESET)
    if getattr(args, "preset_explicit", False) and supplied_preset != requested_preset:
        message = (
            f"--aligner {aligner!r} conflicts with --preset {supplied_preset!r}; "
            f"use --preset {requested_preset!r} or remove one of the two flags."
        )
        if parser is not None:
            parser.error(message)
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="Conflicting aligner and preset were supplied.",
            fix=message,
            details={
                "aligner": aligner,
                "preset": supplied_preset,
                "preset_from_aligner": requested_preset,
            },
        )
    args.preset = requested_preset


def _run_wrapper(args: argparse.Namespace, output_dir: Path) -> int:
    # Always validate/normalize the samplesheet into a temporary staging area and
    # commit it to the output bundle only on a successful execution. This keeps
    # --check and any failed preflight from leaving a half-written
    # reproducibility/ behind (audit H-8), and unifies the resume and non-resume
    # code paths.
    with tempfile.TemporaryDirectory(prefix="clawbio-scrnaseq-") as staging_dir:
        return _run_wrapper_with_staging(
            args, output_dir, staging_dir=Path(staging_dir)
        )


def _run_wrapper_with_staging(
    args: argparse.Namespace, output_dir: Path, *, staging_dir: Path | None
) -> int:
    check_output_dir_available(output_dir, resume=args.resume)
    pipeline_source = resolve_pipeline_source(
        requested_version=args.pipeline_version,
        allow_dirty=bool(getattr(args, "allow_dirty_pipeline", False)),
        require_local=bool(getattr(args, "require_local_pipeline", False)),
    )
    normalized_samplesheet, staged_samplesheet, samplesheet_summary = (
        _prepare_samplesheet(
            args,
            output_dir,
            staging_dir=staging_dir,
        )
    )
    preflight_result = run_preflight(
        args, pipeline_source=pipeline_source, samplesheet_summary=samplesheet_summary
    )
    _print(
        "[preflight] passed "
        f"(warnings: {len(preflight_result.get('warnings', []) or [])})"
    )
    if args.check:
        return _write_check_mode_result(
            output_dir,
            preflight_result=preflight_result,
            pipeline_source=pipeline_source,
        )
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
    # Only warn when the user *explicitly* asked for a non-star preset; the
    # bare default (preset_explicit is False) was never a real request, so
    # claiming "requested: 'standard'" would be misleading.
    if getattr(args, "preset_explicit", False) and args.preset != "star":
        print(
            f"WARNING: --demo forces preset=star (requested: {args.preset!r}). "
            "The nf-core test profile ships STAR-compatible data.",
            file=sys.stderr,
        )
    ignored_flags = _demo_ignored_flags(args)
    if ignored_flags:
        print(
            f"WARNING: --demo ignores pipeline flags ({', '.join(ignored_flags)}); "
            "the nf-core test profile owns every input, reference, protocol and "
            "QC/tuning setting, so these are not written to params.yaml. "
            "Drop --demo to run the pipeline on the supplied inputs.",
            file=sys.stderr,
        )
    args.preset = "star"
    args.skip_cellbender = True

    normalized_samplesheet = _final_samplesheet_path(output_dir, demo=True)
    staged_samplesheet = _staged_samplesheet_path(
        output_dir, staging_dir=staging_dir, demo=True
    )
    _write_demo_samplesheet(staged_samplesheet)
    return (
        normalized_samplesheet,
        staged_samplesheet,
        {
            "normalized_path": normalized_samplesheet,
            "sample_count": 0,
            "sample_names": [],
            "fastq_paths": [],
            "unknown_columns": [],
        },
    )


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
    _warn_about_omitted_unknown_columns(samplesheet_summary)
    return normalized_samplesheet, staged_samplesheet, samplesheet_summary


def _final_samplesheet_path(output_dir: Path, *, demo: bool = False) -> Path:
    filename = "samplesheet.demo.csv" if demo else "samplesheet.valid.csv"
    return output_dir / "reproducibility" / filename


def _staged_samplesheet_path(
    output_dir: Path, *, staging_dir: Path | None, demo: bool = False
) -> Path:
    if staging_dir is not None:
        filename = "samplesheet.demo.csv" if demo else "samplesheet.valid.csv"
        return staging_dir / filename
    return _final_samplesheet_path(output_dir, demo=demo)


def _warn_about_omitted_unknown_columns(
    samplesheet_summary: dict[str, object],
) -> None:
    unknown_columns = samplesheet_summary.get("unknown_columns", [])
    if unknown_columns:
        print(
            "WARNING: samplesheet contains columns outside the nf-core/scrnaseq "
            f"4.1.0 input schema; they are omitted from the normalized upstream sheet: {unknown_columns}",
            file=sys.stderr,
        )


def _write_check_mode_result(
    output_dir: Path,
    *,
    preflight_result: dict[str, Any],
    pipeline_source: dict[str, object],
) -> int:
    payload = {
        "ok": True,
        "skill": SKILL_NAME,
        "preflight": preflight_result,
        "pipeline_source": pipeline_source,
        "parameter_support": _build_parameter_support_summary(),
    }
    write_check_result(output_dir, payload)
    _print("[check] Preflight passed. No pipeline was launched (--check).")
    print(json.dumps(payload, indent=2))
    return 0


def _build_parameter_support_summary() -> dict[str, object]:
    return {
        "nfcore_version": NFCORE_SCRNASEQ_VERSION,
        # The wrapper deliberately exposes a curated, local-first subset of the
        # nf-core surface — not full passthrough (audit finding #7).
        "compatibility_mode": "clawbio_local_first_subset",
        "official_params_total": len(OFFICIAL_PARAMS),
        "supported_upstream_total": len(WRAPPER_SUPPORTED_UPSTREAM_PARAMS),
        "official_params": sorted(OFFICIAL_PARAMS),
        "supported_upstream": sorted(WRAPPER_SUPPORTED_UPSTREAM_PARAMS),
        "intentionally_unsupported": sorted(INTENTIONALLY_UNSUPPORTED_PARAMS),
        "intentionally_unsupported_reasons": dict(
            sorted(INTENTIONALLY_UNSUPPORTED_REASONS.items())
        ),
        "wrapper_policy_params": [
            "preset",
            "check",
            "run_downstream",
            "skip_downstream",
            "expected_cells",
            "timeout_hours",
            "work_dir",
            "allow_dirty_pipeline",
            "require_local_pipeline",
            "allow_conda_cellranger",
            "config",
        ],
    }


def _run_execution_mode(
    args: argparse.Namespace,
    *,
    output_dir: Path,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, Any],
    normalized_samplesheet: Path,
    staged_samplesheet: Path,
    samplesheet_summary: dict[str, object],
) -> int:
    params_payload = build_effective_params(
        args,
        normalized_samplesheet=normalized_samplesheet,
        output_dir=output_dir,
    )
    _check_resume_params_checksum(
        args, output_dir=output_dir, params_payload=params_payload
    )
    _commit_validated_samplesheet(staged_samplesheet, normalized_samplesheet)
    params_path = write_params_yaml(params_payload, output_dir=output_dir)
    command, command_str = _build_nextflow_invocation(
        args, output_dir, pipeline_source, params_path
    )
    nextflow_cwd = _nextflow_execution_cwd(output_dir)
    _print(
        f"[execute] launching Nextflow "
        f"({pipeline_source['source_kind']} → {pipeline_source['resolved_version']})"
    )
    execution_result = execute_nextflow(
        command,
        cwd=nextflow_cwd,
        output_dir=output_dir,
        timeout_seconds=_resolve_timeout_seconds(args),
    )
    _print("[execute] completed")
    _print("[outputs] parsing pipeline outputs")
    parsed_outputs = _parse_outputs_with_effective_aligner(output_dir, args)
    _raise_if_expected_outputs_missing(parsed_outputs, output_dir=output_dir)
    _print("[report] writing report.md, result.json, commands.sh")
    _print("[provenance] writing reproducibility bundle")
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
    _print(f"[done] Wrapper completed successfully. Output: {output_dir}")
    return 0


def _commit_validated_samplesheet(
    staged_samplesheet: Path, normalized_samplesheet: Path
) -> None:
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
    return hashlib.sha256(
        serialize_params_yaml(params_payload).encode("utf-8")
    ).hexdigest()


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
        work_dir=_resolve_nextflow_work_dir(args, output_dir),
        extra_configs=_build_extra_nextflow_configs(args, output_dir),
    )


def _resolve_nextflow_work_dir(args: argparse.Namespace, output_dir: Path) -> Path | str:
    raw_work_dir = getattr(args, "work_dir", None)
    if not raw_work_dir:
        return output_dir / "upstream" / "work"
    work_dir = str(raw_work_dir).strip()
    if "://" in work_dir:
        return work_dir
    return Path(work_dir).expanduser().resolve()


def _build_extra_nextflow_configs(
    args: argparse.Namespace, output_dir: Path
) -> list[Path]:
    configs = [
        Path(value).expanduser().resolve()
        for value in getattr(args, "extra_config", []) or []
    ]
    # The VirtioFS EDEADLK fix is macOS-specific:
    #   Linux native: no VirtioFS; Docker bind-mounts are direct overlayfs.
    #   Windows/WSL2: Nextflow runs inside WSL2 (reports "Linux"); work dir is on ext4.
    #   Windows native: not officially supported by Nextflow; users should use WSL2.
    if profile_includes(args.profile, "docker"):
        if platform.system() == "Darwin":
            # demo gates the test-profile resource ceilings: real runs get only the
            # always-safe VirtioFS / Apple-Silicon / STAR-FIFO workarounds (audit H-1).
            configs.append(
                write_macos_docker_config(
                    output_dir, demo=bool(getattr(args, "demo", False))
                )
            )
        elif not bool(getattr(args, "demo", False)):
            # Non-macOS docker uses the local executor, which aborts a real run when a
            # process's production requirement exceeds this host's RAM (e.g. STARsolo
            # requesting 72 GB on a 63 GB machine). Cap requests to the host per
            # nf-core resourceLimits guidance. --demo relies on -profile test's own
            # (smaller) resourceLimits, so it is deliberately left untouched. Applied
            # to the live run only — never baked into the portable commands.sh bundle,
            # which is rebuilt from params.yaml (see repro_commands.py).
            configs.append(_write_resource_limits_config(output_dir, args=args))
    return configs


_RESOURCE_LIMITS_MEMORY_HEADROOM_GB = 4
_RESOURCE_LIMITS_MEMORY_FLOOR_GB = 8


def _host_memory_gb() -> int:
    """Best-effort total physical RAM in whole GB (portable; falls back to 16)."""
    try:
        return max(1, int(os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / (1024 ** 3)))
    except (ValueError, OSError, AttributeError):
        return 16


def _docker_vm_memory_gb() -> int | None:
    """Docker runtime's total memory in whole GB, or None when unknown.

    On native Linux this reports the host's cgroup memory (≈ physical RAM). Best-effort:
    any failure (docker absent, daemon down, unparseable output) returns None.
    """
    if not shutil.which("docker"):
        return None
    try:
        proc = subprocess.run(
            ["docker", "info", "--format", "{{.MemTotal}}"],
            capture_output=True,
            text=True,
            timeout=10,
        )
    except (OSError, subprocess.SubprocessError):
        return None
    if proc.returncode != 0:
        return None
    raw = proc.stdout.strip()
    if not raw.isdigit():
        return None
    gb = int(int(raw) / (1024 ** 3))
    return gb or None


def _resource_limits_memory_gb() -> int:
    """Per-process memory ceiling (GB) for the non-macOS resourceLimits config.

    nf-core's guidance is that ``resourceLimits`` should match the machine's maximum
    (https://nf-co.re/docs/running/configuration/nextflow-for-your-system). Native
    Linux Docker exposes the full host RAM, so there is no 15 GB macOS-VM ceiling
    here. We take the smaller of the detected physical RAM and (when known) the Docker
    runtime's MemTotal, then reserve a few GB of headroom for the OS, the Nextflow JVM
    and the Docker daemon.
    """
    signals = [gb for gb in (_host_memory_gb(), _docker_vm_memory_gb()) if gb]
    budget = (min(signals) - _RESOURCE_LIMITS_MEMORY_HEADROOM_GB) if signals else 12
    return max(_RESOURCE_LIMITS_MEMORY_FLOOR_GB, budget)


def _write_resource_limits_config(output_dir: Path, *, args: argparse.Namespace) -> Path:
    """Write a host-scaled ``resourceLimits`` config for non-macOS docker runs.

    Carries only the portable resourceLimits block — none of the macOS-only
    workarounds (``--platform``, ``stageInMode = 'copy'``, STAR FIFO overrides).

    Written inside ``reproducibility/`` (not the output root) so it ships with the
    portable bundle: ``commands.sh`` re-applies it on non-macOS replay hosts, letting a
    from-scratch reproduction on the generating machine succeed instead of re-aborting
    with "Process requirement exceeds available memory".
    """
    config_path = output_dir / "reproducibility" / "resource_limits.config"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    cpus = max(1, os.cpu_count() or 1)
    memory_gb = _resource_limits_memory_gb()
    timeout_hours = getattr(args, "timeout_hours", _DEFAULT_TIMEOUT_HOURS)
    time_hours = int(_DEFAULT_TIMEOUT_HOURS) if timeout_hours == 0 else max(1, int(timeout_hours))
    write_text_lf(
        config_path,
        "// Cap process resource requests to this host so Nextflow's local executor\n"
        "// does not abort real runs whose production requirements exceed this machine.\n"
        "// Follows nf-core resourceLimits guidance (values match the machine maximum):\n"
        "// https://nf-co.re/docs/running/configuration/nextflow-for-your-system\n"
        "process {\n"
        "    resourceLimits = [\n"
        f"        cpus: {cpus},\n"
        f"        memory: '{memory_gb}.GB',\n"
        f"        time: '{time_hours}.h'\n"
        "    ]\n"
        "}\n",
    )
    return config_path


def _nextflow_execution_cwd(output_dir: Path) -> Path:
    # params.input is written relative to output_dir so it stays schema-safe even
    # when the absolute workspace path contains spaces.
    return output_dir


def _nextflow_replay_command(command_str: str, cwd: Path) -> str:
    return f"cd {shlex.quote(cwd.as_posix())} && {command_str}"


def _parse_outputs_with_effective_aligner(
    output_dir: Path, args: argparse.Namespace
) -> dict[str, object]:
    parsed_outputs = parse_outputs(output_dir)
    aligner_effective = PRESET_ALIGNERS[args.preset]
    demo = bool(getattr(args, "demo", False))
    return {
        **parsed_outputs,
        "aligner_effective": aligner_effective,
        "output_validation": validate_expected_outputs(
            parsed_outputs,
            aligner=aligner_effective,
            # A hermetic demo never writes skip flags (audit H-01), so the test
            # profile always runs FastQC + MultiQC — validate them as required
            # regardless of any --skip-* the user passed alongside --demo.
            skip_multiqc=False if demo else bool(getattr(args, "skip_multiqc", False)),
            skip_fastqc=False if demo else bool(getattr(args, "skip_fastqc", False)),
        ),
    }


def _raise_if_expected_outputs_missing(
    parsed_outputs: dict[str, object], *, output_dir: Path
) -> None:
    output_validation = parsed_outputs.get("output_validation")
    if not isinstance(output_validation, dict):
        return
    missing_required = output_validation.get("missing_required", [])
    if not missing_required:
        return
    raise SkillError(
        stage="parsing",
        error_code=ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND,
        message="nf-core/scrnaseq completed but required output files or directories were not found.",
        fix=(
            "Inspect Nextflow logs in logs/stdout.txt, logs/stderr.txt, and upstream/work. "
            "If MultiQC was intentionally disabled, rerun the wrapper with --skip-multiqc."
        ),
        details={
            "output_dir": str(output_dir),
            "aligner_effective": parsed_outputs.get("aligner_effective", ""),
            "missing_required": missing_required,
            "missing_optional": output_validation.get("missing_optional", []),
        },
    )


def _write_success_outputs(
    output_dir: Path,
    *,
    args: argparse.Namespace,
    pipeline_source: dict[str, object],
    preflight_result: dict[str, Any],
    params_path: Path,
    params_payload: dict[str, object],
    normalized_samplesheet: Path,
    samplesheet_summary: dict[str, object],
    parsed_outputs: dict[str, object],
    execution_result: dict[str, object],
    command_str: str,
) -> None:
    write_repro_commands(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        nextflow_version=preflight_result.get("nextflow", {}).get("version"),
    )
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


_SCRNA_ORCHESTRATOR_DEFAULT = (
    _SKILL_DIR.parent / "scrna-orchestrator" / "scrna_orchestrator.py"
)


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
    if not getattr(args, "run_downstream", False) or getattr(
        args, "skip_downstream", False
    ):
        return
    preferred_h5ad = str(parsed_outputs.get("preferred_h5ad", "")).strip()
    if not preferred_h5ad:
        h5ad_candidates = parsed_outputs.get("h5ad_candidates", [])
        if h5ad_candidates:
            print(
                "WARNING: downstream handoff was requested, but automatic h5ad selection is ambiguous. "
                "Inspect result.json h5ad_candidates and run downstream analysis manually with the chosen file.",
                file=sys.stderr,
            )
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
        "--input",
        preferred_h5ad,
        "--output",
        str(downstream_output),
    ]
    # Audit trail: the orchestrator script is resolvable via the
    # CLAWBIO_SCRNA_ORCHESTRATOR env var, so record exactly which script (and its
    # checksum) was executed for this handoff (audit Hallazgo 9).
    record: dict[str, object] = {
        "orchestrator": orchestrator.as_posix(),
        "orchestrator_sha256": _sha256_file_safe(orchestrator),
        "orchestrator_source": "env:CLAWBIO_SCRNA_ORCHESTRATOR"
        if os.environ.get("CLAWBIO_SCRNA_ORCHESTRATOR")
        else "default_sibling",
        "preferred_h5ad": preferred_h5ad,
        "downstream_output": str(downstream_output),
    }
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
        record["status"] = "timeout"
        _write_handoff_provenance(output_dir, record)
        print(
            "WARNING: downstream scrna_orchestrator timed out. "
            "The nf-core pipeline completed successfully; run downstream analysis manually if needed.",
            file=sys.stderr,
        )
        return
    record["status"] = "ok" if result.returncode == 0 else "failed"
    record["returncode"] = result.returncode
    _write_handoff_provenance(output_dir, record)
    if result.returncode != 0:
        print(
            f"WARNING: downstream scrna_orchestrator exited with code {result.returncode}. "
            "The nf-core pipeline completed successfully; inspect scrna_analysis/ manually.",
            file=sys.stderr,
        )


def _sha256_file_safe(path: Path) -> str:
    """Best-effort SHA-256 of a script; empty string if it cannot be read."""
    try:
        return hashlib.sha256(path.read_bytes()).hexdigest()
    except OSError:
        return ""


def _write_handoff_provenance(output_dir: Path, record: dict[str, object]) -> None:
    """Persist which downstream orchestrator was executed (path + checksum + outcome)."""
    provenance_dir = output_dir / "provenance"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    write_text_lf(provenance_dir / "handoff.json", json.dumps(record, indent=2))


def _handle_skill_error(
    output_dir: Path, exc: SkillError, *, verbose: bool = False
) -> int:
    payload = exc.to_dict()
    _write_error_result_if_safe(output_dir, payload)
    # Human-readable box on stdout (sarek-style traceability).
    _print("")
    _print("================ SkillError ================")
    _print(f"  stage:   {exc.stage}")
    _print(f"  code:    {exc.error_code}")
    _print(f"  message: {exc.message}")
    _print(f"  fix:     {exc.fix}")
    if exc.details and verbose:
        _print("  details:")
        _print(_indent(json.dumps(exc.details, indent=2, default=str), 4))
    _print("============================================")
    # Machine-readable error on stderr: always available even when result.json
    # cannot be written (e.g. the output dir is a file or inside the repo).
    print(json.dumps(payload, indent=2), file=sys.stderr)
    return 1


def _handle_unexpected_error(
    output_dir: Path, exc: Exception, *, verbose: bool = False
) -> int:
    payload = {
        "status": "error",
        "ok": False,
        "stage": "internal",
        "error_code": ErrorCode.UNEXPECTED_ERROR,
        "message": str(exc),
        "fix": "Report this as a bug. Include the full traceback and your command arguments.",
        "details": {
            "exception_type": type(exc).__name__,
            "traceback": traceback.format_exc(),
        },
    }
    _write_error_result_if_safe(output_dir, payload)
    _print("")
    _print("================ Internal error ================")
    _print(f"  type:    {type(exc).__name__}")
    _print(f"  message: {exc}")
    if verbose:
        _print(_indent(traceback.format_exc(), 4))
    else:
        _print("  Re-run with --verbose for the full traceback.")
    _print("================================================")
    # Machine-readable error on stderr (always available; see _handle_skill_error).
    print(json.dumps(payload, indent=2), file=sys.stderr)
    return 1


if __name__ == "__main__":
    raise SystemExit(main())
