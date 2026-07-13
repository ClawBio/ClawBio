#!/usr/bin/env python3
from __future__ import annotations

import argparse
import json
import os
import platform
import shlex
import shutil
import subprocess
import sys
import tempfile
import time
import traceback
from pathlib import Path


_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

_DOWNSTREAM_HANDOFF_TIMEOUT_SECONDS = 60 * 60 * 2


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules


purge_foreign_bare_modules(
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

from command_builder import build_nextflow_command
from errors import ErrorCode, SkillError
from executor import execute_nextflow
from outputs_parser import parse_outputs
from params_builder import build_effective_params, write_params_yaml
from pipeline_source import resolve_pipeline_source
from preflight import check_resume_params_checksum, check_resume_samplesheet_checksum, run_preflight
from provenance import write_provenance_bundle
from reporting import write_check_result, write_report, write_repro_commands, write_result
from samplesheet_builder import validate_and_normalize_samplesheet
from schemas import (
    DEFAULT_PIPELINE_VERSION,
    DEFAULT_PROFILE,
    DEFAULT_TIMEOUT_HOURS,
    DEFAULT_TRIMMER,
    SKILL_NAME,
    SUPPORTED_ALIGNERS,
    SUPPORTED_BRACKEN_PRECISION,
    SUPPORTED_CONTAMINANT_SCREENING,
    SUPPORTED_PSEUDO_ALIGNERS,
    SUPPORTED_PUBLISH_DIR_MODES,
    SUPPORTED_RIBO_TOOLS,
    SUPPORTED_SALMON_LIB_TYPES,
    SUPPORTED_TRIMMERS,
    SUPPORTED_UMI_EXTRACT_METHODS,
    SUPPORTED_UMI_GROUPING_METHODS,
    SUPPORTED_UMI_TOOLS,
)

from clawbio.common.textio import write_text_lf


_DASH_VALUE_OPTIONS = {
    "--extra-trimgalore-args",
    "--extra-fastp-args",
    "--extra-fqlint-args",
    "--extra-star-align-args",
    "--extra-bowtie2-align-args",
    "--extra-salmon-quant-args",
    "--extra-kallisto-quant-args",
    "--rsem-extra-args",
    # GPU container flags (e.g. '--gpus all', '--nv') begin with a dash and would
    # otherwise be misparsed as wrapper options without the dash-value rewrite.
    "--gpu-container-options",
}


class _RnaseqArgumentParser(argparse.ArgumentParser):
    """Allow upstream extra-args values that begin with dashes."""

    def parse_args(self, args: list[str] | None = None, namespace: argparse.Namespace | None = None):
        return super().parse_args(self._normalise_dash_values(args), namespace)

    def parse_known_args(self, args: list[str] | None = None, namespace: argparse.Namespace | None = None):
        return super().parse_known_args(self._normalise_dash_values(args), namespace)

    def _normalise_dash_values(self, args: list[str] | None) -> list[str] | None:
        if args is None:
            return None
        known_options = {option for action in self._actions for option in action.option_strings}
        normalised: list[str] = []
        i = 0
        while i < len(args):
            token = args[i]
            if token in _DASH_VALUE_OPTIONS and i + 1 < len(args):
                value = args[i + 1]
                if value.startswith("-") and value not in known_options:
                    normalised.append(f"{token}={value}")
                    i += 2
                    continue
            normalised.append(token)
            i += 1
        return normalised


def build_parser() -> argparse.ArgumentParser:
    parser = _RnaseqArgumentParser(description="Run nf-core/rnaseq through the ClawBio wrapper.")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    parser.add_argument("--no-banner", action="store_true", help="Suppress the startup banner")
    parser.add_argument("--input", help="Path to a valid nf-core/rnaseq samplesheet.csv")
    parser.add_argument("--output", required=True, help="Output directory for wrapper results")
    parser.add_argument("--demo", action="store_true", help="Run the upstream test profile without user FASTQs")
    parser.add_argument("--check", action="store_true", help="Run wrapper preflight only and exit before Nextflow")
    parser.add_argument(
        "--profile",
        default=DEFAULT_PROFILE,
        help=(
            "Nextflow execution backend profile. "
            "Common values: docker, singularity, apptainer, podman, conda, mamba, shifter, charliecloud. "
            "Modifier profiles (appended automatically by --prokaryotic / --arm / --rapid-quant): "
            "prokaryotic, arm64, rapid_quant. "
            "nf-core/rnaseq test profiles (test, test_full, test_full_aws) and institutional "
            "profiles (e.g. uppmax, slurm, my_hpc) are passed through unchanged to Nextflow."
        ),
    )
    parser.add_argument("--pipeline-version", default=DEFAULT_PIPELINE_VERSION, help="Pinned remote nf-core/rnaseq tag or commit")
    parser.add_argument(
        "--allow-pipeline-version-override",
        action="store_true",
        help=(
            f"Allow a --pipeline-version other than the pinned nf-core/rnaseq {DEFAULT_PIPELINE_VERSION} "
            "contract. The wrapper's parameter/enum/output validations remain "
            f"{DEFAULT_PIPELINE_VERSION}-specific, so this is at your own risk (warned, recorded)."
        ),
    )
    parser.add_argument("--pipeline-local", default=None, help="Optional local nf-core/rnaseq checkout override")
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
    parser.add_argument("--resume", action="store_true", help="Attempt checksum-validated Nextflow resume")
    parser.add_argument(
        "--timeout-hours",
        type=_timeout_hours,
        default=DEFAULT_TIMEOUT_HOURS,
        help=(
            f"Wall-clock ceiling for the local Nextflow run before it is terminated "
            f"(default: {DEFAULT_TIMEOUT_HOURS}). Raise this for large cohorts so a long "
            "but healthy run is not killed. Use 0 to disable the cap for long HPC/cloud "
            "runs whose walltime is enforced by the scheduler."
        ),
    )

    parser.add_argument("--aligner", default=None, choices=sorted(SUPPORTED_ALIGNERS))
    parser.add_argument("--pseudo-aligner", default=None, choices=sorted(SUPPORTED_PSEUDO_ALIGNERS))
    parser.add_argument(
        "--pseudo-aligner-kmer-size",
        type=int,
        default=None,
        help=(
            "Index k-mer length for the Salmon/Kallisto pseudo-aligner. Must be an odd "
            "integer in 1..31 (their 64-bit k-mer hard cap; pipeline default: 31). "
            "Lower it for short reads (<50 bp)."
        ),
    )
    parser.add_argument("--prokaryotic", action="store_true", help="Compose the prokaryotic nf-core profile")

    parser.add_argument("--genome", default=None)
    parser.add_argument("--igenomes-base", default=None, help="Override iGenomes base path (S3 URI or local directory). Required when using --genome on HPC clusters with a local iGenomes mirror.")
    parser.add_argument("--fasta", default=None)
    parser.add_argument("--gtf", default=None)
    parser.add_argument("--gff", default=None)
    parser.add_argument("--gene-bed", default=None)
    parser.add_argument("--transcript-fasta", default=None)
    parser.add_argument("--additional-fasta", default=None)
    parser.add_argument("--splicesites", default=None)
    parser.add_argument("--star-index", default=None)
    parser.add_argument("--hisat2-index", default=None)
    parser.add_argument("--rsem-index", default=None)
    parser.add_argument("--salmon-index", default=None)
    parser.add_argument("--kallisto-index", default=None)
    parser.add_argument("--bowtie2-index", default=None)
    parser.add_argument("--gencode", action="store_true")

    parser.add_argument("--trimmer", default=DEFAULT_TRIMMER, choices=sorted(SUPPORTED_TRIMMERS))
    parser.add_argument("--extra-trimgalore-args", default=None)
    parser.add_argument("--extra-fastp-args", default=None)
    parser.add_argument("--min-trimmed-reads", type=int, default=None, help="Minimum number of reads a sample must retain after trimming to continue (pipeline default: 10000; must be ≥ 0).")
    parser.add_argument("--sortmerna-index", default=None)
    parser.add_argument("--ribo-database-manifest", default=None)
    parser.add_argument("--remove-ribo-rna", action="store_true")
    parser.add_argument("--ribo-removal-tool", default=None, choices=sorted(SUPPORTED_RIBO_TOOLS))
    parser.add_argument("--with-umi", action="store_true")
    parser.add_argument("--umi-dedup-tool", default=None, choices=sorted(SUPPORTED_UMI_TOOLS))
    parser.add_argument("--umitools-extract-method", default=None, choices=sorted(SUPPORTED_UMI_EXTRACT_METHODS))
    parser.add_argument("--umitools-bc-pattern", default=None)
    parser.add_argument("--umitools-bc-pattern2", default=None)
    parser.add_argument("--umitools-umi-separator", default=None)
    parser.add_argument("--umi-discard-read", type=int, default=None, choices=[1, 2])
    parser.add_argument(
        "--umitools-grouping-method",
        default=None,
        choices=sorted(SUPPORTED_UMI_GROUPING_METHODS),
    )
    parser.add_argument("--skip-umi-extract", action="store_true")
    parser.add_argument("--umitools-dedup-stats", action="store_true")
    parser.add_argument("--umitools-dedup-primary-only", action="store_true")

    parser.add_argument("--seq-center", default=None)
    parser.add_argument("--seq-platform", default=None)
    parser.add_argument("--min-mapped-reads", type=float, default=None)
    parser.add_argument("--star-ignore-sjdbgtf", action="store_true")
    parser.add_argument(
        "--salmon-quant-libtype",
        default=None,
        # Enum values from nf-core/rnaseq 3.26.0 nextflow_schema.json (single source: schemas.py).
        choices=sorted(SUPPORTED_SALMON_LIB_TYPES),
    )
    parser.add_argument("--extra-star-align-args", default=None)
    parser.add_argument("--extra-bowtie2-align-args", default=None)
    parser.add_argument("--extra-salmon-quant-args", default=None)
    parser.add_argument("--extra-kallisto-quant-args", default=None)
    parser.add_argument(
        "--kallisto-quant-fraglen",
        type=int,
        default=None,
        help="Mean fragment length for single-end Kallisto quantification (single-end only; pipeline default: 200).",
    )
    parser.add_argument(
        "--kallisto-quant-fraglen-sd",
        type=int,
        default=None,
        help="Fragment length standard deviation for single-end Kallisto quantification (single-end only; pipeline default: 200).",
    )
    parser.add_argument("--bam-csi-index", action="store_true")
    parser.add_argument("--stringtie-ignore-gtf", action="store_true")
    parser.add_argument("--stranded-threshold", type=float, default=None)
    parser.add_argument("--unstranded-threshold", type=float, default=None)

    parser.add_argument("--contaminant-screening", default=None, choices=sorted(SUPPORTED_CONTAMINANT_SCREENING))
    parser.add_argument("--contaminant-screening-input", default=None, choices=["trimmed", "unmapped"])
    parser.add_argument("--kraken-db", default=None)
    parser.add_argument(
        "--bracken-precision",
        default=None,
        choices=list(SUPPORTED_BRACKEN_PRECISION),
        help="Bracken taxonomic level for classification (D=domain P=phylum C=class O=order F=family G=genus S=species). Pipeline default: S.",
    )
    parser.add_argument("--sylph-db", default=None)
    parser.add_argument("--sylph-taxonomy", default=None)
    parser.add_argument("--bbsplit-fasta-list", default=None)
    parser.add_argument("--bbsplit-index", default=None)
    parser.add_argument("--skip-bbsplit", action="store_true", help="Skip BBSplit even when --bbsplit-fasta-list or --bbsplit-index is provided")
    parser.add_argument("--save-kraken-assignments", action="store_true")
    parser.add_argument("--save-kraken-unassigned", action="store_true")
    parser.add_argument("--save-bbsplit-reads", action="store_true")

    parser.add_argument("--hisat2-build-memory", default=None)
    parser.add_argument("--gffread-transcript-fasta", action="store_true")
    parser.add_argument("--arm", action="store_true", help="Compose arm64 into the Nextflow profile string (e.g. --profile docker --arm → -profile docker,arm64). Use on Apple M-series or AWS Graviton hosts.")
    parser.add_argument("--use-rustqc", action="store_true")
    parser.add_argument("--use-parabricks-star", action="store_true")
    parser.add_argument("--use-sentieon-star", action="store_true")
    parser.add_argument("--use-gpu-ribodetector", action="store_true")
    parser.add_argument("--gpu-container-options", default=None, help="GPU container flag override (e.g. '--gpus all')")
    parser.add_argument("--rapid-quant", action="store_true")
    parser.add_argument("--featurecounts-group-type", default=None)
    parser.add_argument("--featurecounts-feature-type", default=None)
    parser.add_argument("--gtf-extra-attributes", default=None)
    parser.add_argument("--gtf-group-features", default=None)
    parser.add_argument("--deseq2-vst", dest="deseq2_vst", action="store_true", default=None)
    parser.add_argument("--no-deseq2-vst", dest="deseq2_vst", action="store_false")
    parser.add_argument(
        "--rseqc-modules",
        default=None,
        help=(
            "Comma-separated RSeQC modules to run. "
            "Default set: bam_stat,inner_distance,infer_experiment,junction_annotation,"
            "junction_saturation,read_distribution,read_duplication. "
            "Add 'tin' for Transcript Integrity Number analysis (slow on large datasets)."
        ),
    )

    parser.add_argument("--extra-fqlint-args", default=None, help="Additional arguments for fq lint")
    parser.add_argument("--skip-linting", action="store_true")
    parser.add_argument("--skip-gtf-filter", action="store_true")
    parser.add_argument("--skip-gtf-transcript-filter", action="store_true")
    parser.add_argument("--skip-trimming", action="store_true")
    parser.add_argument("--skip-alignment", action="store_true")
    parser.add_argument("--skip-pseudo-alignment", action="store_true")
    parser.add_argument("--skip-quantification-merge", action="store_true")
    parser.add_argument("--skip-markduplicates", action="store_true")
    parser.add_argument("--skip-bigwig", action="store_true")
    parser.add_argument("--skip-stringtie", action="store_true")
    parser.add_argument("--skip-fastqc", action="store_true")
    parser.add_argument("--skip-dupradar", action="store_true")
    parser.add_argument("--skip-qualimap", action="store_true")
    parser.add_argument("--skip-rseqc", action="store_true")
    parser.add_argument("--skip-biotype-qc", action="store_true")
    parser.add_argument("--skip-deseq2-qc", action="store_true")
    parser.add_argument("--skip-multiqc", action="store_true")
    parser.add_argument("--skip-qc", action="store_true")
    parser.add_argument("--enable-preseq", action="store_true")

    parser.add_argument("--save-reference", action="store_true")
    parser.add_argument("--save-trimmed", action="store_true")
    parser.add_argument("--save-align-intermeds", action="store_true")
    parser.add_argument("--save-unaligned", action="store_true")
    parser.add_argument("--save-merged-fastq", action="store_true")
    parser.add_argument("--save-non-ribo-reads", action="store_true")
    parser.add_argument("--save-umi-intermeds", action="store_true")

    parser.add_argument(
        "--publish-dir-mode",
        default=None,
        choices=sorted(SUPPORTED_PUBLISH_DIR_MODES),
        help="Nextflow publishDir mode. Use 'symlink' or 'link' to save disk space on HPC (pipeline default: copy).",
    )
    parser.add_argument("--email", default=None)
    parser.add_argument("--email-on-fail", default=None, help="Send failure notification to this address even when --email is not set.")
    parser.add_argument("--multiqc-title", default=None)
    parser.add_argument("--multiqc-config", default=None)
    parser.add_argument("--multiqc-logo", default=None)
    parser.add_argument("--multiqc-methods-description", default=None)
    parser.add_argument(
        "--rsem-extra-args",
        default=None,
        help=(
            "[No effect] nf-core/rnaseq 3.26.0 does not expose an extra_rsem_quant_args parameter. "
            "A preflight warning is emitted. To customise RSEM args use a Nextflow config: "
            "process { withName: 'RSEM_CALCULATEEXPRESSION' { ext.args = '...' } }"
        ),
    )
    parser.add_argument(
        "--nextflow-config",
        "-c",
        "--config",
        dest="nextflow_config",
        action="append",
        metavar="CONFIG",
        default=None,
        help=(
            "Additional Nextflow config file(s) passed as -c to Nextflow (aliases: -c/--config). "
            "Can be repeated: --nextflow-config hpc.config --nextflow-config rsem.config. "
            "Use this for institution-specific settings, RSEM ext.args, or any non-parametric "
            "configuration not exposed via --params-file."
        ),
    )
    parser.add_argument("--run-downstream", action="store_true", help="Write rnaseq-de handoff template when counts are available")
    parser.add_argument("--skip-downstream", action="store_true", help="Compatibility flag; prevents handoff template creation")
    parser.add_argument("--metadata", default=None, help="Metadata CSV/TSV passed to rnaseq-de when auto-handoff is enabled")
    parser.add_argument("--formula", default=None, help="DESeq2 formula passed to rnaseq-de when auto-handoff is enabled")
    parser.add_argument("--contrast", default=None, help="rnaseq-de contrast factor,numerator,denominator")
    parser.add_argument("--downstream-output", default=None, help="Optional rnaseq-de output directory")
    return parser


def _check_pipeline_version_supported(args: argparse.Namespace) -> None:
    """Keep execution on the pinned nf-core/rnaseq contract unless overridden.

    The wrapper's parameter set, enum validations and output checks are written
    for nf-core/rnaseq ``DEFAULT_PIPELINE_VERSION`` (3.26.0). Running a different
    remote tag/commit would apply 3.26.0 rules to a pipeline that may differ, so
    any non-pinned ``--pipeline-version`` is blocked by default;
    ``--allow-pipeline-version-override`` is a recorded, warned opt-in for advanced
    use (mirrors the scrnaseq wrapper's gate; audit OBS-2).
    """
    requested = str(getattr(args, "pipeline_version", "") or "").strip()
    if requested == DEFAULT_PIPELINE_VERSION:
        return
    if getattr(args, "allow_pipeline_version_override", False):
        print(
            f"WARNING: --pipeline-version {requested!r} differs from the wrapper's pinned "
            f"nf-core/rnaseq {DEFAULT_PIPELINE_VERSION} contract. Parameter, enum and output "
            f"validations remain {DEFAULT_PIPELINE_VERSION}-specific and may not match.",
            file=sys.stderr,
        )
        return
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
        message=(
            f"--pipeline-version must be {DEFAULT_PIPELINE_VERSION} "
            "(the version this wrapper's validations are pinned to)."
        ),
        fix=(
            f"Use --pipeline-version {DEFAULT_PIPELINE_VERSION}, or pass "
            "--allow-pipeline-version-override to run a different version at your own risk "
            f"(validations stay {DEFAULT_PIPELINE_VERSION})."
        ),
        details={"requested": requested, "contract_version": DEFAULT_PIPELINE_VERSION},
    )


def _check_pipeline_source_version(args: argparse.Namespace, pipeline_source: dict[str, object]) -> None:
    """Enforce the pinned-version contract on a *resolved* pipeline source.

    ``_check_pipeline_version_supported`` only inspects ``--pipeline-version`` (the
    remote tag). A sibling ``../rnaseq`` checkout is auto-detected and run regardless
    of which branch/commit it sits on, so its ``manifest.version`` must be checked
    too — otherwise 3.26.0-specific parameter/enum/output validations would silently
    apply to a different pipeline version (audit F-01). An unparseable/absent manifest
    version (``""``) is treated as unknown and warned, never as a hard mismatch.
    """
    if pipeline_source.get("source_kind") != "local_checkout":
        return
    manifest_version = str(pipeline_source.get("manifest_version", "") or "").strip()
    if manifest_version == DEFAULT_PIPELINE_VERSION:
        return
    allow_override = getattr(args, "allow_pipeline_version_override", False)
    if not manifest_version:
        print(
            "WARNING: could not determine the local nf-core/rnaseq checkout version "
            f"(expected {DEFAULT_PIPELINE_VERSION}). The wrapper's parameter, enum and "
            f"output validations remain {DEFAULT_PIPELINE_VERSION}-specific.",
            file=sys.stderr,
        )
        return
    if allow_override:
        print(
            f"WARNING: local nf-core/rnaseq checkout reports version {manifest_version!r}, "
            f"which differs from the wrapper's pinned {DEFAULT_PIPELINE_VERSION} contract. "
            f"Validations remain {DEFAULT_PIPELINE_VERSION}-specific and may not match.",
            file=sys.stderr,
        )
        return
    raise SkillError(
        stage="validation",
        error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
        message=(
            f"Local nf-core/rnaseq checkout is version {manifest_version!r}, but this wrapper's "
            f"validations are pinned to {DEFAULT_PIPELINE_VERSION}."
        ),
        fix=(
            f"Check out nf-core/rnaseq {DEFAULT_PIPELINE_VERSION} in the sibling 'rnaseq' "
            "directory, point --pipeline-local at a 3.26.0 checkout, remove the local "
            "checkout to use the pinned remote pipeline, or pass "
            "--allow-pipeline-version-override to run it at your own risk "
            f"(validations stay {DEFAULT_PIPELINE_VERSION})."
        ),
        details={
            "manifest_version": manifest_version,
            "contract_version": DEFAULT_PIPELINE_VERSION,
            "source_ref": pipeline_source.get("source_ref", ""),
        },
    )


def _raise_if_expected_outputs_missing(
    parsed_outputs: dict[str, object], *, args: argparse.Namespace, output_dir: Path
) -> None:
    """Fail a completed run that produced none of its required outputs.

    nf-core/rnaseq always writes ``pipeline_info/``, and a standard quantifying run
    produces a merged gene-count matrix. Treating their absence as success would
    hide a broken run, so raise EXPECTED_OUTPUTS_NOT_FOUND instead (audit OBS-1;
    mirrors the scrnaseq wrapper's h5ad gate). The count-matrix check is
    conditional — it is NOT applied when the merge was skipped
    (``--skip-quantification-merge`` / rapid_quant) or when no quantifier ran
    (hisat2 without a pseudo-aligner), which legitimately produce no merged matrix.
    """
    missing: list[str] = []
    if not parsed_outputs.get("pipeline_info_dir"):
        missing.append("pipeline_info")

    aligner = str(parsed_outputs.get("aligner_effective", ""))
    pseudo = str(parsed_outputs.get("pseudo_aligner_effective", ""))
    skip_merge = bool(parsed_outputs.get("skip_quantification_merge"))
    skip_alignment = bool(getattr(args, "skip_alignment", False))
    skip_pseudo = bool(getattr(args, "skip_pseudo_alignment", False))
    aligner_quantifies = (
        aligner in {"star_salmon", "star_rsem", "bowtie2_salmon"} and not skip_alignment
    )
    pseudo_quantifies = bool(pseudo) and not skip_pseudo
    counts_expected = (not skip_merge) and (aligner_quantifies or pseudo_quantifies)
    has_counts = bool(
        parsed_outputs.get("preferred_counts_tsv")
        or parsed_outputs.get("raw_counts_tsv")
    )
    if counts_expected and not has_counts:
        missing.append("merged gene-count matrix (<quant>.merged.gene_counts.tsv)")

    if not missing:
        return
    raise SkillError(
        stage="parsing",
        error_code=ErrorCode.EXPECTED_OUTPUTS_NOT_FOUND,
        message="nf-core/rnaseq completed but required output files were not found.",
        fix=(
            "Inspect logs/stdout.txt, logs/stderr.txt and upstream/work for the "
            "failing step. If quantification was intentionally skipped, use "
            "--skip-quantification-merge or a hisat2/--skip-alignment configuration."
            + _skipped_samples_hint(output_dir)
        ),
        details={
            "output_dir": str(output_dir),
            "aligner_effective": aligner,
            "pseudo_aligner_effective": pseudo,
            "missing_required": missing,
        },
    )


def _skipped_samples_hint(output_dir: Path) -> str:
    """Hint appended when nf-core reports it completed with skipped sample(s).

    nf-core/rnaseq drops samples with fewer than ``min_trimmed_reads`` trimmed
    reads (schema default 10000) before quantification, so an all-skipped run
    completes 'successfully' yet produces no merged count matrix — the exact
    EXPECTED_OUTPUTS_NOT_FOUND we raise here. Read from the captured run log so we
    only surface the hint when that specific signal is present. Best-effort:
    empty string on any read error or when the signal is absent.
    """
    for name in ("stdout.txt", "stderr.txt"):
        try:
            text = (output_dir / "logs" / name).read_text(encoding="utf-8", errors="replace")
        except OSError:
            continue
        if "skipped sample" in text.lower():
            return (
                " nf-core reported it completed with skipped sample(s): samples below "
                "--min-trimmed-reads (nf-core/rnaseq default 10000) are dropped before "
                "quantification, so no merged count matrix is produced. Lower "
                "--min-trimmed-reads (e.g. 0 for tiny test data), or check "
                "upstream/results/multiqc and pipeline_info for which samples were skipped."
            )
    return ""


_BANNER = (
    "==============================================\n"
    f"  ClawBio :: {SKILL_NAME}\n"
    f"  nf-core/rnaseq {DEFAULT_PIPELINE_VERSION} orchestrator\n"
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


def _record_aligner_explicit(args: argparse.Namespace) -> bool:
    """Record (once) whether the user explicitly passed --aligner.

    Must run before any default/override mutates ``args.aligner`` (demo coercion,
    prokaryotic/rapid_quant overrides, the star_salmon default). Self-contained
    nf-core test profiles ship their own aligner, so ``params_builder`` omits
    ``aligner`` from ``params.yaml`` for those runs unless the user chose one
    explicitly (audit F7).
    """
    if not hasattr(args, "_aligner_explicit"):
        args._aligner_explicit = getattr(args, "aligner", None) is not None
    return args._aligner_explicit


def _run_wrapper(args: argparse.Namespace, output_dir: Path) -> int:
    _record_aligner_explicit(args)
    _apply_demo_overrides(args)
    _sync_profile_flags(args)
    _apply_noinput_overrides(args)
    _apply_prokaryotic_overrides(args)
    _apply_rapid_quant_overrides(args)
    _apply_aligner_default(args)
    if args.resume:
        with tempfile.TemporaryDirectory(prefix="clawbio-rnaseq-resume-") as staging:
            return _run_wrapper_with_staging(args, output_dir, staging_dir=Path(staging))
    return _run_wrapper_with_staging(args, output_dir, staging_dir=None)


# Reference / index flags that --demo must clear before they reach params.yaml.
# Why: the upstream `test` profile bundles test FASTQs paired with its own reference
# data. params-file values override profile values in Nextflow, so leaving e.g.
# args.fasta set would silently desynchronise the test samples from their reference
# genome and produce garbage counts. Clearing keeps the demo hermetic.
_DEMO_CLEARED_REFERENCE_FIELDS = (
    "genome",
    "igenomes_base",
    "fasta",
    "gtf",
    "gff",
    "transcript_fasta",
    "additional_fasta",
    "gene_bed",
    "splicesites",
    "star_index",
    "rsem_index",
    "hisat2_index",
    "salmon_index",
    "kallisto_index",
    "bowtie2_index",
    "sortmerna_index",
)


# High-signal tuning flags a user is most likely to *believe* took effect under
# --demo but that the hermetic upstream `test` profile silently owns and ignores
# (a -params-file value would override the profile, so build_effective_params
# deliberately writes none of them in demo mode). Warned — never errored — so a
# demo run is not mistaken for e.g. a UMI-deduplicated or contaminant-screened one
# (audit F-02). Reference/index flags are handled separately above.
_DEMO_IGNORED_TUNING_FIELDS = (
    "pseudo_aligner",
    "contaminant_screening",
    "remove_ribo_rna",
    "ribo_removal_tool",
    "with_umi",
    "skip_umi_extract",
    "gencode",
    "bbsplit_fasta_list",
    "bbsplit_index",
)


def _warn_about_ignored_hermetic_tuning_flags(args: argparse.Namespace, *, mechanism: str = "--demo") -> None:
    """Warn that tuning flags are ignored under a hermetic test profile.

    Both ``--demo`` and a self-contained nf-core test profile (``--profile test``/
    ``test_full``/…) run a profile that owns every pipeline parameter, so any tuning
    flag the user set has no effect. ``mechanism`` names which one fired so the
    message is accurate for either path (audit F1).
    """
    set_flags = [
        "--" + field.replace("_", "-")
        for field in _DEMO_IGNORED_TUNING_FIELDS
        if getattr(args, field, None)
    ]
    # --trimmer defaults to a non-None sentinel (DEFAULT_TRIMMER); only a non-default
    # choice is a meaningful, silently-ignored override.
    if getattr(args, "trimmer", DEFAULT_TRIMMER) not in (None, DEFAULT_TRIMMER):
        set_flags.append("--trimmer")
    if set_flags:
        print(
            f"WARNING: {mechanism} ignores tuning flags ({', '.join(set_flags)}); the upstream "
            "test profile is hermetic and owns every pipeline parameter, so they have no effect.",
            file=sys.stderr,
        )


def _apply_demo_overrides(args: argparse.Namespace) -> None:
    """Coerce demo-incompatible flags before any validation runs.

    --demo runs the upstream `test` profile, which ships with a fixed
    samplesheet and STAR-Salmon-compatible test data. The wrapper enforces:

    * `--aligner` is forced to `star_salmon` (warn if the user requested another).
    * `--input` is dropped (the test profile provides its own samplesheet).
    * All reference/index flags are cleared (the test profile bundles its own refs;
      a partial override — e.g. only --fasta — would desync samples from refs).

    `--resume` is deliberately NOT touched. Nextflow's `-resume` is orthogonal to
    `-profile test` (nf-core documents no incompatibility), and a demo run keeps both
    the things a resume needs inside its own output dir: the work tree
    (`<output>/upstream/work`) and the Nextflow session cache (`<output>/.nextflow`).
    Stripping it here used to make a demo bundle's reproducibility/commands.sh
    unreplayable — the replay hit OUTPUT_DIR_NOT_EMPTY and the suggested `--resume`
    fix was then silently discarded.

    Centralising the override here keeps preflight, samplesheet preparation,
    and params construction consistent — every downstream stage sees the same
    coerced ``args``.
    """
    if not getattr(args, "demo", False):
        return
    requested_aligner = getattr(args, "aligner", None)
    if requested_aligner and requested_aligner != "star_salmon":
        print(
            f"WARNING: --demo forces --aligner star_salmon (requested: {requested_aligner!r}).",
            file=sys.stderr,
        )
    args.aligner = "star_salmon"
    if getattr(args, "input", None):
        print(
            "WARNING: --demo ignores --input; the upstream test profile provides its own samplesheet.",
            file=sys.stderr,
        )
        args.input = None
    cleared = [field for field in _DEMO_CLEARED_REFERENCE_FIELDS if getattr(args, field, None)]
    if cleared:
        flags = ", ".join("--" + field.replace("_", "-") for field in cleared)
        print(
            f"WARNING: --demo ignores reference flags ({flags}); "
            "the upstream test profile provides its own bundled references.",
            file=sys.stderr,
        )
        for field in cleared:
            setattr(args, field, None)
    _warn_about_ignored_hermetic_tuning_flags(args, mechanism="--demo")


def _apply_noinput_overrides(args: argparse.Namespace) -> None:
    """Clear reference/tuning overrides for self-contained nf-core test profiles.

    Self-contained test profiles (``test``, ``test_full``, ``test_prokaryotic``, …)
    ship BOTH their own ``params.input`` and bundled reference data. Because a
    ``-params-file`` value overrides profile config in Nextflow, writing a user
    ``--genome``/``--fasta``/index would silently desynchronise the profile's test
    samples from their matched reference and produce garbage with no error (audit F1).

    ``--demo`` already clears these via ``_apply_demo_overrides``; this mirrors that
    for the ``--profile test*`` path. It must run AFTER ``_sync_profile_flags`` has set
    ``args._noinput`` and is a no-op for ``--demo`` (handled separately) and real runs.
    Reuses ``_DEMO_CLEARED_REFERENCE_FIELDS`` so the cleared set has a single source.
    """
    if not getattr(args, "_noinput", False) or getattr(args, "demo", False):
        return
    cleared = [field for field in _DEMO_CLEARED_REFERENCE_FIELDS if getattr(args, field, None)]
    if cleared:
        flags = ", ".join("--" + field.replace("_", "-") for field in cleared)
        print(
            f"WARNING: a self-contained nf-core test profile owns its own references; "
            f"ignoring reference flags ({flags}). The profile bundles samples paired with "
            "their reference, so overriding would desynchronise them.",
            file=sys.stderr,
        )
        for field in cleared:
            setattr(args, field, None)
    _warn_about_ignored_hermetic_tuning_flags(args, mechanism="a self-contained nf-core test profile")


def _apply_prokaryotic_overrides(args: argparse.Namespace) -> None:
    """Coerce aligner to bowtie2_salmon when --prokaryotic is set and no aligner was specified.

    The nf-core/rnaseq prokaryotic profile (conf/prokaryotic.config) sets
    ``aligner = 'bowtie2_salmon'`` because STAR is a splice-aware eukaryotic
    aligner and is scientifically wrong for prokaryotic RNA-seq.  Nextflow
    params-file values override profile defaults, so writing ``aligner:
    star_salmon`` from the argparse default would silently defeat the profile.

    When the user explicitly passes ``--aligner <anything>`` (including
    ``star_salmon``), the sentinel value ``args.aligner`` is non-None and the
    user's choice is honoured without coercion.
    """
    if not getattr(args, "prokaryotic", False):
        return
    if getattr(args, "aligner", None) is None:
        print(
            "WARNING: --prokaryotic forces --aligner bowtie2_salmon "
            "(STAR is splice-aware and not suitable for prokaryotic RNA-seq).",
            file=sys.stderr,
        )
        args.aligner = "bowtie2_salmon"


_PROKARYOTIC_PROFILE_TOKENS = frozenset({"prokaryotic", "test_prokaryotic"})

# nf-core/rnaseq profiles that ship with params.input in their profile config.
# Running these without --input is valid upstream; the wrapper must not block them.
# NOTE: 'debug' is intentionally excluded — it only toggles debug logging flags
# (dumpHashes, cleanup=false) and carries no params.input or reference data.
_SELF_CONTAINED_TEST_PROFILES = frozenset({
    "test",
    "test_full",
    "test_prokaryotic",
    "test_full_aws",
    "test_full_gcp",
    "test_full_azure",
    "test_gpu",
})


def _sync_profile_flags(args: argparse.Namespace) -> None:
    """Detect modifier keywords in the --profile string and sync the corresponding args flags.

    Users may pass ``--profile prokaryotic,docker``, ``--profile test_prokaryotic``,
    or ``--profile rapid_quant`` instead of the dedicated ``--prokaryotic`` /
    ``--rapid-quant`` flags.  When a modifier keyword appears in the profile string,
    the corresponding flag is set to True so that ``_apply_prokaryotic_overrides`` and
    ``_apply_rapid_quant_overrides`` fire correctly — including their aligner coercion
    and params.yaml writes.

    ``test_prokaryotic`` is an official nf-core/rnaseq profile that includes
    ``conf/prokaryotic.config`` and sets ``aligner = 'bowtie2_salmon'``.
    ``arm64`` is detected so that ``params.yaml`` receives ``arm: true`` even when
    the user passes ``--profile docker,arm64`` instead of ``--arm``.

    Self-contained test profiles (test, test_full, test_prokaryotic, etc.) include
    ``params.input`` in their profile config and do not require a user samplesheet.
    When detected, ``args._noinput`` is set to True so that ``_prepare_samplesheet``
    and ``params_builder`` skip the input requirement.
    """
    profile_parts = {p.strip() for p in getattr(args, "profile", "").split(",") if p.strip()}
    if profile_parts & _PROKARYOTIC_PROFILE_TOKENS and not getattr(args, "prokaryotic", False):
        args.prokaryotic = True
    if "rapid_quant" in profile_parts and not getattr(args, "rapid_quant", False):
        args.rapid_quant = True
    if "arm64" in profile_parts and not getattr(args, "arm", False):
        args.arm = True
    if profile_parts & _SELF_CONTAINED_TEST_PROFILES and not getattr(args, "_noinput", False):
        args._noinput = True


def _apply_aligner_default(args: argparse.Namespace) -> None:
    """Apply the star_salmon default after all profile overrides have run.

    ``--aligner`` defaults to None (sentinel) so that ``_apply_prokaryotic_overrides``
    can distinguish "user did not pass --aligner" from "user explicitly passed
    --aligner star_salmon".  After all override functions have had a chance to
    set a non-None value (e.g. bowtie2_salmon for prokaryotic, star_salmon for demo),
    this function fills in the star_salmon default for any run that still has None.
    """
    if getattr(args, "aligner", None) is None:
        args.aligner = "star_salmon"


def _apply_rapid_quant_overrides(args: argparse.Namespace) -> None:
    """Mirror the rapid_quant profile's implied args so outputs_parser finds the right files.

    conf/rapid_quant.config sets skip_alignment=true, pseudo_aligner='salmon', and
    skip_quantification_merge=true.  Because Nextflow profile values are not visible to
    the Python wrapper, the outputs parser would otherwise look for star_salmon merged
    counts instead of per-sample salmon quant dirs.  Mirroring the profile settings
    into args also writes them explicitly to params.yaml, making the reproducibility
    bundle self-describing.
    """
    if not getattr(args, "rapid_quant", False):
        return
    if not getattr(args, "pseudo_aligner", None):
        args.pseudo_aligner = "salmon"
    if not getattr(args, "skip_alignment", False):
        args.skip_alignment = True
    if not getattr(args, "skip_quantification_merge", False):
        args.skip_quantification_merge = True


def _run_wrapper_with_staging(args: argparse.Namespace, output_dir: Path, *, staging_dir: Path | None) -> int:
    pipeline_source = resolve_pipeline_source(
        requested_version=args.pipeline_version,
        local_pipeline_dir=Path(args.pipeline_local).expanduser().resolve() if args.pipeline_local else None,
    )
    _check_pipeline_source_version(args, pipeline_source)
    normalized_samplesheet, staged_samplesheet, samplesheet_summary = _prepare_samplesheet(
        args,
        output_dir,
        staging_dir=staging_dir,
    )
    preflight_result = run_preflight(args, pipeline_source=pipeline_source, samplesheet_summary=samplesheet_summary)
    _print(
        "[preflight] passed "
        f"(warnings: {len(preflight_result.get('warnings', []) or [])})"
    )
    if args.check:
        if staged_samplesheet != normalized_samplesheet:
            _commit_validated_samplesheet(staged_samplesheet, normalized_samplesheet)
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
    if getattr(args, "_noinput", False):
        return _prepare_noinput_samplesheet(args, output_dir, staging_dir=staging_dir)
    return _prepare_user_samplesheet(args, output_dir, staging_dir=staging_dir)


def _prepare_demo_samplesheet(
    args: argparse.Namespace,
    output_dir: Path,
    *,
    staging_dir: Path | None,
) -> tuple[Path, Path, dict[str, object]]:
    # Demo overrides (aligner=star_salmon, input cleared, reference flags cleared) are
    # applied centrally in ``_apply_demo_overrides`` before this point. The samplesheet
    # written here is a fixed, content-stable stub, so its checksum is identical across
    # runs — which is what lets a demo bundle resume itself on an in-place replay.
    normalized_samplesheet = _final_samplesheet_path(output_dir, demo=True)
    staged_samplesheet = _staged_samplesheet_path(output_dir, staging_dir=staging_dir, demo=True)
    _write_demo_samplesheet(staged_samplesheet)
    return normalized_samplesheet, staged_samplesheet, {
        "normalized_path": normalized_samplesheet,
        "sample_count": 0,
        "sample_names": [],
        "fastq_paths": [],
        "bam_paths": [],
        "unknown_columns": [],
        "strandedness_counts": {},
    }


def _prepare_noinput_samplesheet(
    args: argparse.Namespace,
    output_dir: Path,
    *,
    staging_dir: Path | None,
) -> tuple[Path, Path, dict[str, object]]:
    """Write a provenance stub for self-contained nf-core test profiles.

    Profiles like test_full, test_prokaryotic, and test_full_aws ship with
    params.input in their own profile config.  The wrapper does not write an
    input key to params.yaml for these runs.  This stub (samplesheet.noinput.csv)
    lets auditors distinguish a test-profile run from a real --demo run.
    """
    normalized_samplesheet = _final_samplesheet_path(output_dir, noinput=True)
    staged_samplesheet = _staged_samplesheet_path(output_dir, staging_dir=staging_dir, noinput=True)
    _write_noinput_samplesheet_stub(staged_samplesheet)
    return normalized_samplesheet, staged_samplesheet, {
        "normalized_path": normalized_samplesheet,
        "sample_count": 0,
        "sample_names": [],
        "fastq_paths": [],
        "bam_paths": [],
        "unknown_columns": [],
        "strandedness_counts": {},
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
        skip_alignment=getattr(args, "skip_alignment", False),
    )
    samplesheet_summary["normalized_path"] = normalized_samplesheet
    _warn_about_preserved_unknown_columns(samplesheet_summary)
    return normalized_samplesheet, staged_samplesheet, samplesheet_summary


def _write_demo_samplesheet(samplesheet_path: Path) -> None:
    samplesheet_path.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(samplesheet_path, "sample,fastq_1,fastq_2,strandedness\n")


def _write_noinput_samplesheet_stub(samplesheet_path: Path) -> None:
    samplesheet_path.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(
        samplesheet_path,
        "# Input provided by self-contained nf-core test profile (no --input required)\n",
    )


def _final_samplesheet_path(output_dir: Path, *, demo: bool = False, noinput: bool = False) -> Path:
    if demo:
        filename = "samplesheet.demo.csv"
    elif noinput:
        filename = "samplesheet.noinput.csv"
    else:
        filename = "samplesheet.valid.csv"
    return output_dir / "reproducibility" / filename


def _staged_samplesheet_path(
    output_dir: Path, *, staging_dir: Path | None, demo: bool = False, noinput: bool = False
) -> Path:
    if staging_dir is not None:
        if demo:
            filename = "samplesheet.demo.csv"
        elif noinput:
            filename = "samplesheet.noinput.csv"
        else:
            filename = "samplesheet.valid.csv"
        return staging_dir / filename
    return _final_samplesheet_path(output_dir, demo=demo, noinput=noinput)


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
    _print("[check] Preflight passed. No pipeline was launched (--check).")
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
        gencode_autodetected=bool(preflight_result.get("gencode_autodetected", False)),
    )
    _check_resume_params_checksum(args, output_dir=output_dir, params_payload=params_payload)
    _check_resume_samplesheet_checksum(args, output_dir=output_dir, staged_samplesheet=staged_samplesheet)
    _commit_validated_samplesheet(staged_samplesheet, normalized_samplesheet)
    params_path = write_params_yaml(params_payload, output_dir=output_dir)
    command, command_str = _build_nextflow_invocation(args, output_dir, pipeline_source, params_path)
    nextflow_cwd = _nextflow_execution_cwd(output_dir)
    started = time.monotonic()
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
    duration_seconds = round(time.monotonic() - started, 3)
    _print(f"[execute] completed in {duration_seconds:.1f}s")
    _print("[outputs] parsing pipeline outputs")
    parsed_outputs = _parse_outputs_with_effective_aligner(output_dir, args)
    _raise_if_expected_outputs_missing(parsed_outputs, args=args, output_dir=output_dir)
    _print("[report] writing report.md, result.json, commands.sh")
    _print("[provenance] writing reproducibility bundle")
    if preflight_result.get("handoff_available") is False:
        parsed_outputs["handoff_available"] = False
    handoff_result = _run_downstream_handoff(args, parsed_outputs=parsed_outputs, output_dir=output_dir)
    _record_downstream_handoff_result(parsed_outputs, handoff_result)
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
        duration_seconds=duration_seconds,
        command_str=_nextflow_replay_command(command_str, nextflow_cwd),
    )
    _print(f"[done] Wrapper completed successfully. Output: {output_dir}")
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
    if not getattr(args, "resume", False):
        return
    check_resume_params_checksum(params_payload, output_dir)


def _check_resume_samplesheet_checksum(
    args: argparse.Namespace,
    *,
    output_dir: Path,
    staged_samplesheet: Path,
) -> None:
    if not getattr(args, "resume", False):
        return
    check_resume_samplesheet_checksum(staged_samplesheet, output_dir)


def _build_nextflow_invocation(
    args: argparse.Namespace,
    output_dir: Path,
    pipeline_source: dict[str, object],
    params_path: Path,
) -> tuple[list[str], str]:
    return build_nextflow_command(
        pipeline_source=pipeline_source,
        profile=args.profile,
        params_path=params_path,
        resume=args.resume,
        work_dir=_resolve_nextflow_work_dir(args, output_dir),
        extra_configs=_build_extra_nextflow_configs(args, output_dir),
        demo=bool(getattr(args, "demo", False)),
        prokaryotic=bool(getattr(args, "prokaryotic", False)),
        rapid_quant=bool(getattr(args, "rapid_quant", False)),
        arm=bool(getattr(args, "arm", False)),
    )


def _build_extra_nextflow_configs(args: argparse.Namespace, output_dir: Path) -> list[Path]:
    configs: list[Path] = []
    profile_parts = {p.strip() for p in getattr(args, "profile", "").split(",") if p.strip()}
    arm = bool(getattr(args, "arm", False))
    if "docker" in profile_parts:
        if platform.system() == "Darwin":
            configs.append(_write_macos_docker_config(output_dir, arm=arm, timeout_hours=_resolve_timeout_hours(args)))
        elif not bool(getattr(args, "demo", False)):
            # Non-macOS docker runs use the local executor, which aborts a real run
            # when a process's production requirement exceeds this host's RAM (e.g.
            # PREPARE_GENOME requesting 72 GB on a 63 GB machine). Cap requests to the
            # host per nf-core resourceLimits guidance. --demo relies on -profile test's
            # own (smaller) resourceLimits, so it is deliberately left untouched.
            configs.append(_write_resource_limits_config(output_dir, timeout_hours=_resolve_timeout_hours(args)))
    for user_cfg in getattr(args, "nextflow_config", None) or []:
        # --nextflow-config accepts local paths only; Path.resolve() is intentional.
        # Nextflow's own -c flag also expects local file paths in typical use.
        configs.append(Path(user_cfg).expanduser().resolve())
    return configs


_MACOS_DOCKER_MEMORY_CEILING_GB = 15
_MACOS_DOCKER_MEMORY_FLOOR_GB = 8


def _detect_host_memory_gb() -> int | None:
    """Return total physical RAM in whole GB, or None when it cannot be determined."""
    try:
        page_size = os.sysconf("SC_PAGE_SIZE")
        phys_pages = os.sysconf("SC_PHYS_PAGES")
    except (ValueError, OSError, AttributeError):
        return None
    if page_size <= 0 or phys_pages <= 0:
        return None
    return int((page_size * phys_pages) / (1024 ** 3))


def _docker_vm_memory_gb() -> int | None:
    """Return the Docker VM's total memory in whole GB, or None when unknown.

    On macOS, Docker Desktop runs containers in a Linux VM whose RAM is usually
    smaller than the host's. `docker info` reports that VM's MemTotal (bytes), which
    is the true ceiling a container process can use before being OOM-killed
    (audit F-8). Best-effort: any failure (docker absent, daemon down, unparseable
    output) returns None so the caller falls back to the host-derived heuristic.
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


def _macos_docker_memory_gb() -> int:
    """Pick a per-process memory ceiling for the macOS Docker resourceLimits block.

    The ceiling never exceeds the historical 15 GB default (Docker Desktop VMs are
    typically smaller than the host), but it is lowered toward a 75% share on hosts
    with less RAM so the limit cannot exceed what the machine physically has
    (audit F-07). A floor keeps it above the pipeline's per-process minimums.

    When the Docker VM's actual memory is known it overrides everything (including
    the floor): the limit is capped to 90% of the VM so a container process is not
    OOM-killed by requesting more than the VM physically has (audit F-8).
    """
    host_gb = _detect_host_memory_gb()
    if not host_gb:
        budget = _MACOS_DOCKER_MEMORY_CEILING_GB
    else:
        budget = max(_MACOS_DOCKER_MEMORY_FLOOR_GB, min(_MACOS_DOCKER_MEMORY_CEILING_GB, int(host_gb * 0.75)))
    vm_gb = _docker_vm_memory_gb()
    if vm_gb:
        # The VM ceiling overrides the floor — requesting more than the VM has would
        # OOM-kill the process, which is worse than a smaller-but-valid limit.
        budget = min(budget, max(1, int(vm_gb * 0.9)))
    return budget


def _write_macos_docker_config(output_dir: Path, *, arm: bool = False, timeout_hours: float = DEFAULT_TIMEOUT_HOURS) -> Path:
    config_path = output_dir / ".nextflow_macos_docker.config"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    # --platform linux/amd64 forces qemu emulation on Apple Silicon; omit when using
    # ARM64-native containers (--arm / arm64 profile) so they run at native speed.
    platform_opts = "" if arm else "    containerOptions = '--platform linux/amd64'\n"
    docker_block = "" if arm else "docker {\n    runOptions = '--platform linux/amd64'\n}\n"
    # Autodetect CPU ceiling so users on multi-core Macs aren't capped at 4. Floor of
    # 4 preserves prior behaviour on tiny VMs. Memory/time stay conservative because
    # over-allocating memory can trigger OOM-kills (worse than throttling); users on
    # high-RAM hosts can override via --nextflow-config.
    cpus = max(4, os.cpu_count() or 4)
    memory_gb = _macos_docker_memory_gb()
    # The per-process time ceiling tracks --timeout-hours (audit F-4) so a large-cohort
    # run raised above the 12h default is not capped back to 12h by this config. Floored
    # at 1h so a sub-hour --timeout-hours never emits an invalid '0.h' duration.
    # --timeout-hours 0 disables the wall-clock cap, so the per-process ceiling falls back
    # to the pinned default rather than being floored to a misleading 1h.
    time_hours = int(DEFAULT_TIMEOUT_HOURS) if timeout_hours == 0 else max(1, int(timeout_hours))
    write_text_lf(
        config_path,
        "// macOS Docker compatibility for nf-core/rnaseq.\n"
        "process {\n"
        "    stageInMode = 'copy'\n"
        "    resourceLimits = [\n"
        f"        cpus: {cpus},\n"
        f"        memory: '{memory_gb}.GB',\n"
        f"        time: '{time_hours}.h'\n"
        "    ]\n"
        + platform_opts
        + "}\n"
        + docker_block,
    )
    return config_path


_RESOURCE_LIMITS_MEMORY_HEADROOM_GB = 4


def _resource_limits_memory_gb() -> int:
    """Per-process memory ceiling (GB) for the resourceLimits config on non-macOS hosts.

    nf-core's guidance is that ``resourceLimits`` should match the machine's maximum
    (https://nf-co.re/docs/running/configuration/nextflow-for-your-system). Unlike the
    macOS Docker Desktop VM, native Linux Docker exposes the full host RAM, so there is
    no 15 GB ceiling here — clamping a large host to 15 GB would starve genuine
    process_high tasks. We take the smaller of the detected physical RAM and (when
    Docker runs natively) the Docker runtime's MemTotal, then reserve a few GB of
    headroom for the OS, the Nextflow JVM and the Docker daemon so a process using its
    full budget cannot OOM the host. Falls back to the conservative ceiling when the
    host size cannot be determined.
    """
    signals = [gb for gb in (_detect_host_memory_gb(), _docker_vm_memory_gb()) if gb]
    if not signals:
        return _MACOS_DOCKER_MEMORY_CEILING_GB
    budget = min(signals) - _RESOURCE_LIMITS_MEMORY_HEADROOM_GB
    return max(_MACOS_DOCKER_MEMORY_FLOOR_GB, budget)


def _write_resource_limits_config(output_dir: Path, *, timeout_hours: float = DEFAULT_TIMEOUT_HOURS) -> Path:
    """Write a host-scaled ``resourceLimits`` config for non-macOS docker runs.

    Carries only the portable resourceLimits block — none of the macOS-only
    workarounds (``--platform linux/amd64``, ``stageInMode = 'copy'``), which would
    needlessly slow or break a native Linux run.
    """
    config_path = output_dir / ".nextflow_resource_limits.config"
    config_path.parent.mkdir(parents=True, exist_ok=True)
    cpus = max(1, os.cpu_count() or 1)
    memory_gb = _resource_limits_memory_gb()
    time_hours = int(DEFAULT_TIMEOUT_HOURS) if timeout_hours == 0 else max(1, int(timeout_hours))
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


def _resolve_timeout_hours(args: argparse.Namespace) -> float:
    """Return the effective --timeout-hours, falling back to the pinned default."""
    hours = getattr(args, "timeout_hours", DEFAULT_TIMEOUT_HOURS)
    return DEFAULT_TIMEOUT_HOURS if hours is None else float(hours)


def _resolve_nextflow_work_dir(args: argparse.Namespace, output_dir: Path) -> Path | str:
    """Resolve the Nextflow work directory: default <output>/upstream/work, or the
    user's --work-dir (a local path, or an object-store URI kept verbatim for cloud
    executors). Parity with nfcore-scrnaseq/sarek."""
    raw_work_dir = getattr(args, "work_dir", None)
    if not raw_work_dir:
        return output_dir / "upstream" / "work"
    work_dir = str(raw_work_dir).strip()
    if "://" in work_dir:
        return work_dir
    return Path(work_dir).expanduser().resolve()


def _resolve_timeout_seconds(args: argparse.Namespace) -> int | None:
    """Translate --timeout-hours into seconds; 0 disables the cap (returns None).

    Parity with nfcore-scrnaseq/sarek: the default keeps a wall-clock cap, but an
    explicit ``--timeout-hours 0`` opts out (long HPC/cloud runs governed by the
    scheduler's walltime).
    """
    hours = _resolve_timeout_hours(args)
    if hours == 0:
        return None
    return int(hours * 60 * 60)


def _nextflow_execution_cwd(output_dir: Path) -> Path:
    return output_dir


def _nextflow_replay_command(command_str: str, cwd: Path) -> str:
    return f"cd {shlex.quote(cwd.as_posix())} && {command_str}"


def _parse_outputs_with_effective_aligner(output_dir: Path, args: argparse.Namespace) -> dict[str, object]:
    aligner = getattr(args, "aligner", "star_salmon")
    parsed = parse_outputs(
        output_dir,
        aligner=aligner,
        pseudo_aligner=getattr(args, "pseudo_aligner", None),
        skip_alignment=bool(getattr(args, "skip_alignment", False)),
        skip_quantification_merge=bool(getattr(args, "skip_quantification_merge", False)),
    )
    parsed["aligner_effective"] = aligner
    # parse_outputs already sets handoff_available (counts present AND not a skipped
    # merge); the only adjustment here is forcing it off for the HISAT2 alignment-only
    # route, which produces BAMs but no merged quantification.
    hisat2_no_quant = aligner == "hisat2" and not parsed.get("preferred_counts_tsv")
    parsed["hisat2_no_quant"] = hisat2_no_quant
    if hisat2_no_quant:
        parsed["handoff_available"] = False
    return parsed


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
    duration_seconds: float = 0,
) -> None:
    post_run_warnings: list[str] = []
    write_repro_commands(output_dir, args=args)
    try:
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
            duration_seconds=duration_seconds,
            command_str=command_str,
        )
    except Exception as exc:
        post_run_warnings.append(
            f"Provenance bundle could not be written ({type(exc).__name__}: {exc}). "
            "Run --repair-bundle to regenerate missing files."
        )
    write_report(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
        post_run_warnings=post_run_warnings,
    )
    write_result(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
        post_run_warnings=post_run_warnings,
    )


def _run_downstream_handoff(
    args: argparse.Namespace,
    *,
    parsed_outputs: dict[str, object],
    output_dir: Path,
) -> dict[str, object] | None:
    if not getattr(args, "run_downstream", False) or getattr(args, "skip_downstream", False):
        return None
    counts = str(parsed_outputs.get("preferred_counts_tsv", "")).strip()
    if not counts or not parsed_outputs.get("handoff_available", bool(counts)):
        return None
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    handoff_path = repro_dir / "rnaseq_de_handoff.sh"
    write_text_lf(
        handoff_path,
        "#!/usr/bin/env bash\n"
        "# Downstream rnaseq-de handoff — generated by nfcore-rnaseq-wrapper\n"
        "# Usage: CLAWBIO_REPO=/path/to/ClawBio bash rnaseq_de_handoff.sh\n"
        "set -euo pipefail\n\n"
        'if [[ -z "${CLAWBIO_REPO:-}" ]]; then\n'
        '  echo "ERROR: CLAWBIO_REPO is not set." >&2\n'
        '  echo "  CLAWBIO_REPO=/path/to/ClawBio bash ${BASH_SOURCE[0]}" >&2\n'
        "  exit 1\n"
        "fi\n\n"
        '"${PYTHON:-python3}" "${CLAWBIO_REPO}/clawbio.py" run rnaseq \\\n'
        f"  --counts {shlex.quote(counts)} \\\n"
        "  --metadata <your_metadata.csv> \\\n"
        '  --formula "~ batch + condition" \\\n'
        '  --contrast "condition,treated,control" \\\n'
        "  --output <dir>\n",
    )
    if not _downstream_flags_complete(args):
        return _downstream_handoff_result(template_path=handoff_path, status="template_only")
    downstream_output = (
        Path(getattr(args, "downstream_output")).expanduser().resolve()
        if getattr(args, "downstream_output", None)
        else output_dir / "rnaseq_de"
    )
    clawbio_path = _PROJECT_ROOT / "clawbio.py"
    if not clawbio_path.exists():
        print(
            f"WARNING: clawbio.py not found at {clawbio_path}; downstream handoff skipped.",
            file=sys.stderr,
        )
        return _downstream_handoff_result(template_path=handoff_path, status="error")
    command = [
        sys.executable,
        str(clawbio_path),
        "run",
        "rnaseq",
        "--counts",
        counts,
        "--metadata",
        str(getattr(args, "metadata")),
        "--formula",
        str(getattr(args, "formula")),
        "--contrast",
        str(getattr(args, "contrast")),
        "--output",
        str(downstream_output),
    ]
    try:
        result = subprocess.run(
            command,
            capture_output=True,
            text=True,
            timeout=_DOWNSTREAM_HANDOFF_TIMEOUT_SECONDS,
            cwd=str(_PROJECT_ROOT),
        )
    except subprocess.TimeoutExpired:
        print(
            f"WARNING: rnaseq-de downstream handoff timed out after {_DOWNSTREAM_HANDOFF_TIMEOUT_SECONDS} seconds.",
            file=sys.stderr,
        )
        return _downstream_handoff_result(template_path=handoff_path, status="timeout")
    except OSError as exc:
        print(f"WARNING: rnaseq-de downstream handoff could not start: {exc}", file=sys.stderr)
        return _downstream_handoff_result(template_path=handoff_path, status="error")
    else:
        if result.returncode != 0:
            detail = (result.stderr or result.stdout or "").strip()
            suffix = f": {detail}" if detail else ""
            print(f"WARNING: rnaseq-de downstream handoff failed with exit code {result.returncode}{suffix}", file=sys.stderr)
            return _downstream_handoff_result(
                template_path=handoff_path,
                status="failed",
                returncode=result.returncode,
            )
    return _downstream_handoff_result(
        template_path=handoff_path,
        status="completed",
        output_dir=downstream_output,
        returncode=0,
    )


def _downstream_handoff_result(
    *,
    template_path: Path | None,
    status: str,
    output_dir: Path | None = None,
    returncode: int | None = None,
) -> dict[str, object]:
    return {
        "template_path": str(template_path) if template_path is not None else None,
        "downstream_output_dir": str(output_dir) if output_dir is not None else None,
        "downstream_status": status,
        "downstream_returncode": returncode,
    }


def _record_downstream_handoff_result(
    parsed_outputs: dict[str, object],
    handoff_result: dict[str, object] | None,
) -> None:
    if not handoff_result:
        return
    template_path = handoff_result.get("template_path")
    if template_path:
        parsed_outputs["rnaseq_de_handoff"] = str(template_path)
    status = handoff_result.get("downstream_status")
    if status:
        parsed_outputs["rnaseq_de_status"] = str(status)
    returncode = handoff_result.get("downstream_returncode")
    if returncode is not None:
        parsed_outputs["rnaseq_de_returncode"] = returncode
    output_dir = handoff_result.get("downstream_output_dir")
    if status == "completed" and output_dir:
        parsed_outputs["rnaseq_de_output_dir"] = str(output_dir)


def _downstream_flags_complete(args: argparse.Namespace) -> bool:
    return bool(getattr(args, "metadata", None) and getattr(args, "formula", None) and getattr(args, "contrast", None))


def _write_error_result_if_safe(output_dir: Path, payload: dict[str, object]) -> None:
    error_code = payload.get("error_code")
    if error_code in {ErrorCode.OUTPUT_DIR_NOT_EMPTY, ErrorCode.OUTPUT_DIR_NOT_WRITABLE}:
        return
    if error_code == ErrorCode.INVALID_RESUME_STATE and (output_dir / "result.json").exists():
        return
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_text_lf(output_dir / "result.json", json.dumps(payload, indent=2))
    except OSError:
        return


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
        "fix": "Report this as a bug. Include the traceback and command arguments.",
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
