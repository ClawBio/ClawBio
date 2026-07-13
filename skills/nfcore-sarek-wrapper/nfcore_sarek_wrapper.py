#!/usr/bin/env python3
# ruff: noqa: E402
# nfcore-sarek-wrapper / nfcore_sarek_wrapper.py
"""Orchestrator for the nf-core/sarek 3.8.1 ClawBio wrapper.

Ties together all submodules — samplesheet, preflight, params, command, executor,
outputs, provenance, reporting — behind a single CLI.

Run with ``--help`` to see every supported flag.
"""
from __future__ import annotations

import argparse
import json
import os
import shlex
import shutil
import subprocess
import sys
import time
import traceback
from pathlib import Path
from typing import Any, cast
from urllib.parse import urlsplit


_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
# Bootstrap the repo root (which holds the `clawbio` package) onto sys.path so the
# documented direct invocation — `python skills/nfcore-sarek-wrapper/
# nfcore_sarek_wrapper.py --help` — resolves `from clawbio.common...` without the
# caller setting PYTHONPATH. Running via `clawbio.py run sarek-pipeline` already has
# the root on the path, so the guard makes this a no-op there. Mirrors
# nfcore-scrnaseq-wrapper.
_PROJECT_ROOT = _SKILL_DIR.parent.parent
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))


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

from clawbio.common.textio import write_text_lf
from command_builder import build_nextflow_command, compose_profile
from errors import ErrorCode, SkillError
from executor import execute_nextflow
from outputs_parser import parse_outputs
from params_builder import _load_user_params_file, build_effective_params, write_params_yaml
from pipeline_source import resolve_pipeline_source
from preflight import REFERENCE_PATH_PARAMS, run_preflight
from provenance import compute_params_checksum, compute_reference_checksums, load_manifest, write_provenance_bundle
from reporting import write_reports
from samplesheet_builder import validate_and_normalize_samplesheet
from schemas import (
    DEFAULT_ALIGNER,
    DEFAULT_PIPELINE_VERSION,
    DEFAULT_PROFILE,
    DEFAULT_TIMEOUT_SECONDS,
    SKILL_NAME,
    SKILL_VERSION,
    SUPPORTED_ALIGNERS,
    SUPPORTED_GROUP_BY_UMI,
    SUPPORTED_PUBLISH_DIR_MODES,
    SUPPORTED_UMI_LOCATIONS,
    SUPPORTED_VEP_OUT_FORMAT,
)


_DOWNSTREAM_CHOICES = (
    "clinical-variant-reporter",
    "clinical-trial-finder",
    "omics-target-evidence-mapper",
    "wes-clinical-report-en",
    "wes-clinical-report-es",
)

_RESTART_CSV_BY_STEP = {
    "markduplicates": "mapped.csv",
    "prepare_recalibration": "markduplicates_no_table.csv",
    "recalibrate": "markduplicates.csv",
    "variant_calling": "recalibrated.csv",
    "annotate": "variantcalled.csv",
}


# Sarek schema params exposed as wrapper flags. Each row is
# (cli_flag, dest, type, default, help, group).
# The wrapper exposes the flags that the upstream `preflight`, `params_builder`,
# and §5.4 flag-compatibility rules touch directly. Other schema parameters can
# be supplied via ``--extra-param key=value``.
_SAREK_PASSTHROUGH_PARAMS: tuple[tuple[str, str, type, Any, str, str], ...] = (
    # -------- Main --------
    ("--step", "step", str, None, "Pipeline step", "main"),
    ("--tools", "tools", str, None, "Comma-separated tool list (haplotypecaller, mutect2, strelka, ...)", "main"),
    ("--skip-tools", "skip_tools", str, None, "Comma-separated skip-tools list", "main"),
    ("--aligner", "aligner", str, None, "Aligner: " + ",".join(sorted(SUPPORTED_ALIGNERS)), "main"),
    ("--no-intervals", "no_intervals", bool, False, "Disable intervals processing", "main"),
    ("--wes", "wes", bool, False, "Whole-exome/panel mode (provide a target --intervals BED; enforced as BED when given)", "main"),
    ("--joint-germline", "joint_germline", bool, False, "Joint germline variant calling", "main"),
    ("--joint-mutect2", "joint_mutect2", bool, False, "Joint Mutect2 variant calling", "main"),
    ("--only-paired-variant-calling", "only_paired_variant_calling", bool, False, "Only paired variant calling", "main"),
    ("--ignore-soft-clipped-bases", "ignore_soft_clipped_bases", bool, False, "Ignore soft-clipped bases", "main"),
    ("--filter-vcfs", "filter_vcfs", bool, False, "Filter VCFs after calling", "main"),
    ("--normalize-vcfs", "normalize_vcfs", bool, False, "Normalize VCFs after calling", "main"),
    ("--snv-consensus-calling", "snv_consensus_calling", bool, False, "Enable SNV consensus", "main"),
    ("--concatenate-vcfs", "concatenate_vcfs", bool, False, "Concatenate VCFs", "main"),
    ("--build-only-index", "build_only_index", bool, False, "Build reference indices only", "main"),
    ("--download-cache", "download_cache", bool, False, "Download VEP/snpEff caches", "main"),
    ("--use-gatk-spark", "use_gatk_spark", str, None, "Spark tool: baserecalibrator|markduplicates (NOT the spark profile)", "main"),

    # -------- I/O & metadata --------
    ("--seq-platform", "seq_platform", str, None, "Sequencing platform tag", "io"),
    ("--seq-center", "seq_center", str, None, "Sequencing center tag", "io"),
    ("--email", "email", str, None, "Notification email", "io"),
    ("--email-on-fail", "email_on_fail", str, None, "Failure-only email", "io"),
    ("--publish-dir-mode", "publish_dir_mode", str, None,
        "Nextflow publishDir mode: " + ",".join(sorted(SUPPORTED_PUBLISH_DIR_MODES)), "io"),
    ("--outdir-cache", "outdir_cache", str, None, "Override cache outdir", "io"),
    ("--multiqc-title", "multiqc_title", str, None, "MultiQC report title", "io"),
    ("--multiqc-config", "multiqc_config", str, None, "MultiQC config YAML", "io"),
    ("--multiqc-logo", "multiqc_logo", str, None, "MultiQC logo", "io"),
    ("--multiqc-methods-description", "multiqc_methods_description", str, None, "MultiQC methods description", "io"),
    ("--hook-url", "hook_url", str, None, "Slack/Teams hook URL", "io"),
    ("--trace-report-suffix", "trace_report_suffix", str, None, "Custom trace report suffix", "io"),
    ("--max-multiqc-email-size", "max_multiqc_email_size", str, None, "MultiQC email attachment size cap", "io"),
    ("--input-restart", "input_restart", str, None, "Restart CSV to validate and supply as pipeline input", "io"),

    # -------- FASTQ preprocessing --------
    ("--trim-fastq", "trim_fastq", bool, False, "Trim FASTQ reads with fastp", "trim"),
    ("--trim-nextseq", "trim_nextseq", bool, False, "fastp --trim_poly_g for NextSeq/NovaSeq", "trim"),
    ("--clip-r1", "clip_r1", int, None, "Clip N bases from R1 5'", "trim"),
    ("--clip-r2", "clip_r2", int, None, "Clip N bases from R2 5'", "trim"),
    ("--three-prime-clip-r1", "three_prime_clip_r1", int, None, "3' clip R1", "trim"),
    ("--three-prime-clip-r2", "three_prime_clip_r2", int, None, "3' clip R2", "trim"),
    ("--length-required", "length_required", int, None, "Minimum read length after trimming", "trim"),
    ("--split-fastq", "split_fastq", int, None, "Split FASTQ chunks", "trim"),
    ("--save-split-fastqs", "save_split_fastqs", bool, False, "Save split FASTQs", "trim"),
    ("--save-trimmed", "save_trimmed", bool, False, "Save trimmed FASTQs", "trim"),

    # -------- UMI --------
    ("--umi-read-structure", "umi_read_structure", str, None, "fgbio read structure (e.g. 8M+T 8M+T)", "umi"),
    ("--group-by-umi-strategy", "group_by_umi_strategy", str, None,
        "fgbio GroupReadsByUmi strategy: " + ",".join(sorted(SUPPORTED_GROUP_BY_UMI)), "umi"),
    ("--umi-location", "umi_location", str, None,
        "UMI location: " + ",".join(sorted(SUPPORTED_UMI_LOCATIONS)), "umi"),
    ("--umi-tag", "umi_tag", str, None, "BAM tag carrying UMI", "umi"),
    ("--umi-length", "umi_length", int, None, "UMI length", "umi"),
    ("--umi-base-skip", "umi_base_skip", int, None, "Bases to skip after UMI", "umi"),
    ("--umi-in-read-header", "umi_in_read_header", bool, False, "UMI is in read header", "umi"),
    ("--sentieon-consensus", "sentieon_consensus", bool, False, "Sentieon UMI consensus", "umi"),

    # -------- Preprocessing --------
    ("--save-mapped", "save_mapped", bool, False, "Save aligned BAM/CRAM", "preprocessing"),
    ("--save-output-as-bam", "save_output_as_bam", bool, False, "Save preprocessing output as BAM (default: CRAM)", "preprocessing"),
    ("--markduplicates-pixel-distance", "markduplicates_pixel_distance", int, None, "MarkDuplicates pixel distance", "preprocessing"),
    ("--nucleotides-per-second", "nucleotides_per_second", int, None, "Throughput hint", "preprocessing"),

    # -------- Variant calling tuning --------
    ("--ascat-min-base-qual", "ascat_min_base_qual", int, None, "ASCAT min base quality", "varcall"),
    ("--ascat-min-counts", "ascat_min_counts", int, None, "ASCAT min counts", "varcall"),
    ("--ascat-min-map-qual", "ascat_min_map_qual", int, None, "ASCAT min mapping quality", "varcall"),
    ("--ascat-ploidy", "ascat_ploidy", float, None, "ASCAT ploidy", "varcall"),
    ("--ascat-purity", "ascat_purity", float, None, "ASCAT purity", "varcall"),
    ("--ascat-genome", "ascat_genome", str, None, "ASCAT genome (hg19|hg38)", "varcall"),
    ("--cf-coeff", "cf_coeff", float, None, "ControlFreec coeff", "varcall"),
    ("--cf-contamination", "cf_contamination", int, None, "ControlFreec contamination", "varcall"),
    ("--cf-contamination-adjustment", "cf_contamination_adjustment", bool, False, "ControlFreec contamination adjustment", "varcall"),
    ("--cf-minqual", "cf_minqual", int, None, "ControlFreec min qual", "varcall"),
    ("--cf-mincov", "cf_mincov", int, None, "ControlFreec min cov", "varcall"),
    ("--cf-ploidy", "cf_ploidy", str, None, "ControlFreec ploidy", "varcall"),
    ("--cf-window", "cf_window", float, None, "ControlFreec window (upstream type: number)", "varcall"),
    ("--cf-chrom-len", "cf_chrom_len", str, None, "ControlFreec chromosome length file", "varcall"),
    ("--cnvkit-reference", "cnvkit_reference", str, None, "CNVkit pre-built reference", "varcall"),
    ("--freebayes-filter", "freebayes_filter", str, None, "Freebayes vcflib/vcffilter expression (upstream type: string, default '30')", "varcall"),
    ("--sentieon-haplotyper-emit-mode", "sentieon_haplotyper_emit_mode", str, None, "Sentieon Haplotyper emit mode", "varcall"),
    ("--sentieon-dnascope-emit-mode", "sentieon_dnascope_emit_mode", str, None, "Sentieon DNAscope emit mode", "varcall"),
    ("--sentieon-dnascope-pcr-indel-model", "sentieon_dnascope_pcr_indel_model", str, None, "Sentieon DNAscope PCR indel model", "varcall"),
    ("--gatk-pcr-indel-model", "gatk_pcr_indel_model", str, None, "GATK PCR indel model", "varcall"),
    ("--varlociraptor-chunk-size", "varlociraptor_chunk_size", int, None, "Varlociraptor chunk size", "varcall"),
    ("--varlociraptor-scenario-tumor-only", "varlociraptor_scenario_tumor_only", str, None, "Varlociraptor tumor-only scenario YAML", "varcall"),
    ("--varlociraptor-scenario-somatic", "varlociraptor_scenario_somatic", str, None, "Varlociraptor somatic scenario YAML", "varcall"),
    ("--varlociraptor-scenario-germline", "varlociraptor_scenario_germline", str, None, "Varlociraptor germline scenario YAML", "varcall"),
    ("--consensus-min-count", "consensus_min_count", int, None, "Consensus min count for SNV consensus calling", "varcall"),

    # -------- Post-variant calling --------
    ("--bcftools-filter-criteria", "bcftools_filter_criteria", str, None, "bcftools view -e expression", "postvar"),
    ("--bcftools-columns", "bcftools_columns", str, None, "bcftools annotate -c", "postvar"),
    ("--bcftools-header-lines", "bcftools_header_lines", str, None, "Path to extra header lines for bcftools annotate", "postvar"),

    # -------- Annotation --------
    ("--snpeff-db", "snpeff_db", str, None, "snpEff database identifier", "annotation"),
    ("--snpeff-cache", "snpeff_cache", str, None, "Path to snpEff cache", "annotation"),
    ("--vep-cache", "vep_cache", str, None, "Path to VEP cache", "annotation"),
    ("--vep-cache-version", "vep_cache_version", str, None, "VEP cache version", "annotation"),
    ("--vep-genome", "vep_genome", str, None, "VEP genome", "annotation"),
    ("--vep-species", "vep_species", str, None, "VEP species", "annotation"),
    ("--vep-version", "vep_version", str, None, "VEP version", "annotation"),
    ("--vep-out-format", "vep_out_format", str, None,
        "VEP output format: " + ",".join(sorted(SUPPORTED_VEP_OUT_FORMAT)), "annotation"),
    ("--vep-custom-args", "vep_custom_args", str, None, "Extra VEP args", "annotation"),
    ("--vep-include-fasta", "vep_include_fasta", bool, False, "Use FASTA with VEP", "annotation"),
    ("--vep-condel", "vep_condel", bool, False, "VEP Condel plugin", "annotation"),
    ("--vep-dbnsfp", "vep_dbnsfp", bool, False, "VEP dbNSFP plugin", "annotation"),
    ("--vep-loftee", "vep_loftee", bool, False, "VEP LOFTEE plugin", "annotation"),
    ("--vep-mastermind", "vep_mastermind", bool, False, "VEP Mastermind plugin", "annotation"),
    ("--vep-phenotypes", "vep_phenotypes", bool, False, "VEP Phenotypes plugin", "annotation"),
    ("--vep-spliceai", "vep_spliceai", bool, False, "VEP SpliceAI plugin", "annotation"),
    ("--vep-spliceregion", "vep_spliceregion", bool, False, "VEP SpliceRegion plugin", "annotation"),
    ("--dbnsfp", "dbnsfp", str, None, "dbNSFP file path", "annotation"),
    ("--dbnsfp-tbi", "dbnsfp_tbi", str, None, "dbNSFP TBI index", "annotation"),
    ("--dbnsfp-consequence", "dbnsfp_consequence", str, None, "dbNSFP consequence filter", "annotation"),
    ("--dbnsfp-fields", "dbnsfp_fields", str, None, "dbNSFP fields to extract", "annotation"),
    ("--mastermind-file", "mastermind_file", str, None, "Mastermind VCF", "annotation"),
    ("--mastermind-mutations", "mastermind_mutations", bool, False, "Mastermind mutations plugin flag", "annotation"),
    ("--mastermind-var-iden", "mastermind_var_iden", bool, False, "Mastermind variant-identifier plugin flag", "annotation"),
    ("--mastermind-url", "mastermind_url", bool, False, "Mastermind URL plugin flag", "annotation"),
    ("--phenotypes-file", "phenotypes_file", str, None, "VEP Phenotypes file", "annotation"),
    ("--phenotypes-file-tbi", "phenotypes_file_tbi", str, None, "VEP Phenotypes TBI", "annotation"),
    ("--phenotypes-include-types", "phenotypes_include_types", str, None, "Phenotype types to include", "annotation"),
    ("--spliceai-snv", "spliceai_snv", str, None, "SpliceAI SNV VCF", "annotation"),
    ("--spliceai-snv-tbi", "spliceai_snv_tbi", str, None, "SpliceAI SNV TBI", "annotation"),
    ("--spliceai-indel", "spliceai_indel", str, None, "SpliceAI indel VCF", "annotation"),
    ("--spliceai-indel-tbi", "spliceai_indel_tbi", str, None, "SpliceAI indel TBI", "annotation"),
    ("--bcftools-annotations", "bcftools_annotations", str, None, "bcftools annotations VCF", "annotation"),
    ("--bcftools-annotations-tbi", "bcftools_annotations_tbi", str, None, "bcftools annotations TBI", "annotation"),
    ("--condel-config", "condel_config", str, None, "Condel config", "annotation"),
    ("--snpsift-databases", "snpsift_databases", str, None, "SnpSift annotation databases", "annotation"),

    # -------- Reference --------
    ("--genome", "genome", str, None, "iGenomes key (e.g. GATK.GRCh38)", "reference"),
    ("--igenomes-base", "igenomes_base", str, None, "iGenomes base path", "reference"),
    ("--igenomes-ignore", "igenomes_ignore", bool, False, "Ignore iGenomes", "reference"),
    ("--fasta", "fasta", str, None, "Reference FASTA", "reference"),
    ("--fasta-fai", "fasta_fai", str, None, "Reference FASTA .fai", "reference"),
    ("--dict", "dict", str, None, "Reference .dict", "reference"),
    ("--bwa", "bwa", str, None, "BWA index dir", "reference"),
    ("--bwamem2", "bwamem2", str, None, "BWA-MEM2 index dir", "reference"),
    ("--dragmap", "dragmap", str, None, "DragMap hash table dir", "reference"),
    ("--dbsnp", "dbsnp", str, None, "dbSNP VCF", "reference"),
    ("--dbsnp-tbi", "dbsnp_tbi", str, None, "dbSNP TBI", "reference"),
    ("--dbsnp-vqsr", "dbsnp_vqsr", str, None, "dbSNP VQSR resource flag", "reference"),
    ("--known-indels", "known_indels", str, None, "Known indels VCF", "reference"),
    ("--known-indels-tbi", "known_indels_tbi", str, None, "Known indels TBI", "reference"),
    ("--known-indels-vqsr", "known_indels_vqsr", str, None, "Known indels VQSR resource flag", "reference"),
    ("--known-snps", "known_snps", str, None, "Known SNPs VCF", "reference"),
    ("--known-snps-tbi", "known_snps_tbi", str, None, "Known SNPs TBI", "reference"),
    ("--known-snps-vqsr", "known_snps_vqsr", str, None, "Known SNPs VQSR resource flag", "reference"),
    ("--germline-resource", "germline_resource", str, None, "Germline resource VCF", "reference"),
    ("--germline-resource-tbi", "germline_resource_tbi", str, None, "Germline resource TBI", "reference"),
    ("--pon", "pon", str, None, "Panel of Normals VCF", "reference"),
    ("--pon-tbi", "pon_tbi", str, None, "PON TBI", "reference"),
    ("--intervals", "intervals", str, None, "Intervals BED/interval_list", "reference"),
    ("--ascat-alleles", "ascat_alleles", str, None, "ASCAT alleles", "reference"),
    ("--ascat-loci", "ascat_loci", str, None, "ASCAT loci", "reference"),
    ("--ascat-loci-gc", "ascat_loci_gc", str, None, "ASCAT loci GC", "reference"),
    ("--ascat-loci-rt", "ascat_loci_rt", str, None, "ASCAT loci replication timing", "reference"),
    ("--chr-dir", "chr_dir", str, None, "Per-chromosome FASTA dir", "reference"),
    ("--mappability", "mappability", str, None, "Mappability track", "reference"),
    ("--msisensor2-models", "msisensor2_models", str, None, "MSIsensor2 models", "reference"),
    ("--msisensorpro-scan", "msisensorpro_scan", str, None, "MSIsensor-pro scan file", "reference"),
    ("--ngscheckmate-bed", "ngscheckmate_bed", str, None, "NGSCheckMate BED", "reference"),
    ("--sentieon-dnascope-model", "sentieon_dnascope_model", str, None, "Sentieon DNAscope model", "reference"),
    ("--bbsplit-fasta-list", "bbsplit_fasta_list", str, None, "BBSplit FASTA list", "reference"),
    ("--bbsplit-index", "bbsplit_index", str, None, "BBSplit pre-built index", "reference"),
    ("--save-reference", "save_reference", bool, False, "Publish reference indices to outdir", "reference"),
    ("--save-bbsplit-reads", "save_bbsplit_reads", bool, False, "Save BBSplit reads", "reference"),
)

_SAREK_PASSTHROUGH_TYPES = {
    dest: ptype for _flag, dest, ptype, _default, _help, _group in _SAREK_PASSTHROUGH_PARAMS
}
# Native nf-core/institutional parameters intentionally kept behind
# ``--extra-param`` still need schema-correct typing. ``help`` is omitted here
# because upstream accepts either boolean or string for that parameter; keeping
# its supplied text verbatim remains schema-valid.
_GENERIC_SAREK_EXTRA_TYPES: dict[str, type] = {
    "config_profile_contact": str,
    "config_profile_description": str,
    "config_profile_name": str,
    "config_profile_url": str,
    "custom_config_base": str,
    "custom_config_version": str,
    "help_full": bool,
    "modules_testdata_base_path": str,
    "monochrome_logs": bool,
    "pipelines_testdata_base_path": str,
    "plaintext_email": bool,
    "show_hidden": bool,
    "test_data_base": str,
    "validate_params": bool,
    "version": bool,
}
_PROTECTED_EXTRA_PARAMS = {"input", "input_restart", "outdir"}


# Reference / index flags --demo must clear before they reach params.yaml.
_DEMO_CLEARED_REFERENCE_FIELDS = tuple(REFERENCE_PATH_PARAMS) + (
    "genome",
    # igenomes_base is already in REFERENCE_PATH_PARAMS; only the VQSR string
    # params below are extras (they're not local paths, so not in that list).
    "dbsnp_vqsr",
    "known_indels_vqsr",
    "known_snps_vqsr",
)


_BANNER = (
    "==============================================\n"
    f"  ClawBio :: {SKILL_NAME} v{SKILL_VERSION}\n"
    "  nf-core/sarek 3.8.1 orchestrator\n"
    "=============================================="
)


def _print(msg: str) -> None:
    """Log a human-readable status line to stdout (progress/traceability)."""
    print(msg, flush=True)


_DEFAULT_TIMEOUT_HOURS = DEFAULT_TIMEOUT_SECONDS / 3600


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


def _check_pipeline_version_supported(args: argparse.Namespace) -> None:
    """Keep execution on the pinned 3.8.1 contract unless explicitly overridden.

    The wrapper's parameter set, profile matrix and output validations are
    hardcoded for nf-core/sarek 3.8.1. Running a different remote tag/commit would
    apply 3.8.1 rules to a pipeline that may differ, so any non-3.8.1
    ``--pipeline-version`` is blocked by default; ``--allow-pipeline-version-override``
    is a recorded, warned opt-in for advanced use. Parity with
    nfcore-scrnaseq/rnaseq.
    """
    requested = str(getattr(args, "pipeline_version", "") or "").strip()
    if requested == DEFAULT_PIPELINE_VERSION:
        return
    if getattr(args, "allow_pipeline_version_override", False):
        print(
            f"WARNING: --pipeline-version {requested!r} differs from the wrapper's pinned "
            f"nf-core/sarek {DEFAULT_PIPELINE_VERSION} contract. Parameter, profile and "
            f"output validations remain {DEFAULT_PIPELINE_VERSION}-specific and may not match.",
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


def _resolve_nextflow_work_dir(args: argparse.Namespace, output_dir: Path) -> Path | str:
    """Resolve the Nextflow work directory: default <output>/upstream/work, or the
    user's --work-dir (a local path, or an object-store URI kept verbatim for cloud
    executors). Parity with nfcore-scrnaseq/rnaseq."""
    raw_work_dir = getattr(args, "work_dir", None)
    if not raw_work_dir:
        return output_dir / "upstream" / "work"
    work_dir = str(raw_work_dir).strip()
    if "://" in work_dir:
        return work_dir
    return Path(work_dir).expanduser().resolve()


# ---------------------------------------------------------------------------
# Argument parser
# ---------------------------------------------------------------------------


def build_parser() -> argparse.ArgumentParser:
    """Construct the wrapper's argparse parser."""
    parser = argparse.ArgumentParser(
        prog="nfcore_sarek_wrapper",
        description="Run nf-core/sarek 3.8.1 through the ClawBio wrapper.",
    )
    parser.add_argument("--input", "-i", default=None, help="Path to a Sarek samplesheet (.csv/.tsv/.json/.yaml/.yml; required unless --demo)")
    parser.add_argument("--output", "-o", required=True, help="Output directory")
    parser.add_argument("--demo", action="store_true", help="Run the upstream `test` profile")
    parser.add_argument("--check", action="store_true", help="Run preflight + samplesheet validation only; no Nextflow")
    parser.add_argument("--resume", action="store_true", help="Validate manifest drift then run with -resume")
    parser.add_argument("--arm", action="store_true", help="Compose arm64 profile (and skip --platform linux/amd64)")
    parser.add_argument("--gpu", action="store_true", help="Compose gpu profile")
    parser.add_argument("--spark-profile", action="store_true", help="Compose spark profile")
    parser.add_argument("--mutect-profile", action="store_true", help="Compose the mutect TEST-DATA profile (only valid under --demo; hardcodes Mutect2 --normal-sample normal)")
    parser.add_argument("--run-downstream", action="store_true", help="Opt-in: write downstream skill handoff")
    parser.add_argument(
        "--downstream-skill",
        default=None,
        choices=list(_DOWNSTREAM_CHOICES),
        help="Choice of downstream ClawBio skill (required with --run-downstream)",
    )
    parser.add_argument(
        "--profile",
        default=DEFAULT_PROFILE,
        help="Nextflow profile string (default: docker). Wrapper composes test/arm64/gpu/spark/mutect on top.",
    )
    parser.add_argument(
        "--nextflow-config",
        "-c",
        "--config",
        dest="nextflow_config",
        action="append",
        metavar="CONFIG",
        default=None,
        help="Extra Nextflow -c config file(s) (aliases: -c/--config). Can be repeated.",
    )
    parser.add_argument("--pipeline-version", default=DEFAULT_PIPELINE_VERSION, help="Pinned nf-core/sarek tag/ref")
    parser.add_argument(
        "--allow-pipeline-version-override",
        action="store_true",
        help=(
            f"Run a --pipeline-version other than the pinned {DEFAULT_PIPELINE_VERSION} "
            "contract (recorded, warned). Validations stay "
            f"{DEFAULT_PIPELINE_VERSION}-specific, so this is at your own risk."
        ),
    )
    parser.add_argument("--pipeline-local", default=None, help="Local nf-core/sarek checkout")
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
            "runtime warning lists the paths fetched over the network. (The iGenomes "
            "mirror base and object-store --work-dir are not gated.)"
        ),
    )
    parser.add_argument("--params-file", default=None, help="Sarek-native --params-file (advanced)")
    parser.add_argument("--no-banner", action="store_true", help="Suppress console banner")
    parser.add_argument("-v", "--verbose", action="store_true", help="Verbose logging")
    parser.add_argument(
        "--timeout-hours",
        type=_timeout_hours,
        default=_DEFAULT_TIMEOUT_HOURS,
        help=(
            f"Wall-clock cap for the Nextflow run in hours (default {_DEFAULT_TIMEOUT_HOURS:g}). "
            "Use 0 to disable the cap for long HPC/cloud runs whose walltime is enforced "
            "by the scheduler."
        ),
    )
    parser.add_argument(
        "--extra-param",
        action="append",
        metavar="KEY=VALUE",
        default=None,
        help="Pass an arbitrary Sarek param (repeatable). e.g. --extra-param tools=haplotypecaller",
    )

    # ----- passthrough param groups (Sarek schema) -----
    groups: dict[str, argparse._ArgumentGroup] = {}

    def _group(name: str, title: str) -> argparse._ArgumentGroup:
        if name not in groups:
            groups[name] = parser.add_argument_group(title)
        return groups[name]

    _GROUP_TITLES = {
        "main": "Main pipeline options",
        "io": "I/O & metadata",
        "trim": "FASTQ preprocessing",
        "umi": "UMI handling",
        "preprocessing": "Preprocessing",
        "varcall": "Variant calling",
        "postvar": "Post-variant calling",
        "annotation": "Annotation",
        "reference": "Reference & indices",
    }

    for cli_flag, dest, ptype, default, help_text, group_key in _SAREK_PASSTHROUGH_PARAMS:
        g = _group(group_key, _GROUP_TITLES.get(group_key, group_key))
        # Accept the nf-core-native snake_case spelling (``--<dest>``) alongside the
        # hyphenated wrapper flag, so a command copied from the upstream nf-core
        # docs is not rejected. ``dest`` is exactly the nf-core parameter name.
        native = f"--{dest}"
        option_strings = [cli_flag] if native == cli_flag else [cli_flag, native]
        if ptype is bool:
            g.add_argument(*option_strings, dest=dest, action="store_true", default=default, help=help_text)
        else:
            g.add_argument(*option_strings, dest=dest, type=ptype, default=default, help=help_text)

    return parser


# ---------------------------------------------------------------------------
# Entry point
# ---------------------------------------------------------------------------


def main(argv: list[str] | None = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)
    if not args.no_banner:
        _print(_BANNER)
    output_dir = Path(args.output).expanduser().resolve()
    try:
        _check_pipeline_version_supported(args)
        return _run_wrapper(args, output_dir)
    except SkillError as exc:
        return _handle_skill_error(output_dir, exc, verbose=args.verbose)
    except KeyboardInterrupt:
        _print("[abort] Interrupted by user.")
        return 130
    except Exception as exc:
        return _handle_unexpected_error(output_dir, exc, verbose=args.verbose)


def _run_wrapper(args: argparse.Namespace, output_dir: Path) -> int:
    # Materialise config-provided step/tools before input validation. In
    # particular, downstream Sarek steps may omit input and retrieve their CSV
    # handoff, whereas mapping may not.
    _merge_extra_params(args)
    _backfill_samplesheet_args_from_params_file(args)
    # Sync profile modifier tokens (test/arm64/gpu/...) onto the boolean flags
    # BEFORE enforcing the input requirement: a bare `--profile test` is a demo
    # run and must not trip the MISSING_INPUT check.
    _sync_profile_flags(args)
    _validate_wrapper_flags(args)

    # --- demo reference cleanup (after validation; clears input/refs) ---
    _apply_demo_overrides(args)

    composed_profile = compose_profile(
        args.profile or DEFAULT_PROFILE,
        demo=bool(args.demo),
        arm=bool(args.arm),
        gpu=bool(args.gpu),
        spark=bool(args.spark_profile),
        mutect_profile=bool(args.mutect_profile),
    )
    args.profile = composed_profile  # propagate to preflight + provenance

    # Custom --fasta without --genome → disable iGenomes so Sarek does not load its
    # default GATK.GRCh38 bundle and validate ~20 remote paths as local files. Must
    # run after demo-override cleanup and profile composition (see helper docstring).
    igenomes_notice = _auto_disable_igenomes_for_custom_fasta(args)
    if igenomes_notice:
        print(f"WARNING: {igenomes_notice}", file=sys.stderr)

    output_dir.mkdir(parents=True, exist_ok=True)

    # --- pipeline source resolution ---
    pipeline_source = resolve_pipeline_source(
        requested_version=args.pipeline_version,
        local_pipeline_dir=(
            Path(args.pipeline_local).expanduser().resolve()
            if args.pipeline_local
            else None
        ),
    )

    # --- samplesheet preparation ---
    samplesheet_report, normalized_csv = _prepare_samplesheet(args, output_dir)

    # --- preflight ---
    repo_root = _detect_repo_root()
    params_for_preflight = _build_params_for_preflight(args, composed_profile=composed_profile)
    resume_manifest = (
        load_manifest(output_dir / "reproducibility") if args.resume else None
    )
    if args.resume and resume_manifest is None:
        _print("[preflight] WARNING: --resume requested but no prior manifest.json found; continuing as a fresh run.")
    elif args.resume:
        # Compare the complete emitted parameter set and local/remote reference
        # fingerprints before allowing reuse of any cached Nextflow work.
        drift_params = build_effective_params(
            args,
            normalized_samplesheet=normalized_csv,
            output_dir=output_dir,
        )
        params_for_preflight["params_checksum"] = compute_params_checksum(drift_params)
        params_for_preflight["reference_checksums"] = compute_reference_checksums(drift_params)
    preflight_result = run_preflight(
        params=params_for_preflight,
        samplesheet=samplesheet_report,
        pipeline_source=pipeline_source,
        output_dir=output_dir,
        repo_root=repo_root,
        resume=bool(args.resume),
        resume_manifest=resume_manifest,
        allow_remote_inputs=bool(getattr(args, "allow_remote_inputs", False)),
    )
    # Record the auto-igenomes-ignore notice in the report's warning list too, so the
    # behaviour change is provenance-visible (not only on stderr at emit time).
    if igenomes_notice:
        preflight_result.warnings.insert(0, igenomes_notice)
    _print(f"[preflight] passed (warnings: {len(preflight_result.warnings)})")
    if args.verbose:
        for w in preflight_result.warnings:
            _print(f"  · warning: {w}")
        for n in preflight_result.notes:
            _print(f"  · note: {n}")

    # --- check mode: stop here ---
    if args.check:
        return _write_check_report(
            output_dir=output_dir,
            args=args,
            samplesheet_report=samplesheet_report,
            preflight_result=preflight_result,
            pipeline_source=pipeline_source,
        )

    # --- params + command ---
    params = build_effective_params(
        args,
        normalized_samplesheet=normalized_csv,
        output_dir=output_dir,
    )
    params_path = write_params_yaml(params, output_dir=output_dir)

    # On macOS+docker: the macOS workarounds config (also carries resourceLimits).
    # On non-macOS+docker (real runs): a host-scaled resourceLimits cap so the local
    # executor does not abort when a process's production requirement exceeds this
    # host. Exactly one of the two applies per host; both are self-gating.
    auto_cfg = _write_macos_docker_config(output_dir, args=args)
    if auto_cfg is None:
        auto_cfg = _write_resource_limits_config(output_dir, args=args)
    extra_configs = _resolve_extra_configs(args, macos_cfg=auto_cfg)
    work_dir = _resolve_nextflow_work_dir(args, output_dir)

    nextflow_command, command_str = build_nextflow_command(
        pipeline_source=pipeline_source,
        profile=args.profile,
        params_path=params_path,
        resume=bool(args.resume),
        work_dir=work_dir,
        extra_configs=extra_configs,
        demo=bool(args.demo),
        arm=bool(args.arm),
        gpu=bool(args.gpu),
        spark=bool(args.spark_profile),
        mutect_profile=bool(args.mutect_profile),
    )

    # --- execute ---
    _print(f"[execute] launching Nextflow ({pipeline_source['source_kind']} → {pipeline_source['resolved_version']})")
    started = time.monotonic()
    execute_nextflow(
        nextflow_command,
        cwd=output_dir,
        output_dir=output_dir,
        timeout_seconds=_resolve_timeout_seconds(args),
    )
    elapsed = round(time.monotonic() - started, 3)
    _print(f"[execute] completed in {elapsed:.1f}s")

    # --- outputs parsing ---
    # Nextflow runs with cwd=output_dir and params.outdir="upstream/results", so the
    # pipeline writes its output tree under <output_dir>/upstream/results/.
    _print("[outputs] parsing pipeline outputs")
    results_dir = output_dir / "upstream" / "results"
    analysis_mode = str(samplesheet_report.get("analysis_mode") or "germline")
    build_only_index = bool(params.get("build_only_index", False))
    outputs_report = parse_outputs(
        results_dir,
        step=str(params.get("step") or "mapping"),
        # Upstream runs with an empty samplesheet channel when only resources
        # are built, so no per-sample caller output is expected.
        tools=[] if build_only_index else _tool_list(params.get("tools")),
        skip_tools=_tool_list(params.get("skip_tools")),
        analysis_mode=analysis_mode,
        aligner=str(params.get("aligner") or DEFAULT_ALIGNER),
        save_mapped=bool(params.get("save_mapped", False)),
        save_output_as_bam=bool(params.get("save_output_as_bam", False)),
        # PREPARE_GENOME / PREPARE_INTERVALS publish their reference outputs
        # when building indices even without an explicit --save-reference.
        save_reference=bool(params.get("save_reference", False) or params.get("build_only_index", False)),
        validate=not build_only_index,
    )

    # --- reporting ---
    _print("[report] writing report.md, result.json, commands.sh")
    reporting_artifacts = write_reports(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        params=params,
        samplesheet_report=samplesheet_report,
        pipeline_source=pipeline_source,
        outputs_report=outputs_report,
        nextflow_command=nextflow_command,
        pipeline_source_kind=str(pipeline_source["source_kind"]),
        resume_used=bool(args.resume),
        arm=bool(args.arm),
        profile=composed_profile,
        java_version=getattr(preflight_result, "java_version", "") or None,
        nextflow_version=getattr(preflight_result, "nextflow_version", "") or None,
        elapsed_seconds=elapsed,
    )

    # --- provenance ---
    _print("[provenance] writing reproducibility bundle")
    commands_sh_src = reporting_artifacts.commands_sh
    write_provenance_bundle(
        output_dir=output_dir,
        skill_dir=_SKILL_DIR,
        samplesheet_csv_src=normalized_csv,
        params_yaml_src=params_path,
        commands_sh_src=commands_sh_src,
        params=params,
        samplesheet_report=samplesheet_report,
        pipeline_source=pipeline_source,
        outputs_report=outputs_report,
        resume_used=bool(args.resume),
        arm=bool(args.arm),
        spark=bool(args.spark_profile),
        gpu=bool(args.gpu),
        wes=bool(params.get("wes", False)),
        profile=args.profile,
        java_version=_detect_java_version(),
        nextflow_version=_detect_nextflow_version(),
    )

    # --- downstream handoff ---
    if args.run_downstream:
        _emit_downstream_handoff(args, outputs_report=outputs_report, output_dir=output_dir)

    _print(f"[done] Wrapper completed successfully. Output: {output_dir}")
    return 0


# ---------------------------------------------------------------------------
# Helpers — validation, overrides, profile composition
# ---------------------------------------------------------------------------


def _validate_wrapper_flags(args: argparse.Namespace) -> None:
    """Cheap CLI-level invariants that don't need preflight."""
    # The `mutect` profile only pulls in conf/test_mutect2.config, whose sole
    # effect is forcing `--normal-sample normal` on MUTECT2_PAIRED — a value taken
    # from when the upstream test data was generated. On a real somatic paired run
    # whose normal sample is not literally named "normal" it silently mislabels the
    # normal and corrupts/fails the paired calls. Allow it only under --demo (the
    # upstream test profile, where the normal IS "normal"); reject it otherwise.
    if getattr(args, "mutect_profile", False) and not args.demo:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=(
                "The `mutect` profile is only valid for the upstream test dataset; it "
                "hardcodes Mutect2's `--normal-sample normal`, which mislabels the normal "
                "sample on any real paired run."
            ),
            fix=(
                "Remove `--mutect-profile` and any `mutect` token from --profile for real "
                "runs (it is intended only under --demo / the upstream test profile)."
            ),
            details={"profile": getattr(args, "profile", None)},
        )
    params_file_values = _load_user_params_file(getattr(args, "params_file", None))
    if not args.input and params_file_values.get("input") not in (None, "", False, "false", "False"):
        args.input = str(params_file_values["input"])
    if not getattr(args, "input_restart", None) and params_file_values.get("input_restart") not in (None, "", False):
        args.input_restart = str(params_file_values["input_restart"])

    # Sarek uses the sentinel `--input false` for build-only runs. Cache download
    # is often layered onto that mode, but --download_cache alone still consumes input.
    if isinstance(args.input, str) and args.input.strip().lower() in ("false", ""):
        args.input = None
    # Although input_restart is present in Sarek's schema, v3.8.1 overwrites it
    # internally with retrieveInput(...). The wrapper therefore validates a
    # user-supplied restart sheet and submits it as regular input.
    # The official cache-only example combines --build_only_index,
    # --download_cache and --input false; build_only_index is the switch that
    # creates an empty samplesheet channel in Sarek.
    no_input_mode = bool(getattr(args, "build_only_index", False))
    effective_step = (
        getattr(args, "step", None)
        or params_file_values.get("step")
        or "mapping"
    )
    auto_restart = effective_step in _RESTART_CSV_BY_STEP
    if (
        not args.demo
        and not args.input
        and not getattr(args, "input_restart", None)
        and not no_input_mode
        and not auto_restart
    ):
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.MISSING_INPUT,
            message="The mapping step requires an input samplesheet.",
            fix="Provide --input <samplesheet.csv>, pass --demo, or start a downstream step where Sarek can retrieve its prior CSV handoff.",
            details={"step": effective_step},
        )
    if args.run_downstream and not args.downstream_skill:
        raise SkillError(
            stage="validation",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message="--run-downstream requires --downstream-skill <name>.",
            fix=f"Pick one of: {', '.join(_DOWNSTREAM_CHOICES)}.",
            details={"choices": list(_DOWNSTREAM_CHOICES)},
        )


def _apply_demo_overrides(args: argparse.Namespace) -> None:
    """When --demo is set, clear reference flags + force test-friendly defaults.

    Clears --input and every reference/index flag (the upstream `test` profile ships its
    own samplesheet and matched references, and a partial override would desynchronise
    them). Does NOT touch --resume — see the note at the bottom of this function.
    """
    if not args.demo:
        return
    cleared: list[str] = []
    for field in _DEMO_CLEARED_REFERENCE_FIELDS:
        if getattr(args, field, None):
            cleared.append(field)
            setattr(args, field, None)
    # Force genome unset so the test profile owns it. Also clear a user-supplied
    # --igenomes-ignore: test.config defines genome + igenomes_base, and
    # igenomes_ignore=True would override that and break the demo run (it is a
    # boolean flag written by _add_flags, so leaving it set leaks into params.yaml).
    args.genome = None
    if getattr(args, "igenomes_ignore", False):
        cleared.append("igenomes_ignore")
        args.igenomes_ignore = False
    # The --extra-param escape hatch bypasses the per-arg clearing above
    # (params_builder writes every _extras entry verbatim), so prune reference
    # keys from the parsed extras too.
    extras = getattr(args, "_extras", None)
    if isinstance(extras, dict):
        for key in (*_DEMO_CLEARED_REFERENCE_FIELDS, "igenomes_ignore"):
            if extras.pop(key, None) is not None and key not in cleared:
                cleared.append(key)
    # Demo coercions override user-supplied flags, so they are advisories the user
    # should notice → stderr WARNING (parity with nfcore-scrnaseq/rnaseq), not a
    # stdout progress line.
    if cleared:
        flags = ", ".join("--" + f.replace("_", "-") for f in cleared)
        print(
            f"WARNING: --demo ignores reference flags ({flags}); "
            "the upstream test profile provides its own bundled references.",
            file=sys.stderr,
        )
    if args.input:
        print(
            "WARNING: --demo ignores --input; the upstream test profile provides its own samplesheet.",
            file=sys.stderr,
        )
        args.input = None
    # --resume is deliberately NOT cleared here. Nextflow's -resume is orthogonal to
    # -profile test (nf-core documents no incompatibility), and a demo keeps everything a
    # resume needs inside its own output dir. Clearing it left an interrupted demo — a
    # common outcome, since the test profile streams its data over the network — in a dead
    # end: the half-populated output dir makes preflight raise OUTPUT_DIR_NOT_EMPTY, whose
    # fix text points at --resume, which this override then silently threw away.
    # Demo/real drift stays blocked: params["profile"] records the *composed* profile
    # (carrying the `test` token) and is compared via _RESUME_TRACKED_FIELDS.


def _auto_disable_igenomes_for_custom_fasta(args: argparse.Namespace) -> str | None:
    """Disable iGenomes automatically for a custom-FASTA run with no --genome.

    nf-core/sarek 3.8.1 defaults ``genome = 'GATK.GRCh38'`` in ``nextflow.config``.
    The wrapper deliberately omits upstream-default parameters from ``params.yaml``,
    so a custom-reference run that supplies ``--fasta`` but leaves ``--genome`` unset
    would still let Sarek load the full GATK.GRCh38 iGenomes configuration and try
    to validate its ~20 ``s3://ngi-igenomes/...`` reference paths as local files —
    failing with a wall of "does not exist" errors that never mention
    ``--igenomes-ignore``. Setting ``igenomes_ignore=true`` matches the nf-core
    guidance for custom references ("``--igenomes_ignore`` must be set ... The
    ``fasta`` file is the only required input file") and mirrors the equivalent
    automatic behaviour already in nfcore-scrnaseq-wrapper.

    Returns an explanatory notice when it changed anything, else ``None``. Must run
    AFTER ``_apply_demo_overrides`` (which clears references for demo) and AFTER the
    profile is composed, so the upstream ``test`` profile — which owns its own
    genome/igenomes_base — is never touched.
    """
    if getattr(args, "demo", False) or getattr(args, "_noinput", False):
        return None
    profile_parts = {p.strip() for p in str(getattr(args, "profile", "") or "").split(",") if p.strip()}
    if any(p == "test" or p.startswith("test_") for p in profile_parts):
        return None
    if not getattr(args, "fasta", None):
        return None
    if getattr(args, "genome", None):
        # Explicit --genome alongside --fasta is a documented override; respect it.
        return None
    if getattr(args, "igenomes_ignore", False):
        return None
    args.igenomes_ignore = True
    return (
        "--fasta was provided without --genome: iGenomes is disabled automatically "
        "(igenomes_ignore=true) so Sarek does not fall back to its default "
        "GATK.GRCh38 iGenomes reference bundle. Pass --genome to use iGenomes, or "
        "--igenomes-ignore explicitly to silence this notice."
    )


def _sync_profile_flags(args: argparse.Namespace) -> None:
    """If the user passed modifier tokens in --profile, mirror them into the boolean flags."""
    parts = {p.strip() for p in (args.profile or "").split(",") if p.strip()}
    if "arm64" in parts:
        args.arm = True
    if "gpu" in parts:
        args.gpu = True
    if "spark" in parts:
        args.spark_profile = True
    if "mutect" in parts:
        args.mutect_profile = True
    # Treat ANY upstream test/test_full* profile as demo so reference flags get
    # cleared — they all carry their own datasets/refs and a user-supplied
    # --genome / --fasta would conflict (audit pass-2 M4).
    if any(p == "test" or p.startswith("test_") for p in parts):
        args.demo = True


def _merge_extra_params(args: argparse.Namespace) -> None:
    """Merge native Sarek extras while keeping one effective value per param."""
    raw = getattr(args, "extra_param", None) or []
    extras: dict[str, Any] = {}
    for entry in raw:
        if "=" not in entry:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"--extra-param must be key=value, got {entry!r}",
                fix="Use the form --extra-param tools=haplotypecaller",
                details={"entry": entry},
            )
        key, value = entry.split("=", 1)
        key = key.strip()
        value = value.strip()
        if key in _PROTECTED_EXTRA_PARAMS:
            raise SkillError(
                stage="validation",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"--extra-param cannot override wrapper-managed parameter '{key}'.",
                fix=(
                    "Use --input/--input-restart for samplesheets and --output for the "
                    "wrapper output root; these values must be normalized and tracked."
                ),
                details={"key": key},
            )
        effective_value: Any = value
        exposed_type = _SAREK_PASSTHROUGH_TYPES.get(key) or _GENERIC_SAREK_EXTRA_TYPES.get(key)
        if exposed_type is bool:
            if value.lower() not in {"true", "false"}:
                raise SkillError(
                    stage="validation",
                    error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                    message=f"--extra-param {key}=... expects a boolean value.",
                    fix=f"Use --extra-param {key}=true or --extra-param {key}=false.",
                    details={"key": key, "value": value},
                )
            effective_value = value.lower() == "true"
        elif exposed_type in {int, float}:
            try:
                effective_value = exposed_type(value)
            except ValueError as exc:
                raise SkillError(
                    stage="validation",
                    error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                    message=f"--extra-param {key}=... has an invalid numeric value.",
                    fix=f"Provide a valid {exposed_type.__name__} value for {key}.",
                    details={"key": key, "value": value},
                ) from exc

        extras[key] = effective_value
        # Exposed Sarek params affect samplesheet/preflight decisions before
        # params.yaml is built, so apply the same final override here. Unknown
        # native params are retained only for pass-through and can never toggle
        # wrapper-only controls such as demo/check/resume.
        if exposed_type is not None:
            setattr(args, key, effective_value)
    args._extras = extras  # for provenance debug if needed


def _backfill_samplesheet_args_from_params_file(args: argparse.Namespace) -> None:
    """Surface ``step``/``tools`` from a native ``--params-file`` onto ``args``.

    The samplesheet is validated per ``--step`` *before* the full params
    snapshot is assembled. When the user supplies step/tools through a
    ``--params-file`` (instead of CLI flags), those values must reach the
    samplesheet validator too — otherwise an ``annotate`` or restart samplesheet
    is wrongly validated as ``mapping`` with empty tools. CLI flags always win.
    """
    pf = _load_user_params_file(getattr(args, "params_file", None))
    if not pf:
        return
    if not getattr(args, "step", None) and pf.get("step"):
        args.step = str(pf["step"]).strip()
    if getattr(args, "tools", None) in (None, "") and pf.get("tools"):
        tools: Any = pf["tools"]
        args.tools = tools if isinstance(tools, str) else ",".join(str(t) for t in cast("list[Any]", tools))


def _build_params_for_preflight(args: argparse.Namespace, *, composed_profile: str) -> dict[str, Any]:
    """Snapshot a minimal params dict for the preflight check.

    This is *not* the final params.yaml — that comes from params_builder.
    We only need the keys preflight rules inspect.
    """
    params: dict[str, Any] = _load_user_params_file(getattr(args, "params_file", None))
    params["profile"] = composed_profile
    params.setdefault("step", "mapping")

    scalar_fields = (
        "aligner", "use_gatk_spark", "umi_read_structure", "email", "email_on_fail",
        "consensus_min_count", "bcftools_filter_criteria", "snpsift_databases",
        "input", "outdir_cache", "umi_location", "umi_length",
        "group_by_umi_strategy", "vep_out_format", "ascat_genome", "publish_dir_mode",
        "ascat_purity", "ascat_ploidy", "vep_custom_args",
        "sentieon_haplotyper_emit_mode", "sentieon_dnascope_emit_mode",
        "sentieon_dnascope_pcr_indel_model", "max_multiqc_email_size",
        "phenotypes_include_types",
        "bcftools_header_lines", "condel_config", "mastermind_file",
        "spliceai_snv", "spliceai_snv_tbi", "spliceai_indel", "spliceai_indel_tbi",
        "dbnsfp", "dbnsfp_tbi",
    )
    bool_fields = (
        "wes", "no_intervals", "joint_germline", "joint_mutect2",
        "only_paired_variant_calling", "normalize_vcfs",
        "build_only_index", "download_cache", "snv_consensus_calling",
        "filter_vcfs", "concatenate_vcfs", "save_mapped", "save_output_as_bam",
        "umi_in_read_header", "vep_dbnsfp", "vep_loftee",
        "vep_condel", "vep_mastermind", "vep_spliceai", "vep_phenotypes",
    )

    if getattr(args, "step", None):
        params["step"] = args.step
    if getattr(args, "tools", None) is not None:
        params["tools"] = _tool_list(args.tools)
    else:
        params["tools"] = _tool_list(params.get("tools"))
    if getattr(args, "skip_tools", None) is not None:
        params["skip_tools"] = _tool_list(args.skip_tools)
    else:
        params["skip_tools"] = _tool_list(params.get("skip_tools"))

    for field in scalar_fields:
        value = getattr(args, field, None)
        if value not in (None, ""):
            params[field] = value
    for field in bool_fields:
        if bool(getattr(args, field, False)):
            params[field] = True
        else:
            params.setdefault(field, False)

    for name in REFERENCE_PATH_PARAMS:
        v = getattr(args, name, None)
        if v:
            params[name] = v
    # iGenomes/genome
    if getattr(args, "genome", None):
        params["genome"] = args.genome
    if getattr(args, "igenomes_base", None):
        params["igenomes_base"] = args.igenomes_base
    if getattr(args, "igenomes_ignore", False):
        params["igenomes_ignore"] = True
    if getattr(args, "snpeff_db", None):
        params["snpeff_db"] = args.snpeff_db
    if getattr(args, "vep_cache_version", None):
        params["vep_cache_version"] = args.vep_cache_version
    if getattr(args, "vep_genome", None):
        params["vep_genome"] = args.vep_genome
    if getattr(args, "vep_species", None):
        params["vep_species"] = args.vep_species
    # `--extra-param` is the explicit final override of a reusable
    # `--params-file`; preflight must inspect the same effective values that
    # params_builder will emit to params.yaml.
    for key, value in (getattr(args, "_extras", None) or {}).items():
        # Values were schema-typed by `_merge_extra_params`; preserve strings
        # such as custom_config_version="false" rather than silently turning
        # an official string parameter into a boolean.
        params[key] = value
    return params


def _tool_list(raw: Any) -> list[str]:
    if not raw:
        return []
    if isinstance(raw, list):
        return [str(x).strip() for x in raw if str(x).strip()]
    return [t.strip() for t in str(raw).split(",") if t.strip()]


# ---------------------------------------------------------------------------
# Samplesheet preparation
# ---------------------------------------------------------------------------


def _prepare_samplesheet(
    args: argparse.Namespace, output_dir: Path
) -> tuple[dict[str, Any], Path]:
    """Validate the user's samplesheet (or write a demo stub) and return (report, csv_path)."""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    if args.demo:
        normalized_csv = repro_dir / "samplesheet.demo.csv"
        write_text_lf(normalized_csv, "patient,sample,lane,fastq_1,fastq_2\n")
        return {
            "normalized_path": normalized_csv,
            "step": getattr(args, "step", None) or "mapping",
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
            "rows_by_patient": {},
        }, normalized_csv
    if not args.input and getattr(args, "input_restart", None):
        args.input = args.input_restart

    if not args.input and getattr(args, "build_only_index", False):
        # Build-only run: no samplesheet is consumed; it may also download caches.
        normalized_csv = repro_dir / "samplesheet.noinput.csv"
        write_text_lf(normalized_csv, "patient,sample,lane,fastq_1,fastq_2\n")
        return {
            "normalized_path": normalized_csv,
            "step": getattr(args, "step", None) or "mapping",
            "sample_count": 0, "patient_count": 0,
            "sample_names": [], "patient_names": [],
            "fastq_paths": [], "bam_paths": [], "cram_paths": [],
            "vcf_paths": [], "spring_paths": [], "tables": [],
            "unknown_columns": [], "sex_counts": {}, "status_counts": {},
            "pairings": [], "analysis_mode": "germline",
            "rows_by_patient": {},
        }, normalized_csv

    if not args.input:
        step = getattr(args, "step", None) or "mapping"
        restart_name = _RESTART_CSV_BY_STEP.get(step)
        restart_path = (
            output_dir / "upstream" / "results" / "csv" / restart_name
            if restart_name
            else None
        )
        if restart_path is None or not restart_path.exists():
            raise SkillError(
                stage="samplesheet",
                error_code=ErrorCode.MISSING_INPUT,
                message=f"No input was provided for step '{step}' and its Sarek handoff CSV was not found.",
                fix=(
                    f"Provide --input, or retain the prior Sarek output file "
                    f"upstream/results/csv/{restart_name or '<handoff>.csv'} under the wrapper output directory."
                ),
                details={"step": step, "expected_handoff": str(restart_path) if restart_path else None},
            )
        args.input = restart_path.as_posix()
        args._auto_restart_input = True

    normalized_csv = repro_dir / "samplesheet.valid.csv"
    input_path = _materialize_samplesheet_input(str(args.input), repro_dir=repro_dir)
    report = validate_and_normalize_samplesheet(
        input_path,
        normalized_csv,
        step=getattr(args, "step", None) or "mapping",
        tools=_tool_list(getattr(args, "tools", None)),
    )
    return report, normalized_csv


def _materialize_samplesheet_input(value: str, *, repro_dir: Path) -> Path:
    """Return a local samplesheet path, staging Nextflow-supported URIs first."""
    if "://" not in value:
        return Path(value).expanduser().resolve()

    parsed = urlsplit(value)
    suffix = Path(parsed.path).suffix.lower()
    if suffix not in {".csv", ".tsv", ".json", ".yaml", ".yml"}:
        raise SkillError(
            stage="samplesheet",
            error_code=ErrorCode.INVALID_SAMPLESHEET,
            message="Remote samplesheet URL does not end in a supported Sarek samplesheet extension.",
            fix="Provide a .csv, .tsv, .json, .yaml or .yml samplesheet URI.",
            details={"input": value},
        )
    staged_path = repro_dir / f"samplesheet.remote{suffix}"
    staged_path.unlink(missing_ok=True)
    try:
        proc = subprocess.run(
            ["nextflow", "fs", "cp", value, staged_path.as_posix()],
            capture_output=True,
            text=True,
            check=False,
            timeout=120,
        )
    except FileNotFoundError as exc:
        raise SkillError(
            stage="samplesheet",
            error_code=ErrorCode.NEXTFLOW_NOT_FOUND,
            message="Nextflow is required to stage a remote samplesheet.",
            fix="Install Nextflow >=25.10.2 or provide a local samplesheet file.",
            details={"input": value},
        ) from exc
    except subprocess.TimeoutExpired as exc:
        raise SkillError(
            stage="samplesheet",
            error_code=ErrorCode.MISSING_INPUT,
            message="Timed out while retrieving remote samplesheet through Nextflow.",
            fix="Check URI access and credentials, or download the samplesheet locally and pass --input to that file.",
            details={"input": value},
        ) from exc
    if proc.returncode != 0 or not staged_path.exists():
        raise SkillError(
            stage="samplesheet",
            error_code=ErrorCode.MISSING_INPUT,
            message="Nextflow could not retrieve the remote samplesheet.",
            fix="Check the URI, network access and any cloud credentials required by Nextflow.",
            details={"input": value, "stderr": proc.stderr.strip()[-1000:]},
        )
    return staged_path


# ---------------------------------------------------------------------------
# Configs + system probes
# ---------------------------------------------------------------------------


def _host_memory_gb() -> int:
    """Best-effort total physical RAM in GiB (portable; falls back to 16)."""
    try:
        return max(1, int(os.sysconf("SC_PHYS_PAGES") * os.sysconf("SC_PAGE_SIZE") / (1024 ** 3)))
    except (ValueError, OSError, AttributeError):
        return 16


def _docker_vm_memory_gb() -> int | None:
    """Docker runtime's total memory in whole GB, or None when unknown.

    On native Linux this reports the host's cgroup memory (≈ physical RAM); it is
    the true ceiling before the OOM-killer intervenes. Best-effort: any failure
    (docker absent, daemon down, unparseable output) returns None.
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


_RESOURCE_LIMITS_MEMORY_HEADROOM_GB = 4
_RESOURCE_LIMITS_MEMORY_FLOOR_GB = 8


def _resource_limits_memory_gb() -> int:
    """Per-process memory ceiling (GB) for the non-macOS resourceLimits config.

    nf-core's guidance is that ``resourceLimits`` should match the machine's maximum
    (https://nf-co.re/docs/running/configuration/nextflow-for-your-system). Native
    Linux Docker exposes the full host RAM, so — unlike the macOS Docker Desktop VM —
    there is no 15 GB ceiling here. We take the smaller of the detected physical RAM
    and (when known) the Docker runtime's MemTotal, then reserve a few GB of headroom
    for the OS, the Nextflow JVM and the Docker daemon.
    """
    signals = [gb for gb in (_host_memory_gb(), _docker_vm_memory_gb()) if gb]
    budget = (min(signals) - _RESOURCE_LIMITS_MEMORY_HEADROOM_GB) if signals else 12
    return max(_RESOURCE_LIMITS_MEMORY_FLOOR_GB, budget)


def _write_resource_limits_config(output_dir: Path, *, args: argparse.Namespace) -> Path | None:
    """Host-scaled ``resourceLimits`` config for non-macOS docker runs.

    Applied to the live run so a real (non-demo) run whose production process
    requirements exceed this host does not abort with "Process requirement exceeds
    available memory". Written inside ``reproducibility/`` so it ships with the
    portable bundle: ``commands.sh`` re-applies it on non-macOS replay hosts, letting a
    from-scratch reproduction on the generating machine succeed. Carries only the
    portable resourceLimits block — none of the macOS-only workarounds. Returns None
    when not applicable (macOS is handled by ``_write_macos_docker_config``; ``--demo``
    relies on -profile test's own limits).
    """
    if sys.platform == "darwin":
        return None
    profile_parts = {p.strip() for p in (args.profile or "").split(",") if p.strip()}
    if "docker" not in profile_parts or bool(getattr(args, "demo", False)):
        return None
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


def _write_macos_docker_config(output_dir: Path, *, args: argparse.Namespace) -> Path | None:
    """On macOS with docker backend, add platform support without losing profile flags."""
    if sys.platform != "darwin":
        return None
    profile_parts = {p.strip() for p in (args.profile or "").split(",") if p.strip()}
    if "docker" not in profile_parts:
        return None
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    config_path = repro_dir / "macos_docker.config"
    arm = bool(args.arm)
    cpus = max(4, os.cpu_count() or 4)
    # Scale the memory cap to the host (with ~4 GiB headroom for the OS) instead
    # of a hardcoded 15 GiB: on a 16 GiB Mac a 15 GiB cap leaves nothing for the
    # host and the Docker VM, causing memory pressure / container failures.
    mem_gb = min(15, max(4, _host_memory_gb() - 4))
    docker_block = ""
    if not arm:
        # The official docker/gpu/spark profiles set mutually overriding
        # docker.runOptions. Reproduce the final selected option and add only
        # the macOS architecture requirement; a bare overwrite would drop
        # UID/GID mapping or GPU access.
        if bool(getattr(args, "spark_profile", False)):
            run_options = "--platform=linux/amd64"
        elif bool(getattr(args, "gpu", False)):
            run_options = "-u $(id -u):$(id -g) --gpus all --platform=linux/amd64"
        else:
            run_options = "-u $(id -u):$(id -g) --platform=linux/amd64"
        docker_block = f"docker {{\n    runOptions = '{run_options}'\n}}\n"
    write_text_lf(
        config_path,
        "// macOS Docker compatibility for nf-core/sarek.\n"
        "process {\n"
        "    stageInMode = 'copy'\n"
        "    resourceLimits = [\n"
        f"        cpus: {cpus},\n"
        f"        memory: '{mem_gb}.GB',\n"
        "        time: '24.h'\n"
        "    ]\n"
        + "}\n"
        + docker_block,
    )
    return config_path


def _resolve_extra_configs(
    args: argparse.Namespace, *, macos_cfg: Path | None
) -> list[Path]:
    configs: list[Path] = []
    if macos_cfg is not None:
        configs.append(macos_cfg)
    for cfg in getattr(args, "nextflow_config", None) or []:
        configs.append(Path(cfg).expanduser().resolve())
    return configs


def _detect_repo_root() -> Path:
    """Find the ClawBio repo root (where clawbio.py lives), or fall back to skill grandparent."""
    here = _SKILL_DIR
    for parent in (here, *here.parents):
        if (parent / "clawbio.py").exists():
            return parent
    return _SKILL_DIR.parent.parent


def _detect_java_version() -> str | None:
    try:
        out = subprocess.run(
            ["java", "-version"], capture_output=True, text=True, timeout=10
        )
    except (FileNotFoundError, OSError, subprocess.TimeoutExpired):
        return None
    text = (out.stdout or "") + (out.stderr or "")
    return text.strip().splitlines()[0] if text.strip() else None


def _detect_nextflow_version() -> str | None:
    try:
        out = subprocess.run(
            ["nextflow", "-v"], capture_output=True, text=True, timeout=10
        )
    except (FileNotFoundError, OSError, subprocess.TimeoutExpired):
        return None
    text = (out.stdout or out.stderr or "").strip()
    return text or None


# ---------------------------------------------------------------------------
# Check mode
# ---------------------------------------------------------------------------


def _write_check_report(
    *,
    output_dir: Path,
    args: argparse.Namespace,
    samplesheet_report: dict[str, Any],
    preflight_result: Any,
    pipeline_source: dict[str, Any],
) -> int:
    somatic_pairs = [
        pairing for pairing in (samplesheet_report.get("pairings", []) or [])
        if pairing.get("mode") == "somatic_paired"
    ]
    payload = {
        "ok": True,
        "skill": SKILL_NAME,
        "skill_version": SKILL_VERSION,
        "mode": "check",
        "demo": bool(args.demo),
        "profile": args.profile,
        "step": getattr(args, "step", None) or "mapping",
        "tools": _tool_list(getattr(args, "tools", None)),
        "skip_tools": _tool_list(getattr(args, "skip_tools", None)),
        "samplesheet": {
            "sample_count": samplesheet_report.get("sample_count", 0),
            "patient_count": samplesheet_report.get("patient_count", 0),
            "analysis_groups": len(samplesheet_report.get("pairings", []) or []),
            "somatic_pairs": len(somatic_pairs),
            "analysis_mode": samplesheet_report.get("analysis_mode"),
        },
        "preflight": {
            "warnings": list(getattr(preflight_result, "warnings", []) or []),
            "notes": list(getattr(preflight_result, "notes", []) or []),
        },
        "pipeline_source": pipeline_source,
    }
    # check_result.json lands at the output root — parallel to result.json and
    # matching nfcore-rnaseq/scrnaseq. Preflight's _IGNORED_ROOT_NAMES tolerates
    # it so a subsequent real run in the same directory is not blocked.
    output_dir.mkdir(parents=True, exist_ok=True)
    write_text_lf(
        output_dir / "check_result.json", json.dumps(payload, indent=2, default=str)
    )
    _print("[check] Preflight passed.")
    _print(f"  Samples: {payload['samplesheet']['sample_count']}")
    _print(f"  Somatic pairs: {payload['samplesheet']['somatic_pairs']}")
    _print(f"  Tools:   {','.join(payload['tools']) or '<none>'}")
    _print(f"  Warnings: {len(payload['preflight']['warnings'])}")
    return 0


# ---------------------------------------------------------------------------
# Downstream handoff
# ---------------------------------------------------------------------------


def _emit_downstream_handoff(
    args: argparse.Namespace,
    *,
    outputs_report: Any,
    output_dir: Path,
) -> None:
    """Write reproducibility/sarek_downstream_handoff.{sh,json}."""
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)
    skill = args.downstream_skill
    if skill is not None and skill not in _DOWNSTREAM_CHOICES:
        _print(
            f"[done] WARNING: --downstream-skill {skill!r} is not in the supported set "
            f"({', '.join(_DOWNSTREAM_CHOICES)}). Handoff still emitted as a placeholder."
        )

    annotation_dir = output_dir / "upstream" / "results" / "annotation"

    routes: dict[str, dict[str, Any]] = {
        "clinical-variant-reporter": {
            "description": "ACMG/AMP clinical variant interpretation from VEP-annotated VCFs.",
            "example": (
                "python3 skills/clinical-variant-reporter/clinical_variant_reporter.py "
                f"--input {shlex.quote(str(annotation_dir))} --output <dir>"
            ),
        },
        "clinical-trial-finder": {
            "description": "Match patient variants to ClinicalTrials.gov / EUCTR trials.",
            "example": (
                "python3 skills/clinical-trial-finder/clinical_trial_finder.py "
                f"--input {shlex.quote(str(annotation_dir))} --output <dir>"
            ),
        },
        "omics-target-evidence-mapper": {
            "description": "Aggregate target-level evidence across multi-omic sources.",
            "example": (
                "python3 skills/omics-target-evidence-mapper/omics_target_evidence_mapper.py "
                f"--input {shlex.quote(str(annotation_dir))} --output <dir>"
            ),
        },
        "wes-clinical-report-en": {
            "description": "Render a WES clinical PDF report in English.",
            "example": (
                "python3 skills/wes-clinical-report-en/wes_clinical_report_en.py "
                f"--report-dir {shlex.quote(str(output_dir))} --output-dir <pdf_dir>"
            ),
        },
        "wes-clinical-report-es": {
            "description": "Render a WES clinical PDF report in Spanish.",
            "example": (
                "python3 skills/wes-clinical-report-es/wes_clinical_report_es.py "
                f"--report-dir {shlex.quote(str(output_dir))} --output-dir <pdf_dir>"
            ),
        },
    }

    json_payload = {
        "ok": True,
        "selected_skill": skill,
        "routes": routes,
        "annotation_dir_hint": str(annotation_dir),
    }
    write_text_lf(
        repro_dir / "sarek_downstream_handoff.json", json.dumps(json_payload, indent=2)
    )

    lines = [
        "#!/usr/bin/env bash",
        "# Downstream handoff template — generated by nfcore-sarek-wrapper",
        "# Pick a route below, fill in placeholders, then run.",
        "set -euo pipefail",
        "",
    ]
    for name, info in routes.items():
        prefix = "# >>> " if name == skill else "# "
        lines.append(f"{prefix}{name}: {info['description']}")
        lines.append(f"{prefix}{info['example']}")
        lines.append("")
    sh_path = repro_dir / "sarek_downstream_handoff.sh"
    write_text_lf(sh_path, "\n".join(lines) + "\n")
    try:
        sh_path.chmod(sh_path.stat().st_mode | 0o111)
    except OSError:
        pass
    _print(f"[done] downstream handoff: {sh_path}")


# ---------------------------------------------------------------------------
# Error handling
# ---------------------------------------------------------------------------


def _handle_skill_error(output_dir: Path, exc: SkillError, *, verbose: bool) -> int:
    payload = exc.to_dict()
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
    _write_error_result_if_safe(output_dir, payload)
    # Machine-readable error on stderr: always available even when result.json
    # cannot be written (e.g. the output dir is a file or inside the repo).
    print(json.dumps(payload, indent=2), file=sys.stderr)
    # Exit 1 for a SkillError (parity with nfcore-scrnaseq/rnaseq). Exit 2 is
    # reserved for argparse CLI-usage errors, so reusing it for validation
    # failures would make the two indistinguishable to machine consumers.
    return 1


def _handle_unexpected_error(output_dir: Path, exc: Exception, *, verbose: bool) -> int:
    payload = {
        "status": "error",
        "ok": False,
        "stage": "internal",
        "error_code": ErrorCode.UNEXPECTED_ERROR,
        "message": str(exc),
        "fix": "Report this as a bug. Include the traceback and command arguments.",
        "details": {"exception_type": type(exc).__name__},
    }
    _print("")
    _print("================ Internal error ================")
    _print(f"  type:    {type(exc).__name__}")
    _print(f"  message: {exc}")
    if verbose:
        _print(_indent(traceback.format_exc(), 4))
    else:
        _print("  Re-run with --verbose for traceback.")
    _print("================================================")
    _write_error_result_if_safe(output_dir, payload)
    # Machine-readable error on stderr (always available; see _handle_skill_error).
    print(json.dumps(payload, indent=2), file=sys.stderr)
    return 1


def _write_error_result_if_safe(output_dir: Path, payload: dict[str, Any]) -> None:
    code = payload.get("error_code")
    if code in {ErrorCode.OUTPUT_DIR_NOT_EMPTY, ErrorCode.OUTPUT_DIR_NOT_WRITABLE}:
        return
    text = json.dumps(payload, indent=2, default=str)
    # The failure marker lives at the output root, next to where a successful
    # run writes result.json, so consumers read <output>/result.json regardless
    # of outcome (parity with nfcore-rnaseq/scrnaseq).
    try:
        output_dir.mkdir(parents=True, exist_ok=True)
        write_text_lf(output_dir / "result.json", text)
    except OSError:
        return


def _indent(text: str, n: int) -> str:
    pad = " " * n
    return "\n".join(pad + line for line in text.splitlines())


if __name__ == "__main__":
    sys.exit(main())
