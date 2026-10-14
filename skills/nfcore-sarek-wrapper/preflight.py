# ruff: noqa: E402
# nfcore-sarek-wrapper / preflight.py
"""Preflight checks for the nf-core/sarek wrapper.

All hard failures raise ``SkillError`` with ``stage='preflight'``.  Non-fatal
findings are collected into :class:`PreflightResult` as ``warnings`` /
``notes`` and returned to the caller.

The module is laid out as one function per rule so that each rule is
individually testable and easy to reason about.  ``run_preflight()`` orchestrates
the rules, accumulating warnings and notes.
"""
from __future__ import annotations

import hashlib
import os
import re
import csv
import glob
import shutil
import subprocess
import sys
from dataclasses import dataclass, field
from pathlib import Path
from typing import Any, Iterable

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import (
    DEFAULT_ALIGNER,
    DEFAULT_GENOME,
    DEFAULT_STEP,
    DEFAULT_IGENOMES_BASE,
    JAVA_MIN_VERSION,
    NEXTFLOW_MIN_VERSION,
    is_under_tmp,
    SUPPORTED_ALIGNERS,
    SUPPORTED_ASCAT_GENOME,
    SUPPORTED_GATK_PCR_INDEL_MODEL,
    SUPPORTED_GROUP_BY_UMI,
    SUPPORTED_IGENOMES_NAMES,
    SUPPORTED_PROFILES,
    SUPPORTED_PUBLISH_DIR_MODES,
    SUPPORTED_SENTIEON_EMIT_MODE,
    SUPPORTED_SKIP_TOOLS,
    SUPPORTED_STEPS,
    SUPPORTED_TOOLS,
    SUPPORTED_UMI_LOCATIONS,
    SUPPORTED_USE_GATK_SPARK,
    SUPPORTED_VEP_OUT_FORMAT,
)


_SUBPROCESS_TIMEOUT = 30
_EMAIL_RE = re.compile(r"^[A-Za-z0-9_\-.]+@[A-Za-z0-9_\-.]+\.[A-Za-z]{2,}$")
_PROFILE_TOKEN_RE = re.compile(r"^[A-Za-z0-9_.-]+$")
_URI_PREFIXES = ("s3://", "gs://", "az://", "http://", "https://", "ftp://", "ftps://")
_SENTIEON_TOOLS = frozenset({
    "sentieon_dedup",
    "sentieon_dnascope",
    "sentieon_haplotyper",
    "sentieon_tnscope",
})

# fgbio read-structure grammar: a sequence of <length><code> segments where the
# length is a positive integer or `+` (the variable-length remainder, allowed
# only as the final segment) and the code is one of T/B/M/S. Sarek validates
# this early via fgbio's readStructure(); we mirror it with a regex so invalid
# structures (e.g. "8X", "lane 1") fail at preflight instead of late in Nextflow.
_READ_STRUCTURE_TOKEN_RE = re.compile(r"^(?:\d+[TBMS])*(?:\+[TBMS])?$")


def _is_valid_read_structure_token(token: str) -> bool:
    return bool(token) and bool(_READ_STRUCTURE_TOKEN_RE.fullmatch(token)) and any(c in "TBMS" for c in token)

# SNV/indel callers eligible for --snv_consensus_calling (§5.4 rule 13/14).
# Sarek's `--tools mpileup` path emits VCFs tagged `variantcaller: 'bcftools'`,
# which is in post_variantcalling/main.nf `small_variantcallers`. The separate
# internal `samtools` pileup stream used by ControlFREEC is the excluded one.
_SNV_CALLERS = frozenset({
    "strelka", "mutect2", "freebayes", "deepvariant", "mpileup",
    "haplotypecaller", "lofreq", "muse",
    "sentieon_haplotyper", "sentieon_dnascope", "sentieon_tnscope",
})
_GERMLINE_VCF_CALLERS = frozenset({
    "deepvariant", "freebayes", "haplotypecaller", "manta", "mpileup",
    "sentieon_dnascope", "sentieon_haplotyper", "strelka", "tiddit",
})
# Varlociraptor consumes VCF candidate streams produced by another caller in
# POST_VARIANTCALLING; analytical/QC-only tools cannot seed those streams.
_VARLOCIRAPTOR_CANDIDATE_CALLERS = frozenset({
    "deepvariant", "freebayes", "haplotypecaller", "lofreq", "manta",
    "mpileup", "muse", "mutect2", "sentieon_dnascope",
    "sentieon_haplotyper", "sentieon_tnscope", "strelka", "tiddit",
})
_ANNOTATABLE_VARIANT_OUTPUT_TOOLS = _VARLOCIRAPTOR_CANDIDATE_CALLERS | {"varlociraptor"}
_VARIANT_CALLING_TOOLS = frozenset({
    "ascat", "cnvkit", "controlfreec", "deepvariant", "freebayes",
    "haplotypecaller", "indexcov", "lofreq", "manta", "mpileup",
    "msisensor2", "msisensorpro", "muse", "mutect2",
    "sentieon_dnascope", "sentieon_haplotyper", "sentieon_tnscope",
    "strelka", "tiddit", "varlociraptor",
})
_ANNOTATION_TOOLS = frozenset({"bcfann", "merge", "snpeff", "snpsift", "vep"})
# Backward-compatible private spelling used by the wrapper's existing audit
# checks; the canonical source remains schemas.SUPPORTED_SKIP_TOOLS.
_SUPPORTED_SKIP_TOOLS = SUPPORTED_SKIP_TOOLS

# Tool/resource fallback availability in the exact 3.8.1 iGenomes catalogue.
# These sets are derived from conf/igenomes.config; a selected genome is not
# evidence that it supplies every optional tool resource.
_IGENOMES_RESOURCE_BUNDLES: dict[str, frozenset[str]] = {
    "ascat": frozenset({"GATK.GRCh37", "GATK.GRCh38"}),
    "controlfreec": frozenset({"GATK.GRCh37", "GATK.GRCh38", "GRCm38"}),
    "msisensor2": frozenset({"GATK.GRCh37", "GATK.GRCh38"}),
    "msisensorpro": frozenset({"GATK.GRCh37", "GATK.GRCh38"}),
    "pon": frozenset({"GATK.GRCh38"}),
    "germline_resource": frozenset({"GATK.GRCh37", "GATK.GRCh38", "testdata.nf-core.sarek"}),
    "ngscheckmate": frozenset({
        "GATK.GRCh37", "GATK.GRCh38", "Ensembl.GRCh37",
        "NCBI.GRCh38", "testdata.nf-core.sarek",
    }),
    "snpeff": frozenset({
        "GATK.GRCh37", "GATK.GRCh38", "Ensembl.GRCh37", "NCBI.GRCh38",
        "GRCm38", "UMD3.1", "WBcel235", "CanFam3.1", "R64-1-1",
        "hg38", "hg19", "mm10", "testdata.nf-core.sarek",
    }),
    "vep": frozenset({
        "GATK.GRCh37", "GATK.GRCh38", "Ensembl.GRCh37", "NCBI.GRCh38",
        "GRCm38", "UMD3.1", "WBcel235", "CanFam3.1", "R64-1-1",
        "hg38", "hg19", "mm10", "testdata.nf-core.sarek",
    }),
    "bqsr": frozenset({"GATK.GRCh37", "GATK.GRCh38", "GRCm38", "testdata.nf-core.sarek"}),
}

# Backend tokens that map to an executable / daemon we can probe.
_BACKEND_PROFILES = frozenset({
    "docker", "podman", "singularity", "apptainer",
    "shifter", "charliecloud", "conda", "mamba",
})

# Reference path params that must exist on disk if set (and not a URI). §5.5.
REFERENCE_PATH_PARAMS: tuple[str, ...] = (
    "fasta", "fasta_fai", "dict",
    "bwa", "bwamem2", "dragmap",
    "dbsnp", "dbsnp_tbi",
    "known_indels", "known_indels_tbi",
    "known_snps", "known_snps_tbi",
    "germline_resource", "germline_resource_tbi",
    "pon", "pon_tbi",
    "intervals",
    "ascat_alleles", "ascat_loci", "ascat_loci_gc", "ascat_loci_rt",
    "chr_dir", "mappability",
    "msisensor2_models", "msisensorpro_scan",
    "ngscheckmate_bed",
    "sentieon_dnascope_model",
    "snpeff_cache", "vep_cache",
    "dbnsfp", "dbnsfp_tbi",
    "mastermind_file",
    "phenotypes_file", "phenotypes_file_tbi",
    "spliceai_snv", "spliceai_snv_tbi", "spliceai_indel", "spliceai_indel_tbi",
    "bcftools_annotations", "bcftools_annotations_tbi", "bcftools_header_lines",
    "condel_config", "cnvkit_reference", "cf_chrom_len",
    "bbsplit_fasta_list", "bbsplit_index",
    "bcftools_columns", "snpsift_databases",
    "varlociraptor_scenario_tumor_only",
    "varlociraptor_scenario_somatic",
    "varlociraptor_scenario_germline",
    "multiqc_config", "multiqc_logo", "multiqc_methods_description",
    "igenomes_base",
)
_REFERENCE_PATH_PATTERN_PARAMS = frozenset({"known_indels", "known_indels_tbi"})


# --- Public API ------------------------------------------------------------

@dataclass(frozen=True)
class PreflightResult:
    """Outcome of preflight checks.

    Hard failures raise ``SkillError`` and never reach this dataclass; the
    fields contain only non-fatal findings.
    """

    warnings: list[str] = field(default_factory=list)
    notes: list[str] = field(default_factory=list)
    java_version: str = ""
    nextflow_version: str = ""


def _check_remote_inputs(
    params: dict[str, Any],
    samplesheet: dict[str, Any],
    *,
    allow_remote_inputs: bool,
) -> None:
    """Local-first by default: reject remote samplesheet inputs and user-supplied
    reference paths (s3://, gs://, https://, ftp://, …) unless allow_remote_inputs
    is set, in which case warn that genetic data is fetched over the network. The
    object-store ``--work-dir`` and the public iGenomes mirror base (``igenomes_base``,
    default ``s3://ngi-igenomes/``) are intentionally remote and are not gated."""
    remote: list[str] = []
    for key in ("fastq_paths", "bam_paths", "cram_paths"):
        for path in samplesheet.get(key, []) or []:
            if "://" in str(path):
                remote.append(str(path))
    for ref_field in REFERENCE_PATH_PARAMS:
        if ref_field == "igenomes_base":
            continue  # public reference-mirror base, remote by design
        value = params.get(ref_field)
        if value and "://" in str(value):
            remote.append(str(value))
    if not remote:
        return
    unique = sorted(set(remote))
    if allow_remote_inputs:
        print(
            "WARNING: --allow-remote-inputs is set; "
            f"{len(unique)} input/reference path(s) will be fetched over the network, "
            "so genetic data may leave the local machine: " + ", ".join(unique),
            file=sys.stderr,
        )
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.REMOTE_INPUT_NOT_ALLOWED,
        message=(
            "Remote input/reference URIs are not allowed by default (local-first); "
            f"detected {len(unique)} remote path(s)."
        ),
        fix=(
            "Download the data locally and use local paths, or pass "
            "--allow-remote-inputs to fetch them over the network at your own risk "
            "(genetic data will leave the local machine)."
        ),
        details={"remote_paths": unique},
    )


def run_preflight(
    *,
    params: dict[str, Any],
    samplesheet: dict[str, Any],
    pipeline_source: dict[str, Any],
    output_dir: Path,
    repo_root: Path,
    resume: bool = False,
    resume_manifest: dict | None = None,
    allow_remote_inputs: bool = False,
) -> PreflightResult:
    """Run every preflight rule for the Sarek wrapper.

    The function dispatches each rule to a private helper.  Rules raise
    :class:`SkillError` on hard failures; warnings are accumulated into the
    returned :class:`PreflightResult`.
    """
    warnings_acc: list[str] = []
    notes_acc: list[str] = []

    # --demo requires network (nf-core's remote -profile test). Fail fast under NXF_OFFLINE.
    _check_demo_network(params)

    # Local-first gate: reject remote inputs/references unless opted in.
    _check_remote_inputs(params, samplesheet, allow_remote_inputs=allow_remote_inputs)

    # §5.1 — cross-cutting checks ------------------------------------------
    java_version = _check_java()
    nextflow_version = _check_nextflow()
    profile = str(params.get("profile") or "")
    _check_profile_string(profile)

    # §5.2 — step + aligner + enumerated params + tool tokens -------------
    # Validate cheap, environment-independent user input BEFORE probing the
    # container backend, so a typo'd --step/--aligner/--tools/enum surfaces its
    # own actionable error instead of being masked by a missing/stopped backend.
    _check_step(params)
    _check_aligner(params)
    _check_enums(params)
    warnings_acc.extend(_check_genome_known(params))
    _check_tools_required_for_step(params)
    _check_tools_known(_normalize_tools(params.get("tools")))

    # --- environment / IO probes (after input is known good) -------------
    _check_backends_for_profile(profile)
    _check_output_dir(output_dir, repo_root=repo_root, resume=resume)
    warnings_acc.extend(_check_macos_docker_tmp(profile, output_dir))
    _check_email(params.get("email"))
    _check_email(params.get("email_on_fail"), field_name="email_on_fail")
    _check_pipeline_source(pipeline_source)

    # §5.3 — tools × pairing matrix ----------------------------------------
    # Sarek processes each patient independently and never enforces a single
    # global analysis mode, so we evaluate tools against the SET of per-patient
    # modes present (germline / tumor_only / somatic_paired). A tool is accepted
    # when at least one patient matches its required mode; a mixed samplesheet
    # (e.g. some germline-only patients + some tumor/normal pairs) is valid.
    analysis_mode = str(samplesheet.get("analysis_mode") or "germline")
    pairings = samplesheet.get("pairings") or []
    present_modes = {str(p.get("mode")) for p in pairings if p.get("mode")}
    if not present_modes:
        present_modes = {analysis_mode}
    rows_by_patient = _extract_rows_by_patient(samplesheet)
    tools = _normalize_tools(params.get("tools"))
    tool_warnings = _check_tools_against_pairing(
        params=params,
        tools=tools,
        present_modes=present_modes,
        rows_by_patient=rows_by_patient,
    )
    warnings_acc.extend(tool_warnings)

    # §5.4 — flag compatibility -------------------------------------------
    flag_warnings = _check_flag_compatibility(params=params, tools=tools)
    warnings_acc.extend(flag_warnings)

    # §5.5 — reference path existence + documented-suffix sanity ----------
    _check_reference_paths(params)
    warnings_acc.extend(_check_reference_path_suffixes(params))

    # §5.6 — annotation cache ---------------------------------------------
    ann_warnings = _check_annotation_cache(params=params, output_dir=output_dir)
    warnings_acc.extend(ann_warnings)

    # §5.1 — resume manifest drift (after everything else so we can include
    # current values) ------------------------------------------------------
    if resume and resume_manifest is not None:
        _check_resume_drift(
            params=params,
            samplesheet=samplesheet,
            pipeline_source=pipeline_source,
            manifest=resume_manifest,
        )

    return PreflightResult(
        warnings=warnings_acc,
        notes=notes_acc,
        java_version=java_version,
        nextflow_version=nextflow_version,
    )


# --- §5.1: Java -----------------------------------------------------------

def _command_output(args: list[str]) -> str:
    """Return combined stdout/stderr of ``args`` (empty string on error)."""
    try:
        proc = subprocess.run(
            args,
            capture_output=True,
            text=True,
            errors="replace",
            timeout=_SUBPROCESS_TIMEOUT,
        )
    except (subprocess.TimeoutExpired, OSError):
        return ""
    return ((proc.stdout or "") + "\n" + (proc.stderr or "")).strip()


def _parse_version_tuple(text: str) -> tuple[int, ...]:
    m = re.search(r"\b(\d+)\.(\d+)\.(\d+)\b", text)
    if m:
        return (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.search(r"\b(\d+)\.(\d+)\b", text)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    m = re.search(r"\b(\d+)\b", text)
    if m:
        return (int(m.group(1)),)
    return ()


def _pad_version(t: tuple[int, ...], length: int = 3) -> tuple[int, ...]:
    return t + (0,) * max(0, length - len(t))


def _detected_version_string(text: str) -> str:
    """Return the version exactly as reported by the tool (e.g. '26.04.3').

    The parsed tuple is for *comparison* only; reconstructing a string from it
    drops zero-padding (26.04.3 → 26.4.3), which is not a real Nextflow release
    and would break NXF_VER / conda pins on replay. Always store the raw token.
    """
    for pattern in (r"\b(\d+\.\d+\.\d+)\b", r"\b(\d+\.\d+)\b", r"\b(\d+)\b"):
        m = re.search(pattern, text)
        if m:
            return m.group(1)
    return ""


def _check_java() -> str:
    """Java >=17 must be present on PATH. Returns the detected version string."""
    java_path = shutil.which("java")
    if not java_path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_JAVA,
            message="Required executable `java` was not found.",
            fix="Install Java 17 or newer and ensure it is on PATH.",
            details={"executable": "java"},
        )
    version_text = _command_output(["java", "-version"])
    version = _parse_version_tuple(version_text)
    if not version:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.JAVA_VERSION_TOO_OLD,
            message="Java is installed but its version could not be parsed.",
            fix="Install Java 17 or newer and ensure `java -version` works.",
            details={"java_path": java_path, "output": version_text},
        )
    version_str = _detected_version_string(version_text)
    if version[0] < JAVA_MIN_VERSION:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.JAVA_VERSION_TOO_OLD,
            message=f"Java {version_str} is older than the required {JAVA_MIN_VERSION}.",
            fix=f"Install Java {JAVA_MIN_VERSION} or newer.",
            details={"detected_version": version_str, "min": JAVA_MIN_VERSION},
        )
    return version_str


# --- §5.1: Nextflow -------------------------------------------------------

def _check_nextflow() -> str:
    """Nextflow >= minimum must be present. Returns the detected version string."""
    minimum = ".".join(map(str, NEXTFLOW_MIN_VERSION))
    nextflow_path = shutil.which("nextflow")
    if not nextflow_path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_NEXTFLOW,
            message="Required executable `nextflow` was not found.",
            fix=f"Install Nextflow >={minimum} and ensure it is on PATH.",
            details={"executable": "nextflow"},
        )
    version_text = _command_output(["nextflow", "-version"])
    version = _parse_version_tuple(version_text)
    if not version:
        # nextflow IS present but its version cannot be confirmed — this is a
        # version-gate failure, not an absent binary (NEXTFLOW_NOT_FOUND).
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.NEXTFLOW_VERSION_TOO_OLD,
            message="Nextflow is installed but its version could not be parsed.",
            fix=f"Install Nextflow >={minimum} and ensure `nextflow -version` works.",
            details={"nextflow_path": nextflow_path, "output": version_text},
        )
    version_str = _detected_version_string(version_text)
    if _pad_version(version) < _pad_version(NEXTFLOW_MIN_VERSION):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.NEXTFLOW_VERSION_TOO_OLD,
            message=(
                "Nextflow "
                f"{version_str} is older than the required "
                f"{minimum}."
            ),
            fix=f"Upgrade Nextflow to {minimum} or newer.",
            details={
                "detected_version": version_str,
                "min_version": minimum,
            },
        )
    return version_str


# --- §5.1: Profile & backends --------------------------------------------

def _profile_tokens(profile: str) -> list[str]:
    return [p.strip() for p in profile.split(",") if p.strip()]


def _check_profile_string(profile: str) -> None:
    if not profile:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PROFILE,
            message="`profile` is required and must be non-empty.",
            fix=f"Choose one or more comma-separated profiles from: {', '.join(sorted(SUPPORTED_PROFILES))}.",
            details={"profile": profile},
        )
    tokens = _profile_tokens(profile)
    invalid = [t for t in tokens if not _PROFILE_TOKEN_RE.fullmatch(t)]
    if invalid:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PROFILE,
            message=f"Invalid profile token(s): {', '.join(invalid)}.",
            fix="Use comma-separated Nextflow profile tokens containing only letters, numbers, dots, underscores or hyphens.",
            details={
                "profile": profile,
                "invalid_tokens": invalid,
                "supported": sorted(SUPPORTED_PROFILES),
            },
        )


def _check_backends_for_profile(profile: str) -> None:
    """For each backend token in the profile, probe the local binary/daemon."""
    tokens = [t for t in _profile_tokens(profile) if t in _BACKEND_PROFILES]
    for token in tokens:
        if token == "docker":
            _probe_daemon(token, binary="docker")
        elif token == "podman":
            _probe_daemon(token, binary="podman")
        elif token in {"singularity", "apptainer"}:
            _probe_binary_present(token, binary=token)
        elif token in {"conda", "mamba"}:
            _probe_binary_present(token, binary=token)
        elif token == "shifter":
            _probe_binary_present(token, binary="shifter")
        elif token == "charliecloud":
            _probe_binary_present(token, binary="ch-run")


# Runtime-specific backend error codes, mirroring nfcore-scrnaseq-wrapper so the
# two sibling skills classify the same backend failure identically (cross-skill
# consistency; activates the previously-dead specific codes).
_MISSING_BACKEND_CODE = {
    "docker": ErrorCode.MISSING_DOCKER,
    "podman": ErrorCode.MISSING_PODMAN,
    "singularity": ErrorCode.MISSING_SINGULARITY,
    "apptainer": ErrorCode.MISSING_SINGULARITY,
    "conda": ErrorCode.MISSING_CONDA,
    "mamba": ErrorCode.MISSING_CONDA,
    "shifter": ErrorCode.MISSING_HPC_RUNTIME,
    "charliecloud": ErrorCode.MISSING_HPC_RUNTIME,
}
_DAEMON_DOWN_CODE = {
    "docker": ErrorCode.DOCKER_NOT_RUNNING,
    "podman": ErrorCode.PODMAN_NOT_RUNNING,
}


def _probe_daemon(profile: str, *, binary: str) -> None:
    path = shutil.which(binary)
    if not path:
        raise SkillError(
            stage="preflight",
            error_code=_MISSING_BACKEND_CODE.get(profile, ErrorCode.BACKEND_UNAVAILABLE),
            message=f"Profile '{profile}' selected but `{binary}` is not installed.",
            fix=f"Install {binary} or remove the '{profile}' profile token.",
            details={"profile": profile, "binary": binary},
        )
    try:
        info = subprocess.run(
            [binary, "info"],
            capture_output=True,
            text=True,
            timeout=_SUBPROCESS_TIMEOUT,
        )
        ok = info.returncode == 0
    except (subprocess.TimeoutExpired, OSError):
        ok = False
    if not ok:
        raise SkillError(
            stage="preflight",
            error_code=_DAEMON_DOWN_CODE.get(profile, ErrorCode.BACKEND_UNAVAILABLE),
            message=f"{binary} is installed but its service is not available.",
            fix=f"Start {binary} (run `{binary} info` to diagnose).",
            details={"profile": profile, "binary": binary},
        )


def _probe_binary_present(profile: str, *, binary: str) -> None:
    if not shutil.which(binary):
        raise SkillError(
            stage="preflight",
            error_code=_MISSING_BACKEND_CODE.get(profile, ErrorCode.BACKEND_UNAVAILABLE),
            message=f"Profile '{profile}' selected but `{binary}` was not found on PATH.",
            fix=f"Install {binary} or remove the '{profile}' profile token.",
            details={"profile": profile, "binary": binary},
        )


# --- §5.1: Output dir -----------------------------------------------------

def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _check_macos_docker_tmp(profile: str, output_dir: Path) -> list[str]:
    """Warn when a macOS + Docker run writes under /tmp.

    Colima (a common macOS Docker runtime) only shares the user HOME into its VM;
    /tmp and /private/tmp live on the VM's own filesystem and are NOT shared with
    the host, so Nextflow work-dir files are invisible to containers and the run
    fails with a confusing 'No such file or directory'. Mirrors the
    nfcore-scrnaseq-wrapper guard so both skills behave the same on macOS.
    """
    if sys.platform != "darwin" or "docker" not in _profile_tokens(profile):
        return []
    if is_under_tmp(output_dir):
        return [
            "Output directory is under /tmp. On macOS with Colima, Docker "
            "containers cannot see files written to /tmp (the VM uses its own "
            "separate /tmp). Move --output under your home directory to avoid "
            "'No such file or directory' errors."
        ]
    return []


# Root-level entries that must not trip OUTPUT_DIR_NOT_EMPTY: the reproducibility
# bundle (written by --check and staging), OS/VCS scratch files, and a prior
# --check summary (check_result.json). Mirrors nfcore-rnaseq/scrnaseq so the
# three wrappers accept the same set of pre-existing output-dir contents; on
# macOS, Finder writes .DS_Store into any visited folder.
_IGNORED_ROOT_NAMES = frozenset(
    {"reproducibility", ".DS_Store", ".gitkeep", ".gitignore", "Thumbs.db", "check_result.json"}
)


def _check_output_dir(output_dir: Path, *, repo_root: Path, resume: bool) -> None:
    output_dir = output_dir.expanduser().resolve()
    repo_root = repo_root.expanduser().resolve()
    if _is_relative_to(output_dir, repo_root):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.OUTPUT_DIR_INSIDE_REPO,
            message="Output directory cannot be inside the ClawBio source tree.",
            fix="Choose an output directory outside the repository (e.g. under your analysis workspace).",
            details={"output_dir": str(output_dir), "repo_root": str(repo_root)},
        )
    if not output_dir.exists():
        return
    if not output_dir.is_dir() or not os.access(output_dir, os.W_OK):
        raise SkillError(
            stage="preflight", error_code=ErrorCode.OUTPUT_DIR_NOT_WRITABLE,
            message="Output path is not a writable directory.",
            fix="Choose a writable directory path for --output.",
            details={"output_dir": str(output_dir)},
        )
    if resume:
        return
    try:
        entries = list(output_dir.iterdir())
    except FileNotFoundError:
        return
    if entries:
        entries = [entry for entry in entries if entry.name not in _IGNORED_ROOT_NAMES]
    if entries:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.OUTPUT_DIR_NOT_EMPTY,
            message="Output directory already contains files.",
            fix="Choose an empty directory, or re-run with --resume against a compatible manifest.",
            details={
                "output_dir": str(output_dir),
                "entry_count": len(entries),
                "examples": sorted(e.name for e in entries)[:5],
            },
        )


# --- §5.1: Email ----------------------------------------------------------

def _check_email(email: Any, *, field_name: str = "email") -> None:
    if email is None or email == "":
        return
    if not isinstance(email, str) or not _EMAIL_RE.match(email):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_EMAIL,
            message=f"Invalid {field_name} address: {email!r}.",
            fix="Provide a syntactically valid email address (e.g. user@example.com).",
            details={"field": field_name, "value": email},
        )


# --- §5.1: Pipeline source sanity ----------------------------------------

def _check_pipeline_source(pipeline_source: dict[str, Any]) -> None:
    kind = pipeline_source.get("source_kind")
    if kind == "remote_repo" and not pipeline_source.get("resolved_version"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.PIPELINE_SOURCE_INVALID,
            message="Remote pipeline source has no resolved version.",
            fix="Re-run after resolving the pipeline version (e.g. set --pipeline-version).",
            details={"pipeline_source": pipeline_source},
        )


def _nextflow_secret_exists(name: str) -> bool:
    """Return True when a Nextflow secret exists without exposing its value."""
    try:
        proc = subprocess.run(
            ["nextflow", "secrets", "get", name],
            capture_output=True,
            text=True,
            timeout=_SUBPROCESS_TIMEOUT,
        )
    except (subprocess.TimeoutExpired, OSError):
        return False
    return proc.returncode == 0 and bool((proc.stdout or "").strip())


def _sentieon_license_available() -> bool:
    return bool(
        os.environ.get("SENTIEON_LICENSE_BASE64")
        or os.environ.get("SENTIEON_LICENSE")
        or _nextflow_secret_exists("SENTIEON_LICENSE_BASE64")
    )


def _require_sentieon_license(context: str) -> None:
    if _sentieon_license_available():
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.MISSING_SENTIEON_LICENSE,
        message=f"{context} requires a Sentieon license, but none was detected.",
        fix=(
            "Set the Nextflow secret SENTIEON_LICENSE_BASE64 as documented by nf-core/sarek, "
            "or export SENTIEON_LICENSE_BASE64 / SENTIEON_LICENSE before running."
        ),
        details={
            "context": context,
            "checked": ["SENTIEON_LICENSE_BASE64", "SENTIEON_LICENSE", "nextflow secret SENTIEON_LICENSE_BASE64"],
        },
    )


# --- §5.2: step -----------------------------------------------------------

def _check_step(params: dict[str, Any]) -> None:
    step = str(params.get("step") or "").strip()
    if not step:
        # No step provided — samplesheet builder uses 'mapping' default.
        return
    if step not in SUPPORTED_STEPS:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_STEP,
            message=f"Unsupported --step value '{step}'.",
            fix=f"Choose one of: {', '.join(sorted(SUPPORTED_STEPS))}.",
            details={"step": step, "supported": sorted(SUPPORTED_STEPS)},
        )


# --- §5.2b: aligner enumerated -------------------------------------------

def _check_aligner(params: dict[str, Any]) -> None:
    aligner = str(params.get("aligner") or "").strip()
    if aligner and aligner not in SUPPORTED_ALIGNERS:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_ALIGNER,
            message=f"Unsupported --aligner value '{aligner}'.",
            fix=f"Choose one of: {', '.join(sorted(SUPPORTED_ALIGNERS))}.",
            details={"aligner": aligner, "supported": sorted(SUPPORTED_ALIGNERS)},
        )


def _check_enums(params: dict[str, Any]) -> None:
    """Validate enumerated pass-through params against their schema enums so an
    invalid value fails at preflight rather than late inside Nextflow."""
    for name, allowed in (
        ("umi_location", SUPPORTED_UMI_LOCATIONS),
        ("group_by_umi_strategy", SUPPORTED_GROUP_BY_UMI),
        ("vep_out_format", SUPPORTED_VEP_OUT_FORMAT),
        ("ascat_genome", SUPPORTED_ASCAT_GENOME),
        ("publish_dir_mode", SUPPORTED_PUBLISH_DIR_MODES),
    ):
        value = params.get(name)
        if value in (None, ""):
            continue
        if str(value) not in allowed:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"Unsupported --{name.replace('_', '-')} value '{value}'.",
                fix=f"Choose one of: {', '.join(sorted(allowed))}.",
                details={name: value, "supported": sorted(allowed)},
            )
    # use_gatk_spark is a comma/space-separated list of GATK process names.
    spark = params.get("use_gatk_spark")
    if spark not in (None, "", True, False):
        tokens = [t.strip() for t in re.split(r"[,\s]+", str(spark)) if t.strip()]
        bad = [t for t in tokens if t not in SUPPORTED_USE_GATK_SPARK]
        if bad:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"Unsupported --use-gatk-spark value(s): {bad}.",
                fix=f"Choose from: {', '.join(sorted(SUPPORTED_USE_GATK_SPARK))}.",
                details={"use_gatk_spark": spark, "supported": sorted(SUPPORTED_USE_GATK_SPARK)},
            )

    skip_tools = _normalize_tools(params.get("skip_tools"))
    bad_skip_tools = [t for t in skip_tools if t not in SUPPORTED_SKIP_TOOLS]
    if bad_skip_tools:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_TOOLS,
            message=f"Unknown --skip_tools value(s): {', '.join(bad_skip_tools)}.",
            fix=f"Choose skip_tools from: {', '.join(sorted(SUPPORTED_SKIP_TOOLS))}.",
            details={"unknown": bad_skip_tools, "supported": sorted(SUPPORTED_SKIP_TOOLS)},
        )

    for name in ("sentieon_haplotyper_emit_mode", "sentieon_dnascope_emit_mode"):
        value = params.get(name)
        if value not in (None, "") and str(value) not in SUPPORTED_SENTIEON_EMIT_MODE:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"Unsupported --{name.replace('_', '-')} value '{value}'.",
                fix=f"Use one of: {', '.join(sorted(SUPPORTED_SENTIEON_EMIT_MODE))}.",
                details={name: value, "supported": sorted(SUPPORTED_SENTIEON_EMIT_MODE)},
            )

    pcr = params.get("sentieon_dnascope_pcr_indel_model")
    if pcr not in (None, "") and str(pcr) not in SUPPORTED_GATK_PCR_INDEL_MODEL:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=f"Unsupported --sentieon-dnascope-pcr-indel-model value '{pcr}'.",
            fix=f"Use one of: {', '.join(sorted(SUPPORTED_GATK_PCR_INDEL_MODEL))}.",
            details={"sentieon_dnascope_pcr_indel_model": pcr},
        )

    max_email = params.get("max_multiqc_email_size")
    if max_email not in (None, "") and not re.fullmatch(r"\d+(\.\d+)?\.?\s*(K|M|G|T)?B", str(max_email)):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=f"Unsupported --max-multiqc-email-size value '{max_email}'.",
            fix="Use a file-size string matching the Sarek schema, e.g. 25.MB or 1 GB.",
            details={"max_multiqc_email_size": max_email},
        )

    phenotypes_types = params.get("phenotypes_include_types")
    _PT = r"Gene|Variation|QTL|StructuralVariation|SupportingStructuralVariation|RegulatoryFeature"
    phenotypes_re = re.compile(rf"^(?:{_PT})(?:&(?:{_PT}))*$")
    if phenotypes_types not in (None, "") and not phenotypes_re.fullmatch(str(phenotypes_types)):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=f"Unsupported --phenotypes-include-types value '{phenotypes_types}'.",
            fix="Use '&'-separated feature types from the official schema, e.g. Gene&Variation.",
            details={"phenotypes_include_types": phenotypes_types},
        )


def _is_default_igenomes_base(params: dict[str, Any]) -> bool:
    """True when no custom iGenomes mirror is configured.

    The wrapper writes ``igenomes_base`` only when the user passes
    ``--igenomes-base``; an absent key means the upstream default
    (``s3://ngi-igenomes/igenomes/``) is in effect. With the default mirror the
    set of valid ``--genome`` keys is exactly conf/igenomes.config; a custom
    mirror may define others, so genome validation softens to a warning there.
    """
    base = str(params.get("igenomes_base") or "").strip()
    if not base:
        return True
    return base.rstrip("/") == DEFAULT_IGENOMES_BASE.rstrip("/")


def _check_genome_known(params: dict[str, Any]) -> list[str]:
    """Validate an explicitly-provided ``--genome`` against the iGenomes catalogue.

    nf-core/sarek 3.8.1 resolves ``--genome`` against conf/igenomes.config
    (SUPPORTED_IGENOMES_NAMES). The schema declares ``genome`` as a free string
    (default ``GATK.GRCh38``), so an invalid key — e.g. bare ``GRCh38`` instead of
    ``GATK.GRCh38`` — is not caught by nf-schema and only fails late inside
    Nextflow with an opaque reference-resolution error. This gate surfaces the
    mistake at preflight.

    Only an *explicitly set* key is validated: an absent key keeps the upstream
    default (always valid), and the ``null``/``false`` disable sentinel (or
    ``--igenomes_ignore``) marks a custom-FASTA run where the label is
    informational. With the default mirror an unknown key is a hard error; with a
    custom ``--igenomes-base`` it softens to a warning (the mirror may define it).
    """
    raw = params.get("genome")
    if raw in (None, "", False):
        return []
    genome = str(raw).strip()
    if not genome or genome.lower() in ("null", "none", "false"):
        return []
    if _truthy(params.get("igenomes_ignore")):
        return []
    if genome == "testdata.nf-core.sarek":
        # nf-core/sarek's tiny test genome is defined only in the test-datasets mirror
        # (conf/test.config sets igenomes_base to that repo), never in the default
        # s3://ngi-igenomes catalogue. Used in a normal run without that override it
        # fails late during nf-schema reference validation with no warning. -profile
        # test / --demo supply the matching base; otherwise require it explicitly. A
        # custom mirror is assumed to provide the key (matches the softening below).
        if _uses_upstream_test_profile(params) or not _is_default_igenomes_base(params):
            return []
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_GENOME,
            message=(
                "--genome 'testdata.nf-core.sarek' is nf-core/sarek's test-only genome; "
                "it exists only in the nf-core test-datasets mirror, not the default "
                "iGenomes catalogue, so a normal run without a matching --igenomes-base "
                "fails late during Nextflow schema validation."
            ),
            fix=(
                "Use --demo (which runs -profile test with the correct mirror), or pass "
                "--igenomes-base https://raw.githubusercontent.com/nf-core/test-datasets/modules/data/ "
                "to match conf/test.config. For a real analysis use a production genome "
                "(e.g. GATK.GRCh38) or --fasta with --igenomes-ignore."
            ),
            details={"genome": genome},
        )
    if genome in SUPPORTED_IGENOMES_NAMES:
        return []
    # An upstream test/test_* profile supplies its own reference (and may set a
    # non-iGenomes genome label), exactly as the GENOME_REQUIRED check exempts it.
    if _uses_upstream_test_profile(params):
        return []
    if not _is_default_igenomes_base(params):
        return [
            f"--genome '{genome}' is not a standard nf-core/sarek iGenomes key. "
            "A custom --igenomes-base is set, so the wrapper assumes your mirror "
            "defines it; verify the key resolves before a long run."
        ]
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_GENOME,
        message=f"--genome '{genome}' is not a known nf-core/sarek iGenomes key.",
        fix=(
            "Use a valid iGenomes key (e.g. GATK.GRCh38, GATK.GRCh37, "
            "Ensembl.GRCh37, NCBI.GRCh38, hg38, mm10), or set --igenomes-ignore "
            "with --fasta for a custom reference."
        ),
        details={"genome": genome, "supported_count": len(SUPPORTED_IGENOMES_NAMES)},
    )


def _check_tools_required_for_step(params: dict[str, Any]) -> None:
    step = str(params.get("step") or "").strip() or "mapping"
    tools = set(_normalize_tools(params.get("tools")))
    if step in {"variant_calling", "annotate"} and not tools:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_TOOLS,
            message=f"--tools is required when --step is {step}.",
            fix="Provide at least one Sarek tool for variant calling or annotation.",
            details={"step": step, "tools": params.get("tools")},
        )
    if step == "annotate":
        ignored = sorted(tools - _ANNOTATION_TOOLS)
        if ignored:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_TOOLS,
                message="--step annotate accepts only annotation tools; other selected tools are not executed.",
                fix="Remove non-annotator tools and choose from: " + ", ".join(sorted(_ANNOTATION_TOOLS)) + ".",
                details={"step": step, "ignored_tools": ignored},
            )
    elif tools & _ANNOTATION_TOOLS and not tools & _ANNOTATABLE_VARIANT_OUTPUT_TOOLS:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_TOOLS,
            message="Annotation tools have no input VCF because no VCF-producing variant caller was selected.",
            fix=(
                "Add a compatible VCF-producing caller (for example haplotypecaller, mutect2, "
                "or strelka), or start from --step annotate with an input VCF."
            ),
            details={"step": step, "annotation_tools": sorted(tools & _ANNOTATION_TOOLS)},
        )


# --- §5.3: tools × pairing -----------------------------------------------

def _normalize_tools(raw: Any) -> list[str]:
    if not raw:
        return []
    if isinstance(raw, str):
        tokens = [t.strip().lower() for t in raw.split(",")]
    elif isinstance(raw, Iterable):
        tokens = [str(t).strip().lower() for t in raw]
    else:
        tokens = [str(raw).strip().lower()]
    return [t for t in tokens if t]


def _effective_genome(params: dict[str, Any]) -> str:
    """Return the iGenomes genome that Sarek will use.

    Sarek 3.8.1 defines ``genome = 'GATK.GRCh38'`` in ``nextflow.config``.
    The wrapper deliberately omits upstream-default parameters from
    ``params.yaml``, so an absent key means that default remains in effect.
    Conversely, an explicit null/empty sentinel disables iGenomes for custom
    reference runs.
    """
    if "genome" not in params:
        return DEFAULT_GENOME
    genome = str(params.get("genome") or "").strip()
    if not genome or genome.lower() in ("null", "none"):
        return ""
    return genome


def _has_genome(params: dict[str, Any]) -> bool:
    """True when an explicit or upstream-default iGenomes genome is selected."""
    return bool(_effective_genome(params))


def _igenomes_active(params: dict[str, Any]) -> bool:
    """True when an iGenomes ``--genome`` is in effect.

    When an iGenomes genome is selected (and not ignored), the genome config can
    supply tool-specific resources (mappability, msisensorpro_scan, …) the user
    did not pass explicitly.  In that case a missing explicit path is a warning,
    not a hard error — only when ``--igenomes_ignore`` is set (no fallback) does
    the absence become fatal.
    """
    return _has_genome(params) and not _truthy(params.get("igenomes_ignore"))


def _igenomes_supplies(params: dict[str, Any], resource_group: str) -> bool:
    return (
        _igenomes_active(params)
        and _effective_genome(params) in _IGENOMES_RESOURCE_BUNDLES[resource_group]
    )


def _explicitly_disabled(params: dict[str, Any], name: str) -> bool:
    value = params.get(name)
    return value is False or (isinstance(value, str) and value.strip().lower() == "false")


def _resource_available(params: dict[str, Any], name: str, resource_group: str) -> bool:
    """Whether a resource is supplied explicitly or inherited from iGenomes.

    Sarek documents ``--<reference> false`` as the way to disable a resource
    supplied by iGenomes.  An explicit disable therefore must not fall back to
    the selected bundle.
    """
    if _explicitly_disabled(params, name):
        return False
    if params.get(name) not in (None, ""):
        return True
    return _igenomes_supplies(params, resource_group)

def _explicit_resource_supplied(params: dict[str, Any], name: str) -> bool:
    """True only for an explicitly usable value, excluding Sarek's `false` sentinel."""
    return params.get(name) not in (None, "") and not _explicitly_disabled(params, name)


def _extract_rows_by_patient(samplesheet: dict[str, Any]) -> dict[str, list[dict[str, Any]]]:
    """Reconstruct rows-by-patient from samplesheet summary if possible.

    The samplesheet_builder returns inferred pairings but does not surface
    the per-row contamination directly.  We rely on ``pairings``
    for tumor/normal sample names and the cached samplesheet rows when
    provided.
    """
    rows = samplesheet.get("rows_by_patient")
    if isinstance(rows, dict):
        return rows
    return {}


def _check_tools_known(tools: list[str]) -> None:
    """Reject --tools tokens absent from the 3.8.1 allowlist.

    Cheap and environment-independent, so run_preflight calls this before the
    backend probe; `_check_tools_against_pairing` also calls it so the pairing
    check stays self-guarding when invoked directly.
    """
    unknown = [t for t in tools if t not in SUPPORTED_TOOLS]
    if unknown:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_TOOLS,
            message=f"Unknown --tools value(s): {', '.join(unknown)}.",
            fix=f"Choose tools from: {', '.join(sorted(SUPPORTED_TOOLS))}.",
            details={"unknown": unknown, "supported": sorted(SUPPORTED_TOOLS)},
        )


def _check_tools_against_pairing(
    *,
    params: dict[str, Any],
    tools: list[str],
    present_modes: set[str],
    rows_by_patient: dict[str, list[dict[str, Any]]],
) -> list[str]:
    warnings_acc: list[str] = []
    step = str(params.get("step") or "mapping").strip()
    annotation_only = step == "annotate"
    build_only_index = _truthy(params.get("build_only_index"))
    aligner = str(params.get("aligner") or "bwa-mem").strip().lower()
    paired_normal_germline_suppressed = (
        not annotation_only
        and not build_only_index
        and _truthy(params.get("only_paired_variant_calling"))
        and "somatic_paired" in present_modes
    )
    _check_tools_known(tools)

    if not annotation_only and not build_only_index and _truthy(params.get("concatenate_vcfs")):
        germline_input_available = (
            "germline" in present_modes
            or ("somatic_paired" in present_modes and not _truthy(params.get("only_paired_variant_calling")))
        )
        germline_vcf_tools = sorted(set(tools) & _GERMLINE_VCF_CALLERS)
        if not germline_input_available or not germline_vcf_tools:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message="--concatenate_vcfs requires a germline VCF stream, but none will be produced.",
                fix=(
                    "Provide a normal/germline sample and a germline VCF-producing caller "
                    f"({', '.join(sorted(_GERMLINE_VCF_CALLERS))}), or remove --concatenate_vcfs."
                ),
                details={
                    "present_modes": sorted(present_modes),
                    "germline_vcf_tools": germline_vcf_tools,
                    "only_paired_variant_calling": _truthy(params.get("only_paired_variant_calling")),
                },
            )

    def _reject_modes(tool: str, allowed: set[str], message: str, fix: str) -> None:
        if build_only_index or annotation_only:
            # Sarek performs no mode check for an empty build-only channel or
            # for VCF annotation input; global parameter guards below still run.
            return
        # A tool is accepted when at least one patient's mode is compatible.
        if not (present_modes & allowed):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PAIRING,
                message=message,
                fix=fix,
                details={"tool": tool, "present_modes": sorted(present_modes), "allowed_modes": sorted(allowed)},
            )

    def _require_paired(tool: str) -> None:
        _reject_modes(
            tool, {"somatic_paired"},
            f"Tool '{tool}' requires tumor/normal paired samples.",
            f"Provide both normal (status=0) and tumor (status=1) samples for the same patient when using {tool}.",
        )

    def _require_tumor_modes(tool: str) -> None:
        _reject_modes(
            tool, {"tumor_only", "somatic_paired"},
            f"Tool '{tool}' requires at least one tumor sample (status=1).",
            f"Provide a tumor sample when using {tool}, or remove it from --tools.",
        )

    def _require_tumor_only(tool: str) -> None:
        _reject_modes(
            tool, {"tumor_only"},
            f"Tool '{tool}' is tumor-only and requires at least one tumor sample without a matched normal.",
            f"Provide a tumor-only sample (no matched normal) when using {tool}.",
        )

    def _require_normal_sample(tool: str) -> None:
        _reject_modes(
            tool, {"germline", "somatic_paired"},
            f"Tool '{tool}' requires at least one normal/germline sample (status=0).",
            f"Provide a normal sample when using {tool}, or choose a tumor-compatible caller.",
        )

    def _require_germline_or_paired(tool: str) -> None:
        _reject_modes(
            tool, {"germline", "somatic_paired"},
            f"Tool '{tool}' supports germline or tumor/normal paired analyses only.",
            f"Use {tool} with normal/germline or paired tumor-normal samples.",
        )

    def _reject_if_paired_normal_is_only_input(tool: str, *, allow_tumor_only: bool = False) -> None:
        """Reject callers whose only possible input is removed by paired-only routing.

        Sarek filters matched normals out of ``cram_variant_calling_status_normal``
        when ``only_paired_variant_calling`` is enabled.  Germline-only callers
        therefore need an unpaired normal patient; mpileup may alternatively run
        on a tumor-only patient through its dedicated subworkflow.
        """
        if not paired_normal_germline_suppressed:
            return
        remaining_modes = {"germline", "tumor_only"} if allow_tumor_only else {"germline"}
        if present_modes & remaining_modes:
            return
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PAIRING,
            message=(
                f"Tool '{tool}' has no remaining input when --only_paired_variant_calling "
                "suppresses germline calling on matched normals."
            ),
            fix=(
                "Remove --only_paired_variant_calling, include an unpaired normal sample"
                + (" or a tumor-only sample" if allow_tumor_only else "")
                + f", or remove {tool} from --tools."
            ),
            details={
                "tool": tool,
                "present_modes": sorted(present_modes),
                "only_paired_variant_calling": True,
            },
        )

    def _require_sex(tool: str) -> None:
        """Mirror Sarek's global per-input-sample sex validation."""
        if build_only_index:
            return
        missing = []
        for patient, entries in rows_by_patient.items():
            for entry in entries:
                sex = str(entry.get("sex") or "NA")
                if sex == "NA":
                    missing.append({"patient": patient, "sample": entry.get("sample"), "line": entry.get("line")})
        if missing:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_SEX,
                message=f"Tool '{tool}' requires sex to be set for every sample.",
                fix="Set the samplesheet `sex` column to XX or XY for each row.",
                details={"tool": tool, "missing": missing},
            )

    def _require_snpeff_db(tool: str) -> None:
        if _resource_available(params, "snpeff_db", "snpeff"):
            if not params.get("snpeff_db"):
                warnings_acc.append(
                    f"{tool}: --snpeff_db not set; relying on the iGenomes genome config. "
                    "If SnpEff fails, pass --snpeff_db explicitly."
                )
            return
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_SNPEFF_DB,
            message=f"{tool} requires --snpeff_db when no supplying iGenomes genome is selected.",
            fix="Provide --snpeff_db matching the SnpEff cache database name, or select an iGenomes genome that defines it.",
            details={"tool": tool, "genome": _effective_genome(params) or None},
        )

    def _require_vep_cache(tool: str) -> None:
        fields = ("vep_cache_version", "vep_genome", "vep_species")
        inherited = [p for p in fields if not params.get(p)]
        missing = [p for p in fields if not _resource_available(params, p, "vep")]
        if not missing:
            if inherited:
                warnings_acc.append(
                    f"{tool}: VEP cache metadata (" + ", ".join(inherited) + ") not set; relying on the "
                    "iGenomes genome config. If VEP fails, pass them explicitly."
                )
            return
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_VEP_CACHE,
            message="VEP requires cache configuration.",
            fix="Provide: " + ", ".join(missing),
            details={"missing": missing, "tool": tool},
        )

    for tool in tools:
        if tool == "bbsplit":
            # BBSplit is invoked only by GATK preprocessing of mapping input.
            # PREPARE_GENOME otherwise leaves its input channel empty.
            if step == "mapping" and aligner != "parabricks":
                if not _explicit_resource_supplied(params, "bbsplit_index") and not _explicit_resource_supplied(params, "bbsplit_fasta_list"):
                    raise SkillError(
                        stage="preflight",
                        error_code=ErrorCode.MISSING_REFERENCE,
                        message="bbsplit requires an index or FASTA list for mapping preprocessing.",
                        fix="Provide --bbsplit_index, or provide --bbsplit_fasta_list so Sarek can build the index.",
                        details={"tool": tool, "step": step, "aligner": aligner},
                    )
        elif tool == "ngscheckmate":
            # CRAM_SAMPLEQC executes through variant_calling and passes this
            # BED directly into BAM_NGSCHECKMATE.
            if step != "annotate" and not _explicit_resource_supplied(params, "ngscheckmate_bed"):
                if _resource_available(params, "ngscheckmate_bed", "ngscheckmate"):
                    warnings_acc.append(
                        f"ngscheckmate: --ngscheckmate_bed not passed explicitly; relying on the "
                        f"{_effective_genome(params)} iGenomes bundle."
                    )
                else:
                    raise SkillError(
                        stage="preflight",
                        error_code=ErrorCode.MISSING_REFERENCE,
                        message="ngscheckmate requires --ngscheckmate_bed when no supplying iGenomes genome is selected.",
                        fix="Provide --ngscheckmate_bed, or select an iGenomes genome that includes the NGSCheckMate SNP BED.",
                        details={"tool": tool, "step": step, "genome": _effective_genome(params) or None},
                    )
        elif tool == "strelka":
            # The rendered usage matrix lists Strelka as somatic only, but the
            # executable 3.8.1 workflow invokes STRELKA_GERMLINE for normal
            # samples (and the upstream test profile uses that path). Follow
            # executable behaviour: germline and paired, never tumor-only.
            _require_germline_or_paired(tool)
        elif tool == "mutect2":
            _require_tumor_modes(tool)
            if not annotation_only and "somatic_paired" in present_modes and _truthy(params.get("no_intervals")):
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.INVALID_INTERVALS,
                    message="Mutect2 somatic tumor-normal calling cannot be combined with --no_intervals.",
                    fix="Enable intervals for paired somatic Mutect2 calling, or remove mutect2 from --tools.",
                    details={"tool": tool, "present_modes": sorted(present_modes), "no_intervals": True},
                )
            # Sarek warns for any Mutect2 run (paired or tumor-only) lacking a PON.
            if not _resource_available(params, "pon", "pon"):
                warnings_acc.append(
                    f"{tool} without --pon: no Panel-of-Normals was specified. It is highly "
                    "recommended to use one to filter recurrent technical artifacts."
                )
            elif (
                "1000g_pon.hg38.vcf.gz" in str(params.get("pon"))
                or (not params.get("pon") and _effective_genome(params) == "GATK.GRCh38")
            ):
                warnings_acc.append(
                    f"{tool}: the default GATK Panel-of-Normals (1000g_pon.hg38.vcf.gz) is in use. "
                    "It is highly recommended to generate a PON from normals technically similar to your tumors."
                )
            # Sarek also warns when no germline resource is provided (no filtering).
            if not _resource_available(params, "germline_resource", "germline_resource"):
                warnings_acc.append(
                    f"{tool} without --germline_resource: no germline-based filtering will be "
                    "applied. Providing one is recommended."
                )
        elif tool == "muse":
            _require_paired(tool)
        elif tool == "lofreq":
            _require_tumor_only(tool)
        elif tool == "sentieon_tnscope":
            _require_tumor_modes(tool)
        elif tool == "mpileup":
            # Sarek runs germline calling for every normal sample, including
            # the normal member of a tumor/normal cohort, and runs mpileup in
            # the tumor-only workflow as well. With paired-only filtering,
            # however, a purely paired cohort no longer feeds this caller.
            _reject_if_paired_normal_is_only_input(tool, allow_tumor_only=True)
            continue
        elif tool in {"deepvariant", "haplotypecaller", "sentieon_haplotyper", "sentieon_dnascope"}:
            _require_normal_sample(tool)
            _reject_if_paired_normal_is_only_input(tool)
        elif tool in {"manta", "tiddit"}:
            # Any mode OK.
            continue
        elif tool == "ascat":
            _require_sex(tool)
            _require_paired(tool)
            if params.get("ascat_purity") is not None and params.get("ascat_ploidy") is None:
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.INVALID_ASCAT_RESOURCES,
                    message="--ascat_purity requires --ascat_ploidy.",
                    fix="Provide --ascat-ploidy when overriding ASCAT purity, or omit --ascat-purity.",
                    details={"tool": tool, "ascat_purity": params.get("ascat_purity")},
                )
            # samplesheet_to_channel checks the reference archives; additionally,
            # the official schema marks ascat_genome as required to run ASCAT
            # and the ASCAT config passes it as genomeVersion.
            core_fields = ("ascat_genome", "ascat_alleles", "ascat_loci")
            missing_core = [p for p in core_fields if not _resource_available(params, p, "ascat")]
            inherited_core = [p for p in core_fields if not params.get(p) and p not in missing_core]
            if missing_core:
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.INVALID_ASCAT_RESOURCES,
                    message="ASCAT requires --ascat_genome, --ascat_alleles, and --ascat_loci.",
                    fix="Provide: " + ", ".join(missing_core),
                    details={"missing": missing_core, "tool": tool},
                )
            if inherited_core:
                warnings_acc.append(
                    f"ascat: relying on the {_effective_genome(params)} iGenomes bundle for "
                    + ", ".join("--" + p for p in inherited_core) + "."
                )
            if bool(params.get("wes")):
                warnings_acc.append(
                    "ASCAT WES/targeted runs should use custom allele/loci and LogR-correction "
                    "resources; the Sarek iGenomes ASCAT files are intended for WGS."
                )
        elif tool == "cnvkit":
            continue
        elif tool == "indexcov":
            _require_germline_or_paired(tool)
            if _truthy(params.get("wes")):
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.INVALID_TOOLS,
                    message="indexcov is WGS-only and is not executed when --wes is enabled.",
                    fix="Remove indexcov from --tools for WES/targeted runs, or remove --wes for whole-genome analysis.",
                    details={"tool": tool, "wes": True},
                )
        elif tool == "controlfreec":
            _require_sex(tool)
            _require_tumor_modes(tool)
            if annotation_only:
                continue
            fields = ("mappability", "chr_dir")
            missing = [p for p in fields if not _resource_available(params, p, "controlfreec")]
            inherited = [p for p in fields if not params.get(p) and p not in missing]
            if missing:
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.MISSING_REFERENCE,
                    message="controlfreec requires --mappability and --chr_dir.",
                    fix="Provide: " + ", ".join(missing) + " (or select an iGenomes --genome that supplies them).",
                    details={"missing": missing, "tool": tool},
                )
            if inherited:
                warnings_acc.append(
                    f"controlfreec: {', '.join(inherited)} not passed explicitly; relying on the "
                    "iGenomes genome config to supply them. If the controlfreec step fails, "
                    "provide --mappability / --chr_dir."
                )
        elif tool == "msisensor2":
            _require_tumor_only(tool)
            if annotation_only:
                continue
            if not _explicit_resource_supplied(params, "msisensor2_models"):
                if _resource_available(params, "msisensor2_models", "msisensor2"):
                    warnings_acc.append(
                        "msisensor2: --msisensor2_models not passed explicitly; relying on the "
                        "iGenomes genome config. If the msisensor2 step fails, provide "
                        "--msisensor2_models pointing at the trained model directory."
                    )
                else:
                    raise SkillError(
                        stage="preflight",
                        error_code=ErrorCode.MISSING_REFERENCE,
                        message="msisensor2 requires --msisensor2_models.",
                        fix="Provide --msisensor2_models pointing at the trained model directory.",
                        details={"tool": tool},
                    )
        elif tool == "msisensorpro":
            _require_paired(tool)
            if annotation_only:
                continue
            if not _resource_available(params, "msisensorpro_scan", "msisensorpro"):
                warnings_acc.append(
                    "msisensorpro: --msisensorpro_scan not provided; Sarek will generate the "
                    "scan from the reference FASTA during PREPARE_GENOME."
                )
            elif not _explicit_resource_supplied(params, "msisensorpro_scan"):
                warnings_acc.append(
                    "msisensorpro: --msisensorpro_scan not passed explicitly; relying on the "
                    f"{_effective_genome(params)} iGenomes bundle."
                )
        elif tool == "varlociraptor":
            _require_sex(tool)
            if present_modes & {"tumor_only", "somatic_paired"}:
                _enforce_varlociraptor_contamination(rows_by_patient, tools=tools)
            if not annotation_only and not build_only_index and not (
                set(tools) & _VARLOCIRAPTOR_CANDIDATE_CALLERS
            ):
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.INVALID_TOOLS,
                    message="Varlociraptor requires candidate variants from another variant caller.",
                    fix=(
                        "Add at least one compatible candidate caller to --tools, such as "
                        "haplotypecaller for germline or mutect2 for somatic calling."
                    ),
                    details={
                        "tool": tool,
                        "candidate_callers": sorted(_VARLOCIRAPTOR_CANDIDATE_CALLERS),
                    },
                )
        elif tool == "snpeff":
            if _annotation_step_active(params):
                _require_snpeff_db(tool)
        elif tool == "merge":
            # `merge` = combined SnpEff + VEP annotation, so both rule sets apply.
            if _annotation_step_active(params):
                _require_snpeff_db(tool)
                _require_vep_cache(tool)
        elif tool == "bcfann":
            missing = [
                p for p in ("bcftools_annotations", "bcftools_annotations_tbi", "bcftools_header_lines")
                if not _explicit_resource_supplied(params, p)
            ]
            if missing:
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.MISSING_REFERENCE,
                    message="bcfann requires bcftools annotation files and header lines.",
                    fix="Provide: " + ", ".join(missing),
                    details={"missing": missing, "tool": tool},
                )
        elif tool == "vep":
            _require_vep_cache(tool)

    # Require a Sentieon license only after every tool has cleared its mode/pairing
    # validation, so an invalid tool×mode combination (e.g. sentieon_tnscope in a
    # germline run) surfaces as INVALID_PAIRING rather than being masked by a missing
    # license. This keeps the check deterministic regardless of the host environment.
    sentieon_selected = sorted(set(tools) & _SENTIEON_TOOLS)
    if sentieon_selected:
        _require_sentieon_license("Sentieon tool(s): " + ", ".join(sentieon_selected))

    return warnings_acc


def _annotation_step_active(params: dict[str, Any]) -> bool:
    """Annotation can run from any Sarek start step when annotator tools are set."""
    step = str(params.get("step") or "").strip()
    return step in SUPPORTED_STEPS


def _enforce_varlociraptor_contamination(
    rows_by_patient: dict[str, list[dict[str, Any]]],
    *,
    tools: list[str],
) -> None:
    if not rows_by_patient:
        # The samplesheet builder already enforces this when given tools;
        # but if rows aren't surfaced here we cannot re-check.
        return
    missing: list[dict[str, Any]] = []
    for patient, entries in rows_by_patient.items():
        for entry in entries:
            try:
                status = int(entry.get("status", 0))
            except (TypeError, ValueError):
                status = 0
            if status != 1:
                continue
            contamination = entry.get("contamination")
            if contamination in (None, "", "NA"):
                missing.append({
                    "patient": patient,
                    "sample": entry.get("sample"),
                    "line": entry.get("line"),
                })
    if missing:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.VARLOCIRAPTOR_MISSING_CONTAMINATION,
            message="Varlociraptor requires `contamination` for every tumor row.",
            fix="Add a `contamination` value (float in [0, 1]) to each tumor row in the samplesheet.",
            details={"missing": missing, "tools": tools},
        )


# --- §5.4: Flag compatibility --------------------------------------------

def _truthy(value: Any) -> bool:
    if isinstance(value, bool):
        return value
    if value is None:
        return False
    if isinstance(value, (int, float)):
        return bool(value)
    if isinstance(value, str):
        return value.strip().lower() in {"true", "1", "yes", "y"}
    return bool(value)


def _uses_upstream_test_profile(params: dict[str, Any]) -> bool:
    profiles = {
        part.strip()
        for part in str(params.get("profile") or "").split(",")
        if part.strip()
    }
    return any(part == "test" or part.startswith("test_") for part in profiles)


def _check_demo_network(params: dict[str, Any]) -> None:
    """--demo composes nf-core's upstream ``-profile test``, whose inputs (test FASTQs
    and reference FASTA/GTF) are remote GitHub URLs by design. Under ``NXF_OFFLINE`` the
    nf-schema plugin aborts during parameter validation with a confusing
    ``does not exist`` before any work starts, so fail early with an actionable message
    instead. Downloading nf-core's PUBLIC test data does not conflict with the
    local-first guarantee, which governs USER genetic data (never uploaded)."""
    if not _uses_upstream_test_profile(params):
        return
    offline = str(os.environ.get("NXF_OFFLINE", "")).strip().lower()
    if offline in {"true", "1", "yes", "on"}:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.DEMO_REQUIRES_NETWORK,
            message=(
                "--demo runs the upstream nf-core `-profile test`, which downloads its test "
                "datasets and reference files from GitHub, but NXF_OFFLINE is set — Nextflow "
                "would abort during parameter validation because the remote test inputs cannot "
                "be reached."
            ),
            fix=(
                "Unset NXF_OFFLINE to run the demo (it fetches only nf-core's public test data — "
                "no user or genetic data is uploaded), or run a real analysis with your own local "
                "--input samplesheet and references, which is fully offline."
            ),
            details={"NXF_OFFLINE": os.environ.get("NXF_OFFLINE")},
        )


def _check_flag_compatibility(*, params: dict[str, Any], tools: list[str]) -> list[str]:
    warnings_acc: list[str] = []
    step = str(params.get("step") or "").strip() or "mapping"

    # 1. The executable Sarek check rejects Spark only when mapped BAM
    # publication is requested: both flags must be enabled together.
    spark_value = params.get("use_gatk_spark")
    spark_active = bool(spark_value)
    if spark_active and params.get("save_output_as_bam") and params.get("save_mapped"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.SPARK_OUTPUT_INCOMPATIBLE,
            message="`--use_gatk_spark` is incompatible with mapped BAM publication (`--save_mapped --save_output_as_bam`).",
            fix="Disable Spark, or drop either --save_mapped or --save_output_as_bam.",
            details={
                "use_gatk_spark": spark_value,
                "save_output_as_bam": params.get("save_output_as_bam"),
                "save_mapped": params.get("save_mapped"),
            },
        )

    # 2. MarkDuplicatesSpark + header/positional UMI dedup.
    # Sarek's exact rule: MarkDuplicatesSpark cannot do UMI-based deduplication
    # when umi_in_read_header OR umi_location is set. (umi_read_structure runs
    # fgbio consensus upstream, so it is NOT incompatible with Spark dedup.)
    spark_md = "markduplicates" in {t.strip() for t in re.split(r"[,\s]+", str(spark_value or "")) if t.strip()}
    if step == "mapping" and spark_md and (params.get("umi_in_read_header") or params.get("umi_location")):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.UMI_SPARK_INCOMPATIBLE,
            message="MarkDuplicatesSpark cannot perform UMI-based deduplication.",
            fix="Drop `--use_gatk_spark markduplicates`, or remove --umi_in_read_header / --umi_location.",
            details={
                "use_gatk_spark": spark_value,
                "umi_in_read_header": params.get("umi_in_read_header"),
                "umi_location": params.get("umi_location"),
            },
        )

    if step == "mapping" and (aligner := str(params.get("aligner") or "").strip().lower()):
        if aligner in {"parabricks", "sentieon-bwamem"} and params.get("umi_read_structure"):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"--aligner {aligner} cannot be combined with --umi_read_structure.",
                fix="Remove --umi_read_structure for parabricks/sentieon-bwamem, or choose a compatible aligner.",
                details={"aligner": aligner, "umi_read_structure": params.get("umi_read_structure")},
            )
        if aligner == "parabricks" and params.get("umi_in_read_header"):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message="--aligner parabricks cannot be combined with --umi_in_read_header.",
                fix="Remove --umi_in_read_header or choose a non-Parabricks aligner.",
                details={"aligner": aligner, "umi_in_read_header": params.get("umi_in_read_header")},
            )
        if aligner == "parabricks" and params.get("umi_location"):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message="--aligner parabricks cannot be combined with --umi_location.",
                fix="Remove --umi_location or choose a non-Parabricks aligner.",
                details={"aligner": aligner, "umi_location": params.get("umi_location")},
            )

    if step == "mapping" and params.get("umi_read_structure") and params.get("umi_location"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message="--umi_read_structure and --umi_location are mutually exclusive.",
            fix="Use either read-structure based UMI parsing or location/length based UMI parsing.",
            details={"umi_read_structure": params.get("umi_read_structure"), "umi_location": params.get("umi_location")},
        )
    if step == "mapping" and params.get("umi_in_read_header") and params.get("umi_location"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message="--umi_in_read_header and --umi_location are mutually exclusive.",
            fix="Use header-based UMIs or positional UMIs, not both.",
            details={"umi_in_read_header": params.get("umi_in_read_header"), "umi_location": params.get("umi_location")},
        )
    if step == "mapping" and params.get("umi_location") and params.get("umi_length") in (None, ""):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message="--umi_location requires --umi_length.",
            fix="Provide --umi_length when using --umi_location.",
            details={"umi_location": params.get("umi_location"), "umi_length": params.get("umi_length")},
        )
    if step == "mapping" and params.get("umi_in_read_header") and params.get("umi_read_structure"):
        if str(params.get("umi_read_structure")).strip() != "+T +T":
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message="--umi_in_read_header with --umi_read_structure requires '+T +T'.",
                fix="Set --umi_read_structure '+T +T' or remove --umi_in_read_header.",
                details={"umi_in_read_header": params.get("umi_in_read_header")},
            )

    # 2b. umi_read_structure must be valid fgbio read-structure syntax (sarek
    #     validates this early via readStructure(); fail at preflight, not late).
    rs_value = params.get("umi_read_structure")
    if rs_value:
        bad = [tok for tok in str(rs_value).split() if not _is_valid_read_structure_token(tok)]
        if bad or not str(rs_value).split():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FLAG_COMBINATION,
                message=f"--umi_read_structure has invalid fgbio read-structure token(s): {bad or [rs_value]}.",
                fix="Use fgbio read-structure syntax, e.g. '8M2S+T' or '+T +T' (codes T/B/M/S; '+' for the remainder).",
                details={"umi_read_structure": rs_value, "invalid_tokens": bad},
            )

    # 3. joint_germline requires one of the germline callers
    if _truthy(params.get("joint_germline")):
        germline_callers = {"haplotypecaller", "sentieon_haplotyper", "sentieon_dnascope"}
        if not (set(tools) & germline_callers):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.JOINT_GERMLINE_REQUIREMENT,
                message="`--joint_germline` requires at least one germline caller in --tools.",
                fix=f"Add one of: {', '.join(sorted(germline_callers))}.",
                details={"tools": tools, "required_any": sorted(germline_callers)},
            )
        for tool, emit_key in (
            ("sentieon_haplotyper", "sentieon_haplotyper_emit_mode"),
            ("sentieon_dnascope", "sentieon_dnascope_emit_mode"),
        ):
            if tool in tools and "gvcf" not in str(params.get(emit_key) or "").lower().split(","):
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.JOINT_GERMLINE_REQUIREMENT,
                    message=f"--joint_germline with {tool} requires a gvcf emit mode.",
                    fix=f"Set --{emit_key.replace('_', '-')} gvcf or a comma-separated mode containing gvcf.",
                    details={"tool": tool, "emit_mode": params.get(emit_key)},
                )
        missing_vqsr_inputs = [
            name for name in ("dbsnp", "known_indels", "known_snps")
            if not params.get(name)
        ]
        if _truthy(params.get("no_intervals")):
            missing_vqsr_inputs.append("intervals (--no_intervals is true)")
        if missing_vqsr_inputs and not _igenomes_active(params):
            warnings_acc.append(
                "joint_germline without all VQSR/interval inputs: missing "
                + ", ".join(missing_vqsr_inputs)
                + ". Sarek will not perform the corresponding variant recalibration/genotyping path."
            )

    # 4. joint_mutect2 requires mutect2
    if _truthy(params.get("joint_mutect2")) and "mutect2" not in tools:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.JOINT_MUTECT2_REQUIREMENT,
            message="`--joint_mutect2` requires `mutect2` in --tools.",
            fix="Add 'mutect2' to --tools or disable --joint_mutect2.",
            details={"tools": tools},
        )

    # 5. genome empty without igenomes_ignore (except annotate step)
    has_genome = _has_genome(params)
    igenomes_ignore = _truthy(params.get("igenomes_ignore"))
    upstream_test_profile = _uses_upstream_test_profile(params)
    if not has_genome and not igenomes_ignore and step != "annotate" and not upstream_test_profile:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.GENOME_REQUIRED,
            message="--genome is required unless --igenomes_ignore is set.",
            fix="Provide --genome (e.g. GATK.GRCh38) or set --igenomes_ignore true and supply --fasta.",
            details={"genome": params.get("genome"), "igenomes_ignore": igenomes_ignore},
        )

    # 6. genome empty without fasta (except annotate)
    has_fasta = _explicit_resource_supplied(params, "fasta")
    if not has_genome and not has_fasta and step != "annotate" and not upstream_test_profile:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.GENOME_REQUIRED,
            message="Either --genome or --fasta is required.",
            fix="Set --genome or provide a custom --fasta.",
            details={"genome": params.get("genome"), "fasta": params.get("fasta")},
        )

    # 7. build_only_index with --input set to a real samplesheet
    if _truthy(params.get("build_only_index")):
        input_value = params.get("input")
        if input_value not in (None, False, "", "false", "False"):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.BUILD_ONLY_INDEX_CONFLICT,
                message="`--build_only_index` cannot be combined with --input.",
                fix="Set --input false (or remove it) when only building indices.",
                details={"input": str(input_value)},
            )

    # 9. parabricks aligner + conda/mamba in profile
    aligner = str(params.get("aligner") or "").strip().lower()
    if aligner == "sentieon-bwamem":
        _require_sentieon_license("Sentieon aligner sentieon-bwamem")
    profile_tokens = set(_profile_tokens(str(params.get("profile") or "")))
    if aligner == "parabricks" and (profile_tokens & {"conda", "mamba"}):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.PARABRICKS_CONDA_CONFLICT,
            message="`--aligner parabricks` cannot be combined with conda/mamba profiles.",
            fix="Use a container profile (docker, singularity, apptainer) with parabricks.",
            details={"aligner": aligner, "profile_tokens": sorted(profile_tokens)},
        )

    # 10. vep_dbnsfp truthy without dbnsfp + dbnsfp_tbi
    if _truthy(params.get("vep_dbnsfp")) and not (
        _explicit_resource_supplied(params, "dbnsfp")
        and _explicit_resource_supplied(params, "dbnsfp_tbi")
    ):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_DBNSFP,
            message="--vep_dbnsfp requires both --dbnsfp and --dbnsfp_tbi.",
            fix="Provide --dbnsfp and --dbnsfp_tbi, or disable --vep_dbnsfp.",
            details={
                "dbnsfp": params.get("dbnsfp"),
                "dbnsfp_tbi": params.get("dbnsfp_tbi"),
            },
        )

    for flag, required in (
        ("vep_condel", ("condel_config",)),
        ("vep_mastermind", ("mastermind_file",)),
        ("vep_spliceai", ("spliceai_snv", "spliceai_snv_tbi", "spliceai_indel", "spliceai_indel_tbi")),
    ):
        if _truthy(params.get(flag)):
            missing = [name for name in required if not _explicit_resource_supplied(params, name)]
            if missing:
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.MISSING_REFERENCE,
                    message=f"--{flag.replace('_', '-')} requires documented VEP plugin resource path(s).",
                    fix="Provide: " + ", ".join(missing),
                    details={"flag": flag, "missing": missing},
                )

    # Phenotypes is different from the other VEP plugins: the official usage
    # docs mark --phenotypes_file as optional because VEP can download phenotype
    # data automatically. Only require an index when the user explicitly supplies
    # a gzipped local phenotypes file.
    phenotypes_file = str(params.get("phenotypes_file") or "").strip()
    if (
        _truthy(params.get("vep_phenotypes"))
        and _explicit_resource_supplied(params, "phenotypes_file")
        and phenotypes_file.endswith(".gz")
        and not _explicit_resource_supplied(params, "phenotypes_file_tbi")
    ):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_REFERENCE,
            message="--phenotypes_file_tbi is required when --phenotypes_file is gzipped.",
            fix="Provide the matching .tbi/.csi index, or omit --phenotypes_file and let VEP download phenotype data.",
            details={"flag": "vep_phenotypes", "phenotypes_file": phenotypes_file},
        )

    # 11. vep_loftee truthy → WARN
    if _truthy(params.get("vep_loftee")):
        warnings_acc.append(
            "--vep_loftee requires a compatible VEP cache (LoFtee plugin); verify the cache supports it."
        )

    # POST_VARIANTCALLING chooses Varlociraptor before the bcftools
    # filter/normalization/consensus/concatenation branch. Accepting both would
    # silently record requested transformations that upstream never executes.
    post_variant_flags = (
        "filter_vcfs", "normalize_vcfs", "snv_consensus_calling", "concatenate_vcfs",
    )
    conflicting_post_variant_flags = [
        flag for flag in post_variant_flags if _truthy(params.get(flag))
    ]
    if "varlociraptor" in tools and conflicting_post_variant_flags:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=(
                "--tools varlociraptor cannot be combined with bcftools post-variant "
                "processing flags because Sarek executes only the Varlociraptor branch."
            ),
            fix=(
                "Remove varlociraptor, or remove: "
                + ", ".join("--" + flag.replace("_", "-") for flag in conflicting_post_variant_flags)
                + "."
            ),
            details={"tool": "varlociraptor", "conflicting_flags": conflicting_post_variant_flags},
        )

    # 12. filter_vcfs + explicitly malformed bcftools_filter_criteria.
    # If omitted, Sarek supplies its schema default of `-f PASS,.`.
    if _truthy(params.get("filter_vcfs")):
        criteria = params.get("bcftools_filter_criteria")
        if criteria not in (None, False) and not _is_valid_bcftools_criteria(str(criteria)):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_BCFTOOLS_FILTER,
                message="An explicitly supplied --bcftools_filter_criteria value must be a non-empty bcftools view argument string.",
                fix="Omit the parameter to use Sarek's '-f PASS,.' default, or provide valid bcftools view criteria.",
                details={"bcftools_filter_criteria": criteria},
            )

    bcf_transform_flags = [
        flag for flag in ("filter_vcfs", "normalize_vcfs") if _truthy(params.get(flag))
    ]
    if bcf_transform_flags and not (set(tools) & _SNV_CALLERS):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_FLAG_COMBINATION,
            message=(
                "The requested bcftools variant transformation has no eligible "
                "small-variant caller output."
            ),
            fix=(
                "Select at least one supported small-variant caller ("
                + ", ".join(sorted(_SNV_CALLERS))
                + "), or remove: "
                + ", ".join("--" + flag.replace("_", "-") for flag in bcf_transform_flags)
                + "."
            ),
            details={"transform_flags": bcf_transform_flags, "tools": tools},
        )

    # 13/14. snv_consensus_calling
    if _truthy(params.get("snv_consensus_calling")):
        if not _truthy(params.get("normalize_vcfs")):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_CONSENSUS_CALLING,
                message="--snv_consensus_calling requires --normalize_vcfs.",
                fix="Enable --normalize_vcfs before requesting SNV consensus calling.",
                details={"snv_consensus_calling": True, "normalize_vcfs": params.get("normalize_vcfs")},
            )
        snv_in_tools = sorted({t for t in tools if t in _SNV_CALLERS})
        if len(snv_in_tools) < 2:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_CONSENSUS_CALLING,
                message="--snv_consensus_calling requires at least 2 SNV callers in --tools.",
                fix=f"Add at least 2 of: {', '.join(sorted(_SNV_CALLERS))}.",
                details={"snv_callers_in_tools": snv_in_tools, "snv_callers": sorted(_SNV_CALLERS)},
            )
        try:
            min_count = int(params.get("consensus_min_count") or 0)
        except (TypeError, ValueError):
            min_count = 0
        if min_count > len(snv_in_tools):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_CONSENSUS_CALLING,
                message="--consensus_min_count exceeds the number of SNV callers in --tools.",
                fix=f"Set --consensus_min_count <= {len(snv_in_tools)}.",
                details={
                    "consensus_min_count": min_count,
                    "snv_callers_in_tools": snv_in_tools,
                },
            )

    # 16. snpsift in tools + missing databases
    if "snpsift" in tools and not _explicit_resource_supplied(params, "snpsift_databases"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_SNPSIFT_DATABASES,
            message="`snpsift` in --tools requires --snpsift_databases.",
            fix="Provide --snpsift_databases pointing to your SnpSift annotation database(s).",
            details={"tools": tools},
        )
    if "snpsift" in tools and _explicit_resource_supplied(params, "snpsift_databases"):
        _validate_snpsift_databases(str(params["snpsift_databases"]))

    # 17. no_intervals True + intervals non-empty
    if _truthy(params.get("no_intervals")) and _explicit_resource_supplied(params, "intervals"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_INTERVALS,
            message="`--no_intervals true` is incompatible with a non-empty --intervals path.",
            fix="Either drop --intervals or set --no_intervals false.",
            details={"no_intervals": True, "intervals": params.get("intervals")},
        )

    # 18. no_intervals True + wes True → WARN
    if _truthy(params.get("no_intervals")) and _truthy(params.get("wes")):
        warnings_acc.append(
            "--no_intervals true with --wes true bypasses target-region restriction; "
            "results may include off-target regions."
        )

    intervals = str(params.get("intervals") or "").strip() if _explicit_resource_supplied(params, "intervals") else ""
    if intervals and not intervals.lower().endswith((".bed", ".interval_list")):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_INTERVALS,
            message="--intervals must end with .bed or .interval_list.",
            fix="Provide an intervals file accepted by the nf-core/sarek 3.8.1 parameter schema: .bed or .interval_list.",
            details={"intervals": intervals},
        )
    if _truthy(params.get("wes")) and intervals and not intervals.lower().endswith(".bed"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_INTERVALS,
            message="--wes requires --intervals to be a BED file.",
            fix="Use a .bed intervals file for WES/targeted Sarek runs.",
            details={"intervals": intervals, "wes": True},
        )

    bqsr_steps = {"mapping", "markduplicates", "prepare_recalibration", "recalibrate"}
    skip_tools = set(_normalize_tools(params.get("skip_tools")))
    if step in bqsr_steps and "baserecalibrator" not in skip_tools and not upstream_test_profile:
        if not any(_resource_available(params, name, "bqsr") for name in ("dbsnp", "known_indels")):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.MISSING_REFERENCE,
                message=(
                    "BQSR requires --dbsnp or --known_indels when baserecalibrator is not "
                    "skipped. Note: baserecalibrator runs by default in the Sarek mapping "
                    "workflow, independently of the requested variant-calling --tools."
                ),
                fix="Provide --dbsnp/--known_indels resources, select an iGenomes genome, or add baserecalibrator to --skip_tools.",
                details={"step": step, "skip_tools": sorted(skip_tools)},
            )

    germline_filter_callers = {"haplotypecaller", "sentieon_haplotyper", "sentieon_dnascope"}
    if (
        set(tools) & germline_filter_callers
        and not any(_resource_available(params, name, "bqsr") for name in ("dbsnp", "known_indels"))
    ):
        warnings_acc.append(
            "HaplotypeCaller/Sentieon germline calling without --dbsnp or --known_indels: "
            "Sarek will not perform the corresponding variant filtering."
        )

    # DragMap + baserecalibrator: sarek recommends skipping BQSR with DragMap.
    aligner = str(params.get("aligner") or "").strip().lower()
    if step == "mapping" and "dragmap" in aligner and "baserecalibrator" not in skip_tools:
        warnings_acc.append(
            "DragMap was specified as aligner but baserecalibrator is not in --skip_tools. "
            "It is recommended to skip base recalibration when using DragMap."
        )

    return warnings_acc


def _is_valid_bcftools_criteria(criteria: str) -> bool:
    """A non-empty bcftools expression; detailed parsing is delegated to bcftools."""
    text = criteria.strip()
    if not text:
        return False
    return "\n" not in text and "\x00" not in text


def _validate_snpsift_databases(path_text: str) -> None:
    """Validate the SnpSift database CSV shape from Sarek's auxiliary schema."""
    if _is_uri(path_text):
        return
    path = Path(os.path.expanduser(path_text))
    if not path.exists():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_SNPSIFT_DATABASES,
            message="--snpsift_databases file does not exist.",
            fix="Provide a readable CSV with vcf,tbi,fields,prefix,vardb columns.",
            details={"path": path_text},
        )
    try:
        with path.open(newline="", encoding="utf-8-sig") as handle:
            reader = csv.DictReader(handle)
            fieldnames = reader.fieldnames or []
            if "vcf" not in fieldnames:
                raise ValueError("missing required column 'vcf'")
            for line_number, row in enumerate(reader, start=2):
                vcf = (row.get("vcf") or "").strip()
                tbi = (row.get("tbi") or "").strip()
                fields = (row.get("fields") or "").strip()
                vardb = (row.get("vardb") or "").strip()
                if not re.match(r"^\S+\.vcf(\.gz)?$", vcf):
                    raise ValueError(f"line {line_number}: vcf must end with .vcf or .vcf.gz and contain no whitespace")
                if tbi and not re.match(r"^\S+\.tbi$", tbi):
                    raise ValueError(f"line {line_number}: tbi must end with .tbi and contain no whitespace")
                if not fields and not vardb:
                    raise ValueError(f"line {line_number}: fields is required when vardb is not provided")
                if vardb and any(ch.isspace() for ch in vardb):
                    raise ValueError(f"line {line_number}: vardb must not contain whitespace")
    except (OSError, csv.Error, ValueError) as exc:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_SNPSIFT_DATABASES,
            message="--snpsift_databases is not compatible with Sarek's CSV schema.",
            fix="Use columns vcf,tbi,fields,prefix,vardb; vcf is required and fields is required when vardb is empty.",
            details={"path": path_text, "error": str(exc)},
        ) from exc


# --- §5.5: Reference path existence --------------------------------------

def _is_uri(value: str) -> bool:
    return any(value.startswith(prefix) for prefix in _URI_PREFIXES)


def _glob_reference_pattern(pattern: str) -> list[str]:
    """Resolve basic Nextflow-style glob patterns, including `{a,b}` groups."""
    expanded = [pattern]
    brace_re = re.compile(r"\{([^{}]+)\}")
    while any(brace_re.search(item) for item in expanded):
        next_expanded: list[str] = []
        for item in expanded:
            match = brace_re.search(item)
            if match is None:
                next_expanded.append(item)
                continue
            for alternative in match.group(1).split(","):
                next_expanded.append(item[:match.start()] + alternative + item[match.end():])
        expanded = next_expanded
    return sorted({match for item in expanded for match in glob.glob(item)})


# Documented filename patterns (nextflow_schema.json) for reference paths that
# declare one. Mismatches are WARNINGS, not errors: a valid file with a
# non-standard name should not block the run, but a likely mistake is surfaced.
_REFERENCE_PATH_PATTERNS = {
    "pon": re.compile(r"\.vcf\.gz$"),
    "pon_tbi": re.compile(r"\.vcf\.gz\.tbi$"),
    "dbsnp": re.compile(r"\.vcf\.gz$"),
    "dbsnp_tbi": re.compile(r"\.vcf\.gz\.tbi$"),
    "germline_resource": re.compile(r"\.vcf\.gz$"),
    "germline_resource_tbi": re.compile(r"\.vcf\.gz\.tbi$"),
    "known_snps": re.compile(r"\.vcf\.gz$"),
    "known_snps_tbi": re.compile(r"\.vcf\.gz\.tbi$"),
    "ascat_alleles": re.compile(r"\.zip$"),
    "ascat_loci": re.compile(r"\.zip$"),
    "ascat_loci_gc": re.compile(r"\.zip$"),
    "ascat_loci_rt": re.compile(r"\.zip$"),
    "cf_chrom_len": re.compile(r"\.(fai|len)$"),
    "cnvkit_reference": re.compile(r"\.cnn$"),
    "dbnsfp": re.compile(r"\.gz$"),
    "dbnsfp_tbi": re.compile(r"\.gz\.(csi|tbi)$"),
    "mappability": re.compile(r"\.gem$"),
    "mastermind_file": re.compile(r"\.vcf\.gz$"),
    "ngscheckmate_bed": re.compile(r"\.bed$"),
    "phenotypes_file": re.compile(r"\.(gff|gvf)(\.gz)?$"),
    "phenotypes_file_tbi": re.compile(r"\.(gff|gvf)\.gz\.(csi|tbi)$"),
    "sentieon_dnascope_model": re.compile(r"\.model$"),
    "spliceai_snv": re.compile(r"\.vcf\.gz$"),
    "spliceai_snv_tbi": re.compile(r"\.vcf\.gz\.(csi|tbi)$"),
    "spliceai_indel": re.compile(r"\.vcf\.gz$"),
    "spliceai_indel_tbi": re.compile(r"\.vcf\.gz\.(csi|tbi)$"),
    "bcftools_annotations": re.compile(r"\.vcf\.gz$"),
    "bcftools_annotations_tbi": re.compile(r"\.vcf\.gz\.tbi$"),
    "snpsift_databases": re.compile(r"\.csv$"),
    "intervals": re.compile(r"\.(bed|interval_list)$"),
    "dict": re.compile(r"\.dict$"),
    "fasta": re.compile(r"\.fn?a(sta)?(\.gz)?$"),
}


def _check_reference_path_suffixes(params: dict[str, Any]) -> list[str]:
    """Warn when a reference path's filename does not match its documented suffix."""
    warnings_acc: list[str] = []
    for name, pattern in _REFERENCE_PATH_PATTERNS.items():
        raw = params.get(name)
        if raw is None or _explicitly_disabled(params, name):
            continue
        text = str(raw).strip()
        if not text:
            continue
        if any(ch.isspace() for ch in text):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.MISSING_REFERENCE,
                message=f"--{name} value contains whitespace, which does not match the nf-core/sarek schema.",
                fix=f"Move or rename the --{name} path so the complete value contains no spaces.",
                details={"param": name, "path": text},
            )
        basename = text.rsplit("/", 1)[-1]
        if not pattern.search(basename):
            warnings_acc.append(
                f"--{name} value '{basename}' does not match the expected pattern "
                f"'{pattern.pattern}'; verify it is the correct file type."
            )
    return warnings_acc


def _check_reference_paths(params: dict[str, Any]) -> None:
    """Every set, non-URI reference path must exist on disk."""
    for name in REFERENCE_PATH_PARAMS:
        raw = params.get(name)
        if raw is None or _explicitly_disabled(params, name):
            continue
        text = str(raw).strip()
        if not text:
            continue
        if _is_uri(text):
            continue
        if name in _REFERENCE_PATH_PATTERN_PARAMS and (glob.has_magic(text) or ("{" in text and "}" in text)):
            if _glob_reference_pattern(text):
                continue
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.MISSING_REFERENCE,
                message=f"Reference pattern for --{name} did not match any files.",
                fix=f"Correct --{name} or remove it.",
                details={"param": name, "path_pattern": text},
            )
        path = Path(os.path.expanduser(text))
        if not path.exists():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.MISSING_REFERENCE,
                message=f"Reference path for --{name} does not exist.",
                fix=f"Correct --{name} or remove it.",
                details={"param": name, "path": text},
            )


# --- §5.6: Annotation cache ----------------------------------------------

def _check_annotation_cache(*, params: dict[str, Any], output_dir: Path) -> list[str]:
    warnings_acc: list[str] = []

    if _truthy(params.get("download_cache")):
        target = params.get("outdir_cache") or str((output_dir / "cache").as_posix())
        target_text = str(target)
        if not _is_uri(target_text):
            target_path = Path(os.path.expanduser(target_text)).resolve()
            parent = target_path if target_path.exists() else target_path.parent
            if not parent.exists() or not os.access(parent, os.W_OK):
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.MISSING_ANNOTATION_CACHE,
                    message="Annotation cache download target is not writable.",
                    fix="Provide a writable --outdir_cache or ensure the output directory is writable.",
                    details={"target": str(target_path)},
                )

    # SnpEff cache layout check
    snpeff_cache = params.get("snpeff_cache")
    snpeff_db = params.get("snpeff_db")
    if snpeff_cache and not _is_uri(str(snpeff_cache)):
        cache_root = Path(os.path.expanduser(str(snpeff_cache)))
        if cache_root.exists() and snpeff_db:
            child = cache_root / str(snpeff_db)
            if not child.exists():
                warnings_acc.append(
                    f"SnpEff cache root '{cache_root}' exists but expected subdir '{child.name}' was not found; "
                    "verify the snpeff_db naming."
                )

    # VEP cache layout check
    vep_cache = params.get("vep_cache")
    vep_species = params.get("vep_species")
    vep_cache_version = params.get("vep_cache_version")
    vep_genome = params.get("vep_genome")
    if vep_cache and not _is_uri(str(vep_cache)):
        cache_root = Path(os.path.expanduser(str(vep_cache)))
        if cache_root.exists() and vep_species and vep_cache_version and vep_genome:
            custom_args = str(params.get("vep_custom_args") or "")
            species_suffix = "_merged" if "--merged" in custom_args else ("_refseq" if "--refseq" in custom_args else "")
            # The cache initialisation subworkflow and usage guide use:
            # ${vep_species}${suffix}/${vep_cache_version}_${vep_genome}.
            expected = cache_root / f"{vep_species}{species_suffix}" / f"{vep_cache_version}_{vep_genome}"
            if not expected.exists():
                warnings_acc.append(
                    f"VEP cache root '{cache_root}' exists but expected subdir '"
                    f"{vep_species}{species_suffix}/{vep_cache_version}_{vep_genome}' was not found; verify cache layout."
                )

    return warnings_acc


# --- §5.1: Resume manifest drift -----------------------------------------

_RESUME_TRACKED_FIELDS = (
    "step",
    "aligner",
    "analysis_mode",
    "joint_germline",
    "joint_mutect2",
    "wes",
    # "arm", "gpu", "spark" omitted: they are always appended to "profile"
    # (arm64/gpu/spark tokens), so the profile comparison covers them redundantly.
    # Including them here would cause false-positive drift because none of these
    # flags are written to params.yaml (they are wrapper-level profile modifiers).
    "profile",
    # params_checksum and reference_checksums are injected into params_for_preflight
    # by the orchestrator (nfcore_sarek_wrapper.py) before calling run_preflight on
    # resume. They are not present on non-resume runs and are never written to
    # params.yaml — they exist only for drift detection.
    "params_checksum",
    "reference_checksums",
    "samplesheet_checksum",
    "work_dir",
)


def _normalize_set(value: Any) -> list[str]:
    if value is None or value is False:
        return []
    if isinstance(value, str):
        return sorted({t.strip().lower() for t in value.split(",") if t.strip()})
    if isinstance(value, Iterable):
        return sorted({str(t).strip().lower() for t in value if str(t).strip()})
    return [str(value).strip().lower()]


def _check_resume_drift(
    *,
    params: dict[str, Any],
    samplesheet: dict[str, Any],
    pipeline_source: dict[str, Any],
    manifest: dict[str, Any],
) -> None:
    """Compare current run against resume manifest; raise on drift."""
    current: dict[str, Any] = {}
    current["step"] = params.get("step") or DEFAULT_STEP
    current["aligner"] = params.get("aligner") or DEFAULT_ALIGNER
    current["tools"] = _normalize_set(params.get("tools"))
    current["skip_tools"] = _normalize_set(params.get("skip_tools"))
    current["analysis_mode"] = samplesheet.get("analysis_mode")
    current["joint_germline"] = _truthy(params.get("joint_germline"))
    current["joint_mutect2"] = _truthy(params.get("joint_mutect2"))
    current["wes"] = _truthy(params.get("wes"))
    current["profile"] = params.get("profile")
    current["work_dir"] = params.get("work_dir") or "upstream/work"
    current["params_checksum"] = params.get("params_checksum")
    current["reference_checksums"] = params.get("reference_checksums")
    _normalized = samplesheet.get("normalized_path")
    current["samplesheet_checksum"] = None
    if _normalized:
        try:
            current["samplesheet_checksum"] = "sha256:" + hashlib.sha256(
                Path(str(_normalized)).read_bytes()
            ).hexdigest()
        except OSError:
            pass

    diffs: dict[str, dict[str, Any]] = {}

    for key in _RESUME_TRACKED_FIELDS:
        expected = current.get(key)
        observed = manifest.get(key)
        # Legacy manifests: missing field means we don't drift on it.
        if key not in manifest:
            continue
        # Booleans normalize through _truthy on both sides if the type differs.
        if isinstance(expected, bool):
            observed = _truthy(observed)
        if expected != observed:
            diffs[key] = {"expected": expected, "observed": observed}

    # Set-valued fields
    for key in ("tools", "skip_tools"):
        if key not in manifest:
            continue
        observed_set = _normalize_set(manifest.get(key))
        if current[key] != observed_set:
            diffs[key] = {"expected": current[key], "observed": observed_set}

    # pipeline_source: only consider keys that exist in the legacy manifest entry.
    manifest_source = manifest.get("pipeline_source") or {}
    current_source = {
        "source_kind": pipeline_source.get("source_kind"),
        "resolved_version": pipeline_source.get("resolved_version"),
    }
    for k, v in current_source.items():
        if k in manifest_source and manifest_source.get(k) != v:
            diffs[f"pipeline_source.{k}"] = {
                "expected": v,
                "observed": manifest_source.get(k),
            }

    if diffs:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.RESUME_DRIFT,
            message="Resume manifest is incompatible with the current run.",
            fix="Use a fresh output directory, or re-run with the same parameters as the previous run.",
            details={"diffs": diffs},
        )
