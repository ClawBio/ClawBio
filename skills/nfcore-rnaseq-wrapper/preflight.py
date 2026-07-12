from __future__ import annotations

import gzip
import hashlib
import json
import os
import re
import shutil
import subprocess
import sys
import warnings
from pathlib import Path

import yaml

from clawbio.common.checksums import sha256_file

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "schemas")

from errors import ErrorCode, SkillError
from schemas import (
    DEFAULT_PIPELINE_VERSION,
    DEFAULT_REMOTE_PIPELINE,
    DEFAULT_TRIMMER,
    FASTQ_BASENAME_RE as _FASTQ_BASENAME_RE,
    JAVA_MIN_VERSION,
    NEXTFLOW_MIN_VERSION,
    PROJECT_ROOT,
    is_under_tmp,
    SUPPORTED_ALIGNERS,
    SUPPORTED_IGENOMES_NAMES,
    SUPPORTED_PSEUDO_ALIGNERS,
    SUPPORTED_PUBLISH_DIR_MODES,
    SUPPORTED_RIBO_TOOLS,
    SUPPORTED_RSEQC_MODULES,
    SUPPORTED_SALMON_LIB_TYPES,
    SUPPORTED_STRANDEDNESS,
    SUPPORTED_TRIMMERS,
    SUPPORTED_UMI_EXTRACT_METHODS,
    SUPPORTED_UMI_GROUPING_METHODS,
    SUPPORTED_UMI_TOOLS,
    WHITESPACE_RE as _WHITESPACE_RE,
)


_SUBPROCESS_TIMEOUT = 30
# _WHITESPACE_RE is imported from schemas (single source shared with params_builder).
_EMAIL_RE = re.compile(r"^([a-zA-Z0-9_\-\.]+)@([a-zA-Z0-9_\-\.]+)\.([a-zA-Z]{2,5})$")
_GENOME_ID_RE = re.compile(r"^[a-zA-Z0-9_\-\.]+$")
# nf-core/rnaseq accepts .fa, .fasta, .fna (optionally .gz). The earlier pattern
# also matched non-extensions like ".fnasta"; tightened to the canonical set (audit F-06).
_FASTA_EXT_RE = re.compile(r"^\S+\.(fa|fasta|fna)(\.gz)?$")
_GTF_EXT_RE = re.compile(r"^\S+\.gtf(\.gz)?$")
_GFF_EXT_RE = re.compile(r"^\S+\.gff3?(\.gz)?$")
_BED_EXT_RE = re.compile(r"^\S+\.bed(\.gz)?$")
_HISAT2_MEM_RE = re.compile(r"^\d+(\.\d+)?\.?\s*(K|M|G|T)?B$")
# _FASTQ_BASENAME_RE is the shared regex imported from schemas (single source of
# truth, also used by samplesheet_builder); it must NOT be redefined locally.
_REFERENCE_GROUPS = (("genome",), ("fasta", "gtf"), ("fasta", "gff"))
_ANNOTATION_REFERENCE_FIELDS = (
    "fasta",
    "gtf",
    "gff",
    "transcript_fasta",
    "additional_fasta",
    "gene_bed",
    "splicesites",
)
_INDEX_REFERENCE_FIELDS = (
    "star_index",
    "rsem_index",
    "hisat2_index",
    "bowtie2_index",
    "salmon_index",
    "kallisto_index",
)
_EXPLICIT_REFERENCE_FIELDS = (*_ANNOTATION_REFERENCE_FIELDS, *_INDEX_REFERENCE_FIELDS)
# Reference fields that genuinely conflict with --genome because they provide an
# alternative *genome sequence* source. nf-core/rnaseq supports iGenomes together
# with additive annotation/transcriptome overrides (e.g. --genome GRCh38 --gtf
# custom.gtf, or --additional-fasta spike_ins.fa, --transcript-fasta, --gene-bed,
# --splicesites, --salmon-index/--kallisto-index), but providing a second genome
# FASTA (--fasta) or a genome-level index (STAR/RSEM/HISAT2/Bowtie2) alongside a
# named iGenomes genome is ambiguous and is rejected (audit M1).
_GENOME_CONFLICTING_REFERENCE_FIELDS = (
    "fasta",
    "star_index",
    "rsem_index",
    "hisat2_index",
    "bowtie2_index",
)
# sortmerna_index is an rRNA filter database, not a genome index — it is compatible with
# --genome and must not participate in the genome-vs-explicit-reference conflict check.
_CONTAMINANT_PATH_FIELDS = ("kraken_db", "sylph_db", "sylph_taxonomy", "bbsplit_fasta_list", "bbsplit_index", "sortmerna_index")
_COMMA_SEPARATED_CONTAMINANT_FIELDS = frozenset({"sylph_db", "sylph_taxonomy"})
# Catch every Groovy form that assigns into `params`: block (`params {`),
# property (`params.x` / `params.put(...)`), assignment (`params = ...`),
# subscript (`params['x']` / `params ['x']`), and map-merge (`params << [...]`).
# The earlier pattern missed subscript/merge forms, leaving a bypass (audit M2).
_NEXTFLOW_PARAMS_OVERRIDE_RE = re.compile(r"^\s*params\s*(?:\{|\.|=|\[|<<)")
# The documented escape hatch: a custom genome catalogue under `params.genomes`.
# Any `params.genomes...` line (block, property, or subscript) is allowed.
_NEXTFLOW_ALLOWED_PARAMS_RE = re.compile(r"^\s*params\.genomes\b")
# includeConfig directive — its target is audited recursively so a params
# override cannot hide in an included file (audit M2).
_NEXTFLOW_INCLUDE_RE = re.compile(r"""^\s*includeConfig\s+['"]([^'"]+)['"]""")
_MAX_CONFIG_INCLUDE_DEPTH = 16
# Root-level entries that must not trip OUTPUT_DIR_NOT_EMPTY. `reproducibility` is the
# wrapper's own bundle directory (commands.sh, checksums.sha256, environment.yml,
# remap_paths.py, params/samplesheet snapshots, provenance JSON) — never user data — so
# its entire contents are ignored. Matches nfcore-sarek / nfcore-scrnaseq so the three
# wrappers accept the same set of pre-existing output-dir contents, and avoids a fragile
# per-file allowlist that silently rejected a copied bundle (OUTPUT_DIR_NOT_EMPTY).
_IGNORED_ROOT_NAMES = frozenset(
    {"reproducibility", ".DS_Store", ".gitkeep", ".gitignore", "Thumbs.db", "check_result.json"}
)
_GENCODE_GTF_MARKER_RE = re.compile(r'(^|[;\s])(gene_type|havana_gene)\s+"')
_GENCODE_GTF_MAX_FEATURE_RECORDS = 10


def serialize_params_for_checksum(params_payload: dict[str, object]) -> str:
    return yaml.dump(params_payload, allow_unicode=True, sort_keys=False)


def params_payload_checksum(params_payload: dict[str, object]) -> str:
    payload = serialize_params_for_checksum(params_payload).encode("utf-8")
    return "sha256:" + hashlib.sha256(payload).hexdigest()


def check_resume_params_checksum(params_payload: dict[str, object], output_dir: Path) -> str:
    output_dir = output_dir.expanduser().resolve()
    repro_dir = output_dir / "reproducibility"
    manifest_path = repro_dir / "manifest.json"
    current_checksum = params_payload_checksum(params_payload)

    manifest: dict[str, object] = {}
    if manifest_path.exists():
        try:
            manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
        except json.JSONDecodeError as exc:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_RESUME_STATE,
                message="Resume manifest is not valid JSON.",
                fix="Use a valid output directory from a previous wrapper run.",
                details={"manifest": str(manifest_path), "error": str(exc)},
            ) from exc

    manifest_checksum = manifest.get("params_checksum")
    if manifest_checksum:
        if manifest_checksum != current_checksum:
            _raise_resume_mismatch("params_checksum", current_checksum, manifest_checksum)
        return current_checksum

    # Legacy manifest without params_checksum: the old fallback compared file-SHA256
    # against payload-SHA256 — two incompatible formats that could never match.
    # Skip validation and warn instead of blocking a valid resume.
    warnings.warn(
        "Legacy manifest: params_checksum absent — parameter change detection skipped for this resume.",
        UserWarning,
        stacklevel=2,
    )
    return current_checksum


def check_resume_samplesheet_checksum(samplesheet_path: Path, output_dir: Path) -> str:
    output_dir = output_dir.expanduser().resolve()
    manifest_path = output_dir / "reproducibility" / "manifest.json"
    if not manifest_path.exists():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_RESUME_STATE,
            message="--resume requires an existing reproducibility manifest.",
            fix="Run without --resume in a new output directory, or resume a completed compatible run.",
            details={"manifest": str(manifest_path)},
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_RESUME_STATE,
            message="Resume manifest is not valid JSON.",
            fix="Use a valid output directory from a previous wrapper run.",
            details={"manifest": str(manifest_path), "error": str(exc)},
        ) from exc
    manifest_checksum = manifest.get("samplesheet_checksum")
    if not manifest_checksum:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_RESUME_STATE,
            message="Resume manifest is missing the samplesheet checksum.",
            fix="Resume from an output directory created by this wrapper version.",
            details={"manifest": str(manifest_path)},
        )
    current_checksum = sha256_file(samplesheet_path)
    if manifest_checksum != current_checksum:
        _raise_resume_mismatch("samplesheet_checksum", current_checksum, manifest_checksum)
    return current_checksum



def _command_output(args: list[str]) -> str:
    try:
        proc = subprocess.run(args, capture_output=True, text=True, errors="replace", timeout=_SUBPROCESS_TIMEOUT)
    except (subprocess.TimeoutExpired, OSError):
        return ""
    if proc.returncode != 0:
        return ""
    return (proc.stdout or proc.stderr).strip()


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


def _check_java() -> dict[str, str]:
    java_path = shutil.which("java")
    if not java_path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_JAVA,
            message="Required executable `java` was not found.",
            fix="Install Java 17 or newer and ensure it is available on PATH.",
            details={"executable": "java"},
        )
    version_text = _command_output(["java", "-version"])
    version_tuple = _parse_version_tuple(version_text)
    if not version_tuple:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_JAVA,
            message="Java is installed but its version could not be determined.",
            fix="Install Java 17 or newer and ensure `java -version` works.",
            details={"java_path": java_path},
        )
    version_str = _detected_version_string(version_text)
    if version_tuple[0] < JAVA_MIN_VERSION:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.JAVA_VERSION_TOO_OLD,
            message="Java version is too old for nf-core/rnaseq.",
            fix="Install Java 17 or newer.",
            details={"detected_version": version_str},
        )
    return {"path": java_path, "version": version_str}


def _check_nextflow() -> dict[str, str]:
    nextflow_path = shutil.which("nextflow")
    if not nextflow_path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_NEXTFLOW,
            message="Required executable `nextflow` was not found.",
            fix="Install Nextflow and ensure it is available on PATH.",
            details={"executable": "nextflow"},
        )
    version_text = _command_output(["nextflow", "-version"])
    version_tuple = _parse_version_tuple(version_text)
    if not version_tuple:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_NEXTFLOW,
            message="Nextflow is installed but its version could not be determined.",
            fix="Install Nextflow 25.04.3 or newer and ensure `nextflow -version` works.",
            details={"nextflow_path": nextflow_path},
        )
    version_str = _detected_version_string(version_text)
    if _pad_version(version_tuple) < _pad_version(NEXTFLOW_MIN_VERSION):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.NEXTFLOW_VERSION_TOO_OLD,
            message="Nextflow version is too old for nf-core/rnaseq 3.26.0.",
            fix="Upgrade Nextflow to 25.04.3 or newer.",
            details={"detected_version": version_str},
        )
    return {"path": nextflow_path, "version": version_str}


def _check_nextflow_presence() -> dict[str, str | bool]:
    nextflow_path = shutil.which("nextflow")
    if not nextflow_path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_NEXTFLOW,
            message="Required executable `nextflow` was not found.",
            fix="Install Nextflow and ensure it is available on PATH.",
            details={"executable": "nextflow"},
        )
    return {"path": nextflow_path, "version": "", "version_checked": False}


_BACKEND_PROFILES = frozenset(
    {"docker", "podman", "conda", "mamba", "singularity", "apptainer", "shifter", "charliecloud"}
)


def _check_profile(profile: str) -> dict[str, str | bool]:
    """Validate the execution backend in a (possibly composite) profile string.

    nf-core/rnaseq accepts comma-separated profile lists such as
    ``docker,prokaryotic`` and ``singularity,arm64``, as well as undocumented
    profiles like ``test_full``, ``test_gpu``, and institutional HPC profiles.
    The wrapper validates only the first recognised backend component (docker,
    singularity, etc.) and passes all other parts through to Nextflow unchanged.
    """
    parts = [p.strip() for p in profile.split(",") if p.strip()]
    backend = next((p for p in parts if p in _BACKEND_PROFILES), None)
    if backend is None:
        # No recognised backend in the profile string — pass through; Nextflow validates.
        return {"profile": profile, "backend_path": "", "backend_ready": True}
    backend_info = _check_single_backend(backend)
    return {"profile": profile, "backend_path": backend_info["backend_path"], "backend_ready": backend_info["backend_ready"]}


def _check_single_backend(backend: str) -> dict[str, str | bool]:
    if backend == "docker":
        return _check_daemon_profile(backend, "docker", ErrorCode.MISSING_DOCKER, ErrorCode.DOCKER_NOT_RUNNING)
    if backend == "podman":
        return _check_daemon_profile(backend, "podman", ErrorCode.MISSING_PODMAN, ErrorCode.PODMAN_NOT_RUNNING)
    if backend == "conda":
        return _check_conda_like_profile(backend, ("conda", "mamba"))
    if backend == "mamba":
        return _check_conda_like_profile(backend, ("mamba", "conda"))
    if backend in {"singularity", "apptainer"}:
        return _check_singularity_compatible_profile(backend)
    if backend in {"shifter", "charliecloud"}:
        return _check_hpc_profile(backend)
    return {"profile": backend, "backend_path": "", "backend_ready": True}


def _check_daemon_profile(profile: str, binary: str, missing_code: str, not_running_code: str) -> dict[str, str | bool]:
    path = shutil.which(binary)
    if not path:
        raise SkillError(
            stage="preflight",
            error_code=missing_code,
            message=f"{profile} profile was selected but `{binary}` is not installed.",
            fix=f"Install {binary} or choose another supported profile.",
            details={"profile": profile},
        )
    try:
        info = subprocess.run([binary, "info"], capture_output=True, text=True, timeout=_SUBPROCESS_TIMEOUT)
        ok = info.returncode == 0
    except (subprocess.TimeoutExpired, OSError):
        ok = False
    if not ok:
        raise SkillError(
            stage="preflight",
            error_code=not_running_code,
            message=f"{binary} is installed but its service is not available.",
            fix=f"Start {binary} before running this skill.",
            details={"profile": profile},
        )
    return {"profile": profile, "backend_path": path, "backend_ready": True}


def _check_conda_like_profile(profile: str, candidates: tuple[str, ...]) -> dict[str, str | bool]:
    backend = next((path for name in candidates if (path := shutil.which(name))), None)
    if not backend:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_CONDA,
            message=f"{profile} profile was selected but no conda-compatible backend was found.",
            fix="Install conda or mamba, or choose another profile.",
            details={"profile": profile, "tried": list(candidates)},
        )
    return {"profile": profile, "backend_path": backend, "backend_ready": True}


def _check_singularity_compatible_profile(profile: str) -> dict[str, str | bool]:
    primary = "apptainer" if profile == "apptainer" else "singularity"
    fallback = "singularity" if profile == "apptainer" else "apptainer"
    backend = shutil.which(primary) or shutil.which(fallback)
    if not backend:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_SINGULARITY,
            message=f"{profile} profile was selected but no Singularity/Apptainer binary was found.",
            fix="Install Singularity or Apptainer, or choose a different profile.",
            details={"profile": profile, "tried": [primary, fallback]},
        )
    return {"profile": profile, "backend_path": backend, "backend_ready": True}


def _check_hpc_profile(profile: str) -> dict[str, str | bool]:
    binary_map = {"shifter": "shifter", "charliecloud": "ch-run"}
    binary = binary_map[profile]
    path = shutil.which(binary)
    if not path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_HPC_RUNTIME,
            message=f"{profile} profile was selected but `{binary}` was not found on PATH.",
            fix=f"Install {profile} and ensure `{binary}` is available on PATH, or choose a different profile.",
            details={"profile": profile, "binary": binary},
        )
    return {"profile": profile, "backend_path": path, "backend_ready": True}


def check_output_dir_available(output_dir: Path, *, resume: bool) -> None:
    _check_output_dir(output_dir, resume=resume)


def _check_output_dir(output_dir: Path, *, resume: bool) -> None:
    output_dir = output_dir.expanduser().resolve()
    if _is_relative_to(output_dir, PROJECT_ROOT.resolve()):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.OUTPUT_DIR_INSIDE_REPO,
            message="Output directory cannot be inside the ClawBio source tree.",
            fix="Choose an output directory outside the repository, for example under your analysis workspace.",
            details={"output_dir": str(output_dir), "project_root": str(PROJECT_ROOT.resolve())},
        )
    if output_dir.exists() and not output_dir.is_dir():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.OUTPUT_DIR_NOT_WRITABLE,
            message="Output path exists but is not a directory.",
            fix="Choose a directory path for --output, not an existing file.",
            details={"output_dir": str(output_dir)},
        )
    output_dir.mkdir(parents=True, exist_ok=True)
    if not os.access(output_dir, os.W_OK):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.OUTPUT_DIR_NOT_WRITABLE,
            message="Output directory is not writable.",
            fix="Choose a writable location for --output.",
            details={"output_dir": str(output_dir)},
        )
    materialized = []
    for entry in output_dir.iterdir():
        if entry.name in _IGNORED_ROOT_NAMES:
            continue
        materialized.append(entry)
    if materialized and not resume:
        # A prior FAILED run leaves result.json but no reproducibility/manifest.json,
        # so --resume is impossible (it needs a manifest) yet the dir is non-empty — a
        # dead end unless the user knows to delete it. Detect that and guide precisely
        # instead of the generic "use --resume" hint that cannot work here (audit F8).
        prior_run_incomplete = (
            (output_dir / "result.json").exists()
            and not (output_dir / "reproducibility" / "manifest.json").exists()
        )
        if prior_run_incomplete:
            fix = (
                "This directory holds an incomplete previous run (result.json present but no "
                "reproducibility/manifest.json), so --resume cannot work here. Delete the "
                "directory and re-run, or choose a new empty --output."
            )
        else:
            fix = "Choose a new empty output directory, or re-run with --resume when a compatible manifest exists."
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.OUTPUT_DIR_NOT_EMPTY,
            message="Output directory already contains files.",
            fix=fix,
            details={"output_dir": str(output_dir), "prior_run_incomplete": prior_run_incomplete},
        )


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _check_remote_inputs(args, samplesheet_summary: dict[str, object]) -> None:
    """Local-first by default: reject remote samplesheet inputs and reference paths
    (s3://, gs://, https://, ftp://, …) unless ``--allow-remote-inputs`` is set, in
    which case warn that genetic data will be fetched over the network. The
    object-store ``--work-dir`` is a separate, intentionally-remote setting and is
    not gated."""
    remote: list[str] = []
    for key in ("fastq_paths", "bam_paths"):
        for path in samplesheet_summary.get(key, []) or []:
            if "://" in str(path):
                remote.append(str(path))
    for ref_field in _EXPLICIT_REFERENCE_FIELDS:
        value = getattr(args, ref_field, None)
        if value and "://" in str(value):
            remote.append(str(value))
    if not remote:
        return
    unique = sorted(set(remote))
    if getattr(args, "allow_remote_inputs", False):
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


def _check_demo_network(args) -> None:
    """--demo runs nf-core's upstream ``-profile test``, whose inputs (test FASTQs and
    reference FASTA/GTF) are remote GitHub URLs by design. Under ``NXF_OFFLINE`` the
    nf-schema plugin aborts during parameter validation with a confusing
    ``does not exist`` before any work starts, so fail early with an actionable message
    instead. Downloading nf-core's PUBLIC test data does not conflict with the
    local-first guarantee, which governs USER genetic data (never uploaded)."""
    if not getattr(args, "demo", False):
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


def run_preflight(args, *, pipeline_source: dict[str, object], samplesheet_summary: dict[str, object]) -> dict[str, object]:
    warnings: list[str] = []
    _check_demo_network(args)
    _check_remote_inputs(args, samplesheet_summary)
    aligner_effective = _check_supported_options(args)
    _check_igenomes_genome(args, pipeline_source=pipeline_source, warnings=warnings)
    _check_email(getattr(args, "email", None))
    _check_email(getattr(args, "email_on_fail", None), field="email_on_fail")
    _check_threshold_bounds(args)
    java_info = _check_java()
    nextflow_info = _check_nextflow_presence() if getattr(args, "check", False) else _check_nextflow()
    if nextflow_info.get("version_checked") is False:
        # --check verifies Nextflow is present but intentionally defers the
        # >=25.04.3 version gate to the real run (audit F-5). Surface that so a
        # check is not silently assumed to confirm a compatible Nextflow.
        warnings.append(
            "Nextflow version was not verified in --check mode; the >=25.04.3 "
            "requirement is enforced on the actual run."
        )
    profile_info = _check_profile(args.profile)
    output_dir = Path(args.output).expanduser().resolve()
    check_output_dir_available(output_dir, resume=getattr(args, "resume", False))
    _check_contaminant_paths(args)
    _check_ribo_paths(args)
    _check_nextflow_configs(args, warnings)
    _check_acceleration_compatibility(args)
    _check_contaminant_screening(args, warnings)
    _check_mutual_exclusions(args)
    _check_samplesheet_summary(args, samplesheet_summary=samplesheet_summary)
    bam_paths = samplesheet_summary.get("bam_paths", [])
    refs = {} if (getattr(args, "demo", False) or getattr(args, "_noinput", False)) else _check_references(args, bam_reprocessing=bool(bam_paths), warnings=warnings)
    _check_skip_alignment_requirements(args, bam_paths=bam_paths)
    _collect_bam_reprocessing_warning(bam_paths, aligner_effective, warnings)
    _collect_index_aligner_warnings(args, warnings)
    gencode_autodetected = _collect_gencode_autodetect_warning(args, warnings)
    _check_resume_compatibility(args, output_dir=output_dir, pipeline_source=pipeline_source)
    _collect_umi_warnings(args, warnings)
    _collect_demo_warnings(args, warnings)
    _collect_platform_warnings(args, output_dir, warnings)
    handoff_available = _collect_alignment_warnings(args, warnings)

    return {
        "ok": True,
        "java": java_info,
        "nextflow": nextflow_info,
        "profile": profile_info,
        "pipeline_source": pipeline_source,
        "references": refs,
        "aligner_effective": aligner_effective,
        "handoff_available": handoff_available,
        "gencode_autodetected": gencode_autodetected,
        "warnings": warnings,
        "samplesheet": {
            "sample_count": samplesheet_summary.get("sample_count", 0),
            "sample_names": samplesheet_summary.get("sample_names", []),
            "unknown_columns": samplesheet_summary.get("unknown_columns", []),
            "strandedness_counts": samplesheet_summary.get("strandedness_counts", {}),
            "bam_count": len(samplesheet_summary.get("bam_paths", [])),
        },
    }


def _check_supported_options(args) -> str:
    aligner = getattr(args, "aligner", "star_salmon")
    if aligner not in SUPPORTED_ALIGNERS:
        _raise_invalid_enum(ErrorCode.INVALID_ALIGNER, "aligner", aligner, SUPPORTED_ALIGNERS)
    pseudo_aligner = getattr(args, "pseudo_aligner", None)
    if pseudo_aligner and pseudo_aligner not in SUPPORTED_PSEUDO_ALIGNERS:
        _raise_invalid_enum(ErrorCode.INVALID_PSEUDO_ALIGNER, "pseudo_aligner", pseudo_aligner, SUPPORTED_PSEUDO_ALIGNERS)
    trimmer = getattr(args, "trimmer", DEFAULT_TRIMMER)
    if trimmer not in SUPPORTED_TRIMMERS:
        _raise_invalid_enum(ErrorCode.INVALID_TRIMMER, "trimmer", trimmer, SUPPORTED_TRIMMERS)
    ribo_tool = getattr(args, "ribo_removal_tool", None)
    if ribo_tool and ribo_tool not in SUPPORTED_RIBO_TOOLS:
        _raise_invalid_enum(ErrorCode.INVALID_RIBO_TOOL, "ribo_removal_tool", ribo_tool, SUPPORTED_RIBO_TOOLS)
    umi_tool = getattr(args, "umi_dedup_tool", None)
    if umi_tool and umi_tool not in SUPPORTED_UMI_TOOLS:
        _raise_invalid_enum(ErrorCode.INVALID_UMI_TOOL, "umi_dedup_tool", umi_tool, SUPPORTED_UMI_TOOLS)
    umi_extract_method = getattr(args, "umitools_extract_method", None)
    if umi_extract_method and umi_extract_method not in SUPPORTED_UMI_EXTRACT_METHODS:
        _raise_invalid_enum(ErrorCode.INVALID_PRESET_CONFIGURATION, "umitools_extract_method", umi_extract_method, SUPPORTED_UMI_EXTRACT_METHODS)
    umi_grouping_method = getattr(args, "umitools_grouping_method", None)
    if umi_grouping_method and umi_grouping_method not in SUPPORTED_UMI_GROUPING_METHODS:
        _raise_invalid_enum(ErrorCode.INVALID_PRESET_CONFIGURATION, "umitools_grouping_method", umi_grouping_method, SUPPORTED_UMI_GROUPING_METHODS)
    salmon_libtype = getattr(args, "salmon_quant_libtype", None)
    if salmon_libtype and salmon_libtype not in SUPPORTED_SALMON_LIB_TYPES:
        _raise_invalid_enum(ErrorCode.INVALID_PRESET_CONFIGURATION, "salmon_quant_libtype", salmon_libtype, SUPPORTED_SALMON_LIB_TYPES)
    publish_dir_mode = getattr(args, "publish_dir_mode", None)
    if publish_dir_mode and publish_dir_mode not in SUPPORTED_PUBLISH_DIR_MODES:
        _raise_invalid_enum(ErrorCode.INVALID_PRESET_CONFIGURATION, "publish_dir_mode", publish_dir_mode, SUPPORTED_PUBLISH_DIR_MODES)
    _check_rseqc_modules(getattr(args, "rseqc_modules", None))
    # Profile validation is handled by _check_profile (called later in run_preflight).
    # Composite profiles (e.g. "docker,prokaryotic") and undocumented profiles
    # (test_full, institutional HPC) must not be rejected here.
    return "star_salmon" if getattr(args, "demo", False) else aligner


def _check_rseqc_modules(rseqc_modules: str | None) -> None:
    """Reject unknown RSeQC module names locally instead of failing late in Nextflow.

    ``--rseqc-modules`` is a comma-separated list; nf-core/rnaseq 3.26.0 only accepts
    the eight names in SUPPORTED_RSEQC_MODULES. Catching a typo here gives a clear
    error before any container is pulled (audit F-04).
    """
    if not rseqc_modules:
        return
    requested = [m.strip() for m in str(rseqc_modules).split(",") if m.strip()]
    unknown = [m for m in requested if m not in SUPPORTED_RSEQC_MODULES]
    if unknown:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message=f"Unsupported rseqc-modules value(s): {', '.join(unknown)}.",
            fix=f"Choose from: {', '.join(sorted(SUPPORTED_RSEQC_MODULES))}.",
            details={"field": "rseqc_modules", "unknown": unknown, "requested": requested},
        )


def _raise_invalid_enum(code: str, field: str, value: str, supported: set[str]) -> None:
    raise SkillError(
        stage="preflight",
        error_code=code,
        message=f"Unsupported {field.replace('_', '-')} value.",
        fix=f"Choose one of: {', '.join(sorted(supported))}.",
        details={"field": field, "value": value},
    )


def _check_igenomes_genome(args, *, pipeline_source: dict[str, object], warnings: list[str]) -> None:
    if getattr(args, "demo", False) or getattr(args, "_noinput", False):
        return
    genome = getattr(args, "genome", None)
    if not genome:
        return
    if not _GENOME_ID_RE.match(genome):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--genome identifier contains invalid characters (nf-core schema: ^[a-zA-Z0-9_\\-\\.]+$).",
            fix="Use only letters, digits, hyphens, underscores, and dots in the genome identifier (e.g. 'GRCh38').",
            details={"field": "genome", "value": genome},
        )
    if not _uses_default_audited_pipeline_source(pipeline_source):
        warnings.append(
            "iGenomes genome validation was skipped because the pipeline source differs from audited nf-core/rnaseq 3.26.0."
        )
        return
    if genome in SUPPORTED_IGENOMES_NAMES:
        return
    warnings.append(
        f"--genome {genome!r} is not a recognised iGenomes name for nf-core/rnaseq 3.26.0. "
        "If you are using a user-defined genome catalogue (via --nextflow-config my_genomes.config), this is expected. "
        "If you intended an iGenomes entry, check the exact spelling and case "
        f"(e.g. 'GRCh38', 'GRCm38'). "
        "Supported iGenomes keys: " + ", ".join(sorted(SUPPORTED_IGENOMES_NAMES)[:8]) + ", ..."
    )


def _uses_default_audited_pipeline_source(pipeline_source: dict[str, object]) -> bool:
    return (
        pipeline_source.get("source_kind") == "remote_repo"
        and pipeline_source.get("source_ref") == DEFAULT_REMOTE_PIPELINE
        and pipeline_source.get("resolved_version") == DEFAULT_PIPELINE_VERSION
    )


def _check_email(email: str | None, *, field: str = "email") -> None:
    if email and not _EMAIL_RE.match(email):
        flag = f"--{field.replace('_', '-')}"
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message=f"{field} address does not match the nf-core/rnaseq schema.",
            fix=f"Use a simple email address such as user@example.org, or omit {flag}.",
            details={"field": field, "value": email},
        )


def _check_threshold_bounds(args) -> None:
    st = getattr(args, "stranded_threshold", None)
    if st is not None and not (0.5 <= st <= 1.0):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="stranded_threshold is out of range (nf-core schema: minimum 0.5, maximum 1.0).",
            fix="Set --stranded-threshold to a value between 0.5 and 1.0 inclusive.",
            details={"field": "stranded_threshold", "value": st, "minimum": 0.5, "maximum": 1.0},
        )
    ut = getattr(args, "unstranded_threshold", None)
    if ut is not None and not (0.0 <= ut <= 1.0):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="unstranded_threshold is out of range (nf-core schema: minimum 0, maximum 1).",
            fix="Set --unstranded-threshold to a value between 0.0 and 1.0 inclusive.",
            details={"field": "unstranded_threshold", "value": ut, "minimum": 0.0, "maximum": 1.0},
        )
    mtr = getattr(args, "min_trimmed_reads", None)
    if mtr is not None and mtr < 0:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="min_trimmed_reads must be ≥ 0.",
            fix="Set --min-trimmed-reads to 0 or a positive integer.",
            details={"field": "min_trimmed_reads", "value": mtr, "minimum": 0},
        )
    fraglen = getattr(args, "kallisto_quant_fraglen", None)
    if fraglen is not None and fraglen < 1:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="kallisto_quant_fraglen must be ≥ 1 (pipeline default: 200).",
            fix="Set --kallisto-quant-fraglen to 1 or a positive integer (pipeline default: 200).",
            details={"field": "kallisto_quant_fraglen", "value": fraglen, "minimum": 1},
        )
    fraglen_sd = getattr(args, "kallisto_quant_fraglen_sd", None)
    if fraglen_sd is not None and fraglen_sd < 0:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="kallisto_quant_fraglen_sd must be ≥ 0 (pipeline default: 200).",
            fix="Set --kallisto-quant-fraglen-sd to 0 or a positive integer (pipeline default: 200).",
            details={"field": "kallisto_quant_fraglen_sd", "value": fraglen_sd, "minimum": 0},
        )
    mmr = getattr(args, "min_mapped_reads", None)
    # min_mapped_reads is a percentage of uniquely-mapped reads (pipeline default 5).
    # The nf-core schema declares no explicit maximum, but a percentage > 100 is
    # meaningless (it would discard every sample), so the wrapper intentionally bounds
    # it to 0..100 — a guard stricter than the schema, never blocking a usable value
    # (audit F4).
    if mmr is not None and not (0 <= mmr <= 100):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="min_mapped_reads must be between 0 and 100 (it is a percentage, pipeline default: 5).",
            fix="Set --min-mapped-reads to a value between 0 and 100 inclusive.",
            details={"field": "min_mapped_reads", "value": mmr, "minimum": 0, "maximum": 100},
        )
    kmer = getattr(args, "pseudo_aligner_kmer_size", None)
    # Salmon and Kallisto both encode the index k-mer in a 64-bit machine word, so
    # the k-mer must be ODD and ≤ 31 (their shared hard cap; pipeline default 31).
    # The nf-core schema declares no bound, but an even or out-of-range value is a
    # guaranteed indexing crash — catch it early like every other numeric field
    # rather than let the pipeline abort late (audit F9).
    if kmer is not None and (kmer < 1 or kmer > 31 or kmer % 2 == 0):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message=(
                "pseudo_aligner_kmer_size must be an odd integer between 1 and 31 "
                "(Salmon/Kallisto index constraint; pipeline default: 31)."
            ),
            fix="Set --pseudo-aligner-kmer-size to an odd value in 1..31 (e.g. 31 for reads ≥75 bp, lower for short reads).",
            details={"field": "pseudo_aligner_kmer_size", "value": kmer, "minimum": 1, "maximum": 31},
        )
    timeout_hours = getattr(args, "timeout_hours", None)
    if timeout_hours is not None and timeout_hours < 0:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="timeout_hours must be >= 0 (0 disables the wall-clock cap).",
            fix=(
                "Set --timeout-hours to a positive number of hours (default: 12), "
                "or 0 to disable the cap for long HPC/cloud runs whose walltime is "
                "enforced by the scheduler."
            ),
            details={"field": "timeout_hours", "value": timeout_hours, "minimum": 0},
        )
    hbm = getattr(args, "hisat2_build_memory", None)
    if hbm is not None and not _HISAT2_MEM_RE.match(str(hbm)):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="hisat2_build_memory format is invalid (nf-core schema: ^\\d+(\\.\\d+)?(\\.)?\\s*(K|M|G|T)?B$, e.g. '200.GB').",
            fix="Use a value like '200.GB', '8GB', or '1.5TB'.",
            details={"field": "hisat2_build_memory", "value": hbm},
        )
    umi_sep = getattr(args, "umitools_umi_separator", None)
    if umi_sep is not None:
        # Validate against the raw value: nf-core schema requires pattern ^\S+$ (no whitespace
        # anywhere) AND maxLength: 1.  Stripping first would silently accept " :" (space+colon),
        # which fails the schema pattern check.  Always validate the raw string.
        import re as _re
        if not _re.match(r"^\S+$", umi_sep):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message="umitools_umi_separator must be a non-whitespace string (nf-core schema pattern: ^\\S+$).",
                fix="Set --umitools-umi-separator to a non-whitespace character, e.g. ':' or '_'.",
                details={"field": "umitools_umi_separator", "value": repr(umi_sep)},
            )
        if len(umi_sep) > 1:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message=f"umitools_umi_separator must be a single character (nf-core schema maxLength: 1); got {len(umi_sep)} characters.",
                fix="Set --umitools-umi-separator to a single character, e.g. ':' or '_'.",
                details={"field": "umitools_umi_separator", "value": repr(umi_sep)},
            )


def _check_references(args, *, bam_reprocessing: bool = False, warnings: list[str] | None = None) -> dict[str, str]:
    if warnings is None:
        warnings = []
    values = _collect_reference_values(args)
    # --gtf + --gff together: nf-core/rnaseq accepts both and uses --gtf, ignoring
    # --gff (docs: "the latter [GTF] will be used if both are provided"). Match that
    # exactly — keep --gtf, drop --gff (so it is never written to params.yaml or the
    # replay command), and warn. Rejecting here would block a valid upstream config
    # (audit F2).
    if values.get("gtf") and values.get("gff"):
        warnings.append(
            "Both --gtf and --gff were provided; nf-core/rnaseq uses --gtf and ignores "
            "--gff. Dropping --gff to match upstream behaviour."
        )
        values.pop("gff", None)
        args.gff = None
    if values.get("genome"):
        conflicting = [field for field in _GENOME_CONFLICTING_REFERENCE_FIELDS if values.get(field)]
        if conflicting:
            # --genome provides a complete iGenomes reference. A second genome
            # sequence source (--fasta) or a genome-level index (STAR/RSEM/HISAT2/
            # Bowtie2) is genuinely ambiguous, so it is rejected. Additive
            # annotation/transcriptome overrides (--gtf/--gff/--gene-bed/
            # --transcript-fasta/--additional-fasta/--splicesites/--salmon-index/
            # --kallisto-index) ARE allowed alongside --genome — nf-core/rnaseq
            # supports them and they are common (e.g. ERCC spike-ins, custom GTF).
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.CONFLICTING_REFERENCES,
                message="--genome cannot be combined with an explicit genome FASTA or a genome-level index.",
                fix=(
                    "Use --genome <iGenomes_name> (optionally with additive overrides such as "
                    "--gtf, --gff, --additional-fasta, --transcript-fasta, --gene-bed, --splicesites, "
                    "--salmon-index, or --kallisto-index), OR a fully explicit --fasta/--gtf reference. "
                    "Do not provide both a named genome and your own genome FASTA/index."
                ),
                details={"genome": values["genome"], "explicit_set": conflicting},
            )
    # Enforce standard-group uniqueness first (fasta+gtf vs fasta+gff is a real
    # conflict regardless of BAM mode; genome vs explicit is already handled above).
    satisfied = [group for group in _REFERENCE_GROUPS if all(values.get(field) for field in group)]
    if len(satisfied) > 1:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.CONFLICTING_REFERENCES,
            message="Multiple reference groups were provided.",
            fix="Provide only one of: --genome, --fasta/--gtf, or --fasta/--gff.",
            details={"satisfied_groups": [list(group) for group in satisfied]},
        )
    if not satisfied:
        # BAM reprocessing with --skip-alignment: the official nf-core/rnaseq 3.26.0
        # docs show this mode with NO explicit reference (the pipeline uses info from
        # the BAM header / embedded transcriptome). Requiring refs here would block a
        # documented upstream use-case, so we allow it through. Supplying optional
        # refs (--gtf, --salmon-index, --transcript-fasta) is still accepted because
        # they do not conflict with each other in this mode.
        if bam_reprocessing and getattr(args, "skip_alignment", False):
            pass
        elif _is_transcriptome_only_pseudo_route(args, values):
            # Transcriptome-only pseudo-alignment: Salmon/Kallisto can quantify
            # against a transcript FASTA or a prebuilt salmon/kallisto index without
            # a genome FASTA. A GTF/GFF is still required for tximport's tx2gene
            # gene-level summarisation, so this is accepted only with annotation
            # present (audit F-03).
            pass
        elif not getattr(args, "skip_alignment", False) and _alignment_reference_complete(args, values):
            # Fully prebuilt references: a genome index matching the aligner plus a
            # transcript source and annotation needs no --fasta (audit F-12).
            pass
        else:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.MISSING_REFERENCE,
                message="A reference group is required.",
                fix=(
                    "Provide --genome, --fasta plus --gtf, or --fasta plus --gff. "
                    "For pseudo-alignment only, --transcript-fasta or a --salmon-index/"
                    "--kallisto-index plus --gtf/--gff is also accepted."
                ),
                details={"accepted_combinations": [list(group) for group in _REFERENCE_GROUPS]},
            )
    _check_reference_paths_exist(values)
    return values


# Prebuilt genome index that lets each aligner run without --fasta.
_ALIGNER_GENOME_INDEX_FIELDS = {
    "star_salmon": ("star_index",),
    "star_rsem": ("star_index", "rsem_index"),
    "hisat2": ("hisat2_index",),
    "bowtie2_salmon": ("bowtie2_index",),
}


def _alignment_reference_complete(args, values: dict[str, str]) -> bool:
    """True when prebuilt references fully specify an alignment route without --fasta.

    nf-core/rnaseq does not need a genome FASTA when every artifact it would build
    from it is already supplied (audit F-12). A route is complete when:

    * an annotation (``--gtf``/``--gff``) is present — always required for quantification;
    * a genome source for the aligner is present — ``--fasta`` or the aligner's prebuilt
      index (STAR/HISAT2/Bowtie2; RSEM also accepts ``--rsem-index``); and
    * a transcript source is present for quantifying routes — ``--fasta`` (to derive it),
      ``--transcript-fasta`` or ``--salmon-index`` for the Salmon routes, and
      ``--fasta`` or ``--rsem-index`` for RSEM. HISAT2-only alignment needs none.
    """
    aligner = getattr(args, "aligner", "star_salmon")
    has_annotation = bool(values.get("gtf") or values.get("gff"))
    if not has_annotation:
        return False
    has_fasta = bool(values.get("fasta"))
    index_fields = _ALIGNER_GENOME_INDEX_FIELDS.get(aligner, ())
    has_genome_source = has_fasta or any(values.get(f) for f in index_fields)
    if not has_genome_source:
        return False
    if aligner in {"star_salmon", "bowtie2_salmon"}:
        return has_fasta or bool(values.get("transcript_fasta") or values.get("salmon_index"))
    if aligner == "star_rsem":
        return has_fasta or bool(values.get("rsem_index"))
    return True  # hisat2: alignment only, no transcript source required


def _is_transcriptome_only_pseudo_route(args, values: dict[str, str]) -> bool:
    """True when a Salmon/Kallisto pseudo route can run off a transcriptome alone.

    Requires (a) the genome aligner to be skipped (``--skip-alignment``) so no STAR/
    HISAT2/Bowtie2 step demands a genome FASTA — a pseudo-aligner running *alongside*
    a genome aligner still needs the genome and is intentionally NOT covered here;
    (b) a transcript source (``--transcript-fasta`` or a prebuilt ``--salmon-index``/
    ``--kallisto-index``); and (c) an annotation (``--gtf``/``--gff``) for tximport's
    tx2gene gene-level summarisation. Missing any of these falls through to
    MISSING_REFERENCE (audit F-03).
    """
    if not getattr(args, "skip_alignment", False):
        return False
    has_transcript_source = any(values.get(f) for f in ("transcript_fasta", "salmon_index", "kallisto_index"))
    has_annotation = bool(values.get("gtf") or values.get("gff"))
    return has_transcript_source and has_annotation


def _collect_reference_values(args) -> dict[str, str]:
    values: dict[str, str] = {}
    genome = getattr(args, "genome", None) or ""
    if genome:
        values["genome"] = genome
    for field in _EXPLICIT_REFERENCE_FIELDS:
        val = getattr(args, field, None) or ""
        if val:
            values[field] = val
    return values


def _check_reference_paths_exist(values: dict[str, str]) -> None:
    _EXTENSION_CHECKS: dict[str, tuple[re.Pattern[str], str, str]] = {
        "fasta": (
            _FASTA_EXT_RE,
            "does not have a recognised FASTA extension "
            "(nf-core schema: .fa, .fasta, .fna, .fa.gz, .fasta.gz, .fna.gz).",
            "Rename the file to use a standard FASTA extension, e.g. genome.fa or genome.fasta.gz.",
        ),
        "transcript_fasta": (
            _FASTA_EXT_RE,
            "does not have a recognised FASTA extension "
            "(nf-core schema: .fa, .fasta, .fna, .fa.gz, .fasta.gz, .fna.gz).",
            "Rename the file to use a standard FASTA extension, e.g. transcriptome.fa or transcriptome.fasta.gz.",
        ),
        "additional_fasta": (
            _FASTA_EXT_RE,
            "does not have a recognised FASTA extension "
            "(nf-core schema: .fa, .fasta, .fna, .fa.gz, .fasta.gz, .fna.gz).",
            "Rename the file to use a standard FASTA extension, e.g. spike_in.fa or spike_in.fa.gz.",
        ),
        "gtf": (
            _GTF_EXT_RE,
            "does not have a recognised GTF extension (nf-core schema: .gtf or .gtf.gz).",
            "Rename the file to use a standard GTF extension, e.g. genes.gtf or genes.gtf.gz.",
        ),
        "gff": (
            _GFF_EXT_RE,
            "does not have a recognised GFF extension (nf-core schema: .gff or .gff3, optionally .gz).",
            "Rename the file to use a standard GFF extension, e.g. genes.gff3 or genes.gff.gz.",
        ),
        "gene_bed": (
            _BED_EXT_RE,
            "does not have a recognised BED extension (nf-core schema: .bed or .bed.gz).",
            "Rename the file to use a standard BED extension, e.g. genes.bed or genes.bed.gz.",
        ),
    }
    for field, value in values.items():
        if field == "genome" or not value:
            continue
        if "://" in value:
            # URI scheme (https://, s3://, gs://, ftp://) — existence check not applicable;
            # Nextflow resolves remote refs directly at runtime.
            continue
        if field in _EXTENSION_CHECKS:
            # The wrapper resolves reference paths to absolute POSIX before writing
            # params.yaml, so a path that resolves into a directory containing
            # whitespace would fail the nf-core ^\S+ schema pattern at runtime
            # (audit F-1). Catch it here — on the *resolved* path — with a precise,
            # dedicated error, mirroring the samplesheet input guard. This check
            # precedes the extension check because the ^\S+ extension regex also
            # fails on whitespace and would otherwise report a misleading message.
            resolved = Path(value).expanduser().resolve()
            if _WHITESPACE_RE.search(resolved.as_posix()):
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.REFERENCE_PATH_HAS_WHITESPACE,
                    message=(
                        f"--{field.replace('_', '-')} resolves to a path containing whitespace, "
                        "which the nf-core/rnaseq schema rejects (pattern ^\\S+)."
                    ),
                    fix=(
                        "Move the reference into a directory whose full path has no spaces, "
                        "or symlink it to a whitespace-free location, then point the flag there."
                    ),
                    details={"field": field, "value": value, "resolved": resolved.as_posix()},
                )
            pat, msg, fix = _EXTENSION_CHECKS[field]
            if not pat.match(value):
                raise SkillError(
                    stage="preflight",
                    error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                    message=f"--{field.replace('_', '-')} {msg}",
                    fix=fix,
                    details={"field": field, "value": value},
                )
        if not Path(value).expanduser().exists():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.REFERENCE_PATH_NOT_FOUND,
                message="A reference path was not found.",
                fix="Correct the missing reference path and try again.",
                details={"field": field, "path": value},
            )


def _collect_gencode_autodetect_warning(args, warnings: list[str]) -> bool:
    if getattr(args, "demo", False) or getattr(args, "gencode", False):
        return False
    gtf = getattr(args, "gtf", None)
    if not gtf or not _gtf_has_gencode_markers(Path(gtf).expanduser()):
        return False
    warnings.append("GENCODE-style GTF markers were autodetected; setting gencode: true for nf-core/rnaseq.")
    return True


def _gtf_has_gencode_markers(path: Path) -> bool:
    feature_records_seen = 0
    try:
        opener = gzip.open if path.suffix.lower() == ".gz" else Path.open
        with opener(path, mode="rt", encoding="utf-8") as handle:
            for line in handle:
                text = line.strip()
                if not text or text.startswith("#"):
                    continue
                feature_records_seen += 1
                if _GENCODE_GTF_MARKER_RE.search(text):
                    return True
                if feature_records_seen >= _GENCODE_GTF_MAX_FEATURE_RECORDS:
                    break
    except (OSError, UnicodeError):
        return False
    return False


def _check_contaminant_paths(args) -> None:
    for field in _CONTAMINANT_PATH_FIELDS:
        value = getattr(args, field, None)
        if not value:
            continue
        for item in _split_contaminant_paths(field, value):
            _check_existing_contaminant_path(field, item)


def _split_contaminant_paths(field: str, value: str) -> list[str]:
    if field not in _COMMA_SEPARATED_CONTAMINANT_FIELDS:
        return [value]
    return [item.strip() for item in value.split(",") if item.strip()]


def _check_existing_contaminant_path(field: str, value: str) -> None:
    if "://" in value:
        # Remote databases are staged by Nextflow; local existence checks only
        # apply to filesystem paths.
        return
    path = Path(value).expanduser()
    if not path.exists():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.REFERENCE_PATH_NOT_FOUND,
            message="A contaminant screening reference path was not found.",
            fix="Correct the missing contaminant screening path and try again.",
            details={"field": field, "path": value},
        )
    if field == "bbsplit_fasta_list" and not path.is_file():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--bbsplit-fasta-list must be a TSV file path, not a directory.",
            fix=(
                "Provide the path to a tab-separated file listing FASTA genomes for BBSplit "
                "(one genome per line: <name>\t<path>). See the nf-core/rnaseq docs for format."
            ),
            details={"field": field, "path": value, "expected": "file", "got": "directory"},
        )


def _check_ribo_paths(args) -> None:
    manifest = getattr(args, "ribo_database_manifest", None)
    if not manifest:
        return
    _check_existing_ribo_manifest(manifest)


def _check_existing_ribo_manifest(value: str) -> None:
    if "://" in value:
        # Remote URIs are resolved by Nextflow at runtime; local existence checks
        # are only meaningful for filesystem paths.
        return
    path = Path(value).expanduser()
    if not path.exists():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.REFERENCE_PATH_NOT_FOUND,
            message="The rRNA database manifest path was not found.",
            fix="Correct --ribo-database-manifest or remove it.",
            details={"field": "ribo_database_manifest", "path": value},
        )
    if not path.is_file():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--ribo-database-manifest must point to a manifest file, not a directory.",
            fix="Provide the manifest file path expected by the selected rRNA removal tool.",
            details={"field": "ribo_database_manifest", "path": value, "expected": "file", "got": "directory"},
        )


def _check_nextflow_configs(args, warnings: list[str]) -> None:
    for raw_path in getattr(args, "nextflow_config", None) or []:
        path = Path(raw_path).expanduser()
        if not path.exists():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.REFERENCE_PATH_NOT_FOUND,
                message="A Nextflow config path was not found.",
                fix="Correct --nextflow-config or remove it.",
                details={"field": "nextflow_config", "path": raw_path},
            )
        if not path.is_file():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message="--nextflow-config must point to a config file, not a directory.",
                fix="Provide a .config file path.",
                details={"field": "nextflow_config", "path": raw_path, "expected": "file", "got": "directory"},
            )
        _reject_nextflow_params_overrides(path, warnings=warnings)


def _reject_nextflow_params_overrides(
    path: Path, *, warnings: list[str], _seen: set[Path] | None = None, _depth: int = 0
) -> None:
    """Reject params.* overrides in a --nextflow-config file and its includeConfig chain.

    The audited-parameter-surface guarantee is only meaningful if it cannot be
    bypassed. Besides the direct params forms (block/property/assignment/subscript/
    merge — see _NEXTFLOW_PARAMS_OVERRIDE_RE), a params override could be hidden in
    an `includeConfig` target, so every locally-resolvable include is audited
    recursively (with cycle/depth guards). Includes the wrapper cannot read locally
    (remote URIs, interpolated `${...}` paths, or missing files) are surfaced as
    warnings rather than silently trusted — the wrapper is honest about what it
    could not audit (audit M2).
    """
    if _seen is None:
        _seen = set()
    resolved = path.resolve()
    if resolved in _seen:
        return
    _seen.add(resolved)
    try:
        lines = path.read_text(encoding="utf-8").splitlines()
    except (OSError, UnicodeError) as exc:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--nextflow-config could not be read as UTF-8 text.",
            fix="Provide a readable Nextflow config file.",
            details={"field": "nextflow_config", "path": str(path), "error": str(exc)},
        ) from exc
    for line_number, line in enumerate(lines, start=1):
        stripped = line.strip()
        if not stripped or stripped.startswith(("#", "//", "/*", "*")):
            continue
        if _NEXTFLOW_PARAMS_OVERRIDE_RE.match(line) and not _NEXTFLOW_ALLOWED_PARAMS_RE.match(line):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message="--nextflow-config must not override nf-core params.* values.",
                fix=(
                    "Pass pipeline parameters through audited wrapper flags/params.yaml. "
                    "Keep --nextflow-config for runtime settings such as process, executor, profiles, labels, "
                    "and the documented params.genomes custom genome catalogue."
                ),
                details={
                    "field": "nextflow_config",
                    "path": str(path),
                    "line": line_number,
                    "text": stripped[:160],
                },
            )
        include_match = _NEXTFLOW_INCLUDE_RE.match(line)
        if include_match:
            _audit_included_config(
                include_match.group(1),
                parent=path.parent,
                source=path,
                line_number=line_number,
                warnings=warnings,
                _seen=_seen,
                _depth=_depth,
            )


def _audit_included_config(
    target: str,
    *,
    parent: Path,
    source: Path,
    line_number: int,
    warnings: list[str],
    _seen: set[Path],
    _depth: int,
) -> None:
    if "://" in target or "${" in target:
        warnings.append(
            f"--nextflow-config '{source}' line {line_number}: includeConfig '{target}' is a "
            "remote or dynamically-interpolated path the wrapper cannot audit for params overrides. "
            "Ensure it does not set nf-core params.* values."
        )
        return
    if _depth + 1 > _MAX_CONFIG_INCLUDE_DEPTH:
        warnings.append(
            f"--nextflow-config '{source}' line {line_number}: includeConfig nesting exceeded "
            f"{_MAX_CONFIG_INCLUDE_DEPTH} levels and was not fully audited for params overrides."
        )
        return
    included = Path(target).expanduser()
    if not included.is_absolute():
        included = parent / included
    if not included.exists() or not included.is_file():
        warnings.append(
            f"--nextflow-config '{source}' line {line_number}: includeConfig '{target}' could not be "
            "resolved to a local file, so it was not audited for params overrides. "
            "Ensure it does not set nf-core params.* values."
        )
        return
    _reject_nextflow_params_overrides(included, warnings=warnings, _seen=_seen, _depth=_depth + 1)


def _check_acceleration_compatibility(args) -> None:
    aligner = getattr(args, "aligner", "star_salmon")
    if getattr(args, "use_parabricks_star", False) and aligner != "star_salmon":
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--use-parabricks-star requires --aligner star_salmon.",
            fix="Switch to --aligner star_salmon, or remove --use-parabricks-star.",
            details={"field": "use_parabricks_star", "aligner": aligner},
        )
    if getattr(args, "use_sentieon_star", False) and aligner not in {"star_salmon", "star_rsem"}:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--use-sentieon-star requires a STAR-based aligner.",
            fix="Switch to --aligner star_salmon or --aligner star_rsem, or remove --use-sentieon-star.",
            details={"field": "use_sentieon_star", "aligner": aligner},
        )
    if getattr(args, "use_gpu_ribodetector", False):
        remove_ribo = bool(getattr(args, "remove_ribo_rna", False))
        ribo_tool = getattr(args, "ribo_removal_tool", None)
        if not remove_ribo or ribo_tool != "ribodetector":
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message="--use-gpu-ribodetector requires --remove-ribo-rna with --ribo-removal-tool ribodetector.",
                fix=(
                    "Add --remove-ribo-rna --ribo-removal-tool ribodetector, "
                    "or remove --use-gpu-ribodetector."
                ),
                details={
                    "field": "use_gpu_ribodetector",
                    "remove_ribo_rna": remove_ribo,
                    "ribo_removal_tool": ribo_tool,
                },
            )


def _check_contaminant_screening(args, warnings: list[str]) -> None:
    """Ensure the selected contaminant-screening tool has its required database.

    nf-core/rnaseq 3.26.0 needs a Kraken2 database for ``kraken2``/``kraken2_bracken``
    and a Sylph database for ``sylph``. Without it the run fails inside Nextflow after
    container pulls, so catch it at preflight (audit F-10). ``--bracken-precision`` only
    applies to ``kraken2_bracken``; with any other tool it is inert and merely warned.
    """
    screening = getattr(args, "contaminant_screening", None)
    if screening in {"kraken2", "kraken2_bracken"} and not getattr(args, "kraken_db", None):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message=f"--contaminant-screening {screening} requires a Kraken2 database.",
            fix="Provide --kraken-db <path>, or remove --contaminant-screening.",
            details={"field": "kraken_db", "contaminant_screening": screening},
        )
    if screening == "sylph" and not getattr(args, "sylph_db", None):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--contaminant-screening sylph requires a Sylph database.",
            fix="Provide --sylph-db <path[,path...]>, or remove --contaminant-screening.",
            details={"field": "sylph_db", "contaminant_screening": screening},
        )
    if getattr(args, "bracken_precision", None) and screening != "kraken2_bracken":
        warnings.append(
            "--bracken-precision has no effect without --contaminant-screening kraken2_bracken."
        )


def _check_mutual_exclusions(args) -> None:
    if getattr(args, "transcript_fasta", None) and getattr(args, "additional_fasta", None):
        # Documented exception: nf-core/rnaseq allows both when a pre-built index
        # (--salmon-index or --kallisto-index) already incorporates the spike-ins,
        # so the pipeline never needs to merge the two FASTAs itself.
        has_prebuilt_index = any(getattr(args, f, None) for f in ("salmon_index", "kallisto_index"))
        if not has_prebuilt_index:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.CONFLICTING_FASTA_ARGS,
                message="--transcript-fasta and --additional-fasta cannot be used together.",
                fix=(
                    "Use only one transcript FASTA source, or supply a pre-built index "
                    "(--salmon-index or --kallisto-index) that already contains the spike-in sequences."
                ),
                details={"fields": ["transcript_fasta", "additional_fasta"]},
            )
    if getattr(args, "extra_star_align_args", None) and getattr(args, "aligner", "") != "star_salmon":
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INCOMPATIBLE_ALIGNER_ARGS,
            message="--extra-star-align-args is only supported with --aligner star_salmon.",
            fix=(
                "Switch to --aligner star_salmon, or remove --extra-star-align-args. "
                "For star_rsem, pass STAR options via rsem-calculate-expression's --star-options: "
                "process { withName: 'RSEM_CALCULATEEXPRESSION' { ext.args = '--star-options \"...\"' } }"
            ),
            details={"aligner": getattr(args, "aligner", "")},
        )
    if getattr(args, "extra_bowtie2_align_args", None) and getattr(args, "aligner", "") != "bowtie2_salmon":
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INCOMPATIBLE_ALIGNER_ARGS,
            message="--extra-bowtie2-align-args requires --aligner bowtie2_salmon.",
            fix="Switch to --aligner bowtie2_salmon or remove --extra-bowtie2-align-args.",
            details={"aligner": getattr(args, "aligner", "")},
        )
    if (
        getattr(args, "with_umi", False)
        and not getattr(args, "skip_umi_extract", False)
        and not getattr(args, "umitools_bc_pattern", None)
    ):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--with-umi requires --umitools-bc-pattern unless UMI extraction is skipped.",
            fix="Add --umitools-bc-pattern or pass --skip-umi-extract.",
            details={"field": "umitools_bc_pattern"},
        )


def _check_skip_alignment_requirements(args, *, bam_paths: list) -> None:
    if not getattr(args, "skip_alignment", False):
        return
    has_pseudo = bool(getattr(args, "pseudo_aligner", None))
    has_bams = bool(bam_paths)
    # Preprocessing-only mode: --skip_alignment + --skip_pseudo_alignment is a valid
    # documented pipeline mode that runs only FASTQ QC and trimming (no alignment or
    # pseudo-alignment). No pseudo-aligner or BAM paths are required in this case.
    skip_pseudo = bool(getattr(args, "skip_pseudo_alignment", False))
    if not has_pseudo and not has_bams and not skip_pseudo:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--skip-alignment requires either --pseudo-aligner, BAM paths in the samplesheet, or --skip-pseudo-alignment (preprocessing-only mode).",
            fix=(
                "For pseudo-alignment only: add --pseudo-aligner salmon or kallisto. "
                "For BAM reprocessing: provide a samplesheet with genome_bam or transcriptome_bam columns. "
                "For preprocessing only (trimming + QC): also add --skip-pseudo-alignment."
            ),
            details={"field": "pseudo_aligner"},
        )


def _check_samplesheet_summary(args, *, samplesheet_summary: dict[str, object]) -> None:
    for strandedness in samplesheet_summary.get("strandedness_counts", {}):
        if strandedness not in SUPPORTED_STRANDEDNESS:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_STRANDEDNESS,
                message="Samplesheet summary contains invalid strandedness.",
                fix="Use auto, forward, reverse, or unstranded in every samplesheet row.",
                details={"strandedness": strandedness},
            )
    for fastq_path in samplesheet_summary.get("fastq_paths", []):
        raw_path = str(fastq_path)
        path = Path(raw_path)
        if not _FASTQ_BASENAME_RE.match(path.name):
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_FASTQ,
                message="FASTQ basename is not compatible with nf-core/rnaseq.",
                fix="Use .fq, .fastq, .fq.gz, or .fastq.gz and remove whitespace from the basename.",
                details={"filename": path.name},
            )
        if "://" not in raw_path and not path.exists():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.MISSING_FASTQ,
                message="FASTQ path from samplesheet summary was not found.",
                fix="Correct the FASTQ path and rerun samplesheet validation.",
                details={"path": str(path)},
            )
    raw_bam_paths = [str(path) for path in samplesheet_summary.get("bam_paths", [])]
    if raw_bam_paths and not getattr(args, "skip_alignment", False):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.BAM_REPROCESSING_INCOMPLETE,
            message="BAM reprocessing rows require --skip-alignment.",
            fix="Add --skip-alignment, or remove BAM columns from the samplesheet.",
            details={"bam_count": len(raw_bam_paths)},
        )
    for raw_bam_path in raw_bam_paths:
        if "://" in raw_bam_path:
            continue
        bam_path = Path(raw_bam_path)
        if not bam_path.exists():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.REFERENCE_PATH_NOT_FOUND,
                message="BAM path from samplesheet summary was not found.",
                fix="Correct the BAM path in the samplesheet.",
                details={"field": "bam", "path": str(bam_path)},
            )


def _check_resume_compatibility(args, *, output_dir: Path, pipeline_source: dict[str, object]) -> None:
    if not getattr(args, "resume", False):
        return
    manifest_path = output_dir / "reproducibility" / "manifest.json"
    if not manifest_path.exists():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_RESUME_STATE,
            message="--resume requires an existing reproducibility manifest.",
            fix="Run without --resume in a new output directory, or resume a completed compatible run.",
            details={"manifest": str(manifest_path)},
        )
    try:
        manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    except json.JSONDecodeError as exc:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_RESUME_STATE,
            message="Resume manifest is not valid JSON.",
            fix="Use a valid output directory from a previous wrapper run.",
            details={"manifest": str(manifest_path), "error": str(exc)},
        ) from exc
    expected = {
        "aligner": getattr(args, "aligner", "star_salmon"),
        "pseudo_aligner": getattr(args, "pseudo_aligner", None),
        "profile": getattr(args, "profile", "docker"),
        "prokaryotic": bool(getattr(args, "prokaryotic", False)),
    }
    for key, value in expected.items():
        if manifest.get(key) != value:
            _raise_resume_mismatch(key, value, manifest.get(key))
    # `arm` was added to the manifest after v0.1.0; legacy manifests omit the key.
    # Treat absence as False so older runs can still resume; only block on real drift.
    current_arm = bool(getattr(args, "arm", False))
    manifest_arm = bool(manifest.get("arm", False))
    if current_arm != manifest_arm:
        _raise_resume_mismatch("arm", current_arm, manifest_arm)
    # `demo` composes the upstream `test` profile, which brings its own samplesheet AND
    # its own bundled references. Resuming across the demo/real boundary would swap both
    # underneath the run, so block it explicitly instead of letting it surface later as an
    # opaque params_checksum mismatch. Legacy manifests omit the key — treat absence as False.
    current_demo = bool(getattr(args, "demo", False))
    manifest_demo = bool(manifest.get("demo", False))
    if current_demo != manifest_demo:
        _raise_resume_mismatch("demo", current_demo, manifest_demo)
    previous_source = manifest.get("pipeline_source", {})
    for key in ("source_kind", "resolved_version"):
        if previous_source.get(key) != pipeline_source.get(key):
            _raise_resume_mismatch(f"pipeline_source.{key}", pipeline_source.get(key), previous_source.get(key))


def _raise_resume_mismatch(field: str, expected: object, observed: object) -> None:
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_RESUME_STATE,
        message="Resume manifest is incompatible with the requested run.",
        fix="Use a fresh output directory or rerun with matching pipeline/profile/alignment parameters.",
        details={"field": field, "expected": expected, "observed": observed},
    )


_UMI_DEPENDENT_FIELDS = (
    "umi_dedup_tool",
    "umitools_extract_method",
    "umitools_bc_pattern",
    "umitools_bc_pattern2",
    "umitools_umi_separator",
    "umi_discard_read",
    "umitools_grouping_method",
    "skip_umi_extract",
    "umitools_dedup_stats",
    "umitools_dedup_primary_only",
)


def _collect_umi_warnings(args, warnings: list[str]) -> None:
    """Warn when UMI-specific options are set but --with-umi is off.

    Without --with-umi the pipeline performs no UMI extraction/deduplication, so these
    options are silently inert. Surfacing it prevents a run that the user believes is
    UMI-deduplicated but is not (audit F-11).
    """
    if getattr(args, "with_umi", False):
        return
    set_fields = [f for f in _UMI_DEPENDENT_FIELDS if getattr(args, f, None)]
    if set_fields:
        flags = ", ".join("--" + f.replace("_", "-") for f in set_fields)
        warnings.append(
            f"UMI options ({flags}) have no effect without --with-umi; no UMI "
            "deduplication will be performed."
        )


def _collect_demo_warnings(args, warnings: list[str]) -> None:
    """Emit a structured demo warning when reference/aligner overrides are still present.

    In the normal flow ``_apply_demo_overrides`` clears reference flags before preflight
    runs, so this collector adds nothing — the human-readable warning already went to
    stderr from the override step. When preflight is invoked directly (e.g. from tests
    or external orchestrators that bypass the wrapper's main()), the args still carry
    overrides and this surfaces a single structured warning in ``result["warnings"]``.
    """
    if not getattr(args, "demo", False):
        return
    reference_override = any(getattr(args, field, None) for field in _EXPLICIT_REFERENCE_FIELDS) or getattr(args, "genome", None)
    if getattr(args, "aligner", "star_salmon") != "star_salmon" or reference_override:
        warnings.append(
            "Demo mode clears user input/reference overrides and runs the upstream test profile with star_salmon."
        )


def _collect_platform_warnings(args, output_dir: Path, warnings: list[str]) -> None:
    profile_parts = {p.strip() for p in getattr(args, "profile", "").split(",") if p.strip()}
    if sys.platform == "darwin" and "docker" in profile_parts:
        # Colima (a common macOS Docker runtime) only shares the user HOME into
        # its VM; /tmp and /private/tmp live on the VM's own filesystem and are
        # NOT shared with the host, so Nextflow's work-dir files are invisible to
        # containers and the run fails with '.command.run: No such file or
        # directory'. Use a resolve-based check and the accurate message so this
        # matches the nfcore-sarek / nfcore-scrnaseq guards (not the old, wrong
        # "may be slow / VirtioFS" wording, which describes a non-existent failure
        # mode and would not steer the user away from /tmp).
        if is_under_tmp(output_dir):
            warnings.append(
                "Output directory is under /tmp. On macOS with Colima, Docker "
                "containers cannot see files written to /tmp (the VM uses its own "
                "separate /tmp). Move --output under your home directory to avoid "
                "'No such file or directory' errors."
            )
        _warn_macos_star_index_memory(args, warnings)


def _warn_macos_star_index_memory(args, warnings: list[str]) -> None:
    """Warn when a STAR index will be built from --fasta under the macOS memory cap.

    STAR_GENOMEGENERATE for a large genome needs ~30+ GB, but the macOS Docker
    compatibility config caps per-process memory (~15 GB), so the build can be
    OOM-killed. Only relevant when a genome FASTA is provided and the aligner's
    prebuilt genome index is absent (audit L2). Using --genome (iGenomes) or a
    prebuilt index needs no local build, so no warning fires.
    """
    aligner = getattr(args, "aligner", "")
    if aligner not in {"star_salmon", "star_rsem"} or not getattr(args, "fasta", None):
        return
    prebuilt = getattr(args, "star_index", None)
    if aligner == "star_rsem":
        prebuilt = prebuilt or getattr(args, "rsem_index", None)
    if prebuilt:
        return
    warnings.append(
        "On macOS Docker, building the STAR index from --fasta for a large genome can exceed the "
        "wrapper's per-process memory ceiling (~15 GB) and be OOM-killed. Provide a prebuilt "
        "--star-index, use --genome with a local iGenomes mirror, or raise the Docker VM memory."
    )


def _collect_index_aligner_warnings(args, warnings: list[str]) -> None:
    """Warn when an index flag the chosen aligner/pseudo-aligner cannot consume is set.

    The wrapper validates index presence/completeness but not index *type vs aligner*.
    e.g. ``--aligner star_rsem --salmon-index`` is accepted, yet RSEM never uses Salmon,
    so the index is silently ignored by nf-core. Surface it so the user does not believe a
    prebuilt index took effect when it did not (audit F3). A warning, not an error — the
    pipeline simply ignores the unused index.
    """
    aligner = getattr(args, "aligner", "star_salmon")
    pseudo = getattr(args, "pseudo_aligner", None)
    relevant: set[str] = set(_ALIGNER_GENOME_INDEX_FIELDS.get(aligner, ()))
    if aligner in {"star_salmon", "bowtie2_salmon"}:
        relevant.add("salmon_index")
    if pseudo == "salmon":
        relevant.add("salmon_index")
    if pseudo == "kallisto":
        relevant.add("kallisto_index")
    unused = [f for f in _INDEX_REFERENCE_FIELDS if f not in relevant and getattr(args, f, None)]
    if unused:
        flags = ", ".join("--" + f.replace("_", "-") for f in unused)
        pseudo_suffix = f" --pseudo-aligner {pseudo}" if pseudo else ""
        warnings.append(
            f"--aligner {aligner}{pseudo_suffix} does not use the supplied index flag(s) "
            f"({flags}); they will be ignored. Provide an index matching the chosen aligner, "
            "or remove the unused flag(s)."
        )


def _collect_bam_reprocessing_warning(bam_paths: list, aligner_effective: str, warnings: list[str]) -> None:
    """Remind the user that BAM reprocessing must reuse the original aligner.

    nf-core/rnaseq cannot mix quantifier types between BAM generation and reprocessing:
    star_salmon BAMs must be reprocessed with star_salmon, star_rsem with star_rsem.
    The wrapper cannot infer how a BAM was produced, and the default aligner is
    star_salmon, so it warns whenever BAM columns are present to prevent a silent
    quantifier mismatch.
    """
    if not bam_paths:
        return
    warnings.append(
        f"BAM reprocessing detected: ensure --aligner ({aligner_effective}) matches the aligner used to "
        "generate the BAMs. nf-core/rnaseq cannot mix quantifier types — star_salmon BAMs must be "
        "reprocessed with star_salmon, star_rsem BAMs with star_rsem."
    )


def _collect_alignment_warnings(args, warnings: list[str]) -> bool:
    handoff_available = True
    if getattr(args, "prokaryotic", False) and getattr(args, "aligner", "") != "bowtie2_salmon":
        warnings.append("Prokaryotic mode is best paired with --aligner bowtie2_salmon.")
    if (
        getattr(args, "aligner", "") == "hisat2"
        and getattr(args, "run_downstream", False)
        and not getattr(args, "pseudo_aligner", None)
    ):
        warnings.append("HISAT2 route does not produce merged quantification suitable for rnaseq-de handoff.")
        handoff_available = False
    if getattr(args, "rsem_extra_args", None):
        warnings.append(
            "--rsem-extra-args has no effect: nf-core/rnaseq does not expose an extra_rsem_quant_args parameter. "
            "To customise RSEM arguments use a Nextflow config override: "
            "process { withName: 'RSEM_CALCULATEEXPRESSION' { ext.args = '...' } }"
        )
    return handoff_available
