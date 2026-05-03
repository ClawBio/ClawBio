from __future__ import annotations

import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from errors import ErrorCode, SkillError
from schemas import (
    COMMON_GENOMES,
    JAVA_MIN_VERSION,
    NEXTFLOW_MIN_VERSION,
    PRESET_REQUIREMENTS,
    SUPPORTED_PRESETS,
    SUPPORTED_PROFILES,
)


_SUBPROCESS_TIMEOUT = 30
_FASTA_SCHEMA_RE = re.compile(r"^\S+\.fn?a(sta)?(\.gz)?$")
_EMAIL_SCHEMA_RE = re.compile(r"^([a-zA-Z0-9_\-.]+)@([a-zA-Z0-9_\-.]+)\.([a-zA-Z]{2,5})$")


def _command_output(args: list[str]) -> str:
    try:
        proc = subprocess.run(
            args, capture_output=True, text=True, errors="replace", timeout=_SUBPROCESS_TIMEOUT
        )
    except (subprocess.TimeoutExpired, OSError):
        return ""
    if proc.returncode != 0:
        return ""
    return (proc.stdout or proc.stderr).strip()


def _pad_version(t: tuple[int, ...], length: int = 3) -> tuple[int, ...]:
    return t + (0,) * max(0, length - len(t))


def _parse_version_tuple(text: str) -> tuple[int, ...]:
    m = re.search(r'\b(\d+)\.(\d+)\.(\d+)\b', text)
    if m:
        return (int(m.group(1)), int(m.group(2)), int(m.group(3)))
    m = re.search(r'\b(\d+)\.(\d+)\b', text)
    if m:
        return (int(m.group(1)), int(m.group(2)))
    m = re.search(r'\b(\d+)\b', text)
    if m:
        return (int(m.group(1)),)
    return ()


def _check_executable(name: str) -> str:
    path = shutil.which(name)
    if not path:
        raise SkillError(
            stage="preflight",
            error_code=f"MISSING_{name.upper().replace('-', '_')}",
            message=f"Required executable `{name}` was not found.",
            fix=f"Install `{name}` and ensure it is available on PATH.",
            details={"executable": name},
        )
    return path


def _check_java() -> dict[str, str]:
    java_path = _check_executable("java")
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
    if version_tuple[0] < JAVA_MIN_VERSION:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.JAVA_VERSION_TOO_OLD,
            message="Java version is too old for this wrapper.",
            fix="Install Java 17 or newer.",
            details={"detected_version": ".".join(map(str, version_tuple))},
        )
    return {"path": java_path, "version": ".".join(map(str, version_tuple))}


def _check_nextflow() -> dict[str, str]:
    nextflow_path = _check_executable("nextflow")
    version_text = _command_output(["nextflow", "-version"])
    version_tuple = _parse_version_tuple(version_text)
    if not version_tuple:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_NEXTFLOW,
            message="Nextflow is installed but its version could not be determined.",
            fix="Install Nextflow 25.04.0 or newer and ensure `nextflow -version` works.",
            details={"nextflow_path": nextflow_path},
        )
    if _pad_version(version_tuple) < _pad_version(NEXTFLOW_MIN_VERSION):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.NEXTFLOW_VERSION_TOO_OLD,
            message="Nextflow version is too old for this wrapper.",
            fix="Upgrade Nextflow to 25.04.0 or newer.",
            details={"detected_version": ".".join(map(str, version_tuple))},
        )
    return {"path": nextflow_path, "version": ".".join(map(str, version_tuple))}


def _check_profile(profile: str) -> dict[str, str | bool]:
    if profile not in SUPPORTED_PROFILES:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PROFILE,
            message="Unsupported execution profile.",
            fix=f"Choose one of: {', '.join(sorted(SUPPORTED_PROFILES))}.",
            details={"profile": profile},
        )
    if profile == "docker":
        return _check_docker_profile(profile)
    if profile == "conda":
        return _check_conda_profile(profile)
    if profile == "mamba":
        return _check_mamba_profile(profile)
    if profile == "podman":
        return _check_podman_profile(profile)
    if profile in {"shifter", "charliecloud"}:
        return _check_hpc_profile(profile)
    return _check_singularity_compatible_profile(profile)


def _check_docker_profile(profile: str) -> dict[str, str | bool]:
    docker_path = shutil.which("docker")
    if not docker_path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_DOCKER,
            message="Docker profile was selected but Docker is not installed.",
            fix="Install Docker or choose another supported profile.",
            details={"profile": profile},
        )
    try:
        info = subprocess.run(["docker", "info"], capture_output=True, text=True, timeout=_SUBPROCESS_TIMEOUT)
        docker_ok = info.returncode == 0
    except (subprocess.TimeoutExpired, OSError):
        docker_ok = False
    if not docker_ok:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.DOCKER_NOT_RUNNING,
            message="Docker is installed but the daemon is not available.",
            fix="Start Docker Desktop or the Docker daemon before running this skill.",
            details={"profile": profile},
        )
    return {"profile": profile, "backend_path": docker_path, "backend_ready": True}


def _check_conda_profile(profile: str) -> dict[str, str | bool]:
    # Prefer conda (matches the profile name); fall back to mamba.
    # _check_mamba_profile uses the opposite order: mamba-first, conda-fallback.
    backend = shutil.which("conda") or shutil.which("mamba")
    if not backend:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_CONDA,
            message="Conda profile was selected but neither conda nor mamba is installed.",
            fix="Install conda or mamba, or choose another profile.",
            details={"profile": profile},
        )
    return {"profile": profile, "backend_path": backend, "backend_ready": True}


def _check_singularity_compatible_profile(profile: str) -> dict[str, str | bool]:
    # Singularity and Apptainer are API-compatible; accept either binary for either profile.
    # Many modern HPC clusters ship only one of the two, or renamed the binary during the
    # SingularityCE → Apptainer transition.
    primary = "apptainer" if profile == "apptainer" else "singularity"
    fallback = "singularity" if profile == "apptainer" else "apptainer"
    backend = shutil.which(primary) or shutil.which(fallback)
    if not backend:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_SINGULARITY,
            message=(
                f"{profile} profile was selected but neither `{primary}` nor `{fallback}` "
                "was found on PATH."
            ),
            fix=(
                f"Install Singularity or Apptainer and ensure it is available on PATH, "
                f"or choose a different profile (docker, conda)."
            ),
            details={"profile": profile, "tried": [primary, fallback]},
        )
    return {"profile": profile, "backend_path": backend, "backend_ready": True}


def _check_mamba_profile(profile: str) -> dict[str, str | bool]:
    backend = shutil.which("mamba") or shutil.which("conda")
    if not backend:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_CONDA,
            message="Mamba profile was selected but neither mamba nor conda is installed.",
            fix="Install mamba or conda, or choose another profile.",
            details={"profile": profile},
        )
    return {"profile": profile, "backend_path": backend, "backend_ready": True}


def _check_podman_profile(profile: str) -> dict[str, str | bool]:
    podman_path = shutil.which("podman")
    if not podman_path:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.MISSING_PODMAN,
            message="Podman profile was selected but Podman is not installed.",
            fix="Install Podman or choose another supported profile.",
            details={"profile": profile},
        )
    try:
        info = subprocess.run(["podman", "info"], capture_output=True, text=True, timeout=_SUBPROCESS_TIMEOUT)
        podman_ok = info.returncode == 0
    except (subprocess.TimeoutExpired, OSError):
        podman_ok = False
    if not podman_ok:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.PODMAN_NOT_RUNNING,
            message="Podman is installed but the service is not available.",
            fix="Start the Podman service or socket before running this skill.",
            details={"profile": profile},
        )
    return {"profile": profile, "backend_path": podman_path, "backend_ready": True}


def _check_hpc_profile(profile: str) -> dict[str, str | bool]:
    # HPC runtimes (shifter, charliecloud) are user-space tools that don't run a
    # persistent daemon, so binary presence is the only liveness check needed.
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


_IGNORED_ROOT_NAMES = frozenset({".DS_Store", ".gitkeep", ".gitignore", "Thumbs.db", "check_result.json"})
_ALLOWED_REPRO_FILES = frozenset({"samplesheet.valid.csv", "samplesheet.demo.csv", "params.yaml"})


def _check_output_dir(output_dir: Path, *, resume: bool) -> None:
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
    materialized_entries = []
    for entry in output_dir.iterdir():
        if entry.name in _IGNORED_ROOT_NAMES:
            continue
        if entry.name == "reproducibility" and entry.is_dir():
            repro_entries = [
                child for child in entry.iterdir()
                if child.name not in _ALLOWED_REPRO_FILES and child.name not in _IGNORED_ROOT_NAMES
            ]
            if repro_entries:
                materialized_entries.append(entry)
            continue
        materialized_entries.append(entry)
    if materialized_entries and not resume:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.OUTPUT_DIR_NOT_EMPTY,
            message="Output directory already contains files.",
            fix=(
                "Choose a new empty output directory, or re-run with --resume only if the previous run "
                "completed successfully (manifest.json must exist)."
            ),
            details={"output_dir": str(output_dir)},
        )


def check_output_dir_available(output_dir: Path, *, resume: bool) -> None:
    """Validate output-dir reuse policy before writing wrapper artifacts."""
    _check_output_dir(output_dir, resume=resume)


# Fields that are string identifiers, not local paths — never existence-checked.
_SYMBOLIC_REFERENCE_FIELDS = frozenset({"genome", "cellrangerarc_reference"})

# Primary genome/index fields used to satisfy preset reference requirements.
_PRIMARY_REFERENCE_FIELDS = (
    "fasta", "gtf", "transcript_fasta", "txp2gene",
    "simpleaf_index", "kallisto_index", "star_index", "cellranger_index", "barcode_whitelist",
)

# Additional optional file paths — existence-checked when provided, but not part of requires_any.
_OPTIONAL_PATH_FIELDS = (
    "kb_t1c", "kb_t2c",
    "motifs", "cellrangerarc_config",
    "cellranger_vdj_index",
    "gex_frna_probe_set", "gex_target_panel", "gex_cmo_set", "fb_reference",
    "vdj_inner_enrichment_primers", "gex_barcode_sample_assignment", "cellranger_multi_barcodes",
)

_EXPLICIT_REFERENCE_FIELDS = _PRIMARY_REFERENCE_FIELDS + _OPTIONAL_PATH_FIELDS


def _check_references(args) -> dict[str, str]:
    resolved = _collect_reference_values(args)
    _reject_conflicting_reference_styles(resolved)
    satisfied_group = _find_satisfied_reference_group(args.preset, resolved)
    _check_preset_specific_requirements(args.preset, resolved, satisfied_group)
    _check_reference_paths_exist(resolved)
    return resolved


def _collect_reference_values(args) -> dict[str, str]:
    genome = getattr(args, "genome", None) or ""
    cellrangerarc_reference = getattr(args, "cellrangerarc_reference", None) or ""
    resolved = {"genome": genome, "cellrangerarc_reference": cellrangerarc_reference}
    for field in _EXPLICIT_REFERENCE_FIELDS:
        resolved[field] = getattr(args, field, None) or ""
    return resolved


def _reject_conflicting_reference_styles(resolved: dict[str, str]) -> None:
    if resolved["genome"] and any(resolved[f] for f in _PRIMARY_REFERENCE_FIELDS if f != "barcode_whitelist"):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.CONFLICTING_REFERENCES,
            message="--genome (iGenomes shortcut) and explicit reference paths are mutually exclusive.",
            fix="Use either --genome <shortcut> or explicit --fasta/--gtf/--index flags, not both.",
            details={"genome": resolved["genome"]},
        )


def _find_satisfied_reference_group(preset: str, resolved: dict[str, str]) -> tuple[str, ...]:
    requirements = PRESET_REQUIREMENTS[preset]["requires_any"]
    for option_group in requirements:
        if all(resolved.get(name, "") for name in option_group):
            return option_group
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.MISSING_REFERENCE,
        message="The selected preset is missing required references or indexes.",
        fix="Provide one of the supported reference/index combinations for this preset.",
        details={"preset": preset, "accepted_combinations": [list(group) for group in requirements]},
    )


def _check_preset_specific_requirements(
    preset: str,
    resolved: dict[str, str],
    satisfied_group: tuple[str, ...],
) -> None:
    if preset == "cellrangerarc" and satisfied_group != ("cellranger_index",):
        _check_cellrangerarc_config_pairing(resolved)
    requires_also = PRESET_REQUIREMENTS[preset].get("requires_also", [])
    _require_additional_preset_fields(preset, resolved, requires_also, reason="for this preset")


def _check_cellrangerarc_config_pairing(resolved: dict[str, str]) -> None:
    """Mirror upstream ARC rules: motifs are optional, config/reference are paired."""
    has_config = bool(resolved.get("cellrangerarc_config", ""))
    has_reference = bool(resolved.get("cellrangerarc_reference", ""))
    if has_config == has_reference:
        return
    missing_field = "cellrangerarc_reference" if has_config else "cellrangerarc_config"
    flag = "--" + missing_field.replace("_", "-")
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message="CellRanger ARC custom config and reference name must be provided together.",
        fix=(
            f"Add {flag}, or omit both ARC config fields and let the pipeline build the reference "
            "from --fasta/--gtf or --genome."
        ),
        details={"preset": "cellrangerarc", "missing_field": missing_field},
    )


def _require_additional_preset_fields(
    preset: str,
    resolved: dict[str, str],
    fields: list[str],
    *,
    reason: str,
) -> None:
    for field in fields:
        if not resolved.get(field, ""):
            flag = "--" + field.replace("_", "-")
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message=f"The {preset!r} preset requires {flag} {reason}.",
                fix=f"Add {flag} to your command. See the skill documentation for details.",
                details={"preset": preset, "missing_field": field},
            )


def _check_reference_paths_exist(resolved: dict[str, str]) -> None:
    for key, value in resolved.items():
        if key in _SYMBOLIC_REFERENCE_FIELDS:
            continue  # symbolic identifiers, not local paths
        _check_upstream_path_schema_compatibility(key, value)
        if value and not Path(value).expanduser().exists():
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.MISSING_REFERENCE,
                message="A required reference or index path was not found.",
                fix="Correct the missing reference path and try again.",
                details={"field": key, "path": value},
            )


def _check_upstream_path_schema_compatibility(key: str, value: str) -> None:
    if key != "fasta" or not value:
        return
    fasta_path = Path(value).expanduser().as_posix()
    if _FASTA_SCHEMA_RE.match(fasta_path):
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message="FASTA paths must match the nf-core/scrnaseq 4.1.0 schema.",
        fix=(
            "Use a whitespace-free FASTA path ending in .fa, .fna, .fasta, .fa.gz, "
            ".fna.gz, or .fasta.gz, then pass it via --fasta."
        ),
        details={"field": key, "path": value, "schema_pattern": r"^\S+\.fn?a(sta)?(\.gz)?$"},
    )


def _check_parameter_schema_compatibility(args) -> None:
    """Validate non-path parameters whose upstream schema would fail predictably."""
    _check_email_schema_compatibility(getattr(args, "email", None))


def _check_email_schema_compatibility(email: str | None) -> None:
    if not email or _EMAIL_SCHEMA_RE.match(email):
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message="Email address does not match the nf-core/scrnaseq 4.1.0 schema.",
        fix="Use a simple email address such as user@example.org, or omit --email.",
        details={"field": "email", "value": email, "schema_pattern": _EMAIL_SCHEMA_RE.pattern},
    )


def run_preflight(
    args,
    *,
    pipeline_source: dict[str, str | bool],
    samplesheet_summary: dict[str, object],
) -> dict[str, object]:
    _warn_if_native_windows()
    _check_supported_preset(args.preset)
    _check_parameter_schema_compatibility(args)
    _check_protocol_compatibility(args)
    java_info = _check_java()
    nextflow_info = _check_nextflow()
    profile_info = _check_profile(args.profile)
    output_dir = Path(args.output).expanduser().resolve()
    # Belt-and-suspenders: main() already called check_output_dir_available() before
    # samplesheet normalization. This second call keeps run_preflight() self-consistent
    # when invoked directly (e.g., tests). reproducibility/samplesheet.valid.csv written
    # between the two calls is allowlisted in _ALLOWED_REPRO_FILES and does not trigger
    # OUTPUT_DIR_NOT_EMPTY.
    check_output_dir_available(output_dir, resume=args.resume)
    refs = {} if args.demo else _check_references(args)
    _check_samplesheet_driven_preset_requirements(args, samplesheet_summary=samplesheet_summary)
    _check_resume_compatibility(args, output_dir=output_dir, pipeline_source=pipeline_source)
    _warn_if_macos_docker_tmp(args.profile, output_dir)

    return {
        "ok": True,
        "java": java_info,
        "nextflow": nextflow_info,
        "profile": profile_info,
        "pipeline_source": pipeline_source,
        "references": refs,
        "samplesheet": {
            "sample_count": samplesheet_summary["sample_count"],
            "unknown_columns": samplesheet_summary["unknown_columns"],
            "sample_types": samplesheet_summary.get("sample_types", []),
            "feature_types": samplesheet_summary.get("feature_types", []),
        },
    }


def _check_samplesheet_driven_preset_requirements(args, *, samplesheet_summary: dict[str, object]) -> None:
    if args.demo or args.preset != "cellrangermulti":
        return
    feature_types = {str(value).lower() for value in samplesheet_summary.get("feature_types", [])}
    _check_cellrangermulti_feature_references(args, feature_types)
    _check_cellrangermulti_vdj_reference_policy(args, feature_types)
    _check_cellrangermulti_multiplexing_policy(args, feature_types)


def _check_cellrangermulti_feature_references(args, feature_types: set[str]) -> None:
    if "ab" in feature_types and not getattr(args, "fb_reference", None):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti antibody-capture rows require a feature-barcode reference.",
            fix="Provide --fb-reference for feature_type=ab rows.",
            details={"preset": args.preset, "missing_field": "fb_reference", "feature_type": "ab"},
        )


def _check_cellrangermulti_vdj_reference_policy(args, feature_types: set[str]) -> None:
    if "vdj" not in feature_types or not getattr(args, "skip_cellrangermulti_vdjref", False):
        return
    if getattr(args, "cellranger_vdj_index", None):
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message="cellrangermulti VDJ rows need a VDJ reference unless mkvdjref is allowed to run.",
        fix="Remove --skip-cellrangermulti-vdjref, or provide --cellranger-vdj-index for feature_type=vdj rows.",
        details={
            "preset": args.preset,
            "feature_type": "vdj",
            "missing_field": "cellranger_vdj_index",
            "conflicting_flag": "skip_cellrangermulti_vdjref",
        },
    )


def _check_cellrangermulti_multiplexing_policy(args, feature_types: set[str]) -> None:
    has_cmo_rows = "cmo" in feature_types
    has_ffpe_probe_set = bool(getattr(args, "gex_frna_probe_set", None))
    has_cmo_reference = bool(getattr(args, "gex_cmo_set", None))
    has_ocm_assignment = bool(getattr(args, "gex_barcode_sample_assignment", None))
    multiplexing_modes = {
        "cmo": has_cmo_rows or has_cmo_reference,
        "ffpe": has_ffpe_probe_set,
        "ocm": has_ocm_assignment,
    }
    active_modes = [mode for mode, active in multiplexing_modes.items() if active]
    if len(active_modes) > 1:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti multiplexing modes are mutually exclusive.",
            fix="Use only one multiplexing strategy per run: CMO, FFPE probe-set demultiplexing, or OCM assignment.",
            details={"preset": args.preset, "active_modes": active_modes},
        )
    if has_ffpe_probe_set and not getattr(args, "cellranger_multi_barcodes", None):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti FFPE probe-set demultiplexing requires a barcode-to-sample samplesheet.",
            fix="Provide --cellranger-multi-barcodes when using --gex-frna-probe-set.",
            details={"preset": args.preset, "missing_field": "cellranger_multi_barcodes", "field": "gex_frna_probe_set"},
        )
    if "cmo" in feature_types and not getattr(args, "cellranger_multi_barcodes", None):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti CMO multiplexing requires a barcode-to-sample samplesheet.",
            fix="Provide --cellranger-multi-barcodes for feature_type=cmo rows.",
            details={"preset": args.preset, "missing_field": "cellranger_multi_barcodes", "feature_type": "cmo"},
        )


def _check_protocol_compatibility(args) -> None:
    if getattr(args, "demo", False):
        return
    protocol = (getattr(args, "protocol", None) or "auto").strip()
    preset = getattr(args, "preset", "")
    if preset in {"standard", "star", "kallisto"} and protocol == "auto":
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="The selected preset requires an explicit non-auto protocol.",
            fix=(
                "Provide --protocol (for example 10XV2, 10XV3, dropseq, "
                "or a custom chemistry string supported by the aligner)."
            ),
            details={"preset": preset, "protocol": protocol},
        )
    if preset == "standard" and protocol.lower() == "smartseq":
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="The standard preset uses Simpleaf, which does not support smartseq in nf-core/scrnaseq.",
            fix="Use --preset star or --preset kallisto for smartseq, or choose a Simpleaf-supported protocol.",
            details={"preset": preset, "protocol": protocol},
        )


def _prompt_for_genome() -> str | None:
    """Interactively ask the user to pick a genome from COMMON_GENOMES.

    Returns the genome string, or None if stdin is not a TTY or the user
    presses Enter / supplies an unrecognised selection.
    """
    if not sys.stdin.isatty():
        return None
    lines = ["\nNo reference genome specified. Choose a common iGenomes shortcut:"]
    for i, g in enumerate(COMMON_GENOMES, start=1):
        lines.append(f"  {i:2d}. {g}")
    lines.append("\nEnter number or genome name (or press Enter to skip): ")
    try:
        raw = input("".join(lines)).strip()
    except (EOFError, KeyboardInterrupt):
        return None
    if not raw:
        return None
    if raw.isdigit():
        idx = int(raw) - 1
        if 0 <= idx < len(COMMON_GENOMES):
            return COMMON_GENOMES[idx]
        return None
    if raw in COMMON_GENOMES:
        return raw
    return None


def _warn_if_macos_docker_tmp(profile: str, output_dir: Path) -> None:
    # Colima (a common macOS Docker runtime) only 9p-mounts the user HOME directory
    # into its VM.  /tmp and /private/tmp live on the VM's own ext4 and are NOT shared
    # with the host, so Nextflow's work-dir files are invisible to containers.
    # Docker Desktop handles /tmp differently and is unlikely to hit this, but the safe
    # guidance for macOS + any Docker backend is to keep output dirs under HOME.
    if sys.platform != "darwin" or profile != "docker":
        return
    private_tmp = Path("/private/tmp").resolve()
    tmp = Path("/tmp").resolve()
    try:
        out = output_dir.resolve()
        under_tmp = out == tmp or out == private_tmp or tmp in out.parents or private_tmp in out.parents
    except Exception:
        return
    if under_tmp:
        print(
            "WARNING: Output directory is under /tmp. On macOS with Colima, Docker containers "
            "cannot see files written to /tmp (the VM uses its own separate /tmp). "
            "Move --output to a path under your home directory to avoid 'No such file or directory' errors.",
            file=sys.stderr,
        )


def _warn_if_native_windows() -> None:
    if sys.platform != "win32":
        return
    # Nextflow is not officially supported on native Windows.
    # The recommended path is WSL2, which reports platform as Linux.
    print(
        "WARNING: Running on native Windows. Nextflow is not officially supported "
        "outside WSL2 on Windows. If you encounter issues, run from WSL2 instead.",
        file=sys.stderr,
    )


def _check_supported_preset(preset: str) -> None:
    if preset not in SUPPORTED_PRESETS:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.UNSUPPORTED_MODE,
            message="Unsupported preset requested.",
            fix=f"Choose one of: {', '.join(sorted(SUPPORTED_PRESETS))}.",
            details={"preset": preset},
        )


def _check_resume_compatibility(
    args,
    *,
    output_dir: Path,
    pipeline_source: dict[str, str | bool],
) -> None:
    if not args.resume:
        return
    manifest_path = output_dir / "reproducibility" / "manifest.json"
    if not manifest_path.exists():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_RESUME_STATE,
            message="Resume was requested but no previous manifest was found.",
            fix="Remove --resume or point --output to a compatible prior run directory.",
            details={"manifest": str(manifest_path)},
        )
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    _check_resume_preset_and_profile(args, manifest)
    _check_resume_pipeline_source(pipeline_source, manifest)


def _check_resume_preset_and_profile(args, manifest: dict[str, object]) -> None:
    if manifest.get("preset") == args.preset and manifest.get("profile") == args.profile:
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_RESUME_STATE,
        message="Resume state does not match the requested preset/profile.",
        fix="Use the same preset/profile as the original run or start in a new output directory.",
        details={
            "previous_preset": manifest.get("preset"),
            "previous_profile": manifest.get("profile"),
            "requested_preset": args.preset,
            "requested_profile": args.profile,
        },
    )


def _check_resume_pipeline_source(
    pipeline_source: dict[str, str | bool],
    manifest: dict[str, object],
) -> None:
    previous_source = manifest.get("pipeline_source", {})
    if not isinstance(previous_source, dict):
        previous_source = {}
    if (
        previous_source.get("source_kind") == pipeline_source.get("source_kind")
        and previous_source.get("resolved_version") == pipeline_source.get("resolved_version")
    ):
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_RESUME_STATE,
        message="Resume state does not match the requested pipeline source.",
        fix="Resume only with the same pipeline source/ref as the original run.",
        details={
            "previous_source": previous_source,
            "requested_source": pipeline_source,
        },
    )
