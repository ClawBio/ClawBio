from __future__ import annotations

from clawbio.common.nfcore import runtime as _runtime

import csv
import json
import os
import re
import shutil
import subprocess
import sys
from pathlib import Path
from typing import Any

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("errors", "nfcore_4_1_0_contract", "schemas")

from errors import ErrorCode, SkillError
from nfcore_4_1_0_contract import (
    CELLRANGER_FAMILY_PRESETS,
    KNOWN_PROTOCOL_TOKENS,
    POLICY_SOURCE_CLAWBIO,
    POLICY_SOURCE_NFCORE_DOCS,
    PRESETS_REQUIRING_EXPLICIT_PROTOCOL,
    PRESETS_SUPPORTING_CUSTOM_PROTOCOL,
    PRESETS_SUPPORTING_SMARTSEQ_PROTOCOL,
    PROTOCOLS_JSON_4_1_0,
)
from schemas import (
    ALL_REFERENCE_PATH_FIELDS,
    GENOME_REFERENCE_FIELDS,
    JAVA_MIN_VERSION,
    NEXTFLOW_MIN_VERSION,
    NEXTFLOW_MIN_VERSION_DISPLAY,
    PRESET_REQUIREMENTS,
    PROJECT_ROOT,
    SUPPORTED_PRESETS,
    SUPPORTED_PROFILES,
    SYMBOLIC_REFERENCE_FIELDS,
    is_under_tmp,
    profile_components,
    profile_includes,
)


_SUBPROCESS_TIMEOUT = 60
_FASTA_SCHEMA_RE = re.compile(r"^\S+\.fn?a(sta)?(\.gz)?$")
_EMAIL_SCHEMA_RE = re.compile(
    r"^([a-zA-Z0-9_\-.]+)@([a-zA-Z0-9_\-.]+)\.([a-zA-Z]{2,5})$"
)
_PROFILE_COMPONENT_RE = re.compile(r"^[A-Za-z0-9_.-]+$")


def _command_output(args: list[str]) -> str:
    proc = _runtime.run_probe(args, timeout_seconds=_SUBPROCESS_TIMEOUT)
    if proc is None:
        return ""
    if proc.returncode != 0:
        return ""
    return (proc.stdout or proc.stderr).strip()


def _pad_version(t: tuple[int, ...], length: int = 3) -> tuple[int, ...]:
    return _runtime.pad_version(t, length)


def _parse_version_tuple(text: str) -> tuple[int, ...]:
    return _runtime.parse_version_tuple(text)


def _detected_version_string(text: str) -> str:
    """Return the version exactly as reported by the tool (e.g. '25.04.0').

    The parsed tuple is for *comparison* only; reconstructing a string from it
    drops zero-padding (25.04.0 → 25.4.0), which is not a real Nextflow release
    and would break NXF_VER / conda pins on replay. Always store the raw token.
    """
    for pattern in (r"\b(\d+\.\d+\.\d+)\b", r"\b(\d+\.\d+)\b", r"\b(\d+)\b"):
        m = re.search(pattern, text)
        if m:
            return m.group(1)
    return ""


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
        # java IS present but its version can't be confirmed — this is a
        # version-gate failure, not an absent binary (MISSING_JAVA).
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.JAVA_VERSION_TOO_OLD,
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
            details={"detected_version": _detected_version_string(version_text)},
        )
    return {"path": java_path, "version": _detected_version_string(version_text)}


def _check_nextflow() -> dict[str, str]:
    nextflow_path = _check_executable("nextflow")
    version_text = _command_output(["nextflow", "-version"])
    version_tuple = _parse_version_tuple(version_text)
    if not version_tuple:
        # nextflow IS present but its version can't be confirmed — version-gate
        # failure, not an absent binary (MISSING_NEXTFLOW).
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.NEXTFLOW_VERSION_TOO_OLD,
            message="Nextflow is installed but its version could not be determined.",
            fix=f"Install Nextflow {NEXTFLOW_MIN_VERSION_DISPLAY} or newer and ensure `nextflow -version` works.",
            details={"nextflow_path": nextflow_path},
        )
    if _pad_version(version_tuple) < _pad_version(NEXTFLOW_MIN_VERSION):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.NEXTFLOW_VERSION_TOO_OLD,
            message="Nextflow version is too old for this wrapper.",
            fix=f"Upgrade Nextflow to {NEXTFLOW_MIN_VERSION_DISPLAY} or newer.",
            details={"detected_version": _detected_version_string(version_text)},
        )
    return {"path": nextflow_path, "version": _detected_version_string(version_text)}


def _check_profile(profile: str) -> dict[str, object]:
    components = _profile_components(profile)
    if len(components) > 1:
        return _check_composite_profile(profile, components)
    profile = components[0] if components else profile
    if profile not in SUPPORTED_PROFILES:
        return _check_institutional_profile_component(profile, full_profile=profile)
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
    if profile in {"wave", "gpu"}:
        # Wave and GPU are Nextflow-native features, not external runtimes.
        # No binary check is needed; Nextflow handles them internally.
        return {"profile": profile, "backend_path": None, "backend_ready": True}
    if profile in {
        "debug",
        "arm64",
        "emulate_amd64",
        "test",
        "test_full",
        "test_cellrangermulti",
        "test_multiome",
    }:
        return {"profile": profile, "backend_path": None, "backend_ready": True}
    return _check_singularity_compatible_profile(profile)


def _profile_components(profile: str) -> list[str]:
    return profile_components(profile)


def _check_composite_profile(profile: str, components: list[str]) -> dict[str, object]:
    invalid = [
        component
        for component in components
        if component not in SUPPORTED_PROFILES
        and not _is_safe_institutional_profile_component(component)
    ]
    if invalid:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PROFILE,
            message="Unsupported execution profile component.",
            fix=(
                f"Choose comma-separated values from: {', '.join(sorted(SUPPORTED_PROFILES))}, "
                "or use safe institutional profile names containing only letters, numbers, '.', '_', and '-'."
            ),
            details={"profile": profile, "invalid_components": invalid},
        )
    checked = [_check_profile(component) for component in components]
    runtime_checks = [item for item in checked if item.get("backend_path")]
    return {
        "profile": profile,
        "components": components,
        "backend_path": runtime_checks[-1]["backend_path"] if runtime_checks else None,
        "backend_ready": all(
            bool(item.get("backend_ready", False)) for item in checked
        ),
        "component_checks": checked,
    }


def _is_safe_institutional_profile_component(component: str) -> bool:
    return bool(component and _PROFILE_COMPONENT_RE.fullmatch(component))


def _check_institutional_profile_component(
    profile: str, *, full_profile: str
) -> dict[str, object]:
    if _is_safe_institutional_profile_component(profile):
        return {
            "profile": profile,
            "backend_path": None,
            "backend_ready": True,
            "institutional_profile": True,
        }
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PROFILE,
        message="Unsupported execution profile.",
        fix=(
            f"Choose one of: {', '.join(sorted(SUPPORTED_PROFILES))}, "
            "or use a safe institutional profile name containing only letters, numbers, '.', '_', and '-'."
        ),
        details={"profile": full_profile, "invalid_components": [profile]},
    )


def _check_docker_profile(profile: str) -> dict[str, object]:
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
        info = subprocess.run(
            ["docker", "info"],
            capture_output=True,
            text=True,
            timeout=_SUBPROCESS_TIMEOUT,
        )
        docker_ok = info.returncode == 0
    except (subprocess.TimeoutExpired, OSError):
        docker_ok = False
    if not docker_ok:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.DOCKER_NOT_RUNNING,
            message="Docker is installed but its service is not available.",
            fix="Start Docker Desktop or the Docker daemon before running this skill.",
            details={"profile": profile},
        )
    return {"profile": profile, "backend_path": docker_path, "backend_ready": True}


def _check_conda_profile(profile: str) -> dict[str, object]:
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


def _check_singularity_compatible_profile(profile: str) -> dict[str, object]:
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
                "Install Singularity or Apptainer and ensure it is available on PATH, "
                "or choose a different profile (docker, conda)."
            ),
            details={"profile": profile, "tried": [primary, fallback]},
        )
    return {"profile": profile, "backend_path": backend, "backend_ready": True}


def _check_mamba_profile(profile: str) -> dict[str, object]:
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


def _check_podman_profile(profile: str) -> dict[str, object]:
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
        info = subprocess.run(
            ["podman", "info"],
            capture_output=True,
            text=True,
            timeout=_SUBPROCESS_TIMEOUT,
        )
        podman_ok = info.returncode == 0
    except (subprocess.TimeoutExpired, OSError):
        podman_ok = False
    if not podman_ok:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.PODMAN_NOT_RUNNING,
            message="Podman is installed but its service is not available.",
            fix="Start the Podman service or socket before running this skill.",
            details={"profile": profile},
        )
    return {"profile": profile, "backend_path": podman_path, "backend_ready": True}


def _check_hpc_profile(profile: str) -> dict[str, object]:
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


# Root-level entries that must not trip OUTPUT_DIR_NOT_EMPTY. `reproducibility` is the
# wrapper's own bundle directory (commands.sh, checksums.sha256, environment.yml,
# remap_paths.py, resource_limits.config, params/samplesheet snapshots, provenance JSON)
# — never user data — so its entire contents are ignored. Matches nfcore-sarek /
# nfcore-rnaseq so the three wrappers accept the same set of pre-existing output-dir
# contents, and avoids a fragile per-file allowlist that silently rejected a copied
# bundle (OUTPUT_DIR_NOT_EMPTY). Genuine result artifacts (report.md, upstream/,
# provenance/, result.json at the output root) still block a non-resume re-run.
_IGNORED_ROOT_NAMES = frozenset(
    {"reproducibility", ".DS_Store", ".gitkeep", ".gitignore", "Thumbs.db", "check_result.json"}
)


def _is_relative_to(path: Path, parent: Path) -> bool:
    try:
        path.relative_to(parent)
    except ValueError:
        return False
    return True


def _check_output_dir(output_dir: Path, *, resume: bool) -> None:
    # Refuse to write pipeline outputs inside the ClawBio source tree — large
    # upstream artifacts would pollute the repository. Parity with nfcore-rnaseq /
    # nfcore-sarek. Checked before mkdir so no directory is created on rejection.
    if _is_relative_to(output_dir.expanduser().resolve(), PROJECT_ROOT.resolve()):
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
    materialized_entries = []
    for entry in output_dir.iterdir():
        if entry.name in _IGNORED_ROOT_NAMES:
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
_SYMBOLIC_REFERENCE_FIELDS = frozenset(SYMBOLIC_REFERENCE_FIELDS)

# Every local path collected and existence-checked (genome refs + auxiliary
# files). Categories are centralised in schemas.py so preflight and
# params_builder cannot diverge (audit H-3).
_EXPLICIT_REFERENCE_FIELDS = ALL_REFERENCE_PATH_FIELDS


def _check_remote_inputs(args, samplesheet_summary: dict[str, Any]) -> None:
    """Local-first by default: reject remote input/reference URIs (s3://, gs://,
    https://, ftp://, …) unless ``--allow-remote-inputs`` is set, in which case
    warn that genetic data will be fetched over the network. The object-store
    ``--work-dir`` is a separate, intentionally-remote setting and is not gated."""
    remote: list[str] = []
    for key in ("fastq_paths", "bam_paths", "cram_paths"):
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


def _check_references(args) -> dict[str, str]:
    resolved = _collect_reference_values(args)
    _reject_conflicting_reference_styles(resolved)
    satisfied_group = _find_satisfied_reference_group(args.preset, resolved)
    _check_preset_specific_requirements(args.preset, resolved, satisfied_group)
    _check_reference_paths_exist(resolved)
    # Record only references the user actually supplied. The full dict (with empty
    # slots for every unused field) is needed by the checks above, but persisting
    # 20+ empty-string entries into inputs.json/preflight.json would be misleading
    # provenance and noise — provenance must reflect what was actually used.
    return {field: value for field, value in resolved.items() if value}


def _collect_reference_values(args) -> dict[str, str]:
    genome = getattr(args, "genome", None) or ""
    cellrangerarc_reference = getattr(args, "cellrangerarc_reference", None) or ""
    resolved = {"genome": genome, "cellrangerarc_reference": cellrangerarc_reference}
    for field in _EXPLICIT_REFERENCE_FIELDS:
        resolved[field] = getattr(args, field, None) or ""
    return resolved


def _reject_conflicting_reference_styles(resolved: dict[str, str]) -> None:
    # Only GENOME references conflict with --genome. Auxiliary files (barcode
    # whitelists, CMO/probe/feature sets, ...) are compatible with a --genome
    # shortcut and must not be flagged as conflicts (audit H-2).
    if resolved["genome"] and any(resolved[f] for f in GENOME_REFERENCE_FIELDS):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.CONFLICTING_REFERENCES,
            message="--genome (iGenomes shortcut) and explicit reference paths are mutually exclusive.",
            fix="Use either --genome <shortcut> or explicit --fasta/--gtf/--index flags, not both.",
            details={"genome": resolved["genome"]},
        )


def _find_satisfied_reference_group(
    preset: str, resolved: dict[str, str]
) -> tuple[str, ...]:
    requirements = PRESET_REQUIREMENTS[preset]["requires_any"]
    for option_group in requirements:
        if all(resolved.get(name, "") for name in option_group):
            return option_group
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.MISSING_REFERENCE,
        message="The selected preset is missing required references or indexes.",
        fix="Provide one of the supported reference/index combinations for this preset.",
        details={
            "preset": preset,
            "accepted_combinations": [list(group) for group in requirements],
        },
    )


def _check_preset_specific_requirements(
    preset: str,
    resolved: dict[str, str],
    satisfied_group: tuple[str, ...],
) -> None:
    if preset == "cellrangerarc" and satisfied_group != ("cellranger_index",):
        _check_cellrangerarc_config_pairing(resolved)


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
        details={
            "preset": "cellrangerarc",
            "missing_field": missing_field,
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _check_reference_paths_exist(resolved: dict[str, str]) -> None:
    for key, value in resolved.items():
        if key in _SYMBOLIC_REFERENCE_FIELDS:
            continue  # symbolic identifiers, not local paths
        if value and "://" in value:
            # Remote URI (s3://, gs://, https://, …): Nextflow resolves it at runtime,
            # so existence is not pre-checked here. Matches nfcore-rnaseq/sarek and the
            # nf-core schema, which types references as strings Nextflow may fetch.
            continue
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
        details={
            "field": key,
            "path": value,
            "schema_pattern": r"^\S+\.fn?a(sta)?(\.gz)?$",
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _check_parameter_schema_compatibility(args) -> None:
    """Validate non-path parameters whose upstream schema would fail predictably."""
    _check_email_schema_compatibility(getattr(args, "email", None), field="email")
    _check_email_schema_compatibility(
        getattr(args, "email_on_fail", None), field="email_on_fail"
    )
    _check_extra_config_paths(args)
    _warn_if_trusted_nextflow_configs(args)
    _check_config_param_overrides(args)
    _check_work_dir(args)


def _check_extra_config_paths(args) -> None:
    for raw_path in getattr(args, "extra_config", []) or []:
        config_path = Path(str(raw_path)).expanduser()
        if config_path.exists() and config_path.is_file():
            continue
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="A Nextflow config file passed with -c/--config was not found or is not a file.",
            fix="Provide an existing Nextflow config file path, or remove the -c/--config argument.",
            details={"field": "extra_config", "path": str(raw_path)},
        )


def _warn_if_trusted_nextflow_configs(args) -> None:
    configs = getattr(args, "extra_config", []) or []
    if not configs:
        return
    print(
        "WARNING: -c/--config files are trusted Nextflow/Groovy code. "
        "Only pass configs you authored or trust; the wrapper validates paths but does not sandbox contents.",
        file=sys.stderr,
    )


# A `-c` config that assigns params OUTSIDE the audited params.yaml overrides
# pipeline parameters — nf-core explicitly advises against setting params via -c.
# Detect every realistic spelling so such a config cannot silently change
# aligner/outdir/input/etc. (audit finding #2, hardened in audit H-04):
#   * dotted assignment:          params.aligner = 'x'  /  params.genomes.GRCh38.fasta = ...
#   * bracket (index) assignment: params['aligner'] = 'x'  /  params["a"]['b'] = ...
#   * whole-map reassignment:     params = [aligner: 'x']
#   * a params { ... } block, including the `{` placed on the next line.
# ``=(?!=)`` matches ASSIGNMENTS only — not a read such as ``if (params.x == 'y')``
# or ``x = params.foo`` (those do not start the line with ``params``). The segment
# group is ``*`` (not ``+``) so the bare ``params = ...`` whole-map reassignment is
# caught as well as the dotted/bracketed member assignments. The contents of a `-c`
# config are still trusted Groovy (not sandboxed); this lint only enforces that the
# audited params.yaml stays the single parameter source unless --trust-config-params.
_PARAMS_SEGMENT = r"""(?:\.\s*\w+|\[\s*['"][^'"]+['"]\s*\])"""
_PARAMS_ASSIGN_RE = re.compile(r"^params\s*" + _PARAMS_SEGMENT + r"*\s*=(?!=)")
_PARAMS_BLOCK_RE = re.compile(r"^params\s*\{")


def _scan_config_for_param_overrides(config_path: str) -> list[str]:
    try:
        text = Path(config_path).expanduser().read_text(encoding="utf-8")
    except OSError:
        return []  # existence/readability already enforced by _check_extra_config_paths
    name = Path(config_path).name
    hits: list[str] = []
    prev_code_was_bare_params = False
    for line_number, raw_line in enumerate(text.splitlines(), start=1):
        line = raw_line.strip()
        if line.startswith("//") or line.startswith("*") or line.startswith("/*"):
            continue
        # Strip a trailing line comment so `process.cpus = 4 // params.x note` is safe.
        code = line.split("//", 1)[0].strip()
        if not code:
            continue
        # `params` block whose opening brace sits on the following line.
        if prev_code_was_bare_params and code.startswith("{"):
            hits.append(f"{name}:{line_number}: {line}")
        if _PARAMS_ASSIGN_RE.search(code) or _PARAMS_BLOCK_RE.match(code):
            hits.append(f"{name}:{line_number}: {line}")
        prev_code_was_bare_params = code == "params"
    return hits


def _check_config_param_overrides(args) -> None:
    overrides: list[str] = []
    for config_path in getattr(args, "extra_config", []) or []:
        overrides.extend(_scan_config_for_param_overrides(str(config_path)))
    if not overrides:
        return
    if getattr(args, "trust_config_params", False):
        # Trusted: allow, warn, and stash for provenance so the bundle records that
        # params.yaml was not the only parameter source for this run.
        args.config_param_overrides = overrides
        print(
            "WARNING: -c/--config sets params.* outside the audited params.yaml "
            f"({len(overrides)} line(s)); allowed by --trust-config-params and recorded in provenance.",
            file=sys.stderr,
        )
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message=(
            "A -c/--config file sets pipeline parameters (params.*), which would override "
            "the wrapper's audited params.yaml. nf-core advises against setting params via -c."
        ),
        fix=(
            "Remove the params.* assignments (use -c only for infrastructure/process tuning), "
            "or pass --trust-config-params to allow them (recorded in provenance)."
        ),
        details={"overrides": overrides},
    )


def _is_remote_uri(value: str) -> bool:
    return bool(re.match(r"^[A-Za-z][A-Za-z0-9+.-]*://", value))


def _check_work_dir(args) -> None:
    raw_work_dir = getattr(args, "work_dir", None)
    if not raw_work_dir:
        return
    work_dir = str(raw_work_dir).strip()
    if not work_dir:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--work-dir cannot be empty.",
            fix="Provide a local directory path or a remote object-store URI such as s3://bucket/work.",
            details={"field": "work_dir"},
        )
    if _is_remote_uri(work_dir):
        return
    path = Path(work_dir).expanduser()
    if path.exists() and not path.is_dir():
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="--work-dir exists but is not a directory.",
            fix="Choose a directory path for --work-dir, or remove the flag to use <output>/upstream/work.",
            details={"field": "work_dir", "path": work_dir},
        )


def _check_email_schema_compatibility(email: str | None, *, field: str) -> None:
    if not email or _EMAIL_SCHEMA_RE.match(email):
        return
    flag = "--" + field.replace("_", "-")
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message=f"{field} does not match the nf-core/scrnaseq 4.1.0 schema.",
        fix=f"Use a simple email address such as user@example.org, or omit {flag}.",
        details={
            "field": field,
            "value": email,
            "schema_pattern": _EMAIL_SCHEMA_RE.pattern,
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _normalize_protocol_token(protocol: str) -> str:
    return re.sub(r"[\s_-]+", "", protocol.lower())


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


def run_preflight(
    args,
    *,
    pipeline_source: dict[str, object],
    samplesheet_summary: dict[str, Any],
) -> dict[str, object]:
    _check_demo_network(args)
    _warn_if_native_windows()
    _check_remote_inputs(args, samplesheet_summary)
    _check_supported_preset(args.preset)
    _check_parameter_schema_compatibility(args)
    _check_protocol_compatibility(args)
    _check_rna_velocity_pairings(args)
    java_info = _check_java()
    nextflow_info = _check_nextflow()
    profile_info = _check_profile(args.profile)
    output_dir = Path(args.output).expanduser().resolve()
    # Belt-and-suspenders: main() already called check_output_dir_available() before
    # samplesheet normalization. This second call keeps run_preflight() self-consistent
    # when invoked directly (e.g., tests). The reproducibility/ bundle (including the
    # samplesheet.valid.csv written between the two calls) is ignored wholesale by the
    # check, so it never triggers OUTPUT_DIR_NOT_EMPTY.
    check_output_dir_available(output_dir, resume=args.resume)
    refs = {} if args.demo else _check_references(args)
    if not args.demo:
        _check_igenomes_base(args)
        _check_cellranger_runtime_policy(args)
        _warn_if_seq_center_non_star(args)
    _check_samplesheet_driven_preset_requirements(
        args, samplesheet_summary=samplesheet_summary
    )
    _check_resume_compatibility(
        args, output_dir=output_dir, pipeline_source=pipeline_source
    )
    _warn_if_macos_docker_tmp(args.profile, output_dir)
    _warn_if_capped_remote_run(args)

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


def _check_samplesheet_driven_preset_requirements(
    args, *, samplesheet_summary: dict[str, Any]
) -> None:
    if args.demo or args.preset != "cellrangermulti":
        return
    feature_types = {
        str(value).lower() for value in samplesheet_summary.get("feature_types", [])
    }
    _check_cellrangermulti_mutually_exclusive_multiplexing(args, feature_types)
    _check_cellrangermulti_barcode_modes(args)
    _check_cellrangermulti_feature_references(args, feature_types)
    _check_cellrangermulti_vdj_reference_policy(args, feature_types)
    _check_cellrangermulti_multiplexing_policy(args, feature_types)


# Feature-barcode capture types that share the single Cell Ranger Multi feature
# reference (--fb_reference). In nf-core/scrnaseq 4.1.0 both Antibody Capture and
# CRISPR Guide Capture are wired through the same ``params.fb_reference`` channel —
# there is no separate CRISPR reference param — so a feature_type=crispr run without
# --fb-reference has no feature reference and fails inside Cell Ranger. Fail fast in
# preflight for both (audit H-03).
_FEATURE_BARCODE_REFERENCE_TYPES = ("ab", "crispr")


def _check_cellrangermulti_feature_references(args, feature_types: set[str]) -> None:
    needs_fb_reference = sorted(
        ft for ft in _FEATURE_BARCODE_REFERENCE_TYPES if ft in feature_types
    )
    if needs_fb_reference and not getattr(args, "fb_reference", None):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti antibody/CRISPR feature-barcode rows require a feature-barcode reference.",
            fix=f"Provide --fb-reference for feature_type {' and '.join(needs_fb_reference)} rows.",
            details={
                "preset": args.preset,
                "missing_field": "fb_reference",
                "feature_type": needs_fb_reference[0],
                "feature_types": needs_fb_reference,
                "policy_source": POLICY_SOURCE_CLAWBIO,
            },
        )


# nf-core/scrnaseq barcodes-samplesheet columns that select a multiplexing mode.
# The mode is encoded per multiplexed sample by exactly one of these columns
# (usage.md "with CMOs/FFPE/OCMs" examples); they are mutually exclusive.
_MULTIPLEXING_MODE_COLUMNS = {
    "probe_barcode_ids": "FFPE",
    "cmo_ids": "CMO",
    "ocm_ids": "OCM",
}


def _check_cellrangermulti_mutually_exclusive_multiplexing(
    args, feature_types: set[str]
) -> None:
    # Coarse run-level signal: FFPE is requested via --gex-frna-probe-set and CMO
    # via feature_type=cmo (or the --gex-cmo-set override). These are mutually
    # exclusive Cell Ranger Multi modes, so reject the flag-level combination early.
    # (OCM has no run-level flag — it is encoded only in the barcodes samplesheet,
    # validated authoritatively by _check_cellrangermulti_barcode_modes below.)
    has_ffpe_probe_set = bool(getattr(args, "gex_frna_probe_set", None))
    has_cmo = "cmo" in feature_types or bool(getattr(args, "gex_cmo_set", None))
    if not has_ffpe_probe_set or not has_cmo:
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message="cellrangermulti FFPE and CMO multiplexing modes are mutually exclusive in nf-core/scrnaseq.",
        fix="Use either --gex-frna-probe-set for FFPE or CMO inputs, not both, and run separate analyses if needed.",
        details={
            "preset": args.preset,
            "conflicting_fields": ["gex_frna_probe_set", "cmo"],
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _check_cellrangermulti_barcode_modes(args) -> None:
    """Reject a cellrangermulti barcodes samplesheet that mixes multiplexing modes.

    The second samplesheet (``--cellranger-multi-barcodes``) encodes the
    multiplexing mode per multiplexed sample via the ``probe_barcode_ids`` (FFPE),
    ``cmo_ids`` (CMO) and ``ocm_ids`` (OCM) columns. nf-core/scrnaseq usage docs:
    "FFPE; CMO and OCM are mutually exclusive ... Using more than one for a single
    sample will cause the module to fail." This validates that authoritative source
    per physical ``sample`` so the violation fails fast in preflight rather than
    mid-run (audit F-2). It is the correct OCM signal — ``gex_barcode_sample_assignment``
    is a separate cell/tag-calling override, NOT a mode selector. Best-effort: an
    unreadable/unparseable CSV is deferred to upstream nf-schema validation.
    """
    barcodes_path = getattr(args, "cellranger_multi_barcodes", None)
    if not barcodes_path:
        return
    try:
        with Path(barcodes_path).expanduser().open(
            "r", encoding="utf-8-sig", newline=""
        ) as handle:
            rows = list(csv.DictReader(handle))
    except (OSError, csv.Error, UnicodeDecodeError):
        return  # malformed/unreadable here → leave it to upstream nf-schema
    modes_by_sample: dict[str, set[str]] = {}
    for row in rows:
        sample = str(row.get("sample", "") or "").strip()
        for column, mode in _MULTIPLEXING_MODE_COLUMNS.items():
            if str(row.get(column, "") or "").strip():
                modes_by_sample.setdefault(sample, set()).add(mode)
    for sample in sorted(modes_by_sample):
        modes = modes_by_sample[sample]
        if len(modes) > 1:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message=(
                    "cellrangermulti FFPE, CMO and OCM multiplexing modes are "
                    "mutually exclusive per sample in nf-core/scrnaseq."
                ),
                fix=(
                    "In the --cellranger-multi-barcodes samplesheet, populate only "
                    "one of probe_barcode_ids (FFPE), cmo_ids (CMO) or ocm_ids (OCM) "
                    "for each sample; split mixed samples into separate runs."
                ),
                details={
                    "preset": args.preset,
                    "sample": sample,
                    "conflicting_modes": sorted(modes),
                    "policy_source": POLICY_SOURCE_NFCORE_DOCS,
                },
            )


def _check_cellrangermulti_vdj_reference_policy(args, feature_types: set[str]) -> None:
    if "vdj" not in feature_types:
        return
    has_vdj_index = bool(getattr(args, "cellranger_vdj_index", None))
    if getattr(args, "skip_cellrangermulti_vdjref", False):
        # mkvdjref is skipped → a prebuilt VDJ reference is the only valid source.
        if has_vdj_index:
            return
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti VDJ rows need a VDJ reference unless mkvdjref is allowed to run.",
            fix="Remove --skip-cellrangermulti-vdjref, or provide --cellranger-vdj-index for feature_type=vdj rows.",
            details={
                "preset": getattr(args, "preset", ""),
                "feature_type": "vdj",
                "missing_field": "cellranger_vdj_index",
                "conflicting_flag": "skip_cellrangermulti_vdjref",
                "policy_source": POLICY_SOURCE_CLAWBIO,
            },
        )
    # mkvdjref will run: the 4.1.0 usage docs say the VDJ reference is "calculated
    # on the fly given the reference files (--fasta and --gtf)". Those references can
    # come from explicit --fasta/--gtf OR from --genome (iGenomes resolves to
    # fasta+gtf). A GEX-only prebuilt --cellranger-index, however, cannot build a VDJ
    # reference, so require buildable references or a prebuilt --cellranger-vdj-index
    # (audit finding #3).
    if has_vdj_index:
        return
    has_buildable_refs = bool(getattr(args, "genome", None)) or (
        bool(getattr(args, "fasta", None)) and bool(getattr(args, "gtf", None))
    )
    if has_buildable_refs:
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message=(
            "cellrangermulti VDJ rows need a VDJ reference: mkvdjref builds it from "
            "genome references (--genome or --fasta/--gtf), which were not supplied."
        ),
        fix=(
            "Provide --genome or --fasta/--gtf so mkvdjref can build the VDJ reference, or "
            "pass a prebuilt --cellranger-vdj-index (a GEX --cellranger-index alone is not enough)."
        ),
        details={
            "preset": getattr(args, "preset", ""),
            "feature_type": "vdj",
            "missing_any_of": ["genome", "fasta+gtf", "cellranger_vdj_index"],
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _check_cellrangermulti_multiplexing_policy(args, feature_types: set[str]) -> None:
    has_barcodes = bool(getattr(args, "cellranger_multi_barcodes", None))
    # nf-core/scrnaseq usage docs: "When working with multiplexed data (FFPE/CMO/OCM),
    # you'll need a second samplesheet relating the multiplexed samples to the
    # physical sample", passed via --cellranger_multi_barcodes. CMO is detectable
    # from the samplesheet (feature_type=cmo); FFPE is signalled by --gex-frna-probe-set.
    if "cmo" in feature_types and not has_barcodes:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti CMO multiplexing requires a barcode-to-sample samplesheet.",
            fix="Provide --cellranger-multi-barcodes (the second samplesheet with cmo_ids) for feature_type=cmo rows.",
            details={
                "preset": args.preset,
                "missing_field": "cellranger_multi_barcodes",
                "feature_type": "cmo",
                "policy_source": POLICY_SOURCE_NFCORE_DOCS,
            },
        )
    has_ffpe_probe_set = bool(getattr(args, "gex_frna_probe_set", None))
    if has_ffpe_probe_set and not has_barcodes:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="cellrangermulti FFPE probe-set demultiplexing requires a barcode-to-sample samplesheet.",
            fix="Provide --cellranger-multi-barcodes when using --gex-frna-probe-set.",
            details={
                "preset": args.preset,
                "missing_field": "cellranger_multi_barcodes",
                "field": "gex_frna_probe_set",
                "policy_source": POLICY_SOURCE_NFCORE_DOCS,
            },
        )


def _check_protocol_compatibility(args) -> None:
    if getattr(args, "demo", False):
        return
    raw_protocol = getattr(args, "protocol", None)
    preset = getattr(args, "preset", "")
    if raw_protocol is None or not str(raw_protocol).strip():
        if preset in PRESETS_REQUIRING_EXPLICIT_PROTOCOL:
            raise SkillError(
                stage="preflight",
                error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
                message="The selected preset requires an explicit protocol in nf-core/scrnaseq 4.1.0.",
                fix="Provide --protocol with a value supported by the selected aligner, such as 10XV2, 10XV3, 10XV4, dropseq, or smartseq where documented. nf-core/scrnaseq takes the protocol from this flag (a global pipeline parameter); a `protocol` column in the samplesheet is not read for it, so pass the value explicitly here even if your samplesheet has one.",
                details={
                    "preset": preset,
                    "missing_field": "protocol",
                    "policy_source": POLICY_SOURCE_NFCORE_DOCS,
                },
            )
        return
    protocol = str(raw_protocol).strip()
    normalized_protocol = _normalize_protocol_token(protocol)
    # `auto` chemistry detection is defined by protocols.json only for the
    # cellranger family; cellrangermulti is samplesheet-driven (not in the matrix)
    # so its protocol is deferred to upstream rather than rejected here.
    if normalized_protocol == "auto" and preset in PRESETS_REQUIRING_EXPLICIT_PROTOCOL:
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="The selected preset does not support protocol=auto in nf-core/scrnaseq 4.1.0.",
            fix="Use an explicit protocol such as 10XV2, 10XV3, dropseq, smartseq, or a custom string supported by the selected aligner.",
            details={
                "preset": preset,
                "protocol": protocol,
                "policy_source": POLICY_SOURCE_NFCORE_DOCS,
            },
        )
    if normalized_protocol == "smartseq2":
        # nf-core/scrnaseq 4.1.0 docs: 'smartseq' denotes Smart-seq3; Smart-seq2 is
        # explicitly unsupported and routed to nf-core/rnaseq (usage.md). Reject it
        # for every preset rather than passing an unrecognised token to an aligner.
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="Smart-seq2 is not supported by nf-core/scrnaseq 4.1.0.",
            fix=(
                "Process Smart-seq2 data with nf-core/rnaseq instead. In nf-core/scrnaseq, "
                "'smartseq' denotes Smart-seq3 (STARsolo and Kallisto/bustools only)."
            ),
            details={
                "preset": preset,
                "protocol": protocol,
                "policy_source": POLICY_SOURCE_NFCORE_DOCS,
            },
        )
    if (
        normalized_protocol == "smartseq"
        and preset not in PRESETS_SUPPORTING_SMARTSEQ_PROTOCOL
    ):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="Smart-seq (Smart-seq3) is supported only by STARsolo and Kallisto/bustools in nf-core/scrnaseq 4.1.0.",
            fix="Use --preset star or --preset kallisto for protocol=smartseq, or choose a protocol supported by the selected preset.",
            details={
                "preset": preset,
                "protocol": protocol,
                "policy_source": POLICY_SOURCE_NFCORE_DOCS,
            },
        )
    _check_protocol_against_nfcore_matrix(
        preset=preset,
        protocol=protocol,
        normalized_protocol=normalized_protocol,
    )


def _check_rna_velocity_pairings(args) -> None:
    if (
        getattr(args, "preset", "") == "star"
        and getattr(args, "star_feature", None) == "Gene Velocyto"
        and not getattr(args, "star_ignore_sjdbgtf", False)
    ):
        raise SkillError(
            stage="preflight",
            error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
            message="STARsolo RNA velocity requires --star-ignore-sjdbgtf with --star-feature 'Gene Velocyto'.",
            fix="Add --star-ignore-sjdbgtf, or use a non-velocity --star-feature value.",
            details={
                "preset": "star",
                "field": "star_feature",
                "value": "Gene Velocyto",
                "missing_field": "star_ignore_sjdbgtf",
                "policy_source": POLICY_SOURCE_NFCORE_DOCS,
            },
        )
    if getattr(args, "preset", "") != "kallisto":
        return
    if getattr(args, "kb_workflow", None) not in {"lamanno", "nac"}:
        return
    has_t1c = bool(getattr(args, "kb_t1c", None))
    has_t2c = bool(getattr(args, "kb_t2c", None))
    # When building the index from --fasta/--gtf (or --genome), `kb ref` generates
    # the cDNA/intron capture files, so neither is required — this matches the
    # documented `--kb-workflow nac --fasta --gtf` example, and the 4.1.0 parameter
    # docs never mark kb_t1c/kb_t2c as required. They are required only with a
    # prebuilt --kallisto-index (not regenerated). A partial spec (exactly one of
    # the two) is always a mistake regardless of how the index is obtained.
    uses_prebuilt_index = bool(getattr(args, "kallisto_index", None))
    if has_t1c and has_t2c:
        return
    if not uses_prebuilt_index and not has_t1c and not has_t2c:
        return
    missing = [
        name for name in ("kb_t1c", "kb_t2c") if not getattr(args, name, None)
    ]
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message=(
            "Kallisto RNA velocity needs both --kb-t1c and --kb-t2c here."
        ),
        fix=(
            "Provide both capture files (required with a prebuilt --kallisto-index; "
            "generated automatically when building from --fasta/--gtf), or use "
            "--kb-workflow standard."
        ),
        details={
            "preset": "kallisto",
            "kb_workflow": getattr(args, "kb_workflow", None),
            "missing_fields": missing,
            "uses_prebuilt_index": uses_prebuilt_index,
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _check_protocol_against_nfcore_matrix(
    *, preset: str, protocol: str, normalized_protocol: str
) -> None:
    supported_tokens = set(PROTOCOLS_JSON_4_1_0.get(preset, ()))
    if not supported_tokens:
        # cellrangermulti is samplesheet-driven and absent from protocols.json.
        return
    known_but_unsupported = (
        normalized_protocol in KNOWN_PROTOCOL_TOKENS
        and normalized_protocol not in supported_tokens
    )
    custom_not_supported = (
        normalized_protocol not in KNOWN_PROTOCOL_TOKENS
        and preset not in PRESETS_SUPPORTING_CUSTOM_PROTOCOL
    )
    if not known_but_unsupported and not custom_not_supported:
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message="The selected protocol is not supported by the selected preset in nf-core/scrnaseq 4.1.0.",
        fix=(
            "Choose a protocol documented for the selected preset. Cell Ranger supports auto and 10XV1-10XV4; "
            "Cell Ranger ARC supports auto only."
        ),
        details={
            "preset": preset,
            "protocol": protocol,
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _warn_if_macos_docker_tmp(profile: str, output_dir: Path) -> None:
    # Colima (a common macOS Docker runtime) only 9p-mounts the user HOME directory
    # into its VM.  /tmp and /private/tmp live on the VM's own ext4 and are NOT shared
    # with the host, so Nextflow's work-dir files are invisible to containers.
    # Docker Desktop handles /tmp differently and is unlikely to hit this, but the safe
    # guidance for macOS + any Docker backend is to keep output dirs under HOME.
    if sys.platform != "darwin" or not profile_includes(profile, "docker"):
        return
    if is_under_tmp(output_dir):
        print(
            "WARNING: Output directory is under /tmp. On macOS with Colima, Docker containers "
            "cannot see files written to /tmp (the VM uses its own separate /tmp). "
            "Move --output to a path under your home directory to avoid 'No such file or directory' errors.",
            file=sys.stderr,
        )


def _warn_if_capped_remote_run(args) -> None:
    """Hint to disable the wall-clock cap on HPC/cloud runs the scheduler bounds.

    The wrapper applies a default 12 h cap (``--timeout-hours``) that kills the
    Nextflow process tree on expiry. On HPC/cloud — where the batch scheduler
    already enforces walltime — a long but legitimate run can hit that cap. Surface
    a hint to pass ``--timeout-hours 0`` when the run targets an object-store work
    directory or an institutional/site profile (the two clear non-local signals),
    while the cap is still active (audit H-05). Advisory only; never blocks.
    """
    timeout_hours = getattr(args, "timeout_hours", None)
    if not timeout_hours:  # None or 0 → the cap is already disabled
        return
    work_dir = str(getattr(args, "work_dir", "") or "").strip()
    remote_work_dir = bool(work_dir) and _is_remote_uri(work_dir)
    components = profile_components(getattr(args, "profile", "") or "")
    institutional_profile = any(
        component not in SUPPORTED_PROFILES for component in components
    )
    if not remote_work_dir and not institutional_profile:
        return
    target = (
        "an object-store work directory"
        if remote_work_dir
        else "an institutional/site profile"
    )
    print(
        f"WARNING: a {timeout_hours:g} h Nextflow wall-clock cap is active "
        f"(--timeout-hours). This run targets {target}; on HPC/cloud the scheduler "
        "already enforces walltime, so pass --timeout-hours 0 to avoid the wrapper "
        "killing a long but legitimate run.",
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


def _check_igenomes_base(args) -> None:
    """Existence-check a LOCAL ``--igenomes-base`` mirror when iGenomes is used.

    Remote bases (s3://, https://, ...) are resolved by Nextflow at runtime and
    are left untouched. A local mirror, however, is exactly the air-gapped use
    case the docs promote (``--genome GRCh38 --igenomes-base /mnt/local_igenomes``),
    so a typo'd path should fail fast in preflight like every other reference
    path instead of surfacing as a late Nextflow error (audit F-08). iGenomes is
    only consumed when a ``--genome`` shortcut is given, so the check is scoped to
    that case to avoid false positives.
    """
    base = getattr(args, "igenomes_base", None)
    if not base or not getattr(args, "genome", None):
        return
    if "://" in str(base):
        return  # remote base resolved by Nextflow (s3://, https://, gs://, ...)
    base_path = Path(str(base)).expanduser()
    if base_path.exists():
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.MISSING_REFERENCE,
        message="The local --igenomes-base mirror was not found.",
        fix=(
            "Point --igenomes-base at an existing local iGenomes mirror, use the "
            "default S3 base by omitting the flag, or correct the path."
        ),
        details={"field": "igenomes_base", "path": str(base)},
    )


def _warn_if_conda_cellranger(args) -> None:
    """Warn when a Cell Ranger preset is combined with a conda/mamba profile.

    Cell Ranger is not distributed via bioconda/biocontainers (10x Genomics
    licensing); nf-core ships it only inside the docker/singularity containers it
    builds. A conda/mamba profile therefore cannot resolve Cell Ranger and the run
    would fail mid-pipeline with a confusing resolver error. Surface it early as a
    warning (not a hard error, since site configs may inject a custom path).
    """
    if getattr(args, "preset", "") not in CELLRANGER_FAMILY_PRESETS:
        return
    components = profile_components(getattr(args, "profile", "") or "")
    if any(component in {"conda", "mamba"} for component in components):
        print(
            f"WARNING: preset {args.preset!r} uses Cell Ranger, which is not "
            "distributed via bioconda/conda (10x Genomics licensing). The "
            f"{args.profile!r} profile can run it only when a trusted site config "
            "injects a valid Cell Ranger installation. Otherwise use -profile docker "
            "or singularity for Cell Ranger.",
            file=sys.stderr,
        )


def _check_cellranger_runtime_policy(args) -> None:
    if getattr(args, "preset", "") not in CELLRANGER_FAMILY_PRESETS:
        return
    components = profile_components(getattr(args, "profile", "") or "")
    uses_conda = any(component in {"conda", "mamba"} for component in components)
    if not uses_conda:
        return
    if getattr(args, "allow_conda_cellranger", False):
        _warn_if_conda_cellranger(args)
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_PRESET_CONFIGURATION,
        message="Cell Ranger presets are not supported with conda/mamba profiles by default.",
        fix=(
            "Use -profile docker or singularity for Cell Ranger, or pass "
            "--allow-conda-cellranger only when a trusted site config injects a valid Cell Ranger installation."
        ),
        details={
            "preset": getattr(args, "preset", ""),
            "profile": getattr(args, "profile", ""),
            "policy_source": POLICY_SOURCE_NFCORE_DOCS,
        },
    )


def _warn_if_seq_center_non_star(args) -> None:
    """Warn when --seq-center is supplied for a non-STARsolo aligner.

    nf-core/scrnaseq 4.1.0 documents seq_center as affecting STARsolo BAM read
    groups only. It is still forwarded verbatim (a harmless no-op for other
    aligners, which ignore it), but the user should know it will have no effect
    outside STARsolo rather than assume it propagated (audit F-6).
    """
    if not getattr(args, "seq_center", None):
        return
    if getattr(args, "preset", "") == "star":
        return
    print(
        f"WARNING: --seq-center is set but preset {getattr(args, 'preset', '')!r} "
        "is not STARsolo. In nf-core/scrnaseq 4.1.0 seq_center only tags STARsolo "
        "BAM read groups; it is forwarded but will have no effect for this aligner.",
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
    pipeline_source: dict[str, object],
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
    _check_resume_work_dir(args, output_dir=output_dir, manifest=manifest)


def _check_resume_preset_and_profile(args, manifest: dict[str, object]) -> None:
    if (
        manifest.get("preset") == args.preset
        and manifest.get("profile") == args.profile
    ):
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
    pipeline_source: dict[str, object],
    manifest: dict[str, object],
) -> None:
    previous_source = manifest.get("pipeline_source", {})
    if not isinstance(previous_source, dict):
        previous_source = {}
    if previous_source.get("source_kind") == pipeline_source.get(
        "source_kind"
    ) and previous_source.get("resolved_version") == pipeline_source.get(
        "resolved_version"
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


def _check_resume_work_dir(args, *, output_dir: Path, manifest: dict[str, object]) -> None:
    previous_work_dir = str(manifest.get("work_dir") or "upstream/work")
    requested_work_dir = _effective_work_dir_for_resume(args, output_dir=output_dir)
    if previous_work_dir == requested_work_dir:
        return
    raise SkillError(
        stage="preflight",
        error_code=ErrorCode.INVALID_RESUME_STATE,
        message="Resume state does not match the requested Nextflow work directory.",
        fix="Resume only with the same --work-dir as the original run, or start in a new output directory.",
        details={
            "previous_work_dir": previous_work_dir,
            "requested_work_dir": requested_work_dir,
        },
    )


def _effective_work_dir_for_resume(args, *, output_dir: Path) -> str:
    raw_work_dir = getattr(args, "work_dir", None)
    if not raw_work_dir:
        return "upstream/work"
    work_dir = str(raw_work_dir).strip()
    if _is_remote_uri(work_dir):
        return work_dir
    return Path(work_dir).expanduser().resolve().as_posix()
