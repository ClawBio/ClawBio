from __future__ import annotations

from clawbio.common.nfcore.artifacts import write_json_lf, write_json_snapshots

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from clawbio.common.checksums import sha256_file
from clawbio.common.reproducibility import write_checksums, write_environment_yml

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))
sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("schemas")

from clawbio.common.textio import write_text_lf
from schemas import DEFAULT_REMOTE_PIPELINE, SKILL_ALIAS, SKILL_NAME, SKILL_VERSION


def _relativise(path: Path | str, anchor: Path) -> str:
    """Return ``path`` as POSIX relative to ``anchor`` when inside it, else absolute POSIX."""
    p = Path(str(path))
    try:
        return p.resolve().relative_to(Path(anchor).resolve()).as_posix()
    except ValueError:
        return p.as_posix()


def _relativise_command(command_str: str, output_dir: Path) -> str:
    """Replace the absolute output dir (the ``cd`` target) with '.' in a recorded
    command, so runtime.json does not leak the generation environment's path.

    Both the as-given and resolved POSIX forms are replaced (longest first) to
    cover whichever the live invocation embedded.
    """
    candidates = {Path(output_dir).as_posix(), Path(output_dir).resolve().as_posix()}
    for anchor in sorted(candidates, key=len, reverse=True):
        command_str = command_str.replace(anchor, ".")
    return command_str


def write_provenance_bundle(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, Any],
    preflight_result: dict[str, Any],
    params_path: Path,
    params_payload: dict[str, Any],
    normalized_samplesheet: Path,
    samplesheet_summary: dict[str, Any],
    parsed_outputs: dict[str, Any],
    execution_result: dict[str, Any],
    command_str: str,
) -> tuple[Path, Path]:
    provenance_dir = output_dir / "reproducibility"
    provenance_dir.mkdir(parents=True, exist_ok=True)
    timestamp = datetime.now(timezone.utc).isoformat()
    payloads = build_provenance_payloads(
        output_dir,
        args=args,
        pipeline_source=pipeline_source,
        preflight_result=preflight_result,
        params_path=params_path,
        params_payload=params_payload,
        normalized_samplesheet=normalized_samplesheet,
        samplesheet_summary=samplesheet_summary,
        parsed_outputs=parsed_outputs,
        command_str=command_str,
        timestamp=timestamp,
    )
    _write_provenance_payloads(provenance_dir, payloads)
    _copy_policy_files(provenance_dir)
    write_reproducibility_environment(output_dir, preflight_result=preflight_result)
    checksum_file = write_reproducibility_checksums(
        output_dir,
        normalized_samplesheet=normalized_samplesheet,
        params_path=params_path,
        preflight_result=preflight_result,
        parsed_outputs=parsed_outputs,
        execution_result=execution_result,
    )
    manifest_path = write_reproducibility_manifest(
        output_dir,
        args=args,
        upstream=payloads["upstream.json"],
        inputs=payloads["inputs.json"],
        runtime=payloads["runtime.json"],
        checksum_file=checksum_file,
    )
    return provenance_dir, manifest_path


def build_provenance_payloads(
    output_dir: Path,
    *,
    args,
    pipeline_source: dict[str, Any],
    preflight_result: dict[str, Any],
    params_path: Path,
    params_payload: dict[str, Any],
    normalized_samplesheet: Path,
    samplesheet_summary: dict[str, Any],
    parsed_outputs: dict[str, Any],
    command_str: str,
    timestamp: str,
) -> dict[str, dict[str, Any]]:
    return {
        "runtime.json": build_runtime_payload(
            output_dir,
            args=args,
            preflight_result=preflight_result,
            command_str=command_str,
            timestamp=timestamp,
        ),
        "upstream.json": build_upstream_payload(pipeline_source),
        "invocation.json": build_invocation_payload(args, timestamp=timestamp),
        "inputs.json": build_inputs_payload(
            normalized_samplesheet=normalized_samplesheet,
            samplesheet_summary=samplesheet_summary,
            preflight_result=preflight_result,
            params_path=params_path,
            output_dir=output_dir,
        ),
        "outputs.json": build_outputs_payload(parsed_outputs),
        "skill.json": build_skill_payload(params_payload),
        "preflight.json": preflight_result,
    }


def build_runtime_payload(
    output_dir: Path,
    *,
    args,
    preflight_result: dict[str, Any],
    command_str: str,
    timestamp: str,
) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "os": platform.system(),
        "arch": platform.machine(),
        "python_version": platform.python_version(),
        "profile": args.profile,
        "resume_used": bool(args.resume),
        "cwd": ".",
        "work_dir": _runtime_work_dir(args, output_dir),
        "command": _relativise_command(command_str, output_dir),
        "java_version": preflight_result["java"]["version"],
        "nextflow_version": preflight_result["nextflow"]["version"],
    }


def _runtime_work_dir(args, output_dir: Path) -> str:
    raw_work_dir = getattr(args, "work_dir", None)
    if not raw_work_dir:
        return _relativise(output_dir / "upstream" / "work", output_dir)
    work_dir = str(raw_work_dir).strip()
    if "://" in work_dir:
        return work_dir
    return Path(work_dir).expanduser().resolve().as_posix()


def build_upstream_payload(pipeline_source: dict[str, Any]) -> dict[str, Any]:
    return {
        "pipeline": DEFAULT_REMOTE_PIPELINE,
        "source_kind": pipeline_source["source_kind"],
        "source_ref": pipeline_source["source_ref"],
        "resolved_version": pipeline_source["resolved_version"],
        "branch": pipeline_source.get("branch", ""),
        "dirty": pipeline_source.get("dirty", False),
    }


def build_invocation_payload(args, *, timestamp: str) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "preset": args.preset,
        "demo": bool(args.demo),
        "check_only": bool(args.check),
        "profile": args.profile,
        "pipeline_version": args.pipeline_version,
        "allow_dirty_pipeline": bool(getattr(args, "allow_dirty_pipeline", False)),
        "require_local_pipeline": bool(getattr(args, "require_local_pipeline", False)),
        "allow_conda_cellranger": bool(getattr(args, "allow_conda_cellranger", False)),
        "allow_pipeline_version_override": bool(
            getattr(args, "allow_pipeline_version_override", False)
        ),
        "trust_config_params": bool(getattr(args, "trust_config_params", False)),
        "config_param_overrides": list(getattr(args, "config_param_overrides", []) or []),
        "work_dir": getattr(args, "work_dir", None) or "",
        "extra_config": [
            str(path) for path in (getattr(args, "extra_config", []) or [])
        ],
    }


def build_inputs_payload(
    *,
    normalized_samplesheet: Path,
    samplesheet_summary: dict[str, Any],
    preflight_result: dict[str, Any],
    params_path: Path,
    output_dir: Path,
) -> dict[str, Any]:
    return {
        "samplesheet": _relativise(normalized_samplesheet, output_dir),
        "samplesheet_checksum": sha256_file(normalized_samplesheet),
        "sample_count": samplesheet_summary["sample_count"],
        # fastq_paths and reference_paths are deliberately kept absolute: they point
        # to external data outside the bundle and are remapped by remap_paths.py on
        # replay, not relativised against the output dir. Remote URIs (s3://, https://,
        # …) are preserved verbatim — Path().as_posix() would collapse the `//` scheme
        # separator and corrupt them.
        "fastq_paths": [
            str(p) if "://" in str(p) else Path(p).as_posix()
            for p in samplesheet_summary["fastq_paths"]
        ],
        "reference_paths": preflight_result.get("references", {}),
        # Reference *file* digests are recorded here (provenance) rather than in
        # checksums.sha256, because references live outside the bundle: putting
        # bare-basename entries in checksums.sha256 would make `sha256sum -c` fail.
        # Directory indexes (star/simpleaf/cellranger) have no file digest, so they
        # are skipped here and tracked only by path.
        "reference_checksums": _reference_file_checksums(
            preflight_result.get("references", {})
        ),
        "params_path": _relativise(params_path, output_dir),
        "params_checksum": sha256_file(params_path),
    }


def _reference_file_checksums(references: dict[str, Any]) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for field, path in (references or {}).items():
        if not path:
            continue
        p = Path(str(path))
        if p.is_file():
            checksums[field] = sha256_file(p)
    return checksums


def build_outputs_payload(parsed_outputs: dict[str, Any]) -> dict[str, Any]:
    return {
        "preferred_h5ad": parsed_outputs.get("preferred_h5ad", ""),
        "multiqc_report": parsed_outputs.get("multiqc_report", ""),
        "pipeline_info_dir": parsed_outputs.get("pipeline_info_dir", ""),
        "h5ad_candidates": parsed_outputs.get("h5ad_candidates", []),
        "rds_candidates": parsed_outputs.get("rds_candidates", []),
        "handoff_available": parsed_outputs.get("handoff_available", False),
    }


def build_skill_payload(params_payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": SKILL_NAME,
        "cli_alias": SKILL_ALIAS,
        "version": SKILL_VERSION,
        "params": params_payload,
    }


def _write_provenance_payloads(provenance_dir: Path, payloads: dict[str, Any]) -> None:
    write_json_snapshots(provenance_dir, payloads)


def _copy_policy_files(bundle_dir: Path) -> None:
    """Copy compatibility_policy.json + pinned_versions.json from the skill into the bundle."""
    import shutil

    src_dir = _SKILL_DIR / "reproducibility"
    for name in ("compatibility_policy.json", "pinned_versions.json"):
        src = src_dir / name
        if src.is_file():
            shutil.copyfile(src, bundle_dir / name)


def write_reproducibility_environment(
    output_dir: Path, *, preflight_result: dict[str, Any]
) -> None:
    java_version = preflight_result["java"]["version"]
    nextflow_version = preflight_result["nextflow"]["version"]
    # Make environment.yml an installable recipe (manifest.json keeps the exact
    # versions for provenance):
    #  - nextflow is published on bioconda, not conda-forge → add that channel
    #    and pin the exact engine version (also enforced at replay via NXF_VER).
    #  - openjdk is pinned to its major series: conda-forge does not publish every
    #    JDK patch, so an exact patch (e.g. 17.0.10) can be unsatisfiable, and the
    #    Java patch level does not affect containerised nf-core task results.
    java_major = str(java_version).split(".")[0]
    write_environment_yml(
        output_dir,
        env_name="clawbio-nfcore-scrnaseq-wrapper",
        channels=["conda-forge", "bioconda"],
        pip_deps=[],
        conda_deps=[
            f"openjdk={java_major}",
            f"nextflow={nextflow_version}",
        ],
        python_version=f"{platform.python_version_tuple()[0]}.{platform.python_version_tuple()[1]}",
    )


def write_reproducibility_checksums(
    output_dir: Path,
    *,
    normalized_samplesheet: Path,
    params_path: Path,
    preflight_result: dict[str, Any],
    parsed_outputs: dict[str, Any],
    execution_result: dict[str, Any],
) -> Path:
    # checksums.sha256 is a self-verifiable manifest: every entry must resolve
    # relative to output_dir so `sha256sum -c checksums.sha256` succeeds when run
    # from the output directory on any OS. Only in-bundle artifacts are listed;
    # external reference digests live in inputs.json (see _reference_file_checksums)
    # because bare-basename entries for out-of-tree files would break `-c`.
    candidate_paths: list[Path] = [
        normalized_samplesheet,
        params_path,
        Path(execution_result["stdout_path"]),
        Path(execution_result["stderr_path"]),
    ]
    for candidate in parsed_outputs.get("h5ad_candidates", []):
        candidate_paths.append(Path(candidate))
    if parsed_outputs.get("multiqc_report"):
        candidate_paths.append(Path(str(parsed_outputs["multiqc_report"])))

    resolved_output_dir = output_dir.resolve()
    # write_checksums accepts list[Path | str]; declare the same element type so
    # the invariant list matches its parameter without a cast.
    checksum_paths: list[Path | str] = []
    for p in candidate_paths:
        if not p:
            continue
        p = Path(p)
        if not _is_within(p, resolved_output_dir):
            # Defensive: anything outside the bundle would get a bare-basename
            # label and break `sha256sum -c`. In-tree by construction today.
            continue
        checksum_paths.append(p)
    checksum_paths = list(dict.fromkeys(checksum_paths))
    return write_checksums(checksum_paths, output_dir, anchor=output_dir)


def _is_within(path: Path, anchor: Path) -> bool:
    try:
        path.resolve().relative_to(anchor)
        return True
    except ValueError:
        return False


def write_reproducibility_manifest(
    output_dir: Path,
    *,
    args,
    upstream: dict[str, Any],
    inputs: dict[str, Any],
    runtime: dict[str, Any],
    checksum_file: Path,
) -> Path:
    manifest = {
        "skill_name": SKILL_NAME,
        "skill_version": SKILL_VERSION,
        "pipeline_source": upstream,
        "preset": args.preset,
        "profile": args.profile,
        "resume_used": bool(args.resume),
        "generated_at": runtime["timestamp"],
        "java_version": runtime["java_version"],
        "nextflow_version": runtime["nextflow_version"],
        "python_version": runtime["python_version"],
        "environment_yml_mode": "install_recipe",
        "work_dir": runtime["work_dir"],
        "params_checksum": inputs["params_checksum"],
        "samplesheet_checksum": inputs["samplesheet_checksum"],
        "checksums_file": _relativise(checksum_file, output_dir),
    }
    manifest_path = output_dir / "reproducibility" / "manifest.json"
    write_json_lf(manifest_path, manifest, trailing_newline=False)
    return manifest_path
