from __future__ import annotations

import json
import platform
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

from clawbio.common.checksums import sha256_file
from clawbio.common.reproducibility import write_checksums, write_environment_yml

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(_SKILL_DIR))

from schemas import DEFAULT_REMOTE_PIPELINE, JAVA_MIN_VERSION, NEXTFLOW_MIN_VERSION, SKILL_ALIAS, SKILL_NAME, SKILL_VERSION


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
    provenance_dir = output_dir / "provenance"
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
        "cwd": str(Path.cwd()),
        "work_dir": str(output_dir / "upstream" / "work"),
        "command": command_str,
        "java_version": preflight_result["java"]["version"],
        "nextflow_version": preflight_result["nextflow"]["version"],
    }


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
    }


def build_inputs_payload(
    *,
    normalized_samplesheet: Path,
    samplesheet_summary: dict[str, Any],
    preflight_result: dict[str, Any],
    params_path: Path,
) -> dict[str, Any]:
    return {
        "samplesheet": str(normalized_samplesheet),
        "samplesheet_checksum": sha256_file(normalized_samplesheet),
        "sample_count": samplesheet_summary["sample_count"],
        "fastq_paths": [Path(p).as_posix() for p in samplesheet_summary["fastq_paths"]],
        "reference_paths": preflight_result.get("references", {}),
        "params_path": str(params_path),
        "params_checksum": sha256_file(params_path),
    }


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
    for filename, payload in payloads.items():
        (provenance_dir / filename).write_text(json.dumps(payload, indent=2), encoding="utf-8")


def write_reproducibility_environment(output_dir: Path, *, preflight_result: dict[str, Any]) -> None:
    write_environment_yml(
        output_dir,
        env_name="clawbio-nfcore-scrnaseq-wrapper",
        pip_deps=[],
        conda_deps=[
            f"openjdk>={JAVA_MIN_VERSION}",
            f"nextflow>={'.'.join(map(str, NEXTFLOW_MIN_VERSION))}",
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
    checksum_paths: list[Path] = [
        normalized_samplesheet,
        params_path,
        Path(execution_result["stdout_path"]),
        Path(execution_result["stderr_path"]),
    ]
    for ref_path in preflight_result.get("references", {}).values():
        if ref_path:
            p = Path(str(ref_path))
            if p.is_file():
                checksum_paths.append(p)
    for candidate in parsed_outputs.get("h5ad_candidates", []):
        checksum_paths.append(Path(candidate))
    if parsed_outputs.get("multiqc_report"):
        checksum_paths.append(Path(str(parsed_outputs["multiqc_report"])))
    checksum_paths = list(dict.fromkeys(Path(p) for p in checksum_paths if p))
    return write_checksums(checksum_paths, output_dir, anchor=output_dir)


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
        "params_checksum": inputs["params_checksum"],
        "samplesheet_checksum": inputs["samplesheet_checksum"],
        "checksums_file": str(checksum_file),
    }
    manifest_path = output_dir / "reproducibility" / "manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2), encoding="utf-8")
    return manifest_path
