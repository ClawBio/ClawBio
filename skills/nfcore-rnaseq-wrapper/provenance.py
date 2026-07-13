from __future__ import annotations

import json
import os
import platform
import socket
import sys
from datetime import datetime, timedelta, timezone
from pathlib import Path
from typing import Any

from clawbio.common.checksums import sha256_file
from clawbio.common.reproducibility import write_checksums, write_environment_yml
from clawbio.common.textio import write_text_lf

_SKILL_DIR = Path(__file__).resolve().parent
if str(_SKILL_DIR) in sys.path:
    sys.path.remove(str(_SKILL_DIR))
sys.path.insert(0, str(_SKILL_DIR))


sys.modules.pop("_isolated_imports", None)
from _isolated_imports import purge_foreign_bare_modules

purge_foreign_bare_modules("preflight", "schemas")

from preflight import params_payload_checksum
from schemas import DEFAULT_REMOTE_PIPELINE, JAVA_MIN_VERSION, NEXTFLOW_MIN_VERSION, PROJECT_ROOT, SKILL_ALIAS, SKILL_NAME, SKILL_VERSION


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
    duration_seconds: float = 0,
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
        duration_seconds=duration_seconds,
        timestamp=timestamp,
    )
    _write_provenance_payloads(provenance_dir, payloads)
    write_reproducibility_environment(output_dir, preflight_result=preflight_result)
    checksum_file = output_dir / "reproducibility" / "checksums.sha256"
    manifest_path = write_reproducibility_manifest(
        output_dir,
        args=args,
        upstream=payloads["upstream.json"],
        inputs=payloads["inputs.json"],
        runtime=payloads["runtime.json"],
        params_payload=params_payload,
        checksum_file=checksum_file,
    )
    checksum_file = write_reproducibility_checksums(output_dir, execution_result=execution_result)
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
    duration_seconds: float = 0,
) -> dict[str, dict[str, Any]]:
    return {
        "runtime.json": build_runtime_payload(
            output_dir,
            args=args,
            preflight_result=preflight_result,
            command_str=command_str,
            duration_seconds=duration_seconds,
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
    duration_seconds: float = 0,
) -> dict[str, Any]:
    _finished = datetime.fromisoformat(timestamp)
    _started = _finished - timedelta(seconds=duration_seconds)
    return {
        "timestamp": timestamp,
        "started_at": _started.isoformat(),
        "finished_at": _finished.isoformat(),
        "duration_seconds": duration_seconds,
        "os": platform.system(),
        "arch": platform.machine(),
        "host": socket.gethostname(),
        "python_version": platform.python_version(),
        "java_version": preflight_result.get("java", {}).get("version", ""),
        "nextflow_version": preflight_result.get("nextflow", {}).get("version", ""),
        "cpu_count": os.cpu_count() or 0,
        "ram_total_gb": _ram_total_gb(),
        "profile": getattr(args, "profile", ""),
        "resume_used": bool(getattr(args, "resume", False)),
        "cwd": str(Path.cwd()),
        "work_dir": str(output_dir / "upstream" / "work"),
        "command": command_str,
    }


def build_upstream_payload(pipeline_source: dict[str, Any]) -> dict[str, Any]:
    source_kind = pipeline_source.get("source_kind", "")
    source_ref = pipeline_source.get("source_ref", "")
    resolved_version = pipeline_source.get("resolved_version", "")
    is_local = source_kind == "local_checkout"
    payload: dict[str, Any] = {
        "pipeline": DEFAULT_REMOTE_PIPELINE,
        "version": resolved_version,
        "source_kind": source_kind,
        "local_path": source_ref if is_local else "",
        "remote_ref": "" if is_local else source_ref,
        "git_commit": pipeline_source.get("git_commit", pipeline_source.get("commit", resolved_version if is_local else "")),
        "git_branch": pipeline_source.get("git_branch", pipeline_source.get("branch", "")),
        "git_dirty": bool(pipeline_source.get("git_dirty", pipeline_source.get("dirty", False))),
        "source_ref": source_ref,
        "resolved_version": resolved_version,
    }
    # When a local checkout was found but rejected (e.g., whitespace in path),
    # preserve the diagnostic fields so upstream.json tells the full story.
    if pipeline_source.get("local_attempted"):
        payload["local_attempted"] = pipeline_source["local_attempted"]
        payload["local_rejected_reason"] = pipeline_source.get("local_rejected_reason", "")
    return payload


def build_invocation_payload(args, *, timestamp: str) -> dict[str, Any]:
    return {
        "timestamp": timestamp,
        "argv": list(sys.argv),
        "env_filtered": _filtered_env(),
        "working_dir": str(Path.cwd()),
        "clawbio_repo_root": str(PROJECT_ROOT),
        "aligner": getattr(args, "aligner", ""),
        "demo": bool(getattr(args, "demo", False)),
        "check_only": bool(getattr(args, "check", False)),
        "profile": getattr(args, "profile", ""),
        "pipeline_version": getattr(args, "pipeline_version", ""),
    }


_FASTQ_CHECKSUM_SIZE_LIMIT_GB = 50.0
_DIR_CHECKSUM_MAX_GB: float = 5.0


def _checksum_input_file(path_str: str) -> str:
    """SHA256 of a FASTQ/BAM file, or a sentinel if remote/missing/oversized."""
    if "://" in path_str:
        return "<remote-uri>"
    p = Path(path_str)
    if not p.exists():
        return "<missing>"
    size_gb = p.stat().st_size / (1024 ** 3)
    if size_gb > _FASTQ_CHECKSUM_SIZE_LIMIT_GB:
        return f"<skipped:{size_gb:.1f}gb>"
    return sha256_file(p)


def build_inputs_payload(
    *,
    normalized_samplesheet: Path,
    samplesheet_summary: dict[str, Any],
    preflight_result: dict[str, Any],
    params_path: Path,
) -> dict[str, Any]:
    reference_paths = preflight_result.get("references", {})
    fastq_paths = [_path_or_uri_as_posix(p) for p in samplesheet_summary.get("fastq_paths", [])]
    bam_paths = [_path_or_uri_as_posix(p) for p in samplesheet_summary.get("bam_paths", [])]
    return {
        "samplesheet_path": str(normalized_samplesheet),
        "samplesheet_checksum": sha256_file(normalized_samplesheet),
        "samples_count": samplesheet_summary.get("sample_count", 0),
        "fastq_paths": fastq_paths,
        "fastq_checksums": {p: _checksum_input_file(p) for p in fastq_paths},
        "bam_paths": bam_paths,
        "bam_checksums": {p: _checksum_input_file(p) for p in bam_paths},
        "reference_paths": reference_paths,
        "reference_checksums": _reference_checksums(reference_paths),
        "params_path": str(params_path),
        "params_checksum": sha256_file(params_path),
    }


def _path_or_uri_as_posix(value: object) -> str:
    text = str(value)
    if "://" in text:
        return text
    return Path(text).as_posix()


def build_outputs_payload(parsed_outputs: dict[str, Any]) -> dict[str, Any]:
    return dict(parsed_outputs)


def build_skill_payload(params_payload: dict[str, Any]) -> dict[str, Any]:
    return {
        "name": SKILL_NAME,
        "cli_alias": SKILL_ALIAS,
        "version": SKILL_VERSION,
        "params": params_payload,
    }


def write_reproducibility_environment(output_dir: Path, *, preflight_result: dict[str, Any]) -> None:
    # Use exact versions captured at runtime when available so the environment
    # file is a true reproduction spec, not just a minimum-bound recipe.
    java_version = preflight_result.get("java", {}).get("version", "")
    nextflow_version = preflight_result.get("nextflow", {}).get("version", "")
    java_dep = f"openjdk={java_version}" if java_version else f"openjdk>={JAVA_MIN_VERSION}"
    nf_dep = f"nextflow={nextflow_version}" if nextflow_version else f"nextflow>={'.'.join(map(str, NEXTFLOW_MIN_VERSION))}"
    write_environment_yml(
        output_dir,
        env_name="clawbio-nfcore-rnaseq-wrapper",
        pip_deps=[],
        conda_deps=[java_dep, nf_dep],
        python_version=f"{platform.python_version_tuple()[0]}.{platform.python_version_tuple()[1]}",
    )


def write_reproducibility_checksums(output_dir: Path, **_ignored: Any) -> Path:
    checksum_paths: list[Path] = []
    for root in (
        output_dir / "upstream" / "results",
        output_dir / "reproducibility",
        output_dir / "provenance",
        output_dir / "logs",
    ):
        if root.exists():
            checksum_paths.extend(path for path in sorted(root.rglob("*")) if path.is_file())
    checksum_paths = [p for p in checksum_paths if p.name != "checksums.sha256"]
    return write_checksums(list(dict.fromkeys(checksum_paths)), output_dir, anchor=output_dir)


def write_reproducibility_manifest(
    output_dir: Path,
    *,
    args,
    upstream: dict[str, Any],
    inputs: dict[str, Any],
    runtime: dict[str, Any],
    params_payload: dict[str, Any],
    checksum_file: Path,
) -> Path:
    manifest = {
        "skill_name": SKILL_NAME,
        "skill_version": SKILL_VERSION,
        "pipeline_source": upstream,
        "aligner": getattr(args, "aligner", ""),
        "pseudo_aligner": getattr(args, "pseudo_aligner", None),
        "profile": getattr(args, "profile", ""),
        "prokaryotic": bool(getattr(args, "prokaryotic", False)),
        # `demo` composes the upstream `test` profile, which supplies its own samplesheet
        # AND its own bundled references — resuming across the demo/real boundary would
        # swap both underneath the run, so it is part of the resume compatibility contract.
        "demo": bool(getattr(args, "demo", False)),
        # `arm` affects the generated .nextflow_macos_docker.config (platform flag
        # presence) — drift between runs would silently change container architecture.
        "arm": bool(getattr(args, "arm", False)),
        "resume_used": bool(getattr(args, "resume", False)),
        "generated_at": runtime["timestamp"],
        "java_version": runtime["java_version"],
        "nextflow_version": runtime["nextflow_version"],
        "python_version": runtime["python_version"],
        "environment_yml_mode": "exact_versions" if (runtime["java_version"] and runtime["nextflow_version"]) else "install_recipe",
        "params_checksum": params_payload_checksum(params_payload),
        "params_file_sha256": inputs["params_checksum"],
        "samplesheet_checksum": inputs["samplesheet_checksum"],
        "checksums_file": str(checksum_file),
    }
    manifest_path = output_dir / "reproducibility" / "manifest.json"
    manifest_path.parent.mkdir(parents=True, exist_ok=True)
    write_text_lf(manifest_path, json.dumps(manifest, indent=2))
    return manifest_path


def _write_provenance_payloads(provenance_dir: Path, payloads: dict[str, Any]) -> None:
    for filename, payload in payloads.items():
        write_text_lf(provenance_dir / filename, json.dumps(payload, indent=2))


def _dir_checksum(path: Path) -> str:
    """Stable SHA256 over sorted relative file paths + their individual hashes.

    Using relative paths makes the checksum location-independent: copying the
    directory to a different prefix produces the same hash if content is equal.

    Directories larger than _DIR_CHECKSUM_MAX_GB are skipped to avoid blocking
    on large reference indexes (e.g. STAR genome >=30 GB).
    """
    import hashlib as _hashlib
    total_bytes = sum(f.stat().st_size for f in path.rglob("*") if f.is_file())
    total_gb = total_bytes / (1024 ** 3)
    if total_gb > _DIR_CHECKSUM_MAX_GB:
        return f"<skipped:{total_gb:.1f}gb>"
    digest = _hashlib.sha256()
    for f in sorted(path.rglob("*")):
        if f.is_file():
            relative = f.relative_to(path).as_posix()
            digest.update(relative.encode())
            digest.update(sha256_file(f).encode())
    return digest.hexdigest()


def _reference_checksums(reference_paths: dict[str, Any]) -> dict[str, str]:
    checksums: dict[str, str] = {}
    for key, value in reference_paths.items():
        s = str(value)
        if not s:
            # Empty string means the user did not provide this reference field.
            # Path("") resolves to the CWD, so we must skip rather than hash it.
            continue
        if "://" in s:
            # Remote URI — not checksummable with a local hash; record the fact.
            checksums[key] = "<remote-uri>"
            continue
        if s.startswith("$"):
            # Bash variable template (e.g. $REFS/hg38.fa or ${REFS}/hg38.fa).
            # These are replay-time placeholders, not literal filesystem paths.
            checksums[key] = "<bash-variable>"
            continue
        path = Path(s)
        if path.is_file():
            checksums[key] = sha256_file(path)
        elif path.is_dir():
            checksums[key] = _dir_checksum(path)
        else:
            checksums[key] = "<missing>"
    return checksums


def _filtered_env() -> dict[str, str]:
    result: dict[str, str] = {}
    for key, value in os.environ.items():
        if key.startswith("NEXTFLOW_") or key == "JAVA_HOME":
            result[key] = value
        elif key == "PATH":
            result[key] = "<redacted>"
    return result


def _ram_total_gb() -> float:
    if hasattr(os, "sysconf"):
        try:
            pages = os.sysconf("SC_PHYS_PAGES")
            page_size = os.sysconf("SC_PAGE_SIZE")
            return round((pages * page_size) / (1024**3), 2)
        except (OSError, ValueError):
            return 0.0
    return 0.0
