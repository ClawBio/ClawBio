"""Container images a finished Nextflow run used, resolved to immutable digests.

A tag can be re-pushed; a digest cannot. The record lists, per image, the exact
reference Nextflow ran, its digest and where the digest came from, and every
process that used it. When every image is pinned, a Nextflow config pins each
process to ``repository@digest``. Registry-only observations of mutable tags do
not establish which image ran and never produce replay pins.
"""

from __future__ import annotations

import json
import os
import subprocess
from collections.abc import Callable, Mapping, Sequence
from pathlib import Path

from clawbio.common.textio import write_text_lf

from .artifacts import write_json_lf
from .runtime import nxf_offline

SCHEMA_VERSION = 1
RECORD_NAME = "containers.json"
PIN_CONFIG_NAME = "containers.config"
_TASK_LOG = ["nextflow", "log", "last", "-f", "process,container,status"]
_FINISHED = {"COMPLETED", "CACHED"}
# Engines whose images have registry digests the wrappers can pin.
_DIGEST_ENGINES = {"docker", "podman"}
_ENGINES = ("docker", "podman", "singularity", "apptainer", "charliecloud", "shifter")

Runner = Callable[..., "subprocess.CompletedProcess[str] | None"]


def _run(
    command: Sequence[str], *, cwd: Path | str | None = None, env: Mapping[str, str] | None = None,
    timeout_seconds: float,
) -> subprocess.CompletedProcess[str] | None:
    try:
        return subprocess.run(
            list(command), cwd=cwd, env={**os.environ, **env} if env else None, capture_output=True,
            text=True, errors="replace", timeout=timeout_seconds, check=False,
        )
    except (subprocess.TimeoutExpired, OSError):
        return None


def container_engine(profile: str) -> str | None:
    """The container engine a Nextflow profile selects, or None (conda, mamba, none)."""
    parts = [part.strip() for part in str(profile or "").split(",")]
    return next((engine for engine in _ENGINES if engine in parts), None)


def _repository(image: str) -> str:
    if "@" in image:
        return image.split("@", 1)[0]
    name, _, tag = image.rpartition(":")
    return name if name and "/" not in tag else image


def _task_images(output_dir, env, runner) -> tuple[dict[str, set[str]], int, list[str]]:
    proc = runner(_TASK_LOG, cwd=output_dir, env=env, timeout_seconds=180)
    if proc is None or proc.returncode != 0:
        return {}, 0, ["`nextflow log` could not list this run's tasks, so its containers are unrecorded."]
    images: dict[str, set[str]] = {}
    tasks = 0
    for line in proc.stdout.splitlines():
        fields = line.split("\t")
        if len(fields) != 3 or fields[2].strip() not in _FINISHED:
            continue
        tasks += 1
        process, image = fields[0].strip(), fields[1].strip()
        if image:
            images.setdefault(image, set()).add(process)
    return images, tasks, []


def _engine_platform(engine, runner) -> tuple[str | None, str | None]:
    proc = runner([engine, "version", "--format", "{{.Server.Version}} {{.Server.Os}}/{{.Server.Arch}}"],
                  timeout_seconds=30)
    if proc is None or proc.returncode != 0 or not proc.stdout.split():
        return None, None
    fields = proc.stdout.split()
    return fields[0], fields[1] if len(fields) > 1 else None


def _normalized_repository(repository: str) -> str:
    """Match Docker Hub's explicit registry/default-library aliases locally."""
    for prefix in ("docker.io/", "index.docker.io/"):
        if repository.startswith(prefix):
            repository = repository[len(prefix):]
            break
    return repository if "/" in repository else f"library/{repository}"


def _local_image(engine, image, runner) -> dict | None:
    proc = runner([engine, "image", "inspect", "--format", "{{json .}}", image], timeout_seconds=60)
    if proc is None or proc.returncode != 0:
        return None
    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None
    repository = _normalized_repository(_repository(image))
    digest = next(
        (entry.split("@", 1)[1] for entry in data.get("RepoDigests") or [] if _normalized_repository(entry.split("@", 1)[0]) == repository),
        None,
    )
    if digest is None:
        return None
    platform = "/".join(part for part in (data.get("Os"), data.get("Architecture"), data.get("Variant")) if part)
    return {"digest": digest, "digest_source": "local", "image_id": data.get("Id"), "platform": platform or None}


class _RegistryUnreachable(Exception):
    """The registry query timed out or could not start; later images would too."""


def _registry_image(image, platform, runner) -> dict | None:
    # quay.io intermittently answers a valid manifest with an error (seen as
    # `unsupported manifest format`), so a refused query is asked once more.
    for _attempt in range(2):
        proc = runner(["docker", "manifest", "inspect", "-v", image], timeout_seconds=120)
        if proc is None:
            raise _RegistryUnreachable(image)
        if proc.returncode == 0:
            break
    else:
        return None
    try:
        data = json.loads(proc.stdout)
    except json.JSONDecodeError:
        return None
    entries = data if isinstance(data, list) else [data]

    def platform_of(entry):
        found = (entry.get("Descriptor") or {}).get("platform") or {}
        return "/".join(part for part in (found.get("os"), found.get("architecture"), found.get("variant")) if part)

    chosen = next((entry for entry in entries if platform and platform_of(entry) == platform), None)
    if chosen is None and len(entries) == 1:
        chosen = entries[0]
    digest = ((chosen or {}).get("Descriptor") or {}).get("digest")
    if not digest:
        return None
    return {"digest": digest, "digest_source": "registry", "image_id": None, "platform": platform_of(chosen) or None}


def collect_container_record(
    output_dir, *, profile: str, env: Mapping[str, str] | None = None, runner: Runner | None = None,
) -> dict:
    """Images the run in ``output_dir`` used, per ``nextflow log``, with digests.

    Call it right after the run: digests are read from the local image first,
    then from the registry (another job may have removed the image meanwhile).
    Registry-only observations cannot establish historical task identity and do
    not qualify the record for immutable replay pins.
    """
    runner = runner or _run
    engine = container_engine(profile)
    if not (Path(output_dir) / ".nextflow" / "history").is_file():
        return {
            "schema_version": SCHEMA_VERSION, "source": "nextflow log", "engine": engine,
            "engine_version": None, "engine_platform": None, "task_count": 0, "images": [],
            "complete": False,
            "warnings": ["No Nextflow run history in this output directory, so its containers are unrecorded."],
        }
    images, tasks, warnings = _task_images(output_dir, env, runner)
    version, platform = _engine_platform(engine, runner) if engine else (None, None)
    # Only an image with no local digest is looked up in its registry. Offline, or
    # once one query times out, the others are skipped instead of waiting in turn.
    offline = nxf_offline({**os.environ, **(env or {})})
    registry = engine == "docker" and not offline
    entries = []
    for image in sorted(images):
        entry = {"image": image, "digest": None, "pinned_reference": None, "digest_source": None,
                 "image_id": None, "platform": None, "processes": sorted(images[image])}
        if "@sha256:" in image:
            entry.update(digest=image.split("@", 1)[1], pinned_reference=image, digest_source="reference")
        elif engine in _DIGEST_ENGINES:
            resolved = _local_image(engine, image, runner)
            if resolved is None and offline:
                warnings.append(f"NXF_OFFLINE is set; no registry lookup was made for {image} without a local digest.")
            if resolved is None and registry:
                try:
                    resolved = _registry_image(image, platform, runner)
                except _RegistryUnreachable:
                    registry = False
                    warnings.append(
                        f"The registry was unreachable for {image}, so no further images were looked up there."
                    )
            if resolved:
                entry.update(resolved)
                if resolved["digest_source"] == "registry":
                    warnings.append(
                        f"The registry lookup for {image} observes the current tag only; the original "
                        "task image and platform cannot be confirmed, so replay pins are withheld."
                    )
                else:
                    entry["pinned_reference"] = f"{_repository(image)}@{resolved['digest']}"
        if entry["digest"] is None and engine is not None:
            warnings.append(f"No digest recorded for {image}; a replay may run a re-pushed tag.")
        entries.append(entry)
    owners: dict[str, list[str]] = {}
    for image, processes in images.items():
        for process in processes:
            owners.setdefault(process, []).append(image)
    for process in sorted(name for name, used in owners.items() if len(used) > 1):
        warnings.append(f"{process} ran several images ({', '.join(sorted(owners[process]))}); it cannot be pinned by name.")
    return {
        "schema_version": SCHEMA_VERSION,
        "source": "nextflow log",
        "engine": engine,
        "engine_version": version,
        "engine_platform": platform,
        "task_count": tasks,
        "images": entries,
        "complete": not warnings,
        "warnings": warnings,
    }


def pin_config(record: dict) -> str | None:
    """Nextflow config pinning every process to the digest it ran, if all are known."""
    if record.get("engine") not in _DIGEST_ENGINES or not record.get("complete") or not record.get("images"):
        return None
    pins = sorted(
        (process, image["pinned_reference"]) for image in record["images"] for process in image["processes"]
    )
    lines = [
        "// Generated by ClawBio from this run's containers.json: every process runs the",
        "// exact image (repository@digest) it ran originally, even if a tag is re-pushed.",
        "process {",
    ]
    for process, reference in pins:
        lines += [f"    withName: '{process}' {{", f"        container = '{reference}'", "    }"]
    return "\n".join(lines + ["}", ""])


def write_container_provenance(bundle_dir: Path, record: dict) -> list[Path]:
    """Write containers.json and, when every image is pinned, containers.config."""
    bundle_dir = Path(bundle_dir)
    bundle_dir.mkdir(parents=True, exist_ok=True)
    record_path = bundle_dir / RECORD_NAME
    write_json_lf(record_path, record, trailing_newline=True)
    written = [record_path]
    config = pin_config(record)
    pin_path = bundle_dir / PIN_CONFIG_NAME
    if config is None:
        pin_path.unlink(missing_ok=True)
    else:
        write_text_lf(pin_path, config)
        written.append(pin_path)
    return written


def record_run_containers(
    output_dir, *, profile: str, env: Mapping[str, str] | None = None, runner: Runner | None = None,
) -> dict:
    """Collect the run's container record and write it into its reproducibility bundle."""
    record = collect_container_record(output_dir, profile=profile, env=env, runner=runner)
    write_container_provenance(Path(output_dir) / "reproducibility", record)
    return record


def manifest_summary(record: dict) -> dict:
    """What manifest.json says about the run's containers."""
    pinned = pin_config(record) is not None
    return {
        "engine": record.get("engine"),
        "engine_version": record.get("engine_version"),
        "images": len(record.get("images", [])),
        "complete": bool(record.get("complete")),
        "record": f"reproducibility/{RECORD_NAME}",
        "pin_config": f"reproducibility/{PIN_CONFIG_NAME}" if pinned else None,
    }


def replay_pin_guard(*, resume_variable: str, bundle_dir: str = "$SCRIPT_DIR") -> str:
    """Bash lines setting $CONTAINER_PIN for a Nextflow-direct replay.

    On a fresh replay it is ``-c reproducibility/containers.config`` (the replay runs
    from the output directory), so every process runs the image digest it ran
    originally. An in-place resume keeps the cached tasks, which already ran those
    images, so it pins nothing. A plain string keeps macOS's bash 3.2 happy under
    ``set -u``.
    """
    return (
        "# Run every process on the exact image (repository@digest) it ran originally,\n"
        "# as recorded in containers.config, when replaying into a fresh directory.\n"
        'CONTAINER_PIN=""\n'
        f'if [[ -z "${resume_variable}" && -f "{bundle_dir}/{PIN_CONFIG_NAME}" ]]; then\n'
        f'  CONTAINER_PIN="-c reproducibility/{PIN_CONFIG_NAME}"\n'
        "fi\n"
    )
