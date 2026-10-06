"""Preflight: where the pipeline comes from, and whether this machine can run it.

Four questions are answered before anything is launched, because each one is
cheap now and expensive later: is the pipeline source a real checkout of the
pinned workflow, is snakemake new enough, is conda available to build tool
environments, and is the output directory safe to write into.

`network_requirements` exists so the downloads a run will make are stated out
loud before it starts, rather than discovered in a log an hour later.
"""

from __future__ import annotations

import os
import re
import shutil
import subprocess
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from magpipe_errors import ErrorCode, SkillError
from magpipe_schemas import (
    MIN_SNAKEMAKE,
    PINNED_REF,
    PIPELINE_REQUIRED_FILES,
    SNAKEMAKE_PIN,
    UPSTREAM_REPO,
)

_VERSION = re.compile(r"(\d+)\.(\d+)(?:\.(\d+))?")


def default_cache_dir() -> Path:
    """`$CLAWBIO_CACHE_DIR`, else `~/.cache/clawbio/pipelines`."""
    override = os.environ.get("CLAWBIO_CACHE_DIR")
    if override:
        return Path(override).expanduser()
    return Path.home() / ".cache" / "clawbio" / "pipelines"


def _fail(
    stage: str, code: str, message: str, fix: str, **details: object
) -> SkillError:
    return SkillError(
        stage=stage, error_code=code, message=message, fix=fix, details=dict(details)
    )


def _missing_files(root: Path) -> list[str]:
    return [name for name in PIPELINE_REQUIRED_FILES if not (root / name).is_file()]


def _git(root: Path, *args: str) -> str:
    """A git answer, or "" when git is absent or the directory is not a repo."""
    try:
        result = subprocess.run(
            ["git", *args],
            cwd=str(root),
            capture_output=True,
            text=True,
            timeout=30,
            check=False,
        )
    except (OSError, subprocess.SubprocessError):
        return ""
    return result.stdout.strip() if result.returncode == 0 else ""


def _describe(root: Path, source_kind: str, ref: str) -> dict:
    commit = _git(root, "rev-parse", "HEAD")
    status = _git(root, "status", "--porcelain")
    return {
        "path": Path(root),
        "source_kind": source_kind,
        "commit": commit or "unknown",
        "dirty": bool(status),
        "ref": ref,
    }


def resolve_pipeline_source(
    pipeline_dir: Path | None,
    *,
    cache_dir: Path,
    ref: str,
    allow_fetch: bool,
) -> dict:
    """Locate a checkout of the pinned workflow: given, cached, or cloned."""
    if pipeline_dir is not None:
        root = Path(pipeline_dir).expanduser()
        missing = _missing_files(root)
        if missing:
            raise _fail(
                "preflight",
                ErrorCode.PIPELINE_SOURCE_INVALID,
                f"{root} is not a checkout of the upstream workflow; missing: "
                + ", ".join(missing),
                "Pass --pipeline-dir pointing at a clone of "
                f"{UPSTREAM_REPO} checked out at {ref}.",
                path=str(root),
                missing=missing,
            )
        return _describe(root.resolve(), "local_checkout", ref)

    cache = Path(cache_dir).expanduser()
    target = cache / f"metagenomics-workflow-{ref[:12]}"
    if not _missing_files(target):
        return _describe(target.resolve(), "cached_clone", ref)

    if not allow_fetch:
        raise _fail(
            "preflight",
            ErrorCode.PIPELINE_FETCH_FAILED,
            f"no cached checkout of the upstream workflow at {target}",
            "Pass --pipeline-dir <checkout> to use one you already have, or allow "
            "the skill to clone the public workflow (no user data is fetched).",
            cache_dir=str(cache),
            ref=ref,
        )

    cache.mkdir(parents=True, exist_ok=True)
    try:
        subprocess.run(
            ["git", "clone", "--quiet", UPSTREAM_REPO, str(target)],
            capture_output=True,
            text=True,
            timeout=600,
            check=True,
        )
        subprocess.run(
            ["git", "checkout", "--quiet", ref],
            cwd=str(target),
            capture_output=True,
            text=True,
            timeout=300,
            check=True,
        )
    except subprocess.CalledProcessError as error:
        raise _fail(
            "preflight",
            ErrorCode.PIPELINE_FETCH_FAILED,
            f"could not clone {UPSTREAM_REPO} at {ref}: "
            f"{(error.stderr or '').strip() or error}",
            "Check outbound network access to github.com, or clone the workflow "
            "yourself and pass --pipeline-dir.",
            ref=ref,
            target=str(target),
        ) from error
    except subprocess.SubprocessError as error:
        raise _fail(
            "preflight",
            ErrorCode.PIPELINE_FETCH_FAILED,
            f"cloning {UPSTREAM_REPO} did not complete: {error}",
            "Clone the workflow yourself and pass --pipeline-dir.",
            ref=ref,
        ) from error

    missing = _missing_files(target)
    if missing:
        raise _fail(
            "preflight",
            ErrorCode.PIPELINE_SOURCE_INVALID,
            f"the clone at {target} is missing " + ", ".join(missing),
            f"The commit {ref} may not be the expected layout; pass --pipeline-dir "
            "with a checkout that has these files.",
            path=str(target),
            missing=missing,
        )
    return _describe(target.resolve(), "cached_clone", ref)


def check_snakemake() -> str:
    """Snakemake's version, or a SkillError naming how to get one."""
    if shutil.which("snakemake") is None:
        raise _fail(
            "preflight",
            ErrorCode.MISSING_SNAKEMAKE,
            "snakemake is not on PATH",
            f"conda create -n metagenomics -c conda-forge -c bioconda "
            f"snakemake={SNAKEMAKE_PIN} && conda activate metagenomics",
        )
    try:
        result = subprocess.run(
            ["snakemake", "--version"],
            capture_output=True,
            text=True,
            timeout=120,
            check=False,
        )
    except (OSError, subprocess.SubprocessError) as error:
        raise _fail(
            "preflight",
            ErrorCode.MISSING_SNAKEMAKE,
            f"could not run snakemake --version: {error}",
            f"Install it with `conda create -n metagenomics -c conda-forge "
            f"-c bioconda snakemake={SNAKEMAKE_PIN}`.",
        ) from error

    output = f"{result.stdout}\n{result.stderr}"
    match = _VERSION.search(output)
    if match is None:
        raise _fail(
            "preflight",
            ErrorCode.SNAKEMAKE_VERSION_TOO_OLD,
            f"could not read a version from `snakemake --version`: {output.strip()!r}",
            f"Install the pinned version: conda create -n metagenomics "
            f"-c conda-forge -c bioconda snakemake={SNAKEMAKE_PIN}",
        )
    major, minor = int(match.group(1)), int(match.group(2))
    if (major, minor) < MIN_SNAKEMAKE:
        found = ".".join(part for part in match.groups() if part)
        required = ".".join(str(part) for part in MIN_SNAKEMAKE)
        raise _fail(
            "preflight",
            ErrorCode.SNAKEMAKE_VERSION_TOO_OLD,
            f"snakemake {found} is too old; the upstream workflow requires "
            f"{required} or newer (min_version in workflow/Snakefile)",
            f"conda create -n metagenomics -c conda-forge -c bioconda "
            f"snakemake={SNAKEMAKE_PIN} && conda activate metagenomics",
            found=found,
            required=required,
        )
    return match.group(0)


def check_conda() -> str:
    """`conda` or `mamba`, whichever is on PATH."""
    for name in ("conda", "mamba"):
        if shutil.which(name) is not None:
            return name
    raise _fail(
        "preflight",
        ErrorCode.MISSING_CONDA,
        "neither conda nor mamba is on PATH",
        "Snakemake creates one conda environment per tool. Install Miniforge or "
        "Mambaforge and re-run.",
    )


def check_output_dir(output_dir: Path, *, force: bool) -> None:
    """Refuse to overwrite a populated directory, or one that cannot be written."""
    target = Path(output_dir).expanduser()
    if target.exists() and not target.is_dir():
        raise _fail(
            "preflight",
            ErrorCode.OUTPUT_DIR_NOT_WRITABLE,
            f"{target} exists and is not a directory",
            "Pass --output a directory path that does not exist yet, or move the "
            "file aside.",
            path=str(target),
        )
    if target.is_dir():
        entries = sorted(entry.name for entry in target.iterdir())
        if entries and not force:
            raise _fail(
                "preflight",
                ErrorCode.OUTPUT_DIR_NOT_EMPTY,
                f"{target} already holds {len(entries)} entr"
                f"{'y' if len(entries) == 1 else 'ies'}: " + ", ".join(entries[:5]),
                "Choose an empty --output directory, or pass --force to write into "
                "it anyway (existing files may be overwritten).",
                path=str(target),
                entries=entries,
            )
        probe = target
    else:
        probe = target.parent
        while not probe.exists() and probe != probe.parent:
            probe = probe.parent
    if not os.access(probe, os.W_OK | os.X_OK):
        raise _fail(
            "preflight",
            ErrorCode.OUTPUT_DIR_NOT_WRITABLE,
            f"{probe} is not writable by this user",
            "Choose an --output directory you own, or fix the permissions on "
            f"{probe}.",
            path=str(probe),
        )


def network_requirements(config: dict, samples: list) -> list[str]:
    """What this configuration will download, in the words of the upstream docs."""
    requirements = [
        (
            "Conda environments for the pipeline tools are created on first run "
            "from conda-forge and bioconda."
        )
    ]
    fetch = bool(config.get("fetch", {}).get("sra", {}).get("enabled", False))
    named = [sample for sample in samples if getattr(sample, "sra_run", None)]
    if fetch and named:
        accessions = ", ".join(sorted(sample.sra_run for sample in named))
        requirements.append(
            f"Public SRA/ENA/DRA runs are downloaded with sra-tools: {accessions}."
        )
    decontaminate = config.get("qc", {}).get("decontaminate", {})
    if decontaminate.get("enabled") and not decontaminate.get("fasta"):
        name = decontaminate.get("fasta_name", "GRCh38.primary_assembly.genome.fa.gz")
        requirements.append(
            f"The host reference {name} (~3 GB) is downloaded from the release "
            "pinned in qc.decontaminate.fasta_url. Pass --host-fasta to use a local "
            "copy instead."
        )
    profiler = config.get("profile", {}).get("profiler", "none")
    if profiler != "none":
        index = config.get("profile", {}).get(profiler, {}).get("index", "")
        requirements.append(
            f"The MetaPhlAn index {index} is installed from the configured release."
        )
    if config.get("annotation", {}).get("minpath", {}).get("data_dir"):
        requirements.append(
            "The MinPath script is fetched from a pinned commit and checked against "
            "annotation.minpath.script_md5."
        )
    return requirements


def run_preflight(
    *,
    output_dir: Path,
    config: dict,
    samples: list,
    pipeline_dir: Path | None = None,
    cache_dir: Path | None = None,
    ref: str = PINNED_REF,
    allow_fetch: bool = True,
    force: bool = False,
    check_tools: bool = True,
) -> dict:
    """Run every check and return what was found. Raises on the first failure."""
    check_output_dir(output_dir, force=force)
    source = resolve_pipeline_source(
        pipeline_dir,
        cache_dir=cache_dir or default_cache_dir(),
        ref=ref,
        allow_fetch=allow_fetch,
    )
    snakemake = check_snakemake() if check_tools else None
    conda = check_conda() if check_tools else None
    return {
        "output_dir": str(Path(output_dir)),
        "pipeline_source": source,
        "snakemake": snakemake,
        "conda": conda,
        "network": network_requirements(config, samples),
    }
