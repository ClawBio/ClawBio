"""Shared reproducibility bundle integration for WGS-PRS.

The bundle records SHA-256 digests of FASTQ/VCF inputs without storing their
paths or contents. Non-demo replay scripts require caller-supplied environment
variables for those private files and for free-text trait queries.

Score selection is recorded in the same order ``WgsToPrsBridge._run_stage3``
resolves it (``--pgs-id``, then ``--trait``, then ``--panel-id``), except that
``--demo`` short-circuits to a self-contained replay that regenerates the
synthetic VCF. If that order in ``_run_stage3`` changes, change ``_selection``
here in the same commit.
"""

from __future__ import annotations

import hashlib
import json
import re
import shlex
import sys
from collections.abc import Iterable
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

_PROJECT_ROOT = Path(__file__).resolve().parents[2]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from clawbio.common.checksums import sha256_file  # noqa: E402
from clawbio.common.reproducibility import (  # noqa: E402
    ReproCommand,
    ReproPath,
    write_checksums,
    write_environment_yml,
    write_portable_commands_sh,
)
from clawbio.common.textio import write_text_lf  # noqa: E402

SCHEMA_VERSION = 1

# Importing clawbio.common pulls audit (opentelemetry) and scrna_io (numpy,
# pandas) through package __init__. requests is required because stage 3
# invokes gwas-prs. nextflow/docker/bcftools stay undeclared here: they are
# binaries, not pip packages, and SKILL.md already lists them.
REPLAY_PIP_DEPENDENCIES: tuple[str, ...] = (
    "requests>=2.31",
    "opentelemetry-sdk>=1.20,<2",
    "numpy>=1.24",
    "pandas>=2.0",
)

_VERSION_LINE = re.compile(r"^\s*version:\s*['\"]?([0-9]+\.[0-9]+\.[0-9]+)['\"]?\s*$")


def _read_skill_version() -> str:
    """The ``metadata.version`` declared in this skill's SKILL.md frontmatter."""
    text = (Path(__file__).with_name("SKILL.md")).read_text(encoding="utf-8")
    in_frontmatter = False
    for line in text.splitlines():
        if line.strip() == "---":
            if in_frontmatter:
                break
            in_frontmatter = True
            continue
        if not in_frontmatter:
            continue
        match = _VERSION_LINE.match(line)
        if match:
            return match.group(1)
    raise RuntimeError("wgs-prs SKILL.md does not declare metadata.version")


TOOL_VERSION = _read_skill_version()


def create_reproducibility_bundle(
    *,
    output_dir: Path | str,
    config: Any,
    fastq_r1: Path | str | None = None,
    fastq_r2: Path | str | None = None,
    input_vcf: Path | str | None = None,
    samplesheet: Path | str | None = None,
    output_paths: Iterable[Path | str] = (),
) -> dict[str, Path]:
    """Write commands, environment, provenance, and checksum artefacts."""

    output_dir = Path(output_dir)
    fastq_r1 = Path(fastq_r1) if fastq_r1 else None
    fastq_r2 = Path(fastq_r2) if fastq_r2 else None
    input_vcf = Path(input_vcf) if input_vcf else None
    samplesheet = Path(samplesheet) if samplesheet else None

    commands_path = write_portable_commands_sh(
        output_dir,
        _repro_command(
            config,
            output_dir,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            input_vcf=input_vcf,
            samplesheet=samplesheet,
        ),
        repo_root=None,
    )
    environment_path = write_environment_yml(
        output_dir,
        env_name="clawbio-wgs-prs",
        pip_deps=list(REPLAY_PIP_DEPENDENCIES),
        python_version="3.11",
    )

    provenance_path = output_dir / "reproducibility" / "provenance.json"
    provenance = {
        "schema_version": SCHEMA_VERSION,
        "tool": {"name": "ClawBio WGS-PRS", "version": TOOL_VERSION},
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "input": _input_provenance(
            config,
            fastq_r1=fastq_r1,
            fastq_r2=fastq_r2,
            input_vcf=input_vcf,
        ),
        "parameters": _safe_parameters(config),
    }
    write_text_lf(provenance_path, json.dumps(provenance, indent=2) + "\n")

    checksum_targets = [
        *(Path(path) for path in output_paths),
        commands_path,
        environment_path,
        provenance_path,
    ]
    checksums_path = write_checksums(
        checksum_targets,
        output_dir,
        anchor=output_dir,
    )
    return {
        "commands": commands_path,
        "environment": environment_path,
        "provenance": provenance_path,
        "checksums": checksums_path,
    }


def _selection(config: Any) -> tuple[str, str | None]:
    """Which selector ``_run_stage3`` acted on.

    Returns ``(mode, value)``: ``("demo", None)``, ``("pgs_id", id)``,
    ``("trait", query)``, ``("panel_id", id)`` or ``("default", None)``.
    """
    if bool(getattr(config, "demo", False)):
        return "demo", None
    pgs_ids = list(getattr(config, "pgs_ids", []) or [])
    if pgs_ids:
        return "pgs_id", str(pgs_ids[0])
    traits = list(getattr(config, "prs_traits", []) or [])
    if traits:
        return "trait", str(traits[0])
    panel_id = str(getattr(config, "panel_id", "") or "")
    if panel_id:
        return "panel_id", panel_id
    return "default", None


def _repro_command(
    config: Any,
    output_dir: Path,
    *,
    fastq_r1: Path | None,
    fastq_r2: Path | None,
    input_vcf: Path | None,
    samplesheet: Path | None,
) -> ReproCommand:
    command_args: list[str | ReproPath] = []
    preflight: list[str] = []
    mode, value = _selection(config)

    if mode == "demo":
        command_args.append("--demo")
    elif input_vcf is not None:
        preflight.append(
            ': "${INPUT_VCF:?Set INPUT_VCF to the VCF used for this run}"'
        )
        command_args.extend(["--input-vcf", '"${INPUT_VCF}"'])
    else:
        preflight.append(
            ': "${FASTQ_R1:?Set FASTQ_R1 to the forward reads used for this run}"'
        )
        command_args.extend(["--fastq-r1", '"${FASTQ_R1}"'])
        if fastq_r2 is not None:
            preflight.append(
                ': "${FASTQ_R2:?Set FASTQ_R2 to the reverse reads used for this run}"'
            )
            command_args.extend(["--fastq-r2", '"${FASTQ_R2}"'])
        if samplesheet is not None:
            preflight.append(
                ': "${SAMPLESHEET:?Set SAMPLESHEET to the sarek samplesheet used for this run}"'
            )
            command_args.extend(["--samplesheet", '"${SAMPLESHEET}"'])

    if mode == "pgs_id":
        command_args.extend(["--pgs-id", shlex.quote(str(value))])
    elif mode == "trait":
        preflight.append(
            ': "${TRAIT_QUERY:?Set TRAIT_QUERY to the trait used for this run}"'
        )
        command_args.extend(["--trait", '"${TRAIT_QUERY}"'])
    elif mode == "panel_id":
        command_args.extend(["--panel-id", shlex.quote(str(value))])

    command_args.extend(
        [
            "--output-dir",
            ReproPath(output_dir, anchor="output_dir"),
            "--sample-id",
            shlex.quote(str(config.sample_id)),
            "--sex",
            str(config.sex),
            "--genome",
            str(config.sarek.genome),
            "--profile",
            str(config.sarek.profile),
            "--sarek-version",
            str(config.sarek.sarek_version),
            "--min-qual",
            str(config.qc.min_qual),
            "--min-depth",
            str(config.qc.min_depth),
        ]
    )
    if config.sarek.skip_bqsr:
        command_args.append("--skip-bqsr")
    if not config.fail_fast:
        command_args.append("--no-fail-fast")
    if config.dry_run:
        command_args.append("--dry-run")
    ref = str(getattr(config.qc, "reference_fasta", "") or "")
    if ref:
        preflight.append(
            ': "${REFERENCE_FASTA:?Set REFERENCE_FASTA to the FASTA used for bcftools norm}"'
        )
        command_args.extend(["--reference-fasta", '"${REFERENCE_FASTA}"'])

    return ReproCommand(
        script_path=Path("skills/wgs-prs/wgs_prs.py"),
        args=command_args,
        comment="Replay this ClawBio WGS-PRS run",
        preflight=preflight,
    )


def _safe_parameters(config: Any) -> dict[str, Any]:
    mode, value = _selection(config)
    selection: dict[str, Any] = {"mode": mode}
    if mode == "panel_id":
        selection["panel_id"] = value
    elif mode == "pgs_id":
        selection["pgs_id"] = value
    elif mode == "trait":
        selection["query_sha256"] = _hash_text(str(value))
    return {
        "selection": selection,
        "sample_id": str(config.sample_id),
        "sex": str(config.sex),
        "genome": str(config.sarek.genome),
        "profile": str(config.sarek.profile),
        "sarek_version": str(config.sarek.sarek_version),
        "skip_bqsr": bool(config.sarek.skip_bqsr),
        "min_qual": float(config.qc.min_qual),
        "min_depth": int(config.qc.min_depth),
        "fail_fast": bool(config.fail_fast),
        "dry_run": bool(config.dry_run),
    }


def _input_provenance(
    config: Any,
    *,
    fastq_r1: Path | None,
    fastq_r2: Path | None,
    input_vcf: Path | None,
) -> dict[str, Any]:
    if bool(getattr(config, "demo", False)):
        entry = "demo"
    elif input_vcf is not None:
        entry = "vcf"
    else:
        entry = "fastq"
    files = []
    for role, path in (
        ("input_vcf", input_vcf),
        ("fastq_r1", fastq_r1),
        ("fastq_r2", fastq_r2),
    ):
        digest = _file_digest(path)
        if digest is None:
            continue
        files.append({"role": role, "sha256": digest})
    return {"entry": entry, "files": files}


def _file_digest(path: Path | None) -> str | None:
    if path is None or not path.is_file():
        return None
    return sha256_file(path)


def _hash_text(value: str) -> str:
    return f"sha256:{hashlib.sha256(value.encode('utf-8')).hexdigest()}"
