"""Shared reproducibility bundle integration for GWAS-PRS.

The bundle records input and scoring-file digests without storing genotype
paths or contents. Non-demo replay scripts require caller-supplied environment
variables for private input paths and free-text trait queries.

Score selection is recorded in the same order ``main()`` in ``gwas_prs.py``
resolves it (``--demo``, then ``--panel-id``, then ``--pgs-id``, then
``--trait``), so the replay command and ``provenance.json`` always describe the
run that actually happened. The CLI rejects more than one selector, but the
bundle does not rely on that: if the order in ``main()`` changes, change
``_selection`` here in the same commit.
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

# Everything the replay environment must install before ``gwas_prs.py`` can
# import. ``requests`` is imported by the skill itself. The other three arrive
# through ``clawbio.common``: its package ``__init__`` imports ``audit``
# (opentelemetry) and ``scrna_io`` (numpy, pandas) eagerly, so any skill that
# imports ``clawbio.common.checksums`` pays for them even though this skill
# never calls them. Verified by running ``--demo`` in a clean interpreter with
# only ``requests`` and ``opentelemetry-sdk`` installed: it fails on
# ``import numpy``. ``scipy`` and ``matplotlib`` are optional at runtime and
# are deliberately not declared.
REPLAY_PIP_DEPENDENCIES: tuple[str, ...] = (
    "requests>=2.31",
    "opentelemetry-sdk>=1.20,<2",
    "numpy>=1.24",
    "pandas>=2.0",
)

_VERSION_LINE = re.compile(r"^\s*version:\s*['\"]?([0-9]+\.[0-9]+\.[0-9]+)['\"]?\s*$")


def _read_skill_version() -> str:
    """The ``metadata.version`` declared in this skill's SKILL.md frontmatter.

    Matches the first ``version:`` key inside the leading ``---`` block at any
    indentation, so a re-indented frontmatter does not break import.
    """
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
    raise RuntimeError("gwas-prs SKILL.md does not declare metadata.version")


TOOL_VERSION = _read_skill_version()


def create_reproducibility_bundle(
    *,
    output_dir: Path | str,
    input_path: Path | str,
    input_info: dict[str, Any],
    scoring_files: list[dict[str, Any]],
    args: Any,
    output_paths: Iterable[Path | str],
) -> dict[str, Path]:
    """Write commands, environment, provenance, and checksum artefacts."""

    output_dir = Path(output_dir)
    input_path = Path(input_path)

    commands_path = write_portable_commands_sh(
        output_dir,
        _repro_command(args, output_dir),
        repo_root=None,
    )
    environment_path = write_environment_yml(
        output_dir,
        env_name="clawbio-gwas-prs",
        pip_deps=list(REPLAY_PIP_DEPENDENCIES),
        python_version="3.11",
    )

    provenance_path = output_dir / "reproducibility" / "provenance.json"
    provenance = {
        "schema_version": SCHEMA_VERSION,
        "tool": {"name": "ClawBio GWAS-PRS", "version": TOOL_VERSION},
        "generated_at": datetime.now(timezone.utc).strftime("%Y-%m-%dT%H:%M:%SZ"),
        "input": {
            "sha256": sha256_file(input_path),
            "format": str(input_info.get("format", "unknown")),
            "total_snps": int(input_info.get("total_snps", 0)),
        },
        "parameters": _safe_parameters(args),
        "scoring_files": [
            _scoring_file_provenance(item)
            for item in sorted(
                scoring_files,
                key=lambda item: (
                    str(item.get("score_id") or item.get("pgs_id") or ""),
                    Path(item["filepath"]).name,
                ),
            )
        ],
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


def _selection(args: Any) -> tuple[str, str | None]:
    """Which selector ``main()`` acted on, in ``main()``'s own precedence.

    Returns ``(mode, value)``: ``("demo", None)``, ``("panel_id", id)``,
    ``("pgs_id", id)`` or ``("trait", query)``.
    """
    if bool(args.demo):
        return "demo", None
    if args.panel_id:
        return "panel_id", str(args.panel_id)
    if args.pgs_id:
        return "pgs_id", str(args.pgs_id)
    if args.trait:
        return "trait", str(args.trait)
    raise ValueError("no score selector set: expected --demo, --panel-id, --pgs-id or --trait")


def _repro_command(args: Any, output_dir: Path) -> ReproCommand:
    command_args: list[str | ReproPath] = []
    preflight: list[str] = []
    mode, value = _selection(args)
    if mode == "demo":
        command_args.append("--demo")
    else:
        preflight.append(
            ': "${INPUT_FILE:?Set INPUT_FILE to the genotype file used for this run}"'
        )
        command_args.extend(["--input", '"${INPUT_FILE}"'])
        if mode == "panel_id":
            command_args.extend(["--panel-id", shlex.quote(str(value))])
        elif mode == "pgs_id":
            command_args.extend(["--pgs-id", shlex.quote(str(value))])
        else:
            preflight.append(
                ': "${TRAIT_QUERY:?Set TRAIT_QUERY to the trait used for this run}"'
            )
            command_args.extend(["--trait", '"${TRAIT_QUERY}"'])

    command_args.extend(
        [
            "--output",
            ReproPath(output_dir, anchor="output_dir"),
            "--min-overlap",
            str(args.min_overlap),
            "--max-variants",
            str(args.max_variants),
            "--build",
            str(args.build),
            "--cache-dir",
            '"${PGS_CACHE_DIR:-$HOME/.clawbio/pgs_cache}"',
        ]
    )
    if bool(args.no_cache):
        command_args.append("--no-cache")
    return ReproCommand(
        script_path=Path("skills/gwas-prs/gwas_prs.py"),
        args=command_args,
        comment="Replay this ClawBio GWAS-PRS run",
        preflight=preflight,
    )


def _safe_parameters(args: Any) -> dict[str, Any]:
    mode, value = _selection(args)
    selection: dict[str, Any] = {"mode": mode}
    if mode == "panel_id":
        selection["panel_id"] = value
    elif mode == "pgs_id":
        selection["pgs_id"] = value
    elif mode == "trait":
        # A fingerprint, not anonymisation: trait names are a small search
        # space, so this digest only lets a replay confirm it used the same
        # query without the query itself being written into the bundle.
        selection["query_sha256"] = _hash_text(str(value))
    return {
        "selection": selection,
        "build": str(args.build),
        "min_overlap": float(args.min_overlap),
        "max_variants": int(args.max_variants),
        "cache_enabled": not bool(args.no_cache),
    }


def _scoring_file_provenance(item: dict[str, Any]) -> dict[str, Any]:
    path = Path(item["filepath"])
    pgs_id = item.get("pgs_id")
    return {
        "score_id": str(item.get("score_id") or pgs_id or ""),
        "pgs_id": None if pgs_id is None else str(pgs_id),
        "filename": path.name,
        "sha256": sha256_file(path),
        "curated_demo_panel": bool(item.get("curated_demo_panel")),
        "curated_panel_id": item.get("curated_panel_id"),
        "legacy_pgs_id": item.get("legacy_pgs_id"),
        "legacy_pgs_compatibility": bool(item.get("legacy_pgs_compatibility")),
        "pgs_catalog_id": item.get("pgs_catalog_id"),
    }


def _hash_text(value: str) -> str:
    return f"sha256:{hashlib.sha256(value.encode('utf-8')).hexdigest()}"
