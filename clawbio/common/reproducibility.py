"""Reproducibility helpers for ClawBio skills.

Provides write_checksums and write_environment_yml, both writing into
<output_dir>/reproducibility/.
"""

from __future__ import annotations

from pathlib import Path

from clawbio.common.checksums import sha256_file


def write_checksums(
    paths: list[Path | str],
    output_dir: Path | str,
    anchor: Path | str | None = None,
) -> Path:
    """Write sha256sum-compatible checksums for output files.

    Each line: '<sha256>  <label>'
    - If anchor is None, label is the bare filename.
    - If anchor is given, label is the path relative to anchor.

    Files that do not exist are silently skipped.
    Creates reproducibility/ if it doesn't exist.
    Returns the path of the written checksums.sha256 file.
    """
    output_dir = Path(output_dir)
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    anchor_path = Path(anchor) if anchor is not None else None
    lines: list[str] = []
    for p in paths:
        p = Path(p)
        if not p.exists():
            continue
        if anchor_path is not None:
            try:
                label = str(p.relative_to(anchor_path))
            except ValueError:
                label = p.name
        else:
            label = p.name
        lines.append(f"{sha256_file(p)}  {label}")

    checksum_path = repro_dir / "checksums.sha256"
    checksum_path.write_text("\n".join(lines) + ("\n" if lines else ""))
    return checksum_path


def write_environment_yml(
    output_dir: Path | str,
    env_name: str,
    pip_deps: list[str],
    conda_deps: list[str] | None = None,
    python_version: str = "3.10",
) -> Path:
    """Write reproducibility/environment.yml for a ClawBio skill.

    Args:
        output_dir:     Skill output directory.
        env_name:       Conda environment name (e.g. 'clawbio-cell-detection').
        pip_deps:       Packages to install via pip (e.g. ['cellpose>=4.0']).
        conda_deps:     Extra conda packages beyond python (e.g. ['numpy', 'scipy']).
                        Do not include 'python=X.Y' here — use python_version instead.
        python_version: Python version string (default '3.10').

    Returns the path of the written environment.yml file.
    """
    output_dir = Path(output_dir)
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    # Strip any python= entries from conda_deps to avoid duplicating the python line.
    filtered_conda = [d for d in (conda_deps or []) if not d.lower().startswith("python=")]
    conda_lines = "\n".join(f"  - {dep}" for dep in filtered_conda)
    conda_block = f"\n{conda_lines}" if conda_lines else ""

    # Only emit the pip block when there are pip deps — empty pip: is invalid YAML for conda.
    if pip_deps:
        pip_lines = "\n".join(f"      - {dep}" for dep in pip_deps)
        pip_block = f"  - pip:\n{pip_lines}\n"
    else:
        pip_block = ""

    content = f"""name: {env_name}
channels:
  - conda-forge
dependencies:
  - python={python_version}{conda_block}
  - pip
{pip_block}"""
    path = repro_dir / "environment.yml"
    path.write_text(content)
    return path


def write_commands_sh(output_dir: Path | str, command: str) -> Path:
    """Write reproducibility/commands.sh containing the exact command to reproduce a run.

    Args:
        output_dir: Skill output directory.
        command:    The full CLI command string (may be multi-line with continuations).

    Creates reproducibility/ if it doesn't exist.
    Returns the path of the written commands.sh file.
    """
    output_dir = Path(output_dir)
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    content = f"#!/usr/bin/env bash\n{command}\n"
    path = repro_dir / "commands.sh"
    path.write_text(content)
    path.chmod(path.stat().st_mode | 0o111)
    return path
