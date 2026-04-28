"""Reproducibility helpers for ClawBio skills.

Provides write_checksums, write_environment_yml, and portable-bundle helpers
(ReproPath, ReproCommand, write_portable_commands_sh), all writing into
<output_dir>/reproducibility/.
"""

from __future__ import annotations

import subprocess
from dataclasses import dataclass, field
from pathlib import Path
from typing import Union

from clawbio.common.checksums import sha256_file


@dataclass
class ReproPath:
    """A path reference that renders portably in reproducibility scripts.

    anchor controls how the path is rendered in commands.sh:
      "repo_root"  → $CLAWBIO_ROOT/<path relative to repo root>
      "output_dir" → $OUTPUT_DIR/<path relative to output dir>
      "auto"       → absolute path (for user-supplied inputs outside the repo)
    """

    path: Path
    anchor: str = "auto"

    def __post_init__(self) -> None:
        self.path = Path(self.path)


@dataclass
class ReproCommand:
    """Structured representation of a reproducibility command.

    Used by write_portable_commands_sh to generate a portable commands.sh.
    """

    script_path: Path
    args: list[Union[str, ReproPath]] = field(default_factory=list)
    comment: str = ""
    preflight: list[str] = field(default_factory=list)

    def __post_init__(self) -> None:
        self.script_path = Path(self.script_path)


def write_portable_commands_sh(
    output_dir: Path | str,
    command: ReproCommand,
    *,
    repo_root: Path | str | None = None,
) -> Path:
    """Write reproducibility/commands.sh with portable path references.

    The script sets CLAWBIO_ROOT (defaulting to the repo root at bundle
    generation time) and OUTPUT_DIR (auto-detected from the script location),
    then renders ReproPath arguments as shell variable expansions so the bundle
    can be replayed on any machine that has the repo checked out.
    """
    output_dir = Path(output_dir)
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(parents=True, exist_ok=True)

    repo_root_resolved = Path(repo_root).resolve() if repo_root else None
    default_root = str(repo_root_resolved) if repo_root_resolved else "/path/to/ClawBio"

    def render_arg(arg: str | ReproPath) -> str:
        if not isinstance(arg, ReproPath):
            return arg
        p = arg.path
        if arg.anchor == "repo_root" and repo_root_resolved is not None:
            try:
                rel = p.relative_to(repo_root_resolved)
                return f'"$CLAWBIO_ROOT/{rel}"'
            except ValueError:
                pass
        if arg.anchor == "output_dir":
            if p == output_dir:
                return '"$OUTPUT_DIR"'
            try:
                rel = p.relative_to(output_dir)
                return f'"$OUTPUT_DIR/{rel}"'
            except ValueError:
                pass
        return f'"{p}"'

    rendered = [render_arg(a) for a in command.args]
    script_ref = f'"$CLAWBIO_ROOT/{command.script_path}"'

    parts = [f"python {script_ref}"] + rendered
    if len(parts) <= 2:
        cmd_line = " ".join(parts)
    else:
        cmd_line = parts[0] + " \\\n  " + " \\\n  ".join(parts[1:])

    lines = ["#!/usr/bin/env bash"]
    if command.comment:
        lines.append(f"# {command.comment}")
    lines += [
        "set -euo pipefail",
        'SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"',
        'OUTPUT_DIR="$(dirname "$SCRIPT_DIR")"',
        f': "${{CLAWBIO_ROOT:={default_root}}}"',
        "",
    ]
    if command.preflight:
        lines.extend(command.preflight)
        lines.append("")
    lines.append(cmd_line)

    content = "\n".join(lines) + "\n"
    path = repro_dir / "commands.sh"
    path.write_text(content)
    path.chmod(path.stat().st_mode | 0o111)
    return path


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


def write_conda_lock(output_dir: Path | str) -> Path:
    """Write reproducibility/conda-lock.yml from an existing environment.yml.

    Runs ``conda-lock lock`` in the reproducibility directory. conda-lock
    defaults to multi-platform resolution and writes conda-lock.yml.

    Args:
        output_dir: Skill output directory containing reproducibility/environment.yml.

    Raises:
        FileNotFoundError: If reproducibility/environment.yml is missing.
        subprocess.CalledProcessError: If conda-lock exits with a non-zero status.

    Returns the path of the written conda-lock.yml file.
    """
    output_dir = Path(output_dir)
    repro_dir = output_dir / "reproducibility"
    if not (repro_dir / "environment.yml").exists():
        raise FileNotFoundError(repro_dir / "environment.yml")

    try:
        subprocess.run(
            ["conda-lock", "lock", "-f", "environment.yml"],
            cwd=repro_dir,
            check=True,
        )
    except FileNotFoundError:
        raise FileNotFoundError(
            "conda-lock is not installed. Install it with: pip install conda-lock"
        )
    return repro_dir / "conda-lock.yml"
