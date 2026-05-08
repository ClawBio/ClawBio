from __future__ import annotations

import shlex
import subprocess
import sys
from pathlib import Path


def build_nextflow_command(
    *,
    pipeline_source: dict[str, str | bool],
    profile: str,
    params_path: Path,
    resume: bool,
    work_dir: Path | None = None,
    extra_configs: list[Path] | None = None,
) -> tuple[list[str], str]:
    source_kind = str(pipeline_source["source_kind"])
    source_ref = str(pipeline_source["source_ref"])
    resolved_version = str(pipeline_source["resolved_version"])

    command = ["nextflow", "run"]
    if source_kind == "local_checkout":
        # Convert to forward slashes: Nextflow (Java) accepts them on all platforms,
        # and they avoid backslash ambiguity when the path is logged or copy-pasted.
        command.append(Path(source_ref).as_posix())
    else:
        command.extend([source_ref, "-r", resolved_version])

    # Use forward-slash paths for -params-file and -c so the command is portable
    # across macOS, Linux, and Windows (native or via WSL2).
    command.extend(["-profile", profile, "-params-file", params_path.as_posix()])
    if work_dir is not None:
        command.extend(["-work-dir", work_dir.as_posix()])
    for cfg in extra_configs or []:
        command.extend(["-c", cfg.as_posix()])
    if resume:
        command.append("-resume")

    # command_str is the human-readable form stored in result.json / report.md.
    # Use platform-appropriate quoting so it can be pasted directly into a terminal.
    if sys.platform == "win32":
        command_str = subprocess.list2cmdline(command)
    else:
        command_str = " ".join(shlex.quote(part) for part in command)

    return command, command_str
