"""Bounded Git and file-presence probes; callers decide source selection."""

import re
import subprocess
from pathlib import Path


def git_stdout(path: Path, args: list[str], *, timeout_seconds: float = 10) -> str:
    try:
        proc = subprocess.run(
            ["git", *args],
            cwd=str(path),
            capture_output=True,
            text=True,
            check=True,
            timeout=timeout_seconds,
        )
    except (
        subprocess.CalledProcessError,
        FileNotFoundError,
        OSError,
        subprocess.TimeoutExpired,
    ):
        return ""
    return proc.stdout.strip()


def missing_pipeline_files(
    path: Path, required_files: list[str] | tuple[str, ...]
) -> list[str]:
    return [name for name in required_files if not (path / name).exists()]


_CONFIG_TOKEN_RE = re.compile(
    r'/\*.*?\*/|//[^\n]*|\'(?:\\.|[^\'\\])*\'|"(?:\\.|[^"\\])*"|[A-Za-z_][A-Za-z_0-9]*|[^\s]',
    re.DOTALL,
)

def read_manifest_version(path: Path) -> str:
    """Read a literal top-level manifest version without interpreting Groovy."""
    try:
        tokens = [
            token
            for token in _CONFIG_TOKEN_RE.findall(path.read_text(encoding="utf-8"))
            if not token.startswith(("//", "/*"))
        ]
    except (OSError, UnicodeError):
        return ""
    depth, in_manifest = 0, False
    for index, token in enumerate(tokens):
        if (
            token == "manifest"
            and depth == 0
            and tokens[index + 1 : index + 2] == ["{"]
        ):
            in_manifest = True
        elif token == "{":
            depth += 1
        elif token == "}":
            depth -= 1
            if in_manifest and depth == 0:
                return ""
        elif in_manifest and depth == 1 and token == "version":
            if tokens[index - 1 : index] == ["."]:
                continue
            assignment = tokens[index + 1 : index + 3]
            if (
                len(assignment) == 2
                and assignment[0] == "="
                and assignment[1][0] in "'\""
            ):
                return assignment[1][1:-1].strip()
    return ""
