"""Run a generated download script against fake `curl`/`wget`/`unzip`.

The archive-fetch skills emit bash that real users run unattended, so the
properties that matter are only visible at run time: whether a value reached
the transfer tool byte-for-byte, and whether bash expanded anything on the way.
The shims record their argv and write the `-o`/`-O` target, so a script can be
executed end to end with no network.
"""

from __future__ import annotations

import os
import subprocess
from pathlib import Path

SEPARATOR = "--end-of-call--"

# argv is logged one value per line, each call closed by SEPARATOR. A value
# after -o/-O is the destination: the shim writes $SHIM_CONTENT there, so a
# test controls what the "download" produced.
_TRANSFER_SHIM = f"""#!/bin/bash
{{ for a in "$@"; do printf '%s\\n' "$a"; done; echo {SEPARATOR}; }} >> "$SHIM_LOG"
if [ "${{1:-}}" = "--help" ]; then exit 0; fi
prev=""
for a in "$@"; do
  if [ "$prev" = "-o" ] || [ "$prev" = "-O" ]; then
    mkdir -p "$(dirname "$a")"
    printf '%s' "${{SHIM_CONTENT:-data}}" > "$a"
  fi
  prev="$a"
done
"""

_RECORD_SHIM = f"""#!/bin/bash
{{ for a in "$@"; do printf '%s\\n' "$a"; done; echo {SEPARATOR}; }} >> "$SHIM_LOG"
"""


def run_with_shims(script: Path, cwd: Path, *, content: str = "data",
                   ) -> tuple[subprocess.CompletedProcess, list[list[str]]]:
    """Execute `script` with shims first on PATH. Returns (result, calls)."""
    bin_dir = cwd / "_shims"
    bin_dir.mkdir(parents=True, exist_ok=True)
    for name, body in (("curl", _TRANSFER_SHIM), ("wget", _TRANSFER_SHIM),
                       ("unzip", _RECORD_SHIM)):
        shim = bin_dir / name
        shim.write_text(body)
        shim.chmod(0o755)
    log = cwd / "_shim.log"
    log.write_text("")
    env = dict(os.environ, PATH=f"{bin_dir}{os.pathsep}{os.environ['PATH']}",
               SHIM_LOG=str(log), SHIM_CONTENT=content)
    result = subprocess.run(["bash", str(script)], cwd=cwd, env=env,
                            capture_output=True, text=True, timeout=60)
    calls, current = [], []
    for line in log.read_text().splitlines():
        if line == SEPARATOR:
            calls.append(current)
            current = []
        else:
            current.append(line)
    return result, calls
