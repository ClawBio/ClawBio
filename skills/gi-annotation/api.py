"""API entry point for skill registry and orchestrator integration.

Invokes the skill CLI as a subprocess and returns the parsed
``result.json`` so programmatic callers (`ClawBio` runner, the
orchestrator) get the same output shape as a CLI invocation.

Note: this skill always submits annotation asynchronously (the API accepts
either mode); the underlying CLI handles
the submit → poll → render cycle transparently. A demo run typically
takes ~8–25 s depending on GPU warmth.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

_SKILL_DIR = Path(__file__).resolve().parent
_SCRIPT = _SKILL_DIR / "gi_annotation.py"


def run(input_path: Optional[str] = None, output_dir: Optional[str] = None, demo: bool = False) -> Dict[str, Any]:
    """Run gi-annotation and return the parsed ``result.json``.

    Either ``input_path`` or ``demo=True`` must be provided.
    """
    out = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="gi-annotation-"))
    cmd = [sys.executable, str(_SCRIPT), "--output", str(out)]
    if demo or input_path is None:
        cmd.append("--demo")
    else:
        cmd.extend(["--input", str(input_path)])
    subprocess.run(cmd, check=True)
    return json.loads((out / "result.json").read_text())
