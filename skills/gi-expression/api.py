"""API entry point for skill registry and orchestrator integration.

Invokes the skill CLI as a subprocess and returns the parsed
``result.json`` so programmatic callers (`ClawBio` runner, the
orchestrator) get the same output shape as a CLI invocation.
"""

from __future__ import annotations

import json
import subprocess
import sys
import tempfile
from pathlib import Path
from typing import Any, Dict, Optional

_SKILL_DIR = Path(__file__).resolve().parent
_SCRIPT = _SKILL_DIR / "gi_expression.py"


def run(input_path: Optional[str] = None, output_dir: Optional[str] = None, demo: bool = False, description: Optional[str] = None, tss_index: Optional[int] = None) -> Dict[str, Any]:
    """Run gi-expression and return the parsed ``result.json``.

    Either ``input_path`` or ``demo=True`` must be provided. ``description``
    overrides the default K562 cell-type prompt the expression model is
    conditioned on. ``tss_index`` is the 0-based TSS offset into the
    whitespace-stripped sequence — required unless the FASTA holds exactly
    one 9,198 bp window.
    """
    out = Path(output_dir) if output_dir else Path(tempfile.mkdtemp(prefix="gi-expression-"))
    cmd = [sys.executable, str(_SCRIPT), "--output", str(out)]
    if demo or input_path is None:
        cmd.append("--demo")
    else:
        cmd.extend(["--input", str(input_path)])
    if description is not None:
        cmd.extend(["--description", description])
    if tss_index is not None:
        cmd.extend(["--tss-index", str(tss_index)])
    subprocess.run(cmd, check=True)
    return json.loads((out / "result.json").read_text())
