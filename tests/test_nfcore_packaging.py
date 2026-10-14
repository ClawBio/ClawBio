"""Qualified nf-core helpers must not load a skill's bare modules."""

import json
import os
import subprocess
import sys
from pathlib import Path

import pytest


@pytest.mark.parametrize(
    "module",
    [
        "errors",
        "imports",
        "paths",
        "commands",
        "parameters",
        "execution",
        "runtime",
        "source",
        "artifacts",
        "bundle_runtime",
        "bundle_emit",
        "cli",
        "failures",
    ],
)
def test_nfcore_modules_do_not_import_bare_skill_dependencies(module, tmp_path):
    root = Path(__file__).resolve().parent.parent
    env = dict(os.environ, PYTHONPATH=str(root))
    code = "import importlib,json,sys; importlib.import_module(sys.argv[1]); print(json.dumps(sorted({'errors','schemas','preflight','executor','params_builder'} & set(sys.modules))))"
    result = subprocess.run(
        [sys.executable, "-c", code, "clawbio.common.nfcore." + module],
        cwd=tmp_path,
        env=env,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
    assert json.loads(result.stdout) == []
