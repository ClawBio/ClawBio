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
        "launcher",
        "runtime",
        "resources",
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


def test_nfcore_has_one_package_and_no_competing_module():
    """The shared refactor and the audit must expose one unambiguous namespace."""
    from clawbio.common import nfcore

    package = Path(nfcore.__file__).resolve().parent
    assert package.name == "nfcore"
    assert Path(nfcore.__file__).name == "__init__.py"
    assert not package.with_suffix(".py").exists()


def test_wrapper_test_setup_preserves_public_package_api(tmp_path):
    """Collecting wrapper tests cannot replace the package with a partial stub."""
    root = Path(__file__).resolve().parent.parent
    code = """
import runpy, sys
runpy.run_path(sys.argv[1])
from clawbio import __all__, list_skills, run_skill, upload_profile
assert all(callable(value) for value in (list_skills, run_skill, upload_profile))
assert 'run_skill' in __all__
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(root / "skills/nfcore-rnaseq-wrapper/tests/conftest.py")],
        cwd=tmp_path,
        env={**os.environ, "PYTHONPATH": str(root)},
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0, result.stderr
