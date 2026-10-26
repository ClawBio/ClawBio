"""Shared error serialization and explicit module ownership contracts."""

import importlib
import itertools
import subprocess
import sys
import types
from pathlib import Path

import pytest


def test_skill_error_payload_is_unchanged():
    module = importlib.import_module("clawbio.common.nfcore.errors")
    error = module.SkillError("validation", "INVALID_FASTQ", "bad", "fix", {"line": 3})
    assert str(error) == "bad"
    assert error.to_dict() == {
        "status": "error",
        "ok": False,
        "stage": "validation",
        "error_code": "INVALID_FASTQ",
        "message": "bad",
        "fix": "fix",
        "details": {"line": 3},
    }
    assert (
        module.SkillError("s", "c", "m", "f").details
        is not module.SkillError("s", "c", "m", "f").details
    )


def test_purge_respects_explicit_owner(tmp_path, monkeypatch):
    module = importlib.import_module("clawbio.common.nfcore.imports")
    owner = tmp_path / "owner"
    local = types.ModuleType("nfcore_test_owned")
    local.__file__ = str(owner / "nfcore_test_owned.py")
    foreign = types.ModuleType("nfcore_test_foreign")
    foreign.__file__ = str(tmp_path / "other" / "nfcore_test_foreign.py")
    monkeypatch.setitem(sys.modules, local.__name__, local)
    monkeypatch.setitem(sys.modules, foreign.__name__, foreign)
    module.purge_foreign_bare_modules(owner, local.__name__, foreign.__name__)
    assert sys.modules[local.__name__] is local
    assert foreign.__name__ not in sys.modules


@pytest.mark.parametrize(
    "order", list(itertools.permutations(["rnaseq", "sarek", "scrnaseq"]))
)
def test_wrapper_exception_identity_survives_loading_order(order):
    root = Path(__file__).resolve().parents[3]
    code = """
import importlib.util, pathlib, sys
from clawbio.common.nfcore.errors import SkillError
root = pathlib.Path(sys.argv[1])
for name in [*sys.argv[2:], sys.argv[2]]:
    p = root / 'skills' / ('nfcore-' + name + '-wrapper') / ('nfcore_' + name + '_wrapper.py')
    spec = importlib.util.spec_from_file_location('loaded_' + name, p)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    assert module.SkillError is SkillError
    import errors, schemas
    assert errors.SkillError is SkillError
    assert schemas.SKILL_NAME == 'nfcore-' + name + '-wrapper'
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(root), *order],
        cwd=root,
        capture_output=True,
        text=True,
        check=False,
    )
    assert result.returncode == 0, result.stderr
