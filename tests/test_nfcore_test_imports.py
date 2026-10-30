"""Wrapper tests must exercise their owner even when a sibling is cached."""

import subprocess
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


@pytest.mark.parametrize("owner,filename,dependency,function", [
    ("rnaseq", "test_samplesheet_builder.py", "samplesheet_builder", "validate_and_normalize_samplesheet"),
    ("rnaseq", "test_provenance.py", "provenance", "write_provenance_bundle"),
    ("sarek", "test_samplesheet_builder.py", "samplesheet_builder", "validate_and_normalize_samplesheet"),
    ("sarek", "test_preflight.py", "preflight", "run_preflight"),
    ("sarek", "test_provenance.py", "provenance", "write_provenance_bundle"),
    ("scrnaseq", "test_preflight.py", "preflight", "run_preflight"),
    ("scrnaseq", "test_remap_paths.py", "outputs_parser", "find_multiqc_report"),
])
@pytest.mark.parametrize("foreign", ["rnaseq", "sarek", "scrnaseq"])
def test_reviewed_tests_keep_their_own_implementation(owner, filename, dependency, function, foreign):
    if owner == foreign:
        pytest.skip("the other two wrappers supply the foreign module")
    code = r"""
import importlib.util
import inspect
import sys
from pathlib import Path

root, owner, filename, dependency, function, foreign = sys.argv[1:]
root = Path(root)

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    spec.loader.exec_module(module)
    return module

foreign_dir = root / "skills" / f"nfcore-{foreign}-wrapper"
foreign_module = load("foreign_dependency", foreign_dir / f"{dependency}.py")
owner_dir = root / "skills" / f"nfcore-{owner}-wrapper"
# Restore the owner's sibling context, then leave just the colliding primary
# name cached, as can happen between test collection and a runtime import.
load("owner_conftest", owner_dir / "tests/conftest.py")
sys.modules[dependency] = foreign_module
tested = load("reviewed_test", owner_dir / "tests" / filename)
implementation = getattr(tested, function, None)
if implementation is None:
    implementation = getattr(tested.preflight, function)
assert Path(inspect.getfile(implementation)).resolve() == owner_dir / f"{dependency}.py"
"""
    result = subprocess.run(
        [sys.executable, "-c", code, str(ROOT), owner, filename, dependency, function, foreign],
        cwd=ROOT, capture_output=True, text=True, check=False,
    )
    assert result.returncode == 0, result.stderr


@pytest.fixture
def isolated_test_modules():
    previous = sys.modules.copy()
    yield
    for name in set(sys.modules) - previous.keys():
        if name.startswith("_nfcore_test_"):
            sys.modules.pop(name, None)
    if "preflight" in previous:
        sys.modules["preflight"] = previous["preflight"]
    else:
        sys.modules.pop("preflight", None)


def test_explicit_test_import_reuses_the_canonical_module(tmp_path, monkeypatch, isolated_test_modules):
    import importlib.util
    from clawbio.common.tests.nfcore_test_imports import load_tested_module

    source = tmp_path / "preflight.py"
    source.write_text("from dataclasses import dataclass\n@dataclass\nclass Result:\n    value: int\n")
    spec = importlib.util.spec_from_file_location("preflight", source)
    canonical = importlib.util.module_from_spec(spec)
    monkeypatch.setitem(sys.modules, "preflight", canonical)
    spec.loader.exec_module(canonical)
    loaded = load_tested_module(tmp_path, "preflight")
    assert loaded is canonical
    monkeypatch.setattr(loaded, "Result", lambda value: value + 1)
    # A sibling takes the bare name between two imports of the owner.
    import types
    foreign = types.ModuleType("preflight")
    foreign.__file__ = str(tmp_path / "foreign" / "preflight.py")
    sys.modules["preflight"] = foreign
    assert load_tested_module(tmp_path, "preflight").Result(4) == 5
    assert sys.modules["preflight"] is canonical


def test_explicit_test_import_registers_dataclasses_and_caches_by_path(tmp_path, monkeypatch, isolated_test_modules):
    from clawbio.common.tests.nfcore_test_imports import load_tested_module

    for name, value in [("first", 7), ("second", 9)]:
        directory = tmp_path / name
        directory.mkdir()
        (directory / "preflight.py").write_text(
            "from dataclasses import dataclass\n@dataclass\nclass Result:\n"
            f"    value: int = {value}\n"
        )
    first = load_tested_module(tmp_path / "first", "preflight")
    second = load_tested_module(tmp_path / "second", "preflight")
    assert first.Result().value == 7
    assert second.Result().value == 9
    assert first is not second
    monkeypatch.setattr(first, "Result", lambda: "patched")
    assert load_tested_module(tmp_path / "first", "preflight").Result() == "patched"
    assert sys.modules["preflight"] is first
