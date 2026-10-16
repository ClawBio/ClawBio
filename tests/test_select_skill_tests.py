import importlib.util
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def _load():
    path = ROOT / "scripts" / "select_skill_tests.py"
    spec = importlib.util.spec_from_file_location("select_skill_tests", path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


select = _load().select_skills


def _write(root: Path, rel: str, text: str = "") -> None:
    path = root / rel
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(text)


@pytest.fixture(scope="module")
def repo(tmp_path_factory):
    root = tmp_path_factory.mktemp("repo")
    _write(root, "clawbio/__init__.py")
    _write(root, "clawbio/cli.py", "from clawbio.common.report import x\n")
    _write(root, "clawbio/common/__init__.py")
    _write(root, "clawbio/common/report.py", "from .textio import write\n")
    _write(root, "clawbio/common/textio.py")
    _write(root, "clawbio/common/checksums.py")
    _write(root, "clawbio/common/data/table.json", "{}")
    # one skill per import style
    _write(root, "skills/uses-report/run.py", "from clawbio.common.report import x\n")
    _write(root, "skills/uses-checksums/run.py", "from clawbio.common import checksums\n")
    _write(root, "skills/uses-textio/run.py", "import clawbio.common.textio as t\n")
    _write(root, "skills/runs-cli/tests/test_cli.py", "CLI = ROOT / 'clawbio.py'\n")
    _write(root, "skills/standalone/run.py", "import json\n")
    _write(root, "skills/borrows-standalone/run.py", "p = SKILLS / 'standalone' / 'data'\n")
    _write(root, "skills/fakes-standalone/tests/test_x.py", "d = tmp_path / 'skills' / 'standalone'\n")
    _write(root, "skills/router/run.py", "ROUTES = {'standalone': 'fasta'}\n")
    _write(root, "skills/catalog.json", "{}")
    return root


SKILLS = ["borrows-standalone", "fakes-standalone", "router", "runs-cli", "standalone", "uses-checksums", "uses-report", "uses-textio"]


def test_change_inside_one_skill_selects_only_that_skill(repo):
    assert select(repo, ["skills/uses-report/run.py"]) == ["uses-report"]


def test_skill_that_reaches_another_skill_by_path_is_selected_with_it(repo):
    # router names 'standalone' in a routing table, not a path; fakes-standalone
    # builds a throwaway tmp_path/skills/standalone folder: neither is selected
    assert select(repo, ["skills/standalone/run.py"]) == ["borrows-standalone", "standalone"]


def test_files_under_skills_that_are_not_skill_folders_select_nothing(repo):
    assert select(repo, ["skills/catalog.json"]) == []


def test_deleted_skill_selects_nothing(repo):
    assert select(repo, ["skills/gone/run.py"]) == []


def test_repo_level_tests_and_docs_select_no_skills(repo):
    assert select(repo, ["tests/test_x.py", "README.md", "docs/a.md"]) == []


def test_shared_module_selects_its_importers_only(repo):
    assert select(repo, ["clawbio/common/checksums.py"]) == ["uses-checksums"]


def test_shared_module_selects_importers_of_modules_that_import_it(repo):
    # report imports textio (relative import), cli imports report
    assert select(repo, ["clawbio/common/textio.py"]) == ["runs-cli", "uses-report", "uses-textio"]


def test_clawbio_py_counts_as_the_cli(repo):
    assert select(repo, ["clawbio/cli.py"]) == ["runs-cli"]
    assert select(repo, ["clawbio.py"]) == ["runs-cli"]


def test_package_init_selects_only_skills_that_import_the_package(repo):
    # uses-textio imports clawbio.common.textio, which runs clawbio/common/__init__.py
    # implicitly; implicit imports are left to the runtime tests (see module docstring)
    assert select(repo, ["clawbio/common/__init__.py"]) == ["uses-checksums"]


def test_data_file_in_shared_package_selects_that_package_importers(repo):
    expected = ["runs-cli", "uses-checksums", "uses-report", "uses-textio"]
    assert select(repo, ["clawbio/common/data/table.json"]) == expected


@pytest.mark.parametrize("path", ["uv.lock", "pyproject.toml", "conftest.py", "pytest.ini"])
def test_environment_files_select_every_skill(repo, path):
    assert select(repo, [path]) == SKILLS


MIXED = [
    "skills/uses-report/run.py", "skills/standalone/run.py", "skills/catalog.json",
    "clawbio/common/checksums.py", "README.md",
]


def test_selection_does_not_depend_on_input_order(repo):
    assert select(repo, MIXED) == select(repo, list(reversed(MIXED)))


def test_another_change_never_drops_a_selected_skill(repo):
    before = select(repo, ["skills/standalone/run.py"])
    after = select(repo, ["skills/standalone/run.py", "clawbio/common/checksums.py"])
    assert set(before) < set(after)


def test_real_repo_one_skill_change_does_not_select_everything():
    selected = select(ROOT, ["skills/equity-scorer/equity_scorer.py"])
    assert "equity-scorer" in selected
    assert len(selected) < 10


def test_real_repo_shared_helper_selects_importers_not_all_skills():
    every = select(ROOT, ["uv.lock"])
    selected = select(ROOT, ["clawbio/common/checksums.py"])
    assert "pharmgx-reporter" in selected
    assert 0 < len(selected) < len(every)
