"""Tests for scarf-single-cell.

The SKILL.md and CLI checks run everywhere. The end-to-end demo needs Python 3.12+ and
scarf>=1.0.0rc17 (not part of ClawBio's lockfile) and is skipped when Scarf is missing.
"""

from __future__ import annotations

import hashlib
import importlib.util
import json
import os
import re
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
SCRIPT = SKILL_DIR / "scarf_single_cell.py"
SKILL_MD = SKILL_DIR / "SKILL.md"
GENERATOR = SKILL_DIR / "examples" / "make_demo_data.py"
ENV = {**os.environ, "SCARF_MEM_BUDGET": "2G", "SCARF_WORKERS": "2", "MPLBACKEND": "Agg"}


def _scarf_available() -> bool:
    if sys.version_info < (3, 12) or importlib.util.find_spec("scarf") is None:
        return False
    from importlib.metadata import version

    from packaging.version import Version

    return Version(version("scarf")) >= Version("1.0.0rc17")


needs_scarf = pytest.mark.skipif(
    not _scarf_available(), reason="needs Python 3.12+ and scarf[extra]>=1.0.0rc17"
)


def _frontmatter() -> dict:
    import yaml

    text = SKILL_MD.read_text(encoding="utf-8")
    match = re.match(r"^---\s*\n(.*?)\n---", text, re.DOTALL)
    assert match, "SKILL.md has no YAML frontmatter"
    return yaml.safe_load(match.group(1))


def _run(*args: str, cwd: Path = PROJECT_ROOT) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(SCRIPT), *args], capture_output=True, text=True, env=ENV, cwd=cwd
    )


class TestSkillDocument:
    def test_frontmatter_conformance(self):
        meta = _frontmatter()
        assert meta["name"] == SKILL_DIR.name
        assert meta["license"]
        assert re.fullmatch(r"\d+\.\d+\.\d+", meta["metadata"]["version"])
        assert meta["metadata"]["author"]
        assert meta["metadata"]["inputs"] and meta["metadata"]["outputs"]
        assert len(meta["metadata"]["openclaw"]["trigger_keywords"]) >= 3

    def test_install_pins_the_prerelease_floor(self):
        """A bare scarf[extra] resolves to 0.32, whose API the skill does not describe."""
        meta = _frontmatter()
        packages = [item["package"] for item in meta["metadata"]["openclaw"]["install"]]
        assert "scarf[extra]>=1.0.0rc17" in packages
        assert "scarf[extra]>=1.0.0rc17" in meta["metadata"]["dependencies"]["packages"]
        assert "1.0.0rc17" in meta["compatibility"]

    def test_required_sections_present(self):
        text = SKILL_MD.read_text(encoding="utf-8")
        for section in (
            "Trigger", "Scope", "Workflow", "Example Output", "Output Structure", "Gotchas",
            "Safety", "Agent Boundary", "Integration with Bio Orchestrator", "Maintenance",
        ):
            assert re.search(rf"^## {re.escape(section)}\s*$", text, re.M), section
        assert "not a medical device" in text

    def test_under_500_lines(self):
        assert len(SKILL_MD.read_text(encoding="utf-8").splitlines()) < 500

    def test_referenced_modules_exist(self):
        text = SKILL_MD.read_text(encoding="utf-8")
        for rel in set(re.findall(r"(references/[a-z-]+\.md|scripts/[a-z_]+\.py)", text)):
            assert (SKILL_DIR / rel).is_file(), rel


class TestCLI:
    def test_help_works_without_scarf(self):
        result = _run("--help")
        assert result.returncode == 0
        assert "--demo" in result.stdout and "--input" in result.stdout

    def test_no_args_exits_nonzero(self):
        assert _run().returncode != 0

    def test_holdout_must_differ_from_sample_column(self, tmp_path):
        result = _run("--demo", "--output", str(tmp_path), "--sample-column", "x",
                      "--holdout-column", "x")
        assert result.returncode != 0

    @pytest.mark.skipif(_scarf_available(), reason="Scarf is installed")
    def test_missing_scarf_gives_install_hint(self, tmp_path):
        result = _run("--demo", "--output", str(tmp_path / "out"))
        assert result.returncode == 2
        assert "scarf[extra]>=1.0.0rc17" in result.stderr


class TestDemoData:
    def test_generator_is_seeded_raw_counts(self, tmp_path):
        ad = pytest.importorskip("anndata")
        np = pytest.importorskip("numpy")
        spec = importlib.util.spec_from_file_location("scarf_single_cell_demo_data", GENERATOR)
        module = importlib.util.module_from_spec(spec)
        spec.loader.exec_module(module)
        first = ad.read_h5ad(module.make_demo_h5ad(tmp_path / "a.h5ad"))
        second = ad.read_h5ad(module.make_demo_h5ad(tmp_path / "b.h5ad"))
        assert first.shape == (2000, 1500)
        assert (first.X != second.X).nnz == 0
        values = first.X.data
        assert np.all(values >= 0) and np.all(values == np.round(values))
        assert set(first.obs["donor_id"]) == {"D1", "D2"}
        assert first.obs["author_cell_type"].nunique() == 4


@pytest.fixture(scope="module")
def demo_output(tmp_path_factory):
    if not _scarf_available():
        pytest.skip("needs Python 3.12+ and scarf[extra]>=1.0.0rc17")
    out = tmp_path_factory.mktemp("scarf_demo") / "out"
    result = _run("--demo", "--output", str(out))
    assert result.returncode == 0, result.stderr[-3000:]
    return out


@needs_scarf
class TestDemoRun:
    def test_core_outputs(self, demo_output):
        for rel in ("report.md", "result.json", "run_report.md", "handoff/run_cells.h5ad",
                    "figures/umap_clusters.png", "tables/markers.csv"):
            assert (demo_output / rel).is_file(), rel

    def test_result_json(self, demo_output):
        result = json.loads((demo_output / "result.json").read_text())
        assert result["skill"] == "scarf-single-cell"
        summary, data = result["summary"], result["data"]
        assert summary["run_status"] == "completed"
        assert summary["cells_total"] == 2000
        assert 1900 <= summary["cells_in_run"] <= 2000
        assert summary["clusters"] == 4
        assert summary["handoff_raw_count_guard"] == "accepted"
        assert data["resources"] == {"mem_budget": "2G", "workers": 2}
        assert data["holdout"]["adjusted_rand_index"] >= 0.9

    def test_report_has_disclaimer_and_audits(self, demo_output):
        report = (demo_output / "report.md").read_text()
        assert "not a medical device" in report
        assert "QC retention by `donor_id`" in report
        assert "Held-out comparison: `author_cell_type`" in report

    def test_holdout_values_hidden_from_store_profile(self, demo_output):
        profile = json.loads((demo_output / "tables" / "store_profile.json").read_text())
        column = profile["columns"]["author_cell_type"]
        assert "annotation_like" in column["flags"]
        assert column["values"] == "hidden"

    def test_handoff_loads_with_scrna_loader(self, demo_output):
        ad = pytest.importorskip("anndata")
        sys.path.insert(0, str(PROJECT_ROOT))
        from clawbio.common.scrna_io import load_count_adata

        adata, _ = load_count_adata(
            demo_output / "handoff" / "run_cells.h5ad",
            h5ad_loader=ad.read_h5ad,
            expected_input="raw-count .h5ad",
        )
        assert adata.X.dtype.kind == "f"
        assert "clusters" in adata.obs and "X_umap" in adata.obsm

    def test_checksums_verify(self, demo_output):
        lines = (demo_output / "reproducibility" / "checksums.sha256").read_text().splitlines()
        assert lines
        for line in lines:
            digest, label = line.split("  ", 1)
            assert hashlib.sha256((demo_output / label).read_bytes()).hexdigest() == digest, label

    def test_refuses_to_overwrite(self, demo_output):
        result = _run("--demo", "--output", str(demo_output))
        assert result.returncode == 1
        assert "does not overwrite" in result.stderr


def _parse_output_contract(skill_md):
    """Files promised in the SKILL.md '## Output Structure' tree (scaffold_skill.py logic)."""
    text = skill_md.read_text()
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files = []
    parents = {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        if not name:
            continue
        depth = len(prefix) // 4
        if depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(rel + "/" + name if rel else name)
    return files


@needs_scarf
class TestOutputContract:
    """Every artifact promised in SKILL.md '## Output Structure' must be produced."""

    def test_documented_outputs_are_produced(self, demo_output):
        promised = _parse_output_contract(SKILL_MD)
        assert promised, "No parseable '## Output Structure' section in SKILL.md"
        missing = [p for p in promised if not (demo_output / p).exists()]
        assert not missing, f"SKILL.md promises outputs the demo did not write: {missing}"
