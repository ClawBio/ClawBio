"""Tests for the deepspot-m virtual spatial transcriptomics skill.

Every test here runs without torch, without deepspotm and without the gated
model weights, so the suite is green on a base ClawBio checkout.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
MODULE_PATH = SKILL_DIR / "deepspot_m.py"

DISCLAIMER = (
    "ClawBio is a research and educational tool. "
    "It is not a medical device and does not provide clinical diagnoses. "
    "Consult a healthcare professional before making any medical decisions."
)


def load_module():
    spec = importlib.util.spec_from_file_location("deepspot_m", MODULE_PATH)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.fixture(scope="module")
def module():
    return load_module()


@pytest.fixture(scope="module")
def demo_output(tmp_path_factory, module):
    out = tmp_path_factory.mktemp("deepspot_demo")
    module.main(["--demo", "--output", str(out)])
    return out


# ---------------------------------------------------------------------------
# Tile validation
# ---------------------------------------------------------------------------


def test_validate_tile_size_accepts_224(module):
    assert module.validate_tile_size(224, 224) is None


@pytest.mark.parametrize("width,height", [(223, 224), (224, 223), (256, 256), (112, 112)])
def test_validate_tile_size_rejects_other_sizes(module, width, height):
    with pytest.raises(ValueError) as exc:
        module.validate_tile_size(width, height)
    assert "224x224" in str(exc.value)


def test_bundled_demo_tile_is_224_square(module):
    tile = module.load_tile(module.DEMO_TILE)
    assert (tile.width, tile.height) == (224, 224)


# ---------------------------------------------------------------------------
# Gene panel and embedding sources
# ---------------------------------------------------------------------------


def test_parse_genes_normalises_and_deduplicates(module):
    assert module.parse_genes(" braf , cd37 ,BRAF, COL1A1 ") == ["BRAF", "CD37", "COL1A1"]


def test_parse_genes_defaults_to_the_bundled_panel(module):
    assert module.parse_genes(None) == list(module.DEFAULT_GENES)


def test_parse_genes_rejects_an_empty_list(module):
    with pytest.raises(ValueError):
        module.parse_genes(" , , ")


def test_five_embedding_sources_are_supported(module):
    assert module.EMBEDDING_SOURCES == ("evo2", "orthrus", "prott5", "scgpt", "apertus")
    assert module.DEFAULT_SOURCE in module.EMBEDDING_SOURCES


def test_every_embedding_source_is_accepted_by_the_parser(module):
    parser = module.build_parser()
    for source in module.EMBEDDING_SOURCES:
        args = parser.parse_args(["--demo", "--output", "/tmp/x", "--source", source])
        assert args.source == source


def test_demo_mode_refuses_genes_outside_the_fixture(module):
    with pytest.raises(ValueError) as exc:
        module.load_demo_expression(["EPCAM", "NOTAGENE"])
    assert "NOTAGENE" in str(exc.value)


# ---------------------------------------------------------------------------
# Demo output contract
# ---------------------------------------------------------------------------


def test_demo_writes_every_documented_artifact(demo_output):
    expected = [
        "report.md",
        "result.json",
        "tables/gene_expression.csv",
        "reproducibility/commands.sh",
        "reproducibility/environment.yml",
        "reproducibility/checksums.sha256",
    ]
    missing = [name for name in expected if not (demo_output / name).exists()]
    assert missing == []


def test_demo_result_json_reports_log1p_cpm(demo_output, module):
    payload = json.loads((demo_output / "result.json").read_text(encoding="utf-8"))
    assert payload["unit"] == "log1p-CPM"
    assert payload["model"] == "ratschlab/DeepSpotM"
    assert payload["tile_size_px"] == 224
    assert set(payload["expression"]) == set(module.DEFAULT_GENES)
    assert all(row["unit"] == "log1p-CPM" for row in payload["ranked"])


def test_demo_gene_table_has_one_row_per_gene(demo_output, module):
    lines = (demo_output / "tables" / "gene_expression.csv").read_text(
        encoding="utf-8"
    ).strip().splitlines()
    assert lines[0] == "rank,gene,expression,unit"
    assert len(lines) == len(module.DEFAULT_GENES) + 1


def test_demo_report_is_tagged_demo_and_carries_the_disclaimer(demo_output):
    report = (demo_output / "report.md").read_text(encoding="utf-8")
    assert "(demo)" in report
    assert "not from a model run" in report
    assert DISCLAIMER in report
    assert "log1p-CPM" in report


def test_demo_environment_pins_conda_forge_and_nodefaults(demo_output):
    env = (demo_output / "reproducibility" / "environment.yml").read_text(encoding="utf-8")
    assert "conda-forge" in env
    assert "nodefaults" in env
    assert "deepspotm>=1.0" in env


def test_ranking_is_descending(demo_output):
    payload = json.loads((demo_output / "result.json").read_text(encoding="utf-8"))
    values = [row["expression"] for row in payload["ranked"]]
    assert values == sorted(values, reverse=True)


# ---------------------------------------------------------------------------
# Error paths
# ---------------------------------------------------------------------------


def test_a_wrong_sized_tile_exits_with_a_readable_message(tmp_path):
    from PIL import Image

    tile = tmp_path / "too_big.png"
    Image.new("RGB", (256, 256), (220, 170, 200)).save(tile)

    result = subprocess.run(
        [sys.executable, str(MODULE_PATH), "--input", str(tile), "--output", str(tmp_path / "out")],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "224x224" in result.stderr
    assert "256x256" in result.stderr


def test_a_missing_model_stack_points_at_requirements_and_demo(module, monkeypatch):
    monkeypatch.setitem(sys.modules, "deepspotm", None)
    with pytest.raises(RuntimeError) as exc:
        module.predict_expression(object(), ["EPCAM"], "scgpt")
    message = str(exc.value)
    assert "skills/deepspot-m/requirements.txt" in message
    assert "huggingface-cli login" in message
    assert "--demo" in message


# ---------------------------------------------------------------------------
# Lazy imports
# ---------------------------------------------------------------------------


def test_importing_the_module_does_not_import_the_model_stack():
    """deepspotm and torch load inside predict_expression, never at import."""
    probe = (
        "import importlib.util, sys;"
        f"spec = importlib.util.spec_from_file_location('deepspot_m', r'{MODULE_PATH}');"
        "module = importlib.util.module_from_spec(spec);"
        "spec.loader.exec_module(module);"
        "print(sorted(n for n in ('deepspotm', 'torch') if n in sys.modules))"
    )
    result = subprocess.run(
        [sys.executable, "-c", probe],
        capture_output=True,
        text=True,
        check=True,
    )
    assert result.stdout.strip() == "[]"
