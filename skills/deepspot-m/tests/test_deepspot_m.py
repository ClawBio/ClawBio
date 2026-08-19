"""Tests for the deepspot-m virtual spatial transcriptomics skill.

Every test here runs without torch, without deepspotm and without the gated
model weights, so the suite is green on a base ClawBio checkout. The live
prediction path is covered by injecting a stub `deepspotm`, `torch` and
`huggingface_hub` into sys.modules, so the gene resolution, the pinned revision,
the download gate and the tensor handling are all exercised even though CI has
no weights.
"""
from __future__ import annotations

import importlib.util
import json
import subprocess
import sys
import types
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
# A stub model stack, so the live path is exercised without weights
# ---------------------------------------------------------------------------


class _FakeTensor:
    """Minimal stand-in for the (1, n_genes) tensor predict_genes returns."""

    def __init__(self, values):
        self._values = list(values)

    def squeeze(self, _dim):
        return self

    def tolist(self):
        return list(self._values)


class _FakeBatch:
    def __init__(self, tile):
        self.tile = tile

    def unsqueeze(self, _dim):
        return self


class _FakeModel:
    """Mimics the parts of DeepSpotM the skill actually touches."""

    # Deliberately mixed case, including an HGNC symbol with lower-case 'orf'.
    gene_names = ["EPCAM", "COL1A1", "PTPRC", "C9orf72", "MKI67", "TP53"]

    def __init__(self):
        self.eval_called = False
        self.seen_genes = None

    def eval(self):
        self.eval_called = True
        return self

    def predict_genes(self, batch, genes):
        # Upstream raises KeyError for anything outside the panel.
        missing = [g for g in genes if g not in self.gene_names]
        if missing:
            raise KeyError(f"{len(missing)} gene(s) not in the model panel: {missing}")
        self.seen_genes = list(genes)
        return _FakeTensor(float(self.gene_names.index(g)) for g in genes)


def install_stub_stack(monkeypatch, calls: dict, snapshot="/hf-cache/snapshots/pinned"):
    """Put a fake `deepspotm`, `torch` and `huggingface_hub` on sys.modules.

    The fake `hf_hub_download` records the arguments the gate depends on, so the
    tests can check that the skill asked the library to stay local rather than
    that it set an environment variable the library had already read.
    """
    model = _FakeModel()

    class _FakeDeepSpotM:
        @staticmethod
        def from_pretrained(repo, *, source, revision=None, **kwargs):
            calls["repo"] = repo
            calls["source"] = source
            calls["revision"] = revision
            return model, _FakeBatch

    deepspotm = types.ModuleType("deepspotm")
    deepspotm.DeepSpotM = _FakeDeepSpotM

    def _fake_hf_hub_download(repo_id, filename, *, revision=None,
                              local_files_only=False, **kwargs):
        calls["download_repo"] = repo_id
        calls["download_revision"] = revision
        calls["local_files_only"] = local_files_only
        calls.setdefault("downloaded", []).append(filename)
        return f"{snapshot}/{filename}"

    huggingface_hub = types.ModuleType("huggingface_hub")
    huggingface_hub.hf_hub_download = _fake_hf_hub_download
    monkeypatch.setitem(sys.modules, "huggingface_hub", huggingface_hub)

    class _NoGrad:
        def __enter__(self):
            return self

        def __exit__(self, *exc):
            return False

    torch = types.ModuleType("torch")
    torch.no_grad = _NoGrad

    monkeypatch.setitem(sys.modules, "deepspotm", deepspotm)
    monkeypatch.setitem(sys.modules, "torch", torch)
    calls["model"] = model
    return model


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
# Gene symbols: case-fold to look up, emit what the panel calls it
# ---------------------------------------------------------------------------


def test_parse_genes_deduplicates_without_destroying_case(module):
    """HGNC keeps 'orf' lower case, so parse_genes must not upper-case."""
    assert module.parse_genes(" c9orf72 , CD37 ,C9ORF72, COL1A1 ") == [
        "c9orf72",
        "CD37",
        "COL1A1",
    ]


def test_parse_genes_defaults_to_the_bundled_panel(module):
    assert module.parse_genes(None) == list(module.DEFAULT_GENES)


def test_parse_genes_rejects_an_empty_list(module):
    with pytest.raises(ValueError):
        module.parse_genes(" , , ")


def test_resolve_genes_returns_the_panel_spelling(module):
    panel = ["EPCAM", "C9orf72", "COL1A1"]
    assert module.resolve_genes(["c9orf72", "epcam"], panel) == ["C9orf72", "EPCAM"]


def test_resolve_genes_reports_off_panel_symbols(module):
    with pytest.raises(ValueError) as exc:
        module.resolve_genes(["EPCAM", "NOTAGENE"], ["EPCAM", "COL1A1"])
    message = str(exc.value)
    assert "NOTAGENE" in message
    assert "EPCAM" not in message.split("not in the released panel")[1]


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


def test_demo_mode_resolves_genes_case_insensitively(module):
    assert module.load_demo_expression(["epcam"]) == {"EPCAM": pytest.approx(5.82)}


# ---------------------------------------------------------------------------
# The live prediction path, against a stub model stack
# ---------------------------------------------------------------------------


def test_predict_expression_maps_panel_symbols_onto_scores(module, monkeypatch):
    calls: dict = {}
    model = install_stub_stack(monkeypatch, calls)

    scores = module.predict_expression(object(), ["c9orf72", "EPCAM"], "scgpt")

    # Keys come back as the panel spells them, not as the user typed them.
    assert list(scores) == ["C9orf72", "EPCAM"]
    assert scores == {"C9orf72": 3.0, "EPCAM": 0.0}
    assert model.eval_called is True
    assert model.seen_genes == ["C9orf72", "EPCAM"]


def test_predict_expression_pins_the_weight_revision(module, monkeypatch):
    """The pin now rides on the download, which is the call that can reach out."""
    calls: dict = {}
    install_stub_stack(monkeypatch, calls)

    module.predict_expression(object(), ["EPCAM"], "evo2")

    assert calls["download_repo"] == module.MODEL_REPO
    assert calls["source"] == "evo2"
    assert calls["download_revision"] == module.MODEL_REVISION
    assert len(module.MODEL_REVISION) == 40


def test_predict_expression_stays_offline_unless_download_is_allowed(module, monkeypatch):
    """The gate has to be an argument to the library, not an environment variable.

    huggingface_hub reads HF_HUB_OFFLINE once at import, and `deepspotm` imports
    it before the skill gets a say, so setting the variable at call time proved
    nothing about what the library would do.
    """
    calls: dict = {}
    install_stub_stack(monkeypatch, calls)

    module.predict_expression(object(), ["EPCAM"], "scgpt")
    assert calls["local_files_only"] is True

    module.predict_expression(object(), ["EPCAM"], "scgpt", allow_download=True)
    assert calls["local_files_only"] is False


def test_predict_expression_loads_the_resolved_directory_not_the_repo_id(module, monkeypatch):
    """Handing upstream a repo id would let its own loader fetch, ungated."""
    calls: dict = {}
    install_stub_stack(monkeypatch, calls, snapshot="/hf-cache/snapshots/abc123")

    module.predict_expression(object(), ["EPCAM"], "scgpt")

    assert calls["repo"] == "/hf-cache/snapshots/abc123"
    assert calls["repo"] != module.MODEL_REPO
    # A directory carries no revision, so none is passed on.
    assert calls["revision"] is None
    assert calls["downloaded"] == list(module.MODEL_FILES)


def test_resolve_checkpoint_dir_rejects_a_split_snapshot(module, monkeypatch):
    """The directory is what gets loaded, so it has to hold all three files."""
    huggingface_hub = types.ModuleType("huggingface_hub")
    huggingface_hub.hf_hub_download = (
        lambda repo_id, filename, **kwargs: f"/hf-cache/{filename}/{filename}"
    )
    monkeypatch.setitem(sys.modules, "huggingface_hub", huggingface_hub)

    with pytest.raises(RuntimeError) as exc:
        module._resolve_checkpoint_dir(allow_download=False)
    assert "single directory" in str(exc.value)


def test_predict_expression_explains_a_cold_cache_without_downloading(module, monkeypatch):
    """Offline miss is the common first run; the message has to say what to do."""
    calls: dict = {}
    install_stub_stack(monkeypatch, calls)

    def _cold(repo_id, filename, *, revision=None, local_files_only=False, **kwargs):
        raise OSError("not found in the local cache")

    sys.modules["huggingface_hub"].hf_hub_download = _cold

    with pytest.raises(RuntimeError) as exc:
        module.predict_expression(object(), ["EPCAM"], "scgpt")
    message = str(exc.value)
    assert "--allow-download" in message
    assert module.MODEL_REPO in message


def test_predict_expression_turns_an_off_panel_symbol_into_a_readable_error(module, monkeypatch):
    calls: dict = {}
    install_stub_stack(monkeypatch, calls)

    with pytest.raises(ValueError) as exc:
        module.predict_expression(object(), ["NOTAGENE"], "scgpt")
    assert "NOTAGENE" in str(exc.value)


def test_an_off_panel_gene_exits_cleanly_rather_than_tracebacking(module):
    """KeyError from genes_to_indices must be caught like every other input error."""
    assert KeyError in module.INPUT_ERRORS


def test_a_missing_model_stack_points_at_requirements_and_demo(module, monkeypatch):
    monkeypatch.setitem(sys.modules, "deepspotm", None)
    with pytest.raises(RuntimeError) as exc:
        module.predict_expression(object(), ["EPCAM"], "scgpt")
    message = str(exc.value)
    assert "skills/deepspot-m/requirements.txt" in message
    assert "huggingface-cli login" in message
    assert "--demo" in message


# ---------------------------------------------------------------------------
# Microns per pixel is measured or declared, never assumed
# ---------------------------------------------------------------------------


def test_the_module_asserts_no_default_microns_per_pixel(module):
    """0.5 um/px appears nowhere upstream, so the skill must not hardcode it."""
    assert not hasattr(module, "TARGET_MPP")


def test_pixel_size_is_none_without_resolution_metadata(module):
    tile = module.load_tile(module.DEMO_TILE)
    assert module.read_pixel_size_um(tile) is None


def test_pixel_size_is_read_from_tiff_resolution_tags(module, tmp_path):
    from PIL import Image

    path = tmp_path / "tile.tif"
    # 50000 pixels per inch -> 25400 / 50000 = 0.508 microns per pixel.
    Image.new("RGB", (224, 224), (200, 150, 190)).save(path, dpi=(50000, 50000))

    tile = Image.open(path)
    assert module.read_pixel_size_um(tile) == pytest.approx(0.508, abs=1e-3)


def test_an_implausible_resolution_tag_is_ignored(module, tmp_path):
    from PIL import Image

    path = tmp_path / "screen.tif"
    # 72 dpi is a print/screen resolution, not a microscopy one: 353 um/px.
    Image.new("RGB", (224, 224), (200, 150, 190)).save(path, dpi=(72, 72))

    tile = Image.open(path)
    assert module.read_pixel_size_um(tile) is None


def test_declared_mpp_must_be_positive(module):
    parser = module.build_parser()
    args = parser.parse_args(["--demo", "--output", "/tmp/x", "--mpp", "0.5"])
    assert args.mpp == 0.5
    with pytest.raises(ValueError):
        module.validate_mpp(0.0)
    with pytest.raises(ValueError):
        module.validate_mpp(-1.0)


def test_demo_records_an_undeclared_pixel_size_as_null(demo_output):
    payload = json.loads((demo_output / "result.json").read_text(encoding="utf-8"))
    assert payload["microns_per_pixel"] is None
    assert payload["microns_per_pixel_source"] is None

    report = (demo_output / "report.md").read_text(encoding="utf-8")
    assert "not declared" in report
    assert "0.5 microns per pixel" not in report


def test_a_declared_mpp_wins_over_image_metadata_but_says_so(module, tmp_path):
    from PIL import Image

    path = tmp_path / "tile.tif"
    Image.new("RGB", (224, 224), (200, 150, 190)).save(path, dpi=(50000, 50000))
    tile = Image.open(path)

    resolved, source = module.resolve_pixel_size(tile, declared=0.42)
    assert resolved == 0.42
    assert "--mpp" in source
    assert "0.508" in source  # the metadata it disagreed with

    resolved, source = module.resolve_pixel_size(tile, declared=None)
    assert resolved == pytest.approx(0.508, abs=1e-3)
    assert "resolution tags" in source


# ---------------------------------------------------------------------------
# Tiles that do not look like H&E
# ---------------------------------------------------------------------------


def test_a_stained_tile_raises_no_tile_warnings(module):
    tile = module.load_tile(module.DEMO_TILE)
    assessment = module.assess_tile(tile)
    assert assessment["warnings"] == []
    assert assessment["mean_pixel"] < module.WHITE_MEAN_DEFAULT
    assert assessment["mean_saturation"] > module.MIN_SATURATION_DEFAULT


def test_a_near_white_tile_is_flagged_as_background(module):
    from PIL import Image

    tile = Image.new("RGB", (224, 224), (250, 248, 250))
    warnings = module.assess_tile(tile)["warnings"]
    assert any("background" in w for w in warnings)


def test_a_greyscale_tile_is_flagged_as_not_h_and_e(module):
    from PIL import Image

    tile = Image.new("RGB", (224, 224), (120, 120, 120))
    warnings = module.assess_tile(tile)["warnings"]
    assert any("H&E" in w for w in warnings)


def test_skip_background_refuses_to_score_a_blank_tile(tmp_path):
    from PIL import Image

    tile = tmp_path / "blank.png"
    Image.new("RGB", (224, 224), (250, 248, 250)).save(tile)

    result = subprocess.run(
        [
            sys.executable, str(MODULE_PATH),
            "--input", str(tile),
            "--skip-background",
            "--output", str(tmp_path / "out"),
        ],
        capture_output=True,
        text=True,
    )
    assert result.returncode == 1
    assert "background" in result.stderr


def test_tile_warnings_travel_into_the_report(module, tmp_path):
    """A warning is useless if it only ever reaches the terminal."""
    from PIL import Image

    tile = Image.new("RGB", (224, 224), (250, 248, 250))
    rows = [{"gene": "EPCAM", "expression": 5.0, "unit": "log1p-CPM", "rank": 1}]
    meta = module.build_meta(
        demo=False, tile_label="blank.png", source="scgpt", genes=["EPCAM"],
        mpp=None, mpp_source=None, assessment=module.assess_tile(tile),
    )
    out = tmp_path / "out"
    out.mkdir()
    report = module.write_report(out, rows, meta).read_text(encoding="utf-8")
    assert "background" in report

    payload = json.loads(module.write_result_json(out, rows, meta).read_text(encoding="utf-8"))
    assert any("background" in w for w in payload["tile_warnings"])


def test_replay_command_carries_the_thresholds_that_decide_scoring(module, tmp_path):
    """report.md calls this the exact command, so a raised threshold has to survive.

    --white-mean and --min-saturation gate whether a tile is scored at all. A
    bundle that omitted a raised --white-mean would refuse to score the very
    tile the original run scored.
    """
    out = tmp_path / "out"
    out.mkdir()
    table = out / "genes.csv"
    table.write_text("gene,expression\n", encoding="utf-8")

    parser = module.build_parser()
    args = parser.parse_args([
        "--input", "tile.png",
        "--white-mean", "245",
        "--min-saturation", "0.01",
        "--output", str(out),
    ])
    meta = module.build_meta(
        demo=False, tile_label="tile.png", source="scgpt", genes=["EPCAM"],
        mpp=None, mpp_source=None, assessment={"mean_pixel": 200.0,
                                               "mean_saturation": 0.3,
                                               "warnings": []},
    )
    module.write_reproducibility(out, table, meta, str(tmp_path / "tile.png"), args)

    # The bundle puts one argument per line, and argparse parsed these as floats.
    tokens = (out / "reproducibility" / "commands.sh").read_text(encoding="utf-8").split()
    assert tokens[tokens.index("--white-mean") + 2] == "245.0"
    assert tokens[tokens.index("--min-saturation") + 2] == "0.01"


def test_replay_command_stays_quiet_about_default_thresholds(module, tmp_path):
    """Defaults are already the default; spelling them out is just noise."""
    out = tmp_path / "out"
    out.mkdir()
    table = out / "genes.csv"
    table.write_text("gene,expression\n", encoding="utf-8")

    parser = module.build_parser()
    args = parser.parse_args(["--input", "tile.png", "--output", str(out)])
    meta = module.build_meta(
        demo=False, tile_label="tile.png", source="scgpt", genes=["EPCAM"],
        mpp=None, mpp_source=None, assessment={"mean_pixel": 200.0,
                                               "mean_saturation": 0.3,
                                               "warnings": []},
    )
    module.write_reproducibility(out, table, meta, str(tmp_path / "tile.png"), args)

    command = (out / "reproducibility" / "commands.sh").read_text(encoding="utf-8")
    assert "--white-mean" not in command
    assert "--min-saturation" not in command


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
    assert payload["model_revision"] == module.MODEL_REVISION
    assert payload["tile_size_px"] == 224
    assert set(payload["expression"]) == set(module.DEFAULT_GENES)


def test_demo_gene_table_keeps_rank_as_a_secondary_column(demo_output, module):
    lines = (demo_output / "tables" / "gene_expression.csv").read_text(
        encoding="utf-8"
    ).strip().splitlines()
    assert lines[0] == "gene,expression,unit,rank"
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
    assert "deepspotm>=1.0,<2" in env


# ---------------------------------------------------------------------------
# What the report claims
# ---------------------------------------------------------------------------


def test_values_are_reported_in_the_order_they_were_requested(module, tmp_path):
    """No cross-gene ranking in the headline: the model predicts relative
    expression, so ordering genes within one tile mostly recovers each gene's
    training-set mean."""
    out = tmp_path / "ordered"
    module.main(["--demo", "--genes", "MKI67,EPCAM,CD68", "--output", str(out)])

    payload = json.loads((out / "result.json").read_text(encoding="utf-8"))
    assert payload["genes"] == ["MKI67", "EPCAM", "CD68"]
    assert list(payload["expression"]) == ["MKI67", "EPCAM", "CD68"]

    table = (out / "report.md").read_text(encoding="utf-8").split("| Gene | Expression")[1]
    assert table.index("MKI67") < table.index("EPCAM") < table.index("CD68")


def test_the_report_no_longer_headlines_the_highest_gene(demo_output):
    report = (demo_output / "report.md").read_text(encoding="utf-8")
    assert "holds the highest value" not in report
    assert "| Rank |" not in report


def test_the_report_explains_that_cross_gene_ordering_tracks_average_abundance(demo_output):
    report = (demo_output / "report.md").read_text(encoding="utf-8")
    assert "## How to Read These Values" in report
    assert "average abundance" in report
    assert "relative expression rather than absolute counts" in report


@pytest.mark.parametrize(
    "limitation",
    [
        "Trained on a finite set of cancer indications.",
        "Performance on unseen tissue types, stains, scanners or resolutions may degrade.",
        "Predicts relative expression rather than absolute counts.",
        "Under-sequenced genes are predicted less reliably.",
        "Trained on oncology cohorts, so it is not representative of healthy tissue "
        "or non-oncology contexts.",
    ],
)
def test_every_upstream_limitation_reaches_the_report(demo_output, limitation):
    report = (demo_output / "report.md").read_text(encoding="utf-8")
    assert limitation in report


@pytest.mark.parametrize(
    "limitation",
    [
        "Trained on a finite set of cancer indications.",
        "Performance on unseen tissue types, stains, scanners or resolutions may degrade.",
        "Predicts relative expression rather than absolute counts.",
        "Under-sequenced genes are predicted less reliably.",
        "Trained on oncology cohorts, so it is not representative of healthy tissue "
        "or non-oncology contexts.",
    ],
)
def test_every_upstream_limitation_reaches_skill_md(limitation):
    skill_md = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert limitation in skill_md


def test_skill_md_does_not_assert_an_unsourced_pixel_size():
    skill_md = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert "0.5 microns per pixel" not in skill_md
    assert "(source: DeepSpot-M model card)" not in skill_md


def test_skill_md_example_output_is_a_real_demo_run():
    """The block used to show fixture numbers dressed up as a model run."""
    skill_md = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    block = skill_md.split("## Example Output")[1].split("\n## ")[0]
    assert "(demo)" in block
    assert "Fixture Expression" in block
    assert "**Tile**: tile.png" not in block


def test_skill_md_declares_the_upstream_authorship():
    skill_md = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert "Kalin Nonchev" in skill_md


# ---------------------------------------------------------------------------
# The non-commercial term travels with the numbers
# ---------------------------------------------------------------------------


def test_a_real_run_stamps_the_weights_licence_on_its_outputs(module, tmp_path):
    """Upstream puts NonCommercial on 'the weights or their outputs', so a
    forwarded report has to carry the restriction with it."""
    out = tmp_path / "real"
    out.mkdir()
    meta = module.build_meta(
        demo=False, tile_label="tile.png", source="scgpt", genes=["EPCAM"],
        mpp=None, mpp_source=None,
        assessment={"mean_pixel": 180.0, "mean_saturation": 0.3, "warnings": []},
    )
    rows = [{"gene": "EPCAM", "expression": 5.0, "unit": "log1p-CPM", "rank": 1}]

    report = module.write_report(out, rows, meta).read_text(encoding="utf-8")
    payload = json.loads(module.write_result_json(out, rows, meta).read_text(encoding="utf-8"))

    assert "CC-BY-NC-SA-4.0" in report
    assert "non-commercial" in report
    assert payload["weights_license"] == "CC-BY-NC-SA-4.0"
    assert "10.64898" in payload["output_license_note"]


def test_the_demo_run_claims_no_weights_licence(demo_output):
    """Fixture values never touched the weights, so the stamp would be a lie."""
    payload = json.loads((demo_output / "result.json").read_text(encoding="utf-8"))
    assert "weights_license" not in payload
    assert "CC-BY-NC-SA-4.0" not in (demo_output / "report.md").read_text(encoding="utf-8")


def test_skill_md_states_that_the_restriction_covers_outputs():
    skill_md = (SKILL_DIR / "SKILL.md").read_text(encoding="utf-8")
    assert "the weights or their outputs" in skill_md


# ---------------------------------------------------------------------------
# Input paths do not ride along into shareable artefacts
# ---------------------------------------------------------------------------


def test_the_report_records_the_tile_name_not_its_directory(module, tmp_path):
    from PIL import Image

    secret = tmp_path / "PATIENT_12345"
    secret.mkdir()
    tile = secret / "slide.png"
    Image.new("RGB", (224, 224), (200, 150, 190)).save(tile)

    out = tmp_path / "out"
    out.mkdir()
    meta = module.build_meta(
        demo=False, tile_label=module.tile_label(str(tile)), source="scgpt",
        genes=["EPCAM"], mpp=None, mpp_source=None,
        assessment=module.assess_tile(module.load_tile(tile)),
    )
    rows = [{"gene": "EPCAM", "expression": 5.0, "unit": "log1p-CPM", "rank": 1}]

    report = module.write_report(out, rows, meta).read_text(encoding="utf-8")
    payload = module.write_result_json(out, rows, meta).read_text(encoding="utf-8")

    assert "PATIENT_12345" not in report
    assert "PATIENT_12345" not in payload
    assert "slide.png" in report


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
