"""The figure names the exposure's QTL kind from ONE rule.

Before this rule every eQTL Catalogue dataset was titled "Exposure (eQTL)",
including aptamer (protein), leafcutter (splicing) and txrev (transcript usage)
datasets, and both scatter x-axes said "eQTL / exposure" even for a UKB-PPP
protein QTL. A pre-fetched file was titled "... | ge" and recorded as
`exposure_source: eqtl_catalogue` although nothing said it came from there.

Each test here drives the real render path (cassette replay or the offline
example) and reads what was actually drawn or written.
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import matplotlib

matplotlib.use("Agg")

import pytest
import yaml
from matplotlib.figure import Figure

SKILL_DIR = Path(__file__).resolve().parent.parent
TESTS_DIR = Path(__file__).resolve().parent
for _p in (SKILL_DIR, TESTS_DIR):
    if str(_p) not in sys.path:
        sys.path.insert(0, str(_p))

import cli  # noqa: E402
import locuscompare_region_render as composer  # noqa: E402  (puts the sibling skills on sys.path)
import regional_plot  # noqa: E402
from test_golden_parity_sort1_dual_modality import (  # noqa: E402
    _CassetteEQTLClient,
    _CassetteGWASClient,
    _CassetteUKBPPPClient,
    _eqtl_spec,
    _load_eqtl_cassette,
    _load_gwas_cassette,
    _load_pqtl_cassette,
    _pqtl_spec,
)

EXAMPLE_01 = SKILL_DIR / "examples" / "01_synthetic_demo"


@pytest.fixture
def drawn(monkeypatch):
    """Record the renderer input and both scatter x-axis labels on the real path."""
    seen: dict = {}
    real_render = composer.render_full_locuscompare
    real_lc = regional_plot._render_locuscompare_scatter
    real_es = regional_plot._render_effect_size_scatter
    real_gene_track = regional_plot.render_gene_track

    def render_spy(inp, out_path):
        seen["inp"] = inp
        return real_render(inp, out_path)

    def lc_spy(ax, *a, **kw):
        out = real_lc(ax, *a, **kw)
        seen["lc_xlabel"] = ax.get_xlabel()
        return out

    def es_spy(ax, *a, **kw):
        out = real_es(ax, *a, **kw)
        seen["es_xlabel"] = ax.get_xlabel()
        return out

    def gene_track_spy(ax, *a, **kw):
        seen["gene_track_label"] = kw.get("track_label")
        return real_gene_track(ax, *a, **kw)

    monkeypatch.setattr(composer, "render_full_locuscompare", render_spy)
    monkeypatch.setattr(regional_plot, "_render_locuscompare_scatter", lc_spy)
    monkeypatch.setattr(regional_plot, "_render_effect_size_scatter", es_spy)
    monkeypatch.setattr(regional_plot, "render_gene_track", gene_track_spy)
    return seen


# quant_method -> (modality code, tag)
_CASES = {
    "ge": ("eqtl", "eQTL"),
    "microarray": ("eqtl", "eQTL"),
    "exon": ("exon_qtl", "exon eQTL"),
    "tx": ("tx_qtl", "transcript eQTL"),
    "txrev": ("tuqtl", "tuQTL"),
    "leafcutter": ("sqtl", "sQTL"),
    "aptamer": ("pqtl", "pQTL"),
    None: ("qtl", "QTL"),
    "zzz": ("qtl", "QTL"),
}


def test_exposure_modality_rule():
    exposure_modality = composer.exposure_modality
    labels = composer.EXPOSURE_MODALITY_LABELS
    eq = composer.EXPOSURE_KIND_EQTL_CATALOGUE
    ukb = composer.EXPOSURE_KIND_UKB_PPP

    for qm, (code, tag) in _CASES.items():
        assert exposure_modality(eq, qm) == (code, tag, labels[code]), qm
    # UKB-PPP is a protein QTL whatever quantification code rides along.
    assert exposure_modality(ukb, None)[:2] == ("pqtl", "pQTL")
    # The catalogue spells it "aFC"; matching is case-insensitive, so it must map.
    assert exposure_modality(eq, "aFC")[:2] == ("eqtl", "eQTL")
    assert exposure_modality(ukb, "ge")[:2] == ("pqtl", "pQTL")
    # The single-cell hook turns an expression QTL into a single-cell one, and
    # nothing else.
    assert exposure_modality(eq, "ge", "single-cell")[:2] == ("sceqtl", "sceQTL")
    assert exposure_modality(eq, "ge", "single_cell")[:2] == ("sceqtl", "sceQTL")
    assert exposure_modality(eq, "ge", "bulk")[:2] == ("eqtl", "eQTL")
    assert exposure_modality(eq, "leafcutter", "single-cell")[:2] == ("sqtl", "sQTL")


@pytest.mark.parametrize("quant_method", list(_CASES), ids=lambda q: str(q))
def test_cassette_replay_names_the_modality_everywhere(tmp_path, drawn, quant_method):
    code, tag = _CASES[quant_method]
    cassette = _load_eqtl_cassette()
    cassette.release.quant_method = quant_method
    result = composer.render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=_CassetteEQTLClient(cassette),  # type: ignore[arg-type]
        gwas_client=_CassetteGWASClient(_load_gwas_cassette()),  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "lc.png",
    )
    title = drawn["inp"].exposure_short_label
    assert title.startswith(f"Exposure ({tag}): "), title
    assert drawn["lc_xlabel"] == f"−log10 p ({tag} / exposure)"
    assert drawn["es_xlabel"] == f"β ({tag} / exposure)"

    block = result.manifest_block
    assert block["exposure_modality"] == code
    assert block["exposure_modality_label"]
    assert block["exposure_quant_method"] == quant_method
    if quant_method == "ge":
        assert "gene expression (ge)" in title
        assert block["exposure_quant_method_label"] == "gene expression"
    elif quant_method is None:
        assert block["exposure_quant_method_label"] is None
        assert "| None" not in title
    elif quant_method == "zzz":
        # An unknown code is kept as itself, once.
        assert "| zzz (QTD000276)" in title
        assert block["exposure_quant_method_label"] == "zzz"


def test_ukb_ppp_replay_is_a_pqtl_on_every_label(tmp_path, drawn):
    result = composer.render_locuscompare_for_lead(
        spec=_pqtl_spec(),
        eqtl_client=_CassetteEQTLClient(_load_eqtl_cassette()),  # type: ignore[arg-type]
        gwas_client=_CassetteGWASClient(_load_gwas_cassette()),  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "pqtl.png",
        ukb_ppp_client=_CassetteUKBPPPClient(_load_pqtl_cassette()),  # type: ignore[arg-type]
    )
    assert drawn["inp"].exposure_short_label.startswith("Exposure (pQTL): ")
    assert drawn["lc_xlabel"] == "−log10 p (pQTL / exposure)"
    assert drawn["es_xlabel"] == "β (pQTL / exposure)"
    block = result.manifest_block
    assert block["exposure_modality"] == "pqtl"
    assert block["exposure_modality_label"] == "protein QTL"
    assert block["exposure_quant_method"] is None


@pytest.mark.xfail(
    strict=True,
    reason=(
        "the eQTL Catalogue fetcher's release record does not carry the study "
        "type (single-cell vs bulk) yet, so a single-cell study is still titled "
        "as bulk expression"
    ),
)
def test_onek1k_single_cell_dataset_is_titled_sceqtl(tmp_path, drawn):
    cassette = _load_eqtl_cassette()
    cassette.release.study_label = "OneK1K"
    cassette.release.sample_group = "CD14_monocyte"
    cassette.release.quant_method = "ge"
    composer.render_locuscompare_for_lead(
        spec=_eqtl_spec(),
        eqtl_client=_CassetteEQTLClient(cassette),  # type: ignore[arg-type]
        gwas_client=_CassetteGWASClient(_load_gwas_cassette()),  # type: ignore[arg-type]
        ld_client=None,
        out_path=tmp_path / "onek1k.png",
    )
    assert drawn["inp"].exposure_short_label.startswith("Exposure (sceQTL): ")


def _run_example_01(
    tmp_path: Path, monkeypatch, *, extra_exposure: dict | None = None, expect_rc: int = 0,
):
    config = json.loads((EXAMPLE_01 / "config.json").read_text())
    for side in ("exposure", "outcome"):
        config[side]["sumstats_path"] = str(EXAMPLE_01 / config[side]["sumstats_path"])
    config["ld"]["ld_matrix_path"] = str(EXAMPLE_01 / config["ld"]["ld_matrix_path"])
    config["gene_track"]["genes_path"] = str(EXAMPLE_01 / config["gene_track"]["genes_path"])
    config["exposure"].update(extra_exposure or {})
    cfg = tmp_path / "config.json"
    cfg.write_text(json.dumps(config))

    texts: list[str] = []
    real_savefig = Figure.savefig

    def savefig_spy(self, *args, **kwargs):
        texts.extend(t.get_text() for t in self.texts)
        return real_savefig(self, *args, **kwargs)

    monkeypatch.setattr(Figure, "savefig", savefig_spy)
    out = tmp_path / "out"
    rc = cli.main(["--input", str(cfg), "--output", str(out)])
    if expect_rc != 0:
        assert rc == expect_rc
        return None, ""
    assert rc == 0
    manifest = yaml.safe_load((out / "manifest.yaml").read_text())
    return manifest["render_block"], "\n".join(texts)


def test_example_01_prefetched_file_claims_no_source_or_quantification(
    tmp_path, monkeypatch, drawn,
):
    block, fig_text = _run_example_01(tmp_path, monkeypatch)
    title = drawn["inp"].exposure_short_label
    assert title.startswith("Exposure (QTL): ")
    assert "| ge" not in title
    assert drawn["lc_xlabel"] == "−log10 p (QTL / exposure)"
    assert block["exposure_source"] == "prefetched"
    assert block["outcome_source"] == "prefetched"
    assert block["exposure_modality"] == "qtl"
    assert block["exposure_quant_method"] is None
    warnings = block["data_source_warnings"]
    assert warnings and all(w.startswith("prefetched") for w in warnings), warnings
    assert "Exposure: pre-fetched file" in fig_text
    assert "Outcome: pre-fetched file" in fig_text
    assert "eQTL Catalogue" not in fig_text
    assert "GWAS Catalog" not in fig_text
    # Short gene-track label; the caption names the gene source actually used,
    # here the example's made-up TSV, never a GENCODE release it did not come from.
    assert drawn["gene_track_label"] == "Genes"
    assert "Genes: synthetic TSV genes.tsv" in fig_text
    assert "GENCODE" not in fig_text


def test_example_01_with_declared_aptamer_is_a_pqtl(tmp_path, monkeypatch, drawn):
    block, fig_text = _run_example_01(
        tmp_path, monkeypatch, extra_exposure={"quant_method": "aptamer"},
    )
    # A declared quantification names the QTL kind and nothing else: no
    # credible-set (.cc) caveat, which is a claim about eQTL Catalogue files.
    assert "eQTL Catalogue" not in fig_text
    assert "eQTL Catalogue" not in yaml.safe_dump(block)
    assert drawn["inp"].exposure_short_label.startswith("Exposure (pQTL): ")
    assert drawn["es_xlabel"] == "β (pQTL / exposure)"
    assert block["exposure_modality"] == "pqtl"
    assert block["exposure_quant_method"] == "aptamer"
    assert block["exposure_source"] == "prefetched"


def test_example_01_declared_quant_method_ignores_case(tmp_path, monkeypatch, drawn):
    block, _ = _run_example_01(tmp_path, monkeypatch, extra_exposure={"quant_method": "GE"})
    assert drawn["inp"].exposure_short_label.startswith("Exposure (eQTL): ")
    assert block["exposure_quant_method"] == "ge"


def test_example_01_unknown_quant_method_is_refused_not_relabelled(tmp_path, monkeypatch, capsys):
    """A misspelt code used to render a 'QTL' figure with no word about why."""
    _run_example_01(
        tmp_path, monkeypatch, extra_exposure={"quant_method": "gene_expr"}, expect_rc=2,
    )
    err = capsys.readouterr().err
    assert "'gene_expr' is not a known quantification code" in err
    assert not (tmp_path / "out" / "manifest.yaml").exists()
