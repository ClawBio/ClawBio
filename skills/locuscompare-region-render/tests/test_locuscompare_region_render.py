"""Smoke tests for the locuscompare CLI.

Verifies the CLI contract (argument parsing, config loading, output
directory creation, helpful error messages) without exercising the full
fetch / harmonise / render pipeline. End-to-end tests with real network
calls live in test_live_locuscompare_region_render.py (marked @pytest.mark.live).
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import pytest
import yaml


SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import cli as locuscompare  # noqa: E402


def test_cli_requires_input_or_demo(tmp_path: Path):
    """argparse should reject invocations missing both --input and --demo."""
    with pytest.raises(SystemExit) as excinfo:
        locuscompare.main(["--output", str(tmp_path)])
    assert excinfo.value.code == 2


def test_cli_creates_output_directory_on_invalid_config(tmp_path: Path):
    """--output is created up-front, even when the config later fails validation."""
    out = tmp_path / "sub1" / "sub2"
    bad_config = {"schema_version": "1.0"}  # missing lead/exposure/outcome
    config_path = tmp_path / "config.json"
    config_path.write_text(json.dumps(bad_config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(out)])
    assert out.is_dir()
    assert rc == 2


def test_cli_rejects_missing_lead(tmp_path: Path, capsys):
    config = {
        "exposure": {"trait_label": "x", "fetch": {"source": "eqtl_catalogue", "dataset_id": "X"}},
        "outcome": {"trait_label": "y", "fetch": {"source": "gwas_catalog", "accession": "Y"}},
    }
    config_path = tmp_path / "c.json"
    config_path.write_text(json.dumps(config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(tmp_path / "out")])
    assert rc == 2
    captured = capsys.readouterr()
    assert "lead" in captured.err.lower()


def test_cli_rejects_unsupported_exposure_source(tmp_path: Path, capsys):
    config = {
        "lead": {"variant_id": "1_1_A_T", "chromosome": "1", "position_bp": 1, "window_bp": 1000000},
        "exposure": {"trait_label": "x", "fetch": {"source": "made_up_source", "dataset_id": "X"}},
        "outcome": {"trait_label": "y", "fetch": {"source": "gwas_catalog", "accession": "Y"}},
    }
    config_path = tmp_path / "c.json"
    config_path.write_text(json.dumps(config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(tmp_path / "out")])
    assert rc == 2
    captured = capsys.readouterr()
    assert "made_up_source" in captured.err


def test_cli_rejects_missing_exposure_input_vector(tmp_path: Path, capsys):
    """Exposure must declare either `fetch:` (live tabix) or `sumstats_path:`
    (pre-fetched canonical TSV); missing both is a hard error."""
    config = {
        "lead": {"variant_id": "1_1_A_T", "chromosome": "1", "position_bp": 1, "window_bp": 1000000},
        "exposure": {"trait_label": "x"},  # no fetch and no sumstats_path
        "outcome": {"trait_label": "y", "fetch": {"source": "gwas_catalog", "accession": "Y"}},
    }
    config_path = tmp_path / "c.json"
    config_path.write_text(json.dumps(config))
    rc = locuscompare.main(["--input", str(config_path), "--output", str(tmp_path / "out")])
    assert rc == 2
    captured = capsys.readouterr()
    assert "sumstats_path" in captured.err and "fetch" in captured.err


def test_load_config_yaml(tmp_path: Path):
    """_load_config dispatches on file extension."""
    yaml_path = tmp_path / "c.yaml"
    yaml_path.write_text("schema_version: '1.0'\nfoo: bar\n")
    cfg = locuscompare._load_config(yaml_path)
    assert cfg == {"schema_version": "1.0", "foo": "bar"}


def test_load_config_json(tmp_path: Path):
    json_path = tmp_path / "c.json"
    json_path.write_text('{"schema_version": "1.0", "foo": "bar"}')
    cfg = locuscompare._load_config(json_path)
    assert cfg == {"schema_version": "1.0", "foo": "bar"}


def test_load_config_unknown_extension(tmp_path: Path):
    bad = tmp_path / "c.toml"
    bad.write_text("[section]\nfoo = 'bar'\n")
    with pytest.raises(ValueError, match="unsupported config extension"):
        locuscompare._load_config(bad)


def test_format_provenance_prefix_empty():
    assert locuscompare._format_provenance_prefix({}) == ""


def test_format_provenance_prefix_ot_release():
    assert locuscompare._format_provenance_prefix({"provenance": {"ot_release": "26.03"}}) == "OT release: 26.03 | "


def test_format_provenance_prefix_gwas_lookup():
    assert locuscompare._format_provenance_prefix(
        {"provenance": {"gwas_lookup_run_dir": "runs/gl/"}}
    ) == "gwas-lookup chain: runs/gl/ | "


# ----- Pre-fetched TSV input mode (INPUT_SCHEMA.md path) -----


CANONICAL_TSV_HEADER = (
    "variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\tp"
)


def _write_canonical_tsv(path: Path, rows: list[str]) -> None:
    path.write_text(CANONICAL_TSV_HEADER + "\n" + "\n".join(rows) + "\n")


def test_prefetched_load_sumstats_tsv_parses_canonical_schema(tmp_path: Path):
    from _prefetched import load_sumstats_tsv

    tsv = tmp_path / "exposure.tsv"
    _write_canonical_tsv(tsv, [
        "1_500000_A_T\t1\t500000\tA\tT\t0.5\t0.05\t1e-22",
        "1_500100_C_G\t1\t500100\tC\tG\t0.3\t0.05\t1e-9",
    ])
    variants = load_sumstats_tsv(tsv)
    assert len(variants) == 2
    assert variants[0].variant_id == "1_500000_A_T"
    assert variants[0].ref == "A" and variants[0].alt == "T"
    assert variants[0].beta == 0.5 and variants[0].se == 0.05
    assert variants[0].p_value == 1e-22


def test_prefetched_load_sumstats_tsv_rejects_missing_columns(tmp_path: Path):
    from _prefetched import PrefetchedSchemaError, load_sumstats_tsv

    tsv = tmp_path / "bad.tsv"
    # Missing the required `p` column.
    tsv.write_text("variant_id\tchromosome\tposition_bp\tallele_a\tallele_b\tbeta\tse\n")
    with pytest.raises(PrefetchedSchemaError, match="missing required columns"):
        load_sumstats_tsv(tsv)


def test_prefetched_load_sumstats_tsv_drops_rows_with_na(tmp_path: Path):
    from _prefetched import load_sumstats_tsv

    tsv = tmp_path / "exposure.tsv"
    _write_canonical_tsv(tsv, [
        "1_500000_A_T\t1\t500000\tA\tT\t0.5\t0.05\t1e-22",
        "1_500100_C_G\t1\t500100\tC\tG\tNA\tNA\tNA",  # dropped
    ])
    variants = load_sumstats_tsv(tsv)
    assert len(variants) == 1
    assert variants[0].variant_id == "1_500000_A_T"


def test_prefetched_load_synthetic_ld_parses_two_columns(tmp_path: Path):
    from _prefetched import load_synthetic_ld

    tsv = tmp_path / "ld.tsv"
    tsv.write_text("partner_variant_id\tr2\n1_500100_C_G\t0.85\n1_499900_G_A\t0.42\n")
    r2 = load_synthetic_ld(tsv)
    assert r2 == {"1_500100_C_G": 0.85, "1_499900_G_A": 0.42}


def test_prefetched_load_synthetic_ld_rejects_out_of_range(tmp_path: Path):
    from _prefetched import PrefetchedSchemaError, load_synthetic_ld

    tsv = tmp_path / "ld.tsv"
    tsv.write_text("partner_variant_id\tr2\n1_500100_C_G\t1.5\n")
    with pytest.raises(PrefetchedSchemaError, match="out of"):
        load_synthetic_ld(tsv)


def test_prefetched_load_synthetic_gene_track_parses(tmp_path: Path):
    from _prefetched import load_synthetic_gene_track

    tsv = tmp_path / "genes.tsv"
    tsv.write_text(
        "gene_symbol\tstart\tend\tstrand\tbiotype\n"
        "DEMOGENE_A\t100\t200\t+\tprotein_coding\n"
        "DEMOGENE_B\t450\t550\t-\tprotein_coding\n"
    )
    genes = load_synthetic_gene_track(tsv)
    assert [g.gene_symbol for g in genes] == ["DEMOGENE_A", "DEMOGENE_B"]
    assert genes[0].start == 100 and genes[0].end == 200 and genes[0].strand == "+"
    assert genes[1].strand == "-"


def test_prefetched_synthetic_demo_runs_end_to_end_offline(tmp_path: Path):
    """01_synthetic_demo must run without network: TSV-loaded sumstats +
    synthetic LD matrix + synthetic gene track. Asserts a PNG + manifest land
    in the output directory and the expected lead emerges from harmonisation."""
    out = tmp_path / "out"
    demo_config = SKILL_DIR / "examples" / "01_synthetic_demo" / "config.json"
    rc = locuscompare.main(["--input", str(demo_config), "--output", str(out)])
    assert rc == 0, "01_synthetic_demo failed offline"
    plot = out / "1_500000_A_T_full_locuscompare.png"
    manifest = out / "manifest.yaml"
    assert plot.is_file() and plot.stat().st_size > 10_000
    assert manifest.is_file()
    manifest_text = manifest.read_text()
    assert "1_500000_A_T" in manifest_text
    # window_bp is the full width (1,000,000 here); report.md reports the
    # half-window the fetch actually covered, +/-500 kb, never +/-1000 kb. The
    # figure title is asserted in test_regional_plot.py and the caption's
    # `Window:` line in test_composer_passes_the_half_window_label_to_the_renderer.
    report = (out / "report.md").read_text()
    assert "±500 kb" in report and "1000 kb" not in report


def test_composer_passes_the_half_window_label_to_the_renderer(tmp_path: Path, monkeypatch):
    """The caption's `Window:` line is the composer's window_label, so it is
    checked where it is built (the RegionalLocusCompareInput handed to the
    renderer) and where it is emitted (the Figure's text artists at savefig time).
    For the synthetic demo's 1,000,000 bp full window both say +/-500 kb."""
    from matplotlib.figure import Figure

    import locuscompare_region_render as composer

    handed: list = []
    real_render = composer.render_full_locuscompare

    def render_spy(inp, out_path):
        handed.append(inp)
        return real_render(inp, out_path)

    emitted: list[str] = []
    real_savefig = Figure.savefig

    def savefig_spy(self, *args, **kwargs):
        emitted.extend(t.get_text() for t in self.texts)
        return real_savefig(self, *args, **kwargs)

    monkeypatch.setattr(composer, "render_full_locuscompare", render_spy)
    monkeypatch.setattr(Figure, "savefig", savefig_spy)
    rc = locuscompare.main(["--demo", "01_synthetic_demo", "--output", str(tmp_path / "out")])
    assert rc == 0
    (inp,) = handed
    assert inp.window_bp == 1_000_000 and inp.title is None
    assert inp.window_label == "+/-500 kb of lead 1_500000_A_T"
    assert any("Window: +/-500 kb of lead 1_500000_A_T" in text for text in emitted), emitted
    assert not any("1000 kb" in text for text in emitted)


def test_focal_gene_highlight_styles_matching_label(tmp_path: Path):
    """render_gene_track must bold + tint the label of the gene whose symbol
    matches `focal_gene_symbol` and leave all other labels in the default
    grey-black style. The focal-gene highlight is the visual anchor for the
    QTL gene in 10-30-gene windows."""
    import matplotlib

    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    sys.path.insert(0, str(SKILL_DIR))
    from regional_plot import FOCAL_GENE_COLOR, GeneTrackEntry, render_gene_track

    genes = [
        GeneTrackEntry(gene_symbol="GENE_A", start=100, end=200, strand="+"),
        GeneTrackEntry(gene_symbol="SORT1",  start=450, end=550, strand="-"),
        GeneTrackEntry(gene_symbol="GENE_C", start=800, end=900, strand="+"),
    ]
    fig, ax = plt.subplots()
    render_gene_track(
        ax, genes, xlim_bp=(0, 1000),
        focal_gene_symbol="SORT1",
    )
    texts = {t.get_text(): t for t in ax.texts}
    sort1_label = next(t for k, t in texts.items() if "SORT1" in k)
    other_label = next(t for k, t in texts.items() if "GENE_A" in k)
    assert sort1_label.get_fontweight() == "bold"
    assert sort1_label.get_color() == FOCAL_GENE_COLOR
    assert other_label.get_fontweight() == "normal"
    assert other_label.get_color() != FOCAL_GENE_COLOR
    plt.close(fig)


def test_lead_rs_id_propagates_into_manifest_and_report(tmp_path: Path):
    """Optional `lead.rs_id` must flow end-to-end: config -> spec -> manifest +
    report. The orchestrator joins on variant_id, not rs_id; this field is
    human-readability metadata only. ai_scientist's coloc_with_mr workflow
    sets it automatically; manual configs may set it for human clarity."""
    out = tmp_path / "out"
    # Reuse the 01_synthetic_demo bundle but override the lead block with rs_id.
    base_dir = SKILL_DIR / "examples" / "01_synthetic_demo"
    base_config = json.loads((base_dir / "config.json").read_text())
    base_config["lead"]["rs_id"] = "rs99999999"
    # Make paths absolute so the override config can live in tmp_path.
    base_config["exposure"]["sumstats_path"] = str(base_dir / "exposure.tsv")
    base_config["outcome"]["sumstats_path"] = str(base_dir / "outcome.tsv")
    base_config["ld"]["ld_matrix_path"] = str(base_dir / "ld_matrix.tsv")
    base_config["gene_track"]["genes_path"] = str(base_dir / "genes.tsv")
    override_path = tmp_path / "config_with_rsid.json"
    override_path.write_text(json.dumps(base_config))
    rc = locuscompare.main(["--input", str(override_path), "--output", str(out)])
    assert rc == 0
    manifest = yaml.safe_load((out / "manifest.yaml").read_text())
    assert manifest["lead_rs_id"] == "rs99999999"
    report_text = (out / "report.md").read_text()
    assert "rs99999999" in report_text
