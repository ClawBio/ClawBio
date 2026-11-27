"""Unit tests for `cli.py` helpers.

Audit gap (PR #272): `test_locuscompare_region_render.py` covers the argparse
contract + a couple of helpers; this file fills in the remaining helpers the
demo / config / path-resolution surface depends on:

  - `_list_demos`, `_resolve_demo_path`, `_print_available_demos` (demo
    selection)
  - `_resolve_path` (relative-path resolution against the config file's dir)
  - `--list-demos` end-to-end CLI flag

The orchestrator's full `_run` path is exercised by
`test_prefetched_synthetic_demo_runs_end_to_end_offline` in the sibling file,
so we focus on the demo-resolution surface here.
"""
from __future__ import annotations

import sys
from pathlib import Path

import pytest


SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import cli  # noqa: E402


# ----- _list_demos


def test_list_demos_returns_configs_for_each_numbered_subdir():
    paths = cli._list_demos()
    # Each entry must point at a config.{json,yaml,yml} under examples/<dir>/.
    assert paths, "expected at least one bundled demo"
    names = [p.parent.name for p in paths]
    # The synthetic-demo and eqtl-cat-x-gwas-cat must always be present.
    assert "01_synthetic_demo" in names
    assert "02_eqtl_catalogue_x_gwas_catalog" in names
    # recipes/ and chains/ live alongside the numbered demos but are docs, not
    # runnable configs; _list_demos must skip them.
    assert "recipes" not in names
    assert "chains" not in names
    # Each returned path is a real config file.
    for p in paths:
        assert p.is_file()
        assert p.name.startswith("config.")


# ----- _resolve_demo_path


def test_resolve_demo_path_default_picks_eqtl_x_gwas_demo():
    p = cli._resolve_demo_path("__default__")
    assert p.parent.name == "02_eqtl_catalogue_x_gwas_catalog"
    assert p.name.startswith("config.")


def test_resolve_demo_path_named_demo():
    p = cli._resolve_demo_path("01_synthetic_demo")
    assert p.parent.name == "01_synthetic_demo"


def test_resolve_demo_path_short_prefix_match():
    """`--demo 01` should map to `01_synthetic_demo/config.*`."""
    p = cli._resolve_demo_path("01")
    assert p.parent.name == "01_synthetic_demo"


def test_resolve_demo_path_unknown_demo_lists_available():
    with pytest.raises(FileNotFoundError) as excinfo:
        cli._resolve_demo_path("99_does_not_exist")
    msg = str(excinfo.value)
    assert "no bundled demo named" in msg
    assert "01_synthetic_demo" in msg, "error should enumerate available demos"


# ----- _print_available_demos


def test_print_available_demos_lists_each_numbered_demo(capsys):
    cli._print_available_demos()
    out = capsys.readouterr().out
    assert "01_synthetic_demo" in out
    assert "02_eqtl_catalogue_x_gwas_catalog" in out
    assert "(default)" in out
    # Each bundle line includes the config filename in brackets.
    assert "[config." in out


# ----- _resolve_path


def test_resolve_path_returns_absolute_unchanged(tmp_path):
    abs_path = tmp_path / "abs.tsv"
    abs_path.write_text("")
    out = cli._resolve_path(abs_path, tmp_path / "other_config_dir")
    assert out == abs_path


def test_resolve_path_resolves_relative_against_config_dir(tmp_path):
    cfg_dir = tmp_path / "cfg"
    cfg_dir.mkdir()
    target = cfg_dir / "rel.tsv"
    target.write_text("")
    out = cli._resolve_path("rel.tsv", cfg_dir)
    assert out == target.resolve()


# ----- --list-demos end-to-end (returns 0 via main)


def test_cli_list_demos_flag_returns_zero(capsys):
    rc = cli.main(["--list-demos"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "01_synthetic_demo" in out


# ----- _format_provenance_prefix combined block (the empty + ot_release +
# gwas_lookup forms are already in test_locuscompare_region_render.py; here we
# verify the two-bit composite to lock the separator behaviour).


def test_format_provenance_prefix_combines_ot_and_gwas_lookup():
    cfg = {
        "provenance": {
            "ot_release": "26.03",
            "gwas_lookup_run_dir": "runs/sort1_vldl",
        }
    }
    prefix = cli._format_provenance_prefix(cfg)
    assert "OT release: 26.03" in prefix
    assert "gwas-lookup chain: runs/sort1_vldl" in prefix
    assert prefix.endswith(" | "), "prefix must terminate with a separator for downstream join"


# ----- Reproducibility bundle


def test_run_writes_a_verifiable_reproducibility_bundle(tmp_path):
    """Every successful run writes <output>/reproducibility/ via
    clawbio.common.reproducibility, and checksums.sha256 verifies from <output>."""
    _project_root = SKILL_DIR.parents[1]
    if str(_project_root) not in sys.path:
        sys.path.insert(0, str(_project_root))
    from clawbio.common.checksums import sha256_file

    out = tmp_path / "out"
    config = SKILL_DIR / "examples" / "01_synthetic_demo" / "config.json"
    assert cli.main(["--input", str(config), "--output", str(out)]) == 0

    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), name

    lines = (repro / "checksums.sha256").read_text().splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"1_500000_A_T_full_locuscompare.png", "manifest.yaml", "report.md"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert not Path(label).is_absolute()
        assert digest == sha256_file(out / label), label

    commands = (repro / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/skills/locuscompare-region-render/cli.py"' in commands
    assert '"$CLAWBIO_ROOT/skills/locuscompare-region-render/examples/01_synthetic_demo/config.json"' in commands
    assert '"$OUTPUT_DIR"' in commands
    assert str(out) not in commands


def test_failed_run_writes_no_bundle(tmp_path):
    out = tmp_path / "out"
    assert cli.main(["--demo", "no_such_demo", "--output", str(out)]) == 2
    assert not (out / "reproducibility").exists()
