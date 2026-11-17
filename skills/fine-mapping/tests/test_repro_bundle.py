"""Reproducibility bundle tests for fine-mapping.

Every run (SuSiE demo or ABF on user sumstats) must write
<output_dir>/reproducibility/{commands.sh, environment.yml, checksums.sha256}
via clawbio.common.reproducibility, and
`cd <output_dir> && sha256sum -c reproducibility/checksums.sha256` must pass.
"""
import importlib.util
import shlex
import sys
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

import fine_mapping

from clawbio.common.checksums import sha256_file

needs_jax = pytest.mark.skipif(
    importlib.util.find_spec("jax") is None,
    reason="needs jax (uv sync --extra fine-mapping)",
)


def read_checksums(out: Path) -> dict[str, str]:
    lines = (out / "reproducibility" / "checksums.sha256").read_text().strip().splitlines()
    return {label: digest for digest, label in (line.split("  ", 1) for line in lines)}


def assert_bundle_verifies(out: Path) -> dict[str, str]:
    repro = out / "reproducibility"
    for name in ("commands.sh", "environment.yml", "checksums.sha256"):
        assert (repro / name).is_file(), f"missing reproducibility/{name}"
    entries = read_checksums(out)
    assert entries, "checksums.sha256 is empty"
    for label, digest in entries.items():
        assert not Path(label).is_absolute(), label
        assert (out / label).is_file(), f"entry does not resolve from output_dir: {label}"
        assert digest == sha256_file(out / label), f"digest mismatch for {label}"
    return entries


@needs_jax
def test_demo_writes_verifiable_bundle(tmp_path):
    out = tmp_path / "out"
    fine_mapping.run_finemapping(None, None, out, demo=True)
    entries = assert_bundle_verifies(out)
    assert {"report.md", "fine_mapping.json", "tables/pips.tsv",
            "figures/pip_locus_plot.png", "figures/ld_heatmap.png"} <= set(entries)
    commands = (out / "reproducibility" / "commands.sh").read_text()
    assert commands.startswith("#!/usr/bin/env bash")
    assert "--demo" in commands
    assert str(out) not in commands  # output anchored to $OUTPUT_DIR, not a machine path


def test_abf_run_writes_verifiable_bundle(tmp_path):
    df, _ = fine_mapping.make_demo_data()
    sumstats = tmp_path / "locus.tsv"
    df.to_csv(sumstats, sep="\t", index=False)
    out = tmp_path / "out"
    # Stale figure from an earlier LD run must not be certified by this ABF run.
    (out / "figures").mkdir(parents=True)
    (out / "figures" / "ld_heatmap.png").write_bytes(b"stale")

    fine_mapping.run_finemapping(sumstats, None, out, coverage=0.9)
    entries = assert_bundle_verifies(out)
    assert "figures/ld_heatmap.png" not in entries
    assert "figures/pip_locus_plot.png" in entries
    commands = (out / "reproducibility" / "commands.sh").read_text()
    tokens = shlex.split(commands.replace("\\\n", ""), comments=True)
    assert tokens[tokens.index("--coverage") + 1] == "0.9"
    assert tokens[tokens.index("--sumstats") + 1] == str(sumstats.resolve())
