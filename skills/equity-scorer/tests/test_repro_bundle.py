"""
test_repro_bundle.py — Reproducibility bundle tests for Equity Scorer.

Both pipelines (VCF and ancestry CSV) must write a reproducibility bundle
via the shared clawbio.common layer into <output_dir>/reproducibility/:
commands.sh, environment.yml, and checksums.sha256.

The bundle must honour the contract in docs/reproducibility.md:
`cd <output_dir> && sha256sum -c reproducibility/checksums.sha256` passes,
and a byte-identical replay produces an identical checksums.sha256.
"""

import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))
_PROJECT_ROOT = Path(__file__).resolve().parents[3]
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

from equity_scorer import DEFAULT_WEIGHTS, run_csv_pipeline, run_vcf_pipeline
from clawbio.common.checksums import sha256_file

PROJ = Path(__file__).resolve().parents[3]
DEMO_VCF = PROJ / "examples" / "demo_populations.vcf"
DEMO_MAP = PROJ / "examples" / "demo_population_map.csv"
DEMO_CSV = PROJ / "examples" / "sample_ancestry.csv"


def run_vcf(tmp_path, weights=DEFAULT_WEIGHTS):
    out = tmp_path / "out"
    run_vcf_pipeline(DEMO_VCF, DEMO_MAP, out, weights)
    return out


def read_checksums(out):
    lines = (out / "reproducibility" / "checksums.sha256").read_text().strip().splitlines()
    return dict(reversed(line.split("  ", 1)) for line in lines)


def test_vcf_pipeline_writes_bundle(tmp_path):
    out = run_vcf(tmp_path)
    repro = out / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "environment.yml").exists()
    assert (repro / "checksums.sha256").exists()


def test_commands_sh_is_portable_and_faithful(tmp_path):
    out = run_vcf(tmp_path, weights=(0.001, 0.001, 0.001, 0.997))
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert content.startswith("#!/usr/bin/env bash")
    # Script and repo-resident inputs are anchored to $CLAWBIO_ROOT, quoted.
    assert '"$CLAWBIO_ROOT/skills/equity-scorer/equity_scorer.py"' in content
    assert '"$CLAWBIO_ROOT/examples/demo_populations.vcf"' in content
    assert '"$CLAWBIO_ROOT/examples/demo_population_map.csv"' in content
    # Output is anchored to the bundle's own location, not a machine path.
    assert '"$OUTPUT_DIR"' in content
    assert str(out) not in content
    # Weights are recorded at full precision, not rounded to 2 decimals.
    assert "--weights" in content
    assert "0.001,0.001,0.001,0.997" in content
    assert "0.00," not in content


def test_environment_yml_pins_skill_versions(tmp_path):
    out = run_vcf(tmp_path)
    content = (out / "reproducibility" / "environment.yml").read_text()
    assert "name: clawbio-equity-scorer" in content
    # Version bounds from SKILL.md "Dependencies"
    assert "pandas>=2.0" in content
    assert "numpy>=1.24" in content
    assert "scikit-learn>=1.3" in content
    assert "matplotlib>=3.7" in content
    assert "python=3.11" in content


def test_checksums_verify_from_output_dir(tmp_path):
    """Every entry must resolve relative to output_dir and match the file's
    digest — the docs/reproducibility.md `sha256sum -c` contract."""
    out = run_vcf(tmp_path)
    entries = read_checksums(out)
    assert entries, "checksums.sha256 is empty"
    for label, digest in entries.items():
        target = out / label
        assert not Path(label).is_absolute(), "absolute path in manifest: %s" % label
        assert target.is_file(), "manifest entry does not resolve from output_dir: %s" % label
        assert digest == sha256_file(target), "digest mismatch for %s" % label


def test_checksums_cover_derived_artefacts(tmp_path):
    out = run_vcf(tmp_path)
    entries = read_checksums(out)
    for rel in (
        "tables/population_summary.csv",
        "tables/fst_matrix.csv",
        "tables/heterozygosity.csv",
        "tables/heim_score.json",
        "figures/heim_gauge.png",
        "figures/ancestry_bar.png",
        "figures/heterozygosity.png",
        "figures/fst_heatmap.png",
        "figures/pca_plot.png",
    ):
        assert rel in entries, "missing from manifest: %s" % rel
    # Out-of-tree inputs and wall-clock-stamped outputs must not be listed:
    # they either never resolve from output_dir or never re-verify.
    assert DEMO_VCF.name not in entries
    assert DEMO_MAP.name not in entries
    assert "report.md" not in entries
    assert "result.json" not in entries


def test_checksums_reproduce_on_replay(tmp_path):
    """A byte-identical replay must produce a byte-identical manifest."""
    first = run_vcf(tmp_path / "a")
    second = run_vcf(tmp_path / "b")
    assert (
        (first / "reproducibility" / "checksums.sha256").read_text()
        == (second / "reproducibility" / "checksums.sha256").read_text()
    )


def test_csv_pipeline_writes_bundle(tmp_path):
    out = tmp_path / "out"
    run_csv_pipeline(DEMO_CSV, out, DEFAULT_WEIGHTS)
    repro = out / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "environment.yml").exists()
    content = (repro / "commands.sh").read_text()
    assert '"$CLAWBIO_ROOT/examples/sample_ancestry.csv"' in content
    assert "--pop-map" not in content
    entries = read_checksums(out)
    assert "tables/population_summary.csv" in entries
    assert "tables/heim_score.json" in entries
    for label, digest in entries.items():
        assert digest == sha256_file(out / label)
