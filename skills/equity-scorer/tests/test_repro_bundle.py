"""
test_repro_bundle.py — Reproducibility bundle tests for Equity Scorer.

Both pipelines (VCF and ancestry CSV) must write a reproducibility bundle
via the shared clawbio.common layer into <output_dir>/reproducibility/:
commands.sh, environment.yml, and checksums.sha256.
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


def run_vcf(tmp_path):
    out = tmp_path / "out"
    run_vcf_pipeline(DEMO_VCF, DEMO_MAP, out, DEFAULT_WEIGHTS)
    return out


def test_vcf_pipeline_writes_bundle(tmp_path):
    out = run_vcf(tmp_path)
    repro = out / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "environment.yml").exists()
    assert (repro / "checksums.sha256").exists()


def test_commands_sh_reproduces_invocation(tmp_path):
    out = run_vcf(tmp_path)
    content = (out / "reproducibility" / "commands.sh").read_text()
    assert content.startswith("#!/usr/bin/env bash")
    assert "equity_scorer.py" in content
    assert str(DEMO_VCF) in content
    assert str(DEMO_MAP) in content
    assert "--weights" in content


def test_environment_yml_names_skill_env(tmp_path):
    out = run_vcf(tmp_path)
    content = (out / "reproducibility" / "environment.yml").read_text()
    assert "name: clawbio-equity-scorer" in content
    assert "numpy" in content
    assert "scikit-learn" in content


def test_checksums_match_common_sha256(tmp_path):
    out = run_vcf(tmp_path)
    lines = (out / "reproducibility" / "checksums.sha256").read_text().strip().splitlines()
    entries = dict(reversed(line.split("  ", 1)) for line in lines)
    assert entries[DEMO_VCF.name] == sha256_file(DEMO_VCF)
    assert entries[DEMO_MAP.name] == sha256_file(DEMO_MAP)
    assert entries["report.md"] == sha256_file(out / "report.md")
    assert "result.json" in entries


def test_csv_pipeline_writes_bundle(tmp_path):
    out = tmp_path / "out"
    run_csv_pipeline(DEMO_CSV, out, DEFAULT_WEIGHTS)
    repro = out / "reproducibility"
    assert (repro / "commands.sh").exists()
    assert (repro / "environment.yml").exists()
    content = (repro / "commands.sh").read_text()
    assert str(DEMO_CSV) in content
    assert "--pop-map" not in content
    lines = (repro / "checksums.sha256").read_text().strip().splitlines()
    entries = dict(reversed(line.split("  ", 1)) for line in lines)
    assert entries[DEMO_CSV.name] == sha256_file(DEMO_CSV)
