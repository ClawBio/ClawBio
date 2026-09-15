"""Canonical identity tests for ClawBio's bundled PRS demo panels.

Issue #356 established that the bundled files are ClawBio-curated panels,
not PGS Catalog scores. The first remediation pass added honest metadata but
left real PGS accessions as dictionary keys and filenames. These tests pin the
follow-up migration: ``CLAWBIO-*`` identifiers are canonical, PGS accessions
are provenance or a narrowly-scoped compatibility alias only.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

SKILL_DIR = Path(__file__).resolve().parents[1]
PROJECT_ROOT = SKILL_DIR.parents[1]
DATA_DIR = SKILL_DIR / "data"

EXPECTED_PANELS = {
    "CLAWBIO-T2D-8": "PGS000013",
    "CLAWBIO-AF-12": "PGS000011",
    "CLAWBIO-CAD-46": "PGS000004",
    "CLAWBIO-BC-77": "PGS000001",
    "CLAWBIO-PC-147": "PGS000057",
    "CLAWBIO-BMI-97": "PGS000039",
}


def _load_module(name: str, path: Path):
    spec = spec_from_file_location(name, path)
    module = module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


ENGINE = _load_module("gwas_prs_panel_ids", SKILL_DIR / "gwas_prs.py")
API = _load_module("gwas_prs_api_panel_ids", SKILL_DIR / "api.py")


def _run_cli(
    tmp_path: Path, *args: str, env: dict | None = None
) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [
            sys.executable,
            str(SKILL_DIR / "gwas_prs.py"),
            "--input",
            str(SKILL_DIR / "demo_patient_prs.txt"),
            "--output",
            str(tmp_path),
            *args,
        ],
        capture_output=True,
        text=True,
        timeout=120,
        env={**os.environ, **(env or {})},
    )


class TestCanonicalPanelRegistry:
    def test_registry_is_keyed_by_honest_panel_ids(self):
        assert set(ENGINE.CURATED_SCORES) == set(EXPECTED_PANELS)

    @pytest.mark.parametrize("panel_id,legacy_pgs_id", EXPECTED_PANELS.items())
    def test_legacy_accession_is_metadata_not_identity(self, panel_id, legacy_pgs_id):
        meta = ENGINE.CURATED_SCORES[panel_id]
        assert meta["curated_panel_id"] == panel_id
        assert meta["legacy_pgs_id"] == legacy_pgs_id

    @pytest.mark.parametrize("panel_id", EXPECTED_PANELS)
    def test_panel_files_use_clawbio_names(self, panel_id):
        path = ENGINE.curated_panel_path(panel_id, "GRCh37")
        assert path == DATA_DIR / f"{panel_id}_GRCh37.txt"
        assert path.exists()
        assert ENGINE.is_curated_demo_panel(path)

    @pytest.mark.parametrize("legacy_pgs_id", EXPECTED_PANELS.values())
    def test_no_panel_occupies_a_catalog_cache_path(self, legacy_pgs_id):
        assert not (DATA_DIR / f"{legacy_pgs_id}_hmPOS_GRCh37.txt").exists()
        assert not (DATA_DIR / f"{legacy_pgs_id}_hmPOS_GRCh37.txt.gz").exists()

    def test_only_the_pinned_benchmark_alias_remains(self):
        assert ENGINE.LEGACY_PGS_PANEL_COMPAT == {
            "PGS000013": "CLAWBIO-T2D-8",
        }


class TestPanelCli:
    def test_panel_id_scores_the_honest_panel(self, tmp_path):
        proc = _run_cli(tmp_path, "--panel-id", "CLAWBIO-T2D-8")
        assert proc.returncode == 0, proc.stdout + proc.stderr
        results = json.loads((tmp_path / "prs_results.json").read_text())
        assert len(results) == 1
        assert results[0]["score_id"] == "CLAWBIO-T2D-8"
        assert results[0]["curated_panel_id"] == "CLAWBIO-T2D-8"
        assert results[0]["pgs_id"] is None
        assert results[0]["legacy_pgs_id"] == "PGS000013"
        assert results[0]["legacy_pgs_compatibility"] is False

    def test_result_envelope_uses_the_panel_id(self, tmp_path):
        proc = _run_cli(tmp_path, "--panel-id", "CLAWBIO-T2D-8")
        assert proc.returncode == 0, proc.stdout + proc.stderr
        summary = json.loads((tmp_path / "result.json").read_text())["summary"]
        assert summary["score_id"] == "CLAWBIO-T2D-8"
        assert summary["curated_panel_id"] == "CLAWBIO-T2D-8"
        assert summary["pgs_id"] is None

    def test_unknown_panel_id_is_rejected_without_network(self, tmp_path):
        proc = _run_cli(tmp_path, "--panel-id", "CLAWBIO-NOT-REAL")
        assert proc.returncode != 0
        assert "unknown curated panel id" in (proc.stdout + proc.stderr).lower()

    def test_conflicting_score_selectors_are_rejected(self, tmp_path):
        proc = _run_cli(
            tmp_path,
            "--panel-id",
            "CLAWBIO-T2D-8",
            "--pgs-id",
            "PGS000031",
        )
        assert proc.returncode != 0
        assert "choose exactly one" in (proc.stdout + proc.stderr).lower()

    def test_pinned_benchmark_alias_stays_machine_compatible(self, tmp_path):
        # The alias is opt-in since #356's refusal landed; the benchmark
        # workflow sets this variable. See test_legacy_alias_opt_in.py.
        proc = _run_cli(
            tmp_path, "--pgs-id", "PGS000013",
            env={"CLAWBIO_ALLOW_LEGACY_PGS_ALIAS": "1"},
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "benchmark compatibility" in (proc.stdout + proc.stderr).lower()
        result = json.loads((tmp_path / "prs_results.json").read_text())[0]
        assert result["score_id"] == "CLAWBIO-T2D-8"
        assert result["pgs_id"] == "PGS000013"
        assert result["curated_panel_id"] == "CLAWBIO-T2D-8"
        assert result["legacy_pgs_compatibility"] is True

    def test_top_level_runner_forwards_panel_id(self, tmp_path):
        proc = subprocess.run(
            [
                sys.executable,
                str(PROJECT_ROOT / "clawbio.py"),
                "run",
                "prs",
                "--input",
                str(SKILL_DIR / "demo_patient_prs.txt"),
                "--panel-id",
                "CLAWBIO-T2D-8",
                "--output",
                str(tmp_path),
            ],
            capture_output=True,
            text=True,
            timeout=120,
        )
        assert proc.returncode == 0, proc.stdout + proc.stderr
        result = json.loads((tmp_path / "prs_results.json").read_text())[0]
        assert result["score_id"] == "CLAWBIO-T2D-8"


class TestImportableApi:
    def test_api_accepts_curated_panel_id(self):
        _fmt, _count, genotypes = ENGINE.load_genotypes(
            SKILL_DIR / "demo_patient_prs.txt"
        )
        result = API.run(
            genotypes,
            {"curated_panel_id": "CLAWBIO-T2D-8", "min_overlap": 0.0},
        )
        assert result["scores_calculated"] == 1
        assert result["results"][0]["score_id"] == "CLAWBIO-T2D-8"
        assert result["results"][0]["pgs_id"] is None


class TestGeneratedMetadata:
    def test_json_uses_the_same_canonical_keys(self):
        data = json.loads((SKILL_DIR / "curated_scores.json").read_text())
        assert set(data) == set(EXPECTED_PANELS)

    def test_consumer_examples_do_not_label_panels_as_pgs_scores(self):
        ancestry = (
            PROJECT_ROOT / "skills" / "ancestry-risk-profiler"
            / "ancestry_risk_profiler.py"
        ).read_text()
        workshop = (PROJECT_ROOT / "workshop-gwas-slides.html").read_text()
        profile = (
            PROJECT_ROOT / "skills" / "profile-report" / "demo_full_profile.json"
        ).read_text()
        assert "--pgs-id PGS000013` for Type 2 Diabetes" not in ancestry
        for legacy_pgs_id in EXPECTED_PANELS.values():
            assert f"<td>{legacy_pgs_id}</td>" not in workshop
            assert f'"pgs_id": "{legacy_pgs_id}"' not in profile
