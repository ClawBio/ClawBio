"""Issue #356, remaining half: the ``--pgs-id PGS000013`` alias must be opt-in.

PGS000013 is the PGS Catalog accession for Khera 2018 coronary artery disease
(6,630,150 variants). ClawBio ships an 8-variant curated type 2 diabetes panel
that historically answered to that accession. After #357 and #380 the alias
still fired silently on every ``--pgs-id PGS000013`` run, because the pinned
clawbio_bench revision drives eight scoring cases through exactly that path.

These tests pin the resolution: the alias fires only when the caller sets
``CLAWBIO_ALLOW_LEGACY_PGS_ALIAS=1``. The benchmark workflow sets it. Nothing
else does, so an ordinary user asking for PGS000013 is told what the accession
really is and pointed at ``--panel-id CLAWBIO-T2D-8`` for the panel.
"""

from __future__ import annotations

import json
import os
import subprocess
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import yaml

SKILL_DIR = Path(__file__).resolve().parent.parent
PROJECT_ROOT = SKILL_DIR.parent.parent
OPT_IN = "CLAWBIO_ALLOW_LEGACY_PGS_ALIAS"

# A proxy nothing listens on. requests honours HTTPS_PROXY, so any attempt to
# reach the PGS Catalog fails immediately and deterministically instead of
# fetching a 6.6M-variant score in the middle of a unit test.
OFFLINE = {"HTTPS_PROXY": "http://127.0.0.1:9", "HTTP_PROXY": "http://127.0.0.1:9",
           "NO_PROXY": ""}


def _load_engine():
    spec = spec_from_file_location("gwas_prs_alias_opt_in", SKILL_DIR / "gwas_prs.py")
    module = module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    return module


ENGINE = _load_engine()


def _run(tmp_path: Path, *args: str, env: dict | None = None) -> subprocess.CompletedProcess[str]:
    merged = {k: v for k, v in os.environ.items() if k != OPT_IN}
    merged.update(OFFLINE)
    merged.update(env or {})
    return subprocess.run(
        [sys.executable, str(SKILL_DIR / "gwas_prs.py"),
         "--input", str(SKILL_DIR / "demo_patient_prs.txt"),
         "--output", str(tmp_path), *args],
        capture_output=True, text=True, timeout=120, env=merged,
    )


class TestOptInSwitch:
    def test_disabled_when_unset(self):
        assert ENGINE.legacy_alias_enabled({}) is False

    def test_enabled_only_by_an_explicit_true_value(self):
        for value in ("1", "true", "TRUE", "yes", "on"):
            assert ENGINE.legacy_alias_enabled({OPT_IN: value}) is True, value
        for value in ("", "0", "false", "no", "off", "PGS000013"):
            assert ENGINE.legacy_alias_enabled({OPT_IN: value}) is False, value

    def test_env_var_name_is_the_documented_one(self):
        assert ENGINE.LEGACY_ALIAS_ENV == OPT_IN


class TestRefusalWithoutOptIn:
    def test_pgs000013_is_not_substituted(self, tmp_path):
        proc = _run(tmp_path, "--pgs-id", "PGS000013")
        combined = proc.stdout + proc.stderr
        assert "benchmark compatibility" not in combined.lower()
        assert "Traceback" not in proc.stderr
        # Offline, the genuine Catalog fetch cannot succeed, so the run fails
        # rather than quietly scoring the wrong panel.
        assert proc.returncode != 0
        assert "--panel-id CLAWBIO-T2D-8" in combined
        assert "PGS Catalog accession" in combined
        results = tmp_path / "prs_results.json"
        if results.exists():
            for row in json.loads(results.read_text()):
                assert row.get("legacy_pgs_compatibility") is not True

    def test_refusal_names_what_the_accession_really_is(self, tmp_path):
        proc = _run(tmp_path, "--pgs-id", "PGS000013")
        combined = (proc.stdout + proc.stderr).lower()
        assert "coronary artery disease" in combined
        assert "khera" in combined

    def test_a_false_value_does_not_enable_the_alias(self, tmp_path):
        proc = _run(tmp_path, "--pgs-id", "PGS000013", env={OPT_IN: "0"})
        assert "benchmark compatibility" not in (proc.stdout + proc.stderr).lower()
        assert proc.returncode != 0


class TestAliasWithOptIn:
    def test_benchmark_path_still_works_when_opted_in(self, tmp_path):
        proc = _run(tmp_path, "--pgs-id", "PGS000013", env={OPT_IN: "1"})
        assert proc.returncode == 0, proc.stdout + proc.stderr
        assert "benchmark compatibility" in (proc.stdout + proc.stderr).lower()
        result = json.loads((tmp_path / "prs_results.json").read_text())[0]
        assert result["score_id"] == "CLAWBIO-T2D-8"
        assert result["pgs_id"] == "PGS000013"
        assert result["legacy_pgs_compatibility"] is True

    def test_opt_in_does_not_widen_the_alias_table(self):
        assert ENGINE.LEGACY_PGS_PANEL_COMPAT == {"PGS000013": "CLAWBIO-T2D-8"}


class TestTheBenchmarkWorkflowOptsIn:
    """The refusal is only safe to ship if the pinned benchmark keeps working.
    That depends on one line of workflow YAML, so pin it."""

    def test_smoke_run_step_sets_the_opt_in(self):
        wf = yaml.safe_load(
            (PROJECT_ROOT / ".github" / "workflows" / "bench-leaderboard.yml").read_text()
        )
        steps = wf["jobs"]["bench"]["steps"]
        smoke = [s for s in steps if "uv run clawbio-bench" in str(s.get("run", ""))]
        assert len(smoke) == 1, "expected exactly one step that runs clawbio-bench"
        env = smoke[0].get("env") or {}
        assert str(env.get(OPT_IN)) == "1"


class TestDocsAgreeWithCode:
    def test_skill_md_documents_the_opt_in(self):
        text = (SKILL_DIR / "SKILL.md").read_text()
        assert OPT_IN in text
        assert "--panel-id CLAWBIO-T2D-8" in text
