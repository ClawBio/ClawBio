"""Scientific-validity tests for population-equity-auditor.

Red/green TDD on the bundled SYNTHETIC demo cohort (data/demo_cohort.tsv + .meta.json) — no real
patient data. The synthetic cohort is constructed to exercise the documented mechanism deterministically:
  - three LoF-tolerant-gene LoF SNVs, common in the cohort but absent/ultra-rare in gnomAD
    (DEMOLOFTOL1/2/3) -> false *actionable* under naive gnomAD-blind automation;
  - one LoF in a ClinGen haploinsufficient gene (BRCA1) common in the cohort -> the safety case where
    population frequency would MASK a real pathogenic (must be flagged, never auto-benigned);
  - benign / missense controls that must NOT be flagged.

The auditor is deterministic (ACMG/AMP + ClinGen PVS1 gate + Tavtigian points) — no model inference.
"""
from __future__ import annotations

import json
import sys
from importlib.util import module_from_spec, spec_from_file_location
from pathlib import Path

import pytest

_SKILL = Path(__file__).resolve().parents[1]


def _load(name, filename):
    spec = spec_from_file_location(name, _SKILL / filename)
    mod = module_from_spec(spec)
    sys.modules[name] = mod  # register so dataclass string-annotation fields resolve
    spec.loader.exec_module(mod)
    return mod


PEA = _load("population_equity_auditor", "population_equity_auditor.py")

ACTIONABLE = {"Pathogenic", "Likely Pathogenic"}
DEMO = _SKILL / "data" / "demo_cohort.tsv"
SIZES = json.loads((_SKILL / "data" / "demo_cohort.meta.json").read_text())["subpopulation_sizes"]


def _by_gene():
    return {r["gene"]: PEA.audit_variant(r) for r in PEA.parse_cohort_tsv(DEMO, SIZES)}


# ---------------------------------------------------------------------------
# Parsing + cohort scope (cohort-agnostic, metadata-driven)
# ---------------------------------------------------------------------------
def test_parse_cohort_tsv_computes_cohort_af_from_metadata_sizes():
    records = PEA.parse_cohort_tsv(DEMO, SIZES)
    assert len(records) == 6
    v1 = next(r for r in records if r["gene"] == "DEMOLOFTOL1")
    assert v1["cohort_af"] == pytest.approx(9 / 300, abs=1e-6)   # 9 alt alleles / (2*150)
    assert v1["strata"]["POP_G"]["af"] == pytest.approx(5 / 60, abs=1e-6)


def test_scope_is_built_from_metadata_not_hardcoded():
    data = PEA.validate_input(DEMO)
    s = data["scope"].to_dict()
    assert s["variant_classes"] == ["SNV"]
    assert s["build"] == "GRCh38"          # comes from the sidecar, not the code
    assert s["n_samples"] == 150
    assert s["assay"] == "WGS"


def test_comparability_restricts_to_shared_classes_and_flags_build_difference():
    a = PEA.CohortScope(cohort_id="a", assay="WGS", caller="c1",
                        variant_classes=frozenset({"SNV"}), build="GRCh37", n_samples=150)
    b = PEA.CohortScope(cohort_id="b", assay="WGS", caller="c2",
                        variant_classes=frozenset({"SNV", "indel"}), build="GRCh38", n_samples=5)
    shared = a.comparable_with(b)
    assert shared["variant_classes"] == frozenset({"SNV"})
    assert any("indel" in w.lower() for w in shared["warnings"])
    assert any("build" in w.lower() for w in shared["warnings"])


# ---------------------------------------------------------------------------
# The equity mechanism
# ---------------------------------------------------------------------------
def test_naive_blind_manufactures_false_actionable_calls():
    g = _by_gene()
    for gene in ("DEMOLOFTOL1", "DEMOLOFTOL2", "DEMOLOFTOL3"):
        a = g[gene]
        assert a["configs"]["naive_blind"]["acmg_class"] in ACTIONABLE, gene
        assert a["configs"]["hardened_pop"]["acmg_class"] not in ACTIONABLE, gene
        assert a["false_actionable"] is True, gene


def test_each_safeguard_independently_corrects_the_overcall():
    v = _by_gene()["DEMOLOFTOL1"]
    assert v["corrected_by_gate"] is True
    assert v["corrected_by_frequency"] is True
    assert v["corrected_by_pm2_strength"] is True


def test_controls_are_not_flagged():
    g = _by_gene()
    assert g["DEMOCOMMON1"]["false_actionable"] is False   # common in gnomAD -> BA1 both ways
    assert g["DEMOMIS1"]["false_actionable"] is False       # missense, no PVS1


# ---------------------------------------------------------------------------
# Cohort audit, PM2 epoch, and the safety invariant
# ---------------------------------------------------------------------------
def test_actionable_harm_is_specific_to_2015_pm2_moderate_strength():
    report = PEA.audit_cohort(PEA.parse_cohort_tsv(DEMO, SIZES), PEA.validate_input(DEMO)["scope"])
    assert report["false_actionable_naive_blind"] >= 3        # 2015-strength: real actionable harm
    assert report["false_actionable_pm2_supporting"] == 0     # SVI 2020 (PM2=supporting): gone
    fwr = report["referral_frequency_weighted"]
    assert fwr["hardened_pop"] < fwr["naive_blind"]           # population-aware cuts per-genome referral


def test_safety_flags_a_founder_pathogenic_masked_by_population_frequency():
    """BRCA1 (ClinGen haploinsufficient) common in the cohort: actionable from gate-surviving evidence,
    but population frequency would push it to benign -> must be FLAGGED, never auto-benigned."""
    report = PEA.audit_cohort(PEA.parse_cohort_tsv(DEMO, SIZES), PEA.validate_input(DEMO)["scope"])
    assert _by_gene()["BRCA1"]["population_masking_flag"] is True
    assert report["safety"]["unsafe_frequency_downgrades"] == 1
    assert report["safety"]["safety_invariant_holds"] is True
    assert len(report["content_sha256"]) == 64


def test_run_analysis_end_to_end(tmp_path):
    result = PEA.run_analysis(PEA.validate_input(DEMO))
    assert result["status"] == "ok"
    assert result["cohort_scope"]["variant_classes"] == ["SNV"]
    assert result["false_actionable_naive_blind"] >= 3
    PEA.write_report(result, tmp_path)
    assert (tmp_path / "report.md").exists()
    assert (tmp_path / "result.json").exists()
