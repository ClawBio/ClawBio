"""Tests for prs-applicability-gate.

The gate is a pure, deterministic function of (gate input, calibration). These tests cover every rule outcome and
reason code, fail-closed handling of malformed input, determinism, tamper detection, the CLI and output contract,
the synthetic demo, and the ClawBio runner boundary: the shipped calibration cannot be replaced through
`clawbio.py run prs-gate`.
"""

from __future__ import annotations

import copy
import hashlib
import importlib.util
import json
import os
import random
import re
import subprocess
import sys
from pathlib import Path

import pytest
import yaml

SKILL = Path(__file__).resolve().parents[1]
REPO_ROOT = SKILL.parents[1]
SCRIPT = SKILL / "prs_applicability_gate.py"
EXAMPLES = SKILL / "examples"
SHIPPED_CONFIG = SKILL / "config" / "calibration.yaml"
CLAWBIO_DISCLAIMER = (
    "ClawBio is a research and educational tool. It is not a medical device and does not provide clinical "
    "diagnoses. Consult a healthcare professional before making any medical decisions."
)

spec = importlib.util.spec_from_file_location("_prs_gate_under_test", SCRIPT)
gate = importlib.util.module_from_spec(spec)
sys.modules["_prs_gate_under_test"] = gate
spec.loader.exec_module(gate)
CFG = gate.load_config()

DEMO_EXPECTED = {
    "synthetic_A_supported.gate_input.json": ("SUPPORTED", []),
    "synthetic_B_raw_only.gate_input.json": ("RAW_ONLY", ["TARGET_REFERENCE_UNRESOLVED"]),
    "synthetic_C_abstain.gate_input.json": (
        "ABSTAIN", ["LOW_SCOREABILITY", "VARIANTS_MISSING", "TARGET_REFERENCE_UNRESOLVED"]),
}
SUPPORTED_EXAMPLE = EXAMPLES / "synthetic_A_supported.gate_input.json"


def base_input() -> dict:
    """A fully supported case: every rule passes."""
    return {
        "schema": gate.INPUT_SCHEMA,
        "candidate": {"pgs_id": "PGS000000", "pre_rank": 1, "trait_reported": "Breast cancer", "sex_specific": None},
        "score_file": {"name": "x.txt.gz", "sha256": "0" * 64, "build": "GRCh37", "weight_type": "beta",
                       "ratio_weight_type": False, "unsupported_features": [], "n_parse_problems": 0,
                       "variants_interactions": 0},
        "genotype": {"file_sha256": "1" * 64, "format": "vcf", "n_calls": 10000, "n_called": 10000,
                     "build": "GRCh37", "build_evidence": {"method": "empirical"}},
        "person": {"sex": "female"},
        "harmonisation": {"n_variants": 300, "n_matched": 290, "n_located": 295, "status_counts": {},
                          "allele_mismatch_fraction": 0.0, "weight_loss_by_status": {"palindromic_excluded": 0.02},
                          "fraction_matched": 0.97, "fraction_abs_weight_matched": 0.98},
        "scoreability": {"r": 0.98, "method": "reference_panel_correlation"},
        "placement": {"status": "RESOLVED", "placement": "EUR", "nearest_reference": "EUR",
                      "placement_stability": 1.0, "n_sites_used": 6000, "detail": "inside EUR"},
        "catalog_metadata": {"status": "resolved", "failed_checks": [], "detail": "ok", "warnings": []},
        "evaluation": {"reported": True, "unit": "(publication, sample set) pairs", "units": [
            {"pgp_id": "PGP1", "pss_id": "PSS1", "code": "EUR", "pooled": False, "n": 5000, "cases": 900,
             "percent_male": 0.0, "metrics": [{"name": "OR", "estimate": 1.6, "ci_lower": 1.5, "ci_upper": 1.7,
                                               "null": 1.0}]},
            {"pgp_id": "PGP2", "pss_id": "PSS2", "code": "AFR", "pooled": False, "n": 800, "cases": 100,
             "percent_male": 0.0, "metrics": [{"name": "AUROC", "estimate": 0.55, "ci_lower": 0.49,
                                               "ci_upper": 0.61, "null": 0.5}]}]},
        "reference_distribution": {"available": True, "reference_group": "EUR", "reference_n": 502,
                                   "n_intersection": 290, "reference_sensitive": False,
                                   "reference_sensitive_pairs": [], "detail": None},
    }


def run(gi: dict) -> dict:
    return gate.evaluate(gi, CFG)


def rule(res: dict, rid: str) -> dict:
    return next(r for r in res["rule_trace"] if r["rule"] == rid)


def mutate(**changes) -> dict:
    gi = base_input()
    for path, value in changes.items():
        obj = gi
        keys = path.split("__")
        for k in keys[:-1]:
            obj = obj[k]
        obj[keys[-1]] = value
    return gi


def run_cli(*args, cwd: Path | None = None) -> subprocess.CompletedProcess:
    return subprocess.run([sys.executable, str(SCRIPT), *map(str, args)], capture_output=True, text=True,
                          cwd=cwd, timeout=120)


def load_result(out: Path) -> dict:
    return json.loads((out / "result.json").read_text(encoding="utf-8"))


def alt_config(tmp_path: Path, r_min: str = "0.99") -> Path:
    alt = tmp_path / "alt_calibration.yaml"
    alt.write_text(SHIPPED_CONFIG.read_text(encoding="utf-8").replace("r_min: 0.90", f"r_min: {r_min}"),
                   encoding="utf-8")
    return alt


def sha256_of(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


# ---- rules and reason codes ---------------------------------------------------------------------------------------

def test_supported_base_case():
    res = run(base_input())
    assert res["status"] == "SUPPORTED"
    assert res["reason_codes"] == [] and res["primary_reason"] is None
    assert res["allowed_claims"] == {"raw_score": True, "standardized_score": True, "percentile": True,
                                     "absolute_risk": False}
    assert all(r["outcome"] in ("pass", "not_applicable") for r in res["rule_trace"])
    assert [r["rule"] for r in res["rule_trace"]] == [f"G{i}" for i in range(1, 13)]
    assert res["calibration_version"] == CFG["calibration_version"]
    for key in ("status", "primary_reason", "reason_codes", "evidence_used", "evidence_missing",
                "calibration_version", "rule_trace", "what_would_change_result", "provenance"):
        assert key in res


@pytest.mark.parametrize("changes, status, code", [
    ({"score_file__unsupported_features": ["dosage_0_weight"]}, "ABSTAIN", "UNSUPPORTED_SCORE_FORMAT"),
    ({"score_file__ratio_weight_type": True, "score_file__weight_type": "OR"}, "ABSTAIN", "UNSUPPORTED_SCORE_FORMAT"),
    ({"score_file__variants_interactions": 3}, "ABSTAIN", "UNSUPPORTED_SCORE_FORMAT"),
    ({"score_file__n_parse_problems": 2}, "ABSTAIN", "UNSUPPORTED_SCORE_FORMAT"),
    ({"genotype__build": "UNRESOLVED"}, "ABSTAIN", "BUILD_UNRESOLVED"),
    ({"harmonisation__allele_mismatch_fraction": 0.4}, "ABSTAIN", "ALLELE_HARMONIZATION_FAILED"),
    ({"scoreability__r": 0.7}, "ABSTAIN", "LOW_SCOREABILITY"),
    ({"scoreability__r": None, "scoreability__method": None}, "ABSTAIN", "SCOREABILITY_UNVERIFIED"),
    ({"candidate__sex_specific": "female", "person__sex": "male"}, "ABSTAIN", "SEX_POPULATION_MISMATCH"),
    ({"candidate__sex_specific": "female", "person__sex": None}, "RAW_ONLY", "SEX_NOT_PROVIDED"),
    ({"catalog_metadata__status": "contradictory"}, "RAW_ONLY", "METADATA_CONTRADICTION"),
    ({"catalog_metadata__status": "unresolved"}, "RAW_ONLY", "EVALUATION_METADATA_UNAVAILABLE"),
    ({"placement__status": "INTERMEDIATE", "placement__placement": None}, "RAW_ONLY", "TARGET_REFERENCE_UNRESOLVED"),
    ({"placement__status": "UNSTABLE", "placement__placement": None}, "RAW_ONLY", "TARGET_REFERENCE_UNRESOLVED"),
    ({"placement__status": "UNRESOLVED", "placement__placement": None}, "RAW_ONLY", "TARGET_REFERENCE_UNRESOLVED"),
    ({"placement__placement": "SAS"}, "RAW_ONLY", "NO_RELEVANT_EVALUATION"),
    ({"placement__placement": "AFR"}, "RAW_ONLY", "EVALUATION_NOT_INFORMATIVE"),
    ({"person__sex": "male"}, "RAW_ONLY", "SEX_POPULATION_MISMATCH"),
    ({"reference_distribution__available": False}, "RAW_ONLY", "REFERENCE_DISTRIBUTION_UNAVAILABLE"),
    ({"reference_distribution__reference_sensitive": True,
      "reference_distribution__reference_sensitive_pairs": [["FIN", "TSI"]]}, "RAW_ONLY", "REFERENCE_SENSITIVE"),
])
def test_each_reason_code(changes, status, code):
    gi = mutate(**changes)
    if gi["placement"]["placement"] in ("SAS", "AFR"):
        gi["reference_distribution"]["reference_group"] = gi["placement"]["placement"]
    res = run(gi)
    assert res["status"] == status, res["rule_trace"]
    assert code in res["reason_codes"]
    assert res["primary_reason"]["code"] in res["reason_codes"]
    assert res["what_would_change_result"], "every failure must say what would change it"
    assert res["allowed_claims"]["percentile"] is False
    assert res["allowed_claims"]["absolute_risk"] is False
    assert res["allowed_claims"]["raw_score"] is (status == "RAW_ONLY")


# ---- G8: a RESOLVED claim must be backed by enough sites and a stable placement ---------------------------------

def _assert_g8_refuses(res: dict) -> None:
    assert rule(res, "G8")["outcome"] == "fail", rule(res, "G8")
    assert res["status"] == "RAW_ONLY" and res["reason_codes"] == ["TARGET_REFERENCE_UNRESOLVED"]
    assert res["allowed_claims"]["percentile"] is False and res["allowed_claims"]["standardized_score"] is False
    for rid in ("G9", "G10", "G11", "G12"):
        assert rule(res, rid)["outcome"] == "not_applicable", f"{rid} must not use an unverified group"


@pytest.mark.parametrize("changes", [
    {"placement__n_sites_used": 47, "placement__placement_stability": 0.5},  # the over-call class from review
    {"placement__n_sites_used": 47},
    {"placement__n_sites_used": 199},
    {"placement__n_sites_used": 0},
    {"placement__n_sites_used": -5},
    {"placement__n_sites_used": None},
    {"placement__n_sites_used": "6000"},
    {"placement__n_sites_used": True},
    {"placement__n_sites_used": 250.0},
    {"placement__placement_stability": 0.5},
    {"placement__placement_stability": 0.949},
    {"placement__placement_stability": None},
])
def test_g8_resolved_needs_min_sites_and_min_stability(changes):
    _assert_g8_refuses(run(mutate(**changes)))


@pytest.mark.parametrize("field", ["n_sites_used", "placement_stability"])
def test_g8_missing_field_fails_and_is_reported_missing(field):
    gi = base_input()
    del gi["placement"][field]
    res = run(gi)
    _assert_g8_refuses(res)
    assert f"placement.{field}" in res["evidence_missing"]


def test_g8_thresholds_are_inclusive():
    p = CFG["placement"]
    res = run(mutate(placement__n_sites_used=p["min_sites"], placement__placement_stability=p["min_stability"]))
    assert rule(res, "G8")["outcome"] == "pass" and res["status"] == "SUPPORTED"


def test_g8_thresholds_come_from_the_calibration(tmp_path):
    alt = tmp_path / "alt_placement.yaml"
    alt.write_text(SHIPPED_CONFIG.read_text(encoding="utf-8").replace("min_sites: 200", "min_sites: 7000")
                   .replace("min_stability: 0.95", "min_stability: 0.999"), encoding="utf-8")
    cfg = gate.load_config(alt)
    gi = mutate(placement__n_sites_used=6000)
    assert rule(gate.evaluate(gi, cfg), "G8")["outcome"] == "fail"
    assert rule(gate.evaluate(mutate(placement__placement_stability=0.998), cfg), "G8")["outcome"] == "fail"
    assert rule(gate.evaluate(base_input(), CFG), "G8")["outcome"] == "pass"


@pytest.mark.parametrize("bad", ["min_sites: 0", "min_sites: 2.5", "min_sites: true", "min_sites: null"])
def test_config_requires_a_positive_integer_min_sites(tmp_path, bad):
    alt = tmp_path / "bad_sites.yaml"
    alt.write_text(SHIPPED_CONFIG.read_text(encoding="utf-8").replace("min_sites: 200", bad), encoding="utf-8")
    with pytest.raises(ValueError):
        gate.load_config(alt)


@pytest.mark.parametrize("changes", [
    {},
    {"placement__status": "INTERMEDIATE", "placement__placement": None},
    {"placement__status": "UNRESOLVED", "placement__placement": None, "person__sex": None},
])
def test_valid_input_traces_every_rule(changes):
    """SKILL.md promises a complete trace for valid input: G1-G12, each exactly once, in order."""
    assert [r["rule"] for r in run(mutate(**changes))["rule_trace"]] == [f"G{i}" for i in range(1, 13)]


def test_every_reason_code_is_documented_in_skill_md():
    skill_md = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    for code in gate.REASON_CODES:
        assert code in skill_md, f"{code} missing from SKILL.md"


def test_abstain_outranks_raw_only_and_primary_is_most_severe():
    res = run(mutate(scoreability__r=0.5, placement__status="INTERMEDIATE", placement__placement=None))
    assert res["status"] == "ABSTAIN"
    assert res["primary_reason"]["code"] == "LOW_SCOREABILITY"
    assert "TARGET_REFERENCE_UNRESOLVED" in res["reason_codes"]


def test_low_scoreability_names_the_largest_loss():
    res = run(mutate(scoreability__r=0.8, harmonisation__weight_loss_by_status={"palindromic_excluded": 0.3,
                                                                                "missing": 0.1}))
    assert res["reason_codes"][:2] == ["LOW_SCOREABILITY", "PALINDROMIC_VARIANT_UNRESOLVED"]


def test_threshold_boundaries_are_inclusive_as_documented():
    assert run(mutate(scoreability__r=CFG["scoreability"]["r_min"]))["status"] == "SUPPORTED"
    cut = CFG["allele_harmonisation"]["max_mismatch_fraction"]
    assert run(mutate(harmonisation__allele_mismatch_fraction=cut))["status"] == "SUPPORTED"


def test_pooled_units_do_not_count_as_group_evidence():
    gi = base_input()
    for u in gi["evaluation"]["units"]:
        u["pooled"] = True
    assert "NO_RELEVANT_EVALUATION" in run(gi)["reason_codes"]


def test_metric_without_ci_is_not_informative():
    gi = base_input()
    gi["evaluation"]["units"][0]["metrics"] = [{"name": "AUROC", "estimate": 0.65, "ci_lower": None,
                                                "ci_upper": None, "null": 0.5}]
    assert "EVALUATION_NOT_INFORMATIVE" in run(gi)["reason_codes"]


def test_inverse_association_is_not_supporting_evidence():
    """A CI entirely below the null (e.g. a case-only subtype comparison, OR 0.86 [0.82, 0.89]) is not evidence
    for interpreting the score, even though it excludes the null."""
    gi = base_input()
    gi["evaluation"]["units"][0]["metrics"] = [{"name": "OR", "estimate": 0.86, "ci_lower": 0.82, "ci_upper": 0.89,
                                                "null": 1.0}]
    res = run(gi)
    assert "EVALUATION_NOT_INFORMATIVE" in res["reason_codes"] and res["status"] == "RAW_ONLY"


def test_missing_evidence_is_listed_not_inferred():
    res = run(mutate(scoreability__r=None))
    assert "scoreability.r" in res["evidence_missing"]


def test_negative_correlation_is_valid_input_but_low_scoreability():
    res = run(mutate(scoreability__r=-0.5))
    assert res["status"] == "ABSTAIN" and "LOW_SCOREABILITY" in res["reason_codes"]


# ---- unknown evidence never passes a rule ---------------------------------------------------------------------------
# SKILL.md: "Missing stays in evidence_missing and the affected rule fails." Each test below fails if the guard it
# names is deleted or relaxed.

MISSING = object()  # sentinel: delete the key instead of setting a value


def with_field(block: str, key: str, value, gi: dict | None = None) -> dict:
    gi = base_input() if gi is None else gi
    if value is MISSING:
        del gi[block][key]
    else:
        gi[block][key] = value
    return gi


def set_units(gi: dict, *percent_male) -> dict:
    """Make the informative EUR evaluation units carry exactly these percent_male values (MISSING deletes it)."""
    template = gi["evaluation"]["units"][0]
    units = []
    for pm in percent_male:
        u = copy.deepcopy(template)
        if pm is MISSING:
            del u["percent_male"]
        else:
            u["percent_male"] = pm
        units.append(u)
    gi["evaluation"]["units"] = units
    return gi


# G2: ratio_weight_type and unsupported_features must be reported explicitly.

@pytest.mark.parametrize("value, outcome, code", [
    (False, "pass", None),
    (True, "fail", "UNSUPPORTED_SCORE_FORMAT"),
    (MISSING, "fail", "SCORE_FORMAT_UNVERIFIED"),
    (None, "fail", "SCORE_FORMAT_UNVERIFIED"),
    ("false", "fail", "SCORE_FORMAT_UNVERIFIED"),
    (0, "fail", "SCORE_FORMAT_UNVERIFIED"),
    (1, "fail", "SCORE_FORMAT_UNVERIFIED"),
])
def test_g2_ratio_weight_type_must_be_an_explicit_boolean(value, outcome, code):
    res = run(with_field("score_file", "ratio_weight_type", value, mutate(score_file__weight_type="OR")))
    assert rule(res, "G2")["outcome"] == outcome, rule(res, "G2")
    if code is None:
        assert res["status"] == "SUPPORTED"
    else:
        assert res["status"] == "ABSTAIN" and res["primary_reason"]["code"] == code
        assert res["allowed_claims"]["raw_score"] is False
    if code == "SCORE_FORMAT_UNVERIFIED":
        assert "score_file.ratio_weight_type" in res["evidence_missing"]


@pytest.mark.parametrize("value", [MISSING, None, "", 0, {}])
def test_g2_unsupported_features_must_be_a_reported_list(value):
    res = run(with_field("score_file", "unsupported_features", value))
    assert res["status"] == "ABSTAIN" and "SCORE_FORMAT_UNVERIFIED" in res["reason_codes"]


def test_g2_known_unsupported_feature_outranks_unverified_flag():
    gi = with_field("score_file", "ratio_weight_type", MISSING, mutate(score_file__variants_interactions=3))
    res = run(gi)
    assert res["reason_codes"][:2] == ["UNSUPPORTED_SCORE_FORMAT", "SCORE_FORMAT_UNVERIFIED"]


# G3: the build must be one of the supported assemblies, exactly.

def test_supported_builds_are_the_pgs_catalog_harmonised_assemblies():
    assert gate.SUPPORTED_BUILDS == ("GRCh37", "GRCh38")


@pytest.mark.parametrize("build", ["GRCh37", "GRCh38"])
def test_g3_passes_every_supported_build(build):
    res = run(mutate(genotype__build=build))
    assert rule(res, "G3")["outcome"] == "pass" and res["status"] == "SUPPORTED"


@pytest.mark.parametrize("build", [MISSING, None, "", "UNRESOLVED", "GRCh99", "NCBI36", "hg19", "grch38",
                                   " GRCh38", 38, ["GRCh38"], {"build": "GRCh38"}])
def test_g3_fails_anything_but_a_supported_build(build):
    res = run(with_field("genotype", "build", build))
    assert rule(res, "G3")["outcome"] == "fail", rule(res, "G3")
    assert res["status"] == "ABSTAIN" and res["primary_reason"]["code"] == "BUILD_UNRESOLVED"
    assert "genotype.build" in res["evidence_missing"]


# G5: the largest loss is the largest value, not the first key.

@pytest.mark.parametrize("loss", [{"missing": 0.1, "palindromic_excluded": 0.3},
                                  {"palindromic_excluded": 0.3, "missing": 0.1}])
def test_g5_names_the_largest_loss_whatever_the_key_order(loss):
    res = run(mutate(scoreability__r=0.8, harmonisation__weight_loss_by_status=loss))
    assert res["reason_codes"] == ["LOW_SCOREABILITY", "PALINDROMIC_VARIANT_UNRESOLVED"]
    assert rule(res, "G5")["detail"].endswith("largest loss: palindromic_excluded")


def test_g5_tied_largest_losses_are_all_named_independently_of_key_order():
    a = run(mutate(scoreability__r=0.8,
                   harmonisation__weight_loss_by_status={"palindromic_excluded": 0.2, "missing": 0.2}))
    b = run(mutate(scoreability__r=0.8,
                   harmonisation__weight_loss_by_status={"missing": 0.2, "palindromic_excluded": 0.2}))
    assert a["reason_codes"] == b["reason_codes"] == ["LOW_SCOREABILITY", "VARIANTS_MISSING",
                                                      "PALINDROMIC_VARIANT_UNRESOLVED"]
    assert rule(a, "G5")["detail"] == rule(b, "G5")["detail"]


def test_g5_zero_loss_names_no_loss():
    res = run(mutate(scoreability__r=0.8, harmonisation__weight_loss_by_status={"missing": 0.0}))
    assert res["reason_codes"] == ["LOW_SCOREABILITY"]
    assert rule(res, "G5")["detail"].endswith("largest loss: none")


# G9: strict CI above the null, finite numbers, explicitly single-ancestry units.

@pytest.mark.parametrize("ci_lower, ci_upper, null, informative", [
    (1.0, 1.3, 1.0, False),              # lower bound ON the null: not above it (strict >)
    (0.5, 0.6, 0.5, False),
    (1.0001, 1.3, 1.0, True),
    (float("inf"), float("inf"), 1.0, False),
    (float("nan"), 1.3, 1.0, False),
    (1.5, 1.2, 1.0, False),              # lower > upper is not a confidence interval
    (True, 2, False, False),             # booleans are not numbers
])
def test_g9_metric_must_be_a_finite_ci_strictly_above_the_null(ci_lower, ci_upper, null, informative):
    gi = base_input()
    gi["evaluation"]["units"][0]["metrics"] = [{"name": "OR", "ci_lower": ci_lower, "ci_upper": ci_upper,
                                                "null": null}]
    res = run(gi)
    assert (rule(res, "G9")["outcome"] == "pass") is informative, rule(res, "G9")
    if not informative:
        assert res["primary_reason"]["code"] == "EVALUATION_NOT_INFORMATIVE"


@pytest.mark.parametrize("pooled", [MISSING, None, "false", 0])
def test_g9_unit_counts_only_when_explicitly_not_pooled(pooled):
    gi = base_input()
    if pooled is MISSING:
        del gi["evaluation"]["units"][0]["pooled"]
    else:
        gi["evaluation"]["units"][0]["pooled"] = pooled
    assert run(gi)["primary_reason"]["code"] == "NO_RELEVANT_EVALUATION"


# G10: sex composition must be known to include the person's sex.

@pytest.mark.parametrize("sex, percent_male, outcome, code", [
    # known compatible, including the single-sex boundaries
    ("female", [48.0], "pass", None),
    ("male", [48.0], "pass", None),
    ("female", [0.0], "pass", None),
    ("male", [100.0], "pass", None),
    ("female", [99.9], "pass", None),
    ("male", [0.1], "pass", None),
    # known incompatible: only the other sex
    ("female", [100.0], "fail", "SEX_POPULATION_MISMATCH"),
    ("female", [100], "fail", "SEX_POPULATION_MISMATCH"),
    ("male", [0.0], "fail", "SEX_POPULATION_MISMATCH"),
    ("male", [0.0, 0], "fail", "SEX_POPULATION_MISMATCH"),
    # unknown composition is never assumed to include the person
    ("female", [MISSING], "fail", "SEX_EVALUATION_UNVERIFIED"),
    ("male", [None], "fail", "SEX_EVALUATION_UNVERIFIED"),
    ("female", ["50"], "fail", "SEX_EVALUATION_UNVERIFIED"),
    ("female", [float("nan")], "fail", "SEX_EVALUATION_UNVERIFIED"),
    ("male", [150.0], "fail", "SEX_EVALUATION_UNVERIFIED"),
    ("male", [-1.0], "fail", "SEX_EVALUATION_UNVERIFIED"),
    ("male", [0.0, None], "fail", "SEX_EVALUATION_UNVERIFIED"),  # the unknown unit is not a mismatch either
    ("female", [None, 40.0], "pass", None),                      # one known compatible unit suffices
    # unknown target sex: only an evaluation including both sexes covers the person
    (None, [48.0], "pass", None),
    (None, [0.0], "fail", "SEX_EVALUATION_UNVERIFIED"),
    (None, [100.0], "fail", "SEX_EVALUATION_UNVERIFIED"),
    (None, [None], "fail", "SEX_EVALUATION_UNVERIFIED"),
    (None, [0.0, 100.0], "fail", "SEX_EVALUATION_UNVERIFIED"),
])
def test_g10_evaluation_sex_composition(sex, percent_male, outcome, code):
    res = run(set_units(mutate(person__sex=sex), *percent_male))
    g10 = rule(res, "G10")
    assert g10["outcome"] == outcome, g10
    if code is None:
        assert res["status"] == "SUPPORTED"
    else:
        assert res["status"] == "RAW_ONLY" and g10["codes"] == [code]
        assert res["allowed_claims"]["percentile"] is False


def test_g10_not_applicable_without_an_informative_evaluation():
    res = run(mutate(placement__placement="SAS", reference_distribution__reference_group="SAS"))
    assert rule(res, "G10")["outcome"] == "not_applicable"


# G11: a real reference distribution for the placed group, on a non-empty variant intersection.

@pytest.mark.parametrize("key, value", [
    ("reference_n", MISSING), ("reference_n", None), ("reference_n", 0), ("reference_n", -5),
    ("reference_n", "502"), ("reference_n", 502.0), ("reference_n", True),
    ("n_intersection", MISSING), ("n_intersection", None), ("n_intersection", 0), ("n_intersection", -1),
    ("n_intersection", "290"), ("n_intersection", 290.5), ("n_intersection", True),
    ("reference_group", "AFR"), ("reference_group", None), ("reference_group", MISSING),
    ("available", "true"), ("available", 1),
])
def test_g11_rejects_missing_or_invalid_reference_structure(key, value):
    res = run(with_field("reference_distribution", key, value))
    assert rule(res, "G11")["outcome"] == "fail", rule(res, "G11")
    assert res["status"] == "RAW_ONLY" and res["primary_reason"]["code"] == "REFERENCE_DISTRIBUTION_UNAVAILABLE"
    assert rule(res, "G12")["outcome"] == "not_applicable"


@pytest.mark.parametrize("reference_n, n_intersection", [(1, 1), (502, 290)])
def test_g11_passes_positive_integer_sizes(reference_n, n_intersection):
    res = run(mutate(reference_distribution__reference_n=reference_n,
                     reference_distribution__n_intersection=n_intersection))
    assert rule(res, "G11")["outcome"] == "pass" and res["status"] == "SUPPORTED"
    assert {"reference_distribution.reference_n", "reference_distribution.n_intersection"} <= set(res["evidence_used"])


def test_g11_invalid_size_of_an_available_distribution_is_missing_evidence():
    res = run(mutate(reference_distribution__reference_n=None))
    assert "reference_distribution.reference_n" in res["evidence_missing"]


# G12: only an explicit false passes.

@pytest.mark.parametrize("value, outcome, code", [
    (False, "pass", None),
    (True, "fail", "REFERENCE_SENSITIVE"),
    (MISSING, "fail", "REFERENCE_SENSITIVITY_UNVERIFIED"),
    (None, "fail", "REFERENCE_SENSITIVITY_UNVERIFIED"),
    ("false", "fail", "REFERENCE_SENSITIVITY_UNVERIFIED"),
    (0, "fail", "REFERENCE_SENSITIVITY_UNVERIFIED"),
    ([], "fail", "REFERENCE_SENSITIVITY_UNVERIFIED"),
])
def test_g12_reference_sensitive_must_be_explicitly_false(value, outcome, code):
    gi = with_field("reference_distribution", "reference_sensitive", value)
    if value is True:
        gi["reference_distribution"]["reference_sensitive_pairs"] = [["FIN", "TSI"]]
    res = run(gi)
    assert rule(res, "G12")["outcome"] == outcome, rule(res, "G12")
    if code is None:
        assert res["status"] == "SUPPORTED" and res["allowed_claims"]["percentile"] is True
    else:
        assert res["status"] == "RAW_ONLY" and res["primary_reason"]["code"] == code
        assert res["allowed_claims"]["percentile"] is False
    if code == "REFERENCE_SENSITIVITY_UNVERIFIED":
        assert "reference_distribution.reference_sensitive" in res["evidence_missing"]


# Numeric domains: finite numbers only.

@pytest.mark.parametrize("value, expected", [
    (0.0, True), (1.0, True), (0, True), (1, True), (0.5, True),
    (-1e-12, False), (1.0000001, False),
    (float("nan"), False), (float("inf"), False), (float("-inf"), False),
    (10 ** 400, False), (-(10 ** 400), False),
    (True, False), (False, False), ("0.5", False), (None, False),
])
def test_number_in_accepts_only_finite_numbers_in_range(value, expected):
    assert gate._number_in(value, 0.0, 1.0) is expected


@pytest.mark.parametrize("changes", [
    {"scoreability__r": float("inf")}, {"scoreability__r": float("-inf")}, {"scoreability__r": float("nan")},
    {"placement__placement_stability": float("inf")}, {"placement__placement_stability": float("-inf")},
    {"harmonisation__allele_mismatch_fraction": float("inf")},
    {"harmonisation__weight_loss_by_status": {"missing": float("nan")}},
    {"harmonisation__weight_loss_by_status": {"missing": float("inf")}},
])
def test_non_finite_numbers_are_invalid_input(changes):
    assert _invalid(run(mutate(**changes)))


@pytest.mark.parametrize("r", [-1.0, 1.0, -1, 1])
def test_scoreability_r_domain_boundaries_are_valid(r):
    assert "INVALID_GATE_INPUT" not in run(mutate(scoreability__r=r))["reason_codes"]


# ---- fail-closed input validation -----------------------------------------------------------------------------------

def _invalid(res: dict) -> bool:
    return res["status"] == "ABSTAIN" and res["reason_codes"] == ["INVALID_GATE_INPUT"]


@pytest.mark.parametrize("bad", [
    {"schema": "v1"}, {"placement": None}, {"harmonisation": {"n_variants": -1}},
])
def test_invalid_input_abstains(bad):
    gi = base_input()
    for k, v in bad.items():
        if isinstance(v, dict) and isinstance(gi.get(k), dict):
            gi[k].update(v)
        else:
            gi[k] = v
    assert _invalid(run(gi))


def test_out_of_range_numbers_are_rejected():
    assert _invalid(run(mutate(scoreability__r=1.7)))


def test_null_mismatch_fraction_with_located_variants_is_invalid_not_a_crash():
    assert _invalid(run(mutate(harmonisation__allele_mismatch_fraction=None)))


def test_null_mismatch_fraction_is_valid_only_when_nothing_located():
    gi = mutate(harmonisation__allele_mismatch_fraction=None, harmonisation__n_located=0,
                harmonisation__n_matched=0)
    res = run(gi)
    assert "INVALID_GATE_INPUT" not in res["reason_codes"]
    assert rule(res, "G4")["outcome"] == "not_applicable"
    assert _invalid(run(mutate(harmonisation__n_located=0, harmonisation__n_matched=0)))  # 0.0 given, must be null


@pytest.mark.parametrize("changes", [
    {"harmonisation__allele_mismatch_fraction": -0.1},
    {"harmonisation__allele_mismatch_fraction": 1.5},
    {"harmonisation__allele_mismatch_fraction": "0.01"},
    {"harmonisation__allele_mismatch_fraction": float("nan")},
    {"placement__placement_stability": -0.2},
    {"placement__placement_stability": 1.01},
    {"scoreability__r": -1.5},
    {"harmonisation__weight_loss_by_status": {"missing": -0.3}},
    {"score_file__n_parse_problems": -1},
    {"candidate__sex_specific": "both"},
    {"reference_distribution": ["not", "an", "object"]},
    {"evaluation__units": ["not a unit"]},
    {"evaluation__units": [{"code": "EUR", "metrics": "not a list"}]},
])
def test_out_of_domain_or_malformed_inputs_abstain(changes):
    assert _invalid(run(mutate(**changes)))


@pytest.mark.parametrize("not_an_object", [[1, 2], "text", 42, None])
def test_non_object_input_abstains(not_an_object):
    assert _invalid(run(not_an_object))


def test_gate_never_raises_on_mutated_inputs():
    rng = random.Random(20260926)
    weird = [None, -1, 2.5, "x", [], {}, True, float("nan"), float("inf"), 10 ** 400, {"a": [1]}, [None]]
    base = base_input()
    paths = [(blk, k) for blk, v in base.items() if isinstance(v, dict) for k in v]
    for _ in range(400):
        gi = copy.deepcopy(base)
        for blk, k in rng.sample(paths, 3):
            gi[blk][k] = rng.choice(weird)
        res = run(gi)  # must not raise
        assert res["status"] in ("SUPPORTED", "RAW_ONLY", "ABSTAIN")
        assert res["allowed_claims"]["absolute_risk"] is False


def test_invalid_input_says_which_schema_would_be_accepted():
    res = run({"schema": "something else"})
    change = res["what_would_change_result"][0]["change"]
    assert gate.INPUT_SCHEMA in change
    assert "prsguard" not in change.lower(), "any producer of the documented schema is acceptable"


# ---- determinism, tamper evidence, calibration ------------------------------------------------------------------------

def test_tampered_input_digest_is_rejected():
    gi = base_input()
    gi["input_digest"] = gate.canonical_digest(gi)
    assert run(gi)["status"] == "SUPPORTED"
    gi["scoreability"]["r"] = 0.99
    assert _invalid(run(gi))


def test_deterministic_and_order_independent():
    gi = base_input()
    a = run(copy.deepcopy(gi))
    b = run(json.loads(json.dumps(gi, sort_keys=True)))
    assert json.dumps(a, sort_keys=True) == json.dumps(b, sort_keys=True)


def test_config_validation(tmp_path):
    bad = tmp_path / "c.yaml"
    bad.write_text("calibration_version: x\nscoreability: {r_min: 1.5}\n")
    with pytest.raises(ValueError):
        gate.load_config(bad)


def test_shipped_config_is_canonical_and_alternatives_are_labelled(tmp_path):
    assert CFG["_canonical"] is True
    assert gate.evaluate(base_input(), CFG)["provenance"]["config_canonical"] is True
    cfg = gate.load_config(alt_config(tmp_path, "0.80"))
    res = gate.evaluate(base_input(), cfg)
    assert res["provenance"]["config_canonical"] is False
    assert "NON-CANONICAL CALIBRATION" in gate.render_report(res)


def test_gate_has_no_side_channels():
    src = SCRIPT.read_text(encoding="utf-8")
    body = src.split('"""', 2)[2]
    for forbidden in ("requests", "urllib", "random", "datetime.now", "pgscatalog", "socket", "subprocess"):
        assert forbidden not in body, forbidden
    assert not re.search(r"^\s*(import|from)\s+prsguard\b", src, re.M), "no runtime dependency on PRSGuard"


def test_decision_carries_both_disclaimers():
    text = run(base_input())["disclaimer"]
    assert "not a clinical recommendation" in text and "discrimination or calibration" in text
    assert CLAWBIO_DISCLAIMER in text


def test_markdown_rule_trace_escapes_table_breaking_text():
    gi = mutate(placement__status="INTERMEDIATE", placement__placement=None,
                placement__detail="between | two clouds\nsecond line")
    md = gate.render_report(run(gi))
    g8 = [line for line in md.splitlines() if line.startswith("| G8 ")]
    assert len(g8) == 1 and "\\|" in g8[0] and "second line" in g8[0]


# ---- synthetic demo data -------------------------------------------------------------------------------------------------

def test_demo_ships_exactly_the_documented_synthetic_inputs():
    assert sorted(p.name for p in EXAMPLES.glob("*.gate_input.json")) == sorted(DEMO_EXPECTED)


@pytest.mark.parametrize("name", sorted(DEMO_EXPECTED))
def test_demo_inputs_are_labelled_synthetic_and_carry_no_identifiers(name):
    path = EXAMPLES / name
    text = path.read_text(encoding="utf-8")
    data = json.loads(text)
    assert data["demo_label"].startswith("SYNTHETIC"), "demo inputs must say they are synthetic"
    assert not re.search(r"\b(HG|NA)\d{5}\b", text), "no 1000 Genomes (or any real) sample identifiers"
    assert not re.fullmatch(r"PGS\d{6}", data["candidate"]["pgs_id"]), "must not pose as a real PGS Catalog record"
    assert path.stat().st_size < 256 * 1024, "small enough to ship in the wheel"


@pytest.mark.parametrize("name, expected", sorted(DEMO_EXPECTED.items()))
def test_demo_inputs_give_the_documented_decisions(name, expected):
    gi = json.loads((EXAMPLES / name).read_text(encoding="utf-8"))
    assert gate.validate_input(gi) == [], "demo input must be valid with an intact input_digest"
    assert "input_digest" in gi
    res = run(gi)
    assert (res["status"], res["reason_codes"]) == expected


# ---- CLI and output contract --------------------------------------------------------------------------------------------

def test_no_arguments_exits_nonzero():
    assert run_cli().returncode != 0


def test_output_without_input_or_demo_exits_nonzero(tmp_path):
    proc = run_cli("--output", tmp_path / "out")
    assert proc.returncode != 0
    assert "Traceback" not in proc.stderr


def test_missing_input_file_is_a_usage_error_not_a_traceback(tmp_path):
    proc = run_cli("--input", tmp_path / "does_not_exist.json", "--output", tmp_path / "out")
    assert proc.returncode != 0
    assert "Traceback" not in proc.stderr
    assert "does_not_exist.json" in proc.stderr


@pytest.mark.parametrize("content", ["not json at all", "[1, 2, 3]", "", "﻿{\"schema\": 1}"])
def test_unreadable_or_non_object_input_file_fails_closed(tmp_path, content):
    bad = tmp_path / "bad.gate_input.json"
    bad.write_text(content, encoding="utf-8")
    proc = run_cli("--input", bad, "--output", tmp_path / "out")
    assert proc.returncode == 0, proc.stderr
    assert "Traceback" not in proc.stderr
    decision = load_result(tmp_path / "out")["data"]["decisions"][0]
    assert _invalid(decision)
    assert decision["allowed_claims"] == {"raw_score": False, "standardized_score": False, "percentile": False,
                                          "absolute_risk": False}


def test_single_input_writes_the_clawbio_envelope(tmp_path):
    out = tmp_path / "out"
    proc = run_cli("--input", SUPPORTED_EXAMPLE, "--output", out)
    assert proc.returncode == 0, proc.stderr
    result = load_result(out)
    assert result["skill"] == "prs-applicability-gate"
    assert result["version"] == gate.GATE_VERSION
    assert result["input_checksum"] == "sha256:" + sha256_of(SUPPORTED_EXAMPLE)
    summary = result["summary"]
    assert summary["mode"] == "input" and summary["synthetic_demo"] is False
    assert summary["config_canonical"] is True
    assert summary["config_sha256"] == "sha256:" + sha256_of(SHIPPED_CONFIG)
    assert summary["calibration_version"] == CFG["calibration_version"]
    assert [d["status"] for d in summary["decisions"]] == ["SUPPORTED"]
    expected = gate.evaluate(json.loads(SUPPORTED_EXAMPLE.read_text(encoding="utf-8")), CFG)
    assert result["data"]["decisions"] == [expected], "the CLI must return exactly the pure gate decision"
    assert CLAWBIO_DISCLAIMER in result["disclaimer"]


def test_demo_evaluates_the_three_synthetic_cases(tmp_path):
    out = tmp_path / "demo"
    proc = run_cli("--demo", "--output", out)
    assert proc.returncode == 0, proc.stderr
    summary = load_result(out)["summary"]
    assert summary["mode"] == "demo" and summary["synthetic_demo"] is True
    got = {d["input"]: (d["status"], d["reason_codes"]) for d in summary["decisions"]}
    assert got == DEMO_EXPECTED
    assert summary["status_counts"] == {"SUPPORTED": 1, "RAW_ONLY": 1, "ABSTAIN": 1}
    report = (out / "report.md").read_text(encoding="utf-8")
    assert "SYNTHETIC" in report


def test_report_has_rule_trace_claims_and_disclaimers(tmp_path):
    out = tmp_path / "demo"
    assert run_cli("--demo", "--output", out).returncode == 0
    report = (out / "report.md").read_text(encoding="utf-8")
    assert CLAWBIO_DISCLAIMER in report
    assert "not a clinical recommendation" in report
    for status in ("SUPPORTED", "RAW_ONLY", "ABSTAIN"):
        assert status in report
    for rid in range(1, 13):
        assert report.count(f"| G{rid} ") == 3, f"G{rid} must appear once per decision"
    assert report.count("Absolute risk: never") == 3


def _parse_output_contract(skill_md: Path) -> list[str]:
    """Files promised in the SKILL.md '## Output Structure' tree ('(optional)' entries excluded)."""
    text = skill_md.read_text(encoding="utf-8")
    m = re.search(r"##\s*Output Structure\s*\n+```[^\n]*\n(.*?)\n```", text, re.S)
    if not m:
        return []
    files, parents = [], {}
    for raw in m.group(1).splitlines():
        if not raw.strip():
            continue
        parts = re.split(r"\s+#", raw, maxsplit=1)
        entry, comment = parts[0], (parts[1] if len(parts) > 1 else "")
        mm = re.match(r"^([\s│├└─]*)(.*)$", entry)
        prefix, name = mm.group(1), mm.group(2).strip()
        depth = len(prefix) // 4
        if not name or depth == 0:
            continue
        if name.endswith("/"):
            parents[depth] = name.rstrip("/")
            for d in [k for k in parents if k > depth]:
                del parents[d]
            continue
        if "optional" in comment.lower():
            continue
        rel = "/".join(parents[d] for d in sorted(parents) if d < depth)
        files.append(f"{rel}/{name}" if rel else name)
    return files


def test_documented_outputs_are_produced(tmp_path):
    promised = _parse_output_contract(SKILL / "SKILL.md")
    assert "reproducibility/commands.sh" in promised, "Output Structure must list the reproducibility bundle"
    assert run_cli("--demo", "--output", tmp_path).returncode == 0
    missing = [p for p in promised if not (tmp_path / p).exists()]
    assert not missing, f"SKILL.md Output Structure promises files the skill did not write: {missing}"


def test_reproducibility_bundle_verifies(tmp_path):
    out = tmp_path / "out"
    assert run_cli("--input", SUPPORTED_EXAMPLE, "--output", out).returncode == 0
    repro = out / "reproducibility"
    commands = (repro / "commands.sh").read_text(encoding="utf-8")
    assert "prs_applicability_gate.py" in commands and "--input" in commands and "--config" not in commands
    assert os.access(repro / "commands.sh", os.X_OK)
    assert "pyyaml" in (repro / "environment.yml").read_text(encoding="utf-8")
    lines = (repro / "checksums.sha256").read_text(encoding="utf-8").splitlines()
    labels = {line.split("  ", 1)[1] for line in lines}
    assert {"report.md", "result.json"} <= labels
    for line in lines:
        digest, label = line.split("  ", 1)
        assert not label.startswith("/"), "labels are relative to the output directory"
        assert sha256_of(out / label) == digest, label


def test_repeated_runs_are_byte_identical_and_path_free(tmp_path):
    a, b = tmp_path / "a", tmp_path / "b"
    assert run_cli("--demo", "--output", a).returncode == 0
    assert run_cli("--demo", "--output", b).returncode == 0
    for name in ("result.json", "report.md"):
        text = (a / name).read_text(encoding="utf-8")
        assert text == (b / name).read_text(encoding="utf-8"), f"{name} differs between identical runs"
        assert str(tmp_path) not in text and str(SKILL) not in text, f"{name} leaks a local path"


def test_rerun_into_existing_output_warns(tmp_path):
    out = tmp_path / "out"
    assert run_cli("--demo", "--output", out).returncode == 0
    proc = run_cli("--demo", "--output", out)
    assert proc.returncode == 0
    assert "overwrit" in proc.stderr.lower()


def test_direct_cli_non_canonical_config_is_explicit_and_labelled(tmp_path):
    out = tmp_path / "out"
    proc = run_cli("--demo", "--output", out, "--config", alt_config(tmp_path))
    assert proc.returncode == 0, proc.stderr
    assert "NON-CANONICAL CALIBRATION" in proc.stderr
    result = load_result(out)
    assert result["summary"]["config_canonical"] is False
    assert all(d["provenance"]["config_canonical"] is False for d in result["data"]["decisions"])
    assert "NON-CANONICAL CALIBRATION" in (out / "report.md").read_text(encoding="utf-8")
    assert "--config" in (out / "reproducibility" / "commands.sh").read_text(encoding="utf-8")


def test_invalid_config_is_a_clean_error(tmp_path):
    bad = tmp_path / "bad.yaml"
    bad.write_text("calibration_version: x\nscoreability: {r_min: 1.5}\n", encoding="utf-8")
    proc = run_cli("--demo", "--output", tmp_path / "out", "--config", bad)
    assert proc.returncode == 2, "SKILL.md: an unusable --config is a usage error (exit 2)"
    assert "Traceback" not in proc.stderr and "config" in proc.stderr
    assert not (tmp_path / "out" / "result.json").exists(), "no decision is written without a usable config"


@pytest.mark.parametrize("content", ["calibration_version: x\nscoreability: {r_min: 1.5}\n", None])
def test_main_exits_2_on_unusable_config_without_writing(tmp_path, content, capsys):
    config = tmp_path / "calibration.yaml"
    if content is not None:
        config.write_text(content, encoding="utf-8")
    with pytest.raises(SystemExit) as exc:
        gate.main(["--demo", "--output", str(tmp_path / "out"), "--config", str(config)])
    assert exc.value.code == 2
    assert "cannot use calibration config" in capsys.readouterr().err
    assert not (tmp_path / "out").exists()


@pytest.mark.parametrize("token", ["Infinity", "-Infinity", "NaN", "1" + "0" * 400])
def test_non_finite_or_huge_json_numbers_fail_closed_from_the_cli(tmp_path, token):
    text = SUPPORTED_EXAMPLE.read_text(encoding="utf-8").replace('"r": 0.962', f'"r": {token}')
    assert token in text
    bad = tmp_path / "bad.gate_input.json"
    bad.write_text(text, encoding="utf-8")
    proc = run_cli("--input", bad, "--output", tmp_path / "out")
    assert proc.returncode == 0 and "Traceback" not in proc.stderr, proc.stderr
    assert _invalid(load_result(tmp_path / "out")["data"]["decisions"][0])


# ---- ClawBio integration ----------------------------------------------------------------------------------------------

def test_runner_registration_forwards_no_extra_flags():
    from clawbio.cli import SKILLS

    entry = SKILLS["prs-gate"]
    assert Path(entry["script"]).resolve() == SCRIPT
    assert entry["demo_args"] == ["--demo"]
    assert entry["allowed_extra_flags"] == set()
    assert "--config" not in entry.get("allowed_extra_flags_without_values", set())


def test_runner_demo(tmp_path):
    from clawbio.cli import run_skill

    res = run_skill("prs-gate", demo=True, output_dir=str(tmp_path / "out"))
    assert res["success"], res["stderr"]
    statuses = {d["input"]: d["status"] for d in res["skill_result_json"]["summary"]["decisions"]}
    assert statuses == {name: exp[0] for name, exp in DEMO_EXPECTED.items()}


def test_runner_cannot_replace_the_calibration(tmp_path):
    """`clawbio.py run prs-gate` is the agent-facing interface: --config must be dropped, not forwarded."""
    from clawbio.cli import run_skill

    out = tmp_path / "out"
    res = run_skill("prs-gate", input_path=str(SUPPORTED_EXAMPLE), output_dir=str(out),
                    extra_args=["--config", str(alt_config(tmp_path))])
    assert res["success"], res["stderr"]
    summary = load_result(out)["summary"]
    assert summary["config_canonical"] is True
    assert summary["config_sha256"] == "sha256:" + sha256_of(SHIPPED_CONFIG)
    assert [d["status"] for d in summary["decisions"]] == ["SUPPORTED"], "r_min 0.99 would have made this ABSTAIN"


def test_catalog_lists_the_skill_with_its_alias():
    catalog = json.loads((REPO_ROOT / "skills" / "catalog.json").read_text(encoding="utf-8"))
    entry = next(s for s in catalog["skills"] if s["name"] == "prs-applicability-gate")
    assert entry["cli_alias"] == "prs-gate"
    assert entry["demo_command"] == "python clawbio.py run prs-gate --demo"
    assert entry["has_tests"] and entry["has_script"] and entry["license"] == "MIT"
    assert entry["version"] == gate.GATE_VERSION


# ---- SKILL.md conformance (AGENTS.md checklist) ------------------------------------------------------------------------

def _skill_md() -> tuple[dict, str]:
    text = (SKILL / "SKILL.md").read_text(encoding="utf-8")
    m = re.match(r"^---\n(.*?)\n---\n(.*)$", text, re.S)
    assert m, "SKILL.md must start with YAML frontmatter"
    return yaml.safe_load(m.group(1)), m.group(2)


def _section(body: str, title: str) -> str:
    m = re.search(rf"^## {re.escape(title)}\s*\n(.*?)(?=^## |\Z)", body, re.S | re.M)
    assert m, f"missing section: ## {title}"
    return m.group(1)


def test_skill_md_frontmatter_conforms():
    front, _ = _skill_md()
    meta = front["metadata"]
    assert front["name"] == SKILL.name
    assert front["license"] == "MIT"
    assert isinstance(front["description"], str) and "\n" not in front["description"].strip()
    assert re.fullmatch(r"\d+\.\d+\.\d+", meta["version"]) and meta["version"] == gate.GATE_VERSION
    assert meta["author"] == "Rinat Rizvanov"
    assert all(i.get("format") and "required" in i for i in meta["inputs"])
    assert all(o.get("format") for o in meta["outputs"])
    assert len(meta["openclaw"]["trigger_keywords"]) >= 3
    for demo in meta["demo_data"]:
        assert (SKILL / demo["path"]).exists(), demo["path"]


def test_skill_md_sections_conform():
    _, body = _skill_md()
    trigger = _section(body, "Trigger")
    assert "Fire this skill when" in trigger and "Do NOT fire when" in trigger
    assert "one task" in _section(body, "Scope").lower()
    workflow = _section(body, "Workflow")
    assert re.search(r"^1\. ", workflow, re.M) and re.search(r"^4\. ", workflow, re.M)
    assert "```" in _section(body, "Example Output")
    assert len(re.findall(r"^- \*\*", _section(body, "Gotchas"), re.M)) >= 3
    assert "not a medical device" in _section(body, "Safety")
    assert "must NOT" in _section(body, "Agent Boundary")
    assert "Chaining partners" in _section(body, "Integration with Bio Orchestrator")
    assert "Review cadence" in _section(body, "Maintenance")
    assert len((SKILL / "SKILL.md").read_text(encoding="utf-8").splitlines()) < 500
