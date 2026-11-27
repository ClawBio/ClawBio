"""Tests for target-validation-scorer — evidence-grounded target validation."""

import json
import sys
from pathlib import Path


sys.path.insert(0, str(Path(__file__).resolve().parent.parent))

from target_validation_scorer import (
    score_disease_association,
    score_druggability,
    score_chemical_matter,
    score_clinical_precedent,
    score_structural_data,
    compute_safety_penalty,
    decide,
    build_rationale,
    validate_target,
    generate_report,
    DEMO_DATA,
)


# ---------------------------------------------------------------------------
# score_disease_association
# ---------------------------------------------------------------------------

class TestScoreDiseaseAssociation:
    def test_none_input(self):
        r = score_disease_association(None)
        assert r["score"] is None
        assert r["confidence"] == "low"

    def test_high_score(self):
        r = score_disease_association({"overall_score": 0.8})
        assert r["score"] == 20
        assert r["confidence"] == "high"

    def test_medium_score(self):
        r = score_disease_association({"overall_score": 0.5})
        assert r["score"] == 10
        assert r["confidence"] == "medium"

    def test_low_score(self):
        r = score_disease_association({"overall_score": 0.2})
        assert r["score"] == 0
        assert r["confidence"] == "low"

    def test_boundary_0_7(self):
        r = score_disease_association({"overall_score": 0.7})
        assert r["score"] == 20

    def test_boundary_0_4(self):
        r = score_disease_association({"overall_score": 0.4})
        assert r["score"] == 10


# ---------------------------------------------------------------------------
# score_druggability
# ---------------------------------------------------------------------------

class TestScoreDruggability:
    def test_none_input(self):
        r = score_druggability(None)
        assert r["score"] is None

    def test_kinase_class(self):
        r = score_druggability({"target_class": "kinase", "known_ligands": 0})
        assert r["score"] == 10

    def test_unknown_class(self):
        r = score_druggability({"target_class": "transporter", "known_ligands": 0})
        assert r["score"] == 5

    def test_empty_class(self):
        r = score_druggability({"target_class": "", "known_ligands": 0})
        assert r["score"] == 0

    def test_ligand_tiers(self):
        assert score_druggability({"target_class": "", "known_ligands": 5})["score"] == 1
        assert score_druggability({"target_class": "", "known_ligands": 50})["score"] == 3
        assert score_druggability({"target_class": "", "known_ligands": 200})["score"] == 5
        assert score_druggability({"target_class": "", "known_ligands": 600})["score"] == 6

    def test_cap_at_20(self):
        r = score_druggability({"target_class": "kinase", "known_ligands": 1000})
        assert r["score"] == 16  # 10 + 6, under cap

    def test_gpcr_class(self):
        r = score_druggability({"target_class": "GPCR", "known_ligands": 0})
        assert r["score"] == 10


# ---------------------------------------------------------------------------
# score_chemical_matter
# ---------------------------------------------------------------------------

class TestScoreChemicalMatter:
    def test_none_input(self):
        r = score_chemical_matter(None)
        assert r["score"] is None

    def test_high_compounds(self):
        r = score_chemical_matter({"bioactive_compounds": 600, "compounds_sub_100nM": 60})
        assert r["score"] == 15  # 8 + 7

    def test_medium_compounds(self):
        r = score_chemical_matter({"bioactive_compounds": 200, "compounds_sub_100nM": 20})
        assert r["score"] == 11  # 6 + 5

    def test_low_compounds(self):
        r = score_chemical_matter({"bioactive_compounds": 5, "compounds_sub_100nM": 1})
        assert r["score"] == 4  # 1 + 3

    def test_zero_compounds(self):
        r = score_chemical_matter({"bioactive_compounds": 0, "compounds_sub_100nM": 0})
        assert r["score"] == 0


# ---------------------------------------------------------------------------
# score_clinical_precedent
# ---------------------------------------------------------------------------

class TestScoreClinicalPrecedent:
    def test_none_input(self):
        r = score_clinical_precedent(None)
        assert r["score"] is None

    def test_phase_1_plus(self):
        r = score_clinical_precedent({"max_phase": 2})
        assert r["score"] == 20
        assert r["confidence"] == "high"

    def test_trials_but_no_phase(self):
        r = score_clinical_precedent({"max_phase": 0, "trials_active": 3})
        assert r["score"] == 10

    def test_nothing(self):
        r = score_clinical_precedent({"max_phase": 0, "trials_active": 0, "trials_completed": 0})
        assert r["score"] == 0


# ---------------------------------------------------------------------------
# score_structural_data
# ---------------------------------------------------------------------------

class TestScoreStructuralData:
    def test_none_input(self):
        r = score_structural_data(None)
        assert r["score"] is None

    def test_full_structural(self):
        r = score_structural_data({
            "pdb_structures": 5,
            "co_crystal_ligands": True,
            "best_resolution_A": 1.8,
        })
        assert r["score"] == 20
        assert r["confidence"] == "high"

    def test_pdb_only(self):
        r = score_structural_data({
            "pdb_structures": 1,
            "co_crystal_ligands": False,
            "best_resolution_A": 3.5,
        })
        assert r["score"] == 10

    def test_alphafold_only(self):
        r = score_structural_data({
            "pdb_structures": 0,
            "alphafold_available": True,
        })
        assert r["score"] == 10

    def test_nothing(self):
        r = score_structural_data({
            "pdb_structures": 0,
            "alphafold_available": False,
        })
        assert r["score"] == 0


# ---------------------------------------------------------------------------
# compute_safety_penalty
# ---------------------------------------------------------------------------

class TestComputeSafetyPenalty:
    def test_empty_signals(self):
        total, flags = compute_safety_penalty([])
        assert total == 0
        assert flags == []

    def test_none_signals(self):
        total, flags = compute_safety_penalty(None)
        assert total == 0

    def test_single_penalty(self):
        total, flags = compute_safety_penalty([
            {"risk": "hepatotoxicity", "penalty": -10, "tier": "T2", "evidence": "rat study"}
        ])
        assert total == -10
        assert len(flags) == 1
        assert flags[0]["risk"] == "hepatotoxicity"

    def test_multiple_penalties_accumulate(self):
        total, _ = compute_safety_penalty([
            {"risk": "a", "penalty": -5},
            {"risk": "b", "penalty": -8},
        ])
        assert total == -13


# ---------------------------------------------------------------------------
# decide
# ---------------------------------------------------------------------------

class TestDecide:
    def test_go(self):
        decision, _ = decide(80)
        assert decision == "GO"

    def test_conditional_go(self):
        decision, _ = decide(60)
        assert decision == "CONDITIONAL_GO"

    def test_review(self):
        decision, _ = decide(30)
        assert decision == "REVIEW"

    def test_no_go(self):
        decision, _ = decide(10)
        assert decision == "NO_GO"

    def test_boundary_75(self):
        decision, _ = decide(75)
        assert decision == "GO"

    def test_boundary_50(self):
        decision, _ = decide(50)
        assert decision == "CONDITIONAL_GO"

    def test_boundary_25(self):
        decision, _ = decide(25)
        assert decision == "REVIEW"

    def test_zero(self):
        decision, _ = decide(0)
        assert decision == "NO_GO"


# ---------------------------------------------------------------------------
# build_rationale
# ---------------------------------------------------------------------------

class TestBuildRationale:
    def test_includes_decision(self):
        sub = {"dim1": {"score": 18}, "dim2": {"score": 5}, "dim3": {"score": None}}
        text = build_rationale("TGFBR1", "IPF", sub, [], "GO")
        assert "GO" in text

    def test_includes_strong_and_weak(self):
        sub = {"strong_dim": {"score": 16}, "weak_dim": {"score": 3}}
        text = build_rationale("T", "D", sub, [], "REVIEW")
        assert "Strong evidence" in text
        assert "Weak evidence" in text

    def test_includes_safety_concerns(self):
        flags = [{"risk": "cardiotoxicity", "penalty": -15}]
        text = build_rationale("T", "D", {"d": {"score": 10}}, flags, "NO_GO")
        assert "cardiotoxicity" in text


# ---------------------------------------------------------------------------
# validate_target (integration)
# ---------------------------------------------------------------------------

class TestValidateTarget:
    def test_with_demo_data(self):
        data = json.loads(DEMO_DATA.read_text())
        result = validate_target(data)
        assert result["target"] == data["target"]
        assert "decision" in result
        assert result["decision"] in ("GO", "CONDITIONAL_GO", "REVIEW", "NO_GO")
        assert result["adjusted_score"] >= 0
        assert result["raw_score"] == sum(
            v["score"] for v in result["sub_scores"].values() if v["score"] is not None
        )

    def test_all_none_evidence(self):
        result = validate_target({"target": "X", "evidence": {}})
        assert result["raw_score"] == 0
        assert result["decision"] == "NO_GO"


# ---------------------------------------------------------------------------
# generate_report
# ---------------------------------------------------------------------------

class TestGenerateReport:
    def test_report_structure(self, tmp_path):
        data = json.loads(DEMO_DATA.read_text())
        result = validate_target(data)
        report_path = generate_report(result, tmp_path)
        text = report_path.read_text()
        assert "# Target Validation Report" in text
        assert "Scoring summary" in text
        assert "Sub-scores" in text
        assert "Rationale" in text
        assert "research and educational tool" in text

    def test_creates_output_dir(self, tmp_path):
        result = validate_target({"target": "T", "evidence": {}})
        out = tmp_path / "new_dir"
        path = generate_report(result, out)
        assert path.exists()
