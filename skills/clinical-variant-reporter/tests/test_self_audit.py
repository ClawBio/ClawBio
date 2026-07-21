"""Tests for the fail-closed self-audit (hard-abstain on invariant violation).

The anchor case is the CTH wrong-variant regression: a coordinate that resolves to
p.Ser403Ile while the caller asserted p.Thr67Ile must be abstained, not classified.
"""
from __future__ import annotations

import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
sys.path.insert(0, str(SKILL_DIR))

from acmg_engine import ClassifiedVariant, EvidenceCriterion, VariantEvidence
from self_audit import (
    ABSTAIN_LABEL,
    audit_classified,
    expected_from_record,
)
from clinical_variant_reporter import VcfRecord, run_classification


def _crit(code, direction, source="SIFT=deleterious", strength="supporting"):
    return EvidenceCriterion(code=code, triggered=True, strength=strength,
                             direction=direction, source=source, detail="test")


def _cv(gene="CTH", hgvsp="ENSP00000418447.2:p.Ser403Ile", criteria=None,
        is_missense=True):
    ev = VariantEvidence(chrom="chr1", pos=70439117, ref="G", alt="T",
                         gene=gene, hgvsp=hgvsp, consequence="missense_variant",
                         is_missense=is_missense)
    return ClassifiedVariant(evidence=ev, criteria=criteria or [], classification="Benign")


# ---- IDENTITY_MISMATCH (the CTH regression) --------------------------------

def test_identity_mismatch_protein_abstains():
    cv = _cv(hgvsp="ENSP00000418447.2:p.Ser403Ile")
    res = audit_classified(cv, {"gene": "CTH", "hgvsp": "p.Thr67Ile"})
    assert res.passed is False
    assert "IDENTITY_MISMATCH" in res.reason_codes()


def test_identity_match_protein_passes():
    cv = _cv(hgvsp="ENSP0:p.Thr67Ile")
    res = audit_classified(cv, {"gene": "CTH", "hgvsp": "p.Thr67Ile"})
    assert res.passed is True


def test_identity_mismatch_gene_abstains():
    cv = _cv(gene="BRCA2")
    res = audit_classified(cv, {"gene": "CTH"})
    assert res.passed is False
    assert "IDENTITY_MISMATCH" in res.reason_codes()


def test_protein_normalisation_handles_ter_and_star():
    cv = _cv(hgvsp="p.Lys1638Ter")
    assert audit_classified(cv, {"hgvsp": "p.Lys1638*"}).passed is True


def test_no_assertion_no_identity_check():
    """A bare coordinate with no asserted identity is not identity-audited."""
    cv = _cv(gene="CTH", hgvsp="ENSP0:p.Ser403Ile")
    assert audit_classified(cv, {}).passed is True


# ---- CONTRADICTORY_EVIDENCE ------------------------------------------------

def test_pp3_and_bp4_together_abstains():
    cv = _cv(criteria=[_crit("PP3", "pathogenic", "SIFT=deleterious"),
                       _crit("BP4", "benign", "PolyPhen=benign")])
    res = audit_classified(cv, {})
    assert res.passed is False
    assert "CONTRADICTORY_EVIDENCE" in res.reason_codes()


def test_pp3_alone_is_fine():
    cv = _cv(criteria=[_crit("PP3", "pathogenic", "SIFT=deleterious")])
    assert audit_classified(cv, {}).passed is True


# ---- MISSING_PROVENANCE ----------------------------------------------------

def test_triggered_criterion_without_source_abstains():
    cv = _cv(criteria=[_crit("PP3", "pathogenic", source="")])
    res = audit_classified(cv, {})
    assert res.passed is False
    assert "MISSING_PROVENANCE" in res.reason_codes()


def test_no_in_silico_data_marker_abstains():
    cv = _cv(criteria=[_crit("PP3", "pathogenic", source="No in silico data available")])
    assert audit_classified(cv, {}).passed is False


# ---- expected_from_record --------------------------------------------------

def test_expected_from_record_reads_info_and_id():
    rec = VcfRecord(chrom="chr1", pos=70415987, id="rs28941785", ref="C", alt="T",
                    qual=".", filt="PASS",
                    info={"GENE": "CTH", "EXPECTED_HGVSP": "p.Thr67Ile"})
    exp = expected_from_record(rec)
    assert exp["gene"] == "CTH"
    assert exp["hgvsp"] == "p.Thr67Ile"
    assert exp["rsid"] == "rs28941785"


# ---- end-to-end hard-abstain through run_classification --------------------

def test_run_classification_hard_abstains_on_identity_mismatch():
    """A demo record asserting a protein change that does not match the resolved
    evidence must come back abstained, not classified."""
    cache_key = "chr17:43045684:AG:A"  # BRCA1 demo variant (frameshift)
    chrom, pos, ref, alt = "chr17", 43045684, "AG", "A"
    rec = VcfRecord(chrom=chrom, pos=pos, id=".", ref=ref, alt=alt,
                    qual=".", filt="PASS",
                    info={"GENE": "BRCA1", "EXPECTED_HGVSP": "p.Thr67Ile"})
    out = run_classification([rec], demo=True)
    assert len(out) == 1
    assert out[0].classification == ABSTAIN_LABEL
    assert "IDENTITY_MISMATCH" in [v.code for v in out[0].audit_violations]
