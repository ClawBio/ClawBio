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
    out, _ = run_classification([rec], demo=True)
    assert len(out) == 1
    assert out[0].classification == ABSTAIN_LABEL
    assert "IDENTITY_MISMATCH" in [v.code for v in out[0].audit_violations]


# ---- ALLELE_IMBALANCE (low-fraction indel cluster regression) --------------

def _rec(gt="0/1", ad="30,4", alt="GTTTT", ref="G"):
    return VcfRecord(chrom="19", pos=39056319, id=".", ref=ref, alt=alt, qual="69.6",
                     filt="PASS", info={}, genotype=gt,
                     sample={"GT": gt, "AD": ad, "DP": "34"}, alt_index=1)


def test_low_fraction_het_abstains():
    # RYR1 cluster seen on a real exome: 8 calls in 17 bp, each on the same 4/34 reads
    res = audit_classified(_cv(gene="RYR1", hgvsp=""), {}, record=_rec(ad="30,4"))
    assert not res.passed
    assert "ALLELE_IMBALANCE" in res.reason_codes()


def test_balanced_het_passes():
    res = audit_classified(_cv(gene="RYR1", hgvsp=""), {}, record=_rec(ad="144,149"))
    assert "ALLELE_IMBALANCE" not in res.reason_codes()


def test_too_few_alt_reads_abstains_even_if_fraction_ok():
    res = audit_classified(_cv(gene="RYR1", hgvsp=""), {}, record=_rec(ad="3,3"))
    assert "ALLELE_IMBALANCE" in res.reason_codes()


def test_hom_alt_with_low_fraction_abstains():
    res = audit_classified(_cv(gene="RYR1", hgvsp=""), {}, record=_rec(gt="1/1", ad="20,10"))
    assert "ALLELE_IMBALANCE" in res.reason_codes()


def test_no_ad_field_is_not_judged():
    rec = _rec(); rec.sample = {"GT": "0/1"}
    res = audit_classified(_cv(gene="RYR1", hgvsp=""), {}, record=rec)
    assert "ALLELE_IMBALANCE" not in res.reason_codes()


def test_multiallelic_uses_this_alt_index():
    rec = _rec(gt="1/2", ad="2,30,4", alt="T"); rec.alt_index = 2
    res = audit_classified(_cv(gene="RYR1", hgvsp=""), {}, record=rec)
    assert "ALLELE_IMBALANCE" in res.reason_codes()   # second ALT has only 4 reads


def test_parse_vcf_keeps_format_fields(tmp_path):
    from clinical_variant_reporter import parse_vcf
    vcf = tmp_path / "t.vcf"
    vcf.write_text("##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tS1\n"
                   "19\t39056319\t.\tG\tGTTTT\t69.6\tPASS\t.\tGT:AD:DP:GQ\t0/1:30,4:34:77\n")
    rec = parse_vcf(vcf)[0]
    assert rec.sample["AD"] == "30,4" and rec.alt_index == 1
