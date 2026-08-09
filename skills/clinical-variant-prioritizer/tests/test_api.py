"""Red/green tests for the clinical-variant-prioritizer skill.

The skill screens a genotype dict against a curated clinical panel (OMIM-morbid,
ACMG-SF, Hereditary-Cancer loci) and prioritises carried variants by ClinVar
significance, gnomAD frequency, inheritance and zygosity, following the
pathogenicity-screening method of Corpas et al. 2021 (Front Genet 12:535123).
"""
import sys
from pathlib import Path

SKILL_DIR = Path(__file__).resolve().parent.parent
if str(SKILL_DIR) not in sys.path:
    sys.path.insert(0, str(SKILL_DIR))

from api import run  # noqa: E402


# Manuel Corpas (the son, PT00002A) real array calls at panel loci.
MANUEL = {
    "rs28941785": "CT",   # CTH c.200C>T heterozygous (real chip call)
    "rs1800562": "GG",    # HFE C282Y reference
    "rs1799945": "CC",    # HFE H63D reference
    "rs6025": "CC",       # F5 Leiden reference
}


def _by_gene(result, gene):
    return [f for f in result["findings"] if f["gene"] == gene]


def test_returns_expected_shape():
    r = run(MANUEL)
    assert r["skill"] == "clinical-variant-prioritizer"
    assert "summary" in r and "findings" in r and "headline" in r
    for key in ("panel_size", "loci_tested", "loci_carried", "reference"):
        assert key in r["summary"], f"missing summary.{key}"


def test_cth_flagged_as_heterozygous_carrier():
    r = run(MANUEL)
    cth = _by_gene(r, "CTH")
    assert len(cth) == 1
    assert cth[0]["zygosity"] == "het"
    # AR + heterozygous + recessive benign-spectrum trait => carrier, not actionable
    assert cth[0]["category"] == "carrier"
    assert cth[0]["genotype"] == "CT"


def test_hfe_reference_not_carried():
    r = run(MANUEL)
    hfe = [f for f in _by_gene(r, "HFE") if f["rsid"] == "rs1800562"][0]
    assert hfe["zygosity"] == "ref"
    assert hfe["category"] == "reference"


def test_no_actionable_pathogenic_for_manuel():
    r = run(MANUEL)
    assert r["summary"]["actionable"] == 0
    assert r["summary"]["affected"] == 0
    assert r["summary"]["carriers"] == 1  # CTH only


def test_homozygous_recessive_pathogenic_is_affected():
    # Homozygous for HFE C282Y (alt allele A) => affected category
    r = run({"rs1800562": "AA"})
    hfe = [f for f in _by_gene(r, "HFE") if f["rsid"] == "rs1800562"][0]
    assert hfe["zygosity"] == "hom"
    assert hfe["category"] == "affected"


def test_risk_factor_heterozygote_is_actionable():
    # Factor V Leiden heterozygote (risk-factor inheritance) => actionable flag
    r = run({"rs6025": "CT"})
    f5 = _by_gene(r, "F5")[0]
    assert f5["zygosity"] == "het"
    assert f5["category"] == "actionable"


def test_summary_counts_are_consistent():
    r = run(MANUEL)
    s = r["summary"]
    carried_cats = s["actionable"] + s["affected"] + s["carriers"] + s["uncertain"] + s["benign"]
    assert s["loci_carried"] == carried_cats
    assert s["loci_tested"] == s["loci_carried"] + s["reference"]
    # findings list only contains tested loci (carried + reference)
    assert len(r["findings"]) == s["loci_tested"]


def test_findings_ranked_actionable_before_reference():
    r = run({"rs1800562": "AA", "rs28941785": "CT", "rs1799945": "CC"})
    cats = [f["category"] for f in r["findings"]]
    # affected/actionable must come before reference in the ranked list
    assert cats.index("affected") < cats.index("reference")


def test_no_panel_loci_present():
    r = run({"rs9999999": "AA"})
    assert r["summary"]["loci_tested"] == 0
    assert r["summary"]["loci_carried"] == 0
    assert r["findings"] == []
