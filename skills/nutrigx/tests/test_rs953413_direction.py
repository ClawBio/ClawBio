"""
test_rs953413_direction.py — ELOVL2 rs953413: A is the lower-DHA allele.

Tanaka et al. 2009 (PMID 19148276): DHA falls GG > AG > AA in two cohorts, and
"the presence of the minor (A) allele was associated with higher EPA/DPA and
lower DHA". The panel previously scored G as the risk allele. G/A is the plus
strand in GRCh37 and GRCh38, and the SNP is not palindromic, so T/C calls on the
minus strand must resolve to the same counts.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from extract_genotypes import extract_snp_genotypes
from score_variants import compute_nutrient_risk_scores

PANEL = Path(__file__).parent.parent / "data" / "snp_panel.json"


def load_panel():
    with open(PANEL) as f:
        return json.load(f)


def _call(genotype):
    return extract_snp_genotypes({"rs953413": genotype}, load_panel())["rs953413"]


def test_panel_entry_direction_and_citation():
    entry = next(s for s in load_panel() if s["rsid"] == "rs953413")
    assert entry["ref_allele"] == "G"
    assert entry["risk_allele"] == "A"
    assert entry["pmid"] == "19148276"
    assert entry["effect_direction"] == "lower_epa_to_dha_conversion"


def test_plus_strand_counts():
    for genotype, expected in (("AA", 2), ("AG", 1), ("GA", 1), ("GG", 0)):
        call = _call(genotype)
        assert call["status"] == "found", genotype
        assert call["risk_count"] == expected, f"{genotype}: {call['risk_count']} != {expected}"


def test_minus_strand_counts():
    for genotype, expected in (("TT", 2), ("TC", 1), ("CT", 1), ("CC", 0)):
        call = _call(genotype)
        assert call["status"] == "found", genotype
        assert call["risk_count"] == expected, f"{genotype}: {call['risk_count']} != {expected}"


def test_gg_scores_lower_omega3_risk_than_aa():
    panel = load_panel()
    gg = compute_nutrient_risk_scores(extract_snp_genotypes({"rs953413": "GG"}, panel), panel)["omega3"]
    aa = compute_nutrient_risk_scores(extract_snp_genotypes({"rs953413": "AA"}, panel), panel)["omega3"]
    assert gg["score"] < aa["score"]
