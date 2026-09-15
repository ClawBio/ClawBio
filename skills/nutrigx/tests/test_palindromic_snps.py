"""
test_palindromic_snps.py — strand handling for palindromic SNPs in NutriGx Advisor.

Flipping a palindromic (A/T or C/G) genotype yields the other allele of the same
pair, so the flip always succeeds and turns homozygous reference into homozygous
risk. At rs9939609 (FTO, ref T, risk A) a TT call -- no risk alleles -- was
normalised to AA and scored 2. Three panel entries are palindromic.
"""

import json
import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).parent.parent))

from extract_genotypes import extract_snp_genotypes

PANEL = Path(__file__).parent.parent / "data" / "snp_panel.json"


def load_panel():
    with open(PANEL) as f:
        return json.load(f)


def test_palindromic_snps_are_not_strand_flipped():
    panel = load_panel()
    palindromic = [
        s for s in panel
        if frozenset([s["ref_allele"], s["risk_allele"]]) in
        (frozenset(["A", "T"]), frozenset(["C", "G"]))
    ]
    assert palindromic, "expected palindromic entries in the panel"

    for snp in palindromic:
        ref, risk = snp["ref_allele"], snp["risk_allele"]
        for genotype, expected in ((ref * 2, 0), (ref + risk, 1), (risk * 2, 2)):
            call = extract_snp_genotypes({snp["rsid"]: genotype}, panel)[snp["rsid"]]
            assert call["status"] == "found", f"{snp['rsid']} {genotype}: {call['status']}"
            assert call["risk_count"] == expected, (
                f"{snp['rsid']} ({ref}/{risk}) {genotype}: risk_count "
                f"{call['risk_count']}, expected {expected}"
            )


def test_non_palindromic_snps_still_strand_flip():
    """rs4988235 is C/T against G/A, so flipping is unambiguous and still needed."""
    panel = load_panel()
    for genotype in ("CC", "GG"):          # the same call on opposite strands
        call = extract_snp_genotypes({"rs4988235": genotype}, panel)["rs4988235"]
        assert call["risk_count"] == 2, f"{genotype} should give 2 risk alleles"
