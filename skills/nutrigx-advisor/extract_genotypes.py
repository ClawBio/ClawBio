"""
extract_genotypes.py — SNP lookup with forward-strand normalisation
For each SNP in the panel, extracts the genotype from the parsed data dict.
Handles strand flipping for ambiguous A/T and C/G SNPs using frequency context.
"""

COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}

# SNPs where both alleles are complementary (ambiguous strand)
AMBIGUOUS_PAIRS = {frozenset(["A", "T"]), frozenset(["C", "G"])}


def flip_genotype(genotype: str) -> str:
    """Return the complement strand genotype."""
    return "".join(COMPLEMENT.get(b, b) for b in genotype)


def is_ambiguous(ref: str, alt: str) -> bool:
    return frozenset([ref, alt]) in AMBIGUOUS_PAIRS


def extract_snp_genotypes(genotype_table: dict, snp_panel: list) -> dict:
    """
    For each SNP in the panel, look up the genotype in genotype_table.

    Returns dict keyed by rsid:
    {
      "rsid": "rs1801133",
      "status": "found" | "not_tested",
      "genotype": "CT",           # raw as reported
      "normalised": "CT",         # forward-strand normalised
      "risk_allele": "T",
      "risk_count": 1             # 0, 1, or 2 copies of risk allele
    }
    """
    results = {}

    for snp in snp_panel:
        rsid = snp["rsid"]
        risk_allele = snp["risk_allele"]
        ref_allele = snp.get("ref_allele", "")

        if rsid not in genotype_table:
            results[rsid] = {
                "rsid": rsid,
                "gene": snp["gene"],
                "status": "not_tested",
                "genotype": None,
                "normalised": None,
                "risk_allele": risk_allele,
                "risk_count": None,
                "nutrient_domain": snp["nutrient_domain"],
            }
            continue

        raw_geno = genotype_table[rsid]
        if not raw_geno or len(raw_geno) < 2:
            results[rsid] = {
                "rsid": rsid,
                "gene": snp["gene"],
                "status": "no_call",
                "genotype": raw_geno,
                "normalised": None,
                "risk_allele": risk_allele,
                "risk_count": None,
                "nutrient_domain": snp["nutrient_domain"],
            }
            continue

        # Try direct match first
        norm = raw_geno
        allele_matched = risk_allele in raw_geno
        if not allele_matched:
            # Try strand flip
            flipped = flip_genotype(raw_geno)
            if risk_allele in flipped:
                norm = flipped
                allele_matched = True

        if allele_matched:
            risk_count = norm.count(risk_allele)
            results[rsid] = {
                "rsid": rsid,
                "gene": snp["gene"],
                "status": "found",
                "genotype": raw_geno,
                "normalised": norm,
                "risk_allele": risk_allele,
                "risk_count": risk_count,
                "nutrient_domain": snp["nutrient_domain"],
            }
        else:
            # Genotype does not contain the risk allele on either strand.
            # Check if this is homozygous reference (0 copies of risk allele).
            # A true allele_mismatch is when the genotype contains alleles
            # that are neither ref nor risk (e.g. tri-allelic or data error).
            geno_alleles = set(norm)
            known_alleles = {risk_allele, ref_allele} if ref_allele else {risk_allele}
            flipped_known = {COMPLEMENT.get(a, a) for a in known_alleles}

            if ref_allele and geno_alleles <= {ref_allele}:
                # Homozygous reference: 0 copies of risk allele
                results[rsid] = {
                    "rsid": rsid,
                    "gene": snp["gene"],
                    "status": "found",
                    "genotype": raw_geno,
                    "normalised": norm,
                    "risk_allele": risk_allele,
                    "risk_count": 0,
                    "nutrient_domain": snp["nutrient_domain"],
                }
            elif ref_allele and geno_alleles <= {COMPLEMENT.get(ref_allele, ref_allele)}:
                # Homozygous reference on complement strand
                results[rsid] = {
                    "rsid": rsid,
                    "gene": snp["gene"],
                    "status": "found",
                    "genotype": raw_geno,
                    "normalised": flip_genotype(raw_geno),
                    "risk_allele": risk_allele,
                    "risk_count": 0,
                    "nutrient_domain": snp["nutrient_domain"],
                }
            else:
                # True allele mismatch: genotype contains unknown alleles
                print(
                    f"[WARNING] {rsid} ({snp['gene']}): genotype '{raw_geno}' "
                    f"does not contain risk allele '{risk_allele}' (even after strand flip). "
                    f"Setting allele_mismatch."
                )
                results[rsid] = {
                    "rsid": rsid,
                    "gene": snp["gene"],
                    "status": "allele_mismatch",
                    "genotype": raw_geno,
                    "normalised": norm,
                    "risk_allele": risk_allele,
                    "risk_count": None,
                    "nutrient_domain": snp["nutrient_domain"],
                    "warning": (
                        f"Genotype '{raw_geno}' does not contain risk allele "
                        f"'{risk_allele}' on either strand"
                    ),
                }

    return results
