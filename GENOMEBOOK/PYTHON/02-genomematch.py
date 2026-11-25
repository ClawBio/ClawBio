"""
02-genomematch.py — GenomeMatch Compatibility Engine for Genomebook

Purpose: Score genetic compatibility between all M×F pairings in a generation.
Input:  DATA/GENOMES/*.genome.json (for a given generation)
Output: Ranked compatibility list with reproductive risk assessment.

Compatibility scoring:
  - Allelic complementarity: heterozygous offspring are favoured (heterosis advantage)
  - Disease risk: flag pairs where both carry recessive disease alleles
  - Trait diversity: offspring genetic diversity score
  - Only M×F pairings are scored for reproduction.
"""

import json
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"
GENOMES_DIR = DATA / "GENOMES"
DISEASE_REGISTRY = DATA / "disease_registry.json"


def load_genomes(generation=0, genome_dir=None):
    """Load all genomes for a given generation.

    Uses filename prefix to narrow file scanning:
      - Gen 0: files NOT starting with 'g' (e.g. darwin-g0.genome.json)
      - Gen N: files starting with 'gN-' (e.g. g1-001-abc123.genome.json)
    Falls back to full scan if prefix matching finds nothing.
    """
    gdir = genome_dir or GENOMES_DIR
    genomes = {}

    # Fast path: match by filename prefix
    if generation == 0:
        candidates = [gf for gf in sorted(Path(gdir).glob("*.genome.json"))
                       if not gf.name.startswith("g") or gf.name.startswith("g0")]
    else:
        prefix = f"g{generation}-"
        candidates = sorted(Path(gdir).glob(f"{prefix}*.genome.json"))

    for gf in candidates:
        g = json.load(open(gf))
        if g["generation"] == generation:
            genomes[g["id"]] = g

    # Fallback: full scan if prefix matching found nothing
    if not genomes:
        for gf in sorted(Path(gdir).glob("*.genome.json")):
            g = json.load(open(gf))
            if g["generation"] == generation:
                genomes[g["id"]] = g

    return genomes


def load_disease_registry():
    with open(DISEASE_REGISTRY, "r") as f:
        return json.load(f)


def heterozygosity_score(parent_a, parent_b):
    """Score based on how many offspring loci would be heterozygous.

    Higher = more genetic diversity in offspring = better.
    """
    shared_loci = set(parent_a["loci"].keys()) & set(parent_b["loci"].keys())
    if not shared_loci:
        return 0.0

    het_count = 0
    for lid in shared_loci:
        a_alleles = parent_a["loci"][lid]["alleles"]
        b_alleles = parent_b["loci"][lid]["alleles"]
        # Simulate: pick one allele from each parent
        # If parents differ, offspring is more likely heterozygous
        a_unique = set(a_alleles)
        b_unique = set(b_alleles)
        if a_unique != b_unique:
            het_count += 1

    return het_count / len(shared_loci)


def disease_risk_score(parent_a, parent_b, disease_reg):
    """Assess risk of offspring inheriting disease conditions.

    Returns:
        risk_score: 0.0 (no risk) to 1.0 (high risk)
        flagged_diseases: list of disease names with risk details
    """
    diseases = disease_reg.get("diseases", {})
    flagged = []

    for dname, ddef in diseases.items():
        if ddef.get("inheritance") == "autosomal_recessive":
            req = ddef.get("required_genotype", {})
            for locus_id, req_geno in req.items():
                if req_geno != "alt/alt":
                    continue
                # Check if both parents carry at least one ALT
                a_locus = parent_a["loci"].get(locus_id)
                b_locus = parent_b["loci"].get(locus_id)
                if not a_locus or not b_locus:
                    continue

                a_has_alt = a_locus["alt"] in a_locus["alleles"]
                b_has_alt = b_locus["alt"] in b_locus["alleles"]

                if a_has_alt and b_has_alt:
                    # Both carriers — 25% chance of affected offspring
                    flagged.append({
                        "disease": dname,
                        "risk": "25% affected",
                        "mechanism": f"Both parents carry ALT at {locus_id}",
                        "severity": ddef.get("severity", "unknown"),
                    })

    risk_score = min(1.0, len(flagged) * 0.15)
    return risk_score, flagged


def trait_complementarity(parent_a, parent_b):
    """Score how well trait profiles complement each other.

    Penalise identical extremes, reward balanced combinations.
    """
    a_traits = parent_a.get("trait_scores", {})
    b_traits = parent_b.get("trait_scores", {})
    shared = set(a_traits.keys()) & set(b_traits.keys())
    if not shared:
        return 0.5

    complementarity = 0.0
    for t in shared:
        diff = abs(a_traits[t] - b_traits[t])
        avg = (a_traits[t] + b_traits[t]) / 2.0
        # Reward: moderate difference (complementary) + high average (both strong)
        complementarity += (diff * 0.3 + avg * 0.7)

    return complementarity / len(shared)


def compute_compatibility(parent_a, parent_b, disease_reg):
    """Full compatibility score between two agents."""
    het = heterozygosity_score(parent_a, parent_b)
    risk, flagged = disease_risk_score(parent_a, parent_b, disease_reg)
    comp = trait_complementarity(parent_a, parent_b)

    # Weighted final score: heterozygosity(40%) + complementarity(40%) - risk(20%)
    final = (het * 0.40) + (comp * 0.40) - (risk * 0.20)
    final = max(0.0, min(1.0, final))

    return {
        "score": round(final, 4),
        "heterozygosity": round(het, 4),
        "complementarity": round(comp, 4),
        "disease_risk": round(risk, 4),
        "flagged_diseases": flagged,
    }


def match_generation(genomes, disease_reg, top_n=None):
    """Score all valid M×F pairings and return ranked list."""
    males = {gid: g for gid, g in genomes.items() if g["sex"] == "Male"}
    females = {gid: g for gid, g in genomes.items() if g["sex"] == "Female"}

    pairings = []
    for mid, mg in males.items():
        for fid, fg in females.items():
            compat = compute_compatibility(mg, fg, disease_reg)
            pairings.append({
                "male": mid,
                "female": fid,
                "male_name": mg.get("name", mid),
                "female_name": fg.get("name", fid),
                **compat,
            })

    pairings.sort(key=lambda x: x["score"], reverse=True)

    if top_n:
        pairings = pairings[:top_n]

    return pairings


def select_mating_pairs(pairings, max_pairs=10):
    """Select non-overlapping mating pairs from ranked compatibility list.

    Each agent can only mate once per generation. Greedy selection from top.
    """
    used = set()
    selected = []

    for p in pairings:
        if p["male"] in used or p["female"] in used:
            continue
        selected.append(p)
        used.add(p["male"])
        used.add(p["female"])
        if len(selected) >= max_pairs:
            break

    return selected


def main():
    disease_reg = load_disease_registry()
    genomes = load_genomes(generation=0)

    if not genomes:
        print("ERROR: No generation-0 genomes found.")
        return

    males = [g for g in genomes.values() if g["sex"] == "Male"]
    females = [g for g in genomes.values() if g["sex"] == "Female"]
    print(f"Loaded {len(genomes)} genomes ({len(males)}M / {len(females)}F)")

    pairings = match_generation(genomes, disease_reg)
    print(f"\nAll {len(pairings)} M×F pairings scored.\n")

    print(f"{'Rank':>4}  {'Male':>20} × {'Female':<20}  {'Score':>6}  {'Het':>5}  {'Comp':>5}  {'Risk':>5}  Flags")
    print("-" * 100)
    for i, p in enumerate(pairings[:20], 1):
        flags = ", ".join(d["disease"] for d in p["flagged_diseases"]) or "—"
        print(f"{i:4d}  {p['male']:>20} × {p['female']:<20}  {p['score']:6.4f}  {p['heterozygosity']:5.3f}  {p['complementarity']:5.3f}  {p['disease_risk']:5.3f}  {flags}")

    selected = select_mating_pairs(pairings, max_pairs=10)
    print(f"\n{'='*60}")
    print(f"SELECTED MATING PAIRS (generation 0 → 1):")
    print(f"{'='*60}")
    for p in selected:
        print(f"  {p['male_name']} × {p['female_name']}  (compat: {p['score']:.4f})")


if __name__ == "__main__":
    main()
