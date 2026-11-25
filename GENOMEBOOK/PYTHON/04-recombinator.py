"""
04-recombinator.py — Recombinator Engine for Genomebook

Purpose: Produce offspring from two parent genomes via meiotic recombination + mutation.
Input:  Two parent .genome.json files + trait_registry.json + disease_registry.json
Output: Offspring .genome.json with inferred trait scores and clinical history.

Biology model:
  - One allele inherited from each parent (random selection per locus)
  - Crossover: configurable rate (default 1 crossover per chromosome)
  - Mutation: configurable rate per locus per generation (~0.1%)
  - Sex determination: 50/50 coin flip
  - Clinical evaluation: check offspring genotype against disease registry
"""

import json
import random
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"
TRAIT_REGISTRY = DATA / "trait_registry.json"
DISEASE_REGISTRY = DATA / "disease_registry.json"

# Mutation configuration
MUTATION_CONFIG = {
    "rate_per_locus": 0.001,
    "disease_introduction_rate": 0.01,
    "protective_rate": 0.005,
    "neutral_rate": 0.985,
    "hotspot_categories": ["cognitive", "immune", "metabolic"],
    "hotspot_multiplier": 3.0,
}


def load_registries():
    with open(TRAIT_REGISTRY) as f:
        traits = json.load(f)
    with open(DISEASE_REGISTRY) as f:
        diseases = json.load(f)
    return traits, diseases


def inherit_allele(parent_locus):
    """Pick one allele from a diploid parent locus (Mendelian segregation)."""
    return random.choice(parent_locus["alleles"])


def mutate_allele(allele, locus_def, trait_category, config):
    """Possibly mutate an allele. Returns (allele, mutation_record_or_None)."""
    rate = config["rate_per_locus"]
    if trait_category in config["hotspot_categories"]:
        rate *= config["hotspot_multiplier"]

    if random.random() > rate:
        return allele, None

    # Mutation fires — flip allele
    ref = locus_def["ref"]
    alt = locus_def["alt"]
    new_allele = alt if allele == ref else ref

    # Classify mutation type
    roll = random.random()
    if roll < config["disease_introduction_rate"]:
        mut_type = "disease_risk"
    elif roll < config["disease_introduction_rate"] + config["protective_rate"]:
        mut_type = "protective"
    else:
        mut_type = "neutral"

    record = {
        "locus": locus_def["id"] if "id" in locus_def else "unknown",
        "from": allele,
        "to": new_allele,
        "type": mut_type,
        "category": trait_category,
    }
    return new_allele, record


def _short_name(genome):
    """Extract a short display name, stripping recursive 'Offspring of' prefixes."""
    name = genome.get("name", genome["id"])
    # For generation-0 founders, use their real name
    if genome.get("generation", 0) == 0:
        # Take first name + last name (e.g. "Albert Einstein" -> "Einstein")
        parts = name.split()
        return parts[-1] if parts else name
    # For offspring, use their ID which is always short
    return genome["id"]


def _offspring_name(parent_a, parent_b, generation, index):
    """Generate a concise offspring name."""
    a_short = _short_name(parent_a)
    b_short = _short_name(parent_b)
    return f"G{generation}-{index:03d} ({a_short} x {b_short})"


def recombine(parent_a, parent_b, generation, offspring_index, trait_reg, disease_reg):
    """Produce one offspring from two parents via recombination + mutation.

    Args:
        parent_a: Father genome (must be Male)
        parent_b: Mother genome (must be Female)
        generation: Offspring generation number
        offspring_index: Index within this generation (for ID)
        trait_reg: Trait registry
        disease_reg: Disease registry

    Returns:
        Offspring genome dict
    """
    assert parent_a["sex"] == "Male", f"Parent A must be Male, got {parent_a['sex']}"
    assert parent_b["sex"] == "Female", f"Parent B must be Female, got {parent_b['sex']}"

    # Sex determination
    sex = random.choice(["Male", "Female"])
    sex_chr = "XY" if sex == "Male" else "XX"

    # Derive offspring ID — short hash to prevent filename explosion
    import hashlib
    hash_input = f"{parent_a['id']}:{parent_b['id']}:{generation}:{offspring_index}"
    short_hash = hashlib.md5(hash_input.encode()).hexdigest()[:6]
    offspring_id = f"g{generation}-{offspring_index:03d}-{short_hash}"

    # Build locus-to-trait-category map
    locus_category = {}
    locus_defs = {}
    for tname, tdef in trait_reg["traits"].items():
        for locus in tdef["loci"]:
            lid = locus["id"]
            locus_category[lid] = tdef.get("category", "unknown")
            locus_defs[lid] = locus

    # Inherit and possibly mutate each locus
    offspring_loci = {}
    mutations = []
    shared_loci = set(parent_a["loci"].keys()) & set(parent_b["loci"].keys())

    for lid in shared_loci:
        a_locus = parent_a["loci"][lid]
        b_locus = parent_b["loci"][lid]

        # Mendelian: one allele from each parent
        from_father = inherit_allele(a_locus)
        from_mother = inherit_allele(b_locus)

        # Mutation pass
        category = locus_category.get(lid, "unknown")
        ldef = locus_defs.get(lid, a_locus)

        from_father, mut_f = mutate_allele(from_father, ldef, category, MUTATION_CONFIG)
        if mut_f:
            mut_f["parent"] = "father"
            mutations.append(mut_f)

        from_mother, mut_m = mutate_allele(from_mother, ldef, category, MUTATION_CONFIG)
        if mut_m:
            mut_m["parent"] = "mother"
            mutations.append(mut_m)

        offspring_loci[lid] = {
            "chromosome": a_locus["chromosome"],
            "position": a_locus["position"],
            "ref": a_locus["ref"],
            "alt": a_locus["alt"],
            "dominance": a_locus["dominance"],
            "effect_size": a_locus["effect_size"],
            "alleles": [from_father, from_mother],
        }

    # Infer trait scores from genotype
    trait_scores = infer_traits(offspring_loci, trait_reg)

    # Clinical evaluation
    clinical_history = evaluate_clinical(offspring_loci, trait_scores, disease_reg, generation)

    # Compute health score
    total_fitness_cost = sum(c.get("fitness_cost", 0) for c in clinical_history)
    health_score = max(0.0, min(1.0, 1.0 + total_fitness_cost))

    # Derive blended ancestry
    a_anc = parent_a.get("ancestry", "Unknown")
    b_anc = parent_b.get("ancestry", "Unknown")
    if a_anc == b_anc:
        ancestry = a_anc
    else:
        ancestry = f"{a_anc} / {b_anc}"

    # Derive domain (blend of parents)
    a_dom = parent_a.get("domain", "")
    b_dom = parent_b.get("domain", "")
    domain = f"{a_dom} + {b_dom}" if a_dom != b_dom else a_dom

    offspring = {
        "id": offspring_id,
        "name": _offspring_name(parent_a, parent_b, generation, offspring_index),
        "sex": sex,
        "sex_chromosomes": sex_chr,
        "ancestry": ancestry,
        "domain": domain,
        "era": f"Generation {generation}",
        "summary": None,  # Can be generated by LLM later
        "generation": generation,
        "parents": [parent_a["id"], parent_b["id"]],
        "loci": offspring_loci,
        "trait_scores": trait_scores,
        "mutations": mutations,
        "clinical_history": clinical_history,
        "health_score": round(health_score, 4),
    }

    return offspring


def infer_traits(loci, trait_reg):
    """Reverse-map genotype to trait scores."""
    scores = {}
    for tname, tdef in trait_reg["traits"].items():
        total_score = 0.0
        total_weight = 0.0
        for locus_def in tdef["loci"]:
            lid = locus_def["id"]
            if lid not in loci:
                continue
            alleles = loci[lid]["alleles"]
            ref = locus_def["ref"]
            alt = locus_def["alt"]
            dominance = locus_def["dominance"]
            effect = locus_def["effect"]

            # Count ALT alleles
            alt_count = alleles.count(alt)

            if dominance == "additive":
                locus_score = alt_count * 0.5  # 0, 0.5, or 1.0
            elif dominance == "dominant":
                locus_score = 1.0 if alt_count >= 1 else 0.0
            elif dominance == "recessive":
                locus_score = 1.0 if alt_count == 2 else 0.0
            else:
                locus_score = alt_count * 0.5

            total_score += locus_score * effect
            total_weight += effect

        if total_weight > 0:
            scores[tname] = round(total_score / total_weight, 4)
        else:
            scores[tname] = 0.5

    return scores


def evaluate_clinical(loci, trait_scores, disease_reg, generation):
    """Check offspring genotype against disease registry."""
    conditions = []

    for dname, ddef in disease_reg.get("diseases", {}).items():
        req = ddef.get("required_genotype", {})
        threshold_expr = ddef.get("threshold", None)

        # Check genotype requirements
        genotype_met = True
        for locus_id, req_geno in req.items():
            if locus_id not in loci:
                genotype_met = False
                break
            alleles = loci[locus_id]["alleles"]
            ref = loci[locus_id]["ref"]
            alt = loci[locus_id]["alt"]
            alt_count = alleles.count(alt)

            if req_geno == "alt/alt" and alt_count != 2:
                genotype_met = False
                break
            elif req_geno == "alt/?" and alt_count < 1:
                genotype_met = False
                break
            elif req_geno == "ref/ref" and alt_count != 0:
                genotype_met = False
                break

        if not genotype_met:
            continue

        # Check penetrance
        penetrance = ddef.get("penetrance", 1.0)
        if random.random() > penetrance:
            continue

        # Check onset probability
        onset_prob = ddef.get("onset_probability_per_gen", 1.0)
        if random.random() > onset_prob:
            continue

        conditions.append({
            "name": dname,
            "onset_generation": generation,
            "severity": ddef.get("severity", "unknown"),
            "inheritance": ddef.get("inheritance", "unknown"),
            "fitness_cost": ddef.get("fitness_cost", 0),
            "longevity_modifier": ddef.get("longevity_modifier", 0),
        })

    return conditions


def breed_pair(father, mother, generation, num_offspring, trait_reg, disease_reg):
    """Produce multiple offspring from a mating pair."""
    offspring = []
    for i in range(num_offspring):
        child = recombine(father, mother, generation, i + 1, trait_reg, disease_reg)
        offspring.append(child)
    return offspring


if __name__ == "__main__":
    # Demo: breed Einstein × Anning
    trait_reg, disease_reg = load_registries()
    genomes_dir = DATA / "GENOMES"

    father = json.load(open(genomes_dir / "einstein-g0.genome.json"))
    mother = json.load(open(genomes_dir / "anning-g0.genome.json"))

    children = breed_pair(father, mother, generation=1, num_offspring=3,
                          trait_reg=trait_reg, disease_reg=disease_reg)

    for c in children:
        print(f"\n{'='*60}")
        print(f"ID:     {c['id']}")
        print(f"Sex:    {c['sex']} ({c['sex_chromosomes']})")
        print(f"Health: {c['health_score']}")
        print(f"Mutations: {len(c['mutations'])}")
        if c['mutations']:
            for m in c['mutations']:
                print(f"  - {m['locus']}: {m['from']}→{m['to']} ({m['type']}, from {m['parent']})")
        print(f"Conditions: {len(c['clinical_history'])}")
        for cond in c['clinical_history']:
            print(f"  - {cond['name']} ({cond['severity']}, fitness: {cond['fitness_cost']})")
        print(f"Top traits:")
        top = sorted(c['trait_scores'].items(), key=lambda x: x[1], reverse=True)[:5]
        for t, s in top:
            print(f"  - {t}: {s}")
