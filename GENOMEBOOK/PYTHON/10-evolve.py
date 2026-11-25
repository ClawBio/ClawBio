"""
10-evolve.py -- Genomebook Evolution Orchestrator

Runs multi-generation agent evolution: breed pairs, generate offspring
SOUL.md + DNA.md, run Moltbook interactions, log population stats.

Usage:
    python 10-evolve.py --dry-run --generations 100        # Genetics only, zero cost
    python 10-evolve.py --generations 100 --rounds 3       # Full run with Moltbook
    python 10-evolve.py --generations 10 --budget 200      # Budget-capped
    python 10-evolve.py --generations 5 --rounds 1 --model claude-haiku-4-5-20251001
"""

import argparse
import json
import os
import random
import time
from datetime import datetime, timezone
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec

# Paths
BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"
GENOMES_DIR = DATA / "GENOMES"
SOULS_DIR = DATA / "SOULS"
DNA_DIR = DATA / "DNA"
EVOLUTION_LOG = DATA / "evolution_log.jsonl"
PYTHON_DIR = Path(__file__).resolve().parent

DNA_DIR.mkdir(parents=True, exist_ok=True)
SOULS_DIR.mkdir(parents=True, exist_ok=True)

# Environment
ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY", "")
MOLTBOOK_URL = os.environ.get("MOLTBOOK_URL", "http://127.0.0.1:8800")
DEFAULT_MODEL = os.environ.get("MOLTBOOK_MODEL", "claude-sonnet-4-5-20250929")


# ── Module imports (sibling scripts) ─────────────────────────────────────

def _load_module(name, filename):
    spec = spec_from_file_location(name, PYTHON_DIR / filename)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_genomematch = None
_recombinator = None
_dna_compiler = None


def get_genomematch():
    global _genomematch
    if _genomematch is None:
        _genomematch = _load_module("genomematch", "02-genomematch.py")
    return _genomematch


def get_recombinator():
    global _recombinator
    if _recombinator is None:
        _recombinator = _load_module("recombinator", "04-recombinator.py")
    return _recombinator


def get_dna_compiler():
    global _dna_compiler
    if _dna_compiler is None:
        _dna_compiler = _load_module("dna_compiler", "06-dna_compiler.py")
    return _dna_compiler


# ── Offspring SOUL.md generation ─────────────────────────────────────────

def generate_offspring_soul(genome):
    """Generate a SOUL.md for an offspring from its genome data.

    Pure Python, no LLM. Reverse of parse_soul().
    """
    name = genome.get("name", genome["id"])
    sex = genome["sex"]
    era = genome.get("era", f"Generation {genome['generation']}")
    domain = genome.get("domain", "Unknown")
    ancestry = genome.get("ancestry", "Unknown")
    traits = genome.get("trait_scores", {})
    health = genome.get("health_score", 1.0)
    parents = genome.get("parents", [None, None])

    # Build top traits description
    top = sorted(traits.items(), key=lambda x: x[1], reverse=True)[:5]
    top_str = ", ".join(f"{t.replace('_', ' ')} ({s:.2f})" for t, s in top)

    # Build parent names from IDs
    parent_desc = ""
    if parents and parents[0] and parents[1]:
        parent_desc = f"Offspring of {parents[0]} and {parents[1]}."

    lines = []
    lines.append(f"# {name}")
    lines.append(f"## Identity")
    lines.append(f"- **Name:** {name}")
    lines.append(f"- **Sex:** {sex}")
    lines.append(f"- **Era:** {era}")
    lines.append(f"- **Domain:** {domain}")
    lines.append(f"- **Ancestry:** {ancestry}")
    lines.append(f"## Trait Scores")
    for trait_name in sorted(traits.keys()):
        lines.append(f"{trait_name}: {traits[trait_name]:.2f}")
    lines.append(f"## Summary")
    summary = f"{parent_desc} Strongest traits: {top_str}. Health score: {health:.2f}."
    lines.append(summary)

    return "\n".join(lines)


# ── Moltbook interaction ────────────────────────────────────────────────

def run_moltbook_round(agent_ids, model, round_num):
    """Run one Moltbook interaction round for a set of agents.

    Returns number of API calls made.
    """
    # Import agent runner
    agent_mod = _load_module("moltbook_agent", "09-moltbook_agent.py")

    # Override model and URL
    agent_mod.LLM_MODEL = model
    agent_mod.MOLTBOOK_URL = MOLTBOOK_URL

    api_calls = 0
    agents = []
    for aid in agent_ids:
        try:
            agent = agent_mod.GenomebookAgent(aid)
            agents.append(agent)
        except Exception as e:
            print(f"    [skip] {aid}: {e}")

    random.shuffle(agents)
    for agent in agents:
        try:
            agent.decide_and_act()
            api_calls += 1
        except Exception as e:
            print(f"    [{agent.soul_name}] error: {e}")
        time.sleep(1)  # Rate limit spacing

    return api_calls


# ── Generation stats ────────────────────────────────────────────────────

def compute_gen_stats(generation, genomes, pairings, offspring, selected_pairs):
    """Compute population statistics for one generation."""
    all_genomes = list(genomes.values())
    if offspring:
        all_genomes.extend(offspring)

    # Trait means
    trait_means = {}
    trait_stdevs = {}
    if all_genomes:
        all_traits = set()
        for g in all_genomes:
            all_traits.update(g.get("trait_scores", {}).keys())

        for trait in sorted(all_traits):
            values = [g["trait_scores"].get(trait, 0.5) for g in all_genomes if "trait_scores" in g]
            if values:
                mean = sum(values) / len(values)
                variance = sum((v - mean) ** 2 for v in values) / len(values)
                trait_means[trait] = round(mean, 4)
                trait_stdevs[trait] = round(variance ** 0.5, 4)

    # Health
    health_scores = [g.get("health_score", 1.0) for g in all_genomes]
    avg_health = round(sum(health_scores) / len(health_scores), 4) if health_scores else 1.0
    min_health = round(min(health_scores), 4) if health_scores else 1.0
    max_health = round(max(health_scores), 4) if health_scores else 1.0

    # Mutations (offspring only)
    total_mutations = 0
    mutation_types = {"disease_risk": 0, "protective": 0, "neutral": 0}
    for child in (offspring or []):
        muts = child.get("mutations", [])
        total_mutations += len(muts)
        for m in muts:
            mt = m.get("type", "neutral")
            mutation_types[mt] = mutation_types.get(mt, 0) + 1

    # Disease counts
    disease_counts = {}
    for g in all_genomes:
        for cond in g.get("clinical_history", []):
            dname = cond["name"]
            disease_counts[dname] = disease_counts.get(dname, 0) + 1

    # Heterozygosity
    het_values = []
    for g in all_genomes:
        loci = g.get("loci", {})
        if loci:
            het = sum(1 for l in loci.values() if len(set(l["alleles"])) > 1) / len(loci)
            het_values.append(het)
    avg_het = round(sum(het_values) / len(het_values), 4) if het_values else 0.0

    # Sex ratio
    males = sum(1 for g in all_genomes if g["sex"] == "Male")
    females = sum(1 for g in all_genomes if g["sex"] == "Female")

    return {
        "generation": generation,
        "timestamp": datetime.now(timezone.utc).isoformat(),
        "population_size": len(all_genomes),
        "parents_count": len(genomes),
        "offspring_count": len(offspring) if offspring else 0,
        "pairs_selected": len(selected_pairs) if selected_pairs else 0,
        "males": males,
        "females": females,
        "avg_health": avg_health,
        "min_health": min_health,
        "max_health": max_health,
        "avg_heterozygosity": avg_het,
        "total_mutations": total_mutations,
        "mutation_types": mutation_types,
        "trait_means": trait_means,
        "trait_stdevs": trait_stdevs,
        "disease_prevalence": disease_counts,
    }


# ── Main evolution loop ─────────────────────────────────────────────────

def evolve(
    start_gen=0,
    num_gens=100,
    offspring_per_pair=3,
    rounds_per_gen=3,
    max_budget=6000,
    model=DEFAULT_MODEL,
    dry_run=False,
    seed=None,
):
    """Run the full evolution simulation."""

    if seed is not None:
        random.seed(seed)

    gm = get_genomematch()
    rec = get_recombinator()
    dna = get_dna_compiler()
    trait_reg, disease_reg = rec.load_registries()

    total_api_calls = 0
    budget_exhausted = False

    print("=" * 70)
    print("GENOMEBOOK EVOLUTION ORCHESTRATOR")
    print("=" * 70)
    print(f"Generations:      {start_gen} -> {start_gen + num_gens - 1}")
    print(f"Offspring/pair:   {offspring_per_pair}")
    print(f"Rounds/gen:       {rounds_per_gen}")
    print(f"API budget:       {max_budget} calls")
    print(f"Model:            {model}")
    print(f"Dry run:          {dry_run}")
    print(f"Evolution log:    {EVOLUTION_LOG}")
    print()

    for gen in range(start_gen, start_gen + num_gens):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"\n{'='*70}")
        print(f"GENERATION {gen} ({ts})")
        print(f"{'='*70}")

        # 1. Load current generation genomes
        genomes = gm.load_genomes(generation=gen)
        if not genomes:
            if gen == 0:
                print("ERROR: No generation-0 genomes found. Run Soul2DNA first.")
                return
            else:
                print(f"  No genomes for generation {gen}. Evolution complete.")
                break

        males = {gid: g for gid, g in genomes.items() if g["sex"] == "Male"}
        females = {gid: g for gid, g in genomes.items() if g["sex"] == "Female"}
        print(f"  Population: {len(genomes)} ({len(males)}M / {len(females)}F)")

        if not males or not females:
            print(f"  Cannot breed: need both males and females. Stopping.")
            break

        # 2. Score compatibility
        pairings = gm.match_generation(genomes, disease_reg)
        print(f"  Pairings scored: {len(pairings)}")

        # 3. Select mating pairs
        max_pairs = min(len(males), len(females))
        selected = gm.select_mating_pairs(pairings, max_pairs=max_pairs)
        print(f"  Pairs selected: {len(selected)}")

        for p in selected:
            print(f"    {p['male_name']:<30s} x {p['female_name']:<30s}  (compat: {p['score']:.4f})")

        # 4. Breed offspring
        all_offspring = []
        next_gen = gen + 1

        for pair in selected:
            father = genomes[pair["male"]]
            mother = genomes[pair["female"]]

            children = rec.breed_pair(
                father, mother,
                generation=next_gen,
                num_offspring=offspring_per_pair,
                trait_reg=trait_reg,
                disease_reg=disease_reg,
            )

            for child in children:
                all_offspring.append(child)

                # 5. Write genome
                genome_path = GENOMES_DIR / f"{child['id']}.genome.json"
                with open(genome_path, "w") as f:
                    json.dump(child, f, indent=2)

                # 6. Generate and write SOUL.md
                soul_text = generate_offspring_soul(child)
                soul_path = SOULS_DIR / f"{child['id']}.soul.md"
                soul_path.write_text(soul_text)

                # 7. Generate and write DNA.md
                dna_text = dna.compile_dna_md(child, trait_reg, disease_reg)
                dna_path = DNA_DIR / f"{child['id']}.dna.md"
                dna_path.write_text(dna_text)

        print(f"  Offspring born: {len(all_offspring)} (generation {next_gen})")

        # Print offspring summary
        for child in all_offspring:
            mut_count = len(child.get("mutations", []))
            cond_count = len(child.get("clinical_history", []))
            print(f"    {child['id']:20s} | {child['sex']:6s} | health: {child.get('health_score', 1.0):.2f} | {mut_count} mutations | {cond_count} conditions")

        # 8. Moltbook interaction
        if not dry_run and not budget_exhausted:
            # Get agent IDs for current generation (parents who just bred)
            agent_ids = list(genomes.keys())

            # Add some offspring if they have SOUL.md + DNA.md (they all do now)
            # But limit to keep costs reasonable
            offspring_agents = [c["id"] for c in all_offspring[:5]]  # Top 5 offspring
            agent_ids.extend(offspring_agents)

            remaining_budget = max_budget - total_api_calls
            if remaining_budget <= 0:
                print(f"  Budget exhausted ({total_api_calls}/{max_budget}). Skipping Moltbook.")
                budget_exhausted = True
            else:
                for r in range(rounds_per_gen):
                    if total_api_calls >= max_budget:
                        print(f"  Budget hit during round {r+1}. Stopping interactions.")
                        budget_exhausted = True
                        break

                    print(f"\n  Moltbook round {r+1}/{rounds_per_gen} (budget: {total_api_calls}/{max_budget})")
                    calls = run_moltbook_round(agent_ids, model, r)
                    total_api_calls += calls
                    print(f"  API calls this round: {calls} (total: {total_api_calls})")

        # 9. Log generation stats
        stats = compute_gen_stats(gen, genomes, pairings, all_offspring, selected)
        stats["api_calls_total"] = total_api_calls
        stats["api_calls_budget"] = max_budget

        with open(EVOLUTION_LOG, "a") as f:
            f.write(json.dumps(stats) + "\n")

        # 10. Summary
        print(f"\n  Gen {gen} summary: {stats['population_size']} agents, "
              f"avg health {stats['avg_health']:.3f}, "
              f"het {stats['avg_heterozygosity']:.3f}, "
              f"{stats['total_mutations']} new mutations")

    # Final summary
    print(f"\n{'='*70}")
    print(f"EVOLUTION COMPLETE")
    print(f"{'='*70}")
    print(f"Generations:   {start_gen} to {start_gen + num_gens - 1}")
    print(f"API calls:     {total_api_calls} / {max_budget}")
    print(f"Evolution log: {EVOLUTION_LOG}")
    print(f"Genomes:       {GENOMES_DIR}/")
    print(f"Souls:         {SOULS_DIR}/")
    print(f"DNA profiles:  {DNA_DIR}/")


def _set_moltbook_url(url):
    global MOLTBOOK_URL
    MOLTBOOK_URL = url


def main():
    parser = argparse.ArgumentParser(description="Genomebook Evolution Orchestrator")
    parser.add_argument("--generations", type=int, default=100, help="Number of generations (default 100)")
    parser.add_argument("--start-gen", type=int, default=0, help="Starting generation (default 0)")
    parser.add_argument("--offspring", type=int, default=3, help="Offspring per pair (default 3)")
    parser.add_argument("--rounds", type=int, default=3, help="Moltbook rounds per generation (default 3)")
    parser.add_argument("--budget", type=int, default=6000, help="Max API calls (default 6000)")
    parser.add_argument("--model", type=str, default=DEFAULT_MODEL, help="LLM model")
    parser.add_argument("--dry-run", action="store_true", help="Genetics only, no Moltbook (zero cost)")
    parser.add_argument("--seed", type=int, default=None, help="Random seed for reproducibility")
    parser.add_argument("--url", type=str, default=MOLTBOOK_URL, help="Moltbook server URL")
    args = parser.parse_args()

    _set_moltbook_url(args.url)

    evolve(
        start_gen=args.start_gen,
        num_gens=args.generations,
        offspring_per_pair=args.offspring,
        rounds_per_gen=args.rounds,
        max_budget=args.budget,
        model=args.model,
        dry_run=args.dry_run,
        seed=args.seed,
    )


if __name__ == "__main__":
    main()
