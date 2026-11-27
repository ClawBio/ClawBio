"""
20-evolve_v2.py -- Genomebook v2 Evolution Orchestrator

Integrates skill execution + fitness scoring + Mendelian breeding.
Each generation: assign tasks, execute skill-gated research rounds,
score fitness, use fitness-weighted mating selection, breed offspring.

Usage:
    python 20-evolve_v2.py --dry-run --generations 10        # Genetics + skill gates, no API calls
    python 20-evolve_v2.py --generations 10 --task '{"gene":"BRCA1"}'  # With real skills
    python 20-evolve_v2.py --dry-run --generations 50 --seed 42       # Reproducible run
    python 20-evolve_v2.py --random-skills --dry-run --generations 10  # Control: random assignment
"""

import argparse
import json
import random
from datetime import datetime, timezone
from pathlib import Path
from importlib.util import spec_from_file_location, module_from_spec

# Paths
BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"
GENOMES_DIR = DATA / "GENOMES"
SOULS_DIR = DATA / "SOULS"
DNA_DIR = DATA / "DNA"
PYTHON_DIR = Path(__file__).resolve().parent
OUTPUT_DIR = BASE / "OUTPUT"
EVOLUTION_LOG_V2 = DATA / "evolution_log_v2.jsonl"

for d in [GENOMES_DIR, SOULS_DIR, DNA_DIR, OUTPUT_DIR]:
    d.mkdir(parents=True, exist_ok=True)


# -- Module imports (sibling scripts) ----------------------------------------

def _load_module(name, filename):
    spec = spec_from_file_location(name, PYTHON_DIR / filename)
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


_genomematch = None
_recombinator = None
_dna_compiler = None
_skill_executor = None
_fitness_evaluator = None


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


def get_skill_executor():
    global _skill_executor
    if _skill_executor is None:
        _skill_executor = _load_module("skill_executor", "15-skill_executor.py")
    return _skill_executor


def get_fitness_evaluator():
    global _fitness_evaluator
    if _fitness_evaluator is None:
        _fitness_evaluator = _load_module("fitness_evaluator", "16-fitness_evaluator.py")
    return _fitness_evaluator


# -- Offspring SOUL.md generation (from 10-evolve.py) ------------------------

def generate_offspring_soul(genome):
    """Generate a SOUL.md for an offspring from its genome data. Pure Python."""
    name = genome.get("name", genome["id"])
    sex = genome["sex"]
    era = genome.get("era", f"Generation {genome['generation']}")
    domain = genome.get("domain", "Unknown")
    ancestry = genome.get("ancestry", "Unknown")
    traits = genome.get("trait_scores", {})
    health = genome.get("health_score", 1.0)
    parents = genome.get("parents", [None, None])

    top = sorted(traits.items(), key=lambda x: x[1], reverse=True)[:5]
    top_str = ", ".join(f"{t.replace('_', ' ')} ({s:.2f})" for t, s in top)

    parent_desc = ""
    if parents and parents[0] and parents[1]:
        parent_desc = f"Offspring of {parents[0]} and {parents[1]}."

    lines = [
        f"# {name}",
        f"## Identity",
        f"- **Name:** {name}",
        f"- **Sex:** {sex}",
        f"- **Era:** {era}",
        f"- **Domain:** {domain}",
        f"- **Ancestry:** {ancestry}",
        f"## Trait Scores",
    ]
    for trait_name in sorted(traits.keys()):
        lines.append(f"{trait_name}: {traits[trait_name]:.2f}")
    lines.append(f"## Summary")
    lines.append(f"{parent_desc} Strongest traits: {top_str}. Health score: {health:.2f}.")

    return "\n".join(lines)


# -- Fitness-weighted mating -------------------------------------------------

def select_mating_pairs_v2(pairings, fitness_map, max_pairs=10):
    """Select mating pairs weighted by both compatibility and fitness.

    fitness_map: {agent_id: composite_fitness_score}
    """
    # Weight each pairing by compatibility * average fitness of the pair
    for p in pairings:
        male_fit = fitness_map.get(p["male"], 0.5)
        female_fit = fitness_map.get(p["female"], 0.5)
        avg_fitness = (male_fit + female_fit) / 2
        p["fitness_weighted_score"] = p["score"] * (0.5 + 0.5 * avg_fitness)

    pairings.sort(key=lambda x: x["fitness_weighted_score"], reverse=True)

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


# -- Random skill assignment (control condition) -----------------------------

def randomize_skill_gates(agents, skill_registry, seed=None):
    """Randomly assign skill availability, ignoring trait gates.
    Used for Control Condition C (random assignment baseline).
    """
    rng = random.Random(seed)
    skill_names = sorted(skill_registry["skills"].keys())

    for agent in agents:
        # Give each agent a random subset of skills (same distribution as Mendelian)
        n_skills = rng.randint(1, len(skill_names))
        assigned = rng.sample(skill_names, n_skills)
        agent["_random_skills"] = set(assigned)


# -- Main evolution loop -----------------------------------------------------

def evolve_v2(
    start_gen=0,
    num_gens=10,
    offspring_per_pair=2,
    task=None,
    dry_run=False,
    random_skills=False,
    seed=None,
    population_cap=30,
):
    """Run the v2 evolution simulation with skill-gated research rounds."""

    if seed is not None:
        random.seed(seed)

    gm = get_genomematch()
    rec = get_recombinator()
    dna = get_dna_compiler()
    se = get_skill_executor()
    fe = get_fitness_evaluator()

    skill_registry = se.load_skill_registry()
    trait_reg, disease_reg = rec.load_registries()

    # Default task if none provided
    if task is None:
        task = {"gene": "BRCA1"}

    # Use separate log/output files for random-skills control
    suffix = "_random" if random_skills else ""
    evolution_log = DATA / f"evolution_log_v2{suffix}.jsonl"
    summary_csv = OUTPUT_DIR / f"evolution_v2_summary{suffix}.csv"

    print("=" * 72)
    print("GENOMEBOOK v2 EVOLUTION ORCHESTRATOR")
    print("=" * 72)
    print(f"Generations:      {start_gen} -> {start_gen + num_gens - 1}")
    print(f"Offspring/pair:   {offspring_per_pair}")
    print(f"Population cap:   {population_cap}")
    print(f"Task:             {json.dumps(task)}")
    print(f"Dry run:          {dry_run}")
    print(f"Random skills:    {random_skills}")
    print(f"Seed:             {seed}")
    print(f"Log:              {evolution_log}")
    print()

    all_gen_stats = []
    retirement_age = 3  # Agents removed after this many generations

    # Persistent population carried across generations
    current_population = {}

    for gen in range(start_gen, start_gen + num_gens):
        ts = datetime.now().strftime("%H:%M:%S")
        print(f"\n{'='*72}")
        print(f"GENERATION {gen} ({ts})")
        print(f"{'='*72}")

        # 1. Load or maintain population
        if gen == start_gen:
            genomes = gm.load_genomes(generation=gen)
            if not genomes:
                print("ERROR: No generation-0 genomes. Run Soul2DNA first.")
                return
            current_population = dict(genomes)
        else:
            # Retire old agents (born more than retirement_age generations ago)
            to_remove = [gid for gid, g in current_population.items()
                         if gen - g.get("generation", 0) > retirement_age]
            for gid in to_remove:
                del current_population[gid]

        genomes = dict(current_population)
        agents = list(genomes.values())

        # Apply population cap: keep agents sorted by fitness (if available),
        # otherwise by generation (newer first)
        if len(agents) > population_cap:
            agents.sort(key=lambda a: (-a.get("_fitness", 0.5),
                                        -a.get("generation", 0)))
            agents = agents[:population_cap]
            genomes = {a["id"]: a for a in agents}
            current_population = dict(genomes)

        males = {gid: g for gid, g in genomes.items() if g["sex"] == "Male"}
        females = {gid: g for gid, g in genomes.items() if g["sex"] == "Female"}
        print(f"  Population: {len(genomes)} ({len(males)}M / {len(females)}F)")

        if not males or not females:
            print(f"  Cannot breed: need both sexes. Stopping.")
            break

        # 2. Execute research rounds (skill-gated)
        #    For random-skills control: override trait scores with random values
        #    so skill gates activate randomly rather than by inheritance.
        if random_skills:
            skill_names = sorted(skill_registry["skills"].keys())
            for agent in agents:
                n_skills = random.randint(1, len(skill_names))
                assigned = set(random.sample(skill_names, n_skills))
                # Temporarily set trait scores to 1.0 for assigned skills' gates,
                # and 0.0 for all others, so only assigned skills pass gates
                fake_traits = {t: 0.0 for t in agent.get("trait_scores", {})}
                for sname in assigned:
                    for trait in skill_registry["skills"][sname]["gate_traits"]:
                        fake_traits[trait] = 1.0
                agent["_original_traits"] = agent.get("trait_scores", {})
                agent["trait_scores"] = fake_traits

        print(f"\n  --- Research Rounds ---")
        round_results = []
        for agent in agents:
            rr = se.execute_round(agent, task, skill_registry, dry_run=dry_run)
            round_results.append(rr)
            n_used = len(rr["skills_used"])
            n_avail = len(rr["skills_available"])
            print(f"    {agent['id']:20s} | {n_used}/{n_avail} skills used | "
                  f"entropy {rr['skill_entropy']:.3f}")

        # 3. Compute fitness
        # Restore original traits if random-skills mode was used
        if random_skills:
            for agent in agents:
                if "_original_traits" in agent:
                    agent["trait_scores"] = agent.pop("_original_traits")

        print(f"\n  --- Fitness Scoring ---")
        fitness_scores = fe.compute_generation_fitness(round_results)
        fitness_map = {f["agent_id"]: f["composite_fitness"] for f in fitness_scores}

        # Store fitness on agents for population cap sorting
        for agent in agents:
            agent["_fitness"] = fitness_map.get(agent["id"], 0.5)

        # Sort and display
        fitness_scores.sort(key=lambda x: x["composite_fitness"], reverse=True)
        for f in fitness_scores[:5]:
            print(f"    {f['agent_id']:20s} | fitness {f['composite_fitness']:.3f} "
                  f"(F1={f['f1_utilization']:.2f} F2={f['f2_synthesis_quality']:.2f} "
                  f"F3={f['f3_novelty']:.2f} F4={f['f4_skill_chaining']:.2f})")
        if len(fitness_scores) > 5:
            print(f"    ... ({len(fitness_scores) - 5} more)")

        gen_stats = fe.generation_stats(fitness_scores)

        # 4. Score compatibility with fitness weighting
        pairings = gm.match_generation(genomes, disease_reg)
        print(f"\n  --- Mating Selection ---")
        print(f"  Pairings scored: {len(pairings)}")

        max_pairs = min(len(males), len(females))
        selected = select_mating_pairs_v2(pairings, fitness_map, max_pairs=max_pairs)
        print(f"  Pairs selected: {len(selected)}")

        for p in selected:
            mf = fitness_map.get(p["male"], 0.5)
            ff = fitness_map.get(p["female"], 0.5)
            print(f"    {p['male_name']:<25s} (fit={mf:.2f}) x "
                  f"{p['female_name']:<25s} (fit={ff:.2f})")

        # 5. Breed offspring
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

                # Write genome
                genome_path = GENOMES_DIR / f"{child['id']}.genome.json"
                with open(genome_path, "w") as f:
                    json.dump(child, f, indent=2)

                # Write SOUL.md
                soul_text = generate_offspring_soul(child)
                soul_path = SOULS_DIR / f"{child['id']}.soul.md"
                soul_path.write_text(soul_text)

                # Write DNA.md
                dna_text = dna.compile_dna_md(child, trait_reg, disease_reg)
                dna_path = DNA_DIR / f"{child['id']}.dna.md"
                dna_path.write_text(dna_text)

        # Add offspring to persistent population
        for child in all_offspring:
            current_population[child["id"]] = child

        print(f"\n  --- Offspring (gen {next_gen}) ---")
        print(f"  Born: {len(all_offspring)}")

        # Show offspring skill portfolios
        for child in all_offspring:
            child_traits = child.get("trait_scores", {})
            portfolio = se.skill_portfolio_vector(child_traits, skill_registry)
            n_skills = sum(portfolio)
            mut_count = len(child.get("mutations", []))
            print(f"    {child['id']:20s} | {child['sex']:6s} | "
                  f"{n_skills}/7 skills | {mut_count} mutations | "
                  f"health {child.get('health_score', 1.0):.2f}")

        # 6. Compute population skill diversity metrics
        all_portfolios = []
        for a in agents:
            all_portfolios.append(
                se.skill_portfolio_vector(a.get("trait_scores", {}), skill_registry))
        for child in all_offspring:
            all_portfolios.append(
                se.skill_portfolio_vector(child.get("trait_scores", {}), skill_registry))

        # Population skill coverage
        skill_names = sorted(skill_registry["skills"].keys())
        coverage = [0] * len(skill_names)
        for pv in all_portfolios:
            for i, v in enumerate(pv):
                coverage[i] += v

        pop_size = len(all_portfolios)
        mean_skills = sum(sum(pv) for pv in all_portfolios) / pop_size
        mean_entropy = sum(se.skill_entropy(pv) for pv in all_portfolios) / pop_size

        # Unique portfolio count (diversity)
        unique_portfolios = len(set(tuple(pv) for pv in all_portfolios))

        # Skill chaining potential: agents with 3+ skills
        chaining_capable = sum(1 for pv in all_portfolios if sum(pv) >= 3)

        print(f"\n  --- Population Skill Metrics ---")
        print(f"  Mean skills/agent:   {mean_skills:.2f}")
        print(f"  Mean entropy:        {mean_entropy:.4f}")
        print(f"  Unique portfolios:   {unique_portfolios}/{pop_size}")
        print(f"  Chaining capable:    {chaining_capable}/{pop_size} (3+ skills)")
        print(f"  Skill coverage:")
        for i, sname in enumerate(skill_names):
            pct = coverage[i] / pop_size * 100
            print(f"    {sname:25s}: {coverage[i]:3d}/{pop_size} ({pct:.0f}%)")

        # 7. Log generation
        log_entry = {
            "generation": gen,
            "timestamp": datetime.now(timezone.utc).isoformat(),
            "population_size": pop_size,
            "parents_count": len(genomes),
            "offspring_count": len(all_offspring),
            "pairs_selected": len(selected),
            "males": len(males),
            "females": len(females),
            "task": task,
            "dry_run": dry_run,
            "random_skills": random_skills,
            **gen_stats,
            "mean_skills_per_agent": round(mean_skills, 4),
            "mean_skill_entropy": round(mean_entropy, 4),
            "unique_portfolios": unique_portfolios,
            "chaining_capable": chaining_capable,
            "skill_coverage": {skill_names[i]: coverage[i] for i in range(len(skill_names))},
        }
        all_gen_stats.append(log_entry)

        with open(evolution_log, "a") as f:
            f.write(json.dumps(log_entry) + "\n")

        print(f"\n  Gen {gen} summary: fitness={gen_stats.get('mean_fitness', 0):.3f} "
              f"skills={mean_skills:.1f} entropy={mean_entropy:.3f} "
              f"unique={unique_portfolios}")

    # Final summary
    print(f"\n{'='*72}")
    print(f"EVOLUTION v2 COMPLETE")
    print(f"{'='*72}")

    if all_gen_stats:
        print(f"\nFitness trajectory:")
        for gs in all_gen_stats:
            g = gs["generation"]
            mf = gs.get("mean_fitness", 0)
            ms = gs.get("mean_skills_per_agent", 0)
            me = gs.get("mean_skill_entropy", 0)
            cc = gs.get("chaining_capable", 0)
            ps = gs.get("population_size", 0)
            print(f"  Gen {g:3d}: fitness={mf:.3f}  skills={ms:.1f}  "
                  f"entropy={me:.3f}  chaining={cc}/{ps}")

    # Write summary CSV
    summary_path = summary_csv
    with open(summary_path, "w") as f:
        f.write("generation,population,mean_fitness,mean_f1,mean_f2,mean_f3,mean_f4,"
                "mean_skills,mean_entropy,unique_portfolios,chaining_capable\n")
        for gs in all_gen_stats:
            f.write(f"{gs['generation']},{gs['population_size']},"
                    f"{gs.get('mean_fitness',0)},{gs.get('mean_f1',0)},"
                    f"{gs.get('mean_f2',0)},{gs.get('mean_f3',0)},"
                    f"{gs.get('mean_f4',0)},"
                    f"{gs.get('mean_skills_per_agent',0)},"
                    f"{gs.get('mean_skill_entropy',0)},"
                    f"{gs.get('unique_portfolios',0)},"
                    f"{gs.get('chaining_capable',0)}\n")
    print(f"\nSummary written to: {summary_path}")
    print(f"Log written to: {evolution_log}")


# -- Main --------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Genomebook v2 Evolution Orchestrator")
    parser.add_argument("--generations", type=int, default=10)
    parser.add_argument("--start-gen", type=int, default=0)
    parser.add_argument("--offspring", type=int, default=2)
    parser.add_argument("--population-cap", type=int, default=30)
    parser.add_argument("--task", type=str, default=None,
                        help='Task JSON, e.g. \'{"gene":"BRCA1","rsid":"rs80357906"}\'')
    parser.add_argument("--dry-run", action="store_true",
                        help="Skill gates evaluated but no API calls (zero cost)")
    parser.add_argument("--random-skills", action="store_true",
                        help="Control: randomly assign skills instead of trait gates")
    parser.add_argument("--seed", type=int, default=None)
    args = parser.parse_args()

    task = json.loads(args.task) if args.task else None

    evolve_v2(
        start_gen=args.start_gen,
        num_gens=args.generations,
        offspring_per_pair=args.offspring,
        task=task,
        dry_run=args.dry_run,
        random_skills=args.random_skills,
        seed=args.seed,
        population_cap=args.population_cap,
    )


if __name__ == "__main__":
    main()
