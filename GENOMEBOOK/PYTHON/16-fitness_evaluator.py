"""
16-fitness_evaluator.py -- Genomebook v2 Fitness Evaluator

Computes fitness scores from research round outputs. Four components:

  F1: Skill utilization (0.30) - Did the agent use its available skills?
  F2: Synthesis quality  (0.30) - Does the synthesis cite real data?
  F3: Novelty           (0.20) - Unique finding vs other agents this gen?
  F4: Skill chaining    (0.20) - Combined outputs from 2+ skills?

Fitness feeds back into mating selection, creating evolutionary pressure
for agents that use their heritable skills effectively.

Usage:
    # Score a single round result
    python 16-fitness_evaluator.py --round round_result.json

    # Score all agents in a generation
    python 16-fitness_evaluator.py --generation-results gen_results.json

    # Demo: create and score synthetic round results
    python 16-fitness_evaluator.py --demo
"""

import argparse
import json
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"


# -- F1: Skill Utilization ---------------------------------------------------

def score_utilization(round_result):
    """How well did the agent use its available skills for this task?

    Measures: skills_used / skills_that_were_relevant_to_task
    An agent with 5 available skills but only 2 relevant to the task
    gets full marks for using those 2.
    """
    available = set(round_result.get("skills_available", []))
    used = set(round_result.get("skills_used", []))
    task = round_result.get("task", {})

    # Determine which available skills were relevant to this task
    relevant = set()
    for skill in available:
        # Check input_type compatibility
        if skill == "pubmed-search" and "gene" in task:
            relevant.add(skill)
        elif skill == "gwas-variant-lookup" and "rsid" in task:
            relevant.add(skill)
        elif skill == "prs-calculator" and "genotypes" in task:
            relevant.add(skill)
        elif skill == "equity-assessment" and "vcf_path" in task:
            relevant.add(skill)
        elif skill == "clinpgx-lookup" and "gene" in task:
            relevant.add(skill)
        elif skill == "hypothesis-generation" and len(used) > 0:
            relevant.add(skill)
        elif skill == "literature-synthesis" and "pubmed-search" in used:
            relevant.add(skill)

    if not relevant:
        return 0.5  # No relevant skills; neutral score

    used_relevant = relevant & used
    return len(used_relevant) / len(relevant)


# -- F2: Synthesis Quality ---------------------------------------------------

def score_synthesis_quality(round_result):
    """Does the synthesis cite data that actually came from skill outputs?

    For the PoC, this measures:
    - Did skills produce data? (success count)
    - For each successful skill, did it return non-empty data?

    In the full version, this will check whether the LLM synthesis
    references specific data points from skill results.
    """
    skill_results = round_result.get("skill_results", {})
    if not skill_results:
        return 0.0

    successful = 0
    total = 0
    for sname, sresult in skill_results.items():
        data = sresult.get("data")
        # Skip LLM-only skill placeholders
        if isinstance(data, dict) and data.get("type") == "llm_skill":
            continue
        if isinstance(data, dict) and data.get("dry_run"):
            total += 1
            successful += 1  # Dry runs count as successful
            continue

        total += 1
        if sresult.get("success") and sresult.get("data"):
            data = sresult["data"]
            # Check data is non-trivial
            if isinstance(data, list) and len(data) > 0:
                successful += 1
            elif isinstance(data, dict) and len(data) > 0:
                successful += 1

    if total == 0:
        return 0.0
    return successful / total


# -- F3: Novelty -------------------------------------------------------------

def score_novelty(round_result, all_round_results):
    """Did this agent produce a finding no other agent found?

    For PoC: based on skill portfolio uniqueness within this generation.
    Full version: semantic comparison of synthesis claims.
    """
    my_portfolio = tuple(round_result.get("portfolio_vector", []))
    my_skills_used = set(round_result.get("skills_used", []))
    my_id = round_result.get("agent_id")

    if not all_round_results or len(all_round_results) <= 1:
        return 0.5  # Can't compute novelty with < 2 agents

    # Count how many other agents used the exact same skill combination
    same_combo = 0
    for other in all_round_results:
        if other.get("agent_id") == my_id:
            continue
        other_used = set(other.get("skills_used", []))
        if other_used == my_skills_used:
            same_combo += 1

    n_others = len(all_round_results) - 1
    uniqueness = 1.0 - (same_combo / n_others)
    return uniqueness


# -- F4: Skill Chaining ------------------------------------------------------

def score_skill_chaining(round_result):
    """Did the agent combine outputs from 2+ data skills into synthesis?

    This is the key emergent metric. An offspring who inherits both
    pubmed-search (from one parent) and gwas-variant-lookup (from the other)
    can chain them into a finding neither parent could produce alone.
    """
    skills_used = round_result.get("skills_used", [])
    skill_results = round_result.get("skill_results", {})

    # Count data-producing skills (exclude LLM-only skills)
    data_skills = []
    for sname in skills_used:
        sresult = skill_results.get(sname, {})
        data = sresult.get("data", {})
        if isinstance(data, dict) and data.get("type") == "llm_skill":
            continue
        if sresult.get("success"):
            data_skills.append(sname)

    n_data = len(data_skills)

    if n_data == 0:
        return 0.0
    elif n_data == 1:
        return 0.2  # Single skill, minimal chaining
    elif n_data == 2:
        return 0.6  # Two data skills combined
    elif n_data == 3:
        return 0.8  # Three data skills
    else:
        return 1.0  # Four or more data skills chained

    # Bonus: was a synthesis or hypothesis skill also used?
    # (meaning the agent actually combined the data outputs)
    llm_skills = [s for s in skills_used
                  if s in ("literature-synthesis", "hypothesis-generation")]
    if llm_skills and n_data >= 2:
        return min(1.0, score + 0.2)

    return score


# -- Composite Fitness -------------------------------------------------------

WEIGHTS = {
    "utilization": 0.30,
    "synthesis_quality": 0.30,
    "novelty": 0.20,
    "skill_chaining": 0.20,
}


def compute_fitness(round_result, all_round_results=None):
    """Compute composite fitness score for one agent's round.

    Returns dict with component scores and composite.
    """
    f1 = score_utilization(round_result)
    f2 = score_synthesis_quality(round_result)
    f3 = score_novelty(round_result, all_round_results or [])
    f4 = score_skill_chaining(round_result)

    composite = (WEIGHTS["utilization"] * f1 +
                 WEIGHTS["synthesis_quality"] * f2 +
                 WEIGHTS["novelty"] * f3 +
                 WEIGHTS["skill_chaining"] * f4)

    return {
        "agent_id": round_result.get("agent_id"),
        "f1_utilization": round(f1, 4),
        "f2_synthesis_quality": round(f2, 4),
        "f3_novelty": round(f3, 4),
        "f4_skill_chaining": round(f4, 4),
        "composite_fitness": round(composite, 4),
        "num_skills_available": len(round_result.get("skills_available", [])),
        "num_skills_used": len(round_result.get("skills_used", [])),
    }


def compute_generation_fitness(all_round_results):
    """Score all agents in a generation, computing novelty across the group."""
    fitness_scores = []
    for rr in all_round_results:
        fit = compute_fitness(rr, all_round_results)
        fitness_scores.append(fit)
    return fitness_scores


def generation_stats(fitness_scores):
    """Compute population-level fitness statistics."""
    if not fitness_scores:
        return {}

    composites = [f["composite_fitness"] for f in fitness_scores]
    mean = sum(composites) / len(composites)
    variance = sum((c - mean) ** 2 for c in composites) / len(composites)
    stdev = variance ** 0.5

    return {
        "population_size": len(fitness_scores),
        "mean_fitness": round(mean, 4),
        "stdev_fitness": round(stdev, 4),
        "min_fitness": round(min(composites), 4),
        "max_fitness": round(max(composites), 4),
        "mean_f1": round(sum(f["f1_utilization"] for f in fitness_scores) / len(fitness_scores), 4),
        "mean_f2": round(sum(f["f2_synthesis_quality"] for f in fitness_scores) / len(fitness_scores), 4),
        "mean_f3": round(sum(f["f3_novelty"] for f in fitness_scores) / len(fitness_scores), 4),
        "mean_f4": round(sum(f["f4_skill_chaining"] for f in fitness_scores) / len(fitness_scores), 4),
    }


# -- Main --------------------------------------------------------------------

def _demo():
    """Create and score synthetic round results for demonstration."""
    # Simulate 5 agents with different skill portfolios
    agents = [
        {"agent_id": "agent-A", "skills_available": ["pubmed-search", "gwas-variant-lookup", "prs-calculator"],
         "skills_used": ["pubmed-search", "gwas-variant-lookup"],
         "skill_results": {
             "pubmed-search": {"skill": "pubmed-search", "success": True,
                               "data": [{"title": "paper1"}, {"title": "paper2"}], "error": None},
             "gwas-variant-lookup": {"skill": "gwas-variant-lookup", "success": True,
                                     "data": {"merged": {"associations": []}}, "error": None},
         },
         "task": {"gene": "BRCA1", "rsid": "rs80357906"},
         "portfolio_vector": [1, 0, 1, 0, 0, 1, 1], "skill_entropy": 0.86, "num_skills": 4},

        {"agent_id": "agent-B", "skills_available": ["pubmed-search", "equity-assessment"],
         "skills_used": ["pubmed-search"],
         "skill_results": {
             "pubmed-search": {"skill": "pubmed-search", "success": True,
                               "data": [{"title": "paper1"}], "error": None},
         },
         "task": {"gene": "BRCA1", "rsid": "rs80357906"},
         "portfolio_vector": [0, 1, 0, 0, 0, 0, 1], "skill_entropy": 0.59, "num_skills": 2},

        {"agent_id": "agent-C", "skills_available": ["gwas-variant-lookup", "hypothesis-generation"],
         "skills_used": ["gwas-variant-lookup", "hypothesis-generation"],
         "skill_results": {
             "gwas-variant-lookup": {"skill": "gwas-variant-lookup", "success": True,
                                     "data": {"merged": {"associations": [1, 2]}}, "error": None},
             "hypothesis-generation": {"skill": "hypothesis-generation", "success": True,
                                       "data": {"type": "llm_skill", "requires_synthesis": True}, "error": None},
         },
         "task": {"gene": "BRCA1", "rsid": "rs80357906"},
         "portfolio_vector": [0, 0, 1, 1, 0, 0, 0], "skill_entropy": 0.59, "num_skills": 2},

        {"agent_id": "agent-D", "skills_available": ["pubmed-search"],
         "skills_used": ["pubmed-search"],
         "skill_results": {
             "pubmed-search": {"skill": "pubmed-search", "success": True,
                               "data": [{"title": "p1"}, {"title": "p2"}, {"title": "p3"}], "error": None},
         },
         "task": {"gene": "BRCA1"},
         "portfolio_vector": [0, 0, 0, 0, 0, 0, 1], "skill_entropy": 0.41, "num_skills": 1},

        {"agent_id": "agent-E", "skills_available": ["pubmed-search", "gwas-variant-lookup",
                                                      "clinpgx-lookup", "literature-synthesis"],
         "skills_used": ["pubmed-search", "gwas-variant-lookup", "clinpgx-lookup", "literature-synthesis"],
         "skill_results": {
             "pubmed-search": {"skill": "pubmed-search", "success": True,
                               "data": [{"title": "p1"}, {"title": "p2"}], "error": None},
             "gwas-variant-lookup": {"skill": "gwas-variant-lookup", "success": True,
                                     "data": {"merged": {}}, "error": None},
             "clinpgx-lookup": {"skill": "clinpgx-lookup", "success": True,
                                "data": {"gene": "BRCA1", "drugs": []}, "error": None},
             "literature-synthesis": {"skill": "literature-synthesis", "success": True,
                                      "data": {"type": "llm_skill", "requires_synthesis": True}, "error": None},
         },
         "task": {"gene": "BRCA1", "rsid": "rs80357906"},
         "portfolio_vector": [1, 0, 1, 0, 1, 0, 1], "skill_entropy": 0.86, "num_skills": 4},
    ]

    print("Demo: Scoring 5 synthetic agent rounds\n")
    fitness_scores = compute_generation_fitness(agents)

    print(f"{'Agent':10s} | {'F1 Util':>8s} {'F2 Synth':>8s} {'F3 Novel':>8s} {'F4 Chain':>8s} | {'Fitness':>8s}")
    print("-" * 70)
    for f in fitness_scores:
        print(f"{f['agent_id']:10s} | {f['f1_utilization']:8.3f} {f['f2_synthesis_quality']:8.3f} "
              f"{f['f3_novelty']:8.3f} {f['f4_skill_chaining']:8.3f} | {f['composite_fitness']:8.3f}")

    stats = generation_stats(fitness_scores)
    print(f"\nPopulation stats:")
    for k, v in stats.items():
        print(f"  {k}: {v}")


def main():
    parser = argparse.ArgumentParser(description="Genomebook v2 Fitness Evaluator")
    parser.add_argument("--demo", action="store_true", help="Run demo scoring")
    parser.add_argument("--round", type=str, help="Score a single round result JSON")
    parser.add_argument("--generation-results", type=str,
                        help="Score all agents from a generation results JSON")
    args = parser.parse_args()

    if args.demo:
        _demo()
        return

    if args.round:
        with open(args.round) as f:
            rr = json.load(f)
        fit = compute_fitness(rr)
        print(json.dumps(fit, indent=2))
        return

    if args.generation_results:
        with open(args.generation_results) as f:
            all_rr = json.load(f)
        scores = compute_generation_fitness(all_rr)
        stats = generation_stats(scores)
        print(json.dumps({"scores": scores, "stats": stats}, indent=2))
        return

    _demo()


if __name__ == "__main__":
    main()
