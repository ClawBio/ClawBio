"""
15-skill_executor.py -- Genomebook v2 Skill Executor

Evaluates trait gates from skill_registry.json against agent genomes,
dispatches to ClawBio skill APIs, and collects structured results.

The key insight: because trait scores derive from diploid loci via Mendelian
inheritance, skill access is heritable. Offspring of a PubMed-specialist
and a GWAS-specialist may inherit both skills or neither.

Usage:
    python 15-skill_executor.py --matrix              # Print skill matrix for gen-0
    python 15-skill_executor.py --agent einstein-g0   # Show one agent's skills
    python 15-skill_executor.py --agent einstein-g0 --task '{"gene":"BRCA1"}'  # Execute
    python 15-skill_executor.py --generation 0 --task '{"gene":"BRCA1"}' --dry-run
"""

import argparse
import json
import math
import os
import sys
from pathlib import Path

BASE = Path(__file__).resolve().parent.parent
DATA = BASE / "DATA"
SKILL_REGISTRY_PATH = DATA / "skill_registry.json"
GENOMES_DIR = DATA / "GENOMES"
GENERATIONS_DIR = DATA / "GENERATIONS"

CLAWBIO_ROOT = BASE.parent  # 02-APPS/06-CLAWBIO
ANTHROPIC_API_KEY = os.environ.get("ANTHROPIC_API_KEY", "")


# -- Registry loading -------------------------------------------------------

def load_skill_registry():
    with open(SKILL_REGISTRY_PATH) as f:
        return json.load(f)


def load_genome(agent_id):
    path = GENOMES_DIR / f"{agent_id}.genome.json"
    with open(path) as f:
        return json.load(f)


def load_generation_agents(generation):
    gen_path = GENERATIONS_DIR / f"gen-{generation:04d}" / "agents.json"
    if gen_path.exists():
        with open(gen_path) as f:
            return json.load(f)
    # Fallback: scan GENOMES_DIR
    agents = []
    for gf in sorted(GENOMES_DIR.glob("*.genome.json")):
        g = json.load(open(gf))
        if g["generation"] == generation:
            agents.append(g)
    return agents


# -- Gate evaluation ---------------------------------------------------------

def evaluate_gates(trait_scores, skill_registry):
    """Return dict of {skill_name: {available: bool, config: dict, missing: list}}.

    For each skill, check if all gate traits meet thresholds.
    If available, compute configuration from config_traits formulas.
    """
    results = {}
    for skill_name, skill_def in skill_registry["skills"].items():
        gates = skill_def["gate_traits"]
        missing = []
        for trait, threshold in gates.items():
            agent_score = trait_scores.get(trait, 0.0)
            if agent_score < threshold:
                missing.append({
                    "trait": trait,
                    "required": threshold,
                    "actual": round(agent_score, 3),
                })

        available = len(missing) == 0

        # Compute config from trait scores
        config = {}
        if available:
            for param_name, param_def in skill_def.get("config_traits", {}).items():
                trait = param_def.get("trait")
                formula = param_def.get("formula")
                if trait and formula:
                    score = trait_scores.get(trait, 0.5)
                    try:
                        config[param_name] = eval(formula, {"__builtins__": {}},
                                                  {"score": score, "int": int,
                                                   "round": round, "math": math})
                    except Exception:
                        config[param_name] = None
                elif trait:
                    config[param_name] = round(trait_scores.get(trait, 0.5), 3)

        results[skill_name] = {
            "available": available,
            "config": config,
            "missing": missing,
        }
    return results


def available_skills(trait_scores, skill_registry):
    """Return only the skills an agent can use, with their configs."""
    all_gates = evaluate_gates(trait_scores, skill_registry)
    return {k: v for k, v in all_gates.items() if v["available"]}


def skill_portfolio_vector(trait_scores, skill_registry):
    """Return binary vector [0/1] for each skill. Used for diversity metrics."""
    gates = evaluate_gates(trait_scores, skill_registry)
    skill_names = sorted(skill_registry["skills"].keys())
    return [1 if gates[s]["available"] else 0 for s in skill_names]


def skill_entropy(portfolio_vector):
    """Shannon entropy of a binary skill vector. Measures specialization."""
    n = len(portfolio_vector)
    if n == 0:
        return 0.0
    p = sum(portfolio_vector) / n
    if p == 0.0 or p == 1.0:
        return 0.0
    return -(p * math.log2(p) + (1 - p) * math.log2(1 - p))


# -- Skill dispatch ----------------------------------------------------------

def _load_skill_module(module_path):
    """Import a ClawBio skill module by path."""
    full_path = CLAWBIO_ROOT / module_path
    if not full_path.exists():
        return None

    skill_dir = str(full_path.parent)
    if skill_dir not in sys.path:
        sys.path.insert(0, skill_dir)

    # Ensure project root is importable for clawbio.common
    project_root = str(CLAWBIO_ROOT)
    if project_root not in sys.path:
        sys.path.insert(0, project_root)

    from importlib.util import spec_from_file_location, module_from_spec
    mod_name = full_path.stem.replace("-", "_")
    spec = spec_from_file_location(mod_name, str(full_path))
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def execute_skill(skill_name, skill_def, config, task):
    """Execute a single ClawBio skill and return structured results.

    Returns dict with: skill, success, data, error
    """
    module_path = skill_def.get("module_path")
    func_name = skill_def.get("function")

    # LLM-only skills (hypothesis-generation, literature-synthesis)
    if module_path is None:
        return {
            "skill": skill_name,
            "success": True,
            "data": {"type": "llm_skill", "requires_synthesis": True},
            "error": None,
        }

    try:
        mod = _load_skill_module(module_path)
        if mod is None:
            return {
                "skill": skill_name,
                "success": False,
                "data": None,
                "error": f"Module not found: {module_path}",
            }

        func = getattr(mod, func_name, None)
        if func is None:
            return {
                "skill": skill_name,
                "success": False,
                "data": None,
                "error": f"Function {func_name} not found in {module_path}",
            }

        # Dispatch based on skill type
        if skill_name == "pubmed-search":
            gene = task.get("gene", "BRCA1")
            max_results = config.get("max_results", 10)
            data = func(query=gene, max_results=max_results)

        elif skill_name == "gwas-variant-lookup":
            rsid = task.get("rsid")
            if not rsid:
                return {"skill": skill_name, "success": False,
                        "data": None, "error": "No rsid in task"}
            options = {
                "rsid": rsid,
                "skip_apis": config.get("skip_apis", []),
                "max_hits": config.get("max_hits", 100),
                "use_cache": True,
                "demo": task.get("demo", False),
            }
            data = func(genotypes=None, options=options)

        elif skill_name == "prs-calculator":
            genotypes = task.get("genotypes", {})
            if not genotypes:
                return {"skill": skill_name, "success": False,
                        "data": None, "error": "No genotypes in task"}
            options = {"min_overlap": config.get("min_overlap", 0.5)}
            data = func(genotypes=genotypes, options=options)

        elif skill_name == "equity-assessment":
            vcf_path = task.get("vcf_path")
            if not vcf_path:
                return {"skill": skill_name, "success": False,
                        "data": None, "error": "No vcf_path in task"}
            options = {"weights": config.get("weights", "0.35,0.25,0.20,0.20")}
            data = func(input_path=vcf_path, options=options)

        elif skill_name == "clinpgx-lookup":
            gene = task.get("gene")
            if not gene:
                return {"skill": skill_name, "success": False,
                        "data": None, "error": "No gene in task"}
            # clinpgx requires a client instance
            from pathlib import Path as _P
            cache_dir = _P(os.environ.get("CLAWBIO_CACHE",
                                          str(CLAWBIO_ROOT / "cache" / "clinpgx")))
            cache_dir.mkdir(parents=True, exist_ok=True)
            client = mod.ClinPGxClient(cache_dir=cache_dir)
            data = func(client, gene)

        else:
            return {"skill": skill_name, "success": False,
                    "data": None, "error": f"Unknown dispatch for {skill_name}"}

        return {
            "skill": skill_name,
            "success": True,
            "data": data,
            "error": None,
        }

    except Exception as e:
        return {
            "skill": skill_name,
            "success": False,
            "data": None,
            "error": str(e),
        }


def execute_round(agent, task, skill_registry, dry_run=False):
    """Execute a full research round for one agent.

    Returns:
        dict with keys:
            agent_id, skills_available, skills_used, skill_results,
            task, portfolio_vector, skill_entropy
    """
    traits = agent.get("trait_scores", {})
    avail = available_skills(traits, skill_registry)
    portfolio = skill_portfolio_vector(traits, skill_registry)
    entropy = skill_entropy(portfolio)

    skill_results = {}
    skills_used = []

    for skill_name, gate_info in avail.items():
        skill_def = skill_registry["skills"][skill_name]
        input_type = skill_def.get("input_type", "")

        # Skip skills that don't match the task
        if input_type == "rsid" and "rsid" not in task:
            continue
        if input_type == "gene_name" and "gene" not in task:
            continue
        if input_type == "genotypes" and "genotypes" not in task:
            continue
        if input_type == "vcf_path" and "vcf_path" not in task:
            continue
        if input_type == "synthesis_context":
            continue  # LLM skills run after data skills
        if input_type == "paper_list":
            continue  # Runs after pubmed-search

        if dry_run:
            skill_results[skill_name] = {
                "skill": skill_name,
                "success": True,
                "data": {"dry_run": True},
                "error": None,
            }
            skills_used.append(skill_name)
        else:
            result = execute_skill(skill_name, skill_def, gate_info["config"], task)
            skill_results[skill_name] = result
            if result["success"]:
                skills_used.append(skill_name)

    # LLM skills: literature-synthesis runs if we got pubmed results
    if "literature-synthesis" in avail and "pubmed-search" in skill_results:
        pm_result = skill_results["pubmed-search"]
        if pm_result["success"] and pm_result["data"]:
            skill_results["literature-synthesis"] = {
                "skill": "literature-synthesis",
                "success": True,
                "data": {"type": "llm_skill", "requires_synthesis": True,
                         "input_papers": len(pm_result["data"]) if isinstance(pm_result["data"], list) else 0},
                "error": None,
            }
            skills_used.append("literature-synthesis")

    # hypothesis-generation runs if we have any data skill results
    data_skills_used = [s for s in skills_used if s not in ("literature-synthesis", "hypothesis-generation")]
    if "hypothesis-generation" in avail and len(data_skills_used) >= 1:
        skill_results["hypothesis-generation"] = {
            "skill": "hypothesis-generation",
            "success": True,
            "data": {"type": "llm_skill", "requires_synthesis": True,
                     "input_skills": data_skills_used},
            "error": None,
        }
        skills_used.append("hypothesis-generation")

    return {
        "agent_id": agent.get("id", "unknown"),
        "agent_name": agent.get("name", "unknown"),
        "generation": agent.get("generation", 0),
        "skills_available": list(avail.keys()),
        "skills_used": skills_used,
        "skill_results": skill_results,
        "task": task,
        "portfolio_vector": portfolio,
        "skill_entropy": round(entropy, 4),
        "num_skills": sum(portfolio),
    }


# -- Display helpers ---------------------------------------------------------

def print_skill_matrix(agents, skill_registry):
    """Print a matrix of agents x skills showing availability."""
    skill_names = sorted(skill_registry["skills"].keys())
    abbrevs = {
        "pubmed-search": "PM",
        "gwas-variant-lookup": "GW",
        "prs-calculator": "PR",
        "equity-assessment": "EQ",
        "clinpgx-lookup": "CL",
        "hypothesis-generation": "HY",
        "literature-synthesis": "LS",
    }

    header = f"{'Agent':20s} |"
    for s in skill_names:
        header += f" {abbrevs.get(s, s[:2]):>3s}"
    header += " | Total"
    print(header)
    print("-" * len(header))

    pop_counts = {s: 0 for s in skill_names}
    for a in agents:
        traits = a.get("trait_scores", {})
        gates = evaluate_gates(traits, skill_registry)
        row = f"{a['id']:20s} |"
        total = 0
        for s in skill_names:
            if gates[s]["available"]:
                row += "  Y "
                total += 1
                pop_counts[s] += 1
            else:
                row += "  . "
        row += f" | {total}/7"
        print(row)

    print()
    print("Population coverage:")
    for s in skill_names:
        print(f"  {s:25s}: {pop_counts[s]}/{len(agents)}")


# -- Main --------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(description="Genomebook v2 Skill Executor")
    parser.add_argument("--matrix", action="store_true",
                        help="Print skill matrix for a generation")
    parser.add_argument("--generation", type=int, default=0,
                        help="Generation to analyse (default 0)")
    parser.add_argument("--agent", type=str,
                        help="Single agent ID to inspect or execute")
    parser.add_argument("--task", type=str,
                        help='Task JSON, e.g. \'{"gene":"BRCA1"}\'')
    parser.add_argument("--dry-run", action="store_true",
                        help="Evaluate gates and mock execution (no API calls)")
    args = parser.parse_args()

    registry = load_skill_registry()

    if args.matrix:
        agents = load_generation_agents(args.generation)
        if not agents:
            print(f"No agents found for generation {args.generation}")
            return
        print(f"\nSkill Matrix: Generation {args.generation} ({len(agents)} agents)\n")
        print_skill_matrix(agents, registry)
        return

    if args.agent:
        genome = load_genome(args.agent)
        traits = genome.get("trait_scores", {})
        gates = evaluate_gates(traits, registry)

        print(f"\nAgent: {args.agent} ({genome.get('name', '?')})")
        print(f"Generation: {genome.get('generation', '?')}")
        print(f"Sex: {genome.get('sex', '?')}")
        print()

        for skill_name in sorted(gates.keys()):
            info = gates[skill_name]
            status = "AVAILABLE" if info["available"] else "LOCKED"
            print(f"  [{status:9s}] {skill_name}")
            if info["available"] and info["config"]:
                for k, v in info["config"].items():
                    print(f"              {k}: {v}")
            elif info["missing"]:
                for m in info["missing"]:
                    print(f"              needs {m['trait']}: {m['actual']:.2f} < {m['required']:.2f}")

        portfolio = skill_portfolio_vector(traits, registry)
        print(f"\n  Portfolio: {sum(portfolio)}/7 skills")
        print(f"  Entropy:   {skill_entropy(portfolio):.4f}")

        if args.task:
            task = json.loads(args.task)
            print(f"\n  Executing task: {task}")
            result = execute_round(genome, task, registry, dry_run=args.dry_run)
            print(f"  Skills used: {result['skills_used']}")
            for sname, sresult in result["skill_results"].items():
                status = "OK" if sresult["success"] else "FAIL"
                print(f"    [{status}] {sname}")
                if sresult["error"]:
                    print(f"           Error: {sresult['error']}")

        return

    # Default: show matrix for generation 0
    agents = load_generation_agents(0)
    if agents:
        print(f"\nSkill Matrix: Generation 0 ({len(agents)} agents)\n")
        print_skill_matrix(agents, registry)


if __name__ == "__main__":
    main()
