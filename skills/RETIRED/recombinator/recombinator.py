#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ClawBio Recombinator
Produce offspring genomes from parent pairs via meiotic recombination.

Usage:
    python recombinator.py --demo
    python recombinator.py --father einstein-g0 --mother anning-g0 --offspring 3
"""

import argparse
import json
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
_GENOMEBOOK_PYTHON = _PROJECT_ROOT / "GENOMEBOOK" / "PYTHON"

if str(_GENOMEBOOK_PYTHON) not in sys.path:
    sys.path.insert(0, str(_GENOMEBOOK_PYTHON))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))

GENOMES_DIR = _PROJECT_ROOT / "GENOMEBOOK" / "DATA" / "GENOMES"


def _import_recombinator():
    """Import the core recombinator module."""
    from importlib.util import spec_from_file_location, module_from_spec
    spec = spec_from_file_location(
        "recombinator_core",
        _GENOMEBOOK_PYTHON / "04-recombinator.py"
    )
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def load_genome(genome_id):
    """Load a genome by ID from the genomes directory."""
    path = GENOMES_DIR / f"{genome_id}.genome.json"
    if not path.exists():
        print(f"ERROR: Genome not found: {path}")
        sys.exit(1)
    with open(path) as f:
        return json.load(f)


def print_offspring(child):
    """Pretty-print an offspring genome summary."""
    print(f"\n{'=' * 60}")
    print(f"ID:     {child['id']}")
    print(f"Name:   {child['name']}")
    print(f"Sex:    {child['sex']} ({child['sex_chromosomes']})")
    print(f"Health: {child['health_score']}")
    print(f"Mutations: {len(child['mutations'])}")
    if child['mutations']:
        for m in child['mutations']:
            print(f"  - {m['locus']}: {m['from']}->{m['to']} ({m['type']}, from {m['parent']})")
    print(f"Conditions: {len(child['clinical_history'])}")
    for cond in child['clinical_history']:
        print(f"  - {cond['name']} ({cond['severity']}, fitness: {cond['fitness_cost']})")
    print(f"Top traits:")
    top = sorted(child['trait_scores'].items(), key=lambda x: x[1], reverse=True)[:5]
    for t, s in top:
        print(f"  - {t}: {s}")


def run_demo():
    """Demo: breed Einstein x Anning, produce 3 offspring."""
    core = _import_recombinator()
    trait_reg, disease_reg = core.load_registries()

    father = load_genome("einstein-g0")
    mother = load_genome("anning-g0")

    print("=" * 70)
    print("RECOMBINATOR - DEMO MODE")
    print("=" * 70)
    print(f"Father: {father['name']} ({father['id']})")
    print(f"Mother: {mother['name']} ({mother['id']})")
    print(f"Offspring: 3")

    children = core.breed_pair(
        father, mother, generation=1, num_offspring=3,
        trait_reg=trait_reg, disease_reg=disease_reg
    )

    for child in children:
        print_offspring(child)

    print(f"\n[DEMO] {len(children)} offspring generated. No files written.")


def run_breed(father_id, mother_id, num_offspring, generation):
    """Breed a specific pair."""
    core = _import_recombinator()
    trait_reg, disease_reg = core.load_registries()

    father = load_genome(father_id)
    mother = load_genome(mother_id)

    print("=" * 70)
    print("RECOMBINATOR")
    print("=" * 70)
    print(f"Father: {father.get('name', father_id)} ({father_id})")
    print(f"Mother: {mother.get('name', mother_id)} ({mother_id})")
    print(f"Offspring: {num_offspring}, Generation: {generation}")

    children = core.breed_pair(
        father, mother, generation=generation, num_offspring=num_offspring,
        trait_reg=trait_reg, disease_reg=disease_reg
    )

    for child in children:
        print_offspring(child)

        # Write offspring genome
        out_path = GENOMES_DIR / f"{child['id']}.genome.json"
        with open(out_path, "w") as f:
            json.dump(child, f, indent=2)
        print(f"  -> Written: {out_path.name}")

    print(f"\n{len(children)} offspring written to {GENOMES_DIR}/")


def main():
    parser = argparse.ArgumentParser(
        description="Recombinator: Produce offspring via meiotic recombination"
    )
    parser.add_argument("--demo", action="store_true", help="Demo with Einstein x Anning")
    parser.add_argument("--father", type=str, help="Father genome ID (e.g. einstein-g0)")
    parser.add_argument("--mother", type=str, help="Mother genome ID (e.g. anning-g0)")
    parser.add_argument("--offspring", type=int, default=3, help="Number of offspring (default: 3)")
    parser.add_argument("--generation", type=int, default=1, help="Generation number (default: 1)")
    args = parser.parse_args()

    if args.demo:
        run_demo()
    elif args.father and args.mother:
        run_breed(args.father, args.mother, args.offspring, args.generation)
    else:
        parser.print_help()
        print("\nUse --demo for a quick demo, or specify --father and --mother.")


if __name__ == "__main__":
    main()
