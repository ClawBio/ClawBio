#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
ClawBio GenomeMatch
Score genetic compatibility across all M x F pairings in a Genomebook generation.

Usage:
    python genome_match.py              # Score generation 0
    python genome_match.py --demo       # Demo mode
    python genome_match.py --generation 1 --top 10
"""

import argparse
import sys
from pathlib import Path

_PROJECT_ROOT = Path(__file__).resolve().parent.parent.parent
_GENOMEBOOK_PYTHON = _PROJECT_ROOT / "GENOMEBOOK" / "PYTHON"

if str(_GENOMEBOOK_PYTHON) not in sys.path:
    sys.path.insert(0, str(_GENOMEBOOK_PYTHON))
if str(_PROJECT_ROOT) not in sys.path:
    sys.path.insert(0, str(_PROJECT_ROOT))


def _import_genomematch():
    """Import the core genomematch module."""
    from importlib.util import spec_from_file_location, module_from_spec
    spec = spec_from_file_location(
        "genomematch_core",
        _GENOMEBOOK_PYTHON / "02-genomematch.py"
    )
    mod = module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod


def run_demo():
    """Run demo with generation-0 genomes."""
    run_match(generation=0, top_n=20, demo=True)


def run_match(generation=0, top_n=None, demo=False):
    """Run the matching engine."""
    core = _import_genomematch()

    genomes = core.load_genomes(generation=generation)
    if not genomes:
        print(f"ERROR: No generation-{generation} genomes found.")
        sys.exit(1)

    disease_reg = core.load_disease_registry()
    males = [g for g in genomes.values() if g["sex"] == "Male"]
    females = [g for g in genomes.values() if g["sex"] == "Female"]

    label = "DEMO MODE" if demo else f"Generation {generation}"
    print("=" * 100)
    print(f"GENOMEMATCH COMPATIBILITY ENGINE - {label}")
    print("=" * 100)
    print(f"Loaded {len(genomes)} genomes ({len(males)}M / {len(females)}F)")

    pairings = core.match_generation(genomes, disease_reg)
    display = pairings[:top_n] if top_n else pairings[:20]

    print(f"\nAll {len(pairings)} M x F pairings scored.\n")
    print(f"{'Rank':>4}  {'Male':>20} x {'Female':<20}  {'Score':>6}  {'Het':>5}  {'Comp':>5}  {'Risk':>5}  Flags")
    print("-" * 100)

    for i, p in enumerate(display, 1):
        flags = ", ".join(d["disease"] for d in p["flagged_diseases"]) or "--"
        print(
            f"{i:4d}  {p['male']:>20} x {p['female']:<20}  "
            f"{p['score']:6.4f}  {p['heterozygosity']:5.3f}  "
            f"{p['complementarity']:5.3f}  {p['disease_risk']:5.3f}  {flags}"
        )

    selected = core.select_mating_pairs(pairings, max_pairs=10)
    print(f"\n{'=' * 60}")
    print(f"SELECTED MATING PAIRS (generation {generation} -> {generation + 1}):")
    print(f"{'=' * 60}")
    for p in selected:
        print(f"  {p['male_name']} x {p['female_name']}  (compat: {p['score']:.4f})")

    if demo:
        print(f"\n[DEMO] Showing generation-0 results.")


def main():
    parser = argparse.ArgumentParser(
        description="GenomeMatch: Score genetic compatibility across Genomebook pairings"
    )
    parser.add_argument("--demo", action="store_true", help="Run demo with generation-0 genomes")
    parser.add_argument("--generation", type=int, default=0, help="Generation to score (default: 0)")
    parser.add_argument("--top", type=int, default=None, help="Show top N pairings only")
    args = parser.parse_args()

    if args.demo:
        run_demo()
    else:
        run_match(generation=args.generation, top_n=args.top)


if __name__ == "__main__":
    main()
