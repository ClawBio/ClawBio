#!/usr/bin/env python3
"""
generate_patient.py — Random synthetic patient generator for NutriGx Advisor demos.

Generates a realistic 23andMe-format CSV file with randomly assigned genotypes
for all SNPs in the NutriGx panel. Each run produces a different risk profile,
demonstrating the diversity of personalised nutrition reports.

Usage:
    # Generate one random patient
    python examples/generate_patient.py

    # Generate with a specific seed (reproducible)
    python examples/generate_patient.py --seed 42

    # Generate and immediately run the report
    python examples/generate_patient.py --run

    # Generate N patients at once (for batch demos)
    python examples/generate_patient.py --count 5

Output:
    examples/data/patient_<seed>.csv
    (or nutrigx_results_<seed>/ if --run is used)
"""

import argparse
import json
import random
import subprocess
import sys
from datetime import date
from pathlib import Path

# ── Population allele frequencies ─────────────────────────────────────────────
# Approximate minor allele frequencies from gnomAD v4 (global population).
# Used to weight random genotype assignment so the distribution is realistic
# rather than purely uniform.
#
# Format: rsid -> risk_allele_frequency (0–1)
# Where 1.0 means the risk allele is always present (monomorphic in panel),
# and 0.0 means it is never present.
#
# SNPs not listed here default to 0.30 (typical common variant MAF).

RISK_ALLELE_FREQS = {
    "rs1801133":  0.32,   # MTHFR C677T  — very common
    "rs1801131":  0.20,   # MTHFR A1298C
    "rs1805087":  0.17,   # MTR A2756G
    "rs2228570":  0.42,   # VDR FokI
    "rs731236":   0.38,   # VDR TaqI
    "rs4588":     0.28,   # GC Thr436Lys
    "rs174546":   0.47,   # FADS1
    "rs1535":     0.38,   # FADS2
    "rs953413":   0.44,   # ELOVL2
    "rs429358":   0.15,   # APOE ε4 allele
    "rs7412":     0.08,   # APOE ε2 allele
    "rs7501331":  0.24,   # BCMO1
    "rs12934922": 0.42,   # BCMO1
    "rs33972313": 0.05,   # SLC23A1 (rare)
    "rs1256335":  0.35,   # ALPL
    "rs9939609":  0.43,   # FTO
    "rs7903146":  0.29,   # TCF7L2
    "rs1801282":  0.12,   # PPARG Pro12Ala
    "rs662799":   0.09,   # APOA5
    "rs762551":   0.31,   # CYP1A2 slow metaboliser
    "rs4410790":  0.41,   # AHR
    "rs1229984":  0.05,   # ADH1B (mainly East Asian: ~0.70; European: ~0.05)
    "rs671":      0.00,   # ALDH2 (essentially absent in Europeans; ~0.30 East Asian)
    "rs4988235":  0.35,   # MCM6 lactase non-persistence
    "rs4880":     0.47,   # SOD2
    "rs1050450":  0.28,   # GPX1
    "rs1800566":  0.22,   # NQO1
    "rs4680":     0.49,   # COMT Val158Met
}

# Chromosomal positions (GRCh37) — for the 23andMe header format
POSITIONS = {
    "rs1801133":  ("1",  "11854476"),
    "rs1801131":  ("1",  "11856378"),
    "rs1805087":  ("1",  "237048500"),
    "rs2228570":  ("12", "48272895"),
    "rs731236":   ("12", "48239835"),
    "rs4588":     ("4",  "72608790"),
    "rs174546":   ("11", "61327359"),
    "rs1535":     ("11", "61311797"),
    "rs953413":   ("6",  "11044620"),
    "rs429358":   ("19", "45411941"),
    "rs7412":     ("19", "45412079"),
    "rs7501331":  ("16", "81274254"),
    "rs12934922": ("16", "81271111"),
    "rs33972313": ("5",  "110032987"),
    "rs1256335":  ("1",  "21882173"),
    "rs9939609":  ("16", "53820527"),
    "rs7903146":  ("10", "114758349"),
    "rs1801282":  ("3",  "12393125"),
    "rs662799":   ("11", "116700773"),
    "rs762551":   ("15", "75041917"),
    "rs4410790":  ("7",  "17381394"),
    "rs1229984":  ("4",  "100239319"),
    "rs671":      ("12", "111803962"),
    "rs4988235":  ("2",  "136616754"),
    "rs4880":     ("6",  "160113872"),
    "rs1050450":  ("3",  "49394834"),
    "rs1800566":  ("16", "69748869"),
    "rs4680":     ("22", "19951271"),
}

# Reference allele for each rsid (the non-risk base)
REF_ALLELES = {
    "rs1801133":  "C",
    "rs1801131":  "A",
    "rs1805087":  "A",
    "rs2228570":  "T",
    "rs731236":   "T",
    "rs4588":     "C",
    "rs174546":   "T",
    "rs1535":     "G",
    "rs953413":   "A",
    "rs429358":   "T",
    "rs7412":     "C",
    "rs7501331":  "C",
    "rs12934922": "A",
    "rs33972313": "C",
    "rs1256335":  "C",
    "rs9939609":  "T",
    "rs7903146":  "C",
    "rs1801282":  "C",
    "rs662799":   "T",
    "rs762551":   "A",
    "rs4410790":  "C",
    "rs1229984":  "G",
    "rs671":      "G",
    "rs4988235":  "G",
    "rs4880":     "T",
    "rs1050450":  "C",
    "rs1800566":  "C",
    "rs4680":     "G",
}

RISK_ALLELES = {
    "rs1801133":  "T",
    "rs1801131":  "C",
    "rs1805087":  "G",
    "rs2228570":  "C",
    "rs731236":   "C",
    "rs4588":     "A",
    "rs174546":   "C",
    "rs1535":     "A",
    "rs953413":   "G",
    "rs429358":   "C",
    "rs7412":     "T",
    "rs7501331":  "T",
    "rs12934922": "T",
    "rs33972313": "T",
    "rs1256335":  "T",
    "rs9939609":  "A",
    "rs7903146":  "T",
    "rs1801282":  "G",
    "rs662799":   "C",
    "rs762551":   "C",
    "rs4410790":  "T",
    "rs1229984":  "A",
    "rs671":      "A",
    "rs4988235":  "A",
    "rs4880":     "C",
    "rs1050450":  "T",
    "rs1800566":  "T",
    "rs4680":     "A",
}


def assign_genotype(rsid: str, rng: random.Random) -> str:
    """
    Assign a genotype for one SNP using Hardy-Weinberg equilibrium probabilities
    derived from the risk allele frequency.

    Returns a two-character genotype string, e.g. 'CT', 'TT', 'CC'.
    Homozygous risk and homozygous ref are both possible; het is most likely
    at intermediate frequencies.
    """
    q = RISK_ALLELE_FREQS.get(rsid, 0.30)   # risk allele frequency
    p = 1.0 - q                               # ref allele frequency

    # Hardy-Weinberg: p², 2pq, q²
    hom_ref  = p * p
    het      = 2 * p * q
    # hom_risk = q * q  (implicit: the remainder)

    r = rng.random()
    ref = REF_ALLELES.get(rsid, "A")
    risk = RISK_ALLELES.get(rsid, "T")

    if r < hom_ref:
        return ref + ref
    elif r < hom_ref + het:
        # Randomise allele order to mimic real raw data
        alleles = [ref, risk]
        rng.shuffle(alleles)
        return "".join(alleles)
    else:
        return risk + risk


def generate_csv(seed: int, output_path: Path) -> Path:
    """Generate one synthetic patient CSV and write to output_path."""
    rng = random.Random(seed)
    today = date.today().strftime("%a %b %d 00:00:00 %Y")

    rsids = list(POSITIONS.keys())

    lines = [
        f"# This data file generated by 23andMe at {today}",
        "# This file contains raw genotype data, including data that is not used in 23andMe reports.",
        f"# SYNTHETIC DATA — generated by NutriGx generate_patient.py (seed={seed})",
        "# NOT REAL PATIENT DATA.",
        "#",
        "# rsid\tchromosome\tposition\tgenotype",
    ]

    for rsid in rsids:
        chrom, pos = POSITIONS[rsid]
        genotype = assign_genotype(rsid, rng)
        lines.append(f"{rsid}\t{chrom}\t{pos}\t{genotype}")

    output_path.write_text("\n".join(lines) + "\n")
    return output_path


def main():
    parser = argparse.ArgumentParser(
        description="Generate random synthetic patient file(s) for NutriGx demos"
    )
    parser.add_argument(
        "--seed", type=int, default=None,
        help="Random seed for reproducibility (default: random)"
    )
    parser.add_argument(
        "--count", type=int, default=1,
        help="Number of patients to generate (default: 1)"
    )
    parser.add_argument(
        "--run", action="store_true",
        help="Immediately run NutriGx Advisor on each generated file"
    )
    parser.add_argument(
        "--no-figures", action="store_true",
        help="Skip figures when --run is used"
    )
    args = parser.parse_args()

    skill_root = Path(__file__).parent.parent
    examples_data = Path(__file__).parent / "data"
    examples_data.mkdir(exist_ok=True)

    seeds = []
    for i in range(args.count):
        s = (args.seed + i) if args.seed is not None else random.randint(0, 999_999)
        seeds.append(s)

    for seed in seeds:
        out_csv = examples_data / f"patient_{seed}.csv"
        generate_csv(seed, out_csv)
        print(f"[generate_patient] Written: {out_csv}")

        if args.run:
            out_dir = Path(__file__).parent / "output" / f"results_{seed}"
            cmd = [
                sys.executable,
                str(skill_root / "nutrigx_advisor.py"),
                "--input", str(out_csv),
                "--output", str(out_dir),
            ]
            if args.no_figures:
                cmd.append("--no-figures")
            print(f"[generate_patient] Running: {' '.join(cmd)}")
            subprocess.run(cmd, check=True)


if __name__ == "__main__":
    main()
