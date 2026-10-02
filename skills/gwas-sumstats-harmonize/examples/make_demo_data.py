#!/usr/bin/env python3
"""
make_demo_data.py — regenerate the synthetic demo inputs for gwas-sumstats-harmonize.

Everything here is simulated from a fixed seed: no real cohort, sample or
variant data. One synthetic "truth" panel (reference.tsv) is rendered into four
real-world summary-statistics conventions, each seeded with the problems the
harmonizer must handle:

  cohort_ssf.tsv              GWAS-Catalog SSF columns; ~20% alleles swapped vs the
                              reference, some strand-flipped (palindromic A/T, C/G
                              rows included, so only EAF can resolve them), some
                              EAF missing, a duplicate row, invalid p-values
  cohort_plink2.glm.logistic  PLINK2 REF/ALT/A1 trio (A1 is sometimes REF), OR only
  cohort_metal.tbl            METAL: no CHR/BP (MarkerName "chr:pos"), lowercase
                              alleles
  cohort_regenie.regenie      REGENIE ALLELE0/ALLELE1 with LOG10P only, including
                              one p-value far below the float64 minimum
  cohort_saige.txt            SAIGE step 2: BETA and AF_Allele2 refer to Allele2
                              (= ALT), the opposite of METAL/BOLT's Allele1

Example:
    python skills/gwas-sumstats-harmonize/examples/make_demo_data.py
"""
from __future__ import annotations

import math
import random
from pathlib import Path

HERE = Path(__file__).resolve().parent
N_VARIANTS = 200
COMPLEMENT = {"A": "T", "T": "A", "C": "G", "G": "C"}


def _truth(rng: random.Random) -> list[dict]:
    rows = []
    for i in range(N_VARIANTS):
        chrom = 1 + i * 22 // N_VARIANTS
        bp = 1_000_000 + (i % 10) * 250_000 + rng.randint(0, 200_000)
        if i % 25 == 3:  # palindromic
            ref, alt = rng.choice([("A", "T"), ("T", "A"), ("C", "G"), ("G", "C")])
        else:
            ref, alt = rng.sample([a for a in "ACGT"], 2)
            while COMPLEMENT[ref] == alt:
                ref, alt = rng.sample([a for a in "ACGT"], 2)
        af = round(rng.uniform(0.02, 0.98), 4)
        se = round(rng.uniform(0.015, 0.06), 4)
        beta = round(rng.gauss(0, se) + (0.25 if i % 50 == 7 else 0.0), 4)
        rows.append({"CHR": chrom, "BP": bp, "REF": ref, "ALT": alt, "AF": af,
                     "BETA": beta, "SE": se, "rsid": f"synth_rs{i + 1}"})
    rows.sort(key=lambda r: (r["CHR"], r["BP"]))
    return rows


def _p(beta: float, se: float) -> float:
    z = abs(beta / se)
    return math.erfc(z / math.sqrt(2))


def _write(path: Path, header: list[str], rows: list[list]) -> None:
    lines = ["\t".join(header)] + ["\t".join(str(v) for v in r) for r in rows]
    path.write_text("\n".join(lines) + "\n", encoding="utf-8", newline="\n")


def main() -> None:
    rng = random.Random(20261001)
    truth = _truth(rng)

    _write(HERE / "reference.tsv", ["CHR", "BP", "REF", "ALT", "AF"],
           [[t["CHR"], t["BP"], t["REF"], t["ALT"], t["AF"]] for t in truth])

    # GWAS-Catalog SSF: swaps, strand flips, missing EAF, duplicate, bad p
    ssf = []
    for i, t in enumerate(truth):
        ea, nea, beta, eaf = t["ALT"], t["REF"], t["BETA"], t["AF"]
        palin = COMPLEMENT[t["REF"]] == t["ALT"]
        if i % 5 == 1:
            ea, nea, beta, eaf = nea, ea, -beta, round(1 - eaf, 4)
        # Reverse strand. For A/T and C/G this is indistinguishable from a swap
        # by alleles alone; only EAF vs the reference AF can tell them apart.
        if i % 17 == 2 or (palin and i % 2 == 1):
            ea, nea = COMPLEMENT[ea], COMPLEMENT[nea]
        p = f"{_p(beta, t['SE']):.4g}"
        eaf_s = "" if i % 20 == 4 else eaf
        ssf.append([t["CHR"], t["BP"], ea, nea, beta, t["SE"], eaf_s, p, t["rsid"], 10000])
    ssf[10][7] = "0"
    ssf[11][7] = "1.5"
    ssf.append(list(ssf[30]))  # exact duplicate row
    _write(HERE / "cohort_ssf.tsv",
           ["chromosome", "base_pair_location", "effect_allele", "other_allele", "beta",
            "standard_error", "effect_allele_frequency", "p_value", "variant_id", "n"], ssf)

    # PLINK2 glm.logistic: A1 may be REF or ALT; effect only as OR
    plink = []
    for i, t in enumerate(truth):
        a1_is_alt = i % 3 != 0
        a1 = t["ALT"] if a1_is_alt else t["REF"]
        beta = t["BETA"] if a1_is_alt else -t["BETA"]
        freq = t["AF"] if a1_is_alt else round(1 - t["AF"], 4)
        plink.append([t["CHR"], t["BP"], t["rsid"], t["REF"], t["ALT"], a1, freq, 8000,
                      f"{math.exp(beta):.5g}", t["SE"], f"{_p(beta, t['SE']):.4g}"])
    _write(HERE / "cohort_plink2.glm.logistic",
           ["#CHROM", "POS", "ID", "REF", "ALT", "A1", "A1_FREQ", "OBS_CT", "OR",
            "LOG(OR)_SE", "P"], plink)

    # METAL: no CHR/BP columns, lowercase alleles, MarkerName chr:pos
    metal = []
    for i, t in enumerate(truth):
        metal.append([f"{t['CHR']}:{t['BP']}", t["ALT"].lower(), t["REF"].lower(), t["AF"],
                      t["BETA"], t["SE"], f"{_p(t['BETA'], t['SE']):.4g}", "+-+"])
    _write(HERE / "cohort_metal.tbl",
           ["MarkerName", "Allele1", "Allele2", "Freq1", "Effect", "StdErr", "P-value",
            "Direction"], metal)

    # REGENIE: LOG10P only; one extreme signal below float64 range
    regenie = []
    for i, t in enumerate(truth):
        log10p = -math.log10(max(_p(t["BETA"], t["SE"]), 1e-300))
        regenie.append([t["CHR"], t["BP"], t["rsid"], t["REF"], t["ALT"], t["AF"], 12000,
                        t["BETA"], t["SE"], f"{log10p:.4f}"])
    regenie[0][7], regenie[0][8], regenie[0][9] = 0.8, 0.02, "349.4000"
    _write(HERE / "cohort_regenie.regenie",
           ["CHROM", "GENPOS", "ID", "ALLELE0", "ALLELE1", "A1FREQ", "N", "BETA", "SE",
            "LOG10P"], regenie)

    # SAIGE step 2: Allele1 = REF, Allele2 = ALT; BETA and AF_Allele2 are for Allele2
    saige = []
    for t in truth:
        saige.append([t["CHR"], t["BP"], t["rsid"], t["REF"], t["ALT"], round(2 * 9000 * t["AF"]), t["AF"],
                      0, t["BETA"], t["SE"], round(t["BETA"] / t["SE"], 4), round(1 / t["SE"] ** 2, 2),
                      f"{_p(t['BETA'], t['SE']):.4g}", 9000])
    _write(HERE / "cohort_saige.txt",
           ["CHR", "POS", "MarkerID", "Allele1", "Allele2", "AC_Allele2", "AF_Allele2", "MissingRate",
            "BETA", "SE", "Tstat", "var", "p.value", "N"], saige)


if __name__ == "__main__":
    main()
