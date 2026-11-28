#!/usr/bin/env python3
"""Write the synthetic demo summary statistics for coloc-abf-region-run.

Nothing here is real data. Four tables are generated over a 200 kb stretch of
chromosome 1, with association signals injected by hand so each demo has a known
right answer:

  synthetic_eqtl_shared.tsv  + synthetic_gwas_shared.tsv
      one causal variant, shared by both traits, so PP.H4 should dominate.

  synthetic_eqtl_distinct.tsv + synthetic_gwas_distinct.tsv
      a strong signal in each trait, at variants 60 kb apart, so PP.H3 should
      dominate.

The signal is a z-score peak decaying with distance from the causal variant, so
neighbouring variants carry correlated evidence the way a real locus does, plus a
little reproducible noise. There is no genotype data behind it, so these files
demonstrate the pipeline and the expected call, not the behaviour of any real
locus.

Deterministic: a fixed seed and a fixed layout, so re-running reproduces the
committed files byte for byte. `TestDemos.test_the_demo_data_is_reproducible`
checks exactly that, which is what makes the committed tables auditable rather
than something to be taken on trust.

They are written to `examples/sumstats/` rather than `examples/data/` because the
repository's .gitignore excludes `data/`, and committed demo data has to be
tracked. Naming the directory for what is in it avoids a negation rule.

    python skills/coloc-abf-region-run/examples/make_demo_data.py
"""

from __future__ import annotations

import argparse
import math
from pathlib import Path

import numpy as np

# A 200 kb window on chromosome 1, 1 kb apart, in the region of SORT1 on GRCh38.
# The coordinates are real; the statistics attached to them are not.
CHROM = "1"
START = 109_200_000
STEP = 1_000
N_VARIANTS = 200
DECAY_BP = 8_000.0  # how fast the injected peak falls away with distance

HEADER = ["variant_id", "chromosome", "position", "ref", "alt",
          "beta", "se", "p_value", "eaf", "n"]


def _positions() -> np.ndarray:
    return START + STEP * np.arange(N_VARIANTS)


def _alleles(rng: np.random.Generator) -> list[tuple[str, str]]:
    """One ref/alt pair per variant, drawn from the non-palindromic pairs.

    A/T and G/C pairs are left out so that nobody reading the demo has to wonder
    whether the strand matters here. It does not: the Bayes factor reads the
    effect only as z squared, so the posteriors do not change if an allele pair
    is flipped. Keeping them out of the demo keeps that a separate point.
    """
    choices = [("A", "G"), ("G", "A"), ("C", "T"), ("T", "C"),
               ("A", "C"), ("C", "A"), ("G", "T"), ("T", "G")]
    idx = rng.integers(0, len(choices), N_VARIANTS)
    return [choices[i] for i in idx]


def _trait(rng: np.random.Generator, causal_index: int, peak_z: float,
           se_value: float) -> tuple[np.ndarray, np.ndarray, np.ndarray]:
    """Beta, se and p for one trait with a peak at `causal_index`.

    `z = peak_z * exp(-distance / DECAY_BP) + N(0, 1)`: the decay stands in for
    correlation with the causal variant, and the unit-variance term for the
    sampling noise a real z-score carries.
    """
    pos = _positions()
    distance = np.abs(pos - pos[causal_index]).astype(float)
    z = peak_z * np.exp(-distance / DECAY_BP) + rng.normal(0.0, 1.0, N_VARIANTS)
    se = np.full(N_VARIANTS, se_value)
    beta = z * se
    p = np.array([2.0 * 0.5 * math.erfc(abs(v) / math.sqrt(2.0)) for v in z])
    return beta, se, p


def _write(path: Path, alleles, beta, se, p, eaf, n) -> None:
    pos = _positions()
    lines = ["\t".join(HEADER)]
    for i in range(N_VARIANTS):
        ref, alt = alleles[i]
        lines.append("\t".join([
            f"{CHROM}_{pos[i]}_{ref}_{alt}",
            CHROM,
            str(pos[i]),
            ref,
            alt,
            f"{beta[i]:.6f}",
            f"{se[i]:.6f}",
            f"{p[i]:.4e}",
            f"{eaf[i]:.4f}",
            str(n),
        ]))
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def write_demo_data(out_dir: Path) -> list[Path]:
    """Write the four demo tables and return their paths."""
    rng = np.random.default_rng(20261101)
    alleles = _alleles(rng)
    eaf = rng.uniform(0.08, 0.92, N_VARIANTS)

    shared_causal = 100            # position 109,300,000
    distinct_causal_gwas = 160     # 60 kb away, well outside the decay length

    # One shared causal variant: the expression and the lipid trait peak together.
    e_beta, e_se, e_p = _trait(rng, shared_causal, peak_z=10.0, se_value=0.025)
    g_beta, g_se, g_p = _trait(rng, shared_causal, peak_z=9.0, se_value=0.012)

    # Two distinct causal variants: each trait peaks somewhere the other does not.
    d_beta, d_se, d_p = _trait(rng, shared_causal, peak_z=10.0, se_value=0.025)
    h_beta, h_se, h_p = _trait(rng, distinct_causal_gwas, peak_z=9.0, se_value=0.012)

    paths = [
        out_dir / "synthetic_eqtl_shared.tsv",
        out_dir / "synthetic_gwas_shared.tsv",
        out_dir / "synthetic_eqtl_distinct.tsv",
        out_dir / "synthetic_gwas_distinct.tsv",
    ]
    _write(paths[0], alleles, e_beta, e_se, e_p, eaf, 500)
    _write(paths[1], alleles, g_beta, g_se, g_p, eaf, 250_000)
    _write(paths[2], alleles, d_beta, d_se, d_p, eaf, 500)
    _write(paths[3], alleles, h_beta, h_se, h_p, eaf, 250_000)
    return paths


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--out", type=Path,
                        default=Path(__file__).resolve().parent / "sumstats",
                        help="Directory to write the four tables into.")
    args = parser.parse_args(argv)
    for path in write_demo_data(args.out):
        print(f"wrote {path}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
