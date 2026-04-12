"""One-off generator for tune-finemap task data.

Produces 3 synthetic GWAS loci with known causal variants.
Outputs: data/locus_{1,2,3}.tsv, data/locus_{1,2,3}_ld.npy, ground_truth.json
"""
from __future__ import annotations

import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd

_HERE = Path(__file__).resolve().parent
_SKILL = _HERE.parent.parent.parent / "fine-mapping"
sys.path.insert(0, str(_SKILL))
from fine_mapping import _make_block_ld  # noqa: E402


def simulate_locus(seed, n_variants, n_samples, causal_idx_to_beta,
                   block_size, rho, inter_rho, locus_id):
    rng = np.random.default_rng(seed)
    positions = np.linspace(109_000_000, 109_200_000, n_variants, dtype=int)
    rsids = [f"rs_{locus_id}_{i+1:04d}" for i in range(n_variants)]
    R = _make_block_ld(n_variants, block_size=block_size, rho=rho,
                       inter_rho=inter_rho, rng=rng)
    true_betas = np.zeros(n_variants)
    for idx, beta in causal_idx_to_beta.items():
        true_betas[idx] = beta
    signal = np.sqrt(n_samples) * (R @ true_betas)
    noise = rng.normal(0, 1, n_variants)
    z = signal + noise
    from math import erfc
    p = np.array([2 * 0.5 * erfc(abs(zi) / np.sqrt(2)) for zi in z])
    df = pd.DataFrame({
        "rsid": rsids,
        "chr": "1",
        "pos": positions,
        "z": z,
        "se": np.full(n_variants, 1.0 / np.sqrt(n_samples)),
        "n": n_samples,
        "p": p,
    })
    causal_rsids = [rsids[i] for i in causal_idx_to_beta]
    return df, R, causal_rsids


LOCI = [
    # (locus_id, seed, n, n_samples, causal map, block_size, rho, inter_rho)
    # Harder: lower n_samples + effect sizes so signals hover near boundary.
    # High LD around causals makes naive SuSiE spread PIP.
    ("L1", 11, 200, 3500, {50: 0.11, 150: 0.10}, 40, 0.92, 0.03),
    ("L2", 23, 200, 2500, {90: 0.12},             40, 0.90, 0.02),
    ("L3", 37, 200, 3000, {30: 0.10, 110: 0.09, 170: 0.09}, 30, 0.90, 0.02),
]


def main():
    data_dir = _HERE / "data"
    data_dir.mkdir(exist_ok=True)
    truth = {}
    for (lid, seed, n, ns, causal, bs, rho, irho) in LOCI:
        df, R, causal_rsids = simulate_locus(seed, n, ns, causal, bs, rho, irho, lid)
        df.to_csv(data_dir / f"locus_{lid}.tsv", sep="\t", index=False)
        np.save(data_dir / f"locus_{lid}_ld.npy", R)
        truth[lid] = sorted(causal_rsids)
        print(f"  {lid}: {len(df)} variants, causal={causal_rsids}")
    (_HERE / "ground_truth.json").write_text(json.dumps(truth, indent=2))
    print(f"Wrote ground_truth.json with {len(truth)} loci")


if __name__ == "__main__":
    main()
