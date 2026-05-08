#!/usr/bin/env python3
"""DnaSP  -  Population Genetics Analysis of DNA Sequence Alignments.

Python reimplementation of core DnaSP statistics, faithful to the original
Visual Basic source (Rozas et al., J. Hered. 2017, doi:10.1093/jhered/esx062).

Statistical formulas follow:
  - Tajima (1989) Genetics 123:585-595  (Tajima's D)
  - Fu & Li (1993) Genetics 133:693-709  (D, F with outgroup; D*, F* without)
  - Simonsen et al. (1995) Genetics 141:413-429  (variance coefficients A3-A6, Appendix B)
  - Hudson, Kreitman & Aguadé (1987) Genetics 116:153-159  (HKA test)
  - Nei & Tajima (1981) Genetics 97:145-163  (haplotype diversity)
  - Ramos-Onsins & Rozas (2002) Mol. Biol. Evol. 19:2092-2100  (R2)
  - Kelly (1997) Genetics 146:1197-1206  (ZnS)
  - Rozas et al. (2001) Genetics 158:1147-1155  (Za, ZZ)
  - Hill & Robertson (1968) Theor. Appl. Genet. 38:226-231  (R²)
  - Lewontin (1964) Genetics 49:49-67  (D')
  - Hudson & Kaplan (1985) Genetics 111:147-164  (Rm)
  - Rogers & Harpending (1992) Mol. Biol. Evol. 9:552-569  (mismatch)
  - Harpending (1994) Hum. Biol. 66:591-600  (raggedness)
  - Nei (1987) Molecular Evolutionary Genetics  (Dxy, Da, divergence)

Supported input formats:
  FASTA  (.fas, .fa, .fasta)  -  including DnaSP header style >'name'  [comment]
  NEXUS  (.nex, .nexus, .nxs)  -  interleaved or non-interleaved, with MATCHCHAR

Available analyses (--analysis flag):
  polymorphism  : π, k, S, Eta, H, Hd, θ_W, Tajima's D, Fu & Li D*/F*, R2  [default]
  ld            : Linkage disequilibrium (D, D', R², ZnS, Za, ZZ)
  recombination : Minimum recombination events Rm (Hudson & Kaplan 1985)
  popsize       : Mismatch distribution, raggedness r (Harpending 1994)
  indel         : InDel polymorphism statistics
  divergence    : Dxy, Da, fixed/shared/private differences (needs --pop-file or --input2)
  fuliout       : Fu & Li D/F with outgroup (needs --outgroup)
  hka           : HKA neutrality test across loci (needs --hka-file)
  mk            : McDonald-Kreitman test (needs --outgroup; coding alignment)
  kaks          : Ka/Ks (dN/dS) via Nei-Gojobori 1986 (coding alignment)
  fufs          : Fu's Fs neutrality test (Fu 1997; Ewens sampling formula)
  sfs           : Site frequency spectrum  -  folded; unfolded if --outgroup provided
  tstv          : Transition/transversion ratio (Ts/Tv) across all pairs
  codon         : Codon usage bias  -  RSCU (Sharp & Li 1987) + ENC (Wright 1990)
  faywu         : Fay & Wu's H (2000) + Zeng's E (2006)  -  outgroup required
  fst           : Population differentiation Fst  -  Hudson et al. (1992) pairwise
  all           : Run all applicable analyses given the inputs

Usage:
  python dnasp.py --input alignment.fas --output results/
  python dnasp.py --input alignment.fas --analysis all --output results/
  python dnasp.py --input alignment.fas --analysis ld,recombination --output results/
  python dnasp.py --input alignment.fas --window 100 --step 25 --output results/
  python dnasp.py --input pop1.fas --input2 pop2.fas --analysis divergence --output results/
  python dnasp.py --input alignment.fas --pop-file pops.txt --analysis divergence --output results/
  python dnasp.py --demo --output /tmp/dnasp_demo
"""

from __future__ import annotations

__version__ = "0.4.0"
__author__  = "David De Lorenzo"
__credits__ = [
    # Python reimplementation and ClawBio adaptation
    "David De Lorenzo",
    # Original DnaSP 6 algorithms and Visual Basic implementation
    "Julio Rozas", "Albert Ferrer-Mata", "Juan Carlos Sánchez-DelBarrio",
    "Sara Guirao-Rico", "Pablo Librado", "Sebastián Ramos-Onsins",
    "Alejandro Sánchez-Gracia",
]
__licence__ = "MIT"

import argparse
import csv
import hashlib
import math
import re
import shutil
import sys
from collections import Counter
from dataclasses import dataclass, field
from datetime import datetime
from itertools import combinations
from pathlib import Path
from typing import Callable, Optional

# ── Optional heavy deps (graceful degradation) ───────────────────────────────
try:
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    HAS_MPL = True
except ImportError:
    HAS_MPL = False

# ─────────────────────────────────────────────────────────────────────────────
# Constants
# ─────────────────────────────────────────────────────────────────────────────

VALID_ANALYSES = {"polymorphism", "ld", "recombination", "popsize", "indel", "divergence",
                  "fuliout", "hka", "mk", "kaks", "fufs", "sfs", "tstv", "codon",
                  "faywu", "fst"}
_GAP_CHARS = frozenset("-?N")
_NUCLEOTIDES = ('A', 'T', 'C', 'G')
_PURINES     = frozenset('AG')
_PYRIMIDINES = frozenset('CT')

# Standard genetic code (NCBI transl_table=1)
GENETIC_CODE: dict[str, str] = {
    'TTT': 'F', 'TTC': 'F', 'TTA': 'L', 'TTG': 'L',
    'CTT': 'L', 'CTC': 'L', 'CTA': 'L', 'CTG': 'L',
    'ATT': 'I', 'ATC': 'I', 'ATA': 'I', 'ATG': 'M',
    'GTT': 'V', 'GTC': 'V', 'GTA': 'V', 'GTG': 'V',
    'TCT': 'S', 'TCC': 'S', 'TCA': 'S', 'TCG': 'S',
    'CCT': 'P', 'CCC': 'P', 'CCA': 'P', 'CCG': 'P',
    'ACT': 'T', 'ACC': 'T', 'ACA': 'T', 'ACG': 'T',
    'GCT': 'A', 'GCC': 'A', 'GCA': 'A', 'GCG': 'A',
    'TAT': 'Y', 'TAC': 'Y', 'TAA': '*', 'TAG': '*',
    'CAT': 'H', 'CAC': 'H', 'CAA': 'Q', 'CAG': 'Q',
    'AAT': 'N', 'AAC': 'N', 'AAA': 'K', 'AAG': 'K',
    'GAT': 'D', 'GAC': 'D', 'GAA': 'E', 'GAG': 'E',
    'TGT': 'C', 'TGC': 'C', 'TGA': '*', 'TGG': 'W',
    'CGT': 'R', 'CGC': 'R', 'CGA': 'R', 'CGG': 'R',
    'AGT': 'S', 'AGC': 'S', 'AGA': 'R', 'AGG': 'R',
    'GGT': 'G', 'GGC': 'G', 'GGA': 'G', 'GGG': 'G',
}

# Synonymous codon families (Group E: codon usage bias)
# Maps amino acid → list of codons; stop codons ('*') excluded.
_SYNONYMOUS_FAMILIES: dict[str, list[str]] = {}
for _codon, _aa in GENETIC_CODE.items():
    if _aa != '*':
        _SYNONYMOUS_FAMILIES.setdefault(_aa, []).append(_codon)
del _codon, _aa  # clean up loop variables

# ─────────────────────────────────────────────────────────────────────────────
# Data structures
# ─────────────────────────────────────────────────────────────────────────────

@dataclass
class Alignment:
    names: list[str]
    seqs: list[str]
    source: str = ""

    def __post_init__(self) -> None:
        lengths = {len(s) for s in self.seqs}
        if len(lengths) > 1:
            raise ValueError(
                f"Sequences have unequal lengths: {sorted(lengths)}. "
                "Is the alignment correctly formatted?"
            )

    @property
    def n(self) -> int:
        return len(self.seqs)

    @property
    def L(self) -> int:
        return len(self.seqs[0]) if self.seqs else 0


@dataclass
class RegionStats:
    """All polymorphism statistics for one region (whole alignment or sliding window)."""
    region: str = "1-?"
    n: int = 0
    L_total: int = 0
    L_net: int = 0
    S: int = 0
    Eta: int = 0
    H: int = 0
    Hd: float = 0.0
    VarHd: float = 0.0
    Pi: float = 0.0
    k: float = 0.0
    GC: float = 0.0
    ThetaW: float = 0.0
    ThetaW_nuc: float = 0.0
    TajimaD: Optional[float] = None
    FuLiD_star: Optional[float] = None
    FuLiF_star: Optional[float] = None
    R2: Optional[float] = None

    def as_tsv_row(self) -> list:
        def fmt(v: Optional[float]) -> str:
            return "n.a." if v is None else f"{v:.10g}"
        return [
            self.region, self.n, self.L_total, self.L_net,
            self.S, self.Eta, self.H,
            fmt(self.Hd), fmt(self.VarHd),
            fmt(self.Pi), fmt(self.k),
            fmt(self.GC),
            fmt(self.ThetaW_nuc), fmt(self.ThetaW),
            fmt(self.TajimaD),
            fmt(self.FuLiD_star), fmt(self.FuLiF_star),
            fmt(self.R2),
        ]


TSV_HEADER = [
    "Region", "Sample_Size", "Sites", "NetSites",
    "S", "Eta", "Hap",
    "Hd", "VarHd",
    "Pi", "ThetaK",
    "G+Ctot",
    "ThetaWattNuc", "ThetaWatt",
    "TajimaD",
    "FuLiD*", "FuLiF*",
    "Ramos-Onsins_Rozas_R2",
]


@dataclass
class LDPair:
    """LD statistics for one pair of biallelic sites."""
    site1: int          # 1-based position in net alignment
    site2: int
    dist: int           # nucleotide distance
    n_valid: int        # sequences without gaps at either site
    D: float
    D_prime: Optional[float]
    R2: Optional[float]
    chi2: Optional[float]
    p_chi2: Optional[float]


@dataclass
class LDStats:
    """Whole-alignment LD summary."""
    n_biallelic: int = 0
    n_pairs: int = 0
    ZnS: Optional[float] = None
    Za: Optional[float] = None
    ZZ: Optional[float] = None
    pairs: list[LDPair] = field(default_factory=list)


@dataclass
class RecombStats:
    """Recombination statistics."""
    Rm: int = 0
    n_incompatible_pairs: int = 0
    incompatible_pairs: list[tuple[int, int]] = field(default_factory=list)


@dataclass
class MismatchStats:
    """Mismatch distribution statistics (Population Size Changes module)."""
    observed: dict[int, int] = field(default_factory=dict)
    n_pairs: int = 0
    mean: float = 0.0
    variance: float = 0.0
    raggedness: Optional[float] = None
    cv: Optional[float] = None


@dataclass
class InDelEvent:
    """One identified InDel event."""
    start: int          # 0-based alignment position
    end: int            # inclusive
    length: int         # bp
    seq_indices: frozenset = field(default_factory=frozenset)  # indices carrying gap


@dataclass
class InDelStats:
    """InDel polymorphism statistics."""
    n_positions_with_gaps: int = 0
    n_events: int = 0
    mean_event_length: float = 0.0
    n_haplotypes: int = 0
    haplotype_diversity: float = 0.0
    k_indel: float = 0.0
    pi_indel: float = 0.0
    theta_indel: float = 0.0
    tajima_d_indel: Optional[float] = None
    events: list[InDelEvent] = field(default_factory=list)


@dataclass
class DivergenceStats:
    """Between-population divergence statistics (DNA Divergence module)."""
    pop1_name: str = "Pop1"
    pop2_name: str = "Pop2"
    n1: int = 0
    n2: int = 0
    L_net: int = 0
    # Within-population
    Pi1: float = 0.0
    Pi2: float = 0.0
    k1: float = 0.0
    k2: float = 0.0
    # Between populations
    Dxy: float = 0.0   # average pairwise differences between populations per site
    Da: float = 0.0    # net nucleotide differences (Da = Dxy - (Pi1+Pi2)/2)
    n_fixed: int = 0           # fixed differences
    n_shared: int = 0          # shared polymorphisms
    n_private1: int = 0        # private to pop1
    n_private2: int = 0        # private to pop2


@dataclass
class FuLiOutgroupStats:
    """Fu & Li (1993) D and F statistics polarised by an outgroup sequence."""
    n: int = 0           # number of ingroup sequences
    eta: int = 0         # total derived mutations (outgroup-polarised)
    eta_e: int = 0       # derived mutations carried by exactly 1 ingroup seq
    k_bar: float = 0.0   # mean pairwise differences (standard π estimate)
    D: Optional[float] = None  # Fu & Li D (outgroup version)
    F: Optional[float] = None  # Fu & Li F (outgroup version)


@dataclass
class HKALocus:
    """One locus for the HKA test."""
    name: str = ""   # locus identifier
    n: int = 0       # ingroup sample size
    S: int = 0       # segregating sites within ingroup
    D: int = 0       # fixed differences (divergence) to outgroup/sister species


@dataclass
class HKAStats:
    """Results of the HKA test (Hudson, Kreitman & Aguadé 1987)."""
    n_loci: int = 0
    T_hat: float = 0.0          # MLE divergence time (units of N_e generations)
    chi2: float = 0.0
    df: int = 0
    p_value: Optional[float] = None
    loci_results: list = field(default_factory=list)  # per-locus details (list of dicts)


@dataclass
class MKStats:
    """McDonald-Kreitman test (McDonald & Kreitman 1991)."""
    Pn: int = 0               # nonsynonymous polymorphisms in ingroup
    Ps: int = 0               # synonymous polymorphisms in ingroup
    Dn: int = 0               # nonsynonymous fixed differences (ingroup vs outgroup)
    Ds: int = 0               # synonymous fixed differences (ingroup vs outgroup)
    alpha: Optional[float] = None   # 1 − (Ds·Pn)/(Dn·Ps)  -  proportion of adaptive substitutions
    NI: Optional[float] = None      # Neutrality Index = (Pn/Ps)/(Dn/Ds)
    DoS: Optional[float] = None     # Direction of Selection = Dn/(Dn+Ds) − Pn/(Pn+Ps)
    fisher_p: Optional[float] = None  # Two-tailed Fisher's exact test P-value


@dataclass
class KaKsStats:
    """Ka/Ks (dN/dS) estimated by the Nei-Gojobori (1986) method."""
    n_codons: int = 0
    S_sites: float = 0.0     # mean synonymous sites per sequence
    N_sites: float = 0.0     # mean nonsynonymous sites per sequence
    Sd: float = 0.0          # mean synonymous differences per sequence pair
    Nd: float = 0.0          # mean nonsynonymous differences per sequence pair
    Ks: Optional[float] = None    # synonymous substitutions per syn site (JC corrected)
    Ka: Optional[float] = None    # nonsynonymous substitutions per nonsyn site (JC corrected)
    omega: Optional[float] = None  # Ka/Ks; < 1 purifying, ≈ 1 neutral, > 1 positive selection


@dataclass
class FuFsStats:
    """Fu's Fs neutrality test (Fu 1997).

    Fs = ln(S_k / (1 - S_k)) where S_k = P(K_n ≤ H | θ_π, n) under the
    Ewens sampling formula (infinite-alleles model).
    Significant at the 0.02 level (conventional threshold for Fs).
    """
    n: int = 0
    H: int = 0                        # observed number of haplotypes
    theta_pi: float = 0.0             # θ_π = k (mean pairwise differences)
    S_k: Optional[float] = None       # P(K_n ≤ H | θ_π, n)  -  Ewens CDF
    Fs: Optional[float] = None        # ln(S_k / (1 - S_k)); Fs << 0 → expansion/selection


@dataclass
class SFSStats:
    """Site frequency spectrum (folded and optionally unfolded).

    folded[i]   = number of segregating sites where minor allele count = i
                  (i = 1 … n//2); independent of outgroup.
    unfolded[i] = number of segregating sites where derived allele count = i
                  (i = 1 … n-1); requires an outgroup to polarise.
    """
    n: int = 0
    folded: dict[int, int] = field(default_factory=dict)
    unfolded: Optional[dict[int, int]] = None
    has_outgroup: bool = False


@dataclass
class TsTvStats:
    """Transition / transversion ratio across all pairwise comparisons.

    Computed over clean (non-gap, unambiguous) columns only (complete deletion).
    n_transitions and n_transversions are *total* counts summed over all pairs.
    ts_per_site and tv_per_site are means per sequence pair per site.
    ts_tv is None when n_transversions == 0 (all differences are transitions).
    """
    n: int = 0          # number of ingroup sequences
    L_net: int = 0      # number of clean columns used
    n_transitions: int = 0      # total Ts across all pairs
    n_transversions: int = 0    # total Tv across all pairs
    ts_tv: Optional[float] = None       # Ts/Tv ratio; None if Tv == 0
    ts_per_site: float = 0.0            # mean Ts per pair per site
    tv_per_site: float = 0.0            # mean Tv per pair per site


@dataclass
class CodonUsageStats:
    """Codon usage bias: RSCU (Sharp & Li 1987) and ENC (Wright 1990).

    rscu[codon] = observed / expected usage within synonymous family.
      RSCU = 1.0  → uniform usage among synonyms.
      RSCU > 1.0  → preferred codon.
      RSCU < 1.0  → avoided codon.

    ENC (Effective Number of Codons):
      20  → maximum bias (only one codon used per amino acid).
      61  → no bias (all synonymous codons used equally).
      ENC < 35 suggests strong codon usage bias.

    n_codons is the mean number of clean (non-stop, non-gap) codons per sequence.
    """
    n: int = 0                  # number of ingroup sequences
    n_codons: float = 0.0       # mean clean coding codons per sequence
    codon_counts: dict[str, float] = field(default_factory=dict)  # codon → mean count/seq
    rscu: dict[str, float] = field(default_factory=dict)          # codon → RSCU value
    ENC: Optional[float] = None     # 20 (max bias) … 61 (no bias); None if data insufficient


@dataclass
class FayWuStats:
    """Fay & Wu's H and Zeng's E neutrality statistics.

    Both require an outgroup to polarise each segregating site (determine
    which allele is ancestral vs derived).

    theta_pi  = standard nucleotide diversity (π), computed from polarised sites.
    theta_w   = Watterson's θ, computed from polarised sites.
    theta_h   = Fay & Wu's θ_H (weights derived alleles by frequency²).
    theta_l   = Zeng's θ_L (weights derived alleles by frequency).

    H  = θ_π − θ_H  (Fay & Wu 2000; negative → excess high-freq derived alleles).
    E  = θ_L − θ_W  (Zeng et al. 2006; complements H for different distortions).

    All θ values are per-site (divided by L_net).
    None when n_polarised == 0.
    """
    n: int = 0                     # ingroup sample size
    L_net: int = 0                 # sites surviving complete deletion with outgroup
    n_polarised: int = 0           # # polarisable segregating sites (ancestral ∈ ingroup)
    theta_pi: Optional[float] = None
    theta_w: Optional[float] = None
    theta_h: Optional[float] = None
    theta_l: Optional[float] = None
    H: Optional[float] = None      # θ_π − θ_H
    E: Optional[float] = None      # θ_L − θ_W


@dataclass
class FstStats:
    """Population differentiation (Fst) using Hudson et al. (1992) estimator.

    For each pair of populations A and B:
      Fst = 1 − π_s / π_t
    where π_s = (π_A + π_B) / 2  (mean within-population π, per site)
    and   π_t = π_AB              (between-population Dxy, per site).

    fst_pairwise maps (pop1, pop2) → Fst value (None if π_t == 0).
    fst_mean is the unweighted mean across all pairs.
    pi_within[pop] = within-pop nucleotide diversity per site.
    pi_between[(pop1, pop2)] = Dxy per site.
    """
    n_pops: int = 0
    pop_names: list[str] = field(default_factory=list)
    pop_sizes: dict[str, int] = field(default_factory=dict)
    pi_within: dict[str, float] = field(default_factory=dict)
    pi_between: dict[tuple, float] = field(default_factory=dict)
    fst_pairwise: dict[tuple, Optional[float]] = field(default_factory=dict)
    fst_mean: Optional[float] = None


# ─────────────────────────────────────────────────────────────────────────────
# Parsers
# ─────────────────────────────────────────────────────────────────────────────

_DNASP_HEADER = re.compile(r"^>?\s*'?([^'\[]+?)'?\s*(?:\[.*?\])?\s*$")


def parse_fasta(path: Path) -> Alignment:
    """Parse FASTA, including DnaSP's >'name'  [comment] format."""
    names: list[str] = []
    seqs: list[str] = []
    current_seq: list[str] = []

    with open(path, encoding="utf-8", errors="replace") as fh:
        for raw_line in fh:
            line = raw_line.rstrip("\n")
            if line.startswith(">"):
                if names:
                    seqs.append("".join(current_seq).upper())
                    current_seq = []
                m = _DNASP_HEADER.match(line)
                names.append(m.group(1).strip() if m else line[1:].strip())
            else:
                seq_part = line.strip()
                if seq_part and names:
                    current_seq.append(seq_part)

    if current_seq:
        seqs.append("".join(current_seq).upper())

    if not names:
        raise ValueError(f"No sequences found in {path}")
    if len(names) != len(seqs):
        raise ValueError(
            f"Parsed {len(names)} headers but {len(seqs)} sequences in {path}"
        )
    return Alignment(names=names, seqs=seqs, source=str(path))


def parse_nexus(path: Path) -> Alignment:
    """Parse NEXUS (DnaSP style), interleaved or sequential, with MATCHCHAR."""
    text = path.read_text(encoding="utf-8", errors="replace")

    matrix_m = re.search(r"MATRIX\s*(.*?)\s*;", text, re.IGNORECASE | re.DOTALL)
    if not matrix_m:
        raise ValueError(f"No MATRIX block found in {path}")

    matchchar = "."
    mc_m = re.search(r"MATCHCHAR\s*=\s*(\S)", text, re.IGNORECASE)
    if mc_m:
        matchchar = mc_m.group(1)

    matrix_text = matrix_m.group(1)
    seq_dict: dict[str, list[str]] = {}
    order: list[str] = []

    for line in matrix_text.splitlines():
        stripped = line.strip()
        if not stripped or stripped.startswith("["):
            continue
        m = re.match(r"^'?([^'\s]+)'?\s+(\S+)\s*$", stripped)
        if m:
            name, data = m.group(1), m.group(2)
            if name not in seq_dict:
                seq_dict[name] = []
                order.append(name)
            seq_dict[name].append(data.upper())

    if not order:
        raise ValueError(f"Could not parse sequences from MATRIX in {path}")

    raw_seqs = ["".join(seq_dict[n]) for n in order]

    if matchchar and any(matchchar in s for s in raw_seqs[1:]):
        ref = raw_seqs[0]
        expanded = [ref]
        for s in raw_seqs[1:]:
            exp = "".join(ref[i] if c == matchchar else c for i, c in enumerate(s))
            expanded.append(exp)
        raw_seqs = expanded

    return Alignment(names=order, seqs=raw_seqs, source=str(path))


def load_alignment(path: Path) -> Alignment:
    """Auto-detect format from extension."""
    suffix = path.suffix.lower()
    if suffix in {".fas", ".fa", ".fasta"}:
        return parse_fasta(path)
    if suffix in {".nex", ".nexus", ".nxs"}:
        return parse_nexus(path)
    try:
        return parse_fasta(path)
    except Exception:
        return parse_nexus(path)


def load_pop_file(pop_file: Path) -> dict[str, str]:
    """Read a population assignment file (tab-separated: seq_name<TAB>pop_name)."""
    assignments: dict[str, str] = {}
    with open(pop_file, encoding="utf-8") as fh:
        for line in fh:
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) >= 2:
                assignments[parts[0].strip()] = parts[1].strip()
    return assignments


# ─────────────────────────────────────────────────────────────────────────────
# Gap treatment
# ─────────────────────────────────────────────────────────────────────────────

def complete_deletion(seqs: list[str]) -> tuple[list[str], int]:
    """Exclude any column that has a gap/missing in ANY sequence."""
    if not seqs:
        return [], 0
    L = len(seqs[0])
    keep = [pos for pos in range(L) if not any(s[pos] in _GAP_CHARS for s in seqs)]
    cleaned = ["".join(s[i] for i in keep) for s in seqs]
    return cleaned, len(keep)


# ─────────────────────────────────────────────────────────────────────────────
# Core polymorphism functions
# ─────────────────────────────────────────────────────────────────────────────

def _harmonic(n: int, power: int = 1) -> float:
    """Σ 1/i^power for i = 1 .. n-1."""
    return sum(1.0 / i**power for i in range(1, n))


def compute_gc(seqs: list[str]) -> float:
    total = gc = 0
    for s in seqs:
        for c in s:
            if c in "ACGT":
                total += 1
                if c in "GC":
                    gc += 1
    return gc / total if total else 0.0


def compute_haplotypes(seqs: list[str]) -> tuple[int, float, float]:
    """H, Hd, VarHd (Nei & Tajima 1981)."""
    n = len(seqs)
    if n < 2:
        return 1, 0.0, 0.0
    hap_counts: dict[str, int] = {}
    for s in seqs:
        hap_counts[s] = hap_counts.get(s, 0) + 1
    H = len(hap_counts)
    freqs = [c / n for c in hap_counts.values()]
    sum_sq = sum(p * p for p in freqs)
    Hd = (n / (n - 1)) * (1.0 - sum_sq)
    sum_cu = sum(p * p * p for p in freqs)
    VarHd = (2 / (n * (n - 1))) * (
        2 * (n - 2) * (sum_cu - sum_sq**2) + sum_sq - sum_sq**2
    )
    return H, max(0.0, Hd), max(0.0, VarHd)


def compute_k(seqs: list[str]) -> float:
    """Mean pairwise differences (absolute count, not per site)."""
    n = len(seqs)
    if n < 2 or not seqs[0]:
        return 0.0
    total = sum(
        sum(a != b for a, b in zip(seqs[i], seqs[j]))
        for i in range(n) for j in range(i + 1, n)
    )
    return total / (n * (n - 1) / 2)


def _site_states(col: list[str]) -> frozenset[str]:
    return frozenset(c for c in col if c not in _GAP_CHARS)


def compute_segregating(seqs: list[str]) -> tuple[int, int]:
    """S (segregating sites) and Eta (total mutations)."""
    if not seqs:
        return 0, 0
    L = len(seqs[0])
    S = Eta = 0
    for pos in range(L):
        col = [s[pos] for s in seqs]
        states = _site_states(col)
        if len(states) > 1:
            S += 1
            Eta += len(states) - 1
    return S, Eta


def compute_singletons(seqs: list[str]) -> tuple[int, list[float]]:
    """η_s (singleton count) and per-sequence attribution.

    Matches DnaSP BusqSingletones(): counts alleles with count == 1 only.
    """
    n = len(seqs)
    if n < 2 or not seqs:
        return 0, [0.0] * n
    L = len(seqs[0])
    eta_s = 0
    per_seq = [0.0] * n
    n_triallelic = n_quadriallelic = 0

    for pos in range(L):
        col = [s[pos] for s in seqs]
        states = _site_states(col)
        if len(states) <= 1:
            continue
        if len(states) == 3:
            n_triallelic += 1
        elif len(states) == 4:
            n_quadriallelic += 1
        counts: dict[str, list[int]] = {st: [] for st in states}
        for idx, c in enumerate(col):
            if c in states:
                counts[c].append(idx)
        for allele, carriers in counts.items():
            if len(carriers) == 1:
                eta_s += 1
                per_seq[carriers[0]] += 1

    if n == 2:
        total = sum(per_seq)
        per_seq = [total / 2] * 2
    elif n == 3 and n_triallelic > 0:
        corr = n_triallelic / 3
        per_seq = [x - corr for x in per_seq]
    elif n == 4 and n_quadriallelic > 0:
        corr = n_quadriallelic / 4
        per_seq = [x - corr for x in per_seq]

    return eta_s, per_seq


def tajima_d(k: float, S: int, n: int) -> Optional[float]:
    """Tajima's D (Tajima 1989, Genetics 123:585-595)."""
    if S == 0 or n < 3:
        return None
    a1 = _harmonic(n, 1)
    a2 = _harmonic(n, 2)
    b1 = (n + 1) / (3 * (n - 1))
    b2 = (2 * (n**2 + n + 3)) / (9 * n * (n - 1))
    c1 = b1 - 1.0 / a1
    c2 = b2 - (n + 2) / (a1 * n) + a2 / a1**2
    e1 = c1 / a1
    e2 = c2 / (a1**2 + a2)
    vD = e1 * S + e2 * S * (S - 1)
    if vD <= 0:
        return None
    return (k - S / a1) / math.sqrt(vD)


def fu_li_d_star_f_star(
    k: float, S: int, eta_s: int, n: int
) -> tuple[Optional[float], Optional[float]]:
    """Fu & Li D* and F* (Simonsen et al. 1995, equations A3-A6)."""
    if S == 0 or n < 3:
        return None, None
    An = _harmonic(n, 1)
    Bn = _harmonic(n, 2)
    vD = (Bn / An**2 - (2 / n) * (1 + 1 / An - An + An / n) - 1 / n**2)
    vD /= An**2 + Bn
    vD = max(vD, 0.0)
    uD = ((n - 1) / n - 1 / An) / An - vD
    vF = (2 * n**3 + 110 * n**2 - 255 * n + 153) / (9 * n**2 * (n - 1))
    vF += (2 * (n - 1) * An) / n**2
    vF -= (8 * Bn) / n
    vF /= An**2 + Bn
    uF = (4 * n**2 + 19 * n + 3 - 12 * (n + 1) * (An + 1 / n))
    uF /= 3 * n * (n - 1)
    uF = uF / An - vF
    denom_D = uD * S + vD * S**2
    denom_F = uF * S + vF * S**2
    if denom_D <= 0 or denom_F <= 0:
        return None, None
    D_star = (S / An - eta_s * (n - 1) / n) / math.sqrt(denom_D)
    F_star = (k - eta_s * (n - 1) / n) / math.sqrt(denom_F)
    return D_star, F_star


def ramos_onsins_r2(seqs: list[str], k: float, Sw: int) -> Optional[float]:
    """R2 (Ramos-Onsins & Rozas 2002)."""
    if Sw == 0:
        return None
    n = len(seqs)
    _, per_seq = compute_singletons(seqs)
    total = sum((u - k / 2) ** 2 for u in per_seq)
    return math.sqrt(total / n) / Sw


def watterson_theta(S: int, n: int, L_net: int) -> tuple[float, float]:
    a1 = _harmonic(n, 1)
    if a1 == 0:
        return 0.0, 0.0
    theta_abs = S / a1
    theta_nuc = theta_abs / L_net if L_net > 0 else 0.0
    return theta_abs, theta_nuc


# ─────────────────────────────────────────────────────────────────────────────
# Linkage Disequilibrium
# ─────────────────────────────────────────────────────────────────────────────

def _chi2_1df_pvalue(chi2: float) -> float:
    """Two-tailed p-value for chi-square statistic with 1 degree of freedom."""
    if chi2 <= 0:
        return 1.0
    return math.erfc(math.sqrt(chi2 / 2))


def _bisect(f: Callable, a: float, b: float, tol: float = 1e-10, maxiter: int = 300) -> float:
    """Bisection root-finder: return x in [a,b] with f(x) ≈ 0."""
    fa = f(a)
    fb = f(b)
    if fa == 0.0:
        return a
    if fb == 0.0:
        return b
    if fa * fb > 0:
        # f does not change sign  -  return midpoint as fallback
        return (a + b) / 2.0
    for _ in range(maxiter):
        mid = (a + b) / 2.0
        fm = f(mid)
        if abs(fm) < tol or (b - a) / 2.0 < tol:
            return mid
        if fa * fm <= 0:
            b = mid
            fb = fm
        else:
            a = mid
            fa = fm
    return (a + b) / 2.0


def _gammaincl_series(a: float, x: float, max_iter: int = 300) -> float:
    """Lower regularised incomplete gamma P(a, x) via series expansion.

    Converges best when x < a + 1.  Uses the recurrence:
        P(a, x) = e^{-x} x^a / Γ(a+1) * Σ_{n=0}^∞ x^n / [(a+1)(a+2)…(a+n)]
    """
    if x == 0.0:
        return 0.0
    lnpre = a * math.log(x) - x - math.lgamma(a + 1)
    term = 1.0
    total = 1.0
    ap = a
    for _ in range(max_iter):
        ap += 1.0
        term *= x / ap
        total += term
        if term < total * 3e-15:
            break
    return total * math.exp(lnpre)


def _gammaincl_cf(a: float, x: float, max_iter: int = 300) -> float:
    """Upper regularised incomplete gamma Q(a, x) via Lentz continued fraction.

    Converges best when x >= a + 1 (Numerical Recipes §6.2).
    """
    fpmin = 1e-300
    b = x + 1.0 - a
    c = 1.0 / fpmin
    d = 1.0 / b if abs(b) > fpmin else 1.0 / fpmin
    h = d
    for i in range(1, max_iter + 1):
        an = -i * (i - a)
        b += 2.0
        d = an * d + b
        if abs(d) < fpmin:
            d = fpmin
        c = b + an / c
        if abs(c) < fpmin:
            c = fpmin
        d = 1.0 / d
        delta = d * c
        h *= delta
        if abs(delta - 1.0) < 1e-14:
            break
    return math.exp(-x + a * math.log(x) - math.lgamma(a)) * h


def _gammaincc(a: float, x: float) -> float:
    """Regularised upper incomplete gamma Q(a, x) = 1 − P(a, x).

    Uses series expansion for small x and continued fraction for large x.
    """
    if x <= 0.0:
        return 1.0
    if x < a + 1.0:
        return 1.0 - _gammaincl_series(a, x)
    return _gammaincl_cf(a, x)


def _chi2_pvalue(chi2: float, df: int) -> float:
    """Survival function P(X > chi2) for X ~ χ²(df)."""
    if chi2 <= 0.0:
        return 1.0
    if df <= 0:
        return 1.0
    if df == 1:
        return math.erfc(math.sqrt(chi2 / 2.0))
    if df == 2:
        return math.exp(-chi2 / 2.0)
    return _gammaincc(df / 2.0, chi2 / 2.0)


def _get_biallelic_positions(seqs: list[str]) -> list[tuple[int, str, str]]:
    """Return list of (position, minor_allele, major_allele) for biallelic sites.

    Only strictly biallelic sites (exactly 2 nucleotide states, no gaps) are
    returned.  Sites with gaps in any sequence are excluded (complete deletion).
    Multiallelic sites (3+ states) are excluded.
    """
    if not seqs:
        return []
    L = len(seqs[0])
    result = []
    for pos in range(L):
        col = [s[pos] for s in seqs]
        if any(c in _GAP_CHARS for c in col):
            continue
        counts: Counter = Counter(col)
        if len(counts) != 2:
            continue
        alleles = sorted(counts.keys(), key=lambda a: counts[a])
        minor, major = alleles[0], alleles[1]
        result.append((pos, minor, major))
    return result


def _ld_for_pair(
    seqs: list[str], pos_i: int, minor_i: str, pos_j: int, minor_j: str
) -> tuple[float, Optional[float], Optional[float], int]:
    """Compute D, D', R² for a pair of biallelic sites.

    Returns (D, D_prime, R2, n_valid).
    Minor alleles at each site define the "A" and "B" haplotype.
    """
    n11 = n12 = n21 = n22 = 0
    valid = 0
    for s in seqs:
        ci, cj = s[pos_i], s[pos_j]
        if ci in _GAP_CHARS or cj in _GAP_CHARS:
            continue
        valid += 1
        a = (ci == minor_i)
        b = (cj == minor_j)
        if a and b:
            n11 += 1
        elif a and not b:
            n12 += 1
        elif not a and b:
            n21 += 1
        else:
            n22 += 1

    if valid == 0:
        return 0.0, None, None, 0

    p_A = (n11 + n12) / valid   # freq of minor allele at site i
    p_a = 1.0 - p_A
    p_B = (n11 + n21) / valid   # freq of minor allele at site j
    p_b = 1.0 - p_B
    p_AB = n11 / valid

    D = p_AB - p_A * p_B

    if D > 0:
        D_max = min(p_A * p_b, p_a * p_B)
    elif D < 0:
        D_max = min(p_A * p_B, p_a * p_b)
    else:
        D_max = 1.0

    D_prime = (D / D_max) if D_max > 1e-12 else None

    denom_r2 = p_A * p_a * p_B * p_b
    R2 = (D**2 / denom_r2) if denom_r2 > 1e-12 else None

    return D, D_prime, R2, valid


def compute_ld(seqs: list[str], positions: Optional[list[int]] = None) -> LDStats:
    """Full LD analysis: D, D', R², ZnS, Za, ZZ.

    References:
      Lewontin & Kojima 1960 (D), Lewontin 1964 (D'),
      Hill & Robertson 1968 (R²), Kelly 1997 (ZnS), Rozas et al. 2001 (Za, ZZ).

    Args:
        seqs: Sequences after complete deletion (no gap columns).
        positions: Original alignment positions (1-based) for each site in seqs,
                   used to compute nucleotide distance. If None, uses index+1.
    """
    biallelic = _get_biallelic_positions(seqs)
    stats = LDStats(n_biallelic=len(biallelic))

    if len(biallelic) < 2:
        return stats

    if positions is None:
        pos_map = {i: i + 1 for i in range(len(seqs[0]))}
    else:
        pos_map = {i: positions[i] for i in range(len(positions))}

    pairs: list[LDPair] = []
    for (idx_i, (pi, mi, _)), (idx_j, (pj, mj, _)) in combinations(
        enumerate(biallelic), 2
    ):
        D, D_prime, R2, n_valid = _ld_for_pair(seqs, pi, mi, pj, mj)
        dist = abs(pos_map.get(pj, pj + 1) - pos_map.get(pi, pi + 1))
        chi2 = n_valid * R2 if R2 is not None else None
        p_chi2 = _chi2_1df_pvalue(chi2) if chi2 is not None else None
        pairs.append(LDPair(
            site1=pos_map.get(pi, pi + 1),
            site2=pos_map.get(pj, pj + 1),
            dist=dist,
            n_valid=n_valid,
            D=D,
            D_prime=D_prime,
            R2=R2,
            chi2=chi2,
            p_chi2=p_chi2,
        ))

    stats.pairs = pairs
    stats.n_pairs = len(pairs)

    # ZnS: average R² over all pairs (Kelly 1997, eq 3)
    r2_vals = [p.R2 for p in pairs if p.R2 is not None]
    stats.ZnS = sum(r2_vals) / len(r2_vals) if r2_vals else None

    # Za: average R² for adjacent (consecutive) pairs of biallelic sites (Rozas 2001)
    adjacent_r2 = []
    for k in range(len(biallelic) - 1):
        # The pair (biallelic[k], biallelic[k+1])
        pi, mi, _ = biallelic[k]
        pj, mj, _ = biallelic[k + 1]
        _, _, R2, _ = _ld_for_pair(seqs, pi, mi, pj, mj)
        if R2 is not None:
            adjacent_r2.append(R2)

    stats.Za = sum(adjacent_r2) / len(adjacent_r2) if adjacent_r2 else None

    # ZZ = Za - ZnS
    if stats.Za is not None and stats.ZnS is not None:
        stats.ZZ = stats.Za - stats.ZnS

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# Recombination (Rm, four-gamete test)
# ─────────────────────────────────────────────────────────────────────────────

def compute_recombination(seqs: list[str], net_positions: Optional[list[int]] = None) -> RecombStats:
    """Minimum recombination events Rm (Hudson & Kaplan 1985).

    The four-gamete test identifies pairs of biallelic sites that are
    incompatible (all four haplotype combinations observed).  Each incompatible
    pair (i, j) requires at least one recombination event between positions i
    and j.  Rm is the minimum number of recombination events needed to account
    for all incompatible pairs, computed by the interval-stabbing algorithm
    (sort intervals by right endpoint; greedily place events).

    Reference: Hudson RR, Kaplan NL (1985) Genetics 111:147-164.
    """
    biallelic = _get_biallelic_positions(seqs)
    stats = RecombStats()

    if len(biallelic) < 2:
        return stats

    if net_positions is None:
        pos_map = {i: i + 1 for i in range(len(seqs[0]))}
    else:
        pos_map = {i: net_positions[i] for i in range(len(net_positions))}

    incompatible: list[tuple[int, int]] = []
    for (idx_i, (pi, mi, _)), (idx_j, (pj, mj, _)) in combinations(
        enumerate(biallelic), 2
    ):
        gametes: set[tuple[bool, bool]] = set()
        for s in seqs:
            ci, cj = s[pi], s[pj]
            if ci not in _GAP_CHARS and cj not in _GAP_CHARS:
                gametes.add((ci == mi, cj == mj))
        if len(gametes) == 4:  # all four gamete types present → incompatible
            p1 = pos_map.get(pi, pi + 1)
            p2 = pos_map.get(pj, pj + 1)
            incompatible.append((min(p1, p2), max(p1, p2)))

    stats.n_incompatible_pairs = len(incompatible)
    stats.incompatible_pairs = incompatible

    # Rm: minimum number of points to stab all intervals
    # Greedy: sort by right endpoint; if interval not yet stabbed, place point at right endpoint
    if incompatible:
        sorted_intervals = sorted(incompatible, key=lambda x: x[1])
        last_point = -1
        rm = 0
        for left, right in sorted_intervals:
            if left > last_point:
                rm += 1
                last_point = right
        stats.Rm = rm

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# Population Size Changes  -  mismatch distribution
# ─────────────────────────────────────────────────────────────────────────────

def compute_mismatch(seqs: list[str]) -> MismatchStats:
    """Mismatch distribution and raggedness statistic.

    Computes the observed distribution of pairwise nucleotide differences,
    the raggedness statistic r (Harpending 1994, equation 1), and the
    coefficient of variation (Rogers & Harpending 1992).

    Raggedness r quantifies the smoothness of the mismatch distribution.
    Small r → smooth (consistent with population expansion).
    Large r → ragged (consistent with constant size or decline).
    """
    n = len(seqs)
    stats = MismatchStats()

    if n < 2 or not seqs[0]:
        return stats

    # Compute all pairwise differences
    diffs: list[int] = []
    for i in range(n):
        for j in range(i + 1, n):
            d = sum(a != b for a, b in zip(seqs[i], seqs[j]))
            diffs.append(d)

    if not diffs:
        return stats

    stats.n_pairs = len(diffs)
    obs: Counter = Counter(diffs)
    stats.observed = dict(sorted(obs.items()))

    stats.mean = sum(diffs) / len(diffs)
    variance = sum((d - stats.mean) ** 2 for d in diffs) / len(diffs)
    stats.variance = variance
    stats.cv = math.sqrt(variance) / stats.mean if stats.mean > 0 else None

    # Raggedness (Harpending 1994, eq 1)
    # r = Σ (f(i) - f(i-1))² where f(i) = proportion of pairs with i differences
    max_d = max(diffs)
    total_pairs = len(diffs)
    f = [obs.get(d, 0) / total_pairs for d in range(max_d + 2)]
    raggedness = sum((f[i] - f[i - 1]) ** 2 for i in range(1, max_d + 2))
    stats.raggedness = raggedness

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# InDel Polymorphism
# ─────────────────────────────────────────────────────────────────────────────

def _identify_indel_events(seqs: list[str]) -> list[InDelEvent]:
    """Identify InDel events from an alignment.

    An InDel event is a maximal run of alignment columns where the same set
    of sequences carries gaps (and the remaining sequences do not).  When the
    gap-bearing set changes mid-run, a new event begins.

    This implements the 'diallelic' option of DnaSP: overlapping InDel events
    (columns where different subsets of sequences are gapped simultaneously)
    are noted but not analysed further in this version.
    """
    if not seqs:
        return []
    L = len(seqs[0])
    events: list[InDelEvent] = []
    current_gap_set: Optional[frozenset] = None
    event_start: Optional[int] = None

    for pos in range(L):
        gap_set = frozenset(i for i, s in enumerate(seqs) if s[pos] == '-')
        non_gap = [s[pos] for i, s in enumerate(seqs) if s[pos] != '-' and s[pos] not in _GAP_CHARS]

        if gap_set and non_gap:  # some sequences gapped, some not → InDel column
            if gap_set == current_gap_set:
                pass  # continuing same event
            else:
                if current_gap_set is not None and event_start is not None:
                    length = pos - event_start
                    events.append(InDelEvent(
                        start=event_start,
                        end=pos - 1,
                        length=length,
                        seq_indices=current_gap_set,
                    ))
                current_gap_set = gap_set
                event_start = pos
        else:
            if current_gap_set is not None and event_start is not None:
                length = pos - event_start
                events.append(InDelEvent(
                    start=event_start,
                    end=pos - 1,
                    length=length,
                    seq_indices=current_gap_set,
                ))
            current_gap_set = None
            event_start = None

    if current_gap_set is not None and event_start is not None:
        length = L - event_start
        events.append(InDelEvent(
            start=event_start,
            end=L - 1,
            length=length,
            seq_indices=current_gap_set,
        ))

    return events


def compute_indel(seqs: list[str]) -> InDelStats:
    """InDel polymorphism statistics (DnaSP InDel module, diallelic option).

    Reference: DnaSP v6 InDel (Insertion-Deletion) Polymorphism module.
    """
    stats = InDelStats()
    n = len(seqs)

    if n < 2:
        return stats

    # Total alignment positions with gaps in any sequence
    L = len(seqs[0])
    stats.n_positions_with_gaps = sum(
        1 for pos in range(L)
        if any(s[pos] == '-' for s in seqs)
    )

    events = _identify_indel_events(seqs)
    stats.events = events
    stats.n_events = len(events)

    if not events:
        return stats

    stats.mean_event_length = sum(e.length for e in events) / len(events)

    # Number of net positions analysed = positions NOT involved in overlapping events
    # For diallelic option: exclude positions that belong to overlapping events
    # Simple approach: use all event positions, note overlap regions
    net_positions = sum(e.length for e in events)

    # Build binary haplotype: for each sequence, a binary vector over events
    # 1 = sequence carries the gap (InDel), 0 = sequence does not
    binary_haplotypes: list[tuple] = []
    for seq_idx in range(n):
        hap = tuple(1 if seq_idx in e.seq_indices else 0 for e in events)
        binary_haplotypes.append(hap)

    # InDel haplotype diversity
    H_indel, Hd_indel, _ = compute_haplotypes(
        ["".join(str(b) for b in h) for h in binary_haplotypes]
    )
    stats.n_haplotypes = H_indel
    stats.haplotype_diversity = Hd_indel

    # k_indel: average pairwise differences in InDel pattern
    diff_sum = 0
    n_pairs = 0
    for i in range(n):
        for j in range(i + 1, n):
            diff_sum += sum(
                binary_haplotypes[i][e] != binary_haplotypes[j][e]
                for e in range(len(events))
            )
            n_pairs += 1
    k_indel = diff_sum / n_pairs if n_pairs > 0 else 0.0
    stats.k_indel = k_indel

    # pi_indel = k_indel / net_positions
    stats.pi_indel = k_indel / net_positions if net_positions > 0 else 0.0

    # theta_indel from number of events (Watterson)
    S_indel = stats.n_events
    theta_abs, _ = watterson_theta(S_indel, n, net_positions)
    stats.theta_indel = theta_abs

    # Tajima's D on InDel data
    stats.tajima_d_indel = tajima_d(k_indel, S_indel, n)

    return stats


# ─────────────────────────────────────────────────────────────────────────────
# Divergence between populations
# ─────────────────────────────────────────────────────────────────────────────

def compute_divergence(
    seqs1: list[str],
    seqs2: list[str],
    pop1_name: str = "Pop1",
    pop2_name: str = "Pop2",
) -> DivergenceStats:
    """DNA divergence between two populations.

    Computes: π per population, Dxy (average between-population divergence
    per site), Da (net divergence), fixed differences, shared polymorphisms,
    and private polymorphisms.

    References:
      Nei (1987) Molecular Evolutionary Genetics, eq 10.20 (Dxy)
      Hey (1991) Genetics 128:831-840 (fixed differences)
      Tajima (1983) Genetics 105:437-460 (k, eq A3)
    """
    stats = DivergenceStats(pop1_name=pop1_name, pop2_name=pop2_name)
    n1, n2 = len(seqs1), len(seqs2)
    stats.n1, stats.n2 = n1, n2

    if n1 < 1 or n2 < 1 or not seqs1[0]:
        return stats

    # Complete deletion across BOTH populations combined
    combined = seqs1 + seqs2
    _, L_net = complete_deletion(combined)
    if L_net == 0:
        return stats

    # Re-apply gap mask to get cleaned sequences
    L = len(seqs1[0])
    gap_cols = set()
    for pos in range(L):
        if any(s[pos] in _GAP_CHARS for s in combined):
            gap_cols.add(pos)
    keep = [pos for pos in range(L) if pos not in gap_cols]

    c1 = ["".join(s[i] for i in keep) for s in seqs1]
    c2 = ["".join(s[i] for i in keep) for s in seqs2]
    stats.L_net = len(keep)

    if not keep:
        return stats

    # Within-population diversity
    stats.k1 = compute_k(c1)
    stats.k2 = compute_k(c2)
    stats.Pi1 = stats.k1 / stats.L_net if stats.L_net > 0 else 0.0
    stats.Pi2 = stats.k2 / stats.L_net if stats.L_net > 0 else 0.0

    # Dxy: average number of pairwise differences between populations, per site
    # Nei (1987) eq 10.20
    total_between = 0
    for s1 in c1:
        for s2 in c2:
            total_between += sum(a != b for a, b in zip(s1, s2))
    Dxy_abs = total_between / (n1 * n2)  # absolute (per gene)
    stats.Dxy = Dxy_abs / stats.L_net if stats.L_net > 0 else 0.0

    # Da: net divergence (removes within-population diversity)
    stats.Da = stats.Dxy - (stats.Pi1 + stats.Pi2) / 2

    # Fixed differences, shared polymorphisms, private polymorphisms
    # (Hey 1991 classification at each site)
    alleles1_per_site: list[frozenset] = []
    alleles2_per_site: list[frozenset] = []
    for pos in range(stats.L_net):
        a1 = frozenset(s[pos] for s in c1 if s[pos] not in _GAP_CHARS)
        a2 = frozenset(s[pos] for s in c2 if s[pos] not in _GAP_CHARS)
        alleles1_per_site.append(a1)
        alleles2_per_site.append(a2)

    n_fixed = n_shared = n_private1 = n_private2 = 0
    for a1, a2 in zip(alleles1_per_site, alleles2_per_site):
        if not a1 or not a2:
            continue
        seg1 = len(a1) > 1
        seg2 = len(a2) > 1
        disjoint = a1.isdisjoint(a2)
        if disjoint and not seg1 and not seg2:
            n_fixed += 1           # fixed difference: one allele each, different
        elif not disjoint and (seg1 or seg2):
            shared = a1 & a2
            if shared:
                n_shared += 1      # shared polymorphism
        elif seg1 and a1.isdisjoint(a2):
            n_private1 += 1        # private to pop1
        elif seg2 and a2.isdisjoint(a1):
            n_private2 += 1        # private to pop2

    stats.n_fixed = n_fixed
    stats.n_shared = n_shared
    stats.n_private1 = n_private1
    stats.n_private2 = n_private2

    return stats


def split_alignment_by_pop(
    aln: Alignment, pop_assignments: dict[str, str]
) -> dict[str, Alignment]:
    """Split an alignment into sub-alignments by population assignment."""
    groups: dict[str, tuple[list[str], list[str]]] = {}
    for name, seq in zip(aln.names, aln.seqs):
        pop = pop_assignments.get(name)
        if pop is None:
            continue
        if pop not in groups:
            groups[pop] = ([], [])
        groups[pop][0].append(name)
        groups[pop][1].append(seq)
    return {
        pop: Alignment(names=names, seqs=seqs, source=aln.source)
        for pop, (names, seqs) in groups.items()
    }


# ─────────────────────────────────────────────────────────────────────────────
# Fu & Li D / F with outgroup
# ─────────────────────────────────────────────────────────────────────────────

def _count_derived(seqs: list[str], outgroup: str) -> tuple[int, int]:
    """Count outgroup-polarised derived mutations.

    Applies complete deletion: any column with a gap in any ingroup sequence OR
    in the outgroup is skipped.

    Parameters
    ----------
    seqs : list[str]
        Ingroup sequences (all same length, already uppercase).
    outgroup : str
        Outgroup sequence (same length as seqs).

    Returns
    -------
    eta : int
        Total number of derived mutations across all polarisable sites.
    eta_e : int
        Derived mutations carried by exactly one ingroup sequence (external).
    """
    if not seqs:
        return 0, 0
    L = len(seqs[0])
    if len(outgroup) < L:
        raise ValueError(
            f"Outgroup length ({len(outgroup)}) shorter than alignment ({L})."
        )
    eta = 0
    eta_e = 0
    for pos in range(L):
        anc = outgroup[pos]
        if anc in _GAP_CHARS:
            continue
        col = [s[pos] for s in seqs]
        if any(c in _GAP_CHARS for c in col):
            continue
        states = frozenset(col)
        if anc not in states or len(states) < 2:
            continue  # monomorphic or ancestral allele absent
        for derived in states - {anc}:
            n_carriers = sum(1 for c in col if c == derived)
            eta += 1
            if n_carriers == 1:
                eta_e += 1
    return eta, eta_e


def compute_fu_li_outgroup(seqs: list[str], outgroup: str) -> FuLiOutgroupStats:
    """Fu & Li D and F statistics using an outgroup to polarise mutations.

    The outgroup sequence identifies which allele is ancestral at each
    segregating site.  Derived mutations carried by exactly one ingroup
    sequence are 'external' (η_e); all derived mutations count as η.

    Variance coefficients follow Simonsen et al. (1995), Appendix B.
    Note: these are the same u/v structure as the no-outgroup Appendix A
    formulas; the exact Appendix B coefficients (which differ slightly for
    small n) can be substituted once validated against DnaSP 6 output.

    Parameters
    ----------
    seqs : list[str]
        Ingroup sequences (already uppercase, all same length).
    outgroup : str
        Outgroup sequence (same length), uppercase.

    Returns
    -------
    FuLiOutgroupStats
    """
    n = len(seqs)
    if n < 4 or not seqs:
        return FuLiOutgroupStats(n=n)

    # k_bar: mean pairwise differences over ALL clean ingroup sites
    clean, L_net = complete_deletion(seqs)
    k_bar = compute_k(clean)

    # Outgroup-polarised counts
    eta, eta_e = _count_derived(seqs, outgroup)
    if eta == 0:
        return FuLiOutgroupStats(n=n, eta=0, eta_e=0, k_bar=k_bar)

    An = _harmonic(n, 1)   # Σ 1/i for i = 1..n-1
    Bn = _harmonic(n, 2)   # Σ 1/i² for i = 1..n-1

    # Variance coefficients (Simonsen et al. 1995, Appendix B structure)
    v_D = (Bn / An**2 - (2.0 / n) * (1.0 + 1.0 / An - An + An / n) - 1.0 / n**2)
    v_D /= An**2 + Bn
    v_D = max(v_D, 0.0)
    u_D = ((n - 1.0) / n - 1.0 / An) / An - v_D

    v_F = (2 * n**3 + 110 * n**2 - 255 * n + 153) / (9 * n**2 * (n - 1))
    v_F += (2 * (n - 1) * An) / n**2
    v_F -= 8 * Bn / n
    v_F /= An**2 + Bn
    v_F = max(v_F, 0.0)
    u_F = (4 * n**2 + 19 * n + 3 - 12 * (n + 1) * (An + 1.0 / n)) / (3 * n * (n - 1))
    u_F = u_F / An - v_F

    # D statistic: compares η_e with expected η/aₙ under neutrality
    num_D = eta_e - eta / An
    denom_D_sq = u_D * eta + v_D * eta * (eta - 1)
    D: Optional[float] = (num_D / math.sqrt(denom_D_sq)) if denom_D_sq > 0 else None

    # F statistic: compares mean pairwise differences with η_e
    num_F = k_bar - eta_e
    denom_F_sq = u_F * eta + v_F * eta * (eta - 1)
    F: Optional[float] = (num_F / math.sqrt(denom_F_sq)) if denom_F_sq > 0 else None

    return FuLiOutgroupStats(n=n, eta=eta, eta_e=eta_e, k_bar=k_bar, D=D, F=F)


# ─────────────────────────────────────────────────────────────────────────────
# HKA test
# ─────────────────────────────────────────────────────────────────────────────

def load_hka_file(path: Path) -> list[HKALocus]:
    """Parse an HKA locus file.

    Expected format (TSV, header optional, lines starting with # ignored)::

        locus   S   D   n
        ACE     5   10  10
        G6PD    2   8   12

    Columns:
        locus : locus name (any string)
        S     : segregating sites in ingroup (int)
        D     : fixed differences vs outgroup/sister species (int)
        n     : ingroup sample size (int)

    Parameters
    ----------
    path : Path

    Returns
    -------
    list[HKALocus]
    """
    loci: list[HKALocus] = []
    with open(path, encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            if parts[1].lower() in ("s", "seg"):
                continue  # header row
            try:
                loci.append(HKALocus(
                    name=parts[0],
                    S=int(parts[1]),
                    D=int(parts[2]),
                    n=int(parts[3]),
                ))
            except ValueError:
                continue  # skip unparseable rows
    return loci


def compute_hka(loci: list[HKALocus]) -> HKAStats:
    """HKA neutrality test (Hudson, Kreitman & Aguadé 1987).

    Tests whether the ratio of polymorphism to divergence is uniform across
    loci. The null model assumes the neutral model with a single underlying
    θ_i per locus, a shared scaled divergence time T, and Poisson sampling.

    MLE of T is found by bisection on the constraint equation derived from
    ∂log L / ∂T = 0::

        Σ D_i / (1 + 2T) = Σ (S_i + D_i) / (f_i + 1 + 2T)

    where f_i = Σ_{j=1}^{n_i - 1} 1/j (Watterson's harmonic number for
    locus i).

    Parameters
    ----------
    loci : list[HKALocus]

    Returns
    -------
    HKAStats
    """
    k = len(loci)
    if k < 2:
        return HKAStats(n_loci=k)

    # Validate loci
    valid = [loc for loc in loci if loc.n >= 2]
    if len(valid) < 2:
        return HKAStats(n_loci=k)

    fs = [_harmonic(loc.n, 1) for loc in valid]  # f_i per locus

    # MLE constraint: Σ D_i/(1+2T) = Σ (S_i+D_i)/(f_i+1+2T)
    def _g(T: float) -> float:
        left = sum(valid[i].D / (1.0 + 2.0 * T) for i in range(len(valid)))
        right = sum(
            (valid[i].S + valid[i].D) / (fs[i] + 1.0 + 2.0 * T)
            for i in range(len(valid))
        )
        return left - right

    T_hat = _bisect(_g, 0.0, 1000.0)

    # MLE θ̂_i and expected values
    theta_hats = [
        (valid[i].S + valid[i].D) / (fs[i] + 1.0 + 2.0 * T_hat)
        for i in range(len(valid))
    ]
    E_S = [theta_hats[i] * fs[i] for i in range(len(valid))]
    E_D = [theta_hats[i] * (1.0 + 2.0 * T_hat) for i in range(len(valid))]

    # Chi-square (Poisson variance approximation: Var[X] = E[X])
    chi2 = 0.0
    loci_results: list[dict] = []
    for i, loc in enumerate(valid):
        if E_S[i] > 0:
            chi2 += (loc.S - E_S[i]) ** 2 / E_S[i]
        if E_D[i] > 0:
            chi2 += (loc.D - E_D[i]) ** 2 / E_D[i]
        loci_results.append({
            "name": loc.name,
            "n": loc.n,
            "S": loc.S,
            "D": loc.D,
            "theta_hat": theta_hats[i],
            "E_S": E_S[i],
            "E_D": E_D[i],
        })

    df = len(valid) - 1  # one parameter (T) estimated from data
    p_value = _chi2_pvalue(chi2, df)

    return HKAStats(
        n_loci=k,
        T_hat=T_hat,
        chi2=chi2,
        df=df,
        p_value=p_value,
        loci_results=loci_results,
    )


# ─────────────────────────────────────────────────────────────────────────────
# Coding sequence helpers (Group C: mk, kaks)
# ─────────────────────────────────────────────────────────────────────────────

def _log_factorial(n: int) -> float:
    """Compute log(n!) using math.lgamma (numerically stable)."""
    return math.lgamma(n + 1)


def _fisher_exact_2x2(a: int, b: int, c: int, d: int) -> float:
    """Two-tailed Fisher's exact test P-value for a 2×2 contingency table.

    Table layout::

        a  b  (row 1)
        c  d  (row 2)

    Given the observed cell count ``a``, sums probabilities of all tables
    (with the same marginals) that are at least as extreme (≤ observed
    hypergeometric probability) under the null of row/column independence.
    """
    R1 = a + b
    R2 = c + d
    C1 = a + c
    C2 = b + d
    N  = R1 + R2

    if N == 0:
        return 1.0

    lf = _log_factorial
    log_const = lf(C1) + lf(C2) - lf(N)   # normalisation factor

    def _log_prob(k: int) -> float:
        k2 = C1 - k
        if k < 0 or k > R1 or k2 < 0 or k2 > R2:
            return -math.inf
        return (lf(R1) - lf(k) - lf(R1 - k) +
                lf(R2) - lf(k2) - lf(R2 - k2) + log_const)

    p_obs = _log_prob(a)
    k_lo  = max(0, C1 - R2)
    k_hi  = min(R1, C1)

    total = 0.0
    for k in range(k_lo, k_hi + 1):
        lp = _log_prob(k)
        if lp <= p_obs + 1e-10:   # as extreme or more extreme
            total += math.exp(lp)
    return min(total, 1.0)


def _count_syn_sites_codon(codon: str) -> float:
    """Synonymous sites in a codon by the Nei-Gojobori (1986) method.

    For each of the 3 positions, counts the fraction of the 3 possible
    single-nucleotide alternatives that are synonymous (same amino acid).
    Returns the sum across positions (0-3 per codon).

    Codons with gaps, ambiguous bases, or stop codons return 0.0.
    """
    if len(codon) != 3 or any(c not in 'ATCG' for c in codon):
        return 0.0
    aa = GENETIC_CODE.get(codon)
    if aa is None or aa == '*':
        return 0.0
    syn = 0.0
    for pos in range(3):
        n_syn = sum(
            1 for alt in _NUCLEOTIDES
            if alt != codon[pos]
            and GENETIC_CODE.get(codon[:pos] + alt + codon[pos + 1:]) == aa
        )
        syn += n_syn / 3.0   # 3 possible alternatives at each position
    return syn


def _classify_codon_pair(c1: str, c2: str) -> tuple[float, float]:
    """Classify the changes between two codons as synonymous or nonsynonymous.

    Uses the Nei-Gojobori (1986) pathway-averaging method.  For codons
    differing at k positions, all k! pathways through intermediate codons are
    enumerated and averaged.  Paths through stop codons are excluded; if all
    paths go through stops, the entire difference is counted as nonsynonymous.

    Returns
    -------
    (syn_changes, nonsyn_changes) : tuple[float, float]
    """
    if len(c1) != 3 or len(c2) != 3:
        return (0.0, 0.0)
    if any(nt not in 'ATCG' for nt in c1 + c2):
        return (0.0, 0.0)
    aa1 = GENETIC_CODE.get(c1)
    aa2 = GENETIC_CODE.get(c2)
    if aa1 is None or aa1 == '*' or aa2 is None or aa2 == '*':
        return (0.0, 0.0)

    diff_pos = [i for i in range(3) if c1[i] != c2[i]]
    n_diff = len(diff_pos)
    if n_diff == 0:
        return (0.0, 0.0)
    if n_diff == 1:
        return (1.0, 0.0) if aa1 == aa2 else (0.0, 1.0)

    from itertools import permutations as _perms
    total_syn = 0.0
    total_nonsyn = 0.0
    n_valid = 0

    for perm in _perms(diff_pos):
        current = c1
        path_syn = 0.0
        path_nonsyn = 0.0
        path_ok = True
        for pos in perm:
            nxt = current[:pos] + c2[pos] + current[pos + 1:]
            cur_aa = GENETIC_CODE.get(current)
            nxt_aa = GENETIC_CODE.get(nxt)
            if cur_aa is None or nxt_aa is None or cur_aa == '*' or nxt_aa == '*':
                path_ok = False
                break
            if cur_aa == nxt_aa:
                path_syn += 1.0
            else:
                path_nonsyn += 1.0
            current = nxt
        if path_ok:
            total_syn += path_syn
            total_nonsyn += path_nonsyn
            n_valid += 1

    if n_valid == 0:
        return (0.0, float(n_diff))   # all paths through stops → all nonsyn
    return (total_syn / n_valid, total_nonsyn / n_valid)


def _jc_correct(p: float) -> Optional[float]:
    """Jukes-Cantor correction: K = −3/4 × ln(1 − 4p/3).

    Returns None when p ≥ 0.75 (divergence is saturated).
    """
    x = 1.0 - 4.0 * p / 3.0
    if x <= 0.0:
        return None
    return -0.75 * math.log(x)


def compute_mk(seqs: list[str], outgroup: str) -> MKStats:
    """McDonald-Kreitman test (McDonald & Kreitman 1991).

    Classifies each codon of an ingroup coding alignment into one of four
    categories using the ingroup sequences and a single outgroup sequence:

    - **Ps**  -  codon is polymorphic within ingroup; all variant pairs are
      synonymous.
    - **Pn**  -  codon is polymorphic within ingroup; at least one variant pair
      is nonsynonymous.
    - **Ds**  -  all ingroup sequences carry the same codon AND it differs from
      the outgroup by a synonymous change (fixed synonymous difference).
    - **Dn**  -  same as Ds but the fixed difference is nonsynonymous.

    Codons are skipped (complete deletion) when any ingroup sequence or the
    outgroup has a gap, ambiguous base, or stop codon at those three positions.
    Codons that are simultaneously polymorphic in the ingroup *and* divergent
    from the outgroup are excluded (conservative).

    Derived statistics:

    - **α** = 1 − (Ds·Pn)/(Dn·Ps): proportion of adaptive substitutions.
    - **NI** = (Pn/Ps)/(Dn/Ds): Neutrality Index.
    - **DoS** = Dn/(Dn+Ds) − Pn/(Pn+Ps): Direction of Selection.
    - **Fisher P**: two-tailed Fisher's exact test on the 2×2 table.

    Parameters
    ----------
    seqs : list[str]
        Ingroup sequences (pre-aligned, same length, length divisible by 3).
    outgroup : str
        Outgroup sequence (same length as ingroup).

    Returns
    -------
    MKStats
    """
    result = MKStats()
    n = len(seqs)
    if n < 2:
        return result
    L = len(seqs[0])
    if L == 0 or L % 3 != 0 or len(outgroup) != L:
        return result

    n_codons = L // 3
    Pn = Ps = Dn = Ds = 0

    for ci in range(n_codons):
        in_codons = [s[ci * 3:(ci + 1) * 3] for s in seqs]
        out_codon  = outgroup[ci * 3:(ci + 1) * 3]

        # Complete deletion: skip if any base is ambiguous in ingroup or outgroup
        if any(nt not in 'ATCG' for nt in out_codon):
            continue
        if any(nt not in 'ATCG' for codon in in_codons for nt in codon):
            continue
        # Skip if outgroup or any ingroup variant is a stop codon
        if GENETIC_CODE.get(out_codon) == '*':
            continue
        unique_in = set(in_codons)
        if any(GENETIC_CODE.get(co) == '*' for co in unique_in):
            continue

        is_poly  = len(unique_in) > 1
        is_fixed = (len(unique_in) == 1) and (next(iter(unique_in)) != out_codon)

        if is_poly:
            has_nonsyn = False
            has_syn    = False
            for ca, cb in combinations(sorted(unique_in), 2):
                s_d, ns_d = _classify_codon_pair(ca, cb)
                if ns_d > 0:
                    has_nonsyn = True
                if s_d > 0:
                    has_syn = True
            if has_nonsyn:
                Pn += 1
            elif has_syn:
                Ps += 1

        if is_fixed:
            in_codon = next(iter(unique_in))
            s_d, ns_d = _classify_codon_pair(in_codon, out_codon)
            if ns_d > 0:
                Dn += 1
            elif s_d > 0:
                Ds += 1

    result.Pn = Pn
    result.Ps = Ps
    result.Dn = Dn
    result.Ds = Ds

    # α = 1 − (Ds·Pn) / (Dn·Ps)
    if Dn > 0 and Ps > 0:
        result.alpha = 1.0 - (Ds * Pn) / (Dn * Ps)

    # NI = (Pn/Ps) / (Dn/Ds)
    if Ps > 0 and Ds > 0:
        result.NI = (Pn / Ps) / (Dn / Ds)

    # DoS = Dn/(Dn+Ds) − Pn/(Pn+Ps)
    if (Dn + Ds) > 0 and (Pn + Ps) > 0:
        result.DoS = Dn / (Dn + Ds) - Pn / (Pn + Ps)

    # Fisher's exact test on:
    #   Pn  Dn   (nonsynonymous row)
    #   Ps  Ds   (synonymous row)
    result.fisher_p = _fisher_exact_2x2(Pn, Dn, Ps, Ds)

    return result


def compute_ka_ks(seqs: list[str]) -> KaKsStats:
    """Ka/Ks (dN/dS) by the Nei-Gojobori (1986) method.

    Estimates synonymous (Ks) and nonsynonymous (Ka) substitution rates by:

    1. Computing synonymous (S_i) and nonsynonymous (N_i = 3L − S_i) site
       counts for each sequence using ``_count_syn_sites_codon``.
    2. For each pair (i, j), counting synonymous differences (sd) and
       nonsynonymous differences (nd) using ``_classify_codon_pair``.
    3. Computing proportions pS = sd / S_ij and pN = nd / N_ij where
       S_ij = (S_i + S_j)/2.
    4. Applying the Jukes-Cantor correction to obtain Ks_ij and Ka_ij.
    5. Averaging Ks and Ka over all pairs.

    Codons with gaps, ambiguous bases, or stop codons in either sequence are
    excluded (complete deletion at the codon level). Pairs where the
    Jukes-Cantor correction is undefined (pS or pN ≥ 0.75) are excluded from
    the relevant average.

    Parameters
    ----------
    seqs : list[str]
        Ingroup sequences (pre-aligned, same length, length divisible by 3).

    Returns
    -------
    KaKsStats
        Summary statistics averaged over all pairwise comparisons.
        ``omega`` = Ka/Ks; < 1 purifying selection, ≈ 1 neutral, > 1 positive.
    """
    result = KaKsStats()
    n = len(seqs)
    if n < 2:
        return result
    L = len(seqs[0])
    if L == 0 or L % 3 != 0:
        return result

    n_codons = L // 3
    result.n_codons = n_codons

    # Per-sequence synonymous site counts
    S_per_seq = []
    for seq in seqs:
        S = 0.0
        for ci in range(n_codons):
            codon = seq[ci * 3:(ci + 1) * 3]
            if any(nt not in 'ATCG' for nt in codon):
                continue
            if GENETIC_CODE.get(codon) == '*':
                continue
            S += _count_syn_sites_codon(codon)
        S_per_seq.append(S)

    result.S_sites = sum(S_per_seq) / n
    result.N_sites = L - result.S_sites

    # Pairwise computation
    Ks_all: list[float] = []
    Ka_all: list[float] = []
    Sd_all: list[float] = []
    Nd_all: list[float] = []

    for i, j in combinations(range(n), 2):
        total_S  = 0.0
        total_N  = 0.0
        total_sd = 0.0
        total_nd = 0.0
        n_valid  = 0

        for ci in range(n_codons):
            c1 = seqs[i][ci * 3:(ci + 1) * 3]
            c2 = seqs[j][ci * 3:(ci + 1) * 3]
            if any(nt not in 'ATCG' for nt in c1 + c2):
                continue
            if GENETIC_CODE.get(c1) == '*' or GENETIC_CODE.get(c2) == '*':
                continue
            n_valid += 1
            s_sites = (_count_syn_sites_codon(c1) + _count_syn_sites_codon(c2)) / 2.0
            total_S  += s_sites
            total_N  += 3.0 - s_sites
            sd, nd    = _classify_codon_pair(c1, c2)
            total_sd += sd
            total_nd += nd

        if n_valid == 0 or total_S <= 0.0 or total_N <= 0.0:
            continue

        Sd_all.append(total_sd)
        Nd_all.append(total_nd)

        ks = _jc_correct(total_sd / total_S)
        ka = _jc_correct(total_nd / total_N)
        if ks is not None:
            Ks_all.append(ks)
        if ka is not None:
            Ka_all.append(ka)

    if Sd_all:
        result.Sd = sum(Sd_all) / len(Sd_all)
    if Nd_all:
        result.Nd = sum(Nd_all) / len(Nd_all)
    if Ks_all:
        result.Ks = sum(Ks_all) / len(Ks_all)
    if Ka_all:
        result.Ka = sum(Ka_all) / len(Ka_all)
    if result.Ka is not None and result.Ks is not None and result.Ks > 0.0:
        result.omega = result.Ka / result.Ks

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Frequency spectrum helpers (Group D: fufs, sfs)
# ─────────────────────────────────────────────────────────────────────────────

def _stirling1_unsigned(n: int) -> list[int]:
    """Row n of unsigned Stirling numbers of the first kind.

    Returns a list of length n+1 where entry k holds |s(n, k)|.
    Recurrence: |s(n+1, k)| = n × |s(n, k)| + |s(n, k-1)|
    Uses Python arbitrary-precision integers  -  no overflow.
    """
    row: list[int] = [0] * (n + 1)
    row[0] = 1  # |s(0, 0)| = 1 by convention
    for i in range(1, n + 1):
        new_row = [0] * (n + 1)
        for k in range(1, i + 1):
            new_row[k] = (i - 1) * row[k] + (row[k - 1] if k >= 1 else 0)
        row = new_row
    return row


def _ewens_cdf(k_max: int, n: int, theta: float) -> float:
    """P(K_n ≤ k_max) under the Ewens sampling formula.

    K_n is the number of distinct alleles in a sample of n sequences under
    the infinite-alleles model with scaled mutation rate theta (= θ_π).

    P(K_n = k) = |s(n, k)| × θ^k / [θ(θ+1)…(θ+n−1)]
    """
    if n < 1 or theta <= 0.0:
        return 1.0
    k_max = min(k_max, n)
    stirling = _stirling1_unsigned(n)
    # Rising factorial θ^(n) = θ(θ+1)…(θ+n−1)
    rising: float = 1.0
    for i in range(n):
        rising *= theta + i
    if rising == 0.0:
        return 1.0
    prob_sum: float = 0.0
    theta_pow: float = 1.0
    for k in range(1, k_max + 1):
        theta_pow *= theta   # θ^k
        prob_sum += float(stirling[k]) * theta_pow / rising
    return min(1.0, max(0.0, prob_sum))


def compute_fu_fs(seqs: list[str], H: int, k: float) -> FuFsStats:
    """Fu's Fs neutrality test (Fu 1997).

    Uses the Ewens sampling formula to evaluate how extreme the observed
    number of haplotypes H is given the nucleotide diversity estimate θ_π = k
    (mean pairwise differences, which equals π × L_net).

    Parameters
    ----------
    seqs : list[str]   Sequences after complete deletion (same length, clean).
    H    : int         Observed number of haplotypes (distinct sequences).
    k    : float       Mean pairwise differences (θ_π = π × L_net).

    Returns
    -------
    FuFsStats with Fs, S_k, theta_pi filled.
    """
    n = len(seqs)
    result = FuFsStats(n=n, H=H, theta_pi=k)
    if n < 2 or H < 1:
        return result

    S_k = _ewens_cdf(H, n, k)
    result.S_k = S_k

    if S_k <= 0.0:
        result.Fs = -1e308           # effectively -inf
    elif S_k >= 1.0:
        result.Fs = 1e308            # effectively +inf
    else:
        result.Fs = math.log(S_k / (1.0 - S_k))

    return result


def compute_sfs(seqs: list[str], outgroup_seq: Optional[str] = None) -> SFSStats:
    """Site frequency spectrum  -  folded and (with outgroup) unfolded.

    Parameters
    ----------
    seqs         : list[str]          Original (unclean) ingroup sequences.
    outgroup_seq : str | None         Outgroup sequence (same length) for
                                      unfolded SFS.  None → folded only.

    Returns
    -------
    SFSStats with folded dict and, if outgroup_seq is not None, unfolded dict.

    Gap treatment
    -------------
    Folded    -  complete deletion on ingroup only: skip any column where any
               ingroup sequence has a gap/ambiguous base.
    Unfolded  -  additionally requires the outgroup to be a clean ATCG base
               at that column; gap in outgroup → column skipped for unfolded
               but still counted for folded if ingroup is clean.
    """
    n = len(seqs)
    result = SFSStats(n=n, has_outgroup=outgroup_seq is not None)
    if n < 2 or not seqs:
        return result

    L = len(seqs[0])
    folded: dict[int, int] = {}
    unfolded: Optional[dict[int, int]] = {} if outgroup_seq is not None else None

    for col in range(L):
        # Ingroup column
        ing_col = [s[col] for s in seqs]
        # Skip if any ingroup base is not clean ATCG
        if any(c not in 'ATCG' for c in ing_col):
            continue

        counts = Counter(ing_col)
        alleles = [a for a, cnt in counts.items()]
        if len(alleles) < 2:
            continue  # monomorphic  -  not a segregating site

        # ── Folded SFS ──────────────────────────────────────────────────
        # Minor allele = least frequent; fold at n/2
        minor_count = min(counts.values())
        folded_key = minor_count  # i = count of minor allele (1 ≤ i ≤ n//2)
        folded[folded_key] = folded.get(folded_key, 0) + 1

        # ── Unfolded SFS (needs outgroup) ────────────────────────────────
        if outgroup_seq is not None and unfolded is not None and col < len(outgroup_seq):
            og = outgroup_seq[col]
            if og in 'ATCG' and og in counts:
                # Ancestral allele = outgroup; derived = everything else
                derived_count = n - counts[og]
                if 1 <= derived_count <= n - 1:
                    unfolded[derived_count] = unfolded.get(derived_count, 0) + 1

    result.folded = folded
    result.unfolded = unfolded
    return result


# ─────────────────────────────────────────────────────────────────────────────
# Substitution pattern helpers (Group E: tstv, codon)
# ─────────────────────────────────────────────────────────────────────────────

def compute_ts_tv(seqs: list[str]) -> TsTvStats:
    """Compute transition/transversion ratio across all pairwise comparisons.

    A transition (Ts) is a purine↔purine (A↔G) or pyrimidine↔pyrimidine (C↔T)
    substitution.  A transversion (Tv) is a purine↔pyrimidine substitution.
    Ambiguous or gap characters at either position are skipped (complete deletion
    applied per-column: any gap in any sequence drops that column).

    Parameters
    ----------
    seqs : list[str]
        Ingroup sequences (upper-case; equal length).

    Returns
    -------
    TsTvStats
        n_transitions / n_transversions are totals across all pairs.
        ts_per_site / tv_per_site are means per pair per site.
        ts_tv is None when n_transversions == 0.
    """
    result = TsTvStats(n=len(seqs))
    if len(seqs) < 2:
        return result

    L = len(seqs[0])
    clean_cols: list[int] = []
    for col in range(L):
        bases = [s[col] for s in seqs]
        if all(b in _NUCLEOTIDES for b in bases):
            clean_cols.append(col)

    result.L_net = len(clean_cols)
    if not clean_cols:
        return result

    n_seqs = len(seqs)
    total_ts = 0
    total_tv = 0
    n_pairs = n_seqs * (n_seqs - 1) // 2

    for i in range(n_seqs - 1):
        for j in range(i + 1, n_seqs):
            for col in clean_cols:
                a, b = seqs[i][col], seqs[j][col]
                if a == b:
                    continue
                # Both already confirmed ATCG by clean_cols filter
                if (a in _PURINES) == (b in _PURINES):
                    total_ts += 1
                else:
                    total_tv += 1

    result.n_transitions = total_ts
    result.n_transversions = total_tv

    if n_pairs > 0 and result.L_net > 0:
        denom = n_pairs * result.L_net
        result.ts_per_site = total_ts / denom
        result.tv_per_site = total_tv / denom

    if total_tv > 0:
        result.ts_tv = total_ts / total_tv

    return result


def compute_codon_usage(seqs: list[str]) -> CodonUsageStats:
    """Compute codon usage bias: RSCU and ENC.

    Sequences must already be in-frame coding alignments.  Each sequence is
    read in non-overlapping triplets.  Triplets containing gap/ambiguous
    characters or translating to a stop codon are skipped per sequence.

    RSCU (Sharp & Li 1987):
        RSCU_ij = X_ij / (X_i / n_i)
    where X_ij is the count of codon j for amino acid i, X_i is the total
    count for amino acid i, and n_i is the synonymous family size.
    RSCU = 1.0 → uniform usage; > 1.0 → preferred; < 1.0 → avoided.

    ENC (Wright 1990):
        ENC = 2 + 9/F_2 + 1/F_3 + 5/F_4 + 3/F_6
    where F_k is the mean corrected homozygosity for amino acids with k-fold
    degeneracy:
        F_k = (n_aa * Σ p_j² - 1) / (n_aa - 1)   (n_aa = total codons for aa)
    Degeneracy classes (standard genetic code):
        2-fold: 9 amino acids (F, L[2], I, M→skip, V, S[2], P[skip], A, T, C, Y, H, Q, N, K, D, E → refined below)
        Exact mapping used: {2: Cys,Asp,Glu,Phe,His,Lys,Asn,Gln,Tyr}  (9 aa, 2 codons each)
                            {3: Ile}                                      (1 aa, 3 codons)
                            {4: Ala,Gly,Pro,Thr,Val}                     (5 aa, 4 codons each)
                            {6: Arg,Leu,Ser}                              (3 aa, 6 codons each)
        Met (1 codon) and Trp (1 codon) are excluded (no synonymy).
        ENC ranges from 20 (maximum bias) to 61 (no bias).

    Parameters
    ----------
    seqs : list[str]
        Ingroup sequences (upper-case; in-frame coding alignment).

    Returns
    -------
    CodonUsageStats
        codon_counts : mean counts per sequence (pooled then divided by n).
        rscu         : RSCU for all 61 sense codons.
        ENC          : None if any degeneracy class has no data.
    """
    result = CodonUsageStats(n=len(seqs))
    if len(seqs) < 1:
        return result

    # Pool raw codon counts across all sequences
    from collections import Counter
    raw: Counter[str] = Counter()
    total_codons_per_seq: list[int] = []

    for seq in seqs:
        seq_count = 0
        # Strip gaps so we read in-frame (complete-deletion at codon level)
        # We process the gapped alignment triplet-by-triplet; any gap in triplet → skip
        L = len(seq)
        for start in range(0, L - 2, 3):
            triplet = seq[start:start + 3]
            if len(triplet) < 3:
                break
            if any(b not in _NUCLEOTIDES for b in triplet):
                continue  # gap or ambiguous
            aa = GENETIC_CODE.get(triplet, '*')
            if aa == '*':
                continue  # stop codon
            raw[triplet] += 1
            seq_count += 1
        total_codons_per_seq.append(seq_count)

    n = len(seqs)
    result.n_codons = sum(total_codons_per_seq) / n if n > 0 else 0.0

    if not raw:
        return result

    # Mean counts per sequence
    result.codon_counts = {codon: cnt / n for codon, cnt in raw.items()}

    # ── RSCU ─────────────────────────────────────────────────────────────────
    rscu: dict[str, float] = {}
    for aa, codons in _SYNONYMOUS_FAMILIES.items():
        n_syn = len(codons)                       # synonymous family size
        total_aa = sum(raw.get(c, 0) for c in codons)
        expected = total_aa / n_syn if total_aa > 0 else 0.0
        for c in codons:
            if expected > 0:
                rscu[c] = raw.get(c, 0) / expected
            else:
                rscu[c] = 0.0
    result.rscu = rscu

    # ── ENC (Wright 1990) ────────────────────────────────────────────────────
    # Determine degeneracy class for each amino acid (from GENETIC_CODE)
    # Exclude Met (ATG only) and Trp (TGG only)  -  n_i = 1, no synonymy.
    deg_classes: dict[int, list[str]] = {}  # degeneracy → list of amino acids
    for aa, codons in _SYNONYMOUS_FAMILIES.items():
        k = len(codons)
        if k < 2:
            continue  # Met, Trp  -  single codon, excluded from ENC
        deg_classes.setdefault(k, []).append(aa)

    # Corrected homozygosity F_k for each degeneracy class k
    # F_k = mean over amino acids in that class of:
    #       (n_aa * Σ p_j² - 1) / (n_aa - 1)   where n_aa = total codons for aa
    def _mean_F(aa_list: list[str]) -> Optional[float]:
        """Mean corrected homozygosity for a set of amino acids."""
        F_values: list[float] = []
        for aa in aa_list:
            codons = _SYNONYMOUS_FAMILIES[aa]
            n_aa = sum(raw.get(c, 0) for c in codons)
            if n_aa < 2:
                # Insufficient data for this amino acid; skip it
                continue
            sum_p2 = sum((raw.get(c, 0) / n_aa) ** 2 for c in codons)
            F_aa = (n_aa * sum_p2 - 1) / (n_aa - 1)
            F_values.append(F_aa)
        if not F_values:
            return None
        return sum(F_values) / len(F_values)

    # ENC = 2 + 9/F_2 + 1/F_3 + 5/F_4 + 3/F_6
    # Coefficients are the number of amino acids in each class (standard code)
    enc_components: dict[int, tuple[int, float]] = {}  # k → (n_aa_in_class, F_k)
    required_classes = {2: 9, 3: 1, 4: 5, 6: 3}

    enc_ok = True
    for k, n_aa_expected in required_classes.items():
        aa_list = deg_classes.get(k, [])
        F_k = _mean_F(aa_list)
        if F_k is None or F_k <= 0:
            enc_ok = False
            break
        enc_components[k] = (n_aa_expected, F_k)

    if enc_ok:
        enc = 2.0
        for k, (coeff, F_k) in enc_components.items():
            enc += coeff / F_k
        # Clamp to biological range
        result.ENC = max(20.0, min(61.0, enc))

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Group F  -  Fay & Wu's H, Zeng's E, and Fst
# ─────────────────────────────────────────────────────────────────────────────

def compute_fay_wu(seqs: list[str], outgroup: str) -> FayWuStats:
    """Fay & Wu's H (2000) and Zeng's E (2006) from an outgroup-polarised SFS.

    Applies complete deletion: a site is used only if all ingroup sequences AND
    the outgroup have an unambiguous ATCG nucleotide.  A segregating site is
    'polarisable' when the outgroup allele appears in the ingroup (i.e. the
    ancestral allele is present).

    Parameters
    ----------
    seqs : list[str]
        Ingroup sequences (uppercase, equal length).
    outgroup : str
        Outgroup sequence (same length as seqs).

    Returns
    -------
    FayWuStats with theta_pi, theta_w, theta_h, theta_l (all per site), H, E.
    """
    result = FayWuStats(n=len(seqs))
    n = len(seqs)
    if n < 2 or not seqs:
        return result
    L = len(seqs[0])
    if len(outgroup) < L:
        raise ValueError(
            f"Outgroup length ({len(outgroup)}) shorter than alignment ({L})."
        )

    # Step 1: complete deletion with outgroup  -  keep only clean ATCG columns
    clean_cols = [
        col for col in range(L)
        if (outgroup[col] in _NUCLEOTIDES
            and all(s[col] in _NUCLEOTIDES for s in seqs))
    ]
    L_net = len(clean_cols)
    result.L_net = L_net
    if L_net == 0:
        return result

    # Step 2: build unfolded SFS from polarised sites
    # xi[i] = number of sites where derived allele has count i (1 ≤ i ≤ n-1)
    xi: dict[int, int] = {}
    n_polarised = 0

    for col in clean_cols:
        anc = outgroup[col]
        col_bases = [s[col] for s in seqs]
        states = frozenset(col_bases)
        if len(states) < 2:
            continue  # monomorphic site  -  skip
        if anc not in states:
            continue  # can't polarise (ancestral allele absent from ingroup)
        # For each derived allele at this site
        for derived in states - {anc}:
            i = sum(1 for b in col_bases if b == derived)
            if 1 <= i <= n - 1:
                xi[i] = xi.get(i, 0) + 1
                n_polarised += 1

    result.n_polarised = n_polarised
    if n_polarised == 0:
        return result

    # Step 3: compute θ estimates per site
    # θ_π = Σ_i ξ_i * 2*i*(n-i) / [n*(n-1)] / L_net
    # θ_H = Σ_i ξ_i * 2*i²    / [n*(n-1)] / L_net
    # θ_L = Σ_i ξ_i * i        / (n-1)     / L_net
    # θ_W = n_polarised / a1                / L_net
    denom_pi = n * (n - 1)
    a1 = _harmonic(n)  # Σ_{k=1}^{n-1} 1/k  (_harmonic(n) sums i=1..n-1)

    sum_pi = 0.0  # Σ ξ_i * 2*i*(n-i)
    sum_h  = 0.0  # Σ ξ_i * 2*i²
    sum_l  = 0.0  # Σ ξ_i * i

    for i, xi_i in xi.items():
        sum_pi += xi_i * 2 * i * (n - i)
        sum_h  += xi_i * 2 * i * i
        sum_l  += xi_i * i

    theta_pi = sum_pi / denom_pi / L_net
    theta_h  = sum_h  / denom_pi / L_net
    theta_l  = sum_l  / (n - 1)  / L_net
    theta_w  = n_polarised / a1  / L_net

    result.theta_pi = theta_pi
    result.theta_w  = theta_w
    result.theta_h  = theta_h
    result.theta_l  = theta_l
    result.H = theta_pi - theta_h
    result.E = theta_l  - theta_w

    return result


def compute_fst(pop_seqs: dict[str, list[str]]) -> FstStats:
    """Population differentiation using Hudson et al. (1992) pairwise Fst.

    For each pair of populations A and B:
      π_s = (π_A + π_B) / 2   (mean within-population diversity, per site)
      π_t = Dxy_AB             (between-population diversity, per site)
      Fst = 1 − π_s / π_t

    Fst = 0 → populations identical; Fst = 1 → fixed differences only.
    Returns None for a pair when π_t == 0 (no between-pop variation).

    Parameters
    ----------
    pop_seqs : dict[str, list[str]]
        Map of population label → list of sequences (uppercase, equal length).
        All sequences across all pops must be the same length.

    Returns
    -------
    FstStats with per-population π, pairwise Dxy, pairwise Fst, and mean Fst.
    """
    pop_names = sorted(pop_seqs.keys())
    result = FstStats(n_pops=len(pop_names), pop_names=pop_names)

    if len(pop_names) < 2:
        return result

    # Determine L_net via complete deletion across all sequences
    all_seqs = [s for seqs in pop_seqs.values() for s in seqs]
    _, L_net = complete_deletion(all_seqs)

    if L_net == 0:
        for pop in pop_names:
            result.pop_sizes[pop] = len(pop_seqs[pop])
            result.pi_within[pop] = 0.0
        for i, p1 in enumerate(pop_names):
            for p2 in pop_names[i + 1:]:
                result.pi_between[(p1, p2)] = 0.0
                result.fst_pairwise[(p1, p2)] = None
        return result

    # Clean sequences per pop (same complete-deletion mask)
    # Re-derive clean columns from the global mask
    L = len(all_seqs[0])
    clean_cols = [
        col for col in range(L)
        if all(s[col] in _NUCLEOTIDES for s in all_seqs)
    ]
    L_clean = len(clean_cols)

    def _pi_within(seqs_pop: list[str]) -> float:
        """Mean pairwise differences per site within a population."""
        n = len(seqs_pop)
        if n < 2:
            return 0.0
        total = 0.0
        for i in range(n - 1):
            for j in range(i + 1, n):
                total += sum(
                    1 for col in clean_cols
                    if seqs_pop[i][col] != seqs_pop[j][col]
                )
        return total / (n * (n - 1) / 2) / L_clean

    def _dxy(seqs_a: list[str], seqs_b: list[str]) -> float:
        """Mean pairwise differences per site between two populations."""
        na, nb = len(seqs_a), len(seqs_b)
        if na == 0 or nb == 0:
            return 0.0
        total = 0.0
        for sa in seqs_a:
            for sb in seqs_b:
                total += sum(
                    1 for col in clean_cols
                    if sa[col] != sb[col]
                )
        return total / (na * nb) / L_clean

    # Compute within-pop π for each population
    for pop in pop_names:
        result.pop_sizes[pop] = len(pop_seqs[pop])
        result.pi_within[pop] = _pi_within(pop_seqs[pop])

    # Compute pairwise Dxy and Fst
    fst_values: list[float] = []
    for i, p1 in enumerate(pop_names):
        for p2 in pop_names[i + 1:]:
            dxy = _dxy(pop_seqs[p1], pop_seqs[p2])
            result.pi_between[(p1, p2)] = dxy
            pi_s = (result.pi_within[p1] + result.pi_within[p2]) / 2.0
            if dxy > 0:
                fst = 1.0 - pi_s / dxy
                # Clamp to [0, 1]  -  negative values arise with small samples
                fst = max(0.0, min(1.0, fst))
                result.fst_pairwise[(p1, p2)] = fst
                fst_values.append(fst)
            else:
                result.fst_pairwise[(p1, p2)] = None

    if fst_values:
        result.fst_mean = sum(fst_values) / len(fst_values)

    return result


# ─────────────────────────────────────────────────────────────────────────────
# Core analysis pipeline
# ─────────────────────────────────────────────────────────────────────────────

def analyse_region(
    seqs: list[str],
    names: list[str],
    region_label: str,
    L_total: int,
) -> RegionStats:
    """Run polymorphism + neutrality tests on a (possibly windowed) slice."""
    n = len(seqs)
    stats = RegionStats(region=region_label, n=n, L_total=L_total)

    if n < 2:
        return stats

    clean, L_net = complete_deletion(seqs)
    stats.L_net = L_net

    if L_net == 0:
        return stats

    stats.S, stats.Eta = compute_segregating(clean)
    stats.H, stats.Hd, stats.VarHd = compute_haplotypes(clean)
    stats.GC = compute_gc(clean)
    stats.k = compute_k(clean)
    stats.Pi = stats.k / L_net if L_net > 0 else 0.0
    stats.ThetaW, stats.ThetaW_nuc = watterson_theta(stats.S, n, L_net)
    stats.TajimaD = tajima_d(stats.k, stats.S, n)
    eta_s, _ = compute_singletons(clean)
    stats.FuLiD_star, stats.FuLiF_star = fu_li_d_star_f_star(
        stats.k, stats.S, eta_s, n
    )
    stats.R2 = ramos_onsins_r2(clean, stats.k, stats.Eta)

    return stats


def run_analysis(
    aln: Alignment,
    window_size: int = 0,
    step_size: int = 0,
    analyses: Optional[set[str]] = None,
    pop_assignments: Optional[dict[str, str]] = None,
    aln2: Optional[Alignment] = None,
    outgroup: Optional[str] = None,
    hka_loci: Optional[list[HKALocus]] = None,
) -> dict:
    """Run all requested analyses and return a results bundle.

    Parameters
    ----------
    aln : Alignment
    window_size : int
        Sliding window size in bp (0 = whole alignment only).
    step_size : int
        Sliding window step in bp.
    analyses : set[str] | None
        Which analyses to run.  None defaults to {{'polymorphism'}}.
    pop_assignments : dict[str, str] | None
        Sequence-name → population label mapping (for divergence analysis).
    aln2 : Alignment | None
        Second alignment (for divergence analysis).
    outgroup : str | None
        Outgroup sequence string (same length as aln) for fuliout analysis.
    hka_loci : list[HKALocus] | None
        Pre-loaded HKA loci for the hka analysis.

    Returns
    -------
    dict with keys: 'global', 'windows', 'ld', 'recombination',
    'popsize', 'indel', 'divergence', 'fuliout', 'hka'.
    """
    if analyses is None:
        analyses = {"polymorphism"}
    if "all" in analyses:
        analyses = VALID_ANALYSES.copy()

    L = aln.L
    global_label = f"1-{L}"

    # Always run polymorphism
    global_stats = analyse_region(aln.seqs, aln.names, global_label, L)
    window_stats: list[RegionStats] = []

    if window_size > 0 and step_size > 0:
        pos = 0
        while pos + window_size <= L:
            slices = [s[pos: pos + window_size] for s in aln.seqs]
            label = f"{pos + 1}-{pos + window_size}"
            ws = analyse_region(slices, aln.names, label, window_size)
            window_stats.append(ws)
            pos += step_size

    results: dict = {
        "global": global_stats,
        "windows": window_stats,
        "ld": None,
        "recombination": None,
        "popsize": None,
        "indel": None,
        "divergence": None,
        "fuliout": None,
        "hka": None,
        "mk": None,
        "kaks": None,
        "fufs": None,
        "sfs": None,
        "tstv": None,
        "codon": None,
        "faywu": None,
        "fst": None,
    }

    # Get cleaned sequences for non-polymorphism analyses
    clean, _ = complete_deletion(aln.seqs)

    if "ld" in analyses:
        results["ld"] = compute_ld(clean)

    if "recombination" in analyses:
        results["recombination"] = compute_recombination(clean)

    if "popsize" in analyses:
        results["popsize"] = compute_mismatch(clean)

    if "indel" in analyses:
        results["indel"] = compute_indel(aln.seqs)  # use original seqs (keep gaps)

    if "divergence" in analyses:
        # Determine population split
        if aln2 is not None:
            # Two separate files provided
            div_stats = compute_divergence(
                aln.seqs, aln2.seqs,
                pop1_name=Path(aln.source).stem,
                pop2_name=Path(aln2.source).stem,
            )
            results["divergence"] = div_stats
        elif pop_assignments:
            pops = split_alignment_by_pop(aln, pop_assignments)
            pop_names = sorted(pops.keys())
            if len(pop_names) >= 2:
                a = pops[pop_names[0]]
                b = pops[pop_names[1]]
                div_stats = compute_divergence(
                    a.seqs, b.seqs,
                    pop1_name=pop_names[0],
                    pop2_name=pop_names[1],
                )
                results["divergence"] = div_stats
            else:
                print("Warning: divergence analysis requires at least 2 populations in --pop-file", file=sys.stderr)
        else:
            print("Warning: divergence analysis requires --input2 or --pop-file", file=sys.stderr)

    if "fuliout" in analyses:
        if outgroup is not None:
            results["fuliout"] = compute_fu_li_outgroup(aln.seqs, outgroup)
        else:
            print("Warning: fuliout analysis requires --outgroup <seq_name>", file=sys.stderr)

    if "hka" in analyses:
        if hka_loci:
            results["hka"] = compute_hka(hka_loci)
        else:
            print("Warning: hka analysis requires --hka-file <file>", file=sys.stderr)

    if "mk" in analyses:
        if outgroup is not None:
            if aln.L % 3 != 0:
                print(
                    f"Warning: alignment length ({aln.L} bp) is not divisible by 3  -  "
                    "mk analysis requires an in-frame coding sequence.",
                    file=sys.stderr,
                )
            else:
                results["mk"] = compute_mk(aln.seqs, outgroup)
        else:
            print("Warning: mk analysis requires --outgroup <seq_name>", file=sys.stderr)

    if "kaks" in analyses:
        if aln.L % 3 != 0:
            print(
                f"Warning: alignment length ({aln.L} bp) is not divisible by 3  -  "
                "kaks analysis requires an in-frame coding sequence.",
                file=sys.stderr,
            )
        else:
            results["kaks"] = compute_ka_ks(aln.seqs)

    if "fufs" in analyses:
        # Uses polymorphism stats already computed by analyse_region
        results["fufs"] = compute_fu_fs(clean, global_stats.H, global_stats.k)

    if "sfs" in analyses:
        # Uses original (uncleaned) seqs so complete deletion is handled internally
        results["sfs"] = compute_sfs(aln.seqs, outgroup)

    if "tstv" in analyses:
        results["tstv"] = compute_ts_tv(clean)

    if "codon" in analyses:
        results["codon"] = compute_codon_usage(aln.seqs)

    if "faywu" in analyses:
        if outgroup is not None:
            results["faywu"] = compute_fay_wu(aln.seqs, outgroup)
        else:
            print("Warning: faywu analysis requires --outgroup <seq_name>", file=sys.stderr)

    if "fst" in analyses:
        if pop_assignments:
            pops = split_alignment_by_pop(aln, pop_assignments)
            if len(pops) >= 2:
                results["fst"] = compute_fst({k: v.seqs for k, v in pops.items()})
            else:
                print("Warning: fst analysis requires at least 2 populations in --pop-file", file=sys.stderr)
        else:
            print("Warning: fst analysis requires --pop-file", file=sys.stderr)

    return results


# ─────────────────────────────────────────────────────────────────────────────
# Formatting helpers
# ─────────────────────────────────────────────────────────────────────────────

def _fmt(v: Optional[float], digits: int = 6) -> str:
    return "n.a." if v is None else f"{v:.{digits}f}"


def _tajima_interp(d: Optional[float]) -> str:
    if d is None:
        return "n.a."
    if d > 2.0:
        return "Excess intermediate-frequency variants; balancing selection or contraction"
    if d < -2.0:
        return "Excess rare variants; selective sweep or population expansion"
    return "Consistent with neutrality"


def _fu_li_interp(v: Optional[float]) -> str:
    if v is None:
        return "n.a."
    if abs(v) < 2.0:
        return "Consistent with neutrality"
    return "Significant departure from neutrality"


def _r2_interp(v: Optional[float]) -> str:
    if v is None:
        return "n.a."
    return "Consistent with neutrality" if v >= 0.1 else "Low R2 suggests population expansion"


def _ld_significance(p: Optional[float]) -> str:
    if p is None:
        return ""
    if p < 0.001:
        return "***"
    if p < 0.01:
        return "**"
    if p < 0.05:
        return "*"
    return "n.s."


# ─────────────────────────────────────────────────────────────────────────────
# Output writers
# ─────────────────────────────────────────────────────────────────────────────

def write_tsv(
    output_dir: Path,
    source_name: str,
    global_stats: RegionStats,
    window_stats: list[RegionStats],
) -> Path:
    tsv_path = output_dir / "results.tsv"
    with open(tsv_path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["DnaSP-Python", "Source:", source_name, "Date:", datetime.now().strftime("%Y-%m-%d %H:%M")])
        w.writerow([])
        w.writerow(TSV_HEADER)
        w.writerow(global_stats.as_tsv_row())
        for ws in window_stats:
            w.writerow(ws.as_tsv_row())
    return tsv_path


def write_ld_tsv(output_dir: Path, ld: LDStats) -> Path:
    """Write LD pairwise results as TSV."""
    path = output_dir / "ld_pairs.tsv"
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh, delimiter="\t")
        w.writerow(["Site1", "Site2", "Dist", "n", "D", "D_prime", "R2", "Chi2", "P_chi2", "Sig"])
        for p in ld.pairs:
            w.writerow([
                p.site1, p.site2, p.dist, p.n_valid,
                _fmt(p.D, 6),
                _fmt(p.D_prime, 6),
                _fmt(p.R2, 6),
                _fmt(p.chi2, 4),
                _fmt(p.p_chi2, 6),
                _ld_significance(p.p_chi2),
            ])
    return path


def write_report(
    output_dir: Path,
    source_name: str,
    aln: Alignment,
    results: dict,
    figures: list[Path],
) -> Path:
    rs = results["global"]
    window_stats = results["windows"]
    ld: Optional[LDStats] = results.get("ld")
    recomb: Optional[RecombStats] = results.get("recombination")
    popsize: Optional[MismatchStats] = results.get("popsize")
    indel_s: Optional[InDelStats] = results.get("indel")
    div_s: Optional[DivergenceStats] = results.get("divergence")
    fuliout_s: Optional[FuLiOutgroupStats] = results.get("fuliout")
    hka_s: Optional[HKAStats] = results.get("hka")
    mk_s: Optional[MKStats] = results.get("mk")
    kaks_s: Optional[KaKsStats] = results.get("kaks")

    report_path = output_dir / "report.md"
    lines: list[str] = [
        "# DnaSP Population Genetics Report",
        "",
        f"**Input**: {source_name}  ",
        f"**Date**: {datetime.now().strftime('%Y-%m-%d %H:%M')}  ",
        f"**Tool**: DnaSP-Python (ClawBio skill)  ",
        "",
        "---",
        "",
        "## Input Summary",
        "",
        "| Parameter | Value |",
        "|-----------|-------|",
        f"| Number of sequences | {rs.n} |",
        f"| Alignment length (bp) | {rs.L_total} |",
        f"| Net sites (after gap removal) | {rs.L_net} |",
        "",
        "## Polymorphism Statistics",
        "",
        "| Statistic | Value |",
        "|-----------|-------|",
        f"| Segregating sites (S) | {rs.S} |",
        f"| Total mutations (η) | {rs.Eta} |",
        f"| Number of haplotypes (H) | {rs.H} |",
        f"| Haplotype diversity (Hd) | {_fmt(rs.Hd)} |",
        f"| Variance of Hd | {_fmt(rs.VarHd)} |",
        f"| Nucleotide diversity (π) | {_fmt(rs.Pi, 8)} |",
        f"| Mean pairwise differences (k) | {_fmt(rs.k, 6)} |",
        f"| G+C content | {_fmt(rs.GC, 4)} |",
        "",
        "## Theta Estimates",
        "",
        "| Estimator | Per site | Absolute |",
        "|-----------|---------|----------|",
        f"| Watterson's θ_W | {_fmt(rs.ThetaW_nuc, 8)} | {_fmt(rs.ThetaW, 6)} |",
        f"| Tajima's θ_k (= π) | {_fmt(rs.Pi, 8)} | {_fmt(rs.k, 6)} |",
        "",
        "## Neutrality Tests",
        "",
        "| Test | Value | Interpretation |",
        "|------|-------|----------------|",
        f"| Tajima's D | {_fmt(rs.TajimaD, 6)} | {_tajima_interp(rs.TajimaD)} |",
        f"| Fu & Li's D* | {_fmt(rs.FuLiD_star, 6)} | {_fu_li_interp(rs.FuLiD_star)} |",
        f"| Fu & Li's F* | {_fmt(rs.FuLiF_star, 6)} | {_fu_li_interp(rs.FuLiF_star)} |",
        f"| Ramos-Onsins & Rozas R2 | {_fmt(rs.R2, 6)} | {_r2_interp(rs.R2)} |",
        "",
    ]

    # ── Linkage Disequilibrium ───────────────────────────────────────────────
    if ld is not None:
        lines += [
            "## Linkage Disequilibrium",
            "",
            f"Biallelic sites analysed: {ld.n_biallelic}  ",
            f"Pairwise comparisons: {ld.n_pairs}  ",
            "",
            "| Statistic | Value | Reference |",
            "|-----------|-------|-----------|",
            f"| ZnS (average R²) | {_fmt(ld.ZnS, 6)} | Kelly 1997, eq 3 |",
            f"| Za (adjacent R²) | {_fmt(ld.Za, 6)} | Rozas et al. 2001 |",
            f"| ZZ = Za − ZnS | {_fmt(ld.ZZ, 6)} | Rozas et al. 2001 |",
            "",
        ]
        if ld.pairs:
            lines += [
                "### Pairwise LD (first 20 pairs shown; full table in ld_pairs.tsv)",
                "",
                "| Site1 | Site2 | Dist | D | D' | R² | χ² | P | Sig |",
                "|-------|-------|------|---|----|----|-----|---|-----|",
            ]
            for p in ld.pairs[:20]:
                lines.append(
                    f"| {p.site1} | {p.site2} | {p.dist} "
                    f"| {_fmt(p.D, 4)} | {_fmt(p.D_prime, 4)} | {_fmt(p.R2, 4)} "
                    f"| {_fmt(p.chi2, 3)} | {_fmt(p.p_chi2, 4)} | {_ld_significance(p.p_chi2)} |"
                )
            lines.append("")

    # ── Recombination ────────────────────────────────────────────────────────
    if recomb is not None:
        lines += [
            "## Recombination",
            "",
            "| Statistic | Value | Reference |",
            "|-----------|-------|-----------|",
            f"| Minimum recombination events (Rm) | {recomb.Rm} | Hudson & Kaplan 1985 |",
            f"| Incompatible site pairs | {recomb.n_incompatible_pairs} | Four-gamete test |",
            "",
        ]
        if recomb.incompatible_pairs and len(recomb.incompatible_pairs) <= 30:
            lines.append("**Incompatible pairs** (positions 1-based):")
            lines.append("")
            for p1, p2 in recomb.incompatible_pairs:
                lines.append(f"  - Sites {p1} and {p2}")
            lines.append("")

    # ── Population Size Changes ──────────────────────────────────────────────
    if popsize is not None:
        lines += [
            "## Population Size Changes  -  Mismatch Distribution",
            "",
            "| Statistic | Value | Reference |",
            "|-----------|-------|-----------|",
            f"| Number of pairs | {popsize.n_pairs} | |",
            f"| Mean pairwise differences | {_fmt(popsize.mean, 4)} | |",
            f"| Variance | {_fmt(popsize.variance, 4)} | |",
            f"| CV (coefficient of variation) | {_fmt(popsize.cv, 4)} | Rogers & Harpending 1992 |",
            f"| Raggedness r | {_fmt(popsize.raggedness, 6)} | Harpending 1994 |",
            "",
            "**Mismatch distribution** (differences → pair count):",
            "",
        ]
        if popsize.observed:
            lines.append("| Differences | Pairs | Proportion |")
            lines.append("|-------------|-------|------------|")
            for d, count in sorted(popsize.observed.items()):
                prop = count / popsize.n_pairs if popsize.n_pairs > 0 else 0
                lines.append(f"| {d} | {count} | {prop:.4f} |")
            lines.append("")
        r_note = (
            "Low raggedness → smooth distribution (consistent with population expansion).  \n"
            "High raggedness → ragged distribution (consistent with constant size or contraction)."
        )
        lines.append(r_note)
        lines.append("")

    # ── InDel Polymorphism ───────────────────────────────────────────────────
    if indel_s is not None:
        lines += [
            "## InDel Polymorphism",
            "",
            "| Statistic | Value |",
            "|-----------|-------|",
            f"| Alignment positions with gaps | {indel_s.n_positions_with_gaps} |",
            f"| Number of InDel events | {indel_s.n_events} |",
            f"| Mean InDel length (bp) | {_fmt(indel_s.mean_event_length, 2)} |",
            f"| InDel haplotypes | {indel_s.n_haplotypes} |",
            f"| InDel haplotype diversity (Hd) | {_fmt(indel_s.haplotype_diversity, 6)} |",
            f"| InDel diversity k(i) | {_fmt(indel_s.k_indel, 6)} |",
            f"| InDel diversity per site π(i) | {_fmt(indel_s.pi_indel, 6)} |",
            f"| θ from InDel events | {_fmt(indel_s.theta_indel, 6)} |",
            f"| Tajima's D (InDel) | {_fmt(indel_s.tajima_d_indel, 6)} |",
            "",
        ]

    # ── Divergence ───────────────────────────────────────────────────────────
    if div_s is not None:
        lines += [
            "## DNA Divergence Between Populations",
            "",
            f"**Population 1**: {div_s.pop1_name} (n = {div_s.n1})  ",
            f"**Population 2**: {div_s.pop2_name} (n = {div_s.n2})  ",
            f"**Net sites analysed**: {div_s.L_net}  ",
            "",
            "### Within-population diversity",
            "",
            "| Statistic | Pop1 | Pop2 |",
            "|-----------|------|------|",
            f"| π (nucleotide diversity) | {_fmt(div_s.Pi1, 6)} | {_fmt(div_s.Pi2, 6)} |",
            f"| k (mean pairwise diff) | {_fmt(div_s.k1, 4)} | {_fmt(div_s.k2, 4)} |",
            "",
            "### Between-population divergence",
            "",
            "| Statistic | Value | Reference |",
            "|-----------|-------|-----------|",
            f"| Dxy (nucleotide divergence) | {_fmt(div_s.Dxy, 6)} | Nei 1987, eq 10.20 |",
            f"| Da (net divergence) | {_fmt(div_s.Da, 6)} | Nei 1987 |",
            f"| Fixed differences | {div_s.n_fixed} | Hey 1991 |",
            f"| Shared polymorphisms | {div_s.n_shared} | |",
            f"| Private to {div_s.pop1_name} | {div_s.n_private1} | |",
            f"| Private to {div_s.pop2_name} | {div_s.n_private2} | |",
            "",
        ]

    # ── Fu & Li D/F with outgroup ─────────────────────────────────────────────
    if fuliout_s is not None:
        lines += [
            "## Fu & Li D and F (with Outgroup)",
            "",
            f"**Ingroup sequences**: {fuliout_s.n}  ",
            "",
            "| Statistic | Value | Reference |",
            "|-----------|-------|-----------|",
            f"| Total derived mutations (η) | {fuliout_s.eta} | Fu & Li 1993 |",
            f"| External mutations (η_e) | {fuliout_s.eta_e} | Fu & Li 1993 |",
            f"| Mean pairwise differences (k̄) | {_fmt(fuliout_s.k_bar, 4)} | |",
            f"| Fu & Li D | {_fmt(fuliout_s.D, 6)} | {_fu_li_interp(fuliout_s.D)} |",
            f"| Fu & Li F | {_fmt(fuliout_s.F, 6)} | {_fu_li_interp(fuliout_s.F)} |",
            "",
            "> **Note**: D and F with outgroup differ from D\\* and F\\* (no outgroup): "
            "the outgroup polarises mutations as ancestral/derived, enabling counting "
            "of 'external' mutations on terminal branches only (η_e).  "
            "Variance coefficients follow Simonsen et al. (1995) Appendix B.",
            "",
        ]

    # ── HKA test ─────────────────────────────────────────────────────────────
    if hka_s is not None and hka_s.n_loci >= 2:
        sig_str = " (significant)" if (hka_s.p_value is not None and hka_s.p_value < 0.05) else ""
        lines += [
            "## HKA Test (Hudson, Kreitman & Aguadé 1987)",
            "",
            f"**Loci analysed**: {hka_s.n_loci}  ",
            f"**Estimated divergence time** T̂ = {_fmt(hka_s.T_hat, 4)} (coalescent units)  ",
            "",
            "| Statistic | Value |",
            "|-----------|-------|",
            f"| χ² | {_fmt(hka_s.chi2, 4)} |",
            f"| Degrees of freedom | {hka_s.df} |",
            f"| P-value | {_fmt(hka_s.p_value, 6)}{sig_str} |",
            "",
        ]
        if hka_s.loci_results:
            lines += [
                "### Per-locus results",
                "",
                "| Locus | n | S | D | θ̂ | E[S] | E[D] |",
                "|-------|---|---|---|---|------|------|",
            ]
            for lr in hka_s.loci_results:
                lines.append(
                    f"| {lr['name']} | {lr['n']} | {lr['S']} | {lr['D']} "
                    f"| {_fmt(lr['theta_hat'], 4)} | {_fmt(lr['E_S'], 2)} "
                    f"| {_fmt(lr['E_D'], 2)} |"
                )
            lines.append("")
        lines += [
            "> **Interpretation**: A significant P-value (< 0.05) indicates that the "
            "ratio of polymorphism to divergence varies among loci, inconsistent with "
            "the neutral model. This can signal positive selection at specific loci.",
            "",
        ]

    # ── McDonald-Kreitman test ────────────────────────────────────────────────
    if mk_s is not None:
        sig_str = " (significant)" if (mk_s.fisher_p is not None and mk_s.fisher_p < 0.05) else ""
        lines += [
            "## McDonald-Kreitman Test",
            "",
            "| Category | Count |",
            "|----------|-------|",
            f"| Synonymous polymorphisms (Ps) | {mk_s.Ps} |",
            f"| Nonsynonymous polymorphisms (Pn) | {mk_s.Pn} |",
            f"| Synonymous fixed differences (Ds) | {mk_s.Ds} |",
            f"| Nonsynonymous fixed differences (Dn) | {mk_s.Dn} |",
            "",
            "| Statistic | Value | Reference |",
            "|-----------|-------|-----------|",
            f"| α (proportion adaptive) | {_fmt(mk_s.alpha, 4)} | McDonald & Kreitman 1991 |",
            f"| NI (Neutrality Index) | {_fmt(mk_s.NI, 4)} | Rand & Kann 1996 |",
            f"| DoS (Direction of Selection) | {_fmt(mk_s.DoS, 4)} | Stoletzki & Eyre-Walker 2011 |",
            f"| Fisher's exact P | {_fmt(mk_s.fisher_p, 6)}{sig_str} | |",
            "",
            "> **Interpretation**: α > 0 and P < 0.05 → positive (adaptive) selection driving "
            "nonsynonymous fixation. α < 0 → slight excess of nonsynonymous polymorphism "
            "(weak purifying selection or demographic effects). NI > 1 → purifying selection; "
            "NI < 1 → positive selection. DoS > 0 → net positive selection; DoS < 0 → "
            "net negative selection.",
            "",
        ]

    # ── Ka/Ks ─────────────────────────────────────────────────────────────────
    if kaks_s is not None and kaks_s.n_codons > 0:
        lines += [
            "## Ka/Ks (dN/dS)  -  Nei-Gojobori Method",
            "",
            "| Statistic | Value | Reference |",
            "|-----------|-------|-----------|",
            f"| Codons analysed | {kaks_s.n_codons} | |",
            f"| Mean synonymous sites (S) | {_fmt(kaks_s.S_sites, 2)} | Nei & Gojobori 1986 |",
            f"| Mean nonsynonymous sites (N) | {_fmt(kaks_s.N_sites, 2)} | |",
            f"| Mean synonymous differences (Sd) | {_fmt(kaks_s.Sd, 4)} | |",
            f"| Mean nonsynonymous differences (Nd) | {_fmt(kaks_s.Nd, 4)} | |",
            f"| Ks (synonymous substitutions/site) | {_fmt(kaks_s.Ks, 6)} | JC corrected |",
            f"| Ka (nonsynonymous substitutions/site) | {_fmt(kaks_s.Ka, 6)} | JC corrected |",
            f"| ω = Ka/Ks | {_fmt(kaks_s.omega, 4)} | |",
            "",
            "> **Interpretation**: ω < 1 → purifying (negative) selection constrains "
            "amino-acid change. ω ≈ 1 → neutral evolution. ω > 1 → positive (adaptive) "
            "selection driving amino-acid change. Values are averages over all pairwise "
            "comparisons; per-branch estimates require a phylogenetic framework.",
            "",
        ]

    # ── Fu's Fs ──────────────────────────────────────────────────────────────
    fufs_s: Optional[FuFsStats] = results.get("fufs")
    if fufs_s is not None:
        sig = ""
        if fufs_s.Fs is not None:
            if fufs_s.Fs < -1.61:    # S_k < 0.165 → conventional 0.02 threshold
                sig = " (**significant at 0.02 level**  -  fewer haplotypes than expected)"
        lines += [
            "## Fu's Fs Test (Fu 1997)",
            "",
            "| Statistic | Value |",
            "|-----------|-------|",
            f"| n (sequences) | {fufs_s.n} |",
            f"| H (observed haplotypes) | {fufs_s.H} |",
            f"| θ_π (= k, mean pairwise differences) | {_fmt(fufs_s.theta_pi, 4)} |",
            f"| S_k = P(K ≤ H \\| θ_π, n) | {_fmt(fufs_s.S_k, 6)} |",
            f"| Fs | {_fmt(fufs_s.Fs, 4)} |",
            "",
            f"> **Interpretation**: Fs = {_fmt(fufs_s.Fs, 4)}.{sig} "
            "Large negative Fs indicates fewer haplotypes than expected given nucleotide "
            "diversity  -  signature of recent population expansion or positive selection. "
            "Use the conventional significance threshold Fs < 0 with S_k ≤ 0.02. "
            "Fs > 0 (excess haplotypes) is rarely significant and usually arises from "
            "balancing selection or population subdivision.",
            "",
        ]

    # ── Site Frequency Spectrum ───────────────────────────────────────────────
    sfs_s: Optional[SFSStats] = results.get("sfs")
    if sfs_s is not None and sfs_s.folded:
        lines += [
            "## Site Frequency Spectrum",
            "",
        ]
        # Folded SFS
        lines += [
            "### Folded SFS (minor allele frequency)",
            "",
            "| Minor allele count (i) | Frequency (i/n) | Sites |",
            "|------------------------|-----------------|-------|",
        ]
        for i in sorted(sfs_s.folded.keys()):
            freq = i / sfs_s.n if sfs_s.n > 0 else 0.0
            lines.append(f"| {i} | {freq:.3f} | {sfs_s.folded[i]} |")
        lines.append("")
        # Unfolded SFS
        if sfs_s.unfolded is not None and sfs_s.unfolded:
            lines += [
                "### Unfolded SFS (derived allele frequency, outgroup-polarised)",
                "",
                "| Derived allele count (i) | Frequency (i/n) | Sites |",
                "|--------------------------|-----------------|-------|",
            ]
            for i in sorted(sfs_s.unfolded.keys()):
                freq = i / sfs_s.n if sfs_s.n > 0 else 0.0
                lines.append(f"| {i} | {freq:.3f} | {sfs_s.unfolded[i]} |")
            lines.append("")
        elif sfs_s.has_outgroup:
            lines += ["*No unfolded SFS sites  -  outgroup did not polarise any segregating sites.*", ""]
        else:
            lines += ["*Unfolded SFS not available  -  no outgroup provided. Rerun with --outgroup <seq_name>.*", ""]
        lines += [
            "> **Interpretation**: A skew toward singletons (i = 1) → excess rare variants → "
            "population expansion or purifying selection (consistent with negative Tajima's D). "
            "A uniform or U-shaped spectrum → balancing selection. "
            "The unfolded SFS distinguishes ancestral from derived alleles; the folded SFS folds "
            "the spectrum at i = n/2, treating both polarities equivalently.",
            "",
        ]

    # ── Transition / Transversion ratio ──────────────────────────────────────
    tstv_s: Optional[TsTvStats] = results.get("tstv")
    if tstv_s is not None:
        ts_tv_str = _fmt(tstv_s.ts_tv, 4) if tstv_s.ts_tv is not None else "N/A (Tv = 0)"
        lines += [
            "## Transition / Transversion Ratio",
            "",
            f"| Statistic | Value |",
            f"|-----------|-------|",
            f"| Sequences (n) | {tstv_s.n} |",
            f"| Clean sites (L) | {tstv_s.L_net} |",
            f"| Total transitions (Ts) | {tstv_s.n_transitions} |",
            f"| Total transversions (Tv) | {tstv_s.n_transversions} |",
            f"| Ts/Tv ratio | {ts_tv_str} |",
            f"| Mean Ts per pair per site | {_fmt(tstv_s.ts_per_site, 5)} |",
            f"| Mean Tv per pair per site | {_fmt(tstv_s.tv_per_site, 5)} |",
            "",
            "> **Interpretation**: Transitions (purine↔purine: A↔G; pyrimidine↔pyrimidine: C↔T) "
            "are expected to outnumber transversions due to mutational bias. "
            "Ts/Tv ≈ 2 is typical for nuclear DNA; mitochondrial DNA often shows Ts/Tv > 10. "
            "Ts/Tv < 0.5 may indicate saturation or non-neutral evolution. "
            "'N/A' is reported when no transversions are observed.",
            "",
        ]

    # ── Codon usage bias ─────────────────────────────────────────────────────
    codon_s: Optional[CodonUsageStats] = results.get("codon")
    if codon_s is not None:
        enc_str = _fmt(codon_s.ENC, 2) if codon_s.ENC is not None else "N/A (insufficient data)"
        lines += [
            "## Codon Usage Bias",
            "",
            f"| Statistic | Value |",
            f"|-----------|-------|",
            f"| Sequences (n) | {codon_s.n} |",
            f"| Mean coding codons per sequence | {_fmt(codon_s.n_codons, 1)} |",
            f"| Effective number of codons (ENC) | {enc_str} |",
            "",
        ]
        if codon_s.ENC is not None:
            if codon_s.ENC < 35:
                bias_note = "**Strong codon usage bias** (ENC < 35)."
            elif codon_s.ENC < 50:
                bias_note = "**Moderate codon usage bias** (35 ≤ ENC < 50)."
            else:
                bias_note = "**Weak or no codon usage bias** (ENC ≥ 50; close to 61)."
            lines += [f"> {bias_note}", ""]

        # RSCU table  -  group by amino acid family
        if codon_s.rscu:
            lines += [
                "### RSCU Values",
                "",
                "| Amino acid | Codons | RSCU |",
                "|------------|--------|------|",
            ]
            for aa in sorted(_SYNONYMOUS_FAMILIES.keys()):
                codons = sorted(_SYNONYMOUS_FAMILIES[aa])
                for c in codons:
                    rscu_val = codon_s.rscu.get(c)
                    rscu_str = _fmt(rscu_val, 3) if rscu_val is not None else " - "
                    lines.append(f"| {aa} | {c} | {rscu_str} |")
            lines.append("")
            lines += [
                "> RSCU = 1.0: uniform usage among synonyms. "
                "RSCU > 1.0: preferred codon. RSCU < 1.0: avoided codon.",
                "",
            ]

    # ── Fay & Wu ─────────────────────────────────────────────────────────────
    faywu_s: Optional[FayWuStats] = results.get("faywu")
    if faywu_s is not None:
        def _fw(v: Optional[float]) -> str:
            return "N/A" if v is None else _fmt(v, 6)
        lines += [
            "## Fay & Wu's H and Zeng's E",
            "",
            f"| Statistic | Value |",
            f"|-----------|-------|",
            f"| Ingroup sequences (n) | {faywu_s.n} |",
            f"| Sites surviving complete deletion (L_net) | {faywu_s.L_net} |",
            f"| Polarisable segregating sites | {faywu_s.n_polarised} |",
            f"| θ_π (from polarised sites) | {_fw(faywu_s.theta_pi)} |",
            f"| θ_W (Watterson, polarised) | {_fw(faywu_s.theta_w)} |",
            f"| θ_H (Fay & Wu 2000) | {_fw(faywu_s.theta_h)} |",
            f"| θ_L (Zeng et al. 2006) | {_fw(faywu_s.theta_l)} |",
            f"| **H = θ_π − θ_H** | {_fw(faywu_s.H)} |",
            f"| **E = θ_L − θ_W** | {_fw(faywu_s.E)} |",
            "",
        ]
        if faywu_s.H is not None:
            if faywu_s.H < 0:
                h_note = ("**H < 0**: excess of high-frequency derived alleles  -  "
                          "consistent with a selective sweep or background selection.")
            else:
                h_note = ("**H ≥ 0**: no excess of high-frequency derived alleles; "
                          "consistent with neutrality or purifying selection.")
            lines += [f"> {h_note}", ""]
        if faywu_s.E is not None:
            if faywu_s.E < 0:
                e_note = ("**E < 0**: θ_L < θ_W  -  excess of low-frequency derived alleles; "
                          "consistent with a recent bottleneck or purifying selection.")
            else:
                e_note = ("**E ≥ 0**: θ_L ≥ θ_W  -  no excess of low-frequency derived alleles.")
            lines += [f"> {e_note}", ""]

    # ── Fst ──────────────────────────────────────────────────────────────────
    fst_s: Optional[FstStats] = results.get("fst")
    if fst_s is not None and fst_s.n_pops >= 2:
        lines += [
            "## Population Differentiation (Fst)",
            "",
            f"| Population | n | π (within-pop) |",
            f"|------------|---|----------------|",
        ]
        for pop in fst_s.pop_names:
            lines.append(
                f"| {pop} | {fst_s.pop_sizes.get(pop, '?')} "
                f"| {_fmt(fst_s.pi_within.get(pop, 0.0), 6)} |"
            )
        lines.append("")
        lines += [
            "### Pairwise Fst (Hudson et al. 1992)",
            "",
            "| Population pair | Dxy (between) | Fst |",
            "|-----------------|---------------|-----|",
        ]
        for (p1, p2), fst_val in fst_s.fst_pairwise.items():
            dxy = fst_s.pi_between.get((p1, p2), 0.0)
            fst_str = _fmt(fst_val, 4) if fst_val is not None else "N/A (Dxy=0)"
            lines.append(f"| {p1} vs {p2} | {_fmt(dxy, 6)} | {fst_str} |")
        lines.append("")
        if fst_s.fst_mean is not None:
            lines += [f"> **Mean Fst across all pairs: {_fmt(fst_s.fst_mean, 4)}**", ""]
            if fst_s.fst_mean < 0.05:
                lines += ["> Little genetic differentiation (Fst < 0.05; Wright 1978).", ""]
            elif fst_s.fst_mean < 0.15:
                lines += ["> Moderate genetic differentiation (0.05 ≤ Fst < 0.15; Wright 1978).", ""]
            elif fst_s.fst_mean < 0.25:
                lines += ["> Great genetic differentiation (0.15 ≤ Fst < 0.25; Wright 1978).", ""]
            else:
                lines += ["> Very great genetic differentiation (Fst ≥ 0.25; Wright 1978).", ""]

    # ── Sliding window ───────────────────────────────────────────────────────
    if window_stats:
        lines += [
            "## Sliding Window Summary",
            "",
            "| Region | S | π | Tajima D |",
            "|--------|---|---|---------|",
        ]
        for ws in window_stats:
            lines.append(
                f"| {ws.region} | {ws.S} | {_fmt(ws.Pi, 5)} | {_fmt(ws.TajimaD, 4)} |"
            )
        lines.append("")

    # ── Figures ──────────────────────────────────────────────────────────────
    if figures:
        lines += ["## Figures", ""]
        for fig in figures:
            lines.append(f"![{fig.stem}]({fig.name})")
        lines.append("")

    # ── Methods & References ─────────────────────────────────────────────────
    lines += [
        "---",
        "",
        "## Methods",
        "",
        "Statistics computed using DnaSP-Python (ClawBio), reimplementing "
        "core DnaSP 6 algorithms (Rozas et al. 2017). "
        "Gap treatment: complete deletion. Neutrality tests include both intraspecific "
        "(D\\*, F\\*  -  no outgroup) and outgroup-polarised (D, F) formulations.",
        "",
        "## References",
        "",
        "- Rozas et al. (2017) J. Hered. 108:591-593  -  DnaSP v6",
        "- Tajima (1989) Genetics 123:585-595  -  Tajima's D",
        "- Fu & Li (1993) Genetics 133:693-709  -  D, F (outgroup) and D*, F* (no outgroup)",
        "- Simonsen et al. (1995) Genetics 141:413-429  -  variance coefficients",
        "- Nei & Tajima (1981) Genetics 97:145-163  -  haplotype diversity",
        "- Ramos-Onsins & Rozas (2002) Mol. Biol. Evol. 19:2092-2100  -  R2",
        "- Kelly (1997) Genetics 146:1197-1206  -  ZnS",
        "- Rozas et al. (2001) Genetics 158:1147-1155  -  Za, ZZ",
        "- Hill & Robertson (1968) Theor. Appl. Genet. 38:226-231  -  R²",
        "- Lewontin (1964) Genetics 49:49-67  -  D'",
        "- Hudson & Kaplan (1985) Genetics 111:147-164  -  Rm",
        "- Rogers & Harpending (1992) Mol. Biol. Evol. 9:552-569  -  mismatch",
        "- Harpending (1994) Hum. Biol. 66:591-600  -  raggedness",
        "- Nei (1987) Molecular Evolutionary Genetics  -  Dxy, Da",
        "- Hudson, Kreitman & Aguadé (1987) Genetics 116:153-159  -  HKA test",
        "- McDonald & Kreitman (1991) Nature 351:652-654  -  MK test",
        "- Nei & Gojobori (1986) Mol. Biol. Evol. 3:418-426  -  Ka/Ks",
        "- Fu (1997) Genetics 147:915-925  -  Fu's Fs",
        "- Sharp & Li (1987) Nucleic Acids Res. 15:1281-1295  -  RSCU",
        "- Wright (1990) Gene 87:23-29  -  ENC",
        "- Fay & Wu (2000) Genetics 155:1405-1413  -  H statistic",
        "- Zeng et al. (2006) Genetics 174:1431-1439  -  E statistic",
        "- Hudson et al. (1992) Genetics 132:583-589  -  Fst estimator",
        "",
        "---",
        "",
        "*ClawBio is a research and educational tool. It is not a medical device "
        "and does not provide clinical diagnoses.*",
    ]

    report_path.write_text("\n".join(lines), encoding="utf-8")
    return report_path


# ─────────────────────────────────────────────────────────────────────────────
# Figures
# ─────────────────────────────────────────────────────────────────────────────

def make_figures(output_dir: Path, results: dict) -> list[Path]:
    figs_dir = output_dir / "figures"
    figs_dir.mkdir(exist_ok=True)
    paths: list[Path] = []

    if not HAS_MPL:
        return paths

    global_stats = results["global"]
    window_stats = results["windows"]
    ld: Optional[LDStats] = results.get("ld")
    popsize: Optional[MismatchStats] = results.get("popsize")

    # Sliding window: π and Tajima's D
    if window_stats:
        regions = [ws.region for ws in window_stats]
        midpoints = []
        for r in regions:
            parts = r.split("-")
            try:
                midpoints.append((int(parts[0]) + int(parts[1])) / 2)
            except (IndexError, ValueError):
                midpoints.append(0)
        pi_vals = [ws.Pi for ws in window_stats]
        d_vals = [ws.TajimaD if ws.TajimaD is not None else float("nan") for ws in window_stats]
        fig, axes = plt.subplots(2, 1, figsize=(10, 6), sharex=True)
        axes[0].plot(midpoints, pi_vals, color="#2166ac", linewidth=1.5)
        axes[0].set_ylabel("Nucleotide diversity (π)")
        axes[0].set_title("Sliding Window Analysis")
        axes[1].plot(midpoints, d_vals, color="#d6604d", linewidth=1.5)
        axes[1].axhline(0, color="grey", linewidth=0.8, linestyle="--")
        axes[1].set_ylabel("Tajima's D")
        axes[1].set_xlabel("Position (bp)")
        plt.tight_layout()
        fig_path = figs_dir / "sliding_window.png"
        plt.savefig(fig_path, dpi=150)
        plt.close()
        paths.append(fig_path)

    # Summary bar chart
    labels = ["π", "θ_W", "Hd"]
    values = [global_stats.Pi, global_stats.ThetaW_nuc, global_stats.Hd]
    colors = ["#2166ac", "#4dac26", "#d6604d"]
    fig, ax = plt.subplots(figsize=(5, 4))
    bars = ax.bar(labels, values, color=colors, edgecolor="white", linewidth=1.2)
    ax.set_ylabel("Value")
    ax.set_title("Population Genetics Summary")
    for bar, val in zip(bars, values):
        ax.text(
            bar.get_x() + bar.get_width() / 2,
            bar.get_height() + max(values) * 0.01 if max(values) > 0 else 0.001,
            f"{val:.4f}", ha="center", va="bottom", fontsize=9,
        )
    plt.tight_layout()
    fig_path = figs_dir / "summary.png"
    plt.savefig(fig_path, dpi=150)
    plt.close()
    paths.append(fig_path)

    # LD scatter: R² vs distance
    if ld is not None and ld.pairs:
        dists = [p.dist for p in ld.pairs]
        r2s = [p.R2 if p.R2 is not None else 0 for p in ld.pairs]
        fig, ax = plt.subplots(figsize=(6, 4))
        ax.scatter(dists, r2s, alpha=0.5, s=20, color="#762a83")
        ax.set_xlabel("Nucleotide distance (bp)")
        ax.set_ylabel("R²")
        ax.set_title("LD Decay")
        ax.set_ylim(0, 1.05)
        plt.tight_layout()
        fig_path = figs_dir / "ld_decay.png"
        plt.savefig(fig_path, dpi=150)
        plt.close()
        paths.append(fig_path)

    # Mismatch distribution
    if popsize is not None and popsize.observed:
        diffs = sorted(popsize.observed.keys())
        counts = [popsize.observed[d] for d in diffs]
        total = sum(counts)
        props = [c / total for c in counts]
        fig, ax = plt.subplots(figsize=(7, 4))
        ax.bar(diffs, props, color="#1a9850", edgecolor="white", linewidth=0.8)
        ax.set_xlabel("Pairwise nucleotide differences")
        ax.set_ylabel("Proportion of pairs")
        ax.set_title(f"Mismatch Distribution  (r = {_fmt(popsize.raggedness, 4)})")
        plt.tight_layout()
        fig_path = figs_dir / "mismatch.png"
        plt.savefig(fig_path, dpi=150)
        plt.close()
        paths.append(fig_path)

    # SFS bar chart (folded; add unfolded as second panel if available)
    sfs: Optional[SFSStats] = results.get("sfs")
    if sfs is not None and sfs.folded:
        n = sfs.n
        # Determine panels
        has_unfolded = sfs.unfolded is not None and bool(sfs.unfolded)
        n_panels = 2 if has_unfolded else 1
        fig, axes = plt.subplots(1, n_panels, figsize=(5 * n_panels, 4))
        if n_panels == 1:
            axes = [axes]

        # Folded SFS
        max_i = n // 2
        folded_x = list(range(1, max_i + 1))
        folded_y = [sfs.folded.get(i, 0) for i in folded_x]
        axes[0].bar(folded_x, folded_y, color="#4575b4", edgecolor="white", linewidth=0.8)
        axes[0].set_xlabel("Minor allele count (i)")
        axes[0].set_ylabel("Number of sites")
        axes[0].set_title("Folded SFS")
        axes[0].set_xticks(folded_x)

        # Unfolded SFS
        if has_unfolded:
            unfolded_x = list(range(1, n))
            unfolded_y = [sfs.unfolded.get(i, 0) for i in unfolded_x]  # type: ignore[union-attr]
            axes[1].bar(unfolded_x, unfolded_y, color="#d73027", edgecolor="white", linewidth=0.8)
            axes[1].set_xlabel("Derived allele count (i)")
            axes[1].set_ylabel("Number of sites")
            axes[1].set_title("Unfolded SFS (outgroup-polarised)")
            axes[1].set_xticks(unfolded_x)

        plt.suptitle("Site Frequency Spectrum", y=1.01)
        plt.tight_layout()
        fig_path = figs_dir / "sfs.png"
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(fig_path)

    # ── Codon usage RSCU bar chart ────────────────────────────────────────────
    codon_s = results.get("codon")
    if codon_s is not None and codon_s.rscu:
        # Build ordered list: group by amino acid family, sorted by aa
        bar_labels: list[str] = []
        bar_values: list[float] = []
        bar_colors: list[str] = []
        # Cycle colours across amino acid families for visual grouping
        palette = ["#4575b4", "#d73027", "#1a9850", "#fdae61",
                   "#74add1", "#f46d43", "#66bd63", "#fee090",
                   "#313695", "#a50026", "#006837", "#ffffbf"]
        for idx, aa in enumerate(sorted(_SYNONYMOUS_FAMILIES.keys())):
            codons = sorted(_SYNONYMOUS_FAMILIES[aa])
            col = palette[idx % len(palette)]
            for c in codons:
                bar_labels.append(f"{c}\n({aa})")
                bar_values.append(codon_s.rscu.get(c, 0.0))
                bar_colors.append(col)

        fig, ax = plt.subplots(figsize=(max(12, len(bar_labels) * 0.35), 5))
        x_pos = list(range(len(bar_labels)))
        ax.bar(x_pos, bar_values, color=bar_colors, edgecolor="white", linewidth=0.5)
        ax.axhline(1.0, color="black", linewidth=0.8, linestyle="--", label="RSCU = 1 (uniform)")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(bar_labels, fontsize=6, rotation=90)
        ax.set_ylabel("RSCU")
        enc_label = f"  ENC = {codon_s.ENC:.1f}" if codon_s.ENC is not None else ""
        ax.set_title(f"Relative Synonymous Codon Usage (RSCU){enc_label}")
        ax.legend(fontsize=8)
        plt.tight_layout()
        fig_path = figs_dir / "codon_usage.png"
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(fig_path)

    # ── Fst bar chart ────────────────────────────────────────────────────────
    fst_s: Optional[FstStats] = results.get("fst")
    if HAS_MPL and fst_s is not None and fst_s.fst_pairwise:
        pairs = [f"{p1}\nvs\n{p2}" for (p1, p2) in fst_s.fst_pairwise]
        fst_vals = [v if v is not None else 0.0 for v in fst_s.fst_pairwise.values()]
        fig, ax = plt.subplots(figsize=(max(4, len(pairs) * 1.2), 4))
        colours = ["#E05C5C" if v >= 0.25 else "#F5A623" if v >= 0.15 else
                   "#7ED321" if v >= 0.05 else "#4A90D9" for v in fst_vals]
        bars = ax.bar(range(len(pairs)), fst_vals, color=colours, edgecolor="white", linewidth=0.5)
        ax.set_xticks(range(len(pairs)))
        ax.set_xticklabels(pairs, fontsize=8)
        ax.set_ylim(0, max(1.0, max(fst_vals) * 1.15))
        ax.axhline(0.05, color="steelblue", lw=0.8, ls="--", label="Little (0.05)")
        ax.axhline(0.15, color="goldenrod", lw=0.8, ls="--", label="Moderate (0.15)")
        ax.axhline(0.25, color="tomato",    lw=0.8, ls="--", label="Great (0.25)")
        ax.set_ylabel("Fst")
        ax.set_title("Population Differentiation (Fst)  -  Hudson et al. 1992")
        ax.legend(fontsize=7, loc="upper right")
        for bar, val in zip(bars, fst_vals):
            ax.text(bar.get_x() + bar.get_width() / 2, bar.get_height() + 0.01,
                    f"{val:.3f}", ha="center", va="bottom", fontsize=8)
        plt.tight_layout()
        fig_path = figs_dir / "fst.png"
        plt.savefig(fig_path, dpi=150, bbox_inches="tight")
        plt.close()
        paths.append(fig_path)

    return paths


# ─────────────────────────────────────────────────────────────────────────────
# Reproducibility bundle
# ─────────────────────────────────────────────────────────────────────────────

def write_reproducibility(
    output_dir: Path,
    input_path: Optional[Path],
    cli_args: list[str],
    result_files: list[Path],
) -> None:
    repro_dir = output_dir / "reproducibility"
    repro_dir.mkdir(exist_ok=True)
    cmd = " ".join(["python", "skills/dnasp/dnasp.py"] + cli_args)
    (repro_dir / "commands.sh").write_text(
        f"#!/bin/bash\n# DnaSP-Python  -  exact reproduction command\n{cmd}\n",
        encoding="utf-8",
    )
    env_src = Path(__file__).parent / "environment.yml"
    if env_src.exists():
        shutil.copy(env_src, repro_dir / "environment.yml")
    lines = []
    all_files = ([input_path] if input_path else []) + result_files
    for fp in all_files:
        if fp and fp.exists():
            digest = hashlib.sha256(fp.read_bytes()).hexdigest()
            lines.append(f"{digest}  {fp.name}")
    (repro_dir / "checksums.sha256").write_text("\n".join(lines) + "\n", encoding="utf-8")


# ─────────────────────────────────────────────────────────────────────────────
# Demo data
# ─────────────────────────────────────────────────────────────────────────────

DEMO_FASTA = """\
>pop1_seq1
TTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop1_seq2
TTTTATCATGATAATAAAGAACAATGTATTATCATAGCCGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop1_seq3
TTCTATCATGATAATAAAGAACAATGTATTATCATAGCCGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop1_seq4
TTCTATCATAATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop1_seq5
TTTTATCATAATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop2_seq1
TTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop2_seq2
TTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTCTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop2_seq3
TTCTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTCTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop2_seq4
TTTTATCATGATAATAAAGAAAAATGTATTATCATAGCCGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>pop2_seq5
TTTTATCATGATAATAAAGAAAAATGTATTATCATAGCCGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
>outgroup
TTTTATCATGATAATAAAGAACAATGTATTATCATAACTGGTCCCACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTCGTATGTGGTTTTATCATGATAATAAAGAACAATGTATTATCATAGCTGGTCCTACTGTTTTGCTTTCTGGTATGTGGGCTGGTACTGTTTTTAAAGAATCT
"""

DEMO_POP_FILE = """\
pop1_seq1\tPop1
pop1_seq2\tPop1
pop1_seq3\tPop1
pop1_seq4\tPop1
pop1_seq5\tPop1
pop2_seq1\tPop2
pop2_seq2\tPop2
pop2_seq3\tPop2
pop2_seq4\tPop2
pop2_seq5\tPop2
"""

DEMO_POP_ASSIGNMENTS: dict[str, str] = {
    "pop1_seq1": "Pop1", "pop1_seq2": "Pop1", "pop1_seq3": "Pop1",
    "pop1_seq4": "Pop1", "pop1_seq5": "Pop1",
    "pop2_seq1": "Pop2", "pop2_seq2": "Pop2", "pop2_seq3": "Pop2",
    "pop2_seq4": "Pop2", "pop2_seq5": "Pop2",
}

DEMO_DESCRIPTION = """\
Demo alignment: 10 ingroup sequences + 1 outgroup × 300 bp (2 populations, in-frame CDS)
  Pop1: pop1_seq1-5  |  Pop2: pop2_seq1-5  |  Outgroup: outgroup
  Segregating sites S=5, haplotypes H=8, Hd≈0.9556, Tajima's D≈0.6789
  Ts=77, Tv=16, Ts/Tv≈4.81  |  ENC≈23.00 (strong codon-usage bias)
  MK: Pn=2, Ps=3, Dn=1, Ds=1, NI≈0.667, α≈0.333
  KaKs: Ka≈0.00298, Ks≈0.02281, ω≈0.131
"""


def run_demo(output_dir: Path, analyses: set[str]) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)
    demo_path = output_dir / "demo_input.fas"
    demo_path.write_text(DEMO_FASTA, encoding="utf-8")
    pop_path = output_dir / "demo_pop.txt"
    pop_path.write_text(DEMO_POP_FILE, encoding="utf-8")
    print(DEMO_DESCRIPTION)
    print(f"Demo input written to: {demo_path}")
    print(f"Demo pop file written to: {pop_path}")
    _run(demo_path, output_dir, 0, 0, analyses,
         DEMO_POP_ASSIGNMENTS, None, ["--demo"],
         outgroup_name="outgroup")


# ─────────────────────────────────────────────────────────────────────────────
# Main pipeline
# ─────────────────────────────────────────────────────────────────────────────

def _run(
    input_path: Optional[Path],
    output_dir: Path,
    window_size: int,
    step_size: int,
    analyses: set[str],
    pop_assignments: Optional[dict[str, str]],
    aln2: Optional[Alignment],
    cli_args: list[str],
    outgroup_name: Optional[str] = None,
    hka_loci: Optional[list[HKALocus]] = None,
) -> None:
    output_dir.mkdir(parents=True, exist_ok=True)

    if input_path:
        print(f"Loading alignment: {input_path}")
        aln = load_alignment(input_path)
    else:
        # Demo mode: parse inline
        tmp = output_dir / "demo_input.fas"
        aln = load_alignment(tmp)

    print(f"  {aln.n} sequences, {aln.L} sites")

    # Extract outgroup sequence from alignment (removes it from ingroup)
    outgroup_seq: Optional[str] = None
    if outgroup_name is not None:
        if outgroup_name in aln.names:
            idx = aln.names.index(outgroup_name)
            outgroup_seq = aln.seqs[idx]
            # Remove outgroup from ingroup sequences
            aln = Alignment(
                names=[n for i, n in enumerate(aln.names) if i != idx],
                seqs=[s for i, s in enumerate(aln.seqs) if i != idx],
                source=aln.source,
            )
            print(f"  Outgroup '{outgroup_name}' extracted; {aln.n} ingroup sequences remain.")
        else:
            print(
                f"Warning: outgroup sequence '{outgroup_name}' not found in alignment. "
                "fuliout analysis skipped.",
                file=sys.stderr,
            )

    active = ", ".join(sorted(analyses))
    print(f"Running analyses: {active}")

    results = run_analysis(
        aln, window_size, step_size, analyses, pop_assignments, aln2,
        outgroup=outgroup_seq, hka_loci=hka_loci,
    )

    rs = results["global"]
    print(
        f"\n{'─'*50}\n"
        f"  n={rs.n}  L={rs.L_total}  NetSites={rs.L_net}\n"
        f"  S={rs.S}  η={rs.Eta}  H={rs.H}\n"
        f"  Hd={_fmt(rs.Hd,4)}  π={_fmt(rs.Pi,6)}  k={_fmt(rs.k,4)}\n"
        f"  θ_W={_fmt(rs.ThetaW_nuc,6)}  Tajima D={_fmt(rs.TajimaD,4)}\n"
        f"  Fu&Li D*={_fmt(rs.FuLiD_star,4)}  F*={_fmt(rs.FuLiF_star,4)}\n"
        f"  R2={_fmt(rs.R2,6)}\n"
        f"{'─'*50}"
    )

    ld = results.get("ld")
    if ld is not None:
        print(f"  ZnS={_fmt(ld.ZnS,6)}  Za={_fmt(ld.Za,6)}  ZZ={_fmt(ld.ZZ,6)}")

    recomb = results.get("recombination")
    if recomb is not None:
        print(f"  Rm={recomb.Rm}  incompatible_pairs={recomb.n_incompatible_pairs}")

    mism = results.get("popsize")
    if mism is not None:
        print(f"  raggedness_r={_fmt(mism.raggedness,6)}  mean_diff={_fmt(mism.mean,4)}")

    fuliout = results.get("fuliout")
    if fuliout is not None:
        print(f"  Fu&Li D={_fmt(fuliout.D,4)}  F={_fmt(fuliout.F,4)}  eta={fuliout.eta}  eta_e={fuliout.eta_e}")

    hka = results.get("hka")
    if hka is not None and hka.n_loci >= 2:
        print(f"  HKA chi2={_fmt(hka.chi2,4)}  df={hka.df}  p={_fmt(hka.p_value,6)}  T_hat={_fmt(hka.T_hat,4)}")

    mk = results.get("mk")
    if mk is not None:
        print(f"  MK  Pn={mk.Pn}  Ps={mk.Ps}  Dn={mk.Dn}  Ds={mk.Ds}  "
              f"alpha={_fmt(mk.alpha,4)}  DoS={_fmt(mk.DoS,4)}  p={_fmt(mk.fisher_p,6)}")

    kaks = results.get("kaks")
    if kaks is not None and kaks.n_codons > 0:
        print(f"  Ka={_fmt(kaks.Ka,6)}  Ks={_fmt(kaks.Ks,6)}  omega={_fmt(kaks.omega,4)}")

    fufs = results.get("fufs")
    if fufs is not None and fufs.Fs is not None:
        print(f"  Fu's Fs={_fmt(fufs.Fs,4)}  S_k={_fmt(fufs.S_k,6)}  H={fufs.H}  theta_pi={_fmt(fufs.theta_pi,4)}")

    sfs = results.get("sfs")
    if sfs is not None and sfs.folded:
        total_seg = sum(sfs.folded.values())
        singleton_frac = sfs.folded.get(1, 0) / total_seg if total_seg > 0 else 0.0
        unfolded_note = f"  unfolded_sites={sum(sfs.unfolded.values())}" if sfs.unfolded else ""
        print(f"  SFS folded_sites={total_seg}  singletons={sfs.folded.get(1,0)}  "
              f"singleton_frac={singleton_frac:.3f}{unfolded_note}")

    tstv = results.get("tstv")
    if tstv is not None:
        ts_tv_str = f"{tstv.ts_tv:.4f}" if tstv.ts_tv is not None else "N/A"
        print(f"  TsTv Ts={tstv.n_transitions}  Tv={tstv.n_transversions}  Ts/Tv={ts_tv_str}"
              f"  L_net={tstv.L_net}")

    codon_r = results.get("codon")
    if codon_r is not None:
        enc_str = f"{codon_r.ENC:.2f}" if codon_r.ENC is not None else "N/A"
        print(f"  Codon n_codons={codon_r.n_codons:.1f}  ENC={enc_str}")

    faywu_r = results.get("faywu")
    if faywu_r is not None:
        h_str = f"{faywu_r.H:.6f}" if faywu_r.H is not None else "N/A"
        e_str = f"{faywu_r.E:.6f}" if faywu_r.E is not None else "N/A"
        print(f"  FayWu H={h_str}  E={e_str}  n_polarised={faywu_r.n_polarised}")

    fst_r = results.get("fst")
    if fst_r is not None and fst_r.fst_pairwise:
        fst_mean_str = f"{fst_r.fst_mean:.4f}" if fst_r.fst_mean is not None else "N/A"
        pair_strs = "  ".join(
            f"{p1}/{p2}={f'{v:.4f}' if v is not None else 'N/A'}"
            for (p1, p2), v in fst_r.fst_pairwise.items()
        )
        print(f"  Fst mean={fst_mean_str}  pairs: {pair_strs}")

    print(f"\nWriting output to: {output_dir}")

    tsv = write_tsv(output_dir, input_path.name if input_path else "demo", rs, results["windows"])

    result_files = [tsv]
    if ld is not None and ld.pairs:
        ld_tsv = write_ld_tsv(output_dir, ld)
        result_files.append(ld_tsv)

    figs = make_figures(output_dir, results)
    if not HAS_MPL:
        print("  Note: matplotlib not installed  -  figures skipped.")

    report = write_report(output_dir, input_path.name if input_path else "demo", aln, results, figs)
    result_files.append(report)

    write_reproducibility(output_dir, input_path, cli_args, result_files)

    print(f"  Report:  {report}")
    print(f"  TSV:     {tsv}")
    if figs:
        print(f"  Figures: {', '.join(str(f) for f in figs)}")
    print("Done.")


# ─────────────────────────────────────────────────────────────────────────────
# CLI
# ─────────────────────────────────────────────────────────────────────────────

def build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        description="DnaSP-Python  -  population genetics analysis of DNA alignments",
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog=(
            "Analyses available (--analysis):\n"
            "  polymorphism  π, k, S, Hd, θ_W, Tajima D, Fu & Li D*/F*, R2  [always run]\n"
            "  ld            Linkage disequilibrium (D, D', R², ZnS, Za, ZZ)\n"
            "  recombination Minimum recombination events Rm (four-gamete test)\n"
            "  popsize       Mismatch distribution, raggedness r\n"
            "  indel         InDel polymorphism statistics\n"
            "  divergence    Dxy, Da, fixed/shared differences (needs --input2 or --pop-file)\n"
            "  fuliout       Fu & Li D/F with outgroup (needs --outgroup)\n"
            "  hka           HKA neutrality test across loci (needs --hka-file)\n"
            "  mk            McDonald-Kreitman test (needs --outgroup; coding aln)\n"
            "  kaks          Ka/Ks via Nei-Gojobori 1986 (coding alignment)\n"
            "  fufs          Fu's Fs neutrality test (Fu 1997)\n"
            "  sfs           Site frequency spectrum (folded; unfolded with --outgroup)\n"
            "  tstv          Transition/transversion ratio (Ts/Tv)\n"
            "  codon         Codon usage bias: RSCU + ENC (coding alignment)\n"
            "  faywu         Fay & Wu's H + Zeng's E (needs --outgroup)\n"
            "  fst           Population differentiation Fst  -  Hudson 1992 (needs --pop-file)\n"
            "  all           Run all applicable analyses\n"
            "\n"
            "Examples:\n"
            "  python dnasp.py --demo --output /tmp/dnasp_demo\n"
            "  python dnasp.py --input aln.fas --output results/\n"
            "  python dnasp.py --input aln.fas --analysis all --output results/\n"
            "  python dnasp.py --input aln.fas --analysis ld,recombination --output results/\n"
            "  python dnasp.py --input aln.nex --window 100 --step 25 --output results/\n"
            "  python dnasp.py --input pop1.fas --input2 pop2.fas --analysis divergence --output results/\n"
            "  python dnasp.py --input aln.fas --pop-file pops.txt --analysis divergence --output results/\n"
            "  python dnasp.py --input aln.fas --outgroup outgroupSeq --analysis fuliout --output results/\n"
            "  python dnasp.py --input aln.fas --hka-file hka_loci.tsv --analysis hka --output results/\n"
            "  python dnasp.py --input coding.fas --outgroup OutSeq --analysis mk --output results/\n"
            "  python dnasp.py --input coding.fas --analysis kaks --output results/\n"
            "  python dnasp.py --input aln.fas --analysis fufs --output results/\n"
            "  python dnasp.py --input aln.fas --analysis sfs --output results/\n"
            "  python dnasp.py --input aln.fas --outgroup OutSeq --analysis sfs --output results/\n"
            "  python dnasp.py --input aln.fas --analysis tstv --output results/\n"
            "  python dnasp.py --input coding.fas --analysis codon --output results/\n"
            "  python dnasp.py --input coding.fas --analysis codon,tstv,kaks --output results/\n"
            "  python dnasp.py --input aln.fas --outgroup OG --analysis faywu --output results/\n"
            "  python dnasp.py --input aln.fas --pop-file pops.txt --analysis fst --output results/\n"
        ),
    )
    p.add_argument("--version", "-V", action="version", version=f"dnasp {__version__}")
    p.add_argument("--input", "-i", type=Path, help="Input alignment (FASTA or NEXUS)")
    p.add_argument("--input2", type=Path, help="Second population alignment (for --analysis divergence)")
    p.add_argument("--pop-file", type=Path, dest="pop_file",
                   help="Population assignment file (TSV: seq_name<TAB>population)")
    p.add_argument("--outgroup", type=str, default=None,
                   help="Sequence name to use as outgroup (for --analysis fuliout). "
                        "This sequence is removed from the ingroup.")
    p.add_argument("--hka-file", type=Path, dest="hka_file", default=None,
                   help="HKA locus file (TSV: locus<TAB>S<TAB>D<TAB>n) for --analysis hka")
    p.add_argument("--output", "-o", type=Path, default=Path("dnasp_output"),
                   help="Output directory (default: dnasp_output/)")
    p.add_argument("--analysis", "-a", default="polymorphism",
                   help="Comma-separated analyses to run (default: polymorphism). Use 'all' for everything.")
    p.add_argument("--window", "-w", type=int, default=0,
                   help="Sliding window size in bp (0 = whole alignment only)")
    p.add_argument("--step", "-s", type=int, default=0,
                   help="Sliding window step size in bp (default: same as window)")
    p.add_argument("--demo", action="store_true",
                   help="Run on built-in synthetic demo data")
    return p


def _parse_analyses(analysis_str: str) -> set[str]:
    if analysis_str.strip().lower() == "all":
        return VALID_ANALYSES.copy()
    parts = {a.strip().lower() for a in analysis_str.split(",")}
    unknown = parts - VALID_ANALYSES - {"all"}
    if unknown:
        print(
            f"Warning: unknown analysis/analyses ignored: {unknown}\n"
            f"Valid options: {', '.join(sorted(VALID_ANALYSES))}, all",
            file=sys.stderr,
        )
    valid = parts & VALID_ANALYSES
    valid.add("polymorphism")  # always run
    return valid


def main(argv: Optional[list[str]] = None) -> int:
    parser = build_parser()
    args = parser.parse_args(argv)

    analyses = _parse_analyses(args.analysis)

    pop_assignments: Optional[dict[str, str]] = None
    if args.pop_file:
        if not args.pop_file.exists():
            print(f"Error: pop-file not found: {args.pop_file}", file=sys.stderr)
            return 1
        pop_assignments = load_pop_file(args.pop_file)

    aln2: Optional[Alignment] = None
    if args.input2:
        if not args.input2.exists():
            print(f"Error: --input2 file not found: {args.input2}", file=sys.stderr)
            return 1
        aln2 = load_alignment(args.input2)

    hka_loci: Optional[list[HKALocus]] = None
    if args.hka_file:
        if not args.hka_file.exists():
            print(f"Error: --hka-file not found: {args.hka_file}", file=sys.stderr)
            return 1
        hka_loci = load_hka_file(args.hka_file)
        print(f"Loaded {len(hka_loci)} HKA loci from {args.hka_file}")

    step = args.step if args.step > 0 else args.window
    cli_args = sys.argv[1:]

    if args.demo:
        output_dir = args.output
        output_dir.mkdir(parents=True, exist_ok=True)
        demo_path = output_dir / "demo_input.fas"
        demo_path.write_text(DEMO_FASTA, encoding="utf-8")
        pop_path = output_dir / "demo_pop.txt"
        pop_path.write_text(DEMO_POP_FILE, encoding="utf-8")
        print(DEMO_DESCRIPTION)
        print(f"Demo input written to: {demo_path}")
        print(f"Demo pop file written to: {pop_path}")
        # Demo always runs all analyses (hka skipped  -  no --hka-file)
        demo_analyses = VALID_ANALYSES - {"hka"}
        _run(demo_path, output_dir, step, args.window, demo_analyses,
             DEMO_POP_ASSIGNMENTS, None, cli_args,
             outgroup_name="outgroup")
        return 0

    if not args.input:
        parser.error("Provide --input <file> or --demo")

    if not args.input.exists():
        print(f"Error: input file not found: {args.input}", file=sys.stderr)
        return 1

    _run(
        args.input, args.output, args.window, step,
        analyses, pop_assignments, aln2, cli_args,
        outgroup_name=args.outgroup,
        hka_loci=hka_loci,
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
