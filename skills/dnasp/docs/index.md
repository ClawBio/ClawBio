#  DnaSP Skill  -  Documentation

**Version 0.4.0 · ClawBio · MIT Licence**

A Python reimplementation of all major population genetics modules computed by [DnaSP 6](https://www.ub.edu/dnasp/) (Rozas et al. 2017), designed to run from the command line anywhere ClawBio is installed.

> **About this implementation.** This skill is a Python reimplementation of the core statistical modules of DnaSP 6, developed by David De Lorenzo for the ClawBio platform. The original algorithms and Visual Basic source code were authored by Julio Rozas, Albert Ferrer-Mata, Juan Carlos Sánchez-DelBarrio, Sara Guirao-Rico, Pablo Librado, Sebastián Ramos-Onsins, and Alejandro Sánchez-Gracia, and were shared by Julio Rozas. All statistical formulas are transcribed from the primary literature and cross-checked against the DnaSP 6 source; any implementation errors are the responsibility of the adapter, not the original authors.

This document serves two audiences:

- **Users**  -  explains what each module computes, how to run it, and how to interpret results.
- **Claude** (via SKILL.md)  -  provides the statistical detail and interpretation guidance needed to answer questions about DnaSP analyses without hallucinating formulas.

---

## Contents

1. [Overview](#overview)
2. [Quick start](#quick-start)
3. [Input formats](#input-formats)
4. [Running the examples](#running-the-examples)
5. [Running the tests](#running-the-tests)
6. [Statistics reference](#statistics-reference)
   - [Polymorphism & neutrality tests](#polymorphism--neutrality-tests)
   - [Linkage disequilibrium](#linkage-disequilibrium)
   - [Recombination](#recombination)
   - [Population size history  -  mismatch distribution](#population-size-history--mismatch-distribution)
   - [InDel polymorphism](#indel-polymorphism)
   - [Divergence between populations](#divergence-between-populations)
   - [Fu & Li D/F with outgroup](#fu--li-df-with-outgroup)
   - [HKA test](#hka-test)
   - [McDonald-Kreitman test](#mcdonald-kreitman-test)
   - [Ka/Ks (dN/dS)](#kaks-dnds)
   - [Fu's Fs test](#fus-fs-test)
   - [Site frequency spectrum](#site-frequency-spectrum)
   - [Transition/Transversion ratio](#transitiontransversion-ratio)
   - [Codon usage bias](#codon-usage-bias)
   - [Fay & Wu's H and Zeng's E](#fay--wus-h-and-zengs-e)
   - [Population differentiation  -  Fst](#population-differentiation--fst)
7. [Output files](#output-files)
8. [Sliding-window analysis](#sliding-window-analysis)
9. [Implementation notes and gotchas](#implementation-notes-and-gotchas)
10. [Citations](#citations)

---

## Overview

The DnaSP skill computes the standard diversity, neutrality, linkage, and divergence statistics that population geneticists typically generate with the Windows-only DnaSP GUI. Results match DnaSP 6 output to floating-point precision for the standard single-alignment analysis mode.

**Modules implemented:**

| Module (`--analysis` value) | What it computes |
|---|---|
| `polymorphism` (default) | π, k, S, Eta, H, Hd, VarHd, θ_W, Tajima's D, Fu & Li D\*/F\*, R2, GC |
| `ld` | D, D', R² per site pair; ZnS, Za, ZZ genome-wide; LD decay scatter plot |
| `recombination` | Rm (minimum recombination events, four-gamete test) |
| `popsize` | Mismatch distribution, raggedness r, coefficient of variation |
| `indel` | InDel events, InDel haplotypes, k(i), π(i), θ(i), Tajima's D(i) |
| `divergence` | Dxy, Da, fixed differences, shared & private polymorphisms |
| `fuliout` | Fu & Li D/F using an outgroup sequence to polarise mutations |
| `hka` | HKA test  -  multi-locus neutrality test (Hudson, Kreitman & Aguadé 1987) |
| `mk` | McDonald-Kreitman test  -  adaptive evolution via Pn/Ps vs Dn/Ds (McDonald & Kreitman 1991) |
| `kaks` | Ka/Ks (dN/dS)  -  pairwise synonymous and nonsynonymous substitution rates (Nei & Gojobori 1986) |
| `fufs` | Fu's Fs test  -  haplotype-frequency neutrality test (Fu 1997; Ewens sampling formula) |
| `sfs` | Site frequency spectrum  -  folded (always) and unfolded with `--outgroup` |
| `tstv` | Transition/transversion ratio (Ts/Tv)  -  substitution pattern across all pairwise comparisons |
| `codon` | Codon usage bias  -  RSCU (Sharp & Li 1987) and ENC (Wright 1990) for in-frame coding alignments |
| `faywu` | Fay & Wu's H (2000) + Zeng's E (2006)  -  requires `--outgroup`; tests for selective sweeps via high-frequency derived alleles |
| `fst` | Population differentiation  -  pairwise Hudson Fst (Hudson et al. 1992); requires `--pop-file` or `--input2` |

**What it does not do (currently):**

- Sequence alignment (use MAFFT or MUSCLE first  -  sequences must be pre-aligned)
- Variant calling from VCF files
- Parametric simulation for SFS or mismatch P-values

---

## Quick start

```bash
# Polymorphism + neutrality tests (default)
python skills/dnasp/dnasp.py --input alignment.fas --output results/

# Specific analyses
python skills/dnasp/dnasp.py --input alignment.fas --analysis ld,recombination --output results/

# All analyses (alignment only  -  no second population)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis polymorphism,ld,recombination,popsize,indel \
    --output results/

# Divergence  -  two separate FASTA files
python skills/dnasp/dnasp.py \
    --input pop1.fas --input2 pop2.fas \
    --analysis divergence \
    --output results/

# Divergence  -  one alignment + population assignment file
python skills/dnasp/dnasp.py \
    --input combined.fas --pop-file populations.txt \
    --analysis divergence \
    --output results/

# All analyses including divergence
python skills/dnasp/dnasp.py \
    --input pop1.fas --input2 pop2.fas \
    --analysis all \
    --output results/

# Fu & Li D/F with outgroup (outgroup sequence named "OutSeq" in alignment)
python skills/dnasp/dnasp.py \
    --input alignment_with_outgroup.fas \
    --analysis fuliout \
    --outgroup OutSeq \
    --output results/

# HKA test (multi-locus TSV file)
python skills/dnasp/dnasp.py \
    --analysis hka \
    --hka-file loci.tsv \
    --output results/

# McDonald-Kreitman test (outgroup sequence named "OutSeq" in the alignment)
python skills/dnasp/dnasp.py \
    --input coding_aln_with_outgroup.fas \
    --analysis mk \
    --outgroup OutSeq \
    --output results/

# Ka/Ks  -  Nei-Gojobori pairwise dN/dS (in-frame coding alignment)
python skills/dnasp/dnasp.py \
    --input coding_aln.fas \
    --analysis kaks \
    --output results/

# Fu's Fs test
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis fufs \
    --output results/

# Site frequency spectrum (folded only)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis sfs \
    --output results/

# Site frequency spectrum (folded + unfolded with outgroup)
python skills/dnasp/dnasp.py \
    --input aln_with_outgroup.fas \
    --outgroup OutSeq \
    --analysis sfs \
    --output results/

# Transition/transversion ratio (any alignment)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis tstv \
    --output results/

# Codon usage bias (in-frame coding alignment)
python skills/dnasp/dnasp.py \
    --input coding.fas \
    --analysis codon \
    --output results/

# Full coding panel (Ka/Ks + MK + Ts/Tv + Codon usage)
python skills/dnasp/dnasp.py \
    --input coding.fas \
    --outgroup OutSeq \
    --analysis kaks,mk,tstv,codon \
    --output results/

# Fay & Wu's H + Zeng's E (outgroup named "OutSeq" in alignment)
python skills/dnasp/dnasp.py \
    --input aln_with_outgroup.fas \
    --outgroup OutSeq \
    --analysis faywu \
    --output results/

# Population differentiation Fst (pop-file assigns sequences to populations)
python skills/dnasp/dnasp.py \
    --input combined.fas \
    --pop-file populations.txt \
    --analysis fst \
    --output results/

# Complete population panel (polymorphism + divergence + Fst + Fay & Wu)
python skills/dnasp/dnasp.py \
    --input combined_with_outgroup.fas \
    --pop-file populations.txt \
    --outgroup OutSeq \
    --analysis polymorphism,divergence,fst,faywu \
    --output results/

# Sliding window (100 bp window, 25 bp step)
python skills/dnasp/dnasp.py --input alignment.fas --window 100 --step 25 --output results/

# Built-in demo (no input file required)
python skills/dnasp/dnasp.py --demo --output /tmp/dnasp_demo
```

**All CLI options:**

| Option | Description | Default |
|---|---|---|
| `--input FILE` | FASTA or NEXUS alignment |  -  |
| `--input2 FILE` | Second-population FASTA (for `divergence`) |  -  |
| `--pop-file FILE` | Population assignment TSV (alternative to `--input2`) |  -  |
| `--analysis LIST` | Comma-separated analyses or `all` | `polymorphism` |
| `--output DIR` | Output directory (created if absent) | `dnasp_output/` |
| `--window INT` | Sliding window size in bp (0 = whole alignment only) | `0` |
| `--step INT` | Sliding window step in bp | `= window` |
| `--outgroup STR` | Sequence name to use as outgroup (for `fuliout` and `mk`) |  -  |
| `--hka-file FILE` | HKA locus TSV file (for `hka`) |  -  |
| `--demo` | Run on the built-in 11-sequence × 300-bp synthetic alignment (10 ingroup + 1 outgroup, 2 populations) |  -  |

**Population assignment file format** (for `--pop-file`):

```
# Tab-separated: sequence_name<TAB>population_name
# Lines starting with # are ignored
seq1	Pop_Africa
seq2	Pop_Africa
seq3	Pop_Europe
seq4	Pop_Europe
```

---

## Input formats

### FASTA

Standard FASTA files with any extension (`.fas`, `.fa`, `.fasta`). Multi-line sequences are fully supported. DnaSP-style headers (`>'name'  [comment]`) are handled automatically  -  the single-quote delimiters and bracketed comments are stripped, leaving only the sequence name.

```fasta
>Seq1
ATCGATCGATCG
>Seq2
ATCGTTCGATCG
```

**Requirement:** all sequences must be the same length (pre-aligned).

### NEXUS

NEXUS files (`.nex`, `.nexus`, `.nxs`) in both sequential and interleaved formats. Supported features:

- `MATCHCHAR=.`  -  dot characters are expanded relative to the first sequence
- `INTERLEAVE`  -  multi-block interleaved alignments are reassembled correctly
- `CHARSTATELABELS`  -  position labels are read but not required

A typical DnaSP-exported NEXUS header looks like:

```nexus
#NEXUS
BEGIN DATA;
  DIMENSIONS NTAX=6 NCHAR=300;
  FORMAT DATATYPE=DNA MISSING=? GAP=- MATCHCHAR=.;
  MATRIX
    Seq1  ATCGATCG...
    ...
  ;
END;
```

### Gap treatment

DnaSP applies **complete deletion** by default: any alignment column that contains a gap (`-`), missing (`?`), or ambiguous (`N`) character in *any* sequence is excluded from all calculations. The number of net sites used (L\_net) is reported in the output and may be substantially smaller than the total alignment length if the data are gap-rich.

For divergence analyses, complete deletion is applied across *both* populations combined  -  a gap in either population removes the column.

---

## Running the examples

Two example files are provided in `skills/dnasp/examples/`:

### demo\_simple.fas  -  synthetic alignment

A 6-sequence × 10-bp alignment constructed to have exact, verifiable statistics.

```bash
python skills/dnasp/dnasp.py \
    --input skills/dnasp/examples/demo_simple.fas \
    --output /tmp/demo_simple
```

Expected whole-alignment output:

| Statistic | Expected value |
|---|---|
| L\_net | 10 |
| S (segregating sites) | 2 |
| Eta (total mutations) | 2 |
| H (haplotypes) | 3 |
| Hd (haplotype diversity) | 0.7333 |
| π (nucleotide diversity) | 0.08667 |
| θ\_W per site | 0.08759 |
| Tajima's D | −0.047 |
| Fu & Li D\* | +0.062 |
| Fu & Li F\* | +0.040 |
| R2 | ~0.229 |

### demo\_rp49.fas  -  real Drosophila data

17 *Drosophila* sequences of the *rp49* region (~300 bp). This is a real DnaSP example dataset and produces non-trivial statistics suitable for validating the implementation.

```bash
python skills/dnasp/dnasp.py \
    --input skills/dnasp/examples/demo_rp49.fas \
    --analysis all \
    --output /tmp/demo_rp49
```

Add `--window 50 --step 10` for a sliding-window plot:

```bash
python skills/dnasp/dnasp.py \
    --input skills/dnasp/examples/demo_rp49.fas \
    --window 50 --step 10 \
    --output /tmp/demo_rp49_window
```

---

## Running the tests

The test suite uses [pytest](https://docs.pytest.org/) and covers all statistics, both parsers, all new modules, edge cases, and the full analysis pipeline.

```bash
# From the skill root
cd skills/dnasp
pip install pytest --break-system-packages   # if not already installed
python -m pytest tests/ -v
```

Expected output: **355 tests pass**, 0 failures.

Individual test groups:

```bash
python -m pytest tests/ -v -k "parser"        # FASTA and NEXUS parsers
python -m pytest tests/ -v -k "tajima"        # Tajima's D
python -m pytest tests/ -v -k "fu_li"         # Fu & Li D* and F*
python -m pytest tests/ -v -k "ld"            # Linkage disequilibrium
python -m pytest tests/ -v -k "recombination" # Recombination (Rm)
python -m pytest tests/ -v -k "mismatch"      # Mismatch distribution
python -m pytest tests/ -v -k "indel"         # InDel polymorphism
python -m pytest tests/ -v -k "divergence"    # Divergence analysis
python -m pytest tests/ -v -k "pop"           # Population file parsing
python -m pytest tests/ -v -k "fuliout or outgroup or derived" # Fu & Li with outgroup
python -m pytest tests/ -v -k "hka"           # HKA test
python -m pytest tests/ -v -k "mk or fisher"  # McDonald-Kreitman test
python -m pytest tests/ -v -k "kaks or ka_ks or jc" # Ka/Ks (dN/dS)
python -m pytest tests/ -v -k "FuFs or Ewens or Stirling" # Fu's Fs
python -m pytest tests/ -v -k "SFS"           # Site frequency spectrum
python -m pytest tests/ -v -k "TsTv"          # Transition/Transversion ratio
python -m pytest tests/ -v -k "CodonUsage or Synonymous" # Codon usage bias (RSCU + ENC)
python -m pytest tests/ -v -k "FayWu or faywu"          # Fay & Wu's H + Zeng's E
python -m pytest tests/ -v -k "Fst or fst"              # Population differentiation Fst
```

---

## Statistics reference

All formulas are transcribed from the primary literature and cross-checked against the DnaSP 6 Visual Basic source code. Each section gives the definition, formula, and reference for every statistic.

---

### Polymorphism & neutrality tests

#### Sequence summary symbols

| Symbol | Name |
|---|---|
| n | Number of sequences |
| L | Total alignment length |
| L\_net | Net sites after complete deletion |
| GC | GC content (computed on L\_net sites only) |

#### S  -  Number of segregating sites

Count of alignment columns (after complete deletion) at which at least two sequences carry different nucleotides.

#### Eta  -  Total number of mutations

Total nucleotide variants across all segregating sites. For strictly biallelic sites Eta = S; for triallelic or quadriallelic sites Eta > S (each extra allele counted separately). Fu & Li tests use Eta; Tajima's D uses S.

#### k  -  Average pairwise differences (absolute)

```
k = Σ d(i,j) / C(n,2)    over all i < j
```

where d(i,j) is the count of sites differing between sequences i and j. Reference: Tajima (1983), eq. A3.

#### π  -  Nucleotide diversity (per site)

```
π = k / L_net
```

Reference: Nei (1987), eq. 10.5.

#### Haplotypes

- **H**  -  Count of distinct sequence variants.
- **Hd**  -  Haplotype (gene) diversity:
  ```
  Hd = n/(n−1) × (1 − Σ pᵢ²)
  ```
  where pᵢ is the frequency of haplotype i (Nei & Tajima 1981).
- **VarHd**  -  Sampling variance:
  ```
  VarHd = 2/(n(n−1)) × [2(n−2)(Σxᵢ³ − (Σxᵢ²)²) + Σxᵢ² − (Σxᵢ²)²]
  ```

#### θ\_W  -  Watterson's theta (Watterson 1975)

```
a₁ = Σ 1/i    (i = 1 to n−1)

θ_W (per gene) = S / a₁
θ_W (per site) = S / (a₁ × L_net)
```

---

#### Tajima's D (Tajima 1989)

Compares k (π × L_net) with Watterson's estimate S/a₁.

```
D = (k − S/a₁) / √(e₁S + e₂S(S−1))
```

Coefficients e₁ and e₂ are derived from b₁, b₂, c₁, c₂ as in Tajima (1989, eq. 29-38).

**Interpretation:**
- D < 0  -  excess low-frequency variants → population expansion or purifying selection
- D > 0  -  excess intermediate-frequency variants → balancing selection or bottleneck
- D ≈ 0  -  consistent with neutrality + constant population size

**Requirements:** n ≥ 3, S > 0. Returns `n.a.` otherwise.

---

#### Fu & Li D\* (Fu & Li 1993; Simonsen et al. 1995)

Compares singletons (η\_s) with total mutations.

```
D* = (S/Aₙ − η_s(n−1)/n) / √(u_D·S + v_D·S²)
```

Variance coefficients u\_D, v\_D from Simonsen et al. (1995), equations A3-A4.

**Definition of η\_s:** Following DnaSP's `BusqSingletones()`, a singleton is any allele present in exactly 1 sequence. Alleles in n−1 sequences are *not* counted (no folded singletons). This matches DnaSP v6 behaviour precisely.

---

#### Fu & Li F\* (Fu & Li 1993; Simonsen et al. 1995)

Compares k with singletons.

```
F* = (k − η_s(n−1)/n) / √(u_F·S + v_F·S²)
```

Variance coefficients u\_F, v\_F from Simonsen et al. (1995), equations A5-A6.

> **Note:** DnaSP v6 RAD-seq multi-MSA mode uses Achaz (2009) variances for F\*; the standard FASTA/NEXUS mode uses Simonsen (1995). This skill implements the standard mode. D\* agrees exactly regardless.

---

#### R2  -  Ramos-Onsins & Rozas (2002)

```
R2 = √(Σᵢ(Uᵢ − k/2)² / n) / Sw
```

where Uᵢ = singletons attributed to sequence i; Sw = total mutations. Always ≥ 0. Small R2 suggests recent population expansion.

Implementation matches DnaSP's `JulioSebas_R2_Calculo()` (revision April 2002).

---

### Linkage Disequilibrium

Run with `--analysis ld`. Requires ≥ 2 strictly biallelic sites after complete deletion.

**What counts as a biallelic site:** exactly 2 nucleotide states observed, no gaps in any sequence at that position. Multiallelic sites and sites with any gap are excluded from LD analysis.

#### D  -  Linkage disequilibrium coefficient (Lewontin & Kojima 1960)

For a pair of biallelic sites i (alleles A/a) and j (alleles B/b):

```
D = p_AB − p_A × p_B
```

where p_AB is the observed frequency of the AB haplotype, p_A is the frequency of allele A, and p_B is the frequency of allele B. D = 0 indicates gametic equilibrium.

#### D'  -  Normalised D (Lewontin 1964)

```
D' = D / D_max

D_max = min(p_A × p_b, p_a × p_B)   when D > 0
D_max = min(p_A × p_B, p_a × p_b)   when D < 0
```

D' ranges from −1 to +1. |D'| = 1 means no evidence of historical recombination between the two sites (given the data). Note: D' is inflated toward ±1 when sample size is small.

#### R²  -  Squared correlation coefficient (Hill & Robertson 1968)

```
R² = D² / (p_A × p_a × p_B × p_b)
```

R² ranges from 0 (independent) to 1 (complete association). More robust than D' when sample sizes are small. Used to measure LD decay with physical distance.

#### Chi-square test for LD

```
χ² = n × R²    (1 degree of freedom)
p = erfc(√(χ²/2))
```

Uses only the Python standard library (`math.erfc`). Exact p-value for df = 1.

#### ZnS  -  Genome-wide LD statistic (Kelly 1997)

Average R² over all pairs of biallelic sites:

```
ZnS = Σ R²(i,j) / C(S_b, 2)
```

where S_b is the number of biallelic sites and the sum is over all C(S_b, 2) pairs.

#### Za  -  Adjacent-pairs statistic (Rozas et al. 2001)

Average R² restricted to *consecutive* biallelic sites (i.e., adjacent pairs in the biallelic-site ordering):

```
Za = Σ R²(i, i+1) / (S_b − 1)
```

#### ZZ = Za − ZnS (Rozas et al. 2001)

ZZ > 0 indicates that LD is stronger between physically adjacent sites than between non-adjacent ones  -  the signature expected when recombination erodes long-range LD. ZZ ≈ 0 indicates uniform LD across distances.

**Output files produced by `--analysis ld`:**
- `ld_pairs.tsv`  -  one row per site pair: positions, distance, n_valid, D, D', R², χ², p
- `figures/ld_decay.png`  -  R² vs distance scatter plot

---

### Recombination

Run with `--analysis recombination`. Uses the same biallelic site set as the LD module.

#### Four-gamete test (Hudson & Kaplan 1985)

A pair of biallelic sites is **incompatible** when all four possible gamete combinations (AB, Ab, aB, ab) are observed in the data. Under an infinite-sites model without recombination, only three gamete types can exist. Observing all four is therefore evidence of at least one recombination event between those sites.

For sites at positions i and j, define gametes as (allele at i, allele at j). If all four combinations are present → incompatible pair.

#### Rm  -  Minimum number of recombination events (Hudson & Kaplan 1985)

Rm is the minimum number of past recombination events required to explain all observed incompatible pairs. It is computed by the **interval-stabbing greedy algorithm**:

1. Represent each incompatible pair (i, j) as an interval [left, right] on the alignment.
2. Sort intervals by right endpoint.
3. For each interval not yet "stabbed" (i.e., its left endpoint > current rightmost point), place a recombination point at its right endpoint and increment Rm.

This gives the minimum number of points needed to stab all intervals  -  equivalent to Rm.

**Important:** Rm is a **lower bound**. The true number of recombination events is ≥ Rm. Rm = 0 does not mean no recombination occurred; it means the data are consistent with no recombination.

**Interpretation:**
- Rm = 0  -  data are compatible with no recombination (no incompatible pairs)
- Rm ≥ 1  -  at least Rm recombination events are required
- High Rm relative to S suggests a high recombination rate

---

### Population size history  -  Mismatch distribution

Run with `--analysis popsize`.

#### Mismatch distribution

The observed distribution of pairwise nucleotide differences  -  how many sequence pairs differ by exactly 0 sites, 1 site, 2 sites, and so on.

Under a model of sudden population expansion, this distribution is expected to be unimodal and smooth (Rogers & Harpending 1992). Under constant population size, it is irregular (ragged).

The skill computes:
- **Histogram** of pairwise differences: `{d: count}` for each observed difference count d
- **Mean** (= k, consistent with the polymorphism module)
- **Variance** of pairwise differences
- **Coefficient of variation (CV)** = σ/μ

#### Raggedness r (Harpending 1994)

```
r = Σ (f(i) − f(i−1))²     i = 1 to max_d + 1
```

where f(i) is the proportion of pairs with exactly i differences (and f(−1) = 0 by convention).

**Interpretation:**
- r small (< 0.05) → smooth, unimodal distribution → signature of population expansion
- r large → ragged, irregular distribution → consistent with constant population size or decline
- Significance can be assessed by simulation (not currently implemented in this skill)

Reference: Harpending (1994), *Hum. Biol.* 66:591-600.

**Output produced by `--analysis popsize`:**
- `figures/mismatch.png`  -  mismatch distribution bar chart with observed distribution

---

### InDel polymorphism

Run with `--analysis indel`.

#### InDel event identification

A DnaSP InDel event is identified using the **diallelic option**: a maximal run of consecutive alignment columns where the same subset of sequences carries gaps (and the remaining sequences carry nucleotides). When the gap-bearing set changes, a new event begins. All-gap columns (gaps in every sequence) are not counted as InDel events.

This captures the biologically meaningful case: a gap in some sequences but not others at the same alignment region represents a single insertion/deletion event.

Overlapping events (columns where different subsets are simultaneously gapped) are noted but treated as separate events.

#### Statistics computed

- **n\_positions\_with\_gaps**  -  total alignment columns with any gap
- **n\_events**  -  number of InDel events identified
- **mean\_event\_length**  -  average length of InDel events in bp

After identifying events, each sequence is represented as a binary vector over events (1 = carries the gap, 0 = does not). Standard population genetics statistics are then computed on this binary representation:

- **InDel haplotypes (H\_i)**  -  distinct InDel patterns across sequences
- **InDel haplotype diversity (Hd\_i)**  -  haplotype diversity of InDel patterns
- **k\_indel**  -  mean pairwise differences in InDel pattern (count of events that differ)
- **π\_indel** = k\_indel / (total InDel positions)
- **θ\_indel**  -  Watterson's θ applied to InDel data (S = n\_events)
- **Tajima's D (InDel)**  -  computed on InDel data; same formula, interpreted the same way

**Interpretation:** Significant Tajima's D on InDel data independently from nucleotide D can indicate that insertion/deletion variants have a different history from substitutions  -  sometimes a signature of selection acting specifically on insertions or deletions.

---

### Divergence between populations

Run with `--analysis divergence`. Requires either `--input2` (two separate FASTA files, one per population) or `--pop-file` (population assignment TSV used to split a single alignment).

#### Dxy  -  Mean between-population divergence (Nei 1987)

Average number of nucleotide differences per site between a randomly drawn sequence from population 1 and one from population 2:

```
Dxy = [Σᵢ Σⱼ d(i,j) / (n₁ × n₂)] / L_net
```

where the sum is over all n₁ × n₂ cross-population pairs, d(i,j) is the pairwise difference count, and L_net is the net sites after complete deletion across both populations combined. Reference: Nei (1987), eq. 10.20.

#### Da  -  Net divergence

```
Da = Dxy − (π₁ + π₂) / 2
```

Subtracts the average within-population diversity from the gross between-population divergence. Da approximates the divergence that has accumulated *since* the two populations separated, removing the ancestral polymorphism shared by both. Da can be slightly negative by chance (especially when L_net is small).

Reference: Nei (1987), eq. 10.20 (net nucleotide divergence).

#### Fixed differences, shared polymorphisms, private polymorphisms (Hey 1991)

Each site (after complete deletion) is classified into one of the following categories:

- **Fixed difference**  -  each population is monomorphic at this site, but they carry different alleles. These sites accumulated after the populations split.
- **Shared polymorphism**  -  the same alleles are segregating in both populations. These sites were polymorphic in the common ancestor.
- **Private polymorphism (Pop 1)**  -  site is polymorphic in Pop 1 but the Pop 1 alleles do not overlap with Pop 2.
- **Private polymorphism (Pop 2)**  -  site is polymorphic in Pop 2 but the Pop 2 alleles do not overlap with Pop 1.

Reference: Hey (1991), *Genetics* 128:831-840.

**Interpretation:**
- Many fixed differences + few shared → long divergence time, likely little ongoing gene flow
- Many shared + few fixed → recent divergence or ongoing gene flow
- Many private polymorphisms in one population → that population has expanded or has a higher mutation rate

#### Within-population diversity in the divergence module

The divergence module also reports π₁ and π₂ (nucleotide diversity within each population, computed on L_net sites after the joint gap mask), providing a direct comparison of within-population and between-population variation in the same output.

---

### Fu & Li D/F with outgroup

Run with `--analysis fuliout --outgroup SEQNAME`. The named sequence is extracted from the alignment as the outgroup and removed from the ingroup before computation.

#### Biological motivation

The standard Fu & Li D\* and F\* tests use folded (unrooted) data  -  a singleton is any allele present in exactly one sequence, regardless of whether it is ancestral or derived. When an outgroup is available, mutations can be *polarised*: we can determine which allele is ancestral (matches the outgroup) and which is derived (differs from the outgroup). This allows a more powerful test using the *external* (derived) branch of the genealogy.

#### Outgroup polarisation algorithm

For each alignment column, after applying **complete deletion** across all ingroup sequences *and* the outgroup:

1. Identify the ancestral allele as the nucleotide carried by the outgroup sequence.
2. If the ancestral allele is absent from all ingroup sequences at that column, skip the column (uninformative for polarisation).
3. Count derived alleles: any ingroup allele that differs from the outgroup allele.
4. **η (eta)**  -  total derived mutations at this site = number of derived ingroup alleles (counting each sequence carrying a derived allele). Summed across all informative sites.
5. **η_e (eta_e, external mutations)**  -  derived mutations in exactly one ingroup sequence (derived singletons). These occur at the external branches of the genealogy.

#### D and F statistics (Fu & Li 1993)

Using η and η\_e as defined above:

```
D = (η_e − η/Aₙ) / √(u_D·η + v_D·η²)
```

```
F = (k̄ − η_e) / √(u_F·η + v_F·η²)
```

where k̄ is the mean pairwise difference count and Aₙ = Σ 1/i for i = 1 to n−1 (the same harmonic number used in θ_W). Variance coefficients u\_D, v\_D, u\_F, v\_F follow Simonsen et al. (1995) Appendix A. The exact outgroup-specific coefficients from Appendix B differ only slightly for small n and may be substituted once validated against DnaSP 6 output.

**Requirements:** n ≥ 3, η > 0. Returns `n.a.` if the sequence named by `--outgroup` is not found in the alignment, or if no informative sites remain after polarisation.

#### Interpretation

| Statistic | Interpretation |
|---|---|
| D < 0 | Excess derived singletons → expansion or purifying selection on derived alleles |
| D > 0 | Deficit of derived singletons → balancing selection |
| D ≈ 0 | Consistent with neutrality + constant size |
| F < 0 | k̄ smaller than expected given η_e → purifying selection or expansion |
| F > 0 | k̄ larger than expected → balancing selection |

D and F (with outgroup) are generally more sensitive than D\* and F\* when an outgroup can polarise mutations reliably.

#### CLI usage

```bash
# Outgroup sequence named "Dmel" in the alignment
python skills/dnasp/dnasp.py \
    --input dsim_with_dmel.fas \
    --analysis fuliout \
    --outgroup Dmel \
    --output results/

# Combined with polymorphism
python skills/dnasp/dnasp.py \
    --input dsim_with_dmel.fas \
    --analysis polymorphism,fuliout \
    --outgroup Dmel \
    --output results/
```

Reference: Fu & Li (1993); Simonsen et al. (1995).

---

### HKA test

Run with `--analysis hka --hka-file loci.tsv`. The HKA test does not use a sequence alignment  -  it reads pre-computed counts from a locus TSV file.

#### Biological motivation

Hudson, Kreitman & Aguadé (1987) showed that under strict neutrality, the ratio of within-species polymorphism (S) to between-species divergence (D) should be the same across all loci. A locus under positive selection will have elevated D relative to S; a locus under balancing selection will have elevated S relative to D.

#### HKA locus file format

The `--hka-file` argument points to a tab-separated file with one row per locus:

```
# locus_name  n_seqs  S_poly  D_diverg
locus_A       20      12      8
locus_B       18      3       25
locus_C       22      9       11
```

Columns:

| Column | Description |
|---|---|
| `locus_name` | Identifier string (used in output) |
| `n_seqs` | Number of ingroup sequences at this locus |
| `S_poly` | Number of segregating sites (polymorphisms within species) |
| `D_diverg` | Number of fixed differences between species |

Lines beginning with `#` are treated as comments and ignored. A header line beginning with `#` is therefore valid. Lines containing only whitespace are skipped.

#### Statistical model

Under neutrality, for locus i with n_i ingroup sequences:

```
E[Sᵢ] = θ̂ᵢ · fᵢ    where fᵢ = Σ 1/j  (j = 1 to nᵢ−1)   (a₁ for locus i)
E[Dᵢ] = θ̂ᵢ · (1 + 2T̂)
```

The MLE of the relative mutation rates θ̂ᵢ and the divergence time T̂ (in units of 2N generations) are estimated jointly. Because the θ̂ᵢ can be eliminated analytically, T̂ satisfies the constraint:

```
Σᵢ Dᵢ / (1 + 2T̂) = Σᵢ (Sᵢ + Dᵢ) / (fᵢ + 1 + 2T̂)
```

T̂ is found by bisection over [0, 1000]. Once T̂ is known:

```
θ̂ᵢ = (Sᵢ + Dᵢ) / (fᵢ + 1 + 2T̂)
E[Sᵢ] = θ̂ᵢ · fᵢ
E[Dᵢ] = θ̂ᵢ · (1 + 2T̂)
```

#### Chi-square statistic

```
χ² = Σᵢ [ (Sᵢ − E[Sᵢ])² / E[Sᵢ] + (Dᵢ − E[Dᵢ])² / E[Dᵢ] ]
```

with **df = k − 1** where k is the number of loci. The p-value uses the regularised incomplete gamma function Q(df/2, χ²/2), computed without scipy using a series expansion (x < a+1) or Lentz continued fraction (x ≥ a+1).

**Requirements:** k ≥ 2 loci. Returns `n.a.` for p-value if df < 1. Each locus must have n ≥ 2 sequences (for fᵢ > 0).

#### Interpretation

| Result | Interpretation |
|---|---|
| p ≥ 0.05 | Fail to reject neutrality  -  polymorphism/divergence ratios consistent across loci |
| p < 0.05 | At least one locus deviates significantly  -  candidate for selection |
| Large (Sᵢ − E[Sᵢ])² | That locus has unusual polymorphism → possible balancing selection |
| Large (Dᵢ − E[Dᵢ])² | That locus has unusual divergence → possible positive selection |
| T̂ large | Long divergence time between species at typical loci |

Per-locus θ̂ᵢ and expected counts are reported in the output, allowing identification of which locus or loci drive a significant result.

#### CLI usage

```bash
# Basic HKA test (no alignment required)
python skills/dnasp/dnasp.py \
    --analysis hka \
    --hka-file loci.tsv \
    --output results/

# HKA alongside polymorphism (uses same alignment; HKA reads its own file)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis polymorphism,hka \
    --hka-file loci.tsv \
    --output results/
```

Reference: Hudson et al. (1987).

---

### McDonald-Kreitman test

Run with `--analysis mk --outgroup SEQNAME`. The alignment must be an in-frame coding sequence (length divisible by 3). The named outgroup sequence is extracted and used to classify fixed differences; it is **not** included in the ingroup polymorphism counts.

#### Biological motivation

McDonald & Kreitman (1991) observed that under strict neutrality, the ratio of nonsynonymous to synonymous changes should be the same for polymorphisms (within-species) as for fixed differences (between-species). Positive (adaptive) selection accelerates the fixation of nonsynonymous changes, increasing Dn/Ds relative to Pn/Ps.

#### Classification table

The test builds a 2 × 2 contingency table:

|  | Polymorphic (within-species) | Fixed (between-species) |
|---|---|---|
| Nonsynonymous | Pn | Dn |
| Synonymous | Ps | Ds |

A codon position is classified as:

- **Polymorphic**  -  at least two ingroup sequences carry different nucleotides at that position.
- **Fixed**  -  all ingroup sequences agree at that position, but the outgroup carries a different nucleotide.
- **Synonymous / nonsynonymous**  -  determined by the NCBI standard genetic code (transl\_table=1): a change is synonymous if the two codons encode the same amino acid; nonsynonymous otherwise.

**Complete deletion at codon level**: any codon where any nucleotide position in any sequence (ingroup or outgroup) is not in {A, T, C, G} is skipped. Stop codons (in ingroup or outgroup) are also skipped.

#### Derived statistics

```
α  = 1 − (Ds × Pn) / (Dn × Ps)    # proportion of nonsynonymous substitutions driven by positive selection
NI = (Pn / Ps) / (Dn / Ds)         # neutrality index (= 1 under strict neutrality)
DoS = Dn/(Dn+Ds) − Pn/(Pn+Ps)     # direction of selection
```

All three are `None` when any denominator is zero.

**α interpretation:** α > 0 → fraction of nonsynonymous fixations that are adaptive; α < 0 → excess nonsynonymous polymorphism (slightly deleterious alleles segregating).

**NI interpretation:** NI = 1 under neutrality; NI < 1 → positive selection; NI > 1 → slightly deleterious polymorphisms (Ka > expected from Ks).

**DoS interpretation:** Scale-free measure ranging −1 to +1. DoS > 0 → divergence more nonsynonymous than polymorphism (positive selection signature); DoS < 0 → divergence more synonymous (purifying selection removing nonsynonymous variants before fixation).

#### Statistical test

Fisher's exact test on the 2 × 2 table (Pn, Ps, Dn, Ds), two-tailed. Computed via the hypergeometric PMF using log-factorials (`math.lgamma`)  -  no scipy required. P-value < 0.05 → ratio of nonsynonymous to synonymous changes differs significantly between polymorphism and divergence.

**Requirements:** n ≥ 2 ingroup sequences, alignment length % 3 == 0, `--outgroup` provided.

#### CLI usage

```bash
# Basic MK test
python skills/dnasp/dnasp.py \
    --input coding_aln_with_outgroup.fas \
    --analysis mk \
    --outgroup OutSeq \
    --output results/

# MK test alongside polymorphism
python skills/dnasp/dnasp.py \
    --input coding_aln_with_outgroup.fas \
    --analysis polymorphism,mk \
    --outgroup OutSeq \
    --output results/
```

Reference: McDonald & Kreitman (1991).

---

### Ka/Ks (dN/dS)

Run with `--analysis kaks`. The alignment must be an in-frame coding sequence (length divisible by 3). No outgroup is required  -  Ka/Ks is computed pairwise across all ingroup sequences.

#### Biological motivation

The ratio of nonsynonymous substitution rate (Ka, or dN) to synonymous substitution rate (Ks, or dS) measures the selective pressure acting on a protein-coding gene. Synonymous sites serve as an internal neutral reference.

- **ω = Ka/Ks < 1**  -  purifying (negative) selection: most nonsynonymous mutations are deleterious and removed.
- **ω ≈ 1**  -  neutral evolution: nonsynonymous and synonymous rates similar.
- **ω > 1**  -  positive (adaptive) selection: nonsynonymous mutations are preferentially fixed.

#### Nei-Gojobori (1986) method

**Step 1  -  Count synonymous sites (S) per codon:** For each nucleotide position within a codon, compute the fraction of the 3 possible single-nucleotide substitutions that are synonymous (do not change the amino acid). Sum across the 3 positions to get a per-codon synonymous site count (0-3). Accumulate across all clean codons.

For a pairwise comparison, the synonymous site count is averaged: S_ij = (S_i + S_j) / 2; nonsynonymous sites: N_ij = 3 × L_codon − S_ij.

**Step 2  -  Count observed differences (sd, nd):** For codons differing at k positions, enumerate all k! orderings of the differing positions. For each ordering, trace through the intermediate codons and classify each step as synonymous or nonsynonymous. Paths through stop codons are excluded. Average over all valid orderings (pathway averaging). Sum sd (synonymous differences) and nd (nonsynonymous differences) over all clean codons.

**Step 3  -  Proportions and Jukes-Cantor correction:**

```
pS = sd / S_ij       pN = nd / N_ij

Ks = −(3/4) × ln(1 − 4pS/3)    (JC correction for synonymous sites)
Ka = −(3/4) × ln(1 − 4pN/3)    (JC correction for nonsynonymous sites)
```

If pS ≥ 0.75 or pN ≥ 0.75 (JC saturation threshold), that pair is excluded from averages.

**Step 4  -  Average over pairs:** Report mean Ks, mean Ka, and ω = Ka/Ks across all valid pairs.

**Complete deletion at codon level**: any codon where any nucleotide in any sequence is not in {A, T, C, G} is skipped. Stop codons are skipped.

#### Output statistics

| Statistic | Description |
|---|---|
| n_codons | Number of clean codons used |
| S_sites | Total synonymous site count (summed across all sequences, all codons) |
| N_sites | Total nonsynonymous site count |
| Sd | Mean synonymous differences per pair |
| Nd | Mean nonsynonymous differences per pair |
| Ks | Mean JC-corrected synonymous substitution rate |
| Ka | Mean JC-corrected nonsynonymous substitution rate |
| omega (ω) | Ka/Ks; `None` if Ks = 0 or no valid pairs |

#### CLI usage

```bash
# Ka/Ks on a coding alignment
python skills/dnasp/dnasp.py \
    --input coding_aln.fas \
    --analysis kaks \
    --output results/

# Ka/Ks alongside polymorphism and MK
python skills/dnasp/dnasp.py \
    --input coding_aln_with_outgroup.fas \
    --analysis polymorphism,mk,kaks \
    --outgroup OutSeq \
    --output results/
```

**Caution:** This implementation uses the original Nei-Gojobori (1986) method without Jukes-Cantor correction for multiple hits at individual sites within the codon (only the overall proportion is corrected). For highly divergent sequences (> ~15% difference), the Yang & Nielsen (2000) maximum likelihood method is preferred.

Reference: Nei & Gojobori (1986).

---

### Fu's Fs test

Run with `--analysis fufs`. No extra flags needed  -  Fu's Fs reuses nucleotide diversity (π) and haplotype count (H) already computed by the polymorphism module, which always runs.

#### Biological motivation

Fu (1997) proposed Fs as a test of the neutral model based on the *number* of haplotypes rather than segregating sites alone. Under positive selection or population expansion, the number of haplotypes is lower than expected for the observed nucleotide diversity  -  reflected as a large negative Fs value.

#### Statistical model

Given observed haplotype count H and scaled mutation rate θ_π = k (mean pairwise differences), Fs asks: how probable is it to observe H or fewer haplotypes?

**Ewens sampling formula** (Ewens 1972): under the infinite-alleles model with mutation rate θ, the probability of observing exactly k distinct alleles in a sample of n sequences is:

```
P(K_n = k) = |s(n, k)| × θ^k / θ^(n)
```

where |s(n, k)| = unsigned Stirling number of the first kind, and θ^(n) = θ(θ+1)…(θ+n−1) is the rising factorial (Pochhammer symbol). Unsigned Stirling numbers are computed by the recurrence |s(n+1, k)| = n × |s(n, k)| + |s(n, k−1)| using Python arbitrary-precision integers.

The CDF:

```
S_k = P(K_n ≤ H | θ_π, n) = Σ_{k=1}^{H} P(K_n = k)
```

**Fs statistic:**

```
Fs = ln( S_k / (1 − S_k) )
```

This is a logit transformation of S_k. Fs << 0 means the observed haplotype count is improbably low given the nucleotide diversity.

#### Interpretation

| Result | Interpretation |
|---|---|
| Fs << 0, S_k ≤ 0.02 | Significant  -  far fewer haplotypes than expected; population expansion or positive selection |
| Fs < 0, S_k > 0.02 | Non-significant departure; consistent with neutrality |
| Fs ≈ 0 | Observed haplotypes match neutral expectation |
| Fs > 0 | More haplotypes than expected; balancing selection or subdivision |

**Significance threshold:** Fu (1997) recommended S_k ≤ 0.02 (not the usual 0.05) because Fs is anti-conservative.

**Requirements:** n ≥ 2, H ≥ 1, k > 0. If k = 0 (all sequences identical), θ_π = 0 and the result is degenerate.

#### CLI usage

```bash
# Fu's Fs alongside polymorphism (polymorphism always runs)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis fufs \
    --output results/

# Combined with all standard neutrality tests
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis polymorphism,fufs \
    --output results/
```

References: Fu (1997); Ewens (1972).

---

### Site frequency spectrum

Run with `--analysis sfs`. Folded SFS is always produced; unfolded SFS additionally requires `--outgroup <seq_name>`.

#### Biological motivation

The site frequency spectrum (SFS) records the distribution of allele frequencies across all segregating sites. Its shape is a sensitive indicator of selection and demographic history:

- **Singleton excess (i=1 dominant)** → negative Tajima's D, consistent with population expansion or purifying selection removing non-neutral alleles.
- **Intermediate-frequency excess (flat or hump-shaped)** → positive Tajima's D, consistent with balancing selection.
- **High-frequency derived allele excess (unfolded SFS skewed toward n−1)** → directional selection driving alleles toward fixation.

#### Folded SFS

For each segregating site (after complete deletion of the ingroup), count the number of sequences carrying the **minor allele** (the less-common allele):

```
folded[i] = number of sites where i sequences carry the minor allele
            (i = 1, 2, …, n//2)
```

The folded SFS does not distinguish ancestral from derived alleles  -  it treats both polarities equivalently. It can be computed from ingroup sequences alone.

#### Unfolded SFS (requires outgroup)

With an outgroup sequence, each allele can be polarised:

- **Ancestral allele** = the allele matching the outgroup.
- **Derived allele** = any ingroup allele that differs from the outgroup.

```
unfolded[i] = number of sites where i sequences carry the derived allele
              (i = 1, 2, …, n−1)
```

Gap treatment: complete deletion on ingroup only for folded; both ingroup and outgroup must have clean ATCG bases for the unfolded count.

#### Output

- `report.md`  -  tables of folded and (if available) unfolded SFS counts.
- `figures/sfs.png`  -  bar chart with one panel (folded) or two panels (folded + unfolded).

#### CLI usage

```bash
# Folded SFS only
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis sfs \
    --output results/

# Folded + unfolded (outgroup named "Dmel" in alignment)
python skills/dnasp/dnasp.py \
    --input aln_with_outgroup.fas \
    --outgroup Dmel \
    --analysis sfs \
    --output results/

# SFS + Fu's Fs + polymorphism in one run
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis polymorphism,fufs,sfs \
    --output results/
```

---

### Transition/Transversion ratio

Run with `--analysis tstv`. Works on any alignment  -  coding or non-coding. Ts/Tv is computed across all n(n−1)/2 pairwise comparisons at every clean (non-gap, unambiguous ATCG) column.

#### Substitution classification

| Change | Type | Reason |
|--------|------|--------|
| A ↔ G | Transition (Ts) | Purine ↔ Purine |
| C ↔ T | Transition (Ts) | Pyrimidine ↔ Pyrimidine |
| A ↔ C | Transversion (Tv) | Purine ↔ Pyrimidine |
| A ↔ T | Transversion (Tv) | Purine ↔ Pyrimidine |
| G ↔ C | Transversion (Tv) | Purine ↔ Pyrimidine |
| G ↔ T | Transversion (Tv) | Purine ↔ Pyrimidine |

#### Formulas

```
n_pairs = n(n−1)/2

n_transitions   = Σ over all pairs and all clean columns: 1 if (same chemical class and a ≠ b)
n_transversions = Σ over all pairs and all clean columns: 1 if (different chemical class)

ts_tv       = n_transitions / n_transversions   (None if n_transversions = 0)
ts_per_site = n_transitions  / (n_pairs × L_net)
tv_per_site = n_transversions / (n_pairs × L_net)
```

#### Interpretation

| Ts/Tv | Interpretation |
|-------|----------------|
| ≈ 2 | Typical nuclear DNA  -  transitions are inherently more mutable |
| > 10 | Mitochondrial DNA or highly conserved regions |
| < 0.5 | Possible saturation of transitions (high divergence) or non-neutral substitution patterns |
| None | No transversions observed (all differences are transitions)  -  common for closely related sequences |

#### Output

- `report.md`  -  table of Ts, Tv, Ts/Tv, ts_per_site, tv_per_site, L_net.
- No separate figure; pair with `ld` or `divergence` for a comprehensive substitution analysis.

#### Example

```bash
# Transition/transversion ratio only
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis tstv \
    --output results/

# Combined coding evolution panel
python skills/dnasp/dnasp.py \
    --input coding.fas \
    --outgroup OutSeq \
    --analysis kaks,tstv,codon \
    --output results/
```

---

### Codon usage bias

Run with `--analysis codon`. The alignment must be an in-frame coding sequence (position 0 must be the first codon position). Stop codons and triplets containing gap or ambiguous characters are skipped automatically.

#### RSCU  -  Relative Synonymous Codon Usage (Sharp & Li 1987)

RSCU measures observed vs expected codon frequency within each synonymous family:

```
RSCU_ij = X_ij / (X_i / n_i)
```

where:
- X_ij = count of codon j for amino acid i (pooled across all sequences)
- X_i  = total codon count for amino acid i = Σⱼ X_ij
- n_i  = synonymous family size (number of codons encoding amino acid i)

| RSCU value | Meaning |
|------------|---------|
| 1.0 | Uniform usage  -  this codon is used at the expected frequency |
| > 1.0 | Preferred codon  -  used more often than expected by chance |
| < 1.0 | Avoided codon  -  used less often than expected |
| 0.0 | Codon not observed in the alignment |

The sum of RSCU values across all codons encoding a given amino acid equals the synonymous family size (e.g. for a 4-fold degenerate amino acid, ΣRSCU = 4).

#### ENC  -  Effective Number of Codons (Wright 1990)

ENC measures the overall strength of codon usage bias across the entire gene:

```
ENC = 2 + 9/F̄₂ + 1/F̄₃ + 5/F̄₄ + 3/F̄₆
```

where F̄_k is the mean corrected homozygosity over all amino acids with k-fold degeneracy:

```
F̄_k = mean over amino acids with k synonymous codons of:
       F_aa = (n_aa × Σⱼ pⱼ² − 1) / (n_aa − 1)
```

with pⱼ = X_ij / X_i (observed codon frequency within the family) and n_aa = X_i (total codon observations for that amino acid). Amino acids with n_aa < 2 are excluded from the average for that degeneracy class.

Coefficients in the ENC formula (9, 1, 5, 3) equal the number of amino acids in each degeneracy class under the standard genetic code:

| Degeneracy class | Amino acids | Count |
|-----------------|-------------|-------|
| 2-fold | Cys, Asp, Glu, Phe, His, Lys, Asn, Gln, Tyr | 9 |
| 3-fold | Ile | 1 |
| 4-fold | Ala, Gly, Pro, Thr, Val | 5 |
| 6-fold | Arg, Leu, Ser | 3 |

Met (ATG only) and Trp (TGG only) have no synonymous codons and are excluded.

ENC is clamped to the range [20, 61]:

| ENC | Interpretation |
|-----|----------------|
| 61 | No codon usage bias  -  all synonymous codons used equally |
| 50-60 | Weak bias |
| 35-50 | Moderate bias |
| < 35 | Strong codon usage bias  -  typical of highly expressed genes |
| 20 | Maximum bias  -  only one codon used per amino acid |
| None | Insufficient data; one or more degeneracy classes has no usable amino acids |

#### Output

- `report.md`  -  ENC value, n_codons, and full RSCU table grouped by amino acid.
- `figures/codon_usage.png`  -  bar chart of RSCU values for all 61 sense codons, colour-coded by amino acid family, with RSCU = 1.0 reference line and ENC annotation in the title.

#### Requirements

- Alignment length must be divisible by 3 and in-frame from position 0.
- ENC requires n_aa ≥ 2 observations for at least one amino acid in each of the four degeneracy classes. For short genes (< 300 bp) ENC = None is common.
- For reliable ENC estimates, use > 300 bp of clean coding sequence.

#### Example

```bash
# Codon usage bias alone
python skills/dnasp/dnasp.py \
    --input coding.fas \
    --analysis codon \
    --output results/

# Full coding evolution panel (Ka/Ks + MK + Ts/Tv + Codon usage)
python skills/dnasp/dnasp.py \
    --input coding.fas \
    --outgroup OutSeq \
    --analysis kaks,mk,tstv,codon \
    --output results/
```

---

### Fay & Wu's H and Zeng's E

Run with `--analysis faywu --outgroup SEQNAME`. The named outgroup sequence is extracted from the alignment and used to polarise mutations (determine which allele is ancestral vs derived). It is not included in ingroup calculations.

#### Biological motivation

Fay & Wu (2000) proposed the H statistic as a test for the signature of a selective sweep. A hard sweep drives a positively selected allele to high frequency, dragging linked variants along with it. In the genealogy, this creates an excess of *high-frequency derived alleles*  -  variants where the derived allele is present in n−1 or n−2 sequences. H < 0 signals this pattern. Zeng et al. (2006) extended the framework with E, which detects *low-frequency derived alleles* (n=1 or n=2 copies) and is more sensitive to certain demographic effects.

#### Unfolded SFS

Both H and E are computed from the **unfolded site frequency spectrum** (ξ), which requires an outgroup to polarise each segregating site. For each alignment column (after complete deletion across all ingroup sequences *and* the outgroup):

1. Identify the ancestral allele as the nucleotide carried by the outgroup.
2. Count *i* = number of ingroup sequences carrying the derived allele at that site.
3. Increment ξ[i] by 1.

Sites where the outgroup allele is not in {A, T, C, G}, or where any ingroup sequence is not in {A, T, C, G}, or where the outgroup allele is not present in the ingroup at all, are skipped.

#### Estimators

All four θ estimators are per-site (divided by L_net = number of clean columns):

```
a₁ = Σ_{k=1}^{n-1} 1/k   (a₁ = Watterson normalisation constant)

θ_W = n_polarised / (a₁ × L_net)     (Watterson  -  uses count of polarised sites)

θ_π = [Σ_{i=1}^{n-1} ξ_i × i(n−i)] / [n(n−1)/2 × L_net]    (nucleotide diversity from SFS)

θ_H = [Σ_{i=1}^{n-1} ξ_i × 2i²] / [n(n−1) × L_net]          (Fay & Wu 2000, eq. 2)

θ_L = [Σ_{i=1}^{n-1} ξ_i × i] / [(n−1) × L_net]              (Zeng et al. 2006)
```

#### Test statistics

```
H = θ_π − θ_H      (Fay & Wu 2000)
E = θ_L − θ_W      (Zeng et al. 2006)
```

#### Interpretation

| Statistic | Sign | Interpretation |
|---|---|---|
| H < 0 | | Excess high-frequency derived alleles → hard selective sweep signature |
| H > 0 | | Deficit of high-frequency derived alleles  -  uncommon in practice |
| H ≈ 0 | | Consistent with neutrality |
| E < 0 | | Excess low-frequency derived alleles → rapid expansion or selection against recessive alleles |
| E > 0 | | Deficit of low-frequency variants → balancing selection or bottleneck |
| E ≈ 0 | | Consistent with neutrality |

H and E are complementary: a selective sweep produces H << 0 but may produce E close to 0, whereas population expansion produces E << 0 but H close to 0. Used together they help distinguish selective from demographic explanations.

**Requirements:** n ≥ 2, `--outgroup` provided, at least one polarisable segregating site. Returns `None` values if no polarisable sites are found.

#### CLI usage

```bash
# Fay & Wu H + Zeng E (outgroup named "Dmel" in alignment)
python skills/dnasp/dnasp.py \
    --input dsim_with_dmel.fas \
    --outgroup Dmel \
    --analysis faywu \
    --output results/

# Combined with polymorphism and SFS
python skills/dnasp/dnasp.py \
    --input dsim_with_dmel.fas \
    --outgroup Dmel \
    --analysis polymorphism,sfs,faywu \
    --output results/
```

References: Fay & Wu (2000); Zeng et al. (2006).

---

### Population differentiation  -  Fst

Run with `--analysis fst`. Requires either `--pop-file` (population assignment TSV) or `--input2` (two separate FASTA files). Computes the Hudson et al. (1992) pairwise Fst estimator for all pairs of populations.

#### Biological motivation

Fst (fixation index) measures the degree to which populations differ in allele frequencies relative to total genetic diversity. High Fst between two populations indicates strong differentiation  -  either due to a long history of isolation, local adaptation, or genetic drift in a small population. Low Fst indicates similarity, consistent with ongoing gene flow or recent divergence.

#### Hudson et al. (1992) estimator

For a pair of populations with sequences **a** and **b**:

```
π_a   = mean pairwise nucleotide differences per site within population a
π_b   = mean pairwise nucleotide differences per site within population b
π_s   = mean within-population diversity = (π_a + π_b) / 2
π_t   = Dxy = mean between-population divergence per site (cross-population pairs)

Fst = 1 − π_s / π_t
```

This formulation (Hudson et al. 1992, eq. 3) is analogous to the classical Fst = (H_T − H_S) / H_T but uses nucleotide diversity as the measure of heterozygosity.

All calculations use **complete deletion**: any column where any sequence (in any population) carries a non-ATCG character is excluded. L_net is reported in the output.

Fst is clamped to [0, 1]: negative values (which can arise when π_s > π_t due to sampling with small n) are set to 0. If π_t = 0 (no between-population variation at any site), Fst is reported as `None`.

For analyses with more than 2 populations, all pairwise Fst values are computed and a mean Fst across all pairs is reported.

#### Interpretation (Wright 1978 thresholds)

| Fst | Level of differentiation |
|-----|--------------------------|
| 0 - 0.05 | Little or negligible differentiation |
| 0.05 - 0.15 | Moderate differentiation |
| 0.15 - 0.25 | Great differentiation |
| > 0.25 | Very great differentiation |

These thresholds are commonly used guidelines, not hard boundaries. Fst interpretation depends heavily on the organism, marker type, and geographic scale.

#### Output statistics (FstStats dataclass)

| Field | Description |
|-------|-------------|
| `n_pops` | Number of populations |
| `pop_names` | List of population names |
| `pop_sizes` | Dict: population → number of sequences |
| `pi_within` | Dict: population → π per site |
| `pi_between` | Dict: (pop1, pop2) → Dxy |
| `fst_pairwise` | Dict: (pop1, pop2) → Fst (or None) |
| `fst_mean` | Mean Fst across all pairwise comparisons |

#### CLI usage

```bash
# Fst from a single alignment + pop-file
python skills/dnasp/dnasp.py \
    --input combined.fas \
    --pop-file populations.txt \
    --analysis fst \
    --output results/

# Fst from two separate FASTA files
python skills/dnasp/dnasp.py \
    --input pop1.fas --input2 pop2.fas \
    --analysis fst \
    --output results/

# Full population panel: divergence + Fst + Fay & Wu
python skills/dnasp/dnasp.py \
    --input combined_with_outgroup.fas \
    --pop-file populations.txt \
    --outgroup OutSeq \
    --analysis divergence,fst,faywu \
    --output results/
```

**Output figure:** `figures/fst.png`  -  bar chart of pairwise Fst values with horizontal reference lines at Wright's 0.05, 0.15, and 0.25 thresholds.

References: Hudson et al. (1992); Wright (1978).

---

## Output files

Running the skill produces the following files in the output directory:

```
output_directory/
├── report.md                  # Full narrative Markdown report
├── results.tsv                # DnaSP-compatible tab-delimited table
├── ld_pairs.tsv               # Pairwise LD table (only when --analysis ld)
├── figures/
│   ├── summary.png            # Bar chart: π, θ_W, Hd  (requires matplotlib)
│   ├── sliding_window.png     # π and Tajima's D per window (if --window)
│   ├── ld_decay.png           # R² vs distance scatter (if --analysis ld)
│   ├── mismatch.png           # Mismatch distribution bar chart (if --analysis popsize)
│   ├── sfs.png                # Site frequency spectrum (if --analysis sfs)
│   ├── codon_usage.png        # RSCU bar chart per codon (if --analysis codon)
│   └── fst.png                # Pairwise Fst bar chart with Wright thresholds (if --analysis fst)
└── reproducibility/
    ├── commands.sh            # Exact command used (rerun-ready)
    ├── environment.yml        # Python + package versions
    └── checksums.sha256       # SHA-256 of input file and all outputs
```

### report.md

Markdown narrative including:
- Input file(s), date, sequence summary (n, L, L\_net, GC)
- Polymorphism table (always)
- Neutrality tests table (always)
- LD summary + top pairs (if `--analysis ld`)
- Recombination summary (if `--analysis recombination`)
- Mismatch distribution (if `--analysis popsize`)
- InDel summary (if `--analysis indel`)
- Divergence summary (if `--analysis divergence`)
- Fu & Li D/F with outgroup (if `--analysis fuliout`)
- HKA test summary with per-locus breakdown (if `--analysis hka`)
- McDonald-Kreitman table (Pn, Ps, Dn, Ds) and derived statistics (if `--analysis mk`)
- Ka/Ks summary (S_sites, N_sites, Ks, Ka, ω) (if `--analysis kaks`)
- Fu's Fs table (θ_π, H, S_k, Fs) with significance note (if `--analysis fufs`)
- SFS tables  -  folded and unfolded (if `--analysis sfs`)
- Transition/Transversion table (Ts, Tv, Ts/Tv, per-site rates, L_net) (if `--analysis tstv`)
- Codon Usage section  -  ENC, n_codons, and full RSCU table grouped by amino acid (if `--analysis codon`)
- Fay & Wu section  -  θ_W, θ_π, θ_H, θ_L, H, E values; interpretation of H and E sign (if `--analysis faywu`)
- Population Fst section  -  within-population π, pairwise Dxy, pairwise Fst table, Wright-threshold interpretation, mean Fst (if `--analysis fst`)

### results.tsv

Tab-delimited file with the same column order as DnaSP 6's `.MF.out` output. Readable in Excel or R:

```r
df <- read.table("results.tsv", header=TRUE, sep="\t")
```

### ld\_pairs.tsv

One row per biallelic site pair, columns: `site1`, `site2`, `dist`, `n_valid`, `D`, `D_prime`, `R2`, `chi2`, `p_chi2`. Produced only when `--analysis ld` is active.

### figures/

- `summary.png`  -  bar chart: π, θ\_W, Hd. Requires `matplotlib ≥ 3.7`.
- `sliding_window.png`  -  π (blue) and Tajima's D (orange) per window. Only when `--window` is specified.
- `ld_decay.png`  -  R² vs physical distance (bp) scatter plot. Only when `--analysis ld`.
- `mismatch.png`  -  mismatch distribution histogram. Only when `--analysis popsize`.
- `sfs.png`  -  folded SFS bar chart (+ unfolded panel if `--outgroup`). Only when `--analysis sfs`.
- `codon_usage.png`  -  RSCU values for all 61 sense codons, colour-coded by amino acid family. Only when `--analysis codon`.
- `fst.png`  -  pairwise Fst bar chart with Wright (1978) threshold lines at 0.05, 0.15, 0.25. Only when `--analysis fst`.

### reproducibility/

- `commands.sh`  -  the exact invocation used, so the analysis can be re-run identically
- `environment.yml`  -  Python version and package versions at analysis time
- `checksums.sha256`  -  SHA-256 digests of all output files and the input alignment, for archiving

---

## Sliding-window analysis

Use `--window` and (optionally) `--step` to compute polymorphism statistics in overlapping windows across the alignment. Each window uses the same complete-deletion gap treatment as the whole-alignment analysis. Windows with fewer than 3 sequences after gap removal, or with S = 0, report `n.a.` for tests that require variation.

```bash
# 100 bp windows, 25 bp step
python skills/dnasp/dnasp.py \
  --input alignment.fas \
  --window 100 --step 25 \
  --output results_window/
```

The `results.tsv` in sliding-window mode contains one row per window, with columns for window start, end, midpoint, L\_net, S, π, θ\_W, Tajima's D, D\*, F\*, R2. The LD, recombination, mismatch, InDel, and divergence modules are not currently computed per window (whole-alignment only).

---

## Implementation notes and gotchas

### Complete deletion vs. pairwise deletion

DnaSP applies complete deletion by default: any site with a gap or missing base in *any* sequence is excluded from the entire analysis. This skill matches that behaviour. If your alignment is very gap-rich, L\_net may be much smaller than L\_total  -  always check both values in the report. For divergence, the gap mask is applied across both populations combined.

### Singletons: count == 1 only (no folded counting)

DnaSP's `BusqSingletones()` counts a singleton as any allele present in exactly 1 sequence. It does **not** count alleles present in n−1 sequences as additional singletons (folded definition). For a biallelic site with counts (1, n−1) only one singleton is recorded. Using the folded definition would approximately double η\_s and produce incorrect D\* and F\* values.

### n < 3 limitation

Tajima's D, Fu & Li D\* and F\* are mathematically undefined for n < 3 (variance formulas involve `n − 2`). The script returns `n.a.` and prints a warning.

### Eta vs. S for Fu & Li tests

DnaSP uses Eta in the D\* numerator and S in the variance denominator. For strictly biallelic datasets Eta = S. For triallelic or quadriallelic sites Eta > S. The implementation follows DnaSP's convention exactly.

### F\* discrepancy with DnaSP multi-MSA output

DnaSP's RAD-seq module uses Achaz (2009) variances for F\*. The standard FASTA/NEXUS analysis uses Simonsen (1995) A5-A6. This skill implements the standard mode. D\* matches exactly regardless of mode.

### D' inflation with small samples

D' is known to be upwardly biased (inflated toward ±1) when sample size is small. Always report D' alongside n. With n < 10, treat D' values cautiously; prefer R² for comparisons.

### Rm is a lower bound

The four-gamete test gives the *minimum* number of past recombination events consistent with the data. The true count is ≥ Rm. Rm = 0 does not mean no recombination occurred; it means the data do not require any.

### Mismatch raggedness significance

The skill computes the raggedness statistic r but does not currently provide a p-value (which requires parametric simulation under a specific demographic model). r is best interpreted qualitatively: values < 0.05 are conventionally taken as consistent with expansion.

### Da can be negative

Net divergence Da = Dxy − (π₁ + π₂)/2 can be slightly negative due to sampling error, especially when L_net is small. Treat negative Da as effectively zero; do not interpret it as evidence for anything.

### FASTA wrapping and DnaSP headers

DnaSP exports sequences wrapped across multiple lines with `>'name'  [comment]` headers. Standard parsers may fail on these. Always use `dnasp.py` to parse DnaSP-generated files.

### NEXUS MATCHCHAR expansion

NEXUS files from DnaSP use `.` as MATCHCHAR. The parser expands dots relative to the **first** sequence. If your first sequence is an outgroup, re-order before analysis.

---

## Roadmap

### v0.5.0 (planned)

**VCF support.** A new `parse_vcf_population` function for `clawbio/common/parsers.py` will accept multi-sample population VCF files and return aligned haplotype sequences compatible with all DnaSP analyses. Design proposal filed with the ClawBio team. This will add `--vcf` as an alternative input format alongside `--fasta`.

**Simulation-based p-values.** Kingman coalescent simulations (pure Python, no new dependencies) will provide empirical p-values for Tajima's D (two-tailed) and R2 (left-tail). Opt-in via `--n-sim INT` (default 0); reproducible with `--sim-seed INT`. Design proposal filed with the ClawBio team.

---

## Citations

If you use results from this skill in a publication, cite the relevant references for each module used.

**DnaSP 6 (the software this skill reimplements):**
> Rozas J, Ferrer-Mata A, Sánchez-DelBarrio JC, Guirao-Rico S, Librado P, Ramos-Onsins SE, Sánchez-Gracia A (2017). DnaSP 6: DNA Sequence Polymorphism Analysis of Large Data Sets. *Mol. Biol. Evol.* 34:3299-3302. https://doi.org/10.1093/molbev/msx248

**Tajima's D:**
> Tajima F (1989). Statistical method for testing the neutral mutation hypothesis by DNA polymorphism. *Genetics* 123:585-595.

**Fu & Li D\* and F\*:**
> Fu YX, Li WH (1993). Statistical tests of neutrality of mutations. *Genetics* 133:693-709.

**Variance coefficients for D\* and F\* (Simonsen et al.):**
> Simonsen KL, Churchill GA, Aquadro CF (1995). Properties of statistical tests of neutrality for DNA polymorphism data. *Genetics* 141:413-429.

**Haplotype diversity:**
> Nei M, Tajima F (1981). DNA polymorphism detectable by restriction endonucleases. *Genetics* 97:145-163.

**Nucleotide diversity:**
> Nei M (1987). *Molecular Evolutionary Genetics*. Columbia University Press, New York.

**Watterson's theta:**
> Watterson GA (1975). On the number of segregating sites in genetical models without recombination. *Theor. Popul. Biol.* 7:256-276.

**Ramos-Onsins & Rozas R2:**
> Ramos-Onsins SE, Rozas J (2002). Statistical properties of new neutrality tests against population growth. *Mol. Biol. Evol.* 19:2092-2100. https://doi.org/10.1093/oxfordjournals.molbev.a004068

**Linkage disequilibrium (D):**
> Lewontin RC, Kojima K (1960). The evolutionary dynamics of complex polymorphisms. *Evolution* 14:458-472.

**D' (normalised LD):**
> Lewontin RC (1964). The interaction of selection and linkage. I. General considerations; heterotic models. *Genetics* 49:49-67.

**R² (LD):**
> Hill WG, Robertson A (1968). Linkage disequilibrium in finite populations. *Theor. Appl. Genet.* 38:226-231.

**ZnS:**
> Kelly JK (1997). A test of neutrality based on interlocus associations. *Genetics* 146:1197-1206.

**Za, ZZ:**
> Rozas J, Gullaud M, Blandin G, Aguadé M (2001). DNA variation at the *rp49* gene region of *Drosophila simulans*: evolutionary inferences from an unusual haplotype structure. *Genetics* 158:1321-1330.

**Rm (minimum recombination events):**
> Hudson RR, Kaplan NL (1985). Statistical properties of the number of recombination events in the history of a sample of DNA sequences. *Genetics* 111:147-164.

**Mismatch distribution:**
> Rogers AR, Harpending H (1992). Population growth makes waves in the distribution of pairwise genetic differences. *Mol. Biol. Evol.* 9:552-569.

**Raggedness statistic:**
> Harpending HC (1994). Signature of ancient population growth in a low-resolution mitochondrial DNA mismatch distribution. *Hum. Biol.* 66:591-600.

**Dxy, Da (between-population divergence):**
> Nei M (1987). *Molecular Evolutionary Genetics*. Columbia University Press, New York. Equation 10.20.

**Fixed/shared/private classification:**
> Hey J (1991). A multi-dimensional coalescent process applied to multi-allelic selection models and migration models. *Theor. Popul. Biol.* 39:30-48.

> Hey J, Wakeley J (1997). A coalescent estimator of the population recombination rate. *Genetics* 145:833-846.

**Fu & Li D/F with outgroup:**
> Fu YX, Li WH (1993). Statistical tests of neutrality of mutations. *Genetics* 133:693-709.

> Simonsen KL, Churchill GA, Aquadro CF (1995). Properties of statistical tests of neutrality for DNA polymorphism data. *Genetics* 141:413-429.

**HKA test:**
> Hudson RR, Kreitman M, Aguadé M (1987). A test of neutral molecular evolution based on nucleotide data. *Genetics* 116:153-159.

**McDonald-Kreitman test:**
> McDonald JH, Kreitman M (1991). Adaptive protein evolution at the Adh locus in Drosophila. *Nature* 351:652-654. https://doi.org/10.1038/351652a0

**Fu's Fs test:**
> Fu YX (1997). Statistical tests of neutrality of mutations against population growth, hitchhiking and background selection. *Genetics* 147:915-925.

**Ewens sampling formula (basis for Fu's Fs):**
> Ewens WJ (1972). The sampling theory of selectively neutral alleles. *Theor. Popul. Biol.* 3:87-112. https://doi.org/10.1016/0040-5809(72)90035-4

**Ka/Ks  -  synonymous sites and substitution rates:**
> Nei M, Gojobori T (1986). Simple methods for estimating the numbers of synonymous and nonsynonymous nucleotide substitutions. *Mol. Biol. Evol.* 3:418-426. https://doi.org/10.1093/oxfordjournals.molbev.a040410

**Transition/transversion ratio:**
> (Classic biochemistry; no single definitive citation  -  standard practice in molecular evolution textbooks, e.g. Li WH (1997) *Molecular Evolution*. Sinauer Associates.)

**Codon usage bias  -  RSCU:**
> Sharp PM, Li WH (1987). The codon adaptation index  -  a measure of directional synonymous codon usage bias, and its potential applications. *Nucleic Acids Res.* 15:1281-1295. https://doi.org/10.1093/nar/15.3.1281

**Codon usage bias  -  ENC (Effective Number of Codons):**
> Wright F (1990). The 'effective number of codons' used in a gene. *Gene* 87:23-29. https://doi.org/10.1016/0378-1119(90)90491-9

**Fay & Wu's H:**
> Fay JC, Wu CI (2000). Hitchhiking under positive Darwinian selection. *Genetics* 155:1405-1413. https://doi.org/10.1093/genetics/155.3.1405

**Zeng's E (θ_L-based test):**
> Zeng K, Fu YX, Shi S, Wu CI (2006). Statistical tests for detecting positive selection by utilizing high-frequency variants. *Genetics* 174:1431-1439. https://doi.org/10.1534/genetics.106.061432

**Population differentiation Fst  -  Hudson et al. (1992):**
> Hudson RR, Slatkin M, Maddison WP (1992). Estimation of levels of gene flow from DNA sequence data. *Genetics* 132:583-589.

**Wright (1978) Fst interpretation thresholds:**
> Wright S (1978). *Evolution and the Genetics of Populations*, Vol. 4: Variability Within and Among Natural Populations. University of Chicago Press, Chicago.

**Achaz (2009) variance (DnaSP v6 RAD-seq mode, for reference only):**
> Achaz G (2009). Frequency spectrum neutrality tests: one for all and all for one. *Genetics* 183:249-258.
