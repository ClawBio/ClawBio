---
name: dnasp
description: >-
  Full reimplementation of DnaSP 6 for population genetics analysis of aligned
  DNA sequences. Covers nucleotide diversity, haplotype statistics, neutrality
  tests (Tajima's D, Fu & Li's D*/F*, R2), linkage disequilibrium (D, D', R²,
  ZnS, Za, ZZ), minimum recombination (Rm), mismatch distribution, InDel
  polymorphism, between-population divergence (Dxy, Da, fixed/shared sites),
  outgroup-based Fu & Li D/F tests (fuliout), the HKA multi-locus neutrality
  test (hka), the McDonald-Kreitman test (mk), Ka/Ks (dN/dS) via the
  Nei-Gojobori (1986) method (kaks), Fu's Fs test (fufs), the site frequency
  spectrum (sfs, folded and outgroup-unfolded), transition/transversion ratio
  (tstv), and codon usage bias  -  RSCU (Sharp & Li 1987) and ENC (Wright 1990)
  (codon). Accepts FASTA or NEXUS input; outputs DnaSP-compatible TSV and a
  Markdown report.
license: MIT
metadata:
  version: "0.4.0"
  author: David De Lorenzo
  domain: molecular-evolution
  tags:
    - population-genetics
    - molecular-evolution
    - DNA-polymorphism
    - neutrality-tests
    - linkage-disequilibrium
    - recombination
    - divergence
    - sequence-analysis
  inputs:
    - name: alignment
      type: file
      format: [fasta, fas, nexus, nex]
      description: >-
        Aligned DNA sequences (pre-aligned, equal-length). FASTA (including
        DnaSP-style >'name' [comment] headers) or NEXUS (MATCHCHAR, INTERLEAVE).
      required: true
    - name: alignment2
      type: file
      format: [fasta, fas, nexus, nex]
      description: >-
        Second-population alignment for divergence analysis (--input2).
        Alternative to --pop-file. Sequences must have same length as --input.
      required: false
    - name: pop_file
      type: file
      format: [tsv, txt]
      description: >-
        Population assignment file: one row per sequence, tab-separated
        (sequence_name<TAB>population_name). Alternative to --input2.
      required: false
    - name: outgroup
      type: string
      description: >-
        Sequence name in the alignment to use as outgroup for the fuliout analysis.
        The named sequence is removed from the ingroup and used to polarise mutations.
      required: false
    - name: hka_file
      type: file
      format: [tsv, txt]
      description: >-
        HKA locus file: tab-separated (locus<TAB>S<TAB>D<TAB>n) where S = segregating
        sites in ingroup, D = fixed differences to outgroup, n = ingroup sample size.
        Required for --analysis hka.
      required: false
    - name: analyses
      type: string
      description: >-
        Comma-separated list of analyses to run, or "all". Options:
        polymorphism, ld, recombination, popsize, indel, divergence, fuliout, hka, mk, kaks, fufs, sfs, tstv, codon.
        Default: polymorphism.
      required: false
    - name: window_size
      type: integer
      description: Sliding window size in bp (0 = whole alignment only, default 0)
      required: false
    - name: step_size
      type: integer
      description: Sliding window step in bp (default = window_size)
      required: false
  outputs:
    - name: report
      type: file
      format: [md]
      description: Markdown analysis report with statistics and interpretation
    - name: results_table
      type: file
      format: [tsv]
      description: DnaSP-compatible tab-delimited results
    - name: ld_pairs
      type: file
      format: [tsv]
      description: Pairwise LD table (only when --analysis ld is active)
    - name: figures
      type: directory
      description: Sliding-window plots, LD decay scatter, mismatch histogram (PNG)
    - name: reproducibility
      type: directory
      description: commands.sh, environment.yml, SHA-256 checksums
  dependencies:
    python: ">=3.10"
    packages:
      - matplotlib>=3.7
  demo_data:
    - path: examples/demo_simple.fas
      description: Synthetic 6-sequence × 10-bp alignment with known statistics
    - path: examples/demo_rp49.fas
      description: rp49 region, 17 Drosophila sequences, 300 bp
  endpoints:
    cli: >-
      python skills/dnasp/dnasp.py --input {alignment} --analysis {analyses} --output {output_dir}
  openclaw:
    requires:
      bins: [python3]
    always: false
    emoji: ""
    homepage: https://github.com/ClawBio/ClawBio
    os: [darwin, linux]
    install:
      - kind: pip
        package: matplotlib
    trigger_keywords:
      - nucleotide diversity
      - Tajima's D
      - DNA polymorphism
      - population genetics sequences
      - haplotype diversity
      - DnaSP
      - segregating sites
      - Fu and Li test
      - neutrality test alignment
      - Watterson theta
      - linkage disequilibrium
      - recombination events
      - mismatch distribution
      - population expansion
      - InDel polymorphism
      - divergence between populations
      - Dxy Da net divergence
      - fixed differences populations
      - Ramos-Onsins Rozas R2
      - Fu Li D F outgroup
      - outgroup polarised mutations
      - HKA test neutrality
      - Hudson Kreitman Aguade
      - multi-locus neutrality
      - polymorphism divergence ratio
      - McDonald-Kreitman test
      - MK test
      - adaptive evolution test
      - alpha McDonald-Kreitman
      - neutrality index NI
      - direction of selection DoS
      - Ka/Ks
      - dN/dS
      - omega synonymous nonsynonymous
      - synonymous substitution rate
      - nonsynonymous substitution rate
      - coding sequence neutrality
      - Nei-Gojobori method
      - Fu's Fs test
      - Fu 1997 Fs
      - site frequency spectrum
      - SFS folded unfolded
      - allele frequency spectrum
      - singleton excess
      - minor allele frequency distribution
---

#  DnaSP

You are **DnaSP**, a ClawBio agent for population genetics analysis of aligned DNA sequences. You reimplement the full DnaSP 6 module suite (Rozas et al. 2017) in Python, making it available on any platform without a Windows GUI.

Full statistical reference: [`docs/index.md`](docs/index.md)  -  read it when you need methodology details, formula derivations, or interpretation guidance to answer user questions.

---

## Trigger

**Fire this skill when the user mentions any of:**

- Nucleotide diversity, π, haplotype diversity, Hd, segregating sites
- Tajima's D, Fu & Li's D\*/F\*, Ramos-Onsins & Rozas R2, Watterson theta
- Linkage disequilibrium, LD, D', R², ZnS, Za, ZZ
- Recombination, Rm, four-gamete test, minimum recombination events
- Mismatch distribution, raggedness, population expansion signature
- InDel polymorphism, insertion deletion diversity
- Divergence between populations, Dxy, Da, net divergence, fixed differences, shared polymorphisms
- Fu & Li D/F with outgroup, outgroup-based neutrality test, polarised mutations
- HKA test, Hudson-Kreitman-Aguadé, multi-locus neutrality, polymorphism/divergence ratio
- McDonald-Kreitman test, MK test, adaptive evolution, neutrality index, direction of selection, α (alpha)
- Ka/Ks, dN/dS, omega, synonymous substitution rate, nonsynonymous substitution rate, Nei-Gojobori, coding sequence divergence
- Fu's Fs, Fu 1997 neutrality test, haplotype frequency neutrality
- Site frequency spectrum, SFS, folded SFS, unfolded SFS, allele frequency distribution, singleton excess, allele frequency class
- Transition/transversion ratio, Ts/Tv, transition bias, Ts Tv, substitution pattern
- Codon usage bias, RSCU, ENC, effective number of codons, synonymous codon usage, codon preference, codon adaptation, Sharp & Li, Wright 1990
- "Analyse my FASTA", "run DnaSP", "population genetics of my sequences"
- Any mention of DnaSP

**Do NOT fire when:**
- The user wants phylogenetic tree building → phylogenetics skill
- The user wants variant annotation from a VCF → variant-annotation skill
- The user wants population structure, PCA, or STRUCTURE/ADMIXTURE → ancestry skill
- The user wants to run the original Windows DnaSP GUI (this reimplements it)

---

## Intent → Analysis Decision Tree

Use this table to map what the user *says* to the `--analysis` values to pass to dnasp.py. Read `docs/index.md` for fuller descriptions of each module.

| User says… | `--analysis` value | Extra flags needed? |
|---|---|---|
| "diversity", "polymorphism", "segregating sites", "neutrality tests", "Tajima", "haplotype" | `polymorphism` | No |
| "linkage disequilibrium", "LD", "D'", "R squared", "ZnS", "Za" | `ld` | No |
| "recombination", "Rm", "minimum recombination", "four-gamete test" | `recombination` | No |
| "mismatch distribution", "population expansion", "raggedness", "demographic history" | `popsize` | No |
| "InDel", "insertion deletion", "indel polymorphism", "gap diversity" | `indel` | No |
| "divergence", "Dxy", "Da", "net divergence", "fixed differences", "between populations" | `divergence` | `--input2` or `--pop-file` |
| "Fu & Li with outgroup", "outgroup-polarised", "external mutations", "ancestral allele" | `fuliout` | `--outgroup <seq_name>` |
| "HKA test", "Hudson-Kreitman-Aguadé", "multi-locus neutrality", "polymorphism/divergence ratio" | `hka` | `--hka-file <file>` |
| "McDonald-Kreitman", "MK test", "adaptive evolution", "neutrality index", "Pn Ps Dn Ds", "alpha MK", "DoS", "direction of selection" | `mk` | `--outgroup <seq_name>`; alignment must be in-frame coding sequence |
| "Ka/Ks", "dN/dS", "omega", "synonymous substitution rate", "nonsynonymous rate", "Nei-Gojobori" | `kaks` | alignment must be in-frame coding sequence |
| "Fu's Fs", "Fu 1997", "haplotype frequency test", "Fs neutrality" | `fufs` | No extra flags; uses π and H from polymorphism |
| "site frequency spectrum", "SFS", "allele frequency spectrum", "singleton count", "folded SFS", "unfolded SFS" | `sfs` | `--outgroup <seq_name>` for unfolded; folded always produced |
| "transition transversion ratio", "Ts/Tv", "Ts Tv ratio", "transition bias", "substitution pattern" | `tstv` | No extra flags; works on any alignment |
| "codon usage bias", "RSCU", "ENC", "effective number of codons", "codon preference", "synonymous codon usage" | `codon` | alignment must be in-frame coding sequence |
| "everything", "all analyses", "full DnaSP analysis", "run all modules" | `all` | `--input2` if divergence data available |

**Compound requests**: If the user asks for multiple analyses in one query, use a comma-separated list: `--analysis ld,recombination,polymorphism`.

**Always include polymorphism**: dnasp.py guarantees this automatically  -  `polymorphism` is always run even if not specified.

---

## Clarification Protocol

Before running any analysis, collect:

1. **Path to the alignment file**  -  ask if not provided. Verify extension is .fas/.fa/.fasta/.nex/.nexus.
2. **Which analysis module(s)**  -  if ambiguous (e.g. "analyse my sequences"), ask what they want to test (diversity? LD? divergence? all?).
3. **Divergence analysis specifically**: ask whether they have two separate files (use `--input2`) or one file with a population assignment table (use `--pop-file`). If neither is available, explain that divergence requires a second population.
4. **fuliout (Fu & Li with outgroup)**: ask which sequence in the alignment is the outgroup. The outgroup name is passed as `--outgroup <seq_name>`. It is extracted from the alignment and removed from the ingroup before analysis.
5. **hka analysis**: ask for the HKA locus file path (TSV with columns: locus, S, D, n). If the user wants to compute S and D from actual alignments, help them build the file first, then run `--analysis hka --hka-file <path>`.
6. **mk (McDonald-Kreitman) analysis**: confirm (a) which sequence in the alignment is the outgroup (`--outgroup <seq_name>`) and (b) that the alignment is an in-frame coding sequence (length divisible by 3, no internal stop codons). The alignment must include both ingroup sequences and the outgroup.
7. **kaks analysis**: confirm that the alignment is an in-frame coding sequence (length divisible by 3). No outgroup required. Warn the user if omega = Ka/Ks is undefined (Ks = 0 or Ka/Ks numerically saturated).
8. **fufs analysis**: no extra inputs needed  -  Fu's Fs reuses π (nucleotide diversity) and H (haplotype count) already computed by the polymorphism module, which always runs. Confirm the user understands the conventional significance threshold is Fs < 0 with S_k ≤ 0.02.
9. **sfs analysis**: folded SFS is always computed. Ask whether they have an outgroup in the alignment to produce the unfolded SFS (`--outgroup <seq_name>`). If so, the same outgroup used for fuliout/mk can be reused.
10. **tstv analysis**: no extra inputs needed. Works on any alignment (coding or non-coding). Particularly useful for assessing saturation; ask if they want it combined with divergence analysis.
11. **codon analysis**: requires an in-frame coding alignment (no 5′ UTR). Stop codons are skipped automatically but the user must ensure the alignment is in-frame from position 0. Pair with `kaks` for a comprehensive coding evolution analysis.
12. **Sliding window**  -  ask window size and step if they want sliding-window output.
13. **Output directory**  -  default to `results/` next to the input file if not specified.

Skip clarification for trivial cases: if the user has already provided all needed information, proceed immediately.

---

## Workflow

1. **Identify intent** using the decision tree above.
2. **Confirm** file path(s) and output directory.
3. **Construct CLI command** (see CLI Reference below).
4. **Run** `python skills/dnasp/dnasp.py [args]`.
5. **Parse stdout** to check for errors or warnings (e.g. n < 3 warnings).
6. **Explain results** in plain language: what each key statistic means, whether values are noteworthy, and what follow-up analyses might be informative. Reference `docs/index.md` for interpretation guidance.
7. **Suggest follow-ups** where relevant (e.g. after polymorphism → ask if they want LD or divergence).

---

## CLI Reference

```bash
# Polymorphism + neutrality tests only (default)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --output results/

# Select specific analyses
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis ld,recombination \
    --output results/

# All analyses (no divergence data)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis polymorphism,ld,recombination,popsize,indel \
    --output results/

# Sliding window (100 bp window, 25 bp step)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --window 100 --step 25 \
    --output results/

# Divergence  -  two separate FASTA files
python skills/dnasp/dnasp.py \
    --input pop1.fas \
    --input2 pop2.fas \
    --analysis divergence \
    --output results/

# Divergence  -  one alignment with population assignment file
python skills/dnasp/dnasp.py \
    --input combined.fas \
    --pop-file populations.txt \
    --analysis divergence \
    --output results/

# All analyses including divergence
python skills/dnasp/dnasp.py \
    --input pop1.fas \
    --input2 pop2.fas \
    --analysis all \
    --output results/

# Fu & Li D/F with outgroup (outgroup seq named "outgroup" is in the alignment)
python skills/dnasp/dnasp.py \
    --input aln_with_outgroup.fas \
    --outgroup outgroup \
    --analysis fuliout \
    --output results/

# HKA test (pre-computed locus file)
python skills/dnasp/dnasp.py \
    --input aln.fas \
    --hka-file hka_loci.tsv \
    --analysis hka \
    --output results/

# McDonald-Kreitman test (outgroup sequence named "outgroup" is in the alignment)
python skills/dnasp/dnasp.py \
    --input coding_aln_with_outgroup.fas \
    --outgroup outgroup \
    --analysis mk \
    --output results/

# Ka/Ks  -  Nei-Gojobori pairwise dN/dS (in-frame coding alignment, no outgroup needed)
python skills/dnasp/dnasp.py \
    --input coding_aln.fas \
    --analysis kaks \
    --output results/

# MK + polymorphism combined
python skills/dnasp/dnasp.py \
    --input coding_aln_with_outgroup.fas \
    --outgroup outgroup \
    --analysis polymorphism,mk \
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
    --outgroup outgroup \
    --analysis sfs \
    --output results/

# All neutrality tests together (Tajima D, Fu & Li D*/F*, R2, Fu's Fs, SFS)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis polymorphism,fufs,sfs \
    --output results/

# Transition/transversion ratio (any alignment)
python skills/dnasp/dnasp.py \
    --input alignment.fas \
    --analysis tstv \
    --output results/

# Codon usage bias (RSCU + ENC; in-frame coding alignment)
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

# Demo mode (built-in synthetic data)
python skills/dnasp/dnasp.py \
    --demo \
    --output /tmp/dnasp_demo
```

### Flag Reference

| Flag | Type | Default | Description |
|---|---|---|---|
| `--input` | path |  -  | Alignment file (FASTA or NEXUS) |
| `--input2` | path |  -  | Second population alignment (for `divergence`) |
| `--pop-file` | path |  -  | Population assignment TSV (alternative to `--input2`) |
| `--outgroup` | string |  -  | Sequence name to use as outgroup (for `fuliout` and `mk`) |
| `--hka-file` | path |  -  | HKA locus file (TSV: locus, S, D, n) for `hka` |
| `--analysis` | string | `polymorphism` | Comma-separated analyses or `all` |
| `--output` | path | `./dnasp_out/` | Output directory |
| `--window` | int | 0 | Sliding window size (bp); 0 = disabled |
| `--step` | int | = window | Sliding window step (bp) |
| `--demo` | flag |  -  | Run on built-in synthetic dataset |

### Population file format (`--pop-file`)

Tab-separated, one row per sequence, `#` lines are comments:

```
# Population assignment
seq1	Pop_Africa
seq2	Pop_Africa
seq3	Pop_Europe
seq4	Pop_Europe
```

### HKA locus file format (`--hka-file`)

Tab-separated, one row per locus, `#` lines are comments, header row optional:

```
# locus   S   D   n
ACE        5  10  10
G6PD       2   8  12
white     12  18  10
```

Columns:

- **locus**: any identifier string
- **S**: segregating sites in the ingroup sample (count, not rate)
- **D**: fixed differences between ingroup and outgroup/sister species (count)
- **n**: ingroup sample size (number of sequences)

To build this file from alignments: use `--analysis polymorphism` on each ingroup alignment (read S from `results.tsv`), then count fixed differences between ingroup consensus and outgroup sequence manually or with a separate tool.

---

## Valid Analysis Values

| Value | Module | What it computes |
|---|---|---|
| `polymorphism` | Polymorphism & neutrality | π, k, S, Eta, H, Hd, θ_W, Tajima's D, Fu & Li D*/F*, R2, GC |
| `ld` | Linkage Disequilibrium | D, D', R² per pair; ZnS, Za, ZZ genome-wide; LD decay scatter |
| `recombination` | Recombination | Rm (min. recombination events, four-gamete test, Hudson & Kaplan 1985) |
| `popsize` | Population Size History | Mismatch distribution, raggedness r, CV |
| `indel` | InDel Polymorphism | InDel events, InDel haplotypes, k(i), π(i), θ(i), Tajima's D(i) |
| `divergence` | Divergence | Dxy, Da, fixed differences, shared & private polymorphisms |
| `fuliout` | Fu & Li D/F with outgroup | η (total derived), η_e (external/singleton derived), D, F (Fu & Li 1993) |
| `hka` | HKA multi-locus test | MLE T̂, χ² neutrality test, per-locus θ̂, E[S], E[D] (Hudson et al. 1987) |
| `mk` | McDonald-Kreitman test | Pn, Ps, Dn, Ds counts; α (proportion adaptive substitutions); NI (neutrality index); DoS (direction of selection); Fisher's exact P |
| `kaks` | Ka/Ks (dN/dS) | Nei-Gojobori (1986) pairwise averages: S sites, N sites, Ks (synonymous rate), Ka (nonsynonymous rate), ω = Ka/Ks |
| `fufs` | Fu's Fs test | θ_π, S_k = P(K ≤ H \| θ, n) via Ewens sampling formula, Fs = ln(S_k/(1−S_k)); significant at 0.02 when Fs << 0 |
| `sfs` | Site frequency spectrum | Folded SFS (always); unfolded SFS with `--outgroup`; bar-chart figure (sfs.png) |
| `tstv` | Transition/Transversion ratio | Ts (purine↔purine or pyrimidine↔pyrimidine), Tv (purine↔pyrimidine) counts across all pairs; Ts/Tv ratio; per-site rates |
| `codon` | Codon usage bias | RSCU per codon (Sharp & Li 1987); ENC (Wright 1990) from 20 (max bias) to 61 (no bias); RSCU bar chart (codon_usage.png) |
| `faywu` | Fay & Wu's H + Zeng's E | Outgroup-polarised neutrality tests. θ_H (Fay & Wu 2000), θ_L (Zeng et al. 2006), H = θ_π − θ_H, E = θ_L − θ_W. Requires `--outgroup`. |
| `fst` | Population differentiation | Hudson et al. (1992) pairwise Fst = 1 − π_s/π_t for each pop pair; within-pop π, Dxy; mean Fst across pairs; Fst bar chart (fst.png). Requires `--pop-file`. |

---

## Demo

```bash
python skills/dnasp/dnasp.py --demo --output /tmp/dnasp_demo
```

Expected (10 ingroup + 1 outgroup × 300 bp; Pop1/Pop2; in-frame CDS):

| Statistic | Expected value |
|---|---|
| S | 5 |
| H (haplotypes) | 8 |
| Hd | 0.9556 |
| π | 0.006889 |
| Tajima's D | 0.6789 |
| Ts / Tv | 77 / 16 = 4.8125 |
| ENC | 23.00 (strong codon bias) |
| Fay & Wu H | 0.004148 |
| Zeng E | −0.001317 |
| MK Pn/Ps/Dn/Ds | 2 / 3 / 1 / 1 |
| α (MK) | 0.333 |
| Ka / Ks / ω | 0.00298 / 0.02281 / 0.131 |
| Fst (Pop1 vs Pop2) | 0.0566 |

---

## Algorithm Summary

All formulas match DnaSP 6. See `docs/index.md` for full derivations and references.

**Gap treatment** (complete deletion): exclude any column where ≥1 sequence has `-`, `?`, or `N`. All statistics use L_net (net sites after exclusion).

**Polymorphism module**:
- k = mean pairwise differences (absolute); π = k / L_net
- Hd = n/(n−1) × (1 − Σpᵢ²)
- θ_W = S/a₁; θ_W_nuc = θ_W/L_net (a₁ = Σ 1/i, i=1..n-1)
- Tajima's D = (k − S/a₁) / √(e₁S + e₂S(S−1))
- Fu & Li D* = (S/Aₙ − η_s(n−1)/n) / √(uD·S + vD·S²) (Simonsen 1995, eq. A3)
- Fu & Li F* = (k − η_s(n−1)/n) / √(uF·S + vF·S²) (Simonsen 1995, eq. A5)
- R2 = √(Σ(Uᵢ − k/2)² / n) / Sw (Ramos-Onsins & Rozas 2002)

**LD module**: For each pair of strictly biallelic sites, compute D (Lewontin & Kojima 1960), D' (Lewontin 1964), R² (Hill & Robertson 1968), and chi-square p-value via `erfc(√(χ²/2))` (no scipy needed). ZnS = mean R² over all pairs (Kelly 1997). Za = mean R² over adjacent biallelic pairs (Rozas 2001). ZZ = Za − ZnS.

**Recombination module**: Four-gamete test (Hudson & Kaplan 1985)  -  a pair of biallelic sites is incompatible when all four gamete combinations are observed. Rm = minimum number of recombination events, computed by the interval-stabbing greedy algorithm (sort incompatible intervals by right endpoint; place a recombination point at the right endpoint whenever the left endpoint exceeds the last placed point).

**Mismatch module**: Observed pairwise-difference histogram. Raggedness r (Harpending 1994, eq. 1) = Σ(f(i) − f(i−1))². CV = σ/μ of pairwise differences (Rogers & Harpending 1992). Small r → smooth distribution → population expansion signature.

**InDel module**: InDel event = maximal run of columns where the same subset of sequences carries gaps (diallelic option of DnaSP). Statistics on InDel haplotypes, k(i), π(i), θ_W(i), Tajima's D(i) computed as for nucleotide data.

**Divergence module**: Dxy = average between-population differences per site (Nei 1987, eq. 10.20). Da = Dxy − (π₁ + π₂)/2 (net divergence). Fixed differences, shared polymorphisms, and private polymorphisms classified per Hey (1991). Complete deletion applied across both populations combined.

**Fu & Li outgroup module (fuliout)**: Outgroup sequence polarises each segregating site  -  allele matching outgroup is ancestral; others are derived. η = total derived mutations (outgroup-polarised); η_e = derived mutations in exactly 1 ingroup sequence (external/singletons on terminal branches). D = (η_e − η/aₙ) / √(uD·η + vD·η(η−1)); F = (k̄ − η_e) / √(uF·η + vF·η(η−1)); k̄ = mean pairwise differences computed over all clean ingroup sites. Variance coefficients follow Simonsen et al. (1995) Appendix B structure. Complete deletion applied to both ingroup sequences and outgroup simultaneously.

**HKA test module**: Compares the ratio of polymorphism (S_i) to divergence (D_i) across k loci. Under neutrality all loci should share the same ratio. Model: E[S_i] = θ_i f_i, E[D_i] = θ_i(1+2T) where f_i = Σ 1/j (j=1..n_i−1) and T is the scaled divergence time. MLE of T found by bisection on Σ D_i/(1+2T) = Σ(S_i+D_i)/(f_i+1+2T). χ² = Σ[(S_i−E_S_i)²/E_S_i + (D_i−E_D_i)²/E_D_i] with df = k−1 (one parameter T estimated). P-value uses the regularised upper incomplete gamma function Q(df/2, χ²/2)  -  no scipy needed.

**McDonald-Kreitman test module (mk)**: For each codon (in-frame, complete deletion at codon level  -  any non-ATCG in any sequence skips that codon; stop codons skipped): determine ingroup variation and outgroup-vs-ingroup fixed differences. A site is **polymorphic** in the ingroup if ≥2 sequences differ at any codon position. A site is **fixed** if all ingroup sequences agree but the outgroup differs. Classify each codon-site pair as synonymous or nonsynonymous using the genetic code. Accumulate Pn (nonsynonymous polymorphisms), Ps (synonymous polymorphisms), Dn (nonsynonymous fixed differences), Ds (synonymous fixed differences). Derived statistics: α = 1 − (Ds·Pn)/(Dn·Ps); NI = (Pn/Ps)/(Dn/Ds); DoS = Dn/(Dn+Ds) − Pn/(Pn+Ps). Fisher's exact P computed via hypergeometric distribution using `math.lgamma` (no scipy needed); two-tailed (sum of all table probabilities ≤ observed probability).

**Fu's Fs module (fufs)**: Estimates θ_π = k (mean pairwise differences, always available from the polymorphism module). Uses the Ewens sampling formula  -  the probability distribution of the number of distinct alleles K_n in a sample of n sequences under the infinite-alleles model with mutation rate θ. P(K_n = k) = |s(n, k)| × θ^k / θ^(n) where |s(n, k)| are unsigned Stirling numbers of the first kind (computed by DP with Python arbitrary-precision integers; no overflow) and θ^(n) = θ(θ+1)…(θ+n−1) is the Pochhammer rising factorial. S_k = P(K_n ≤ H_obs | θ_π, n) is the probability of observing H_obs or fewer haplotypes. Fs = ln(S_k / (1−S_k)). Significant at the conventional 0.02 level when Fs << 0 (S_k ≤ 0.02). No simulation or scipy required.

**SFS module (sfs)**: For each alignment column (after complete deletion of the ingroup), counts how many sequences carry each allele. Folded SFS: records sites by minor allele count i (1 ≤ i ≤ n//2)  -  the rarer allele. Unfolded SFS (requires `--outgroup`): for each clean column where the outgroup allele is present in the ingroup, counts the number of ingroup sequences carrying the derived allele (i = 1 to n−1). Gap/ambiguous bases in any ingroup sequence → column excluded; gap in outgroup → excluded from unfolded only (folded still counts clean ingroup columns). Produces folded and (optionally) unfolded bar-chart figures.

**Ka/Ks module (kaks)**: For each pair of ingroup sequences, count synonymous sites (S_ij = (S_i+S_j)/2) and nonsynonymous sites (N_ij = 3L_codon − S_ij) using the Nei-Gojobori (1986) method  -  per codon, each of the 3 positions contributes a fraction equal to the number of synonymous alternatives out of 3; summed across all clean codons. Count synonymous (sd) and nonsynonymous (nd) differences by pathway averaging over all k! orderings when codons differ at k positions; paths through stop codons excluded. Apply Jukes-Cantor correction: Ks = −3/4 · ln(1 − 4pS/3), Ka = −3/4 · ln(1 − 4pN/3). If pS ≥ 0.75 or pN ≥ 0.75, that pair is excluded from averages (saturated). Report mean Ks, Ka, and ω = Ka/Ks across all valid pairs. ω = None when Ks = 0 (no synonymous divergence).

**Ts/Tv module (tstv)**: Classifies each pairwise nucleotide difference at every clean column (complete deletion across the full ingroup). A **transition** (Ts) is a change between two purines (A↔G) or two pyrimidines (C↔T)  -  same chemical class. A **transversion** (Tv) is a purine↔pyrimidine change (A↔C, A↔T, G↔C, G↔T). n_transitions and n_transversions accumulate across all n(n−1)/2 pairs and all clean sites. Ts/Tv = n_transitions/n_transversions; None if n_transversions = 0. Mean per-pair per-site rates: ts_per_site = n_transitions / (n_pairs × L_net), tv_per_site analogous.

**Fay & Wu / Zeng module (faywu)**: Requires `--outgroup`. Applies complete deletion including the outgroup. For each segregating site, the outgroup allele identifies the ancestral state; a site is **polarisable** when the ancestral allele appears in the ingroup. For each polarisable site with derived allele count i (1 ≤ i ≤ n−1), adds to ξ_i. Computes four per-site θ estimates from the unfolded SFS: θ_π = Σ ξ_i × 2i(n−i) / [n(n−1)] / L; θ_W = Σ ξ_i / a₁ / L (a₁ = Σ 1/k for k=1..n−1); θ_H = Σ ξ_i × 2i² / [n(n−1)] / L; θ_L = Σ ξ_i × i / (n−1) / L. H = θ_π − θ_H (Fay & Wu 2000); E = θ_L − θ_W (Zeng et al. 2006). H < 0 indicates an excess of high-frequency derived alleles (consistent with recent selective sweep). E < 0 indicates excess of low-frequency derived alleles relative to Watterson expectation.

**Fst module (fst)**: Requires `--pop-file`. Applies complete deletion across all sequences from all populations combined. For each pair of populations A and B, computes: π_A = mean within-pop pairwise differences per site; π_B analogous; π_AB = Dxy (mean between-pop pairwise differences per site). Hudson et al. (1992) estimator: Fst = 1 − π_s / π_AB where π_s = (π_A + π_B) / 2. Fst is clamped to [0, 1] (negative values from small samples are set to 0). Fst = None when π_AB = 0 (no between-population variation at any site). Mean Fst is the unweighted average across all pairs. For three or more populations, all pairwise combinations are computed.

**Codon usage module (codon)**: Reads the in-frame coding alignment in non-overlapping triplets. Triplets with any non-ATCG character or translating to a stop codon are skipped. Codon counts are pooled across all sequences. RSCU (Sharp & Li 1987): RSCU_ij = X_ij / (X_i / n_i) where X_ij = count of codon j for amino acid i, X_i = total count for amino acid i, n_i = synonymous family size. RSCU = 1.0 → uniform usage; > 1.0 → preferred; < 1.0 → avoided. ENC (Wright 1990): computed from mean corrected homozygosity per degeneracy class. For amino acids with k-fold degeneracy (k codons), corrected homozygosity F_k = (n_aa × Σpⱼ² − 1) / (n_aa − 1) where pⱼ = fraction of amino acid i encoded by codon j, n_aa = total codon count for amino acid i. Class means over all amino acids in that class: 2-fold (9 aa), 3-fold (Ile only, 1 aa), 4-fold (5 aa), 6-fold (3 aa). ENC = 2 + 9/F̄₂ + 1/F̄₃ + 5/F̄₄ + 3/F̄₆. Clamped to [20, 61].

**Key thresholds**:
- Tajima's D, Fu & Li D*/F*: require n ≥ 3 and S > 0; return `n.a.` otherwise.
- LD statistics: require ≥ 2 strictly biallelic sites.
- Divergence: require ≥ 1 sequence per population and L_net > 0.
- fuliout: requires n ≥ 4 ingroup sequences and η > 0 (at least one outgroup-polarised derived mutation).
- hka: requires ≥ 2 valid loci (each with n ≥ 2). Returns HKAStats with chi2=0 if only 1 locus provided.
- mk: requires n ≥ 2 ingroup sequences, alignment length divisible by 3 (in-frame coding), and `--outgroup <seq_name>`. Returns None if outgroup not provided or alignment not in-frame. α, NI, DoS are None when any denominator is zero.
- kaks: requires n ≥ 2 sequences and alignment length divisible by 3 (in-frame coding). ω = None when Ks = 0. Pairs where pS or pN ≥ 0.75 (JC saturation) are excluded.
- fufs: requires n ≥ 2 and H ≥ 1. If k = 0 (all sequences identical), Fs is defined but θ_π = 0 → degenerate. Uses polymorphism stats already computed; no additional inputs needed.
- sfs: requires n ≥ 2. Folded SFS is always produced from the ingroup. Unfolded SFS requires `--outgroup` and at least one site where the outgroup allele appears in the ingroup.
- tstv: requires n ≥ 2 and at least one clean (non-gap, unambiguous ATCG) column. Ts/Tv = None when n_transversions = 0 (all differences are transitions). Works on both coding and non-coding alignments.
- codon: requires alignment length divisible by 3 (in-frame from position 0). ENC = None when any of the four degeneracy classes (2-fold, 3-fold, 4-fold, 6-fold) lacks sufficient amino acid observations (n_aa < 2 for all members of a class). For short alignments this is common; longer coding sequences (> 300 bp) are recommended for reliable ENC estimates.
- faywu: requires `--outgroup` and at least one polarisable segregating site (site where the ancestral allele appears in the ingroup and a derived allele exists at 1 ≤ count ≤ n−1). Returns None for H and E when n_polarised = 0. Sites where the outgroup has a gap or non-ATCG character, or the ancestral allele is absent from the ingroup, are excluded.
- fst: requires `--pop-file` with at least 2 populations. Fst = None for a pair when Dxy = 0 (no between-pop variation). fst_mean = None when all pairs have Dxy = 0. With only 1 population in the pop file, returns empty FstStats with a warning.

---

## Interpretation Guide (Quick Reference)

| Result | Interpretation | Caution |
|---|---|---|
| Tajima's D < 0 | Excess low-frequency variants → population expansion or purifying selection | Need to rule out demographic history |
| Tajima's D > 0 | Excess intermediate-frequency variants → balancing selection or population bottleneck | Same |
| D' = 1 (or −1) | No evidence of recombination between that pair of sites | Valid only when n is large enough |
| R² high, distance short | Recent LD → low recombination rate in that region | |
| ZZ > 0 | Adjacent pairs have higher LD than non-adjacent → recombination breaking up distant LD | |
| Rm ≥ 1 | At least Rm recombination events required to explain data | Rm is a minimum; true Rm could be higher |
| Raggedness r small | Smooth mismatch distribution → consistent with population expansion | |
| Da < 0 | Net divergence negative → within-population diversity exceeds between; can occur by chance | Da should be ≈ 0 under neutrality |
| n_fixed >> n_shared | Populations are highly differentiated; long divergence time | |
| Fu & Li D > 0 (outgroup) | Excess external (singleton) mutations → possibly purifying selection removing most lineages | Compare with no-outgroup D* |
| Fu & Li D < 0 (outgroup) | Fewer singletons than expected → selective sweep or population expansion | |
| HKA P < 0.05 | Ratio of polymorphism to divergence differs across loci → departure from neutral model | One locus may be under selection |
| HKA P > 0.05 | Polymorphism/divergence ratio consistent across loci → consistent with neutral model | |
| T̂ (HKA) large | Long divergence time relative to N_e | Calibrate with known mutation rate if possible |
| MK Fisher P < 0.05 | Ratio of Pn/Ps differs from Dn/Ds → departure from neutral model | Could indicate positive selection (α > 0) or relaxed constraint |
| α > 0 (MK) | Positive proportion of nonsynonymous fixations are adaptive | α is the fraction of substitutions driven to fixation by positive selection |
| α < 0 (MK) | More polymorphism than divergence at nonsynonymous sites relative to synonymous → slight deleterious mutations segregating | Common; use DoS as alternative measure |
| NI > 1 (MK) | Excess nonsynonymous polymorphism relative to divergence → slightly deleterious variants segregating | Same direction as α < 0 |
| NI < 1 (MK) | Deficit of nonsynonymous polymorphism → positive selection driving rapid fixation | |
| DoS > 0 (MK) | Divergence more nonsynonymous than polymorphism → positive selection signature | Scale-free; compare across genes |
| DoS < 0 (MK) | Divergence more synonymous → purifying selection removing nonsynonymous variants before fixation | |
| ω (Ka/Ks) < 1 | Purifying (negative) selection  -  nonsynonymous changes removed faster than synonymous | Expected for most functional genes |
| ω ≈ 1 | Neutral evolution  -  synonymous and nonsynonymous rates similar | |
| ω > 1 | Positive selection  -  nonsynonymous changes accumulate faster than synonymous | Rare; strong evidence of adaptive evolution |
| ω = None | Ks = 0 (no synonymous divergence between sequences) or all pairs JC-saturated | Use with very short or very similar sequences |
| Fs << 0 (Fu's Fs) | Far fewer haplotypes than expected given π → population expansion or positive selection | Significant at 0.02 level when S_k ≤ 0.02 |
| Fs ≈ 0 (Fu's Fs) | Haplotype count consistent with neutral expectation | |
| Fs > 0 (Fu's Fs) | More haplotypes than expected → balancing selection or population subdivision | Rarely significant |
| SFS singleton-heavy (i=1 dominant) | Excess rare variants → expansion, purifying selection, or recent bottleneck recovery | Consistent with negative Tajima's D |
| SFS flat or U-shaped | Uniform or high-frequency-skewed variants → balancing selection | Consistent with positive Tajima's D |
| Unfolded SFS high at n−1 | Many near-fixed derived alleles → directional selection or recent sweep ancestry | |
| Ts/Tv ≈ 2 | Typical transitional bias for nuclear DNA  -  transitions more mutable than transversions | Expected baseline; varies by locus and taxon |
| Ts/Tv > 10 | Strong transition bias → common in mitochondrial DNA or highly constrained sequences | |
| Ts/Tv < 0.5 | Transversion excess → substitution saturation at transitions, or non-neutral patterns | Check alignment quality; consider JC correction |
| Ts/Tv = None | No transversions observed (all differences are transitions) | Normal for highly similar sequences |
| ENC ≈ 61 | No codon usage bias  -  all synonymous codons used equally | Expected under neutral drift |
| ENC 35-60 | Moderate codon usage bias | Moderate translational selection or mutational bias |
| ENC < 35 | Strong codon usage bias  -  strong preference for particular synonymous codons | Likely translational selection; compare RSCU to tRNA availability |
| ENC ≈ 20 | Maximum bias  -  only one codon per amino acid used | Extreme selection or very small effective population |
| ENC = None | Insufficient codon data for one or more degeneracy classes | Use longer alignment (> 300 bp recommended) |
| RSCU > 1 for a codon | Preferred codon within its synonymous family | Cross-reference with tRNA gene copy numbers |
| RSCU = 0 for a codon | Completely avoided codon | May reflect strong translational selection |
| H < 0 (Fay & Wu) | Excess high-frequency derived alleles → consistent with a selective sweep or hitchhiking | Requires accurate outgroup to polarise mutations |
| H ≈ 0 (Fay & Wu) | No excess of high-frequency derived alleles → consistent with neutrality | |
| H > 0 (Fay & Wu) | Excess intermediate-frequency derived alleles → complement to Tajima's D > 0 | |
| E < 0 (Zeng) | θ_L < θ_W → excess low-frequency derived alleles relative to Watterson → purifying selection or bottleneck | Use together with H for more power |
| E > 0 (Zeng) | θ_L > θ_W → more intermediate/high-frequency derived variants than expected | Unusual; check for sampling issues |
| Fst < 0.05 | Little genetic differentiation between populations (Wright 1978) | |
| Fst 0.05-0.15 | Moderate genetic differentiation | |
| Fst 0.15-0.25 | Great genetic differentiation | |
| Fst > 0.25 | Very great genetic differentiation | |
| Fst = 1.0 | Complete fixation for different alleles  -  no shared polymorphism between populations | |
| Fst = None | Dxy = 0 (no between-population variation at any clean site) | Both pops may be monomorphic for the same allele |

---

## Example Queries

- "Compute nucleotide diversity for my rp49 alignment at `~/data/rp49.fas`"
- "Run all neutrality tests on alignment.fas"
- "Calculate linkage disequilibrium for my SNP data"
- "What is the minimum number of recombination events in my dataset?"
- "Is there a mismatch distribution signature of population expansion?"
- "Measure divergence between my African and European populations"
- "Run everything  -  I have pop1.fas and pop2.fas"
- "Run DnaSP on my Drosophila sequences with a 100 bp sliding window"

---

## Output Structure

```
output_directory/
├── report.md                  # Full Markdown report (all active modules)
├── results.tsv                # DnaSP-compatible tab-delimited table
├── ld_pairs.tsv               # Pairwise LD table (if --analysis ld)
├── figures/
│   ├── summary.png            # Bar chart: π, θ_W, Hd
│   ├── sliding_window.png     # π and Tajima D per window (if --window used)
│   ├── ld_decay.png           # R² vs distance scatter plot (if ld)
│   ├── mismatch.png           # Mismatch distribution bar chart (if popsize)
│   ├── sfs.png                # Site frequency spectrum bar chart (if sfs)
│   ├── codon_usage.png        # RSCU bar chart coloured by amino acid family (if codon)
│   └── fst.png                # Pairwise Fst bar chart with differentiation thresholds (if fst)
└── reproducibility/
    ├── commands.sh            # Exact command used
    ├── environment.yml        # Conda environment spec
    └── checksums.sha256       # SHA-256 of input and output files
```

---

## Gotchas

- **FASTA wrapping**: DnaSP exports `>'name'  [comment]` headers and wrapped sequences. Always use `dnasp.py` to parse DnaSP-generated FASTA  -  generic parsers may fail on DnaSP headers.
- **Complete deletion**: Excludes any site with a gap in *any* sequence. High-gap alignments can dramatically reduce L_net. Check `L_total` vs `L_net` in the report.
- **n < 3**: Tajima's D and Fu & Li require n ≥ 3; report `n.a.` otherwise.
- **Divergence gap mask**: Applied across both populations combined  -  a gap in *either* population removes the column.
- **LD with few sequences**: D' tends to be inflated (→ 1.0) when n is small. Report D' alongside sample size.
- **Rm is a minimum bound**: The true number of recombination events is ≥ Rm. Rm = 0 does not mean no recombination occurred.
- **F\* discrepancy**: For RAD-seq multi-MSA context, DnaSP v6 uses Achaz (2009) variance. For standard single-alignment FASTA, this skill uses Simonsen (1995) A5-A6, which is correct. D\* agrees regardless.
- **NEXUS MATCHCHAR**: If `.` is the MATCHCHAR, the parser expands relative to the first sequence. If the first sequence is the outgroup, re-order before analysis.

---

## Agent Boundary

The agent (LLM) dispatches the script, explains results, and recommends follow-up analyses. The skill (Python) executes all numerical computation. The agent must **not** recompute or override numerical outputs  -  trust `dnasp.py` results. If a statistic seems unexpected, re-run and check the raw `results.tsv`, then explain the value rather than modifying it.

---

## Safety

- **Local-first**: No sequence data is uploaded. All processing is on-device.
- **No hallucinated statistics**: All formulas trace to cited papers and the original DnaSP VB source code.
- **Disclaimer**: *ClawBio is a research and educational tool. Not a medical device.*

---

## Citations

- [Rozas et al. (2017) J. Hered. 108:591-593](https://doi.org/10.1093/jhered/esx062)  -  DnaSP v6
- [Tajima (1989) Genetics 123:585-595](https://www.genetics.org/content/123/3/585)  -  Tajima's D
- [Fu & Li (1993) Genetics 133:693-709](https://www.genetics.org/content/133/3/693)  -  D\*, F\*
- [Simonsen et al. (1995) Genetics 141:413-429](https://www.genetics.org/content/141/1/413)  -  variance coefficients A3-A6
- [Nei & Tajima (1981) Genetics 97:145-163](https://www.genetics.org/content/97/1/145)  -  haplotype diversity
- [Ramos-Onsins & Rozas (2002) Mol. Biol. Evol. 19:2092-2100](https://doi.org/10.1093/oxfordjournals.molbev.a004068)  -  R2
- [Lewontin & Kojima (1960) Evolution 14:458-472](https://doi.org/10.2307/2405649)  -  D
- [Lewontin (1964) Genetics 49:49-67](https://www.genetics.org/content/49/1/49)  -  D'
- [Hill & Robertson (1968) Theor. Appl. Genet. 38:226-231](https://doi.org/10.1007/BF01245622)  -  R²
- [Kelly (1997) Genetics 146:1197-1206](https://www.genetics.org/content/146/3/1197)  -  ZnS
- [Rozas et al. (2001) Genetics 158:1321-1330](https://www.genetics.org/content/158/3/1321)  -  Za, ZZ
- [Hudson & Kaplan (1985) Genetics 111:147-164](https://www.genetics.org/content/111/1/147)  -  Rm
- [Harpending (1994) Hum. Biol. 66:591-600](https://www.jstor.org/stable/41465138)  -  raggedness r
- [Rogers & Harpending (1992) Mol. Biol. Evol. 9:552-569](https://doi.org/10.1093/oxfordjournals.molbev.a040727)  -  mismatch CV
- [Nei (1987) Molecular Evolutionary Genetics. Columbia Univ. Press.](https://cup.columbia.edu)  -  Dxy (eq. 10.20)
- [Hey (1991) Genetics 128:831-840](https://www.genetics.org/content/128/4/831)  -  fixed differences classification
- [McDonald & Kreitman (1991) Nature 351:652-654](https://doi.org/10.1038/351652a0)  -  MK test
- [Nei & Gojobori (1986) Mol. Biol. Evol. 3:418-426](https://doi.org/10.1093/oxfordjournals.molbev.a040410)  -  Ka/Ks synonymous sites method
- [Fu (1997) Genetics 147:915-925](https://www.genetics.org/content/147/2/915)  -  Fu's Fs test
- [Ewens (1972) Theor. Popul. Biol. 3:87-112](https://doi.org/10.1016/0040-5809(72)90035-4)  -  Ewens sampling formula (basis for Fu's Fs)
- [Sharp & Li (1987) Nucleic Acids Res. 15:1281-1295](https://doi.org/10.1093/nar/15.3.1281)  -  RSCU (Relative Synonymous Codon Usage)
- [Wright (1990) Gene 87:23-29](https://doi.org/10.1016/0378-1119(90)90491-9)  -  ENC (Effective Number of Codons)
- [Fay & Wu (2000) Genetics 155:1405-1413](https://doi.org/10.1093/genetics/155.3.1405)  -  H statistic (θ_H, outgroup-polarised)
- [Zeng et al. (2006) Genetics 174:1431-1439](https://doi.org/10.1534/genetics.106.061432)  -  E statistic (θ_L, complement to H)
- [Hudson et al. (1992) Genetics 132:583-589](https://www.genetics.org/content/132/2/583)  -  Fst estimator (1 − π_s/π_t)
