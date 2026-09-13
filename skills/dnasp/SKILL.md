---
name: dnasp
description: >-
  Full reimplementation of DnaSP 6 for population genetics analysis of aligned
  DNA sequences. Covers nucleotide diversity, haplotype statistics, neutrality
  tests (Tajima's D, Fu & Li's D*/F*, R2), linkage disequilibrium (D, D', R²,
  ZnS, Za, ZZ), minimum recombination (Rm), mismatch distribution, InDel
  polymorphism, between-population divergence (Dxy, Da, fixed/shared sites),
  outgroup-based Fu & Li D/F tests (fuliout), the HKA two-locus neutrality
  test (hka), the McDonald-Kreitman test (mk), Ka/Ks (dN/dS) via the
  Nei-Gojobori (1986) method (kaks), Fu's Fs test (fufs), the site frequency
  spectrum (sfs, folded and outgroup-unfolded), transition/transversion ratio
  (tstv), and codon usage bias  -  RSCU (Sharp & Li 1987) and ENC (Wright 1990)
  (codon). Accepts pre-aligned FASTA/NEXUS or a multi-sample VCF (one MSA per
  CHROM); outputs DnaSP-compatible TSV and a Markdown report.
license: MIT
metadata:
  version: "0.5.0"
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
      format:
        - fasta
        - fas
        - nexus
        - nex
      description: >-
        Aligned DNA sequences (pre-aligned, equal-length). FASTA (including
        DnaSP-style >'name' [comment] headers) or NEXUS (MATCHCHAR, INTERLEAVE).
      required: true
    - name: vcf
      type: file
      format:
        - vcf
      description: >-
        Multi-sample VCF (--vcf). Converted to one aligned MSA per CHROM
        (biallelic SNPs only; phased -> haplotype rows). Alternative to
        --input. Optional --region CHROM, --vcf-merge to pool all CHROMs.
      required: false
    - name: alignment2
      type: file
      format:
        - fasta
        - fas
        - nexus
        - nex
      description: >-
        Second-population alignment for divergence analysis (--input2).
        Alternative to --pop-file. Sequences must have same length as --input.
      required: false
    - name: pop_file
      type: file
      format:
        - tsv
        - txt
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
      format:
        - tsv
        - txt
      description: >-
        HKA locus file: whitespace-separated, exactly two loci, columns
        locus n S L_poly D [L_div] [chrom]. Required for --analysis hka.
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
    - name: genetic_code
      type: string
      description: >-
        Codon table for mk/kaks/codon: "standard" or "vertebrate-mitochondrial"
        (TGA=Trp, AGA/AGG=stop, ATA=Met; for COII/cytb/ND-type loci). Default: standard.
      required: false
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Markdown analysis report with statistics and interpretation
    - name: results_table
      type: file
      format:
        - tsv
      description: DnaSP-compatible tab-delimited results
    - name: ld_pairs
      type: file
      format:
        - tsv
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
      bins:
        - python3
    always: false
    emoji: ""
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
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
      - two-locus neutrality
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
      - VCF population genetics
      - multi-sample VCF nucleotide diversity
      - VCF to haplotypes
      - population genomics from VCF
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
- Divergence between populations, Dxy, Da, net divergence, fixed / shared / exclusive sites
- Fu & Li D/F with outgroup, outgroup-based neutrality test, polarised mutations
- HKA test, Hudson-Kreitman-Aguadé, two-locus neutrality, polymorphism/divergence ratio
- McDonald-Kreitman test, MK test, adaptive evolution, neutrality index, direction of selection, α (alpha)
- Ka/Ks, dN/dS, omega, synonymous substitution rate, nonsynonymous substitution rate, Nei-Gojobori, coding sequence divergence
- Fu's Fs, Fu 1997 neutrality test, haplotype frequency neutrality
- Population genetics from a VCF, multi-sample VCF, VCF to haplotypes, per-CHROM diversity from variant calls
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
| "HKA test", "Hudson-Kreitman-Aguadé", "two-locus neutrality", "polymorphism/divergence ratio" | `hka` | `--hka-file <file>` |
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
5. **hka analysis**: ask for the HKA locus file path (whitespace-separated, exactly two loci, columns `locus n S L_poly D [L_div] [chrom]`). If the user needs to compute S and D from alignments, help them build the file first, then run `--analysis hka --hka-file <path>`.
6. **mk (McDonald-Kreitman) analysis**: confirm (a) which sequence in the alignment is the outgroup (`--outgroup <seq_name>`) and (b) that the alignment is an in-frame coding sequence (length divisible by 3, no internal stop codons under the genetic code in use). The alignment must include both ingroup sequences and the outgroup. If the locus is mitochondrial (e.g. COII/COX2, cytb, ND genes), ask and pass `--genetic-code vertebrate-mitochondrial` -- under the standard code TGA is misread as a stop, and a stop inside the coding region is analysed as a 21st amino acid (as DnaSP does), so every TGA/TGG change would be scored as a replacement; the report counts such codons and the CLI warns.
7. **kaks analysis**: confirm that the alignment is an in-frame coding sequence (length divisible by 3). No outgroup required for a plain ingroup-vs-ingroup diversity comparison; if `--outgroup <seq_name>` is set (e.g. already in use for mk/fuliout in the same run), kaks switches to ingroup-vs-outgroup divergence pairs instead, matching DnaSP's own behaviour once an outgroup is defined. Warn the user if omega = Ka/Ks is undefined (Ks = 0 or Ka/Ks numerically saturated).
8. **fufs analysis**: no extra inputs needed  -  Fu's Fs reuses π (nucleotide diversity) and H (haplotype count) already computed by the polymorphism module. It reports Fs and S′ but no significance level; a formal test needs coalescent simulation (planned for a future release).
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

# Ka/Ks  -  Nei-Gojobori pairwise dN/dS (in-frame coding alignment; --outgroup optional, switches to ingroup-vs-outgroup divergence)
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
| `--vcf` | path |  -  | Multi-sample VCF; analyses run once per CHROM (one MSA each), as in DnaSP |
| `--region` | string |  -  | With `--vcf`: restrict to this CHROM only |
| `--vcf-merge` | flag |  -  | With `--vcf`: pool every CHROM into one MSA (genome-wide summary only  -  not valid for per-site statistics) |
| `--input2` | path |  -  | Second population alignment (for `divergence`) |
| `--pop-file` | path |  -  | Population assignment TSV (alternative to `--input2`) |
| `--outgroup` | string |  -  | Sequence name to use as outgroup (for `fuliout` and `mk`) |
| `--hka-file` | path |  -  | HKA locus file (2 loci: `locus n S L_poly D [L_div] [chrom]`) for `hka` |
| `--genetic-code` | `standard` \| `vertebrate-mitochondrial` | `standard` | Codon table for `mk`/`kaks`/`codon`; use `vertebrate-mitochondrial` for mtDNA loci (e.g. COII), where TGA is Trp, AGA/AGG are stop, ATA is Met |
| `--analysis` | string | `polymorphism` | Comma-separated analyses or `all` |
| `--output` | path | `./dnasp_out/` | Output directory |
| `--window` | int | 0 | Sliding window size (bp); 0 = disabled |
| `--step` | int | = window | Sliding window step (bp) |
| `--demo` | flag |  -  | Run on built-in synthetic dataset |

### VCF input (`--vcf`)

Multi-sample VCF converted to aligned haplotype sequences, following DnaSP 6
(`multifilefrmvcf.vb::readvcf`):

- **One MSA per `CHROM`**   -   analyses are run once per CHROM, exactly as DnaSP's
  VCF/RAD mode does. `--region CHROM` restricts to one; `--vcf-merge` pools all
  CHROMs into a single MSA (a deliberate genome-wide summary only   -   it mixes
  unlinked regions and is not valid for π, Tajima's D or the SFS).
- **Biallelic SNPs only**: indels, multi-base REF/ALT, and multiallelic sites
  are skipped (reported in the run summary).
- **FORMAT** must start with `GT`. **FILTER is ignored** (as in DnaSP).
- **Diploid**: phased `|` → two haplotype rows per sample (`<sample>_h1`,
  `<sample>_h2`); unphased `/` → two rows only when homozygous, otherwise both
  become gaps (DnaSP cannot phase them); `.` → gaps. **On unphased data this is
  lossy**: any SNP with one or more heterozygous calls has a gap in the column,
  and complete deletion then removes the whole column, so on a heterozygous
  unphased VCF most sites drop. The run summary reports the count
  (`n_unphased_het_sites`); phase the VCF for a full site set.
- **Haploid** (`GT` = `0`/`1`): one row per sample. Ploidy is validated per
  genotype  -  a VCF that mixes haploid and diploid calls, or contains a
  polyploid call, is rejected with an error rather than silently decoded.
- **Per-CHROM sample set**: a sample whose GT starts with `.` at the first
  retained variant of a CHROM is dropped from that whole MSA.
- Population split: pass `--pop-file` keyed by base sample ID (DnaSP's
  `.SG.txt` files work directly   -   `sample<space>population`).
- **Per-site diversity is per variant site, not per base.** A VCF carries only
  variant rows, so the MSA has one column per retained SNP and `NetSites` is the
  SNP count. π and θ_W per site are therefore computed over variant sites, as in
  DnaSP 6 (`multifilefrmvcf.vb`), which prints the same caveat. The report and
  TSV carry a note; rescale with (variant sites / callable sites) if you need
  per-base values. Counts (S, η, H) and scale-free statistics (Hd, Tajima's D,
  Fu & Li D\*/F\*, R2) are unaffected.
- **`--window` on a VCF slides over SNP index, not base pairs** — VCF POS is not
  used to place columns, so a window is a range of consecutive retained SNPs.
  The report labels this.

This is standard-mode analysis on a VCF-derived alignment; it is **not** a full
port of DnaSP's RAD engine (no Achaz F\* variances, no per-MSA Mean row, no
`.Hetz`/`.Btw`/`.GFlow` outputs, and DnaSP's downstream MNP/multiallelic site
handling is not reproduced   -   a CHROM with those can differ from DnaSP by a site
or two).

### What "DnaSP parity" means here

The per-module methodology notes above map each estimator to the DnaSP 6 Visual
Basic routine it was ported from, and the unit tests assert values derived from
the primary literature and from those routines (analytically verifiable inputs:
monomorphic alignments, single-segregating-site alignments, hand-computable
cases). The unit suite is **not** a bulk regression against DnaSP 6 GUI output,
and its fixtures are synthetic.

The comparison against real DnaSP output is kept with the associated manuscript,
not in this repo: a supplementary table compares every shared statistic on the
classic example alignments (rp49, COII, DmelOSRegion) against DnaSP 6.12 for
Windows, and against DnaSP's own shipped results for the RAD example
(`rp49_5regions`) and the `Data_Example_*.vcf` files, for which DnaSP 6.0.60
distributes its computed output. The `Adh` / `5flank` rows in the `--hka-file`
example above are placeholders showing the column layout, not a reproduction of a
published HKA result.

### Population file format (`--pop-file`)

Tab- or space-separated, one row per sequence, `#` lines are comments:

```
# Population assignment
seq1	Pop_Africa
seq2	Pop_Africa
seq3	Pop_Europe
seq4	Pop_Europe
```

### HKA locus file format (`--hka-file`)

Whitespace-separated, **exactly two loci**, `#` lines are comments, header row optional:

```
# locus   n    S    L_poly   D    L_div   chrom
Adh       81   9    4052     210  4052    A
5flank    81   8    3200     78   3200    A
```

Columns:

- **locus**: identifier
- **n**: ingroup sample size (number of sequences)
- **S**: segregating sites within the ingroup
- **L_poly**: sites analysed within the ingroup
- **D**: differences to the sister species (divergence)
- **L_div**: sites analysed for divergence (optional; defaults to L_poly)
- **chrom**: optional  -  `A` autosomal (default), `X`/`Z` (factor 0.75), `Y`/`W` (0.25)

To build this file: run `--analysis polymorphism` on each ingroup alignment (S, sites from `results.tsv`), and count divergent sites to the sister species separately.

---

## Valid Analysis Values

| Value | Module | What it computes |
|---|---|---|
| `polymorphism` | Polymorphism & neutrality | π, k, S, Eta, H, Hd, θ_W, Tajima's D, Fu & Li D*/F*, R2, GC |
| `ld` | Linkage Disequilibrium | D, D', R² per pair; ZnS, Za, ZZ genome-wide; LD decay scatter |
| `recombination` | Recombination | Rm (min. recombination events, four-gamete test, Hudson & Kaplan 1985) |
| `popsize` | Population Size History | Mismatch distribution, raggedness r, CV |
| `indel` | InDel Polymorphism | InDel events, InDel haplotypes, k(i), π(i), θ(i), Tajima's D(i) |
| `divergence` | Divergence | Dxy, Da, fixed-difference sites (Sf), shared (Ss) & exclusive (Sx) mutations |
| `fuliout` | Fu & Li D/F with outgroup | S (orientable segregating sites), η (total derived, reported), η_e (external/singleton derived sites), D, F (Fu & Li 1993) |
| `hka` | HKA two-locus test | closed-form θ̂₁, θ̂₂, T̂; χ² (df 1) with HKA (1987) variances; error if no positive-θ solution |
| `mk` | McDonald-Kreitman test | Pn, Ps, Dn, Ds counts; α (proportion adaptive substitutions); NI (neutrality index); DoS (direction of selection); Fisher's exact P |
| `kaks` | Ka/Ks (dN/dS) | Nei-Gojobori (1986) counts as DnaSP 6 applies them: mean S and N sites (stop-excluded denominators), mean Sd and Nd differences, Ks and Ka as one JC-corrected ratio, ω = Ka/Ks |
| `fufs` | Fu's Fs test | θ_π; S′ = P(K ≥ H \| θ, n) via Ewens sampling formula (log space); Fs = ln(S′/(1−S′)); no significance level (needs coalescent simulation) |
| `sfs` | Site frequency spectrum | Folded SFS (always); unfolded SFS with `--outgroup`; bar-chart figure (sfs.png) |
| `tstv` | Transition/Transversion ratio | one change per biallelic segregating site; Ts (A↔G, C↔T), Tv (purine↔pyrimidine); Ts/Tv ratio; optional outgroup polarisation |
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
- R2 = √(Σ(Uᵢ − k/2)² / n) / S (Ramos-Onsins & Rozas 2002; divides by S, the segregating-site count, matching `Dnasp_51.vb::JulioSebas_R2_CalculoAdaptado`'s `stot`, not η)

**LD module**: For each pair of strictly biallelic sites, compute D (Lewontin & Kojima 1960), D' (Lewontin 1964), R² (Hill & Robertson 1968), and chi-square p-value via `erfc(√(χ²/2))` (no scipy needed). ZnS = mean R² over all pairs (Kelly 1997). Za = mean R² over adjacent biallelic pairs (Rozas 2001). ZZ = Za − ZnS.

**Recombination module**: Four-gamete test (Hudson & Kaplan 1985)  -  a pair of biallelic sites is incompatible when all four gamete combinations are observed. Rm = minimum number of recombination events, computed by the interval-stabbing greedy algorithm (sort incompatible intervals by right endpoint; place a recombination point at the right endpoint whenever the left endpoint exceeds the last placed point).

**Mismatch module**: Observed pairwise-difference histogram. Raggedness r (Harpending 1994, eq. 1) = Σ(f(i) − f(i−1))². CV = σ/μ of pairwise differences (Rogers & Harpending 1992). Small r → smooth distribution → population expansion signature.

**InDel module**: InDel event = maximal run of columns where the same subset of sequences carries gaps (diallelic option of DnaSP). Statistics on InDel haplotypes, k(i), π(i), θ_W(i), Tajima's D(i) computed as for nucleotide data.

**Divergence module**: Dxy = average between-population differences per site (Nei 1987, eq. 10.20). Da = Dxy − (π₁ + π₂)/2 (net divergence). Fixed / shared / exclusive sites classified as in DnaSP 6 (`Divergencia.vb::Mod3BuscaShareFixDifferences`): a site with no shared allele between the populations is a **fixed-difference site** (Sf, a site count), and any within-population variation on top of it is exclusive; a site polymorphic in exactly one population is **exclusive** to that population (Sx, counted as mutations = n_alleles − 1); a site polymorphic in both is split into shared (Ss) and exclusive mutations by a lookup on (alleles in pop1, alleles in pop2, alleles in the union). Ss / Sx1 / Sx2 are mutation counts; the total mutations in a population = Sx + Ss. Complete deletion applied across both populations combined.

**Fu & Li outgroup module (fuliout)**: Outgroup sequence polarises each segregating site  -  allele matching outgroup is ancestral; others derived. S = number of orientable segregating sites (one per site, DnaSP's `pv1`); η = total derived mutations (outgroup-polarised, extra per site for a tri-/quadri-allelic site; reported but not used to scale D/F); η_e = derived-mutation SITES carried by exactly 1 ingroup sequence (external branches), capped at 1 per site even when more than one derived allele there is a singleton. Formulas match DnaSP 6's default mode (`FULI.vb::Mod12FuLiOutgroupNew`, `SSOrMutations=1`): D = (S − aₙ·η_e) / √(u_D·S + v_D·S²); F = (k̄ − η_e) / √(u_F·S + v_F·S²); c_n = 2(n·aₙ − 2(n−1))/((n−1)(n−2)); v_D = 1 + (aₙ²/(bₙ+aₙ²))(c_n − (n+1)/(n−1)); u_D = aₙ − 1 − v_D; v_F, u_F per Fu & Li (1993) / Simonsen et al. (1995). k̄ is the mean pairwise difference over the **orientable-site set only** (outgroup clean, ingroup clean, ancestral allele present)   -   the same column mask as S/η/η_e, matching DnaSP's `SitioIesInformativo` gate. Negative D or F → excess of external (singleton) mutations.

**HKA test module (hka)**: The DnaSP 6 two-locus model (`HKA.vb::HKAResolEcuacion` case 1). Given one species' polymorphism (S_i, sample size n_i, L_i sites) and divergence to a sister species (D_i, L_div_i sites) at **exactly two loci**, the neutral model has parameters θ₁, θ₂ (per site) and a scaled divergence time T. θ₁ is the positive root of a quadratic in θ₁; θ₂ and T follow; roots giving a negative θ are rejected. Goodness of fit is χ² with **df = 1** using the HKA (1987) variances (`HKAJiCuadrado`): E[S_i] = aₙᵢ·sexᵢ·θᵢ·Lᵢ, Var[S_i] = E[S_i] + sexᵢ²·bₙᵢ·θᵢ²·Lᵢ² (the bₙ term); E[D_i] = (T + sexᵢ)·θᵢ·Ldivᵢ, Var[D_i] = E[D_i] + sexᵢ²·θᵢ²·Ldivᵢ². P-value via the regularised upper incomplete gamma Q(1/2, χ²/2)  -  no scipy. If the equations have no positive-θ solution the test is not run (HKAStats.error is set); no divergence time is fabricated.

**McDonald-Kreitman test module (mk)**: `--genetic-code` selects the codon table (default `standard`; `vertebrate-mitochondrial` for mtDNA loci, where TGA is Trp, ATA is Met, AGA/AGG are stop). Counts are per nucleotide site, not per codon, following DnaSP 6's own routines (`Módulos/McDonaldK.vb`). For each codon (in-frame, complete deletion at codon level -- any non-ATCG in any sequence skips that codon; a stop codon under the selected table does NOT skip it: as in DnaSP it is analysed as a 21st amino acid, synonymous with another stop and a replacement against any sense codon, and counted in `n_internal_stop_codons`; a final codon that is a stop in every sequence is the terminal stop and is left out, the equivalent of DnaSP's noncoding annotation of it) two independent tallies are made. Within species (Pn/Ps): the ingroup's distinct codons are classified by a port of `BuscaSitiosReemplazamientoMK`: a site segregating for k bases contributes k - 1 changes, each labelled synonymous or replacement by DnaSP's path rules (most synonymous ordering; for two codons two steps apart, an intermediate codon observed in the outgroup decides the path; two codons for the same amino acid are synonymous at every differing site; the four-codon circular-path case is resolved only under nuclear codes). Codons DnaSP does not analyse (three codons differing at all three sites, four codons with a three-base site, five or more codons, or a circular path under a mitochondrial code) are "complex codons": excluded from all four counts and reported. Between species (Dn/Ds): a site is a fixed difference whenever no ingroup sequence carries the outgroup's base (`EntrePobsMod.vb`'s `ht3` = 1), even if the ingroup itself segregates there; its label comes from the ingroup codon(s) closest to the outgroup codon, fewest replacements on ties. The two tallies do not exclude each other: a site that segregates within the ingroup and whose bases are all absent from the outgroup counts once as a polymorphism and once as a fixed difference (DnaSP's help file, codon 13-15 worked example). Derived statistics: α = 1 − (Ds·Pn)/(Dn·Ps); NI = (Pn·Ds)/(Dn·Ps); DoS = Dn/(Dn+Ds) − Pn/(Pn+Ps); α and NI are reported only when Dn·Ps > 0, as in DnaSP. Fisher's exact P computed via hypergeometric distribution using `math.lgamma` (no scipy needed); two-tailed (sum of all table probabilities ≤ observed probability).

**Fu's Fs module (fufs)**: Estimates θ_π = k (mean pairwise differences, from the polymorphism module). Uses the Ewens sampling formula for the number of distinct alleles K_n in a sample of n under the infinite-alleles model: P(K_n = k) = |s(n, k)| × θ^k / θ^(n), where |s(n, k)| are unsigned Stirling numbers of the first kind and θ^(n) is the rising factorial. Evaluated in **log space** (`math.log` of the exact Stirling integer) so the central coefficients, which exceed float range for n ≳ 171, do not overflow. S′ = P(K_n ≥ H_obs | θ_π, n)   -   the **upper** tail (Fu 1997). Fs = ln(S′ / (1 − S′)); large negative Fs → more haplotypes than expected → population expansion or hitchhiking. No significance level is reported: S′ is not a P-value because θ_π is estimated; a formal test needs coalescent simulation of the null (planned for a future release).

**SFS module (sfs)**: For each alignment column (after complete deletion of the ingroup), counts how many sequences carry each allele. **Only biallelic columns contribute** (DnaSP `FULI.vb` gates on exactly two states); multiallelic columns are excluded and tallied in `n_multiallelic_excluded`. Folded SFS: records sites by minor allele count i (1 ≤ i ≤ n//2). Unfolded SFS (requires `--outgroup`): for each biallelic column where the outgroup allele is present in the ingroup, counts the ingroup sequences carrying the derived allele (i = 1 to n−1). Gap/ambiguous in any ingroup sequence → column excluded; gap in outgroup → excluded from unfolded only. Produces folded and (optionally) unfolded bar-chart figures.

**Ka/Ks module (kaks)**: also honours `--genetic-code`. Without `--outgroup`, pairs are every ingroup-vs-ingroup combination (DnaSP's own manual: "for any pair of sequences", the plain no-population-structure case). With `--outgroup <seq_name>`, pairs become each ingroup sequence vs the outgroup only -- this matches DnaSP's own behaviour once an outgroup is defined via Define Sequence Sets (confirmed on a real DnaSP 6.12 GUI comparison: ingroup-vs-ingroup averaging was 3-3.6x off from DnaSP's own Ka/Ks figures). The outgroup is also folded into the per-sequence site averages in that case (DnaSP: "the total number of synonymous and nonsynonymous sites... is estimated as the average... of all sequences"). Counting follows DnaSP 6's own routines (`SINONIMO.vb`, `EntrePobsMod.vb`), read line by line for the COII calibration. Sites: per sequence, each codon position contributes the fraction of its single-base alternatives that are synonymous, with alternatives that would create a stop codon left out of the denominator (DnaSP's `ComputeFoldPos`: one stop path and one synonymous of the two remaining = 1/2 a site, TGT being DnaSP's own example; two stop paths and a synonymous third = a whole site) -- under the vertebrate mitochondrial code this makes every AGY (Ser) and TAY (Tyr) third position a whole synonymous site; nonsynonymous sites are 3 x codons analysed - S; S and N are means over sequences. Differences: per pair, synonymous (sd) and nonsynonymous (nd) differences by pathway averaging over all k! orderings when codons differ at k positions, paths through stop codons excluded (DnaSP's `NumSynonEntreCodons`); Sd and Nd are means over pairs. Rates: Ks = -3/4 ln(1 - 4pS/3) with pS = Sd/S and Ka likewise with pN = Nd/N -- ONE Jukes-Cantor correction of the ratio of mean differences to mean sites, as DnaSP's `PolDivergenceOut` does, not a mean of per-pair corrected distances (which the convex correction inflates, by about 8% on DnaSP's COII example). A rate is undefined (None) if its p >= 0.75. omega = Ka/Ks, None when Ks = 0. Stop codons follow DnaSP's 21st-amino-acid rule: a stop inside the coding region is analysed (zero synonymous sites, three nonsynonymous; synonymous with another stop, a replacement against a sense codon) and counted in `n_internal_stop_codons`; a final codon that is a stop in every sequence is the terminal stop and is left out. Calibration: on DnaSP's COII example under the mitochondrial code, S 168.222, Ks 0.83011 and Ka 0.02561 match DnaSP 6.12 exactly.

**Ts/Tv module (tstv)**: Counts **one change per biallelic segregating column**, as in DnaSP 6 (`Mutational.vb::Mod31Compute_1` / `RellenoMatrizCambios` / `CalculaTransitionTransversionRatio`)   -   not summed over sequence pairs, so the ratio does not depend on sample size. A **transition** (Ts) is A↔G or C↔T; a **transversion** (Tv) is any purine↔pyrimidine change. Multiallelic columns are excluded (`n_multiallelic_excluded`). With `--outgroup`, DnaSP's polarised mode also drops columns the outgroup cannot orient   -   outgroup gap, or outgroup allele not among the ingroup alleles (`n_unpolarisable_excluded`); the Ts/Tv classification is polarity-independent so the ratio is unchanged for the surviving columns. Ts/Tv = n_transitions / n_transversions; None when n_transversions = 0.

**Fay & Wu / Zeng module (faywu)**: Requires `--outgroup`. Applies complete deletion including the outgroup. For each segregating site, the outgroup allele identifies the ancestral state; a site is **polarisable** when the ancestral allele appears in the ingroup. For each polarisable site with derived allele count i (1 ≤ i ≤ n−1), adds to ξ_i. Computes four per-site θ estimates from the unfolded SFS: θ_π = Σ ξ_i × 2i(n−i) / [n(n−1)] / L; θ_W = Σ ξ_i / a₁ / L (a₁ = Σ 1/k for k=1..n−1); θ_H = Σ ξ_i × 2i² / [n(n−1)] / L; θ_L = Σ ξ_i × i / (n−1) / L. H = θ_π − θ_H (Fay & Wu 2000); E = θ_L − θ_W (Zeng et al. 2006). H < 0 indicates an excess of high-frequency derived alleles (consistent with recent selective sweep). E < 0 indicates excess of low-frequency derived alleles relative to Watterson expectation.

**Fst module (fst)**: Requires `--pop-file`. Applies complete deletion across all sequences from all populations combined. For each pair of populations A and B, computes: π_A = mean within-pop pairwise differences per site; π_B analogous; π_AB = Dxy (mean between-pop pairwise differences per site). Hudson et al. (1992) estimator: Fst = 1 − π_s / π_AB where π_s = (π_A + π_B) / 2. Fst is clamped to [0, 1] (negative values from small samples are set to 0). Fst = None when π_AB = 0 (no between-population variation at any site). Mean Fst is the unweighted average across all pairs. For three or more populations, all pairwise combinations are computed.

**Codon usage module (codon)**: also honours `--genetic-code`; the RSCU family groupings and ENC degeneracy classes follow the selected table (a mitochondrial locus regroups Ile/Met/Trp/Arg). ENC reports n.a. under a non-standard code whose class sizes don't match Wright's (1990) standard-code coefficients (9/1/5/3) rather than compute a wrong number  -  RSCU is unaffected. Reads the in-frame coding alignment in non-overlapping triplets. Triplets with any non-ATCG character or translating to a stop codon are skipped. Codon counts are pooled across all sequences. RSCU (Sharp & Li 1987): RSCU_ij = X_ij / (X_i / n_i) where X_ij = count of codon j for amino acid i, X_i = total count for amino acid i, n_i = synonymous family size. RSCU = 1.0 → uniform usage; > 1.0 → preferred; < 1.0 → avoided. ENC (Wright 1990): computed from mean corrected homozygosity per degeneracy class. For amino acids with k-fold degeneracy (k codons), corrected homozygosity F_k = (n_aa × Σpⱼ² − 1) / (n_aa − 1) where pⱼ = fraction of amino acid i encoded by codon j, n_aa = total codon count for amino acid i. Class means over all amino acids in that class: 2-fold (9 aa), 3-fold (Ile only, 1 aa), 4-fold (5 aa), 6-fold (3 aa). ENC = 2 + 9/F̄₂ + 1/F̄₃ + 5/F̄₄ + 3/F̄₆. Clamped to [20, 61].

**Key thresholds**:
- Tajima's D, Fu & Li D*/F*: require n ≥ 3 and S > 0; return `n.a.` otherwise.
- LD statistics: require ≥ 2 strictly biallelic sites.
- Divergence: require ≥ 1 sequence per population and L_net > 0.
- fuliout: requires n ≥ 4 ingroup sequences and S > 0 (at least one orientable segregating site).
- hka: requires exactly 2 loci, each with n ≥ 2 and positive site counts (L_poly, L_div). HKAStats.error is set (test not run) for the wrong number of loci, bad inputs, or equations with no positive-θ solution.
- mk: requires n ≥ 2 ingroup sequences, alignment length divisible by 3 (in-frame coding), and `--outgroup <seq_name>`. Returns None if outgroup not provided or alignment not in-frame. α, NI, DoS are None when any denominator is zero.
- kaks: requires n ≥ 2 sequences and alignment length divisible by 3 (in-frame coding). ω = None when Ks = 0. Ks or Ka is None when its mean ratio pS or pN ≥ 0.75 (JC saturation).
- fufs: requires n ≥ 2 and H ≥ 1. If k = 0 (all sequences identical), Fs is defined but θ_π = 0 → degenerate. Uses polymorphism stats already computed; no additional inputs needed.
- sfs: requires n ≥ 2. Folded SFS is always produced from the ingroup. Unfolded SFS requires `--outgroup` and at least one site where the outgroup allele appears in the ingroup.
- tstv: requires n ≥ 2 and at least one biallelic segregating column. Multiallelic columns are excluded; with `--outgroup`, columns the outgroup cannot orient are excluded. Ts/Tv = None when n_transversions = 0.
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
| Fu & Li D or F < 0 (outgroup) | Excess external (singleton) mutations → selective sweep or population expansion | Compare with no-outgroup D*/F* |
| Fu & Li D or F > 0 (outgroup) | Fewer singletons than expected → balancing selection or population subdivision | |
| HKA P < 0.05 | Polymorphism/divergence ratio differs between the two loci → departure from neutrality | One locus may be under selection |
| HKA P > 0.05 | Ratio consistent between the two loci → consistent with the neutral model | |
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
| Fs << 0 (Fu's Fs) | More haplotypes than expected given π → population expansion or genetic hitchhiking | No significance level reported (needs coalescent simulation) |
| Fs ≈ 0 (Fu's Fs) | Haplotype count consistent with neutral expectation | |
| Fs > 0 (Fu's Fs) | Fewer haplotypes than expected → balancing selection, population subdivision, or recent bottleneck | Rarely significant |
| SFS singleton-heavy (i=1 dominant) | Excess rare variants → expansion, purifying selection, or recent bottleneck recovery | Consistent with negative Tajima's D |
| SFS flat or U-shaped | Uniform or high-frequency-skewed variants → balancing selection | Consistent with positive Tajima's D |
| Unfolded SFS high at n−1 | Many near-fixed derived alleles → directional selection or recent sweep ancestry | |
| Ts/Tv ≈ 2 | Typical transitional bias for nuclear DNA  -  transitions more mutable than transversions | Expected baseline; varies by locus and taxon |
| Ts/Tv > 10 | Strong transition bias → common in mitochondrial DNA or highly constrained sequences | |
| Ts/Tv < 1 | Transversion excess → can be genuine for AT-rich / fast-evolving non-coding DNA, or transitional saturation | Check base composition and alignment quality |
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
