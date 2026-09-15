---
name: dnasp
description: Population genetics of pre-aligned DNA sequences or multi-sample VCFs
  using selected DnaSP 6 methods. Use for diversity, neutrality statistics, linkage
  disequilibrium, InDel polymorphism, divergence, MK, Ka/Ks and codon usage; not alignment,
  phasing or clinical interpretation.
license: MIT
metadata:
  version: 0.5.2
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
    description: Aligned DNA sequences (pre-aligned, equal-length). FASTA (including
      DnaSP-style >'name' [comment] headers) or NEXUS (MATCHCHAR, INTERLEAVE).
    required: false
  - name: vcf
    type: file
    format:
    - vcf
    description: Multi-sample VCF (--vcf). Converted to one aligned MSA per CHROM
      (biallelic SNPs only; phased -> haplotype rows). Alternative to --input. Optional
      --region CHROM, --vcf-merge to pool all CHROMs.
    required: false
  - name: alignment2
    type: file
    format:
    - fasta
    - fas
    - nexus
    - nex
    description: Second-population alignment for divergence analysis (--input2). Alternative
      to --pop-file. Sequences must have same length as --input.
    required: false
  - name: pop_file
    type: file
    format:
    - tsv
    - txt
    description: 'Population assignment file: one row per sequence, tab-separated
      (sequence_name<TAB>population_name). Required for Fst; alternative to --input2 for divergence.'
    required: false
  - name: outgroup
    type: string
    description: Unique sequence identifier removed from the ingroup. Required for
      fuliout, mk and faywu; also polarises SFS/TsTv and selects ingroup-versus-outgroup
      Ka/Ks.
    required: false
  - name: hka_file
    type: file
    format:
    - tsv
    - txt
    description: 'HKA locus file: whitespace-separated, exactly two loci, columns
      locus n S L_poly D [L_div] [chrom]. Required for --analysis hka.'
    required: false
  - name: analyses
    type: string
    description: Comma-separated polymorphism, ld, recombination, popsize, indel,
      divergence, fuliout, hka, mk, kaks, fufs, sfs, tstv, codon, faywu, fst; or all.
      Polymorphism runs with alignment input.
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
    description: 'Codon table for mk/kaks/codon: "standard" or "vertebrate-mitochondrial"
      (TGA=Trp, AGA/AGG=stop, ATA=Met; for COII/cytb/ND-type loci). Default: standard.'
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
    description: Polymorphism and sliding-window TSV only; other modules appear in
      report.md.
  - name: result
    type: file
    format:
    - json
    description: ClawBio result.json envelope (skill, version, input checksum,
      headline summary, the summary.json payload plus artifact list,
      chat_summary_lines and preferred_artifacts for the runner). A multi-CHROM
      VCF run writes one envelope per CHROM directory and a root envelope that
      summarises them with relative paths.
  - name: summary
    type: file
    format:
    - json
    description: Module summaries (global and per-window statistics with DnaSP
      midpoints, LD/recombination/mismatch/InDel/divergence/MK/Ka-Ks/codon results
      including named per-sequence ENC); LD pair grids are in ld_pairs, not here.
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
    description: Archived inputs, commands.sh, environment.yml, manifest.json and
      output-relative checksums.sha256.
  dependencies:
    python: '>=3.10'
    packages:
    - matplotlib>=3.7
    - numpy>=1.24
    - pandas>=2.0
    - opentelemetry-sdk>=1.20
  demo_data:
  - path: examples/demo_simple.fas
    description: Synthetic 6-sequence × 10-bp alignment with known statistics
  - path: examples/demo_rp49.fas
    description: rp49 region, 17 Drosophila sequences, 300 bp
  endpoints:
    cli: python skills/dnasp/dnasp.py --input {alignment} --analysis {analyses} --output
      {output_dir}
  openclaw:
    requires:
      bins:
      - python3
    always: false
    emoji: ''
    homepage: https://github.com/ClawBio/ClawBio
    os:
    - darwin
    - linux
    - win32
    install:
    - kind: pip
      package: matplotlib>=3.7
    - kind: pip
      package: numpy>=1.24
    - kind: pip
      package: pandas>=2.0
    - kind: pip
      package: opentelemetry-sdk>=1.20
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

# DnaSP

## Trigger

Fire when a user requests population-genetic analysis of aligned DNA, a supported
VCF, or a DnaSP-compatible statistic listed below. Do NOT fire for sequence
alignment, read mapping, haplotype phasing, clinical advice or unsupported
coalescent significance tests.

## Scope

Analyse genetic variation in supplied alignments using 16 selected DnaSP methods.
This is not a complete replacement for the DnaSP GUI or all its analysis modes.
Read [the statistical reference](docs/index.md) for definitions, exclusions,
source conventions, file formats, examples and release validation evidence.

## Core Capabilities

- Polymorphism: S, Eta, haplotypes, Hd, VarHd, pi, k, theta-W, G+C,
  Tajima's D, Fu and Li D*/F*, Ramos-Onsins and Rozas R2.
- LD: D, D', R2, ZnS, Za, ZZ and original-column pair labels; recombination Rm.
- Mismatch distribution with unbiased variance over unordered pairs and Sokal &
  Rohlf corrected CV; Harpending (1994) raggedness; Model 1 diallelic InDel diversity.
- Divergence and Hudson Fst between populations; outgroup Fu and Li D/F.
- Two-locus HKA, McDonald-Kreitman, Nei-Gojobori Ka/Ks, Fu's Fs and SFS.
- Ts/Tv, codon counts/RSCU including stops, named per-sequence ENC and its
  synonymous-codon-weighted summary.
- Sliding windows mirror DnaSP's **Gaps in Sliding Window = considered** mode:
  starts advance from 1 by the step, capped at alignment length; ends are capped
  there too, and the first window reaching the end terminates the loop. Empty
  windows are retained. The not-considered mode is not implemented.
- Window midpoints are the alignment position of the ceil(net/2)-th gap-free
  column, or the window start if none exists; exported in TSV, JSON and report
  and used for the plot. VCF windows use retained SNP indices, including the
  final partial window.
- Raw per-site Fay-Wu H and theta-L minus theta-W; normalised Hn/ZE are absent.

Fu and Li D*/F* and outgroup D/F mirror **Data > Segregating Sites/Mutations =
Segregating sites**, using DnaSP's v5-style panels. The rp49 Eta-setting D*/F*
figures can be reproduced by substituting eta for S, but this is not a general
conversion: `FULI.vb` also changes singleton and external-mutation capping
(`SingleMut` and `ExternaMut` subtraction). No Eta-mode switch is implemented.

## Workflow

1. Confirm the input is pre-aligned and identify the scientific comparison:
   ingroup, explicit outgroup, populations, coding interval and genetic code.
   Do not infer an outgroup from record order or silently select a code.
2. Use the full repository checkout or supplied validation code snapshot so the
   shared reproducibility helpers are available. Install the declared dependencies.
3. Choose an empty output directory. For alignment input use `--input`; for VCF
   use `--vcf`; for two-locus HKA alone use `--hka-file --analysis hka`.
4. Select actual implemented names with `--analysis`. Supply `--outgroup`,
   `--input2` or `--pop-file` when needed. Use `--genetic-code
   vertebrate-mitochondrial` for the matching mitochondrial table. The default is
   standard. Coding intervals must be preselected and divisible by three.
5. Execute the CLI. An explicit analysis that cannot run returns a non-zero code.
   `--analysis all` is opportunistic: inspect its completed/skipped manifest.
6. Read stderr diagnostics, `report.md` and `reproducibility/manifest.json`.
   Distinguish a failed analysis, an undefined statistic and an excluded site.
7. Interpret the chosen statistic within its documented assumptions. Do not turn
   a signed value or an uncalibrated threshold into a significance claim.
8. Retain the input archive, settings, hashes and environment with the report.
   `commands.sh` replays on the recorded host/code path into a new output folder;
   moving a run requires the code and dependencies as well as its input archive.
9. For Windows GUI comparison, follow the separate validation checklist and
   capture raw DnaSP output. Source-derived expectations are not GUI observations.

## Example Output

```bash
python skills/dnasp/dnasp.py --demo --output new_demo_run
python skills/dnasp/dnasp.py --input alignment.fas --analysis polymorphism,ld --output new_run
python skills/dnasp/dnasp.py --input coding.fas --outgroup OutSeq --genetic-code vertebrate-mitochondrial --analysis mk,kaks,codon --output new_coding_run
python skills/dnasp/dnasp.py --vcf samples.vcf --analysis polymorphism,sfs,fufs --output new_vcf_run
```

The synthetic demo contains 10 ingroup sequences, one outgroup and 300 sites:

| Quantity | Demo value |
|---|---:|
| S / haplotypes | 5 / 8 |
| Hd / Tajima D | 0.9556 / 0.6789 |
| MK Pn / Ps / Dn / Ds | 2 / 3 / 2 / 1 |
| MK alpha | 0.6667 |
| Ka / Ks / omega | 0.010239 / 0.030291 / 0.3380 |
| Ts / Tv | 4 / 1 |

The demo runs 15 modules; HKA uses a separate two-locus input. The regression
suite, rather than the printed banner alone, checks the expected figures.

## Output Structure

```text
output/
  report.md
  results.tsv                  # polymorphism and windows
  summary.json                 # module summaries, window midpoints, named ENC/null
  result.json                  # ClawBio envelope: headline summary, summary.json payload, artifacts
  ld_pairs.tsv                 # when LD pairs exist
  figures/                     # when matplotlib is available
  reproducibility/
    inputs/
    commands.sh
    environment.yml
    manifest.json
    checksums.sha256
```

## Dependencies

Python 3.10 or later. Core estimators use the standard library. Plotting uses
matplotlib; the repository's shared reproducibility package also imports NumPy,
pandas and OpenTelemetry. Use `environment.yml` or the validation package's
pinned requirements. The Windows validation runner targets Python 3.12.

## Gotchas

- The model will want to count IUPAC symbols as alleles. Do not. Supported
  ambiguity symbols are missing data and the nucleotide mask excludes their columns.
- The model will want to analyse an unknown outgroup as an ingroup-only run.
  Do not. Identifiers must be unique and an explicit outgroup must resolve.
- The model will want to call VCF diversity per-base diversity. Do not. VCF
  columns represent retained variant records; invariant callable bases are absent.
- The model will want to use any coding annotation in a NEXUS file. Do not.
  This implementation requires a preselected coding alignment and does not read
  CHARSET coding annotations. In COII, use positions 1-681 for codon usage.
- The model will want to apply one stop-codon rule everywhere. Do not. Selected
  stops count in RSCU and as family 21 in MK/KaKs; ENC and coding G+C use sense codons.
- The model will want to read raw H/E as normalised DnaSP Hn/ZE. Do not.
  These outputs are explicitly different, and no significance test is supplied.
- The model will want to overwrite an earlier output folder. Do not. Choose a
  fresh folder; runs reject non-empty destinations to preserve prior artefacts.
- The model will want to infer genomic-scale performance from the tiny bundled
  VCF examples. Do not. LD enumerates all biallelic-site pairs and uses quadratic
  storage; read the measured limits in the reference before a large run.

## Safety

All sequence analysis and output remain local. ClawBio is a research and
educational tool. It is not a medical device and does not provide clinical
diagnoses. Consult a healthcare professional before making any medical decisions.

## Agent Boundary

The agent selects documented inputs/options, executes the skill and explains
reported results. The code computes the statistics. Neither the agent nor the
skill may fabricate GUI validation, P-values or missing results.

## Integration with Bio Orchestrator

CLI alias: `dnasp`. Use `--analysis` for module selection, not invented flags
such as `--pi`, `--kaks` or `--n-sim`. The repository dispatcher permits the
implemented options, including VCF, populations and genetic-code selection.

Route here when the user asks for a statistic this skill computes (for example
nucleotide or haplotype diversity, Tajima's D, Fu and Li's tests, Fu's Fs,
McDonald-Kreitman, Ka/Ks, HKA, the mismatch distribution, InDel polymorphism or
codon usage bias) on an aligned FASTA or NEXUS file or a multi-sample VCF.
`INTENTS.json` publishes the `dnasp` aliases to ClawBio's intent planner and
plans the demo only when the user explicitly asks for a demo; a real analysis
needs the user's file, and `--analysis` selects modules other than the default
polymorphism summary.

## Chaining Partners

- **Alignment upstream.** Sequences must be aligned before this skill.
  `phylogenetics-builder` aligns with MAFFT, MUSCLE and other aligners. Versions
  of that skill that save their alignment write `alignment/aligned.fasta`; use
  that untrimmed file rather than the trimAl output, because trimming removes
  gapped or poorly aligned columns and so changes S, pi and the InDel results.
  Otherwise align the sequences separately and pass the aligned FASTA.
- **`fastreer`.** Builds distance trees from the same multi-sample VCF: run this
  skill for diversity and neutrality statistics, then fastreeR for sample
  relationships.
- **`claw-ancestry-pca`.** Gives population-structure context before samples are
  grouped in a `--pop-file` for Fst or divergence.
- **`equity-scorer`.** Reports HEIM heterozygosity and Fst from VCF or ancestry
  data with its own estimators; this skill's Hudson Fst and diversity statistics
  follow DnaSP 6, so the two sets of values are not interchangeable.
- **VCF preparation.** A filtering or phasing workflow may prepare input, but
  every filter and phasing choice must be recorded; unphased heterozygotes stay
  unresolved here.
- **Reporting downstream.** Use `report.md`, `summary.json`, `result.json` and
  the TSV files; the TSV does not contain every module's results.

## Maintenance

Recheck source/help-derived regression cases and the 170 historical comparison
fixtures after changes to formulas, masks or parsers. Review GUI differences
when the target DnaSP build or analysis mode changes. Keep this file, the method
reference, CLI metadata, version and catalogue consistent. New Windows GUI
observations must be reviewed before changing published concordance counts.
