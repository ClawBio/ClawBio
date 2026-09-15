# DnaSP ClawBio statistical reference

Version 0.5.2 implements 16 selected population-genetic analysis modules. Its
methods combine published estimators with documented DnaSP 6 implementation
conventions. Agreement on the validation datasets does not establish equivalence
for every DnaSP mode, input or statistic.

## Inputs and analysis scope

FASTA and DNA NEXUS inputs must be pre-aligned, have unique identifiers and equal
sequence lengths. Lower-case bases are normalised to upper case. A, C, G and T
are nucleotides; supported IUPAC ambiguities, N and ? are missing data. Unsupported
symbols are rejected. Nucleotide complete deletion removes a column if any
sequence in the relevant sample has a gap or missing base. Coding, outgroup and
population analyses use their documented masks rather than one universal mask.

The NEXUS reader supports a single DNA/nucleotide MATRIX, quoted labels,
sequential or interleaved data, MATCHCHAR and declared missing/gap symbols. It
checks NTAX/NCHAR when present and validates reconstructed lengths. Unlabelled
continuation blocks need an unambiguous preceding taxon block. Transposed,
tokenised and NOLABELS formats are rejected. Coding-region annotations such as
CHARSET are not interpreted. MATCHCHAR expands against the first sequence before
outgroup removal; do not reorder records to put the outgroup first.

An explicit `--outgroup` must match exactly one sequence. That record is removed
from the ingroup for the complete run. It supplies ancestral states for SFS,
Ts/Tv, outgroup Fu-Li and Fay-Wu, the comparison sequence for MK, and the
outgroup comparisons for Ka/Ks. Without it, Ka/Ks averages ingroup pairs.

Population input is either `--input2` or `--pop-file`. Two alignments must cover
the same columns. The population file maps sequence identifiers to population
labels. Divergence uses the first two population labels in sorted order; Fst
reports all population pairs. Unassigned records do not contribute to the
population comparisons. Fay-Wu requires an outgroup, not a population file.

Coding analyses require a preselected alignment starting at codon position one,
with length divisible by three. Codon usage processes the selected triplets;
ambiguous/gapped triplets are omitted per sequence. MK and Ka/Ks recognise a
shared terminal stop according to their existing coding-selection rules.
NEXUS annotations are not a substitute for selecting the coding interval.
For the COII example, positions 1-681 are coding, whereas the whole nucleotide
alignment has 684 sites. Ribosomal 12S/16S regions are not protein-coding examples.

The default genetic code is standard. `--genetic-code vertebrate-mitochondrial`
selects the DnaSP mtDNA Mammals convention: TGA is Trp, ATA is Met, AGA and AGG are
stops. Selected internal stops are treated as a 21st amino-acid family by MK and
Ka/Ks, with a warning. RSCU includes the stop family. ENC excludes stops.

## Module definitions

### Polymorphism

`polymorphism` runs on alignment input and is the default. After ingroup complete
deletion, S counts segregating columns and Eta counts the sum of (number of
nucleotide states minus one) over those columns. They differ at multiallelic sites.

Haplotype number H counts distinct cleaned sequences. Haplotype diversity is
n/(n-1) times (1 minus the sum of squared haplotype frequencies), with the
Nei/Tajima sampling variance used by the implementation. k is the mean number of
pairwise nucleotide differences; pi = k/L_net. Watterson theta is S/a1, with
per-site theta-W = S/(a1 L_net), where a1 = sum(1/i), i = 1,...,n-1. G+C is the
fraction of clean nucleotide characters that are G or C.

Tajima's D uses k and S with the Tajima (1989) variance. Fu and Li D*/F* use the
standard Simonsen et al. (1995) variance convention, with the DnaSP singleton-site
count capped at one per site. D*/F* and outgroup D/F mirror **Data > Segregating
Sites/Mutations = Segregating sites**: compare the v5-style Segregating-sites
panels (the second panels in `A_rp49_all.out`). Substituting eta=92 for S=89
reproduces that capture's rp49 Eta D*/F*, but is not a general Eta-mode rule:
`FULI.vb` also changes singleton and external-mutation capping through the
`SingleMut` and `ExternaMut` subtractions. The default is unchanged and there is
no Eta-mode switch. R2 follows the DnaSP caller's S denominator, not
Eta. Small-sample, monomorphic and zero-denominator cases can be undefined.
Signed neutrality statistics are descriptive; uncalibrated cutoffs do not supply
P-values. `theta-L` belongs to Fay-Wu, not this module.

### Linkage disequilibrium

`ld` uses strictly biallelic columns after complete deletion and enumerates all
pairs. It reports D, normalised D', R2, chi-square and its one-degree-of-freedom
P-value. ZnS is mean R2 over all pairs; Za is mean R2 over adjacent eligible
sites; ZZ = Za - ZnS. These LD association P-values are distinct from neutrality
P-values.

The sign of D and D' follows DnaSP's allele ordering (`CODIGO2.vb::calculo_mas_freq1`):
the major allele is the reference, and when the two alleles are equally frequent
the first analysed sequence's allele takes that role. |D|, |D'| and R2 do not
depend on the ordering.

Site1/Site2 retain original, one-based alignment columns. DnaSP's inline LD loop
in `CODIGO2.vb` calculates Dist as CInt(j-i-gap_count/n), with gaps summed across
the inclusive original interval in the selected sequences. CInt rounds ties to
even. Missing bases are not gaps. The linkage-disequilibrium help example has
sites 1 and 18, Dist 13; full rp49 has sites 2 and 69, Dist 64. The unrelated
`DistanciaNucleotidica` routine is not the LD-grid distance calculation.

VCF columns are retained variant records, not genomic positions; distance and
sliding-window labels on that input must be interpreted in variant-column units.
Cross-chromosome merged VCF output has no meaningful physical pair distance.
The grid does not implement LD R, Fisher's exact P or the Bonferroni flag.
Some small-sample DnaSP modes suppress ZnS; this implementation does not impose
an unverified cutoff.

### Recombination

`recombination` uses biallelic four-gamete incompatibilities. Rm follows the
two-pass interval reduction in `CODIGO2.vb::RecombinacionRM`, including its rule
for intervals that meet at a boundary. It is not a generic greedy interval-
stabbing solution. Pair labels retain original alignment columns. Restoring
coordinates does not change the incompatibility counts or Rm.

### Population size changes

`popsize` reports the pairwise mismatch distribution, mean, variance,
coefficient of variation and Harpending (1994, equation 1) raggedness. Its
distribution is computed from the cleaned ingroup. Variance is the sum of
squared deviations from the mean divided by (number of unordered pairs minus
1), or zero for two sequences. CV is `(1 + 1/(4n)) * sqrt(variance) / mean`,
where n is the number of sequences, following DnaSP's Sokal & Rohlf correction
(`PairwiseDiff.vb`, lines 565-567 and 732); CV is undefined for a zero mean.
It does not fit a demographic model or compute a
coalescent P-value. A smooth or ragged curve alone does not identify a population
history.

### InDel polymorphism

`indel` implements DnaSP Model 1, diallelic non-overlapping events. Gap runs are
followed within sequences and events sharing a start and length share carriers.
Distinct overlapping events cause the affected fragment to be excluded;
adjacent non-overlapping events remain separate. Missing columns are ignored
inside the source-defined fragment traversal. Fixed-gap and missing sites do not
contribute to analysed positions.

The report gives accepted and excluded events, mean event length, carrier-
weighted mean deletion length, InDel haplotypes/diversity, k(i), pi(i), theta(i)
and Tajima's D(i). k(i) is the sum of binary event diversities. Its per-site
normalisation uses non-InDel sites plus accepted InDel sites, not gap sites alone.
Theta(i) is I/a1 per sequence. InDel Tajima's D requires at least four sequences.

Source: `DNAPolymorphism.vb::Mod34_Compute`,
`Mod34_BuscoNumEventosEnInDel`, `mod34_HacerAnalisis` and
`mod34_ComputeInDelStatistics`. The 13-sequence help example gives I=2,
net sites=11, k(i)=0.435897, pi(i)=0.039627 and D(i)=-0.909202.
Other DnaSP InDel models are not implemented.

### Divergence and population differentiation

`divergence` reports within-population pi, Dxy, net divergence Da, fixed-
difference sites, shared mutations and population-exclusive mutations.
The DnaSP shared/private quantities are mutation counts, not simply site counts;
multiallelic columns can contribute more than one mutation.

`fst` implements the Hudson pairwise estimator 1 - mean(pi_A, pi_B)/Dxy_AB.
The mask spans all selected populations. If Dxy is zero, Fst is undefined;
reports and plots preserve that missing value rather than displaying zero.
Negative estimates are possible. Pairwise values are averaged for the mean Fst.

### Outgroup Fu and Li

`fuliout` reports D and F using orientable segregating sites and the DnaSP
singleton-site convention. The total polarised mutation count is also reported;
it is not substituted for the S-based test denominator. The standard-mode
comparison and biallelic/Achaz mode are different reference targets. No
significance is asserted from an arbitrary absolute-value threshold.

### HKA

`hka` follows DnaSP's two-locus model in `HKA.vb`: one species' polymorphism and
between-species divergence, solving for locus-specific theta and shared scaled
divergence time. The chi-square has one degree of freedom. The input format is:

```text
locus n S L_poly D [L_div] [chrom]
a 11 10 1000 20 1000 A
b 11 5 1000 10 1000 A
```

A header and # comments are optional. n must be at least two; S and D must be
non-negative integers; lengths must be finite and positive. An omitted L_div
means L_poly; zero is not an omission. Chromosome values are A/autosomal,
X/Z/XL or Y/W, with scaling factors 1, 0.75 and 0.25. Malformed rows and unknown
chromosome tokens are errors. Unidentifiable models are reported as not computed.

### McDonald-Kreitman

`mk` uses an ingroup and outgroup coding alignment. Pn/Ps/Dn/Ds count
nonsynonymous/synonymous polymorphic and fixed nucleotide changes. Within-species
and fixed-difference counts are independent: a site can contribute to both.
The port follows `McDonaldK.vb` per-site pathway classification and multiallelic
rules, not whole-codon fixed-versus-polymorphic heuristics.

Reported summaries include neutrality index, alpha, direction of selection
(DoS) and two-sided Fisher's exact P. Undefined ratios remain missing. The
single-outgroup and count-equivalent pathway tie-break limitations remain.

### Ka/Ks

`kaks` uses Nei-Gojobori synonymous/nonsynonymous site and difference counts with
DnaSP's `SINONIMO.vb` and `EntrePobsMod.vb` conventions. It averages differences
and site counts before applying Jukes-Cantor, rather than averaging corrected
pairwise rates. An explicit outgroup selects ingroup-versus-outgroup comparisons.
Saturation makes corrected rates undefined. omega = Ka/Ks uses full precision.

The implementation's per-sequence/per-pair treatment of gapped codons differs
from alignment-wide deletion used in some DnaSP comparisons. Use a common coding
mask for a directly comparable run; this release does not claim universal
coding-mask equivalence.

### Fu's Fs

`fufs` uses the Ewens probability of at least the observed number of haplotypes,
conditional on theta estimated by k. Fs is log(upper tail) minus log(lower tail).
Both tails are accumulated in log space from unsigned Stirling numbers so extreme
probabilities do not produce arbitrary finite sentinels. S_k is the upper-tail
probability, not a calibrated neutrality P-value. Negative Fs describes an
excess of haplotypes relative to this model. A formal test requires simulation.

### Site frequency spectrum and Ts/Tv

`sfs` counts strictly biallelic segregating sites. The folded spectrum uses
minor-allele counts. The unfolded spectrum additionally requires a clean outgroup
allele present in the ingroup. Multiallelic exclusions are counted.

`tstv` classifies one transition or transversion per eligible biallelic site,
with polarisation when an outgroup is supplied. It does not sum substitutions
across every sequence pair. The validation GUI build did not expose a dedicated
Ts/Tv panel, so these rows are excluded from the historical concordance count.

### Codon usage

`codon` reports mean triplet counts and pooled RSCU across the selected sequences.
RSCU = observed codon count / (amino-acid-family count / family size), with all
64 codons displayed. Stops belong to family 21 in
`CodonUsage.vb::MuestraRSCU`; they are not silently dropped.

ENC is computed per sequence using `M23ENC`, then averaged using synonymous-
codon counts as weights. Degeneracy classes follow the chosen genetic code.
Unusable family homozygosities are omitted; a missing three-fold estimate may
be interpolated from the two- and four-fold classes. The source cap at 61 is
retained. Stop codons are excluded from ENC. Missing required class estimates
remain undefined; pooling sequences before estimating ENC is not equivalent.
`summary.json` exports `codon.per_sequence_enc`, mapping each analysed sequence
name to ENC or JSON null; the report includes the same values at three decimals.
An explicitly selected outgroup is excluded. The summary ENC remains the
synonymous-codon-weighted mean.

For the 11-sequence COII ingroup restricted to 1-681, independent source-derived
weighted means are 46.977002 (standard) and 46.780546 (vertebrate mitochondrial).
The round-5 `M_ENC_standard.out` and `M_ENC_vertebrate-mitochondrial.out`
summaries agree at three decimals. The 11 + 11 per-sequence values are checked
against the labelled `M_ENC_*_GUI_tables_transcribed.txt` files; those tables
are operator-screen transcriptions, not DnaSP file exports.

### Fay-Wu and Zeng summaries

`faywu` reports raw per-site H = theta-pi - theta-H and the raw numerator
theta-L - theta-W. These are not normalised Hn or ZE. Eligible columns must be
clean in the ingroup and outgroup and cannot be multiallelic; polymorphic columns
must contain the ancestral state. Counts, diversity and per-site denominators
use this same eligible set. The module does not supply normalised-test P-values.

## Windows and reproducibility

`--window` and `--step` mirror **Gaps in Sliding Window = considered**. Coordinates
are one-based inclusive: start at 1, advance by the step with the next start
capped at alignment length L, and use `end = min(start + window - 1, L)`.
Stop after the first window reaching L, retaining that truncated window and
its actual length. Windows with no analysable sites are kept. This follows
`CODIGO2.vb` and `CONTROLE.vb::BuscaDesplazamientoSW/BuscaTamanyoVentanaSW`.
The not-considered mode, which discounts gaps in placement, is not implemented.

`Midpoint` in `results.tsv`, `midpoint` in each `summary.json` window and the
report table give the alignment position of the ceil(net/2)-th gap-free column,
or the window start when none exists (`CONTROLE.vb::BuscaPuntoMedioWithSynSW`).
The window figure uses these midpoints on its x axis. On VCF, windows contain
retained variant records, not genomic bases, and also retain the final partial
window. Tajima's D remains undefined for S=0 windows even when DnaSP prints 0.0000.

Explicitly requested analyses fail with a non-zero exit when prerequisites are
missing. `all` runs applicable analyses and records skipped modules/reasons in
both report and manifest. Mathematically undefined estimates on otherwise valid
data are not the same as failed input validation. Non-empty output destinations
are rejected to prevent stale figures and overwritten reports.

`results.tsv` contains polymorphism and windows only. Other modules are reported
in Markdown and `summary.json`; LD also has `ld_pairs.tsv` (pair grids are
omitted from the JSON summary). `result.json` is ClawBio's structured envelope
(skill, version, input checksum, a headline summary, the `summary.json` payload
with the artifact list, and the chat lines and preferred artifacts the ClawBio
runner promotes); it is written before the reproducibility bundle so the
checksums cover it. A multi-CHROM VCF run (no `--region`, no `--vcf-merge`)
writes one envelope per CHROM directory and a root `result.json` listing each
run's headline statistics and artifacts by relative path. Names taken from
input files or CHROM identifiers are stripped of control characters and markup
before they appear in chat lines. The reproducibility folder uses ClawBio's
shared writers and records the actual arguments, versions, archived input bytes,
analysis status and output-relative hashes. `commands.sh` quotes arguments and
writes to a fresh replay destination, but assumes the recorded interpreter and
code path. A portable PC package must also include code and dependencies.

The first two TSV lines are metadata and a blank separator. In R, for example,
use `read.delim("results.tsv", skip=2, check.names=FALSE)` and inspect headers.

## Validation evidence and limits

The checked-in historical fixture maps all 170 recorded S1 comparisons to input
files, settings and expected values. The recorded assessment is 163 matches,
three last-displayed-digit differences and four F* mode differences. Unit tests
protect those observations; they do not constitute a new DnaSP GUI run. The
four historical F* Python cells were padded to six decimal places after rounding
to four; their regression precision is explicitly recorded as four decimals.

Additional regressions cover the InDel help example, original LD coordinates,
RSCU stops, ENC missing classes, extreme Fs tails, ambiguous/invalid inputs,
NEXUS labels and dimensions, output protection, figure semantics and replay.
The next Windows round must capture new GUI observations before these are added
to concordance totals. A passing test suite is not proof of all-mode equivalence.

LD requires n times the square of the biallelic-site count in comparisons and
quadratic result storage. On the review machine, n=200 and 500 biallelic sites
produced 124,750 pairs in 3.349 seconds. At 5,000 biallelic sites there are
12,497,500 pairs; scaling that timing suggests roughly 335 seconds and several GB
of Python object storage, not a measured benchmark at that size. No pair sampling
or distance restriction is silently applied. The tiny bundled VCF examples do
not establish genomic-scale performance.

## References

Primary source verification for the corrected bibliographic details:

- Rozas et al. (2017). DnaSP 6: DNA Sequence Polymorphism Analysis of Large Data
  Sets. Molecular Biology and Evolution 34:3299-3302.
  https://doi.org/10.1093/molbev/msx248
- Rozas, Gullaud, Blandin and Aguade (2001). DNA Variation at the rp49 Gene Region
  of Drosophila simulans: Evolutionary Inferences From an Unusual Haplotype
  Structure. Genetics 158:1147-1155. https://doi.org/10.1093/genetics/158.3.1147
- DnaSP 6.12 documentation: https://www.ub.edu/dnasp/DnaSP6_Documentation_6.12.pdf

The generated report cites the original estimator literature: Tajima (1989),
Fu and Li (1993), Simonsen et al. (1995), Nei and Tajima (1981), Ramos-Onsins and
Rozas (2002), Kelly (1997), Hill and Robertson (1968), Lewontin (1964), Hudson
and Kaplan (1985), Rogers and Harpending (1992), Harpending (1994), Nei (1987),
Hudson et al. (1987, 1992), McDonald and Kreitman (1991), Nei and Gojobori (1986),
Fu (1997), Sharp and Li (1987), Wright (1990), Fay and Wu (2000) and Zeng et al.
(2006). Targeted VB routines and help-derived expectations are recorded beside
the corresponding code and regression cases.
