# DnaSP ClawBio version history

Every change listed alters results unless marked otherwise.

**0.6.0** (17 September 2026)

- New, opt-in: coalescent-simulation P-values for Tajima's D (two-tailed),
  Ramos-Onsins and Rozas R2 and Fu's Fs (lower tail) over the whole region, via
  `--n-sim`, `--sim-given S|theta` and `--sim-seed`. The null is the Kingman
  coalescent with constant size, infinite sites and no recombination, conditioned
  by default on the observed number of segregating sites. Without `--n-sim` no
  simulation runs and no statistic changes; the report's interpretation text does
  change, as the next entry describes.
- Not result-changing for any statistic: the report no longer turns a value into
  a significance claim. Tajima's D and R2 were labelled from fixed thresholds
  ("Consistent with neutrality", "suggests population expansion"); they are now
  described as not assessed unless `--n-sim` supplies a P-value.
- Not result-changing: the version history moved from SKILL.md to
  `docs/version_history.md` to keep SKILL.md within the repository's length limit.
- The ClawBio dispatcher forwards `--n-sim`, `--sim-given` and `--sim-seed`; without
  them in its allowlist, `clawbio run dnasp` would have dropped the options silently.
- Two sequences now give R2 and Fu's Fs P-values; only Tajima's D needs three.
- Each P-value is (b + 1)/(N + 1) for b of N valid replicates at least as extreme,
  so none is zero (Phipson and Smyth 2010); summary.json also stores the counts b,
  from which DnaSP's proportion b/N follows.

**0.5.3** (16 September 2026)

- Not result-changing: every input in the skill metadata now names the
  command-line flag it maps to (`cli_flag`). Without it an agent reading only the
  registered metadata had to guess, and two independent models turned
  `window_size` into the non-existent `--window-size`. Numerical behaviour is
  unchanged from 0.5.2, which is the version compared with DnaSP 6.12.03.

**0.5.2** (15 September 2026, compared with DnaSP 6.12.03)

- Sliding windows follow DnaSP's placement: the final window, truncated at the
  alignment end, is kept, and each window reports DnaSP's midpoint. Window
  counts and the last window's values change.
- Mismatch distribution: the observed variance of k is the unbiased variance over
  sequence pairs, and its coefficient of variation carries the (1 + 1/(4n))
  correction. rp49 moves from 40.7054 and 0.3954 to 40.7780 and 0.3987.
- LD: where two alleles are tied, the first sequence's allele is the reference,
  so the sign of D can change; |D|, |D'| and r2 do not.
- Not result-changing: per-sequence ENC export, `summary.json` and `result.json`,
  a root envelope and bundle for VCF runs split by CHROM, and safe CHROM
  directory names.

**0.5.1** (candidate, never released on its own; included in 0.5.2)

- InDel polymorphism follows DnaSP's Model 1 (diallelic) event rules and
  denominator, the model the skill implements; in DnaSP it must be selected, as
  the dialog defaults to Model 2. InDel events, lengths, haplotypes and diversity
  change.
- Fay and Wu H and Zeng E use only clean columns that are at most biallelic in the
  ingroup and, when polymorphic, carry the outgroup allele in the ingroup. The
  eligible-site count and the theta estimates behind H and E change.
- ENC is the synonymous-codon-weighted mean of per-sequence ENC, as in DnaSP,
  instead of a single ENC from pooled codon counts.
- RSCU counts stop codons and reports the stop family.
- LD distances use original alignment coordinates with DnaSP's integer gap
  adjustment.
- Fu's Fs computes both tails in log space without finite sentinel values, which
  changes extreme values.
- IUPAC ambiguity symbols R, Y, S, W, K, M, B, D, H and V are missing data, like
  N; any other symbol now stops the run with an error.
