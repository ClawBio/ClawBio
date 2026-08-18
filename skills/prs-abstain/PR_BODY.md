## What this adds

`prs-abstain` decides whether a polygenic score percentile may be reported for a given
individual, and withholds it when it cannot be interpreted. It does not infer ancestry
(`ancestry-risk-profiler`, `claw-ancestry-pca`) and does not compute scores (`gwas-prs`).
It sits between them.

Built for Challenge 3, ClawBio + Nebius Berlin, 18 August 2026.

## Why

Every curated score in the library discloses `Reference population: EUR` and still returns a
confident percentile for any genome. Disclosure is not abstention.

While building this we found the reference mean is exactly `Sum(2 * AF_EUR(i) * w(i))` for all
six curated scores (PGS000013 1.1186 vs 1.12; PGS000004 2.8428 vs 2.84; PGS000057 7.1080 vs
7.11). The percentile is centred on a European allele-frequency calculation, so the error
decomposes per variant rather than being a property of the whole score.

## Findings from the existing demo data

- `gwas-prs --demo` returns a breast cancer **and** a prostate cancer percentile for the same
  individual. At least one is inapplicable.
- PGS000001 and PGS000039 each carry two different rsIDs at one coordinate
  (`chr1:114448389` = rs11552449 + rs7072776). Coordinate error or a locus counted twice;
  either corrupts the sum. This blocks the score.
- PGS000013's two largest weights (rs7903146, rs12255372) are one locus, TCF7L2, 50 kb apart
  and 48% of the score. Effective n is 3.5, not 6.
- PGS000057 looks like 147 independent contributions and behaves like 35.9.

## Checks, in order

1. Trait applicability (sex-specific scores)
2. Score integrity (weight coverage, concentration, strand-ambiguous, duplicate positions)
3. Ancestry: markers >= 30 (Kosoy et al. 2009, PMID 18683858), placeable, distance <= threshold
4. Reference-population match

Raw scores are retained in every branch. Only the percentile is withheld: summing weights is
valid for anyone, converting that sum to a percentile needs a reference distribution that fits.

## Validation

- 59 tests, written red/green before implementation
- 17/17 SKILL.md conformance; 10/10 on the submit-page validator
- 120 tests pass across `prs-abstain`, `gwas-prs`, `equity-scorer` (no regressions)
- Stress-tested 10 ways; three regressions found and fixed, all now Gotchas
- `--k-sd 20` is detected and stamped THRESHOLD OVERREACH rather than silently passing
  non-European individuals

## Demo

```bash
python skills/prs-abstain/prs_abstain.py --demo --output /tmp/prs-abstain
```

Three individuals, three verdicts, two of them deliberate refusals. Writes a jargon-free
clinician report and a technical report, both also as PDF.

## Honest limitations

LD is approximated by physical distance, not measured, so the grouped effective_n remains an
overestimate. Admixed individuals are the principal unhandled case. The bundled allele-frequency
table is synthetic and labelled as such in the file. Twelve limitations are listed in the
technical report and in SKILL.md.
