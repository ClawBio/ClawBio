# DEMO RUNBOOK — prs-abstain

**3 minutes. One terminal, one PDF open. Everything runs live.**

Pre-stage: `cd ClawBio`, terminal font large, `/tmp/demo` deleted, `report_clinician.pdf` open
in a second window. Nothing else on screen.

---

## The line to open with

> "Every polygenic score in this library is labelled *Reference population: EUR*. That's honest.
> And every one of them still hands back a confident percentile for anybody's genome.
> We measured the gap between disclosure and abstention."

---

## Beat 1 — Show the problem (25 s)

```bash
uv run python clawbio.py run prs --demo
```

Six percentiles come back. Point at two of them:

> "Prostate cancer, 100th percentile. Breast cancer, 94th. **Same person.** At least one of
> those is meaningless, and the tool said neither."

Then:

```bash
grep -c "Reference population.*EUR" output/*/prs_report.md
```

> "Six of six, EUR."

## Beat 2 — The mechanism (35 s)

```bash
uv run python -c "
import json,glob,os
cur=json.load(open('skills/gwas-prs/curated_scores.json'))
for f in sorted(glob.glob('skills/gwas-prs/data/*.txt')):
    pid=os.path.basename(f).split('_')[0]; exp=0
    for l in open(f):
        if l.startswith(('#','rsID')): continue
        p=l.split('\t')
        if len(p)>6 and p[6].strip(): exp+=2*float(p[6])*float(p[5])
    print(f'{pid}  curated mean {cur[pid][\"reference_distribution\"][\"mean\"]:.2f}   sum 2*AF*w {exp:.4f}')
"
```

> "The reference mean isn't *derived from* Europeans. It **is** a European allele-frequency sum:
> two times frequency times weight, added up. That's the centre of every percentile on the
> previous screen. Which means the error decomposes — SNP by SNP."

**This is the beat that wins the round. Do not rush it.**

## Beat 3 — Run the gate (40 s)

```bash
uv run python skills/prs-abstain/prs_abstain.py --demo --output /tmp/demo
```

```
EUR_001         REPORT                   distance 1.74
AFR_001         REFUSE_DISTANT           distance 10.38
DEMO_PATIENT    REFUSE_UNDETERMINABLE    distance n/a
```

> "Three people. One released. Two refused, for two different reasons.
> The threshold is 3.47 — mean plus three SD of the European spread. Zero false refusals on
> 22 Europeans, 28 of 28 correct on the rest. Nobody sits between 3.14 and 7.38."

## Beat 4 — The deliberate failure (35 s)

```bash
awk '/^### DEMO_PATIENT/,0' /tmp/demo/report_technical.md | head -8
```

> "This one is the point. Zero markers shared with the reference panel. We never computed a
> distance, because we were never entitled to. The score still returns six confident numbers —
> we withhold all six and say why. And we say the thing tools never say: **a withheld
> percentile is not evidence of low risk.**"

## Beat 5 — SNP by SNP (30 s)

```bash
grep "^PGS000013" /tmp/demo/tables/af_shift_per_variant.csv | head -3
grep "PGS000013" /tmp/demo/report_technical.md | head -2
```

> "T2D score, eight variants. Re-centred on African frequencies it shifts 1.28 standard
> deviations — a median person reported at the 90th percentile. **rs7903146 alone is 92% of
> that shift.** And its neighbour rs12255372, 50 kb away, is the same TCF7L2 locus. Two of the
> eight variants, 48% of the score, one signal. Effective n isn't 6, it's 3.5."

## Beat 6 — Land it (15 s)

Switch to the PDF.

> "Two reports out of one run. Clinician version has no jargon — a test in CI fails the build if
> the word 'centroid' appears in it. Technical version carries the mechanism and twelve named
> limitations, including the one we can't fix: we approximate LD with physical distance, so our
> effective-n is still an overestimate. That's in the report, not hidden."

**Closing line:**

> "Disclosure says *this score was built in Europeans*. Abstention says *so I won't give you a
> number*. 59 tests, a pull request is open. Thank you."

---

## Answers to the questions you will get

**"How did you pick 3.47?"**
Mean plus three SD of within-European distance. It sits in an empty gap: Europeans reach 3.14,
the nearest non-European is 7.38. It's a policy choice inside a measured gap, and the report
says so. Push k to 20 and the tool prints THRESHOLD OVERREACH and voids its own REPORT verdicts.

**"Isn't your demo data synthetic?"**
The scores are real PGS Catalog files with real PMIDs. The allele-frequency table is synthetic
and labelled synthetic inside the file — real arithmetic on illustrative numbers. Fabricating
gnomAD frequencies would have broken the repo's own rule against hallucinated parameters.

**"Doesn't ancestry-risk-profiler already do this?"**
It infers ancestry and abstains below 30 markers. We reuse that exact bound rather than
inventing a second one. It doesn't gate a percentile, doesn't check trait applicability, and
doesn't decompose the error per SNP.

**"What about admixed people?"**
Our principal unhandled case, and it's limitation three in the technical report. One centroid
and one radius can't express partial membership. The demo panel has no admixed samples to
calibrate against, so we name it instead of pretending.

**"Why not just show the raw score?"**
We do. Every branch keeps the raw score. Summing weights is valid for anyone; converting that
sum to a percentile needs a reference distribution that fits the person. We withhold the
interpretation, not the arithmetic.

---

## If something breaks

Reports are already generated at `/mnt/user-data/outputs/`. Open the PDF and walk Beats 4–6
from it. The numbers on the slide are the numbers in the file.
