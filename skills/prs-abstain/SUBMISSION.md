# SUBMISSION PACK — Challenge 3, "Whose genome does this fail?"

## Post this in #berlin-demos before 16:10

> **prs-abstain** — https://github.com/<YOUR-FORK>/ClawBio/tree/skill/prs-abstain
> Every PRS in the library discloses "Reference population: EUR" and still returns a confident
> percentile for anyone. We built the skill that withholds it instead, and traces the error to
> the individual SNPs that cause it.

(Replace `<YOUR-FORK>`. One line is what they asked for; the link carries the rest.)

## Open the PR (Route A)

```bash
gh repo fork ClawBio/ClawBio --clone && cd ClawBio
git checkout -b skill/prs-abstain
# copy skills/prs-abstain/ in, plus the cli.py alias block from HOWTO.md
git add skills/prs-abstain/ clawbio/cli.py skills/catalog.json
git commit -m "Add prs-abstain: withhold PRS percentiles that cannot be interpreted"
git push -u origin skill/prs-abstain
gh pr create --title "Add skill: prs-abstain" --body-file skills/prs-abstain/PR_BODY.md
```

No `gh`? Fork on the web, upload the `skills/prs-abstain/` folder, commit to a new branch,
click "Compare & pull request", paste `PR_BODY.md`.

## Against the four "done" criteria

| Requirement | Status |
|---|---|
| Runs live on real data | Yes. Real PGS Catalog scoring files, real PMIDs. The AF table is synthetic and says so in the file. |
| Every claim traces to a source that resolves | Kosoy 2009 verified live: PMID 18683858, Hum Mutat 30(1):69-78, DOI 10.1002/humu.20822. Thresholds computed from supplied data, not asserted. |
| At least one input where it correctly says it cannot answer | Two. DEMO_PATIENT (0 shared markers, undeterminable) and AFR_001 (distance 10.38 > 3.47). |
| Does something not obvious before you built it | Reference mean = Sum(2·AF·w) exactly. Breast + prostate percentiles for one person. Duplicate coordinates in two scoring files. TCF7L2 pair = 48% of the T2D score. |

## Against the three judging criteria

**Originality.** The room will build ancestry plots. We built the refusal, then found the
mechanism underneath it and attributed the error per SNP. The deliverable is the abstention,
which is what the brief asked for and what most teams skip.

**Impact.** Any PRS report generator can adopt this gate. It is the difference between a
disclosure footnote nobody reads and a number that is not printed. It also caught two real
defects in the library's own demo data.

**ClawBio implementation.** A conforming skill, registered in the runner, chaining from
`claw-ancestry-pca` and `gwas-prs` without reformatting, reusing `ancestry-risk-profiler`'s
marker bound rather than inventing a second one. 59 tests, red/green TDD, PR open.

## Files

| File | Use |
|---|---|
| `DEMO.md` | Stage runbook, 3 min, every command tested |
| `PR_BODY.md` | Paste into the pull request |
| `HOWTO.md` | Install, flags, real-data usage |
| `SKILL.md` | The skill definition |
| `report_clinician.pdf` | Backup if the terminal fails |
| `report_technical.pdf` | For the questions afterwards |

## Do not say on stage

Don't claim the AF numbers are gnomAD. Don't claim LD is measured. Don't claim the clean
3.14–7.38 gap generalises. All three are in the reports as limitations, and being the team that
volunteers them is worth more than being the team that gets caught.
