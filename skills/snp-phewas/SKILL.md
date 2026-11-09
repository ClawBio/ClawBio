---
name: snp-phewas
description: >-
  Variant PheWAS across three curated biobanks, Pan-UK Biobank, FinnGen R12 and the Million
  Veteran Program: one variant in, beta / SE / p / allele frequency for every outcome each
  biobank tested (nulls included), on one effect allele, every untested or unreadable
  outcome kept as a row with a reason, outcomes harmonized across the biobanks, and a
  replication verdict per shared outcome for every biobank pair and across all three.
license: MIT
metadata:
  version: "0.3.0"
  author: Aviv Madar
  domain: genomics
  tags:
    - phewas
    - pan-ukbb
    - finngen
    - million-veteran-program
    - replication
    - summary-statistics
    - tabix
  data_license: >-
    CC-BY-4.0 (Pan-UK Biobank outcome manifest and demo slices); FinnGen R12 public
    results terms, free download with citation and acknowledgment required (FinnGen
    outcome catalog, outcome crosswalk, demo slices); Million Veteran Program summary
    statistics through the NHGRI-EBI GWAS Catalog, EMBL-EBI terms of use with no
    restrictions beyond the data owner's, citation given (MVP outcome catalog and its
    source snapshots, demo slices)
  inputs:
    - name: variant
      type: string
      format:
        - rsid
        - chr_pos_ref_alt
      description: One variant, an rsID (rs12740374) or chr_pos_ref_alt (1_109274968_G_T, GRCh38 unless --build GRCh37), given as --variant or as the first line of an --input text file.
      required: true
  outputs:
    - name: report
      type: file
      format:
        - md
      description: Plain-language report with scope, strongest associations, replication for every biobank pair and all three, what was read, and data citations
    - name: result
      type: file
      format:
        - json
      description: Summary counts, per-dataset status histogram, replication verdict counts per pair and across three, provenance
    - name: rows
      type: file
      format:
        - tsv
      description: One row per catalog outcome per dataset with status and plain-language label
  dependencies:
    python: ">=3.11"
    packages:
      - requests>=2.28
      - matplotlib>=3.7
      - numpy>=1.24
      - pyyaml>=6.0
  demo_data:
    - path: demo/
      description: Real +/-10 kb slices around rs12740374 of 3 Pan-UK Biobank, 3 FinnGen R12 and 2 Million Veteran Program (GWAS Catalog harmonized) outcome files, re-indexed, plus a trimmed FinnGen API answer
  endpoints:
    cli: python skills/snp-phewas/snp_phewas.py --variant {variant} --output {output_dir}
  openclaw:
    requires:
      bins:
        - python3
    always: false
    emoji: "🧭"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install:
      - kind: pip
        package: requests
    trigger_keywords:
      - snp phewas
      - variant phewas across biobanks
      - pan-ukbb phewas
      - finngen phewas
      - million veteran program phewas
      - mvp phewas
      - phewas replication
      - phenome-wide association for rs
---

# 🧭 SNP PheWAS (Pan-UK Biobank + FinnGen + Million Veteran Program)

You are **SNP PheWAS**, a specialised ClawBio agent for phenome-wide association lookups.
Your role is to read one variant's association with every outcome of three curated biobanks,
put them on one effect allele and one outcome key, and say which associations replicate.

## Trigger

**Fire this skill when the user says any of:**
- "PheWAS for rs12740374 across UK Biobank, FinnGen and the Million Veteran Program"
- "Does rs12740374 replicate in FinnGen?", "Does it replicate in MVP?"
- "every phenotype this variant was tested for", "all outcomes for this SNP, including nulls"
- "Pan-UKBB PheWAS", "FinnGen PheWAS", "MVP PheWAS", "phenome-wide scan of a variant in biobanks"
- "is this association seen in all three biobanks?", "the same PheWAS in African ancestry"

**Do NOT fire when:**
- The user wants reported associations across the GWAS Catalog / Open Targets credible
  sets, eQTLs (GTEx), or BioBank Japan / UKB-TOPMed PheWeb: route to `gwas-lookup`.
- The user wants summary statistics for a REGION (many variants) from one study: route to
  `gwas-catalog-region-fetch` or `eqtl-catalogue-region-fetch`.
- The user gives a gene, a disease or a trait rather than a variant: this skill takes one
  variant; resolve the variant first.
- The user wants colocalization or Mendelian randomization: route to `mr-region-run` /
  `locuscompare-region-render`.

## Why This Exists

- **Without it**: a PheWAS lookup is a list of reported associations (significance-selected,
  from mixed studies), or one biobank's web page at a time. Nulls are invisible, each source
  uses its own effect allele, a FinnGen endpoint, a UK Biobank outcome and an MVP phecode for
  the same disease cannot be lined up, and an outage looks like "no association".
- **With it**: every outcome of each biobank is a row (tested or not, with the reason), on
  the Open Targets alternative allele, with FinnGen endpoints, Pan-UK Biobank outcomes and MVP
  studies on one join key and a verdict per shared outcome for every pair of biobanks and
  across all three. The default three-biobank read took 159 s and 613 s in two runs with an
  empty index cache, and 57 s and 66 s with the indexes cached (measured 2026-10-19 on a Mac
  mini, 16 workers per biobank; the cold spread is the GWAS Catalog's download speed).
- **Why ClawBio**: values are read from the biobanks' own files; nothing is estimated. The
  skill `gwas-lookup` calls the same FinnGen web endpoint for its FinnGen leg. This skill
  differs in five ways: it reads Pan-UK Biobank and the Million Veteran Program (not in
  `gwas-lookup`); it reads FinnGen's files where the web API hides values (the API reports
  values only below p = 0.05); it accounts for every outcome with a status; it harmonizes
  outcomes across the three biobanks; and it calls replication. Use `gwas-lookup` for a fast
  federated overview first and this skill for the full-phenome, null-complete read of three
  biobanks.

## Core Capabilities

1. **Every outcome, nulls included**: one row per outcome in each biobank's published
   catalog (Pan-UK Biobank 7,160 European-ancestry outcomes, FinnGen R12 2,469 endpoints,
   Million Veteran Program 2,070 European-ancestry GWAS Catalog studies), each with a
   twelve-value `status`; rows equal the catalog as SETS (asserted).
2. **One effect allele**: the Open Targets alternative allele; swapped records are
   re-signed; a strand-ambiguous (A/T, C/G) variant in swapped order has its direction
   withheld, never guessed. An rsID that names several alleles at one position reads the one
   with the most GWAS credible sets in Open Targets, then the most frequent, and says why.
3. **Both builds**: Pan-UK Biobank is GRCh37, read at the Ensembl-mapped position; FinnGen
   and the Million Veteran Program are GRCh38.
4. **Ancestry groups**: Pan-UK Biobank's six groups and MVP's four groups plus its
   multi-population meta-analysis; `--ancestry` reads the same group in every biobank that
   has it and names the biobank that does not.
5. **Outcome harmonization**: 622 of 2,469 FinnGen endpoints map onto 589 Pan-UK Biobank
   outcomes (523 by FinnGen's own published mapping, 99 by a single whole ICD-10 block); MVP
   phecode studies join Pan-UK Biobank phecodes on the code, and MVP lab and vital-sign means
   join Pan-UK Biobank measurements through a reviewed name table.
6. **Replication verdicts**: per shared outcome for every pair of biobanks, at
   p < 0.05 / (outcomes both share), and a three-way class at p < 0.05 / (outcomes all three
   share), with "not tested", "not read" and "low confidence" kept apart from a tested null.
   Two datasets that share participants (MVP's meta-analysis and a group it pools) are never
   compared, and the report says which pair and why.
7. **Three figures**: outcome Manhattan per biobank, replication matrix (one strip per pair
   plus all three), paired effects (one row of panels per pair with shared outcomes). The
   Manhattan groups the outcomes read into 20 body-system categories taken from each
   source's own grouping, colors only the outcomes past the significance line (red raises,
   blue lowers, open circle for no direction), log-compresses the y axis above 20 when the
   top hit is past 40, and labels every outcome past the line on a background patch; one
   that fits nowhere is left off and the panel says how many.

## Scope

**One skill, one task.** This skill looks up one variant across the outcomes of curated
biobanks and nothing else. It does not fine-map, colocalize, run MR, or search literature.

## Input Formats

| Format | Example | Notes |
|--------|---------|-------|
| rsID | `rs12740374` | Resolved through the Open Targets variant search. When the rsID names several alleles at one position, the one with the most GWAS credible sets is read, then the highest gnomAD frequency, then the search's order; the report names the others and what separated the chosen allele from the runner-up |
| GRCh38 variant id | `1_109274968_G_T` | `_`, `:` or `-` separators; `chr` prefix allowed |
| GRCh37 variant id | `1:109817590:G:T` with `--build GRCh37` | Mapped to GRCh38 with the Ensembl assembly mapper |
| Text file | `--input variant.txt` | First non-comment line; one variant per run |

## Workflow

1. **Resolve** (prescriptive): rsID to GRCh38 `chr_pos_ref_alt` with Open Targets; the
   effect allele is that id's alternative allele. Several alleles at one position: rank by
   GWAS credible sets, then gnomAD frequency, then search order. If the weight query fails,
   the search's first id is read and the report says the weights could not be read.
2. **Place** (prescriptive): map to GRCh37 with the Ensembl assembly mapper for Pan-UK
   Biobank. Zero or several placements, a length change, or an indel on the opposite
   strand is `build_unmappable`; an Ensembl outage is `dataset_unreachable`.
3. **Scope** (prescriptive): with two or more biobanks, each is read over its outcomes with
   a readable counterpart in another biobank of the lookup (the same rule in each);
   `--all-outcomes` reads everything. The datasets are read at once.
4. **Read Pan-UK Biobank** (prescriptive): for each outcome in scope, download its tabix
   index once into the cache and read the one position with the skill's own index walker.
5. **Read FinnGen** (prescriptive): one web-API request for all endpoints, then a file read
   for every endpoint the API censored (route `hybrid`).
6. **Read the Million Veteran Program** (prescriptive): a study the GWAS Catalog has not
   harmonized (or holds without summary statistics) is a `source_file_not_readable` row
   from the catalog, with no network call; every other study is read from its harmonized
   file with the same index walker. An X / Y / MT variant is absent by rule (MVP released
   autosomal results only).
7. **Orient and account** (prescriptive): re-express each record on the effect allele;
   assert that each biobank's rows equal its catalog.
8. **Harmonize and compare**: attach the outcome crosswalk; compute replication verdicts
   per pair and across three.
9. **Report**: write `report.md` (scope first), `result.json`, `tables/`, `figures/`,
   `reproducibility/`, and print the report. Interpretation of the biology is the agent's
   job, not the skill's.

## CLI Reference

```bash
# Default: the outcomes matched across the three biobanks (1,156 Pan-UKBB + 622 FinnGen + 727 MVP)
# (first run downloads 4.53 GB of index files into the cache, summed from Content-Length)
python skills/snp-phewas/snp_phewas.py --variant rs12740374 --output /tmp/snp_phewas

# Every outcome in every biobank (7,160 + 2,469 + the 1,221 readable MVP studies; 20.0 GB of indexes)
python skills/snp-phewas/snp_phewas.py --variant rs12740374 --all-outcomes --output /tmp/x

# One ancestry group in every biobank that has it (AFR, AMR, CSA, EAS, EUR, MID, META)
python skills/snp-phewas/snp_phewas.py --variant rs12740374 --ancestry AFR --output /tmp/x

# Quick look: FinnGen web API only (one request; values only where p < 0.05)
python skills/snp-phewas/snp_phewas.py --variant rs12740374 --quick --output /tmp/x

# Other options
#   --dataset pan_ukbb_eur,finngen_r12,mvp_eur | all | pan_ukbb_afr,...   (see --list-datasets)
#   --build GRCh37        --outcomes pan_ukbb_eur:icd10-E78-both_sexes,mvp_eur:GCST90475929
#   --finngen-route hybrid|api|tabix   --workers 16   --index-cache DIR
python skills/snp-phewas/snp_phewas.py --list-datasets

# Via the ClawBio runner
python clawbio.py run snp-phewas --variant rs12740374 --output /tmp/x
python clawbio.py run snp-phewas --demo
```

```python
import sys; sys.path.insert(0, "skills/snp-phewas")
from snp_phewas_lookup import lookup_variant
from snp_phewas_replication import build_replication, replication_lines
result = lookup_variant("rs12740374")      # rows, uncatalogued, provenance, notes
table = build_replication(result)          # verdicts per harmonized outcome per pair, and three-way
print("\n".join(replication_lines(table)))
```

## Demo

```bash
python clawbio.py run snp-phewas --demo
```

Offline, about 2 seconds (1.98 s wall time through the runner, measured 2026-10-19 on a
laptop): rs12740374 over real +/-10 kb slices of 3 Pan-UK Biobank outcome files (LDL
direct, ICD-10 E78, phecode 411), 3 FinnGen endpoints (I9_CHD, E4_LIPOPROT, I9_IHD) and 2
Million Veteran Program studies (GCST90475929 ischemic heart disease, an odds-ratio file;
GCST90475416 LDL cholesterol mean, whose p is 0 in the file). Ischemic heart disease
replicates in all three biobanks; E78 replicates in Pan-UK Biobank and FinnGen, LDL in
Pan-UK Biobank and MVP. Every other catalog outcome is a `not_queried` row (or, for an MVP
study with no readable file, a `source_file_not_readable` row), and the report says so at
the top. The rsID and GRCh37 position use answers recorded from Open Targets and Ensembl
(`snp_phewas_demo.py`). The slices are real source rows, not synthetic (provenance and
acknowledgments in `demo/README.md`).

## Algorithm / Methodology

**Datasets** (`datasets.yaml`, `--list-datasets`):

| dataset_id | outcomes | build | access |
|---|---|---|---|
| `pan_ukbb_eur` (default) | 7,160 | GRCh37 | per-outcome tabix files on AWS Open Data |
| `pan_ukbb_afr` / `_amr` / `_csa` / `_eas` / `_mid` | 2,445 / 982 / 2,739 / 1,533 / 1,248 | GRCh37 | same files |
| `finngen_r12` (default) | 2,469 | GRCh38 | PheWeb variant API + per-endpoint tabix files |
| `mvp_eur` (default) | 2,070 (1,221 readable) | GRCh38, autosomes | the GWAS Catalog's harmonized per-study tabix files |
| `mvp_afr` / `_amr` / `_eas` / `_meta` | 1,757 / 1,479 / 500 / 1,731 | same | same |

Pan-UK Biobank counts are the 2024-06-26 manifest, FinnGen's the R12 manifest (md5 pinned).
**Million Veteran Program** (Verma A et al. Science 2024;385:eadj1182) is served by the GWAS
Catalog as one study per (trait, population): 7,537 studies on the 2026-10-19 snapshot. The
catalog (`snp_phewas_mvp_catalog.py`, `catalogs/mvp_outcomes.tsv.gz`) is built offline from
three md5-pinned snapshots in `catalogs/mvp_source/` and lists EVERY study with its access
status: 4,442 harmonized (sorted, indexed, readable), 1,580 held only as an unsorted file
the Catalog has not harmonized, 1,515 deposited without summary statistics. A study with no
readable file is a row with status `source_file_not_readable`, never read and never a null.
BioBank Japan is listed as not included, with the reason.

**Status** (every catalog outcome is exactly one row): `ok`, `low_confidence` (tested;
flagged by Pan-UK Biobank for too few carriers), `not_tested_in_outcome`,
`absent_in_dataset`, `allele_mismatch`, `build_unmappable` (permanent for this variant);
`dataset_unreachable` (an outage, never an absence), `source_file_missing`,
`record_unparsable`, `above_api_reporting_threshold` (FinnGen API: tested, p >= 0.05),
`not_queried` (outside this lookup), `source_file_not_readable` (an MVP study the GWAS
Catalog holds with no position-indexed file; outside every significance denominator) (not
yet: a re-run or a later release can change them); and `outcome_not_in_catalog`, kept apart
(the FinnGen API returns one endpoint its R12 manifest lacks, `K11_GINGIVITIS_PERIODONTAL`).

**Effect allele**: all sources report reference / alternative alleles on the forward strand
of their build, so allele ORDER decides orientation. Same order: kept. Swapped: beta
negated, EAF 1 - EAF, `allele_flip_applied`. Swapped and strand-ambiguous: direction
withheld (`direction_withheld_reason`), p and SE kept. Other alleles: `allele_mismatch`.
EAF is `af_<ancestry>` (Pan-UK Biobank continuous), the controls' `af_controls_<ancestry>`
(Pan-UK Biobank case/control), `af_alt` (FinnGen) or `effect_allele_frequency` (MVP).

**Reading the GWAS Catalog harmonized files** (MVP): lines end in CRLF (stripped from the
header and each record). The file is reference-anchored with ref = `other_allele`, alt =
`effect_allele`. Case/control files carry an odds ratio with `standard_error = NA`: beta =
ln OR, and SE is the one that reproduces the file's own p, because the odds ratio and its
interval are rounded to about three significant digits (on the bundled ischemic heart
disease slice, an interval-width SE missed the file's p on 85 of 150 records by more than
1e-3 in p or 0.05 in log10 p); the interval is used only where p cannot give an SE (p = 0 or
1, a null effect). `p_value` underflows to 0 at strong hits: -log10 p is then derived from
beta / SE in log space (within 2e-11 of scipy over |z| 0.1 to 300) and `status_detail`
says so. A record the Catalog's harmonizer could not orient (`hm_code` 9 or 14-18) keeps
its p and SE; beta and EAF are withheld (`source_not_orientated`). In the two MVP files the
demo slices come from (GCST90475416 and GCST90475929, both harmonized with palindromic mode
`forward`), every record carries code 5, 6 or 10 to 13 (all 19,681,119 and 19,669,670
records read, 2026-10-24); the 38,609 and 38,598 sites their own logs, bundled beside the
demo slices, count as failed (codes 15 and 16) are not in the files. The withhold rule
covers the codes these files do not exercise, such as 9 (palindromic mode `drop`); a code
outside the Catalog's code table
(https://www.ebi.ac.uk/gwas/docs/methods/summary-statistics), or a harmonized file with no
`hm_code` column, is a `record_unparsable` row.

**Harmonization** (`catalogs/biobank_outcome_crosswalk.tsv.gz`, one row per outcome of each
biobank, mapped or not; 2,469 + 7,223 + 7,537 rows): join key `ICD10:<block>` for a Pan-UK
Biobank ICD-10 outcome, `PHECODE:<code>` for a Pan-UK Biobank both-sexes phecode outcome and
the MVP studies on the same code, `PANUKBB:<outcome_id>` for another Pan-UK Biobank outcome
FinnGen or an MVP measurement matched; a FinnGen endpoint takes its match's key. Bases:
`finngen_published_match` (FinnGen's R12 mapping files), `icd10_exact` (an ICD-10-only
definition naming one block whole), `icd10_code`, `phecode_code`, `phecode_same_code`,
`curated_name_match` (the reviewed MVP measurement table,
`catalogs/mvp_measurement_matches.tsv`), `summary_statistic_variant` (an MVP maximum /
minimum: keeps the id, never compared), `icd10_broad_only` (only a broader or narrower
counterpart: NOT compared), `measurement_near_equivalent`, `sex_stratum_only`, `unmapped`.
When several outcomes of one dataset share a key (25 of the 589 Pan-UK Biobank x FinnGen
keys), one representative is compared (one the source tested first, then strongest basis,
then most cases); the others keep their rows.

**Replication**: per pair, `n_shared` = harmonized outcomes BOTH catalogs carry with a
readable file, whatever the lookup could read (an outage cannot loosen the threshold; an
unreadable MVP study is known from the catalog before any read). Default three biobanks:
Pan-UK Biobank x FinnGen 589, Pan-UK Biobank x MVP 724, FinnGen x MVP 157. Replicates = both
p < 0.05 / n_shared, same sign. Other verdicts: same direction significant in one only;
opposite directions significant in both; significant in one only, direction not shared;
significant in both with a direction withheld; tested in both, significant in neither; and
three no-comparison verdicts (not tested in one, not read in one, low confidence in one). An
outcome one dataset of a pair does not carry is `no_counterpart` for that pair: drawn blank,
not counted. With three datasets, outcomes all three carry (157 by default) get a
`three_way` class at 0.05 / n_shared_3: replicates in all three, in two (the pair named),
discordant, a direction unavailable, significant in one only, null in all three, or the
worst missing cell's no-comparison class. Fewer outcomes are shared by all three than by a
pair, so that line is looser than the pairs'; the figure says so. 5e-8 is drawn on the
figures for reference and decides nothing.

**Sample overlap**: each dataset in `datasets.yaml` declares the participant groups it is
computed from (`samples:`; `mvp_meta` lists the four MVP groups). Two datasets whose groups
intersect share participants, and agreement between them is not independent replication
(the bias participant overlap causes is documented for two-sample Mendelian randomization:
Burgess S, Davies NM, Thompson SG. Genet Epidemiol. 2016;40(7):597-608.
doi:10.1002/gepi.21998. PMID 27625185). Such a pair is
never compared: it has no `n_shared`, threshold, verdicts, counts, strip or panel; it is
listed in `result.json` under `replication.excluded_pairs` with the reason `not_independent`,
and the report names it and the shared participants; the matrix title names it. With three
datasets, any such pair withholds the three-way class. No dataset's rows are dropped.
When an MVP group and the MVP meta-analysis are both read with a third biobank, each can
replicate against it; those two verdicts are not independent of each other (the group's
participants are in both), so they are one piece of evidence, not two.
Overlap is declared rather than derived from ancestry labels so that a meta-analysis pooling
two biobanks (FinnGen publishes FinnGen + UK Biobank meta-analyses) is caught against either.
That UK Biobank, FinnGen and MVP share no participants rests on study design, not
measurement; `datasets.yaml` marks it `assumed_by_design (unverified)`. A pair from one
biobank and two ancestry groups (separate participants, shared outcome definitions) is
compared and headed "same biobank, different ancestry group".

**Cross-ancestry pairs**: a compared pair of different ancestry groups (FinnGen's Finnish
participants count as European), or one with a multi-population meta-analysis, carries
"Effect sizes are not compared across ancestry groups (allele frequency and linkage with the
causal variant differ); the verdict uses sign and significance only." The verdicts are
unchanged; that pair's paired-effects panel gives no correlation.

**Exit status**: when no source answered for any outcome (every read failed), the run is
still written, `result.json` has status `error`, and the script exits 3. A null FinnGen API
answer or a variant absent from a dataset is an answer and exits 0.

**Two views**: by default every biobank is read over its outcomes with a READABLE
counterpart in another biobank of the same lookup, the ones the replication view compares,
one rule for all three (never more of one because its data is quicker to get: FinnGen's
one API request returns every endpoint, and the ones outside the set are not kept). With
the default three that is 1,156 Pan-UK Biobank, 622 FinnGen and 727 MVP outcomes; with
Pan-UK Biobank and FinnGen alone, 589 and 622. The others are `not_queried` rows, and the
report states this first, per biobank, counted from the rows. `--all-outcomes` reads every
outcome in every biobank. A lookup of one biobank alone reads its whole catalog.

**Ancestry**: `--ancestry AFR` reads Pan-UK Biobank AFR and MVP AFR and says that FinnGen
(one population) has no such group; `--ancestry META` reads MVP's multi-population
meta-analysis beside the other biobanks' default groups and says so. A default lookup ends
its report with the other groups and an example command.

**Cost and cache**: each outcome file needs its tabix index, downloaded once into
`--index-cache` (default `~/.clawbio/snp_phewas/index`, or `SNP_PHEWAS_INDEX_CACHE`), keyed
by release. Index sizes, summed from the sources' Content-Length (2026-10-18 for Pan-UK
Biobank and FinnGen in full, 2026-10-19 for the rest): 4.53 GB for the default matched
files (2.15 GB Pan-UK Biobank, 1.10 GB FinnGen, 1.27 GB MVP), 20.0 GB for `--all-outcomes`
(13.54 + 4.37 + 2.12 GB). The CLI says what it will read and download on stderr before it
starts. The `--all-outcomes` wall time is not measured.

## Example Queries

- "PheWAS for rs12740374 across UK Biobank, FinnGen and the Million Veteran Program"
- "Which diseases is rs7412 associated with in FinnGen, Pan-UKBB and MVP, and do they replicate?"
- "Show every Pan-UKBB outcome rs12740374 was tested for, including nulls"
- "The same PheWAS for rs7903146 in African ancestry"
- "Quick FinnGen PheWAS for 19_44908822_C_T"

## Example Output

From `--demo` (abridged; the full file is `examples/expected_output.md`):

```markdown
## Summary
- **Million Veteran Program, European genetic ancestry**: tested in 2 of 2 outcomes read
  (2,070 in the catalog, 849 of them with no readable file); 2 below the threshold 0.05 / 2
  (over the 2 outcomes looked up).
- **Replication, all three** (threshold 0.05 / 157 outcomes all three share): 1 replicate
  in all three, 0 in two.

## Replication across biobanks
**Pan-UK Biobank (European ancestry) and Million Veteran Program (European genetic ancestry)**:
724 shared outcomes, threshold 0.05 / 724.
- replicates (significant in both, same direction): 2 (LDL direct, Ischemic Heart Disease)

**All three biobanks** (Pan-UKBB EUR, FinnGen and MVP EUR): 157 outcomes with a readable
file are shared by all three; each is compared at 0.05 / 157.
- replicates in all three (significant in each, same direction): 1 (Ischemic Heart Disease)
```

Live, Pan-UK Biobank and FinnGen only (`--dataset pan_ukbb_eur,finngen_r12`), rs12740374
(measured 2026-10-18 on a laptop, 16 workers, empty index cache): 87 s wall time; 589 Pan-UK
Biobank files read in 40 s; 1 FinnGen API request plus 564 FinnGen files read in 42 s; 1.9
GiB (2.08 GB) of indexes downloaded. Every one of the 589 + 622 outcomes was tested.
Replication over the 589 shared outcomes (threshold 0.05 / 589): 4 replicate (E78
lipoprotein disorders, ischemic heart disease, angina pectoris, abdominal aortic aneurysm),
5 are same-direction significant in one only, 580 are null in both. Live, the default three
biobanks (measured 2026-10-19 on a Mac mini, 16 workers each), rs12740374 with an empty
index cache: 159 s in one run and 613 s in another, the second with the MVP leg at 605 s
(Pan-UK Biobank 92 s, FinnGen 50 s), so the first-run time is mostly the GWAS Catalog's
download speed; rs2228145 with the indexes cached: 66 s (57 s in an earlier warm run). Read
one after another instead of at once, the first cold run's legs took 308 s (104 + 47 + 153 s).

## Output Structure

```
output_directory/
├── report.md                      # Scope, summary, strongest associations, replication, sources
├── result.json                    # Counts, status histograms, verdict counts, provenance
├── tables/
│   ├── snp_phewas_rows.tsv        # One row per catalog outcome per dataset, every status
│   ├── excluded.tsv               # Rows with no value, each with its reason
│   ├── replication.tsv            # Harmonized outcomes, a verdict per pair, three-way class
│   └── uncataloged.tsv            # (optional) source outcomes outside the published catalog
├── figures/
│   ├── manhattan.png              # Outcome Manhattan, one panel per dataset
│   ├── replication_matrix.png     # Signed z per harmonized outcome x dataset, pair strips
│   └── paired_effects.png         # (optional) beta vs beta with 95% intervals, per pair
└── reproducibility/
    ├── commands.sh                # The exact invocation
    ├── environment.yml            # Python and package versions
    └── checksums.sha256           # SHA-256 of every output, relative to output_directory
```

`replication.tsv` and `replication_matrix.png` need two datasets; `--quick` (FinnGen only)
writes neither.

## Dependencies

**Required**: `requests` (HTTP and range reads), `matplotlib` and `numpy` (figures), `pyyaml`
(dataset catalog). No `pysam`, no `tabix` binary: the index walker is pure Python (`zlib`).

## Gotchas

- **You will want to use `pysam` / `tabix` for the Pan-UK Biobank reads. Do not.** An
  exact-position htslib query on these files misses records: measured 2026-10-17, it
  returned nothing for rs12740374 in 32 of 32 randomly sampled European-ancestry outcome
  files, and for 16 of 250 sampled positions in the LDL file. The index offsets point into
  an empty bgzip block; this skill's reader decompresses the stream as one concatenation,
  in which that block contributes nothing. Against a +/-500 kb htslib fetch as truth it
  agreed on all 32 files and on every sampled position (250 + 62 in the LDL file).
- **You will want to call a FinnGen API result with no value "not associated". Do not.**
  The API fills beta / SE / p only where p < 0.05 (measured 2026-10-17 on four variants:
  154 to 249 of 2,470 results carry values, the largest p 0.0499). `--quick` labels the rest
  `above_api_reporting_threshold`; the default route reads them from the files.
- **You will want to report a missing value as a null. Do not.** `dataset_unreachable`,
  `source_file_missing`, `record_unparsable`, `not_queried` and `source_file_not_readable`
  are reads that failed or did not happen; the replication view calls them "not read",
  never "both null".
- **You will want to count an MVP study the GWAS Catalog has not harmonized. Do not.** 849
  of the 2,070 European-ancestry studies have no readable file on the 2026-10-19 snapshot;
  they are rows with their reason and sit outside every significance denominator.
- **You will want to read an MVP odds ratio's SE from its interval. Do not.** The interval
  is rounded; the SE that reproduces the file's own p is used (see Methodology).
- **You will want to say "no association" for an outcome outside the default scope. Do
  not.** Without `--all-outcomes`, 6,004 Pan-UK Biobank, 1,847 FinnGen and 494 readable MVP
  outcomes are not read; say so and offer the full run.
- **You will want to report the first pair alone. Do not.** With three biobanks there are
  three pairs and an all-three class; the report and stdout print them all.
- **You will want to equate an ICD-10 or phecode match with the same disease definition.
  Do not.** FinnGen endpoints add register sources and control exclusions; phecodes exclude
  related codes from controls; MVP phecodes come from VA health records.
  `icd10_broad_only` and `measurement_near_equivalent` pairs are deliberately not compared.
- **You will want to flip a palindromic (A/T, C/G) variant whose alleles are swapped.
  Do not.** The skill withholds its direction; report the p-value only.
- **You will want to read the first allele the rsID search lists. Do not.** The search
  lists rs7903146's rare C>G (no GWAS credible set) before its common C>T (888); reading the
  first made every outcome an allele mismatch.

## Data access

FinnGen asks users to fill its online form (https://www.finngen.fi/en/access_results)
before downloading summary statistics; instructions arrive by e-mail. The R12 files this
skill reads sit in a public bucket, so nothing enforces the form, but users should fill it.
FinnGen R13 (DF13) was released on June 2, 2026; this skill reads R12. Million Veteran
Program summary statistics are read from the GWAS Catalog's public FTP
(https://ftp.ebi.ac.uk/pub/databases/gwas/summary_statistics/); the Catalog's terms
(https://www.ebi.ac.uk/about/terms-of-use/) place no restrictions beyond the data owner's,
and MVP deposited no license text.

## Safety

- **Local-first**: no user data is uploaded. The variant you query (an rsID or a position)
  is sent to Open Targets, Ensembl and FinnGen; see `docs/data-handling.md`.
- **Disclaimer**: every report ends with the ClawBio medical disclaimer.
- **Audit trail**: `reproducibility/` records the command, the environment and a checksum of
  every output; `result.json` records the dataset releases read.
- **No hallucinated science**: every value is read from the source file or API; thresholds
  are stated with their basis.

## Agent Boundary

The agent (LLM) dispatches and explains. The skill (Python) executes. The agent must NOT
invent associations, override the replication threshold, call an unread or unreadable
outcome a null, report one biobank pair as the whole replication view, or present the
default matched scope as the whole of a biobank.

## Integration with Bio Orchestrator

**Trigger conditions**: an rsID or `chr_pos_ref_alt` together with "PheWAS across biobanks",
"every outcome", "Pan-UKBB", "FinnGen", "Million Veteran Program", "replicate in FinnGen" or
"replicate in MVP".

## Chaining Partners

- `gwas-lookup`: run first for a fast federated overview (GWAS Catalog, Open Targets, GTEx,
  three PheWebs); then this skill for the full, null-complete read of three biobanks.
- `mr-region-run` / `locuscompare-region-render`: the variant usually arrives as a
  colocalization lead; this skill shows what else the lead is associated with.
- `ld-1000g-region-compute`: check whether two associations share a signal before reading
  them as pleiotropy.

## Scaling up

Each lookup reads every outcome file live, so cost grows with outcomes x variants. For many
variants, the same rows can be precomputed: stream each outcome file once and keep the
records at a fixed variant list (for example the lead variants of Open Targets GWAS credible
sets), storing values as float64 so they equal the live read exactly. That is a cluster job
in the files' cloud region, not a laptop one: Pan-UK Biobank European-ancestry files total
about 7.4 TB, and a laptop read the bucket at 7.3 MB/s (measured 2026-10-17), about 12 days
per full pass (an estimate from those two numbers). A reader that returns the same records
from such a table plugs into `lookup_variant(readers=...)`, and every downstream step
(orientation, status accounting, harmonization, replication, figures) is unchanged.

## Maintenance

- **Review cadence**: on each Pan-UK Biobank, FinnGen or GWAS Catalog MVP release.
- **Staleness signals**: FinnGen R13 (released June 2, 2026; this skill reads R12); a changed
  manifest md5 (the catalog rebuild script `snp_phewas_build_catalogs.py` fails loudly on a
  mismatch); the GWAS Catalog harmonizing more MVP studies (refresh the snapshots with
  `snp_phewas_mvp_catalog.refresh_snapshot`, re-pin, rebuild with
  `snp_phewas_build_catalogs.py --families mvp`); the FinnGen API changing its reporting
  threshold; Open Targets changing its variant id format.
- **On a release bump**: rebuild the outcome catalogs, the crosswalk, the demo slices and
  the golden test values together; the index cache is keyed by release, so old indexes are
  never reused.
- **Deprecation**: if every biobank ships a full-null PheWAS API, replace the file reads.

## Citations

- Karczewski KJ, Gupta R, Kanai M, Lu W, Tsuo K, Wang Y, et al. Pan-UK Biobank genome-wide
  association analyses enhance discovery and resolution of ancestry-enriched effects. Nat
  Genet. 2025;57(10):2408-2417. doi:10.1038/s41588-025-02335-7. PMID: 40968291. This research
  has been conducted using the UK Biobank Resource.
- Kurki MI, Karjalainen J, Palta P, Sipilä TP, Kristiansson K, Donner KM, et al. FinnGen
  provides genetic insights from a well-phenotyped isolated population. Nature.
  2023;613(7944):508-518. doi:10.1038/s41586-022-05473-8. PMID: 36653562. We want to
  acknowledge the participants and investigators of the FinnGen study.
- Verma A, Huffman JE, Rodriguez A, Conery M, Liu M, Ho YL, et al. Diversity and scale:
  genetic architecture of 2068 traits in the VA Million Veteran Program. Science.
  2024;385(6706):eadj1182. doi:10.1126/science.adj1182. PMID: 39024449. Summary statistics
  via the NHGRI-EBI GWAS Catalog.
- Bastarache L. Using Phecodes for Research with the Electronic Health Record: From PheWAS
  to PheRS. Annu Rev Biomed Data Sci. 2021;4:1-19.
  doi:10.1146/annurev-biodatasci-122320-112352. PMID: 34465180. (The phecode join key)
- Buniello A, Suveges D, Cruz-Castillo C, Llinares MB, Cornu H, Lopez I, et al. Open Targets
  Platform: facilitating therapeutic hypotheses building in drug discovery. Nucleic Acids
  Res. 2025;53(D1):D1467-D1475. doi:10.1093/nar/gkae1128. PMID: 39657122. (rsID resolution
  and allele weights)
- Dyer SC, Austine-Orimoloye O, Azov AG, Barba M, Barnes I, Barrera-Enriquez VP, et al.
  Ensembl 2025. Nucleic Acids Res. 2025;53(D1):D948-D957. doi:10.1093/nar/gkae1071.
  PMID: 39656687. (GRCh38 to GRCh37 assembly mapping)
- Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files.
  Bioinformatics. 2011;27(5):718-719. doi:10.1093/bioinformatics/btq671. PMID: 21208982.
  (The index format the reader walks)
