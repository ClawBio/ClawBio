---
name: snp-phewas
description: >-
  Variant PheWAS across two curated biobanks, Pan-UK Biobank and FinnGen R12: one variant
  in, beta / SE / p / allele frequency for every outcome each biobank tested (nulls
  included), on one effect allele, every untested outcome kept as a row with a reason,
  outcomes harmonized across the biobanks, and a replication verdict per shared outcome.
license: MIT
metadata:
  version: "0.1.0"
  author: Aviv Madar
  domain: genomics
  tags:
    - phewas
    - pan-ukbb
    - finngen
    - replication
    - summary-statistics
    - tabix
  data_license: >-
    CC-BY-4.0 (Pan-UK Biobank outcome manifest and demo slices); FinnGen R12 public
    results terms, free download with citation and acknowledgment required (FinnGen
    outcome catalog, outcome crosswalk, demo slices)
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
      description: Plain-language report with scope, strongest associations, replication summary, what was read, and data citations
    - name: result
      type: file
      format:
        - json
      description: Summary counts, per-dataset status histogram, replication verdict counts, provenance
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
      - pyyaml>=6.0
  demo_data:
    - path: demo/
      description: Real +/-10 kb slices around rs12740374 of 3 Pan-UK Biobank and 3 FinnGen R12 outcome files, re-indexed, plus a trimmed FinnGen API answer
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
      - phewas replication
      - phenome-wide association for rs
---

# 🧭 SNP PheWAS (Pan-UK Biobank + FinnGen)

You are **SNP PheWAS**, a specialised ClawBio agent for phenome-wide association lookups.
Your role is to read one variant's association with every outcome of two curated biobanks,
put both on one effect allele and one outcome key, and say which associations replicate.

## Trigger

**Fire this skill when the user says any of:**
- "PheWAS for rs12740374 across UK Biobank and FinnGen"
- "Does rs12740374 replicate in FinnGen?"
- "every phenotype this variant was tested for", "all outcomes for this SNP, including nulls"
- "Pan-UKBB PheWAS", "FinnGen PheWAS", "phenome-wide scan of a variant in biobanks"
- "is this association seen in both biobanks?", "replicate this hit in Pan-UK Biobank"

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
  uses its own effect allele, a FinnGen endpoint and a UK Biobank outcome for the same
  disease cannot be lined up, and an outage looks like "no association".
- **With it**: every outcome of each biobank is a row (tested or not, with the reason), on
  the Open Targets alternative allele, with FinnGen endpoints and Pan-UK Biobank outcomes on
  one join key and a verdict per shared outcome, in about four minutes on an empty index
  cache (229 s measured for rs12740374; a warm-cache run was not measured).
- **Why ClawBio**: values are read from the biobanks' own files; nothing is estimated. The
  skill `gwas-lookup` calls the same FinnGen web endpoint for its FinnGen leg. This skill
  differs in four ways: it reads Pan-UK Biobank (not in `gwas-lookup`); it reads FinnGen's
  files where the web API hides values (the API reports values only below p = 0.05); it
  accounts for every outcome with a status; and it harmonizes outcomes across the two
  biobanks and calls replication. Use `gwas-lookup` for a fast federated overview first and
  this skill for the full-phenome, null-complete read of two biobanks.

## Core Capabilities

1. **Every outcome, nulls included**: one row per outcome in each biobank's published
   catalog (Pan-UK Biobank 7,160 European-ancestry outcomes, FinnGen R12 2,469 endpoints),
   each with an eleven-value `status`; rows equal the catalog as SETS (asserted).
2. **One effect allele**: the Open Targets alternative allele; swapped records are
   re-signed; a strand-ambiguous (A/T, C/G) variant in swapped order has its direction
   withheld, never guessed.
3. **Both builds**: Pan-UK Biobank is GRCh37, read at the Ensembl-mapped position; FinnGen
   is GRCh38.
4. **Outcome harmonization**: 622 of 2,469 FinnGen endpoints map onto 589 Pan-UK Biobank
   outcomes (523 by FinnGen's own published mapping, 99 by a single whole ICD-10 block).
5. **Replication verdicts**: per shared outcome, at p < 0.05 / (outcomes both biobanks
   share), with "not tested", "not read" and "low confidence" kept apart from a tested null.
6. **Three figures**: outcome Manhattan per biobank, replication matrix, paired effects.

## Scope

**One skill, one task.** This skill looks up one variant across the outcomes of curated
biobanks and nothing else. It does not fine-map, colocalize, run MR, or search literature.

## Input Formats

| Format | Example | Notes |
|--------|---------|-------|
| rsID | `rs12740374` | Resolved through the Open Targets variant search; the top hit is used and other alleles at a multi-allelic site are named |
| GRCh38 variant id | `1_109274968_G_T` | `_`, `:` or `-` separators; `chr` prefix allowed |
| GRCh37 variant id | `1:109817590:G:T` with `--build GRCh37` | Mapped to GRCh38 with the Ensembl assembly mapper |
| Text file | `--input variant.txt` | First non-comment line; one variant per run |

## Workflow

1. **Resolve** (prescriptive): rsID to GRCh38 `chr_pos_ref_alt` with Open Targets; the
   effect allele is that id's alternative allele.
2. **Place** (prescriptive): map to GRCh37 with the Ensembl assembly mapper for Pan-UK
   Biobank. Zero or several placements, a length change, or an indel on the opposite
   strand is `build_unmappable`; an Ensembl outage is `dataset_unreachable`.
3. **Read Pan-UK Biobank** (prescriptive): for each outcome in scope, download its tabix
   index once into the cache and read the one position with the skill's own index walker.
4. **Read FinnGen** (prescriptive): one web-API request for all endpoints, then a file read
   for every endpoint the API censored (route `hybrid`).
5. **Orient and account** (prescriptive): re-express each record on the effect allele;
   assert that each biobank's rows equal its catalog.
6. **Harmonize and compare**: attach the outcome crosswalk; compute replication verdicts.
7. **Report**: write `report.md` (scope first), `result.json`, `tables/`, `figures/`,
   `reproducibility/`. Interpretation of the biology is the agent's job, not the skill's.

## CLI Reference

```bash
# Default: Pan-UKBB EUR over the 589 outcomes with a FinnGen counterpart + all FinnGen
python skills/snp-phewas/snp_phewas.py --variant rs12740374 --output /tmp/snp_phewas

# Every Pan-UKBB outcome (7,160; an estimated 13 GB of Pan-UKBB indexes on a cold cache)
python skills/snp-phewas/snp_phewas.py --variant rs12740374 --all-outcomes --output /tmp/x

# Quick look: FinnGen web API only (one request; values only where p < 0.05)
python skills/snp-phewas/snp_phewas.py --variant rs12740374 --quick --output /tmp/x

# Other options
#   --dataset pan_ukbb_eur,finngen_r12 | all | pan_ukbb_afr,...   (see --list-datasets)
#   --build GRCh37        --outcomes pan_ukbb_eur:icd10-E78-both_sexes,finngen_r12:I9_CHD
#   --finngen-route hybrid|api|tabix   --workers 16   --index-cache DIR
python skills/snp-phewas/snp_phewas.py --list-datasets

# Via the ClawBio runner
python clawbio.py run snp-phewas --variant rs12740374 --output /tmp/x
python clawbio.py run snp-phewas --demo
```

```python
import sys; sys.path.insert(0, "skills/snp-phewas")
from snp_phewas_lookup import lookup_variant
from snp_phewas_replication import build_replication
result = lookup_variant("rs12740374")      # rows, uncatalogued, provenance, notes
table = build_replication(result)          # verdicts per harmonized outcome
```

## Demo

```bash
python clawbio.py run snp-phewas --demo
```

Offline, about 2 seconds: rs12740374 over real +/-10 kb slices of 3 Pan-UK Biobank outcome
files (LDL direct, ICD-10 E78, phecode 411) and 3 FinnGen endpoints (I9_CHD, E4_LIPOPROT,
I9_IHD). Both harmonized pairs replicate. Every other catalog outcome is a
`not_queried` row, and the report says so at the top. The rsID and GRCh37 position use
answers recorded from Open Targets and Ensembl (`snp_phewas_demo.py`). The slices are real
source rows, not synthetic (provenance and acknowledgments in `demo/README.md`).

## Algorithm / Methodology

**Datasets** (`datasets.yaml`, `--list-datasets`): `pan_ukbb_eur` (default; 7,160 outcomes)
and the five other Pan-UK Biobank ancestry groups `pan_ukbb_afr` / `_amr` / `_csa` / `_eas` /
`_mid` (2,445 / 982 / 2,739 / 1,533 / 1,248), all from the 2024-06-26 manifest (md5 pinned);
`finngen_r12` (default; 2,469 endpoints, R12 manifest, md5 pinned). Million Veteran Program
and BioBank Japan are listed as not included, with the reason.

**Status** (every catalog outcome is exactly one row): `ok`, `low_confidence` (tested;
flagged by Pan-UK Biobank for too few carriers), `not_tested_in_outcome`,
`absent_in_dataset`, `allele_mismatch`, `build_unmappable` (permanent for this variant);
`dataset_unreachable` (an outage, never an absence), `source_file_missing`,
`record_unparsable`, `above_api_reporting_threshold` (FinnGen API: tested, p >= 0.05),
`not_queried` (outside this lookup) (not yet: a re-run can change them); and
`outcome_not_in_catalog`, kept apart (the FinnGen API returns one endpoint its R12 manifest
lacks, `K11_GINGIVITIS_PERIODONTAL`).

**Effect allele**: all sources report reference / alternative alleles on the forward strand
of their build, so allele ORDER decides orientation. Same order: kept. Swapped: beta
negated, EAF 1 - EAF, `allele_flip_applied`. Swapped and strand-ambiguous: direction
withheld (`direction_withheld_reason`), p and SE kept. Other alleles: `allele_mismatch`.
EAF is `af_<ancestry>` (Pan-UK Biobank continuous), the controls' `af_controls_<ancestry>`
(Pan-UK Biobank case/control), or `af_alt` (FinnGen).

**Harmonization** (`catalogs/biobank_outcome_crosswalk.tsv.gz`, one row per outcome of each
biobank, mapped or not): join key `ICD10:<block>` for a Pan-UK Biobank ICD-10 outcome,
`PANUKBB:<outcome_id>` for another Pan-UK Biobank outcome FinnGen matched; a FinnGen
endpoint takes its match's key. Bases: `finngen_published_match` (FinnGen's R12 mapping
files), `icd10_exact` (an ICD-10-only definition naming one block whole), `icd10_code`,
`icd10_broad_only` (only a broader or narrower counterpart: NOT compared), `unmapped`.
When several FinnGen endpoints share a key (25 of the 589 keys), one representative is
compared (strongest basis, then most cases); the others keep their rows.

**Replication**: `n_shared` = harmonized outcomes BOTH catalogs carry, whatever the lookup
could read (an outage cannot loosen the threshold). Replicates = both p < 0.05 / n_shared,
same sign. Other verdicts: same direction significant in one only; opposite directions
significant in both; significant in one only, direction not shared; significant in both
with a direction withheld; tested in both, significant in neither; and three no-comparison
verdicts (not tested in one, not read in one, low confidence in one). 5e-8 is drawn on the
figures for reference and decides nothing.

**Default scope**: Pan-UK Biobank is read over the 589 outcomes that have a FinnGen
counterpart (the outcomes the replication view can compare); the other 6,571 are
`not_queried` rows, and the report states this first. `--all-outcomes` reads all 7,160.

**Cost and cache**: each outcome file needs its tabix index (about 1.8 MB each, measured as
2,863 indexes = 5.1 GB for the default scope), downloaded once into `--index-cache`
(default `~/.clawbio/snp_phewas/index`, or `SNP_PHEWAS_INDEX_CACHE`), keyed by release.
The full `--all-outcomes` read adds about 12 GB more of Pan-UK Biobank indexes and takes an
estimated 12 minutes cold (not measured: 43 s x 7,160 / 589 for Pan-UK Biobank,
plus the 175 s FinnGen read).

## Example Queries

- "PheWAS for rs12740374 across UK Biobank and FinnGen"
- "Which diseases is rs7412 associated with in FinnGen and Pan-UKBB, and do they replicate?"
- "Show every Pan-UKBB outcome rs12740374 was tested for, including nulls"
- "Quick FinnGen PheWAS for 19_44908822_C_T"

## Example Output

From `--demo` (abridged; the full file is `examples/expected_output.md`):

```markdown
## What this answer covers
- **Demo data.** This run reads bundled slices of 3 Pan-UK Biobank outcome files and
  3 FinnGen endpoint files around rs12740374, offline; ...

## Strongest associations: Pan-UK Biobank, European ancestry
| Outcome | Outcome id | beta | SE | p | EAF | Status |
|---|---|---|---|---|---|---|
| LDL direct | `biomarkers-30780-both_sexes-irnt` | -0.1345 | 0.002913 | 1.58e-465 | 0.222 | tested |
| E78 Disorders of lipoprotein metabolism and other lipidaemias | `icd10-E78-both_sexes` | -0.1545 | 0.009285 | 3.31e-62 | 0.224 | tested |

## Replication across biobanks (589 harmonized outcomes)
| Harmonized outcome | Pan-UK Biobank p | FinnGen R12 p | Verdict |
|---|---|---|---|
| E78 Disorders of lipoprotein metabolism ... | 3.31e-62 | 3.23e-57 | replicates (significant in both, same direction) |
| Ischemic Heart Disease | 2.14e-28 | 3.68e-24 | replicates (significant in both, same direction) |
```

Live, default scope, rs12740374 (measured 2026-10-18 on a laptop, 16 workers, empty index
cache): 229 s wall time; 589 Pan-UK Biobank files read in 43 s; 1 FinnGen API request plus
2,274 FinnGen files read in 175 s; 2,863 indexes downloaded (5.1 GB). Every one of the 589 +
2,469 outcomes was tested. Replication over the 589 shared outcomes (threshold 0.05 / 589):
4 replicate (E78 lipoprotein disorders, ischemic heart disease, angina pectoris, abdominal
aortic aneurysm), 5 are same-direction significant in one only, 580 are null in both.

## Output Structure

```
output_directory/
├── report.md                      # Scope, summary, strongest associations, replication, sources
├── result.json                    # Counts, status histograms, verdict counts, provenance
├── tables/
│   ├── snp_phewas_rows.tsv        # One row per catalog outcome per dataset, every status
│   ├── excluded.tsv               # Rows with no value, each with its reason
│   ├── replication.tsv            # Harmonized outcomes, a verdict per dataset pair
│   └── uncatalogued.tsv           # (optional) source outcomes outside the published catalog
├── figures/
│   ├── manhattan.png              # Outcome Manhattan, one panel per dataset
│   ├── replication_matrix.png     # Signed z per harmonized outcome x dataset, verdict strip
│   └── paired_effects.png         # (optional) beta vs beta with 95% intervals
└── reproducibility/
    ├── commands.sh                # The exact invocation
    ├── environment.yml            # Python and package versions
    └── checksums.sha256           # SHA-256 of every output, relative to output_directory
```

`replication.tsv` and `replication_matrix.png` need two datasets; `--quick` (FinnGen only)
writes neither.

## Dependencies

**Required**: `requests` (HTTP and range reads), `matplotlib` (figures), `pyyaml` (dataset
catalog). No `pysam`, no `tabix` binary: the index walker is pure Python (`zlib`).

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
  `source_file_missing`, `record_unparsable` and `not_queried` are reads that failed or did
  not happen; the replication view calls them "not read", never "both null".
- **You will want to say "no Pan-UK Biobank association" for an outcome outside the default
  scope. Do not.** Without `--all-outcomes`, 6,571 Pan-UK Biobank outcomes are not read; say
  so and offer the full run.
- **You will want to equate an ICD-10 match with the same disease definition. Do not.**
  FinnGen endpoints add register sources and control exclusions; phecodes exclude related
  codes from controls. `icd10_broad_only` pairs are deliberately not compared.
- **You will want to flip a palindromic (A/T, C/G) variant whose alleles are swapped.
  Do not.** The skill withholds its direction; report the p-value only.

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
invent associations, override the replication threshold, call an unread outcome a null, or
present the default 589-outcome scope as the whole of Pan-UK Biobank.

## Integration with Bio Orchestrator

**Trigger conditions**: an rsID or `chr_pos_ref_alt` together with "PheWAS across biobanks",
"every outcome", "Pan-UKBB", "FinnGen", or "replicate in FinnGen".

## Chaining Partners

- `gwas-lookup`: run first for a fast federated overview (GWAS Catalog, Open Targets, GTEx,
  three PheWebs); then this skill for the full, null-complete read of two biobanks.
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

- **Review cadence**: on each Pan-UK Biobank or FinnGen release.
- **Staleness signals**: a new FinnGen release (R13); a changed manifest md5 (the catalog
  rebuild script `snp_phewas_build_catalogs.py` fails loudly on a mismatch); the FinnGen API
  changing its reporting threshold; Open Targets changing its variant id format.
- **On a release bump**: rebuild the outcome catalogs, the crosswalk, the demo slices and
  the golden test values together; the index cache is keyed by release, so old indexes are
  never reused.
- **Deprecation**: if both biobanks ship a full-null PheWAS API, replace the file reads.

## Citations

- Karczewski KJ, Gupta R, Kanai M, Lu W, Tsuo K, Wang Y, et al. Pan-UK Biobank GWAS
  improves discovery, analysis of genetic architecture, and resolution into
  ancestry-enriched effects. medRxiv [Preprint]. 2024. doi:10.1101/2024.03.13.24303864.
  (Preprint; no PMID. Europe PMC PPR822111.) This research has been conducted using the UK
  Biobank Resource.
- Kurki MI, Karjalainen J, Palta P, Sipilä TP, Kristiansson K, Donner KM, et al. FinnGen
  provides genetic insights from a well-phenotyped isolated population. Nature.
  2023;613(7944):508-518. doi:10.1038/s41586-022-05473-8. PMID: 36653562. We want to
  acknowledge the participants and investigators of the FinnGen study.
- Buniello A, Suveges D, Cruz-Castillo C, Llinares MB, Cornu H, Lopez I, et al. Open Targets
  Platform: facilitating therapeutic hypotheses building in drug discovery. Nucleic Acids
  Res. 2025;53(D1):D1467-D1475. doi:10.1093/nar/gkae1128. PMID: 39657122. (rsID resolution)
- Dyer SC, Austine-Orimoloye O, Azov AG, Barba M, Barnes I, Barrera-Enriquez VP, et al.
  Ensembl 2025. Nucleic Acids Res. 2025;53(D1):D948-D957. doi:10.1093/nar/gkae1071.
  PMID: 39656687. (GRCh38 to GRCh37 assembly mapping)
- Li H. Tabix: fast retrieval of sequence features from generic TAB-delimited files.
  Bioinformatics. 2011;27(5):718-719. doi:10.1093/bioinformatics/btq671. PMID: 21208982.
  (The index format the reader walks)
