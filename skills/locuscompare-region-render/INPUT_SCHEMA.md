# Canonical sumstats-slice TSV schema

The canonical input format for `locuscompare` — and the contract any data-fetcher
skill must emit if it wants to plug into `locuscompare` (or `fine-mapping`, or
any downstream coloc-analysis skill that consumes regional sumstats).

## Why a canonical schema

`locuscompare` is plot + harmonisation only; it does not know how to fetch
from arbitrary sources. Bundled convenience fetchers cover eQTL Catalogue and
GWAS Catalog harmonised. For everything else (Pan-UKBB, FinnGen direct,
UKB-PPP pQTL, GTEx v10 direct, custom cohorts), the user supplies a
pre-fetched harmonised TSV in this format.

Fetcher skills (separate ClawBio skills, one per source) emit this format too,
so they can chain into `locuscompare` without per-source plumbing.

## File format

- **Encoding:** UTF-8.
- **Line endings:** Unix (`\n`).
- **Compression:** `bgzip`-compressed (`.tsv.bgz` or `.tsv.gz`) or uncompressed (`.tsv`).
  bgzip is preferred for slices > 10 MB.
- **Field separator:** literal TAB (`\t`).
- **Header:** first line is the column header. Required.
- **Sort order:** ascending by `chromosome` then `position_bp`. Required if a tabix index ships alongside; optional otherwise.
- **Tabix index:** if compressed with bgzip and sorted, a `.tbi` sibling is permitted but not required.

## Required columns

These must be present in every row. Missing values: empty string OR explicit `NA` (`NaN` reads the same). A row with a missing required field is invalid and the skill will reject the file.

| Column | Type | Description |
|---|---|---|
| `variant_id` | string | Canonical id in `chrom_pos_ref_alt` form, GRCh38, no `chr` prefix (e.g. `1_109274968_G_T`). Used as the join key between exposure and outcome slices. |
| `chromosome` | string | Chromosome name without the `chr` prefix (`1`, `2`, …, `22`, `X`, `Y`, `MT`). |
| `position_bp` | integer | 1-based GRCh38 position of the variant. |
| `allele_a` | string | The non-effect allele (the reference allele in most sources). Single-base for SNVs; multi-base for indels. |
| `allele_b` | string | The effect allele. `beta` is interpreted as the change in trait per copy of `allele_b`. |
| `beta` | float | Effect size estimate (regression coefficient). |
| `se` | float | Standard error of `beta`. Must be > 0. |
| `p` | float | Two-sided p-value for the variant–trait association. Range (0, 1]. |

## Optional columns

Note: the *config-level* `lead.rs_id` field (orchestrator-only; propagates into the manifest and report) is separate from the TSV `rsid` column below. See `SKILL.md` for the `lead.rs_id` config knob; the column below is the per-variant rs-ID inside a sumstats slice.

Present if the source emits them; ignored if absent. Skills SHOULD produce these when the source has them, for downstream filtering / weighting.

| Column | Type | Description |
|---|---|---|
| `maf` | float | Minor-allele frequency in the study cohort. Range [0, 0.5]. |
| `eaf` | float | Effect-allele frequency (frequency of `allele_b`). Range [0, 1]. |
| `n` | integer | Effective sample size for the variant. For binary traits this is typically the case-equivalent N or the harmonic-mean N per study convention. |
| `rsid` | string | dbSNP rs identifier if known. The join is on `variant_id`, not `rsid`: the two alleles of a multi-allelic site can share one rsID. A variant with several rsIDs may appear once per rsID (the eQTL Catalogue writes its files this way); in an exposure file those rows must have identical statistics and are collapsed to one, keeping every rsID, and rows of one `variant_id` that disagree stop the render. Outcome files are not checked for repeats; give one row per `variant_id`. |
| `info` | float | Imputation INFO score, [0, 1]. |
| `molecular_trait_id` | string | For QTL data: the ENSG / ENST / molecular feature id the variant is tested against. Required by eQTL Cat–style datasets where one TSV bundles multiple traits. |
| `study_id` | string | Source-native study identifier (e.g. eQTL Cat `dataset_id`, GWAS Cat `accession`, FinnGen phenotype code). For self-documentation; not used for join. |
| `chromosome_grch37` | string | Optional GRCh37 chromosome (for cross-build provenance). Same value as `chromosome` if the source is already GRCh38. |
| `position_bp_grch37` | integer | Optional GRCh37 position. Same as `position_bp` if the source is already GRCh38. |

## Conventions

### Coordinate system
GRCh38 throughout. Sources that ship in GRCh37 (e.g. legacy GWAS Catalog studies, BBJ PheWeb) must be lifted over before producing this format. Document the liftover tool and version in the file's accompanying `provenance.json` if shipping a fetcher that does liftover; flag any liftover failures (multi-mapping, unmappable) by dropping the variant.

### Effect-allele convention
`allele_b` is the effect allele. `beta` is per-copy of `allele_b`. This matches:
- the Open Targets coloc graph convention,
- the eQTL Catalogue summary-stats convention (`alt` = effect allele),
- the GWAS Catalog harmonised convention (`hm_effect_allele` = effect allele after harmonisation),
- ClawBio `fine-mapping`'s expected input.

If a source's native convention is the opposite (`allele_a` = effect), the fetcher MUST flip `beta` → `-beta` and swap `allele_a` ↔ `allele_b` before emitting this format.

### Palindromic variants (A/T, G/C)
Palindromic variants have ambiguous strand alignment. `locuscompare`'s harmoniser excludes them by default (they're flagged in the output's `n_palindromic_excluded` count, not silently dropped). Producers should NOT pre-filter palindromic variants; let the consumer decide. Including them is informative even if they're excluded from the LD/coloc analysis.

### Multi-allelic loci
The `variant_id` `chrom_pos_ref_alt` convention naturally splits multi-allelic loci into per-allele rows. Producers MUST emit one row per (chrom, pos, ref, alt) tuple; consumers join on `variant_id`. Do not concatenate multiple alts in a single row.

### INDEL representation
Use the VCF normalisation convention: `allele_a` is the longest common prefix-removed reference, `allele_b` is the longest common prefix-removed alt. `position_bp` is the 1-based position of the leftmost base. This matches `bcftools norm` output. Sources that emit pre-normalised TSVs (eQTL Cat, GWAS Cat harmonised) already follow this; producers should `bcftools norm` first if working from raw VCF.

### Missing data
- Missing `beta`/`se`/`p` for a single variant: row is dropped before harmonisation.
- Missing optional columns: tolerated, downstream features that depend on them gracefully degrade (e.g. MAF-weighted plots fall back to unweighted; `n`-aware diagnostics fall back to N-agnostic).
- `p = 0` is replaced with the underflow floor (`5e-324`) before plotting on a `-log10` axis. Producers MAY emit `p = 0` for variants where the reported p is below floating-point precision; consumers handle the substitution.

### One trait per exposure file
The plot draws one molecular trait. A source that bundles several in one TSV (an eQTL Catalogue or GTEx file) must be filtered to one `molecular_trait_id` before it is used as an exposure; an exposure file with more than one stops the render and names them.

### Sign-checked vs unchecked sumstats
Some sources emit summary stats that have NOT been sign-checked against a reference (e.g. raw GTEx v8 cis-eQTL bgz files). Effect-allele direction in those is undefined relative to the reference allele. Producers MUST document sign-checking status in `provenance.json` (see below) and ideally do the sign check themselves before emitting this format.

## Provenance sidecar (optional but encouraged)

Producers SHOULD ship a `<file>.provenance.json` next to the TSV with at least:

```json
{
  "source": "eqtl_catalogue",
  "source_release": "v7+",
  "study_id": "QTD000276",
  "molecular_trait_id": "ENSG00000134243",
  "fetched_at": "2026-05-06T14:32:17Z",
  "fetched_via": "tabix-on-FTP at https://ftp.ebi.ac.uk/pub/databases/spot/eQTL/sumstats/...",
  "region": {"chromosome": "1", "start_bp": 108774968, "end_bp": 109774968},
  "n_variants_emitted": 3127,
  "build": "GRCh38",
  "liftover": null,
  "sign_check": "passed",
  "harmonisation_version": "",
  "effect_allele_convention": "allele_b is effect allele",
  "source_native_format_link": "https://www.ebi.ac.uk/eqtl/Studies/"
}
```

`locuscompare` reads this sidecar when present and folds it into the rendered plot's manifest block. Absence is non-fatal but reduces reproducibility.

## Validation

`locuscompare` validates input files at the start of every run:

1. **Header check:** all required columns present.
2. **Type check:** numeric columns parse as floats/ints; categorical columns are non-empty.
3. **Range check:** `p ∈ (0, 1]`, `se > 0`, `maf ∈ [0, 0.5]`, `eaf ∈ [0, 1]`, `info ∈ [0, 1]`.
4. **Coordinate check:** `chromosome` ∈ canonical set, `position_bp > 0`.
5. **Allele check:** `allele_a` and `allele_b` both non-empty, both alphanumeric, not equal.
6. **Variant-id check:** `variant_id` matches `^[1-9XYM][0-9]?_[0-9]+_[A-Z]+_[A-Z]+$` and is consistent with `chromosome_position_alleleA_alleleB`.
7. **Sort check:** if a `.tbi` sibling is present, file is sorted as required.

Any check failure is reported with the file path, line number, and field name; the run aborts.

## Reference fetcher implementations

These bundled / planned fetcher skills emit this exact format:

| Skill | Source | Status |
|---|---|---|
| `locuscompare/_fetchers/eqtl_catalogue` | EBI eQTL Catalogue | bundled in v1 |
| `locuscompare/_fetchers/gwas_catalog` | NHGRI-EBI GWAS Catalog harmonised | bundled in v1 |
| `pan-ukbb-region` (separate skill) | Pan-UKBB FTP | planned |
| `finngen-region` (separate skill) | FinnGen direct (registration-aware) | planned |
| `ukb-ppp-region` (separate skill) | UKB-PPP pQTL (registration-gated, Yellow license) | planned |
| `gtex-v10-region` (separate skill) | GTEx v10 cis-eQTL bgz | planned |

User-supplied TSVs from any other source are first-class citizens: provide a path that satisfies this schema, point the config at it, and the skill renders.

## Worked harmonisation recipes

`examples/recipes/` ships small scripts (~20 lines each) that convert a non-bundled source's native format into this schema. v1 covers FinnGen, Pan-UKBB, UKB-PPP, GTEx v10. New recipes can be added by community contribution — they're independent of `locuscompare` core and don't require a release.

## Versioning

The schema itself is versioned in the file header via a `# locuscompare-schema-version: 1.0` comment line. v1 is locked once the upstream PR lands; future versions bump the integer for breaking changes (column rename, semantics change), the minor for additive changes (new optional columns).
