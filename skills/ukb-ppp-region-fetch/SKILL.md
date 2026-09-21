---
name: ukb-ppp-region-fetch
description: |
  Fetch a regional slice of plasma pQTL summary statistics from the UK
  Biobank Pharma Proteomics Project (UKB-PPP; Sun 2023 Nature) for a
  specific (protein, ancestry) measurement. Use when an agent needs
  per-variant beta / SE / p-value around a coloc-lead variant for
  downstream colocalisation, Mendelian randomisation, or regional
  plotting against a pQTL exposure. The canonical use case is the
  cis-window around the protein's coding gene TSS, but UKB-PPP releases
  full-genome summary stats per protein so any GRCh38 window (including
  trans loci) is supported when the user supplies an explicit
  (chromosome, start_bp, end_bp). Input: protein_label (HGNC or
  UniProt), ancestry, chromosome, start_bp, end_bp. Output: harmonised
  TSV slice + manifest + human-readable report.
license: MIT
metadata:
  skill-author: Aviv Madar
  version: 0.1.0
  domain: bioinformatics
  tags:
    - pqtl
    - ukb-ppp
    - region-fetch
    - regenie
    - summary-statistics
    - proteomics
    - olink
  inputs:
    - name: protein_label
      type: string
      description: HGNC symbol (e.g. SORT1) or UniProt accession (e.g. Q99523). Resolves to the canonical UKB-PPP per-protein file via the Synapse listing.
      required: true
    - name: ancestry
      type: string
      description: One of EUR (European discovery; N=46,673), AFR (African; N=931), CSA / SAS (Central / South Asian; N=920), EAS (East Asian; N=262), MID (Middle East; N=124), AMR (American Hispanic; N=60), ALL (Combined multi-ancestry meta; N=47,970). Per Sun 2023 Nature Table 1.
      required: true
    - name: chromosome
      type: string
      description: Chromosome without `chr` prefix (1, 2, ..., X). UKB-PPP per-protein archive has one file per chromosome.
      required: true
    - name: start_bp
      type: integer
      description: Region start, 1-based GRCh38.
      required: true
    - name: end_bp
      type: integer
      description: Region end, 1-based GRCh38 (inclusive).
      required: true
  outputs:
    - name: variants
      type: list
      description: Per-variant rows with variant_id (chr_pos_ref_alt, OT convention), chromosome, position, ref, alt, beta, se, p_value (linear; converted from LOG10P), maf, effect_allele_frequency, molecular_trait_id (Olink reagent id), study_id (Synapse fileID).
    - name: release
      type: object
      description: UKBPPPRelease with protein_hgnc, protein_uniprot, olink_reagent_id, olink_panel, ancestry, ancestry_label, n_samples, synapse_id, source_url, release_label, fetched_at_utc.
  dependencies:
    - python>=3.10
    - synapseclient>=3.0
    - requests>=2.28
  env:
    - name: SYNAPSE_AUTH_TOKEN
      description: Free Synapse personal access token. Required for the file-download path; not required for the listing-only smoke test. Obtain via https://www.synapse.org/Profile:settings ("Personal Access Tokens"). No UKB Application required for the summary-stats layer.
      required: true
  demo_data:
    - examples/sort1_ukb_ppp_eur.json
  endpoints:
    - https://repo-prod.prod.sagebase.org/repo/v1/                # Synapse REST (listing, file download)
    - s3://ukbiobank.opendata.sagebase.org                        # AWS Open Data Registry mirror (gated)
  openclaw:
    requires:
      bins:
        - python3
      env:
        - SYNAPSE_AUTH_TOKEN
      config:
    always: false
    emoji: "🧪"
    homepage: https://github.com/ClawBio/ClawBio
    os:
      - darwin
      - linux
    install: |
      pip install synapseclient requests
    trigger_keywords:
      - ukb-ppp pqtl region fetch
      - pqtl regional summary stats
      - sun 2023 ukbppp
      - protein qtl summary stats
      - plasma pqtl regional fetch
---

# 🧪 UKB-PPP Region Fetch

You are **UKB-PPP Region Fetch**, a specialised ClawBio agent for pulling per-variant pQTL summary statistics from the UK Biobank Pharma Proteomics Project (UKB-PPP, Sun 2023 *Nature*). Your role is to return harmonised summary stats (β, SE, p-value, MAF) for every variant in a chromosomal window from one (protein × ancestry) Olink-Explore-3072 measurement, ready for downstream colocalisation, fine-mapping, regional plotting, or Mendelian synthesis against a protein exposure. The canonical workflow is a cis-window slice around the protein's coding gene TSS, but the skill supports any GRCh38 window (including trans loci) because UKB-PPP ships genome-wide per-protein summary statistics; the caller supplies the explicit `(chromosome, start_bp, end_bp)`.

## First-time setup (IMPORTANT)

The skill ships with two fetch paths. **Most users only need the first**:

1. **Bundled-slice path (no auth, no setup).** Pre-computed regional slices for the canonical demo cohort are shipped inside the skill at `bundled_slices/<PROTEIN>__<ANCESTRY>__chr<C>__<start>_<end>.json.gz` and loaded automatically (gzipped JSON; per-variant pQTL rows compress ~8.5x, so a 5,000-variant slice is ~430 KB on disk vs ~3.5 MB raw). v0.1.0 ships the SORT1 / EUR / OID20213 slice (chr1:108,774,968-109,774,968, the 1p13.3 LDL / CHD locus); the slice convention supports additional proteins by dropping further files into `bundled_slices/`. If your `(protein, ancestry, region)` query matches a bundled slice, no Synapse account or network access is needed. Redistribution is permitted under CC-BY 4.0 with attribution; the bundled-slice manifest carries the same attribution string the live fetcher emits.

2. **Live Synapse fetch (free PAT required).** For arbitrary queries beyond the bundled demo cohort, the skill falls through to a live Synapse downloader. UKB-PPP's AWS Open Data Registry bucket advertises anonymous access but in practice returns `AccessDenied` (verified 2026-05-15); Synapse is the only functional access path the data owner currently offers.

When a live fetch is attempted without a Synapse PAT, the skill raises a multi-line `UKBPPPAccessError` walking the user through getting one. Summary of the steps:

   1. Register a free account at <https://www.synapse.org>.
   2. Accept the Synapse Terms of Use (one-time click-through).
   3. Open <https://www.synapse.org/Profile:settings> → "Personal Access Tokens".
   4. Click "Create new token". Tick scopes `view` and `download`. Copy the token immediately (Synapse shows it once).
   5. Export it: `export SYNAPSE_AUTH_TOKEN=<token>` and re-run.

**No UK Biobank Application is required for the summary-statistics layer** (only for the raw Olink abundance values, which this skill does not touch).

## Overview

UKB-PPP (Sun et al. 2023 *Nature*) is the largest open-access plasma proteomic GWAS resource, profiling 2,923 Olink Explore 3072 proteins across 54,219 UK Biobank participants stratified into European discovery (N=46,673) plus six smaller ancestry breakouts (African, Central/South Asian, East Asian, Middle East, American Hispanic, and a Combined multi-ancestry meta-analysis). Summary statistics are released per (protein × ancestry) as REGENIE step-2 outputs and packaged as `<HGNC>_<UniProt>_<OlinkID>_v1_<Panel>.tar` archives on Synapse (`syn51364943`). Each tar contains one gzipped REGENIE file per autosome + X. This skill resolves a protein label (HGNC or UniProt) to the canonical Synapse fileID, downloads the protein's tar to a local cache, extracts the per-chromosome file, filters to a `(chr, start, end)` window, and emits a harmonised TSV slice plus a provenance manifest.

## Trigger

**Fire when** the user (or upstream agent step) wants:

- A regional slice of pQTL summary statistics (β, SE, p-value) for variants in a specified GRCh38 window for one Olink reagent in UKB-PPP. Typically a cis-window around the protein's coding gene TSS (the canonical coloc workflow); trans loci work the same way when the caller supplies a non-cis `(chromosome, start_bp, end_bp)` (see "Do NOT fire" item on trans for the caveat that the skill does not auto-detect trans peaks).
- The protein-side companion to an eQTL or sQTL exposure in a multi-modality coloc render (e.g. SORT1 eQTL liver × pQTL plasma × LDL-C GWAS).
- Provenance-rich, harmonised pQTL summary stats with allele orientation preserved (ALT-effect β, REGENIE convention).

**Do NOT fire when** the user wants:

- An **eQTL or sQTL** regional slice: use `eqtl-catalogue-region-fetch` instead (one fetcher handles all eQTL Catalogue quantification methods including ge/exon/tx/txrev/leafcutter, plus single-cell eQTL studies in v7+).
- A **single-variant pQTL lookup**: UKB-PPP's full archive is large; a per-variant lookup against the Open Targets pQTL coloc table is cheaper for one-point queries.
- The **deCODE pQTL** panel (Ferkingstad 2021): a different upstream cohort with separate access terms; this skill targets UKB-PPP only. Choose the upstream-source skill at the orchestrator level.
- The **raw Olink abundance values** linked to phenotypes: those are gated behind a UK Biobank Application via Synapse `syn52364558` and are out of scope for the public locuscompare render path.
- **trans-pQTL signals** at distant loci: UKB-PPP releases full-genome summary stats per protein, so trans signals are present in the data, but the `(chromosome, start_bp, end_bp)` window must be supplied explicitly; the skill does not auto-detect trans peaks.

## Scope

**One skill, one task.** This skill fetches one `(protein × ancestry)` pair's regional summary statistics from UKB-PPP and writes them as a harmonised TSV plus a provenance manifest. It does NOT iterate proteins, ancestries, or windows; it does NOT do pQTL fine-mapping or coloc directly; it does NOT fetch eQTL / sQTL / sceQTL (use `eqtl-catalogue-region-fetch`); it does NOT fetch the deCODE pQTL panel. The caller composes those workflows on top.

## Workflow

When an agent asks for a regional pQTL slice from UKB-PPP:

1. **Resolve the protein.** The skill lists the requested ancestry folder on Synapse (`syn51365303` for EUR, etc.) and parses `<HGNC>_<UniProt>_<OlinkID>_v1_<Panel>.tar` filenames into a (HGNC, UniProt) -> Synapse fileID index. Lookup tolerates both keys; HGNC is the default surface. The listing call is auth-free; only the subsequent download requires a Synapse PAT.
2. **Resolve the ancestry.** EUR is the default for most coloc renders given its 50x cohort size advantage; use AFR / EAS / CSA / MID / AMR for ancestry-specific renders. The 1000 Genomes super-population for the LD reference panel (1000G Phase 3) maps EUR -> EUR, AFR -> AFR, EAS -> EAS, CSA/SAS -> SAS, MID -> EUR proxy with caption caveat, AMR -> AMR.
3. **Download the tar (cached).** First fetch downloads the protein tar via `synapseclient` to the local cache (`UKB_PPP_CACHE_DIR` env or `~/.clawbio/ukb_ppp_region_fetch_cache/`). Repeat fetches across regions on the same protein reuse the cached tar.
4. **Extract the chromosome.** The tar contains one gzipped REGENIE file per autosome + X. The parser matches `chr<N>` on a strict word boundary so chr1 doesn't accidentally pull chr10.
5. **Stream-filter to the window.** REGENIE files are plain gzip (not BGZ / tabix), so regional fetches scan one chromosome's file linearly. ~1M rows per chromosome; ~2 s on a modern laptop per window. Rows are normalised to OT `chr_pos_ref_alt` ALT-effect convention; LOG10P is converted to a linear p-value; A1FREQ above 0.5 is folded to MAF.
6. **Write outputs** to `--output <dir>/`: a flat `variants.tsv` (effect-allele-aligned, GRCh38, ALT-effect β), a `manifest.yaml` with provenance (study_label, release_label, protein_hgnc, protein_uniprot, olink_reagent_id, olink_panel, ancestry, ancestry_label, n_samples, synapse_id, source_url, fetched-at UTC timestamp, attribution string), and a `report.md` human-readable summary.

## CLI Reference

```bash
# Standard usage with a config file (Synapse PAT in env)
SYNAPSE_AUTH_TOKEN=... python skills/ukb-ppp-region-fetch/ukb_ppp_region_fetch.py \
    --input <config.json> --output <output_dir>

# Bundled demo (SORT1 plasma pQTL in EUR; the canonical 1p13.3 LDL/CHD locus)
SYNAPSE_AUTH_TOKEN=... python skills/ukb-ppp-region-fetch/ukb_ppp_region_fetch.py \
    --demo sort1_ukb_ppp_eur --output /tmp/sort1_ukbppp_demo

# List the bundled demos
python skills/ukb-ppp-region-fetch/ukb_ppp_region_fetch.py --list-demos

# Via ClawBio runner
SYNAPSE_AUTH_TOKEN=... python clawbio.py run ukb-ppp-region-fetch --input <config.json>
```

Config schema (JSON or YAML):

```json
{
  "protein_label": "SORT1",
  "ancestry": "EUR",
  "chromosome": "1",
  "start_bp": 108774968,
  "end_bp": 109774968
}
```

## Example Output

Running `--demo sort1_ukb_ppp_eur` (see `examples/expected_output.md` for the full reproduction):

```
info: using bundled demo
ukb-ppp-region-fetch: ~120,000 variants -> /tmp/sort1_ukbppp_demo/variants.tsv
  source: UKB-PPP | SORT1 (Q99523, OID20213) | European (discovery) (EUR)
```

## Gotchas

1. **Live fetch requires a free Synapse PAT, not anonymous AWS Open Data.** The AWS Open Data Registry page advertises `arn:aws:s3:::ukbiobank.opendata.sagebase.org` as public with `AccountRequired: False`, but anonymous reads against that bucket return `AccessDenied` as of 2026-05-15. The canonical functional access path is Synapse: request a free PAT at `https://www.synapse.org/Profile:settings` and export it as `SYNAPSE_AUTH_TOKEN`. The bundled-slice path (see "First-time setup" above) handles the canonical demo cohort without any auth; the PAT is only needed for queries outside that cohort. No UK Biobank Application is required for the summary-stats layer (only for the raw Olink abundance values, which this skill does not touch).

2. **One protein, one tar, full-genome.** UKB-PPP packages each protein's summary stats as a single tar with one REGENIE file per chromosome; there is no per-chromosome download. First-fetch for a protein downloads ~100–500 MB. Subsequent regional fetches on the same protein reuse the cached tar.

3. **REGENIE LOG10P, not -log10(p).** The REGENIE column reports `|log10(p)|` (always positive). The skill converts to linear `p_value = 10^-LOG10P` at the row boundary. Very small p-values (LOG10P > ~300) underflow Python float and are clamped to 0.0 rather than raising.

4. **A1FREQ is the ALLELE1 (ALT, effect) frequency, not MAF.** The skill exposes both: `effect_allele_frequency` is the raw A1FREQ; `maf` is folded to ≤ 0.5. Downstream code (e.g. palindromic-variant excluder) reads `effect_allele_frequency`.

5. **β is on the ALT allele.** Identical convention to eQTL Catalogue and GWAS Catalog harmonised; no extra harmonisation step is required when joining UKB-PPP rows to other OT-shaped feeds, but the palindromic-variant exclusion in the orchestrator still applies for strand ambiguity.

6. **Some HGNC symbols map to >1 Olink reagent.** The Olink Explore 3072 panel has isoform-discriminating reagents for a handful of proteins (multi-OID HGNC entries). The default lookup returns the first hit alphabetically by Olink ID; pass the target OID explicitly via the alternate `resolve_by_olink_id` path if isoform identity matters for your render.

7. **Per-chromosome file names vary slightly across the release.** The parser matches `chr<N>` on a strict word boundary inside `.tar` members, accepting names like `discovery_chr1_<protein>_*.regenie.gz` or `chr1_*.tsv.gz`. The strict boundary prevents chr1 from accidentally matching chr10 / chr11, a class of bug that would silently return the wrong chromosome's data.

## Safety

**Not for clinical decisions.** This skill returns research-grade summary statistics from a public proteomic GWAS. Do not use the output for direct clinical decision-making, diagnosis, or treatment selection without independent validation by a qualified clinician.

**Effect estimates may not generalise across populations.** UKB-PPP's discovery cohort is overwhelmingly European (N=46,673 vs N=931 for African, the next-largest stratum). Effect sizes from EUR-discovery analyses should not be assumed to apply uniformly across other ancestries; the orchestrator's caption layer flags this when the ancestry side of an LD reference panel mismatches the source study.

**Plasma vs tissue.** UKB-PPP measures circulating plasma proteins, which is biologically distinct from cell- or tissue-level protein abundance. Downstream interpretation should not assume a plasma cis-pQTL implies an identical effect on intra-cellular abundance for the same protein.

## Agent Boundary

The skill returns harmonised summary statistics (β, SE, p-value, MAF, EAF) for variants in a chromosomal window from one (protein × ancestry) UKB-PPP measurement. The agent should:

- **Use the output as input to colocalisation, fine-mapping, or Mendelian randomisation tooling.** These are the appropriate downstream methods for inferring causal effects.
- **NOT make causal-effect claims directly from a single pQTL p-value.** Statistical association ≠ causation; instrumental-variable assumptions must be satisfied for MR.
- **NOT cherry-pick variants by p-value alone.** Statistical inference requires the full credible set / window context.
- **NOT compare effect sizes across ancestries without acknowledging cohort N.** EUR's 50x sample-size advantage means EUR effect estimates have much tighter SEs; "no effect" in a small-N ancestry stratum may reflect power, not biology.
- **Surface protein identity, Olink reagent ID, and ancestry in the user-facing reply** alongside any β / p-value the agent quotes. Per the user-friendly enum-expansion rule (`AGENTS.md`), expand all three fields: `protein = SORT1 (Q99523, OID20213); ancestry = European (discovery) (EUR); N = 46,673`.
- **NOT silently swap Olink reagents.** If the user asked for `SORT1` and the dataset is `SORT1-AOH2` (a different isoform reagent), the agent must say so explicitly.

## Citations

- Sun, B.B., Chiou, J., Traylor, M. et al. (2023). *Plasma proteomic associations with genetics and health in the UK Biobank.* Nature 622, 329–338. doi:10.1038/s41586-023-06592-6 (PMID 37794186).
- UKB-PPP data release: <https://www.synapse.org/Synapse:syn51364943> (Sage Bionetworks, CC-BY 4.0).
- REGENIE: Mbatchou, J. et al. (2021). *Computationally efficient whole-genome regression for quantitative and binary traits.* Nat Genet 53, 1097–1103.
- Olink Explore 3072 panel content: <https://olink.com/products-services/explore/>.
- `synapseclient` Python library: Sage Bionetworks (Apache-2.0). Used by the live-fetch path; not invoked when serving a bundled slice.
