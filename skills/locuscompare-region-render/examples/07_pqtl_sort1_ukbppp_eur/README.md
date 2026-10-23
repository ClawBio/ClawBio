# 08: pQTL × GWAS: SORT1 plasma sortilin in UKB-PPP EUR × cholesterol-VLDL

**v1.3 multi-modality demo, family 4 of 4**: the canonical pQTL render
the v1.3 work was built around. Shows that the locuscompare orchestrator
now dispatches by OT studyId prefix: `UKB_PPP_*` routes to the new
UKB-PPP fetcher (added in v1.3), everything else stays on eQTL Catalogue.

This is the **protein-side companion to demo `02_*`** (which renders
the SORT1 ge-eQTL in GTEx minor salivary gland against the same VLDL
GWAS). Together, demos `02_*` + `08_*` form the canonical SORT1
dual-modality pair (eQTL + pQTL × LDL-C) that anchors v1.3.

## What this demo proves

> The OT studyId prefix `UKB_PPP_*` dispatches correctly to the
> UKB-PPP fetcher; the bundled regional slice for the canonical
> SORT1 1p13.3 window means downstream ClawBio + K-Dense users
> render this demo with NO Synapse account; auth is opt-in
> per skill design (see SKILL.md "First-time setup").

## OT coloc row resolved

| Field | Value |
|-------|-------|
| OT studyId (exposure) | `UKB_PPP_EUR_SORT1` |
| OT studyId (outcome) | `GCST90269602` |
| Pattern | `UKB_PPP_<ancestry>_<protein>` (the canonical OT pQTL studyId form) |
| Mapping row | `study_id_mappings.yaml` (caller-supplied, see `load_study_id_mappings`) |

## Upstream studies

- **Exposure**: UKB-PPP r1 (Sun 2023 *Nature*, PMID 37794186) plasma
  cis-pQTL for SORT1 (UniProt Q99523, Olink reagent OID20213 on the
  Cardiometabolic panel) in the European-ancestry discovery cohort
  (published n = 34,557; per-variant N in this slice 33,628).
  Verified live 2026-05-15: rs12740374 at chr1:109,274,968 shows
  β=0.120, SE=0.0089, p=4.95e-41, MAF=0.221: the textbook plasma-
  sortilin pQTL signal, protein-side counterpart to the Musunuru 2010
  *Nature* LDL/CHD pairing.
- **Outcome**: GWAS Catalog harmonised `GCST90269602`
  (cholesterol in medium VLDL).
- **LD**: 1000G Phase 3 GRCh38, super-pop EUR.

## How to run

**No Synapse PAT needed** (bundled slice ships with the skill):

```bash
python skills/locuscompare-region-render/cli.py \
    --input skills/locuscompare-region-render/examples/07_pqtl_sort1_ukbppp_eur/config.yaml \
    --output runs/pqtl_sort1_ukbppp/
```

The bundled slice (`skills/ukb-ppp-region-fetch/bundled_slices/
SORT1__EUR__chr1__108774968_109774968.json`, 5,205 variants, 3.5 MB) is
loaded directly with no network round-trip on the exposure side. The
outcome side (GWAS Catalog harmonised) and LD (1000G) still touch the
EBI tabix-on-FTP service.

## Caveats

- This demo's bundled slice is CC-BY 4.0 redistribution under the
  cited Sun 2023 attribution. The slice carries the attribution
  string in its release manifest block. Downstream redistribution
  is also permitted with the same attribution.
- The pQTL signal direction (β=+0.12 with the T allele) means the
  T allele increases plasma sortilin abundance. The same T allele
  DECREASES hepatic SORT1 mRNA in eQTL data (β<0 in the GTEx liver
  ge-eQTL). The mRNA-vs-protein discordance is a biological feature
  of SORT1, not an analysis artifact. Musunuru 2010 documented the
  same direction-flip between hepatic and plasma compartments.
- For arbitrary `(protein, ancestry, region)` queries outside the
  bundled cohort, a free Synapse PAT is required. The skill prints
  a multi-line walkthrough (PAT_REQUIRED_MESSAGE) on first attempt
  if SYNAPSE_AUTH_TOKEN is unset.
