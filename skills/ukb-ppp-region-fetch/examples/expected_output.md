# Expected output: SORT1 ge-pQTL in UKB-PPP European (discovery)

Running:

```bash
# Bundled-slice path (no auth needed; slice ships in this repo)
python skills/ukb-ppp-region-fetch/ukb_ppp_region_fetch.py \
    --demo sort1_ukb_ppp_eur --output /tmp/sort1_ukbppp_demo

# Or live Synapse fetch (PAT required only if you ask for a region the
# bundle doesn't cover)
SYNAPSE_AUTH_TOKEN=eyJ... python skills/ukb-ppp-region-fetch/ukb_ppp_region_fetch.py \
    --demo sort1_ukb_ppp_eur --output /tmp/sort1_ukbppp_demo
```

produces (numbers below verified live 2026-05-15 against UKB-PPP r1):

```
info: using bundled slice SORT1__EUR__chr1__108774968_109774968.json.gz (no Synapse auth needed)
ukb-ppp-region-fetch: 5205 variants -> /tmp/sort1_ukbppp_demo/variants.tsv
  source: UKB-PPP | SORT1 (Q99523, OID20213) | European (discovery) (EUR) | N = 33,628
```

`<output_dir>/manifest.yaml`:

```yaml
skill: ukb-ppp-region-fetch
version: 0.1.2
protein_label: SORT1
ancestry: EUR
region:
  chromosome: '1'
  start_bp: 108774968
  end_bp: 109774968
n_variants: 5205
release:
  study_label: UKB-PPP
  release_label: UKB-PPP r1 2023 (Sun 2023)
  protein_hgnc: SORT1
  protein_uniprot: Q99523
  protein_label: SORT1 (Q99523, OID20213)
  olink_reagent_id: OID20213
  olink_panel: Cardiometabolic
  ancestry: EUR
  ancestry_label: European (discovery)
  n_samples: 33628
  n_samples_min: 33628
  n_samples_label: 33,628
  cohort_n: 34557
  cohort_n_source: 'Sun BB et al. Plasma proteomic associations with genetics and
    health in the UK Biobank. Nature. 2023;622(7982):329-338. doi:10.1038/s41586-023-06592-6.
    PMID 37794186 (Results: discovery cohort and non-European pQTL mapping).'
  synapse_id: syn51469328
  source_url: https://www.synapse.org/Synapse:syn51469328
  fetched_at_utc: '2026-05-15T18:00:00Z'
attribution: 'UK Biobank Pharma Proteomics Project (Sun et al. 2023 Nature; PMID 37794186); accessed via Synapse syn51364943.'
outputs:
  variants_tsv: variants.tsv
```

`n_samples` is the per-variant N carried by the slice's own rows (all 5,205
rows state N = 33,628); `cohort_n` is the published EUR discovery cohort size
(Sun et al. 2023). They differ, and the report states both on separate lines:

```
- **Ancestry:** European (discovery) (EUR)
- **Sample size (per-variant N in the summary statistics):** 33,628
- **Published cohort size:** 34,557 (European (discovery); Sun et al. 2023). This is the cohort as published, not the N of this slice's rows.
```

`<output_dir>/variants.tsv`: the canonical 1p13.3 lead variant rs12740374
sits at chr1:109274968 (T/G in REGENIE convention) and shows the textbook
SORT1 plasma sortilin pQTL signal in EUR:

```
variant_id          chromosome  position_bp  allele_a  allele_b  beta       se        p          maf       molecular_trait_id  study_id
1_109274968_G_T     1           109274968    G         T         0.119857   0.00893   4.95e-41   0.220962  OID20213            syn51469328
```

That p-value (5e-41) is the protein-side counterpart to the classic SORT1 LDL-C / CHD pairing
(Musunuru 2010 *Nature*: minor allele -> increased hepatic SORT1 expression -> reduced plasma
LDL/VLDL). MAF ~0.22 matches the EUR-cohort literature.
