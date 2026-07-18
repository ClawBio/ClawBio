# population-equity-auditor / data

## Bundled files

| File | Provenance | Notes |
|------|-----------|-------|
| `clingen/ClinGen_haploinsufficiency_GRCh38.tsv` | ClinGen Dosage Sensitivity (public) | gene-level haploinsufficiency scores for the PVS1 gate |
| `clingen/PROVENANCE.json` | — | source URL + retrieval note |
| `gnomad/gnomad_constraint_by_gene.tsv` | gnomAD (public) | per-gene LOEUF / pLI for the LoF-constraint gate |
| `demo_cohort.tsv` + `demo_cohort.meta.json` | **SYNTHETIC** | a fabricated demonstration cohort — **not real patient data** |

## No patient data

This skill ships **no real cohort or patient data**. The demo is synthetic (fabricated allele counts,
placeholder population labels `POP_A`…`POP_G`, mostly placeholder gene names). It exists only to exercise
the classification mechanism. Only public gene-level reference tables (ClinGen, gnomAD constraint) are bundled.

## Bring your own cohort

The skill is cohort-agnostic. Provide two files:

1. `mycohort.tsv` — one row per variant, tab-separated, columns:
   `LOCATION REF ALT SYMBOL CADD_PHRED Consequence CLIN_SIG <one column per subpopulation of ALT allele counts> GNOMAD_AF GNOMAD_AF_AMR GNOMAD_AN GNOMAD_AN_AMR`
2. `mycohort.meta.json` — declares the cohort:
   ```json
   {
     "cohort_id": "mycohort",
     "assay": "WGS",
     "caller": "your-caller",
     "build": "GRCh38",
     "variant_classes": ["SNV"],
     "subpopulation_sizes": {"POP_A": 30, "POP_B": 16},
     "notes": "assay/build caveats a reviewer should know"
   }
   ```

Cohort allele frequency is computed as `sum(subpopulation ALT allele counts) / (2 * total N)`. gnomAD +
consequence + CADD annotation is expected to be produced upstream (e.g. by `variant-annotation` /
`clinical-variant-reporter`); for a sealed clinical deployment, run that annotation **locally** so no
cohort-derived coordinates leave the environment.
