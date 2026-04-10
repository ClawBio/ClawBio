# Demo Data Provenance

## demo_olink_npx.csv.gz

- **Type**: Fully synthetic — no real patient data
- **Generated**: 2026-04-10
- **Samples**: 20 synthetic individuals (DEMO_000 to DEMO_019)
- **Proteins**: 26 Olink proteins (Heart: 7, Kidney: 9, Brain subset: 10)
- **Values**: NPX values drawn from N(3.0, 1.5), rounded to 3 decimal places
- **Metadata**: sample_id, age (N(60,12) clipped 30-85), sex (random M/F)
- **Purpose**: Quick demo and test execution — not for biological interpretation
- **Seed**: numpy random seed 42 for reproducibility
