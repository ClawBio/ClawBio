# Synthetic Demo Fixture Provenance

The bundled demo is `synthetic_beta_demo.csv.gz`: two entirely simulated rows
with 485,512 public CpG/probe identifiers. It contains **no donor methylation
measurements**. The identifiers in `probe_ids.txt.gz` were extracted from the
header of the public GSE139307 example previously used by this skill; only the
column names were retained.

`generate_synthetic_demo.py` hashes each public probe identifier with seed
`20261019` and applies an integer-only bimodal transform. Each sample gets a
separately hashed offset of at most `0.04`; values are clipped to `[0, 1]` and
written to six decimal places. This makes regeneration byte-stable across
Python and NumPy versions. The fixture labels
both rows as `blood` and supplies placeholder ages (35 and 65 years) and sexes
so the five default clocks can be exercised. Those labels are **not encoded in
the simulated methylation values**. Clock predictions from this fixture test
software operation and must not be interpreted biologically.

Rebuild from the repository root:

```bash
python skills/methylation-clock/data/generate_synthetic_demo.py \
  --output skills/methylation-clock/data/synthetic_beta_demo.csv.gz
```

SHA-256 of `synthetic_beta_demo.csv.gz`:
`c6fe0fa1574ec9cd002128f53eaf0b67fce6374e1dc32821a3c24ecd3246ffca`

SHA-256 of `probe_ids.txt.gz`:
`cd72d91c20878ac2dff24f4b2a1c28a5eb03541b4b5339f8abe84699036b0c01`

ClawBio is a research and educational tool. It is not a medical device and does
not provide clinical diagnoses. Consult a healthcare professional before making
any medical decisions.
