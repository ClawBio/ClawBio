# prs-abstain v0.3.0 — how to run it

## Install

```bash
git clone https://github.com/ClawBio/ClawBio.git && cd ClawBio
curl -LsSf https://astral.sh/uv/install.sh | sh
uv sync && uv pip install matplotlib reportlab pillow pypdf pytest
```

Then drop `skills/prs-abstain/` into place (from the tarball) and register the alias by adding
this block to the `SKILLS` dict in `clawbio/cli.py`, just above `"just-prs"`:

```python
"prs-abstain": {
    "script": SKILLS_DIR / "prs-abstain" / "prs_abstain.py",
    "demo_args": ["--demo"],
    "description": "Gate PRS percentiles on ancestry transferability; refuse when uninterpretable",
    "allowed_extra_flags": {
        "--reference-panel", "--individuals", "--prs-results", "--ref-pop", "--k-sd",
        "--min-markers", "--pcs", "--no-figures", "--scores", "--genotype",
        "--population-af", "--af-population", "--min-weight-coverage", "--min-effective-n",
        "--ld-window-kb", "--no-pdf",
    },
},
```

## Run the demo

```bash
uv run python skills/prs-abstain/prs_abstain.py --demo --output /tmp/prs-abstain
# or through the ClawBio runner, once the alias is registered
uv run python clawbio.py run prs-abstain --demo
```

Takes about 15 seconds. Produces three individuals, three verdicts, six scores each.

## In the hosted BioNeMo agent

```
Use your ClawBio skill-listing tool to find prs-abstain, describe its contract,
then run it in demo mode and show me the clinician report and the technical report.
```

The demo needs no input files, so it works even where `input_path` is refused.

## Real data

```bash
uv run python skills/prs-abstain/prs_abstain.py \
  --reference-panel /tmp/pca/tables/pc_coordinates.csv \
  --individuals my_cohort.csv \
  --prs-results /tmp/prs/prs_results.json \
  --scores skills/gwas-prs/data \
  --genotype my_genotype.txt \
  --population-af gnomad_af.tsv --af-population AFR \
  --output /tmp/run
```

`pc_coordinates.csv` comes straight from `claw-ancestry-pca`. `prs_results.json` comes straight
from `gwas-prs`. Neither needs reformatting.

**Individuals CSV** (`sex` is optional but sex-specific scores are withheld without it):

```csv
sample_id,population,PC1,PC2,PC3,PC4,n_markers_shared,sex
P001,,-3.01,-2.39,-2.00,0.35,480,female
```

**Population AF TSV** — the bundled demo table is synthetic and must not be used for
interpretation. Supply real gnomAD or 1000 Genomes frequencies:

```
rsid	population	effect_allele_frequency
rs7903146	AFR	0.2900
```

## Useful flags

| Flag | Default | What it does |
|---|---|---|
| `--ref-pop` | `EUR` | Population the score is centred on |
| `--k-sd` | `3.0` | Threshold = mean + k·sd of within-reference spread |
| `--min-markers` | `30` | Below this, ancestry is undeterminable (Kosoy 2009) |
| `--min-weight-coverage` | `0.90` | Fraction of a score's total weight that must be genotyped |
| `--min-effective-n` | `10` | Warn below this many independent contributions |
| `--af-population` | `AFR` | Population column to re-centre on |
| `--ld-window-kb` | `250` | Group variants this close as potentially correlated |
| `--no-figures` | off | Skip plots (runs without matplotlib) |
| `--no-pdf` | off | Skip PDF rendering |

## Output

```
report_clinician.md    plain language, no jargon
report_technical.md    mechanism, per-variant audit, LD, limitations
report_clinician.pdf   typeset, 4 pages
report_technical.pdf   typeset, 8 pages
result.json            machine-readable verdicts and provenance
tables/                decisions, gated_scores, variant_audit, af_shift_per_variant
figures/               gate overview + one panel per individual
reproducibility/       commands.sh, environment.yml
```

## Tests

```bash
uv run python -m pytest skills/prs-abstain/tests/ -v   # 59 tests
```

## Reading the LD numbers

The technical report gives two effective-n figures per score. The first assumes every variant
is an independent piece of evidence; the second groups variants that sit within
`--ld-window-kb` of each other. Quote the second. On the bundled scores PGS000057 drops from
122.5 to 35.9, and PGS000013's two largest weights turn out to be one locus.

Physical distance is a proxy, not a measurement. It finds correlated variants that sit close
together and misses correlated variants that do not, so the grouped figure is still an
overestimate. Real r-squared needs a haplotype reference panel.

## Two things to check before trusting a run

`--k-sd` above about 5 makes the threshold wider than the distance to the nearest non-reference
individual. The tool detects this and stamps **THRESHOLD OVERREACH** on the report. If you see
that, the REPORT verdicts are meaningless.

The bundled demo allele frequencies are synthetic. The `shift (sd)` column shows real
arithmetic on fake numbers, so use it to demonstrate the method and never to interpret a person.
