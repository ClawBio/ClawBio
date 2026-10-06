# metagenome-mag-pipeline

A ClawBio wrapper around the upstream Snakemake workflow at
[mobashirrahman/metagenomics-workflow](https://github.com/mobashirrahman/metagenomics-workflow),
pinned to commit `7351a702d29857801963cae7aeff6600d59ebe77` (tag `v0.2.0`, MIT).

It validates a samplesheet, renders the pipeline config from a preset, checks the
machine can run it, launches snakemake, and turns the pipeline's `final/` tables
into `report.md` and `result.json`.

**Read [SKILL.md](SKILL.md)** for the trigger conditions, the CLI reference, the
output contract, the gotchas and the safety boundaries. This file is a pointer.

```bash
python clawbio.py run mag-pipeline --demo --check --output /tmp/plan
python clawbio.py run mag-pipeline --input data/reads --output out --preset assembly
```
