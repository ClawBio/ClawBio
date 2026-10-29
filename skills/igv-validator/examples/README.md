# igv-validator examples (synthetic demo data only)

- `samples.csv`: a samplesheet for the demo files, one row per sample. Run from the repository root:
  ```bash
  python skills/igv-validator/igv_validator.py --samplesheet skills/igv-validator/examples/samples.csv \
    --reference skills/igv-validator/demo/demo_ref.fa --regions skills/igv-validator/examples/genes.bed \
    --curated-calls skills/igv-validator/examples/curated_calls.csv --no-igv --reports-dir /tmp/igv_reports
  ```
- `curated_calls.csv`: a curated table (`sample, gene, alteration`) to compare with; WRONGDEL is a planted wrong call.
- `demo_report.md`: the report `--demo --no-igv` writes for one check.
- `genes.bed`: the demo genes (chrom, 0-based start, end, name).
