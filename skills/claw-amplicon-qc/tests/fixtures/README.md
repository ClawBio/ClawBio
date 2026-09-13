# Test Fixtures

Bundled minimal FASTQ pair for the `--demo` mode and end-to-end test.

## Provenance

Derived from a subset of sample B4 from astrobiomike's Deep Sea Rock 16S
tutorial dataset (https://astrobiomike.github.io/amplicon/dada2_workflow_ex),
subsetted to the first 500 read pairs.

Original dataset: CC-BY 4.0 (astrobiomike / Michael D. Lee).

## Specifications

- Platform: Illumina MiSeq 2×300 paired-end
- Region: 16S rRNA V4
- Primers: 515F (GTGYCAGCMGCCGCGGTAA) and 806R (GGACTACNVGGGTWTCTAAT), attached
- Reads: 500 pairs
- Size: ~50KB per file

## Usage

Automatically used by `--demo` mode:

    Rscript amplicon_qc.R --demo --output /tmp/demo_output

Also used by `tests/run_test.sh` for CI.
