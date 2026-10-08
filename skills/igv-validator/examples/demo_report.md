# IGV Validator Report

**Input**: `demo_calls.vcf` (6 variant(s) checked) · tumor: demo_tumor · normal: demo_normal
**Date**: 2026-10-08  
**Skill**: igv-validator 0.1.0

**3 supported, 2 flagged, 1 insufficient**.

IGV screenshots skipped: --no-igv given.

Read counts come from the BAMs (MAPQ >= 20, base quality >= 20, duplicate/secondary/supplementary reads removed, each DNA molecule counted once), never from the screenshots. **Status** is a rule-based summary of the flags, not a verdict on whether the variant is real.

| ID | Gene | Variant | Tumor support | Normal support | Caller reported (VCF) | Flags | Status |
|---|---|---|---|---|---|---|---|
| V1 | DEMO1 | demo1:3,000 C>A | 18/64 (28.1%) | 0/60 (0.0%) | 18/64 (28.1%) | none | **supported** |
| V2 | DEMO2 | demo1:5,000 A>T | 6/53 (11.3%) | 0/55 (0.0%) | 6/53 (11.3%) | strand_bias | **flagged** |
| V3 | DEMO3 | demo1:7,000 A>T | 2/46 (4.3%) | 0/60 (0.0%) | 12/56 (21.4%) | low_support; low_vaf; caller_disagrees | **insufficient** |
| V4 | DEMO4 | demo1:9,000 12 bp deletion | 15/60 (25.0%) | 0/63 (0.0%) | 15/60 (25.0%) | none | **supported** |
| V5_1 | DEMO5 | demo1:12,000 <-> demo2:5,001 breakend | 14 (6 split, 8 pairs) / 96 | 0 (0 split, 0 pairs) / 99 | - | none | **supported** |
| V6 | DEMO6 | demo1:15,000 A>G | 31/63 (49.2%) | 25/50 (50.0%) | 31/63 (49.2%) | normal_support | **flagged** |

## Flags

- `normal_support`: supporting reads in the normal
- `low_support`: fewer than 3 supporting reads
- `low_vaf`: under 5% of reads
- `strand_bias`: all supporting reads on one strand
- `read_end`: alt bases mostly within 10 bp of read ends
- `low_depth`: under 10 reads in tumor or normal
- `germline_site`: the normal carries another allele here (>= 20% of reads)
- `caller_disagrees`: the caller's reported counts (VCF AD) differ from the BAM by over 10 VAF points
- `no_coverage`: no reads at this position in either BAM (region not in the BAM, or wrong BAM/contig)
- `high_depth`: depth over 2.5x the median of this run (possible repeat or mismapped reads)

## Variants

### V1 DEMO1: supported

- demo_tumor: 18/64 (28.1%), 10 forward / 8 reverse
- demo_normal: 0/60 (0.0%)
- Caller reported (VCF): 18/64 (28.1%)
- Flags: none


### V2 DEMO2: flagged

- demo_tumor: 6/53 (11.3%), 0 forward / 6 reverse
- demo_normal: 0/55 (0.0%)
- Caller reported (VCF): 6/53 (11.3%)
- Flags: all supporting reads on one strand


### V3 DEMO3: insufficient

- demo_tumor: 2/46 (4.3%), 0 forward / 2 reverse
- demo_normal: 0/60 (0.0%)
- Caller reported (VCF): 12/56 (21.4%)
- Flags: fewer than 3 supporting reads, under 5% of reads, the caller's reported counts (VCF AD) differ from the BAM by over 10 VAF points


### V4 DEMO4: supported

- demo_tumor: 15/60 (25.0%), 8 forward / 7 reverse
- demo_normal: 0/63 (0.0%)
- Caller reported (VCF): 15/60 (25.0%)
- Flags: none


### V5_1 DEMO5: supported

- demo_tumor: 14 (6 split, 8 pairs) / 96
- demo_normal: 0 (0 split, 0 pairs) / 99
- Caller reported (VCF): -
- Flags: none


### V6 DEMO6: flagged

- demo_tumor: 31/63 (49.2%), 12 forward / 19 reverse
- demo_normal: 25/50 (50.0%), 17 forward / 8 reverse
- Caller reported (VCF): 31/63 (49.2%)
- Flags: supporting reads in the normal


---

## Disclaimer

*ClawBio is a research and educational tool. It is not a medical device and does not provide clinical diagnoses. Consult a healthcare professional before making any medical decisions.*
