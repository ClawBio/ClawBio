#!/usr/bin/env Rscript
################################################################################
## B2 fix verification — does fixed="subject" behave correctly?
##
## Expected behaviour:
##   - N in read is literal (no wildcard) → all-N reads should NOT match either primer
##   - IUPAC in PRIMER still expands (Y matches C/T, etc.) → real forward reads
##     with GTGCCAG... at position 1 should still match FWD (GTGYCAG...)
################################################################################

suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
})

FWD <- "GTGYCAGCMGCCGCGGTAA"
REV <- "GGACTACHVGGGTWTCTAAT"

# The FIXED helper: fixed = "subject" instead of fixed = FALSE
check_orientation_fixed <- function(reads, fwd_primer, rev_primer) {
  n_total <- length(reads)
  if (n_total == 0) return(list(n_total = 0, fwd = 0, rev = 0))
  fwd_len <- nchar(fwd_primer)
  rev_len <- nchar(rev_primer)
  starts_fwd <- subseq(reads, 1, pmin(width(reads), fwd_len))
  starts_rev <- subseq(reads, 1, pmin(width(reads), rev_len))
  fwd_matches <- vcountPattern(fwd_primer, starts_fwd,
                                max.mismatch = 1, fixed = "subject")
  rev_matches <- vcountPattern(rev_primer, starts_rev,
                                max.mismatch = 1, fixed = "subject")
  list(n_total = n_total,
       fwd     = sum(fwd_matches > 0),
       rev     = sum(rev_matches > 0))
}

cat("─────────────────────────────────────────────────────────\n")
cat("  B2 FIX verification — fixed = \"subject\"\n")
cat("─────────────────────────────────────────────────────────\n\n")

# 1. All-N reads (previously matched BOTH primers)
cat("Test 1: all-N reads at start.\n")
n_reads <- DNAStringSet(rep(paste0(strrep("N", 25),
                                    "AAAAAAAAAAAAAA"), 10))
r1 <- check_orientation_fixed(n_reads, FWD, REV)
cat("  fwd matches:", r1$fwd, "/", r1$n_total, "\n")
cat("  rev matches:", r1$rev, "/", r1$n_total, "\n")
if (r1$fwd == 0 && r1$rev == 0) {
  cat("  PASS — N reads no longer inflating either counter.\n\n")
} else {
  cat("  FAIL — N reads still match. Fix is not correct.\n\n")
}

# 2. Real forward-primered reads (should still be detected — IUPAC in primer expands)
cat("Test 2: real forward primer reads (GTGCCAG… resolves 515F GTGYCAG…).\n")
real_fwd <- DNAStringSet(rep(paste0("GTGCCAGCAGCCGCGGTAA",
                                     "TACGAAGGGGGCTAGCGTTGTT"), 10))
r2 <- check_orientation_fixed(real_fwd, FWD, REV)
cat("  fwd matches:", r2$fwd, "/", r2$n_total, "\n")
cat("  rev matches:", r2$rev, "/", r2$n_total, "\n")
if (r2$fwd == 10 && r2$rev == 0) {
  cat("  PASS — real fwd reads still detected, IUPAC still expands in primer.\n\n")
} else {
  cat("  FAIL — fwd detection broken by the fix.\n\n")
}

# 3. Mixed (the scenario Manuel would care about most)
cat("Test 3: mixed batch (10 N-start + 10 real fwd).\n")
mixed <- c(n_reads, real_fwd)
r3 <- check_orientation_fixed(mixed, FWD, REV)
cat("  fwd matches:", r3$fwd, "/", r3$n_total, "\n")
cat("  rev matches:", r3$rev, "/", r3$n_total, "\n")
if (r3$fwd == 10 && r3$rev == 0) {
  cat("  PASS — only the 10 real fwd reads counted.\n")
  cat("  Verdict on this data: 10/20 = 50% fwd, 0% rev.\n")
  cat("  Would still fail the 80% consistent-forward threshold,\n")
  cat("  correctly triggering 'mixed' or 'primers_not_detected'.\n\n")
} else {
  cat("  UNEXPECTED — investigate.\n\n")
}

# 4. Real reverse-primered reads (verify all-reversed detection still works)
cat("Test 4: real reverse primer reads.\n")
real_rev <- DNAStringSet(rep(paste0("GGACTACCAGGGTATCTAAT",  # resolves 806R
                                     "AATTGCTCGTCGCTAGCTAG"), 10))
r4 <- check_orientation_fixed(real_rev, FWD, REV)
cat("  fwd matches:", r4$fwd, "/", r4$n_total, "\n")
cat("  rev matches:", r4$rev, "/", r4$n_total, "\n")
if (r4$fwd == 0 && r4$rev == 10) {
  cat("  PASS — reverse detection still works.\n\n")
} else {
  cat("  FAIL — reverse detection broken by the fix.\n\n")
}

cat("─────────────────────────────────────────────────────────\n")
cat("  If all four tests PASS, the fix is safe to apply.\n")
cat("─────────────────────────────────────────────────────────\n")
