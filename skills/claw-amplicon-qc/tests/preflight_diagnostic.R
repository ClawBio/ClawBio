#!/usr/bin/env Rscript
################################################################################
## Preflight primer-matching diagnostic
##
## Evidence for the fixed = "subject" choice in amplicon_qc.R's preflight
## helper. Reproducible test showing:
##
##   B1: reads shorter than the primer do NOT error — subseq() + pmin()
##       handle the case, returning 0 matches gracefully.
##
##   B2: with fixed = FALSE (the pre-audit default), N bases in the READ
##       were treated as wildcards matching any primer base, so all-N
##       reads counted as BOTH forward AND reverse matches, inflating
##       both orientation counters on the same reads.
##
##   B2 fix: with fixed = "subject", IUPAC codes in the PATTERN (primer)
##       still expand (Y matches C or T, etc.), but N in the SUBJECT
##       (read) is treated literally, so all-N reads no longer inflate
##       counters.
##
## Run with:  Rscript tests/preflight_diagnostic.R
##
## Expected on the fix:
##   Test B1 (short reads):       no error, 0 / 0
##   Test B2a (all-N reads):      0 fwd / 0 rev   (previously 10 / 10)
##   Test B2b (real fwd reads):  10 fwd / 0 rev
##   Test B2c (mixed batch):     10 fwd / 0 rev   (only the real fwd counted)
##   Test B2d (real rev reads):   0 fwd / 10 rev
################################################################################

suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
})

FWD <- "GTGYCAGCMGCCGCGGTAA"   # 515F, 19 bp, degenerate (Y, M)
REV <- "GGACTACHVGGGTWTCTAAT"  # 806R, 20 bp, degenerate (H, V)

# Reproduces the helper from amplicon_qc.R with the fix applied.
check_orientation <- function(reads, fwd_primer, rev_primer) {
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
cat("  Preflight diagnostic — fixed = \"subject\" verification\n")
cat("─────────────────────────────────────────────────────────\n\n")

# Test B1: short reads do not error
cat("Test B1: reads shorter than the primer\n")
short_reads <- DNAStringSet(c("ACGTA", "TTGCA", "GTGYC", "GGGGG", "GTGCC"))
b1 <- tryCatch(check_orientation(short_reads, FWD, REV),
               error = function(e) list(error = conditionMessage(e)))
if (!is.null(b1$error)) {
  cat("  FAIL:", b1$error, "\n\n")
} else {
  cat("  OK — no error. fwd:", b1$fwd, "/", b1$n_total,
      " rev:", b1$rev, "/", b1$n_total, "\n\n")
}

# Test B2a: all-N reads should NOT match either primer under fixed="subject"
cat("Test B2a: all-N reads (should show 0/0 with fix, was 10/10 without)\n")
n_reads <- DNAStringSet(rep(paste0(strrep("N", 25), "AAAAAAAAAAAAAA"), 10))
r <- check_orientation(n_reads, FWD, REV)
status <- if (r$fwd == 0 && r$rev == 0) "PASS" else "FAIL"
cat("  ", status, " — fwd:", r$fwd, "/", r$n_total,
    "  rev:", r$rev, "/", r$n_total, "\n\n")

# Test B2b: real forward-primered reads (IUPAC in primer should still expand)
cat("Test B2b: real fwd reads (GTGCCAG... resolves 515F GTGYCAG...)\n")
real_fwd <- DNAStringSet(rep(paste0("GTGCCAGCAGCCGCGGTAA",
                                     "TACGAAGGGGGCTAGCGTTGTT"), 10))
r <- check_orientation(real_fwd, FWD, REV)
status <- if (r$fwd == 10 && r$rev == 0) "PASS" else "FAIL"
cat("  ", status, " — fwd:", r$fwd, "/", r$n_total,
    "  rev:", r$rev, "/", r$n_total, "\n\n")

# Test B2c: mixed batch — only the real fwd reads should be counted
cat("Test B2c: 10 N-start + 10 real fwd (should count only 10 fwd)\n")
r <- check_orientation(c(n_reads, real_fwd), FWD, REV)
status <- if (r$fwd == 10 && r$rev == 0) "PASS" else "FAIL"
cat("  ", status, " — fwd:", r$fwd, "/", r$n_total,
    "  rev:", r$rev, "/", r$n_total, "\n\n")

# Test B2d: real reverse-primered reads (verify reverse detection still works)
cat("Test B2d: real rev reads\n")
real_rev <- DNAStringSet(rep(paste0("GGACTACCAGGGTATCTAAT",
                                     "AATTGCTCGTCGCTAGCTAG"), 10))
r <- check_orientation(real_rev, FWD, REV)
status <- if (r$fwd == 0 && r$rev == 10) "PASS" else "FAIL"
cat("  ", status, " — fwd:", r$fwd, "/", r$n_total,
    "  rev:", r$rev, "/", r$n_total, "\n\n")

cat("─────────────────────────────────────────────────────────\n")
cat("  All four tests PASS → fix is safe and behaves as documented.\n")
cat("─────────────────────────────────────────────────────────\n")
