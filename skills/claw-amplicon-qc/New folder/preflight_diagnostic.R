#!/usr/bin/env Rscript
################################################################################
## Diagnostic script for Manuel's re-audit questions B1 and B2
##
## B1: does vcountPattern error rather than return 0 when a read is shorter
##     than the primer?
## B2: does fixed = FALSE apply IUPAC ambiguity to the SUBJECT (reads) too,
##     inflating both orientation counters when raw reads have N bases?
##
## Run with:  Rscript /tmp/preflight_diagnostic.R
################################################################################

suppressPackageStartupMessages({
  library(Biostrings)
  library(ShortRead)
})

FWD <- "GTGYCAGCMGCCGCGGTAA"   # 515F, 19 bp
REV <- "GGACTACHVGGGTWTCTAAT"  # 806R, 20 bp

cat("─────────────────────────────────────────────────────────\n")
cat("  Diagnostic — reproduces the preflight helper function\n")
cat("─────────────────────────────────────────────────────────\n\n")

# --- reproduce the exact helper from amplicon_qc.R ---
check_orientation <- function(reads, fwd_primer, rev_primer) {
  n_total <- length(reads)
  if (n_total == 0) return(list(n_total = 0, fwd = 0, rev = 0))
  fwd_len <- nchar(fwd_primer)
  rev_len <- nchar(rev_primer)
  starts_fwd <- subseq(reads, 1, pmin(width(reads), fwd_len))
  starts_rev <- subseq(reads, 1, pmin(width(reads), rev_len))
  fwd_matches <- vcountPattern(fwd_primer, starts_fwd,
                                max.mismatch = 1, fixed = FALSE)
  rev_matches <- vcountPattern(rev_primer, starts_rev,
                                max.mismatch = 1, fixed = FALSE)
  list(n_total = n_total,
       fwd     = sum(fwd_matches > 0),
       rev     = sum(rev_matches > 0))
}

################################################################################
## TEST B1 — reads shorter than the primer
################################################################################

cat("### TEST B1: reads shorter than the primer\n")
cat("Creating 5 reads of 5 bp each (shorter than 19-bp FWD primer).\n\n")

short_reads <- DNAStringSet(c(
  "ACGTA",   # 5 bp
  "TTGCA",   # 5 bp
  "GTGYC",   # 5 bp — starts with FWD but only 5 bases
  "GGGGG",   # 5 bp
  "GTGCC"    # 5 bp — matches first 5 of 515F if degeneracy resolves
))

result_b1 <- tryCatch(
  check_orientation(short_reads, FWD, REV),
  error = function(e) list(error = conditionMessage(e))
)

if (!is.null(result_b1$error)) {
  cat("  Result: ERROR — Manuel's B1 concern is confirmed.\n")
  cat("  Error message:", result_b1$error, "\n\n")
} else {
  cat("  Result: NO error.\n")
  cat("  fwd matches:", result_b1$fwd, "/", result_b1$n_total, "\n")
  cat("  rev matches:", result_b1$rev, "/", result_b1$n_total, "\n")
  cat("  Verdict: preflight tolerates short reads — B1 concern is unfounded.\n\n")
}

################################################################################
## TEST B2 — N in raw reads with fixed = FALSE
################################################################################

cat("### TEST B2: N-in-read behaviour with fixed = FALSE\n\n")

# Case B2a: reads where the entire primer region is Ns
cat("Case B2a: reads with all-N start (19 Ns at position 1).\n")
n_reads <- DNAStringSet(rep(paste0(strrep("N", 25),
                                    "AAAAAAAAAAAAAA"), 10))
result_b2a <- check_orientation(n_reads, FWD, REV)
cat("  fwd matches:", result_b2a$fwd, "/", result_b2a$n_total, "\n")
cat("  rev matches:", result_b2a$rev, "/", result_b2a$n_total, "\n")
b2a_inflated <- (result_b2a$fwd > 0 && result_b2a$rev > 0)
if (b2a_inflated) {
  cat("  BOTH counters non-zero on all-N reads → B2 concern CONFIRMED.\n")
  cat("  N-in-read is matching any primer base.\n\n")
} else {
  cat("  Not both counters non-zero → B2 concern NOT confirmed for all-N reads.\n\n")
}

# Case B2b: reads with real forward primer + Ns after
cat("Case B2b: reads with real forward primer at pos 1 (unambiguous case).\n")
real_fwd_reads <- DNAStringSet(rep(paste0("GTGCCAGCAGCCGCGGTAA",  # resolves 515F
                                           "TACGAAGGGGGCTAGCGTTGTT"), 10))
result_b2b <- check_orientation(real_fwd_reads, FWD, REV)
cat("  fwd matches:", result_b2b$fwd, "/", result_b2b$n_total, "\n")
cat("  rev matches:", result_b2b$rev, "/", result_b2b$n_total, "\n")
if (result_b2b$fwd == result_b2b$n_total && result_b2b$rev == 0) {
  cat("  Only fwd matched — clean orientation detection on unambiguous reads.\n\n")
} else {
  cat("  Unexpected: both counters fired even without Ns. Look at input.\n\n")
}

# Case B2c: mixed — half N-starting, half real forward
cat("Case B2c: mixed batch, half N-starting + half real forward-oriented.\n")
mixed <- c(n_reads, real_fwd_reads)  # 20 reads: 10 N-start, 10 real FWD
result_b2c <- check_orientation(mixed, FWD, REV)
cat("  fwd matches:", result_b2c$fwd, "/", result_b2c$n_total, "\n")
cat("  rev matches:", result_b2c$rev, "/", result_b2c$n_total, "\n")
if (result_b2c$fwd > 10 || result_b2c$rev > 0) {
  cat("  N-reads are inflating counters — orientation verdict on real data\n")
  cat("  could be wrong for datasets with N-heavy reads.\n\n")
} else {
  cat("  Only real fwd reads counted — mixed case behaves correctly.\n\n")
}

################################################################################
## SUMMARY
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Summary for reply to Manuel\n")
cat("─────────────────────────────────────────────────────────\n")

b1_status <- if (!is.null(result_b1$error)) "CONFIRMED (needs fix)" else "not confirmed"
b2_status <- if (b2a_inflated) "CONFIRMED (needs fix)" else "not confirmed"

cat(sprintf("  B1 (short reads error):        %s\n", b1_status))
cat(sprintf("  B2 (N-in-read inflates fwd+rev): %s\n", b2_status))
