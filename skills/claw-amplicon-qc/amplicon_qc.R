#!/usr/bin/env Rscript
################################################################################
##                    CLAW-AMPLICON-QC — v0.1                                 ##
##                                                                            ##
##   Non-interactive amplicon preprocessing skill for the ClawBio ecosystem. ##
##   Measures raw reads with seqkit, removes reads containing N bases,       ##
##   measures again, trims primers with Cutadapt, measures a third time,     ##
##   and writes a machine-readable JSON summary plus a human-readable        ##
##   markdown report.                                                        ##
##                                                                            ##
##   Author: Zabiulla                                                         ##
##   License: MIT                                                             ##
################################################################################

suppressPackageStartupMessages({
  required_pkgs <- c("dada2", "ShortRead", "Biostrings", "jsonlite", "optparse")
  missing_pkgs  <- required_pkgs[!vapply(required_pkgs, requireNamespace,
                                         logical(1), quietly = TRUE)]
  if (length(missing_pkgs) > 0) {
    stop(
      "Missing required R packages: ", paste(missing_pkgs, collapse = ", "), "\n",
      "Install with: BiocManager::install(c('dada2','ShortRead','Biostrings')); ",
      "install.packages(c('jsonlite','optparse'))",
      call. = FALSE
    )
  }
  library(dada2)
  library(ShortRead)
  library(Biostrings)
  library(jsonlite)
  library(optparse)
})

################################################################################
## PART 1: PARSE COMMAND-LINE ARGUMENTS
################################################################################

option_list <- list(
  make_option(c("--raw"), type = "character", default = NULL,
              help = "Path to folder with raw paired-end FASTQ files"),
  make_option(c("--output"), type = "character", default = NULL,
              help = "Path to output folder (created if missing)"),
  make_option(c("--fwd-primer"), type = "character", default = NULL,
              dest = "fwd_primer",
              help = "Forward primer sequence (IUPAC nucleotide codes)"),
  make_option(c("--rev-primer"), type = "character", default = NULL,
              dest = "rev_primer",
              help = "Reverse primer sequence (IUPAC nucleotide codes)"),
  make_option(c("--min-length"), type = "integer", default = NULL,
              dest = "min_length",
              help = "Minimum read length in bp after primer trimming")
)

opt <- parse_args(OptionParser(
  option_list = option_list,
  description = "claw-amplicon-qc — 16S amplicon preprocessing skill"
))

required_args <- c("raw", "output", "fwd_primer", "rev_primer", "min_length")
missing_args  <- required_args[vapply(required_args,
                                      function(a) is.null(opt[[a]]),
                                      logical(1))]
if (length(missing_args) > 0) {
  stop("Missing required arguments: --",
       paste(gsub("_", "-", missing_args), collapse = " --"),
       "\nRun with --help for usage.", call. = FALSE)
}
if (!dir.exists(opt$raw)) {
  stop("Raw folder does not exist: ", opt$raw, call. = FALSE)
}

iupac_pattern <- "^[ACGTMRWSYKVHDBN]+$"
if (!grepl(iupac_pattern, toupper(opt$fwd_primer))) {
  stop("Forward primer contains invalid characters (must be IUPAC codes): ",
       opt$fwd_primer, call. = FALSE)
}
if (!grepl(iupac_pattern, toupper(opt$rev_primer))) {
  stop("Reverse primer contains invalid characters (must be IUPAC codes): ",
       opt$rev_primer, call. = FALSE)
}
opt$fwd_primer <- toupper(opt$fwd_primer)
opt$rev_primer <- toupper(opt$rev_primer)

for (tool in c("seqkit", "cutadapt")) {
  # CRITICAL 3: use Sys.which() rather than the previous try-error check.
  # system2() with stdout=TRUE never returns a try-error object, so
  # inherits(chk, "try-error") could never fire — the check was decorative.
  # Sys.which() returns "" when a tool is not on PATH.
  if (Sys.which(tool) == "") {
    stop(tool, " not found on PATH.\n",
         "  Install with: conda install -c bioconda ", tool, "\n",
         "  If already installed, activate the conda env first:\n",
         "    conda activate amplicon-qc",
         call. = FALSE)
  }
}

if (!dir.exists(opt$output)) dir.create(opt$output, recursive = TRUE)

start_time <- Sys.time()

cat("─────────────────────────────────────────────────────────\n")
cat("  claw-amplicon-qc v0.1.0\n")
cat("─────────────────────────────────────────────────────────\n")
cat("  Raw folder:    ", opt$raw,        "\n")
cat("  Output folder: ", opt$output,     "\n")
cat("  Fwd primer:    ", opt$fwd_primer, "\n")
cat("  Rev primer:    ", opt$rev_primer, "\n")
cat("  Min length:    ", opt$min_length, "\n")
cat("─────────────────────────────────────────────────────────\n\n")

################################################################################
## PART 2: DISCOVER FASTQ FILE PAIRS
##
## Supports two common paired-end FASTQ naming conventions:
##   Pattern A — Illumina default: SAMPLE_R1_...fastq.gz / SAMPLE_R2_...fastq.gz
##   Pattern B — R1/R2 prefix:     R1_SAMPLE.fastq.gz     / R2_SAMPLE.fastq.gz
## The script auto-detects which pattern the raw folder uses.
################################################################################

# All FASTQ files (accept .fastq.gz, .fq.gz, .fastq, .fq)
all_fastq <- list.files(opt$raw,
                        pattern     = "\\.(fastq|fq)(\\.gz)?$",
                        full.names  = TRUE,
                        ignore.case = TRUE)

if (length(all_fastq) == 0) {
  stop("No FASTQ files found in ", opt$raw,
       "\nExpected files ending in .fastq.gz, .fq.gz, .fastq, or .fq",
       call. = FALSE)
}

# CRITICAL 2: pair R1 and R2 by sample name, not by sort position.
#
# The previous approach sorted R1 and R2 vectors independently and zipped
# them positionally, relying on the two sorts producing aligned orderings.
# That silently fails if:
#   - A sample is missing R1 or R2 (counts stay equal if a different sample
#     is also missing its other partner — sample A's R1 pairs with sample B's R2)
#   - Filenames sort in different orders on the two sides (inconsistent
#     capitalisation, padding, or locale collation)
#   - A stray non-sample file matches the R1 or R2 glob
#
# The audit noted this matters when the skill is pointed at data the user
# did not generate (collaborator tarballs, SRA dumps, provider deliveries).
# On healthy data (matched R1/R2 partners), this fix produces identical
# outputs to the previous approach — verified across three datasets during
# the audit fixes.

# Detect naming pattern and set up per-pattern sample-name extraction
# Pattern A (Illumina default): SAMPLE_R1_..., SAMPLE_R2_...
# Pattern B (prefix):           R1_SAMPLE,    R2_SAMPLE
r1_A_idx <- grep("_R1[_.]", basename(all_fastq))
r2_A_idx <- grep("_R2[_.]", basename(all_fastq))
r1_B_idx <- grep("^R1_",    basename(all_fastq))
r2_B_idx <- grep("^R2_",    basename(all_fastq))

if (length(r1_A_idx) > 0 && length(r2_A_idx) > 0) {
  r1_paths <- all_fastq[r1_A_idx]
  r2_paths <- all_fastq[r2_A_idx]
  naming_pattern <- "Illumina default (SAMPLE_R1_..., SAMPLE_R2_...)"
  # Sample name = strip _R1/_R2 and everything after
  extract_sample_name <- function(paths, which_read) {
    pat <- sprintf("_R%d[_.].*", which_read)
    gsub(pat, "", basename(paths))
  }
} else if (length(r1_B_idx) > 0 && length(r2_B_idx) > 0) {
  r1_paths <- all_fastq[r1_B_idx]
  r2_paths <- all_fastq[r2_B_idx]
  naming_pattern <- "Prefix (R1_SAMPLE, R2_SAMPLE)"
  # Sample name = strip R1_/R2_ prefix and file extension
  extract_sample_name <- function(paths, which_read) {
    pat <- sprintf("^R%d_", which_read)
    n <- gsub(pat, "", basename(paths))
    gsub("\\.(fastq|fq)(\\.gz)?$", "", n, ignore.case = TRUE)
  }
} else {
  stop("No paired-end FASTQ files found in ", opt$raw,
       "\nExpected either:",
       "\n  Pattern A: SAMPLE_R1_...fastq.gz + SAMPLE_R2_...fastq.gz",
       "\n  Pattern B: R1_SAMPLE.fastq.gz + R2_SAMPLE.fastq.gz",
       call. = FALSE)
}

# Extract sample names for R1 and R2 files, independently
r1_names <- extract_sample_name(r1_paths, 1)
r2_names <- extract_sample_name(r2_paths, 2)

# ASSERTION 1: no duplicate sample names on the R1 side
r1_dupes <- unique(r1_names[duplicated(r1_names)])
if (length(r1_dupes) > 0) {
  stop("Duplicate R1 sample names detected — same sample appears twice:\n  ",
       paste(r1_dupes, collapse = ", "),
       "\nThis usually means leftover files from a previous run. ",
       "Clean up the raw folder and try again.",
       call. = FALSE)
}

# ASSERTION 2: no duplicate sample names on the R2 side
r2_dupes <- unique(r2_names[duplicated(r2_names)])
if (length(r2_dupes) > 0) {
  stop("Duplicate R2 sample names detected — same sample appears twice:\n  ",
       paste(r2_dupes, collapse = ", "),
       "\nThis usually means leftover files from a previous run. ",
       "Clean up the raw folder and try again.",
       call. = FALSE)
}

# ASSERTION 3: every R1 sample has an R2 partner, and vice versa
r1_without_r2 <- setdiff(r1_names, r2_names)
r2_without_r1 <- setdiff(r2_names, r1_names)
if (length(r1_without_r2) > 0 || length(r2_without_r1) > 0) {
  msg <- "R1/R2 pairing failed: some samples have no partner.\n"
  if (length(r1_without_r2) > 0) {
    msg <- paste0(msg, "  R1 files with no matching R2:\n    ",
                  paste(r1_without_r2, collapse = "\n    "), "\n")
  }
  if (length(r2_without_r1) > 0) {
    msg <- paste0(msg, "  R2 files with no matching R1:\n    ",
                  paste(r2_without_r1, collapse = "\n    "), "\n")
  }
  msg <- paste0(msg,
                "Every sample must have both an R1 and R2 file. ",
                "Fix the raw folder and try again.")
  stop(msg, call. = FALSE)
}

# Build fnFs/fnRs by name-join — alignment guaranteed by construction,
# not by hope that two independent sorts agree.
sample_names <- sort(r1_names)
fnFs <- r1_paths[match(sample_names, r1_names)]
fnRs <- r2_paths[match(sample_names, r2_names)]

# Sanity check: verify alignment holds after the join. If any of these fail,
# the pairing logic itself has a bug (should be impossible given the above
# assertions, but defensive checking is cheap).
stopifnot(
  length(fnFs) == length(fnRs),
  length(fnFs) == length(sample_names),
  identical(extract_sample_name(fnFs, 1), sample_names),
  identical(extract_sample_name(fnRs, 2), sample_names)
)

cat("Detected naming pattern:", naming_pattern, "\n")
cat("Discovered", length(fnFs), "sample pairs (paired by sample name).\n\n")

################################################################################
## PART 3: HELPER — run seqkit stats and parse the output
################################################################################

run_seqkit_stats <- function(files, output_file) {
  # CRITICAL 1: shQuote() every path passed as an argument to system2().
  # system2 pastes args into a shell string, so unquoted paths with spaces
  # or shell metacharacters break the call (or execute embedded commands).
  # stdout=output_file is quoted internally by system2, so it is safe as-is.
  #
  # CRITICAL 3: capture stderr into a temp file so the message is visible
  # if seqkit fails, and check the exit status. Previously stderr="" threw
  # all error messages in the bin, and the integer return value was ignored
  # — a failed seqkit call would produce an empty output file, downstream
  # code would parse an empty table, and the run would continue silently.
  stderr_file <- tempfile()
  on.exit(unlink(stderr_file), add = TRUE)
  status <- system2("seqkit", args = c("stats", "-a", "-T", shQuote(files)),
                    stdout = output_file, stderr = stderr_file)
  if (status != 0) {
    err_msg <- paste(readLines(stderr_file, warn = FALSE), collapse = "\n")
    stop("seqkit failed (exit code ", status, "):\n", err_msg,
         call. = FALSE)
  }
  read.table(output_file, sep = "\t", header = TRUE,
             stringsAsFactors = FALSE, check.names = FALSE)
}

seqkit_reads_by_file <- function(stats_df) {
  keys   <- basename(stats_df$file)
  counts <- as.numeric(gsub(",", "", stats_df$num_seqs))
  names(counts) <- keys
  counts
}

################################################################################
## PART 4: STAGE 1 — BASELINE SEQKIT ON RAW READS
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Stage 1/5 — Baseline seqkit stats on raw reads\n")
cat("─────────────────────────────────────────────────────────\n")

raw_stats_file <- file.path(opt$output, "01_raw_stats.txt")
raw_stats      <- run_seqkit_stats(c(fnFs, fnRs), raw_stats_file)
raw_counts     <- seqkit_reads_by_file(raw_stats)

cat("  Total raw reads:", format(sum(raw_counts), big.mark = ","), "\n")
cat("  Stats saved to: ", raw_stats_file, "\n\n")

################################################################################
## PART 5: STAGE 2 — N-FILTER (filterAndTrim, maxN = 0 ONLY)
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Stage 2/5 — Removing reads containing N bases\n")
cat("─────────────────────────────────────────────────────────\n")

filtN_fp <- file.path(opt$output, "filtN")
if (!dir.exists(filtN_fp)) dir.create(filtN_fp, recursive = TRUE)

fnFs.filtN <- file.path(filtN_fp, basename(fnFs))
fnRs.filtN <- file.path(filtN_fp, basename(fnRs))

is_windows  <- Sys.info()[["sysname"]] == "Windows"
multithread <- FALSE

filterAndTrim(
  fnFs, fnFs.filtN,
  fnRs, fnRs.filtN,
  maxN        = 0,
  multithread = multithread,
  compress    = TRUE,
  verbose     = TRUE
)
cat("\n")

################################################################################
## PART 6: STAGE 3 — SEQKIT AFTER N-FILTER
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Stage 3/5 — seqkit stats on N-filtered reads\n")
cat("─────────────────────────────────────────────────────────\n")

filtN_stats_file <- file.path(opt$output, "02_filtN_stats.txt")
filtN_stats      <- run_seqkit_stats(c(fnFs.filtN, fnRs.filtN), filtN_stats_file)
filtN_counts     <- seqkit_reads_by_file(filtN_stats)

cat("  Total reads after N removal:",
    format(sum(filtN_counts), big.mark = ","), "\n")
cat("  Stats saved to: ", filtN_stats_file, "\n\n")

################################################################################
## PART 7: STAGE 4 — CUTADAPT PRIMER TRIMMING
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Stage 4/5 — Primer trimming with Cutadapt\n")
cat("─────────────────────────────────────────────────────────\n")

FWD    <- opt$fwd_primer
REV    <- opt$rev_primer
FWD.RC <- dada2::rc(FWD)
REV.RC <- dada2::rc(REV)

cat("  Fwd:", FWD, "  Fwd RC:", FWD.RC, "\n")
cat("  Rev:", REV, "  Rev RC:", REV.RC, "\n")
cat("  Min length after trimming:", opt$min_length, "bp\n\n")

trimmed_fp <- file.path(opt$output, "cutadapt_trimmed")
if (!dir.exists(trimmed_fp)) dir.create(trimmed_fp, recursive = TRUE)

fnFs.cut <- file.path(trimmed_fp, basename(fnFs))
fnRs.cut <- file.path(trimmed_fp, basename(fnRs))

cutadapt_log <- file.path(opt$output, "04_cutadapt_log.txt")
if (file.exists(cutadapt_log)) file.remove(cutadapt_log)

# CRITICAL 3: track cutadapt success/failure per sample. Previously the return
# value of system2() was discarded — a failed cutadapt run wrote nothing but
# the loop continued to the next sample and the report showed NA/0 for that
# sample without any indication that cutadapt itself had errored.
# Behavior (b): on failure, warn loudly, skip this sample from Stage 5 seqkit
# and downstream stats, mark the sample in the report, continue.
failed_samples  <- character(0)
successful_idx  <- integer(0)

for (i in seq_along(fnFs)) {
  cat(sprintf("  [%d/%d] %s\n", i, length(fnFs), sample_names[i]))
  # CRITICAL 1: shQuote() every path argument. Primer sequences are already
  # validated by the IUPAC regex allowlist upstream; --minimum-length is an
  # integer coerced by optparse. Only the four paths carry shell-injection
  # risk (or, more commonly, spaces in output dirs like "/mnt/d/Research Data/").
  args_vec <- c(
    "-g", FWD, "-a", REV.RC,
    "-G", REV, "-A", FWD.RC,
    "--minimum-length", as.character(opt$min_length),
    "--nextseq-trim=20",
    "-n", "2",
    "-j", "0",
    "-o", shQuote(fnFs.cut[i]),
    "-p", shQuote(fnRs.cut[i]),
    shQuote(fnFs.filtN[i]), shQuote(fnRs.filtN[i])
  )
  cat("\n\n===== Sample:", sample_names[i], "=====\n",
      file = cutadapt_log, append = TRUE)
  # CRITICAL 3: capture the integer exit status. Cutadapt returns 0 on success,
  # non-zero on error. stdout+stderr already go to cutadapt_log so any error
  # message is preserved in the log for post-mortem investigation.
  status <- system2("cutadapt", args = args_vec,
                    stdout = cutadapt_log, stderr = cutadapt_log, wait = TRUE)
  if (status != 0) {
    warning("Cutadapt failed on sample '", sample_names[i],
            "' (exit code ", status, "). Skipping this sample. ",
            "See ", cutadapt_log, " for the per-sample cutadapt output.",
            call. = FALSE, immediate. = TRUE)
    failed_samples <- c(failed_samples, sample_names[i])
  } else {
    successful_idx <- c(successful_idx, i)
  }
}
cat("\n")

if (length(failed_samples) > 0) {
  cat("  ⚠ Cutadapt failed on", length(failed_samples), "sample(s):",
      paste(failed_samples, collapse = ", "), "\n")
  cat("  These samples are excluded from downstream stats but appear in the report.\n\n")
}

################################################################################
## PART 8: STAGE 5 — SEQKIT AFTER CUTADAPT
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Stage 5/5 — seqkit stats on primer-trimmed reads\n")
cat("─────────────────────────────────────────────────────────\n")

trimmed_stats_file <- file.path(opt$output, "03_trimmed_stats.txt")

# CRITICAL 3: only feed seqkit files that actually exist. If cutadapt failed
# on a sample above (behavior b: skip and continue), that sample's trimmed
# output was never written, and passing a missing path to seqkit would abort
# the whole Stage 5 with a confusing error. Filtering by file.exists() makes
# Stage 5 robust to any reason a file might be absent, not just cutadapt
# failure.
existing_trimmed <- c(fnFs.cut, fnRs.cut)
existing_trimmed <- existing_trimmed[file.exists(existing_trimmed)]

if (length(existing_trimmed) == 0) {
  stop("Cutadapt failed on every sample — no trimmed FASTQ files exist. ",
       "See ", cutadapt_log, " for per-sample diagnostics.",
       call. = FALSE)
}

trimmed_stats  <- run_seqkit_stats(existing_trimmed, trimmed_stats_file)
trimmed_counts <- seqkit_reads_by_file(trimmed_stats)

cat("  Total reads after primer trimming:",
    format(sum(trimmed_counts), big.mark = ","), "\n")
cat("  Stats saved to: ", trimmed_stats_file, "\n\n")

# CRITICAL 3: 100%-data-loss backstop. If cutadapt "succeeded" on every sample
# (exit code 0) but discarded >99% of reads across the whole run, this is a
# configuration error (wrong primers, wrong --min-length for the read length,
# wrong --nextseq-trim for the platform), not a per-sample data-quality issue.
# Aborting here beats writing a report full of zeros — which is exactly what
# happened on the pre-trimmed raw2 dataset before this check existed.
if (length(successful_idx) > 0) {
  successful_bn <- basename(fnFs)[successful_idx]
  total_in  <- sum(filtN_counts  [successful_bn], na.rm = TRUE)
  total_out <- sum(trimmed_counts[successful_bn], na.rm = TRUE)
  if (total_in > 0 && (total_out / total_in) < 0.01) {
    stop(
      "CONFIGURATION ERROR: >99% of reads discarded at the primer-trim stage ",
      "across all successful samples.\n",
      "  Reads in (post-N-filter): ", format(total_in,  big.mark = ","), "\n",
      "  Reads out (post-cutadapt): ", format(total_out, big.mark = ","), "\n",
      "  Retention: ", sprintf("%.2f%%", 100 * total_out / total_in), "\n\n",
      "This is almost certainly one of:\n",
      "  1. Wrong primers specified. Verify against your library prep records.\n",
      "  2. Data already primer-trimmed by the sequencing facility.\n",
      "     Check with: zcat <one_R1.fastq.gz> | head -2\n",
      "     If reads don't start with your forward primer, this skill is not the\n",
      "     right tool — go directly to DADA2 quality filtering.\n",
      "  3. --min-length (", opt$min_length,
      " bp) too high for the actual read length after quality trimming.\n\n",
      "See ", cutadapt_log, " for per-sample cutadapt diagnostics.",
      call. = FALSE
    )
  }
}

################################################################################
## PART 9: RETENTION ANALYSIS + AUTOMATIC FLAGS
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Retention analysis + anomaly flagging\n")
cat("─────────────────────────────────────────────────────────\n")

per_sample <- data.frame(
  sample     = sample_names,
  raw        = raw_counts    [basename(fnFs)],
  after_N    = filtN_counts  [basename(fnFs)],
  after_cut  = trimmed_counts[basename(fnFs)],
  stringsAsFactors = FALSE
)

# CRITICAL 3: cutadapt_failed is a first-class flag. For samples where
# cutadapt errored out (behavior b: skipped, no trimmed file produced), the
# lookup above returned NA for after_cut. Force to 0 so retention shows 0%
# (matching reality — zero reads survived) rather than NA propagating through
# the flag logic. Failed samples will naturally also trip extreme_drop /
# low_read_count / high_overall_loss because 0% retention crosses all three
# thresholds — that's correct behaviour, the sample really did drop out.
per_sample$flag_cutadapt_failed <- per_sample$sample %in% failed_samples
per_sample$after_cut[per_sample$flag_cutadapt_failed] <- 0

per_sample$pct_after_N   <- 100 * per_sample$after_N   / pmax(per_sample$raw,     1)
per_sample$pct_after_cut <- 100 * per_sample$after_cut / pmax(per_sample$after_N, 1)
per_sample$pct_overall   <- 100 * per_sample$after_cut / pmax(per_sample$raw,     1)

per_sample$flag_extreme_drop     <-
  per_sample$pct_after_N < 50 | per_sample$pct_after_cut < 50
per_sample$flag_low_read_count   <- per_sample$after_cut < 1000
per_sample$flag_high_overall_loss <- per_sample$pct_overall < 70

flagged_samples <- per_sample$sample[
  per_sample$flag_cutadapt_failed |
  per_sample$flag_extreme_drop |
  per_sample$flag_low_read_count |
  per_sample$flag_high_overall_loss
]

cat("  Samples processed successfully: ",
    nrow(per_sample) - length(failed_samples), "\n")
cat("  Samples failed (cutadapt):      ", length(failed_samples), "\n")
cat("  Samples flagged:                ", length(flagged_samples), "\n")
if (length(flagged_samples) > 0) {
  cat("  Flagged:                        ",
      paste(flagged_samples, collapse = ", "), "\n")
}
cat("\n")

################################################################################
## PART 10: WRITE qc_summary.json AND report.md
################################################################################

sample_records <- lapply(seq_len(nrow(per_sample)), function(i) {
  s <- per_sample[i, ]
  flags <- c()
  # CRITICAL 3: emit cutadapt_failed first so downstream consumers of
  # qc_summary.json see it before the other flags. A sample that failed
  # cutadapt will typically also carry extreme_drop / low_read_count /
  # high_overall_loss, but the cutadapt_failed flag is the actionable one.
  if (s$flag_cutadapt_failed)   flags <- c(flags, "cutadapt_failed")
  if (s$flag_extreme_drop)      flags <- c(flags, "extreme_drop")
  if (s$flag_low_read_count)    flags <- c(flags, "low_read_count")
  if (s$flag_high_overall_loss) flags <- c(flags, "high_overall_loss")
  list(
    sample                  = s$sample,
    raw_reads               = as.integer(s$raw),
    reads_after_N_filter    = as.integer(s$after_N),
    reads_after_primer_trim = as.integer(s$after_cut),
    pct_retained_N_filter   = round(s$pct_after_N,   2),
    pct_retained_cutadapt   = round(s$pct_after_cut, 2),
    pct_retained_overall    = round(s$pct_overall,   2),
    flags                   = if (length(flags) > 0) flags else list()
  )
})

end_time <- Sys.time()
runtime  <- as.numeric(difftime(end_time, start_time, units = "secs"))

summary_json <- list(
  skill         = "claw-amplicon-qc",
  version       = "0.1.0",
  timestamp_utc = format(end_time, "%Y-%m-%dT%H:%M:%SZ", tz = "UTC"),
  runtime_secs  = round(runtime, 1),
  inputs = list(
    raw_folder    = opt$raw,
    output_folder = opt$output,
    fwd_primer    = opt$fwd_primer,
    rev_primer    = opt$rev_primer,
    min_length    = opt$min_length
  ),
  totals = list(
    samples              = nrow(per_sample),
    samples_flagged      = length(flagged_samples),
    raw_reads_total      = as.integer(sum(per_sample$raw)),
    after_N_total        = as.integer(sum(per_sample$after_N)),
    after_cutadapt_total = as.integer(sum(per_sample$after_cut))
  ),
  samples = sample_records
)

summary_json_file <- file.path(opt$output, "qc_summary.json")
write_json(summary_json, summary_json_file, auto_unbox = TRUE, pretty = TRUE)

# Human-readable report.md
# CRITICAL 3: Summary now distinguishes successful from failed samples, and a
# 'Failed Samples' section always renders — showing 'None' on healthy runs, so
# readers get in the habit of looking at that section and don't miss failures
# the first time they occur.
n_successful <- nrow(per_sample) - length(failed_samples)

report_lines <- c(
  "# claw-amplicon-qc — QC Report",
  "",
  sprintf("**Run timestamp:** %s", summary_json$timestamp_utc),
  sprintf("**Runtime:** %.1f seconds", runtime),
  sprintf("**Skill version:** %s", summary_json$version),
  "",
  "## Inputs",
  "",
  sprintf("- Raw folder: `%s`", opt$raw),
  sprintf("- Output folder: `%s`", opt$output),
  sprintf("- Forward primer: `%s`", opt$fwd_primer),
  sprintf("- Reverse primer: `%s`", opt$rev_primer),
  sprintf("- Minimum length after trimming: %d bp", opt$min_length),
  "",
  "## Summary",
  "",
  sprintf("- Samples processed successfully: **%d**", n_successful),
  sprintf("- Samples failed (cutadapt): **%d**", length(failed_samples)),
  sprintf("- Samples flagged: **%d**", length(flagged_samples)),
  sprintf("- Total raw reads: %s", format(sum(per_sample$raw), big.mark = ",")),
  sprintf("- Total reads after N-filter: %s", format(sum(per_sample$after_N), big.mark = ",")),
  sprintf("- Total reads after primer trim: %s", format(sum(per_sample$after_cut), big.mark = ",")),
  "",
  "## Failed Samples",
  "",
  if (length(failed_samples) == 0) {
    "None. All samples processed successfully through cutadapt."
  } else {
    paste(c(
      sprintf("The following %d sample(s) failed at the cutadapt stage and are",
              length(failed_samples)),
      "excluded from downstream stats. See `04_cutadapt_log.txt` for the",
      "per-sample cutadapt output.",
      "",
      paste(sprintf("- `%s`", failed_samples), collapse = "\n")
    ), collapse = "\n")
  },
  "",
  "## Per-sample retention",
  "",
  "| Sample | Raw | After N | After Cut | %N | %Cut | %Overall | Flags |",
  "|--------|-----|---------|-----------|----|------|---------|-------|"
)

for (i in seq_len(nrow(per_sample))) {
  s <- per_sample[i, ]
  flags <- c()
  # CRITICAL 3: emit cutadapt_failed first (same ordering as JSON)
  if (s$flag_cutadapt_failed)   flags <- c(flags, "cutadapt_failed")
  if (s$flag_extreme_drop)      flags <- c(flags, "extreme_drop")
  if (s$flag_low_read_count)    flags <- c(flags, "low_read_count")
  if (s$flag_high_overall_loss) flags <- c(flags, "high_overall_loss")
  flag_str <- if (length(flags) > 0) paste(flags, collapse = ", ") else "—"
  report_lines <- c(report_lines, sprintf(
    "| %s | %s | %s | %s | %.1f | %.1f | %.1f | %s |",
    s$sample,
    format(s$raw,       big.mark = ","),
    format(s$after_N,   big.mark = ","),
    format(s$after_cut, big.mark = ","),
    s$pct_after_N, s$pct_after_cut, s$pct_overall,
    flag_str
  ))
}

report_lines <- c(report_lines,
  "",
  "## Files produced",
  "",
  "- `01_raw_stats.txt` — seqkit stats on raw reads",
  "- `filtN/` — reads after N removal",
  "- `02_filtN_stats.txt` — seqkit stats after N removal",
  "- `cutadapt_trimmed/` — reads after primer trimming",
  "- `04_cutadapt_log.txt` — Cutadapt per-sample log",
  "- `03_trimmed_stats.txt` — seqkit stats after primer trimming",
  "- `qc_summary.json` — machine-readable summary",
  "- `report.md` — this file",
  "",
  "## Next step",
  "",
  "This skill deliberately stops before quality filtering and DADA2 denoising.",
  "The forthcoming `claw-amplicon-dada2` skill will handle quality profile review,",
  "truncation-length recommendation with a human approval gate, and denoising.",
  ""
)

report_file <- file.path(opt$output, "report.md")
writeLines(report_lines, report_file)

################################################################################
## PART 11: FINAL SUMMARY TO STDOUT
################################################################################

cat("═════════════════════════════════════════════════════════\n")
cat("  claw-amplicon-qc — complete\n")
cat("═════════════════════════════════════════════════════════\n")
cat("  Samples processed: ", nrow(per_sample), "\n")
cat("  Samples flagged:   ", length(flagged_samples), "\n")
cat("  Runtime:           ", sprintf("%.1f s", runtime), "\n")
cat("  Output folder:     ", opt$output, "\n\n")
cat("  Machine-readable:  ", summary_json_file, "\n")
cat("  Human-readable:    ", report_file, "\n")
cat("═════════════════════════════════════════════════════════\n")
