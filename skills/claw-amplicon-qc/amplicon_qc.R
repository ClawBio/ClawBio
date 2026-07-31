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
  chk <- suppressWarnings(system2(tool, args = "--version",
                                  stdout = TRUE, stderr = TRUE))
  if (inherits(chk, "try-error") || length(chk) == 0) {
    stop(tool, " not found on PATH. Install with: conda install -c bioconda ",
         tool, call. = FALSE)
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

# Try Pattern A first (Illumina default: SAMPLE_R1_..., SAMPLE_R2_...)
fnFs_A <- grep("_R1[_.]", basename(all_fastq), value = FALSE)
fnRs_A <- grep("_R2[_.]", basename(all_fastq), value = FALSE)

# Then Pattern B (prefix: R1_SAMPLE, R2_SAMPLE)
fnFs_B <- grep("^R1_",    basename(all_fastq), value = FALSE)
fnRs_B <- grep("^R2_",    basename(all_fastq), value = FALSE)

if (length(fnFs_A) > 0 && length(fnRs_A) > 0) {
  fnFs <- sort(all_fastq[fnFs_A])
  fnRs <- sort(all_fastq[fnRs_A])
  naming_pattern <- "Illumina default (SAMPLE_R1_..., SAMPLE_R2_...)"
  # Sample name = strip everything from _R1 onwards
  sample_names <- gsub("_R1[_.].*", "", basename(fnFs))
} else if (length(fnFs_B) > 0 && length(fnRs_B) > 0) {
  fnFs <- sort(all_fastq[fnFs_B])
  fnRs <- sort(all_fastq[fnRs_B])
  naming_pattern <- "Prefix (R1_SAMPLE, R2_SAMPLE)"
  # Sample name = strip R1_ prefix and file extension
  sample_names <- gsub("^R1_", "", basename(fnFs))
  sample_names <- gsub("\\.(fastq|fq)(\\.gz)?$", "", sample_names, ignore.case = TRUE)
} else {
  stop("No paired-end FASTQ files found in ", opt$raw,
       "\nExpected either:",
       "\n  Pattern A: SAMPLE_R1_...fastq.gz + SAMPLE_R2_...fastq.gz",
       "\n  Pattern B: R1_SAMPLE.fastq.gz + R2_SAMPLE.fastq.gz",
       call. = FALSE)
}

if (length(fnFs) != length(fnRs)) {
  stop("Unequal number of R1 and R2 files: ",
       length(fnFs), " R1 vs ", length(fnRs), " R2", call. = FALSE)
}

cat("Detected naming pattern:", naming_pattern, "\n")
cat("Discovered", length(fnFs), "sample pairs.\n\n")

################################################################################
## PART 3: HELPER — run seqkit stats and parse the output
################################################################################

run_seqkit_stats <- function(files, output_file) {
  system2("seqkit", args = c("stats", "-a", "-T", files),
          stdout = output_file, stderr = "")
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

for (i in seq_along(fnFs)) {
  cat(sprintf("  [%d/%d] %s\n", i, length(fnFs), sample_names[i]))
  args_vec <- c(
    "-g", FWD, "-a", REV.RC,
    "-G", REV, "-A", FWD.RC,
    "--minimum-length", as.character(opt$min_length),
    "--nextseq-trim=20",
    "-n", "2",
    "-j", "0",
    "-o", fnFs.cut[i],
    "-p", fnRs.cut[i],
    fnFs.filtN[i], fnRs.filtN[i]
  )
  cat("\n\n===== Sample:", sample_names[i], "=====\n",
      file = cutadapt_log, append = TRUE)
  system2("cutadapt", args = args_vec,
          stdout = cutadapt_log, stderr = cutadapt_log, wait = TRUE)
}
cat("\n")

################################################################################
## PART 8: STAGE 5 — SEQKIT AFTER CUTADAPT
################################################################################

cat("─────────────────────────────────────────────────────────\n")
cat("  Stage 5/5 — seqkit stats on primer-trimmed reads\n")
cat("─────────────────────────────────────────────────────────\n")

trimmed_stats_file <- file.path(opt$output, "03_trimmed_stats.txt")
trimmed_stats      <- run_seqkit_stats(c(fnFs.cut, fnRs.cut), trimmed_stats_file)
trimmed_counts     <- seqkit_reads_by_file(trimmed_stats)

cat("  Total reads after primer trimming:",
    format(sum(trimmed_counts), big.mark = ","), "\n")
cat("  Stats saved to: ", trimmed_stats_file, "\n\n")

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

per_sample$pct_after_N   <- 100 * per_sample$after_N   / pmax(per_sample$raw,     1)
per_sample$pct_after_cut <- 100 * per_sample$after_cut / pmax(per_sample$after_N, 1)
per_sample$pct_overall   <- 100 * per_sample$after_cut / pmax(per_sample$raw,     1)

per_sample$flag_extreme_drop     <-
  per_sample$pct_after_N < 50 | per_sample$pct_after_cut < 50
per_sample$flag_low_read_count   <- per_sample$after_cut < 1000
per_sample$flag_high_overall_loss <- per_sample$pct_overall < 70

flagged_samples <- per_sample$sample[
  per_sample$flag_extreme_drop |
  per_sample$flag_low_read_count |
  per_sample$flag_high_overall_loss
]

cat("  Samples processed: ", nrow(per_sample), "\n")
cat("  Samples flagged:   ", length(flagged_samples), "\n")
if (length(flagged_samples) > 0) {
  cat("  Flagged:           ", paste(flagged_samples, collapse = ", "), "\n")
}
cat("\n")

################################################################################
## PART 10: WRITE qc_summary.json AND report.md
################################################################################

sample_records <- lapply(seq_len(nrow(per_sample)), function(i) {
  s <- per_sample[i, ]
  flags <- c()
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
  sprintf("- Samples processed: **%d**", nrow(per_sample)),
  sprintf("- Samples flagged: **%d**", length(flagged_samples)),
  sprintf("- Total raw reads: %s", format(sum(per_sample$raw), big.mark = ",")),
  sprintf("- Total reads after N-filter: %s", format(sum(per_sample$after_N), big.mark = ",")),
  sprintf("- Total reads after primer trim: %s", format(sum(per_sample$after_cut), big.mark = ",")),
  "",
  "## Per-sample retention",
  "",
  "| Sample | Raw | After N | After Cut | %N | %Cut | %Overall | Flags |",
  "|--------|-----|---------|-----------|----|------|---------|-------|"
)

for (i in seq_len(nrow(per_sample))) {
  s <- per_sample[i, ]
  flags <- c()
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
