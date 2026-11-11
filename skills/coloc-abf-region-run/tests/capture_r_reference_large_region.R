# Capture R coloc.abf reference numbers on regions large enough for coloc's
# prior adjustment to fire, for tests/fixtures/coloc_abf_large_region_r_reference.json.
#
# Since CRAN coloc 5.2.2, coloc.abf lowers p1, p2 and p12 to 1/(n + 1) whenever
# n * p >= 1, each against its OWN count: p1 against dataset 1's rows, p2 against
# dataset 2's rows, p12 against the shared rows. The 12-variant fixture
# (capture_r_reference.R) is far too small to trigger it. These three
# configurations do:
#
#   A  12,000 shared variants, default priors   -> p1 and p2 lowered
#   B  10,500 vs 9,999 variants, 9,999 shared    -> p1 only. p2 sits one variant
#                                                   under its threshold, and an
#                                                   implementation that checked p1
#                                                   against the SHARED count would
#                                                   not lower it
#   C  configuration A with p12 = 1e-4, the top  -> all three lowered
#      of the sensitivity sweep
#
# The data are a closed-form signal, not random draws, so the Python test
# rebuilds them exactly without bundling 12,000-row inputs:
#   z_i = 0.9 sin(a i + ph) + 0.6 cos(b i), with the lead positions overwritten;
#   beta = 0.01 z, se = 0.01, sdY = 1, type = "quant", snp = "1_<1000000 + i>_A_G".
# Trait 1 leads at i = 6000 (z 4.5); trait 2 at i = 3000 (z 4.5) and 6000 (z 3.5),
# so H3 and H4 both carry mass and a prior change moves all five posteriors.
#
# As with capture_r_reference.R, the GPL R package is never installed, imported
# or vendored into this repository, and this script is NOT run by the test
# suite. Run it in an environment of your own; it prints the "r_reference"
# block of the fixture as JSON.
#
#   Rscript capture_r_reference_large_region.R   # needs install.packages("coloc")
#
# Captured with R 4.6.1 and coloc 5.2.3.

library(coloc)

gen <- function(n, start, leads, a, b, ph) {
  i <- start:(start + n - 1)
  z <- 0.9 * sin(a * i + ph) + 0.6 * cos(b * i)
  for (L in leads) z[i == L[1]] <- L[2]
  list(beta = 0.01 * z, varbeta = rep(0.01, n)^2,
       snp = paste0("1_", 1000000 + i, "_A_G"), type = "quant", sdY = 1.0)
}
L1 <- list(c(6000, 4.5))
L2 <- list(c(3000, 4.5), c(6000, 3.5))

configs <- list(
  A = list(d1 = gen(12000, 0, L1, 0.37, 1.13, 0.0), d2 = gen(12000, 0, L2, 0.53, 0.71, 1.0), p12 = 1e-5),
  B = list(d1 = gen(10500, 0, L1, 0.37, 1.13, 0.0), d2 = gen(9999, 501, L2, 0.53, 0.71, 1.0), p12 = 1e-5),
  C = list(d1 = gen(12000, 0, L1, 0.37, 1.13, 0.0), d2 = gen(12000, 0, L2, 0.53, 0.71, 1.0), p12 = 1e-4)
)

out <- c()
for (nm in names(configs)) {
  cf <- configs[[nm]]
  res <- withCallingHandlers(
    coloc.abf(cf$d1, cf$d2, p1 = 1e-4, p2 = 1e-4, p12 = cf$p12),
    warning = function(w) { message(nm, ": ", conditionMessage(w)); invokeRestart("muffleWarning") }
  )
  s <- res$summary; p <- res$priors
  out <- c(out, sprintf(
    '"%s": {"nsnps": %d, "PP.H0": %.17g, "PP.H1": %.17g, "PP.H2": %.17g, "PP.H3": %.17g, "PP.H4": %.17g, "p1": %.17g, "p2": %.17g, "p12": %.17g}',
    nm, as.integer(s["nsnps"]), s["PP.H0.abf"], s["PP.H1.abf"], s["PP.H2.abf"],
    s["PP.H3.abf"], s["PP.H4.abf"], p["p1"], p["p2"], p["p12"]))
}
cat("{", paste(out, collapse = ",\n "), "}\n", sep = "")
